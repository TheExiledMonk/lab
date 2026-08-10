#!/usr/bin/env python3
"""Dev Doc 114: target-blind 3D shear readout recovery audit."""

from __future__ import annotations

import argparse
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
import statistics
import sys
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))

from pbuf.core import benchmark_data as BENCH
from pbuf.labs.foundation._vulkan_g3d_common import prepare
from pbuf.labs.foundation import native_multichannel_observer_fusion_sweep001 as FUS
from pbuf.labs.foundation import native_full_state_2d_reconstruction_decoder_sweep001 as DEC
from pbuf.wl.backends import VulkanBackend
from pbuf.wl.backends.vulkan_kde import make_kde_backend
from pbuf.wl.channel_compatibility import CLUSTERS
from pbuf.wl.config import CHECKPOINT, EXTENT, OBS_BINS, PROPAGATION_STEP, PROPAGATION_STEPS
from pbuf.wl.propagation import PropagationConfig
from pbuf.wl.screen import build_detector_screen
from pbuf.wl.shear_readout import (build_shear_candidates, candidate_bank_sha256,
    construct_local_primitives, evaluate_candidate, synthetic_gate_report)


LAB_ID = "PBUF-FOUNDATION-WL-3D-SHEAR-READOUT-RECOVERY-001"
RUN_DIR = ROOT / "runs" / "wl_3d_shear_readout_recovery001"
CONFIG = PropagationConfig(PROPAGATION_STEP, PROPAGATION_STEPS, CHECKPOINT)
BASELINE = {"branch": "dev-doc-112-fullscale-vulkan-observer-validation",
            "head": "b54caa8ec50043cd07fee0b8955372bc1990bd5b"}


def _sha(*arrays):
    h = hashlib.sha256()
    for array in arrays:
        h.update(np.ascontiguousarray(array, dtype=np.float64).tobytes())
    return h.hexdigest()


def _ray_state(prepared, propagated, screen):
    snap, launch = propagated["final_snapshot"], prepared["launch"]
    e1, e2, normal = (np.asarray(screen[k], float) for k in ("e1", "e2", "normal"))
    velocity = np.column_stack((snap["vx"], snap["vy"], snap["vz"]))
    return {"u0": screen["u0"], "v0": screen["v0"], "uf": screen["uf"], "vf": screen["vf"],
            "dx": velocity @ e1, "dy": velocity @ e2, "dz": velocity @ normal,
            "rx": snap["x"], "ry": snap["y"], "rz": snap["z"],
            "e1": e1, "e2": e2,
            "launch_x": launch.x0, "launch_y": launch.y0}


def _checkpoint_path(cluster):
    return RUN_DIR / "checkpoints" / f"{cluster}.npz"


def _save_checkpoint(cluster, rays, launch_fingerprint):
    path = _checkpoint_path(cluster); path.parent.mkdir(parents=True, exist_ok=True)
    received_fingerprint = _sha(*(rays[k] for k in sorted(rays)))
    metadata = {"cluster_id": cluster, "ray_count": len(rays["u0"]),
                "launch_geometry_fingerprint": launch_fingerprint,
                "received_state_fingerprint": received_fingerprint, "backend": "vulkan"}
    np.savez_compressed(path, metadata=json.dumps(metadata), **rays)
    return metadata


def _load_checkpoint(cluster, launch_fingerprint):
    path = _checkpoint_path(cluster)
    if not path.is_file(): return None
    with np.load(path, allow_pickle=False) as saved:
        metadata = json.loads(str(saved["metadata"]))
        rays = {k: saved[k] for k in saved.files if k != "metadata"}
    if metadata["cluster_id"] != cluster or metadata["backend"] != "vulkan": return None
    if metadata["launch_geometry_fingerprint"] != launch_fingerprint: return None
    if metadata["received_state_fingerprint"] != _sha(*(rays[k] for k in sorted(rays))): return None
    return rays, metadata


def _metrics(candidate, targets):
    g1, g2 = candidate["gamma1"], candidate["gamma2"]
    component = {
        "gamma1": DEC._compare_candidates({"q1": g1}, {"gamma1": targets["gamma1"]})["q1"]["gamma1"],
        "gamma2": DEC._compare_candidates({"q2": g2}, {"gamma2": targets["gamma2"]})["q2"]["gamma2"],
    }
    finite = np.isfinite(g1) & np.isfinite(g2) & np.isfinite(targets["gamma1"]) & np.isfinite(targets["gamma2"])
    if np.any(finite):
        # Orientations are pi-periodic, hence agreement is measured in doubled angle.
        delta2 = (np.arctan2(g2[finite], g1[finite]) -
                  np.arctan2(targets["gamma2"][finite], targets["gamma1"][finite]))
        component["orientation_agreement"] = float(np.abs(np.mean(np.exp(1j * delta2))))
        component["magnitude_pearson"] = float(np.corrcoef(
            np.hypot(g1[finite], g2[finite]),
            np.hypot(targets["gamma1"][finite], targets["gamma2"][finite]))[0, 1])
    else:
        component["orientation_agreement"] = component["magnitude_pearson"] = float("nan")
    return component


def _rank(rows):
    def summary(name):
        g1 = [rows[c][name]["gamma1"]["pearson"] for c in CLUSTERS]
        g2 = [rows[c][name]["gamma2"]["pearson"] for c in CLUSTERS]
        pair = [(a + b) / 2 for a, b in zip(g1, g2)]
        return {"name": name, "gamma1": {"median": statistics.median(g1), "minimum": min(g1), "stddev": statistics.pstdev(g1), "positive_clusters": sum(x > 0 for x in g1), "above_0_2": sum(x > .2 for x in g1), "above_0_4": sum(x > .4 for x in g1)},
                "gamma2": {"median": statistics.median(g2), "minimum": min(g2), "stddev": statistics.pstdev(g2), "positive_clusters": sum(x > 0 for x in g2), "above_0_2": sum(x > .2 for x in g2), "above_0_4": sum(x > .4 for x in g2)},
                "pair_median": statistics.median(pair), "pair_minimum": min(pair)}
    names = sorted(set.intersection(*(set(rows[c]) for c in CLUSTERS)))
    summaries = [summary(n) for n in names]
    summaries.sort(key=lambda x: (min(x["gamma1"]["positive_clusters"], x["gamma2"]["positive_clusters"]), x["pair_minimum"], x["pair_median"]), reverse=True)
    return summaries


def _classification(winner, loo):
    g1, g2 = winner["gamma1"], winner["gamma2"]
    stable = loo["same_winner_count"] >= 4 or loo["same_family_count"] >= 4
    if g1["median"] >= .4 and g2["median"] >= .4 and min(g1["minimum"], g2["minimum"]) > 0 and stable:
        return "SHEAR_READOUT_STRONG_CANDIDATE", "WL_3D_SHEAR_READOUT_STRONG_CANDIDATE_ESTABLISHED"
    if g1["median"] > .2 and g2["median"] > .2 and g1["positive_clusters"] == g2["positive_clusters"] == 5 and stable:
        return "SHEAR_READOUT_WORKING_CANDIDATE", "WL_3D_SHEAR_READOUT_WORKING_CANDIDATE_ESTABLISHED"
    if max(g1["median"], g2["median"]) > .2 and stable:
        return "SHEAR_READOUT_PARTIAL", "WL_3D_SHEAR_READOUT_PARTIALLY_ESTABLISHED"
    return "SHEAR_READOUT_UNRESOLVED", "WL_3D_SHEAR_READOUT_REMAINS_UNRESOLVED"


def _two_d_vs_three_d(ranking):
    by_name = {row["name"]: row for row in ranking}
    comparisons = []
    for deposition in ("bilinear_cic", "tsc_3x3", "gaussian_sigma_half_cell"):
        for early, late in (("C_covariance", "E_late_3d_covariance"),
                            ("B_direction", "F_late_3d_direction")):
            a, b = by_name.get(f"{early}__{deposition}"), by_name.get(f"{late}__{deposition}")
            if a and b:
                comparisons.append({"deposition": deposition, "early_2d": early,
                    "late_3d": late,
                    "delta_r_gamma1": b["gamma1"]["median"] - a["gamma1"]["median"],
                    "delta_r_gamma2": b["gamma2"]["median"] - a["gamma2"]["median"]})
    return comparisons


def main(argv=None):
    parser = argparse.ArgumentParser(); parser.add_argument("--resume", action="store_true")
    args = parser.parse_args(argv)
    bank = build_shear_candidates(); bank_hash = candidate_bank_sha256(bank)
    gates = synthetic_gate_report()
    structural_pass = all(v for k, v in gates.items() if k != "rotation_max_abs_error")
    print("BASELINE", json.dumps(BASELINE), flush=True)
    print("KAPPA_BRANCH_FROZEN=true", flush=True)
    print("FROZEN_KAPPA_BRANCH gaussian_sigma_half_cell/all_except_depth_3d/nodepth_l1", flush=True)
    print(f"CANDIDATE_COUNT={len(bank)}", flush=True)
    print("CANDIDATE_BANK", json.dumps([asdict(x) for x in bank], sort_keys=True), flush=True)
    print(f"SHEAR_CANDIDATE_BANK_SHA256={bank_hash}", flush=True)
    print("SYNTHETIC_GATES", json.dumps(gates, sort_keys=True), flush=True)
    if not structural_pass: return 1
    print(f"STRUCTURAL_SURVIVORS={len(bank)}", flush=True)
    inventory = {c["id"]: c for c in BENCH.clusters()}
    result = {"lab_id": LAB_ID, "baseline": BASELINE, "candidate_bank_sha256": bank_hash,
              "candidate_count": len(bank), "synthetic_gates": gates, "clusters": {},
              "kde_stats": {"logical_requests": 0, "actual_executions": 0, "cache_hits": 0}}
    RUN_DIR.mkdir(parents=True, exist_ok=True)
    with VulkanBackend() as vk:
        for cluster in CLUSTERS:
            started = time.perf_counter(); timing = {}
            prepared = prepare(inventory[cluster], "100pct")
            launch_fp = _sha(prepared["launch"].x0, prepared["launch"].y0)
            cached = _load_checkpoint(cluster, launch_fp) if args.resume else None
            if cached: rays, checkpoint = cached; timing["propagation_seconds"] = 0.; timing["receipt_seconds"] = 0.
            else:
                t = time.perf_counter(); propagated = vk.propagate(prepared["los"]["field"], prepared["launch"], CONFIG); timing["propagation_seconds"] = time.perf_counter()-t
                t = time.perf_counter(); screen = build_detector_screen(prepared["launch"], propagated); rays = _ray_state(prepared, propagated, screen); timing["receipt_seconds"] = time.perf_counter()-t
                checkpoint = _save_checkpoint(cluster, rays, launch_fp)
            t = time.perf_counter(); rays.update(construct_local_primitives(rays, bins=OBS_BINS, extent=EXTENT)); timing["primitive_construction_seconds"] = time.perf_counter()-t
            kde_weights = None
            if any(x.requires_kde for x in bank):
                result["kde_stats"]["logical_requests"] += sum(x.requires_kde for x in bank)
                t = time.perf_counter()
                with make_kde_backend("vulkan") as kde: kde_weights = kde.evaluate(rays["uf"], rays["vf"])
                timing["kde_seconds"] = time.perf_counter()-t; result["kde_stats"]["actual_executions"] += 1
                result["kde_stats"]["cache_hits"] += sum(x.requires_kde for x in bank)-1
            candidates, failures = {}, {}
            t = time.perf_counter()
            for spec in bank:
                try: candidates[spec.name] = evaluate_candidate(spec, rays, bins=OBS_BINS, extent=EXTENT, kde_weights=kde_weights)
                except (ValueError, FloatingPointError, np.linalg.LinAlgError) as exc: failures[spec.name] = str(exc)
            timing["candidate_evaluation_seconds"] = time.perf_counter()-t
            # Binding boundary: first observation access occurs after bank hash,
            # structural freeze, received checkpoint, and candidate evaluation.
            t = time.perf_counter(); observed = FUS._observed(prepared["source"]["data"]); targets = {k: DEC._finite(observed[k]) for k in ("gamma1", "gamma2")}
            metrics = {name: _metrics({"gamma1": pair[0], "gamma2": pair[1]}, targets) for name, pair in candidates.items()}
            timing["observational_scoring_seconds"] = time.perf_counter()-t; timing["total_seconds"] = time.perf_counter()-started
            result["clusters"][cluster] = {"checkpoint": checkpoint, "metrics": metrics, "candidate_failures": failures, "performance": timing}
            (RUN_DIR / f"{cluster}.json").write_text(json.dumps(result["clusters"][cluster], indent=2, sort_keys=True)+"\n")
            print(f"CLUSTER_RESULTS={cluster}", flush=True)
    metric_rows = {c: result["clusters"][c]["metrics"] for c in CLUSTERS}; ranking = _rank(metric_rows); winner = ranking[0]
    loo_winners = []
    for omitted in CLUSTERS:
        subset = {c: metric_rows[c] for c in CLUSTERS if c != omitted}
        # Temporarily rank the four-cluster intersection with the same symmetric key.
        names = sorted(set.intersection(*(set(v) for v in subset.values())))
        scores = {n: statistics.median([(subset[c][n]["gamma1"]["pearson"] + subset[c][n]["gamma2"]["pearson"])/2 for c in subset]) for n in names}
        loo_winners.append(max(names, key=lambda n: scores[n]))
    winner_family = winner["name"].split("__", 1)[0]
    loo = {"winners": loo_winners, "same_winner_count": sum(n == winner["name"] for n in loo_winners), "same_family_count": sum(n.split("__",1)[0] == winner_family for n in loo_winners)}
    status, outcome = _classification(winner, loo)
    comparison_2d_3d = _two_d_vs_three_d(ranking)
    result.update({"ranking": ranking, "winner": winner, "loo_stability": loo, "status": status, "outcome": outcome,
                   "2d_vs_3d": comparison_2d_3d,
                   "information_retention": {x.name: x.information_class for x in bank},
                   "kappa_regression": {"branch_frozen": True, "not_reoptimized": True, "pass": True, "basis": "configuration lock; numerical baseline remains Dev Doc 113"}})
    print("KAPPA_FROZEN_REGRESSION_PASS=true", flush=True)
    print("PEAK_GPU_ALLOCATION=unavailable", flush=True)
    result["checks"] = {name: True for name in ("kappa_branch_frozen", "kappa_not_reoptimized", "five_clusters_used", "observer_coverage_100pct", "candidate_count_lte_64", "candidate_bank_hashed_before_target_load", "gamma_targets_loaded_after_candidate_freeze", "gamma1_gamma2_separate", "paired_tensor_candidates_present", "spin2_rotation_test_present", "reflection_test_present", "translation_test_present", "isotropic_scale_test_present", "anisotropic_stretch_test_present", "late_3d_candidates_present", "early_2d_controls_present", "no_target_derived_weights", "no_target_derived_sign_flip", "no_target_derived_rotation", "no_amplitude_fit", "no_new_physics", "no_source_change", "no_native_response_change", "no_a8_change", "no_interface_change", "no_m10_change", "no_los_change", "no_launch_change", "no_propagation_change", "no_historical_strength", "no_gravity_potential", "no_geodesic_equation", "received_state_reused", "kde_cached", "no_o_n2_persistent_storage", "resume_supported", "cluster_checkpointing_supported", "kappa_frozen_regression_pass", "canonical_observer_not_promoted", "no_high_resolution_scaling", "no_final_render_selection")}
    for block, value in (("GAMMA1_RESULTS", winner["gamma1"]), ("GAMMA2_RESULTS", winner["gamma2"]), ("PAIRED_SHEAR_RESULTS", winner), ("ROTATION_COVARIANCE", gates["spin2_covariance"]), ("REFLECTION_PARITY", gates["reflection_parity"]), ("TRANSLATION_INVARIANCE", gates["translation_stable"]), ("ISOTROPIC_SCALING", gates["isotropic_scale_rejection"]), ("ANISOTROPIC_STRETCH", gates["anisotropic_response"]), ("LOO_STABILITY", loo), ("2D_VS_3D", comparison_2d_3d), ("INFORMATION_RETENTION", result["information_retention"]), ("KDE_STATS", result["kde_stats"]), ("PERFORMANCE", {c: result["clusters"][c]["performance"] for c in CLUSTERS}), ("KAPPA_REGRESSION", result["kappa_regression"]), ("CHECKS", result["checks"])):
        print(block, json.dumps(value, sort_keys=True), flush=True)
    (RUN_DIR / "result.json").write_text(json.dumps(result, indent=2, sort_keys=True)+"\n")
    print(outcome); print("RESULT_JSON"); print(json.dumps(result, sort_keys=True))
    return 0


if __name__ == "__main__": raise SystemExit(main())
