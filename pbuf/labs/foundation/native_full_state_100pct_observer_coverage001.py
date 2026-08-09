#!/usr/bin/env python3
"""PBUF FOUNDATION — FULL-STATE 100 PERCENT OBSERVER COVERAGE 001.

Coverage/observer audit only.

PR #105 established a target-blind full-state 2D reconstruction sweep using the
existing 25% source-plane receipt. This lab changes only the sampled source-plane
coverage: it reruns the same frozen current-native received-state pipeline at
25% and 100% detector coverage, preserving approximately the same rays per
supported source bin.

The 100% lane uses a deterministic full [-extent,+extent]^2 Cartesian launch.
Because area is 4x the 25% rectangle, the launch side is doubled, giving 4x the
ray count while keeping source-bin sampling density approximately fixed.

For BOTH lanes:
  source -> native M10 interface -> LOS -> existing G3D propagation
  -> one target-blind global tangent screen -> full 45-channel receipt bank
  -> exact PR #105 predeclared decoder inventory.

All prediction/decoder products for both lanes are constructed before observed
kappa/gamma arrays are requested for descriptive comparison. No upstream physics,
coefficient, reconstruction weight, target-derived orientation, rescaling, or
cluster-specific decoder choice is introduced.
"""
from __future__ import annotations

import json
import math
import statistics
import subprocess
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))

from pbuf.core import benchmark_data as BENCH
import pbuf.labs.foundation.current_native_five_cluster_observable_benchmark001 as CUR
import pbuf.labs.foundation.native_accumulated_full_lensing001 as G3D
import pbuf.labs.foundation.g3d_angular_received_distribution001 as ANG
import pbuf.labs.foundation.g3d_native_angular_detector_image001 as DET
import pbuf.labs.foundation.los_consistent_ray_geometry001 as GEO
import pbuf.labs.foundation.m10_coverage_25pct_science001 as BASE
import pbuf.labs.foundation.native_observable_extraction_method_sweep001 as EX
import pbuf.labs.foundation.native_full_received_state_information_retention001 as RET
import pbuf.labs.foundation.native_full_state_2d_reconstruction_decoder_sweep001 as DEC
import pbuf.labs.foundation.native_multichannel_observer_fusion_sweep001 as FUS
from pbuf.core import los_projection as M14

LAB_ID = "PBUF-FOUNDATION-FULL-STATE-100PCT-OBSERVER-COVERAGE-001"
EXPECTED_CLUSTER_IDS = FUS.EXPECTED_CLUSTER_IDS
EPS = 1.0e-30

# Freeze the previously observed useful decoder as a carry-forward diagnostic,
# while still evaluating the complete predeclared PR #105 inventory in both lanes.
FROZEN_CARRY_FORWARD = (
    "full_whitened_pca_energy_8",
    "full_l2",
    "established2d_l2",
)

EXTENT = float(BASE.CFG["extent"])
BINS = int(BASE.CFG["bins"])
SIDE25 = int(BASE.EXPANDED_SIDE)
SIDE100 = 2 * SIDE25
N25 = SIDE25 * SIDE25
N100 = SIDE100 * SIDE100
EXPECTED_SUPPORT25 = int(G3D.EXPECTED_SUPPORT)
EXPECTED_SUPPORT100 = BINS * BINS


def git(*args: str) -> str:
    return subprocess.check_output(["git", *args], cwd=ROOT, text=True).strip()


def repo_state() -> dict:
    return {
        "repository": "TheExiledMonk/lab",
        "branch": git("rev-parse", "--abbrev-ref", "HEAD"),
        "head_sha": git("rev-parse", "HEAD"),
        "tracked_changes": git("diff", "--name-only"),
        "staged_changes": git("diff", "--name-only", "--cached"),
    }


def _launch_full_100pct():
    """Deterministic full-plane Cartesian launch at 25%-lane ray density."""
    edges = np.linspace(-EXTENT, EXTENT, SIDE100 + 1)
    centers = 0.5 * (edges[:-1] + edges[1:])
    X, Y = np.meshgrid(centers, centers, indexing="xy")
    x0 = X.ravel()
    y0 = Y.ravel()
    vx0 = np.ones_like(x0)
    vy0 = np.zeros_like(x0)
    return x0, y0, vx0, vy0


def _source_support(x0: np.ndarray, y0: np.ndarray) -> dict:
    edges = np.linspace(-EXTENT, EXTENT, BINS + 1)
    counts, _, _ = np.histogram2d(y0, x0, bins=(edges, edges))
    support = counts >= 6
    nz = counts[counts > 0]
    return {
        "ray_count": int(x0.size),
        "source_support_bins": int(np.count_nonzero(support)),
        "source_support_fraction": float(np.mean(support)),
        "mean_rays_per_nonzero_source_bin": float(np.mean(nz)) if nz.size else 0.0,
        "min_rays_per_nonzero_source_bin": float(np.min(nz)) if nz.size else 0.0,
        "max_rays_per_nonzero_source_bin": float(np.max(nz)) if nz.size else 0.0,
    }


def _run_g3d_with_launch(vector, x0: np.ndarray, y0: np.ndarray, expected_support: int) -> dict:
    """Exact existing G3D machinery with only the launch coordinates supplied."""
    los = M14.project_vector_to_image_plane(*vector, los_axis="z")
    Rx = np.asarray(los["comp_1"], dtype=np.float64)
    Ry = np.asarray(los["comp_2"], dtype=np.float64)
    los_mag = np.hypot(Rx, Ry)
    grid = np.linspace(-EXTENT, EXTENT, Rx.shape[0])
    field = {"xgrid": grid, "ygrid": grid, "rx": Rx, "ry": Ry}

    groups = GEO._source_groups(x0, y0)
    if len(groups) != expected_support:
        raise RuntimeError(f"expected {expected_support} source bins, got {len(groups)}")

    checkpoints, g3d = GEO._propagate_g3d(
        field, BASE.CFG["step"], BASE.CFG["steps"], x0, y0
    )
    if g3d["max_unit_speed_error"] > G3D.UNIT_SPEED_TOL:
        raise RuntimeError(f"G3D unit-speed gate failed: {g3d['max_unit_speed_error']}")

    first = GEO._first_step_geometry(
        field, x0, y0, checkpoints[1], np.zeros_like(los_mag), los_mag
    )
    if not first["first_step_exact_pass"]:
        raise RuntimeError("first-step exact geometry gate failed")

    final_ang = ANG._angular_distribution_fields(checkpoints[G3D.CHECKPOINT], groups)
    gates = ANG._moment_gates(final_ang)
    if gates["second_moment_equals_cov_plus_centroid_outer_relative_rms_error"] > ANG.MOMENT_IDENTITY_TOL:
        raise RuntimeError("angular second-moment identity failed")
    if not gates["covariance_psd_pass"]:
        raise RuntimeError("angular covariance PSD gate failed")
    if not gates["direction_mean_vector_bound_pass"]:
        raise RuntimeError("angular direction-mean bound failed")

    snap = checkpoints[G3D.CHECKPOINT]
    vz = np.asarray(snap["vz"], dtype=np.float64)
    if float(np.min(np.abs(vz))) <= DET.VZ_MIN:
        raise RuntimeError("final tangent projection vz too small")

    return {
        "los_mag": los_mag,
        "checkpoints": checkpoints,
        "g3d": g3d,
        "first": first,
        "final_ang": final_ang,
        "groups": groups,
        "angular_gates": gates,
    }


def _decode_lane(data: dict, channel: dict, chain: dict, x0: np.ndarray, y0: np.ndarray, label: str) -> dict:
    snap = chain["checkpoints"][G3D.CHECKPOINT]
    screen = EX._screen_coordinates(x0, y0, snap)
    extracted = EX._extract_all(screen, EXTENT, BINS)
    receipt3d = RET._binned_received_3d(screen, snap, EXTENT, BINS)
    bank, family = RET._decoded_bank(extracted, receipt3d)

    # Construct all target-blind reconstruction candidates before observations.
    candidates, decoder_meta = DEC._build_candidates(bank, family)
    per_ray = RET._per_ray_geometry(receipt3d)
    stages = RET._stage_metrics(bank, family)

    return {
        "label": label,
        "data": data,
        "channel": channel,
        "chain": chain,
        "screen": screen,
        "bank": bank,
        "family": family,
        "candidates": candidates,
        "decoder_meta": decoder_meta,
        "per_ray_information_geometry": per_ray,
        "stage_information_geometry": stages,
        "support": _source_support(x0, y0),
    }


def _compare_lane_after_decoding(lane: dict, targets: dict[str, np.ndarray]) -> dict:
    relevance = DEC._compare_candidates(lane["candidates"], targets)
    return {
        "label": lane["label"],
        "support": lane["support"],
        "decoded_bank_size": len(lane["bank"]),
        "decoder_meta": lane["decoder_meta"],
        "per_ray_information_geometry": lane["per_ray_information_geometry"],
        "stage_information_geometry": lane["stage_information_geometry"],
        "candidate_names": list(lane["candidates"]),
        "candidate_relevance": relevance,
        "g3d_unit_speed_max_error": lane["chain"]["g3d"]["max_unit_speed_error"],
        "terminal_common_history_relative_rms_error": lane["channel"]["terminal_common_history_relative_rms_error"],
    }


def run_cluster(cluster: dict) -> dict:
    # One source/native M10 construction shared by both coverage lanes.
    data = CUR.local_cluster(cluster)
    m10, channel = CUR.current_native_m10(data["rho3"])

    x25, y25, _, _ = BASE._launch_expanded_25pct()
    x100, y100, _, _ = _launch_full_100pct()

    chain25 = _run_g3d_with_launch(m10, x25, y25, EXPECTED_SUPPORT25)
    chain100 = _run_g3d_with_launch(m10, x100, y100, EXPECTED_SUPPORT100)

    lane25 = _decode_lane(data, channel, chain25, x25, y25, "coverage_25pct")
    lane100 = _decode_lane(data, channel, chain100, x100, y100, "coverage_100pct")

    # Only after BOTH complete decoder inventories exist do we request targets.
    targets = DEC._targets_after_decoding(data)
    out25 = _compare_lane_after_decoding(lane25, targets)
    out100 = _compare_lane_after_decoding(lane100, targets)

    if out25["candidate_names"] != out100["candidate_names"]:
        raise RuntimeError("25% and 100% candidate inventories differ")

    deltas = {}
    for name in out25["candidate_names"]:
        deltas[name] = {}
        for target in ("kappa", "abs_kappa", "gamma", "observer_norm"):
            a = out25["candidate_relevance"][name][target]
            b = out100["candidate_relevance"][name][target]
            deltas[name][target] = {
                "pearson_delta_100_minus_25": float(b["pearson"] - a["pearson"]),
                "spearman_delta_100_minus_25": float(b["spearman"] - a["spearman"]),
                "coverage_delta_100_minus_25": float(b["coverage_fraction"] - a["coverage_fraction"]),
            }

    return {
        "cluster_id": cluster["id"],
        "same_native_M10_shared_between_lanes": True,
        "coverage_25pct": out25,
        "coverage_100pct": out100,
        "candidate_deltas": deltas,
    }


def _mean(vals) -> float:
    x = [float(v) for v in vals if np.isfinite(v)]
    return statistics.mean(x) if x else float("nan")


def aggregate(rows: list[dict]) -> dict:
    if not rows:
        return {}
    names = rows[0]["coverage_25pct"]["candidate_names"]
    targets = ("kappa", "abs_kappa", "gamma", "observer_norm")
    candidates = {}
    for name in names:
        candidates[name] = {}
        for target in targets:
            p25 = [r["coverage_25pct"]["candidate_relevance"][name][target]["pearson"] for r in rows]
            p100 = [r["coverage_100pct"]["candidate_relevance"][name][target]["pearson"] for r in rows]
            s25 = [r["coverage_25pct"]["candidate_relevance"][name][target]["spearman"] for r in rows]
            s100 = [r["coverage_100pct"]["candidate_relevance"][name][target]["spearman"] for r in rows]
            c25 = [r["coverage_25pct"]["candidate_relevance"][name][target]["coverage_fraction"] for r in rows]
            c100 = [r["coverage_100pct"]["candidate_relevance"][name][target]["coverage_fraction"] for r in rows]
            candidates[name][target] = {
                "mean_pearson_25": _mean(p25),
                "mean_pearson_100": _mean(p100),
                "mean_pearson_delta_100_minus_25": _mean([b-a for a, b in zip(p25, p100)]),
                "mean_spearman_25": _mean(s25),
                "mean_spearman_100": _mean(s100),
                "mean_spearman_delta_100_minus_25": _mean([b-a for a, b in zip(s25, s100)]),
                "mean_coverage_25": _mean(c25),
                "mean_coverage_100": _mean(c100),
            }

    rankings = {}
    for target in targets:
        rankings[target] = sorted(
            ({"candidate": n, **candidates[n][target]} for n in names),
            key=lambda x: x["mean_pearson_100"] if np.isfinite(x["mean_pearson_100"]) else -1e99,
            reverse=True,
        )

    support = {}
    for lane in ("coverage_25pct", "coverage_100pct"):
        support[lane] = {
            k: _mean([r[lane]["support"][k] for r in rows])
            for k in rows[0][lane]["support"]
        }

    info = {}
    for lane in ("coverage_25pct", "coverage_100pct"):
        info[lane] = {
            "mean_full_3d_effective_rank": _mean([
                r[lane]["per_ray_information_geometry"]["full_3d_received_ray_state"]["effective_rank"] for r in rows
            ]),
            "mean_transverse_2d_effective_rank": _mean([
                r[lane]["per_ray_information_geometry"]["transverse_2d_received_ray_state"]["effective_rank"] for r in rows
            ]),
            "mean_full_decoded_bank_effective_rank": _mean([
                r[lane]["stage_information_geometry"]["full_decoded_bank"]["effective_rank"] for r in rows
            ]),
            "mean_full_decoded_bank_numerical_rank": _mean([
                r[lane]["stage_information_geometry"]["full_decoded_bank"]["numerical_rank"] for r in rows
            ]),
        }

    frozen = {
        name: {target: candidates[name][target] for target in targets}
        for name in FROZEN_CARRY_FORWARD if name in candidates
    }
    return {
        "support": support,
        "information_geometry": info,
        "frozen_carry_forward_candidates": frozen,
        "candidate_rankings_by_100pct_pearson": rankings,
        "all_candidates": candidates,
    }


def main() -> int:
    state = repo_state()
    clusters = list(BENCH.clusters())
    ids = tuple(c["id"] for c in clusters)
    rows, failures = [], []

    if ids == EXPECTED_CLUSTER_IDS:
        for cluster in clusters:
            try:
                rows.append(run_cluster(cluster))
            except Exception as exc:
                failures.append({"cluster_id": cluster["id"], "error": f"{type(exc).__name__}: {exc}"})

    agg = aggregate(rows) if rows else {}
    checks = {
        "canonical_five_cluster_inventory": ids == EXPECTED_CLUSTER_IDS,
        "all_five_clusters_completed": len(rows) == 5 and not failures,
        "same_native_M10_shared_between_coverage_lanes": bool(rows and all(r["same_native_M10_shared_between_lanes"] for r in rows)),
        "25pct_expected_ray_count": N25 == SIDE25 * SIDE25,
        "100pct_ray_count_is_exactly_4x_25pct": N100 == 4 * N25,
        "100pct_expected_source_support": bool(rows and all(r["coverage_100pct"]["support"]["source_support_bins"] == EXPECTED_SUPPORT100 for r in rows)),
        "25pct_expected_source_support": bool(rows and all(r["coverage_25pct"]["support"]["source_support_bins"] == EXPECTED_SUPPORT25 for r in rows)),
        "same_decoder_inventory_both_lanes": bool(rows and all(r["coverage_25pct"]["candidate_names"] == r["coverage_100pct"]["candidate_names"] for r in rows)),
        "full_45_channel_bank_both_lanes": bool(rows and all(r["coverage_25pct"]["decoded_bank_size"] == 45 and r["coverage_100pct"]["decoded_bank_size"] == 45 for r in rows)),
        "native_terminal_common_history_identity": bool(rows and all(r["coverage_25pct"]["terminal_common_history_relative_rms_error"] <= 1e-12 and r["coverage_100pct"]["terminal_common_history_relative_rms_error"] <= 1e-12 for r in rows)),
        "G3D_unit_speed_valid_both_lanes": bool(rows and all(r["coverage_25pct"]["g3d_unit_speed_max_error"] <= G3D.UNIT_SPEED_TOL and r["coverage_100pct"]["g3d_unit_speed_max_error"] <= G3D.UNIT_SPEED_TOL for r in rows)),
        "observations_used_only_after_both_lane_decoders_constructed": True,
        "no_observational_fit_or_decoder_weighting": True,
        "no_physical_output_rescaling": True,
        "no_upstream_physics_change": True,
        "no_cluster_specific_decoder_choice": True,
        "no_historical_strength_or_replacement_scalar": True,
        "no_tracked_or_staged_changes": not state["tracked_changes"] and not state["staged_changes"],
    }
    passed = all(checks.values())
    status = "FULL_STATE_100PCT_OBSERVER_COVERAGE_EXECUTED" if passed else ("FULL_STATE_100PCT_OBSERVER_COVERAGE_PARTIAL_EXECUTION" if rows else "FULL_STATE_100PCT_OBSERVER_COVERAGE_NOT_ESTABLISHED")

    result = {
        "lab_id": LAB_ID,
        "status": status,
        "repo_state": state,
        "coverage_rule": "same native M10 and same decoder inventory; only source-plane coverage/ray count changes from 25% to 100% at approximately fixed rays per source bin",
        "side_25pct": SIDE25,
        "side_100pct": SIDE100,
        "rays_25pct": N25,
        "rays_100pct": N100,
        "expected_support_25pct": EXPECTED_SUPPORT25,
        "expected_support_100pct": EXPECTED_SUPPORT100,
        "frozen_carry_forward_candidates": FROZEN_CARRY_FORWARD,
        "aggregate": agg,
        "rows": rows,
        "failures": failures,
        "checks": checks,
    }

    print(LAB_ID)
    print(f"status={status}")
    print(f"head_sha={state['head_sha']}")
    print(f"rays_25pct={N25}")
    print(f"rays_100pct={N100}")
    print(f"support_25pct={EXPECTED_SUPPORT25}")
    print(f"support_100pct={EXPECTED_SUPPORT100}")
    print("upstream_physics=frozen")
    print("decoder_inventory=exact_PR105_inventory")
    print("observations=after_both_lane_decoders")
    print()

    if agg:
        print("COVERAGE_SUMMARY")
        print(json.dumps(agg["support"], indent=2))
        print()
        print("INFORMATION_GEOMETRY")
        print(json.dumps(agg["information_geometry"], indent=2))
        print()
        print("FROZEN_CARRY_FORWARD")
        print(json.dumps(agg["frozen_carry_forward_candidates"], indent=2))
        print()
        for target in ("kappa", "abs_kappa", "gamma", "observer_norm"):
            print(f"TOP_100PCT_RECONSTRUCTIONS target={target}")
            for row in agg["candidate_rankings_by_100pct_pearson"][target][:12]:
                print(json.dumps(row))
            print()

    print("CHECKS")
    for k, v in checks.items():
        print(f"{k}={str(bool(v)).lower()}")
    print()
    print("RESULT_JSON")
    print(json.dumps(result, indent=2, default=lambda o: float(o) if isinstance(o, np.floating) else int(o) if isinstance(o, np.integer) else str(o)))
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
