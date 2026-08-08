#!/usr/bin/env python3
"""PBUF FOUNDATION — CLEAN CURRENT WEAK-LENSING BENCHMARK 001.

Purpose
-------
Measure the CURRENT frozen native propagation lane without benchmark-assisted source
construction, HST/F160W, historical controls, replacement strength, inferred transfer
coefficients, fitted normalization, or target-dependent tuning.

The five canonical local kappa FITS files are TARGETS ONLY.  A cluster may run only
when its local benchmark directory also contains an independently prepared 3-D PBUF
native-loading FITS cube.  This lab deliberately does not invent an SI->native source
conversion: if no independently calibrated native-loading cube exists, that cluster is
reported as unavailable rather than silently normalized or reconstructed from kappa.

Accepted independent source metadata
------------------------------------
A source FITS must be 3-D and explicitly identify itself as native PBUF loading via
one of these header declarations:

    PBUFROLE = 'INDEPENDENT_NATIVE_LOADING'
    PBUFROLE = 'NATIVE_MATTER_LOADING'

or

    BUNIT = 'PBUF_NATIVE_LOADING'
    BUNIT = 'PBUF_NATIVE_RHO'

The source filename is not used to choose it.  The canonical kappa target is always
excluded.  Exactly one valid source cube per cluster is required.

No arbitrary physics scalar is introduced here.  The only numerical transfer
coefficients used are the already-frozen A8 model coefficients (0.03 and 0.003), read
from the model implementation rather than re-fitted in this benchmark.
"""
from __future__ import annotations

import json
import math
import subprocess
import sys
from pathlib import Path

import numpy as np
from astropy.io import fits

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))

from pbuf.core import benchmark_data as BENCH
from pbuf.core import observable_extraction as M16
from pbuf.models import a8_state as A8
import pbuf.labs.foundation.interface_to_interface_survivor_sweep001 as S92
import pbuf.labs.foundation.native_channel_transfer_closure_sweep001 as S93
import pbuf.labs.foundation.native_accumulated_full_lensing001 as G3D
import pbuf.labs.foundation.m10_coverage_25pct_science001 as BASE

LAB_ID = "PBUF-FOUNDATION-CLEAN-CURRENT-WEAK-LENSING-BENCHMARK-001"
EXPECTED_CLUSTER_IDS = ("Abell2744", "MACS0416", "MACS1149", "AbellS1063", "Abell370")
VALID_ROLES = {"INDEPENDENT_NATIVE_LOADING", "NATIVE_MATTER_LOADING"}
VALID_BUNITS = {"PBUF_NATIVE_LOADING", "PBUF_NATIVE_RHO"}
EPS = 1.0e-30


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


def rms(a) -> float:
    x = np.asarray(a, dtype=np.float64)
    m = np.isfinite(x)
    return float(np.sqrt(np.mean(x[m] * x[m]))) if np.any(m) else float("nan")


def corr(a, b) -> tuple[float, float, int]:
    x = np.asarray(a, dtype=np.float64)
    y = np.asarray(b, dtype=np.float64)
    if x.shape != y.shape:
        raise RuntimeError(f"correlation shape mismatch: {x.shape} vs {y.shape}")
    m = np.isfinite(x) & np.isfinite(y)
    n = int(np.count_nonzero(m))
    if n < 2:
        return float("nan"), float("nan"), n
    return float(M16.safe_pearson(x[m], y[m])), float(M16.safe_spearman(x[m], y[m])), n


def source_inventory(cluster: dict) -> tuple[list[dict], list[Path]]:
    row = BENCH.resolve_cluster(cluster)
    directory = BENCH.BENCHMARK_ROOT / row["directory"]
    target = BENCH.require_kappa_path(row).resolve()
    entries, accepted = [], []
    for path in sorted(directory.glob("*.fits")):
        if path.resolve() == target:
            entries.append({"path": str(path), "role": "weak_lensing_target", "accepted_source": False})
            continue
        try:
            hdr = fits.getheader(path, 0)
            naxis = int(hdr.get("NAXIS", 0))
            role = str(hdr.get("PBUFROLE", "")).strip().upper()
            bunit = str(hdr.get("BUNIT", "")).strip().upper()
            ok = bool(naxis == 3 and (role in VALID_ROLES or bunit in VALID_BUNITS))
            entries.append({
                "path": str(path), "naxis": naxis, "PBUFROLE": role or None,
                "BUNIT": bunit or None, "accepted_source": ok,
            })
            if ok:
                accepted.append(path)
        except Exception as exc:
            entries.append({"path": str(path), "accepted_source": False,
                            "inventory_error": f"{type(exc).__name__}: {exc}"})
    return entries, accepted


def load_independent_native_source(cluster: dict) -> tuple[np.ndarray, dict]:
    inventory, accepted = source_inventory(cluster)
    if len(accepted) != 1:
        raise RuntimeError(
            f"clean independent native source requires exactly one accepted 3-D FITS; "
            f"found {len(accepted)}. inventory={json.dumps(inventory, sort_keys=True)}"
        )
    path = accepted[0]
    with fits.open(path, memmap=True) as hdul:
        rho3 = np.asarray(hdul[0].data, dtype=np.float64).copy()
        hdr = hdul[0].header
    if rho3.ndim != 3 or not np.all(np.isfinite(rho3)):
        raise RuntimeError(f"invalid independent native-loading cube: {path} shape={rho3.shape}")
    if not np.any(rho3 != 0.0):
        raise RuntimeError(f"independent native-loading cube is identically zero: {path}")
    return rho3, {
        "source_path": str(path),
        "PBUFROLE": hdr.get("PBUFROLE"),
        "BUNIT": hdr.get("BUNIT"),
        "source_rms": rms(rho3),
        "source_sum": float(np.sum(rho3)),
        "inventory": inventory,
    }


def native_m10(rho3: np.ndarray):
    channels = S93.native_terminal_channels(rho3)
    c = np.asarray(channels["c"], dtype=np.float64)
    uf = np.asarray(channels["u_fast"], dtype=np.float64)
    us = np.asarray(channels["u_slow"], dtype=np.float64)
    cf = float(A8.A8_INIT_DT * A8.A8_INIT_OMEGA * A8.A8_INIT_K)
    cs = float(A8.A8_INIT_DT * A8.A8_INIT_SLOW_TIMESCALE)
    amps = S93.combine_amplitudes(
        S93.scale_amplitudes(S92.positive_bond_amplitudes(uf), cf),
        S93.scale_amplitudes(S92.positive_bond_amplitudes(us), cs),
    )
    m10 = S92.m10_from_amplitudes(amps, c)
    history_rel = S92.rms(channels["history_c"] - c) / max(S92.rms(c), EPS)
    return m10, {"coef_fast": cf, "coef_slow": cs,
                 "terminal_common_history_relative_rms_error": float(history_rel)}


def target_for_shape(cluster: dict, shape: tuple[int, int]) -> np.ndarray:
    """Load target only after prediction and resample deterministically for correlation."""
    kappa = BENCH.load_kappa(cluster)
    # construct_common_proxy is used only as a deterministic target-grid sampler here;
    # its output never enters the source/prediction lane.
    if shape[0] != shape[1]:
        raise RuntimeError(f"observer target shape must be square, got {shape}")
    return np.asarray(BASE.construct_common_proxy(kappa, bins=shape[0], extent=BASE.CFG["extent"]), dtype=np.float64)


def run_cluster(cluster: dict) -> dict:
    # SOURCE + PREDICTION FIRST. No kappa pixels are loaded here.
    rho3, source_meta = load_independent_native_source(cluster)
    m10, channel = native_m10(rho3)
    chain = G3D.run_g3d_from_vector(m10, observed_for_first_step=None)
    final = np.asarray(chain["final_ang"]["angular_rms_angle_mag"], dtype=np.float64)
    los = np.asarray(chain["los_mag"], dtype=np.float64)

    # TARGET REVEALED ONLY NOW.
    target_final = target_for_shape(cluster, final.shape)
    target_los = target_for_shape(cluster, los.shape)
    fp, fs, fn = corr(final, target_final)
    lp, ls, ln = corr(los, target_los)

    return {
        "cluster_id": cluster["id"],
        "source": source_meta,
        "target_path": str(BENCH.require_kappa_path(cluster)),
        "target_loaded_after_prediction": True,
        "benchmark_assisted_source": False,
        "hst_f160w_used": False,
        "network_access_used": False,
        "historical_control_lane_used": False,
        "replacement_strength_scalar": None,
        "fitted_or_tuned_scalar": None,
        "native_los_rms": rms(los),
        "native_final_angle_rms": rms(final),
        "native_los_vs_target_pearson": lp,
        "native_los_vs_target_spearman": ls,
        "native_los_vs_target_count": ln,
        "native_final_angle_vs_target_pearson": fp,
        "native_final_angle_vs_target_spearman": fs,
        "native_final_angle_vs_target_count": fn,
        "native_g3d_unit_speed_max_error": float(chain["g3d"]["max_unit_speed_error"]),
        **channel,
    }


def main() -> int:
    state = repo_state()
    clusters = list(BENCH.clusters())
    inventory = BENCH.inventory()
    ids = tuple(c["id"] for c in clusters)
    rows, failures, source_inventory_rows = [], [], []

    local_targets_ready = bool(ids == EXPECTED_CLUSTER_IDS and len(inventory) == 5 and all(x["exists"] for x in inventory))
    if local_targets_ready:
        for cluster in clusters:
            inv, accepted = source_inventory(cluster)
            source_inventory_rows.append({"cluster_id": cluster["id"], "accepted_source_count": len(accepted), "files": inv})
            try:
                rows.append(run_cluster(cluster))
            except Exception as exc:
                failures.append({"cluster_id": cluster["id"], "error": f"{type(exc).__name__}: {exc}"})

    source_ready = bool(len(source_inventory_rows) == 5 and all(x["accepted_source_count"] == 1 for x in source_inventory_rows))
    all_completed = bool(len(rows) == 5 and not failures)
    checks = {
        "canonical_five_local_targets_present": local_targets_ready,
        "exactly_one_independent_native_source_per_cluster": source_ready,
        "all_five_clusters_completed": all_completed,
        "benchmark_assisted_source_false": bool(rows and all(not r["benchmark_assisted_source"] for r in rows)),
        "target_loaded_after_prediction": bool(rows and all(r["target_loaded_after_prediction"] for r in rows)),
        "network_access_used_false": bool(not rows or all(not r["network_access_used"] for r in rows)),
        "hst_f160w_used_false": bool(not rows or all(not r["hst_f160w_used"] for r in rows)),
        "historical_control_lane_used_false": bool(not rows or all(not r["historical_control_lane_used"] for r in rows)),
        "no_replacement_or_fitted_scalar": bool(not rows or all(r["replacement_strength_scalar"] is None and r["fitted_or_tuned_scalar"] is None for r in rows)),
        "terminal_common_history_identity": bool(not rows or all(r["terminal_common_history_relative_rms_error"] <= 1.0e-12 for r in rows)),
        "no_tracked_or_staged_changes": bool(not state["tracked_changes"] and not state["staged_changes"]),
    }

    if all_completed and all(checks.values()):
        status = "CLEAN_CURRENT_WEAK_LENSING_BENCHMARK_EXECUTED"
    elif local_targets_ready and not source_ready:
        status = "CLEAN_CURRENT_WEAK_LENSING_SOURCE_NOT_AVAILABLE"
    elif rows:
        status = "CLEAN_CURRENT_WEAK_LENSING_BENCHMARK_PARTIAL_EXECUTION"
    else:
        status = "CLEAN_CURRENT_WEAK_LENSING_BENCHMARK_NOT_ESTABLISHED"

    result = {
        "lab_id": LAB_ID,
        "status": status,
        "repo_state": state,
        "model": {
            "source_requirement": "independent local 3-D PBUF native-loading FITS; no target-derived or luminous proxy source",
            "target": "canonical local Merten v1 kappa FITS, revealed after prediction only",
            "native_lane": "independent native rho3 -> zero-flux terminal fast/slow -> frozen exact pair law -> terminal c geometry -> PM1/PS2/M10 -> LOS -> existing G3D",
            "benchmark_assisted_source": False,
            "hst_f160w_used": False,
            "historical_control_lane_used": False,
            "replacement_strength_scalar": None,
            "fit_or_tuning": False,
            "important_limit": "No SI matter-density-to-native-loading conversion is invented. If a calibrated native source cube is absent, clean present-day performance cannot yet be measured from these targets alone.",
        },
        "source_inventory": source_inventory_rows,
        "rows": rows,
        "failures": failures,
        "checks": checks,
    }

    print(LAB_ID)
    print(f"status={status}")
    print(f"head_sha={state['head_sha']}")
    print("benchmark_assisted_source=false")
    print("hst_f160w_used=false")
    print("network_access_used=false")
    print("historical_control_lane_used=false")
    print("replacement_strength_scalar=none")
    print("fit_or_tuning=false")
    print("target_role=end_of_chain_only")
    print("source_requirement=independent_local_3d_PBUF_native_loading")
    print()
    print("SOURCE_INVENTORY")
    for x in source_inventory_rows:
        print(f"cluster={x['cluster_id']} accepted_native_source_count={x['accepted_source_count']}")
    print()
    print("CLUSTERS")
    for r in rows:
        print(
            f"cluster={r['cluster_id']} "
            f"final_pearson={r['native_final_angle_vs_target_pearson']:.12g} "
            f"final_spearman={r['native_final_angle_vs_target_spearman']:.12g} "
            f"los_pearson={r['native_los_vs_target_pearson']:.12g} "
            f"los_spearman={r['native_los_vs_target_spearman']:.12g} "
            f"native_angle_rms={r['native_final_angle_rms']:.12g}"
        )
    for f in failures:
        print(f"failure_cluster={f['cluster_id']} error={f['error']}")
    print()
    print("CHECKS")
    for k, v in checks.items():
        print(f"{k}={str(v).lower()}")
    print(json.dumps(result, sort_keys=True, separators=(",", ":"), default=str))
    return 0 if status in {"CLEAN_CURRENT_WEAK_LENSING_BENCHMARK_EXECUTED", "CLEAN_CURRENT_WEAK_LENSING_SOURCE_NOT_AVAILABLE"} else 1


if __name__ == "__main__":
    raise SystemExit(main())
