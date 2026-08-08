#!/usr/bin/env python3
"""PBUF FOUNDATION — NATIVE CHANNEL FULL WEAK LENSING FIVE CLUSTER 001.

Purpose
-------
Run the newly closed native terminal fast+slow transfer lane through the existing
M10 -> LOS -> G3D -> observer stack on all five canonical local weak-lensing
benchmarks, then reveal observed kappa only after every prediction lane is complete.

Native prediction lane:
    independent rho3 source proxy
    -> frozen zero-flux A8 terminal u_fast/u_slow
    -> exact frozen pair law 0.03*Delta u_fast + 0.003*Delta u_slow
    -> native c_state geometry
    -> frozen PM1/PS2/M10
    -> frozen LOS
    -> existing G3D/observer

No replacement strength scalar, inferred transfer coefficient, normalization,
rescaling, or observational fitting is permitted.
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

import pbuf.labs.foundation.native_channel_transfer_closure_sweep001 as C93
import pbuf.labs.foundation.interface_to_interface_survivor_sweep001 as S92
import pbuf.labs.foundation.native_accumulated_full_lensing001 as FULL
import pbuf.labs.foundation.independent_source_training_wheels_off001_common_footprint_fix as SRC
import pbuf.labs.foundation.independent_source_training_wheels_off001 as IND
import pbuf.labs.foundation.m10_coverage_25pct_science001 as BASE
from pbuf.core import los_projection as M14

LAB_ID = "PBUF-FOUNDATION-NATIVE-CHANNEL-FULL-WEAK-LENSING-FIVE-CLUSTER-001"
EXPECTED_CLUSTER_IDS = ("Abell2744", "MACS0416", "MACS1149", "AbellS1063", "Abell370")
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


def native_full_pair_amplitudes(rho3: np.ndarray) -> tuple[dict, np.ndarray, dict]:
    channels = C93.native_terminal_channels(rho3)
    uf = np.asarray(channels["u_fast"], dtype=np.float64)
    us = np.asarray(channels["u_slow"], dtype=np.float64)
    c = np.asarray(channels["c"], dtype=np.float64)
    cf = 0.03
    cs = 0.003
    amps = C93.combine_amplitudes(
        C93.scale_amplitudes(S92.positive_bond_amplitudes(uf), cf),
        C93.scale_amplitudes(S92.positive_bond_amplitudes(us), cs),
    )
    return amps, c, channels


def native_prediction_chain(rho3: np.ndarray) -> dict:
    amps, cgeom, channels = native_full_pair_amplitudes(rho3)
    m10 = S92.m10_from_amplitudes(amps, cgeom)
    g3d = FULL.run_g3d_from_vector(m10, observed_for_first_step=None)
    history_terminal = np.asarray(channels["history_c"], dtype=np.float64)
    terminal_common_relerr = rms(history_terminal - cgeom) / max(rms(cgeom), EPS)
    return {
        **g3d,
        "m10": m10,
        "c_state": cgeom,
        "terminal_common_history_relative_rms_error": terminal_common_relerr,
    }


def compare_prediction(prefix: str, rho2: np.ndarray, chain: dict, observed: np.ndarray) -> dict:
    return FULL.compare_lane(prefix, rho2, chain, observed)


def run_cluster(cluster: dict) -> dict:
    # Prediction source only. Observed kappa values are intentionally unavailable
    # until all prediction lanes below have completed.
    source = SRC._independent_source(cluster)
    rho3 = np.asarray(source["rho3"], dtype=np.float64)

    native = native_prediction_chain(rho3)
    legacy = FULL.legacy_chain(rho3)
    unit = FULL.unit_control_chain(rho3)

    # Only now reveal the benchmark weak-lensing map.
    kpath = IND._kappa_path(cluster)
    with fits.open(kpath) as hdul:
        kappa_native = np.asarray(hdul[0].data, dtype=np.float64)
    observed, _ = SRC._benchmark_on_common_grid(kappa_native, source["geometry"])

    native_metrics = compare_prediction("native_channel", source["rho2"], native, observed)
    legacy_metrics = compare_prediction("legacy", source["rho2"], legacy, observed)
    unit_metrics = compare_prediction("unit_control", source["rho2"], unit, observed)

    return {
        "cluster_id": cluster["id"],
        "source_role": source["source_role"],
        "source_limit": "HST_F160W_luminous_structure_proxy_not_absolute_baryonic_mass_map",
        "benchmark_values_loaded_after_all_prediction_lanes_complete": True,
        "native_role": "zero_flux_terminal_fast_plus_slow_pair_transfer_to_native_c_geometry_to_PM1_PS2_M10_LOS_G3D",
        "terminal_common_history_relative_rms_error": native["terminal_common_history_relative_rms_error"],
        "native_los_rms": rms(native["los_mag"]),
        "legacy_los_rms": rms(legacy["los_mag"]),
        "unit_los_rms": rms(unit["los_mag"]),
        "native_over_unit_los_rms": rms(native["los_mag"]) / max(rms(unit["los_mag"]), EPS),
        "native_over_legacy_los_rms": rms(native["los_mag"]) / max(rms(legacy["los_mag"]), EPS),
        "native_final_angle_rms": rms(native["final_ang"]["angular_rms_angle_mag"]),
        "legacy_final_angle_rms": rms(legacy["final_ang"]["angular_rms_angle_mag"]),
        "unit_final_angle_rms": rms(unit["final_ang"]["angular_rms_angle_mag"]),
        "native_g3d_unit_speed_max_error": native["g3d"]["max_unit_speed_error"],
        **native_metrics,
        **legacy_metrics,
        **unit_metrics,
    }


def main() -> int:
    state = repo_state()
    clusters = list(BASE.CLUSTERS)
    ids = tuple(c["id"] for c in clusters)
    rows, failures = [], []

    if ids == EXPECTED_CLUSTER_IDS:
        for cluster in clusters:
            try:
                rows.append(run_cluster(cluster))
            except Exception as exc:
                failures.append({"cluster_id": cluster["id"], "error": f"{type(exc).__name__}: {exc}"})

    checks = {
        "canonical_five_cluster_inventory": ids == EXPECTED_CLUSTER_IDS,
        "all_five_clusters_completed": bool(len(rows) == 5 and not failures),
        "terminal_common_history_identity": bool(rows and all(r["terminal_common_history_relative_rms_error"] <= 1.0e-12 for r in rows)),
        "native_G3D_valid": bool(rows and all(r["native_g3d_unit_speed_max_error"] <= FULL.UNIT_SPEED_TOL for r in rows)),
        "native_signal_nonzero": bool(rows and all(r["native_los_rms"] > 0.0 and r["native_final_angle_rms"] > 0.0 for r in rows)),
        "benchmark_revealed_after_predictions": bool(rows and all(r["benchmark_values_loaded_after_all_prediction_lanes_complete"] for r in rows)),
        "all_primary_observation_metrics_finite": bool(rows and all(math.isfinite(r["native_channel_final_angular_rms_angle_mag_vs_observed_pearson"]) and math.isfinite(r["native_channel_final_angular_rms_angle_mag_vs_observed_spearman"]) for r in rows)),
        "no_tracked_or_staged_changes": bool(not state["tracked_changes"] and not state["staged_changes"]),
    }
    execution_gate_pass = bool(all(checks.values()))
    status = (
        "NATIVE_CHANNEL_FULL_WEAK_LENSING_FIVE_CLUSTER_EXECUTED" if execution_gate_pass
        else "NATIVE_CHANNEL_FULL_WEAK_LENSING_FIVE_CLUSTER_PARTIAL_EXECUTION" if rows
        else "NATIVE_CHANNEL_FULL_WEAK_LENSING_FIVE_CLUSTER_NOT_ESTABLISHED"
    )

    native_better_than_legacy_pearson = sum(
        1 for r in rows
        if r["native_channel_final_angular_rms_angle_mag_vs_observed_pearson"]
        > r["legacy_final_angular_rms_angle_mag_vs_observed_pearson"]
    )
    native_better_than_unit_pearson = sum(
        1 for r in rows
        if r["native_channel_final_angular_rms_angle_mag_vs_observed_pearson"]
        > r["unit_control_final_angular_rms_angle_mag_vs_observed_pearson"]
    )

    result = {
        "lab_id": LAB_ID,
        "status": status,
        "repo_state": state,
        "model": {
            "native_prediction": "rho3 -> zero-flux terminal fast/slow -> exact frozen pair law -> native c geometry -> PM1/PS2/M10 -> LOS -> G3D",
            "pair_law": "A_ij=0.03*Delta_u_fast+0.003*Delta_u_slow",
            "replacement_strength_scalar": None,
            "inferred_effective_coefficient_applied": False,
            "normalization_or_rescaling": False,
            "fit_or_tuning": False,
            "observed_kappa_role": "revealed only after all prediction lanes complete",
            "source_limit": "normalized_HST_F160W_luminous_structure_proxy_not_absolute_mass_map",
        },
        "rows": rows,
        "failures": failures,
        "checks": checks,
        "execution_gate_pass": execution_gate_pass,
        "native_better_than_legacy_pearson_count_of_5": native_better_than_legacy_pearson,
        "native_better_than_unit_pearson_count_of_5": native_better_than_unit_pearson,
        "interpretation_rule": "Report end-to-end weak-lensing performance without fitting or rescaling. Do not alter the native lane based on benchmark agreement in this run.",
    }

    print(LAB_ID)
    print(f"status={status}")
    print(f"head_sha={state['head_sha']}")
    print("network_access_used=false")
    print("observed_lensing_values_used=true_end_of_chain_only")
    print("replacement_strength_scalar=none")
    print("inferred_effective_coefficient_applied=false")
    print("native_response_rescaled=false")
    print("fit_or_tuning=false")
    print("existing_G3D_ray_tracker_reused=true")
    print("pair_law=A_ij=0.03*Delta_u_fast+0.003*Delta_u_slow")
    print()
    print("CLUSTERS")
    for r in rows:
        print(
            f"cluster={r['cluster_id']} "
            f"native_over_unit_los={r['native_over_unit_los_rms']:.12g} "
            f"native_over_legacy_los={r['native_over_legacy_los_rms']:.12g} "
            f"native_obs_pearson={r['native_channel_final_angular_rms_angle_mag_vs_observed_pearson']:.12g} "
            f"native_obs_spearman={r['native_channel_final_angular_rms_angle_mag_vs_observed_spearman']:.12g} "
            f"legacy_obs_pearson={r['legacy_final_angular_rms_angle_mag_vs_observed_pearson']:.12g} "
            f"unit_obs_pearson={r['unit_control_final_angular_rms_angle_mag_vs_observed_pearson']:.12g} "
            f"native_angle_rms={r['native_final_angle_rms']:.12g}"
        )
    for f in failures:
        print(f"failure_cluster={f['cluster_id']} error={f['error']}")
    print()
    print("SUMMARY")
    print(f"native_better_than_legacy_pearson_count_of_5={native_better_than_legacy_pearson}")
    print(f"native_better_than_unit_pearson_count_of_5={native_better_than_unit_pearson}")
    print("CHECKS")
    for k, v in checks.items():
        print(f"{k}={str(v).lower()}")
    print(f"execution_gate_pass={str(execution_gate_pass).lower()}")
    print(json.dumps(result, sort_keys=True, separators=(",", ":"), default=str))
    return 0 if execution_gate_pass else 1


if __name__ == "__main__":
    raise SystemExit(main())
