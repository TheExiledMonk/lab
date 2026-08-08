#!/usr/bin/env python3
"""PBUF FOUNDATION — LOCAL BENCHMARK NATIVE CHANNEL FULL LENSING 001.

Purpose
-------
Run the frozen native terminal fast+slow propagation interface end-to-end through
PM1/PS2/M10 -> LOS -> existing G3D/observer on all five canonical local weak-
lensing benchmark FITS files, using exactly the same benchmark-assisted local
source construction employed by the recent native bridge audits (#85-#93).

This is deliberately a benchmark-assisted end-to-end structural test, NOT an
independent source prediction: the local kappa morphology is used to construct the
frozen rho2/rho3 benchmark proxy and is also the comparison morphology. The point
of this lab is to test the repaired native propagation lane consistently across
all five local files without any HST/F160W/network code.

There is NO URL discovery, NO download, NO HST source helper, and NO network
fallback anywhere in this lab. Missing local benchmark files are a hard failure.
"""
from __future__ import annotations

import json
import math
import subprocess
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))

from a8_three_dimensional_projection_lab001 import construct_rho_3d
from pbuf.core import benchmark_data as BENCH
from pbuf.core import observable_extraction as M16
from pbuf.core import los_projection as M14
import pbuf.labs.foundation.m10_coverage_25pct_science001 as BASE
import pbuf.labs.foundation.interface_to_interface_survivor_sweep001 as S92
import pbuf.labs.foundation.native_channel_transfer_closure_sweep001 as S93
import pbuf.labs.foundation.native_accumulated_full_lensing001 as G3D

LAB_ID = "PBUF-FOUNDATION-LOCAL-BENCHMARK-NATIVE-CHANNEL-FULL-LENSING-001"
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


def local_source(cluster: dict) -> dict:
    """Exact recent-audit source construction: local kappa -> common proxy -> rho3."""
    path = BENCH.require_kappa_path(cluster)
    kappa = BENCH.load_kappa(cluster)
    rho2 = BASE.construct_common_proxy(kappa, bins=BASE.OBS_BINS, extent=BASE.CFG["extent"])
    rho3 = construct_rho_3d(rho2, BASE.NZ, profile=BASE.PROFILE)
    return {
        "path": str(path),
        "kappa": np.asarray(kappa, dtype=np.float64),
        "rho2": np.asarray(rho2, dtype=np.float64),
        "rho3": np.asarray(rho3, dtype=np.float64),
    }


def native_m10(rho3: np.ndarray) -> tuple[tuple[np.ndarray, np.ndarray, np.ndarray], dict]:
    """Frozen PR #93 full native terminal fast+slow transfer with native c geometry."""
    channels = S93.native_terminal_channels(rho3)
    c = np.asarray(channels["c"], dtype=np.float64)
    uf = np.asarray(channels["u_fast"], dtype=np.float64)
    us = np.asarray(channels["u_slow"], dtype=np.float64)
    cf = float(BASE.A8.A8_INIT_DT * BASE.A8.A8_INIT_OMEGA * BASE.A8.A8_INIT_K) if hasattr(BASE, "A8") else 0.03
    cs = 0.003
    # Guard against silent coefficient drift.
    if abs(cf - 0.03) > 1.0e-15 or abs(cs - 0.003) > 1.0e-15:
        raise RuntimeError(f"frozen pair coefficients changed: fast={cf}, slow={cs}")
    amps = S93.combine_amplitudes(
        S93.scale_amplitudes(S92.positive_bond_amplitudes(uf), cf),
        S93.scale_amplitudes(S92.positive_bond_amplitudes(us), cs),
    )
    m10 = S92.m10_from_amplitudes(amps, c)
    history_rel = S92.rms(channels["history_c"] - c) / max(S92.rms(c), EPS)
    return m10, {
        "terminal_common_history_relative_rms_error": float(history_rel),
        "coef_fast": cf,
        "coef_slow": cs,
    }


def unit_m10(rho3: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    hist = S92.historical_state_and_interface(rho3)
    return tuple(np.asarray(v, dtype=np.float64) for v in hist["m10"])


def run_lane(vector) -> dict:
    return G3D.run_g3d_from_vector(vector, observed_for_first_step=None)


def lane_metrics(prefix: str, chain: dict, benchmark_proxy: np.ndarray) -> dict:
    fields = {
        "los_mag": chain["los_mag"],
        "final_angular_centroid_mag": chain["final_ang"]["angular_centroid_mag"],
        "final_angular_spread_rms": chain["final_ang"]["angular_spread_rms"],
        "final_angular_rms_angle_mag": chain["final_ang"]["angular_rms_angle_mag"],
    }
    out = {}
    for name, field in fields.items():
        if np.asarray(field).shape != np.asarray(benchmark_proxy).shape:
            raise RuntimeError(
                f"{prefix}:{name} shape {np.asarray(field).shape} != benchmark proxy {benchmark_proxy.shape}"
            )
        p, s, n = corr(field, benchmark_proxy)
        out[f"{prefix}_{name}_rms"] = rms(field)
        out[f"{prefix}_{name}_vs_benchmark_pearson"] = p
        out[f"{prefix}_{name}_vs_benchmark_spearman"] = s
        out[f"{prefix}_{name}_vs_benchmark_count"] = n
    return out


def run_cluster(cluster: dict) -> dict:
    source = local_source(cluster)
    rho3 = source["rho3"]
    benchmark_proxy = source["rho2"]

    # Build both prediction/response lanes solely from the local benchmark-assisted
    # source. No HST, URL, download, or remote data path exists here.
    nm10, channel = native_m10(rho3)
    um10 = unit_m10(rho3)
    native = run_lane(nm10)
    unit = run_lane(um10)

    native_metrics = lane_metrics("native", native, benchmark_proxy)
    unit_metrics = lane_metrics("unit_control", unit, benchmark_proxy)

    return {
        "cluster_id": cluster["id"],
        "local_benchmark_path": source["path"],
        "source_role": "local_kappa_derived_frozen_benchmark_proxy_structural_control_not_independent_prediction",
        "network_access_used": False,
        "hst_f160w_used": False,
        "benchmark_assisted": True,
        "rho2_rms": rms(source["rho2"]),
        "rho3_rms": rms(source["rho3"]),
        "native_los_rms": rms(native["los_mag"]),
        "unit_los_rms": rms(unit["los_mag"]),
        "native_over_unit_los_rms": rms(native["los_mag"]) / max(rms(unit["los_mag"]), EPS),
        "native_final_angle_rms": rms(native["final_ang"]["angular_rms_angle_mag"]),
        "unit_final_angle_rms": rms(unit["final_ang"]["angular_rms_angle_mag"]),
        "native_over_unit_final_angle_rms": rms(native["final_ang"]["angular_rms_angle_mag"]) / max(rms(unit["final_ang"]["angular_rms_angle_mag"]), EPS),
        "native_g3d_unit_speed_max_error": float(native["g3d"]["max_unit_speed_error"]),
        "unit_g3d_unit_speed_max_error": float(unit["g3d"]["max_unit_speed_error"]),
        **channel,
        **native_metrics,
        **unit_metrics,
    }


def main() -> int:
    state = repo_state()
    clusters = list(BENCH.clusters())
    inventory = BENCH.inventory()
    ids = tuple(c["id"] for c in clusters)
    rows, failures = [], []

    ready = bool(ids == EXPECTED_CLUSTER_IDS and len(inventory) == 5 and all(x["exists"] for x in inventory))
    if ready:
        for cluster in clusters:
            try:
                rows.append(run_cluster(cluster))
            except Exception as exc:
                failures.append({"cluster_id": cluster["id"], "error": f"{type(exc).__name__}: {exc}"})

    checks = {
        "canonical_five_cluster_inventory": ids == EXPECTED_CLUSTER_IDS,
        "all_five_local_fits_present": ready,
        "all_five_clusters_completed": bool(len(rows) == 5 and not failures),
        "terminal_common_history_identity": bool(rows and all(r["terminal_common_history_relative_rms_error"] <= 1.0e-12 for r in rows)),
        "native_G3D_valid": bool(rows and all(r["native_g3d_unit_speed_max_error"] <= G3D.UNIT_SPEED_TOL for r in rows)),
        "unit_G3D_valid": bool(rows and all(r["unit_g3d_unit_speed_max_error"] <= G3D.UNIT_SPEED_TOL for r in rows)),
        "all_primary_metrics_finite": bool(rows and all(
            math.isfinite(float(r[k])) for r in rows for k in (
                "native_over_unit_los_rms",
                "native_over_unit_final_angle_rms",
                "native_final_angular_rms_angle_mag_vs_benchmark_pearson",
                "native_final_angular_rms_angle_mag_vs_benchmark_spearman",
            )
        )),
        "network_access_used_false": bool(rows and all(not r["network_access_used"] for r in rows)),
        "hst_f160w_used_false": bool(rows and all(not r["hst_f160w_used"] for r in rows)),
        "no_tracked_or_staged_changes": bool(not state["tracked_changes"] and not state["staged_changes"]),
    }
    passed = bool(all(checks.values()))
    status = (
        "LOCAL_BENCHMARK_NATIVE_CHANNEL_FULL_LENSING_EXECUTED" if passed
        else ("LOCAL_BENCHMARK_NATIVE_CHANNEL_FULL_LENSING_PARTIAL_EXECUTION" if rows
              else "LOCAL_BENCHMARK_NATIVE_CHANNEL_FULL_LENSING_NOT_ESTABLISHED")
    )

    result = {
        "lab_id": LAB_ID,
        "status": status,
        "repo_state": state,
        "model": {
            "native_lane": "local benchmark proxy -> zero-flux terminal fast/slow -> exact frozen pair law -> native c geometry -> PM1/PS2/M10 -> LOS -> existing G3D",
            "unit_control": "same local benchmark proxy -> historical unit-loading A8/PM1/PS2/M10 -> LOS -> existing G3D",
            "source_role": "benchmark-assisted structural response test; not an independent source prediction",
            "network_access_used": False,
            "hst_f160w_used": False,
            "replacement_strength_scalar": None,
            "normalization_or_rescaling": False,
            "fit_or_tuning": False,
        },
        "rows": rows,
        "failures": failures,
        "checks": checks,
        "execution_gate_pass": passed,
        "interpretation_rule": "Report the five-cluster local benchmark-assisted end-to-end response without fitting. Do not call this an independent weak-lensing prediction because the benchmark kappa morphology constructs the source proxy.",
    }

    print(LAB_ID)
    print(f"status={status}")
    print(f"head_sha={state['head_sha']}")
    print("network_access_used=false")
    print("hst_f160w_used=false")
    print("source_mode=canonical_local_PBUF_benchmark_only")
    print("benchmark_assisted=true")
    print("independent_source_prediction=false")
    print("replacement_strength_scalar=none")
    print("fit_or_tuning=false")
    print("pair_law=A_ij=0.03*Delta_u_fast+0.003*Delta_u_slow")
    print()
    print("CLUSTERS")
    for r in rows:
        print(
            f"cluster={r['cluster_id']} "
            f"native_over_unit_los={r['native_over_unit_los_rms']:.12g} "
            f"native_over_unit_angle={r['native_over_unit_final_angle_rms']:.12g} "
            f"native_benchmark_pearson={r['native_final_angular_rms_angle_mag_vs_benchmark_pearson']:.12g} "
            f"native_benchmark_spearman={r['native_final_angular_rms_angle_mag_vs_benchmark_spearman']:.12g} "
            f"unit_benchmark_pearson={r['unit_control_final_angular_rms_angle_mag_vs_benchmark_pearson']:.12g} "
            f"unit_benchmark_spearman={r['unit_control_final_angular_rms_angle_mag_vs_benchmark_spearman']:.12g}"
        )
    for f in failures:
        print(f"failure_cluster={f['cluster_id']} error={f['error']}")
    print()
    print("CHECKS")
    for k, v in checks.items():
        print(f"{k}={str(v).lower()}")
    print(f"execution_gate_pass={str(passed).lower()}")
    print(json.dumps(result, sort_keys=True, separators=(",", ":"), default=str))
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
