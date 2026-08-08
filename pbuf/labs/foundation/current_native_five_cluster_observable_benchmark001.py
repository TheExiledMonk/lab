#!/usr/bin/env python3
"""PBUF FOUNDATION — CURRENT NATIVE FIVE-CLUSTER OBSERVABLE BENCHMARK 001.

Measure the current frozen native propagation lane on the established five-cluster
local benchmark using the canonical loader already used by the project.

Benchmark source construction is unchanged from the existing local benchmark labs:
    local kappa -> construct_common_proxy -> construct_rho_3d

Current native response lane:
    rho3 -> zero-flux terminal fast/slow channels
         -> exact A8 pair law using coefficients derived from A8 constants
         -> terminal c_state geometry
         -> PM1/PS2 -> M10 -> LOS -> photon propagation -> Jacobian
         -> predicted kappa/gamma1/gamma2

Observed local kappa/gamma/gamma1/gamma2 products are loaded only through
pbuf.core.benchmark_data. No HST path, legacy strength scalar, unit-control lane,
replacement scale, inferred transfer coefficient, output normalization, fitting,
or cluster-specific tuning is used.

This is a current-model benchmark using the long-standing kappa-derived benchmark
source construction. It is NOT labeled an independent matter-source prediction.
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
from weak_lensing_observation001 import resample_to_grid
from pbuf.core import benchmark_data as BENCH
from pbuf.core import observable_extraction as M16
from pbuf.models import a8_state as A8
import pbuf.labs.foundation.m10_coverage_25pct_science001 as BASE
import pbuf.labs.foundation.interface_to_interface_survivor_sweep001 as S92
import pbuf.labs.foundation.native_channel_transfer_closure_sweep001 as S93
import pbuf.labs.foundation.primary_candidate_science_rerun_m10_001 as OBS

LAB_ID = "PBUF-FOUNDATION-CURRENT-NATIVE-FIVE-CLUSTER-OBSERVABLE-BENCHMARK-001"
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
        raise RuntimeError(f"shape mismatch {x.shape} vs {y.shape}")
    m = np.isfinite(x) & np.isfinite(y)
    n = int(np.count_nonzero(m))
    if n < 2:
        return float("nan"), float("nan"), n
    return float(M16.safe_pearson(x[m], y[m])), float(M16.safe_spearman(x[m], y[m])), n


def local_cluster(cluster: dict) -> dict:
    kappa = BENCH.load_kappa(cluster)
    gamma = BENCH.load_gamma(cluster)
    gamma1 = BENCH.load_gamma1(cluster)
    gamma2 = BENCH.load_gamma2(cluster)
    rho2 = BASE.construct_common_proxy(kappa, bins=BASE.OBS_BINS, extent=BASE.CFG["extent"])
    rho3 = construct_rho_3d(rho2, BASE.NZ, profile=BASE.PROFILE)
    return {
        "kappa": np.asarray(kappa, dtype=np.float64),
        "gamma": np.asarray(gamma, dtype=np.float64),
        "gamma1": np.asarray(gamma1, dtype=np.float64),
        "gamma2": np.asarray(gamma2, dtype=np.float64),
        "rho2": np.asarray(rho2, dtype=np.float64),
        "rho3": np.asarray(rho3, dtype=np.float64),
    }


def current_native_m10(rho3: np.ndarray):
    channels = S93.native_terminal_channels(rho3)
    c = np.asarray(channels["c"], dtype=np.float64)
    uf = np.asarray(channels["u_fast"], dtype=np.float64)
    us = np.asarray(channels["u_slow"], dtype=np.float64)

    coef_fast = float(A8.A8_INIT_DT * A8.A8_INIT_OMEGA * A8.A8_INIT_K)
    coef_slow = float(A8.A8_INIT_DT * A8.A8_INIT_SLOW_TIMESCALE)

    fast = S93.scale_amplitudes(S92.positive_bond_amplitudes(uf), coef_fast)
    slow = S93.scale_amplitudes(S92.positive_bond_amplitudes(us), coef_slow)
    amps = S93.combine_amplitudes(fast, slow)
    m10 = S92.m10_from_amplitudes(amps, c)

    identity = S92.rms(np.asarray(channels["history_c"]) - c) / max(S92.rms(c), EPS)
    return m10, {
        "pair_fast_coefficient_from_A8": coef_fast,
        "pair_slow_coefficient_from_A8": coef_slow,
        "terminal_common_history_relative_rms_error": float(identity),
    }


def compare_field(pred, obs) -> dict:
    p, s, n = corr(pred, obs)
    pr = rms(pred)
    orms = rms(obs)
    return {
        "pearson": p,
        "spearman": s,
        "count": n,
        "pred_rms": pr,
        "obs_rms": orms,
        "rms_ratio_pred_over_obs": pr / max(orms, EPS),
    }


def run_cluster(cluster: dict) -> dict:
    data = local_cluster(cluster)
    m10, channel = current_native_m10(data["rho3"])

    # Reuse the existing frozen ray/Jacobian observable extraction. The only
    # reference passed here is the already-loaded canonical local kappa target;
    # it is used by package_lensing_observables for reporting, never to alter rays.
    ray = OBS._ray_and_observable(cluster["id"], m10, {"kappa": data["kappa"]})
    pred_k = np.asarray(ray["observable"]["kappa"], dtype=np.float64)
    pred_g1 = np.asarray(ray["observable"]["gamma1"], dtype=np.float64)
    pred_g2 = np.asarray(ray["observable"]["gamma2"], dtype=np.float64)
    pred_g = np.hypot(pred_g1, pred_g2)

    bins = int(BASE.CFG["bins"])
    extent = float(BASE.CFG["extent"])
    obs_k = resample_to_grid(data["kappa"], bins, extent)
    obs_g = resample_to_grid(data["gamma"], bins, extent)
    obs_g1 = resample_to_grid(data["gamma1"], bins, extent)
    obs_g2 = resample_to_grid(data["gamma2"], bins, extent)

    return {
        "cluster_id": cluster["id"],
        "source_mode": "canonical_local_kappa_to_existing_common_proxy_to_rho3",
        "benchmark_assisted_source": True,
        "legacy_strength_used": False,
        "unit_control_lane_used": False,
        "replacement_scale_used": False,
        "fit_or_tuning": False,
        **channel,
        "kappa": compare_field(pred_k, obs_k),
        "gamma": compare_field(pred_g, obs_g),
        "gamma1": compare_field(pred_g1, obs_g1),
        "gamma2": compare_field(pred_g2, obs_g2),
        "ray_metrics": ray["metrics"],
    }


def main() -> int:
    state = repo_state()
    clusters = list(BENCH.clusters())
    ids = tuple(c["id"] for c in clusters)
    inventory = BENCH.inventory()

    rows, failures = [], []
    ready = bool(
        ids == EXPECTED_CLUSTER_IDS
        and len(inventory) == 5
        and all(row["all_products_exist"] for row in inventory)
    )

    if ready:
        for cluster in clusters:
            try:
                rows.append(run_cluster(cluster))
            except Exception as exc:
                failures.append({"cluster_id": cluster["id"], "error": f"{type(exc).__name__}: {exc}"})

    finite_metrics = bool(rows and all(
        math.isfinite(float(row[obs][metric]))
        for row in rows
        for obs in ("kappa", "gamma", "gamma1", "gamma2")
        for metric in ("pearson", "spearman", "rms_ratio_pred_over_obs")
    ))

    checks = {
        "canonical_five_cluster_inventory": ids == EXPECTED_CLUSTER_IDS,
        "all_four_local_lensing_products_present_for_all_five": ready,
        "all_five_clusters_completed": len(rows) == 5 and not failures,
        "all_primary_metrics_finite": finite_metrics,
        "terminal_common_history_identity": bool(rows and all(row["terminal_common_history_relative_rms_error"] <= 1.0e-12 for row in rows)),
        "legacy_strength_used_false": bool(rows and all(not row["legacy_strength_used"] for row in rows)),
        "unit_control_lane_used_false": bool(rows and all(not row["unit_control_lane_used"] for row in rows)),
        "replacement_scale_used_false": bool(rows and all(not row["replacement_scale_used"] for row in rows)),
        "fit_or_tuning_false": bool(rows and all(not row["fit_or_tuning"] for row in rows)),
        "no_tracked_or_staged_changes": not state["tracked_changes"] and not state["staged_changes"],
    }
    passed = bool(all(checks.values()))
    status = (
        "CURRENT_NATIVE_FIVE_CLUSTER_OBSERVABLE_BENCHMARK_EXECUTED" if passed
        else ("CURRENT_NATIVE_FIVE_CLUSTER_OBSERVABLE_BENCHMARK_PARTIAL_EXECUTION" if rows
              else "CURRENT_NATIVE_FIVE_CLUSTER_OBSERVABLE_BENCHMARK_NOT_ESTABLISHED")
    )

    result = {
        "lab_id": LAB_ID,
        "status": status,
        "repo_state": state,
        "source_rule": "reuse established canonical local benchmark loader and existing kappa->proxy->rho3 construction",
        "model_rule": "current native terminal fast+slow transfer; A8 coefficients derived from model constants; no legacy strength/unit lane/replacement scale/fit",
        "interpretation_rule": "This measures current-model performance on the established benchmark. Because kappa constructs the benchmark source proxy, do not label it an independent matter-source prediction. Gamma/gamma1/gamma2 comparisons are separately reported as observable targets.",
        "rows": rows,
        "failures": failures,
        "checks": checks,
    }

    print(LAB_ID)
    print(f"status={status}")
    print(f"head_sha={state['head_sha']}")
    print("source_loader=pbuf.core.benchmark_data")
    print("source_construction=existing_kappa_to_common_proxy_to_rho3")
    print("legacy_strength_used=false")
    print("unit_control_lane_used=false")
    print("replacement_scale_used=false")
    print("fit_or_tuning=false")
    if rows:
        print(f"pair_fast_coefficient_from_A8={rows[0]['pair_fast_coefficient_from_A8']:.17g}")
        print(f"pair_slow_coefficient_from_A8={rows[0]['pair_slow_coefficient_from_A8']:.17g}")
    print()
    print("CLUSTERS")
    for row in rows:
        print(
            f"cluster={row['cluster_id']} "
            f"kappa_r={row['kappa']['pearson']:.12g} kappa_rho={row['kappa']['spearman']:.12g} kappa_amp={row['kappa']['rms_ratio_pred_over_obs']:.12g} "
            f"gamma_r={row['gamma']['pearson']:.12g} gamma_rho={row['gamma']['spearman']:.12g} gamma_amp={row['gamma']['rms_ratio_pred_over_obs']:.12g} "
            f"gamma1_r={row['gamma1']['pearson']:.12g} gamma1_rho={row['gamma1']['spearman']:.12g} gamma1_amp={row['gamma1']['rms_ratio_pred_over_obs']:.12g} "
            f"gamma2_r={row['gamma2']['pearson']:.12g} gamma2_rho={row['gamma2']['spearman']:.12g} gamma2_amp={row['gamma2']['rms_ratio_pred_over_obs']:.12g}"
        )
    for failure in failures:
        print(f"failure_cluster={failure['cluster_id']} error={failure['error']}")
    print()
    print("CHECKS")
    for key, value in checks.items():
        print(f"{key}={str(value).lower()}")
    print(json.dumps(result, sort_keys=True, separators=(",", ":"), default=str))
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
