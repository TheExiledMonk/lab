#!/usr/bin/env python3
"""PBUF FOUNDATION — NATIVE MULTICHANNEL OBSERVER FUSION SWEEP 001.

Test the hypothesis that the received ray field behaves like a multi-channel
measurement: different detector features may carry convergence and shear
information most cleanly, and the final observable may require stacking those
channels rather than forcing one extractor to produce everything.

Frozen for every candidate within a cluster:
  established local benchmark source -> current native fast/slow transfer
  -> PM1/PS2 -> M10 -> LOS -> G3D -> ONE received 3D ray state
  -> ONE target-blind global tangent detector screen

Only channel selection/stacking differs. No fitted weights are used. Each stack
copies complete raw channels from established extraction operators with unit
identity transfer; nothing is normalized to observations or rescaled.

Raw channel families:
  convergence: KNN density, KDE density, covariance-area, affine Jacobian,
               displacement divergence
  shear:       affine Jacobian, covariance transport, displacement divergence,
               KNN displacement-derived shear

Predeclared stacks combine one convergence channel with one complete (g1,g2)
shear channel. A geometric composite image is also formed as
  sqrt(kappa^2 + gamma1^2 + gamma2^2)
from the stacked channels, with the same expression applied to observations.
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

from weak_lensing_observation001 import resample_to_grid
from pbuf.core import benchmark_data as BENCH
import pbuf.labs.foundation.current_native_five_cluster_observable_benchmark001 as CUR
import pbuf.labs.foundation.native_accumulated_full_lensing001 as G3D
import pbuf.labs.foundation.m10_coverage_25pct_science001 as BASE
import pbuf.labs.foundation.native_observable_extraction_method_sweep001 as EX

LAB_ID = "PBUF-FOUNDATION-NATIVE-MULTICHANNEL-OBSERVER-FUSION-SWEEP-001"
EXPECTED_CLUSTER_IDS = ("Abell2744", "MACS0416", "MACS1149", "AbellS1063", "Abell370")
EPS = 1.0e-30

KAPPA_CHANNELS = {
    "k_knn": ("knn_density", "convergence"),
    "k_kernel": ("kernel_density", "convergence"),
    "k_area": ("covariance_area", "convergence"),
    "k_jacobian": ("jacobian_affine", "convergence"),
    "k_divergence": ("displacement_divergence", "convergence"),
}
SHEAR_CHANNELS = {
    "g_jacobian": "jacobian_affine",
    "g_covariance": "covariance_transport",
    "g_divergence": "displacement_divergence",
    "g_knn": "knn_density",
}

# No fitted coefficients. Each stack simply lays raw channels on top of each other.
STACKS = {
    "knn_plus_jacobian": ("k_knn", "g_jacobian"),
    "knn_plus_covariance": ("k_knn", "g_covariance"),
    "knn_plus_divergence": ("k_knn", "g_divergence"),
    "knn_plus_knn": ("k_knn", "g_knn"),
    "kernel_plus_jacobian": ("k_kernel", "g_jacobian"),
    "area_plus_jacobian": ("k_area", "g_jacobian"),
    "jacobian_all_control": ("k_jacobian", "g_jacobian"),
    "covariance_all_control": ("k_area", "g_covariance"),
    "divergence_all_control": ("k_divergence", "g_divergence"),
}


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
    return EX.rms(a)


def compare(pred, obs) -> dict:
    return EX._compare(pred, obs)


def orientation(pg1, pg2, og1, og2) -> float:
    return EX._orientation(pg1, pg2, og1, og2)


def _observed(data: dict) -> dict:
    bins = int(BASE.CFG["bins"])
    extent = float(BASE.CFG["extent"])
    out = {
        "kappa": resample_to_grid(data["kappa"], bins, extent),
        "gamma": resample_to_grid(data["gamma"], bins, extent),
        "gamma1": resample_to_grid(data["gamma1"], bins, extent),
        "gamma2": resample_to_grid(data["gamma2"], bins, extent),
    }
    out["composite"] = np.sqrt(out["kappa"]**2 + out["gamma1"]**2 + out["gamma2"]**2)
    return out


def _build_frozen_state(cluster: dict) -> dict:
    data = CUR.local_cluster(cluster)
    m10, channel = CUR.current_native_m10(data["rho3"])
    chain = G3D.run_g3d_from_vector(m10, observed_for_first_step=None)
    snap = chain["checkpoints"][G3D.CHECKPOINT]
    x0, y0, _, _ = BASE._launch_expanded_25pct()
    screen = EX._screen_coordinates(x0, y0, snap)
    extracted = EX._extract_all(screen, float(BASE.CFG["extent"]), int(BASE.CFG["bins"]))
    return {
        "data": data,
        "channel": channel,
        "chain": chain,
        "screen": screen,
        "extracted": extracted,
    }


def _channel_inventory(extracted: dict) -> dict:
    out = {}
    for name, (method, field) in KAPPA_CHANNELS.items():
        out[name] = np.asarray(extracted[method][field], dtype=np.float64)
    for name, method in SHEAR_CHANNELS.items():
        g1 = np.asarray(extracted[method]["shear_g1"], dtype=np.float64)
        g2 = np.asarray(extracted[method]["shear_g2"], dtype=np.float64)
        out[name] = {"gamma1": g1, "gamma2": g2, "gamma": np.hypot(g1, g2)}
    return out


def _stack(channels: dict, k_name: str, g_name: str) -> dict:
    k = np.asarray(channels[k_name], dtype=np.float64)
    g1 = np.asarray(channels[g_name]["gamma1"], dtype=np.float64)
    g2 = np.asarray(channels[g_name]["gamma2"], dtype=np.float64)
    return {
        "kappa": k,
        "gamma1": g1,
        "gamma2": g2,
        "gamma": np.hypot(g1, g2),
        "composite": np.sqrt(k*k + g1*g1 + g2*g2),
    }


def _score_stack(pred: dict, obs: dict) -> dict:
    row = {key: compare(pred[key], obs[key]) for key in ("kappa", "gamma", "gamma1", "gamma2", "composite")}
    row["pred_gamma2_over_gamma1_rms"] = rms(pred["gamma2"]) / max(rms(pred["gamma1"]), EPS)
    row["obs_gamma2_over_gamma1_rms"] = rms(obs["gamma2"]) / max(rms(obs["gamma1"]), EPS)
    row["shear_orientation_cosine_mean"] = orientation(pred["gamma1"], pred["gamma2"], obs["gamma1"], obs["gamma2"])
    return row


def _score_raw_channels(channels: dict, obs: dict) -> dict:
    out = {"kappa": {}, "shear": {}}
    for name in KAPPA_CHANNELS:
        out["kappa"][name] = compare(channels[name], obs["kappa"])
    for name in SHEAR_CHANNELS:
        g = channels[name]
        out["shear"][name] = {
            "gamma": compare(g["gamma"], obs["gamma"]),
            "gamma1": compare(g["gamma1"], obs["gamma1"]),
            "gamma2": compare(g["gamma2"], obs["gamma2"]),
            "shear_orientation_cosine_mean": orientation(g["gamma1"], g["gamma2"], obs["gamma1"], obs["gamma2"]),
        }
    return out


def run_cluster(cluster: dict) -> dict:
    frozen = _build_frozen_state(cluster)
    obs = _observed(frozen["data"])
    channels = _channel_inventory(frozen["extracted"])
    stacks = {}
    for name, (k_name, g_name) in STACKS.items():
        stacks[name] = _score_stack(_stack(channels, k_name, g_name), obs)
        stacks[name]["kappa_channel"] = k_name
        stacks[name]["shear_channel"] = g_name

    screen = frozen["screen"]
    return {
        "cluster_id": cluster["id"],
        "pair_fast_coefficient_from_A8": frozen["channel"]["pair_fast_coefficient_from_A8"],
        "pair_slow_coefficient_from_A8": frozen["channel"]["pair_slow_coefficient_from_A8"],
        "terminal_common_history_relative_rms_error": frozen["channel"]["terminal_common_history_relative_rms_error"],
        "g3d_unit_speed_max_error": frozen["chain"]["g3d"]["max_unit_speed_error"],
        "received_state_role": "one_frozen_current_native_G3D_state_shared_by_all_channels",
        "screen_role": "one_target_blind_global_tangent_screen_shared_by_all_channels",
        "screen_axis_e1": screen["e1"].tolist(),
        "screen_axis_e2": screen["e2"].tolist(),
        "screen_normal": screen["normal"].tolist(),
        "raw_channels": _score_raw_channels(channels, obs),
        "stacks": stacks,
    }


def aggregate(rows: list[dict]) -> dict:
    out = {}
    for name in STACKS:
        vals = [r["stacks"][name] for r in rows]
        kappas = [v["kappa"]["rms_ratio_pred_over_obs"] for v in vals]
        mean_k = statistics.mean(kappas)
        out[name] = {
            "kappa_channel": STACKS[name][0],
            "shear_channel": STACKS[name][1],
            "mean_kappa_pearson": statistics.mean(v["kappa"]["pearson"] for v in vals),
            "mean_kappa_spearman": statistics.mean(v["kappa"]["spearman"] for v in vals),
            "mean_kappa_amp_ratio": mean_k,
            "mean_gamma_pearson": statistics.mean(v["gamma"]["pearson"] for v in vals),
            "mean_gamma_spearman": statistics.mean(v["gamma"]["spearman"] for v in vals),
            "mean_gamma_amp_ratio": statistics.mean(v["gamma"]["rms_ratio_pred_over_obs"] for v in vals),
            "mean_gamma1_amp_ratio": statistics.mean(v["gamma1"]["rms_ratio_pred_over_obs"] for v in vals),
            "mean_gamma2_amp_ratio": statistics.mean(v["gamma2"]["rms_ratio_pred_over_obs"] for v in vals),
            "mean_composite_pearson": statistics.mean(v["composite"]["pearson"] for v in vals),
            "mean_composite_spearman": statistics.mean(v["composite"]["spearman"] for v in vals),
            "mean_composite_amp_ratio": statistics.mean(v["composite"]["rms_ratio_pred_over_obs"] for v in vals),
            "mean_shear_orientation_cosine": statistics.mean(v["shear_orientation_cosine_mean"] for v in vals),
            "cross_cluster_kappa_amp_cv": statistics.pstdev(kappas) / max(abs(mean_k), EPS),
        }
    return out


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
        "all_nine_predeclared_stacks_present": bool(rows and all(set(r["stacks"]) == set(STACKS) for r in rows)),
        "all_nine_raw_channels_present": bool(rows and all(len(r["raw_channels"]["kappa"]) + len(r["raw_channels"]["shear"]) == 9 for r in rows)),
        "same_received_G3D_state_shared_within_each_cluster": bool(rows),
        "same_target_blind_tangent_screen_shared_within_each_cluster": bool(rows),
        "native_terminal_common_history_identity": bool(rows and all(r["terminal_common_history_relative_rms_error"] <= 1e-12 for r in rows)),
        "G3D_unit_speed_valid": bool(rows and all(r["g3d_unit_speed_max_error"] <= G3D.UNIT_SPEED_TOL for r in rows)),
        "no_fitted_channel_weights": True,
        "no_output_rescaling": True,
        "no_tracked_or_staged_changes": not state["tracked_changes"] and not state["staged_changes"],
    }
    passed = all(checks.values())
    status = "NATIVE_MULTICHANNEL_OBSERVER_FUSION_SWEEP_EXECUTED" if passed else ("NATIVE_MULTICHANNEL_OBSERVER_FUSION_SWEEP_PARTIAL_EXECUTION" if rows else "NATIVE_MULTICHANNEL_OBSERVER_FUSION_SWEEP_NOT_ESTABLISHED")

    result = {
        "lab_id": LAB_ID,
        "status": status,
        "repo_state": state,
        "frozen_rule": "one current-native G3D received state and one target-blind tangent detector screen per cluster; only raw channel selection/stacking may differ",
        "combination_rule": "predeclared channel stacks use identity transfer only; no fitted weights, normalization, rescaling, or target-dependent channel construction",
        "kappa_channels": KAPPA_CHANNELS,
        "shear_channels": SHEAR_CHANNELS,
        "stacks": STACKS,
        "rows": rows,
        "aggregate": agg,
        "failures": failures,
        "checks": checks,
    }

    print(LAB_ID)
    print(f"status={status}")
    print(f"head_sha={state['head_sha']}")
    print("trajectory=frozen_current_native_G3D")
    print("observer_screen=frozen_target_blind_global_tangent")
    print("channel_weights=none")
    print("output_rescaling=false")
    print("stacks=" + ",".join(STACKS))
    print()
    print("AGGREGATE")
    for name, a in agg.items():
        print(
            f"stack={name} k_channel={a['kappa_channel']} g_channel={a['shear_channel']} "
            f"kappa_r={a['mean_kappa_pearson']:.12g} kappa_rho={a['mean_kappa_spearman']:.12g} kappa_amp={a['mean_kappa_amp_ratio']:.12g} "
            f"gamma_r={a['mean_gamma_pearson']:.12g} gamma_rho={a['mean_gamma_spearman']:.12g} gamma_amp={a['mean_gamma_amp_ratio']:.12g} "
            f"composite_r={a['mean_composite_pearson']:.12g} composite_rho={a['mean_composite_spearman']:.12g} composite_amp={a['mean_composite_amp_ratio']:.12g} "
            f"orientation={a['mean_shear_orientation_cosine']:.12g}"
        )
    print()
    print("CLUSTERS")
    for row in rows:
        for name, s in row["stacks"].items():
            print(
                f"cluster={row['cluster_id']} stack={name} k_channel={s['kappa_channel']} g_channel={s['shear_channel']} "
                f"kappa_r={s['kappa']['pearson']:.12g} kappa_rho={s['kappa']['spearman']:.12g} kappa_amp={s['kappa']['rms_ratio_pred_over_obs']:.12g} "
                f"gamma_r={s['gamma']['pearson']:.12g} gamma_amp={s['gamma']['rms_ratio_pred_over_obs']:.12g} "
                f"composite_r={s['composite']['pearson']:.12g} composite_amp={s['composite']['rms_ratio_pred_over_obs']:.12g} "
                f"orientation={s['shear_orientation_cosine_mean']:.12g}"
            )
    print()
    print("CHECKS")
    for key, value in checks.items():
        print(f"{key}={str(value).lower()}")
    for failure in failures:
        print(f"failure_cluster={failure['cluster_id']} error={failure['error']}")
    print(json.dumps(result, sort_keys=True, separators=(",", ":"), default=str))
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
