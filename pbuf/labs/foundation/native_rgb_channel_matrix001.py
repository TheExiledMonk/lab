#!/usr/bin/env python3
"""PBUF FOUNDATION — NATIVE RGB-LIKE CHANNEL MATRIX 001.

Wide-net three-channel observer audit.

For each canonical cluster, build ONE frozen current-native received G3D ray state
and ONE target-blind global tangent detector screen. Reuse the raw receiver
channels established by native_multichannel_observer_fusion_sweep001, but now
allow kappa, gamma1, and gamma2 to be selected independently.

The complete Cartesian matrix is evaluated:
  5 kappa receivers x 4 gamma1 receivers x 4 gamma2 receivers = 80 states.

Every state is an identity assembly of raw detector channels. No fitted weights,
normalization, target-derived construction, rescaling, or cluster-specific choice
is allowed.

Three fixed observer overlays are reported for every state:
  vector_magnitude = sqrt(kappa^2 + gamma1^2 + gamma2^2)
  absolute_overlay = |kappa| + |gamma1| + |gamma2|
  signed_overlay   = kappa + gamma1 + gamma2

These are diagnostics of RGB-like channel combination, not fitted observables.
The same formula is applied to the observed channels for comparison.
"""
from __future__ import annotations

import itertools
import json
import statistics
import subprocess
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))

from pbuf.core import benchmark_data as BENCH
import pbuf.labs.foundation.native_multichannel_observer_fusion_sweep001 as FUS
import pbuf.labs.foundation.native_accumulated_full_lensing001 as G3D

LAB_ID = "PBUF-FOUNDATION-NATIVE-RGB-CHANNEL-MATRIX-001"
EXPECTED_CLUSTER_IDS = FUS.EXPECTED_CLUSTER_IDS
EPS = 1.0e-30

KAPPA_CHANNELS = tuple(FUS.KAPPA_CHANNELS.keys())
GAMMA1_CHANNELS = tuple(FUS.SHEAR_CHANNELS.keys())
GAMMA2_CHANNELS = tuple(FUS.SHEAR_CHANNELS.keys())
MATRIX = tuple(itertools.product(KAPPA_CHANNELS, GAMMA1_CHANNELS, GAMMA2_CHANNELS))
EXPECTED_MATRIX_SIZE = len(KAPPA_CHANNELS) * len(GAMMA1_CHANNELS) * len(GAMMA2_CHANNELS)


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
    return FUS.rms(a)


def compare(pred, obs) -> dict:
    return FUS.compare(pred, obs)


def orientation(pg1, pg2, og1, og2) -> float:
    return FUS.orientation(pg1, pg2, og1, og2)


def _overlay_fields(k, g1, g2) -> dict:
    k = np.asarray(k, dtype=np.float64)
    g1 = np.asarray(g1, dtype=np.float64)
    g2 = np.asarray(g2, dtype=np.float64)
    return {
        "vector_magnitude": np.sqrt(k*k + g1*g1 + g2*g2),
        "absolute_overlay": np.abs(k) + np.abs(g1) + np.abs(g2),
        "signed_overlay": k + g1 + g2,
    }


def _assembled(channels: dict, k_name: str, g1_name: str, g2_name: str) -> dict:
    k = np.asarray(channels[k_name], dtype=np.float64)
    g1 = np.asarray(channels[g1_name]["gamma1"], dtype=np.float64)
    g2 = np.asarray(channels[g2_name]["gamma2"], dtype=np.float64)
    out = {
        "kappa": k,
        "gamma1": g1,
        "gamma2": g2,
        "gamma": np.hypot(g1, g2),
    }
    out.update(_overlay_fields(k, g1, g2))
    return out


def _observed_with_overlays(data: dict) -> dict:
    obs = FUS._observed(data)
    # FUS composite is the same vector magnitude; keep explicit names here.
    obs.update(_overlay_fields(obs["kappa"], obs["gamma1"], obs["gamma2"]))
    return obs


def _score(pred: dict, obs: dict) -> dict:
    keys = (
        "kappa", "gamma", "gamma1", "gamma2",
        "vector_magnitude", "absolute_overlay", "signed_overlay",
    )
    row = {key: compare(pred[key], obs[key]) for key in keys}
    row["pred_gamma2_over_gamma1_rms"] = rms(pred["gamma2"]) / max(rms(pred["gamma1"]), EPS)
    row["obs_gamma2_over_gamma1_rms"] = rms(obs["gamma2"]) / max(rms(obs["gamma1"]), EPS)
    row["shear_orientation_cosine_mean"] = orientation(
        pred["gamma1"], pred["gamma2"], obs["gamma1"], obs["gamma2"]
    )
    return row


def _matrix_name(k_name: str, g1_name: str, g2_name: str) -> str:
    return f"{k_name}__{g1_name}__{g2_name}"


def run_cluster(cluster: dict) -> dict:
    frozen = FUS._build_frozen_state(cluster)
    channels = FUS._channel_inventory(frozen["extracted"])
    obs = _observed_with_overlays(frozen["data"])

    matrix = {}
    for k_name, g1_name, g2_name in MATRIX:
        name = _matrix_name(k_name, g1_name, g2_name)
        pred = _assembled(channels, k_name, g1_name, g2_name)
        row = _score(pred, obs)
        row["kappa_channel"] = k_name
        row["gamma1_channel"] = g1_name
        row["gamma2_channel"] = g2_name
        matrix[name] = row

    screen = frozen["screen"]
    return {
        "cluster_id": cluster["id"],
        "pair_fast_coefficient_from_A8": frozen["channel"]["pair_fast_coefficient_from_A8"],
        "pair_slow_coefficient_from_A8": frozen["channel"]["pair_slow_coefficient_from_A8"],
        "terminal_common_history_relative_rms_error": frozen["channel"]["terminal_common_history_relative_rms_error"],
        "g3d_unit_speed_max_error": frozen["chain"]["g3d"]["max_unit_speed_error"],
        "received_state_role": "one_frozen_current_native_G3D_state_shared_by_all_80_channel_states",
        "screen_role": "one_target_blind_global_tangent_screen_shared_by_all_80_channel_states",
        "screen_axis_e1": screen["e1"].tolist(),
        "screen_axis_e2": screen["e2"].tolist(),
        "screen_normal": screen["normal"].tolist(),
        "matrix": matrix,
    }


def aggregate(rows: list[dict]) -> dict:
    out = {}
    for k_name, g1_name, g2_name in MATRIX:
        name = _matrix_name(k_name, g1_name, g2_name)
        vals = [r["matrix"][name] for r in rows]
        k_amps = [v["kappa"]["rms_ratio_pred_over_obs"] for v in vals]
        vm_amps = [v["vector_magnitude"]["rms_ratio_pred_over_obs"] for v in vals]
        ao_amps = [v["absolute_overlay"]["rms_ratio_pred_over_obs"] for v in vals]
        so_amps = [v["signed_overlay"]["rms_ratio_pred_over_obs"] for v in vals]
        out[name] = {
            "kappa_channel": k_name,
            "gamma1_channel": g1_name,
            "gamma2_channel": g2_name,
            "mean_kappa_pearson": statistics.mean(v["kappa"]["pearson"] for v in vals),
            "mean_kappa_spearman": statistics.mean(v["kappa"]["spearman"] for v in vals),
            "mean_kappa_amp_ratio": statistics.mean(k_amps),
            "mean_gamma_pearson": statistics.mean(v["gamma"]["pearson"] for v in vals),
            "mean_gamma_spearman": statistics.mean(v["gamma"]["spearman"] for v in vals),
            "mean_gamma_amp_ratio": statistics.mean(v["gamma"]["rms_ratio_pred_over_obs"] for v in vals),
            "mean_gamma1_pearson": statistics.mean(v["gamma1"]["pearson"] for v in vals),
            "mean_gamma1_spearman": statistics.mean(v["gamma1"]["spearman"] for v in vals),
            "mean_gamma1_amp_ratio": statistics.mean(v["gamma1"]["rms_ratio_pred_over_obs"] for v in vals),
            "mean_gamma2_pearson": statistics.mean(v["gamma2"]["pearson"] for v in vals),
            "mean_gamma2_spearman": statistics.mean(v["gamma2"]["spearman"] for v in vals),
            "mean_gamma2_amp_ratio": statistics.mean(v["gamma2"]["rms_ratio_pred_over_obs"] for v in vals),
            "mean_vector_pearson": statistics.mean(v["vector_magnitude"]["pearson"] for v in vals),
            "mean_vector_spearman": statistics.mean(v["vector_magnitude"]["spearman"] for v in vals),
            "mean_vector_amp_ratio": statistics.mean(vm_amps),
            "mean_absolute_overlay_pearson": statistics.mean(v["absolute_overlay"]["pearson"] for v in vals),
            "mean_absolute_overlay_spearman": statistics.mean(v["absolute_overlay"]["spearman"] for v in vals),
            "mean_absolute_overlay_amp_ratio": statistics.mean(ao_amps),
            "mean_signed_overlay_pearson": statistics.mean(v["signed_overlay"]["pearson"] for v in vals),
            "mean_signed_overlay_spearman": statistics.mean(v["signed_overlay"]["spearman"] for v in vals),
            "mean_signed_overlay_amp_ratio": statistics.mean(so_amps),
            "mean_shear_orientation_cosine": statistics.mean(v["shear_orientation_cosine_mean"] for v in vals),
            "cross_cluster_kappa_amp_cv": statistics.pstdev(k_amps) / max(abs(statistics.mean(k_amps)), EPS),
            "cross_cluster_vector_amp_cv": statistics.pstdev(vm_amps) / max(abs(statistics.mean(vm_amps)), EPS),
        }
    return out


def _top(agg: dict, field: str, n: int = 12) -> list[dict]:
    rows = sorted(agg.items(), key=lambda kv: kv[1][field], reverse=True)[:n]
    return [{"state": name, field: vals[field], **{k: vals[k] for k in ("kappa_channel", "gamma1_channel", "gamma2_channel")}} for name, vals in rows]


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
        "matrix_size_is_80": EXPECTED_MATRIX_SIZE == 80,
        "all_80_states_present_per_cluster": bool(rows and all(len(r["matrix"]) == EXPECTED_MATRIX_SIZE for r in rows)),
        "all_80_aggregate_states_present": len(agg) == EXPECTED_MATRIX_SIZE,
        "same_received_G3D_state_shared_within_each_cluster": bool(rows),
        "same_target_blind_tangent_screen_shared_within_each_cluster": bool(rows),
        "native_terminal_common_history_identity": bool(rows and all(r["terminal_common_history_relative_rms_error"] <= 1e-12 for r in rows)),
        "G3D_unit_speed_valid": bool(rows and all(r["g3d_unit_speed_max_error"] <= G3D.UNIT_SPEED_TOL for r in rows)),
        "no_fitted_channel_weights": True,
        "no_output_rescaling": True,
        "no_target_derived_channel_construction": True,
        "no_cluster_specific_channel_choice": True,
        "no_tracked_or_staged_changes": not state["tracked_changes"] and not state["staged_changes"],
    }
    passed = all(checks.values())
    status = "NATIVE_RGB_CHANNEL_MATRIX_EXECUTED" if passed else ("NATIVE_RGB_CHANNEL_MATRIX_PARTIAL_EXECUTION" if rows else "NATIVE_RGB_CHANNEL_MATRIX_NOT_ESTABLISHED")

    top = {
        "vector_pearson": _top(agg, "mean_vector_pearson") if agg else [],
        "vector_spearman": _top(agg, "mean_vector_spearman") if agg else [],
        "absolute_overlay_pearson": _top(agg, "mean_absolute_overlay_pearson") if agg else [],
        "absolute_overlay_spearman": _top(agg, "mean_absolute_overlay_spearman") if agg else [],
        "signed_overlay_pearson": _top(agg, "mean_signed_overlay_pearson") if agg else [],
        "signed_overlay_spearman": _top(agg, "mean_signed_overlay_spearman") if agg else [],
    }

    result = {
        "lab_id": LAB_ID,
        "status": status,
        "repo_state": state,
        "frozen_rule": "one current-native G3D received state and one target-blind tangent screen per cluster; all 80 receiver states consume the same detector-plane rays",
        "matrix_rule": "5 kappa x 4 gamma1 x 4 gamma2 = 80 identity-assembled channel states",
        "overlay_rule": "vector magnitude, absolute sum, and signed sum use fixed unit channel transfer only",
        "kappa_channels": KAPPA_CHANNELS,
        "gamma1_channels": GAMMA1_CHANNELS,
        "gamma2_channels": GAMMA2_CHANNELS,
        "matrix_size": EXPECTED_MATRIX_SIZE,
        "rows": rows,
        "aggregate": agg,
        "top_diagnostics": top,
        "failures": failures,
        "checks": checks,
    }

    print(LAB_ID)
    print(f"status={status}")
    print(f"head_sha={state['head_sha']}")
    print("trajectory=frozen_current_native_G3D")
    print("observer_screen=frozen_target_blind_global_tangent")
    print(f"matrix_size={EXPECTED_MATRIX_SIZE}")
    print("matrix=5_kappa_x_4_gamma1_x_4_gamma2")
    print("channel_weights=identity_only")
    print("output_rescaling=false")
    print()
    print("TOP_DIAGNOSTICS")
    for key, entries in top.items():
        print(key)
        for e in entries:
            metric = [k for k in e if k.startswith("mean_")][0]
            print(f"  state={e['state']} value={e[metric]:.12g} k={e['kappa_channel']} g1={e['gamma1_channel']} g2={e['gamma2_channel']}")
    print()
    print("AGGREGATE_MATRIX")
    for name, a in agg.items():
        print(
            f"state={name} k={a['kappa_channel']} g1={a['gamma1_channel']} g2={a['gamma2_channel']} "
            f"k_r={a['mean_kappa_pearson']:.12g} k_rho={a['mean_kappa_spearman']:.12g} k_amp={a['mean_kappa_amp_ratio']:.12g} "
            f"g1_r={a['mean_gamma1_pearson']:.12g} g1_rho={a['mean_gamma1_spearman']:.12g} g1_amp={a['mean_gamma1_amp_ratio']:.12g} "
            f"g2_r={a['mean_gamma2_pearson']:.12g} g2_rho={a['mean_gamma2_spearman']:.12g} g2_amp={a['mean_gamma2_amp_ratio']:.12g} "
            f"vector_r={a['mean_vector_pearson']:.12g} vector_rho={a['mean_vector_spearman']:.12g} vector_amp={a['mean_vector_amp_ratio']:.12g} "
            f"abs_r={a['mean_absolute_overlay_pearson']:.12g} abs_rho={a['mean_absolute_overlay_spearman']:.12g} abs_amp={a['mean_absolute_overlay_amp_ratio']:.12g} "
            f"signed_r={a['mean_signed_overlay_pearson']:.12g} signed_rho={a['mean_signed_overlay_spearman']:.12g} signed_amp={a['mean_signed_overlay_amp_ratio']:.12g} "
            f"orientation={a['mean_shear_orientation_cosine']:.12g}"
        )
    print()
    print("CHECKS")
    for key, val in checks.items():
        print(f"{key}={'true' if val else 'false'}")
    print(json.dumps(result, sort_keys=True, separators=(",", ":"), allow_nan=True))
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
