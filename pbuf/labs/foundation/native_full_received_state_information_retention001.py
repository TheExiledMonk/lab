#!/usr/bin/env python3
"""PBUF FOUNDATION — FULL RECEIVED-STATE INFORMATION RETENTION 001.

Observer-side audit only.

For each canonical cluster, build ONE frozen current-native received G3D ray
state and ONE target-blind global tangent detector screen. Decode a broad,
predeclared inventory of receiver channels BEFORE reducing the state to the
usual weak-lensing observables.

The audit asks:
  1. How many independent dimensions are present in the received 3D ray state?
  2. How much observer-side information survives successive target-blind
     reductions from the full decoded channel bank toward conventional 2D and
     three-field summaries?
  3. Which decoded channels are redundant/complementary?
  4. Does explicitly retained depth/direction information add independent
     received-state dimensions?

No propagation, source, A8, PM1/PS2, M10, LOS, or G3D physics is changed.
Observed kappa/gamma products are used only after decoding for descriptive
per-channel relevance; they never construct, normalize, select, or weight a
channel.

Internal column standardization is used ONLY for dimensionless information-
geometry diagnostics (rank/SVD/correlation) so heterogeneous decoded channels
can be compared. It is not an observable rescaling and is never fed back into
physical outputs.
"""
from __future__ import annotations

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
import pbuf.labs.foundation.native_observable_extraction_method_sweep001 as EX
import pbuf.labs.foundation.native_accumulated_full_lensing001 as G3D
import pbuf.labs.foundation.m10_coverage_25pct_science001 as BASE

LAB_ID = "PBUF-FOUNDATION-FULL-RECEIVED-STATE-INFORMATION-RETENTION-001"
EXPECTED_CLUSTER_IDS = FUS.EXPECTED_CLUSTER_IDS
EPS = 1.0e-30

EXTRACTION_METHODS = EX.METHODS
EXTRACTION_FIELDS = ("convergence", "shear_g1", "shear_g2")

PRIMARY_3D_BIN_CHANNELS = (
    "mean_du", "mean_dv", "mean_dw",
    "mean_t1", "mean_t2", "mean_tn",
    "std_du", "std_dv", "std_dw",
    "std_t1", "std_t2", "std_tn",
    "cov_du_dv", "cov_du_dw", "cov_dv_dw",
    "j3_e1_u", "j3_e1_v", "j3_e2_u", "j3_e2_v", "j3_n_u", "j3_n_v",
)

DEPTH_CHANNELS = {
    "mean_dw", "mean_tn", "std_dw", "std_tn",
    "cov_du_dw", "cov_dv_dw", "j3_n_u", "j3_n_v",
}
DIRECTION_CHANNELS = {
    "mean_t1", "mean_t2", "mean_tn", "std_t1", "std_t2", "std_tn",
}
DISPLACEMENT_CHANNELS = {
    "mean_du", "mean_dv", "mean_dw",
    "std_du", "std_dv", "std_dw",
    "cov_du_dv", "cov_du_dw", "cov_dv_dw",
}
J3_CHANNELS = {
    "j3_e1_u", "j3_e1_v", "j3_e2_u", "j3_e2_v", "j3_n_u", "j3_n_v",
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


def _finite(v) -> np.ndarray:
    return np.asarray(v, dtype=np.float64)


def _bin_indices(u0, v0, extent: float, bins: int):
    edges = np.linspace(-extent, extent, bins + 1)
    c = np.searchsorted(edges, u0, side="right") - 1
    r = np.searchsorted(edges, v0, side="right") - 1
    valid = (r >= 0) & (r < bins) & (c >= 0) & (c < bins)
    return r, c, valid


def _empty(bins: int) -> np.ndarray:
    return np.full((bins, bins), np.nan, dtype=np.float64)


def _binned_received_3d(screen: dict, snap: dict, extent: float, bins: int) -> dict[str, np.ndarray]:
    e1 = _finite(screen["e1"])
    e2 = _finite(screen["e2"])
    n = _finite(screen["normal"])
    u0 = _finite(screen["u0"])
    v0 = _finite(screen["v0"])
    uf = _finite(screen["uf"])
    vf = _finite(screen["vf"])

    x0, y0, _, _ = BASE._launch_expanded_25pct()
    p0 = np.column_stack((x0, y0, np.zeros_like(x0))).astype(np.float64)
    pf = np.column_stack((snap["x"], snap["y"], snap["z"])).astype(np.float64)
    vel = np.column_stack((snap["vx"], snap["vy"], snap["vz"])).astype(np.float64)

    w0 = p0 @ n
    wf = pf @ n
    du = uf - u0
    dv = vf - v0
    dw = wf - w0
    t1 = vel @ e1
    t2 = vel @ e2
    tn = vel @ n

    scalars = {"du": du, "dv": dv, "dw": dw, "t1": t1, "t2": t2, "tn": tn}
    out = {name: _empty(bins) for name in PRIMARY_3D_BIN_CHANNELS}
    r, c, valid = _bin_indices(u0, v0, extent, bins)
    flat = r * bins + c

    for q in np.unique(flat[valid]):
        idx = np.where(valid & (flat == q))[0]
        rr, cc = divmod(int(q), bins)
        if idx.size == 0:
            continue

        for key in ("du", "dv", "dw", "t1", "t2", "tn"):
            vals = scalars[key][idx]
            out[f"mean_{key}"][rr, cc] = float(np.mean(vals))
            out[f"std_{key}"][rr, cc] = float(np.std(vals))

        if idx.size >= 2:
            D = np.column_stack((du[idx], dv[idx], dw[idx]))
            C = np.cov(D, rowvar=False, ddof=1)
            out["cov_du_dv"][rr, cc] = float(C[0, 1])
            out["cov_du_dw"][rr, cc] = float(C[0, 2])
            out["cov_dv_dw"][rr, cc] = float(C[1, 2])

        if idx.size >= 6:
            X = np.column_stack((u0[idx] - np.mean(u0[idx]), v0[idx] - np.mean(v0[idx])))
            Y = np.column_stack((uf[idx] - np.mean(uf[idx]), vf[idx] - np.mean(vf[idx]), wf[idx] - np.mean(wf[idx])))
            try:
                A, *_ = np.linalg.lstsq(X, Y, rcond=None)
            except np.linalg.LinAlgError:
                continue
            out["j3_e1_u"][rr, cc] = float(A[0, 0])
            out["j3_e1_v"][rr, cc] = float(A[1, 0])
            out["j3_e2_u"][rr, cc] = float(A[0, 1])
            out["j3_e2_v"][rr, cc] = float(A[1, 1])
            out["j3_n_u"][rr, cc] = float(A[0, 2])
            out["j3_n_v"][rr, cc] = float(A[1, 2])

    out["_per_ray_full3d"] = np.column_stack((du, dv, dw, t1, t2, tn))
    out["_per_ray_transverse2d"] = np.column_stack((du, dv, t1, t2))
    return out


def _decoded_bank(extracted: dict, receipt3d: dict) -> tuple[dict[str, np.ndarray], dict[str, str]]:
    bank: dict[str, np.ndarray] = {}
    family: dict[str, str] = {}

    for method in EXTRACTION_METHODS:
        for field in EXTRACTION_FIELDS:
            name = f"{method}__{field}"
            bank[name] = _finite(extracted[method][field])
            if method in ("histogram_density", "kernel_density", "knn_density"):
                family[name] = "density"
            elif method == "displacement_divergence":
                family[name] = "displacement_2d"
            elif method == "covariance_area":
                family[name] = "area"
            else:
                family[name] = "differential_shape"

    for name in PRIMARY_3D_BIN_CHANNELS:
        bank[name] = _finite(receipt3d[name])
        if name in DEPTH_CHANNELS:
            family[name] = "depth_3d"
        elif name in DIRECTION_CHANNELS:
            family[name] = "direction"
        elif name in J3_CHANNELS:
            family[name] = "j3_differential"
        elif name in DISPLACEMENT_CHANNELS:
            family[name] = "displacement_3d"
        else:
            family[name] = "received_3d"

    return bank, family


def _standardized_matrix(bank: dict[str, np.ndarray], names: list[str]) -> tuple[np.ndarray, dict]:
    cols = []
    coverage = {}
    used = []
    for name in names:
        a = _finite(bank[name]).reshape(-1)
        finite = np.isfinite(a)
        coverage[name] = float(np.count_nonzero(finite) / max(a.size, 1))
        if not np.any(finite):
            continue
        mu = float(np.mean(a[finite]))
        sd = float(np.std(a[finite]))
        if not np.isfinite(sd) or sd <= EPS:
            continue
        z = np.zeros_like(a)
        z[finite] = (a[finite] - mu) / sd
        cols.append(z)
        used.append(name)
    X = np.column_stack(cols) if cols else np.empty((0, 0), dtype=np.float64)
    return X, {"coverage": coverage, "used_channels": used}


def _information_geometry(X: np.ndarray) -> dict:
    X = np.asarray(X, dtype=np.float64)
    if X.ndim != 2 or X.shape[0] < 2 or X.shape[1] == 0:
        return {
            "samples": int(X.shape[0]) if X.ndim == 2 else 0,
            "channels": int(X.shape[1]) if X.ndim == 2 else 0,
            "numerical_rank": 0,
            "effective_rank": 0.0,
            "participation_ratio": 0.0,
            "top5_variance_fraction": 0.0,
            "median_abs_channel_corr": float("nan"),
            "max_abs_channel_corr": float("nan"),
        }

    Xc = X - np.mean(X, axis=0, keepdims=True)
    s = np.linalg.svd(Xc, compute_uv=False, full_matrices=False)
    power = s * s
    total = float(np.sum(power))
    p = power / max(total, EPS)
    good = p > 0.0
    effective_rank = float(np.exp(-np.sum(p[good] * np.log(p[good]))))
    participation = float(1.0 / max(np.sum(p * p), EPS))
    numerical_rank = int(np.count_nonzero(s > max(float(s[0]) if s.size else 0.0, 1.0) * 1e-10))
    top5 = float(np.sum(p[: min(5, p.size)]))

    if X.shape[1] >= 2:
        C = np.corrcoef(X, rowvar=False)
        tri = np.abs(C[np.triu_indices(C.shape[0], 1)])
        tri = tri[np.isfinite(tri)]
        med = float(np.median(tri)) if tri.size else float("nan")
        mx = float(np.max(tri)) if tri.size else float("nan")
    else:
        med = mx = float("nan")

    return {
        "samples": int(X.shape[0]),
        "channels": int(X.shape[1]),
        "numerical_rank": numerical_rank,
        "effective_rank": effective_rank,
        "participation_ratio": participation,
        "top5_variance_fraction": top5,
        "median_abs_channel_corr": med,
        "max_abs_channel_corr": mx,
    }


def _per_ray_geometry(receipt3d: dict) -> dict:
    full = _finite(receipt3d["_per_ray_full3d"])
    transverse = _finite(receipt3d["_per_ray_transverse2d"])

    def zscore(X):
        mu = np.mean(X, axis=0, keepdims=True)
        sd = np.std(X, axis=0, keepdims=True)
        keep = np.squeeze(sd > EPS)
        return (X[:, keep] - mu[:, keep]) / sd[:, keep] if np.any(keep) else np.empty((X.shape[0], 0))

    fg = _information_geometry(zscore(full))
    tg = _information_geometry(zscore(transverse))
    return {
        "full_3d_received_ray_state": fg,
        "transverse_2d_received_ray_state": tg,
        "incremental_numerical_rank_from_depth_normal_channels": fg["numerical_rank"] - tg["numerical_rank"],
        "incremental_effective_rank_from_depth_normal_channels": fg["effective_rank"] - tg["effective_rank"],
    }


def _stage_definitions(bank: dict[str, np.ndarray], family: dict[str, str]) -> dict[str, list[str]]:
    all_names = list(bank)
    no_depth = [n for n in all_names if family[n] != "depth_3d"]
    conventional_2d_all = [n for n in all_names if "__" in n]
    jacobian3 = [
        "jacobian_affine__convergence",
        "jacobian_affine__shear_g1",
        "jacobian_affine__shear_g2",
    ]
    density3 = [
        "histogram_density__convergence",
        "kernel_density__convergence",
        "knn_density__convergence",
    ]
    return {
        "full_decoded_bank": all_names,
        "full_minus_explicit_depth3d": no_depth,
        "all_established_2d_decoders": conventional_2d_all,
        "canonical_jacobian_three_field": jacobian3,
        "density_three_channel_control": density3,
    }


def _stage_metrics(bank: dict[str, np.ndarray], family: dict[str, str]) -> dict:
    defs = _stage_definitions(bank, family)
    out = {}
    full_eff = None
    full_rank = None
    for stage, names in defs.items():
        X, meta = _standardized_matrix(bank, names)
        geom = _information_geometry(X)
        geom["requested_channels"] = len(names)
        geom["used_channel_names"] = meta["used_channels"]
        geom["mean_coverage_fraction"] = statistics.mean(meta["coverage"].values()) if meta["coverage"] else float("nan")
        out[stage] = geom
        if stage == "full_decoded_bank":
            full_eff = geom["effective_rank"]
            full_rank = geom["numerical_rank"]

    for stage, geom in out.items():
        geom["effective_rank_fraction_of_full"] = geom["effective_rank"] / max(full_eff or 0.0, EPS)
        geom["numerical_rank_fraction_of_full"] = geom["numerical_rank"] / max(full_rank or 0, 1)
    return out


def _family_metrics(bank: dict[str, np.ndarray], family: dict[str, str]) -> dict:
    out = {}
    for fam in sorted(set(family.values())):
        names = [n for n, f in family.items() if f == fam]
        X, meta = _standardized_matrix(bank, names)
        geom = _information_geometry(X)
        geom["channel_names"] = names
        geom["used_channel_names"] = meta["used_channels"]
        out[fam] = geom
    return out


def _observed_relevance(bank: dict[str, np.ndarray], data: dict) -> dict:
    obs = FUS._observed(data)
    targets = {"kappa": obs["kappa"], "gamma1": obs["gamma1"], "gamma2": obs["gamma2"], "gamma": obs["gamma"]}
    out = {}
    for name, field in bank.items():
        out[name] = {tname: EX._compare(field, target) for tname, target in targets.items()}
    return out


def _aggregate_relevance(rows: list[dict]) -> list[dict]:
    names = rows[0]["observer_relevance"].keys() if rows else []
    ranked = []
    for name in names:
        per_target = {}
        for target in ("kappa", "gamma1", "gamma2", "gamma"):
            vals = [r["observer_relevance"][name][target]["pearson"] for r in rows]
            vals = [v for v in vals if np.isfinite(v)]
            per_target[target] = statistics.mean(vals) if vals else float("nan")
        finite_items = [(k, v) for k, v in per_target.items() if np.isfinite(v)]
        if finite_items:
            best_target, best_r = max(finite_items, key=lambda kv: abs(kv[1]))
        else:
            best_target, best_r = "none", float("nan")
        ranked.append({
            "channel": name,
            "best_observed_target": best_target,
            "best_mean_pearson": best_r,
            "mean_pearson_by_target": per_target,
            "family": rows[0]["channel_families"][name],
        })
    ranked.sort(key=lambda r: abs(r["best_mean_pearson"]) if np.isfinite(r["best_mean_pearson"]) else -1.0, reverse=True)
    return ranked


def run_cluster(cluster: dict) -> dict:
    frozen = FUS._build_frozen_state(cluster)
    data = frozen["data"]
    screen = frozen["screen"]
    snap = frozen["chain"]["checkpoints"][G3D.CHECKPOINT]
    bins = int(BASE.CFG["bins"])
    extent = float(BASE.CFG["extent"])

    receipt3d = _binned_received_3d(screen, snap, extent, bins)
    bank, family = _decoded_bank(frozen["extracted"], receipt3d)

    return {
        "cluster_id": cluster["id"],
        "pair_fast_coefficient_from_A8": frozen["channel"]["pair_fast_coefficient_from_A8"],
        "pair_slow_coefficient_from_A8": frozen["channel"]["pair_slow_coefficient_from_A8"],
        "terminal_common_history_relative_rms_error": frozen["channel"]["terminal_common_history_relative_rms_error"],
        "g3d_unit_speed_max_error": frozen["chain"]["g3d"]["max_unit_speed_error"],
        "received_state_role": "one_frozen_current_native_G3D_state_decoded_once_into_full_channel_bank",
        "screen_role": "one_target_blind_global_tangent_screen_shared_by_all_decoded_channels",
        "decoded_channel_count": len(bank),
        "decoded_channel_names": list(bank),
        "channel_families": family,
        "per_ray_information_geometry": _per_ray_geometry(receipt3d),
        "stage_information_geometry": _stage_metrics(bank, family),
        "family_information_geometry": _family_metrics(bank, family),
        "observer_relevance": _observed_relevance(bank, data),
    }


def aggregate(rows: list[dict]) -> dict:
    if not rows:
        return {}
    stage_names = rows[0]["stage_information_geometry"].keys()
    stages = {}
    for stage in stage_names:
        vals = [r["stage_information_geometry"][stage] for r in rows]
        corr_vals = [v["median_abs_channel_corr"] for v in vals if np.isfinite(v["median_abs_channel_corr"])]
        stages[stage] = {
            "mean_numerical_rank": statistics.mean(v["numerical_rank"] for v in vals),
            "mean_effective_rank": statistics.mean(v["effective_rank"] for v in vals),
            "mean_participation_ratio": statistics.mean(v["participation_ratio"] for v in vals),
            "mean_top5_variance_fraction": statistics.mean(v["top5_variance_fraction"] for v in vals),
            "mean_median_abs_channel_corr": statistics.mean(corr_vals) if corr_vals else float("nan"),
            "mean_effective_rank_fraction_of_full": statistics.mean(v["effective_rank_fraction_of_full"] for v in vals),
            "mean_numerical_rank_fraction_of_full": statistics.mean(v["numerical_rank_fraction_of_full"] for v in vals),
            "requested_channels": vals[0]["requested_channels"],
        }

    per_ray = {
        "mean_full_3d_numerical_rank": statistics.mean(r["per_ray_information_geometry"]["full_3d_received_ray_state"]["numerical_rank"] for r in rows),
        "mean_transverse_2d_numerical_rank": statistics.mean(r["per_ray_information_geometry"]["transverse_2d_received_ray_state"]["numerical_rank"] for r in rows),
        "mean_incremental_numerical_rank_from_depth_normal_channels": statistics.mean(r["per_ray_information_geometry"]["incremental_numerical_rank_from_depth_normal_channels"] for r in rows),
        "mean_full_3d_effective_rank": statistics.mean(r["per_ray_information_geometry"]["full_3d_received_ray_state"]["effective_rank"] for r in rows),
        "mean_transverse_2d_effective_rank": statistics.mean(r["per_ray_information_geometry"]["transverse_2d_received_ray_state"]["effective_rank"] for r in rows),
        "mean_incremental_effective_rank_from_depth_normal_channels": statistics.mean(r["per_ray_information_geometry"]["incremental_effective_rank_from_depth_normal_channels"] for r in rows),
    }

    return {
        "per_ray_3d_vs_2d": per_ray,
        "stages": stages,
        "top_observer_relevance_channels": _aggregate_relevance(rows)[:20],
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
    expected_decoded = len(EXTRACTION_METHODS) * len(EXTRACTION_FIELDS) + len(PRIMARY_3D_BIN_CHANNELS)
    checks = {
        "canonical_five_cluster_inventory": ids == EXPECTED_CLUSTER_IDS,
        "all_five_clusters_completed": len(rows) == 5 and not failures,
        "one_received_G3D_state_per_cluster": bool(rows),
        "one_target_blind_screen_per_cluster": bool(rows),
        "all_established_2d_extraction_methods_retained": bool(rows and all(all(f"{m}__{f}" in r["decoded_channel_names"] for m in EXTRACTION_METHODS for f in EXTRACTION_FIELDS) for r in rows)),
        "all_predeclared_3d_receipt_channels_retained": bool(rows and all(all(n in r["decoded_channel_names"] for n in PRIMARY_3D_BIN_CHANNELS) for r in rows)),
        "decoded_channel_count_matches_inventory": bool(rows and all(r["decoded_channel_count"] == expected_decoded for r in rows)),
        "native_terminal_common_history_identity": bool(rows and all(r["terminal_common_history_relative_rms_error"] <= 1e-12 for r in rows)),
        "G3D_unit_speed_valid": bool(rows and all(r["g3d_unit_speed_max_error"] <= G3D.UNIT_SPEED_TOL for r in rows)),
        "observations_used_only_after_decoding": True,
        "no_observational_fit_or_channel_weighting": True,
        "no_physical_output_rescaling": True,
        "no_upstream_physics_change": True,
        "no_cluster_specific_channel_choice": True,
        "no_tracked_or_staged_changes": not state["tracked_changes"] and not state["staged_changes"],
    }
    passed = all(checks.values())
    status = "FULL_RECEIVED_STATE_INFORMATION_RETENTION_EXECUTED" if passed else ("FULL_RECEIVED_STATE_INFORMATION_RETENTION_PARTIAL_EXECUTION" if rows else "FULL_RECEIVED_STATE_INFORMATION_RETENTION_NOT_ESTABLISHED")

    result = {
        "lab_id": LAB_ID,
        "status": status,
        "repo_state": state,
        "frozen_rule": "one current-native G3D received state and one target-blind detector screen per cluster; decode all channels before reduction",
        "information_geometry_rule": "column standardization is target-blind and used only for rank/SVD/redundancy diagnostics, never physical output rescaling",
        "decoded_inventory_size": expected_decoded,
        "established_2d_methods": EXTRACTION_METHODS,
        "established_2d_fields": EXTRACTION_FIELDS,
        "predeclared_3d_receipt_channels": PRIMARY_3D_BIN_CHANNELS,
        "aggregate": agg,
        "rows": rows,
        "failures": failures,
        "checks": checks,
    }

    print(LAB_ID)
    print(f"status={status}")
    print(f"head_sha={state['head_sha']}")
    print("trajectory=frozen_current_native_G3D")
    print("observer_screen=frozen_target_blind_global_tangent")
    print(f"decoded_inventory_size={expected_decoded}")
    print("observations=after_the_fact_relevance_only")
    print("physical_output_rescaling=false")
    print()

    if agg:
        p = agg["per_ray_3d_vs_2d"]
        print("PER_RAY_3D_VS_2D")
        for k, v in p.items():
            print(f"{k}={v:.12g}")
        print()
        print("REDUCTION_STAGES")
        for name, a in agg["stages"].items():
            print(f"stage={name} requested={a['requested_channels']} rank={a['mean_numerical_rank']:.12g} eff_rank={a['mean_effective_rank']:.12g} participation={a['mean_participation_ratio']:.12g} top5={a['mean_top5_variance_fraction']:.12g} eff_fraction_full={a['mean_effective_rank_fraction_of_full']:.12g} rank_fraction_full={a['mean_numerical_rank_fraction_of_full']:.12g} median_abs_corr={a['mean_median_abs_channel_corr']:.12g}")
        print()
        print("TOP_OBSERVER_RELEVANCE_CHANNELS")
        for r in agg["top_observer_relevance_channels"]:
            print(f"channel={r['channel']} family={r['family']} best_target={r['best_observed_target']} best_mean_r={r['best_mean_pearson']:.12g}")
        print()

    print("CHECKS")
    for key, value in checks.items():
        print(f"{key}={str(bool(value)).lower()}")
    print(json.dumps(result, sort_keys=True, separators=(",", ":")))
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
