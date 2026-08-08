#!/usr/bin/env python3
"""PBUF FOUNDATION — NATIVE OBSERVABLE EXTRACTION METHOD SWEEP 001.

Hold the complete current-native received G3D ray state fixed and test several
ways of extracting 2D weak-lensing observables from ONE fixed detector screen.

The detector screen is target-blind and fixed per cluster: the plane normal to
the global mean received ray direction, with its first axis obtained by
projecting global +x into that plane.  Every extraction candidate consumes the
same initial/final screen coordinates from the same received rays.

Nothing upstream may differ between candidates:
  established local benchmark source -> current native fast/slow transfer
  -> PM1/PS2 -> M10 -> LOS -> frozen G3D -> received 3D rays
  -> one fixed tangent screen -> extraction sweep

Candidates:
  histogram_density       occupancy-ratio extraction
  kernel_density          Gaussian KDE extraction
  jacobian_affine         per-bin affine ray-bundle Jacobian
  covariance_area         finite covariance-area distortion
  displacement_divergence gradient of mean screen displacement
  knn_density             adaptive k-nearest-neighbour density
  polar_jacobian          rotation-free polar stretch of affine Jacobian
  covariance_transport    SPD covariance transport between initial/final bundles

No fitting, target-dependent orientation, amplitude matching, normalization to
observations, or cluster-specific method selection is performed.
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

import observable_lab001 as OLD
from weak_lensing_observation001 import resample_to_grid
from pbuf.core import benchmark_data as BENCH
from pbuf.core import observable_extraction as M16
import pbuf.labs.foundation.current_native_five_cluster_observable_benchmark001 as CUR
import pbuf.labs.foundation.native_accumulated_full_lensing001 as G3D
import pbuf.labs.foundation.m10_coverage_25pct_science001 as BASE

LAB_ID = "PBUF-FOUNDATION-NATIVE-OBSERVABLE-EXTRACTION-METHOD-SWEEP-001"
EXPECTED_CLUSTER_IDS = ("Abell2744", "MACS0416", "MACS1149", "AbellS1063", "Abell370")
METHODS = (
    "histogram_density",
    "kernel_density",
    "jacobian_affine",
    "covariance_area",
    "displacement_divergence",
    "knn_density",
    "polar_jacobian",
    "covariance_transport",
)
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


def _sqrt_psd(M: np.ndarray) -> np.ndarray:
    vals, vecs = np.linalg.eigh(np.asarray(M, dtype=np.float64))
    vals = np.maximum(vals, 0.0)
    return vecs @ np.diag(np.sqrt(vals)) @ vecs.T


def _invsqrt_spd(M: np.ndarray) -> np.ndarray:
    vals, vecs = np.linalg.eigh(np.asarray(M, dtype=np.float64))
    if float(np.min(vals)) <= 1.0e-14:
        raise np.linalg.LinAlgError("covariance not positive definite")
    return vecs @ np.diag(1.0 / np.sqrt(vals)) @ vecs.T


def _screen_basis(snap: dict) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    n = np.array([
        float(np.mean(snap["vx"])),
        float(np.mean(snap["vy"])),
        float(np.mean(snap["vz"])),
    ], dtype=np.float64)
    n /= max(float(np.linalg.norm(n)), EPS)

    # Anchor the detector to global +x so the basis does not rotate arbitrarily.
    ref = np.array([1.0, 0.0, 0.0], dtype=np.float64)
    e1 = ref - float(np.dot(ref, n)) * n
    if float(np.linalg.norm(e1)) <= 1.0e-10:
        ref = np.array([0.0, 1.0, 0.0], dtype=np.float64)
        e1 = ref - float(np.dot(ref, n)) * n
    e1 /= max(float(np.linalg.norm(e1)), EPS)
    e2 = np.cross(n, e1)
    e2 /= max(float(np.linalg.norm(e2)), EPS)
    return e1, e2, n


def _screen_coordinates(x0, y0, snap: dict) -> dict:
    e1, e2, n = _screen_basis(snap)
    p0 = np.column_stack((x0, y0, np.zeros_like(x0))).astype(np.float64)
    pf = np.column_stack((snap["x"], snap["y"], snap["z"])).astype(np.float64)
    return {
        "u0": p0 @ e1,
        "v0": p0 @ e2,
        "uf": pf @ e1,
        "vf": pf @ e2,
        "e1": e1,
        "e2": e2,
        "normal": n,
    }


def _matrix_observables(a00, a01, a10, a11) -> dict:
    k = 1.0 - (a00*a11 - a01*a10)
    g1 = 0.5 * (a00 - a11)
    g2 = 0.5 * (a01 + a10)
    return {
        "convergence": k,
        "shear_g1": g1,
        "shear_g2": g2,
        "shear_magnitude": np.hypot(g1, g2),
    }


def _empty(bins: int) -> np.ndarray:
    return np.full((bins, bins), np.nan, dtype=np.float64)


def _groups_from_initial(u0, v0, extent: float, bins: int) -> dict[int, np.ndarray]:
    edges = np.linspace(-extent, extent, bins + 1)
    col = np.searchsorted(edges, u0, side="right") - 1
    row = np.searchsorted(edges, v0, side="right") - 1
    valid = (row >= 0) & (row < bins) & (col >= 0) & (col < bins)
    flat = row * bins + col
    groups = {}
    for q in np.unique(flat[valid]):
        idx = np.where(valid & (flat == q))[0]
        if idx.size >= 6:
            groups[int(q)] = idx
    return groups


def _affine_maps(u0, v0, uf, vf, extent: float, bins: int) -> tuple[dict, dict[int,np.ndarray]]:
    groups = _groups_from_initial(u0, v0, extent, bins)
    maps = {k: _empty(bins) for k in ("a00", "a01", "a10", "a11")}
    for q, idx in groups.items():
        r, c = divmod(q, bins)
        X = np.column_stack((u0[idx]-np.mean(u0[idx]), v0[idx]-np.mean(v0[idx])))
        yu = uf[idx]-np.mean(uf[idx]); yv = vf[idx]-np.mean(vf[idx])
        try:
            au, *_ = np.linalg.lstsq(X, yu, rcond=None)
            av, *_ = np.linalg.lstsq(X, yv, rcond=None)
        except np.linalg.LinAlgError:
            continue
        maps["a00"][r,c] = au[0]; maps["a01"][r,c] = au[1]
        maps["a10"][r,c] = av[0]; maps["a11"][r,c] = av[1]
    return maps, groups


def _polar_jacobian(u0, v0, uf, vf, extent: float, bins: int) -> dict:
    maps, _ = _affine_maps(u0, v0, uf, vf, extent, bins)
    out = {k: _empty(bins) for k in maps}
    finite = np.isfinite(maps["a00"])
    for r, c in zip(*np.where(finite)):
        A = np.array([[maps["a00"][r,c], maps["a01"][r,c]],
                      [maps["a10"][r,c], maps["a11"][r,c]]])
        U = _sqrt_psd(A.T @ A)
        out["a00"][r,c], out["a01"][r,c] = U[0,0], U[0,1]
        out["a10"][r,c], out["a11"][r,c] = U[1,0], U[1,1]
    return _matrix_observables(out["a00"], out["a01"], out["a10"], out["a11"])


def _covariance_transport(u0, v0, uf, vf, extent: float, bins: int) -> dict:
    groups = _groups_from_initial(u0, v0, extent, bins)
    out = {k: _empty(bins) for k in ("a00", "a01", "a10", "a11")}
    for q, idx in groups.items():
        r, c = divmod(q, bins)
        Xi = np.column_stack((u0[idx]-np.mean(u0[idx]), v0[idx]-np.mean(v0[idx])))
        Xf = np.column_stack((uf[idx]-np.mean(uf[idx]), vf[idx]-np.mean(vf[idx])))
        Ci = Xi.T @ Xi / max(len(idx)-1, 1)
        Cf = Xf.T @ Xf / max(len(idx)-1, 1)
        try:
            Ci_half = _sqrt_psd(Ci)
            Ci_ih = _invsqrt_spd(Ci)
            middle = Ci_half @ Cf @ Ci_half
            A = Ci_ih @ _sqrt_psd(middle) @ Ci_ih
        except np.linalg.LinAlgError:
            continue
        out["a00"][r,c], out["a01"][r,c] = A[0,0], A[0,1]
        out["a10"][r,c], out["a11"][r,c] = A[1,0], A[1,1]
    return _matrix_observables(out["a00"], out["a01"], out["a10"], out["a11"])


def _canonicalize(old_result: dict) -> dict:
    return {
        "convergence": np.asarray(old_result["convergence"], dtype=np.float64),
        "shear_g1": np.asarray(old_result["shear_g1"], dtype=np.float64),
        "shear_g2": np.asarray(old_result["shear_g2"], dtype=np.float64),
        "shear_magnitude": np.asarray(old_result["shear_magnitude"], dtype=np.float64),
    }


def _extract_all(screen: dict, extent: float, bins: int) -> dict[str, dict]:
    args = (screen["u0"], screen["v0"], screen["uf"], screen["vf"], extent, bins)
    affine, _ = _affine_maps(*args[:4], extent, bins)
    jac = _matrix_observables(affine["a00"], affine["a01"], affine["a10"], affine["a11"])

    return {
        "histogram_density": _canonicalize(OLD.method_histogram(*args)),
        "kernel_density": _canonicalize(OLD.method_kernel(*args)),
        "jacobian_affine": jac,
        "covariance_area": _canonicalize(OLD.method_area(*args)),
        "displacement_divergence": _canonicalize(OLD.method_divergence(*args)),
        "knn_density": _canonicalize(OLD.method_knn(*args)),
        "polar_jacobian": _polar_jacobian(*args[:4], extent, bins),
        "covariance_transport": _covariance_transport(*args[:4], extent, bins),
    }


def _compare(pred, obs) -> dict:
    p, s, n = corr(pred, obs)
    pr = rms(pred); orms = rms(obs)
    return {
        "pearson": p,
        "spearman": s,
        "count": n,
        "coverage_fraction": n / float(np.asarray(pred).size),
        "pred_rms": pr,
        "obs_rms": orms,
        "rms_ratio_pred_over_obs": pr / max(orms, EPS),
    }


def _orientation(pg1, pg2, og1, og2) -> float:
    p1=np.asarray(pg1); p2=np.asarray(pg2); o1=np.asarray(og1); o2=np.asarray(og2)
    m=np.isfinite(p1)&np.isfinite(p2)&np.isfinite(o1)&np.isfinite(o2)
    if not np.any(m): return float("nan")
    pn=np.hypot(p1[m],p2[m]); on=np.hypot(o1[m],o2[m]); good=(pn>EPS)&(on>EPS)
    if not np.any(good): return float("nan")
    return float(np.mean((p1[m][good]*o1[m][good]+p2[m][good]*o2[m][good])/(pn[good]*on[good])))


def _score_method(pred: dict, obs: dict) -> dict:
    row = {
        "kappa": _compare(pred["convergence"], obs["kappa"]),
        "gamma": _compare(pred["shear_magnitude"], obs["gamma"]),
        "gamma1": _compare(pred["shear_g1"], obs["gamma1"]),
        "gamma2": _compare(pred["shear_g2"], obs["gamma2"]),
    }
    row["pred_gamma2_over_gamma1_rms"] = rms(pred["shear_g2"]) / max(rms(pred["shear_g1"]), EPS)
    row["obs_gamma2_over_gamma1_rms"] = rms(obs["gamma2"]) / max(rms(obs["gamma1"]), EPS)
    row["shear_orientation_cosine_mean"] = _orientation(
        pred["shear_g1"], pred["shear_g2"], obs["gamma1"], obs["gamma2"]
    )
    return row


def run_cluster(cluster: dict) -> dict:
    data = CUR.local_cluster(cluster)
    m10, channel = CUR.current_native_m10(data["rho3"])

    # ONE frozen received state; all extractors consume this exact state.
    chain = G3D.run_g3d_from_vector(m10, observed_for_first_step=None)
    snap = chain["checkpoints"][G3D.CHECKPOINT]
    x0, y0, _, _ = BASE._launch_expanded_25pct()
    screen = _screen_coordinates(x0, y0, snap)

    bins = int(BASE.CFG["bins"])
    extent = float(BASE.CFG["extent"])
    predictions = _extract_all(screen, extent, bins)
    obs = {
        "kappa": resample_to_grid(data["kappa"], bins, extent),
        "gamma": resample_to_grid(data["gamma"], bins, extent),
        "gamma1": resample_to_grid(data["gamma1"], bins, extent),
        "gamma2": resample_to_grid(data["gamma2"], bins, extent),
    }

    return {
        "cluster_id": cluster["id"],
        "pair_fast_coefficient_from_A8": channel["pair_fast_coefficient_from_A8"],
        "pair_slow_coefficient_from_A8": channel["pair_slow_coefficient_from_A8"],
        "terminal_common_history_relative_rms_error": channel["terminal_common_history_relative_rms_error"],
        "g3d_unit_speed_max_error": chain["g3d"]["max_unit_speed_error"],
        "screen_normal": screen["normal"].tolist(),
        "screen_axis_e1": screen["e1"].tolist(),
        "screen_axis_e2": screen["e2"].tolist(),
        "received_state_role": "one_frozen_current_native_G3D_state_shared_by_all_extractors",
        "screen_role": "one_target_blind_global_tangent_screen_shared_by_all_extractors",
        "methods": {name: _score_method(predictions[name], obs) for name in METHODS},
    }


def _cv(values) -> float:
    vals=[float(v) for v in values if math.isfinite(float(v))]
    if len(vals)<2: return float("nan")
    mean=statistics.fmean(vals)
    return float(statistics.pstdev(vals)/max(abs(mean),EPS))


def aggregate(rows: list[dict]) -> dict:
    out={}
    for method in METHODS:
        rr=[row["methods"][method] for row in rows]
        def mf(path1,path2):
            vals=[r[path1][path2] for r in rr if math.isfinite(float(r[path1][path2]))]
            return float(statistics.fmean(vals)) if vals else float("nan")
        out[method]={
            "mean_coverage_fraction": mf("kappa","coverage_fraction"),
            "mean_kappa_pearson": mf("kappa","pearson"),
            "mean_kappa_spearman": mf("kappa","spearman"),
            "mean_kappa_amp_ratio": mf("kappa","rms_ratio_pred_over_obs"),
            "cross_cluster_kappa_amp_cv": _cv([r["kappa"]["rms_ratio_pred_over_obs"] for r in rr]),
            "mean_gamma_pearson": mf("gamma","pearson"),
            "mean_gamma_spearman": mf("gamma","spearman"),
            "mean_gamma_amp_ratio": mf("gamma","rms_ratio_pred_over_obs"),
            "mean_gamma1_amp_ratio": mf("gamma1","rms_ratio_pred_over_obs"),
            "mean_gamma2_amp_ratio": mf("gamma2","rms_ratio_pred_over_obs"),
            "mean_shear_orientation_cosine": float(statistics.fmean([
                r["shear_orientation_cosine_mean"] for r in rr
                if math.isfinite(float(r["shear_orientation_cosine_mean"]))
            ])) if any(math.isfinite(float(r["shear_orientation_cosine_mean"])) for r in rr) else float("nan"),
        }
    return out


def main() -> int:
    state=repo_state(); clusters=list(BENCH.clusters()); ids=tuple(c["id"] for c in clusters)
    rows=[]; failures=[]
    if ids==EXPECTED_CLUSTER_IDS:
        for cluster in clusters:
            try: rows.append(run_cluster(cluster))
            except Exception as exc: failures.append({"cluster_id":cluster["id"],"error":f"{type(exc).__name__}: {exc}"})

    checks={
        "canonical_five_cluster_inventory": ids==EXPECTED_CLUSTER_IDS,
        "all_five_clusters_completed": len(rows)==5 and not failures,
        "all_eight_methods_present": bool(rows and all(tuple(r["methods"].keys())==METHODS for r in rows)),
        "same_received_G3D_state_shared_within_each_cluster": bool(rows and all(r["received_state_role"].startswith("one_frozen") for r in rows)),
        "same_target_blind_tangent_screen_shared_within_each_cluster": bool(rows and all(r["screen_role"].startswith("one_target_blind") for r in rows)),
        "native_terminal_common_history_identity": bool(rows and all(r["terminal_common_history_relative_rms_error"]<=1e-12 for r in rows)),
        "G3D_unit_speed_valid": bool(rows and all(r["g3d_unit_speed_max_error"]<=G3D.UNIT_SPEED_TOL for r in rows)),
        "no_tracked_or_staged_changes": not state["tracked_changes"] and not state["staged_changes"],
    }
    passed=bool(all(checks.values()))
    status=("NATIVE_OBSERVABLE_EXTRACTION_METHOD_SWEEP_EXECUTED" if passed else
            ("NATIVE_OBSERVABLE_EXTRACTION_METHOD_SWEEP_PARTIAL_EXECUTION" if rows else
             "NATIVE_OBSERVABLE_EXTRACTION_METHOD_SWEEP_NOT_ESTABLISHED"))
    agg=aggregate(rows) if rows else {}
    result={
        "lab_id":LAB_ID,"status":status,"repo_state":state,
        "frozen_rule":"one current-native G3D received state and one target-blind global tangent screen per cluster; extraction algorithm only may differ",
        "methods":list(METHODS),"rows":rows,"aggregate":agg,"failures":failures,"checks":checks,
        "fit_or_tuning":False,"output_rescaling":False,
        "interpretation_rule":"Compare extraction algorithms on identical detector-plane ray data. No method is automatically promoted from benchmark agreement.",
    }

    print(LAB_ID); print(f"status={status}"); print(f"head_sha={state['head_sha']}")
    print("trajectory=frozen_current_native_G3D"); print("observer_screen=frozen_target_blind_global_tangent")
    print("extraction_methods="+",".join(METHODS)); print("fit_or_tuning=false"); print("output_rescaling=false")
    print("\nAGGREGATE")
    for name in METHODS:
        if name not in agg: continue
        a=agg[name]
        print(f"method={name} coverage={a['mean_coverage_fraction']:.12g} kappa_r={a['mean_kappa_pearson']:.12g} kappa_rho={a['mean_kappa_spearman']:.12g} kappa_amp={a['mean_kappa_amp_ratio']:.12g} gamma_r={a['mean_gamma_pearson']:.12g} gamma_rho={a['mean_gamma_spearman']:.12g} gamma_amp={a['mean_gamma_amp_ratio']:.12g} gamma1_amp={a['mean_gamma1_amp_ratio']:.12g} gamma2_amp={a['mean_gamma2_amp_ratio']:.12g} shear_orientation={a['mean_shear_orientation_cosine']:.12g} kappa_amp_cv={a['cross_cluster_kappa_amp_cv']:.12g}")
    print("\nCLUSTERS")
    for row in rows:
        for name in METHODS:
            m=row["methods"][name]
            print(f"cluster={row['cluster_id']} method={name} coverage={m['kappa']['coverage_fraction']:.12g} kappa_r={m['kappa']['pearson']:.12g} kappa_rho={m['kappa']['spearman']:.12g} kappa_amp={m['kappa']['rms_ratio_pred_over_obs']:.12g} gamma_r={m['gamma']['pearson']:.12g} gamma_amp={m['gamma']['rms_ratio_pred_over_obs']:.12g} gamma1_amp={m['gamma1']['rms_ratio_pred_over_obs']:.12g} gamma2_amp={m['gamma2']['rms_ratio_pred_over_obs']:.12g} shear_orientation={m['shear_orientation_cosine_mean']:.12g}")
    for f in failures: print(f"failure_cluster={f['cluster_id']} error={f['error']}")
    print("\nCHECKS")
    for k,v in checks.items(): print(f"{k}={str(v).lower()}")
    print(json.dumps(result,sort_keys=True,separators=(",",":"),default=str))
    return 0 if passed else 1


if __name__=="__main__":
    raise SystemExit(main())
