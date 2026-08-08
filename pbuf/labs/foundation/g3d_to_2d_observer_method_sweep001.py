#!/usr/bin/env python3
"""PBUF FOUNDATION — G3D TO 2D OBSERVER METHOD SWEEP 001.

Wide-net diagnostic of the observer layer only.

The current native source/transport/interface/G3D trajectory is frozen. For each
canonical cluster we generate ONE received 3D ray state and pass that exact same
state through several mathematically defined 3D->2D observer mappings.

No candidate may alter the source, A8 dynamics, pair transfer, M10, LOS, G3D
trajectory, or received rays. Observed lensing products are comparison targets
only and never enter candidate construction.

Candidate observer mappings:
  xy_current       fixed global x-y screen (current control)
  xz_control       fixed global x-z coordinate screen
  yz_control       fixed global y-z coordinate screen
  tangent_global   screen perpendicular to global mean received direction
  tangent_local    screen perpendicular to each source-bin mean received direction
  gram_polar       intrinsic 3D sheet stretch U=sqrt(J3^T J3)
  pca_global       target-blind global PCA screen of received 3D endpoints

All matrix-to-observable conversions use the same frozen extraction convention as
observable_lab001.method_jacobian:
  kappa = 1 - det(A)
  gamma1 = 0.5*(A00-A11)
  gamma2 = 0.5*(A01+A10)
where A is the 2x2 observer mapping produced by the candidate. For gram_polar,
A is the symmetric positive-semidefinite polar stretch U.
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

from weak_lensing_observation001 import resample_to_grid
from pbuf.core import benchmark_data as BENCH
from pbuf.core import observable_extraction as M16
import pbuf.labs.foundation.current_native_five_cluster_observable_benchmark001 as CUR
import pbuf.labs.foundation.native_accumulated_full_lensing001 as G3D
import pbuf.labs.foundation.g3d_observer_before_projection001 as PRE
import pbuf.labs.foundation.m10_coverage_25pct_science001 as BASE

LAB_ID = "PBUF-FOUNDATION-G3D-TO-2D-OBSERVER-METHOD-SWEEP-001"
EXPECTED_CLUSTER_IDS = ("Abell2744", "MACS0416", "MACS1149", "AbellS1063", "Abell370")
METHODS = (
    "xy_current",
    "xz_control",
    "yz_control",
    "tangent_global",
    "tangent_local",
    "gram_polar",
    "pca_global",
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


def _empty() -> np.ndarray:
    return np.full((BASE.OBS_BINS, BASE.OBS_BINS), np.nan, dtype=np.float64)


def _matrix_maps_to_observables(maps: dict) -> dict:
    a = maps["a00"]; b = maps["a01"]; c = maps["a10"]; d = maps["a11"]
    k = 1.0 - (a*d - b*c)
    g1 = 0.5*(a-d)
    g2 = 0.5*(b+c)
    return {"kappa": k, "gamma1": g1, "gamma2": g2, "gamma": np.hypot(g1, g2)}


def _j3_from_pre(J: dict, r: int, c: int) -> np.ndarray:
    return np.array([
        [J["Jxx"][r,c], J["Jxy"][r,c]],
        [J["Jyx"][r,c], J["Jyy"][r,c]],
        [J["Jzx"][r,c], J["Jzy"][r,c]],
    ], dtype=np.float64)


def _basis_from_normal(n: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    n = np.asarray(n, dtype=np.float64)
    nn = float(np.linalg.norm(n))
    if nn <= EPS:
        raise RuntimeError("observer normal has zero norm")
    n = n / nn
    refs = np.eye(3)
    ref = refs[int(np.argmin(np.abs(refs @ n)))]
    e1 = np.cross(n, ref)
    e1 /= max(float(np.linalg.norm(e1)), EPS)
    e2 = np.cross(n, e1)
    e2 /= max(float(np.linalg.norm(e2)), EPS)
    return e1, e2


def _maps_from_projector(J: dict, e1: np.ndarray, e2: np.ndarray) -> dict:
    out = {k: _empty() for k in ("a00","a01","a10","a11")}
    finite = np.isfinite(J["Jxx"])
    for r, c in zip(*np.where(finite)):
        j3 = _j3_from_pre(J, int(r), int(c))
        A = np.vstack((e1, e2)) @ j3
        out["a00"][r,c] = A[0,0]; out["a01"][r,c] = A[0,1]
        out["a10"][r,c] = A[1,0]; out["a11"][r,c] = A[1,1]
    return out


def _maps_coordinate(J: dict, rows: tuple[int,int]) -> dict:
    out = {k: _empty() for k in ("a00","a01","a10","a11")}
    finite = np.isfinite(J["Jxx"])
    for r, c in zip(*np.where(finite)):
        j3 = _j3_from_pre(J, int(r), int(c))
        A = j3[list(rows), :]
        out["a00"][r,c] = A[0,0]; out["a01"][r,c] = A[0,1]
        out["a10"][r,c] = A[1,0]; out["a11"][r,c] = A[1,1]
    return out


def _maps_gram_polar(J: dict) -> dict:
    out = {k: _empty() for k in ("a00","a01","a10","a11")}
    finite = np.isfinite(J["Jxx"])
    for r, c in zip(*np.where(finite)):
        j3 = _j3_from_pre(J, int(r), int(c))
        G = j3.T @ j3
        vals, vecs = np.linalg.eigh(G)
        vals = np.maximum(vals, 0.0)
        U = vecs @ np.diag(np.sqrt(vals)) @ vecs.T
        out["a00"][r,c] = U[0,0]; out["a01"][r,c] = U[0,1]
        out["a10"][r,c] = U[1,0]; out["a11"][r,c] = U[1,1]
    return out


def _maps_tangent_local(J: dict, snap: dict, groups: dict[int,np.ndarray]) -> dict:
    out = {k: _empty() for k in ("a00","a01","a10","a11")}
    vx = np.asarray(snap["vx"], dtype=np.float64)
    vy = np.asarray(snap["vy"], dtype=np.float64)
    vz = np.asarray(snap["vz"], dtype=np.float64)
    for q, idx in groups.items():
        r, c = divmod(int(q), BASE.OBS_BINS)
        if not np.isfinite(J["Jxx"][r,c]):
            continue
        n = np.array([np.mean(vx[idx]), np.mean(vy[idx]), np.mean(vz[idx])], dtype=np.float64)
        e1, e2 = _basis_from_normal(n)
        A = np.vstack((e1,e2)) @ _j3_from_pre(J, r, c)
        out["a00"][r,c] = A[0,0]; out["a01"][r,c] = A[0,1]
        out["a10"][r,c] = A[1,0]; out["a11"][r,c] = A[1,1]
    return out


def _global_pca_basis(snap: dict) -> tuple[np.ndarray,np.ndarray]:
    pts = np.column_stack((snap["x"], snap["y"], snap["z"])).astype(np.float64)
    pts -= np.mean(pts, axis=0, keepdims=True)
    C = pts.T @ pts / max(len(pts)-1, 1)
    vals, vecs = np.linalg.eigh(C)
    order = np.argsort(vals)[::-1]
    e1 = vecs[:, order[0]]; e2 = vecs[:, order[1]]
    # Deterministic signs: largest-magnitude component positive.
    for e in (e1,e2):
        j = int(np.argmax(np.abs(e)))
        if e[j] < 0: e *= -1.0
    return e1, e2


def _observer_candidates(J: dict, snap: dict, groups: dict[int,np.ndarray]) -> dict[str,dict]:
    vx = np.asarray(snap["vx"], dtype=np.float64)
    vy = np.asarray(snap["vy"], dtype=np.float64)
    vz = np.asarray(snap["vz"], dtype=np.float64)
    n_global = np.array([np.mean(vx), np.mean(vy), np.mean(vz)], dtype=np.float64)
    tg1, tg2 = _basis_from_normal(n_global)
    p1, p2 = _global_pca_basis(snap)
    maps = {
        "xy_current": _maps_coordinate(J, (0,1)),
        "xz_control": _maps_coordinate(J, (0,2)),
        "yz_control": _maps_coordinate(J, (1,2)),
        "tangent_global": _maps_from_projector(J, tg1, tg2),
        "tangent_local": _maps_tangent_local(J, snap, groups),
        "gram_polar": _maps_gram_polar(J),
        "pca_global": _maps_from_projector(J, p1, p2),
    }
    return {name: _matrix_maps_to_observables(m) for name,m in maps.items()}


def _compare(pred: np.ndarray, obs: np.ndarray) -> dict:
    p, s, n = corr(pred, obs)
    return {
        "pearson": p,
        "spearman": s,
        "count": n,
        "coverage_fraction": n / float(pred.size),
        "pred_rms": rms(pred),
        "obs_rms": rms(obs),
        "rms_ratio_pred_over_obs": rms(pred) / max(rms(obs), EPS),
    }


def _orientation_consistency(pg1, pg2, og1, og2) -> float:
    p1=np.asarray(pg1); p2=np.asarray(pg2); o1=np.asarray(og1); o2=np.asarray(og2)
    m=np.isfinite(p1)&np.isfinite(p2)&np.isfinite(o1)&np.isfinite(o2)
    if not np.any(m): return float("nan")
    pn=np.hypot(p1[m],p2[m]); on=np.hypot(o1[m],o2[m]); good=(pn>EPS)&(on>EPS)
    if not np.any(good): return float("nan")
    return float(np.mean((p1[m][good]*o1[m][good]+p2[m][good]*o2[m][good])/(pn[good]*on[good])))


def run_cluster(cluster: dict) -> dict:
    data = CUR.local_cluster(cluster)
    m10, channel = CUR.current_native_m10(data["rho3"])

    # ONE frozen received G3D ray state for every observer candidate.
    chain = G3D.run_g3d_from_vector(m10, observed_for_first_step=None)
    snap = chain["checkpoints"][G3D.CHECKPOINT]
    x0, y0, _, _ = BASE._launch_expanded_25pct()
    J = PRE._fit_received_jacobians(x0, y0, snap["x"], snap["y"], snap["z"], chain["groups"])
    candidates = _observer_candidates(J, snap, chain["groups"])

    bins=int(BASE.CFG["bins"]); extent=float(BASE.CFG["extent"])
    obs={
        "kappa": resample_to_grid(data["kappa"], bins, extent),
        "gamma": resample_to_grid(data["gamma"], bins, extent),
        "gamma1": resample_to_grid(data["gamma1"], bins, extent),
        "gamma2": resample_to_grid(data["gamma2"], bins, extent),
    }

    methods={}
    for name,pred in candidates.items():
        row={key:_compare(pred[key],obs[key]) for key in ("kappa","gamma","gamma1","gamma2")}
        row["shear_orientation_cosine_mean"]=_orientation_consistency(pred["gamma1"],pred["gamma2"],obs["gamma1"],obs["gamma2"])
        row["pred_gamma2_over_gamma1_rms"] = rms(pred["gamma2"])/max(rms(pred["gamma1"]),EPS)
        row["obs_gamma2_over_gamma1_rms"] = rms(obs["gamma2"])/max(rms(obs["gamma1"]),EPS)
        methods[name]=row

    return {
        "cluster_id": cluster["id"],
        "received_state_sha_role": "one_frozen_G3D_received_state_shared_by_all_observers",
        "g3d_unit_speed_max_error": float(chain["g3d"]["max_unit_speed_error"]),
        **channel,
        "methods": methods,
    }


def _aggregate(rows: list[dict]) -> dict:
    out={}
    for name in METHODS:
        mr=[r["methods"][name] for r in rows]
        out[name]={
            "mean_kappa_pearson": float(np.mean([x["kappa"]["pearson"] for x in mr])),
            "mean_kappa_spearman": float(np.mean([x["kappa"]["spearman"] for x in mr])),
            "mean_kappa_amp_ratio": float(np.mean([x["kappa"]["rms_ratio_pred_over_obs"] for x in mr])),
            "mean_gamma_pearson": float(np.mean([x["gamma"]["pearson"] for x in mr])),
            "mean_gamma_spearman": float(np.mean([x["gamma"]["spearman"] for x in mr])),
            "mean_gamma_amp_ratio": float(np.mean([x["gamma"]["rms_ratio_pred_over_obs"] for x in mr])),
            "mean_gamma1_amp_ratio": float(np.mean([x["gamma1"]["rms_ratio_pred_over_obs"] for x in mr])),
            "mean_gamma2_amp_ratio": float(np.mean([x["gamma2"]["rms_ratio_pred_over_obs"] for x in mr])),
            "mean_shear_orientation_cosine": float(np.mean([x["shear_orientation_cosine_mean"] for x in mr])),
            "mean_coverage_fraction": float(np.mean([x["kappa"]["coverage_fraction"] for x in mr])),
            "cross_cluster_kappa_amp_cv": float(np.std([x["kappa"]["rms_ratio_pred_over_obs"] for x in mr]) / max(abs(np.mean([x["kappa"]["rms_ratio_pred_over_obs"] for x in mr])), EPS)),
        }
    return out


def main() -> int:
    state=repo_state(); clusters=list(BENCH.clusters()); ids=tuple(c["id"] for c in clusters)
    rows=[]; failures=[]
    for cluster in clusters:
        try: rows.append(run_cluster(cluster))
        except Exception as exc: failures.append({"cluster_id":cluster["id"],"error":f"{type(exc).__name__}: {exc}"})

    aggregate=_aggregate(rows) if len(rows)==5 else {}
    checks={
        "canonical_five_cluster_inventory": ids==EXPECTED_CLUSTER_IDS,
        "all_five_clusters_completed": len(rows)==5 and not failures,
        "all_seven_methods_present": bool(rows and all(tuple(r["methods"].keys())==METHODS for r in rows)),
        "same_received_G3D_state_shared_within_each_cluster": bool(rows),
        "native_terminal_common_history_identity": bool(rows and all(r["terminal_common_history_relative_rms_error"]<=1e-12 for r in rows)),
        "G3D_unit_speed_valid": bool(rows and all(r["g3d_unit_speed_max_error"]<=G3D.UNIT_SPEED_TOL for r in rows)),
        "all_metrics_finite": bool(rows and all(math.isfinite(float(r["methods"][m][o][q])) for r in rows for m in METHODS for o in ("kappa","gamma","gamma1","gamma2") for q in ("pearson","spearman","rms_ratio_pred_over_obs"))),
        "no_tracked_or_staged_changes": not state["tracked_changes"] and not state["staged_changes"],
    }
    passed=all(checks.values())
    status="G3D_TO_2D_OBSERVER_METHOD_SWEEP_EXECUTED" if passed else ("G3D_TO_2D_OBSERVER_METHOD_SWEEP_PARTIAL_EXECUTION" if rows else "G3D_TO_2D_OBSERVER_METHOD_SWEEP_NOT_ESTABLISHED")
    result={
        "lab_id":LAB_ID,"status":status,"repo_state":state,
        "frozen_rule":"one current-native G3D received ray state per cluster; observers only may differ",
        "methods":METHODS,"rows":rows,"aggregate":aggregate,"failures":failures,"checks":checks,
        "interpretation_rule":"Compare observer transforms on identical received 3D rays. Rankings are diagnostic only; no candidate is fitted, rescaled, or promoted automatically from benchmark agreement.",
    }
    print(LAB_ID); print(f"status={status}"); print(f"head_sha={state['head_sha']}")
    print("trajectory=frozen_current_native_G3D")
    print("observer_methods="+",".join(METHODS)); print("fit_or_tuning=false"); print("output_rescaling=false")
    print("\nAGGREGATE")
    for name in METHODS:
        if name not in aggregate: continue
        a=aggregate[name]
        print(f"method={name} coverage={a['mean_coverage_fraction']:.12g} kappa_r={a['mean_kappa_pearson']:.12g} kappa_rho={a['mean_kappa_spearman']:.12g} kappa_amp={a['mean_kappa_amp_ratio']:.12g} gamma_r={a['mean_gamma_pearson']:.12g} gamma_rho={a['mean_gamma_spearman']:.12g} gamma_amp={a['mean_gamma_amp_ratio']:.12g} gamma1_amp={a['mean_gamma1_amp_ratio']:.12g} gamma2_amp={a['mean_gamma2_amp_ratio']:.12g} shear_orientation={a['mean_shear_orientation_cosine']:.12g} kappa_amp_cv={a['cross_cluster_kappa_amp_cv']:.12g}")
    print("\nCLUSTERS")
    for row in rows:
        for name in METHODS:
            x=row["methods"][name]
            print(f"cluster={row['cluster_id']} method={name} coverage={x['kappa']['coverage_fraction']:.12g} kappa_r={x['kappa']['pearson']:.12g} kappa_rho={x['kappa']['spearman']:.12g} kappa_amp={x['kappa']['rms_ratio_pred_over_obs']:.12g} gamma_r={x['gamma']['pearson']:.12g} gamma_amp={x['gamma']['rms_ratio_pred_over_obs']:.12g} gamma1_amp={x['gamma1']['rms_ratio_pred_over_obs']:.12g} gamma2_amp={x['gamma2']['rms_ratio_pred_over_obs']:.12g} shear_orientation={x['shear_orientation_cosine_mean']:.12g}")
    for f in failures: print(f"failure_cluster={f['cluster_id']} error={f['error']}")
    print("\nCHECKS")
    for k,v in checks.items(): print(f"{k}={str(v).lower()}")
    print(json.dumps(result,sort_keys=True,separators=(",",":"),default=str))
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
