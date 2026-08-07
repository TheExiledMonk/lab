#!/usr/bin/env python3
"""PBUF FOUNDATION — G3D PRINCIPAL-AXIS ORIENTATION 001.

Diagnostic-only follow-up to G3D displacement topology 001.

The PBUF candidate, M10 interface field, C25 source geometry, photon count,
normalized propagation law, step size/count, and G3D LOS-consistent diagnostic
geometry are frozen. No conventional gravitational law is introduced into the
PBUF pipeline. The public observed-kappa map remains only an external morphology
benchmark.

This lab diagonalizes the symmetric part of the PBUF-produced transverse
 displacement-gradient tensor D=grad(a), with a=(ax,ay), at the same frozen G3D
checkpoints. It measures principal eigenvalues, principal-axis orientation,
headless/nematic orientation coherence, and internal alignment between the
principal axis and (a) the PBUF displacement direction and (b) the sampled M10
response direction.

No topology mode, orientation, sign, geometry, physics change, or production
replacement is selected by this lab.
"""
from __future__ import annotations

import csv
import json
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))

import pbuf.labs.foundation.m10_coverage_25pct_science001 as BASE
import pbuf.labs.foundation.los_consistent_ray_geometry001 as GEO
import pbuf.labs.foundation.g3d_displacement_topology001 as TOP
from pbuf.core import los_projection as M14
from pbuf.core import observable_extraction as M16

LAB_ID = "PBUF-FOUNDATION-G3D-PRINCIPAL-AXIS-ORIENTATION-001"
OUT = ROOT / "runs" / "g3d_principal_axis_orientation001"
CHECKPOINTS = GEO.CHECKPOINTS
EXPECTED_SUPPORT = GEO.EXPECTED_SUPPORT
UNIT_SPEED_TOL = GEO.UNIT_SPEED_TOL
FIRST_STEP_TOL = GEO.FIRST_STEP_TOL
EIGEN_IDENTITY_TOL = 1e-12
RECONSTRUCTION_TOL = 1e-12
ANGLE_EPS = 1e-30


def _git(*args: str) -> str:
    return subprocess.check_output(["git", *args], cwd=str(ROOT), text=True).strip()


def _repo_state() -> dict:
    return {
        "repository": "TheExiledMonk/lab",
        "branch": _git("rev-parse", "--abbrev-ref", "HEAD"),
        "head_sha": _git("rev-parse", "HEAD"),
        "tracked_changes": _git("diff", "--name-only"),
        "staged_changes": _git("diff", "--name-only", "--cached"),
    }


def _json_default(o):
    if isinstance(o, np.floating): return float(o)
    if isinstance(o, np.integer): return int(o)
    if isinstance(o, np.bool_): return bool(o)
    if isinstance(o, Path): return str(o)
    return str(o)


def _write_json(path: Path, obj) -> None:
    path.write_text(json.dumps(obj, indent=2, default=_json_default))


def _write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        path.write_text("")
        return
    keys: list[str] = []
    for row in rows:
        for key in row:
            if key not in keys:
                keys.append(key)
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)


def _corr(a, b):
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    if a.shape != b.shape:
        raise RuntimeError(f"correlation shape mismatch {a.shape} vs {b.shape}")
    mask = np.isfinite(a) & np.isfinite(b)
    n = int(mask.sum())
    if n < 2:
        return float("nan"), float("nan"), n
    return (
        float(M16.safe_pearson(a[mask], b[mask])),
        float(M16.safe_spearman(a[mask], b[mask])),
        n,
    )


def _rms(a) -> float:
    x = np.asarray(a, dtype=np.float64)
    return float(np.sqrt(np.mean(x*x))) if x.size else float("nan")


def _safe_rel_rms(diff, reference) -> float:
    d = np.asarray(diff, dtype=np.float64)
    r = np.asarray(reference, dtype=np.float64)
    mask = np.isfinite(d) & np.isfinite(r)
    if not np.any(mask):
        return float("nan")
    den = _rms(r[mask])
    num = _rms(d[mask])
    if den <= 1e-30:
        return 0.0 if num <= 1e-30 else float("inf")
    return float(num / den)


def _nematic_order(c2, s2) -> tuple[float, int]:
    c2 = np.asarray(c2, dtype=np.float64)
    s2 = np.asarray(s2, dtype=np.float64)
    mask = np.isfinite(c2) & np.isfinite(s2)
    n = int(mask.sum())
    if n == 0:
        return float("nan"), 0
    mc = float(np.mean(c2[mask]))
    ms = float(np.mean(s2[mask]))
    return float(np.hypot(mc, ms)), n


def _neighbor_nematic_coherence(c2, s2) -> tuple[np.ndarray, float, int]:
    """Headless nearest-neighbor orientation coherence, rotation invariant."""
    c2 = np.asarray(c2, dtype=np.float64)
    s2 = np.asarray(s2, dtype=np.float64)
    if c2.shape != s2.shape:
        raise RuntimeError("nematic coherence shape mismatch")
    valid = np.isfinite(c2) & np.isfinite(s2)
    local_sum = np.zeros_like(c2, dtype=np.float64)
    local_n = np.zeros_like(c2, dtype=np.int64)
    pair_values: list[np.ndarray] = []

    # right neighbors
    m = valid[:, :-1] & valid[:, 1:]
    if np.any(m):
        q = c2[:, :-1] * c2[:, 1:] + s2[:, :-1] * s2[:, 1:]
        vals = q[m]
        pair_values.append(vals)
        rr, cc = np.where(m)
        local_sum[rr, cc] += vals
        local_n[rr, cc] += 1
        local_sum[rr, cc + 1] += vals
        local_n[rr, cc + 1] += 1

    # down neighbors
    m = valid[:-1, :] & valid[1:, :]
    if np.any(m):
        q = c2[:-1, :] * c2[1:, :] + s2[:-1, :] * s2[1:, :]
        vals = q[m]
        pair_values.append(vals)
        rr, cc = np.where(m)
        local_sum[rr, cc] += vals
        local_n[rr, cc] += 1
        local_sum[rr + 1, cc] += vals
        local_n[rr + 1, cc] += 1

    local = np.full_like(c2, np.nan, dtype=np.float64)
    use = local_n > 0
    local[use] = local_sum[use] / local_n[use]
    if not pair_values:
        return local, float("nan"), 0
    all_pairs = np.concatenate(pair_values)
    return local, float(np.mean(all_pairs)), int(all_pairs.size)


def _principal_fields(topo: dict) -> dict:
    """Diagonalize the symmetric part of D using invariant 2x2 formulas."""
    Dxx = np.asarray(topo["Dxx"], dtype=np.float64)
    Dxy = np.asarray(topo["Dxy"], dtype=np.float64)
    Dyx = np.asarray(topo["Dyx"], dtype=np.float64)
    Dyy = np.asarray(topo["Dyy"], dtype=np.float64)

    Sxx = Dxx
    Syy = Dyy
    Sxy = 0.5 * (Dxy + Dyx)
    trace = Sxx + Syy
    plus = Sxx - Syy
    cross = 2.0 * Sxy
    gap = np.hypot(plus, cross)

    lam_hi = 0.5 * (trace + gap)
    lam_lo = 0.5 * (trace - gap)
    mean_lam = 0.5 * trace
    abs_hi = np.abs(lam_hi)
    abs_lo = np.abs(lam_lo)
    spectral = np.maximum(abs_hi, abs_lo)
    min_abs = np.minimum(abs_hi, abs_lo)
    symmetric_frobenius = np.sqrt(Sxx*Sxx + 2.0*Sxy*Sxy + Syy*Syy)

    phi = np.full_like(trace, np.nan, dtype=np.float64)
    c2 = np.full_like(trace, np.nan, dtype=np.float64)
    s2 = np.full_like(trace, np.nan, dtype=np.float64)
    valid = np.isfinite(trace) & np.isfinite(gap) & (gap > ANGLE_EPS)
    phi[valid] = 0.5 * np.arctan2(cross[valid], plus[valid])
    c2[valid] = plus[valid] / gap[valid]
    s2[valid] = cross[valid] / gap[valid]

    # Reconstruct symmetric tensor from eigen invariants and orientation.
    rec_xx = np.full_like(trace, np.nan, dtype=np.float64)
    rec_yy = np.full_like(trace, np.nan, dtype=np.float64)
    rec_xy = np.full_like(trace, np.nan, dtype=np.float64)
    rec_xx[valid] = mean_lam[valid] + 0.5*gap[valid]*c2[valid]
    rec_yy[valid] = mean_lam[valid] - 0.5*gap[valid]*c2[valid]
    rec_xy[valid] = 0.5*gap[valid]*s2[valid]

    # Degenerate symmetric tensors have undefined orientation but reconstruct
    # exactly from mean eigenvalue alone.
    deg = np.isfinite(trace) & np.isfinite(gap) & (gap <= ANGLE_EPS)
    rec_xx[deg] = mean_lam[deg]
    rec_yy[deg] = mean_lam[deg]
    rec_xy[deg] = 0.0

    return {
        "symmetric_xx": Sxx,
        "symmetric_xy": Sxy,
        "symmetric_yy": Syy,
        "lambda_high": lam_hi,
        "lambda_low": lam_lo,
        "abs_lambda_high": abs_hi,
        "abs_lambda_low": abs_lo,
        "mean_eigenvalue": mean_lam,
        "eigenvalue_gap": gap,
        "spectral_magnitude": spectral,
        "minimum_abs_eigenvalue": min_abs,
        "symmetric_frobenius_mag": symmetric_frobenius,
        "principal_axis_angle": phi,
        "principal_axis_cos2": c2,
        "principal_axis_sin2": s2,
        "reconstructed_symmetric_xx": rec_xx,
        "reconstructed_symmetric_xy": rec_xy,
        "reconstructed_symmetric_yy": rec_yy,
    }


def _alignment_map(phi, vx, vy) -> np.ndarray:
    """cos(2(phi-theta)); invariant to principal-axis sign reversal."""
    phi = np.asarray(phi, dtype=np.float64)
    vx = np.asarray(vx, dtype=np.float64)
    vy = np.asarray(vy, dtype=np.float64)
    if phi.shape != vx.shape or phi.shape != vy.shape:
        raise RuntimeError("alignment shape mismatch")
    mag = np.hypot(vx, vy)
    out = np.full_like(phi, np.nan, dtype=np.float64)
    valid = np.isfinite(phi) & np.isfinite(vx) & np.isfinite(vy) & (mag > ANGLE_EPS)
    theta = np.arctan2(vy[valid], vx[valid])
    out[valid] = np.cos(2.0 * (phi[valid] - theta))
    return out


def _compare_fields(row: dict, fields: dict, observed, los_mag) -> None:
    names = (
        "lambda_high",
        "lambda_low",
        "abs_lambda_high",
        "abs_lambda_low",
        "mean_eigenvalue",
        "eigenvalue_gap",
        "spectral_magnitude",
        "minimum_abs_eigenvalue",
        "symmetric_frobenius_mag",
        "neighbor_nematic_coherence",
        "principal_vs_displacement_nematic_alignment",
        "principal_vs_response_nematic_alignment",
    )
    for name in names:
        arr = fields[name]
        p, s, n = _corr(arr, observed)
        p2, s2, n2 = _corr(arr, los_mag)
        finite = arr[np.isfinite(arr)]
        row[f"{name}_vs_observed_pearson"] = p
        row[f"{name}_vs_observed_spearman"] = s
        row[f"{name}_vs_observed_count"] = n
        row[f"{name}_vs_los_mag_pearson"] = p2
        row[f"{name}_vs_los_mag_spearman"] = s2
        row[f"{name}_vs_los_mag_count"] = n2
        row[f"{name}_rms"] = _rms(finite) if finite.size else float("nan")


def _checkpoint_orientation(cid, k, snap, x0, y0, groups, observed, los_mag):
    topo = TOP._gradient_tensor_maps(x0, y0, snap["x"], snap["y"], groups)
    disp = TOP._displacement_maps(x0, y0, snap["x"], snap["y"], groups)
    principal = _principal_fields(topo)

    response_x = GEO._mean_map(np.asarray(snap["rx_sample"], dtype=np.float64), groups)
    response_y = GEO._mean_map(np.asarray(snap["ry_sample"], dtype=np.float64), groups)
    principal["sampled_response_x"] = response_x
    principal["sampled_response_y"] = response_y
    principal["displacement_x"] = disp["displacement_x"]
    principal["displacement_y"] = disp["displacement_y"]
    principal["displacement_mag"] = disp["displacement_mag"]

    neighbor_map, neighbor_mean, neighbor_pairs = _neighbor_nematic_coherence(
        principal["principal_axis_cos2"], principal["principal_axis_sin2"]
    )
    principal["neighbor_nematic_coherence"] = neighbor_map
    principal["principal_vs_displacement_nematic_alignment"] = _alignment_map(
        principal["principal_axis_angle"], disp["displacement_x"], disp["displacement_y"]
    )
    principal["principal_vs_response_nematic_alignment"] = _alignment_map(
        principal["principal_axis_angle"], response_x, response_y
    )

    nematic_order, orientation_count = _nematic_order(
        principal["principal_axis_cos2"], principal["principal_axis_sin2"]
    )

    trace_identity = _safe_rel_rms(
        (principal["lambda_high"] + principal["lambda_low"]) - topo["divergence_like"],
        topo["divergence_like"],
    )
    gap_identity = _safe_rel_rms(
        (principal["lambda_high"] - principal["lambda_low"]) - topo["symmetric_traceless_mag"],
        topo["symmetric_traceless_mag"],
    )

    rec_diff_sq = (
        (principal["reconstructed_symmetric_xx"] - principal["symmetric_xx"])**2
        + 2.0*(principal["reconstructed_symmetric_xy"] - principal["symmetric_xy"])**2
        + (principal["reconstructed_symmetric_yy"] - principal["symmetric_yy"])**2
    )
    rec_ref_sq = (
        principal["symmetric_xx"]**2
        + 2.0*principal["symmetric_xy"]**2
        + principal["symmetric_yy"]**2
    )
    reconstruction_rel = _safe_rel_rms(np.sqrt(rec_diff_sq), np.sqrt(rec_ref_sq))

    row = {
        "cluster_id": cid,
        "step_index": int(k),
        "propagation_distance": float(k * BASE.CFG["step"]),
        "eigen_trace_identity_relative_rms_error": trace_identity,
        "eigen_gap_identity_relative_rms_error": gap_identity,
        "symmetric_reconstruction_relative_rms_error": reconstruction_rel,
        "principal_axis_global_nematic_order": nematic_order,
        "principal_axis_orientation_count": orientation_count,
        "neighbor_nematic_coherence_mean": neighbor_mean,
        "neighbor_nematic_pair_count": neighbor_pairs,
    }
    _compare_fields(row, principal, observed, los_mag)

    for name in (
        "principal_vs_displacement_nematic_alignment",
        "principal_vs_response_nematic_alignment",
    ):
        arr = principal[name]
        finite = arr[np.isfinite(arr)]
        row[f"{name}_mean"] = float(np.mean(finite)) if finite.size else float("nan")
        row[f"{name}_mean_abs"] = float(np.mean(np.abs(finite))) if finite.size else float("nan")

    fields = {**topo, **disp, **principal}
    return row, fields


def _run_cluster(cluster):
    cid = cluster["id"]
    real = BASE._load_cluster(cluster)
    state = BASE._evolve(BASE._initial_state(real["rho3"]))
    candidate = BASE._candidate(state)
    vector = BASE._interface_vector(candidate)
    los = M14.project_vector_to_image_plane(*vector, los_axis="z")
    Rx = np.asarray(los["comp_1"], dtype=np.float64)
    Ry = np.asarray(los["comp_2"], dtype=np.float64)
    los_mag = np.hypot(Rx, Ry)
    observed = np.asarray(real["observed_kappa"], dtype=np.float64)
    grid = np.linspace(-BASE.CFG["extent"], BASE.CFG["extent"], Rx.shape[0])
    field = {"xgrid": grid, "ygrid": grid, "rx": Rx, "ry": Ry}

    x0, y0, _, _ = BASE._launch_expanded_25pct()
    groups = GEO._source_groups(x0, y0)
    if len(groups) != EXPECTED_SUPPORT:
        raise RuntimeError(f"{cid}: expected {EXPECTED_SUPPORT} source bins, got {len(groups)}")

    checkpoints, g3d = GEO._propagate_g3d(
        field, BASE.CFG["step"], BASE.CFG["steps"], x0, y0
    )
    if g3d["max_unit_speed_error"] > UNIT_SPEED_TOL:
        raise RuntimeError(f"{cid}: G3D unit-speed gate failed: {g3d['max_unit_speed_error']}")

    first = GEO._first_step_geometry(field, x0, y0, checkpoints[1], observed, los_mag)
    if not first["first_step_exact_pass"]:
        raise RuntimeError(f"{cid}: G3D first-step exact geometry gate failed")

    rows = []
    all_fields = {}
    for k in CHECKPOINTS:
        row, fields = _checkpoint_orientation(
            cid, k, checkpoints[k], x0, y0, groups, observed, los_mag
        )
        for key in (
            "eigen_trace_identity_relative_rms_error",
            "eigen_gap_identity_relative_rms_error",
        ):
            value = row[key]
            if np.isfinite(value) and value > EIGEN_IDENTITY_TOL:
                raise RuntimeError(f"{cid}: {key} failed at step {k}: {value}")
        value = row["symmetric_reconstruction_relative_rms_error"]
        if np.isfinite(value) and value > RECONSTRUCTION_TOL:
            raise RuntimeError(
                f"{cid}: symmetric reconstruction failed at step {k}: {value}"
            )
        rows.append(row)
        all_fields[f"step_{k}"] = fields

    final = rows[-1]
    summary = {
        "cluster_id": cid,
        "candidate_id": BASE.CANDIDATE_ID,
        "physical_source_representation": BASE.PHYSICAL_SOURCE,
        "geometry_lane": "G3D_LOS_consistent_diagnostic",
        "benchmark_role": "external_morphology_comparison_only",
        "n_photons": int(len(x0)),
        "source_supported_bins": int(len(groups)),
        "checkpoint_count": len(rows),
        "g3d_unit_speed_max_error": float(g3d["max_unit_speed_error"]),
        "g3d_unit_speed_pass": bool(g3d["max_unit_speed_error"] <= UNIT_SPEED_TOL),
        "first_step_exact_max_vector_error": first["first_step_exact_max_vector_error"],
        "first_step_exact_pass": first["first_step_exact_pass"],
        "los_mag_vs_observed_pearson": _corr(los_mag, observed)[0],
        "los_mag_vs_observed_spearman": _corr(los_mag, observed)[1],
        "final_eigen_trace_identity_relative_rms_error": final["eigen_trace_identity_relative_rms_error"],
        "final_eigen_gap_identity_relative_rms_error": final["eigen_gap_identity_relative_rms_error"],
        "final_symmetric_reconstruction_relative_rms_error": final["symmetric_reconstruction_relative_rms_error"],
        "final_principal_axis_global_nematic_order": final["principal_axis_global_nematic_order"],
        "final_principal_axis_orientation_count": final["principal_axis_orientation_count"],
        "final_neighbor_nematic_coherence_mean": final["neighbor_nematic_coherence_mean"],
        "final_neighbor_nematic_pair_count": final["neighbor_nematic_pair_count"],
        "final_lambda_high_vs_observed_pearson": final["lambda_high_vs_observed_pearson"],
        "final_lambda_low_vs_observed_pearson": final["lambda_low_vs_observed_pearson"],
        "final_abs_lambda_high_vs_observed_pearson": final["abs_lambda_high_vs_observed_pearson"],
        "final_abs_lambda_low_vs_observed_pearson": final["abs_lambda_low_vs_observed_pearson"],
        "final_eigenvalue_gap_vs_observed_pearson": final["eigenvalue_gap_vs_observed_pearson"],
        "final_spectral_magnitude_vs_observed_pearson": final["spectral_magnitude_vs_observed_pearson"],
        "final_symmetric_frobenius_mag_vs_observed_pearson": final["symmetric_frobenius_mag_vs_observed_pearson"],
        "final_neighbor_nematic_coherence_vs_observed_pearson": final["neighbor_nematic_coherence_vs_observed_pearson"],
        "final_principal_vs_displacement_alignment_mean": final["principal_vs_displacement_nematic_alignment_mean"],
        "final_principal_vs_displacement_alignment_mean_abs": final["principal_vs_displacement_nematic_alignment_mean_abs"],
        "final_principal_vs_displacement_alignment_vs_observed_pearson": final["principal_vs_displacement_nematic_alignment_vs_observed_pearson"],
        "final_principal_vs_response_alignment_mean": final["principal_vs_response_nematic_alignment_mean"],
        "final_principal_vs_response_alignment_mean_abs": final["principal_vs_response_nematic_alignment_mean_abs"],
        "final_principal_vs_response_alignment_vs_observed_pearson": final["principal_vs_response_nematic_alignment_vs_observed_pearson"],
    }
    return summary, rows, all_fields, {
        "los_mag": los_mag,
        "observed_benchmark": observed,
        "g3d_x_final": g3d["x"],
        "g3d_y_final": g3d["y"],
        "g3d_z_final": g3d["z"],
    }


def main() -> int:
    started = time.perf_counter()
    OUT.mkdir(parents=True, exist_ok=True)
    repo = _repo_state()
    _write_json(OUT / "repository_state.json", repo)
    if repo["branch"] != "main" or repo["tracked_changes"] or repo["staged_changes"]:
        validation = {
            "lab_id": LAB_ID,
            "outcome": "REPOSITORY_GATE_FAILURE",
            "head_sha": repo["head_sha"],
        }
        _write_json(OUT / "validation.json", validation)
        print(json.dumps(validation, indent=2))
        return 2

    summaries = []
    checkpoint_rows = []
    failures = []
    for cluster in BASE.CLUSTERS:
        cid = cluster["id"]
        print(f"[{cid}] G3D principal-axis/eigenmode orientation audit")
        try:
            summary, rows, fields, final_arrays = _run_cluster(cluster)
            summaries.append(summary)
            checkpoint_rows.extend(rows)
            cdir = OUT / "clusters" / cid
            cdir.mkdir(parents=True, exist_ok=True)
            _write_json(cdir / "orientation_summary.json", summary)
            _write_csv(cdir / "orientation_checkpoints.csv", rows)
            npz = {}
            for step_name, fd in fields.items():
                for name, arr in fd.items():
                    npz[f"{step_name}__{name}"] = arr
            npz.update(final_arrays)
            np.savez_compressed(cdir / "orientation_checkpoint_fields.npz", **npz)
        except Exception as exc:
            failures.append({"cluster_id": cid, "error": repr(exc)})
            _write_json(OUT / "cluster_failures.json", failures)
            raise

    _write_csv(OUT / "orientation_summary.csv", summaries)
    _write_csv(OUT / "orientation_checkpoint_summary.csv", checkpoint_rows)
    _write_json(OUT / "cluster_failures.json", failures)

    validation = {
        "lab_id": LAB_ID,
        "outcome": "Outcome A — G3D PRINCIPAL-AXIS ORIENTATION AUDIT COMPLETE",
        "head_sha": repo["head_sha"],
        "candidate_id": BASE.CANDIDATE_ID,
        "physical_source_representation": BASE.PHYSICAL_SOURCE,
        "geometry_lane": "G3D_LOS_consistent_diagnostic",
        "benchmark_role": "external_morphology_comparison_only",
        "cluster_count_expected": len(BASE.CLUSTERS),
        "cluster_count_completed": len(summaries),
        "checkpoint_steps": list(CHECKPOINTS),
        "eigen_identity_tolerance": EIGEN_IDENTITY_TOL,
        "symmetric_reconstruction_tolerance": RECONSTRUCTION_TOL,
        "all_cluster_g3d_unit_speed_pass": bool(all(r["g3d_unit_speed_pass"] for r in summaries)),
        "all_cluster_first_step_exact_pass": bool(all(r["first_step_exact_pass"] for r in summaries)),
        "all_cluster_eigen_trace_identity_pass": bool(all(r["final_eigen_trace_identity_relative_rms_error"] <= EIGEN_IDENTITY_TOL for r in summaries)),
        "all_cluster_eigen_gap_identity_pass": bool(all(r["final_eigen_gap_identity_relative_rms_error"] <= EIGEN_IDENTITY_TOL for r in summaries)),
        "all_cluster_symmetric_reconstruction_pass": bool(all(r["final_symmetric_reconstruction_relative_rms_error"] <= RECONSTRUCTION_TOL for r in summaries)),
        "mean_los_mag_vs_observed_pearson": float(np.mean([r["los_mag_vs_observed_pearson"] for r in summaries])),
        "mean_final_abs_lambda_high_vs_observed_pearson": float(np.mean([r["final_abs_lambda_high_vs_observed_pearson"] for r in summaries])),
        "mean_final_abs_lambda_low_vs_observed_pearson": float(np.mean([r["final_abs_lambda_low_vs_observed_pearson"] for r in summaries])),
        "mean_final_eigenvalue_gap_vs_observed_pearson": float(np.mean([r["final_eigenvalue_gap_vs_observed_pearson"] for r in summaries])),
        "mean_final_spectral_magnitude_vs_observed_pearson": float(np.mean([r["final_spectral_magnitude_vs_observed_pearson"] for r in summaries])),
        "mean_final_symmetric_frobenius_mag_vs_observed_pearson": float(np.mean([r["final_symmetric_frobenius_mag_vs_observed_pearson"] for r in summaries])),
        "mean_final_principal_axis_global_nematic_order": float(np.mean([r["final_principal_axis_global_nematic_order"] for r in summaries])),
        "mean_final_neighbor_nematic_coherence": float(np.mean([r["final_neighbor_nematic_coherence_mean"] for r in summaries])),
        "mean_final_principal_vs_displacement_alignment": float(np.mean([r["final_principal_vs_displacement_alignment_mean"] for r in summaries])),
        "mean_final_principal_vs_response_alignment": float(np.mean([r["final_principal_vs_response_alignment_mean"] for r in summaries])),
        "physics_change_authorized": False,
        "candidate_change_authorized": False,
        "production_geometry_change_authorized": False,
        "orientation_selection_authorized": False,
        "eigenmode_selection_authorized": False,
        "next_experiment_authorized": False,
        "science_interpretation_required": True,
        "duration_seconds": float(time.perf_counter() - started),
    }
    _write_json(OUT / "validation.json", validation)
    _write_json(OUT / "run.json", {
        "lab_id": LAB_ID,
        "head_sha": repo["head_sha"],
        "candidate_id": BASE.CANDIDATE_ID,
        "physical_source_representation": BASE.PHYSICAL_SOURCE,
        "geometry_lane": "G3D_LOS_consistent_diagnostic",
        "benchmark_role": "external_morphology_comparison_only",
        "checkpoint_steps": list(CHECKPOINTS),
        "note": "Pure PBUF tensor/eigenorientation diagnostic; no conventional gravitational law injected.",
    })
    print(json.dumps(validation, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
