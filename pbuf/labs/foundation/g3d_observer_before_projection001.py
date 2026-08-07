#!/usr/bin/env python3
"""PBUF FOUNDATION — G3D OBSERVER-BEFORE-PROJECTION 001.

Diagnostic-only follow-up to the G3D principal-axis orientation audit.

The PBUF candidate, M10 interface field, C25 source geometry, photon count,
normalized G3D propagation law, step size/count, and all coefficients are frozen.
No conventional gravitational law is introduced into the PBUF pipeline.

Question tested here:

  Does the full received 3D ray-sheet geometry retain structure that is weakened
  when the received state is collapsed to the transverse 2D observer plane?

The incoming ray sheet is parameterized only by the two source coordinates
(x0,y0), so the natural received differential object is the 3x2 Jacobian

    J3 = d(xf,yf,zf) / d(x0,y0).

Before any observer projection, this lab forms the induced 2x2 metric

    G3 = J3^T J3,

whose eigenvalues are the squared principal stretch factors of the received
2D ray sheet embedded in 3D.  These quantities are intrinsic to the received
3D sheet and do not require choosing a 2D screen.

Only after those 3D quantities are measured does the lab apply the fixed
observer projection normal n=(0,0,1), giving

    J2 = P_xy J3,
    G2 = J2^T J2.

The lab then compares 3D-sheet invariants with their post-projection 2D
counterparts.  Observed kappa remains an external morphology benchmark only;
correlation is never an execution gate and no observable is selected by fit.
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
from pbuf.core import los_projection as M14
from pbuf.core import observable_extraction as M16

LAB_ID = "PBUF-FOUNDATION-G3D-OBSERVER-BEFORE-PROJECTION-001"
OUT = ROOT / "runs" / "g3d_observer_before_projection001"
CHECKPOINTS = GEO.CHECKPOINTS
EXPECTED_SUPPORT = GEO.EXPECTED_SUPPORT
UNIT_SPEED_TOL = GEO.UNIT_SPEED_TOL
FIRST_STEP_TOL = GEO.FIRST_STEP_TOL
METRIC_IDENTITY_TOL = 1e-12
PROJECTION_IDENTITY_TOL = 1e-12
PSD_TOL = 1e-12


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
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)


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
    num = _rms(d[mask])
    den = _rms(r[mask])
    if den <= 1e-30:
        return 0.0 if num <= 1e-30 else float("inf")
    return float(num / den)


def _empty_map() -> np.ndarray:
    return np.full((BASE.OBS_BINS, BASE.OBS_BINS), np.nan, dtype=np.float64)


def _fit_received_jacobians(x0, y0, xf, yf, zf, groups):
    """Fit J3=d(xf,yf,zf)/d(x0,y0) independently in every source bin."""
    names = ("Jxx", "Jxy", "Jyx", "Jyy", "Jzx", "Jzy")
    out = {name: _empty_map() for name in names}

    x0 = np.asarray(x0, dtype=np.float64)
    y0 = np.asarray(y0, dtype=np.float64)
    xf = np.asarray(xf, dtype=np.float64)
    yf = np.asarray(yf, dtype=np.float64)
    zf = np.asarray(zf, dtype=np.float64)

    for q, idx in groups.items():
        r, c = divmod(q, BASE.OBS_BINS)
        xi = x0[idx]
        yi = y0[idx]
        A = np.column_stack([xi-xi.mean(), yi-yi.mean()])
        try:
            gx, *_ = np.linalg.lstsq(A, xf[idx]-xf[idx].mean(), rcond=None)
            gy, *_ = np.linalg.lstsq(A, yf[idx]-yf[idx].mean(), rcond=None)
            gz, *_ = np.linalg.lstsq(A, zf[idx]-zf[idx].mean(), rcond=None)
        except np.linalg.LinAlgError:
            continue
        out["Jxx"][r,c] = float(gx[0])
        out["Jxy"][r,c] = float(gx[1])
        out["Jyx"][r,c] = float(gy[0])
        out["Jyy"][r,c] = float(gy[1])
        out["Jzx"][r,c] = float(gz[0])
        out["Jzy"][r,c] = float(gz[1])

    return out


def _metric_fields(J: dict) -> dict:
    """Construct full-3D received-sheet metric and post-projection 2D metric."""
    Jxx, Jxy = J["Jxx"], J["Jxy"]
    Jyx, Jyy = J["Jyx"], J["Jyy"]
    Jzx, Jzy = J["Jzx"], J["Jzy"]

    # Full 3D induced metric G3=J3^T J3.
    g3_11 = Jxx*Jxx + Jyx*Jyx + Jzx*Jzx
    g3_12 = Jxx*Jxy + Jyx*Jyy + Jzx*Jzy
    g3_22 = Jxy*Jxy + Jyy*Jyy + Jzy*Jzy

    # Fixed observer collapse after 3D observation: P_xy removes only z.
    g2_11 = Jxx*Jxx + Jyx*Jyx
    g2_12 = Jxx*Jxy + Jyx*Jyy
    g2_22 = Jxy*Jxy + Jyy*Jyy

    def invariants(a, b, c, prefix):
        tr = a + c
        det = a*c - b*b
        disc2 = np.maximum((a-c)*(a-c) + 4.0*b*b, 0.0)
        disc = np.sqrt(disc2)
        lam_hi = 0.5*(tr + disc)
        lam_lo = 0.5*(tr - disc)
        lam_hi_clip = np.maximum(lam_hi, 0.0)
        lam_lo_clip = np.maximum(lam_lo, 0.0)
        s_hi = np.sqrt(lam_hi_clip)
        s_lo = np.sqrt(lam_lo_clip)
        area = np.sqrt(np.maximum(det, 0.0))
        stretch_gap = s_hi - s_lo
        anis_ratio = np.full_like(s_hi, np.nan)
        good = s_hi > 1e-30
        anis_ratio[good] = s_lo[good] / s_hi[good]
        metric_frob = np.sqrt(a*a + 2.0*b*b + c*c)
        return {
            f"{prefix}_g11": a,
            f"{prefix}_g12": b,
            f"{prefix}_g22": c,
            f"{prefix}_trace": tr,
            f"{prefix}_det": det,
            f"{prefix}_lambda_high": lam_hi,
            f"{prefix}_lambda_low": lam_lo,
            f"{prefix}_principal_stretch_high": s_hi,
            f"{prefix}_principal_stretch_low": s_lo,
            f"{prefix}_area_scale": area,
            f"{prefix}_stretch_gap": stretch_gap,
            f"{prefix}_anisotropy_ratio_low_over_high": anis_ratio,
            f"{prefix}_metric_frobenius_mag": metric_frob,
        }

    fields = {**J}
    fields.update(invariants(g3_11, g3_12, g3_22, "full3d"))
    fields.update(invariants(g2_11, g2_12, g2_22, "projected2d"))

    # z-gradient contribution removed by the fixed observer projection.
    z_metric_11 = Jzx*Jzx
    z_metric_12 = Jzx*Jzy
    z_metric_22 = Jzy*Jzy
    fields.update({
        "z_gradient_x": Jzx,
        "z_gradient_y": Jzy,
        "z_gradient_mag": np.hypot(Jzx, Jzy),
        "z_metric_g11": z_metric_11,
        "z_metric_g12": z_metric_12,
        "z_metric_g22": z_metric_22,
        "projection_metric_loss_frobenius_mag": np.sqrt(
            z_metric_11*z_metric_11 + 2.0*z_metric_12*z_metric_12 + z_metric_22*z_metric_22
        ),
        "area_scale_loss_full3d_minus_projected2d":
            fields["full3d_area_scale"] - fields["projected2d_area_scale"],
        "principal_stretch_high_loss_full3d_minus_projected2d":
            fields["full3d_principal_stretch_high"] - fields["projected2d_principal_stretch_high"],
        "principal_stretch_low_loss_full3d_minus_projected2d":
            fields["full3d_principal_stretch_low"] - fields["projected2d_principal_stretch_low"],
    })

    return fields


def _metric_gates(fields: dict) -> dict:
    # Identity: G3-G2 must equal outer product of z-gradient row.
    d11 = fields["full3d_g11"] - fields["projected2d_g11"] - fields["z_metric_g11"]
    d12 = fields["full3d_g12"] - fields["projected2d_g12"] - fields["z_metric_g12"]
    d22 = fields["full3d_g22"] - fields["projected2d_g22"] - fields["z_metric_g22"]
    ref = np.sqrt(
        fields["full3d_g11"]**2 + 2.0*fields["full3d_g12"]**2 + fields["full3d_g22"]**2
    )
    diff = np.sqrt(d11*d11 + 2.0*d12*d12 + d22*d22)
    projection_identity = _safe_rel_rms(diff, ref)

    # Trace/determinant eigen identities for full3d metric.
    trace_diff = (
        fields["full3d_lambda_high"] + fields["full3d_lambda_low"] - fields["full3d_trace"]
    )
    det_diff = (
        fields["full3d_lambda_high"] * fields["full3d_lambda_low"] - fields["full3d_det"]
    )
    trace_identity = _safe_rel_rms(trace_diff, fields["full3d_trace"])
    det_identity = _safe_rel_rms(det_diff, fields["full3d_det"])

    mask3 = np.isfinite(fields["full3d_lambda_low"])
    mask2 = np.isfinite(fields["projected2d_lambda_low"])
    min_lam3 = float(np.nanmin(fields["full3d_lambda_low"][mask3])) if np.any(mask3) else float("nan")
    min_lam2 = float(np.nanmin(fields["projected2d_lambda_low"][mask2])) if np.any(mask2) else float("nan")

    return {
        "full3d_metric_trace_identity_relative_rms_error": trace_identity,
        "full3d_metric_det_identity_relative_rms_error": det_identity,
        "projection_metric_identity_relative_rms_error": projection_identity,
        "full3d_metric_min_eigenvalue": min_lam3,
        "projected2d_metric_min_eigenvalue": min_lam2,
        "full3d_metric_psd_pass": bool(min_lam3 >= -PSD_TOL),
        "projected2d_metric_psd_pass": bool(min_lam2 >= -PSD_TOL),
    }


def _add_correlations(row: dict, fields: dict, observed, los_mag) -> None:
    names = (
        "full3d_principal_stretch_high",
        "full3d_principal_stretch_low",
        "full3d_area_scale",
        "full3d_stretch_gap",
        "full3d_anisotropy_ratio_low_over_high",
        "full3d_metric_frobenius_mag",
        "projected2d_principal_stretch_high",
        "projected2d_principal_stretch_low",
        "projected2d_area_scale",
        "projected2d_stretch_gap",
        "projected2d_anisotropy_ratio_low_over_high",
        "projected2d_metric_frobenius_mag",
        "z_gradient_mag",
        "projection_metric_loss_frobenius_mag",
        "area_scale_loss_full3d_minus_projected2d",
        "principal_stretch_high_loss_full3d_minus_projected2d",
        "principal_stretch_low_loss_full3d_minus_projected2d",
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


def _checkpoint(cid, k, snap, x0, y0, groups, observed, los_mag):
    J = _fit_received_jacobians(x0, y0, snap["x"], snap["y"], snap["z"], groups)
    fields = _metric_fields(J)
    gates = _metric_gates(fields)
    row = {
        "cluster_id": cid,
        "step_index": int(k),
        "propagation_distance": float(k * BASE.CFG["step"]),
        **gates,
    }
    _add_correlations(row, fields, observed, los_mag)
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

    checkpoints, g3d = GEO._propagate_g3d(field, BASE.CFG["step"], BASE.CFG["steps"], x0, y0)
    if g3d["max_unit_speed_error"] > UNIT_SPEED_TOL:
        raise RuntimeError(f"{cid}: G3D unit-speed gate failed: {g3d['max_unit_speed_error']}")

    first = GEO._first_step_geometry(field, x0, y0, checkpoints[1], observed, los_mag)
    if not first["first_step_exact_pass"]:
        raise RuntimeError(f"{cid}: first-step exact G3D geometry gate failed")

    rows = []
    all_fields = {}
    for k in CHECKPOINTS:
        row, fields = _checkpoint(cid, k, checkpoints[k], x0, y0, groups, observed, los_mag)
        if np.isfinite(row["full3d_metric_trace_identity_relative_rms_error"]) and row["full3d_metric_trace_identity_relative_rms_error"] > METRIC_IDENTITY_TOL:
            raise RuntimeError(f"{cid}: full3d trace identity failed at step {k}")
        if np.isfinite(row["full3d_metric_det_identity_relative_rms_error"]) and row["full3d_metric_det_identity_relative_rms_error"] > METRIC_IDENTITY_TOL:
            raise RuntimeError(f"{cid}: full3d determinant identity failed at step {k}")
        if np.isfinite(row["projection_metric_identity_relative_rms_error"]) and row["projection_metric_identity_relative_rms_error"] > PROJECTION_IDENTITY_TOL:
            raise RuntimeError(f"{cid}: projection metric identity failed at step {k}")
        if not row["full3d_metric_psd_pass"]:
            raise RuntimeError(f"{cid}: full3d metric PSD gate failed at step {k}")
        if not row["projected2d_metric_psd_pass"]:
            raise RuntimeError(f"{cid}: projected2d metric PSD gate failed at step {k}")
        rows.append(row)
        all_fields[f"step_{k}"] = fields

    final = rows[-1]
    summary = {
        "cluster_id": cid,
        "candidate_id": BASE.CANDIDATE_ID,
        "physical_source_representation": BASE.PHYSICAL_SOURCE,
        "geometry_lane": "G3D_LOS_consistent_diagnostic",
        "observer_order": "full_3d_received_sheet_then_fixed_xy_projection",
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
        "final_full3d_metric_trace_identity_relative_rms_error": final["full3d_metric_trace_identity_relative_rms_error"],
        "final_full3d_metric_det_identity_relative_rms_error": final["full3d_metric_det_identity_relative_rms_error"],
        "final_projection_metric_identity_relative_rms_error": final["projection_metric_identity_relative_rms_error"],
        "final_full3d_metric_psd_pass": final["full3d_metric_psd_pass"],
        "final_projected2d_metric_psd_pass": final["projected2d_metric_psd_pass"],
    }
    copy_keys = (
        "full3d_principal_stretch_high_vs_observed_pearson",
        "full3d_principal_stretch_low_vs_observed_pearson",
        "full3d_area_scale_vs_observed_pearson",
        "full3d_stretch_gap_vs_observed_pearson",
        "full3d_metric_frobenius_mag_vs_observed_pearson",
        "projected2d_principal_stretch_high_vs_observed_pearson",
        "projected2d_principal_stretch_low_vs_observed_pearson",
        "projected2d_area_scale_vs_observed_pearson",
        "projected2d_stretch_gap_vs_observed_pearson",
        "projected2d_metric_frobenius_mag_vs_observed_pearson",
        "z_gradient_mag_vs_observed_pearson",
        "projection_metric_loss_frobenius_mag_vs_observed_pearson",
        "area_scale_loss_full3d_minus_projected2d_vs_observed_pearson",
    )
    for key in copy_keys:
        summary[f"final_{key}"] = final[key]

    return summary, rows, all_fields, {
        "observed_benchmark": observed,
        "los_mag": los_mag,
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
        print(f"[{cid}] full 3D received ray-sheet observation before fixed 2D projection")
        try:
            summary, rows, fields, finals = _run_cluster(cluster)
            summaries.append(summary)
            checkpoint_rows.extend(rows)
            cdir = OUT / "clusters" / cid
            cdir.mkdir(parents=True, exist_ok=True)
            _write_json(cdir / "observer3d_summary.json", summary)
            _write_csv(cdir / "observer3d_checkpoints.csv", rows)
            npz = {}
            for step_name, fd in fields.items():
                for name, arr in fd.items():
                    npz[f"{step_name}__{name}"] = arr
            npz.update(finals)
            np.savez_compressed(cdir / "observer3d_checkpoint_fields.npz", **npz)
        except Exception as exc:
            failures.append({"cluster_id": cid, "error": repr(exc)})
            _write_json(OUT / "cluster_failures.json", failures)
            raise

    _write_csv(OUT / "observer3d_summary.csv", summaries)
    _write_csv(OUT / "observer3d_checkpoint_summary.csv", checkpoint_rows)
    _write_json(OUT / "cluster_failures.json", failures)

    def mean_key(key):
        vals = [float(r[key]) for r in summaries]
        return float(np.mean(vals))

    validation = {
        "lab_id": LAB_ID,
        "outcome": "Outcome A — G3D OBSERVER-BEFORE-PROJECTION AUDIT COMPLETE",
        "head_sha": repo["head_sha"],
        "candidate_id": BASE.CANDIDATE_ID,
        "physical_source_representation": BASE.PHYSICAL_SOURCE,
        "geometry_lane": "G3D_LOS_consistent_diagnostic",
        "observer_order": "full_3d_received_sheet_then_fixed_xy_projection",
        "benchmark_role": "external_morphology_comparison_only",
        "cluster_count_expected": len(BASE.CLUSTERS),
        "cluster_count_completed": len(summaries),
        "checkpoint_steps": list(CHECKPOINTS),
        "metric_identity_tolerance": METRIC_IDENTITY_TOL,
        "projection_identity_tolerance": PROJECTION_IDENTITY_TOL,
        "psd_tolerance": PSD_TOL,
        "all_cluster_g3d_unit_speed_pass": bool(all(r["g3d_unit_speed_pass"] for r in summaries)),
        "all_cluster_first_step_exact_pass": bool(all(r["first_step_exact_pass"] for r in summaries)),
        "all_cluster_full3d_metric_identity_pass": bool(all(
            r["final_full3d_metric_trace_identity_relative_rms_error"] <= METRIC_IDENTITY_TOL and
            r["final_full3d_metric_det_identity_relative_rms_error"] <= METRIC_IDENTITY_TOL
            for r in summaries
        )),
        "all_cluster_projection_metric_identity_pass": bool(all(
            r["final_projection_metric_identity_relative_rms_error"] <= PROJECTION_IDENTITY_TOL
            for r in summaries
        )),
        "all_cluster_full3d_metric_psd_pass": bool(all(r["final_full3d_metric_psd_pass"] for r in summaries)),
        "all_cluster_projected2d_metric_psd_pass": bool(all(r["final_projected2d_metric_psd_pass"] for r in summaries)),
        "mean_los_mag_vs_observed_pearson": mean_key("los_mag_vs_observed_pearson"),
        "mean_final_full3d_principal_stretch_high_vs_observed_pearson": mean_key("final_full3d_principal_stretch_high_vs_observed_pearson"),
        "mean_final_full3d_principal_stretch_low_vs_observed_pearson": mean_key("final_full3d_principal_stretch_low_vs_observed_pearson"),
        "mean_final_full3d_area_scale_vs_observed_pearson": mean_key("final_full3d_area_scale_vs_observed_pearson"),
        "mean_final_full3d_stretch_gap_vs_observed_pearson": mean_key("final_full3d_stretch_gap_vs_observed_pearson"),
        "mean_final_full3d_metric_frobenius_mag_vs_observed_pearson": mean_key("final_full3d_metric_frobenius_mag_vs_observed_pearson"),
        "mean_final_projected2d_principal_stretch_high_vs_observed_pearson": mean_key("final_projected2d_principal_stretch_high_vs_observed_pearson"),
        "mean_final_projected2d_principal_stretch_low_vs_observed_pearson": mean_key("final_projected2d_principal_stretch_low_vs_observed_pearson"),
        "mean_final_projected2d_area_scale_vs_observed_pearson": mean_key("final_projected2d_area_scale_vs_observed_pearson"),
        "mean_final_projected2d_stretch_gap_vs_observed_pearson": mean_key("final_projected2d_stretch_gap_vs_observed_pearson"),
        "mean_final_projected2d_metric_frobenius_mag_vs_observed_pearson": mean_key("final_projected2d_metric_frobenius_mag_vs_observed_pearson"),
        "mean_final_z_gradient_mag_vs_observed_pearson": mean_key("final_z_gradient_mag_vs_observed_pearson"),
        "mean_final_projection_metric_loss_frobenius_mag_vs_observed_pearson": mean_key("final_projection_metric_loss_frobenius_mag_vs_observed_pearson"),
        "mean_final_area_scale_loss_full3d_minus_projected2d_vs_observed_pearson": mean_key("final_area_scale_loss_full3d_minus_projected2d_vs_observed_pearson"),
        "physics_change_authorized": False,
        "candidate_change_authorized": False,
        "production_geometry_change_authorized": False,
        "observer_projection_change_authorized": False,
        "observable_selection_authorized": False,
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
        "observer_order": "full_3d_received_sheet_then_fixed_xy_projection",
        "observer_projection_normal": [0.0, 0.0, 1.0],
        "benchmark_role": "external_morphology_comparison_only",
        "checkpoint_steps": list(CHECKPOINTS),
        "note": "Full received 3D ray-sheet geometry is measured before any fixed 2D observer projection; no conventional gravitational law is injected.",
    })
    print(json.dumps(validation, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
