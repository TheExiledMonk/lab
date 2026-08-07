#!/usr/bin/env python3
"""PBUF FOUNDATION — G3D FULL RECEIVED STATE 001.

Diagnostic-only follow-up to observer-before-projection 001.

The PBUF candidate, M10 interface field, C25 source geometry, photon count,
normalized G3D propagation law, step size/count, coefficients, and fixed
observer normal n=(0,0,1) are frozen. No conventional gravitational law is
introduced into the PBUF pipeline.

Question tested here:

  Does morphology live in the received ray DIRECTION state even when the
  received POSITION sheet looks almost identical before and after xy projection?

The observer receives two independent channels and this lab keeps them
separate so no arbitrary position/direction weighting is introduced:

  position channel:  r=(x,y,z)
  direction channel: v=(vx,vy,vz), |v|=1

For each channel the natural full-3D differential object is a 3x2 Jacobian
with respect to the two source-sheet coordinates (x0,y0). The induced metric
G=J^T J is measured in full 3D first. Only afterward is the fixed xy observer
projection applied. For the direction channel the lab additionally constructs
observer tangent-plane direction coordinates

  tx=vx/vz,  ty=vy/vz

only after the full 3D direction state has been measured. Position and
direction metrics are never added together. Observed kappa is an external
morphology benchmark only; correlation is never an execution gate.
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
import pbuf.labs.foundation.g3d_observer_before_projection001 as OBS3
from pbuf.core import los_projection as M14
from pbuf.core import observable_extraction as M16

LAB_ID = "PBUF-FOUNDATION-G3D-FULL-RECEIVED-STATE-001"
OUT = ROOT / "runs" / "g3d_full_received_state001"
CHECKPOINTS = GEO.CHECKPOINTS
EXPECTED_SUPPORT = GEO.EXPECTED_SUPPORT
UNIT_SPEED_TOL = GEO.UNIT_SPEED_TOL
FIRST_STEP_TOL = GEO.FIRST_STEP_TOL
METRIC_IDENTITY_TOL = 1e-12
PROJECTION_IDENTITY_TOL = 1e-12
PSD_TOL = 1e-12
VZ_MIN = 1e-12


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


def _empty_map() -> np.ndarray:
    return np.full((BASE.OBS_BINS, BASE.OBS_BINS), np.nan, dtype=np.float64)


def _prefixed_metric(J: dict, prefix: str):
    """Reuse the validated 3x2 metric construction, then namespace its fields."""
    raw = OBS3._metric_fields(J)
    out = {}
    for key, value in raw.items():
        if key.startswith("full3d_"):
            out[f"{prefix}_full3d_{key[len('full3d_'):]}"] = value
        elif key.startswith("projected2d_"):
            out[f"{prefix}_projected2d_{key[len('projected2d_'):]}"] = value
        elif key.startswith("z_"):
            out[f"{prefix}_{key}"] = value
        elif key.startswith("projection_metric_loss_"):
            out[f"{prefix}_{key}"] = value
        elif key.startswith("area_scale_loss_"):
            out[f"{prefix}_{key}"] = value
        elif key.startswith("principal_stretch_") and "_loss_" in key:
            out[f"{prefix}_{key}"] = value
        elif key.startswith("J"):
            out[f"{prefix}_{key}"] = value
    gates = OBS3._metric_gates(raw)
    gates = {f"{prefix}_{k}": v for k, v in gates.items()}
    return out, gates


def _fit_2d_jacobian(x0, y0, u, w, groups):
    names = ("Jux", "Juy", "Jwx", "Jwy")
    out = {name: _empty_map() for name in names}
    x0 = np.asarray(x0, dtype=np.float64)
    y0 = np.asarray(y0, dtype=np.float64)
    u = np.asarray(u, dtype=np.float64)
    w = np.asarray(w, dtype=np.float64)
    for q, idx in groups.items():
        r, c = divmod(q, BASE.OBS_BINS)
        A = np.column_stack([x0[idx]-x0[idx].mean(), y0[idx]-y0[idx].mean()])
        try:
            gu, *_ = np.linalg.lstsq(A, u[idx]-u[idx].mean(), rcond=None)
            gw, *_ = np.linalg.lstsq(A, w[idx]-w[idx].mean(), rcond=None)
        except np.linalg.LinAlgError:
            continue
        out["Jux"][r,c] = float(gu[0]); out["Juy"][r,c] = float(gu[1])
        out["Jwx"][r,c] = float(gw[0]); out["Jwy"][r,c] = float(gw[1])
    return out


def _metric2(J: dict, prefix: str):
    a = J["Jux"]*J["Jux"] + J["Jwx"]*J["Jwx"]
    b = J["Jux"]*J["Juy"] + J["Jwx"]*J["Jwy"]
    c = J["Juy"]*J["Juy"] + J["Jwy"]*J["Jwy"]
    tr = a+c
    det = a*c-b*b
    disc = np.sqrt(np.maximum((a-c)*(a-c)+4.0*b*b, 0.0))
    hi = 0.5*(tr+disc); lo = 0.5*(tr-disc)
    shi = np.sqrt(np.maximum(hi, 0.0)); slo = np.sqrt(np.maximum(lo, 0.0))
    area = np.sqrt(np.maximum(det, 0.0))
    ratio = np.full_like(shi, np.nan)
    good = shi > 1e-30
    ratio[good] = slo[good]/shi[good]
    frob = np.sqrt(a*a+2.0*b*b+c*c)
    fields = {
        f"{prefix}_g11": a, f"{prefix}_g12": b, f"{prefix}_g22": c,
        f"{prefix}_trace": tr, f"{prefix}_det": det,
        f"{prefix}_lambda_high": hi, f"{prefix}_lambda_low": lo,
        f"{prefix}_principal_stretch_high": shi,
        f"{prefix}_principal_stretch_low": slo,
        f"{prefix}_area_scale": area,
        f"{prefix}_stretch_gap": shi-slo,
        f"{prefix}_anisotropy_ratio_low_over_high": ratio,
        f"{prefix}_metric_frobenius_mag": frob,
    }
    mask = np.isfinite(lo)
    min_eig = float(np.nanmin(lo[mask])) if np.any(mask) else float("nan")
    trace_err = OBS3._safe_rel_rms(hi+lo-tr, tr)
    det_err = OBS3._safe_rel_rms(hi*lo-det, det)
    gates = {
        f"{prefix}_metric_trace_identity_relative_rms_error": trace_err,
        f"{prefix}_metric_det_identity_relative_rms_error": det_err,
        f"{prefix}_metric_min_eigenvalue": min_eig,
        f"{prefix}_metric_psd_pass": bool(min_eig >= -PSD_TOL),
    }
    return fields, gates


def _mean_maps(snap, groups):
    vx = np.asarray(snap["vx"], dtype=np.float64)
    vy = np.asarray(snap["vy"], dtype=np.float64)
    vz = np.asarray(snap["vz"], dtype=np.float64)
    if np.min(np.abs(vz)) <= VZ_MIN:
        raise RuntimeError(f"observer tangent projection vz too small: {np.min(np.abs(vz))}")
    tx = vx/vz; ty = vy/vz
    return {
        "direction_transverse_mag": GEO._mean_map(np.hypot(vx,vy), groups),
        "direction_los_component": GEO._mean_map(vz, groups),
        "observer_tangent_angle_mag": GEO._mean_map(np.hypot(tx,ty), groups),
    }, tx, ty


def _add_correlations(row: dict, fields: dict, observed, los_mag, names) -> None:
    for name in names:
        arr = fields[name]
        p,s,n = _corr(arr, observed)
        p2,s2,n2 = _corr(arr, los_mag)
        finite = arr[np.isfinite(arr)]
        row[f"{name}_vs_observed_pearson"] = p
        row[f"{name}_vs_observed_spearman"] = s
        row[f"{name}_vs_observed_count"] = n
        row[f"{name}_vs_los_mag_pearson"] = p2
        row[f"{name}_vs_los_mag_spearman"] = s2
        row[f"{name}_vs_los_mag_count"] = n2
        row[f"{name}_rms"] = _rms(finite) if finite.size else float("nan")


def _checkpoint(cid, k, snap, x0, y0, groups, observed, los_mag):
    # Position channel: full received 3D position first, xy projection second.
    Jp = OBS3._fit_received_jacobians(x0,y0,snap["x"],snap["y"],snap["z"],groups)
    pos_fields, pos_gates = _prefixed_metric(Jp, "position")

    # Direction channel: full received 3D unit direction first, xy projection second.
    Jv = OBS3._fit_received_jacobians(x0,y0,snap["vx"],snap["vy"],snap["vz"],groups)
    dir_fields, dir_gates = _prefixed_metric(Jv, "direction")

    direct_maps, tx, ty = _mean_maps(snap, groups)
    Jt = _fit_2d_jacobian(x0,y0,tx,ty,groups)
    tangent_fields, tangent_gates = _metric2(Jt, "observer_tangent")

    fields = {**pos_fields, **dir_fields, **direct_maps, **tangent_fields}
    row = {
        "cluster_id": cid,
        "step_index": int(k),
        "propagation_distance": float(k*BASE.CFG["step"]),
        **pos_gates,
        **dir_gates,
        **tangent_gates,
    }

    names = (
        "position_full3d_principal_stretch_high",
        "position_full3d_principal_stretch_low",
        "position_full3d_area_scale",
        "position_full3d_stretch_gap",
        "position_full3d_metric_frobenius_mag",
        "position_projected2d_principal_stretch_high",
        "position_projected2d_principal_stretch_low",
        "position_projected2d_area_scale",
        "position_projected2d_stretch_gap",
        "position_projected2d_metric_frobenius_mag",
        "direction_full3d_principal_stretch_high",
        "direction_full3d_principal_stretch_low",
        "direction_full3d_area_scale",
        "direction_full3d_stretch_gap",
        "direction_full3d_metric_frobenius_mag",
        "direction_projected2d_principal_stretch_high",
        "direction_projected2d_principal_stretch_low",
        "direction_projected2d_area_scale",
        "direction_projected2d_stretch_gap",
        "direction_projected2d_metric_frobenius_mag",
        "direction_z_gradient_mag",
        "direction_projection_metric_loss_frobenius_mag",
        "direction_transverse_mag",
        "direction_los_component",
        "observer_tangent_angle_mag",
        "observer_tangent_principal_stretch_high",
        "observer_tangent_principal_stretch_low",
        "observer_tangent_area_scale",
        "observer_tangent_stretch_gap",
        "observer_tangent_metric_frobenius_mag",
    )
    _add_correlations(row, fields, observed, los_mag, names)
    return row, fields


def _require_checkpoint_gates(cid, k, row):
    for prefix in ("position", "direction"):
        if row[f"{prefix}_full3d_metric_trace_identity_relative_rms_error"] > METRIC_IDENTITY_TOL:
            raise RuntimeError(f"{cid}: {prefix} full3d trace identity failed at step {k}")
        if row[f"{prefix}_full3d_metric_det_identity_relative_rms_error"] > METRIC_IDENTITY_TOL:
            raise RuntimeError(f"{cid}: {prefix} full3d determinant identity failed at step {k}")
        if row[f"{prefix}_projection_metric_identity_relative_rms_error"] > PROJECTION_IDENTITY_TOL:
            raise RuntimeError(f"{cid}: {prefix} projection identity failed at step {k}")
        if not row[f"{prefix}_full3d_metric_psd_pass"]:
            raise RuntimeError(f"{cid}: {prefix} full3d metric PSD failed at step {k}")
        if not row[f"{prefix}_projected2d_metric_psd_pass"]:
            raise RuntimeError(f"{cid}: {prefix} projected2d metric PSD failed at step {k}")
    if row["observer_tangent_metric_trace_identity_relative_rms_error"] > METRIC_IDENTITY_TOL:
        raise RuntimeError(f"{cid}: tangent metric trace identity failed at step {k}")
    if row["observer_tangent_metric_det_identity_relative_rms_error"] > METRIC_IDENTITY_TOL:
        raise RuntimeError(f"{cid}: tangent metric determinant identity failed at step {k}")
    if not row["observer_tangent_metric_psd_pass"]:
        raise RuntimeError(f"{cid}: tangent metric PSD failed at step {k}")


def _run_cluster(cluster):
    cid = cluster["id"]
    real = BASE._load_cluster(cluster)
    state = BASE._evolve(BASE._initial_state(real["rho3"]))
    candidate = BASE._candidate(state)
    vector = BASE._interface_vector(candidate)
    los = M14.project_vector_to_image_plane(*vector, los_axis="z")
    Rx = np.asarray(los["comp_1"], dtype=np.float64)
    Ry = np.asarray(los["comp_2"], dtype=np.float64)
    los_mag = np.hypot(Rx,Ry)
    observed = np.asarray(real["observed_kappa"], dtype=np.float64)
    grid = np.linspace(-BASE.CFG["extent"], BASE.CFG["extent"], Rx.shape[0])
    field = {"xgrid":grid,"ygrid":grid,"rx":Rx,"ry":Ry}

    x0,y0,_,_ = BASE._launch_expanded_25pct()
    groups = GEO._source_groups(x0,y0)
    if len(groups) != EXPECTED_SUPPORT:
        raise RuntimeError(f"{cid}: expected {EXPECTED_SUPPORT} source bins, got {len(groups)}")

    checkpoints,g3d = GEO._propagate_g3d(field,BASE.CFG["step"],BASE.CFG["steps"],x0,y0)
    if g3d["max_unit_speed_error"] > UNIT_SPEED_TOL:
        raise RuntimeError(f"{cid}: G3D unit-speed gate failed: {g3d['max_unit_speed_error']}")
    first = GEO._first_step_geometry(field,x0,y0,checkpoints[1],observed,los_mag)
    if not first["first_step_exact_pass"]:
        raise RuntimeError(f"{cid}: first-step exact geometry gate failed")

    rows=[]; all_fields={}
    for k in CHECKPOINTS:
        row,fields = _checkpoint(cid,k,checkpoints[k],x0,y0,groups,observed,los_mag)
        _require_checkpoint_gates(cid,k,row)
        rows.append(row); all_fields[f"step_{k}"]=fields

    final = rows[-1]
    summary = {
        "cluster_id": cid,
        "candidate_id": BASE.CANDIDATE_ID,
        "physical_source_representation": BASE.PHYSICAL_SOURCE,
        "geometry_lane": "G3D_LOS_consistent_diagnostic",
        "observer_order": "full_3d_position_plus_direction_then_fixed_xy_projection",
        "observer_channels": "position_and_direction_kept_separate_no_arbitrary_weighting",
        "benchmark_role": "external_morphology_comparison_only",
        "n_photons": int(len(x0)),
        "source_supported_bins": int(len(groups)),
        "checkpoint_count": len(rows),
        "g3d_unit_speed_max_error": float(g3d["max_unit_speed_error"]),
        "g3d_unit_speed_pass": bool(g3d["max_unit_speed_error"] <= UNIT_SPEED_TOL),
        "first_step_exact_max_vector_error": first["first_step_exact_max_vector_error"],
        "first_step_exact_pass": first["first_step_exact_pass"],
        "los_mag_vs_observed_pearson": _corr(los_mag,observed)[0],
        "los_mag_vs_observed_spearman": _corr(los_mag,observed)[1],
    }
    keep = (
        "position_full3d_metric_trace_identity_relative_rms_error",
        "position_projection_metric_identity_relative_rms_error",
        "direction_full3d_metric_trace_identity_relative_rms_error",
        "direction_projection_metric_identity_relative_rms_error",
        "observer_tangent_metric_trace_identity_relative_rms_error",
        "position_full3d_stretch_gap_vs_observed_pearson",
        "position_projected2d_stretch_gap_vs_observed_pearson",
        "direction_full3d_principal_stretch_high_vs_observed_pearson",
        "direction_full3d_principal_stretch_low_vs_observed_pearson",
        "direction_full3d_area_scale_vs_observed_pearson",
        "direction_full3d_stretch_gap_vs_observed_pearson",
        "direction_full3d_metric_frobenius_mag_vs_observed_pearson",
        "direction_projected2d_principal_stretch_high_vs_observed_pearson",
        "direction_projected2d_principal_stretch_low_vs_observed_pearson",
        "direction_projected2d_area_scale_vs_observed_pearson",
        "direction_projected2d_stretch_gap_vs_observed_pearson",
        "direction_projected2d_metric_frobenius_mag_vs_observed_pearson",
        "direction_z_gradient_mag_vs_observed_pearson",
        "direction_projection_metric_loss_frobenius_mag_vs_observed_pearson",
        "direction_transverse_mag_vs_observed_pearson",
        "observer_tangent_angle_mag_vs_observed_pearson",
        "observer_tangent_principal_stretch_high_vs_observed_pearson",
        "observer_tangent_principal_stretch_low_vs_observed_pearson",
        "observer_tangent_area_scale_vs_observed_pearson",
        "observer_tangent_stretch_gap_vs_observed_pearson",
        "observer_tangent_metric_frobenius_mag_vs_observed_pearson",
    )
    for key in keep:
        summary[f"final_{key}"] = final[key]
    return summary,rows,all_fields,{"los_mag":los_mag,"observed_benchmark":observed}


def _nanmean_key(summaries,key):
    vals=np.asarray([r[key] for r in summaries],dtype=np.float64)
    return float(np.nanmean(vals)) if np.any(np.isfinite(vals)) else float("nan")


def main() -> int:
    started=time.perf_counter()
    OUT.mkdir(parents=True,exist_ok=True)
    repo=_repo_state()
    _write_json(OUT/"repository_state.json",repo)
    if repo["branch"] != "main" or repo["tracked_changes"] or repo["staged_changes"]:
        validation={"lab_id":LAB_ID,"outcome":"REPOSITORY_GATE_FAILURE","head_sha":repo["head_sha"]}
        _write_json(OUT/"validation.json",validation)
        print(json.dumps(validation,indent=2)); return 2

    summaries=[]; checkpoint_rows=[]; failures=[]
    for cluster in BASE.CLUSTERS:
        cid=cluster["id"]
        print(f"[{cid}] full received position+direction state before observer projection")
        try:
            summary,rows,fields,final_arrays=_run_cluster(cluster)
            summaries.append(summary); checkpoint_rows.extend(rows)
            cdir=OUT/"clusters"/cid; cdir.mkdir(parents=True,exist_ok=True)
            _write_json(cdir/"received_state_summary.json",summary)
            _write_csv(cdir/"received_state_checkpoints.csv",rows)
            npz={}
            for step_name,fd in fields.items():
                for name,arr in fd.items(): npz[f"{step_name}__{name}"]=arr
            npz.update(final_arrays)
            np.savez_compressed(cdir/"received_state_checkpoint_fields.npz",**npz)
        except Exception as exc:
            failures.append({"cluster_id":cid,"error":repr(exc)})
            _write_json(OUT/"cluster_failures.json",failures)
            raise

    _write_csv(OUT/"received_state_summary.csv",summaries)
    _write_csv(OUT/"received_state_checkpoint_summary.csv",checkpoint_rows)
    _write_json(OUT/"cluster_failures.json",failures)

    validation={
        "lab_id":LAB_ID,
        "outcome":"Outcome A — G3D FULL RECEIVED POSITION+DIRECTION STATE AUDIT COMPLETE",
        "head_sha":repo["head_sha"],
        "candidate_id":BASE.CANDIDATE_ID,
        "physical_source_representation":BASE.PHYSICAL_SOURCE,
        "geometry_lane":"G3D_LOS_consistent_diagnostic",
        "observer_order":"full_3d_position_plus_direction_then_fixed_xy_projection",
        "observer_channels":"position_and_direction_kept_separate_no_arbitrary_weighting",
        "benchmark_role":"external_morphology_comparison_only",
        "cluster_count_expected":len(BASE.CLUSTERS),
        "cluster_count_completed":len(summaries),
        "checkpoint_steps":list(CHECKPOINTS),
        "metric_identity_tolerance":METRIC_IDENTITY_TOL,
        "projection_identity_tolerance":PROJECTION_IDENTITY_TOL,
        "psd_tolerance":PSD_TOL,
        "all_cluster_g3d_unit_speed_pass":bool(all(r["g3d_unit_speed_pass"] for r in summaries)),
        "all_cluster_first_step_exact_pass":bool(all(r["first_step_exact_pass"] for r in summaries)),
        "mean_los_mag_vs_observed_pearson":float(np.mean([r["los_mag_vs_observed_pearson"] for r in summaries])),
        "mean_final_position_full3d_stretch_gap_vs_observed_pearson":_nanmean_key(summaries,"final_position_full3d_stretch_gap_vs_observed_pearson"),
        "mean_final_direction_full3d_principal_stretch_high_vs_observed_pearson":_nanmean_key(summaries,"final_direction_full3d_principal_stretch_high_vs_observed_pearson"),
        "mean_final_direction_full3d_principal_stretch_low_vs_observed_pearson":_nanmean_key(summaries,"final_direction_full3d_principal_stretch_low_vs_observed_pearson"),
        "mean_final_direction_full3d_area_scale_vs_observed_pearson":_nanmean_key(summaries,"final_direction_full3d_area_scale_vs_observed_pearson"),
        "mean_final_direction_full3d_stretch_gap_vs_observed_pearson":_nanmean_key(summaries,"final_direction_full3d_stretch_gap_vs_observed_pearson"),
        "mean_final_direction_full3d_metric_frobenius_mag_vs_observed_pearson":_nanmean_key(summaries,"final_direction_full3d_metric_frobenius_mag_vs_observed_pearson"),
        "mean_final_direction_projected2d_stretch_gap_vs_observed_pearson":_nanmean_key(summaries,"final_direction_projected2d_stretch_gap_vs_observed_pearson"),
        "mean_final_direction_z_gradient_mag_vs_observed_pearson":_nanmean_key(summaries,"final_direction_z_gradient_mag_vs_observed_pearson"),
        "mean_final_direction_transverse_mag_vs_observed_pearson":_nanmean_key(summaries,"final_direction_transverse_mag_vs_observed_pearson"),
        "mean_final_observer_tangent_angle_mag_vs_observed_pearson":_nanmean_key(summaries,"final_observer_tangent_angle_mag_vs_observed_pearson"),
        "mean_final_observer_tangent_stretch_gap_vs_observed_pearson":_nanmean_key(summaries,"final_observer_tangent_stretch_gap_vs_observed_pearson"),
        "physics_change_authorized":False,
        "candidate_change_authorized":False,
        "production_geometry_change_authorized":False,
        "observer_projection_change_authorized":False,
        "observer_channel_combination_authorized":False,
        "observable_selection_authorized":False,
        "next_experiment_authorized":False,
        "science_interpretation_required":True,
        "duration_seconds":float(time.perf_counter()-started),
    }
    _write_json(OUT/"validation.json",validation)
    _write_json(OUT/"run.json",{
        "lab_id":LAB_ID,"head_sha":repo["head_sha"],"candidate_id":BASE.CANDIDATE_ID,
        "physical_source_representation":BASE.PHYSICAL_SOURCE,
        "geometry_lane":"G3D_LOS_consistent_diagnostic",
        "observer_order":"full_3d_position_plus_direction_then_fixed_xy_projection",
        "note":"Position and direction are measured as separate received channels before projection; no arbitrary weighting or conventional gravitational law is introduced.",
    })
    print(json.dumps(validation,indent=2)); return 0


if __name__ == "__main__":
    raise SystemExit(main())
