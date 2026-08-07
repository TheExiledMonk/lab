#!/usr/bin/env python3
"""PBUF FOUNDATION — G3D ANGULAR MOMENT DECOMPOSITION 001.

Diagnostic-only follow-up to g3d_angular_received_distribution001.

Frozen science:
- candidate PL1_PM1_PS2
- physical source M10_interface_field
- C25 source geometry and photon count
- normalized G3D LOS-consistent propagation
- fixed observer normal n=(0,0,1)
- checkpoints 0,1,5,10,20,40,80,120,159

No conventional gravitational law is introduced into the PBUF pipeline.
Observed kappa remains an external morphology benchmark only.

Question tested here:

  Does the observer's natural angular second moment require both the centroid
  (coherent mean deflection) and the central spread (ray-direction dispersion),
  without introducing any fitted or arbitrary weighting between them?

For each supported source bin, the previous validated angular distribution
provides tangent-plane first and second moments. This lab uses the exact identity

    <|t|^2> = |<t>|^2 + tr(C)

where t=(tx,ty), <t> is the angular centroid, and C is the central angular
covariance. Thus the full angular RMS readout is the unique quadrature sum

    theta_rms = sqrt(theta_centroid^2 + theta_spread^2)

with no fitted coefficient. The full raw second-moment tensor is likewise

    M = <t t^T> = <t><t>^T + C.

The lab measures the centroid contribution, spread contribution, their native
fractions, and the tensor contribution of each term. It does not select a final
observable by benchmark correlation and does not combine position with direction.
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
import pbuf.labs.foundation.g3d_angular_received_distribution001 as ANG
from pbuf.core import los_projection as M14
from pbuf.core import observable_extraction as M16

LAB_ID = "PBUF-FOUNDATION-G3D-ANGULAR-MOMENT-DECOMPOSITION-001"
OUT = ROOT / "runs" / "g3d_angular_moment_decomposition001"
CHECKPOINTS = GEO.CHECKPOINTS
EXPECTED_SUPPORT = GEO.EXPECTED_SUPPORT
UNIT_SPEED_TOL = GEO.UNIT_SPEED_TOL
FIRST_STEP_TOL = GEO.FIRST_STEP_TOL
IDENTITY_TOL = 1e-12
FRACTION_TOL = 1e-12
PARITY_TOL = 1e-14


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
    if isinstance(o, np.floating):
        return float(o)
    if isinstance(o, np.integer):
        return int(o)
    if isinstance(o, np.bool_):
        return bool(o)
    if isinstance(o, Path):
        return str(o)
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
    return ANG._safe_rel_rms(diff, reference)


def _decomposition_fields(base: dict) -> dict:
    mux = np.asarray(base["angular_centroid_tx"], dtype=np.float64)
    muy = np.asarray(base["angular_centroid_ty"], dtype=np.float64)
    cxx = np.asarray(base["angular_cov_xx"], dtype=np.float64)
    cxy = np.asarray(base["angular_cov_xy"], dtype=np.float64)
    cyy = np.asarray(base["angular_cov_yy"], dtype=np.float64)
    mxx = np.asarray(base["angular_second_moment_xx"], dtype=np.float64)
    mxy = np.asarray(base["angular_second_moment_xy"], dtype=np.float64)
    myy = np.asarray(base["angular_second_moment_yy"], dtype=np.float64)

    centroid_energy = mux*mux + muy*muy
    spread_energy = cxx + cyy
    total_energy = mxx + myy
    quadrature = np.sqrt(np.maximum(centroid_energy + spread_energy, 0.0))
    centroid_angle = np.sqrt(np.maximum(centroid_energy, 0.0))
    spread_angle = np.sqrt(np.maximum(spread_energy, 0.0))

    centroid_fraction = np.full_like(total_energy, np.nan)
    spread_fraction = np.full_like(total_energy, np.nan)
    good = np.isfinite(total_energy) & (np.abs(total_energy) > 1e-30)
    centroid_fraction[good] = centroid_energy[good] / total_energy[good]
    spread_fraction[good] = spread_energy[good] / total_energy[good]
    balance = centroid_fraction - spread_fraction

    # Centroid outer-product tensor B = mu mu^T.
    bxx = mux*mux
    bxy = mux*muy
    byy = muy*muy

    centroid_tensor_frob = np.sqrt(bxx*bxx + 2.0*bxy*bxy + byy*byy)
    spread_tensor_frob = np.sqrt(cxx*cxx + 2.0*cxy*cxy + cyy*cyy)
    full_tensor_frob = np.sqrt(mxx*mxx + 2.0*mxy*mxy + myy*myy)

    # Exact missing pieces if one keeps only one part of M.
    centroid_only_residual_frob = spread_tensor_frob
    spread_only_residual_frob = centroid_tensor_frob

    centroid_only_relative_residual = np.full_like(full_tensor_frob, np.nan)
    spread_only_relative_residual = np.full_like(full_tensor_frob, np.nan)
    goodf = np.isfinite(full_tensor_frob) & (full_tensor_frob > 1e-30)
    centroid_only_relative_residual[goodf] = centroid_only_residual_frob[goodf] / full_tensor_frob[goodf]
    spread_only_relative_residual[goodf] = spread_only_residual_frob[goodf] / full_tensor_frob[goodf]

    # Frobenius cross term <B,C> is retained explicitly; tensor norms do not add linearly.
    centroid_spread_tensor_inner = bxx*cxx + 2.0*bxy*cxy + byy*cyy
    centroid_spread_tensor_cosine = np.full_like(full_tensor_frob, np.nan)
    den = centroid_tensor_frob * spread_tensor_frob
    gd = np.isfinite(den) & (den > 1e-30)
    centroid_spread_tensor_cosine[gd] = centroid_spread_tensor_inner[gd] / den[gd]

    return {
        "centroid_energy": centroid_energy,
        "spread_energy": spread_energy,
        "total_angular_second_moment_energy": total_energy,
        "centroid_angle_mag": centroid_angle,
        "spread_angle_rms": spread_angle,
        "quadrature_total_angle_mag": quadrature,
        "centroid_energy_fraction": centroid_fraction,
        "spread_energy_fraction": spread_fraction,
        "centroid_minus_spread_fraction_balance": balance,
        "centroid_outer_xx": bxx,
        "centroid_outer_xy": bxy,
        "centroid_outer_yy": byy,
        "centroid_tensor_frobenius_mag": centroid_tensor_frob,
        "spread_tensor_frobenius_mag": spread_tensor_frob,
        "full_second_moment_tensor_frobenius_mag": full_tensor_frob,
        "centroid_only_residual_tensor_frobenius_mag": centroid_only_residual_frob,
        "spread_only_residual_tensor_frobenius_mag": spread_only_residual_frob,
        "centroid_only_relative_tensor_residual": centroid_only_relative_residual,
        "spread_only_relative_tensor_residual": spread_only_relative_residual,
        "centroid_spread_tensor_inner_product": centroid_spread_tensor_inner,
        "centroid_spread_tensor_cosine": centroid_spread_tensor_cosine,
    }


def _decomposition_gates(base: dict, f: dict) -> dict:
    total = f["total_angular_second_moment_energy"]
    energy_identity = _safe_rel_rms(f["centroid_energy"] + f["spread_energy"] - total, total)
    quadrature_identity = _safe_rel_rms(
        f["quadrature_total_angle_mag"] - base["angular_rms_angle_mag"],
        base["angular_rms_angle_mag"],
    )

    fraction_sum = f["centroid_energy_fraction"] + f["spread_energy_fraction"]
    finite_frac = np.isfinite(fraction_sum)
    fraction_sum_err = float(np.max(np.abs(fraction_sum[finite_frac] - 1.0))) if np.any(finite_frac) else float("nan")

    frac_values = np.concatenate([
        f["centroid_energy_fraction"][np.isfinite(f["centroid_energy_fraction"])],
        f["spread_energy_fraction"][np.isfinite(f["spread_energy_fraction"])],
    ])
    fraction_min = float(np.min(frac_values)) if frac_values.size else float("nan")
    fraction_max = float(np.max(frac_values)) if frac_values.size else float("nan")
    fraction_bounds_pass = bool(
        np.isfinite(fraction_min) and np.isfinite(fraction_max)
        and fraction_min >= -FRACTION_TOL and fraction_max <= 1.0 + FRACTION_TOL
    )

    mxx = base["angular_second_moment_xx"]
    mxy = base["angular_second_moment_xy"]
    myy = base["angular_second_moment_yy"]
    rxx = f["centroid_outer_xx"] + base["angular_cov_xx"]
    rxy = f["centroid_outer_xy"] + base["angular_cov_xy"]
    ryy = f["centroid_outer_yy"] + base["angular_cov_yy"]
    tensor_diff = np.sqrt((rxx-mxx)**2 + 2.0*(rxy-mxy)**2 + (ryy-myy)**2)
    tensor_ref = f["full_second_moment_tensor_frobenius_mag"]
    tensor_identity = _safe_rel_rms(tensor_diff, tensor_ref)

    fro_sq_rhs = (
        f["centroid_tensor_frobenius_mag"]**2
        + f["spread_tensor_frobenius_mag"]**2
        + 2.0*f["centroid_spread_tensor_inner_product"]
    )
    fro_sq_lhs = f["full_second_moment_tensor_frobenius_mag"]**2
    fro_identity = _safe_rel_rms(fro_sq_lhs - fro_sq_rhs, fro_sq_lhs)

    return {
        "centroid_plus_spread_energy_identity_relative_rms_error": energy_identity,
        "quadrature_total_angle_vs_prior_rms_angle_relative_rms_error": quadrature_identity,
        "centroid_plus_spread_fraction_max_abs_error": fraction_sum_err,
        "energy_fraction_min": fraction_min,
        "energy_fraction_max": fraction_max,
        "energy_fraction_bounds_pass": fraction_bounds_pass,
        "second_moment_tensor_decomposition_relative_rms_error": tensor_identity,
        "second_moment_tensor_frobenius_identity_relative_rms_error": fro_identity,
    }


def _require_gates(cid: str, k: int, gates: dict) -> None:
    if gates["centroid_plus_spread_energy_identity_relative_rms_error"] > IDENTITY_TOL:
        raise RuntimeError(f"{cid}: angular energy decomposition identity failed at step {k}")
    if gates["quadrature_total_angle_vs_prior_rms_angle_relative_rms_error"] > PARITY_TOL:
        raise RuntimeError(f"{cid}: quadrature/RMS parity failed at step {k}")
    if gates["centroid_plus_spread_fraction_max_abs_error"] > FRACTION_TOL:
        raise RuntimeError(f"{cid}: centroid/spread fraction sum failed at step {k}")
    if not gates["energy_fraction_bounds_pass"]:
        raise RuntimeError(f"{cid}: centroid/spread fraction bounds failed at step {k}")
    if gates["second_moment_tensor_decomposition_relative_rms_error"] > IDENTITY_TOL:
        raise RuntimeError(f"{cid}: second-moment tensor decomposition failed at step {k}")
    if gates["second_moment_tensor_frobenius_identity_relative_rms_error"] > IDENTITY_TOL:
        raise RuntimeError(f"{cid}: second-moment Frobenius identity failed at step {k}")


def _add_correlations(row: dict, fields: dict, observed, los_mag) -> None:
    names = (
        "centroid_energy",
        "spread_energy",
        "total_angular_second_moment_energy",
        "centroid_angle_mag",
        "spread_angle_rms",
        "quadrature_total_angle_mag",
        "centroid_energy_fraction",
        "spread_energy_fraction",
        "centroid_minus_spread_fraction_balance",
        "centroid_tensor_frobenius_mag",
        "spread_tensor_frobenius_mag",
        "full_second_moment_tensor_frobenius_mag",
        "centroid_only_relative_tensor_residual",
        "spread_only_relative_tensor_residual",
        "centroid_spread_tensor_cosine",
    )
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

    x0,y0,_,_ = BASE._launch_expanded_25pct()
    groups = GEO._source_groups(x0,y0)
    if len(groups) != EXPECTED_SUPPORT:
        raise RuntimeError(f"{cid}: expected {EXPECTED_SUPPORT} source bins, got {len(groups)}")

    checkpoints,g3d = GEO._propagate_g3d(field, BASE.CFG["step"], BASE.CFG["steps"], x0, y0)
    if g3d["max_unit_speed_error"] > UNIT_SPEED_TOL:
        raise RuntimeError(f"{cid}: G3D unit-speed gate failed: {g3d['max_unit_speed_error']}")
    first = GEO._first_step_geometry(field, x0, y0, checkpoints[1], observed, los_mag)
    if not first["first_step_exact_pass"]:
        raise RuntimeError(f"{cid}: first-step exact geometry gate failed")

    rows=[]
    all_fields={}
    for k in CHECKPOINTS:
        base = ANG._angular_distribution_fields(checkpoints[k], groups)
        prior_gates = ANG._moment_gates(base)
        if prior_gates["second_moment_equals_cov_plus_centroid_outer_relative_rms_error"] > ANG.MOMENT_IDENTITY_TOL:
            raise RuntimeError(f"{cid}: inherited angular moment identity failed at step {k}")
        if not prior_gates["covariance_psd_pass"]:
            raise RuntimeError(f"{cid}: inherited covariance PSD failed at step {k}")
        if not prior_gates["direction_mean_vector_bound_pass"]:
            raise RuntimeError(f"{cid}: inherited mean-direction bound failed at step {k}")

        fields = _decomposition_fields(base)
        gates = _decomposition_gates(base, fields)
        _require_gates(cid, k, gates)
        row = {
            "cluster_id": cid,
            "step_index": int(k),
            "propagation_distance": float(k*BASE.CFG["step"]),
            **gates,
        }
        _add_correlations(row, fields, observed, los_mag)
        rows.append(row)
        all_fields[f"step_{k}"] = {**fields}

    final = rows[-1]
    summary = {
        "cluster_id": cid,
        "candidate_id": BASE.CANDIDATE_ID,
        "physical_source_representation": BASE.PHYSICAL_SOURCE,
        "geometry_lane": "G3D_LOS_consistent_diagnostic",
        "observer_readout_lane": "native_angular_first_plus_second_moment_decomposition_no_fitted_weighting",
        "observer_coordinates": "fixed_tangent_plane_tx_vx_over_vz_ty_vy_over_vz_after_full_3d_arrival_state",
        "observer_information_used": "centroid_outer_product_plus_central_covariance_equals_raw_second_moment",
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
    }
    for key,value in final.items():
        if key not in ("cluster_id", "step_index", "propagation_distance"):
            summary[f"final_{key}"] = value

    return summary, rows, all_fields, {"los_mag":los_mag, "observed_benchmark":observed}


def _nanmean_key(summaries, key):
    vals=np.asarray([r[key] for r in summaries], dtype=np.float64)
    return float(np.nanmean(vals)) if np.any(np.isfinite(vals)) else float("nan")


def main() -> int:
    started=time.perf_counter()
    OUT.mkdir(parents=True, exist_ok=True)
    repo=_repo_state()
    _write_json(OUT/"repository_state.json", repo)
    if repo["branch"] != "main" or repo["tracked_changes"] or repo["staged_changes"]:
        validation={"lab_id":LAB_ID,"outcome":"REPOSITORY_GATE_FAILURE","head_sha":repo["head_sha"]}
        _write_json(OUT/"validation.json", validation)
        print(json.dumps(validation, indent=2))
        return 2

    summaries=[]
    checkpoint_rows=[]
    failures=[]
    for cluster in BASE.CLUSTERS:
        cid=cluster["id"]
        print(f"[{cid}] native angular centroid+spread moment decomposition")
        try:
            summary,rows,fields,final_arrays=_run_cluster(cluster)
            summaries.append(summary)
            checkpoint_rows.extend(rows)
            cdir=OUT/"clusters"/cid
            cdir.mkdir(parents=True, exist_ok=True)
            _write_json(cdir/"angular_moment_decomposition_summary.json", summary)
            _write_csv(cdir/"angular_moment_decomposition_checkpoints.csv", rows)
            npz={}
            for step_name,fd in fields.items():
                for name,arr in fd.items():
                    npz[f"{step_name}__{name}"]=arr
            npz.update(final_arrays)
            np.savez_compressed(cdir/"angular_moment_decomposition_checkpoint_fields.npz", **npz)
        except Exception as exc:
            failures.append({"cluster_id":cid,"error":repr(exc)})
            _write_json(OUT/"cluster_failures.json", failures)
            raise

    _write_csv(OUT/"angular_moment_decomposition_summary.csv", summaries)
    _write_csv(OUT/"angular_moment_decomposition_checkpoint_summary.csv", checkpoint_rows)
    _write_json(OUT/"cluster_failures.json", failures)

    validation={
        "lab_id":LAB_ID,
        "outcome":"Outcome A — G3D ANGULAR MOMENT DECOMPOSITION AUDIT COMPLETE",
        "head_sha":repo["head_sha"],
        "candidate_id":BASE.CANDIDATE_ID,
        "physical_source_representation":BASE.PHYSICAL_SOURCE,
        "geometry_lane":"G3D_LOS_consistent_diagnostic",
        "observer_readout_lane":"native_angular_first_plus_second_moment_decomposition_no_fitted_weighting",
        "observer_coordinates":"fixed_tangent_plane_tx_vx_over_vz_ty_vy_over_vz_after_full_3d_arrival_state",
        "observer_information_used":"centroid_outer_product_plus_central_covariance_equals_raw_second_moment",
        "benchmark_role":"external_morphology_comparison_only",
        "cluster_count_expected":len(BASE.CLUSTERS),
        "cluster_count_completed":len(summaries),
        "checkpoint_steps":list(CHECKPOINTS),
        "identity_tolerance":IDENTITY_TOL,
        "fraction_tolerance":FRACTION_TOL,
        "parity_tolerance":PARITY_TOL,
        "all_cluster_g3d_unit_speed_pass":bool(all(r["g3d_unit_speed_pass"] for r in summaries)),
        "all_cluster_first_step_exact_pass":bool(all(r["first_step_exact_pass"] for r in summaries)),
        "all_cluster_final_energy_fraction_bounds_pass":bool(all(r["final_energy_fraction_bounds_pass"] for r in summaries)),
        "mean_los_mag_vs_observed_pearson":float(np.mean([r["los_mag_vs_observed_pearson"] for r in summaries])),
        "mean_final_centroid_angle_mag_vs_observed_pearson":_nanmean_key(summaries,"final_centroid_angle_mag_vs_observed_pearson"),
        "mean_final_spread_angle_rms_vs_observed_pearson":_nanmean_key(summaries,"final_spread_angle_rms_vs_observed_pearson"),
        "mean_final_quadrature_total_angle_mag_vs_observed_pearson":_nanmean_key(summaries,"final_quadrature_total_angle_mag_vs_observed_pearson"),
        "mean_final_centroid_energy_vs_observed_pearson":_nanmean_key(summaries,"final_centroid_energy_vs_observed_pearson"),
        "mean_final_spread_energy_vs_observed_pearson":_nanmean_key(summaries,"final_spread_energy_vs_observed_pearson"),
        "mean_final_total_angular_second_moment_energy_vs_observed_pearson":_nanmean_key(summaries,"final_total_angular_second_moment_energy_vs_observed_pearson"),
        "mean_final_centroid_energy_fraction":float(np.nanmean([r["final_centroid_energy_fraction_rms"] for r in summaries])),
        "mean_final_spread_energy_fraction":float(np.nanmean([r["final_spread_energy_fraction_rms"] for r in summaries])),
        "mean_final_centroid_tensor_frobenius_mag_vs_observed_pearson":_nanmean_key(summaries,"final_centroid_tensor_frobenius_mag_vs_observed_pearson"),
        "mean_final_spread_tensor_frobenius_mag_vs_observed_pearson":_nanmean_key(summaries,"final_spread_tensor_frobenius_mag_vs_observed_pearson"),
        "mean_final_full_second_moment_tensor_frobenius_mag_vs_observed_pearson":_nanmean_key(summaries,"final_full_second_moment_tensor_frobenius_mag_vs_observed_pearson"),
        "mean_final_centroid_only_relative_tensor_residual_rms":float(np.nanmean([r["final_centroid_only_relative_tensor_residual_rms"] for r in summaries])),
        "mean_final_spread_only_relative_tensor_residual_rms":float(np.nanmean([r["final_spread_only_relative_tensor_residual_rms"] for r in summaries])),
        "physics_change_authorized":False,
        "candidate_change_authorized":False,
        "production_geometry_change_authorized":False,
        "observer_projection_change_authorized":False,
        "observer_channel_combination_authorized":False,
        "observable_selection_authorized":False,
        "angular_moment_interpretation_required":True,
        "next_experiment_authorized":False,
        "duration_seconds":float(time.perf_counter()-started),
    }
    _write_json(OUT/"validation.json", validation)
    _write_json(OUT/"run.json", {
        "lab_id":LAB_ID,
        "head_sha":repo["head_sha"],
        "output_directory":str(OUT.relative_to(ROOT)),
        "duration_seconds":validation["duration_seconds"],
    })
    print(json.dumps(validation, indent=2, default=_json_default))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
