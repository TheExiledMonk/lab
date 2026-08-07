#!/usr/bin/env python3
"""PBUF FOUNDATION — G3D ANGULAR RECEIVED DISTRIBUTION 001.

Diagnostic-only follow-up to G3D full received state 001.

Frozen science:
- candidate PL1_PM1_PS2
- physical source M10_interface_field
- C25 source geometry and photon count
- normalized G3D LOS-consistent propagation
- fixed observer normal n=(0,0,1)
- checkpoints 0,1,5,10,20,40,80,120,159

No conventional gravitational law is introduced into the PBUF pipeline.
Observed kappa is an external morphology benchmark only.

Question tested here:

  Does the observer-received angular distribution itself retain structured
  morphology when we gather the per-ray arrival directions directly, before
  taking spatial derivatives or forcing them into a convergence-like scalar?

For each supported source bin, every received unit direction v=(vx,vy,vz) is
converted only after the full 3D state exists into the fixed observer tangent
coordinates

    tx = vx/vz,
    ty = vy/vz.

The lab then measures the complete first/second moment content of the per-ray
angular distribution in each bin:

- angular centroid (mean tx, mean ty)
- mean and RMS angle magnitude
- central angular covariance and its eigenstructure
- raw angular second-moment tensor
- full-3D directional coherence |mean(v)| and dispersion 1-|mean(v)|

No spatial derivative of the angular field is used in this lab. No position
and direction channels are combined. No observable is selected by benchmark
correlation.
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

LAB_ID = "PBUF-FOUNDATION-G3D-ANGULAR-RECEIVED-DISTRIBUTION-001"
OUT = ROOT / "runs" / "g3d_angular_received_distribution001"
CHECKPOINTS = GEO.CHECKPOINTS
EXPECTED_SUPPORT = GEO.EXPECTED_SUPPORT
UNIT_SPEED_TOL = GEO.UNIT_SPEED_TOL
FIRST_STEP_TOL = GEO.FIRST_STEP_TOL
MOMENT_IDENTITY_TOL = 1e-12
PSD_TOL = 1e-12
PARITY_TOL = 1e-14
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
    return float(np.sqrt(np.mean(x * x))) if x.size else float("nan")


def _safe_rel_rms(diff, reference) -> float:
    return OBS3._safe_rel_rms(diff, reference)


def _empty_map() -> np.ndarray:
    return np.full((BASE.OBS_BINS, BASE.OBS_BINS), np.nan, dtype=np.float64)


def _angular_distribution_fields(snap: dict, groups: dict[int, np.ndarray]) -> dict:
    vx = np.asarray(snap["vx"], dtype=np.float64)
    vy = np.asarray(snap["vy"], dtype=np.float64)
    vz = np.asarray(snap["vz"], dtype=np.float64)

    min_abs_vz = float(np.min(np.abs(vz)))
    if min_abs_vz <= VZ_MIN:
        raise RuntimeError(f"observer tangent projection vz too small: {min_abs_vz}")

    tx = vx / vz
    ty = vy / vz
    angle_mag = np.hypot(tx, ty)
    transverse_mag = np.hypot(vx, vy)

    names = (
        "angular_centroid_tx",
        "angular_centroid_ty",
        "angular_centroid_mag",
        "angular_mean_angle_mag",
        "angular_rms_angle_mag",
        "angular_cov_xx",
        "angular_cov_xy",
        "angular_cov_yy",
        "angular_cov_lambda_high",
        "angular_cov_lambda_low",
        "angular_sigma_high",
        "angular_sigma_low",
        "angular_cov_trace",
        "angular_cov_det",
        "angular_cov_frobenius_mag",
        "angular_spread_rms",
        "angular_cov_anisotropy_ratio_low_over_high",
        "angular_second_moment_xx",
        "angular_second_moment_xy",
        "angular_second_moment_yy",
        "angular_second_moment_trace",
        "direction_mean_vx",
        "direction_mean_vy",
        "direction_mean_vz",
        "direction_mean_vector_mag",
        "direction_dispersion_one_minus_mean_vector_mag",
        "direction_mean_transverse_mag",
        "photon_count_per_bin",
    )
    out = {name: _empty_map() for name in names}

    for q, idx in groups.items():
        r, c = divmod(q, BASE.OBS_BINS)
        txi = tx[idx]
        tyi = ty[idx]
        ami = angle_mag[idx]

        mux = float(np.mean(txi))
        muy = float(np.mean(tyi))
        dx = txi - mux
        dy = tyi - muy

        cxx = float(np.mean(dx * dx))
        cxy = float(np.mean(dx * dy))
        cyy = float(np.mean(dy * dy))
        tr = cxx + cyy
        det = cxx * cyy - cxy * cxy
        disc = float(np.sqrt(max((cxx - cyy) ** 2 + 4.0 * cxy * cxy, 0.0)))
        lam_hi = 0.5 * (tr + disc)
        lam_lo = 0.5 * (tr - disc)
        sig_hi = float(np.sqrt(max(lam_hi, 0.0)))
        sig_lo = float(np.sqrt(max(lam_lo, 0.0)))
        anis = float(sig_lo / sig_hi) if sig_hi > 1e-30 else float("nan")

        mxx = float(np.mean(txi * txi))
        mxy = float(np.mean(txi * tyi))
        myy = float(np.mean(tyi * tyi))

        mvx = float(np.mean(vx[idx]))
        mvy = float(np.mean(vy[idx]))
        mvz = float(np.mean(vz[idx]))
        mean_v_mag = float(np.sqrt(mvx * mvx + mvy * mvy + mvz * mvz))

        out["angular_centroid_tx"][r, c] = mux
        out["angular_centroid_ty"][r, c] = muy
        out["angular_centroid_mag"][r, c] = float(np.hypot(mux, muy))
        out["angular_mean_angle_mag"][r, c] = float(np.mean(ami))
        out["angular_rms_angle_mag"][r, c] = float(np.sqrt(np.mean(ami * ami)))
        out["angular_cov_xx"][r, c] = cxx
        out["angular_cov_xy"][r, c] = cxy
        out["angular_cov_yy"][r, c] = cyy
        out["angular_cov_lambda_high"][r, c] = lam_hi
        out["angular_cov_lambda_low"][r, c] = lam_lo
        out["angular_sigma_high"][r, c] = sig_hi
        out["angular_sigma_low"][r, c] = sig_lo
        out["angular_cov_trace"][r, c] = tr
        out["angular_cov_det"][r, c] = det
        out["angular_cov_frobenius_mag"][r, c] = float(np.sqrt(cxx*cxx + 2.0*cxy*cxy + cyy*cyy))
        out["angular_spread_rms"][r, c] = float(np.sqrt(max(tr, 0.0)))
        out["angular_cov_anisotropy_ratio_low_over_high"][r, c] = anis
        out["angular_second_moment_xx"][r, c] = mxx
        out["angular_second_moment_xy"][r, c] = mxy
        out["angular_second_moment_yy"][r, c] = myy
        out["angular_second_moment_trace"][r, c] = mxx + myy
        out["direction_mean_vx"][r, c] = mvx
        out["direction_mean_vy"][r, c] = mvy
        out["direction_mean_vz"][r, c] = mvz
        out["direction_mean_vector_mag"][r, c] = mean_v_mag
        out["direction_dispersion_one_minus_mean_vector_mag"][r, c] = 1.0 - mean_v_mag
        out["direction_mean_transverse_mag"][r, c] = float(np.mean(transverse_mag[idx]))
        out["photon_count_per_bin"][r, c] = float(len(idx))

    return out


def _moment_gates(fields: dict) -> dict:
    mux = fields["angular_centroid_tx"]
    muy = fields["angular_centroid_ty"]
    cxx = fields["angular_cov_xx"]
    cxy = fields["angular_cov_xy"]
    cyy = fields["angular_cov_yy"]
    mxx = fields["angular_second_moment_xx"]
    mxy = fields["angular_second_moment_xy"]
    myy = fields["angular_second_moment_yy"]

    dxx = mxx - (cxx + mux * mux)
    dxy = mxy - (cxy + mux * muy)
    dyy = myy - (cyy + muy * muy)
    diff = np.sqrt(dxx*dxx + 2.0*dxy*dxy + dyy*dyy)
    ref = np.sqrt(mxx*mxx + 2.0*mxy*mxy + myy*myy)
    second_central_identity = _safe_rel_rms(diff, ref)

    rms_sq = fields["angular_rms_angle_mag"] ** 2
    second_trace = fields["angular_second_moment_trace"]
    rms_identity = _safe_rel_rms(rms_sq - second_trace, second_trace)

    lam_hi = fields["angular_cov_lambda_high"]
    lam_lo = fields["angular_cov_lambda_low"]
    tr = fields["angular_cov_trace"]
    det = fields["angular_cov_det"]
    eig_trace_identity = _safe_rel_rms(lam_hi + lam_lo - tr, tr)
    eig_det_identity = _safe_rel_rms(lam_hi * lam_lo - det, det)

    mask = np.isfinite(lam_lo)
    min_eig = float(np.nanmin(lam_lo[mask])) if np.any(mask) else float("nan")

    mean_v = fields["direction_mean_vector_mag"]
    finite_mv = mean_v[np.isfinite(mean_v)]
    max_mean_vector_mag = float(np.max(finite_mv)) if finite_mv.size else float("nan")
    min_mean_vector_mag = float(np.min(finite_mv)) if finite_mv.size else float("nan")

    return {
        "second_moment_equals_cov_plus_centroid_outer_relative_rms_error": second_central_identity,
        "rms_angle_squared_equals_second_moment_trace_relative_rms_error": rms_identity,
        "covariance_eigen_trace_identity_relative_rms_error": eig_trace_identity,
        "covariance_eigen_det_identity_relative_rms_error": eig_det_identity,
        "covariance_min_eigenvalue": min_eig,
        "covariance_psd_pass": bool(min_eig >= -PSD_TOL),
        "direction_mean_vector_mag_min": min_mean_vector_mag,
        "direction_mean_vector_mag_max": max_mean_vector_mag,
        "direction_mean_vector_bound_pass": bool(max_mean_vector_mag <= 1.0 + PSD_TOL and min_mean_vector_mag >= -PSD_TOL),
    }


def _add_correlations(row: dict, fields: dict, observed, los_mag) -> None:
    names = (
        "angular_centroid_mag",
        "angular_mean_angle_mag",
        "angular_rms_angle_mag",
        "angular_sigma_high",
        "angular_sigma_low",
        "angular_cov_trace",
        "angular_cov_det",
        "angular_cov_frobenius_mag",
        "angular_spread_rms",
        "angular_cov_anisotropy_ratio_low_over_high",
        "angular_second_moment_trace",
        "direction_mean_vector_mag",
        "direction_dispersion_one_minus_mean_vector_mag",
        "direction_mean_transverse_mag",
        "photon_count_per_bin",
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


def _checkpoint(cid, k, snap, groups, observed, los_mag):
    fields = _angular_distribution_fields(snap, groups)
    gates = _moment_gates(fields)

    # Exact parity against the already-used direct per-bin means. This is not a
    # science comparison; it only protects the aggregation implementation.
    vx = np.asarray(snap["vx"], dtype=np.float64)
    vy = np.asarray(snap["vy"], dtype=np.float64)
    vz = np.asarray(snap["vz"], dtype=np.float64)
    tx = vx / vz
    ty = vy / vz
    ref_trans = GEO._mean_map(np.hypot(vx, vy), groups)
    ref_angle = GEO._mean_map(np.hypot(tx, ty), groups)
    trans_parity = _safe_rel_rms(fields["direction_mean_transverse_mag"] - ref_trans, ref_trans)
    angle_parity = _safe_rel_rms(fields["angular_mean_angle_mag"] - ref_angle, ref_angle)

    row = {
        "cluster_id": cid,
        "step_index": int(k),
        "propagation_distance": float(k * BASE.CFG["step"]),
        **gates,
        "direction_mean_transverse_aggregation_parity_relative_rms_error": trans_parity,
        "angular_mean_angle_aggregation_parity_relative_rms_error": angle_parity,
    }
    _add_correlations(row, fields, observed, los_mag)
    return row, fields


def _require_checkpoint_gates(cid: str, k: int, row: dict) -> None:
    for key in (
        "second_moment_equals_cov_plus_centroid_outer_relative_rms_error",
        "rms_angle_squared_equals_second_moment_trace_relative_rms_error",
        "covariance_eigen_trace_identity_relative_rms_error",
        "covariance_eigen_det_identity_relative_rms_error",
    ):
        value = row[key]
        if not np.isfinite(value) or value > MOMENT_IDENTITY_TOL:
            raise RuntimeError(f"{cid}: {key} failed at step {k}: {value}")
    if not row["covariance_psd_pass"]:
        raise RuntimeError(f"{cid}: angular covariance PSD failed at step {k}")
    if not row["direction_mean_vector_bound_pass"]:
        raise RuntimeError(f"{cid}: mean direction-vector norm bound failed at step {k}")
    for key in (
        "direction_mean_transverse_aggregation_parity_relative_rms_error",
        "angular_mean_angle_aggregation_parity_relative_rms_error",
    ):
        value = row[key]
        if not np.isfinite(value) or value > PARITY_TOL:
            raise RuntimeError(f"{cid}: {key} failed at step {k}: {value}")


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
        raise RuntimeError(f"{cid}: first-step exact geometry gate failed")

    rows: list[dict] = []
    all_fields: dict[str, dict] = {}
    for k in CHECKPOINTS:
        row, fields = _checkpoint(cid, k, checkpoints[k], groups, observed, los_mag)
        _require_checkpoint_gates(cid, k, row)
        rows.append(row)
        all_fields[f"step_{k}"] = fields

    final = rows[-1]
    summary = {
        "cluster_id": cid,
        "candidate_id": BASE.CANDIDATE_ID,
        "physical_source_representation": BASE.PHYSICAL_SOURCE,
        "geometry_lane": "G3D_LOS_consistent_diagnostic",
        "observer_readout_lane": "per_ray_arrival_direction_distribution_before_spatial_differentiation",
        "observer_coordinates": "fixed_tangent_plane_tx_vx_over_vz_ty_vy_over_vz_after_full_3d_arrival_state",
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

    gate_keys = (
        "second_moment_equals_cov_plus_centroid_outer_relative_rms_error",
        "rms_angle_squared_equals_second_moment_trace_relative_rms_error",
        "covariance_eigen_trace_identity_relative_rms_error",
        "covariance_eigen_det_identity_relative_rms_error",
        "covariance_min_eigenvalue",
        "covariance_psd_pass",
        "direction_mean_vector_mag_min",
        "direction_mean_vector_mag_max",
        "direction_mean_vector_bound_pass",
        "direction_mean_transverse_aggregation_parity_relative_rms_error",
        "angular_mean_angle_aggregation_parity_relative_rms_error",
    )
    for key in gate_keys:
        summary[f"final_{key}"] = final[key]

    metric_names = (
        "angular_centroid_mag",
        "angular_mean_angle_mag",
        "angular_rms_angle_mag",
        "angular_sigma_high",
        "angular_sigma_low",
        "angular_cov_trace",
        "angular_cov_det",
        "angular_cov_frobenius_mag",
        "angular_spread_rms",
        "angular_cov_anisotropy_ratio_low_over_high",
        "angular_second_moment_trace",
        "direction_mean_vector_mag",
        "direction_dispersion_one_minus_mean_vector_mag",
        "direction_mean_transverse_mag",
        "photon_count_per_bin",
    )
    for name in metric_names:
        for suffix in (
            "vs_observed_pearson", "vs_observed_spearman", "vs_observed_count",
            "vs_los_mag_pearson", "vs_los_mag_spearman", "vs_los_mag_count", "rms",
        ):
            summary[f"final_{name}_{suffix}"] = final[f"{name}_{suffix}"]

    return summary, rows, all_fields, {"los_mag": los_mag, "observed_benchmark": observed}


def _nanmean_key(summaries: list[dict], key: str) -> float:
    vals = np.asarray([r[key] for r in summaries], dtype=np.float64)
    return float(np.nanmean(vals)) if np.any(np.isfinite(vals)) else float("nan")


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

    summaries: list[dict] = []
    checkpoint_rows: list[dict] = []
    failures: list[dict] = []

    for cluster in BASE.CLUSTERS:
        cid = cluster["id"]
        print(f"[{cid}] received angular distribution moments")
        try:
            summary, rows, fields, final_arrays = _run_cluster(cluster)
            summaries.append(summary)
            checkpoint_rows.extend(rows)
            cdir = OUT / "clusters" / cid
            cdir.mkdir(parents=True, exist_ok=True)
            _write_json(cdir / "angular_distribution_summary.json", summary)
            _write_csv(cdir / "angular_distribution_checkpoints.csv", rows)
            npz = {}
            for step_name, fd in fields.items():
                for name, arr in fd.items():
                    npz[f"{step_name}__{name}"] = arr
            npz.update(final_arrays)
            np.savez_compressed(cdir / "angular_distribution_checkpoint_fields.npz", **npz)
        except Exception as exc:
            failures.append({"cluster_id": cid, "error": repr(exc)})
            _write_json(OUT / "cluster_failures.json", failures)
            raise

    _write_csv(OUT / "angular_distribution_summary.csv", summaries)
    _write_csv(OUT / "angular_distribution_checkpoint_summary.csv", checkpoint_rows)
    _write_json(OUT / "cluster_failures.json", failures)

    validation = {
        "lab_id": LAB_ID,
        "outcome": "Outcome A — G3D ANGULAR RECEIVED DISTRIBUTION AUDIT COMPLETE",
        "head_sha": repo["head_sha"],
        "candidate_id": BASE.CANDIDATE_ID,
        "physical_source_representation": BASE.PHYSICAL_SOURCE,
        "geometry_lane": "G3D_LOS_consistent_diagnostic",
        "observer_readout_lane": "per_ray_arrival_direction_distribution_before_spatial_differentiation",
        "observer_coordinates": "fixed_tangent_plane_tx_vx_over_vz_ty_vy_over_vz_after_full_3d_arrival_state",
        "observer_information_used": "first_and_second_moments_of_full_per_ray_direction_distribution_per_source_bin",
        "benchmark_role": "external_morphology_comparison_only",
        "cluster_count_expected": len(BASE.CLUSTERS),
        "cluster_count_completed": len(summaries),
        "checkpoint_steps": list(CHECKPOINTS),
        "moment_identity_tolerance": MOMENT_IDENTITY_TOL,
        "psd_tolerance": PSD_TOL,
        "aggregation_parity_tolerance": PARITY_TOL,
        "all_cluster_g3d_unit_speed_pass": bool(all(r["g3d_unit_speed_pass"] for r in summaries)),
        "all_cluster_first_step_exact_pass": bool(all(r["first_step_exact_pass"] for r in summaries)),
        "all_cluster_final_covariance_psd_pass": bool(all(r["final_covariance_psd_pass"] for r in summaries)),
        "all_cluster_final_direction_mean_vector_bound_pass": bool(all(r["final_direction_mean_vector_bound_pass"] for r in summaries)),
        "mean_los_mag_vs_observed_pearson": float(np.mean([r["los_mag_vs_observed_pearson"] for r in summaries])),
        "mean_final_angular_centroid_mag_vs_observed_pearson": _nanmean_key(summaries, "final_angular_centroid_mag_vs_observed_pearson"),
        "mean_final_angular_mean_angle_mag_vs_observed_pearson": _nanmean_key(summaries, "final_angular_mean_angle_mag_vs_observed_pearson"),
        "mean_final_angular_rms_angle_mag_vs_observed_pearson": _nanmean_key(summaries, "final_angular_rms_angle_mag_vs_observed_pearson"),
        "mean_final_angular_sigma_high_vs_observed_pearson": _nanmean_key(summaries, "final_angular_sigma_high_vs_observed_pearson"),
        "mean_final_angular_sigma_low_vs_observed_pearson": _nanmean_key(summaries, "final_angular_sigma_low_vs_observed_pearson"),
        "mean_final_angular_spread_rms_vs_observed_pearson": _nanmean_key(summaries, "final_angular_spread_rms_vs_observed_pearson"),
        "mean_final_angular_cov_frobenius_mag_vs_observed_pearson": _nanmean_key(summaries, "final_angular_cov_frobenius_mag_vs_observed_pearson"),
        "mean_final_angular_second_moment_trace_vs_observed_pearson": _nanmean_key(summaries, "final_angular_second_moment_trace_vs_observed_pearson"),
        "mean_final_direction_mean_vector_mag_vs_observed_pearson": _nanmean_key(summaries, "final_direction_mean_vector_mag_vs_observed_pearson"),
        "mean_final_direction_dispersion_vs_observed_pearson": _nanmean_key(summaries, "final_direction_dispersion_one_minus_mean_vector_mag_vs_observed_pearson"),
        "mean_final_direction_mean_transverse_mag_vs_observed_pearson": _nanmean_key(summaries, "final_direction_mean_transverse_mag_vs_observed_pearson"),
        "physics_change_authorized": False,
        "candidate_change_authorized": False,
        "production_geometry_change_authorized": False,
        "observer_projection_change_authorized": False,
        "observer_channel_combination_authorized": False,
        "observable_selection_authorized": False,
        "angular_distribution_interpretation_required": True,
        "next_experiment_authorized": False,
        "duration_seconds": float(time.perf_counter() - started),
    }

    _write_json(OUT / "validation.json", validation)
    _write_json(OUT / "run.json", {
        "lab_id": LAB_ID,
        "head_sha": repo["head_sha"],
        "output_directory": str(OUT.relative_to(ROOT)),
        "cluster_count": len(summaries),
        "duration_seconds": validation["duration_seconds"],
    })

    print(json.dumps(validation, indent=2, default=_json_default))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
