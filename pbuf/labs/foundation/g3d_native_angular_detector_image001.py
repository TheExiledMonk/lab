#!/usr/bin/env python3
"""PBUF FOUNDATION — G3D NATIVE ANGULAR DETECTOR IMAGE 001.

Diagnostic-only follow-up to the validated angular moment and higher-moment
observer audits.

Frozen science:
- candidate PL1_PM1_PS2
- physical source M10_interface_field
- exact C25 source geometry / photon count
- normalized G3D LOS-consistent propagation
- fixed observer normal n=(0,0,1)
- tangent coordinates only after full 3D arrival: tx=vx/vz, ty=vy/vz
- checkpoints 0,1,5,10,20,40,80,120,159

No conventional gravitational law is introduced. Observed kappa is never used
to construct, tune, center, scale, weight, or select the detector image.

Question:

  If the observer forms an image directly from the received angular rays, does
  replacing each empirical per-source-bin angular distribution by a deterministic
  control that preserves only its measured centroid mu and covariance C materially
  change the resulting detector image?

Native detector operation:

  The observer image is simply a 2-D histogram of received tangent-plane ray
  directions (tx,ty). Every empirical ray contributes unit flux. The moment-only
  control uses exactly the same number of rays in every source bin and exactly
  preserves that bin's centroid and covariance before both lanes pass through the
  same histogram operation.

Moment-only control:

  For a source bin containing N rays with empirical centroid mu and covariance C,
  construct N deterministic points on a uniformly sampled unit circle and map them
  through sqrt(2 C):

      t_j = mu + sqrt(2 C) u_j,
      u_j = (cos(2 pi j/N), sin(2 pi j/N)).

  For N>=3 the discrete circle has exact zero mean and covariance I/2, so the
  transformed points reproduce mu and C to floating-point precision. This is a
  control distribution only. It is NOT asserted to be Gaussian or physical.

Detector range:

  For each cluster, one fixed square detector range centered on the observer optical
  axis is determined without benchmark input from the maximum absolute tx/ty value
  occurring in either lane over every committed checkpoint. That same range is then
  used at 16x16, 32x32, and 64x64 resolution for all checkpoints and both lanes.

No smoothing, fitting, rescaling, adaptive benchmark alignment, or observable
selection is performed. Similarity metrics are descriptive and never execution
criteria.
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

LAB_ID = "PBUF-FOUNDATION-G3D-NATIVE-ANGULAR-DETECTOR-IMAGE-001"
OUT = ROOT / "runs" / "g3d_native_angular_detector_image001"
CHECKPOINTS = GEO.CHECKPOINTS
EXPECTED_SUPPORT = GEO.EXPECTED_SUPPORT
UNIT_SPEED_TOL = GEO.UNIT_SPEED_TOL
FIRST_STEP_TOL = GEO.FIRST_STEP_TOL
MOMENT_PARITY_TOL = 1e-12
FLUX_PARITY_TOL = 1e-12
VZ_MIN = 1e-12
DETECTOR_RESOLUTIONS = (16, 32, 64)
PRIMARY_RESOLUTION = 64
RANGE_PAD = 1.0 + 1e-12


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


def _cov2(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    mx = float(np.mean(x))
    my = float(np.mean(y))
    dx = x - mx
    dy = y - my
    return np.array(
        [
            [float(np.mean(dx * dx)), float(np.mean(dx * dy))],
            [float(np.mean(dx * dy)), float(np.mean(dy * dy))],
        ],
        dtype=np.float64,
    )


def _moment_only_surrogate(
    tx: np.ndarray,
    ty: np.ndarray,
    groups: dict[int, np.ndarray],
) -> tuple[np.ndarray, np.ndarray, dict]:
    tx = np.asarray(tx, dtype=np.float64)
    ty = np.asarray(ty, dtype=np.float64)
    sx = np.empty_like(tx)
    sy = np.empty_like(ty)

    max_centroid_error = 0.0
    max_cov_frobenius_error = 0.0
    max_cov_negative_eigenvalue = 0.0
    min_group_size = 10**9
    max_group_size = 0

    for q, idx in groups.items():
        del q
        n = int(len(idx))
        if n < 3:
            raise RuntimeError(f"moment-only surrogate requires >=3 rays per source bin, got {n}")
        min_group_size = min(min_group_size, n)
        max_group_size = max(max_group_size, n)

        ex = tx[idx]
        ey = ty[idx]
        mu = np.array([float(np.mean(ex)), float(np.mean(ey))], dtype=np.float64)
        C = _cov2(ex, ey)
        evals, evecs = np.linalg.eigh(C)
        min_eval = float(np.min(evals))
        if min_eval < max_cov_negative_eigenvalue:
            max_cov_negative_eigenvalue = min_eval
        if min_eval < -ANG.PSD_TOL:
            raise RuntimeError(f"empirical angular covariance not PSD: min eigenvalue={min_eval}")
        evals = np.maximum(evals, 0.0)

        # A A^T / 2 = C, hence A = evecs diag(sqrt(2 lambda)).
        A = evecs @ np.diag(np.sqrt(2.0 * evals))
        theta = (2.0 * np.pi / n) * np.arange(n, dtype=np.float64)
        unit = np.column_stack((np.cos(theta), np.sin(theta)))
        pts = mu[None, :] + unit @ A.T
        sx[idx] = pts[:, 0]
        sy[idx] = pts[:, 1]

        smu = np.array([float(np.mean(pts[:, 0])), float(np.mean(pts[:, 1]))])
        SC = _cov2(pts[:, 0], pts[:, 1])
        max_centroid_error = max(max_centroid_error, float(np.linalg.norm(smu - mu)))
        max_cov_frobenius_error = max(
            max_cov_frobenius_error,
            float(np.linalg.norm(SC - C, ord="fro")),
        )

    empirical_mu = np.array([float(np.mean(tx)), float(np.mean(ty))])
    surrogate_mu = np.array([float(np.mean(sx)), float(np.mean(sy))])
    empirical_cov = _cov2(tx, ty)
    surrogate_cov = _cov2(sx, sy)

    controls = {
        "source_bin_count": int(len(groups)),
        "source_group_size_min": int(min_group_size),
        "source_group_size_max": int(max_group_size),
        "source_max_centroid_vector_error": max_centroid_error,
        "source_max_covariance_frobenius_error": max_cov_frobenius_error,
        "empirical_covariance_most_negative_eigenvalue": max_cov_negative_eigenvalue,
        "global_centroid_vector_error": float(np.linalg.norm(surrogate_mu - empirical_mu)),
        "global_covariance_frobenius_error": float(np.linalg.norm(surrogate_cov - empirical_cov, ord="fro")),
        "surrogate_moment_parity_pass": bool(
            max_centroid_error <= MOMENT_PARITY_TOL
            and max_cov_frobenius_error <= MOMENT_PARITY_TOL
        ),
    }
    return sx, sy, controls


def _detector_hist(tx, ty, bins: int, half_range: float) -> np.ndarray:
    h, _, _ = np.histogram2d(
        np.asarray(tx, dtype=np.float64),
        np.asarray(ty, dtype=np.float64),
        bins=bins,
        range=[[-half_range, half_range], [-half_range, half_range]],
    )
    return np.asarray(h, dtype=np.float64)


def _probability_image(h: np.ndarray) -> np.ndarray:
    h = np.asarray(h, dtype=np.float64)
    total = float(np.sum(h))
    if total <= 0.0:
        raise RuntimeError("detector histogram has zero total flux")
    return h / total


def _js_divergence(p: np.ndarray, q: np.ndarray) -> float:
    p = np.asarray(p, dtype=np.float64).ravel()
    q = np.asarray(q, dtype=np.float64).ravel()
    m = 0.5 * (p + q)
    out = 0.0
    mp = p > 0.0
    mq = q > 0.0
    out += 0.5 * float(np.sum(p[mp] * np.log(p[mp] / m[mp])))
    out += 0.5 * float(np.sum(q[mq] * np.log(q[mq] / m[mq])))
    return out


def _image_metrics(emp: np.ndarray, sur: np.ndarray) -> dict:
    emp = np.asarray(emp, dtype=np.float64)
    sur = np.asarray(sur, dtype=np.float64)
    flux_emp = float(np.sum(emp))
    flux_sur = float(np.sum(sur))
    flux_rel = abs(flux_sur - flux_emp) / max(abs(flux_emp), 1.0)
    p = _probability_image(emp)
    q = _probability_image(sur)
    diff = p - q
    tv = 0.5 * float(np.sum(np.abs(diff)))
    js = _js_divergence(p, q)
    l2 = float(np.linalg.norm(diff.ravel()))
    denom = float(np.linalg.norm(p.ravel()))
    rel_l2 = l2 / denom if denom > 0.0 else float("nan")
    pearson, spearman, count = _corr(p, q)
    return {
        "empirical_flux": flux_emp,
        "moment_control_flux": flux_sur,
        "flux_relative_error": flux_rel,
        "flux_parity_pass": bool(flux_rel <= FLUX_PARITY_TOL),
        "total_variation_distance": tv,
        "jensen_shannon_divergence_nats": js,
        "normalized_l2_distance": rel_l2,
        "image_pearson": pearson,
        "image_spearman": spearman,
        "image_correlation_count": count,
        "empirical_occupied_bins": int(np.count_nonzero(emp)),
        "moment_control_occupied_bins": int(np.count_nonzero(sur)),
        "occupied_bin_intersection": int(np.count_nonzero((emp > 0) & (sur > 0))),
    }


def _prepare_checkpoint_states(checkpoints: dict, groups: dict[int, np.ndarray]) -> tuple[dict, float]:
    prepared = {}
    max_abs = 0.0
    for k in CHECKPOINTS:
        snap = checkpoints[k]
        vx = np.asarray(snap["vx"], dtype=np.float64)
        vy = np.asarray(snap["vy"], dtype=np.float64)
        vz = np.asarray(snap["vz"], dtype=np.float64)
        min_abs_vz = float(np.min(np.abs(vz)))
        if min_abs_vz <= VZ_MIN:
            raise RuntimeError(f"checkpoint {k}: observer tangent projection vz too small: {min_abs_vz}")
        tx = vx / vz
        ty = vy / vz
        sx, sy, controls = _moment_only_surrogate(tx, ty, groups)
        local = float(np.max(np.abs(np.concatenate((tx, ty, sx, sy)))))
        max_abs = max(max_abs, local)
        prepared[k] = {
            "tx": tx,
            "ty": ty,
            "sx": sx,
            "sy": sy,
            "controls": controls,
        }
    if max_abs <= 0.0:
        raise RuntimeError("all committed checkpoints have zero angular detector range")
    return prepared, max_abs * RANGE_PAD


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

    prepared, half_range = _prepare_checkpoint_states(checkpoints, groups)
    rows = []
    npz = {}

    for k in CHECKPOINTS:
        st = prepared[k]
        controls = st["controls"]
        if not controls["surrogate_moment_parity_pass"]:
            raise RuntimeError(f"{cid}: moment-only surrogate parity failed at step {k}: {controls}")

        base_ang = ANG._angular_distribution_fields(checkpoints[k], groups)
        prior_gates = ANG._moment_gates(base_ang)
        if prior_gates["second_moment_equals_cov_plus_centroid_outer_relative_rms_error"] > ANG.MOMENT_IDENTITY_TOL:
            raise RuntimeError(f"{cid}: inherited angular moment identity failed at step {k}")
        if not prior_gates["covariance_psd_pass"]:
            raise RuntimeError(f"{cid}: inherited covariance PSD failed at step {k}")
        if not prior_gates["direction_mean_vector_bound_pass"]:
            raise RuntimeError(f"{cid}: inherited mean-direction bound failed at step {k}")

        row = {
            "cluster_id": cid,
            "step_index": int(k),
            "propagation_distance": float(k * BASE.CFG["step"]),
            "detector_half_range": half_range,
            **controls,
        }

        for res in DETECTOR_RESOLUTIONS:
            eh = _detector_hist(st["tx"], st["ty"], res, half_range)
            sh = _detector_hist(st["sx"], st["sy"], res, half_range)
            metrics = _image_metrics(eh, sh)
            for key, value in metrics.items():
                row[f"r{res}_{key}"] = value
            if not metrics["flux_parity_pass"]:
                raise RuntimeError(f"{cid}: detector flux parity failed at step {k}, resolution {res}")
            npz[f"step_{k}__r{res}__empirical"] = eh
            npz[f"step_{k}__r{res}__moment_control"] = sh

        rows.append(row)

    final_ang = ANG._angular_distribution_fields(checkpoints[CHECKPOINTS[-1]], groups)
    prior_rms_p, prior_rms_s, prior_rms_n = _corr(final_ang["angular_rms_angle_mag"], observed)
    los_p, los_s, _ = _corr(los_mag, observed)

    final = rows[-1]
    summary = {
        "cluster_id": cid,
        "candidate_id": BASE.CANDIDATE_ID,
        "physical_source_representation": BASE.PHYSICAL_SOURCE,
        "geometry_lane": "G3D_LOS_consistent_diagnostic",
        "observer_readout_lane": "native_tangent_direction_histogram_empirical_vs_mu_cov_moment_control",
        "observer_coordinates": "fixed_tangent_plane_tx_vx_over_vz_ty_vy_over_vz_after_full_3d_arrival_state",
        "observer_information_full_lane": "all_empirical_received_ray_directions",
        "observer_information_control_lane": "per_source_bin_centroid_and_covariance_only_deterministic_equal_count_circle_control",
        "detector_operation": "unit_flux_2d_histogram_of_received_tangent_directions_no_smoothing",
        "detector_range_rule": "per_cluster_fixed_symmetric_zero_centered_range_from_union_of_both_lanes_all_committed_checkpoints_no_benchmark_input",
        "detector_resolutions": list(DETECTOR_RESOLUTIONS),
        "primary_detector_resolution": PRIMARY_RESOLUTION,
        "benchmark_role": "external_morphology_comparison_only_not_used_for_detector_construction_or_lane_similarity",
        "moment_control_role": "diagnostic_mu_cov_sufficiency_control_not_physical_distribution_model",
        "n_photons": int(len(x0)),
        "source_supported_bins": int(len(groups)),
        "checkpoint_count": len(rows),
        "detector_half_range": half_range,
        "g3d_unit_speed_max_error": float(g3d["max_unit_speed_error"]),
        "g3d_unit_speed_pass": bool(g3d["max_unit_speed_error"] <= UNIT_SPEED_TOL),
        "first_step_exact_max_vector_error": first["first_step_exact_max_vector_error"],
        "first_step_exact_pass": first["first_step_exact_pass"],
        "los_mag_vs_observed_pearson": los_p,
        "los_mag_vs_observed_spearman": los_s,
        "prior_final_angular_rms_angle_mag_vs_observed_pearson": prior_rms_p,
        "prior_final_angular_rms_angle_mag_vs_observed_spearman": prior_rms_s,
        "prior_final_angular_rms_angle_mag_vs_observed_count": prior_rms_n,
    }
    for key, value in final.items():
        if key not in ("cluster_id", "step_index", "propagation_distance"):
            summary[f"final_{key}"] = value

    npz["observed_benchmark_reference_only"] = observed
    npz["los_magnitude_reference_only"] = los_mag
    return summary, rows, npz


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

    summaries = []
    checkpoint_rows = []
    failures = []

    for cluster in BASE.CLUSTERS:
        cid = cluster["id"]
        print(f"[{cid}] native angular detector image: empirical vs mu,C moment control")
        try:
            summary, rows, npz = _run_cluster(cluster)
            summaries.append(summary)
            checkpoint_rows.extend(rows)
            cdir = OUT / "clusters" / cid
            cdir.mkdir(parents=True, exist_ok=True)
            _write_json(cdir / "native_angular_detector_summary.json", summary)
            _write_csv(cdir / "native_angular_detector_checkpoints.csv", rows)
            np.savez_compressed(cdir / "native_angular_detector_checkpoint_images.npz", **npz)
        except Exception as exc:
            failures.append({"cluster_id": cid, "error": repr(exc)})
            _write_json(OUT / "cluster_failures.json", failures)
            raise

    _write_csv(OUT / "native_angular_detector_summary.csv", summaries)
    _write_csv(OUT / "native_angular_detector_checkpoint_summary.csv", checkpoint_rows)
    _write_json(OUT / "cluster_failures.json", failures)

    validation = {
        "lab_id": LAB_ID,
        "outcome": "Outcome A — G3D NATIVE ANGULAR DETECTOR IMAGE AUDIT COMPLETE",
        "head_sha": repo["head_sha"],
        "candidate_id": BASE.CANDIDATE_ID,
        "physical_source_representation": BASE.PHYSICAL_SOURCE,
        "geometry_lane": "G3D_LOS_consistent_diagnostic",
        "observer_readout_lane": "native_tangent_direction_histogram_empirical_vs_mu_cov_moment_control",
        "observer_coordinates": "fixed_tangent_plane_tx_vx_over_vz_ty_vy_over_vz_after_full_3d_arrival_state",
        "observer_information_full_lane": "all_empirical_received_ray_directions",
        "observer_information_control_lane": "per_source_bin_centroid_and_covariance_only_deterministic_equal_count_circle_control",
        "detector_operation": "unit_flux_2d_histogram_of_received_tangent_directions_no_smoothing",
        "detector_range_rule": "per_cluster_fixed_symmetric_zero_centered_range_from_union_of_both_lanes_all_committed_checkpoints_no_benchmark_input",
        "detector_resolutions": list(DETECTOR_RESOLUTIONS),
        "primary_detector_resolution": PRIMARY_RESOLUTION,
        "benchmark_role": "external_morphology_comparison_only_not_used_for_detector_construction_or_lane_similarity",
        "moment_control_role": "diagnostic_mu_cov_sufficiency_control_not_physical_distribution_model",
        "cluster_count_expected": len(BASE.CLUSTERS),
        "cluster_count_completed": len(summaries),
        "checkpoint_steps": list(CHECKPOINTS),
        "moment_parity_tolerance": MOMENT_PARITY_TOL,
        "flux_parity_tolerance": FLUX_PARITY_TOL,
        "all_cluster_g3d_unit_speed_pass": bool(all(r["g3d_unit_speed_pass"] for r in summaries)),
        "all_cluster_first_step_exact_pass": bool(all(r["first_step_exact_pass"] for r in summaries)),
        "all_cluster_final_surrogate_moment_parity_pass": bool(all(r["final_surrogate_moment_parity_pass"] for r in summaries)),
        "mean_los_mag_vs_observed_pearson": float(np.mean([r["los_mag_vs_observed_pearson"] for r in summaries])),
        "mean_prior_final_angular_rms_angle_mag_vs_observed_pearson": float(np.mean([r["prior_final_angular_rms_angle_mag_vs_observed_pearson"] for r in summaries])),
        "mean_final_r16_total_variation_distance": _nanmean_key(summaries, "final_r16_total_variation_distance"),
        "mean_final_r16_jensen_shannon_divergence_nats": _nanmean_key(summaries, "final_r16_jensen_shannon_divergence_nats"),
        "mean_final_r16_image_pearson": _nanmean_key(summaries, "final_r16_image_pearson"),
        "mean_final_r32_total_variation_distance": _nanmean_key(summaries, "final_r32_total_variation_distance"),
        "mean_final_r32_jensen_shannon_divergence_nats": _nanmean_key(summaries, "final_r32_jensen_shannon_divergence_nats"),
        "mean_final_r32_image_pearson": _nanmean_key(summaries, "final_r32_image_pearson"),
        "mean_final_r64_total_variation_distance": _nanmean_key(summaries, "final_r64_total_variation_distance"),
        "mean_final_r64_jensen_shannon_divergence_nats": _nanmean_key(summaries, "final_r64_jensen_shannon_divergence_nats"),
        "mean_final_r64_image_pearson": _nanmean_key(summaries, "final_r64_image_pearson"),
        "mean_final_source_max_centroid_vector_error": _nanmean_key(summaries, "final_source_max_centroid_vector_error"),
        "mean_final_source_max_covariance_frobenius_error": _nanmean_key(summaries, "final_source_max_covariance_frobenius_error"),
        "physics_change_authorized": False,
        "candidate_change_authorized": False,
        "production_geometry_change_authorized": False,
        "observer_projection_change_authorized": False,
        "observer_distribution_truncation_authorized": False,
        "moment_control_as_physical_model_authorized": False,
        "native_detector_operation_authorized_for_production": False,
        "observable_selection_authorized": False,
        "native_image_interpretation_required": True,
        "next_experiment_authorized": False,
        "duration_seconds": float(time.perf_counter() - started),
    }

    _write_json(OUT / "validation.json", validation)
    _write_json(
        OUT / "run.json",
        {
            "lab_id": LAB_ID,
            "head_sha": repo["head_sha"],
            "output_directory": str(OUT.relative_to(ROOT)),
            "duration_seconds": validation["duration_seconds"],
        },
    )
    print(json.dumps(validation, indent=2, default=_json_default))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
