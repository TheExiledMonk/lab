#!/usr/bin/env python3
"""PBUF FOUNDATION — INDEPENDENT MASS/RESPONSE AMPLITUDE AUDIT 001.

Purpose
-------
The first training-wheels-off run showed that independent HST/F160W morphology
survives the PBUF chain but much more weakly than the historical kappa-derived
source. This audit asks a narrower physics question before changing any law:

    Is absolute lens/source amplitude being erased or suppressed before/during
    the G3D trajectory calculation?

The audit uses ONLY the independent HST/F160W common-footprint source until every
amplitude experiment has completed. Observed kappa pixel values are loaded only
at the end for measurement. They never select an amplitude, lane, coefficient,
or verdict.

Two diagnostic amplitude families are run with the same morphology:

  SOURCE family
      rho3 -> a * rho3 -> existing PBUF evolution -> M10 -> G3D observer

  RESPONSE family
      fixed baseline rho3 -> fixed baseline PBUF/M10 field -> a * (Rx,Ry)
      -> existing G3D observer

with a in {0.25, 0.5, 1, 2, 4}.

These multipliers are NOT masses and NOT fitted parameters. They are dimensionless
sensitivity probes. Their purpose is to localize amplitude loss:

- if source amplitude barely changes M10, the source->response bridge is suppressing
  absolute scale;
- if M10 scales but observer amplitude barely changes, trajectory normalization is
  suppressing it;
- if both carry scale approximately linearly, the missing piece is likely the absent
  physical mass/energy calibration rather than numerical amplitude loss.

No conventional gravitational law is introduced. No assisted/kappa-fed source lane
is run. No amplitude is selected by benchmark correlation.
"""
from __future__ import annotations

import csv
import json
import math
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
from astropy.io import fits

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))

import pbuf.labs.foundation.independent_source_training_wheels_off001 as LAB
import pbuf.labs.foundation.independent_source_training_wheels_off001_common_footprint_fix as FIX
import pbuf.labs.foundation.m10_coverage_25pct_science001 as BASE
import pbuf.labs.foundation.los_consistent_ray_geometry001 as GEO
import pbuf.labs.foundation.g3d_angular_received_distribution001 as ANG
import pbuf.labs.foundation.g3d_native_angular_detector_image001 as DET
from pbuf.core import los_projection as M14

LAB_ID = "PBUF-FOUNDATION-INDEPENDENT-MASS-RESPONSE-AMPLITUDE-AUDIT-001"
OUT = ROOT / "runs" / "independent_mass_response_amplitude_audit001"
DOWNLOADS = OUT / "downloads"
SCALES = (0.25, 0.5, 1.0, 2.0, 4.0)
BASELINE_SCALE = 1.0
EXPECTED_SUPPORT = GEO.EXPECTED_SUPPORT
UNIT_SPEED_TOL = GEO.UNIT_SPEED_TOL
FIRST_STEP_TOL = GEO.FIRST_STEP_TOL
CHECKPOINT = GEO.CHECKPOINTS[-1]


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
        w.writeheader(); w.writerows(rows)


def _rms(a) -> float:
    x = np.asarray(a, dtype=np.float64)
    m = np.isfinite(x)
    return float(np.sqrt(np.mean(x[m] * x[m]))) if np.any(m) else float("nan")


def _corr(a, b):
    return LAB._corr(np.asarray(a, dtype=np.float64), np.asarray(b, dtype=np.float64))


def _field_rms_nonzero(a) -> float:
    x = np.asarray(a, dtype=np.float64)
    x = x[np.isfinite(x)]
    return float(np.sqrt(np.mean(x*x))) if x.size else float("nan")


def _elasticity(scales, values) -> dict:
    """Log-log slope of positive finite output amplitude versus probe scale."""
    s = np.asarray(scales, dtype=np.float64)
    v = np.asarray(values, dtype=np.float64)
    m = np.isfinite(s) & np.isfinite(v) & (s > 0) & (v > 0)
    if int(m.sum()) < 2:
        return {"loglog_slope": float("nan"), "loglog_r2": float("nan"), "count": int(m.sum())}
    x = np.log(s[m]); y = np.log(v[m])
    A = np.column_stack((x, np.ones_like(x)))
    beta, *_ = np.linalg.lstsq(A, y, rcond=None)
    pred = A @ beta
    ss_res = float(np.sum((y-pred)**2))
    ss_tot = float(np.sum((y-y.mean())**2))
    r2 = 1.0 - ss_res/ss_tot if ss_tot > 0 else float("nan")
    return {"loglog_slope": float(beta[0]), "loglog_r2": r2, "count": int(m.sum())}


def _trajectory_from_los(Rx: np.ndarray, Ry: np.ndarray) -> dict:
    Rx = np.asarray(Rx, dtype=np.float64)
    Ry = np.asarray(Ry, dtype=np.float64)
    los_mag = np.hypot(Rx, Ry)
    grid = np.linspace(-BASE.CFG["extent"], BASE.CFG["extent"], Rx.shape[0])
    field = {"xgrid": grid, "ygrid": grid, "rx": Rx, "ry": Ry}
    x0, y0, _, _ = BASE._launch_expanded_25pct()
    groups = GEO._source_groups(x0, y0)
    if len(groups) != EXPECTED_SUPPORT:
        raise RuntimeError(f"expected {EXPECTED_SUPPORT} source bins, got {len(groups)}")

    checkpoints, g3d = GEO._propagate_g3d(field, BASE.CFG["step"], BASE.CFG["steps"], x0, y0)
    if g3d["max_unit_speed_error"] > UNIT_SPEED_TOL:
        raise RuntimeError(f"G3D unit-speed gate failed: {g3d['max_unit_speed_error']}")
    zero_obs = np.zeros_like(los_mag)
    first = GEO._first_step_geometry(field, x0, y0, checkpoints[1], zero_obs, los_mag)
    if not first["first_step_exact_pass"]:
        raise RuntimeError("first-step exact geometry gate failed")

    final_ang = ANG._angular_distribution_fields(checkpoints[CHECKPOINT], groups)
    gates = ANG._moment_gates(final_ang)
    if gates["second_moment_equals_cov_plus_centroid_outer_relative_rms_error"] > ANG.MOMENT_IDENTITY_TOL:
        raise RuntimeError("angular second-moment identity failed")
    if not gates["covariance_psd_pass"]:
        raise RuntimeError("angular covariance PSD failed")
    if not gates["direction_mean_vector_bound_pass"]:
        raise RuntimeError("angular mean-direction bound failed")

    snap = checkpoints[CHECKPOINT]
    vz = np.asarray(snap["vz"], dtype=np.float64)
    if float(np.min(np.abs(vz))) <= DET.VZ_MIN:
        raise RuntimeError("final tangent projection vz too small")
    tx = np.asarray(snap["vx"], dtype=np.float64) / vz
    ty = np.asarray(snap["vy"], dtype=np.float64) / vz
    trans_v = np.hypot(np.asarray(snap["vx"]), np.asarray(snap["vy"]))
    trans_disp = np.hypot(np.asarray(snap["x"])-x0, np.asarray(snap["y"])-y0)

    return {
        "los_mag": los_mag,
        "final_ang": final_ang,
        "tx": tx, "ty": ty,
        "final_transverse_velocity_rms": _rms(trans_v),
        "final_transverse_displacement_rms": _rms(trans_disp),
        "g3d_unit_speed_max_error": float(g3d["max_unit_speed_error"]),
        "first_step_exact_max_vector_error": first["first_step_exact_max_vector_error"],
        "angular_identity_error": gates["second_moment_equals_cov_plus_centroid_outer_relative_rms_error"],
        "angular_covariance_psd_pass": gates["covariance_psd_pass"],
        "angular_direction_mean_bound_pass": gates["direction_mean_vector_bound_pass"],
    }


def _source_scaled_chain(rho3_base: np.ndarray, scale: float) -> dict:
    # Scale before the existing initial-state strength factor. No parameter inside
    # production code is changed.
    rho3 = np.asarray(rho3_base, dtype=np.float64) * float(scale)
    chain = LAB._run_chain_from_rho3(rho3, observed_for_first_step=None)
    return {
        "los_mag": chain["los_mag"],
        "final_ang": chain["final_ang"],
        "tx": chain["tx"], "ty": chain["ty"],
        "final_transverse_velocity_rms": _rms(np.hypot(
            chain["checkpoints"][CHECKPOINT]["vx"], chain["checkpoints"][CHECKPOINT]["vy"])),
        "final_transverse_displacement_rms": _rms(np.hypot(
            chain["checkpoints"][CHECKPOINT]["x"] - BASE._launch_expanded_25pct()[0],
            chain["checkpoints"][CHECKPOINT]["y"] - BASE._launch_expanded_25pct()[1])),
        "g3d_unit_speed_max_error": float(chain["g3d"]["max_unit_speed_error"]),
        "first_step_exact_max_vector_error": chain["first"]["first_step_exact_max_vector_error"],
        "angular_identity_error": chain["angular_gates"]["second_moment_equals_cov_plus_centroid_outer_relative_rms_error"],
        "angular_covariance_psd_pass": chain["angular_gates"]["covariance_psd_pass"],
        "angular_direction_mean_bound_pass": chain["angular_gates"]["direction_mean_vector_bound_pass"],
        "candidate": chain["candidate"],
    }


def _independent_source(cluster: dict) -> dict:
    # Reuse the corrected metadata-only common-footprint source construction, but
    # direct downloads into this lab's own output tree.
    old = FIX.DOWNLOADS
    FIX.DOWNLOADS = DOWNLOADS
    try:
        src = FIX._independent_source(cluster)
    finally:
        FIX.DOWNLOADS = old
    return src


def _amplitude_metrics(chain: dict) -> dict:
    ang = chain["final_ang"]
    return {
        "m10_los_mag_rms": _field_rms_nonzero(chain["los_mag"]),
        "final_angular_centroid_mag_rms": _field_rms_nonzero(ang["angular_centroid_mag"]),
        "final_angular_spread_rms_rms": _field_rms_nonzero(ang["angular_spread_rms"]),
        "final_angular_rms_angle_mag_rms": _field_rms_nonzero(ang["angular_rms_angle_mag"]),
        "final_transverse_velocity_rms": chain["final_transverse_velocity_rms"],
        "final_transverse_displacement_rms": chain["final_transverse_displacement_rms"],
    }


def _row(cid: str, family: str, scale: float, chain: dict) -> dict:
    return {
        "cluster_id": cid,
        "family": family,
        "scale": float(scale),
        **_amplitude_metrics(chain),
        "g3d_unit_speed_max_error": chain["g3d_unit_speed_max_error"],
        "g3d_unit_speed_pass": bool(chain["g3d_unit_speed_max_error"] <= UNIT_SPEED_TOL),
        "first_step_exact_max_vector_error": chain["first_step_exact_max_vector_error"],
        "first_step_exact_pass": bool(chain["first_step_exact_max_vector_error"] <= FIRST_STEP_TOL),
        "angular_identity_error": chain["angular_identity_error"],
        "angular_covariance_psd_pass": bool(chain["angular_covariance_psd_pass"]),
        "angular_direction_mean_bound_pass": bool(chain["angular_direction_mean_bound_pass"]),
    }


def _run_cluster_independent_only(cluster: dict) -> dict:
    cid = cluster["id"]
    src = _independent_source(cluster)
    rho2 = np.asarray(src["rho2"], dtype=np.float64)
    rho3 = np.asarray(src["rho3"], dtype=np.float64)
    luminous = np.asarray(src["luminous_common"], dtype=np.float64)
    valid = np.asarray(src["geometry"]["valid_mask"], dtype=bool)

    # Explicitly expose how much absolute photometric amplitude the current source
    # normalization removes before PBUF sees rho2.
    lum = luminous[valid & np.isfinite(luminous)]
    positive = lum[lum > 0]
    source_diag = {
        "luminous_positive_sum": float(np.sum(positive)) if positive.size else 0.0,
        "luminous_positive_mean": float(np.mean(positive)) if positive.size else 0.0,
        "luminous_positive_rms": _rms(positive),
        "luminous_positive_max": float(np.max(positive)) if positive.size else 0.0,
        "normalized_rho2_max": float(np.nanmax(rho2)),
        "normalized_rho2_rms": _rms(rho2),
        "normalization_divisor": float(np.max(positive)) if positive.size else float("nan"),
        "normalization_erases_absolute_photometric_scale": True,
        "existing_initial_state_strength": float(BASE.STRENGTH),
    }

    source_chains: dict[float, dict] = {}
    rows: list[dict] = []
    for scale in SCALES:
        print(f"[{cid}] SOURCE amplitude probe a={scale:g}")
        ch = _source_scaled_chain(rho3, scale)
        source_chains[scale] = ch
        rows.append(_row(cid, "SOURCE", scale, ch))

    baseline = source_chains[BASELINE_SCALE]
    vector = BASE._interface_vector(baseline["candidate"])
    los = M14.project_vector_to_image_plane(*vector, los_axis="z")
    Rx0 = np.asarray(los["comp_1"], dtype=np.float64)
    Ry0 = np.asarray(los["comp_2"], dtype=np.float64)

    response_chains: dict[float, dict] = {}
    for scale in SCALES:
        print(f"[{cid}] RESPONSE amplitude probe a={scale:g}")
        ch = _trajectory_from_los(Rx0 * scale, Ry0 * scale)
        response_chains[scale] = ch
        rows.append(_row(cid, "RESPONSE", scale, ch))

    return {
        "cluster": cluster,
        "source": src,
        "source_diag": source_diag,
        "source_chains": source_chains,
        "response_chains": response_chains,
        "rows": rows,
    }


def _add_benchmark_metrics(row: dict, chain: dict, observed: np.ndarray) -> None:
    fields = {
        "m10_los_mag": chain["los_mag"],
        "final_angular_centroid_mag": chain["final_ang"]["angular_centroid_mag"],
        "final_angular_spread_rms": chain["final_ang"]["angular_spread_rms"],
        "final_angular_rms_angle_mag": chain["final_ang"]["angular_rms_angle_mag"],
    }
    for name, arr in fields.items():
        p, s, n = _corr(arr, observed)
        row[f"{name}_vs_observed_pearson"] = p
        row[f"{name}_vs_observed_spearman"] = s
        row[f"{name}_vs_observed_count"] = n


def _family_elasticities(rows: list[dict], family: str) -> dict:
    rr = sorted([r for r in rows if r["family"] == family], key=lambda r: r["scale"])
    out = {}
    for key in (
        "m10_los_mag_rms",
        "final_angular_centroid_mag_rms",
        "final_angular_spread_rms_rms",
        "final_angular_rms_angle_mag_rms",
        "final_transverse_velocity_rms",
        "final_transverse_displacement_rms",
    ):
        e = _elasticity([r["scale"] for r in rr], [r[key] for r in rr])
        out[f"{family.lower()}_{key}_elasticity_loglog_slope"] = e["loglog_slope"]
        out[f"{family.lower()}_{key}_elasticity_loglog_r2"] = e["loglog_r2"]
    return out


def _run_all_without_benchmark() -> list[dict]:
    # This function deliberately completes every source/response amplitude probe for
    # every cluster before any kappa pixel value is loaded anywhere in main().
    results = []
    for cluster in BASE.CLUSTERS:
        print(f"[{cluster['id']}] INDEPENDENT mass/response amplitude audit — benchmark withheld")
        results.append(_run_cluster_independent_only(cluster))
    return results


def main() -> int:
    started = time.perf_counter()
    OUT.mkdir(parents=True, exist_ok=True); DOWNLOADS.mkdir(parents=True, exist_ok=True)
    repo = _repo_state(); _write_json(OUT / "repository_state.json", repo)
    if repo["branch"] != "main" or repo["tracked_changes"] or repo["staged_changes"]:
        v = {"lab_id": LAB_ID, "outcome": "REPOSITORY_GATE_FAILURE", "head_sha": repo["head_sha"]}
        _write_json(OUT / "validation.json", v); print(json.dumps(v, indent=2)); return 2

    # ================================================================
    # STRICT INDEPENDENT PHASE — observed kappa pixels are not loaded.
    # ================================================================
    results = _run_all_without_benchmark()
    benchmark_loaded_after_all_independent_amplitude_runs_complete = True

    summaries = []
    all_rows = []
    for result in results:
        cluster = result["cluster"]
        cid = cluster["id"]

        # ONLY NOW load benchmark pixels, after every independent amplitude run for
        # every cluster has completed. Metadata/WCS had been allowed upstream.
        kpath = LAB._kappa_path(cluster)
        with fits.open(kpath) as hdul:
            kappa_native = np.asarray(hdul[0].data, dtype=np.float64)
        observed, _unused_assisted = FIX._benchmark_on_common_grid(
            kappa_native, result["source"]["geometry"]
        )
        del _unused_assisted

        # Attach benchmark measurements only; never select a scale.
        for row in result["rows"]:
            family = row["family"]; scale = float(row["scale"])
            chain = result["source_chains"][scale] if family == "SOURCE" else result["response_chains"][scale]
            _add_benchmark_metrics(row, chain, observed)
            all_rows.append(row)

        source_base = result["source_chains"][BASELINE_SCALE]
        psrc, ssrc, nsrc = _corr(result["source"]["rho2"], observed)
        source_el = _family_elasticities(result["rows"], "SOURCE")
        response_el = _family_elasticities(result["rows"], "RESPONSE")
        base_row = next(r for r in result["rows"] if r["family"] == "SOURCE" and r["scale"] == BASELINE_SCALE)

        summary = {
            "cluster_id": cid,
            "candidate_id": BASE.CANDIDATE_ID,
            "physical_source_representation": BASE.PHYSICAL_SOURCE,
            "independent_source_role": "HST_F160W_positive_luminous_structure_proxy_common_footprint_no_kappa_pixel_values",
            "independent_source_limit": "luminous_structure_proxy_not_mass_map_no_absolute_mass_calibration",
            "amplitude_probe_role": "dimensionless_sensitivity_only_not_mass_not_fit",
            "source_amplitude_family": "scale_rho3_before_existing_initial_state_strength",
            "response_amplitude_family": "scale_fixed_baseline_M10_LOS_Rx_Ry_before_G3D",
            "probe_scales": list(SCALES),
            "benchmark_role": "measurement_only_loaded_after_all_independent_amplitude_runs_complete",
            "benchmark_values_loaded_after_all_independent_amplitude_runs_complete": benchmark_loaded_after_all_independent_amplitude_runs_complete,
            "observed_kappa_used_to_select_amplitude": False,
            "observed_kappa_used_in_independent_source": False,
            "source_rho2_vs_observed_pearson": psrc,
            "source_rho2_vs_observed_spearman": ssrc,
            "source_rho2_vs_observed_count": nsrc,
            **result["source_diag"],
            **source_el,
            **response_el,
            "baseline_m10_los_mag_vs_observed_pearson": base_row["m10_los_mag_vs_observed_pearson"],
            "baseline_final_angular_rms_angle_mag_vs_observed_pearson": base_row["final_angular_rms_angle_mag_vs_observed_pearson"],
            "all_g3d_unit_speed_pass": bool(all(r["g3d_unit_speed_pass"] for r in result["rows"])),
            "all_first_step_exact_pass": bool(all(r["first_step_exact_pass"] for r in result["rows"])),
            "all_angular_covariance_psd_pass": bool(all(r["angular_covariance_psd_pass"] for r in result["rows"])),
            "all_angular_direction_mean_bound_pass": bool(all(r["angular_direction_mean_bound_pass"] for r in result["rows"])),
        }
        summaries.append(summary)

        cdir = OUT / "clusters" / cid; cdir.mkdir(parents=True, exist_ok=True)
        _write_json(cdir / "mass_response_amplitude_summary.json", summary)
        _write_csv(cdir / "mass_response_amplitude_scales.csv", result["rows"])
        np.savez_compressed(
            cdir / "mass_response_amplitude_fields.npz",
            independent_source_rho2=result["source"]["rho2"],
            observed_kappa_reference_only=observed,
            **{
                f"source_a{str(a).replace('.', 'p')}__m10_los_mag": result["source_chains"][a]["los_mag"]
                for a in SCALES
            },
            **{
                f"source_a{str(a).replace('.', 'p')}__angular_rms": result["source_chains"][a]["final_ang"]["angular_rms_angle_mag"]
                for a in SCALES
            },
            **{
                f"response_a{str(a).replace('.', 'p')}__angular_rms": result["response_chains"][a]["final_ang"]["angular_rms_angle_mag"]
                for a in SCALES
            },
        )

    _write_csv(OUT / "mass_response_amplitude_summary.csv", summaries)
    _write_csv(OUT / "mass_response_amplitude_scales.csv", all_rows)

    def mean_key(k):
        vals = np.asarray([s[k] for s in summaries], dtype=np.float64)
        return float(np.nanmean(vals)) if np.any(np.isfinite(vals)) else float("nan")

    validation = {
        "lab_id": LAB_ID,
        "outcome": "Outcome A — INDEPENDENT MASS/RESPONSE AMPLITUDE AUDIT COMPLETE",
        "head_sha": repo["head_sha"],
        "candidate_id": BASE.CANDIDATE_ID,
        "physical_source_representation": BASE.PHYSICAL_SOURCE,
        "independent_source_role": "HST_F160W_positive_luminous_structure_proxy_common_footprint_no_kappa_pixel_values",
        "independent_source_limit": "luminous_structure_proxy_not_mass_map_no_absolute_mass_calibration",
        "existing_initial_state_strength": float(BASE.STRENGTH),
        "probe_scales": list(SCALES),
        "probe_scales_are_physical_masses": False,
        "probe_scales_are_fit_parameters": False,
        "benchmark_role": "measurement_only_after_all_independent_amplitude_runs_complete",
        "benchmark_values_loaded_after_all_independent_amplitude_runs_complete": True,
        "observed_kappa_used_to_select_amplitude": False,
        "cluster_count_expected": len(BASE.CLUSTERS),
        "cluster_count_completed": len(summaries),
        "all_source_normalization_erases_absolute_photometric_scale": bool(all(s["normalization_erases_absolute_photometric_scale"] for s in summaries)),
        "all_g3d_unit_speed_pass": bool(all(s["all_g3d_unit_speed_pass"] for s in summaries)),
        "all_first_step_exact_pass": bool(all(s["all_first_step_exact_pass"] for s in summaries)),
        "mean_source_rho2_vs_observed_pearson": mean_key("source_rho2_vs_observed_pearson"),
        "mean_baseline_m10_los_mag_vs_observed_pearson": mean_key("baseline_m10_los_mag_vs_observed_pearson"),
        "mean_baseline_final_angular_rms_angle_mag_vs_observed_pearson": mean_key("baseline_final_angular_rms_angle_mag_vs_observed_pearson"),
        "mean_source_m10_los_mag_rms_elasticity_loglog_slope": mean_key("source_m10_los_mag_rms_elasticity_loglog_slope"),
        "mean_source_final_angular_rms_angle_mag_rms_elasticity_loglog_slope": mean_key("source_final_angular_rms_angle_mag_rms_elasticity_loglog_slope"),
        "mean_response_final_angular_rms_angle_mag_rms_elasticity_loglog_slope": mean_key("response_final_angular_rms_angle_mag_rms_elasticity_loglog_slope"),
        "mean_response_final_transverse_velocity_rms_elasticity_loglog_slope": mean_key("response_final_transverse_velocity_rms_elasticity_loglog_slope"),
        "mean_response_final_transverse_displacement_rms_elasticity_loglog_slope": mean_key("response_final_transverse_displacement_rms_elasticity_loglog_slope"),
        "physical_mass_bridge_present": False,
        "absolute_mass_calibration_authorized": False,
        "physics_change_authorized": False,
        "candidate_change_authorized": False,
        "trajectory_change_authorized": False,
        "amplitude_selection_authorized": False,
        "mass_response_interpretation_required": True,
        "next_experiment_authorized": False,
        "duration_seconds": float(time.perf_counter() - started),
    }
    _write_json(OUT / "validation.json", validation)
    _write_json(OUT / "run.json", {"validation": validation, "summaries": summaries})
    print(json.dumps(validation, indent=2, default=_json_default))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
