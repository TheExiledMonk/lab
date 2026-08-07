#!/usr/bin/env python3
"""PBUF FOUNDATION — STRENGTH FACTORIZATION / PHYSICAL BRIDGE 001.

Purpose
-------
The current weak-lensing foundation uses a dimensionless constant STRENGTH=0.18
inside the initial-state construction:

    u_slow0 = STRENGTH * rho3
    u_fast0 = STRENGTH * rho3 + scaled injection noise

The original scalar Version-A constitutive equation used the same structure:

    u = deformation_strength * normalized(matter)

No physical derivation of 0.18 is present in the current repository.  The latest
independent amplitude audit also showed that the downstream M10 -> G3D -> observer
chain transmits amplitude almost linearly while the independent F160W source has
no absolute mass calibration.

This audit does NOT choose a new strength and does NOT fit anything.  Instead it
factorizes the historical coefficient from the independent source pipeline and
measures exactly how much of the current result is attributable to that scalar.
It constructs two source-amplitude lanes from the SAME independent HST/F160W
source for each cluster:

  LEGACY : s = 0.18
  UNIT   : s = 1.0

The UNIT lane is not a physical prediction.  It is a response-kernel diagnostic:
it asks what the existing PBUF machinery produces per unit dimensionless source
loading.  If the chain were perfectly homogeneous in s, LEGACY would equal 0.18
of UNIT at each amplitude-bearing stage.  Deviations quantify nonlinearities.

No observed kappa pixel values are loaded anywhere in this lab.  Benchmark FITS
header/shape may be used indirectly by the existing common-footprint helper for
WCS geometry only.  There is no assisted lane, no morphology score, and no
benchmark-dependent selection.

The scientific output is a bridge specification:

    physical source -> A_phys(source, medium, geometry) -> unit response kernel

where A_phys is currently MISSING.  The lab records what physical information is
available from the HST product, what is erased by normalization, and what units a
future bridge must supply before 0.18 can be removed scientifically.
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

import pbuf.labs.foundation.independent_source_training_wheels_off001_common_footprint_fix as SRC
import pbuf.labs.foundation.independent_source_training_wheels_off001 as IND
import pbuf.labs.foundation.los_consistent_ray_geometry001 as GEO
import pbuf.labs.foundation.g3d_angular_received_distribution001 as ANG
import pbuf.labs.foundation.g3d_native_angular_detector_image001 as DET
import pbuf.labs.foundation.m10_coverage_25pct_science001 as BASE
from pbuf.models import a8_state as M06_state
from pbuf.core import los_projection as M14

LAB_ID = "PBUF-FOUNDATION-STRENGTH-FACTORIZATION-PHYSICAL-BRIDGE-001"
OUT = ROOT / "runs" / "strength_factorization_physical_bridge001"
DOWNLOADS = OUT / "downloads"
LEGACY_STRENGTH = float(BASE.STRENGTH)
UNIT_STRENGTH = 1.0
EXPECTED_LEGACY_STRENGTH = 0.18
EXPECTED_SUPPORT = GEO.EXPECTED_SUPPORT
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


def _rel_rms(a, b) -> float:
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    m = np.isfinite(a) & np.isfinite(b)
    if not np.any(m): return float("nan")
    num = float(np.sqrt(np.mean((a[m] - b[m]) ** 2)))
    den = max(float(np.sqrt(np.mean(b[m] ** 2))), 1e-300)
    return num / den


def _ratio(a: float, b: float) -> float:
    return float(a / b) if np.isfinite(a) and np.isfinite(b) and abs(b) > 1e-300 else float("nan")


def _effective_exponent(y_legacy: float, y_unit: float) -> float:
    """p where y(s)/y(1) ~= s^p for the single legacy-vs-unit comparison."""
    if not (np.isfinite(y_legacy) and np.isfinite(y_unit)) or y_legacy <= 0 or y_unit <= 0:
        return float("nan")
    return float(math.log(y_legacy / y_unit) / math.log(LEGACY_STRENGTH))


def _initial_state(rho3: np.ndarray, strength: float) -> dict:
    rng = np.random.RandomState(BASE.SEED)
    eq = float(strength) * np.asarray(rho3, dtype=np.float64)
    noise = M06_state.A8_INIT_INJECTION_NOISE * float(strength) * rng.randn(*rho3.shape)
    return {
        "rho_3d": np.asarray(rho3, dtype=np.float64).copy(),
        "u_slow0": eq.copy(),
        "u_fast0": eq + noise,
    }


def _run_chain(rho3: np.ndarray, strength: float) -> dict:
    initial = _initial_state(rho3, strength)
    state = BASE._evolve(initial)
    candidate = BASE._candidate(state)
    vector = BASE._interface_vector(candidate)
    los = M14.project_vector_to_image_plane(*vector, los_axis="z")
    Rx = np.asarray(los["comp_1"], dtype=np.float64)
    Ry = np.asarray(los["comp_2"], dtype=np.float64)
    los_mag = np.hypot(Rx, Ry)

    grid = np.linspace(-BASE.CFG["extent"], BASE.CFG["extent"], Rx.shape[0])
    field = {"xgrid": grid, "ygrid": grid, "rx": Rx, "ry": Ry}
    x0, y0, _, _ = BASE._launch_expanded_25pct()
    groups = GEO._source_groups(x0, y0)
    if len(groups) != EXPECTED_SUPPORT:
        raise RuntimeError(f"expected {EXPECTED_SUPPORT} source bins, got {len(groups)}")

    checkpoints, g3d = GEO._propagate_g3d(
        field, BASE.CFG["step"], BASE.CFG["steps"], x0, y0
    )
    if g3d["max_unit_speed_error"] > GEO.UNIT_SPEED_TOL:
        raise RuntimeError(f"G3D unit-speed gate failed: {g3d['max_unit_speed_error']}")
    first = GEO._first_step_geometry(
        field, x0, y0, checkpoints[1], np.zeros_like(los_mag), los_mag
    )
    if not first["first_step_exact_pass"]:
        raise RuntimeError("first-step exact geometry gate failed")

    final_ang = ANG._angular_distribution_fields(checkpoints[CHECKPOINT], groups)
    gates = ANG._moment_gates(final_ang)
    if gates["second_moment_equals_cov_plus_centroid_outer_relative_rms_error"] > ANG.MOMENT_IDENTITY_TOL:
        raise RuntimeError("angular second-moment identity failed")
    if not gates["covariance_psd_pass"]:
        raise RuntimeError("angular covariance PSD gate failed")
    if not gates["direction_mean_vector_bound_pass"]:
        raise RuntimeError("angular direction-mean bound failed")

    snap = checkpoints[CHECKPOINT]
    vz = np.asarray(snap["vz"], dtype=np.float64)
    if float(np.min(np.abs(vz))) <= DET.VZ_MIN:
        raise RuntimeError("final tangent projection vz too small")
    tx = np.asarray(snap["vx"], dtype=np.float64) / vz
    ty = np.asarray(snap["vy"], dtype=np.float64) / vz

    return {
        "initial": initial,
        "state": state,
        "candidate": candidate,
        "Rx": Rx,
        "Ry": Ry,
        "los_mag": los_mag,
        "checkpoints": checkpoints,
        "g3d": g3d,
        "first": first,
        "final_ang": final_ang,
        "tx": tx,
        "ty": ty,
        "groups": groups,
        "angular_gates": gates,
    }


def _photometry_header(path: Path) -> dict:
    with fits.open(path, memmap=True) as hdul:
        hdr = hdul[0].header
        keys = [
            "BUNIT", "EXPTIME", "PHOTFLAM", "PHOTPLAM", "PHOTBW", "ZPTMAG",
            "ABMAGZP", "VEGAMAG", "FILTER", "INSTRUME", "DETECTOR",
        ]
        return {k: hdr.get(k) for k in keys}


def _run_cluster(cluster: dict) -> tuple[dict, dict]:
    # Independent source construction may use benchmark HEADER/SHAPE for WCS geometry,
    # but no benchmark pixel values are loaded by this lab.
    old_downloads = SRC.DOWNLOADS
    try:
        SRC.DOWNLOADS = DOWNLOADS
        source = SRC._independent_source(cluster)
    finally:
        SRC.DOWNLOADS = old_downloads

    rho2 = np.asarray(source["rho2"], dtype=np.float64)
    rho3 = np.asarray(source["rho3"], dtype=np.float64)
    luminous = np.asarray(source["luminous_common"], dtype=np.float64)
    valid = np.asarray(source["geometry"]["valid_mask"], dtype=bool)

    legacy = _run_chain(rho3, LEGACY_STRENGTH)
    unit = _run_chain(rho3, UNIT_STRENGTH)

    # Exact factorization gate at initialization: both equilibrium and injected noise
    # are explicitly multiplied by the same scalar in the current foundation.
    init_slow_err = _rel_rms(legacy["initial"]["u_slow0"], LEGACY_STRENGTH * unit["initial"]["u_slow0"])
    init_fast_err = _rel_rms(legacy["initial"]["u_fast0"], LEGACY_STRENGTH * unit["initial"]["u_fast0"])

    fields = {
        "m10_los_mag": (legacy["los_mag"], unit["los_mag"]),
        "final_angular_centroid_mag": (
            legacy["final_ang"]["angular_centroid_mag"],
            unit["final_ang"]["angular_centroid_mag"],
        ),
        "final_angular_spread_rms": (
            legacy["final_ang"]["angular_spread_rms"],
            unit["final_ang"]["angular_spread_rms"],
        ),
        "final_angular_rms_angle_mag": (
            legacy["final_ang"]["angular_rms_angle_mag"],
            unit["final_ang"]["angular_rms_angle_mag"],
        ),
    }
    amp = {}
    for name, (a, b) in fields.items():
        ar, br = _rms(a), _rms(b)
        amp[f"legacy_{name}_rms"] = ar
        amp[f"unit_{name}_rms"] = br
        amp[f"legacy_over_unit_{name}_rms_ratio"] = _ratio(ar, br)
        amp[f"{name}_effective_strength_exponent"] = _effective_exponent(ar, br)
        amp[f"legacy_vs_0p18_unit_{name}_relative_rms_error"] = _rel_rms(a, LEGACY_STRENGTH * b)

    ls = luminous[valid & np.isfinite(luminous)]
    positive = ls[ls > 0.0]
    phot = _photometry_header(Path(source["hst_local_path"]))

    summary = {
        "cluster_id": cluster["id"],
        "candidate_id": BASE.CANDIDATE_ID,
        "physical_source_representation": BASE.PHYSICAL_SOURCE,
        "independent_source_role": source["source_role"],
        "benchmark_pixel_values_loaded": False,
        "legacy_strength": LEGACY_STRENGTH,
        "unit_strength": UNIT_STRENGTH,
        "legacy_strength_role": "dimensionless_initial_state_amplitude_coefficient",
        "legacy_strength_physical_derivation_present_in_current_foundation": False,
        "unit_lane_role": "dimensionless_unit_response_kernel_diagnostic_not_physical_prediction",
        "source_absolute_mass_calibration_present": False,
        "source_normalization_role": "positive_F160W_common_footprint_divided_by_its_own_positive_max",
        "source_normalization_erases_absolute_photometric_scale": True,
        "future_bridge_required_form": "u0_physical_field_or_dimensionless_loading_derived_from_physical_mass_energy_density_and_medium_response",
        "future_bridge_must_not_use_observed_kappa": True,
        "future_bridge_must_replace_not_refit_legacy_strength": True,
        "rho2_max": float(np.max(rho2)),
        "rho2_rms": _rms(rho2),
        "rho3_rms": _rms(rho3),
        "luminous_positive_sum": float(np.sum(positive)) if positive.size else 0.0,
        "luminous_positive_mean": float(np.mean(positive)) if positive.size else 0.0,
        "luminous_positive_rms": _rms(positive),
        "luminous_positive_max": float(np.max(positive)) if positive.size else 0.0,
        "normalization_divisor": float(source["alignment"]["positive_luminous_common_max"]),
        "hst_url": source["hst_url"],
        "hst_sha256": source["hst_sha256"],
        "hst_BUNIT": phot.get("BUNIT"),
        "hst_EXPTIME": phot.get("EXPTIME"),
        "hst_PHOTFLAM": phot.get("PHOTFLAM"),
        "hst_PHOTPLAM": phot.get("PHOTPLAM"),
        "hst_PHOTBW": phot.get("PHOTBW"),
        "hst_FILTER": phot.get("FILTER"),
        "hst_INSTRUME": phot.get("INSTRUME"),
        "hst_DETECTOR": phot.get("DETECTOR"),
        "initial_u_slow_exact_0p18_factorization_relative_rms_error": init_slow_err,
        "initial_u_fast_exact_0p18_factorization_relative_rms_error": init_fast_err,
        "legacy_g3d_unit_speed_pass": bool(legacy["g3d"]["max_unit_speed_error"] <= GEO.UNIT_SPEED_TOL),
        "unit_g3d_unit_speed_pass": bool(unit["g3d"]["max_unit_speed_error"] <= GEO.UNIT_SPEED_TOL),
        "legacy_first_step_exact_pass": legacy["first"]["first_step_exact_pass"],
        "unit_first_step_exact_pass": unit["first"]["first_step_exact_pass"],
        "legacy_angular_covariance_psd_pass": legacy["angular_gates"]["covariance_psd_pass"],
        "unit_angular_covariance_psd_pass": unit["angular_gates"]["covariance_psd_pass"],
        "legacy_angular_direction_mean_bound_pass": legacy["angular_gates"]["direction_mean_vector_bound_pass"],
        "unit_angular_direction_mean_bound_pass": unit["angular_gates"]["direction_mean_vector_bound_pass"],
        **amp,
    }

    arrays = {
        "rho2_dimensionless_normalized": rho2,
        "rho3_dimensionless_normalized": rho3,
        "luminous_common_pre_normalization": luminous,
        "valid_mask": valid,
        "legacy_m10_los_mag": legacy["los_mag"],
        "unit_m10_los_mag": unit["los_mag"],
        "legacy_final_angular_centroid_mag": legacy["final_ang"]["angular_centroid_mag"],
        "unit_final_angular_centroid_mag": unit["final_ang"]["angular_centroid_mag"],
        "legacy_final_angular_spread_rms": legacy["final_ang"]["angular_spread_rms"],
        "unit_final_angular_spread_rms": unit["final_ang"]["angular_spread_rms"],
        "legacy_final_angular_rms_angle_mag": legacy["final_ang"]["angular_rms_angle_mag"],
        "unit_final_angular_rms_angle_mag": unit["final_ang"]["angular_rms_angle_mag"],
    }
    return summary, arrays


def _mean(rows: list[dict], key: str) -> float:
    v = np.asarray([r[key] for r in rows], dtype=np.float64)
    return float(np.nanmean(v)) if np.any(np.isfinite(v)) else float("nan")


def main() -> int:
    started = time.perf_counter()
    OUT.mkdir(parents=True, exist_ok=True)
    DOWNLOADS.mkdir(parents=True, exist_ok=True)
    repo = _repo_state()
    _write_json(OUT / "repository_state.json", repo)
    if repo["branch"] != "main" or repo["tracked_changes"] or repo["staged_changes"]:
        v = {"lab_id": LAB_ID, "outcome": "REPOSITORY_GATE_FAILURE", "head_sha": repo["head_sha"]}
        _write_json(OUT / "validation.json", v); print(json.dumps(v, indent=2)); return 2
    if abs(LEGACY_STRENGTH - EXPECTED_LEGACY_STRENGTH) > 1e-15:
        raise RuntimeError(f"legacy STRENGTH changed unexpectedly: {LEGACY_STRENGTH}")

    summaries: list[dict] = []
    failures: list[dict] = []
    for cluster in BASE.CLUSTERS:
        cid = cluster["id"]
        print(f"[{cid}] factorize legacy strength 0.18 against unit response kernel")
        try:
            summary, arrays = _run_cluster(cluster)
            summaries.append(summary)
            cdir = OUT / "clusters" / cid
            cdir.mkdir(parents=True, exist_ok=True)
            _write_json(cdir / "strength_factorization_summary.json", summary)
            np.savez_compressed(cdir / "strength_factorization_fields.npz", **arrays)
        except Exception as exc:
            failures.append({"cluster_id": cid, "error": repr(exc)})
            _write_json(OUT / "cluster_failures.json", failures)
            raise

    _write_csv(OUT / "strength_factorization_summary.csv", summaries)
    _write_json(OUT / "cluster_failures.json", failures)

    validation = {
        "lab_id": LAB_ID,
        "outcome": "Outcome A — STRENGTH FACTORIZATION / PHYSICAL BRIDGE AUDIT COMPLETE",
        "head_sha": repo["head_sha"],
        "candidate_id": BASE.CANDIDATE_ID,
        "physical_source_representation": BASE.PHYSICAL_SOURCE,
        "cluster_count_expected": len(BASE.CLUSTERS),
        "cluster_count_completed": len(summaries),
        "benchmark_pixel_values_loaded": False,
        "legacy_strength": LEGACY_STRENGTH,
        "legacy_strength_role": "dimensionless_initial_state_amplitude_coefficient",
        "legacy_strength_physical_derivation_present_in_current_foundation": False,
        "source_absolute_mass_calibration_present": False,
        "all_source_normalization_erases_absolute_photometric_scale": bool(all(r["source_normalization_erases_absolute_photometric_scale"] for r in summaries)),
        "unit_response_kernel_is_physical_prediction": False,
        "physical_bridge_present": False,
        "physical_bridge_required": True,
        "physical_bridge_target": "replace_legacy_strength_with_source_and_medium_derived_dimensionless_loading_or_initial_deformation_field",
        "mean_initial_u_slow_exact_0p18_factorization_relative_rms_error": _mean(summaries, "initial_u_slow_exact_0p18_factorization_relative_rms_error"),
        "mean_initial_u_fast_exact_0p18_factorization_relative_rms_error": _mean(summaries, "initial_u_fast_exact_0p18_factorization_relative_rms_error"),
        "mean_legacy_over_unit_m10_los_mag_rms_ratio": _mean(summaries, "legacy_over_unit_m10_los_mag_rms_ratio"),
        "mean_m10_los_mag_effective_strength_exponent": _mean(summaries, "m10_los_mag_effective_strength_exponent"),
        "mean_legacy_over_unit_final_angular_rms_angle_mag_rms_ratio": _mean(summaries, "legacy_over_unit_final_angular_rms_angle_mag_rms_ratio"),
        "mean_final_angular_rms_angle_mag_effective_strength_exponent": _mean(summaries, "final_angular_rms_angle_mag_effective_strength_exponent"),
        "all_legacy_g3d_unit_speed_pass": bool(all(r["legacy_g3d_unit_speed_pass"] for r in summaries)),
        "all_unit_g3d_unit_speed_pass": bool(all(r["unit_g3d_unit_speed_pass"] for r in summaries)),
        "all_legacy_first_step_exact_pass": bool(all(r["legacy_first_step_exact_pass"] for r in summaries)),
        "all_unit_first_step_exact_pass": bool(all(r["unit_first_step_exact_pass"] for r in summaries)),
        "legacy_strength_removal_authorized": False,
        "physical_mass_conversion_authorized": False,
        "physics_change_authorized": False,
        "candidate_change_authorized": False,
        "trajectory_change_authorized": False,
        "benchmark_fit_authorized": False,
        "next_experiment_authorized": False,
        "interpretation_required": True,
        "duration_seconds": float(time.perf_counter() - started),
    }
    _write_json(OUT / "validation.json", validation)
    _write_json(OUT / "run.json", {"validation": validation, "summaries": summaries})
    print(json.dumps(validation, indent=2, default=_json_default))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
