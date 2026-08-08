#!/usr/bin/env python3
"""PBUF FOUNDATION — NATIVE RESPONSE FINGERPRINT MATCH 001.

Compare existing native PBUF fields against the frozen effective weak-field
response fingerprint derived in known_source_inverse_response_fingerprint001.

This lab does not rescale native amplitudes to the macroscopic benchmark and
does not claim a microscopic substrate interpretation. It asks only whether the
existing native variables reproduce the required *scaling structure*:

  local loading       ~ rho
  surface response    ~ M/R
  far response        ~ M/r
  weak-field response ~ additive

The comparison is deliberately inverse/structural. A match identifies a viable
native response role; it does not establish SI amplitude, derive G, or identify
what spacetime is made of.
"""
from __future__ import annotations

import json
import math
import subprocess
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]

import pbuf.labs.foundation.native_field_curvature_dimension_audit001 as NATIVE
import pbuf.labs.foundation.m10_coverage_25pct_science001 as BASE

LAB_ID = "PBUF-FOUNDATION-NATIVE-RESPONSE-FINGERPRINT-MATCH-001"

# Frozen effective fingerprint. These are comparison targets, not fitted values.
TARGET_LOCAL_DENSITY_EXPONENT = 1.0
TARGET_SURFACE_MASS_EXPONENT = 1.0
TARGET_SURFACE_RADIUS_FIXED_M_EXPONENT = -1.0
TARGET_FAR_MASS_EXPONENT = 1.0
TARGET_FAR_RADIUS_EXPONENT = -1.0
TARGET_ADDITIVITY = 0.0  # relative residual

EXPONENT_WINDOW = 0.20
ADDITIVITY_REL_TOL = 1.0e-10
ALG_TOL = 1.0e-12

# Native loading amplitudes are dimensionless implementation units. They are
# never mapped to kg or to the empirical G benchmark in this lab.
FIXED_RADIUS = 4.5
DENSITY_LADDER = (0.0125, 0.025, 0.05, 0.10, 0.20)
MASS_LADDER = (2.0, 4.0, 8.0, 16.0, 32.0)
RADIUS_LADDER = (2.5, 3.5, 4.5, 5.5, 6.5)
FIXED_NATIVE_MASS = 10.0
FAR_SOURCE_RADIUS = 3.5
FAR_SOURCE_DENSITY = 0.05
FAR_PROBE_RADII = (7.0, 8.0, 10.0, 12.0, 14.0)


def git(*args: str) -> str:
    return subprocess.check_output(["git", *args], cwd=ROOT, text=True).strip()


def repo_state() -> dict:
    return {
        "branch": git("rev-parse", "--abbrev-ref", "HEAD"),
        "head_sha": git("rev-parse", "HEAD"),
        "tracked_changes": git("diff", "--name-only"),
        "staged_changes": git("diff", "--name-only", "--cached"),
    }


def loglog_fit(xs, ys) -> dict:
    x = np.asarray(xs, dtype=np.float64)
    y = np.asarray(ys, dtype=np.float64)
    m = np.isfinite(x) & np.isfinite(y) & (x > 0.0) & (y > 0.0)
    if int(np.count_nonzero(m)) < 3:
        return {"slope": float("nan"), "r2": float("nan"), "count": int(np.count_nonzero(m))}
    lx = np.log(x[m])
    ly = np.log(y[m])
    A = np.column_stack((lx, np.ones_like(lx)))
    beta, *_ = np.linalg.lstsq(A, ly, rcond=None)
    pred = A @ beta
    ss_res = float(np.sum((ly - pred) ** 2))
    ss_tot = float(np.sum((ly - np.mean(ly)) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0.0 else float("nan")
    return {"slope": float(beta[0]), "r2": r2, "count": int(np.count_nonzero(m))}


def exponent_match(value: float, target: float) -> bool:
    return bool(math.isfinite(value) and abs(value - target) <= EXPONENT_WINDOW)


def native_fields(rho: np.ndarray) -> dict:
    state = NATIVE._noise_free_state(np.asarray(rho, dtype=np.float64))
    candidate = BASE._candidate(state)
    vx, vy, vz = BASE._interface_vector(candidate)
    vx = np.asarray(vx, dtype=np.float64)
    vy = np.asarray(vy, dtype=np.float64)
    vz = np.asarray(vz, dtype=np.float64)
    c = np.asarray(state["c_state"], dtype=np.float64)
    vmag = np.sqrt(vx * vx + vy * vy + vz * vz)
    return {
        "c_state": c,
        "vx": vx,
        "vy": vy,
        "vz": vz,
        "m10_mag": vmag,
        "all_finite": bool(
            np.all(np.isfinite(c)) and np.all(np.isfinite(vx)) and
            np.all(np.isfinite(vy)) and np.all(np.isfinite(vz))
        ),
        "max_abs_state": float(state["max_abs_state"]),
    }


def shell_mean(field: np.ndarray, source: dict) -> float:
    shell = NATIVE._shell_mask(source["effective_radius"])
    return float(np.mean(np.abs(np.asarray(field, dtype=np.float64)[shell])))


def center_abs(field: np.ndarray) -> float:
    return float(abs(np.asarray(field, dtype=np.float64)[NATIVE.CENTER_INDEX]))


def sphere_with_native_mass(radius: float, native_mass: float) -> dict:
    # Geometry is fixed first. Density is then chosen only to hold the specified
    # integrated native source constant across the radius ladder. This is an
    # orthogonal scaling experiment, not a fit to a response target.
    unit = NATIVE._sphere(radius, 1.0)
    density = float(native_mass) / float(unit["volume"])
    return NATIVE._sphere(radius, density)


def density_scaling_test() -> dict:
    rows = []
    for d in DENSITY_LADDER:
        src = NATIVE._sphere(FIXED_RADIUS, d)
        f = native_fields(src["rho"])
        rows.append({
            "density": d,
            "c_center": center_abs(f["c_state"]),
            "c_shell": shell_mean(f["c_state"], src),
            "m10_shell": shell_mean(f["m10_mag"], src),
            "finite": f["all_finite"],
            "max_abs_state": f["max_abs_state"],
        })
    return {
        "rows": rows,
        "c_center": loglog_fit([r["density"] for r in rows], [r["c_center"] for r in rows]),
        "c_shell": loglog_fit([r["density"] for r in rows], [r["c_shell"] for r in rows]),
        "m10_shell": loglog_fit([r["density"] for r in rows], [r["m10_shell"] for r in rows]),
    }


def mass_scaling_fixed_radius_test() -> dict:
    # Vary integrated native source at fixed radius. Density is implied by M/V.
    unit = NATIVE._sphere(FIXED_RADIUS, 1.0)
    rows = []
    for m in MASS_LADDER:
        d = m / unit["volume"]
        src = NATIVE._sphere(FIXED_RADIUS, d)
        f = native_fields(src["rho"])
        rows.append({
            "native_mass": src["integrated_source"],
            "density": d,
            "c_shell": shell_mean(f["c_state"], src),
            "m10_shell": shell_mean(f["m10_mag"], src),
            "finite": f["all_finite"],
        })
    return {
        "rows": rows,
        "c_shell": loglog_fit([r["native_mass"] for r in rows], [r["c_shell"] for r in rows]),
        "m10_shell": loglog_fit([r["native_mass"] for r in rows], [r["m10_shell"] for r in rows]),
    }


def radius_scaling_fixed_mass_test() -> dict:
    rows = []
    for radius in RADIUS_LADDER:
        src = sphere_with_native_mass(radius, FIXED_NATIVE_MASS)
        f = native_fields(src["rho"])
        rows.append({
            "R_eff": src["effective_radius"],
            "native_mass": src["integrated_source"],
            "density": src["density"],
            "c_shell": shell_mean(f["c_state"], src),
            "m10_shell": shell_mean(f["m10_mag"], src),
            "finite": f["all_finite"],
        })
    return {
        "rows": rows,
        "mass_constancy_rel_span": (
            max(r["native_mass"] for r in rows) - min(r["native_mass"] for r in rows)
        ) / FIXED_NATIVE_MASS,
        "c_shell": loglog_fit([r["R_eff"] for r in rows], [r["c_shell"] for r in rows]),
        "m10_shell": loglog_fit([r["R_eff"] for r in rows], [r["m10_shell"] for r in rows]),
    }


def radial_probe(field: np.ndarray, radius: float) -> float:
    # Probe along +x from the source center. FAR_PROBE_RADII are integral grid
    # positions so this introduces no interpolation or scale fitting.
    zi, yi, xi = NATIVE.CENTER_INDEX
    offset = int(round(radius / NATIVE.DX))
    if abs(offset * NATIVE.DX - radius) > ALG_TOL:
        raise RuntimeError("far probe radius is not an exact grid location")
    return float(abs(np.asarray(field, dtype=np.float64)[zi, yi, xi + offset]))


def far_radius_test() -> dict:
    src = NATIVE._sphere(FAR_SOURCE_RADIUS, FAR_SOURCE_DENSITY)
    f = native_fields(src["rho"])
    cvals = [radial_probe(f["c_state"], r) for r in FAR_PROBE_RADII]
    mvals = [radial_probe(f["m10_mag"], r) for r in FAR_PROBE_RADII]
    return {
        "native_mass": src["integrated_source"],
        "R_eff": src["effective_radius"],
        "probe_radii": FAR_PROBE_RADII,
        "c_values": cvals,
        "m10_values": mvals,
        "c_state": loglog_fit(FAR_PROBE_RADII, cvals),
        "m10_mag": loglog_fit(FAR_PROBE_RADII, mvals),
        "finite": f["all_finite"],
    }


def far_mass_test() -> dict:
    # Fixed source radius and fixed far probe. Only native source mass varies.
    probe_r = 12.0
    unit = NATIVE._sphere(FAR_SOURCE_RADIUS, 1.0)
    rows = []
    for m in MASS_LADDER:
        d = m / unit["volume"]
        src = NATIVE._sphere(FAR_SOURCE_RADIUS, d)
        f = native_fields(src["rho"])
        rows.append({
            "native_mass": src["integrated_source"],
            "c_far": radial_probe(f["c_state"], probe_r),
            "m10_far": radial_probe(f["m10_mag"], probe_r),
            "finite": f["all_finite"],
        })
    return {
        "probe_radius": probe_r,
        "rows": rows,
        "c_state": loglog_fit([r["native_mass"] for r in rows], [r["c_far"] for r in rows]),
        "m10_mag": loglog_fit([r["native_mass"] for r in rows], [r["m10_far"] for r in rows]),
    }


def shifted_sphere(center_offset_x: int, radius: float, density: float) -> np.ndarray:
    x0 = float(center_offset_x) * NATIVE.DX
    rr = np.sqrt((NATIVE.XGRID - x0) ** 2 + NATIVE.YGRID ** 2 + NATIVE.ZGRID ** 2)
    rho = np.zeros(NATIVE.SHAPE, dtype=np.float64)
    rho[rr <= radius] = density
    return rho


def relative_rms(residual: np.ndarray, reference: np.ndarray) -> float:
    a = np.asarray(residual, dtype=np.float64)
    b = np.asarray(reference, dtype=np.float64)
    den = float(np.sqrt(np.mean(b * b)))
    num = float(np.sqrt(np.mean(a * a)))
    return num / den if den > 0.0 else float("nan")


def superposition_test() -> dict:
    rho1 = shifted_sphere(-5, 2.5, 0.025)
    rho2 = shifted_sphere(+5, 2.5, 0.05)
    f1 = native_fields(rho1)
    f2 = native_fields(rho2)
    f12 = native_fields(rho1 + rho2)

    csum = f1["c_state"] + f2["c_state"]
    c_rel = relative_rms(f12["c_state"] - csum, csum)

    # Compare signed M10 components, not magnitudes, because magnitudes do not
    # obey vector superposition even when the underlying vector field does.
    component_residuals = {}
    for key in ("vx", "vy", "vz"):
        ref = f1[key] + f2[key]
        component_residuals[key] = relative_rms(f12[key] - ref, ref)
    m10_rel = max(component_residuals.values())

    return {
        "c_state_relative_rms_residual": c_rel,
        "m10_component_relative_rms_residuals": component_residuals,
        "m10_max_component_relative_rms_residual": m10_rel,
        "all_finite": f1["all_finite"] and f2["all_finite"] and f12["all_finite"],
    }


def candidate_summary(name: str, density: dict, mass: dict, radius: dict, far_r: dict, far_m: dict, superpos: dict) -> dict:
    if name == "c_state":
        local_slope = density["c_center"]["slope"]
        surface_mass_slope = mass["c_shell"]["slope"]
        surface_radius_slope = radius["c_shell"]["slope"]
        far_radius_slope = far_r["c_state"]["slope"]
        far_mass_slope = far_m["c_state"]["slope"]
        add_resid = superpos["c_state_relative_rms_residual"]
    elif name == "m10_vector":
        local_slope = density["m10_shell"]["slope"]
        surface_mass_slope = mass["m10_shell"]["slope"]
        surface_radius_slope = radius["m10_shell"]["slope"]
        far_radius_slope = far_r["m10_mag"]["slope"]
        far_mass_slope = far_m["m10_mag"]["slope"]
        add_resid = superpos["m10_max_component_relative_rms_residual"]
    else:
        raise ValueError(name)

    checks = {
        "local_density_exponent_matches_1": exponent_match(local_slope, TARGET_LOCAL_DENSITY_EXPONENT),
        "surface_mass_exponent_matches_1": exponent_match(surface_mass_slope, TARGET_SURFACE_MASS_EXPONENT),
        "surface_radius_fixed_mass_exponent_matches_minus1": exponent_match(surface_radius_slope, TARGET_SURFACE_RADIUS_FIXED_M_EXPONENT),
        "far_mass_exponent_matches_1": exponent_match(far_mass_slope, TARGET_FAR_MASS_EXPONENT),
        "far_radius_exponent_matches_minus1": exponent_match(far_radius_slope, TARGET_FAR_RADIUS_EXPONENT),
        "weak_field_additivity_matches": bool(math.isfinite(add_resid) and add_resid <= ADDITIVITY_REL_TOL),
    }
    return {
        "candidate": name,
        "measured": {
            "local_density_exponent": local_slope,
            "surface_mass_exponent_fixed_R": surface_mass_slope,
            "surface_radius_exponent_fixed_M": surface_radius_slope,
            "far_mass_exponent_fixed_probe_r": far_mass_slope,
            "far_radius_exponent_fixed_source": far_radius_slope,
            "additivity_relative_residual": add_resid,
        },
        "checks": checks,
        "full_fingerprint_match": all(checks.values()),
        "matched_check_count": sum(bool(v) for v in checks.values()),
    }


def main() -> int:
    repo = repo_state()
    if repo["tracked_changes"] or repo["staged_changes"]:
        raise RuntimeError("tracked or staged repository changes present")

    density = density_scaling_test()
    mass = mass_scaling_fixed_radius_test()
    radius = radius_scaling_fixed_mass_test()
    far_r = far_radius_test()
    far_m = far_mass_test()
    superpos = superposition_test()

    candidates = {
        "c_state": candidate_summary("c_state", density, mass, radius, far_r, far_m, superpos),
        "m10_vector": candidate_summary("m10_vector", density, mass, radius, far_r, far_m, superpos),
    }

    full_matches = [k for k, v in candidates.items() if v["full_fingerprint_match"]]
    best_count = max(v["matched_check_count"] for v in candidates.values())
    best_candidates = [k for k, v in candidates.items() if v["matched_check_count"] == best_count]

    if full_matches:
        status = "NATIVE_EFFECTIVE_RESPONSE_FINGERPRINT_MATCH_FOUND"
    else:
        status = "NATIVE_EFFECTIVE_RESPONSE_FINGERPRINT_PARTIAL_MATCH_ONLY"

    all_rows_finite = all(r["finite"] for r in density["rows"] + mass["rows"] + radius["rows"] + far_m["rows"])
    all_finite = all_rows_finite and far_r["finite"] and superpos["all_finite"]

    result = {
        "lab_id": LAB_ID,
        "status": status,
        "repo_state": repo,
        "frozen_effective_fingerprint": {
            "local_density_exponent": TARGET_LOCAL_DENSITY_EXPONENT,
            "surface_mass_exponent_fixed_R": TARGET_SURFACE_MASS_EXPONENT,
            "surface_radius_exponent_fixed_M": TARGET_SURFACE_RADIUS_FIXED_M_EXPONENT,
            "far_mass_exponent_fixed_probe_r": TARGET_FAR_MASS_EXPONENT,
            "far_radius_exponent_fixed_source": TARGET_FAR_RADIUS_EXPONENT,
            "weak_field_additivity_relative_residual": TARGET_ADDITIVITY,
            "exponent_acceptance_window": EXPONENT_WINDOW,
            "additivity_relative_tolerance": ADDITIVITY_REL_TOL,
        },
        "native_tests": {
            "density_scaling": density,
            "mass_scaling_fixed_radius": mass,
            "radius_scaling_fixed_mass": radius,
            "far_radius_scaling": far_r,
            "far_mass_scaling": far_m,
            "superposition": superpos,
        },
        "candidate_comparison": candidates,
        "summary": {
            "full_matches": full_matches,
            "best_candidates_by_matched_checks": best_candidates,
            "best_matched_check_count_of_6": best_count,
            "interpretation": (
                "A structural match licenses only an effective native response role. "
                "It does not establish absolute SI amplitude or identify the microscopic substrate."
            ),
            "safe_next": (
                "If a native candidate matches the frozen scaling fingerprint, audit the absolute amplitude bridge separately without rescaling to the benchmark. "
                "If only partial matches occur, identify which stage of the native pipeline supplies the missing accumulation/propagation behavior before proposing new substrate physics."
            ),
        },
        "policy": {
            "gravity_fundamental_in_PBUF": False,
            "G_used": False,
            "macroscopic_benchmark_amplitude_used": False,
            "native_amplitude_rescaled": False,
            "fit_or_tuning_used": False,
            "microscopic_substrate_claimed": False,
            "Rmax_used": False,
            "cosmology_used": False,
            "lensing_target_used": False,
            "legacy_0p18_used": False,
            "quantum_engine_used": False,
            "planck_input_used": False,
        },
        "checks": {
            "all_native_arrays_finite": all_finite,
            "fixed_native_mass_control": radius["mass_constancy_rel_span"] <= 1.0e-12,
            "no_G": True,
            "no_macroscopic_amplitude_benchmark": True,
            "no_native_rescaling": True,
            "no_fit_or_tuning": True,
            "no_microscopic_substrate_claim": True,
            "no_Rmax": True,
            "no_cosmology": True,
            "no_lensing_target": True,
            "no_legacy_0p18": True,
            "no_quantum_engine": True,
            "no_planck_input": True,
            "no_tracked_or_staged_changes": True,
            "stdout_only_no_run_directory_created": True,
        },
    }

    print(LAB_ID)
    print(f"status={status}")
    print(f"head_sha={repo['head_sha']}")
    print("G_used=false")
    print("macroscopic_benchmark_amplitude_used=false")
    print("native_amplitude_rescaled=false")
    print("fit_or_tuning_used=false")
    print("microscopic_substrate_claimed=false")
    print()

    print("FROZEN_EFFECTIVE_FINGERPRINT")
    print("local_density_exponent=1")
    print("surface_mass_exponent_fixed_R=1")
    print("surface_radius_exponent_fixed_M=-1")
    print("far_mass_exponent_fixed_probe_r=1")
    print("far_radius_exponent_fixed_source=-1")
    print("weak_field_additivity_relative_residual=0")
    print()

    print("NATIVE_CANDIDATE_COMPARISON")
    for name, row in candidates.items():
        m = row["measured"]
        print(
            f"{name} | local_rho={m['local_density_exponent']:.12g} | "
            f"surface_M={m['surface_mass_exponent_fixed_R']:.12g} | "
            f"surface_R_fixed_M={m['surface_radius_exponent_fixed_M']:.12g} | "
            f"far_M={m['far_mass_exponent_fixed_probe_r']:.12g} | "
            f"far_r={m['far_radius_exponent_fixed_source']:.12g} | "
            f"additivity_resid={m['additivity_relative_residual']:.12g} | "
            f"matched={row['matched_check_count']}/6 | full_match={str(row['full_fingerprint_match']).lower()}"
        )
    print()
    print(f"full_matches={','.join(full_matches) if full_matches else 'NONE'}")
    print(f"best_candidates={','.join(best_candidates)}")
    print(f"best_matched_check_count_of_6={best_count}")
    print(f"safe_next={result['summary']['safe_next']}")
    print()

    print("CHECKS")
    for k, v in result["checks"].items():
        print(f"{k}={str(v).lower()}")
    print(json.dumps(result, sort_keys=True, separators=(",", ":")))

    if not all(result["checks"].values()):
        raise RuntimeError("native-response fingerprint lab integrity gate failed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
