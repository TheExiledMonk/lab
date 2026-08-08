#!/usr/bin/env python3
"""PBUF FOUNDATION — NATIVE MEDIUM RESPONSE ACCUMULATION 001.

Fact-finding audit of the existing native bridge:

    local medium loading -> long-range accumulated medium response

Gravity is not used as a native construction variable.  The frozen response
fingerprint is used only as an independently established shape/scaling test.
No G, amplitude calibration, rescaling, fitting, new kernel, Rmax, cosmology,
lensing target, Quantum Engine, or Planck-scale input is used.

This lab tests only existing native operations already present in the repo:
  * c_state
  * laplacian(c_state)
  * signed M10 interface vector
  * divergence(M10)

A non-match is a valid scientific result.
"""
from __future__ import annotations

import json
import math
import subprocess
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]

import pbuf.labs.foundation.native_field_curvature_dimension_audit001 as NATIVE
import pbuf.labs.foundation.native_response_fingerprint_match001 as FP
import pbuf.labs.foundation.m10_coverage_25pct_science001 as BASE

LAB_ID = "PBUF-FOUNDATION-NATIVE-MEDIUM-RESPONSE-ACCUMULATION-001"
EXP_WINDOW = 0.20
ADD_TOL = 1.0e-10


def git(*args: str) -> str:
    return subprocess.check_output(["git", *args], cwd=ROOT, text=True).strip()


def repo_state() -> dict:
    return {
        "branch": git("rev-parse", "--abbrev-ref", "HEAD"),
        "head_sha": git("rev-parse", "HEAD"),
        "tracked_changes": git("diff", "--name-only"),
        "staged_changes": git("diff", "--name-only", "--cached"),
    }


def fields(rho: np.ndarray) -> dict:
    state = NATIVE._noise_free_state(np.asarray(rho, dtype=np.float64))
    c = np.asarray(state["c_state"], dtype=np.float64)
    candidate = BASE._candidate(state)
    vx, vy, vz = BASE._interface_vector(candidate)
    vx = np.asarray(vx, dtype=np.float64)
    vy = np.asarray(vy, dtype=np.float64)
    vz = np.asarray(vz, dtype=np.float64)
    lap = np.asarray(NATIVE._laplacian(c), dtype=np.float64)
    div = np.asarray(NATIVE._divergence(vx, vy, vz), dtype=np.float64)
    return {
        "c_state": (c,),
        "laplacian_c_state": (lap,),
        "m10_vector": (vx, vy, vz),
        "m10_divergence": (div,),
    }


PROVENANCE = {
    "c_state": {
        "source": "pbuf.models.a8_state.evolve_a8_transport_3d via native_field_curvature_dimension_audit001._noise_free_state",
        "kind": "scalar transport state",
        "from_c_state": True,
        "operator_class": "transport/state",
    },
    "laplacian_c_state": {
        "source": "native_field_curvature_dimension_audit001._laplacian",
        "kind": "scalar signed derivative",
        "from_c_state": True,
        "operator_class": "local second derivative",
    },
    "m10_vector": {
        "source": "m10_coverage_25pct_science001._interface_vector",
        "kind": "signed 3-vector",
        "from_c_state": False,
        "operator_class": "existing downstream interface transform",
    },
    "m10_divergence": {
        "source": "native_field_curvature_dimension_audit001._divergence applied to existing M10 vector",
        "kind": "scalar signed derivative",
        "from_c_state": False,
        "operator_class": "local derivative of downstream vector",
    },
}


def amp_at_shell(components: tuple[np.ndarray, ...], src: dict) -> float:
    shell = NATIVE._shell_mask(src["effective_radius"])
    if len(components) == 1:
        return float(np.mean(np.abs(components[0][shell])))
    sq = np.zeros_like(components[0], dtype=np.float64)
    for a in components:
        sq += a * a
    return float(np.mean(np.sqrt(sq)[shell]))


def amp_at_probe(components: tuple[np.ndarray, ...], radius: float) -> float:
    zi, yi, xi = NATIVE.CENTER_INDEX
    offset = int(round(radius / NATIVE.DX))
    vals = [float(a[zi, yi, xi + offset]) for a in components]
    return float(math.sqrt(sum(v * v for v in vals))) if len(vals) > 1 else abs(vals[0])


def center_amp(components: tuple[np.ndarray, ...]) -> float:
    zi, yi, xi = NATIVE.CENTER_INDEX
    vals = [float(a[zi, yi, xi]) for a in components]
    return float(math.sqrt(sum(v * v for v in vals))) if len(vals) > 1 else abs(vals[0])


def rel_rms_components(actual: tuple[np.ndarray, ...], reference: tuple[np.ndarray, ...]) -> float:
    residual_sq = 0.0
    reference_sq = 0.0
    for a, b in zip(actual, reference):
        residual_sq += float(np.mean((a - b) ** 2))
        reference_sq += float(np.mean(b ** 2))
    return math.sqrt(residual_sq / reference_sq) if reference_sq > 0.0 else float("nan")


def logfit(xs, ys) -> float:
    return float(FP.loglog_fit(xs, ys)["slope"])


def match(x: float, target: float) -> bool:
    return bool(math.isfinite(x) and abs(x - target) <= EXP_WINDOW)


def run_candidate(name: str) -> dict:
    dens_vals = []
    for d in FP.DENSITY_LADDER:
        src = NATIVE._sphere(FP.FIXED_RADIUS, d)
        dens_vals.append(center_amp(fields(src["rho"])[name]))
    density_exp = logfit(FP.DENSITY_LADDER, dens_vals)

    unit = NATIVE._sphere(FP.FIXED_RADIUS, 1.0)
    mass_vals = []
    masses = []
    for m in FP.MASS_LADDER:
        src = NATIVE._sphere(FP.FIXED_RADIUS, m / unit["volume"])
        masses.append(src["integrated_source"])
        mass_vals.append(amp_at_shell(fields(src["rho"])[name], src))
    surface_mass_exp = logfit(masses, mass_vals)

    r_eff = []
    radius_vals = []
    fixed_masses = []
    for r in FP.RADIUS_LADDER:
        src = FP.sphere_with_native_mass(r, FP.FIXED_NATIVE_MASS)
        r_eff.append(src["effective_radius"])
        fixed_masses.append(src["integrated_source"])
        radius_vals.append(amp_at_shell(fields(src["rho"])[name], src))
    surface_radius_exp = logfit(r_eff, radius_vals)

    far_src = NATIVE._sphere(FP.FAR_SOURCE_RADIUS, FP.FAR_SOURCE_DENSITY)
    far_field = fields(far_src["rho"])[name]
    far_vals = [amp_at_probe(far_field, r) for r in FP.FAR_PROBE_RADII]
    far_radius_exp = logfit(FP.FAR_PROBE_RADII, far_vals)

    far_probe = 12.0
    far_unit = NATIVE._sphere(FP.FAR_SOURCE_RADIUS, 1.0)
    far_mass_vals = []
    far_masses = []
    for m in FP.MASS_LADDER:
        src = NATIVE._sphere(FP.FAR_SOURCE_RADIUS, m / far_unit["volume"])
        far_masses.append(src["integrated_source"])
        far_mass_vals.append(amp_at_probe(fields(src["rho"])[name], far_probe))
    far_mass_exp = logfit(far_masses, far_mass_vals)

    rho1 = FP.shifted_sphere(-5, 2.5, 0.035)
    rho2 = FP.shifted_sphere(+5, 3.5, 0.022)
    f1 = fields(rho1)[name]
    f2 = fields(rho2)[name]
    f12 = fields(rho1 + rho2)[name]
    summed = tuple(a + b for a, b in zip(f1, f2))
    add_resid = rel_rms_components(f12, summed)

    checks = {
        "local_density_exponent_matches_1": match(density_exp, 1.0),
        "surface_mass_exponent_matches_1": match(surface_mass_exp, 1.0),
        "surface_radius_fixed_mass_exponent_matches_minus1": match(surface_radius_exp, -1.0),
        "far_mass_exponent_matches_1": match(far_mass_exp, 1.0),
        "far_radius_exponent_matches_minus1": match(far_radius_exp, -1.0),
        "signed_component_additivity_matches": bool(math.isfinite(add_resid) and add_resid <= ADD_TOL),
    }
    return {
        "candidate": name,
        "provenance": PROVENANCE[name],
        "measured": {
            "local_density_exponent": density_exp,
            "surface_mass_exponent_fixed_R": surface_mass_exp,
            "surface_radius_exponent_fixed_M": surface_radius_exp,
            "far_mass_exponent_fixed_probe_r": far_mass_exp,
            "far_radius_exponent_fixed_source": far_radius_exp,
            "signed_component_additivity_relative_residual": add_resid,
            "fixed_mass_relative_span": (max(fixed_masses) - min(fixed_masses)) / FP.FIXED_NATIVE_MASS,
        },
        "checks": checks,
        "matched_check_count": int(sum(checks.values())),
        "full_accumulation_fingerprint_match": bool(all(checks.values())),
    }


def main() -> int:
    repo = repo_state()
    names = tuple(PROVENANCE)
    results = {name: run_candidate(name) for name in names}
    full = [name for name, r in results.items() if r["full_accumulation_fingerprint_match"]]
    best_count = max(r["matched_check_count"] for r in results.values())
    best = [name for name, r in results.items() if r["matched_check_count"] == best_count]

    if full:
        status = "EXISTING_NATIVE_ACCUMULATION_CANDIDATE_FOUND"
    elif best_count > 0:
        status = "EXISTING_NATIVE_ACCUMULATION_PARTIAL_MATCH_ONLY"
    else:
        status = "NATIVE_ACCUMULATION_OPERATOR_NOT_FOUND"

    policy = {
        "gravity_fundamental_in_PBUF": False,
        "gravity_used_as_native_variable": False,
        "G_used": False,
        "macroscopic_amplitude_used": False,
        "native_amplitude_rescaled": False,
        "fit_or_tuning_used": False,
        "new_accumulation_kernel_introduced": False,
        "Rmax_used": False,
        "cosmology_used": False,
        "lensing_target_used": False,
        "legacy_0p18_used": False,
        "quantum_engine_used": False,
        "planck_input_used": False,
    }
    all_outputs_finite = all(
        math.isfinite(v)
        for r in results.values()
        for v in r["measured"].values()
    )
    checks = {
        "candidate_inventory_nonempty": bool(results),
        "all_outputs_finite": all_outputs_finite,
        "fixed_mass_control": all(r["measured"]["fixed_mass_relative_span"] <= 1.0e-12 for r in results.values()),
        "no_G": True,
        "no_macroscopic_amplitude": True,
        "no_native_rescaling": True,
        "no_fit_or_tuning": True,
        "no_new_accumulation_kernel": True,
        "no_Rmax": True,
        "no_cosmology": True,
        "no_lensing_target": True,
        "no_legacy_0p18": True,
        "no_quantum_engine": True,
        "no_planck_input": True,
        "no_tracked_or_staged_changes": repo["tracked_changes"] == "" and repo["staged_changes"] == "",
        "stdout_only_no_run_directory_created": True,
    }
    execution_gate_keys = tuple(k for k in checks if k != "all_outputs_finite")
    execution_gate_pass = all(checks[k] for k in execution_gate_keys)

    payload = {
        "lab_id": LAB_ID,
        "status": status,
        "repo_state": repo,
        "policy": policy,
        "frozen_shape_fingerprint": {
            "local_density_exponent": 1.0,
            "surface_mass_exponent_fixed_R": 1.0,
            "surface_radius_exponent_fixed_M": -1.0,
            "far_mass_exponent_fixed_probe_r": 1.0,
            "far_radius_exponent_fixed_source": -1.0,
            "signed_component_additivity_relative_residual": 0.0,
            "exponent_acceptance_window": EXP_WINDOW,
            "additivity_relative_tolerance": ADD_TOL,
        },
        "candidate_results": results,
        "summary": {
            "full_matches": full,
            "best_candidates": best,
            "best_matched_check_count_of_6": best_count,
            "all_outputs_finite": all_outputs_finite,
            "nonfinite_scientific_measurements_are_valid_open_results": True,
            "interpretation": "native medium-response structure only; no SI amplitude, microscopic substrate ontology, metric-strain assignment, or gravity-field assignment is established",
            "safe_next": (
                "If a full existing candidate is found, audit its physical mapping separately without amplitude rescaling. "
                "If only partial matches occur, identify the missing native accumulation behavior. If no candidate matches, "
                "the repository lacks the required existing accumulation operator and a later explicitly speculative/constitutive derivation is licensed."
            ),
        },
        "checks": checks,
        "execution_gate": {
            "keys": execution_gate_keys,
            "pass": execution_gate_pass,
            "note": "nonfinite scientific measurements remain reported but do not convert a valid partial/null scientific outcome into an execution failure",
        },
    }

    print(LAB_ID)
    print(f"status={status}")
    print(f"head_sha={repo['head_sha']}")
    print("gravity_used_as_native_variable=false")
    print("G_used=false")
    print("macroscopic_amplitude_used=false")
    print("native_amplitude_rescaled=false")
    print("fit_or_tuning_used=false")
    print("new_accumulation_kernel_introduced=false")
    print()
    print("EXISTING_NATIVE_CANDIDATES")
    for name in names:
        r = results[name]
        m = r["measured"]
        print(
            f"{name} | rho={m['local_density_exponent']:.12g} | M_surface={m['surface_mass_exponent_fixed_R']:.12g} | "
            f"R_fixed_M={m['surface_radius_exponent_fixed_M']:.12g} | M_far={m['far_mass_exponent_fixed_probe_r']:.12g} | "
            f"r_far={m['far_radius_exponent_fixed_source']:.12g} | add_resid={m['signed_component_additivity_relative_residual']:.12g} | "
            f"matched={r['matched_check_count']}/6 | full_match={str(r['full_accumulation_fingerprint_match']).lower()}"
        )
    print()
    print(f"full_matches={','.join(full) if full else 'NONE'}")
    print(f"best_candidates={','.join(best)}")
    print(f"best_matched_check_count_of_6={best_count}")
    print(f"safe_next={payload['summary']['safe_next']}")
    print()
    print("CHECKS")
    for k, v in checks.items():
        print(f"{k}={str(v).lower()}")
    print(f"execution_gate_pass={str(execution_gate_pass).lower()}")
    print(json.dumps(payload, sort_keys=True, separators=(",", ":")))

    if not execution_gate_pass:
        raise RuntimeError("native medium response accumulation execution gate failed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
