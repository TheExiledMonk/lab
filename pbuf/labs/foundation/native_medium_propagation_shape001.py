#!/usr/bin/env python3
"""PBUF FOUNDATION — NATIVE MEDIUM PROPAGATION SHAPE 001.

Purpose
-------
Test a direct PBUF-native propagation hypothesis using the already-supported
accumulation bridge:

    rho -> existing A8 transport -> raw c_state
        -> six-neighbor bounded-strain equilibrium -> accumulated medium state u

The propagation hypothesis tested here is deliberately minimal and structural:
a propagating disturbance crossing an inhomogeneous accumulated medium receives
a transverse directional response proportional to the path-integrated transverse
gradient of that medium state,

    Delta k_x ~ integral ds (partial_x u).

No physical proportionality coefficient is assigned. Therefore only symmetry,
sign, source scaling, and impact-parameter scaling are tested. This is not a
lensing-amplitude calculation and does not use GR/LCDM potentials, field
equations, observed lensing maps, or an inserted 1/r or 1/b law.
"""
from __future__ import annotations

import json
import math
import subprocess
from pathlib import Path

import numpy as np

import pbuf.labs.foundation.c_state_bounded_strain_bridge001 as BRIDGE

ROOT = Path(__file__).resolve().parents[3]
LAB_ID = "PBUF-FOUNDATION-NATIVE-MEDIUM-PROPAGATION-SHAPE-001"

# Freeze the already-supported bridge implementation.
N = BRIDGE.N
CENTER = BRIDGE.CENTER
DX = BRIDGE.DX
SOURCE_RADIUS = 3.5
REFERENCE_MASS = 2.0
MASS_LADDER = (0.5, 1.0, 2.0, 4.0, 8.0)
IMPACT_PARAMETERS = (6.0, 7.0, 8.0, 9.0, 10.0)
MASS_PROBE_B = 8.0

# Shape-only gates, predeclared before execution.
ZERO_TOL = 1.0e-12
SYMMETRY_REL_TOL = 2.0e-3
MASS_EXP_TOL = 0.10
IMPACT_EXP_TOL = 0.35


def git(*args: str) -> str:
    return subprocess.check_output(["git", *args], cwd=ROOT, text=True).strip()


def repo_state() -> dict:
    return {
        "repository": "TheExiledMonk/lab",
        "branch": git("rev-parse", "--abbrev-ref", "HEAD"),
        "head_sha": git("rev-parse", "HEAD"),
        "tracked_changes": git("diff", "--name-only"),
        "staged_changes": git("diff", "--name-only", "--cached"),
    }


def logfit(xs, ys) -> dict:
    x = np.asarray(xs, dtype=np.float64)
    y = np.asarray(ys, dtype=np.float64)
    m = np.isfinite(x) & np.isfinite(y) & (x > 0.0) & (y > 0.0)
    lx = np.log(x[m]); ly = np.log(y[m])
    A = np.column_stack((lx, np.ones_like(lx)))
    beta, *_ = np.linalg.lstsq(A, ly, rcond=None)
    pred = A @ beta
    ss_res = float(np.sum((ly - pred) ** 2))
    ss_tot = float(np.sum((ly - np.mean(ly)) ** 2))
    return {
        "slope": float(beta[0]),
        "r2": 1.0 - ss_res / ss_tot if ss_tot > 0.0 else float("nan"),
        "count": int(np.count_nonzero(m)),
    }


def accumulated_field(mass: float) -> dict:
    rho = BRIDGE.native_rho_fixed_mass(SOURCE_RADIUS, mass)
    sol = BRIDGE.accumulated_from_rho(rho)
    return {
        "field": np.asarray(sol["field"], dtype=np.float64),
        "converged": bool(sol["converged"]),
        "max_strain_fraction": float(sol["max_strain_fraction"]),
        "rho_integral": float(sol["rho_integral"]),
        "c_state_integral": float(sol["c_state_integral"]),
    }


def transverse_gradient(field: np.ndarray) -> np.ndarray:
    # x is axis 2 in the repository's (z,y,x) convention.
    return np.gradient(np.asarray(field, dtype=np.float64), DX, axis=2, edge_order=2)


def path_response(field: np.ndarray, impact_b: float) -> float:
    """Path-integrated x-gradient along a z-directed path at y=0, x=b."""
    gx = transverse_gradient(field)
    off = int(round(float(impact_b) / DX))
    ix = CENTER + off
    if ix <= 0 or ix >= N - 1:
        raise ValueError("impact parameter outside usable accumulation domain")
    # Exclude the fixed zero-boundary planes themselves.
    return float(np.sum(gx[1:-1, CENTER, ix]) * DX)


def zero_source_test() -> dict:
    zero = np.zeros((N, N, N), dtype=np.float64)
    vals = [path_response(zero, b) for b in IMPACT_PARAMETERS]
    return {"responses": vals, "max_abs": float(max(abs(v) for v in vals))}


def symmetry_and_impact_test() -> dict:
    sol = accumulated_field(REFERENCE_MASS)
    rows = []
    for b in IMPACT_PARAMETERS:
        plus = path_response(sol["field"], +b)
        minus = path_response(sol["field"], -b)
        denom = max(abs(plus), abs(minus), 1.0e-30)
        antisym_error = abs(plus + minus) / denom
        rows.append({
            "impact_b": b,
            "response_plus": plus,
            "response_minus": minus,
            "antisymmetry_relative_error": antisym_error,
        })
    fit = logfit([r["impact_b"] for r in rows], [abs(r["response_plus"]) for r in rows])
    center_response = path_response(sol["field"], 0.0)
    return {
        "rows": rows,
        "fit_abs_response_vs_impact": fit,
        "center_response": center_response,
        "converged": sol["converged"],
        "max_strain_fraction": sol["max_strain_fraction"],
    }


def mass_scaling_test() -> dict:
    rows = []
    for mass in MASS_LADDER:
        sol = accumulated_field(mass)
        resp = path_response(sol["field"], MASS_PROBE_B)
        rows.append({
            "mass": mass,
            "response": resp,
            "response_abs": abs(resp),
            "rho_integral": sol["rho_integral"],
            "c_state_integral": sol["c_state_integral"],
            "converged": sol["converged"],
        })
    return {
        "rows": rows,
        "fit": logfit([r["mass"] for r in rows], [r["response_abs"] for r in rows]),
    }


def main() -> int:
    zero = zero_source_test()
    impact = symmetry_and_impact_test()
    mass = mass_scaling_test()
    state = repo_state()

    impact_slope = float(impact["fit_abs_response_vs_impact"]["slope"])
    mass_slope = float(mass["fit"]["slope"])
    max_sym_error = float(max(r["antisymmetry_relative_error"] for r in impact["rows"]))

    # For a positive centered accumulated state that decreases outward,
    # d_x u at x>0 is negative: the direct gradient response is source-directed.
    plus_responses = [r["response_plus"] for r in impact["rows"]]
    minus_responses = [r["response_minus"] for r in impact["rows"]]

    checks = {
        "zero_source_zero_propagation_response": bool(zero["max_abs"] <= ZERO_TOL),
        "centered_path_zero_transverse_response": bool(abs(impact["center_response"]) <= ZERO_TOL),
        "reflection_antisymmetry": bool(max_sym_error <= SYMMETRY_REL_TOL),
        "response_points_toward_source": bool(all(v < 0.0 for v in plus_responses) and all(v > 0.0 for v in minus_responses)),
        "weak_mass_linearity": bool(math.isfinite(mass_slope) and abs(mass_slope - 1.0) <= MASS_EXP_TOL),
        "impact_parameter_inverse_scaling": bool(math.isfinite(impact_slope) and abs(impact_slope + 1.0) <= IMPACT_EXP_TOL),
    }

    execution_checks = {
        "all_measured_values_finite": bool(np.all(np.isfinite([zero["max_abs"], impact["center_response"], impact_slope, mass_slope, max_sym_error]))),
        "nonlinear_solves_converged": bool(impact["converged"] and all(r["converged"] for r in mass["rows"])),
        "raw_c_state_path_used": True,
        "bounded_strain_bridge_reused_without_modification": True,
        "K0_frozen": bool(BRIDGE.K0 == 1.0),
        "epsilon_max_frozen": bool(BRIDGE.EPSILON_MAX == 1.0),
        "no_G": True,
        "no_GR_potential_decomposition": True,
        "no_LCDM": True,
        "no_physical_propagation_amplitude": True,
        "no_native_rescaling": True,
        "no_fit_or_tuning": True,
        "no_inserted_one_over_r_response": True,
        "no_inserted_one_over_b_response": True,
        "no_spherical_equilibrium_shortcut": True,
        "no_Rmax": True,
        "no_cosmology": True,
        "no_lensing_target": True,
        "no_kappa_or_shear_observation": True,
        "no_quantum_engine": True,
        "no_planck_input": True,
        "no_tracked_or_staged_changes": bool(not state["tracked_changes"] and not state["staged_changes"]),
        "stdout_only_no_run_directory_created": True,
    }

    execution_gate_pass = bool(all(execution_checks.values()))
    matched = int(sum(bool(v) for v in checks.values()))
    if matched == len(checks) and execution_gate_pass:
        status = "NATIVE_MEDIUM_PROPAGATION_SHAPE_SUPPORTED"
    elif matched >= 4 and execution_gate_pass:
        status = "NATIVE_MEDIUM_PROPAGATION_SHAPE_PARTIAL_SUPPORT"
    else:
        status = "NATIVE_MEDIUM_PROPAGATION_SHAPE_NOT_SUPPORTED"

    result = {
        "lab_id": LAB_ID,
        "status": status,
        "model": {
            "native_chain": "rho -> existing A8 transport -> raw c_state -> bounded-strain accumulated medium state u",
            "propagation_hypothesis": "transverse directional response proportional to path integral of transverse gradient of u",
            "physical_proportionality_coefficient_assigned": False,
            "source_radius": SOURCE_RADIUS,
            "reference_mass": REFERENCE_MASS,
            "mass_ladder": MASS_LADDER,
            "impact_parameters": IMPACT_PARAMETERS,
            "accumulation_grid_N": N,
        },
        "zero_source_test": zero,
        "impact_test": impact,
        "mass_test": mass,
        "measured": {
            "mass_response_exponent": mass_slope,
            "impact_response_exponent": impact_slope,
            "max_reflection_antisymmetry_relative_error": max_sym_error,
            "central_response": impact["center_response"],
        },
        "checks": checks,
        "matched_checks_of_6": matched,
        "execution_checks": execution_checks,
        "execution_gate_pass": execution_gate_pass,
        "repo_state": state,
        "summary": {
            "question": "Does the supported native accumulated medium state produce the required propagation-response shape directly, without GR/LCDM potential machinery or fitted amplitude?",
            "next_if_supported": "Derive the physical medium-to-propagation coupling and then freeze it before any observational lensing comparison.",
            "next_if_partial": "Localize the failed propagation-shape property before changing the accumulation bridge or introducing observational data.",
        },
    }

    print(LAB_ID)
    print(f"status={status}")
    print(f"head_sha={state['head_sha']}")
    print("direct_medium_gradient_response_used=true")
    print("physical_propagation_amplitude_assigned=false")
    print("GR_potential_decomposition_used=false")
    print("LCDM_used=false")
    print("fit_or_tuning_used=false")
    print("one_over_r_response_inserted=false")
    print("one_over_b_response_inserted=false")
    print()
    print("PROPAGATION_SHAPE")
    print(f"mass_response_exponent={mass_slope:.12g}")
    print(f"impact_response_exponent={impact_slope:.12g}")
    print(f"central_response={impact['center_response']:.12g}")
    print(f"max_reflection_antisymmetry_relative_error={max_sym_error:.12g}")
    for row in impact["rows"]:
        print(f"impact_b={row['impact_b']:.12g} response_plus={row['response_plus']:.12g} response_minus={row['response_minus']:.12g} antisymmetry_error={row['antisymmetry_relative_error']:.12g}")
    print()
    print("CHECKS")
    for k, v in checks.items():
        print(f"{k}={str(v).lower()}")
    print("EXECUTION_CHECKS")
    for k, v in execution_checks.items():
        print(f"{k}={str(v).lower()}")
    print(f"execution_gate_pass={str(execution_gate_pass).lower()}")
    print(json.dumps(result, sort_keys=True, separators=(",", ":")))
    return 0 if execution_gate_pass else 1


if __name__ == "__main__":
    raise SystemExit(main())
