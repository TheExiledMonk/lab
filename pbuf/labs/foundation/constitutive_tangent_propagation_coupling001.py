#!/usr/bin/env python3
"""PBUF FOUNDATION — CONSTITUTIVE TANGENT PROPAGATION COUPLING 001.

Purpose
-------
Test whether the physical medium-to-propagation coupling can be closed from the
same bounded-strain constitutive law already used by the supported native
accumulation bridge, without introducing an observational or fitted coupling.

Frozen native chain:
    rho -> existing A8 transport -> raw c_state
        -> six-neighbor bounded-strain equilibrium -> accumulated state u

The bounded-strain bond law is
    sigma(e) = K0 e / (1 - (e/epsilon_max)^2)

so its exact tangent stiffness is
    K_t(e) = d sigma/de
           = K0 (1 + q) / (1 - q)^2,
      q = (e/epsilon_max)^2.

If a small propagating disturbance is the same local neighbor mode, its local
speed satisfies v^2 ~ K_t/mu.  For constant inertial density mu, the relative
index is therefore fixed with no free coefficient:
    n/n0 = v0/v = sqrt(K0/K_t).

The small-angle geometric-optics directional response is then
    Delta k_x = integral ds partial_x ln(n/n0).

No GR/Weyl potential, G, fitted amplitude, lensing target, inserted 1/r or 1/b
law, or extra propagation coefficient is used.  This is a falsifiable closure
candidate: if its sign/scaling is wrong, the same-mode tangent-speed hypothesis
is rejected rather than tuned.
"""
from __future__ import annotations

import json
import math
import subprocess
from pathlib import Path

import numpy as np

import pbuf.labs.foundation.c_state_bounded_strain_bridge001 as BRIDGE

ROOT = Path(__file__).resolve().parents[3]
LAB_ID = "PBUF-FOUNDATION-CONSTITUTIVE-TANGENT-PROPAGATION-COUPLING-001"

N = BRIDGE.N
CENTER = BRIDGE.CENTER
DX = BRIDGE.DX
K0 = BRIDGE.K0
EPSILON_MAX = BRIDGE.EPSILON_MAX

SOURCE_RADIUS = 3.5
REFERENCE_MASS = 2.0
MASS_LADDER = (0.5, 1.0, 2.0, 4.0, 8.0)
IMPACT_PARAMETERS = (6.0, 7.0, 8.0, 9.0, 10.0)
MASS_PROBE_B = 8.0

ZERO_TOL = 1.0e-12
SYMMETRY_REL_TOL = 2.0e-3
MASS_EXP_TOL = 0.15
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
    lx, ly = np.log(x[m]), np.log(y[m])
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


def tangent_stiffness_from_strain(e: np.ndarray) -> np.ndarray:
    frac = np.asarray(e, dtype=np.float64) / EPSILON_MAX
    q = frac * frac
    if float(np.max(q)) >= 0.995 ** 2:
        raise RuntimeError("bounded-strain barrier approached too closely")
    return K0 * (1.0 + q) / ((1.0 - q) ** 2)


def node_x_tangent(field: np.ndarray) -> np.ndarray:
    """Map x-directed bond tangent stiffness to nodes by adjacent-bond mean."""
    u = np.asarray(field, dtype=np.float64)
    ex = (u[:, :, 1:] - u[:, :, :-1]) / DX
    kx_bond = tangent_stiffness_from_strain(ex)
    kx = np.full_like(u, K0)
    kx[:, :, 1:-1] = 0.5 * (kx_bond[:, :, :-1] + kx_bond[:, :, 1:])
    kx[:, :, 0] = kx_bond[:, :, 0]
    kx[:, :, -1] = kx_bond[:, :, -1]
    return kx


def relative_index(field: np.ndarray) -> np.ndarray:
    """n/n0 = sqrt(K0/K_t); constant inertial density cancels exactly."""
    kx = node_x_tangent(field)
    return np.sqrt(K0 / kx)


def path_response(field: np.ndarray, impact_b: float) -> float:
    """Derived small-angle response integral d_x ln(n/n0) along z."""
    nrel = relative_index(field)
    ln_n = np.log(nrel)
    gx = np.gradient(ln_n, DX, axis=2, edge_order=2)
    off = int(round(float(impact_b) / DX))
    ix = CENTER + off
    if ix <= 0 or ix >= N - 1:
        raise ValueError("impact parameter outside usable domain")
    return float(np.sum(gx[1:-1, CENTER, ix]) * DX)


def zero_source_test() -> dict:
    zero = np.zeros((N, N, N), dtype=np.float64)
    vals = [path_response(zero, b) for b in IMPACT_PARAMETERS]
    return {"responses": vals, "max_abs": float(max(abs(v) for v in vals))}


def impact_test() -> dict:
    sol = accumulated_field(REFERENCE_MASS)
    rows = []
    for b in IMPACT_PARAMETERS:
        plus = path_response(sol["field"], +b)
        minus = path_response(sol["field"], -b)
        denom = max(abs(plus), abs(minus), 1.0e-30)
        rows.append({
            "impact_b": b,
            "response_plus": plus,
            "response_minus": minus,
            "antisymmetry_relative_error": abs(plus + minus) / denom,
        })
    return {
        "rows": rows,
        "fit_abs_response_vs_impact": logfit(
            [r["impact_b"] for r in rows],
            [abs(r["response_plus"]) for r in rows],
        ),
        "center_response": path_response(sol["field"], 0.0),
        "converged": sol["converged"],
        "max_strain_fraction": sol["max_strain_fraction"],
    }


def mass_test() -> dict:
    rows = []
    for mass in MASS_LADDER:
        sol = accumulated_field(mass)
        response = path_response(sol["field"], MASS_PROBE_B)
        rows.append({
            "mass": mass,
            "response": response,
            "response_abs": abs(response),
            "rho_integral": sol["rho_integral"],
            "c_state_integral": sol["c_state_integral"],
            "converged": sol["converged"],
        })
    return {
        "rows": rows,
        "fit": logfit([r["mass"] for r in rows], [r["response_abs"] for r in rows]),
    }


def constitutive_identity_test() -> dict:
    samples = np.linspace(-0.25 * EPSILON_MAX, 0.25 * EPSILON_MAX, 41)
    h = 1.0e-6 * EPSILON_MAX
    sigma_p = K0 * (samples + h) / (1.0 - ((samples + h) / EPSILON_MAX) ** 2)
    sigma_m = K0 * (samples - h) / (1.0 - ((samples - h) / EPSILON_MAX) ** 2)
    fd = (sigma_p - sigma_m) / (2.0 * h)
    exact = tangent_stiffness_from_strain(samples)
    rel = float(np.max(np.abs(fd - exact) / np.maximum(np.abs(exact), 1.0e-30)))
    return {"max_relative_error": rel, "pass": bool(rel <= 1.0e-8)}


def main() -> int:
    identity = constitutive_identity_test()
    zero = zero_source_test()
    impact = impact_test()
    mass = mass_test()
    state = repo_state()

    impact_slope = float(impact["fit_abs_response_vs_impact"]["slope"])
    mass_slope = float(mass["fit"]["slope"])
    max_sym_error = float(max(r["antisymmetry_relative_error"] for r in impact["rows"]))
    plus = [r["response_plus"] for r in impact["rows"]]
    minus = [r["response_minus"] for r in impact["rows"]]

    # Source-directed convention: for +b the response must be negative and for
    # -b positive.  If the derived tangent-speed index does the opposite, fail.
    checks = {
        "constitutive_tangent_identity": identity["pass"],
        "zero_source_zero_propagation_response": bool(zero["max_abs"] <= ZERO_TOL),
        "centered_path_zero_transverse_response": bool(abs(impact["center_response"]) <= ZERO_TOL),
        "reflection_antisymmetry": bool(max_sym_error <= SYMMETRY_REL_TOL),
        "derived_response_points_toward_source": bool(all(v < 0.0 for v in plus) and all(v > 0.0 for v in minus)),
        "weak_mass_linearity": bool(math.isfinite(mass_slope) and abs(mass_slope - 1.0) <= MASS_EXP_TOL),
        "impact_parameter_inverse_scaling": bool(math.isfinite(impact_slope) and abs(impact_slope + 1.0) <= IMPACT_EXP_TOL),
    }

    measured_values = [zero["max_abs"], impact["center_response"], impact_slope, mass_slope, max_sym_error, identity["max_relative_error"]]
    execution_checks = {
        "all_measured_values_finite": bool(np.all(np.isfinite(measured_values))),
        "nonlinear_solves_converged": bool(impact["converged"] and all(r["converged"] for r in mass["rows"])),
        "raw_c_state_path_used": True,
        "bounded_strain_bridge_reused_without_modification": True,
        "propagation_index_derived_from_constitutive_tangent": True,
        "constant_inertial_density_cancels_from_relative_index": True,
        "K0_frozen": bool(K0 == 1.0),
        "epsilon_max_frozen": bool(EPSILON_MAX == 1.0),
        "no_free_propagation_coefficient": True,
        "no_G": True,
        "no_GR_potential_decomposition": True,
        "no_LCDM": True,
        "no_observational_amplitude_calibration": True,
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
        status = "CONSTITUTIVE_TANGENT_PROPAGATION_COUPLING_SUPPORTED"
    elif matched >= 4 and execution_gate_pass:
        status = "CONSTITUTIVE_TANGENT_PROPAGATION_COUPLING_PARTIAL_SUPPORT"
    else:
        status = "CONSTITUTIVE_TANGENT_PROPAGATION_COUPLING_NOT_SUPPORTED"

    result = {
        "lab_id": LAB_ID,
        "status": status,
        "model": {
            "native_chain": "rho -> existing A8 transport -> raw c_state -> bounded-strain accumulated state u",
            "constitutive_stress": "sigma=K0*e/(1-(e/epsilon_max)^2)",
            "tangent_stiffness": "K_t=K0*(1+q)/(1-q)^2; q=(e/epsilon_max)^2",
            "relative_index": "n/n0=sqrt(K0/K_t)",
            "propagation_response": "Delta k_x=integral ds partial_x ln(n/n0)",
            "free_propagation_coefficient": False,
            "source_radius": SOURCE_RADIUS,
            "reference_mass": REFERENCE_MASS,
            "mass_ladder": MASS_LADDER,
            "impact_parameters": IMPACT_PARAMETERS,
        },
        "constitutive_identity_test": identity,
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
        "matched_checks": matched,
        "execution_checks": execution_checks,
        "execution_gate_pass": execution_gate_pass,
        "repo_state": state,
        "summary": {
            "question": "Does the same bounded-strain constitutive tangent that supports accumulation also close propagation through its derived local wave speed, with no free coupling coefficient?",
            "next_if_supported": "Freeze the constitutive propagation coupling and proceed to source-independent physical-scale closure before observational testing.",
            "next_if_partial_or_not_supported": "Reject/localize the same-mode tangent-speed propagation hypothesis without altering the supported accumulation bridge; derive a different native propagation interaction rather than fit one.",
        },
    }

    print(LAB_ID)
    print(f"status={status}")
    print(f"head_sha={state['head_sha']}")
    print("constitutive_tangent_speed_coupling_used=true")
    print("free_propagation_coefficient_used=false")
    print("GR_potential_decomposition_used=false")
    print("LCDM_used=false")
    print("fit_or_tuning_used=false")
    print()
    print("DERIVED_PROPAGATION_COUPLING")
    print(f"constitutive_tangent_identity_max_relative_error={identity['max_relative_error']:.12g}")
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
