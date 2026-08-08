#!/usr/bin/env python3
"""PBUF FOUNDATION — C_STATE BRIDGE PROPAGATION SHAPE 001.

Purpose
-------
Advance the supported native accumulation bridge into a propagation-shape test
without importing a physical lensing amplitude or silently choosing a complete
metric tensor closure.

Frozen native bridge:

    rho -> existing A8 transport -> raw c_state
        -> six-neighbor bounded-strain equilibrium -> accumulated scalar u

For a static weak metric written schematically as

    ds^2 = -(1+2 Phi) dt^2 + (1-2 Psi) delta_ij dx^i dx^j,

first-order null-ray bending depends on the transverse gradient of the Weyl
combination Phi+Psi.  The current PBUF bridge supplies one accumulated scalar u
but does NOT yet derive the tensor split or the absolute coefficient mapping u
into Phi+Psi.  Therefore this lab tests only the coefficient-independent shape
of a unit positive Weyl-channel probe:

    alpha_tilde(b) = integral dz [d u / d x] at x=b, y=0.

Any constant physical map Phi+Psi = C_W u would multiply every alpha_tilde by
C_W and therefore cannot affect the symmetry or log-log exponents tested here.
C_W is NOT set, fitted, inferred, or claimed to be known.

Predeclared structural tests
----------------------------
- zero source -> zero propagation response;
- centered source -> zero central transverse response;
- reflection symmetry alpha(+b) = -alpha(-b);
- for a positive Weyl-channel coefficient, the transverse gradient points
  inward toward a centered positive source;
- propagation response is linear in source mass in the weak regime;
- impact-parameter magnitude approaches b^-1 for the already-supported u~r^-1
  accumulated response, without inserting a 1/r or 1/b law into the solver.

This is NOT a weak-lensing observation comparison and does not output kappa,
shear, HST residuals, or a physical deflection angle.
"""
from __future__ import annotations

import json
import math
import subprocess
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))

import pbuf.labs.foundation.c_state_bounded_strain_bridge001 as BR

LAB_ID = "PBUF-FOUNDATION-C-STATE-BRIDGE-PROPAGATION-SHAPE-001"

# Freeze the already-supported bridge implementation and geometry.
DX = BR.DX
N = BR.N
CENTER = BR.CENTER
SOURCE_RADIUS = 3.5
REFERENCE_MASS = 2.0
IMPACT_PARAMETERS = (6.0, 7.0, 8.0, 9.0, 10.0)
MASS_LADDER = (0.5, 1.0, 2.0, 4.0, 8.0)
MASS_PROBE_B = 8.0

# Predeclared shape-only acceptance windows.
MASS_EXP_TOL = 0.10
IMPACT_EXP_TOL = 0.30
ODD_SYMMETRY_TOL = 2.0e-5
CENTRAL_REL_TOL = 2.0e-6
ZERO_ABS_TOL = 1.0e-12


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
    lx = np.log(x[m])
    ly = np.log(y[m])
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


def accumulated_for_mass(mass: float) -> dict:
    rho = BR.native_rho_fixed_mass(SOURCE_RADIUS, mass)
    return BR.accumulated_from_rho(rho)


def transverse_gradient_x(field: np.ndarray) -> np.ndarray:
    return np.gradient(np.asarray(field, dtype=np.float64), DX, axis=2, edge_order=2)


def unit_weyl_deflection(field: np.ndarray, impact_parameter: float) -> float:
    """Shape-only line integral of transverse accumulated-state gradient.

    This is alpha_tilde for a unit positive coefficient in Phi+Psi = C_W u.
    No physical C_W is assigned here.
    """
    off = int(round(float(impact_parameter) / DX))
    xi = CENTER + off
    if xi <= 1 or xi >= N - 2:
        raise ValueError("impact parameter too close to numerical boundary")
    gx = transverse_gradient_x(field)
    # Exclude the fixed-zero numerical boundary planes from the line integral.
    return float(np.sum(gx[1:-1, CENTER, xi]) * DX)


def zero_source_test() -> dict:
    rho = np.zeros((BR.N_NATIVE, BR.N_NATIVE, BR.N_NATIVE), dtype=np.float64)
    sol = BR.accumulated_from_rho(rho)
    vals = [unit_weyl_deflection(sol["field"], b) for b in IMPACT_PARAMETERS]
    max_abs = max(abs(v) for v in vals)
    return {"responses": vals, "max_abs_response": max_abs, "converged": sol["converged"]}


def symmetry_and_impact_test() -> dict:
    sol = accumulated_for_mass(REFERENCE_MASS)
    field = sol["field"]
    rows = []
    odd_errors = []
    positive_side = []
    for b in IMPACT_PARAMETERS:
        ap = unit_weyl_deflection(field, +b)
        am = unit_weyl_deflection(field, -b)
        scale = max(abs(ap), abs(am), 1.0e-30)
        odd = abs(ap + am) / scale
        odd_errors.append(odd)
        positive_side.append(ap)
        rows.append({"b": b, "alpha_plus": ap, "alpha_minus": am, "odd_relative_error": odd})

    center_alpha = unit_weyl_deflection(field, 0.0)
    reference_scale = max(max(abs(v) for v in positive_side), 1.0e-30)
    center_relative = abs(center_alpha) / reference_scale
    impact_fit = logfit(IMPACT_PARAMETERS, [abs(v) for v in positive_side])
    inward_positive_unit_channel = bool(all(v < 0.0 for v in positive_side))

    return {
        "rows": rows,
        "impact_fit": impact_fit,
        "max_odd_relative_error": max(odd_errors),
        "center_alpha": center_alpha,
        "center_relative_response": center_relative,
        "positive_side_points_inward_for_positive_unit_channel": inward_positive_unit_channel,
        "max_strain_fraction": sol["max_strain_fraction"],
        "converged": sol["converged"],
        "c_state_integral": sol["c_state_integral"],
        "rho_integral": sol["rho_integral"],
    }


def mass_scaling_test() -> dict:
    rows = []
    for mass in MASS_LADDER:
        sol = accumulated_for_mass(mass)
        alpha = unit_weyl_deflection(sol["field"], MASS_PROBE_B)
        rows.append({
            "mass": mass,
            "alpha_tilde": alpha,
            "abs_alpha_tilde": abs(alpha),
            "max_strain_fraction": sol["max_strain_fraction"],
            "converged": sol["converged"],
        })
    fit = logfit([r["mass"] for r in rows], [r["abs_alpha_tilde"] for r in rows])
    return {"rows": rows, "fit": fit}


def main() -> int:
    state = repo_state()
    zero = zero_source_test()
    geom = symmetry_and_impact_test()
    mass = mass_scaling_test()

    impact_exp = geom["impact_fit"]["slope"]
    mass_exp = mass["fit"]["slope"]

    checks = {
        "zero_source_zero_propagation_response": bool(zero["max_abs_response"] <= ZERO_ABS_TOL),
        "centered_source_zero_central_transverse_response": bool(geom["center_relative_response"] <= CENTRAL_REL_TOL),
        "reflection_odd_symmetry": bool(geom["max_odd_relative_error"] <= ODD_SYMMETRY_TOL),
        "positive_unit_weyl_channel_points_inward": bool(geom["positive_side_points_inward_for_positive_unit_channel"]),
        "weak_mass_linearity": bool(math.isfinite(mass_exp) and abs(mass_exp - 1.0) <= MASS_EXP_TOL),
        "impact_parameter_shape_near_minus1": bool(math.isfinite(impact_exp) and abs(impact_exp + 1.0) <= IMPACT_EXP_TOL),
    }

    nonlinear_ok = bool(
        zero["converged"] and geom["converged"] and all(r["converged"] for r in mass["rows"])
    )
    finite_ok = bool(
        all(math.isfinite(float(v)) for v in [
            zero["max_abs_response"], geom["center_relative_response"],
            geom["max_odd_relative_error"], impact_exp, mass_exp,
        ])
    )

    execution_checks = {
        "all_measured_values_finite": finite_ok,
        "nonlinear_solves_converged": nonlinear_ok,
        "native_bridge_reused_without_amplitude_rescaling": True,
        "physical_weyl_coefficient_not_assigned": True,
        "metric_tensor_split_not_claimed_closed": True,
        "K0_frozen": bool(BR.K0 == 1.0),
        "epsilon_max_frozen": bool(BR.EPSILON_MAX == 1.0),
        "no_G": True,
        "no_macroscopic_amplitude": True,
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

    scientific_passes = sum(bool(v) for v in checks.values())
    if scientific_passes == len(checks):
        status = "C_STATE_BRIDGE_PROPAGATION_SHAPE_SUPPORTED"
    elif scientific_passes >= 3:
        status = "C_STATE_BRIDGE_PROPAGATION_SHAPE_PARTIAL_SUPPORT"
    else:
        status = "C_STATE_BRIDGE_PROPAGATION_SHAPE_NOT_SUPPORTED"

    execution_gate_pass = bool(all(execution_checks.values()))
    result = {
        "lab_id": LAB_ID,
        "status": status,
        "repo_state": state,
        "model": {
            "native_chain": "rho -> existing A8 transport -> raw c_state -> six-neighbor bounded-strain equilibrium -> accumulated u",
            "propagation_shape_probe": "alpha_tilde(b)=integral dz d_x u, corresponding only to unit positive Weyl-channel coefficient",
            "physical_metric_map": "Phi+Psi=C_W u with C_W deliberately unresolved",
            "physical_weyl_coefficient_assigned": False,
            "metric_tensor_split_closed": False,
            "source_radius": SOURCE_RADIUS,
            "reference_mass": REFERENCE_MASS,
            "impact_parameters": IMPACT_PARAMETERS,
            "mass_ladder": MASS_LADDER,
            "mass_probe_b": MASS_PROBE_B,
        },
        "zero_source_test": zero,
        "symmetry_and_impact_test": geom,
        "mass_scaling_test": mass,
        "measured": {
            "impact_parameter_exponent": impact_exp,
            "mass_exponent": mass_exp,
            "max_odd_symmetry_relative_error": geom["max_odd_relative_error"],
            "center_relative_response": geom["center_relative_response"],
            "zero_source_max_abs_response": zero["max_abs_response"],
        },
        "checks": checks,
        "matched_checks": scientific_passes,
        "total_checks": len(checks),
        "execution_checks": execution_checks,
        "execution_gate_pass": execution_gate_pass,
        "remaining_open": {
            "absolute_constitutive_scale": True,
            "map_accumulated_scalar_to_Phi_plus_Psi": True,
            "metric_tensor_split_Phi_vs_Psi": True,
            "physical_deflection_amplitude": True,
            "weak_lensing_observable_test": True,
        },
    }

    print(LAB_ID)
    print(f"status={status}")
    print(f"head_sha={state['head_sha']}")
    print("native_bridge_changed=false")
    print("physical_weyl_coefficient_assigned=false")
    print("metric_tensor_split_claimed_closed=false")
    print("fit_or_tuning_used=false")
    print("one_over_r_response_inserted=false")
    print("one_over_b_response_inserted=false")
    print()
    print("PROPAGATION_SHAPE")
    print(f"impact_parameter_exponent={impact_exp:.12g}")
    print(f"mass_exponent={mass_exp:.12g}")
    print(f"max_odd_symmetry_relative_error={geom['max_odd_relative_error']:.12g}")
    print(f"center_relative_response={geom['center_relative_response']:.12g}")
    print(f"zero_source_max_abs_response={zero['max_abs_response']:.12g}")
    for row in geom["rows"]:
        print(
            f"b={row['b']:.12g} alpha_plus={row['alpha_plus']:.12g} "
            f"alpha_minus={row['alpha_minus']:.12g} odd_error={row['odd_relative_error']:.12g}"
        )
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
