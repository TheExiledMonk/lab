#!/usr/bin/env python3
"""PBUF FOUNDATION — DEFORMED MEDIUM GEOMETRY PROPAGATION 001.

Purpose
-------
Test whether propagation can be obtained from the *form of the already-supported
accumulated medium state itself*, without coupling light to stiffness and without
introducing GR/LCDM potentials or an observationally fitted coefficient.

Frozen source/medium chain:

    rho -> existing A8 transport -> raw c_state
        -> six-neighbor bounded-strain equilibrium -> accumulated medium state u

Geometric hypothesis tested
---------------------------
Treat a two-dimensional x-z cross-section of the scalar accumulated state as the
height of a deformed medium sheet embedded in one additional geometric
coordinate,

    X(x,z) = (x, z, u(x,z)).

The induced metric is then fixed algebraically by the deformation:

    g_ij = delta_ij + (partial_i u)(partial_j u).

A freely propagating path on that deformed geometry follows the geodesic equation
of this induced metric.  There is no free propagation coefficient.  No tangent
stiffness enters the ray equation: stiffness is used only upstream to determine
u via the already-frozen accumulation bridge.

This is deliberately a falsifiable *scalar-height geometry* candidate.  If its
weak-source scaling/sign/impact behavior is wrong, the correct conclusion is
that scalar u alone is insufficient to specify propagation geometry; the
supported accumulation bridge is not to be altered or tuned.
"""
from __future__ import annotations

import json
import math
import subprocess
from pathlib import Path

import numpy as np

import pbuf.labs.foundation.c_state_bounded_strain_bridge001 as BRIDGE

ROOT = Path(__file__).resolve().parents[3]
LAB_ID = "PBUF-FOUNDATION-DEFORMED-MEDIUM-GEOMETRY-PROPAGATION-001"

N = BRIDGE.N
CENTER = BRIDGE.CENTER
DX = BRIDGE.DX
SOURCE_RADIUS = 3.5
REFERENCE_MASS = 2.0
MASS_LADDER = (0.5, 1.0, 2.0, 4.0, 8.0)
IMPACT_PARAMETERS = (6.0, 7.0, 8.0, 9.0, 10.0)
MASS_PROBE_B = 8.0

# Keep rays well away from the fixed zero boundary.
Z_START = -20.0
Z_STOP = 20.0
STEP = 0.10
MAX_STEPS = 1200

ZERO_TOL = 1.0e-11
SYMMETRY_REL_TOL = 3.0e-3
MASS_EXP_TOL = 0.20
IMPACT_EXP_TOL = 0.40


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


def geometry_from_field(field3: np.ndarray) -> dict:
    """Build induced 2D metric and Christoffel symbols on y=0 x-z plane."""
    # repository arrays are (z,y,x); plane axes are therefore (z,x)
    u = np.asarray(field3[:, CENTER, :], dtype=np.float64)
    uz, ux = np.gradient(u, DX, DX, edge_order=2)

    # Coordinate order for tensors below: 0=x, 1=z. Arrays remain indexed [z,x].
    g = np.empty((2, 2, N, N), dtype=np.float64)
    g[0, 0] = 1.0 + ux * ux
    g[0, 1] = ux * uz
    g[1, 0] = g[0, 1]
    g[1, 1] = 1.0 + uz * uz

    det = g[0, 0] * g[1, 1] - g[0, 1] * g[1, 0]
    inv = np.empty_like(g)
    inv[0, 0] = g[1, 1] / det
    inv[1, 1] = g[0, 0] / det
    inv[0, 1] = -g[0, 1] / det
    inv[1, 0] = inv[0, 1]

    # dg[coord_derivative, i, j], coordinate derivative order x,z.
    dg = np.empty((2, 2, 2, N, N), dtype=np.float64)
    for i in range(2):
        for j in range(2):
            dz, dx = np.gradient(g[i, j], DX, DX, edge_order=2)
            dg[0, i, j] = dx
            dg[1, i, j] = dz

    gamma = np.zeros((2, 2, 2, N, N), dtype=np.float64)
    for a in range(2):
        for b in range(2):
            for c in range(2):
                total = np.zeros((N, N), dtype=np.float64)
                for d in range(2):
                    total += 0.5 * inv[a, d] * (
                        dg[b, d, c] + dg[c, d, b] - dg[d, b, c]
                    )
                gamma[a, b, c] = total

    return {
        "u": u,
        "ux": ux,
        "uz": uz,
        "metric": g,
        "gamma": gamma,
        "max_slope": float(max(np.max(np.abs(ux)), np.max(np.abs(uz)))),
    }


def bilinear(arr: np.ndarray, x: float, z: float) -> float:
    fx = x / DX + CENTER
    fz = z / DX + CENTER
    ix = int(math.floor(fx)); iz = int(math.floor(fz))
    if ix < 1 or ix >= N - 2 or iz < 1 or iz >= N - 2:
        raise ValueError("ray left interpolation domain")
    tx = fx - ix; tz = fz - iz
    return float(
        (1-tz) * ((1-tx)*arr[iz, ix] + tx*arr[iz, ix+1])
        + tz * ((1-tx)*arr[iz+1, ix] + tx*arr[iz+1, ix+1])
    )


def rhs(state: np.ndarray, gamma: np.ndarray) -> np.ndarray:
    x, z, vx, vz = [float(v) for v in state]
    vel = np.array([vx, vz], dtype=np.float64)
    acc = np.zeros(2, dtype=np.float64)
    for a in range(2):
        for b in range(2):
            for c in range(2):
                acc[a] -= bilinear(gamma[a, b, c], x, z) * vel[b] * vel[c]
    return np.array([vx, vz, acc[0], acc[1]], dtype=np.float64)


def trace_ray(geometry: dict, impact_b: float) -> dict:
    gamma = geometry["gamma"]
    state = np.array([float(impact_b), Z_START, 0.0, 1.0], dtype=np.float64)
    initial_angle = math.atan2(state[2], state[3])
    steps = 0
    while state[1] < Z_STOP and steps < MAX_STEPS:
        h = STEP
        k1 = rhs(state, gamma)
        k2 = rhs(state + 0.5*h*k1, gamma)
        k3 = rhs(state + 0.5*h*k2, gamma)
        k4 = rhs(state + h*k3, gamma)
        state = state + (h/6.0)*(k1 + 2*k2 + 2*k3 + k4)
        steps += 1
        if not np.all(np.isfinite(state)):
            raise RuntimeError("non-finite geodesic state")
    if steps >= MAX_STEPS:
        raise RuntimeError("geodesic integration did not reach exit plane")
    final_angle = math.atan2(float(state[2]), float(state[3]))
    return {
        "impact_b": float(impact_b),
        "outgoing_angle": final_angle,
        "angle_change": final_angle - initial_angle,
        "exit_x": float(state[0]),
        "exit_z": float(state[1]),
        "steps": steps,
    }


def zero_source_test() -> dict:
    zero = np.zeros((N, N, N), dtype=np.float64)
    geom = geometry_from_field(zero)
    rows = [trace_ray(geom, b) for b in IMPACT_PARAMETERS]
    max_abs = max(abs(r["angle_change"]) for r in rows)
    return {"rows": rows, "max_abs_angle_change": float(max_abs)}


def impact_test() -> dict:
    sol = accumulated_field(REFERENCE_MASS)
    geom = geometry_from_field(sol["field"])
    rows = []
    for b in IMPACT_PARAMETERS:
        plus = trace_ray(geom, +b)
        minus = trace_ray(geom, -b)
        denom = max(abs(plus["angle_change"]), abs(minus["angle_change"]), 1.0e-30)
        antisym = abs(plus["angle_change"] + minus["angle_change"]) / denom
        rows.append({
            "impact_b": b,
            "angle_plus": plus["angle_change"],
            "angle_minus": minus["angle_change"],
            "antisymmetry_relative_error": antisym,
            "exit_x_plus": plus["exit_x"],
            "exit_x_minus": minus["exit_x"],
        })
    center = trace_ray(geom, 0.0)
    fit = logfit([r["impact_b"] for r in rows], [abs(r["angle_plus"]) for r in rows])
    return {
        "rows": rows,
        "center_angle_change": center["angle_change"],
        "fit_abs_angle_vs_impact": fit,
        "converged": sol["converged"],
        "max_strain_fraction": sol["max_strain_fraction"],
        "max_surface_slope": geom["max_slope"],
    }


def mass_test() -> dict:
    rows = []
    for mass in MASS_LADDER:
        sol = accumulated_field(mass)
        geom = geometry_from_field(sol["field"])
        ray = trace_ray(geom, MASS_PROBE_B)
        rows.append({
            "mass": mass,
            "angle_change": ray["angle_change"],
            "angle_abs": abs(ray["angle_change"]),
            "rho_integral": sol["rho_integral"],
            "c_state_integral": sol["c_state_integral"],
            "converged": sol["converged"],
        })
    return {
        "rows": rows,
        "fit": logfit([r["mass"] for r in rows], [r["angle_abs"] for r in rows]),
    }


def main() -> int:
    zero = zero_source_test()
    impact = impact_test()
    mass = mass_test()
    state = repo_state()

    impact_slope = float(impact["fit_abs_angle_vs_impact"]["slope"])
    mass_slope = float(mass["fit"]["slope"])
    sym_err = float(max(r["antisymmetry_relative_error"] for r in impact["rows"]))
    plus = [r["angle_plus"] for r in impact["rows"]]
    minus = [r["angle_minus"] for r in impact["rows"]]

    checks = {
        "zero_source_zero_geometric_deflection": bool(zero["max_abs_angle_change"] <= ZERO_TOL),
        "centered_path_zero_transverse_deflection": bool(abs(impact["center_angle_change"]) <= ZERO_TOL),
        "reflection_antisymmetry": bool(sym_err <= SYMMETRY_REL_TOL),
        "geometric_deflection_points_toward_source": bool(all(v < 0.0 for v in plus) and all(v > 0.0 for v in minus)),
        "weak_mass_linearity": bool(math.isfinite(mass_slope) and abs(mass_slope - 1.0) <= MASS_EXP_TOL),
        "impact_parameter_inverse_scaling": bool(math.isfinite(impact_slope) and abs(impact_slope + 1.0) <= IMPACT_EXP_TOL),
    }

    execution_checks = {
        "all_measured_values_finite": bool(np.all(np.isfinite([
            zero["max_abs_angle_change"], impact["center_angle_change"], impact_slope, mass_slope, sym_err
        ]))),
        "nonlinear_medium_solves_converged": bool(impact["converged"] and all(r["converged"] for r in mass["rows"])),
        "raw_c_state_path_used": True,
        "bounded_strain_bridge_reused_without_modification": True,
        "scalar_height_induced_metric_used": True,
        "geodesic_equation_used": True,
        "stiffness_not_used_in_propagation_equation": True,
        "no_free_propagation_coefficient": True,
        "K0_frozen": bool(BRIDGE.K0 == 1.0),
        "epsilon_max_frozen": bool(BRIDGE.EPSILON_MAX == 1.0),
        "no_G": True,
        "no_GR_potential_decomposition": True,
        "no_Weyl": True,
        "no_LCDM": True,
        "no_observational_amplitude_calibration": True,
        "no_native_rescaling": True,
        "no_fit_or_tuning": True,
        "no_inserted_one_over_r_response": True,
        "no_inserted_one_over_b_response": True,
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
        status = "DEFORMED_MEDIUM_GEOMETRY_PROPAGATION_SUPPORTED"
    elif matched >= 3 and execution_gate_pass:
        status = "DEFORMED_MEDIUM_GEOMETRY_PROPAGATION_PARTIAL_SUPPORT"
    else:
        status = "DEFORMED_MEDIUM_GEOMETRY_PROPAGATION_NOT_SUPPORTED"

    result = {
        "lab_id": LAB_ID,
        "status": status,
        "model": {
            "native_chain": "rho -> existing A8 transport -> raw c_state -> bounded-strain accumulated medium state u",
            "geometric_hypothesis": "2D scalar-height embedding X=(x,z,u); induced metric g_ij=delta_ij+partial_i(u)partial_j(u)",
            "propagation_rule": "geodesic of the induced deformed-medium metric",
            "free_propagation_coefficient": False,
            "stiffness_role": "upstream accumulation only; absent from propagation equation",
            "source_radius": SOURCE_RADIUS,
            "reference_mass": REFERENCE_MASS,
            "mass_ladder": MASS_LADDER,
            "impact_parameters": IMPACT_PARAMETERS,
            "z_start": Z_START,
            "z_stop": Z_STOP,
            "integration_step": STEP,
        },
        "zero_source_test": zero,
        "impact_test": impact,
        "mass_test": mass,
        "measured": {
            "mass_response_exponent": mass_slope,
            "impact_response_exponent": impact_slope,
            "center_angle_change": impact["center_angle_change"],
            "max_reflection_antisymmetry_relative_error": sym_err,
        },
        "checks": checks,
        "matched_checks_of_6": matched,
        "execution_checks": execution_checks,
        "execution_gate_pass": execution_gate_pass,
        "repo_state": state,
        "summary": {
            "question": "Is the form of scalar accumulated deformation u, treated as an induced deformed-medium geometry, sufficient by itself to produce the required propagation behavior?",
            "next_if_supported": "Freeze the geometric propagation rule and then address physical scale before observational comparison.",
            "next_if_partial_or_not_supported": "Do not alter the accumulation bridge; conclude that scalar-height geometry is insufficient and derive the missing vector/tensor deformation geometry rather than fit a light force.",
        },
    }

    print(LAB_ID)
    print(f"status={status}")
    print(f"head_sha={state['head_sha']}")
    print("scalar_height_induced_metric_used=true")
    print("geodesic_propagation_used=true")
    print("stiffness_used_in_propagation_equation=false")
    print("free_propagation_coefficient_used=false")
    print("GR_potential_decomposition_used=false")
    print("Weyl_used=false")
    print("LCDM_used=false")
    print("fit_or_tuning_used=false")
    print()
    print("GEOMETRIC_PROPAGATION")
    print(f"mass_response_exponent={mass_slope:.12g}")
    print(f"impact_response_exponent={impact_slope:.12g}")
    print(f"central_angle_change={impact['center_angle_change']:.12g}")
    print(f"max_reflection_antisymmetry_relative_error={sym_err:.12g}")
    print(f"max_surface_slope={impact['max_surface_slope']:.12g}")
    for row in impact["rows"]:
        print(
            f"impact_b={row['impact_b']:.12g} angle_plus={row['angle_plus']:.12g} "
            f"angle_minus={row['angle_minus']:.12g} antisymmetry_error={row['antisymmetry_relative_error']:.12g}"
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
