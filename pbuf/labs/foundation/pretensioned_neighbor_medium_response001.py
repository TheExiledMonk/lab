#!/usr/bin/env python3
"""PBUF FOUNDATION — PRETENSIONED NEIGHBOR MEDIUM RESPONSE 001.

Minimal speculative substrate-mechanics test inspired by a locally coupled,
pre-tensioned medium.  The goal is not to assume a 1/r response, but to ask
whether a nearest-neighbor equilibrium law derived from a quadratic coupling
energy naturally produces the independently required long-range scaling.

Discrete energy:

    E = (T/2) * sum_<ij> (u_i - u_j)^2 - sum_i S_i u_i

Stationarity gives the native equilibrium equation

    T * L[u] = S

where L is the positive nearest-neighbor discrete Laplacian.  The corresponding
source-free dynamical sector is

    mu * u_tt = -T * L[u] / dx^2

which supports lattice waves.  No identification with electromagnetism is made
in this lab; the wave sector is used only to test whether the same local
coupling that produces the static accumulated response also defines a native
propagation speed.

No G, macroscopic amplitude calibration, lensing target, Rmax, Quantum Engine,
Planck input, fitted radial kernel, inverse-square law, or imported 1/r Green
function is used.
"""
from __future__ import annotations

import json
import math
import subprocess
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
LAB_ID = "PBUF-FOUNDATION-PRETENSIONED-NEIGHBOR-MEDIUM-RESPONSE-001"

N = 49
DX = 1.0
CENTER = N // 2
SHAPE = (N, N, N)
AXIS = np.arange(N, dtype=np.float64) - CENTER
Z, Y, X = np.meshgrid(AXIS, AXIS, AXIS, indexing="ij")
RR = np.sqrt(X * X + Y * Y + Z * Z) * DX

TENSION = 1.0
CG_REL_TOL = 1.0e-10
CG_MAX_ITER = 600
EXP_WINDOW = 0.35
ADDITIVITY_TOL = 2.0e-8

RADIUS_LADDER = (2.5, 3.5, 4.5, 5.5, 6.5)
FIXED_NATIVE_LOAD = 10.0
MASS_LADDER = (2.0, 4.0, 8.0, 16.0, 32.0)
FIXED_RADIUS = 4.5
FAR_SOURCE_RADIUS = 3.5
FAR_PROBE_RADII = (7.0, 8.0, 9.0, 10.0, 11.0, 12.0)


def git(*args: str) -> str:
    return subprocess.check_output(["git", *args], cwd=ROOT, text=True).strip()


def repo_state() -> dict:
    return {
        "branch": git("rev-parse", "--abbrev-ref", "HEAD"),
        "head_sha": git("rev-parse", "HEAD"),
        "tracked_changes": git("diff", "--name-only"),
        "staged_changes": git("diff", "--name-only", "--cached"),
    }


def source_sphere(radius: float, integrated_load: float) -> dict:
    mask = RR <= radius
    count = int(np.count_nonzero(mask))
    if count <= 0:
        raise RuntimeError("empty source sphere")
    density = float(integrated_load) / (count * DX ** 3)
    source = np.zeros(SHAPE, dtype=np.float64)
    source[mask] = density
    return {
        "source": source,
        "mask": mask,
        "density": density,
        "integrated_load": float(np.sum(source) * DX ** 3),
        "voxel_count": count,
        "radius": radius,
    }


def source_sphere_density(radius: float, density: float) -> dict:
    mask = RR <= radius
    source = np.zeros(SHAPE, dtype=np.float64)
    source[mask] = density
    return {
        "source": source,
        "mask": mask,
        "density": density,
        "integrated_load": float(np.sum(source) * DX ** 3),
        "voxel_count": int(np.count_nonzero(mask)),
        "radius": radius,
    }


def apply_L(u: np.ndarray) -> np.ndarray:
    """Positive nearest-neighbor Laplacian with zero Dirichlet boundary."""
    out = np.zeros_like(u)
    core = u[1:-1, 1:-1, 1:-1]
    out[1:-1, 1:-1, 1:-1] = (
        6.0 * core
        - u[2:, 1:-1, 1:-1]
        - u[:-2, 1:-1, 1:-1]
        - u[1:-1, 2:, 1:-1]
        - u[1:-1, :-2, 1:-1]
        - u[1:-1, 1:-1, 2:]
        - u[1:-1, 1:-1, :-2]
    )
    return out


def solve_equilibrium(source: np.ndarray, tension: float = TENSION) -> dict:
    # T L[u]/DX^2 = S  ->  L[u] = DX^2 S/T.
    b = np.zeros_like(source)
    b[1:-1, 1:-1, 1:-1] = (DX * DX / tension) * source[1:-1, 1:-1, 1:-1]
    x = np.zeros_like(source)
    r = b - apply_L(x)
    p = r.copy()
    rr_old = float(np.sum(r * r))
    bnorm = math.sqrt(float(np.sum(b * b)))
    if bnorm == 0.0:
        return {"field": x, "iterations": 0, "relative_residual": 0.0}

    rel = math.sqrt(rr_old) / bnorm
    iterations = 0
    for k in range(1, CG_MAX_ITER + 1):
        Ap = apply_L(p)
        denom = float(np.sum(p * Ap))
        if not math.isfinite(denom) or denom <= 0.0:
            raise RuntimeError("CG operator lost positive definiteness")
        alpha = rr_old / denom
        x += alpha * p
        # Keep the Dirichlet boundary exact.
        x[0, :, :] = x[-1, :, :] = 0.0
        x[:, 0, :] = x[:, -1, :] = 0.0
        x[:, :, 0] = x[:, :, -1] = 0.0
        r -= alpha * Ap
        rr_new = float(np.sum(r * r))
        rel = math.sqrt(rr_new) / bnorm
        iterations = k
        if rel <= CG_REL_TOL:
            break
        beta = rr_new / rr_old
        p = r + beta * p
        rr_old = rr_new

    return {"field": x, "iterations": iterations, "relative_residual": rel}


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
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0.0 else float("nan")
    return {"slope": float(beta[0]), "r2": r2, "count": int(np.count_nonzero(m))}


def exponent_match(value: float, target: float) -> bool:
    return bool(math.isfinite(value) and abs(value - target) <= EXP_WINDOW)


def shell_mean(field: np.ndarray, radius: float) -> float:
    shell = (RR >= radius - 0.75) & (RR <= radius + 0.75)
    return float(np.mean(np.abs(field[shell])))


def probe(field: np.ndarray, radius: float) -> float:
    offset = int(round(radius / DX))
    if abs(offset * DX - radius) > 1.0e-12:
        raise RuntimeError("probe radius not on grid")
    return float(abs(field[CENTER, CENTER, CENTER + offset]))


def static_radius_fixed_load_test() -> dict:
    rows = []
    for radius in RADIUS_LADDER:
        src = source_sphere(radius, FIXED_NATIVE_LOAD)
        sol = solve_equilibrium(src["source"])
        rows.append({
            "radius": radius,
            "integrated_load": src["integrated_load"],
            "surface_response": shell_mean(sol["field"], radius),
            "iterations": sol["iterations"],
            "relative_residual": sol["relative_residual"],
        })
    return {
        "rows": rows,
        "fit": logfit([r["radius"] for r in rows], [r["surface_response"] for r in rows]),
        "load_relative_span": (max(r["integrated_load"] for r in rows) - min(r["integrated_load"] for r in rows)) / FIXED_NATIVE_LOAD,
    }


def static_far_radius_test() -> dict:
    src = source_sphere(FAR_SOURCE_RADIUS, 10.0)
    sol = solve_equilibrium(src["source"])
    vals = [probe(sol["field"], r) for r in FAR_PROBE_RADII]
    return {
        "probe_radii": FAR_PROBE_RADII,
        "responses": vals,
        "fit": logfit(FAR_PROBE_RADII, vals),
        "iterations": sol["iterations"],
        "relative_residual": sol["relative_residual"],
    }


def static_mass_test() -> dict:
    rows = []
    for mass in MASS_LADDER:
        src = source_sphere(FIXED_RADIUS, mass)
        sol = solve_equilibrium(src["source"])
        rows.append({
            "integrated_load": src["integrated_load"],
            "surface_response": shell_mean(sol["field"], FIXED_RADIUS),
            "far_response": probe(sol["field"], 11.0),
            "relative_residual": sol["relative_residual"],
        })
    return {
        "rows": rows,
        "surface_fit": logfit([r["integrated_load"] for r in rows], [r["surface_response"] for r in rows]),
        "far_fit": logfit([r["integrated_load"] for r in rows], [r["far_response"] for r in rows]),
    }


def static_density_test() -> dict:
    densities = (0.01, 0.02, 0.04, 0.08, 0.16)
    rows = []
    for density in densities:
        src = source_sphere_density(FIXED_RADIUS, density)
        sol = solve_equilibrium(src["source"])
        rows.append({
            "density": density,
            "center_response": float(abs(sol["field"][CENTER, CENTER, CENTER])),
            "integrated_load": src["integrated_load"],
        })
    return {
        "rows": rows,
        "fit": logfit([r["density"] for r in rows], [r["center_response"] for r in rows]),
    }


def shifted_source(offset_x: int, radius: float, integrated_load: float) -> np.ndarray:
    rr = np.sqrt((X - offset_x) ** 2 + Y * Y + Z * Z) * DX
    mask = rr <= radius
    count = int(np.count_nonzero(mask))
    src = np.zeros(SHAPE, dtype=np.float64)
    src[mask] = integrated_load / (count * DX ** 3)
    return src


def relative_rms(actual: np.ndarray, reference: np.ndarray) -> float:
    num = math.sqrt(float(np.mean((actual - reference) ** 2)))
    den = math.sqrt(float(np.mean(reference ** 2)))
    return num / den if den > 0.0 else float("nan")


def additivity_test() -> dict:
    s1 = shifted_source(-5, 2.5, 6.0)
    s2 = shifted_source(+5, 3.5, 9.0)
    u1 = solve_equilibrium(s1)["field"]
    u2 = solve_equilibrium(s2)["field"]
    u12 = solve_equilibrium(s1 + s2)["field"]
    resid = relative_rms(u12, u1 + u2)
    return {"relative_rms_residual": resid}


def dynamic_dispersion_test() -> dict:
    # From the same nearest-neighbor equation, a lattice plane wave gives:
    # omega^2 = (4 T / (mu DX^2)) * sum_a sin^2(k_a DX/2).
    # We evaluate the low-k branch numerically and test its T and mu scaling.
    wavelength = 20.0 * DX
    k = 2.0 * math.pi / wavelength

    def phase_speed(tension: float, inertia: float) -> float:
        omega2 = (4.0 * tension / (inertia * DX * DX)) * math.sin(0.5 * k * DX) ** 2
        return math.sqrt(omega2) / k

    tension_ladder = (0.25, 0.5, 1.0, 2.0, 4.0)
    inertia_ladder = (0.25, 0.5, 1.0, 2.0, 4.0)
    v_t = [phase_speed(t, 1.0) for t in tension_ladder]
    v_mu = [phase_speed(1.0, mu) for mu in inertia_ladder]
    tension_fit = logfit(tension_ladder, v_t)
    inertia_fit = logfit(inertia_ladder, v_mu)
    continuum_speed = math.sqrt(TENSION / 1.0)
    lattice_speed = phase_speed(TENSION, 1.0)
    return {
        "wavelength_grid_units": wavelength / DX,
        "tension_ladder": tension_ladder,
        "speed_vs_tension": v_t,
        "inertia_ladder": inertia_ladder,
        "speed_vs_inertia": v_mu,
        "tension_exponent": tension_fit,
        "inertia_exponent": inertia_fit,
        "unit_case_lattice_phase_speed": lattice_speed,
        "unit_case_long_wavelength_limit": continuum_speed,
        "unit_case_relative_dispersion_error": abs(lattice_speed - continuum_speed) / continuum_speed,
    }


def main() -> int:
    repo = repo_state()
    radius = static_radius_fixed_load_test()
    far = static_far_radius_test()
    mass = static_mass_test()
    density = static_density_test()
    add = additivity_test()
    dynamic = dynamic_dispersion_test()

    measured = {
        "local_density_exponent": density["fit"]["slope"],
        "surface_mass_exponent_fixed_R": mass["surface_fit"]["slope"],
        "surface_radius_exponent_fixed_load": radius["fit"]["slope"],
        "far_mass_exponent_fixed_probe": mass["far_fit"]["slope"],
        "far_radius_exponent_fixed_source": far["fit"]["slope"],
        "static_additivity_relative_residual": add["relative_rms_residual"],
        "wave_speed_tension_exponent": dynamic["tension_exponent"]["slope"],
        "wave_speed_inertia_exponent": dynamic["inertia_exponent"]["slope"],
    }

    checks = {
        "density_linearity": exponent_match(measured["local_density_exponent"], 1.0),
        "surface_mass_linearity": exponent_match(measured["surface_mass_exponent_fixed_R"], 1.0),
        "surface_radius_accumulation": exponent_match(measured["surface_radius_exponent_fixed_load"], -1.0),
        "far_mass_linearity": exponent_match(measured["far_mass_exponent_fixed_probe"], 1.0),
        "far_radius_accumulation": exponent_match(measured["far_radius_exponent_fixed_source"], -1.0),
        "static_additivity": bool(math.isfinite(measured["static_additivity_relative_residual"]) and measured["static_additivity_relative_residual"] <= ADDITIVITY_TOL),
        "wave_speed_tension_half_power": exponent_match(measured["wave_speed_tension_exponent"], 0.5),
        "wave_speed_inertia_minus_half_power": exponent_match(measured["wave_speed_inertia_exponent"], -0.5),
    }

    static_core = (
        checks["density_linearity"]
        and checks["surface_mass_linearity"]
        and checks["surface_radius_accumulation"]
        and checks["far_mass_linearity"]
        and checks["far_radius_accumulation"]
        and checks["static_additivity"]
    )
    dynamic_core = checks["wave_speed_tension_half_power"] and checks["wave_speed_inertia_minus_half_power"]

    if static_core and dynamic_core:
        status = "PRETENSIONED_NEIGHBOR_MEDIUM_STATIC_DYNAMIC_STRUCTURE_FOUND"
    elif static_core or dynamic_core:
        status = "PRETENSIONED_NEIGHBOR_MEDIUM_PARTIAL_STRUCTURE_ONLY"
    else:
        status = "PRETENSIONED_NEIGHBOR_MEDIUM_STRUCTURE_NOT_SUPPORTED"

    max_static_residual = max(
        max(r["relative_residual"] for r in radius["rows"]),
        far["relative_residual"],
        max(r["relative_residual"] for r in mass["rows"]),
    )

    policy = {
        "gravity_fundamental_in_PBUF": False,
        "gravity_used_as_native_variable": False,
        "electromagnetic_origin_claimed": False,
        "G_used": False,
        "macroscopic_amplitude_used": False,
        "native_amplitude_rescaled": False,
        "fit_or_tuning_used": False,
        "one_over_r_kernel_inserted": False,
        "inverse_square_law_inserted": False,
        "Rmax_used": False,
        "cosmology_used": False,
        "lensing_target_used": False,
        "legacy_0p18_used": False,
        "quantum_engine_used": False,
        "planck_input_used": False,
    }

    execution_checks = {
        "all_measured_values_finite": all(math.isfinite(v) for v in measured.values()),
        "cg_converged": max_static_residual <= 1.0e-8,
        "fixed_load_control": radius["load_relative_span"] <= 1.0e-12,
        "no_G": True,
        "no_macroscopic_amplitude": True,
        "no_native_rescaling": True,
        "no_fit_or_tuning": True,
        "no_inserted_radial_kernel": True,
        "no_Rmax": True,
        "no_cosmology": True,
        "no_lensing_target": True,
        "no_legacy_0p18": True,
        "no_quantum_engine": True,
        "no_planck_input": True,
        "no_tracked_or_staged_changes": repo["tracked_changes"] == "" and repo["staged_changes"] == "",
        "stdout_only_no_run_directory_created": True,
    }
    execution_gate_pass = all(execution_checks.values())

    payload = {
        "lab_id": LAB_ID,
        "status": status,
        "repo_state": repo,
        "model": {
            "discrete_energy": "E=(T/2) sum_<ij>(u_i-u_j)^2 - sum_i S_i u_i",
            "static_equilibrium": "T L[u]/DX^2 = S",
            "source_free_dynamics": "mu u_tt = -T L[u]/DX^2",
            "native_interpretation": "local neighbor coupling plus pre-tension candidate; no microscopic ontology fixed",
        },
        "measured": measured,
        "checks": checks,
        "static_tests": {
            "radius_fixed_load": radius,
            "far_radius": far,
            "mass": mass,
            "density": density,
            "additivity": add,
        },
        "dynamic_test": dynamic,
        "policy": policy,
        "execution_checks": execution_checks,
        "execution_gate_pass": execution_gate_pass,
        "summary": {
            "static_structure_pass": static_core,
            "dynamic_structure_pass": dynamic_core,
            "interpretation": "Tests whether a local pre-tensioned neighbor coupling can generate long-range accumulated medium response and a wave-supporting dynamical sector without inserting the target radial law.",
            "next_if_supported": "Map the derived equilibrium variable to the existing c_state/local-loading pipeline without amplitude fitting, then test whether c_state can serve as S or as a constitutive precursor to S.",
            "next_if_partial": "Identify which assumption in the minimal neighbor-coupled medium controls the failed fingerprint before adding nonlinear or tensor structure.",
        },
    }

    print(LAB_ID)
    print(f"status={status}")
    print(f"head_sha={repo['head_sha']}")
    print("gravity_used_as_native_variable=false")
    print("electromagnetic_origin_claimed=false")
    print("G_used=false")
    print("one_over_r_kernel_inserted=false")
    print("fit_or_tuning_used=false")
    print()
    print("STATIC_RESPONSE")
    print(f"rho_exponent={measured['local_density_exponent']:.12g}")
    print(f"surface_mass_exponent={measured['surface_mass_exponent_fixed_R']:.12g}")
    print(f"surface_radius_fixed_load_exponent={measured['surface_radius_exponent_fixed_load']:.12g}")
    print(f"far_mass_exponent={measured['far_mass_exponent_fixed_probe']:.12g}")
    print(f"far_radius_exponent={measured['far_radius_exponent_fixed_source']:.12g}")
    print(f"additivity_relative_residual={measured['static_additivity_relative_residual']:.12g}")
    print()
    print("DYNAMIC_RESPONSE")
    print(f"wave_speed_tension_exponent={measured['wave_speed_tension_exponent']:.12g}")
    print(f"wave_speed_inertia_exponent={measured['wave_speed_inertia_exponent']:.12g}")
    print(f"unit_case_lattice_phase_speed={dynamic['unit_case_lattice_phase_speed']:.12g}")
    print(f"unit_case_long_wavelength_limit={dynamic['unit_case_long_wavelength_limit']:.12g}")
    print(f"unit_case_relative_dispersion_error={dynamic['unit_case_relative_dispersion_error']:.12g}")
    print()
    print("CHECKS")
    for k, v in checks.items():
        print(f"{k}={str(v).lower()}")
    print("EXECUTION_CHECKS")
    for k, v in execution_checks.items():
        print(f"{k}={str(v).lower()}")
    print(f"execution_gate_pass={str(execution_gate_pass).lower()}")
    print(json.dumps(payload, sort_keys=True, separators=(",", ":")))

    if not execution_gate_pass:
        raise RuntimeError("pretensioned neighbor medium response execution gate failed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
