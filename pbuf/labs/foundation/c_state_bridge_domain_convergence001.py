#!/usr/bin/env python3
"""PBUF FOUNDATION — C_STATE BRIDGE DOMAIN CONVERGENCE 001.

Tests whether the native c_state -> bounded-strain 3D accumulation bridge keeps
its source physics fixed while the numerical accumulation boundary is moved
outward.

Frozen native chain:
    rho -> existing A8 transport -> raw c_state

Frozen accumulation chain:
    raw c_state -> six-neighbor bounded-strain equilibrium -> accumulated field

Only the accumulation grid size changes.  The native 33^3 c_state generation,
raw amplitude, source geometry, source mass, constitutive law, K0, epsilon_max,
probe radii, nonlinear tolerances, and damping are fixed.

No G, macroscopic amplitude calibration, fitted K, inserted 1/r law, spherical
shortcut, Rmax, cosmology, lensing target, Quantum Engine, or Planck input.
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

from pbuf.models import a8_state as M06_state

LAB_ID = "PBUF-FOUNDATION-C-STATE-BRIDGE-DOMAIN-CONVERGENCE-001"

N_NATIVE = 33
DX = 1.0
CENTER_NATIVE = N_NATIVE // 2
AXN = np.arange(N_NATIVE, dtype=np.float64) - CENTER_NATIVE
ZN, YN, XN = np.meshgrid(AXN, AXN, AXN, indexing="ij")
RR_NATIVE = np.sqrt(XN * XN + YN * YN + ZN * ZN) * DX

K0 = 1.0
EPSILON_MAX = 1.0
PICARD_TOL = 2.0e-7
PICARD_MAX_ITER = 30
PICARD_DAMP = 0.65
CG_REL_TOL = 2.0e-9
CG_MAX_ITER = 900

GRID_LADDER = (49, 65, 81, 97)
SOURCE_RADIUS = 3.5
SOURCE_MASS = 2.0
PROBE_RADII = (6.0, 7.0, 8.0, 9.0, 10.0)
LOCAL_PROBE = 8.0
FIT_NEAR_MINUS1_TOL = 0.20
STRAIN_STABILITY_REL_TOL = 5.0e-4
CSTATE_INTEGRAL_REL_TOL = 1.0e-12


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


def native_rho_fixed_mass() -> np.ndarray:
    mask = RR_NATIVE <= SOURCE_RADIUS
    count = int(np.count_nonzero(mask))
    rho = np.zeros((N_NATIVE, N_NATIVE, N_NATIVE), dtype=np.float64)
    rho[mask] = SOURCE_MASS / (count * DX**3)
    return rho


def native_c_state(rho: np.ndarray) -> np.ndarray:
    u_slow0 = np.asarray(rho, dtype=np.float64).copy()
    u_fast0 = np.asarray(rho, dtype=np.float64).copy()
    us, uf, history = M06_state.evolve_a8_transport_3d(
        u_slow0, u_fast0, stencil="N6", boundary="reflective"
    )
    max_abs = max(float(np.max(np.abs(us))), float(np.max(np.abs(uf))))
    if max_abs >= M06_state.A8_INIT_CLIP - 1.0e-12:
        raise RuntimeError("native c_state clipping gate failed")
    return np.asarray(history[-1], dtype=np.float64)


def geometry(n: int):
    c = n // 2
    ax = np.arange(n, dtype=np.float64) - c
    z, y, x = np.meshgrid(ax, ax, ax, indexing="ij")
    rr = np.sqrt(x * x + y * y + z * z) * DX
    return c, rr


def embed_raw(native: np.ndarray, n: int) -> np.ndarray:
    out = np.zeros((n, n, n), dtype=np.float64)
    c = n // 2
    lo = c - CENTER_NATIVE
    hi = lo + N_NATIVE
    if lo <= 0 or hi >= n:
        raise RuntimeError("accumulation grid too small for native embedding")
    out[lo:hi, lo:hi, lo:hi] = np.asarray(native, dtype=np.float64)
    return out


def zero_boundary(a: np.ndarray) -> None:
    a[0, :, :] = a[-1, :, :] = 0.0
    a[:, 0, :] = a[:, -1, :] = 0.0
    a[:, :, 0] = a[:, :, -1] = 0.0


def bond_strains(u: np.ndarray):
    return (
        (u[1:, :, :] - u[:-1, :, :]) / DX,
        (u[:, 1:, :] - u[:, :-1, :]) / DX,
        (u[:, :, 1:] - u[:, :, :-1]) / DX,
    )


def secant_weights(u: np.ndarray):
    weights = []
    max_frac = 0.0
    for e in bond_strains(u):
        frac = np.abs(e) / EPSILON_MAX
        max_frac = max(max_frac, float(np.max(frac)))
        if max_frac >= 0.995:
            raise RuntimeError("bounded-strain barrier approached too closely")
        weights.append(K0 / (1.0 - frac * frac))
    return tuple(weights), max_frac


def apply_A(u: np.ndarray, weights) -> np.ndarray:
    wz, wy, wx = weights
    out = np.zeros_like(u)
    dz = (u[1:, :, :] - u[:-1, :, :]) * wz
    dy = (u[:, 1:, :] - u[:, :-1, :]) * wy
    dx = (u[:, :, 1:] - u[:, :, :-1]) * wx
    out[:-1, :, :] -= dz
    out[1:, :, :] += dz
    out[:, :-1, :] -= dy
    out[:, 1:, :] += dy
    out[:, :, :-1] -= dx
    out[:, :, 1:] += dx
    zero_boundary(out)
    return out / (DX * DX)


def cg_solve(source: np.ndarray, weights, x0=None) -> dict:
    b = source.copy()
    zero_boundary(b)
    x = np.zeros_like(b) if x0 is None else x0.copy()
    zero_boundary(x)
    r = b - apply_A(x, weights)
    zero_boundary(r)
    p = r.copy()
    rr_old = float(np.sum(r * r))
    bnorm = math.sqrt(float(np.sum(b * b)))
    if bnorm == 0.0:
        return {"field": x, "iterations": 0, "relative_residual": 0.0}
    rel = math.sqrt(rr_old) / bnorm
    it = 0
    for k in range(1, CG_MAX_ITER + 1):
        Ap = apply_A(p, weights)
        denom = float(np.sum(p * Ap))
        if not math.isfinite(denom) or denom <= 0.0:
            raise RuntimeError("CG operator lost positive definiteness")
        alpha = rr_old / denom
        x += alpha * p
        zero_boundary(x)
        r -= alpha * Ap
        zero_boundary(r)
        rr_new = float(np.sum(r * r))
        rel = math.sqrt(rr_new) / bnorm
        it = k
        if rel <= CG_REL_TOL:
            break
        beta = rr_new / rr_old
        p = r + beta * p
        zero_boundary(p)
        rr_old = rr_new
    return {"field": x, "iterations": it, "relative_residual": rel}


def solve_nonlinear(source: np.ndarray) -> dict:
    n = source.shape[0]
    ones = (
        np.ones((n - 1, n, n), dtype=np.float64) * K0,
        np.ones((n, n - 1, n), dtype=np.float64) * K0,
        np.ones((n, n, n - 1), dtype=np.float64) * K0,
    )
    first = cg_solve(source, ones)
    u = first["field"]
    history = []
    converged = False
    for it in range(1, PICARD_MAX_ITER + 1):
        weights, _ = secant_weights(u)
        sol = cg_solve(source, weights, x0=u)
        candidate = sol["field"]
        new_u = PICARD_DAMP * candidate + (1.0 - PICARD_DAMP) * u
        zero_boundary(new_u)
        scale = max(float(np.sqrt(np.mean(new_u * new_u))), 1.0e-14)
        change = float(np.sqrt(np.mean((new_u - u) ** 2))) / scale
        u = new_u
        _, max_frac = secant_weights(u)
        history.append({
            "picard": it,
            "relative_change": change,
            "max_strain_fraction": max_frac,
            "cg_iterations": sol["iterations"],
            "cg_relative_residual": sol["relative_residual"],
        })
        if change <= PICARD_TOL:
            converged = True
            break
    weights, max_frac = secant_weights(u)
    residual = source - apply_A(u, weights)
    zero_boundary(residual)
    src_norm = math.sqrt(float(np.sum(source[1:-1, 1:-1, 1:-1] ** 2)))
    nl_rel = math.sqrt(float(np.sum(residual * residual))) / src_norm if src_norm > 0 else 0.0
    return {
        "field": u,
        "converged": converged,
        "picard_iterations": len(history),
        "max_strain_fraction": max_frac,
        "nonlinear_relative_residual": nl_rel,
        "history_tail": history[-3:],
    }


def probe(field: np.ndarray, center: int, radius: float) -> float:
    off = int(round(radius / DX))
    return float(abs(field[center, center, center + off]))


def shell_mean(field: np.ndarray, rr: np.ndarray, radius: float) -> float:
    mask = (rr >= radius - 0.75) & (rr <= radius + 0.75)
    return float(np.mean(np.abs(field[mask])))


def logfit(xs, ys) -> dict:
    x = np.asarray(xs, dtype=np.float64)
    y = np.asarray(ys, dtype=np.float64)
    m = np.isfinite(x) & np.isfinite(y) & (x > 0) & (y > 0)
    lx = np.log(x[m])
    ly = np.log(y[m])
    A = np.column_stack((lx, np.ones_like(lx)))
    beta, *_ = np.linalg.lstsq(A, ly, rcond=None)
    pred = A @ beta
    ss_res = float(np.sum((ly - pred) ** 2))
    ss_tot = float(np.sum((ly - np.mean(ly)) ** 2))
    return {
        "slope": float(beta[0]),
        "r2": 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan"),
        "count": int(np.count_nonzero(m)),
    }


def local_exponent(radii, responses, target_radius: float) -> float:
    r = np.asarray(radii, dtype=np.float64)
    y = np.asarray(responses, dtype=np.float64)
    idx = int(np.argmin(np.abs(r - target_radius)))
    if idx == 0 or idx == len(r) - 1:
        raise RuntimeError("local exponent target must have neighbors")
    return float((math.log(y[idx + 1]) - math.log(y[idx - 1])) /
                 (math.log(r[idx + 1]) - math.log(r[idx - 1])))


def monotonic_toward_minus1(values) -> bool:
    errs = [abs(float(v) + 1.0) for v in values]
    return all(errs[i + 1] < errs[i] for i in range(len(errs) - 1))


def main() -> int:
    state = repo_state()
    rho = native_rho_fixed_mass()
    c_state = native_c_state(rho)
    rho_integral = float(np.sum(rho) * DX**3)
    c_state_integral = float(np.sum(c_state) * DX**3)

    rows = []
    for n in GRID_LADDER:
        center, rr = geometry(n)
        source = embed_raw(c_state, n)
        sol = solve_nonlinear(source)
        responses = [probe(sol["field"], center, r) for r in PROBE_RADII]
        fit = logfit(PROBE_RADII, responses)
        local_p = local_exponent(PROBE_RADII, responses, LOCAL_PROBE)
        rows.append({
            "grid_N": n,
            "half_width": float(center * DX),
            "far_fit": fit,
            "local_exponent_r8": local_p,
            "surface_response": shell_mean(sol["field"], rr, SOURCE_RADIUS),
            "responses": responses,
            "probe_radii": PROBE_RADII,
            "max_strain_fraction": sol["max_strain_fraction"],
            "converged": sol["converged"],
            "picard_iterations": sol["picard_iterations"],
            "nonlinear_relative_residual": sol["nonlinear_relative_residual"],
            "history_tail": sol["history_tail"],
        })

    far_slopes = [r["far_fit"]["slope"] for r in rows]
    local_ps = [r["local_exponent_r8"] for r in rows]
    strains = np.asarray([r["max_strain_fraction"] for r in rows], dtype=np.float64)
    surface = np.asarray([r["surface_response"] for r in rows], dtype=np.float64)

    strain_rel_span = float((np.max(strains) - np.min(strains)) / max(np.mean(strains), 1.0e-30))
    surface_rel_span = float((np.max(surface) - np.min(surface)) / max(np.mean(surface), 1.0e-30))
    cstate_integral_rel_error = abs(c_state_integral - rho_integral) / max(abs(rho_integral), 1.0e-30)

    checks = {
        "far_fit_moves_toward_minus1_with_domain": monotonic_toward_minus1(far_slopes),
        "fixed_r8_local_exponent_moves_toward_minus1": monotonic_toward_minus1(local_ps),
        "largest_domain_far_fit_near_minus1": abs(far_slopes[-1] + 1.0) <= FIT_NEAR_MINUS1_TOL,
        "largest_domain_local_exponent_near_minus1": abs(local_ps[-1] + 1.0) <= FIT_NEAR_MINUS1_TOL,
        "source_strain_stable_across_domain": strain_rel_span <= STRAIN_STABILITY_REL_TOL,
        "native_c_state_integral_preserved": cstate_integral_rel_error <= CSTATE_INTEGRAL_REL_TOL,
    }

    execution_checks = {
        "all_measured_values_finite": bool(all(np.isfinite(far_slopes)) and all(np.isfinite(local_ps)) and np.all(np.isfinite(strains))),
        "nonlinear_solves_converged": bool(all(r["converged"] for r in rows)),
        "raw_c_state_used_without_amplitude_rescaling": True,
        "native_grid_frozen": True,
        "source_geometry_frozen": True,
        "source_mass_frozen": True,
        "probe_radii_frozen": True,
        "constitutive_law_frozen": True,
        "K0_frozen": True,
        "epsilon_max_frozen": True,
        "no_G": True,
        "no_macroscopic_amplitude": True,
        "no_native_rescaling": True,
        "no_fit_or_tuning": True,
        "no_inserted_one_over_r_response": True,
        "no_spherical_equilibrium_shortcut": True,
        "no_Rmax": True,
        "no_cosmology": True,
        "no_lensing_target": True,
        "no_quantum_engine": True,
        "no_planck_input": True,
        "no_tracked_or_staged_changes": state["tracked_changes"] == "" and state["staged_changes"] == "",
        "stdout_only_no_run_directory_created": True,
    }

    execution_gate_pass = bool(all(checks.values()) and all(execution_checks.values()))
    if execution_gate_pass:
        status = "C_STATE_BRIDGE_DOMAIN_CONVERGENCE_SUPPORTED"
    elif (checks["far_fit_moves_toward_minus1_with_domain"] or checks["fixed_r8_local_exponent_moves_toward_minus1"]):
        status = "C_STATE_BRIDGE_DOMAIN_CONVERGENCE_PARTIAL_SUPPORT"
    else:
        status = "C_STATE_BRIDGE_DOMAIN_CONVERGENCE_NOT_SUPPORTED"

    report = {
        "lab_id": LAB_ID,
        "status": status,
        "repo_state": state,
        "model": {
            "native_grid_N": N_NATIVE,
            "accumulation_grid_ladder": GRID_LADDER,
            "source_radius": SOURCE_RADIUS,
            "source_mass": SOURCE_MASS,
            "probe_radii": PROBE_RADII,
            "K0": K0,
            "epsilon_max": EPSILON_MAX,
            "native_source": "rho -> existing A8 transport -> raw c_state",
            "accumulation": "raw c_state -> six-neighbor bounded-strain equilibrium",
            "amplitude_rescaling": False,
        },
        "native_source_metrics": {
            "rho_integral": rho_integral,
            "c_state_integral": c_state_integral,
            "c_state_integral_relative_error": cstate_integral_rel_error,
        },
        "rows": rows,
        "checks": checks,
        "execution_checks": execution_checks,
        "execution_gate_pass": execution_gate_pass,
        "summary": {
            "far_slopes": far_slopes,
            "local_r8_exponents": local_ps,
            "strain_relative_span": strain_rel_span,
            "surface_response_relative_span": surface_rel_span,
            "question": "Does the frozen native c_state -> bounded-strain accumulation bridge converge toward a -1 long-range radial exponent when only the accumulation domain is enlarged?",
        },
    }

    print(LAB_ID)
    print(f"status={status}")
    print(f"head_sha={state['head_sha']}")
    print("raw_c_state_used=true")
    print("native_amplitude_rescaled=false")
    print("K0_changed=false")
    print("fit_or_tuning_used=false")
    print("one_over_r_response_inserted=false")
    print("spherical_equilibrium_shortcut_used=false")
    print()
    print("C_STATE_BRIDGE_DOMAIN_CONVERGENCE")
    for r in rows:
        print(
            f"grid_N={r['grid_N']} half_width={r['half_width']:.12g} "
            f"far_fit={r['far_fit']['slope']:.12g} "
            f"local_p_r8={r['local_exponent_r8']:.12g} "
            f"surface_response={r['surface_response']:.12g} "
            f"max_strain_fraction={r['max_strain_fraction']:.12g} "
            f"converged={str(r['converged']).lower()}"
        )
    print()
    print("NATIVE_SOURCE")
    print(f"rho_integral={rho_integral:.12g}")
    print(f"c_state_integral={c_state_integral:.12g}")
    print(f"c_state_integral_relative_error={cstate_integral_rel_error:.12g}")
    print()
    print("CHECKS")
    for k, v in checks.items():
        print(f"{k}={str(v).lower()}")
    print("EXECUTION_CHECKS")
    for k, v in execution_checks.items():
        print(f"{k}={str(v).lower()}")
    print(f"execution_gate_pass={str(execution_gate_pass).lower()}")
    print(json.dumps(report, sort_keys=True, separators=(",", ":")))
    return 0 if execution_gate_pass else 1


if __name__ == "__main__":
    raise SystemExit(main())
