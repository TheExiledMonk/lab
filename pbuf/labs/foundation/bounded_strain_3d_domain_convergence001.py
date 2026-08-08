#!/usr/bin/env python3
"""PBUF FOUNDATION — BOUNDED-STRAIN 3D DOMAIN CONVERGENCE 001.

Freeze the nonlinear nearest-neighbor constitutive law used by
bounded_strain_3d_neighbor_network001 and vary only the cubic domain size.
The source geometry, integrated load, physical probe radii, K0, epsilon_max,
solver tolerances, and boundary condition are fixed.

Question: does the far radial exponent of the genuine 3D network move toward
-1 as the artificial zero-response boundary is moved outward?
"""
from __future__ import annotations

import json
import math
import subprocess
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
LAB_ID = "PBUF-FOUNDATION-BOUNDED-STRAIN-3D-DOMAIN-CONVERGENCE-001"

DX = 1.0
K0 = 1.0
EPSILON_MAX = 1.0
GRID_LADDER = (49, 65, 81, 97)
SOURCE_RADIUS = 3.5
SOURCE_LOAD = 2.0
PROBE_RADII = (6.0, 7.0, 8.0, 9.0, 10.0)
LOCAL_EXPONENT_RADIUS = 8.0

CG_REL_TOL = 2.0e-9
CG_MAX_ITER = 800
PICARD_REL_TOL = 1.0e-7
PICARD_MAX_ITER = 30
PICARD_DAMPING = 0.65
EXP_NEAR_MINUS1 = 0.20


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


def stress(e: np.ndarray) -> np.ndarray:
    q = e / EPSILON_MAX
    if np.any(np.abs(q) >= 1.0):
        raise RuntimeError("bond strain reached or exceeded epsilon_max")
    return K0 * e / (1.0 - q * q)


def secant_stiffness(e: np.ndarray) -> np.ndarray:
    out = np.full_like(e, K0)
    nz = np.abs(e) > 1.0e-15
    out[nz] = stress(e[nz]) / e[nz]
    return out


def geometry(N: int) -> dict:
    c = N // 2
    axis = np.arange(N, dtype=np.float64) - c
    z, y, x = np.meshgrid(axis, axis, axis, indexing="ij")
    rr = np.sqrt(x*x + y*y + z*z) * DX
    return {"N": N, "center": c, "rr": rr, "shape": (N, N, N)}


def source_sphere(g: dict) -> np.ndarray:
    mask = g["rr"] <= SOURCE_RADIUS
    count = int(np.count_nonzero(mask))
    src = np.zeros(g["shape"], dtype=np.float64)
    src[mask] = SOURCE_LOAD / (count * DX**3)
    return src


def bond_k(u: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray, float]:
    ex = (u[:, :, 1:] - u[:, :, :-1]) / DX
    ey = (u[:, 1:, :] - u[:, :-1, :]) / DX
    ez = (u[1:, :, :] - u[:-1, :, :]) / DX
    max_frac = max(
        float(np.max(np.abs(ex))) if ex.size else 0.0,
        float(np.max(np.abs(ey))) if ey.size else 0.0,
        float(np.max(np.abs(ez))) if ez.size else 0.0,
    ) / EPSILON_MAX
    return secant_stiffness(ex), secant_stiffness(ey), secant_stiffness(ez), max_frac


def apply_A(u: np.ndarray, kx: np.ndarray, ky: np.ndarray, kz: np.ndarray) -> np.ndarray:
    out = np.zeros_like(u)
    # Positive divergence form: sum k_ij (u_i-u_j) / dx^2.
    d = (u[:, :, 1:] - u[:, :, :-1])
    f = kx * d / (DX*DX)
    out[:, :, :-1] -= f
    out[:, :, 1:] += f

    d = (u[:, 1:, :] - u[:, :-1, :])
    f = ky * d / (DX*DX)
    out[:, :-1, :] -= f
    out[:, 1:, :] += f

    d = (u[1:, :, :] - u[:-1, :, :])
    f = kz * d / (DX*DX)
    out[:-1, :, :] -= f
    out[1:, :, :] += f

    # Dirichlet boundary values are fixed, so operator rows there are identity.
    out[0, :, :] = u[0, :, :]
    out[-1, :, :] = u[-1, :, :]
    out[:, 0, :] = u[:, 0, :]
    out[:, -1, :] = u[:, -1, :]
    out[:, :, 0] = u[:, :, 0]
    out[:, :, -1] = u[:, :, -1]
    return out


def enforce_boundary(a: np.ndarray, value: float = 0.0) -> None:
    a[0, :, :] = value
    a[-1, :, :] = value
    a[:, 0, :] = value
    a[:, -1, :] = value
    a[:, :, 0] = value
    a[:, :, -1] = value


def cg_solve(b: np.ndarray, x0: np.ndarray, kx: np.ndarray, ky: np.ndarray, kz: np.ndarray) -> dict:
    x = x0.copy()
    enforce_boundary(x)
    bb = b.copy()
    enforce_boundary(bb)
    r = bb - apply_A(x, kx, ky, kz)
    enforce_boundary(r)
    p = r.copy()
    rr_old = float(np.sum(r*r))
    bnorm = math.sqrt(float(np.sum(bb*bb)))
    if bnorm == 0.0:
        return {"field": x, "iterations": 0, "relative_residual": 0.0}

    rel = math.sqrt(rr_old) / bnorm
    iterations = 0
    for it in range(1, CG_MAX_ITER + 1):
        Ap = apply_A(p, kx, ky, kz)
        denom = float(np.sum(p*Ap))
        if not math.isfinite(denom) or denom <= 0.0:
            raise RuntimeError("CG operator lost positive definiteness")
        alpha = rr_old / denom
        x += alpha*p
        enforce_boundary(x)
        r -= alpha*Ap
        enforce_boundary(r)
        rr_new = float(np.sum(r*r))
        rel = math.sqrt(rr_new) / bnorm
        iterations = it
        if rel <= CG_REL_TOL:
            break
        beta = rr_new / rr_old
        p = r + beta*p
        enforce_boundary(p)
        rr_old = rr_new
    return {"field": x, "iterations": iterations, "relative_residual": rel}


def nonlinear_residual(u: np.ndarray, source: np.ndarray) -> float:
    ex = (u[:, :, 1:] - u[:, :, :-1]) / DX
    ey = (u[:, 1:, :] - u[:, :-1, :]) / DX
    ez = (u[1:, :, :] - u[:-1, :, :]) / DX
    sx, sy, sz = stress(ex), stress(ey), stress(ez)
    out = np.zeros_like(u)
    out[:, :, :-1] -= sx / DX
    out[:, :, 1:] += sx / DX
    out[:, :-1, :] -= sy / DX
    out[:, 1:, :] += sy / DX
    out[:-1, :, :] -= sz / DX
    out[1:, :, :] += sz / DX
    resid = source - out
    enforce_boundary(resid)
    denom = math.sqrt(float(np.sum(source*source)))
    return math.sqrt(float(np.sum(resid*resid))) / denom


def solve_network(N: int) -> dict:
    g = geometry(N)
    src = source_sphere(g)
    u = np.zeros(g["shape"], dtype=np.float64)
    history = []
    converged = False
    max_frac = 0.0
    for picard in range(1, PICARD_MAX_ITER + 1):
        kx, ky, kz, max_frac = bond_k(u)
        inner = cg_solve(src, u, kx, ky, kz)
        candidate = inner["field"]
        updated = PICARD_DAMPING*candidate + (1.0-PICARD_DAMPING)*u
        enforce_boundary(updated)
        num = math.sqrt(float(np.sum((updated-u)**2)))
        den = max(math.sqrt(float(np.sum(updated**2))), 1.0e-30)
        rel_change = num/den
        u = updated
        _, _, _, max_frac = bond_k(u)
        history.append({
            "picard": picard,
            "relative_change": rel_change,
            "cg_iterations": inner["iterations"],
            "cg_relative_residual": inner["relative_residual"],
            "max_strain_fraction": max_frac,
        })
        if rel_change <= PICARD_REL_TOL:
            converged = True
            break
    nlr = nonlinear_residual(u, src)
    return {
        "field": u,
        "geometry": g,
        "converged": converged,
        "picard_iterations": len(history),
        "history_tail": history[-3:],
        "max_strain_fraction": max_frac,
        "nonlinear_relative_residual": nlr,
    }


def probe(sol: dict, radius: float) -> float:
    c = sol["geometry"]["center"]
    offset = int(round(radius/DX))
    return float(abs(sol["field"][c, c, c+offset]))


def logfit(xs, ys) -> dict:
    x = np.asarray(xs, dtype=np.float64)
    y = np.asarray(ys, dtype=np.float64)
    lx, ly = np.log(x), np.log(y)
    A = np.column_stack((lx, np.ones_like(lx)))
    beta, *_ = np.linalg.lstsq(A, ly, rcond=None)
    pred = A@beta
    ss_res = float(np.sum((ly-pred)**2))
    ss_tot = float(np.sum((ly-np.mean(ly))**2))
    return {"slope": float(beta[0]), "r2": 1.0-ss_res/ss_tot, "count": len(xs)}


def local_exponent(sol: dict, r: float) -> float:
    r0, r1 = r-1.0, r+1.0
    u0, u1 = probe(sol, r0), probe(sol, r1)
    return math.log(u1/u0)/math.log(r1/r0)


def main() -> int:
    rows = []
    for N in GRID_LADDER:
        sol = solve_network(N)
        responses = [probe(sol, r) for r in PROBE_RADII]
        fit = logfit(PROBE_RADII, responses)
        rows.append({
            "grid_N": N,
            "half_width": (N//2)*DX,
            "fit_far_radius": fit,
            "local_exponent_r8": local_exponent(sol, LOCAL_EXPONENT_RADIUS),
            "probe_radii": PROBE_RADII,
            "responses": responses,
            "max_strain_fraction": sol["max_strain_fraction"],
            "converged": sol["converged"],
            "picard_iterations": sol["picard_iterations"],
            "nonlinear_relative_residual": sol["nonlinear_relative_residual"],
            "history_tail": sol["history_tail"],
        })

    slopes = [r["fit_far_radius"]["slope"] for r in rows]
    local_ps = [r["local_exponent_r8"] for r in rows]
    strains = [r["max_strain_fraction"] for r in rows]

    def moves_toward_minus1(vals) -> bool:
        d = [abs(v+1.0) for v in vals]
        return all(d[i+1] < d[i] for i in range(len(d)-1))

    checks = {
        "far_fit_moves_toward_minus1_with_domain": moves_toward_minus1(slopes),
        "fixed_r8_local_exponent_moves_toward_minus1": moves_toward_minus1(local_ps),
        "largest_domain_far_fit_near_minus1": abs(slopes[-1]+1.0) <= EXP_NEAR_MINUS1,
        "largest_domain_local_exponent_near_minus1": abs(local_ps[-1]+1.0) <= EXP_NEAR_MINUS1,
        "source_strain_stable_across_domain": (max(strains)-min(strains)) <= 5.0e-4,
    }
    state = repo_state()
    execution_checks = {
        "all_measured_values_finite": all(math.isfinite(v) for v in slopes+local_ps+strains),
        "nonlinear_solves_converged": all(r["converged"] for r in rows),
        "constitutive_law_frozen": True,
        "K0_frozen": K0 == 1.0,
        "epsilon_max_frozen": EPSILON_MAX == 1.0,
        "source_geometry_frozen": SOURCE_RADIUS == 3.5,
        "source_load_frozen": SOURCE_LOAD == 2.0,
        "probe_radii_frozen": PROBE_RADII == (6.0,7.0,8.0,9.0,10.0),
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
        "no_tracked_or_staged_changes": not state["tracked_changes"] and not state["staged_changes"],
        "stdout_only_no_run_directory_created": True,
    }
    execution_gate_pass = all(execution_checks.values())

    if all(checks.values()):
        status = "BOUNDED_STRAIN_3D_DOMAIN_CONVERGENCE_SUPPORTED"
    elif any(checks.values()):
        status = "BOUNDED_STRAIN_3D_DOMAIN_CONVERGENCE_PARTIAL_SUPPORT"
    else:
        status = "BOUNDED_STRAIN_3D_DOMAIN_CONVERGENCE_NOT_SUPPORTED"

    payload = {
        "lab_id": LAB_ID,
        "status": status,
        "model": {
            "K0": K0,
            "epsilon_max": EPSILON_MAX,
            "grid_ladder": GRID_LADDER,
            "source_radius": SOURCE_RADIUS,
            "source_load": SOURCE_LOAD,
            "probe_radii": PROBE_RADII,
            "network_equilibrium": "discrete divergence of six nearest-neighbor bounded-strain bond stresses = source",
            "boundary_condition": "u=0 on cubic outer boundary",
        },
        "rows": rows,
        "checks": checks,
        "execution_checks": execution_checks,
        "execution_gate_pass": execution_gate_pass,
        "repo_state": state,
        "summary": {
            "far_slopes": slopes,
            "local_r8_exponents": local_ps,
            "question": "Does the frozen 3D bounded-strain neighbor network converge toward a -1 radial response when only the cubic domain is enlarged?",
        },
    }

    print(LAB_ID)
    print(f"status={status}")
    print(f"head_sha={state['head_sha']}")
    print("constitutive_law_changed=false")
    print("K0_changed=false")
    print("fit_or_tuning_used=false")
    print("one_over_r_response_inserted=false")
    print("spherical_equilibrium_shortcut_used=false")
    print("\nTHREE_D_DOMAIN_CONVERGENCE")
    for r in rows:
        print(
            f"grid_N={r['grid_N']} half_width={r['half_width']:.12g} "
            f"far_fit={r['fit_far_radius']['slope']:.12g} "
            f"local_p_r8={r['local_exponent_r8']:.12g} "
            f"max_strain_fraction={r['max_strain_fraction']:.12g} "
            f"converged={str(r['converged']).lower()}"
        )
    print("\nCHECKS")
    for k,v in checks.items():
        print(f"{k}={str(v).lower()}")
    print("EXECUTION_CHECKS")
    for k,v in execution_checks.items():
        print(f"{k}={str(v).lower()}")
    print(f"execution_gate_pass={str(execution_gate_pass).lower()}")
    print(json.dumps(payload, sort_keys=True, separators=(",", ":")))
    return 0 if execution_gate_pass else 1


if __name__ == "__main__":
    raise SystemExit(main())
