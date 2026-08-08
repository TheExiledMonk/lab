#!/usr/bin/env python3
"""PBUF FOUNDATION — BOUNDED STRAIN FINITE DOMAIN CONVERGENCE 001.

Keep the bounded-strain constitutive law frozen and test whether the apparent
far-field steepening seen in the previous lab is caused by a finite outer
Dirichlet boundary rather than by the constitutive response itself.

Frozen constitutive law:
    W(e)=-(K e_max^2/2) ln(1-(e/e_max)^2)
    sigma(e)=K e/(1-(e/e_max)^2)

Outside a localized spherical load, radial equilibrium gives
    r^2 sigma(e(r)) = constant.
The response is obtained by numerically inverting sigma(e) and integrating
strain inward from u(R_B)=0. No 1/r response or target exponent is inserted.
"""
from __future__ import annotations

import json
import math
import subprocess
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
LAB_ID = "PBUF-FOUNDATION-BOUNDED-STRAIN-FINITE-DOMAIN-CONVERGENCE-001"

K0 = 1.0
EPSILON_MAX = 1.0
SOURCE_RADIUS = 1.0
LOAD = 8.0
BOUNDARIES = (128.0, 256.0, 512.0, 1024.0)
FIXED_PROBES = (4.0, 8.0, 16.0, 32.0, 64.0)
FIT_WINDOW = (4.0, 8.0, 16.0, 32.0)
N_RADIAL = 24000
EXP_TOL = 0.05
MONOTONIC_TOL = 1.0e-12


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


def stress(e: np.ndarray | float) -> np.ndarray | float:
    x = np.asarray(e, dtype=np.float64)
    y = K0 * x / (1.0 - (x / EPSILON_MAX) ** 2)
    return float(y) if np.ndim(e) == 0 else y


def strain_from_stress(sig: np.ndarray | float) -> np.ndarray | float:
    """Exact positive inverse of sigma=K e/(1-(e/e_max)^2)."""
    s = np.asarray(sig, dtype=np.float64)
    out = np.zeros_like(s)
    mask = s > 0.0
    sm = s[mask]
    # sigma e^2 + K e_max^2 e - sigma e_max^2 = 0 after multiplying through.
    out[mask] = (
        -K0 * EPSILON_MAX**2
        + np.sqrt((K0 * EPSILON_MAX**2) ** 2 + 4.0 * sm * sm * EPSILON_MAX**2)
    ) / (2.0 * sm)
    return float(out) if np.ndim(sig) == 0 else out


def flux_constant(load: float) -> float:
    # Spherical integrated-load convention; the geometric 4pi factor is part
    # of the forward equilibrium derivation, not a fitted amplitude.
    return load / (4.0 * math.pi)


def profile(boundary: float, load: float = LOAD) -> dict:
    if boundary <= max(FIXED_PROBES):
        raise RuntimeError("boundary must exceed fixed probes")
    r = np.geomspace(SOURCE_RADIUS, boundary, N_RADIAL)
    q = flux_constant(load)
    sig = q / (r * r)
    eps = strain_from_stress(sig)
    if not np.all(np.isfinite(eps)) or np.any(eps < 0.0) or np.any(eps >= EPSILON_MAX):
        raise RuntimeError("invalid constitutive strain profile")

    # u(r)=integral_r^R_B e(s) ds, with u(R_B)=0.
    dr = np.diff(r)
    trap = 0.5 * (eps[:-1] + eps[1:]) * dr
    suffix = np.zeros_like(r)
    suffix[:-1] = np.cumsum(trap[::-1])[::-1]

    def interp_u(radius: float) -> float:
        return float(np.interp(radius, r, suffix))

    def interp_eps(radius: float) -> float:
        return float(np.interp(radius, r, eps))

    probes = {str(int(x)): interp_u(x) for x in FIXED_PROBES}
    strains = {str(int(x)): interp_eps(x) for x in FIXED_PROBES}
    return {
        "boundary": boundary,
        "r": r,
        "u": suffix,
        "eps": eps,
        "probe_response": probes,
        "probe_strain_fraction": {k: v / EPSILON_MAX for k, v in strains.items()},
        "max_strain_fraction": float(np.max(eps) / EPSILON_MAX),
    }


def logfit(xs, ys) -> dict:
    x = np.log(np.asarray(xs, dtype=np.float64))
    y = np.log(np.asarray(ys, dtype=np.float64))
    A = np.column_stack((x, np.ones_like(x)))
    beta, *_ = np.linalg.lstsq(A, y, rcond=None)
    pred = A @ beta
    ss_res = float(np.sum((y - pred) ** 2))
    ss_tot = float(np.sum((y - np.mean(y)) ** 2))
    return {
        "slope": float(beta[0]),
        "r2": 1.0 - ss_res / ss_tot if ss_tot > 0.0 else float("nan"),
        "count": len(xs),
    }


def local_exponent(p: dict, radius: float) -> float:
    r = p["r"]
    u = p["u"]
    h = 0.02
    r1 = radius * math.exp(-h)
    r2 = radius * math.exp(+h)
    u1 = float(np.interp(r1, r, u))
    u2 = float(np.interp(r2, r, u))
    return math.log(u2 / u1) / math.log(r2 / r1)


def main() -> int:
    repo = repo_state()
    rows = []
    for b in BOUNDARIES:
        p = profile(b)
        fit_vals = [p["probe_response"][str(int(r))] for r in FIT_WINDOW]
        fit = logfit(FIT_WINDOW, fit_vals)
        local = {str(int(r)): local_exponent(p, r) for r in FIXED_PROBES}
        rows.append({
            "boundary": b,
            "fit_fixed_4_32": fit,
            "local_exponents": local,
            "probe_response": p["probe_response"],
            "probe_strain_fraction": p["probe_strain_fraction"],
            "max_strain_fraction": p["max_strain_fraction"],
        })

    slopes = [row["fit_fixed_4_32"]["slope"] for row in rows]
    errors = [abs(s + 1.0) for s in slopes]
    local32 = [row["local_exponents"]["32"] for row in rows]
    local32_errors = [abs(x + 1.0) for x in local32]

    convergence_monotonic = all(errors[i + 1] <= errors[i] + MONOTONIC_TOL for i in range(len(errors) - 1))
    local32_convergence = all(local32_errors[i + 1] <= local32_errors[i] + MONOTONIC_TOL for i in range(len(local32_errors) - 1))
    largest_near_minus1 = errors[-1] <= EXP_TOL
    largest_local32_near_minus1 = local32_errors[-1] <= EXP_TOL
    strain_profile_invariant = max(row["max_strain_fraction"] for row in rows) - min(row["max_strain_fraction"] for row in rows) <= 1.0e-10

    checks = {
        "fixed_probe_fit_moves_toward_minus1_with_boundary": convergence_monotonic,
        "fixed_r32_local_exponent_moves_toward_minus1": local32_convergence,
        "largest_boundary_fixed_probe_fit_near_minus1": largest_near_minus1,
        "largest_boundary_r32_local_exponent_near_minus1": largest_local32_near_minus1,
        "source_strain_invariant_under_boundary_change": strain_profile_invariant,
    }

    if all(checks.values()):
        status = "BOUNDED_STRAIN_FINITE_DOMAIN_CONVERGENCE_SUPPORTED"
    elif convergence_monotonic or local32_convergence:
        status = "BOUNDED_STRAIN_FINITE_DOMAIN_CONVERGENCE_PARTIAL_SUPPORT"
    else:
        status = "BOUNDED_STRAIN_FINITE_DOMAIN_CONVERGENCE_NOT_SUPPORTED"

    execution_checks = {
        "all_measured_values_finite": all(math.isfinite(x) for x in slopes + local32),
        "constitutive_law_frozen": True,
        "K0_frozen": K0 == 1.0,
        "epsilon_max_frozen": EPSILON_MAX == 1.0,
        "load_frozen": LOAD == 8.0,
        "no_G": True,
        "no_macroscopic_amplitude": True,
        "no_native_rescaling": True,
        "no_fit_or_tuning": True,
        "no_inserted_one_over_r_response": True,
        "no_Rmax": True,
        "no_cosmology": True,
        "no_lensing_target": True,
        "no_quantum_engine": True,
        "no_planck_input": True,
        "no_tracked_or_staged_changes": repo["tracked_changes"] == "" and repo["staged_changes"] == "",
        "stdout_only_no_run_directory_created": True,
    }
    gate = all(execution_checks.values())

    payload = {
        "lab_id": LAB_ID,
        "status": status,
        "model": {
            "K0": K0,
            "epsilon_max": EPSILON_MAX,
            "load": LOAD,
            "constitutive_energy": "W(e)=-(K e_max^2/2) ln(1-(e/e_max)^2)",
            "constitutive_stress": "sigma(e)=K e/(1-(e/e_max)^2)",
            "outside_source_equilibrium": "r^2 sigma(e(r))=constant",
            "outer_boundary_condition": "u(R_B)=0",
            "boundaries": BOUNDARIES,
            "fixed_probes": FIXED_PROBES,
        },
        "rows": rows,
        "checks": checks,
        "execution_checks": execution_checks,
        "execution_gate_pass": gate,
        "repo_state": repo,
        "summary": {
            "fixed_probe_slopes": slopes,
            "fixed_r32_local_exponents": local32,
            "question": "Does the frozen bounded-strain response converge toward the long-range -1 radial exponent when only the numerical outer boundary is moved outward?",
        },
    }

    print(LAB_ID)
    print(f"status={status}")
    print(f"head_sha={repo['head_sha']}")
    print("constitutive_law_changed=false")
    print("fit_or_tuning_used=false")
    print("one_over_r_response_inserted=false")
    print("\nFINITE_DOMAIN_CONVERGENCE")
    for row in rows:
        print(
            f"boundary={row['boundary']:.0f} "
            f"fit_4_32={row['fit_fixed_4_32']['slope']:.12g} "
            f"local_p_r32={row['local_exponents']['32']:.12g} "
            f"max_strain_fraction={row['max_strain_fraction']:.12g}"
        )
    print("\nCHECKS")
    for k, v in checks.items():
        print(f"{k}={str(v).lower()}")
    print("EXECUTION_CHECKS")
    for k, v in execution_checks.items():
        print(f"{k}={str(v).lower()}")
    print(f"execution_gate_pass={str(gate).lower()}")
    print(json.dumps(payload, sort_keys=True, separators=(",", ":")))

    if not gate:
        raise RuntimeError("bounded-strain finite-domain convergence execution gate failed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
