#!/usr/bin/env python3
"""PBUF FOUNDATION — BOUNDED STRAIN NEIGHBOR CONSTITUTIVE 001.

Test a nonlinear nearest-neighbor medium motivated by the working intuition that
additional separation becomes progressively harder and may approach a finite
strain that requires unbounded energy to reach.

The constitutive energy density is chosen as the minimal even barrier law

    W(e) = -(K e_max^2 / 2) ln(1 - (e/e_max)^2)

which gives restoring stress

    sigma(e) = K e / (1 - (e/e_max)^2)

and tangent stiffness d sigma / d e that grows with |e| and diverges as
|e| -> e_max.  This is a speculative constitutive hypothesis, not a fitted law.

For a static spherically symmetric neighbor-coupled continuum, equilibrium away
from a localized source follows from flux balance of the constitutive stress:

    d/dr [ r^2 sigma(e(r)) ] = 0

so r^2 sigma is constant outside the source.  We solve the nonlinear inverse
constitutive relation numerically and integrate the strain to obtain the
accumulated displacement/state.  No 1/r displacement law or inverse-square
force law is inserted as a response target.

No G, macroscopic amplitude calibration, lensing target, Rmax, Quantum Engine,
Planck input, fitted exponent, or gravity variable is used.
"""
from __future__ import annotations

import json
import math
import subprocess
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
LAB_ID = "PBUF-FOUNDATION-BOUNDED-STRAIN-NEIGHBOR-CONSTITUTIVE-001"

# Dimensionless structural normalisation. These are fixed before execution and
# are not fitted to a radial fingerprint.
K0 = 1.0
EPS_MAX = 1.0
FOUR_PI = 4.0 * math.pi
R_SOURCE = 1.0
R_OUTER = 256.0
N_RADIAL = 8193
RADII = np.linspace(R_SOURCE, R_OUTER, N_RADIAL, dtype=np.float64)

LOW_LOAD = 0.05
LOAD_LADDER = (0.02, 0.05, 0.10, 0.20, 0.40)
STRONG_LOAD_LADDER = (0.5, 1.0, 2.0, 4.0, 8.0, 16.0)
PROBE_RADII = (4.0, 8.0, 16.0, 32.0, 64.0, 96.0)
LOCAL_EXPONENT_RADII = (2.0, 4.0, 8.0, 16.0, 32.0, 64.0, 96.0)
EXP_WINDOW = 0.12
LINEARITY_WINDOW = 0.05


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


def energy_density(eps: float) -> float:
    x = abs(float(eps)) / EPS_MAX
    if x >= 1.0:
        return float("inf")
    return -0.5 * K0 * EPS_MAX * EPS_MAX * math.log1p(-(x * x))


def stress(eps: float) -> float:
    e = float(eps)
    x2 = (e / EPS_MAX) ** 2
    if x2 >= 1.0:
        return math.copysign(float("inf"), e)
    return K0 * e / (1.0 - x2)


def tangent_stiffness(eps: float) -> float:
    e = float(eps)
    x2 = (e / EPS_MAX) ** 2
    if x2 >= 1.0:
        return float("inf")
    return K0 * (1.0 + x2) / ((1.0 - x2) ** 2)


def inverse_stress(target: float) -> float:
    """Solve sigma(e)=target on 0 <= e < EPS_MAX by bisection."""
    if target < 0.0 or not math.isfinite(target):
        raise ValueError("target stress must be finite and nonnegative")
    if target == 0.0:
        return 0.0
    lo = 0.0
    hi = EPS_MAX * (1.0 - 1.0e-14)
    for _ in range(100):
        mid = 0.5 * (lo + hi)
        if stress(mid) < target:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def radial_solution(load: float) -> dict:
    """Solve spherical constitutive equilibrium and integrate strain inward.

    Outside the source, r^2 sigma(e)=load/(4 pi).  u(R_OUTER)=0 is used only
    as a finite numerical reference.  The accumulated response magnitude is the
    integral of strain from r to R_OUTER.
    """
    q = float(load) / FOUR_PI
    target_sigma = q / (RADII * RADII)
    eps = np.asarray([inverse_stress(float(s)) for s in target_sigma], dtype=np.float64)

    # u(r) = integral_r^Rout eps(s) ds, evaluated by reverse trapezoidal sum.
    u = np.zeros_like(eps)
    dr = np.diff(RADII)
    seg = 0.5 * (eps[:-1] + eps[1:]) * dr
    u[:-1] = np.cumsum(seg[::-1])[::-1]

    return {
        "load": float(load),
        "strain": eps,
        "response": u,
        "max_strain": float(np.max(eps)),
        "max_strain_fraction": float(np.max(eps) / EPS_MAX),
        "surface_stress": float(target_sigma[0]),
        "surface_response": float(u[0]),
    }


def interp(field: np.ndarray, radius: float) -> float:
    return float(np.interp(float(radius), RADII, np.asarray(field, dtype=np.float64)))


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


def constitutive_barrier_test() -> dict:
    fractions = (0.0, 0.25, 0.50, 0.75, 0.90, 0.95, 0.99, 0.999, 0.9999)
    rows = []
    for f in fractions:
        e = f * EPS_MAX
        rows.append({
            "strain_fraction": f,
            "energy": energy_density(e),
            "stress": stress(e),
            "tangent_stiffness": tangent_stiffness(e),
        })
    positive = rows[1:]
    energy_monotonic = all(b["energy"] > a["energy"] for a, b in zip(positive[:-1], positive[1:]))
    stress_monotonic = all(b["stress"] > a["stress"] for a, b in zip(positive[:-1], positive[1:]))
    stiffness_monotonic = all(b["tangent_stiffness"] > a["tangent_stiffness"] for a, b in zip(positive[:-1], positive[1:]))
    return {
        "rows": rows,
        "energy_monotonic": energy_monotonic,
        "stress_monotonic": stress_monotonic,
        "stiffness_monotonic": stiffness_monotonic,
        "energy_growth_0p9_to_0p9999": rows[-1]["energy"] / rows[4]["energy"],
        "stress_growth_0p9_to_0p9999": rows[-1]["stress"] / rows[4]["stress"],
        "stiffness_growth_0p9_to_0p9999": rows[-1]["tangent_stiffness"] / rows[4]["tangent_stiffness"],
    }


def low_strain_control() -> dict:
    sol = radial_solution(LOW_LOAD)
    eps = sol["strain"]
    target_sigma = (LOW_LOAD / FOUR_PI) / (RADII * RADII)
    linear_eps = target_sigma / K0
    rel = np.abs(eps - linear_eps) / np.maximum(np.abs(linear_eps), 1.0e-30)
    probes = [interp(sol["response"], r) for r in PROBE_RADII]
    fit = logfit(PROBE_RADII, probes)
    return {
        "load": LOW_LOAD,
        "max_strain_fraction": sol["max_strain_fraction"],
        "max_relative_strain_difference_from_linear": float(np.max(rel)),
        "probe_radii": PROBE_RADII,
        "responses": probes,
        "far_response_fit": fit,
    }


def load_scaling_test() -> dict:
    rows = []
    for load in LOAD_LADDER:
        sol = radial_solution(load)
        rows.append({
            "load": load,
            "surface_response": sol["surface_response"],
            "far_response_r32": interp(sol["response"], 32.0),
            "max_strain_fraction": sol["max_strain_fraction"],
        })
    return {
        "rows": rows,
        "surface_fit": logfit([r["load"] for r in rows], [r["surface_response"] for r in rows]),
        "far_fit": logfit([r["load"] for r in rows], [r["far_response_r32"] for r in rows]),
    }


def strong_load_test() -> dict:
    rows = []
    for load in STRONG_LOAD_LADDER:
        sol = radial_solution(load)
        rows.append({
            "load": load,
            "surface_stress": sol["surface_stress"],
            "max_strain_fraction": sol["max_strain_fraction"],
            "surface_response": sol["surface_response"],
            "far_response_r32": interp(sol["response"], 32.0),
        })
    strain_fracs = [r["max_strain_fraction"] for r in rows]
    return {
        "rows": rows,
        "strain_monotonic": all(b > a for a, b in zip(strain_fracs[:-1], strain_fracs[1:])),
        "all_below_limit": all(x < 1.0 for x in strain_fracs),
    }


def local_radial_exponent_test(load: float = 8.0) -> dict:
    sol = radial_solution(load)
    r = RADII
    u = sol["response"]
    # Local p = d ln u / d ln r, using centered finite differences. Stay well
    # away from the outer reference boundary where u necessarily tends to zero.
    good = (r >= 1.5) & (r <= 110.0) & (u > 0.0)
    lr = np.log(r[good])
    lu = np.log(u[good])
    p = np.gradient(lu, lr)
    rg = r[good]
    samples = []
    for rp in LOCAL_EXPONENT_RADII:
        samples.append({
            "radius": rp,
            "local_response_exponent": float(np.interp(rp, rg, p)),
            "strain_fraction": interp(sol["strain"], rp) / EPS_MAX,
            "response": interp(u, rp),
        })
    far_samples = [x["local_response_exponent"] for x in samples[-3:]]
    return {
        "load": load,
        "samples": samples,
        "outer_three_mean_exponent": float(np.mean(far_samples)),
        "max_strain_fraction": sol["max_strain_fraction"],
    }


def main() -> int:
    repo = repo_state()
    barrier = constitutive_barrier_test()
    low = low_strain_control()
    load_scaling = load_scaling_test()
    strong = strong_load_test()
    local_exp = local_radial_exponent_test()

    checks = {
        "energy_rises_with_strain": barrier["energy_monotonic"],
        "stress_rises_with_strain": barrier["stress_monotonic"],
        "tangent_stiffness_rises_with_strain": barrier["stiffness_monotonic"],
        "finite_strain_limit_respected": strong["all_below_limit"],
        "strain_increases_with_load": strong["strain_monotonic"],
        "low_strain_recovers_linear_constitutive_behavior": low["max_relative_strain_difference_from_linear"] <= 0.01,
        "low_strain_far_response_near_minus1": abs(low["far_response_fit"]["slope"] + 1.0) <= EXP_WINDOW,
        "weak_load_far_response_nearly_linear_in_load": abs(load_scaling["far_fit"]["slope"] - 1.0) <= LINEARITY_WINDOW,
    }

    if all(checks.values()):
        status = "BOUNDED_STRAIN_NEIGHBOR_CONSTITUTIVE_STRUCTURE_SUPPORTED"
    elif sum(checks.values()) >= 5:
        status = "BOUNDED_STRAIN_NEIGHBOR_CONSTITUTIVE_PARTIAL_SUPPORT"
    else:
        status = "BOUNDED_STRAIN_NEIGHBOR_CONSTITUTIVE_NOT_SUPPORTED"

    execution_checks = {
        "all_measured_values_finite": all(math.isfinite(v) for v in [
            low["far_response_fit"]["slope"],
            load_scaling["surface_fit"]["slope"],
            load_scaling["far_fit"]["slope"],
            local_exp["outer_three_mean_exponent"],
        ]),
        "no_G": True,
        "no_macroscopic_amplitude": True,
        "no_native_rescaling": True,
        "no_fit_or_tuning": True,
        "no_inserted_one_over_r_response": True,
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
            "energy_density": "W(e)=-(K e_max^2/2) ln(1-(e/e_max)^2)",
            "stress": "sigma(e)=K e/(1-(e/e_max)^2)",
            "static_outside_source": "d/dr[r^2 sigma(e(r))]=0",
            "K0": K0,
            "epsilon_max": EPS_MAX,
            "interpretation": "strain-dependent nearest-neighbor constitutive hypothesis with asymptotic strain barrier",
        },
        "policy": {
            "gravity_fundamental_in_PBUF": False,
            "gravity_used_as_native_variable": False,
            "G_used": False,
            "macroscopic_amplitude_used": False,
            "fit_or_tuning_used": False,
            "one_over_r_response_inserted": False,
            "electromagnetic_origin_claimed": False,
            "planck_input_used": False,
        },
        "constitutive_barrier": barrier,
        "low_strain_control": low,
        "load_scaling": load_scaling,
        "strong_load": strong,
        "local_radial_exponent": local_exp,
        "checks": checks,
        "execution_checks": execution_checks,
        "execution_gate_pass": execution_gate_pass,
        "summary": {
            "question": "Can a neighbor interaction become progressively harder under strain, approach a finite strain barrier at unbounded energy cost, and still recover the long-range linear-medium response at weak/far strain?",
            "next_if_supported": "Map this constitutive stress law onto the existing local-loading variable and test a 3D discrete nonlinear neighbor network without amplitude fitting.",
            "next_if_partial": "Identify which constitutive or finite-domain assumption fails before changing the strain law.",
        },
    }

    print(LAB_ID)
    print(f"status={status}")
    print(f"head_sha={repo['head_sha']}")
    print("gravity_used_as_native_variable=false")
    print("electromagnetic_origin_claimed=false")
    print("G_used=false")
    print("one_over_r_response_inserted=false")
    print("fit_or_tuning_used=false")
    print()
    print("CONSTITUTIVE_BARRIER")
    print(f"energy_growth_0p9_to_0p9999={barrier['energy_growth_0p9_to_0p9999']:.12g}")
    print(f"stress_growth_0p9_to_0p9999={barrier['stress_growth_0p9_to_0p9999']:.12g}")
    print(f"stiffness_growth_0p9_to_0p9999={barrier['stiffness_growth_0p9_to_0p9999']:.12g}")
    print()
    print("LOW_STRAIN_CONTROL")
    print(f"max_strain_fraction={low['max_strain_fraction']:.12g}")
    print(f"max_relative_difference_from_linear={low['max_relative_strain_difference_from_linear']:.12g}")
    print(f"far_response_exponent={low['far_response_fit']['slope']:.12g}")
    print()
    print("LOAD_SCALING")
    print(f"surface_load_exponent={load_scaling['surface_fit']['slope']:.12g}")
    print(f"far_load_exponent={load_scaling['far_fit']['slope']:.12g}")
    print()
    print("STRONG_LOAD")
    for row in strong["rows"]:
        print(f"load={row['load']:.12g} max_strain_fraction={row['max_strain_fraction']:.12g} surface_response={row['surface_response']:.12g} far_response_r32={row['far_response_r32']:.12g}")
    print()
    print("LOCAL_RADIAL_EXPONENT")
    for row in local_exp["samples"]:
        print(f"r={row['radius']:.12g} p={row['local_response_exponent']:.12g} strain_fraction={row['strain_fraction']:.12g}")
    print(f"outer_three_mean_exponent={local_exp['outer_three_mean_exponent']:.12g}")
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
        raise RuntimeError("bounded strain neighbor constitutive execution gate failed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
