#!/usr/bin/env python3
"""PBUF FOUNDATION — BOUNDED STRAIN 3D NEIGHBOR NETWORK 001.

Tests the bounded-strain nearest-neighbor constitutive hypothesis in a genuine
3D discrete network, without using the spherical equilibrium shortcut
r^2 sigma(e)=constant and without inserting a 1/r kernel.

Each nearest-neighbor bond has strain e = u_i-u_j and frozen constitutive law

    W(e)=-(K e_max^2/2) ln(1-(e/e_max)^2)
    sigma(e)=K e/(1-(e/e_max)^2)

The node equilibrium is the discrete divergence of bond stress = source.
A Picard iteration freezes bond secant stiffness, solves the resulting positive
variable-coefficient network by matrix-free conjugate gradient, then updates
bond stiffness from the new strain.  Zero Dirichlet boundary is numerical only.

No G, macroscopic amplitude calibration, fitted radial law, Rmax, cosmology,
lensing target, Quantum Engine, Planck input, or spherical-response equation is
used.
"""
from __future__ import annotations

import json
import math
import subprocess
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
LAB_ID = "PBUF-FOUNDATION-BOUNDED-STRAIN-3D-NEIGHBOR-NETWORK-001"

N = 49
DX = 1.0
CENTER = N // 2
SHAPE = (N, N, N)
AXIS = np.arange(N, dtype=np.float64) - CENTER
Z, Y, X = np.meshgrid(AXIS, AXIS, AXIS, indexing="ij")
RR = np.sqrt(X * X + Y * Y + Z * Z) * DX

K0 = 1.0
EPSILON_MAX = 1.0
PICARD_TOL = 2.0e-7
PICARD_MAX_ITER = 30
PICARD_DAMP = 0.65
CG_REL_TOL = 2.0e-9
CG_MAX_ITER = 500
EXP_WINDOW = 0.35
ADDITIVITY_TOL = 2.0e-3

RADIUS_LADDER = (2.5, 3.5, 4.5, 5.5, 6.5)
FIXED_LOAD = 2.0
MASS_LADDER = (0.5, 1.0, 2.0, 4.0, 8.0)
FIXED_RADIUS = 4.5
DENSITY_LADDER = (0.0025, 0.005, 0.01, 0.02, 0.04)
FAR_SOURCE_RADIUS = 3.5
FAR_PROBES = (6.0, 7.0, 8.0, 9.0, 10.0)
STRONG_LOAD = 12.0


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


def source_sphere(radius: float, integrated_load: float) -> dict:
    mask = RR <= radius
    count = int(np.count_nonzero(mask))
    source = np.zeros(SHAPE, dtype=np.float64)
    source[mask] = integrated_load / (count * DX ** 3)
    return {
        "source": source,
        "mask": mask,
        "density": float(source[mask][0]),
        "integrated_load": float(np.sum(source) * DX ** 3),
        "radius": radius,
    }


def source_density(radius: float, density: float) -> dict:
    mask = RR <= radius
    source = np.zeros(SHAPE, dtype=np.float64)
    source[mask] = density
    return {
        "source": source,
        "mask": mask,
        "density": density,
        "integrated_load": float(np.sum(source) * DX ** 3),
        "radius": radius,
    }


def shifted_source(offset_x: int, radius: float, integrated_load: float) -> np.ndarray:
    rr = np.sqrt((X-offset_x)**2 + Y*Y + Z*Z) * DX
    mask = rr <= radius
    out = np.zeros(SHAPE, dtype=np.float64)
    out[mask] = integrated_load / (np.count_nonzero(mask) * DX ** 3)
    return out


def zero_boundary(a: np.ndarray) -> None:
    a[0,:,:] = a[-1,:,:] = 0.0
    a[:,0,:] = a[:,-1,:] = 0.0
    a[:,:,0] = a[:,:,-1] = 0.0


def bond_strains(u: np.ndarray):
    return (
        (u[1:,:,:] - u[:-1,:,:]) / DX,
        (u[:,1:,:] - u[:,:-1,:]) / DX,
        (u[:,:,1:] - u[:,:,:-1]) / DX,
    )


def secant_weights(u: np.ndarray):
    es = bond_strains(u)
    out = []
    max_frac = 0.0
    for e in es:
        frac = np.abs(e) / EPSILON_MAX
        max_frac = max(max_frac, float(np.max(frac)))
        if max_frac >= 0.995:
            raise RuntimeError("bounded-strain barrier approached too closely")
        out.append(K0 / (1.0 - frac*frac))
    return tuple(out), max_frac


def apply_A(u: np.ndarray, weights) -> np.ndarray:
    """Positive weighted graph Laplacian on the 3D nearest-neighbor lattice."""
    wz, wy, wx = weights
    out = np.zeros_like(u)
    dz = (u[1:,:,:] - u[:-1,:,:]) * wz
    dy = (u[:,1:,:] - u[:,:-1,:]) * wy
    dx = (u[:,:,1:] - u[:,:,:-1]) * wx
    out[:-1,:,:] -= dz
    out[1:,:,:] += dz
    out[:,:-1,:] -= dy
    out[:,1:,:] += dy
    out[:,:,:-1] -= dx
    out[:,:,1:] += dx
    zero_boundary(out)
    return out / (DX*DX)


def cg_solve(source: np.ndarray, weights, x0: np.ndarray | None = None) -> dict:
    b = source.copy()
    zero_boundary(b)
    x = np.zeros_like(b) if x0 is None else x0.copy()
    zero_boundary(x)
    r = b - apply_A(x, weights)
    zero_boundary(r)
    p = r.copy()
    rr_old = float(np.sum(r*r))
    bnorm = math.sqrt(float(np.sum(b*b)))
    if bnorm == 0.0:
        return {"field": x, "iterations": 0, "relative_residual": 0.0}
    rel = math.sqrt(rr_old) / bnorm
    it = 0
    for k in range(1, CG_MAX_ITER+1):
        Ap = apply_A(p, weights)
        denom = float(np.sum(p*Ap))
        if not math.isfinite(denom) or denom <= 0.0:
            raise RuntimeError("CG operator lost positive definiteness")
        alpha = rr_old / denom
        x += alpha*p
        zero_boundary(x)
        r -= alpha*Ap
        zero_boundary(r)
        rr_new = float(np.sum(r*r))
        rel = math.sqrt(rr_new) / bnorm
        it = k
        if rel <= CG_REL_TOL:
            break
        beta = rr_new / rr_old
        p = r + beta*p
        zero_boundary(p)
        rr_old = rr_new
    return {"field": x, "iterations": it, "relative_residual": rel}


def solve_nonlinear(source: np.ndarray) -> dict:
    # Begin from the weak-strain K0 network; this is an initial iterate, not a fit.
    ones = (
        np.ones((N-1,N,N), dtype=np.float64)*K0,
        np.ones((N,N-1,N), dtype=np.float64)*K0,
        np.ones((N,N,N-1), dtype=np.float64)*K0,
    )
    first = cg_solve(source, ones)
    u = first["field"]
    converged = False
    history = []
    last_cg = first
    for it in range(1, PICARD_MAX_ITER+1):
        weights, before_frac = secant_weights(u)
        sol = cg_solve(source, weights, x0=u)
        candidate = sol["field"]
        new_u = PICARD_DAMP*candidate + (1.0-PICARD_DAMP)*u
        zero_boundary(new_u)
        scale = max(float(np.sqrt(np.mean(new_u*new_u))), 1.0e-14)
        change = float(np.sqrt(np.mean((new_u-u)**2))) / scale
        u = new_u
        _, after_frac = secant_weights(u)
        history.append({"picard": it, "relative_change": change, "max_strain_fraction": after_frac, "cg_iterations": sol["iterations"], "cg_relative_residual": sol["relative_residual"]})
        last_cg = sol
        if change <= PICARD_TOL:
            converged = True
            break
    weights, max_frac = secant_weights(u)
    residual = source - apply_A(u, weights)
    zero_boundary(residual)
    src_norm = math.sqrt(float(np.sum(source[1:-1,1:-1,1:-1]**2)))
    nl_rel = math.sqrt(float(np.sum(residual*residual))) / src_norm if src_norm > 0 else 0.0
    return {
        "field": u,
        "converged": converged,
        "picard_iterations": len(history),
        "max_strain_fraction": max_frac,
        "nonlinear_relative_residual": nl_rel,
        "last_cg_relative_residual": last_cg["relative_residual"],
        "history_tail": history[-3:],
    }


def shell_mean(field: np.ndarray, radius: float) -> float:
    m = (RR >= radius-0.75) & (RR <= radius+0.75)
    return float(np.mean(np.abs(field[m])))


def probe(field: np.ndarray, radius: float) -> float:
    off = int(round(radius/DX))
    return float(abs(field[CENTER,CENTER,CENTER+off]))


def logfit(xs, ys) -> dict:
    x = np.asarray(xs, dtype=np.float64)
    y = np.asarray(ys, dtype=np.float64)
    m = np.isfinite(x)&np.isfinite(y)&(x>0)&(y>0)
    lx, ly = np.log(x[m]), np.log(y[m])
    A = np.column_stack((lx, np.ones_like(lx)))
    beta, *_ = np.linalg.lstsq(A, ly, rcond=None)
    pred = A@beta
    ss_res = float(np.sum((ly-pred)**2))
    ss_tot = float(np.sum((ly-np.mean(ly))**2))
    return {"slope": float(beta[0]), "r2": 1.0-ss_res/ss_tot if ss_tot>0 else float("nan"), "count": int(np.count_nonzero(m))}


def match(v: float, target: float) -> bool:
    return bool(math.isfinite(v) and abs(v-target) <= EXP_WINDOW)


def relative_rms(a: np.ndarray, b: np.ndarray) -> float:
    num = math.sqrt(float(np.mean((a-b)**2)))
    den = math.sqrt(float(np.mean(b*b)))
    return num/den if den>0 else float("nan")


def density_test() -> dict:
    rows=[]
    for d in DENSITY_LADDER:
        src=source_density(FIXED_RADIUS,d)
        sol=solve_nonlinear(src["source"])
        rows.append({"density":d,"integrated_load":src["integrated_load"],"center_response":float(abs(sol["field"][CENTER,CENTER,CENTER])),"max_strain_fraction":sol["max_strain_fraction"],"converged":sol["converged"]})
    return {"rows":rows,"fit":logfit([r["density"] for r in rows],[r["center_response"] for r in rows])}


def mass_test() -> dict:
    rows=[]
    for mass in MASS_LADDER:
        src=source_sphere(FIXED_RADIUS,mass)
        sol=solve_nonlinear(src["source"])
        rows.append({"load":src["integrated_load"],"surface_response":shell_mean(sol["field"],FIXED_RADIUS),"far_response":probe(sol["field"],10.0),"max_strain_fraction":sol["max_strain_fraction"],"converged":sol["converged"]})
    return {"rows":rows,"surface_fit":logfit([r["load"] for r in rows],[r["surface_response"] for r in rows]),"far_fit":logfit([r["load"] for r in rows],[r["far_response"] for r in rows])}


def radius_test() -> dict:
    rows=[]
    for radius in RADIUS_LADDER:
        src=source_sphere(radius,FIXED_LOAD)
        sol=solve_nonlinear(src["source"])
        rows.append({"radius":radius,"integrated_load":src["integrated_load"],"surface_response":shell_mean(sol["field"],radius),"max_strain_fraction":sol["max_strain_fraction"],"converged":sol["converged"]})
    return {"rows":rows,"fit":logfit([r["radius"] for r in rows],[r["surface_response"] for r in rows])}


def far_radius_test() -> dict:
    src=source_sphere(FAR_SOURCE_RADIUS,FIXED_LOAD)
    sol=solve_nonlinear(src["source"])
    vals=[probe(sol["field"],r) for r in FAR_PROBES]
    return {"probe_radii":FAR_PROBES,"responses":vals,"fit":logfit(FAR_PROBES,vals),"max_strain_fraction":sol["max_strain_fraction"],"converged":sol["converged"]}


def additivity_test() -> dict:
    s1=shifted_source(-5,2.5,0.8)
    s2=shifted_source(+5,3.5,1.2)
    u1=solve_nonlinear(s1)["field"]
    u2=solve_nonlinear(s2)["field"]
    u12=solve_nonlinear(s1+s2)["field"]
    return {"relative_rms_residual":relative_rms(u12,u1+u2)}


def strong_load_test() -> dict:
    src=source_sphere(3.5,STRONG_LOAD)
    sol=solve_nonlinear(src["source"])
    return {"load":STRONG_LOAD,"max_strain_fraction":sol["max_strain_fraction"],"surface_response":shell_mean(sol["field"],3.5),"far_response":probe(sol["field"],10.0),"converged":sol["converged"],"picard_iterations":sol["picard_iterations"],"nonlinear_relative_residual":sol["nonlinear_relative_residual"],"history_tail":sol["history_tail"]}


def main() -> int:
    density=density_test()
    mass=mass_test()
    radius=radius_test()
    far=far_radius_test()
    add=additivity_test()
    strong=strong_load_test()
    state=repo_state()

    measured={
        "density_exponent":density["fit"]["slope"],
        "surface_mass_exponent":mass["surface_fit"]["slope"],
        "surface_radius_fixed_load_exponent":radius["fit"]["slope"],
        "far_mass_exponent":mass["far_fit"]["slope"],
        "far_radius_exponent":far["fit"]["slope"],
        "weak_additivity_relative_residual":add["relative_rms_residual"],
        "strong_load_max_strain_fraction":strong["max_strain_fraction"],
    }
    checks={
        "density_linearity":match(measured["density_exponent"],1.0),
        "surface_mass_linearity":match(measured["surface_mass_exponent"],1.0),
        "surface_radius_accumulation":match(measured["surface_radius_fixed_load_exponent"],-1.0),
        "far_mass_linearity":match(measured["far_mass_exponent"],1.0),
        "far_radius_accumulation":match(measured["far_radius_exponent"],-1.0),
        "weak_regime_additivity":bool(math.isfinite(measured["weak_additivity_relative_residual"]) and measured["weak_additivity_relative_residual"] <= ADDITIVITY_TOL),
        "strong_load_below_strain_limit":bool(strong["max_strain_fraction"] < 1.0),
    }
    six_keys=("density_linearity","surface_mass_linearity","surface_radius_accumulation","far_mass_linearity","far_radius_accumulation","weak_regime_additivity")
    six_count=sum(bool(checks[k]) for k in six_keys)
    all_converged=all(r["converged"] for r in density["rows"]+mass["rows"]+radius["rows"]) and far["converged"] and strong["converged"]
    execution_checks={
        "all_measured_values_finite":all(math.isfinite(v) for v in measured.values()),
        "nonlinear_solves_converged":all_converged,
        "constitutive_law_frozen":True,
        "K0_frozen":K0==1.0,
        "epsilon_max_frozen":EPSILON_MAX==1.0,
        "spherical_equilibrium_shortcut_used":False,
        "no_G":True,
        "no_macroscopic_amplitude":True,
        "no_native_rescaling":True,
        "no_fit_or_tuning":True,
        "no_inserted_one_over_r_response":True,
        "no_Rmax":True,
        "no_cosmology":True,
        "no_lensing_target":True,
        "no_quantum_engine":True,
        "no_planck_input":True,
        "no_tracked_or_staged_changes":state["tracked_changes"]=="" and state["staged_changes"]=="",
        "stdout_only_no_run_directory_created":True,
    }
    execution_gate_pass=all(execution_checks.values())
    if six_count==6 and checks["strong_load_below_strain_limit"]:
        status="BOUNDED_STRAIN_3D_NEIGHBOR_NETWORK_STRUCTURE_SUPPORTED"
    elif six_count>=4 and checks["strong_load_below_strain_limit"]:
        status="BOUNDED_STRAIN_3D_NEIGHBOR_NETWORK_PARTIAL_SUPPORT"
    else:
        status="BOUNDED_STRAIN_3D_NEIGHBOR_NETWORK_NOT_SUPPORTED"

    payload={
        "lab_id":LAB_ID,
        "status":status,
        "repo_state":state,
        "model":{
            "grid_N":N,"dx":DX,"K0":K0,"epsilon_max":EPSILON_MAX,
            "constitutive_energy":"W(e)=-(K e_max^2/2) ln(1-(e/e_max)^2)",
            "constitutive_stress":"sigma(e)=K e/(1-(e/e_max)^2)",
            "network_equilibrium":"discrete divergence of six nearest-neighbor bond stresses = source",
            "spherical_equilibrium_shortcut":False,
        },
        "measured":measured,"checks":checks,"matched_frozen_checks_of_6":six_count,
        "density_test":density,"mass_test":mass,"radius_fixed_load_test":radius,"far_radius_test":far,"additivity_test":add,"strong_load_test":strong,
        "execution_checks":execution_checks,"execution_gate_pass":execution_gate_pass,
        "policy":{"gravity_used_as_native_variable":False,"G_used":False,"one_over_r_response_inserted":False,"fit_or_tuning_used":False,"spherical_response_equation_used":False},
        "summary":{
            "question":"Does the frozen bounded-strain bond law reproduce the required accumulation fingerprint in a true 3D nearest-neighbor network without the spherical equilibrium shortcut?",
            "next_if_supported":"Connect the existing local-loading candidate to the 3D network source term without amplitude fitting and repeat the frozen fingerprint audit.",
            "next_if_partial":"Localize which 3D fingerprint fails before changing the constitutive law or network topology.",
        },
    }

    print(LAB_ID)
    print(f"status={status}")
    print(f"head_sha={state['head_sha']}")
    print("spherical_equilibrium_shortcut_used=false")
    print("constitutive_law_changed=false")
    print("G_used=false")
    print("one_over_r_response_inserted=false")
    print("fit_or_tuning_used=false")
    print("\nFROZEN_RESPONSE_FINGERPRINT")
    for k,v in measured.items(): print(f"{k}={v:.12g}")
    print(f"matched_frozen_checks_of_6={six_count}")
    print("\nCHECKS")
    for k,v in checks.items(): print(f"{k}={str(v).lower()}")
    print("\nSTRONG_LOAD")
    print(f"load={STRONG_LOAD:g} max_strain_fraction={strong['max_strain_fraction']:.12g} converged={str(strong['converged']).lower()} nonlinear_relative_residual={strong['nonlinear_relative_residual']:.12g}")
    print("\nEXECUTION_CHECKS")
    for k,v in execution_checks.items(): print(f"{k}={str(v).lower()}")
    print(f"execution_gate_pass={str(execution_gate_pass).lower()}")
    print(json.dumps(payload,sort_keys=True,separators=(",",":"),allow_nan=False))
    return 0 if execution_gate_pass else 1


if __name__ == "__main__":
    raise SystemExit(main())
