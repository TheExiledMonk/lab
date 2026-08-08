#!/usr/bin/env python3
"""PBUF FOUNDATION — C_STATE BOUNDED-STRAIN BRIDGE 001.

Tests whether the existing native c_state local-loading candidate can serve as
the source term for the bounded-strain 3D nearest-neighbor accumulation network.

The native source chain is:
    rho -> existing A8 transport -> raw c_state

The accumulated-response chain is:
    raw c_state -> six-neighbor bounded-strain equilibrium -> accumulated field

No amplitude normalization, G, fitted K, inserted 1/r law, spherical-response
shortcut, Rmax, cosmology, lensing target, Quantum Engine, or Planck input is
used. K0=1 and epsilon_max=1 remain dimensionless structural normalizations.
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

LAB_ID = "PBUF-FOUNDATION-C-STATE-BOUNDED-STRAIN-BRIDGE-001"

# Native c_state is generated on the existing 33^3 synthetic geometry and then
# embedded without rescaling into a larger accumulation domain so the network
# boundary is farther from the fixed physical probes.
N_NATIVE = 33
N = 65
DX = 1.0
CENTER_NATIVE = N_NATIVE // 2
CENTER = N // 2

AXN = np.arange(N_NATIVE, dtype=np.float64) - CENTER_NATIVE
ZN, YN, XN = np.meshgrid(AXN, AXN, AXN, indexing="ij")
RR_NATIVE = np.sqrt(XN*XN + YN*YN + ZN*ZN) * DX

AX = np.arange(N, dtype=np.float64) - CENTER
Z, Y, X = np.meshgrid(AX, AX, AX, indexing="ij")
RR = np.sqrt(X*X + Y*Y + Z*Z) * DX

K0 = 1.0
EPSILON_MAX = 1.0
PICARD_TOL = 2.0e-7
PICARD_MAX_ITER = 30
PICARD_DAMP = 0.65
CG_REL_TOL = 2.0e-9
CG_MAX_ITER = 700
EXP_WINDOW = 0.35
ADDITIVITY_TOL = 2.0e-3

RADIUS_LADDER = (2.5, 3.5, 4.5, 5.5, 6.5)
FIXED_NATIVE_MASS = 2.0
MASS_LADDER = (0.5, 1.0, 2.0, 4.0, 8.0)
FIXED_RADIUS = 4.5
DENSITY_LADDER = (0.0025, 0.005, 0.01, 0.02, 0.04)
FAR_SOURCE_RADIUS = 3.5
FAR_PROBES = (6.0, 7.0, 8.0, 9.0, 10.0)
STRONG_NATIVE_MASS = 12.0


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


def native_rho_sphere(radius: float, density: float) -> np.ndarray:
    rho = np.zeros((N_NATIVE, N_NATIVE, N_NATIVE), dtype=np.float64)
    rho[RR_NATIVE <= radius] = float(density)
    return rho


def native_rho_fixed_mass(radius: float, mass: float) -> np.ndarray:
    mask = RR_NATIVE <= radius
    count = int(np.count_nonzero(mask))
    if count <= 0:
        raise RuntimeError("empty native source sphere")
    rho = np.zeros((N_NATIVE, N_NATIVE, N_NATIVE), dtype=np.float64)
    rho[mask] = float(mass) / (count * DX**3)
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


def embed_raw(native: np.ndarray) -> np.ndarray:
    """Embed raw native c_state into accumulation grid; no amplitude rescaling."""
    out = np.zeros((N, N, N), dtype=np.float64)
    lo = CENTER - CENTER_NATIVE
    hi = lo + N_NATIVE
    out[lo:hi, lo:hi, lo:hi] = np.asarray(native, dtype=np.float64)
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
    strains = bond_strains(u)
    weights = []
    max_frac = 0.0
    for e in strains:
        frac = np.abs(e) / EPSILON_MAX
        max_frac = max(max_frac, float(np.max(frac)))
        if max_frac >= 0.995:
            raise RuntimeError("bounded-strain barrier approached too closely")
        weights.append(K0 / (1.0 - frac*frac))
    return tuple(weights), max_frac


def apply_A(u: np.ndarray, weights) -> np.ndarray:
    wz, wy, wx = weights
    out = np.zeros_like(u)
    dz = (u[1:,:,:] - u[:-1,:,:]) * wz
    dy = (u[:,1:,:] - u[:,:-1,:]) * wy
    dx = (u[:,:,1:] - u[:,:,:-1]) * wx
    out[:-1,:,:] -= dz; out[1:,:,:] += dz
    out[:,:-1,:] -= dy; out[:,1:,:] += dy
    out[:,:,:-1] -= dx; out[:,:,1:] += dx
    zero_boundary(out)
    return out / (DX*DX)


def cg_solve(source: np.ndarray, weights, x0=None) -> dict:
    b = source.copy(); zero_boundary(b)
    x = np.zeros_like(b) if x0 is None else x0.copy(); zero_boundary(x)
    r = b - apply_A(x, weights); zero_boundary(r)
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
        x += alpha*p; zero_boundary(x)
        r -= alpha*Ap; zero_boundary(r)
        rr_new = float(np.sum(r*r))
        rel = math.sqrt(rr_new) / bnorm
        it = k
        if rel <= CG_REL_TOL:
            break
        beta = rr_new / rr_old
        p = r + beta*p; zero_boundary(p)
        rr_old = rr_new
    return {"field": x, "iterations": it, "relative_residual": rel}


def solve_nonlinear(source: np.ndarray) -> dict:
    ones = (
        np.ones((N-1,N,N), dtype=np.float64)*K0,
        np.ones((N,N-1,N), dtype=np.float64)*K0,
        np.ones((N,N,N-1), dtype=np.float64)*K0,
    )
    first = cg_solve(source, ones)
    u = first["field"]
    history = []
    converged = False
    for it in range(1, PICARD_MAX_ITER+1):
        weights, _ = secant_weights(u)
        sol = cg_solve(source, weights, x0=u)
        candidate = sol["field"]
        new_u = PICARD_DAMP*candidate + (1.0-PICARD_DAMP)*u
        zero_boundary(new_u)
        scale = max(float(np.sqrt(np.mean(new_u*new_u))), 1.0e-14)
        change = float(np.sqrt(np.mean((new_u-u)**2))) / scale
        u = new_u
        _, max_frac = secant_weights(u)
        history.append({"picard":it,"relative_change":change,"max_strain_fraction":max_frac,"cg_iterations":sol["iterations"],"cg_relative_residual":sol["relative_residual"]})
        if change <= PICARD_TOL:
            converged = True
            break
    weights, max_frac = secant_weights(u)
    residual = source - apply_A(u, weights); zero_boundary(residual)
    src_norm = math.sqrt(float(np.sum(source[1:-1,1:-1,1:-1]**2)))
    nl_rel = math.sqrt(float(np.sum(residual*residual))) / src_norm if src_norm > 0 else 0.0
    return {"field":u,"converged":converged,"picard_iterations":len(history),"max_strain_fraction":max_frac,"nonlinear_relative_residual":nl_rel,"history_tail":history[-3:]}


def accumulated_from_rho(rho: np.ndarray) -> dict:
    c = native_c_state(rho)
    source = embed_raw(c)
    sol = solve_nonlinear(source)
    return {
        **sol,
        "c_state": c,
        "c_state_integral": float(np.sum(c)*DX**3),
        "rho_integral": float(np.sum(rho)*DX**3),
        "c_state_max_abs": float(np.max(np.abs(c))),
    }


def shell_mean(field: np.ndarray, radius: float) -> float:
    mask = (RR >= radius-0.75) & (RR <= radius+0.75)
    return float(np.mean(np.abs(field[mask])))


def probe(field: np.ndarray, radius: float) -> float:
    off = int(round(radius/DX))
    return float(abs(field[CENTER,CENTER,CENTER+off]))


def logfit(xs, ys) -> dict:
    x=np.asarray(xs,dtype=np.float64); y=np.asarray(ys,dtype=np.float64)
    m=np.isfinite(x)&np.isfinite(y)&(x>0)&(y>0)
    lx=np.log(x[m]); ly=np.log(y[m])
    A=np.column_stack((lx,np.ones_like(lx)))
    beta,*_=np.linalg.lstsq(A,ly,rcond=None)
    pred=A@beta
    ss_res=float(np.sum((ly-pred)**2)); ss_tot=float(np.sum((ly-np.mean(ly))**2))
    return {"slope":float(beta[0]),"r2":1.0-ss_res/ss_tot if ss_tot>0 else float("nan"),"count":int(np.count_nonzero(m))}


def match(v,target):
    return bool(math.isfinite(v) and abs(v-target)<=EXP_WINDOW)


def relative_rms(a,b):
    num=math.sqrt(float(np.mean((a-b)**2))); den=math.sqrt(float(np.mean(b*b)))
    return num/den if den>0 else float("nan")


def density_test():
    rows=[]
    for d in DENSITY_LADDER:
        rho=native_rho_sphere(FIXED_RADIUS,d)
        sol=accumulated_from_rho(rho)
        rows.append({"density":d,"rho_integral":sol["rho_integral"],"c_state_integral":sol["c_state_integral"],"center_response":float(abs(sol["field"][CENTER,CENTER,CENTER])),"converged":sol["converged"]})
    return {"rows":rows,"fit":logfit([r["density"] for r in rows],[r["center_response"] for r in rows])}


def mass_test():
    rows=[]
    for m in MASS_LADDER:
        rho=native_rho_fixed_mass(FIXED_RADIUS,m)
        sol=accumulated_from_rho(rho)
        rows.append({"rho_mass":m,"c_state_integral":sol["c_state_integral"],"surface_response":shell_mean(sol["field"],FIXED_RADIUS),"far_response":probe(sol["field"],10.0),"converged":sol["converged"]})
    return {"rows":rows,"surface_fit":logfit([r["rho_mass"] for r in rows],[r["surface_response"] for r in rows]),"far_fit":logfit([r["rho_mass"] for r in rows],[r["far_response"] for r in rows])}


def radius_test():
    rows=[]
    for radius in RADIUS_LADDER:
        rho=native_rho_fixed_mass(radius,FIXED_NATIVE_MASS)
        sol=accumulated_from_rho(rho)
        rows.append({"radius":radius,"rho_mass":sol["rho_integral"],"c_state_integral":sol["c_state_integral"],"surface_response":shell_mean(sol["field"],radius),"converged":sol["converged"]})
    return {"rows":rows,"fit":logfit([r["radius"] for r in rows],[r["surface_response"] for r in rows])}


def far_radius_test():
    rho=native_rho_fixed_mass(FAR_SOURCE_RADIUS,FIXED_NATIVE_MASS)
    sol=accumulated_from_rho(rho)
    vals=[probe(sol["field"],r) for r in FAR_PROBES]
    return {"probe_radii":FAR_PROBES,"responses":vals,"fit":logfit(FAR_PROBES,vals),"converged":sol["converged"],"max_strain_fraction":sol["max_strain_fraction"]}


def shifted_native_rho(offset_x: int, radius: float, mass: float) -> np.ndarray:
    rr=np.sqrt((XN-offset_x)**2+YN*YN+ZN*ZN)*DX
    mask=rr<=radius
    out=np.zeros((N_NATIVE,N_NATIVE,N_NATIVE),dtype=np.float64)
    out[mask]=mass/(np.count_nonzero(mask)*DX**3)
    return out


def additivity_test():
    r1=shifted_native_rho(-5,2.5,0.8)
    r2=shifted_native_rho(+5,3.5,1.2)
    u1=accumulated_from_rho(r1)["field"]
    u2=accumulated_from_rho(r2)["field"]
    u12=accumulated_from_rho(r1+r2)["field"]
    return {"relative_rms_residual":relative_rms(u12,u1+u2)}


def strong_test():
    rho=native_rho_fixed_mass(3.5,STRONG_NATIVE_MASS)
    sol=accumulated_from_rho(rho)
    return {"rho_mass":STRONG_NATIVE_MASS,"max_strain_fraction":sol["max_strain_fraction"],"converged":sol["converged"],"nonlinear_relative_residual":sol["nonlinear_relative_residual"],"history_tail":sol["history_tail"]}


def main() -> int:
    density=density_test(); mass=mass_test(); radius=radius_test(); far=far_radius_test(); add=additivity_test(); strong=strong_test(); state=repo_state()
    measured={
        "density_exponent":density["fit"]["slope"],
        "surface_mass_exponent":mass["surface_fit"]["slope"],
        "surface_radius_fixed_mass_exponent":radius["fit"]["slope"],
        "far_mass_exponent":mass["far_fit"]["slope"],
        "far_radius_exponent":far["fit"]["slope"],
        "weak_additivity_relative_residual":add["relative_rms_residual"],
        "strong_load_max_strain_fraction":strong["max_strain_fraction"],
    }
    checks={
        "density_linearity":match(measured["density_exponent"],1.0),
        "surface_mass_linearity":match(measured["surface_mass_exponent"],1.0),
        "surface_radius_accumulation":match(measured["surface_radius_fixed_mass_exponent"],-1.0),
        "far_mass_linearity":match(measured["far_mass_exponent"],1.0),
        "far_radius_accumulation":match(measured["far_radius_exponent"],-1.0),
        "weak_regime_additivity":bool(math.isfinite(measured["weak_additivity_relative_residual"]) and measured["weak_additivity_relative_residual"]<=ADDITIVITY_TOL),
    }
    matched=sum(bool(v) for v in checks.values())
    nonlinear_converged=all(r["converged"] for r in density["rows"]+mass["rows"]+radius["rows"]) and far["converged"] and strong["converged"]
    execution_checks={
        "all_measured_values_finite":all(math.isfinite(v) for v in measured.values()),
        "nonlinear_solves_converged":nonlinear_converged,
        "raw_c_state_used_without_amplitude_rescaling":True,
        "K0_frozen":K0==1.0,
        "epsilon_max_frozen":EPSILON_MAX==1.0,
        "no_G":True,"no_macroscopic_amplitude":True,"no_native_rescaling":True,"no_fit_or_tuning":True,
        "no_inserted_one_over_r_response":True,"no_spherical_equilibrium_shortcut":True,"no_Rmax":True,
        "no_cosmology":True,"no_lensing_target":True,"no_quantum_engine":True,"no_planck_input":True,
        "no_tracked_or_staged_changes":state["tracked_changes"]=="" and state["staged_changes"]=="",
        "stdout_only_no_run_directory_created":True,
    }
    execution_gate_pass=all(execution_checks.values())
    if matched==6:
        status="C_STATE_BOUNDED_STRAIN_BRIDGE_STRUCTURE_SUPPORTED"
    elif matched>=4:
        status="C_STATE_BOUNDED_STRAIN_BRIDGE_PARTIAL_SUPPORT"
    else:
        status="C_STATE_BOUNDED_STRAIN_BRIDGE_NOT_SUPPORTED"
    payload={
        "lab_id":LAB_ID,"status":status,"repo_state":state,
        "model":{"native_source":"rho -> existing A8 transport -> raw c_state","accumulation":"raw c_state -> six-neighbor bounded-strain equilibrium","K0":K0,"epsilon_max":EPSILON_MAX,"native_grid_N":N_NATIVE,"accumulation_grid_N":N,"amplitude_rescaling":False},
        "density_test":density,"mass_test":mass,"radius_fixed_mass_test":radius,"far_radius_test":far,"additivity_test":add,"strong_load_test":strong,
        "measured":measured,"checks":checks,"matched_frozen_checks_of_6":matched,
        "execution_checks":execution_checks,"execution_gate_pass":execution_gate_pass,
        "summary":{"question":"Can the existing raw c_state local-loading candidate drive the bounded-strain 3D network and reproduce the frozen long-range response fingerprint without amplitude fitting?","next_if_supported":"Promote the c_state -> bounded-strain network path as the leading native accumulation bridge candidate and audit absolute constitutive scale separately.","next_if_partial":"Localize the remaining fingerprint failure before altering c_state or the constitutive law."}
    }
    print(LAB_ID)
    print(f"status={status}")
    print(f"head_sha={state['head_sha']}")
    print("raw_c_state_used=true")
    print("native_amplitude_rescaled=false")
    print("G_used=false")
    print("fit_or_tuning_used=false")
    print("one_over_r_response_inserted=false")
    print("\nFROZEN_RESPONSE_FINGERPRINT")
    for k,v in measured.items(): print(f"{k}={v:.12g}")
    print(f"matched_frozen_checks_of_6={matched}")
    print("\nCHECKS")
    for k,v in checks.items(): print(f"{k}={str(v).lower()}")
    print("\nEXECUTION_CHECKS")
    for k,v in execution_checks.items(): print(f"{k}={str(v).lower()}")
    print(f"execution_gate_pass={str(execution_gate_pass).lower()}")
    print(json.dumps(payload,sort_keys=True,separators=(",",":")))
    return 0 if execution_gate_pass else 1


if __name__ == "__main__":
    raise SystemExit(main())
