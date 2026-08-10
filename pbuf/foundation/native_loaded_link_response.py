"""Dev153 candidate registry and coefficient-free response execution."""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from .native_transverse_transfer_capacity import orthogonal_transport

CANDIDATE_NAMES = ["actual neighbor separation", "normalized longitudinal strain", "native stress/tension",
 "tangent stiffness", "bounded constitutive curvature", "local constitutive energy density", "geometric stretch ratio",
 "normalized link-length ratio", "local longitudinal-state gradient", "local stress gradient",
 "local tangent-stiffness gradient", "link-state curvature along propagation", "pairwise loading contrast",
 "local loading anisotropy", "longitudinal equilibrium offset", "bond-response derivative",
 "local state-space capacity", "joint link norm constraint", "geometry + constitutive-state combination",
 "no native longitudinal dependence"]
RESPONSE_NAMES = ["orientation only", "progression allocation", "reversible amplitude redistribution",
 "local wavelength adjustment", "mode redistribution", "polarization rotation", "handedness rotation",
 "neighbor-transfer preference", "state-space capacity modification", "reversible excitation/loading exchange",
 "transverse-to-longitudinal conversion", "no response"]


def candidate_registry():
    return [{"id": f"T{i:02d}", "name": name, "gradient_only": 9 <= i <= 14} for i, name in enumerate(CANDIDATE_NAMES, 1)]


def response_registry(): return [{"id": f"RSP{i:02d}", "name": n} for i,n in enumerate(RESPONSE_NAMES, 1)]


def load_profile(load_id, n=64, amplitude=None):
    x = np.linspace(-1, 1, n); idx = int(load_id[-2:]); a = ([0,.05,.25,.75,.05,.25,.25,.25,.25,.25,.25][idx] if amplitude is None else amplitude)
    if idx == 0: p=np.zeros(n)
    elif idx <= 3: p=np.full(n,a)
    elif idx in (4,5): p=a*np.exp(-(x/.22)**2)
    elif idx == 6: p=a*np.exp(-(x/.55)**2)
    elif idx == 7: p=a*(x+1)/2
    elif idx == 8: p=a*np.exp(-((x-.25)/.3)**2)
    elif idx == 9: p=a*(np.exp(-((x-.4)/.2)**2)+np.exp(-((x+.4)/.2)**2))
    else: p=a*np.exp(-((np.abs(x)-.5)/.12)**2)
    return x,p


def excitation(ex_id, n=64, amplitude=1.0, wavelength=12.0):
    x=np.arange(n); k=2*np.pi/wavelength; phase=k*(x-n/3); env=np.exp(-((x-n/3)/(n/9))**2); i=int(ex_id[-2:])
    if i==1: return amplitude*np.stack([env*np.cos(phase),np.zeros(n)],1)
    if i==2: return amplitude*np.stack([np.zeros(n),env*np.cos(phase)],1)
    if i==3: return amplitude*np.stack([env*np.cos(phase),env*np.cos(phase)],1)/np.sqrt(2)
    if i in (4,5): return amplitude*np.stack([env*np.cos(phase),(1 if i==4 else -1)*env*np.sin(phase)],1)/np.sqrt(2)
    width=n/(14 if i==6 else 6); env=np.exp(-((x-n/3)/width)**2)
    if i==8: env += np.exp(-((x-2*n/3)/width)**2)
    carrier=np.cos(phase)+(0.35*np.cos(2*phase) if i==9 else 0)
    return amplitude*np.stack([env*carrier,np.zeros(n)],1)


def execute(candidate_id, load_id, excitation_id, n=64, amplitude=1.0, wavelength=12.0):
    _, load=load_profile(load_id,n); x=excitation(excitation_id,n,amplitude,wavelength); R=orthogonal_transport(0.17)
    # The frozen Dev151/152 energy contains no L-X mixed Hessian or derived
    # transverse metric. Consequently no T01-T19 defines a transfer operator.
    out=x @ R.T
    derived = candidate_id == "T20"
    local_claim = int(candidate_id[1:]) <= 8 or 15 <= int(candidate_id[1:]) <= 19
    gradient_claim = 9 <= int(candidate_id[1:]) <= 14
    return {"candidate_id":candidate_id,"load_id":load_id,"excitation_id":excitation_id,
      "derived_transfer":derived,"classification":"NO_EFFECT" if derived else "UNDERDETERMINED",
      "dependence":"NEITHER","claimed_dependence":"LOCAL_STATE_DEPENDENT" if local_claim else ("GRADIENT_DEPENDENT" if gradient_claim else "NEITHER"),
      "transfer_fraction":1.0,"progression_ratio":1.0,"norm_in":float(np.sum(x*x)),"norm_out":float(np.sum(out*out)),
      "polarization_rotation":0.17,"handedness_change":0.0,"longitudinal_leakage":0.0,"load_mean":float(np.mean(load))}
