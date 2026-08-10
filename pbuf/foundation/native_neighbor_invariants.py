"""Invariants and basis audits for unified neighbor states."""
from __future__ import annotations
import numpy as np
from .native_neighbor_state import rotate_transverse
from .native_neighbor_constitutive_law import pair_energy

def excitation_norm(state): return float(np.sum(state.transverse**2))
def basis_invariance(values,angles=(.1,.7,1.9),atol=1e-12):
    x=np.asarray(values,float); n=float(np.sum(x*x)); errors=[abs(float(np.sum(rotate_transverse(x,a)**2))-n) for a in angles]
    return {"status":"ESTABLISHED" if max(errors,default=0)<=atol*max(1,n) else "BASIS_DEPENDENT","max_error":max(errors,default=0.)}
def joint_invariant_audit(state,K=1.,epsilon_max=1.):
    n=excitation_norm(state); u=float(np.sum(pair_energy(state.longitudinal,state.transverse,K,epsilon_max)))
    return {"J01":"UNDERDETERMINED","J02":"UNDERDETERMINED","J03":{"value":u,"status":"STRUCTURALLY_SUPPORTED"},
            "J04":"UNDERDETERMINED","J05":"UNDERDETERMINED","J06":{"value":n+float(np.sum(state.longitudinal**2)),"status":"STRUCTURALLY_SUPPORTED"},
            "J07":"UNDERDETERMINED","J08":"UNDERDETERMINED","J09":"UNDERDETERMINED","J10":"UNDERDETERMINED"}
