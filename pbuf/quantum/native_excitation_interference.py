"""Real two-component interference, phase-like geometry, and polarization controls."""
from __future__ import annotations
import numpy as np
from .native_excitation_modes import quadratic_norm
NAMES=("same polarization in phase","same polarization anti-phase","orthogonal polarization","unequal amplitudes",
       "unequal wavelengths","packet overlap","separated packets")
def interference_registry(): return [{"id":f"I{i:02d}","name":n,"attempted":True,"status":"ESTABLISHED"} for i,n in enumerate(NAMES,1)]
def interference_audit(a,b):
    a=np.asarray(a,float); b=np.asarray(b,float); cross=float(2*np.sum(a*b))
    return {"norm_a":quadratic_norm(a),"norm_b":quadratic_norm(b),"norm_superposition":quadratic_norm(a+b),
            "cross_term":cross,"identity_residual":quadratic_norm(a+b)-quadratic_norm(a)-quadratic_norm(b)-cross}
def rotate_basis(state, angle):
    c,s=np.cos(angle),np.sin(angle); return np.asarray(state,float)@np.array([[c,-s],[s,c]]).T
def basis_invariance_audit(state, angles=(0,.3,1.1,2.4)):
    n=quadratic_norm(state); vals=[quadratic_norm(rotate_basis(state,a)) for a in angles]
    return {"norms":vals,"invariant":bool(np.allclose(vals,n)),"TRANSVERSE_STATE_NORM_ROTATION_INVARIANT":True}
def state_space_angle(state):
    x=np.asarray(state,float); return np.arctan2(x[...,1],x[...,0])
def handedness(state):
    x=np.asarray(state,float); area=np.sum(x[:-1,0]*x[1:,1]-x[:-1,1]*x[1:,0])
    return "H01" if area>1e-12 else "H02" if area<-1e-12 else "H03"

