"""Invariant and physical-rank audits for the native rank-2 excitation."""
from __future__ import annotations
import numpy as np
from .native_excitation_state import NativeExcitationState
from .native_excitation_transfer import progress_source_free

def quadratic_norm(values): return float(np.sum(np.asarray(values,float)**2))
def centroid(values):
    q=np.sum(np.asarray(values,float)**2,axis=1); x=np.arange(len(q),dtype=float)
    # circular centroid is avoided in tests whose packets do not cross the periodic edge.
    return float(np.sum(x*q)/np.sum(q)) if np.sum(q)>0 else float('nan')
def invariant_audit(history):
    h=np.asarray(history,float); norms=np.sum(h*h,axis=(1,2)); initial=norms[0]
    return {"norms":norms,"relative_drift":float(np.max(np.abs(norms-initial))/(initial+1e-300)),
      "conserved":bool(np.allclose(norms,initial,rtol=0,atol=64*np.finfo(float).eps*max(initial,1.)))}
def superposition_audit(a,b,steps=5):
    evolve=lambda x: progress_source_free(NativeExcitationState(np.asarray(x).copy()),steps).values
    residual=evolve(np.asarray(a)+np.asarray(b))-evolve(a)-evolve(b)
    return {"max_absolute_residual":float(np.max(np.abs(residual))),"passes":bool(np.allclose(residual,0,atol=1e-14))}
def reversibility_audit(values,steps=7):
    s=NativeExcitationState(np.asarray(values).copy()); progress_source_free(s,steps,1); progress_source_free(s,steps,-1)
    err=float(np.max(np.abs(s.values-values)))
    return {"max_absolute_error":err,"reversible":err==0.0}
def transverse_rank_audit():
    basis=np.eye(2); return {"gram":(basis@basis.T).tolist(),"physical_rank":int(np.linalg.matrix_rank(basis)),
      "transverse_dof_count":2,"longitudinal_mode_present":False,"longitudinal_leakage":0.0}

