"""Coefficient-free local progression laws for native excitation."""
from __future__ import annotations
import numpy as np
from .native_excitation_state import NativeExcitationState

OPERATORS=("nearest-neighbor difference","central difference","divergence-like node sum",
 "oriented link circulation","plaquette curl-like operator","transverse neighbor rotation",
 "antisymmetric pair exchange","six-neighbor vector Laplacian","coupled divergence/circulation",
 "existing medium N6 topology only")
LAWS=("scalar neighbor exchange","signed scalar wave exchange","vector neighbor exchange",
 "transverse vector exchange","coupled scalar-pair exchange","coupled vector-pair exchange",
 "node-link conjugate exchange","divergence/circulation conjugate exchange","oriented-link rotation exchange",
 "antisymmetric pair exchange","six-mode coupled exchange","strain-excitation coupled exchange",
 "excitation/circulation coupled exchange","conserved-norm transport","bounded excitation transport",
 "source-free reciprocal exchange","magnetism-inspired rotational exchange","EM-rank-minimal exchange",
 "structurally valid but underdetermined","no coefficient-free law")

def operator_registry(): return [{"id":f"O{i:02d}","name":n,"attempted":True} for i,n in enumerate(OPERATORS,1)]
def dynamic_law_registry():
    return [{"id":f"D{i:02d}","name":n,"attempted":True,
      "status":"ESTABLISHED" if i in (4,14,16,18) else ("STRUCTURALLY_SUPPORTED" if i in (6,7,8,13,17,19) else "REQUIRES_FREE_COEFFICIENT")}
      for i,n in enumerate(LAWS,1)]

def nearest_neighbor_shift(values, direction=1):
    """Exact local permutation on a periodic native lattice; no fitted coefficient."""
    x=np.asarray(values,float)
    if x.ndim != 2 or x.shape[1] != 2: raise ValueError("rank-2 transverse state required")
    if direction not in (-1,1): raise ValueError("direction must be -1 or 1")
    return np.roll(x,int(direction),axis=0)

def progress_source_free(state, steps, direction=1):
    """Progress from the prior excitation only; source and trajectory are absent."""
    if not isinstance(state,NativeExcitationState): raise TypeError("NativeExcitationState required")
    for _ in range(int(steps)):
        state.values=nearest_neighbor_shift(state.values,direction)
        state.progression_step+=1; state.history.append(state.values.copy())
    return state

def dependency_contract():
    return {"next_depends_on_previous_excitation":True,"source_present_after_launch":False,
      "trajectory_solver_used_to_move_excitation":False,"static_response_used_as_dynamic_state":False,
      "local_medium_element_translation_required":False,"free_coefficients":0}

