"""Coefficient-free Dev150 loading/excitation interaction audit."""
from __future__ import annotations
import numpy as np

INTERACTION_NAMES=("local loading modifies excitation transfer geometry","local loading modifies excitation orientation",
 "loading modifies local excitation support","loading modifies link-transfer admissibility","loading binds standing excitation mode",
 "loading binds circulating excitation mode","loading/excitation norm exchange","loading-gradient interaction",
 "strain-gradient interaction","deformation-energy interaction","bounded-state interaction","mode-conversion interaction",
 "transverse-mode mixing","pair/link interaction","local excitation threshold interaction",
 "nonlinear excitation self-consistency","loading-induced boundary condition","topology-induced localization",
 "no coefficient-free interaction","current state insufficient")

def interaction_registry():
    return [{"id":f"I{i:02d}","name":n,"attempted":True,
             "status":"ESTABLISHED" if i in (19,20) else "MISSING_INTERACTION_LAW",
             "post_hoc_coefficients":0} for i,n in enumerate(INTERACTION_NAMES,1)]

def coefficient_free_coupling_audit():
    return {"families":interaction_registry(),"coefficient_free_binding_law_found":False,
            "classification":"MISSING_COUPLING","outcome":"LOADING_EXCITATION_BINDING_LAW_UNRESOLVED",
            "arbitrary_trapping_potential_used":False,"post_hoc_interaction_coefficients":0}

def norm_exchange(bound_before, bound_after, free_norm, *, atol=1e-12):
    a,b,f=map(float,(bound_before,bound_after,free_norm)); delta=a-b
    denom=max(abs(delta),np.finfo(float).eps); error=abs(delta-f)/denom
    return {"bound_norm_before":a,"bound_norm_after":b,"delta_bound":delta,"free_norm":f,
            "relative_error":error,"conserved":bool(error<=atol)}

def interaction_locality():
    return {"local_loading_region":False,"boundary_shell":False,"distributed_overlap":False,
            "link_state_interaction":False,"instantaneous_global_update":False,"status":"MISSING_INTERACTION_LAW"}

def backreaction_audit():
    return {"classification":"UNRESOLVED","internal_excitation_modifies_effective_loading":False,
            "excitation_bound_identified_with_strain_bound":False}
