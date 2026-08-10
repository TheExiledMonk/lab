"""Dev145 passive loaded-propagation ontology; contains no SR construction."""
from __future__ import annotations
from dataclasses import asdict, dataclass
from typing import Any
import numpy as np

Q0_LEVELS=(0.25,0.5,1.0,2.0,4.0,8.0)
LOADING_LEVELS=(0.0,0.1,0.25,0.5,0.75,0.9)

@dataclass(frozen=True)
class LoadedPropagationState:
    rest_loading_state: Any
    excitation_state: float
    propagation_direction: tuple[float, ...]
    propagation_fraction_beta: float
    zero_mass_limit: bool = False
    super_c_allowed: bool = False
    time_fundamental: bool = False
    def __post_init__(self):
        if not np.isfinite(self.excitation_state) or not 0 <= self.propagation_fraction_beta <= 1:
            raise ValueError("finite excitation and 0 <= beta <= 1 required")
        if self.super_c_allowed or self.time_fundamental: raise ValueError("forbidden Dev145 state")
        if self.zero_mass_limit and self.propagation_fraction_beta != 1: raise ValueError("zero loading requires beta=1")
    def contract(self): return {"state":"PBUF_MASS_LOADED_PROPAGATION_STATE_V1",**asdict(self)}

P_NAMES=("local loading","accumulated loading","strain","stress","elastic deformation energy",
 "normalized loading ratio","loading / excitation","excitable fraction","unloaded medium fraction",
 "local constitutive stiffness state","loading compactness","integrated deformation","response fingerprint",
 "bound-strain proximity","conservation relation","orthogonal loading/excitation decomposition",
 "energy partition","momentum-loading relation","invariant norm","no existing loading variable sufficient")
R_NAMES=("local loading resistance","integrated loading resistance","strain-fraction resistance",
 "stress-state resistance","deformation-energy resistance","loading/excitation ratio",
 "excitation-minus-loading capacity","normalized remaining-capacity","additive state partition",
 "quadratic state partition","orthogonal norm partition","bounded-state geometric partition",
 "local stiffness modulation","compactness-dependent propagation","accumulated-response-dependent propagation",
 "excitation-to-loading ratio","conserved total-state norm","conserved pair-state norm",
 "rest-loading invariant + variable excitation","no derivable relation in current PBUF")

def candidate_manifest(prefix: str):
    names=P_NAMES if prefix=="P" else R_NAMES
    return [{"id":f"{prefix}{i:02d}","name":n,"attempted":True} for i,n in enumerate(names,1)]

def mechanism_results():
    out=[]
    for i,n in enumerate(R_NAMES,1):
        if i in (1,2,3,4,5,13,14,15): status="MISSING_CONSTITUTIVE_LAW"
        elif i in (6,7,8,9,10,12,16,19): status="RELATION_ONLY"
        elif i in (11,17,18): status="MISSING_EXCITATION_DEFINITION"
        else: status="ESTABLISHED"
        out.append({"id":f"R{i:02d}","name":n,"status":status,"native_beta_law":False,
                    "same_mass_variable_speed":"REPRESENTABLE" if i in (6,7,8,9,10,11,12,16,17,18,19) else False})
    return out

def loading_only_audit():
    return {"status":"FAILS_VARIABLE_SPEED_SAME_MASS","reason":"rest loading alone fixes at most one beta",
            "beta_law_accepted":False}

def diagnostic_surface():
    """Unresolved surface: NaN denotes no native prediction, not missing numerics."""
    l=np.asarray(LOADING_LEVELS); q=np.asarray(Q0_LEVELS); b=np.full((l.size,q.size),np.nan); b[0,:]=1
    return l,q,b

def propagation_fraction_contract():
    return {"symbol":"beta_PBUF","definition":"fraction of maximum available medium propagation/change rate",
            "range":[0,1],"SI_velocity_created":False,"time_required":False,"zero_loading_beta":1,
            "loaded_beta_law":"UNRESOLVED","super_c_allowed":False}

def loaded_speed_contract():
    return {"contract":"PBUF_LOADED_PROPAGATION_SPEED_V1","beta_law_established":False,"beta_law":None,
            "loading_input":"persistent native loading proxy","excitation_input":"neutral q control (not identified as energy)",
            "zero_load_beta":1.0,"beta_upper_bound":1.0,"same_mass_variable_speed_supported":"STATE_ONTOLOGY_ONLY",
            "high_excitation_limit":"UNRESOLVED","low_excitation_limit":"REPRESENTABLE_NOT_DERIVED",
            "free_parameters":0,"SR_used_in_derivation":False,"SR_benchmark_compatible":False,
            "time_required":False,"L0_required":False}

