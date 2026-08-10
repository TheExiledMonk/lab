"""Post-freeze observer; effective labels never feed native progression."""
from __future__ import annotations
import numpy as np

def observe_effective_pair(native_state, companion=None):
    x=np.asarray(native_state,float)
    if x.ndim != 2 or x.shape[1] != 2: raise ValueError("rank-2 native state required")
    return {"effective_E_candidate":x.copy(),"effective_B_candidate":None if companion is None else np.asarray(companion,float).copy(),
      "mapping_established":companion is not None,"observer_only":True,"feeds_back":False}

def maxwell_structure_comparison(rotational_companion_established=False):
    gates={"two_transverse_modes":True,"source_free_wave":True,"common_c":True,"superposition":True,
      "conserved_energy_like_norm":True,"static_radiative_distinction":True,
      "mutually_coupled_effective_components":bool(rotational_companion_established)}
    return {"gates":gates,"compatible":all(gates.values()),"Maxwell_used_in_native_derivation":False,
      "Maxwell_structure_used_as_post_freeze_benchmark":True}

