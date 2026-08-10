"""Read-only static-medium exposure and coefficient-free coupling gate."""
from __future__ import annotations
import numpy as np

def expose_static_medium(medium):
    x=np.asarray(medium,float).copy(); x.flags.writeable=False; return x
def coupling_audit(medium):
    frozen=expose_static_medium(medium)
    return {"static_medium_available":bool(frozen.size),"coefficient_free_coupling_available":False,
      "coupling_definition":None,"classification":"STATIC_MEDIUM_EXCITATION_COUPLING_UNRESOLVED",
      "static_medium_modified":False,"trajectory_solver_dependency":False}
def coupling_contract():
    return {"contract":"PBUF_DYNAMIC_EXCITATION_STATIC_MEDIUM_COUPLING_V1",
      "dynamic_excitation_available":True,"static_medium_available":True,
      "coefficient_free_coupling_available":False,"coupling_definition":None,
      "unloaded_propagation":True,"nonuniform_path_emergence":False,
      "frozen_ray_path_comparison_available":False,"trajectory_parity":False,
      "mass_loading_coupling_preliminary":False,"loaded_progression_available":False,
      "zero_mass_propagation_changed":False}
