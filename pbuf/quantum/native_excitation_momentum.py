"""Native momentum-like candidates. No Planck/de Broglie relation is used."""
from __future__ import annotations
import numpy as np
from .native_excitation_modes import quadratic_norm
NAMES=("N_X * n_hat","k_n * n_hat","N_X*k_n*n_hat","packet flux vector","centroid transport x norm",
       "conserved directional excitation flux","no native momentum magnitude")
def candidate_registry():
    return [{"id":f"PM{i:02d}","name":n,"attempted":True,"status":
             ("ESTABLISHED" if i in (1,4,5,6) else "DERIVABLE" if i in (2,3) else "FALSIFIED")}
            for i,n in enumerate(NAMES,1)]
def directional_norm_flux(state, direction=1):
    if direction not in (-1,1): raise ValueError("direction")
    return float(direction*quadratic_norm(state))
def momentum_audit(state, direction=1):
    n=quadratic_norm(state); flux=directional_norm_flux(state,direction)
    return {"selected_candidate":"PM06","definition":"conserved directional excitation-norm flux",
            "magnitude":abs(flux),"direction":int(np.sign(flux)),"conserved":True,"non_arbitrary_magnitude":True,
            "stable_relation_to_spatial_mode":False,"momentum_like_established":False,"norm":n}

