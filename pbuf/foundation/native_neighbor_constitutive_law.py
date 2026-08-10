"""Candidate unified local constitutive laws for Dev151."""
from __future__ import annotations
import numpy as np
from pbuf.wl.native_incremental_elastic_energy import bounded_strain_energy,bounded_strain_stress,bounded_strain_tangent

LAW_NAMES=("magnitude-only bond energy","longitudinal + transverse quadratic norm","bounded longitudinal + linear transverse",
"bounded full-vector bond deformation","separation/orientation energy","symmetric/antisymmetric state energy",
"tension-derived transverse stiffness","local tangent-stiffness state","bond-length-dependent transverse transfer",
"geometry-only neighbor interaction","conserved pair-state norm","constrained state-manifold law","orientation-compatibility law",
"reciprocal push/pull law","finite-capacity neighbor-state law","nonlinear geometric law","pre-tensioned vector-bond law",
"shared state with separate static/dynamic limits","underdetermined family","no coefficient-free unified law")
MEC_NAMES=("longitudinal equilibrium + transverse perturbation","vector bond displacement","pre-tensioned transverse wave",
"strain-dependent transverse stiffness","tangent-stiffness propagation","orientation-transport propagation","link-frame rotation",
"pair-state conserved norm","bounded vector-bond constitutive law","symmetric/antisymmetric decomposition",
"full displacement-vector elasticity","local bond-energy curvature","finite-capacity state manifold","state-space orthogonality",
"geometric constraint-only coupling","constitutive constraint-only coupling","geometry + constitutive coupling",
"unified law reproduces sectors but no cross-coupling","multiple equivalent unified laws","no tested unified law")

def _registry(names,prefix,survivors):
    return [{"id":f"{prefix}{i:02d}","name":n,"attempted":True,
             "status":"STRUCTURALLY_SUPPORTED" if i in survivors else ("UNDERDETERMINED" if i<20 else "FALSIFIED")}
            for i,n in enumerate(names,1)]
def law_registry(): return _registry(LAW_NAMES,"C",{8,10,12,13,16,18})
def mechanism_registry(): return _registry(MEC_NAMES,"MEC",{1,2,6,7,12,15,17,18})

def pair_energy(longitudinal,transverse,K=1.,epsilon_max=1.):
    """One expression on the full state; transverse curvature inherits K."""
    e=np.asarray(longitudinal,float); x=np.asarray(transverse,float)
    if x.shape != e.shape+(2,): raise ValueError("state shape mismatch")
    return bounded_strain_energy(e,K,epsilon_max)+.5*K*np.sum(x*x,axis=-1)

def pair_response(longitudinal,transverse,K=1.,epsilon_max=1.):
    e=np.asarray(longitudinal,float); x=np.asarray(transverse,float)
    return {"longitudinal":bounded_strain_stress(e,K,epsilon_max),"transverse":K*x,
            "tangent":bounded_strain_tangent(e,K,epsilon_max)}

def coefficient_inventory():
    return {"K":"frozen bounded-strain parameter","epsilon_max":"frozen bounded-strain parameter",
            "new_interaction_coefficients":0,"post_hoc_loading_excitation_coefficients":0}
