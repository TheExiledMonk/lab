"""Post-freeze external SR fingerprints.  Native modules must not import this file."""
from __future__ import annotations
import numpy as np

BENCHMARK_BETA=np.array([0,.1,.25,.5,.75,.9,.99,.999,.9999],dtype=float)

def gamma_from_beta(beta):
    b=np.asarray(beta,dtype=float)
    if np.any(~np.isfinite(b)) or np.any(np.abs(b)>=1): raise ValueError("gamma requires finite abs(beta) < 1")
    return 1/np.sqrt(1-b*b)

def required_loading(beta):
    b=np.asarray(beta,dtype=float)
    if np.any(~np.isfinite(b)) or np.any(np.abs(b)>1): raise ValueError("requires abs(beta) <= 1")
    return np.sqrt(1-b*b)

def orthogonal_benchmark_mapping(native_loading):
    ell=np.asarray(native_loading,dtype=float)
    if np.any(~np.isfinite(ell)) or np.any((ell<0)|(ell>1)): raise ValueError("requires 0 <= loading <= 1")
    return np.sqrt(1-ell*ell)

def benchmark_contract():
    return {"lane":"POST_FREEZE_EXTERNAL_BENCHMARK","SR_USED_TO_CONSTRUCT_PBUF_LAW":False,
            "SR_USED_AS_POST_FREEZE_BENCHMARK":True,"mapping_label":"ORTHOGONAL_BENCHMARK_MAPPING",
            "native_orthogonal_derivation":False,"outcome":"ORTHOGONAL_LOADING_RELATION_UNJUSTIFIED"}
