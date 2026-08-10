"""Dev150 passive localized loading/excitation state diagnostics.

Construction is deliberately not dynamics: no trapping potential or loading to
excitation coupling is introduced here.
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from pbuf.quantum.native_excitation_modes import estimate_wavelengths, native_k, quadratic_norm

STATE_NAMES = ("loading only","weak internal excitation","moderate internal excitation","strong internal excitation",
 "compact excitation","diffuse excitation","symmetric excitation","asymmetric excitation","single transverse mode",
 "orthogonal transverse mode","mixed transverse state","rotating transverse state","standing internal excitation",
 "circulating internal excitation","node-centered excitation","shell-like excitation","multi-lobed excitation",
 "loading-bound spatial mode","naturally stable localized survivor","no stable localized excitation")
LOADING_NAMES = ("weak compact","moderate compact","strong compact unsaturated","weak diffuse","strong diffuse",
                 "asymmetric","shell-like","two-center")
TOPOLOGY_NAMES = ("no topology","sign-domain topology","transverse orientation winding","circulation winding",
 "shell topology","node parity","link circulation","closed-loop state","mixed loading/excitation topology","no stable invariant")

def state_registry():
    return [{"id":f"S{i:02d}","name":n,"attempted":True,
             "status":"MISSING_COUPLING" if i not in (1,20) else "ESTABLISHED"}
            for i,n in enumerate(STATE_NAMES,1)]

def loading_registry():
    return [{"id":f"LOAD{i:02d}","name":n,"proxy":["L06 |epsilon|","L08 W(epsilon)","L16 volume-integrated loading"],
             "synthetic":True} for i,n in enumerate(LOADING_NAMES,1)]

def topology_registry():
    return [{"id":f"T{i:02d}","name":n,"attempted":True,"status":"MISSING_LOCALIZED_STATE"}
            for i,n in enumerate(TOPOLOGY_NAMES,1)]

@dataclass(frozen=True)
class LocalizedComposite:
    loading: np.ndarray
    excitation: np.ndarray
    def __post_init__(self):
        l=np.asarray(self.loading,float); x=np.asarray(self.excitation,float)
        if l.ndim != 1 or x.shape != (len(l),2) or not np.isfinite(l).all() or not np.isfinite(x).all():
            raise ValueError("requires finite loading[n] and excitation[n,2]")
        object.__setattr__(self,"loading",l); object.__setattr__(self,"excitation",x)
    @property
    def internal_norm(self): return quadratic_norm(self.excitation)

def construct_composite(loading, excitation): return LocalizedComposite(loading,excitation)

def state_observables(state: LocalizedComposite):
    amp=np.linalg.norm(state.excitation,axis=1); support=amp > max(float(amp.max(initial=0))*1e-3,1e-15)
    est=estimate_wavelengths(state.excitation) if np.any(amp) else {k:None for k in ("L01","L02","L03","L04","L05","L08")}
    lam=est.get("L04")
    return {"total_internal_excitation_norm":state.internal_norm,"localized_excitation_norm":None,
            "escaped_excitation_norm":None,"spatial_support":int(support.sum()),"internal_wavelength":lam,
            "internal_k_n":native_k(lam) if lam else None,"loading_integral":float(np.sum(np.abs(state.loading))),
            "classification":"SYNTHETIC_OVERLAP_NOT_BOUND"}

def stability_analysis(history, loading_history=None):
    h=np.asarray(history,float)
    finite=bool(h.ndim==4 and h.shape[-1]==2 and np.isfinite(h).all())
    return {"persistent_loading_preserved":None if loading_history is None else bool(np.allclose(loading_history,loading_history[0])),
            "internal_excitation_remains_localized":False,"native_excitation_norm_finite":finite,
            "source_free":True,"reproducible":True,"resolution_survivor":False,"perturbation_survivor":False,
            "classification":"FREE" if finite else "UNSTABLE","stable_localized_state":False}
