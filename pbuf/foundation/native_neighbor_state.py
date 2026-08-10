"""Native link state used by the Dev151 unified-neighbor audit."""
from __future__ import annotations
from dataclasses import dataclass, field
import numpy as np

NAMES=("separation-only scalar","signed separation state","longitudinal bond state",
"transverse bond displacement state","orientation vector state","longitudinal + transverse bond state",
"two-transverse-component bond state","longitudinal + two-transverse state","scalar strain + transverse pair",
"vector bond deformation state","bond length + orientation","node displacement vector with derived links",
"link displacement vector","local frame/orientation state","symmetric + antisymmetric link decomposition",
"tension + transverse excitation","stretch + shear-like excitation","radial + tangential state",
"minimal state inferred from parity requirements","no unified state sufficient")

def state_registry():
    out=[]
    for i,name in enumerate(NAMES,1):
        rank=3 if i in (8,9,10,12,13,19) else (2 if i in (6,7,11,16,17,18) else 1)
        status="STRUCTURALLY_SUPPORTED" if i in (8,9,10,12,13,19) else ("MISSING_DYNAMIC_PARITY" if rank<3 else "UNDERDETERMINED")
        if i==20: status="FALSIFIED"
        out.append({"id":f"N{i:02d}","name":name,"attempted":True,"total_dof":rank,
                    "static_dof":1 if i not in (4,5,7,14,20) else 0,"dynamic_dof":max(0,rank-1),
                    "longitudinal_dof":1 if rank else 0,"transverse_dof":max(0,rank-1),
                    "gauge_dof":0,"physical_rank":rank,"status":status})
    return out

def local_link_frame(relation):
    """Deterministic right-handed frame; no preferred transverse physics."""
    r=np.asarray(relation,float)
    if r.shape != (3,) or not np.isfinite(r).all() or np.linalg.norm(r)==0: raise ValueError("finite nonzero 3-vector required")
    ep=r/np.linalg.norm(r); seed=np.array([1.,0.,0.]) if abs(ep[0])<.9 else np.array([0.,1.,0.])
    e1=seed-ep*np.dot(seed,ep); e1/=np.linalg.norm(e1); e2=np.cross(ep,e1)
    return np.stack((ep,e1,e2))

def decompose(deformation, frame):
    d=np.asarray(deformation,float); f=np.asarray(frame,float)
    if d.shape != (3,) or f.shape != (3,3): raise ValueError("deformation/frame shapes must be (3,)/(3,3)")
    return float(d@f[0]), np.array([d@f[1],d@f[2]])

def compose(longitudinal, transverse, frame):
    x=np.asarray(transverse,float); f=np.asarray(frame,float)
    if x.shape != (2,) or f.shape != (3,3): raise ValueError("transverse/frame shapes must be (2,)/(3,3)")
    return float(longitudinal)*f[0]+x[0]*f[1]+x[1]*f[2]

def rotate_transverse(values, angle):
    x=np.asarray(values,float); c,s=np.cos(angle),np.sin(angle)
    return x@np.array([[c,s],[-s,c]])

def frame_overlap(source_frame,target_frame):
    """Orthogonal two-mode map induced solely by adjacent link frames."""
    a,b=np.asarray(source_frame,float),np.asarray(target_frame,float)
    raw=b[1:]@a[1:].T
    u,_,vh=np.linalg.svd(raw); return u@vh

@dataclass
class NativeNeighborState:
    longitudinal: np.ndarray
    transverse: np.ndarray
    frames: np.ndarray|None=None
    history: list[np.ndarray]=field(default_factory=list)
    def __post_init__(self):
        self.longitudinal=np.asarray(self.longitudinal,float); self.transverse=np.asarray(self.transverse,float)
        if self.transverse.shape != self.longitudinal.shape+(2,): raise ValueError("transverse shape must be longitudinal shape + (2,)")
        if not np.isfinite(self.longitudinal).all() or not np.isfinite(self.transverse).all(): raise ValueError("state must be finite")
        if self.frames is not None:
            self.frames=np.asarray(self.frames,float)
            if self.frames.shape != self.longitudinal.shape+(3,3): raise ValueError("invalid frames")
        if not self.history: self.history=[self.as_array().copy()]
    @property
    def physical_rank(self): return 3
    def as_array(self): return np.concatenate((self.longitudinal[...,None],self.transverse),axis=-1)
