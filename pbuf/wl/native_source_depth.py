"""Native-coordinate source-depth diagnostics (Dev138).

This module is diagnostic only: it neither changes nor mathematically reverses
the PBUF propagation law. Candidate paths are supplied by frozen transport or a
forward-constrained inverse.
"""
from __future__ import annotations
from dataclasses import dataclass,field,asdict
from typing import Any,Callable
import numpy as np

ESTIMATOR_CLASSES=("POSITION","DIRECTION","POSITION_DIRECTION","BUNDLE","SECOND_ORDER","TOPOLOGY","MULTIPATH","TRAJECTORY_HISTORY","COMPLEXITY","ROUNDTRIP","GEOMETRIC_CONTROL")
_IDS=tuple(f"D{i:02d}" for i in range(1,36))
_CLASS=("POSITION","POSITION","POSITION","BUNDLE","POSITION","DIRECTION","POSITION_DIRECTION","POSITION_DIRECTION","BUNDLE","SECOND_ORDER","TOPOLOGY","TOPOLOGY","TOPOLOGY","TOPOLOGY","TOPOLOGY","POSITION_DIRECTION","GEOMETRIC_CONTROL","GEOMETRIC_CONTROL","MULTIPATH","MULTIPATH","MULTIPATH","MULTIPATH","MULTIPATH","BUNDLE","POSITION","BUNDLE","TRAJECTORY_HISTORY","SECOND_ORDER","COMPLEXITY","COMPLEXITY","ROUNDTRIP","ROUNDTRIP","GEOMETRIC_CONTROL","GEOMETRIC_CONTROL","ROUNDTRIP")

@dataclass(frozen=True)
class NativeDepthCandidate:
    depth_native: float; score: float; estimator_id: str; estimator_class: str
    support_width: float=0.; local_curvature: float=0.; ambiguity_rank: int=1
    source_component_count: int=1; metadata: dict[str,Any]=field(default_factory=dict)

def estimator_registry():
    return [{"estimator_id":i,"estimator_class":c,"attempted":True} for i,c in zip(_IDS,_CLASS)]

def deterministic_depth_grid(z_min,z_max,coarse=256,refine=4,best=4,max_values=2048,score=None):
    if not np.isfinite([z_min,z_max]).all() or z_max<=z_min: raise ValueError("invalid native depth domain")
    grid=np.linspace(z_min,z_max,int(coarse))
    if score is None:return grid
    q=np.asarray(score(grid),float); step=grid[1]-grid[0]
    centers=grid[np.argsort(q,kind="stable")[:best]]
    extra=np.concatenate([np.linspace(max(z_min,c-step),min(z_max,c+step),2*refine+1) for c in centers])
    return np.unique(np.r_[grid,extra])[:max_values]

def _cov(x):
    x=np.asarray(x,float); return np.cov(x.T,bias=True) if len(x)>1 else np.zeros((x.shape[1],x.shape[1]))

def primitive_score_curves(depths,positions,directions=None):
    """Compute interpretable D01-D08 primitives from candidate event paths.

    ``positions`` has shape (depth,event,2); directions is analogous.
    """
    z=np.asarray(depths,float); p=np.asarray(positions,float)
    if p.shape[:1]!=(len(z),) or p.ndim!=3 or p.shape[2]!=2: raise ValueError("positions must be (depth,event,2)")
    traces=np.array([np.trace(_cov(x)) for x in p]); dets=np.array([max(0.,np.linalg.det(_cov(x))) for x in p])
    area=np.sqrt(dets); sep=np.array([np.var(np.linalg.norm(x[:,None]-x[None,:],axis=2)) for x in p])
    out={"D01":area,"D02":dets,"D03":traces,"D04":sep,"D05":1/(traces+np.finfo(float).eps)}
    if directions is not None:
      d=np.asarray(directions,float); ds=np.array([np.trace(_cov(x)) for x in d])
      phase=np.array([max(0.,np.linalg.det(_cov(np.c_[p[i],d[i]]))) for i in range(len(z))])
      centered=p-p.mean(1,keepdims=True); angular=1-np.mean(np.abs(np.sum(centered*d,axis=2))/(np.linalg.norm(centered,axis=2)*np.linalg.norm(d,axis=2)+1e-15),axis=1)
      out.update(D06=ds,D07=phase,D08=angular)
    return out

def candidates_from_curve(depths,scores,estimator_id,tolerance=.05):
    z=np.asarray(depths,float); q=np.asarray(scores,float)
    if len(z)!=len(q) or not len(z) or np.any(~np.isfinite(q)): return []
    span=float(np.ptp(q)); cut=float(q.min()+tolerance*span); idx=np.flatnonzero(q<=cut)
    groups=np.split(idx,np.flatnonzero(np.diff(idx)>1)+1) if len(idx) else []
    out=[]
    for rank,g in enumerate(sorted(groups,key=lambda a:(q[a].min(),z[a[np.argmin(q[a])]])),1):
      j=int(g[np.argmin(q[g])]); curv=0. if j in (0,len(z)-1) else float((q[j-1]-2*q[j]+q[j+1])/((z[j+1]-z[j])**2))
      out.append(NativeDepthCandidate(float(z[j]),float(q[j]),estimator_id,_CLASS[int(estimator_id[1:])-1],float(z[g[-1]]-z[g[0]]),curv,rank))
    return out

def depth_consensus(candidates,tolerance=.1):
    """Equal-class-weight consensus; raw estimator multiplicity gives no vote."""
    usable=[c for c in candidates if np.isfinite(c.depth_native)]
    by={k:np.median([c.depth_native for c in usable if c.estimator_class==k]) for k in set(c.estimator_class for c in usable)}
    if not by:return {"outcome":"NO_DEPTH_INFORMATION","depth_native":None,"class_estimates":{}}
    vals=np.array(list(by.values())); center=float(np.median(vals)); scale=max(abs(center),np.finfo(float).eps)
    required=any(k in by for k in ("DIRECTION","POSITION_DIRECTION")) and any(k in by for k in ("BUNDLE","SECOND_ORDER")) and any(k in by for k in ("TOPOLOGY","MULTIPATH","ROUNDTRIP"))
    agree=np.abs(vals-center)<=tolerance*scale
    outcome="STRONG_UNIQUE_DEPTH" if len(by)>=4 and required and agree.all() else "MODERATE_UNIQUE_DEPTH" if len(by)>=3 and agree.all() else "ESTIMATORS_INCONSISTENT"
    return {"outcome":outcome,"depth_native":center,"class_estimates":by,"independent_class_count":len(by)}

def scale_free_ratios(z_observer,z_lens,z_source,r_source=None,r_lens=None):
    dol=abs(float(z_lens)-float(z_observer)); dls=abs(float(z_source)-float(z_lens)); dos=abs(float(z_source)-float(z_observer))
    if min(dol,dos)<=0: raise ValueError("distance denominators must be nonzero")
    out={"D_LS_over_D_OL":dls/dol,"D_OS_over_D_OL":dos/dol,"D_LS_over_D_OS":dls/dos,
      "proof":{"distance":"(L0*D_a)/(L0*D_b)=D_a/D_b"}}
    if r_source is not None and r_lens is not None:
      if r_lens<=0:raise ValueError("lens radius must be positive")
      out.update(R_source_over_R_lens=float(r_source)/float(r_lens));out["proof"]["size"]="(L0*R_s)/(L0*R_l)=R_s/R_l"
    return out

def localize_interaction_region(depths,diagnostics,tolerance=.1):
    z=np.asarray(depths,float); rows=[]
    for name,values in diagnostics.items():
      q=np.asarray(values,float); j=int(np.nanargmax(np.abs(q))); rows.append({"estimator":name,"z_l_candidate":float(z[j]),"width":0.,"event_fraction_supporting":1.})
    vals=np.array([r["z_l_candidate"] for r in rows]); center=float(np.median(vals)) if len(vals) else None
    agree=int(np.sum(np.abs(vals-center)<=tolerance*max(abs(center),1e-15))) if len(vals) else 0
    return {"status":"DOMINANT_INTERACTION_REGION_LOCALIZED" if agree>=3 else "DOMINANT_INTERACTION_REGION_BROAD_OR_MULTIPLE","interaction_centroid_depth":center,"estimators":rows}

