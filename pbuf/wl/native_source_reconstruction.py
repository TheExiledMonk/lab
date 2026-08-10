"""Forward-constrained native source reconstruction for Dev138."""
from __future__ import annotations
from dataclasses import dataclass,field
from typing import Any,Callable
import numpy as np
from .native_source_depth import candidates_from_curve,depth_consensus

@dataclass(frozen=True)
class NativeSourceReconstruction:
    source_depth_native: float|None
    source_depth_candidates: tuple
    interaction_region_depth: float|None
    source_positions_native: np.ndarray
    source_size_metrics: dict[str,float]
    scale_free_distance_ratios: dict=field(default_factory=dict)
    scale_free_size_ratios: dict=field(default_factory=dict)
    depth_consensus: dict=field(default_factory=dict)
    depth_ambiguity: str="UNRESOLVED"
    morphology_metrics: dict=field(default_factory=dict)
    roundtrip_score: float=float("inf")
    state_information_used: tuple=()
    deleted_information_controls: dict=field(default_factory=dict)
    scale_classification: str="NATIVE_SCALE_DEPENDENT"
    provenance: dict=field(default_factory=dict)

def source_size_metrics(points):
    p=np.asarray(points,float); c=p.mean(0); x=p-c; cov=np.cov(x.T,bias=True); eig=np.sort(np.maximum(np.linalg.eigvalsh(cov),0))[::-1]
    rms=float(np.sqrt(np.mean(np.sum(x*x,axis=1)))); major,minor=np.sqrt(eig)*2
    return {"R_rms_native":rms,"R_area_native":float(np.sqrt(major*minor)),"major_axis_native":float(major),"minor_axis_native":float(minor),"axis_ratio":float(minor/major) if major else 0.,"source_area_native2":float(np.pi*major*minor)}

def forward_constrained_inverse(depths,received_positions,candidate_builder,forward,*,received_directions=None,direction_weight=1.):
    """Search depth without exposing truth; candidates are re-forwarded."""
    obs=np.asarray(received_positions,float); scores=[]; clouds=[]
    for z in np.asarray(depths,float):
      src=np.asarray(candidate_builder(float(z)),float); pred=forward(src,float(z));
      pp,pd=(pred if isinstance(pred,tuple) else (pred,None)); q=np.mean((np.asarray(pp)-obs)**2)
      if received_directions is not None and pd is not None:q+=direction_weight*np.mean((np.asarray(pd)-received_directions)**2)
      scores.append(float(q));clouds.append(src)
    scores=np.asarray(scores); cands=candidates_from_curve(depths,scores,"D32"); j=int(np.argmin(scores)); consensus=depth_consensus(cands)
    return NativeSourceReconstruction(float(depths[j]),tuple(cands),None,clouds[j],source_size_metrics(clouds[j]),depth_consensus=consensus,
      depth_ambiguity="MULTIPLE_DEPTH_CANDIDATES" if len(cands)>1 else "UNIQUE_CANDIDATE",roundtrip_score=float(scores[j]),
      state_information_used=("arrival_position",)+(('arrival_direction',) if received_directions is not None else ()),provenance={"method":"forward_constrained_inverse","truth_access":False}),scores

def ambiguity_area(surface,threshold_fraction=.05):
    q=np.asarray(surface,float); cut=np.nanmin(q)+threshold_fraction*(np.nanmax(q)-np.nanmin(q)); return float(np.mean(q<=cut))

