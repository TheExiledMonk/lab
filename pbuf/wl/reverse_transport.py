"""Reverse metadata and structural inversion diagnostics for Dev131."""
from __future__ import annotations
from dataclasses import dataclass, field
from typing import Any, Mapping, Sequence
import numpy as np

@dataclass(frozen=True)
class ReverseCandidateSet:
    observation_identity: str
    candidate_launch_identities: tuple[int, ...]
    candidate_states: tuple[Any, ...] = ()
    uniqueness_classification: str = "UNRESOLVED"
    conditioning: Mapping[str, Any] = field(default_factory=dict)
    constraint_provenance: tuple[str, ...] = ()

    @classmethod
    def from_candidates(cls, observation_identity, launch_ids, states=(), **kw):
        ids=tuple(int(x) for x in launch_ids)
        label="NO_LAUNCH_CANDIDATE" if not ids else ("UNIQUE_INVERSE" if len(ids)==1 else "MULTIPLE_LAUNCH_CANDIDATES")
        return cls(str(observation_identity),ids,tuple(states),label,**kw)

def affine_inverse(arrival, matrix, offset):
    x=np.asarray(arrival,float);a=np.asarray(matrix,float);b=np.asarray(offset,float)
    if a.shape != (2,2) or abs(np.linalg.det(a)) <= 128*np.finfo(float).eps:
        return None,"NON_UNIQUE_INVERSE"
    return (x-b)@np.linalg.inv(a).T,"UNIQUE_INVERSE"

def transport_diagnostics(matrices, *, singular_ratio=1e-12, poor_ratio=1e-4, moderate_ratio=1e-2):
    a=np.asarray(matrices,float)
    if a.ndim != 3 or a.shape[1:] != (2,2): raise ValueError("matrices must have shape (N,2,2)")
    finite=np.all(np.isfinite(a),axis=(1,2));s=np.full((len(a),2),np.nan);s[finite]=np.linalg.svd(a[finite],compute_uv=False)
    smax=s[:,0];smin=s[:,1];ratio=np.divide(smin,smax,out=np.full(len(a),np.nan),where=smax>0)
    det=np.full(len(a),np.nan);det[finite]=np.linalg.det(a[finite]);labels=np.full(len(a),"UNRESOLVED",dtype="U28")
    labels[finite&(ratio>=moderate_ratio)]="LOCALLY_WELL_CONDITIONED"
    labels[finite&(ratio<moderate_ratio)&(ratio>=poor_ratio)]="LOCALLY_ANISOTROPIC"
    labels[finite&(ratio<poor_ratio)&(ratio>=singular_ratio)]="LOCALLY_ILL_CONDITIONED"
    labels[finite&(ratio<singular_ratio)]="LOCALLY_SINGULAR_CANDIDATE"
    amplification=np.divide(1.,smin,out=np.full(len(a),np.inf),where=smin>0)
    return {"detJ":det,"sigma_min":smin,"sigma_max":smax,"transport_condition_number":ratio,
            "reverse_sensitivity":amplification,"classification":labels}

def reconstruct_receiver(arrival_point, arrival_direction, intersection_t):
    p=np.asarray(arrival_point,float);d=np.asarray(arrival_direction,float);t=np.asarray(intersection_t,float)
    return p-t[...,None]*d

def roundtrip_errors(reconstructed, original):
    delta=np.asarray(reconstructed,float)-np.asarray(original,float);flat=np.abs(delta).ravel()
    scale=max(1.,float(np.nanmax(np.abs(original))))
    mx=float(np.nanmax(flat));tol=512*np.finfo(float).eps*scale
    label="EXACT_ROUNDTRIP" if np.array_equal(reconstructed,original) else ("NUMERICALLY_EXACT_ROUNDTRIP" if mx<=tol else "APPROXIMATE_ROUNDTRIP")
    return {"classification":label,"tolerance":tol,"bitwise_equal_rows":int(np.all(delta==0,axis=1).sum()),
            "max_abs":mx,"rms":float(np.sqrt(np.nanmean(delta*delta))),
            "p95":float(np.nanpercentile(flat,95)),"p99":float(np.nanpercentile(flat,99)),
            "per_coordinate": [{"max_abs":float(np.nanmax(np.abs(delta[:,i]))),"rms":float(np.sqrt(np.nanmean(delta[:,i]**2))),
                                "p95":float(np.nanpercentile(np.abs(delta[:,i]),95)),"p99":float(np.nanpercentile(np.abs(delta[:,i]),99))} for i in range(delta.shape[1])]}

def correspondence_index(launch_ids, receiver_ids):
    launch=np.asarray(launch_ids,dtype=np.int64);receiver=np.asarray(receiver_ids,dtype=np.int64);events=np.arange(len(launch),dtype=np.int64)
    order=np.lexsort((events,launch));keys,counts=np.unique(launch[order],return_counts=True);offset=np.r_[0,np.cumsum(counts)].astype(np.int64)
    rorder=np.lexsort((events,receiver));rkeys,rcounts=np.unique(receiver[rorder],return_counts=True);roffset=np.r_[0,np.cumsum(rcounts)].astype(np.int64)
    return {"event_to_launch":launch,"launch_keys":keys,"launch_offsets":offset,"launch_event_indices":events[order],
            "receiver_keys":rkeys,"receiver_offsets":roffset,"receiver_event_indices":events[rorder]}

def arrival_knn(points, ks=(4,8,16,32,64)):
    """Deterministic KNN: distance then canonical event-index tie break."""
    from scipy.spatial import cKDTree
    p=np.asarray(points,float);maxk=max(ks);tree=cKDTree(p);dist,idx=tree.query(p,k=maxk+1,workers=1)
    out={}
    for row in range(len(p)):
        keep=idx[row]!=row; pairs=sorted(zip(dist[row][keep],idx[row][keep]),key=lambda q:(q[0],q[1]))[:maxk]
        idx[row,:maxk]=[q[1] for q in pairs]
    for k in ks:out[f"k{k}"]=idx[:,:k].astype(np.int64)
    return out

def reverse_free_propagation(final_position, direction, distance):
    """Algebraic inverse of ``x_f = x_i + distance * direction``."""
    p=np.asarray(final_position,float);d=np.asarray(direction,float);t=np.asarray(distance,float)
    return p-t[...,None]*d

def reverse_optical_record(record):
    """Return the previous state represented by an immutable optical record."""
    status=getattr(record,"aperture_status","")
    if status == "BLOCKED_BY_APERTURE":
        return ReverseCandidateSet(str(record.event_uid),(),(),
            "BLOCKED_INFORMATION_NOT_PRESENT_DOWNSTREAM",
            {"reverse_classification":"INFORMATION_LOSSY_DOWNSTREAM"},
            (str(record.surface_id),))
    state={"position":np.asarray(record.input_position,float),
           "direction":np.asarray(record.input_direction,float),
           "event_uid":str(record.event_uid),"surface_id":str(record.surface_id)}
    return ReverseCandidateSet(str(record.event_uid),(0,),(state,),"UNIQUE_INVERSE",
        {"reverse_classification":str(record.reverse_classification),"metadata_used":True},
        (str(record.surface_id),"stored input position","stored incoming direction"))

def reverse_optical_history(records):
    """Traverse one event's history from output to input without index guessing."""
    ordered=sorted(records,key=lambda r:r.interaction_index,reverse=True)
    candidates=[]
    for record in ordered:
        candidate=reverse_optical_record(record);candidates.append(candidate)
        if candidate.uniqueness_classification == "BLOCKED_INFORMATION_NOT_PRESENT_DOWNSTREAM": break
    return tuple(candidates)
