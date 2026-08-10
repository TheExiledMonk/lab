"""Deterministic information-preservation audits for Dev129 receiver banks."""
from __future__ import annotations
import numpy as np
from .receiver_state import FAMILIES, EPS, receiver_matrix


def standardized_matrix(X):
    X=np.asarray(X,np.float64)
    if X.shape[0] == 0: return np.empty_like(X),np.zeros(X.shape[1],bool)
    mean=X.mean(0);std=X.std(0)
    keep=np.isfinite(std)&(std>0)
    return (X[:,keep]-mean[keep])/std[keep],keep


def rank_metrics(X):
    Z,_=standardized_matrix(X)
    if not Z.size:return {"channel_count":0,"effective_rank":0.0,"variance_dimension_95":0,"variance_dimension_99":0,"condition_number":None}
    s=np.linalg.svd(Z,compute_uv=False);power=s*s
    if power.sum()==0:return {"channel_count":Z.shape[1],"effective_rank":0.0,"variance_dimension_95":0,"variance_dimension_99":0,"condition_number":None}
    p=power/power.sum();nz=p[p>0];tol=s[0]*max(Z.shape)*np.finfo(float).eps;positive=s[s>tol]
    return {"channel_count":int(Z.shape[1]),"effective_rank":float(np.exp(-np.sum(nz*np.log(nz)))),
        "variance_dimension_95":int(np.searchsorted(np.cumsum(p),.95)+1),"variance_dimension_99":int(np.searchsorted(np.cumsum(p),.99)+1),
        "condition_number":float(positive[0]/positive[-1]) if len(positive) else None,"numerical_rank":int(len(positive))}


def rank_ladder(state):
    rows=[];previous=0.0
    for k in range(10):
        X,names,mask,audit=receiver_matrix(state,FAMILIES[:k+1]);metrics=rank_metrics(X[mask]);metrics.update(stage=f"R{k}",retained_sample_count=int(mask.sum()),incremental_effective_rank=metrics["effective_rank"]-previous)
        previous=metrics["effective_rank"];rows.append(metrics)
    return rows


def linear_reconstruction_r2(target,predictors):
    y=np.asarray(target,float);x=np.asarray(predictors,float)
    if y.ndim==1:y=y[:,None]
    if x.ndim==1:x=x[:,None]
    if x.shape[1]==0 or y.shape[1]==0 or len(x)==0:return np.zeros(y.shape[1])
    x=(x-x.mean(0))/(x.std(0)+EPS);ymean=y.mean(0);coef=np.linalg.lstsq(np.column_stack((np.ones(len(x)),x)),y,rcond=None)[0]
    pred=np.column_stack((np.ones(len(x)),x))@coef;den=np.sum((y-ymean)**2,axis=0)
    return np.divide(den-np.sum((y-pred)**2,axis=0),den,out=np.zeros_like(den),where=den>0)


def classify_r2(median):
    return "STRONGLY_INDEPENDENT" if median<.5 else ("PARTIALLY_INDEPENDENT" if median<.9 else ("MOSTLY_REDUNDANT" if median<.99 else "REDUNDANT"))


def _summary(values):
    q=np.asarray(values,float)
    if q.size==0:return {"median":0.0,"mean":0.0,"p10":0.0,"p90":0.0,"minimum":0.0,"maximum":0.0}
    return {"median":float(np.median(q)),"mean":float(np.mean(q)),"p10":float(np.quantile(q,.1)),"p90":float(np.quantile(q,.9)),"minimum":float(q.min()),"maximum":float(q.max())}


def family_reconstruction(state):
    out=[]
    for k in range(1,10):
        prior,pnames,pmask,_=receiver_matrix(state,FAMILIES[:k]);new,nnames,nmask,_=receiver_matrix(state,(FAMILIES[k],));mask=pmask&nmask
        forward=linear_reconstruction_r2(new[mask],prior[mask]);reverse=linear_reconstruction_r2(prior[mask],new[mask])
        summary=_summary(forward);out.append({"family":FAMILIES[k],"per_channel_r2":dict(zip(nnames,forward.tolist())),"summary":summary,"classification":classify_r2(summary["median"]),"reverse_summary":_summary(reverse),"retained_sample_count":int(mask.sum())})
    return out


def leave_one_family_out(state):
    full,_,fm,_=receiver_matrix(state);rfull=rank_metrics(full[fm])["effective_rank"];rows=[]
    for family in FAMILIES:
        X,_,mask,_=receiver_matrix(state,tuple(f for f in FAMILIES if f!=family));r=rank_metrics(X[mask])["effective_rank"]
        rows.append({"family":family,"full_effective_rank":rfull,"without_effective_rank":r,"rank_loss":rfull-r})
    return rows


def packet_preservation(state, packet):
    from .receiver_state import _arrays
    a=_arrays(packet);direct={}
    for family in ("C0","C1","C3","C4","C5","C6"):
        direct.update(state.channel_bank[family])
    results={};lost=[]
    endpoint_map={"endpoint_receive_position":("receive_u","receive_v","receive_w"),"endpoint_final_direction":("final_dir_u","final_dir_v","final_dir_w")}
    for key,value in a.items():
        v=np.asarray(value)
        candidates=[]
        if key in endpoint_map:candidates=[direct.get(q) for q in endpoint_map[key]]
        elif key.startswith("path_"):
            name=key[5:]; candidates=[direct.get(name),direct.get(name.replace("path_curvature","curvature"))]
        elif key.startswith("native_"):candidates=[direct.get(key[7:])]
        candidates=[q for q in candidates if q is not None]
        exact=False
        if v.ndim==2 and len(candidates)==v.shape[1]: exact=np.array_equal(v,np.column_stack(candidates))
        elif v.ndim==1 and candidates: exact=np.array_equal(v,np.asarray(candidates[0]))
        # Launch position and initial direction live in arrival_state, exactly, but are control fields rather than channels.
        if key=="endpoint_launch_position": exact=np.array_equal(v,state.arrival_state["global_launch_position"])
        if key=="endpoint_initial_direction": exact=np.array_equal(v,np.column_stack([direct[f"initial_dir_{q}"] for q in "uvw"]))
        results[key]={"r2":1.0 if exact else 0.0,"exact":bool(exact)}
        if not exact:lost.append(key)
    return {"DEV128_PACKET_RECONSTRUCTION_R2":results,"DIRECT_TRAJECTORY_FIELDS_LOST":len(lost),"lost_fields":lost,
            "status":"PASS" if not lost else "RECEIVER_INFORMATION_PRESERVATION_FAILURE"}
