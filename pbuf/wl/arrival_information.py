"""Target-blind information audits for Dev130 arrival events."""
from __future__ import annotations
import hashlib
import numpy as np
from .receiver_information import rank_metrics, linear_reconstruction_r2, classify_r2


def _summary(q):
    q=np.asarray(q,float)
    return {"median":float(np.median(q)) if q.size else 0.,"mean":float(np.mean(q)) if q.size else 0.,"minimum":float(np.min(q)) if q.size else 0.,"maximum":float(np.max(q)) if q.size else 0.,"p10":float(np.quantile(q,.1)) if q.size else 0.,"p90":float(np.quantile(q,.9)) if q.size else 0.}


def family_arrays(events):
    g=events.event_geometry;l=events.local_relations
    return {"A0":{k:g[k] for k in ("arrival_u","arrival_v","intersection_t")},
      "A1":{k:g[k] for k in ("arrival_dir_u","arrival_dir_v","arrival_dir_n")},
      "A2":{k:g[k] for k in ("receiver_incidence_cosine","receiver_incidence_angle","forward_distance","surface_residual")},
      "A3":{k:v for k,v in l.items() if any(k.endswith(q) for q in ("cov_u_du","cov_u_dv","cov_v_du","cov_v_dv"))},
      "A4":{k:v for k,v in l.items() if k not in {q for q in l if any(q.endswith(x) for x in ("cov_u_du","cov_u_dv","cov_v_du","cov_v_dv"))}}}


def matrix(events,families):
    fam=family_arrays(events);cols=[];names=[];mask=np.asarray(events.validity_masks["primary"],bool).copy()
    for f in families:
        for k,v in sorted(fam[f].items()):
            q=np.asarray(v,float);cols.append(q);names.append(f+":"+k);mask &= np.isfinite(q)
    return np.column_stack(cols) if cols else np.empty((events.ray_count,0)),names,mask


def geometry_rank_ladder(events):
    rows=[];prev=0.
    for i in range(5):
        x,_,m=matrix(events,[f"A{k}" for k in range(i+1)]);r=rank_metrics(x[m]);r.update(stage=f"G{i}",retained_sample_count=int(m.sum()),incremental_effective_rank=r["effective_rank"]-prev);prev=r["effective_rank"];rows.append(r)
    return rows


def reconstruction_audits(events):
    tests=(("arrival_direction","A0_to_A1",("A0",),("A1",)),("incidence","A0_A1_to_A2",("A0","A1"),("A2",)),("position_direction_coupling","A0_A1_A2_to_A3",("A0","A1","A2"),("A3",)),("arrival_relations","A0_A1_A2_A3_to_A4",("A0","A1","A2","A3"),("A4",)))
    out=[]
    for label,name,pf,tf in tests:
        x,_,xm=matrix(events,pf);y,yn,ym=matrix(events,tf);m=xm&ym;r2=linear_reconstruction_r2(y[m],x[m]);s=_summary(r2);out.append({"family_test":label,"mapping":name,"per_channel_r2":dict(zip(yn,r2.tolist())),"summary":s,"classification":classify_r2(s["median"]),"retained_sample_count":int(m.sum())})
    return out


def endpoint_comparison(events, receive_position, plane):
    x=np.asarray(receive_position,float);du=events.event_geometry["arrival_u"]-(x-plane.origin)@plane.e_u;dv=events.event_geometry["arrival_v"]-(x-plane.origin)@plane.e_v;valid=np.isfinite(du)&np.isfinite(dv)
    def stats(q):
        a=np.abs(q[valid]);return {"median_abs":float(np.median(a)),"rms":float(np.sqrt(np.mean(q[valid]**2))),"p95":float(np.quantile(a,.95)),"p99":float(np.quantile(a,.99)),"maximum":float(a.max())}
    radial=np.hypot(du[valid],dv[valid]);mx=float(radial.max())
    classification="ENDPOINT_ALREADY_ON_RECEIVER_SURFACE" if mx<=1e-12 else ("ENDPOINT_NEAR_RECEIVER_SURFACE" if mx<=1e-8 else "EXPLICIT_INTERSECTION_REQUIRED")
    return {"delta_u":stats(du),"delta_v":stats(dv),"radial_maximum":mx,"classification":classification,"thresholds":{"already_surface":1e-12,"near_surface":1e-8}}


def intersection_statistics(events):
    t=np.asarray(events.event_geometry["intersection_t"],float);finite=np.isfinite(t);q=t[finite];n=len(t)
    return {"minimum":float(q.min()),"median":float(np.median(q)),"mean":float(q.mean()),"p95":float(np.quantile(q,.95)),"p99":float(np.quantile(q,.99)),"maximum":float(q.max()),"fraction_near_zero":float(np.mean(np.abs(t[finite])<=1e-12)),"fraction_positive":float(np.sum(t>0)/n),"fraction_negative":float(np.sum(t<0)/n),"fraction_parallel":float(np.mean(events.event_geometry["intersection_status"]=="PARALLEL"))}


def information_preservation(events, dev129_primary, relational=None, bundle=None):
    n=events.ray_count;row=np.asarray(events.receiver_reference["receiver_row_index"]);identity=np.array_equal(row,np.arange(n)) and np.array_equal(events.receiver_reference["ray_index"],np.asarray(dev129_primary["ray_index"]))
    groups={"primary":dev129_primary,"relational":relational or {},"bundle":bundle or {}};bad=[];count=0;digests={}
    for group,data in groups.items():
        for k,v in data.items():
            a=np.asarray(v);count+=1
            if a.shape[:1]!=(n,):bad.append(f"{group}:{k}")
            else:digests[f"{group}:{k}"]=hashlib.sha256(np.ascontiguousarray(a).tobytes()).hexdigest()
    if not identity:bad.append("receiver_row_identity")
    return {"DEV129_RECEIVER_FIELDS_LOST":len(bad),"field_count":count,"lost_fields":bad,"foreign_key_exact":identity,"field_sha256":digests,"status":"PASS" if not bad else "WL_RECEIVER_TO_ARRIVAL_INFORMATION_PRESERVATION_FAILURE"}


def survival_audit(events, channels, label):
    x,_,xm=matrix(events,("A0","A1","A2","A3","A4"));names=[];cols=[]
    for k,v in sorted(channels.items()):
        q=np.asarray(v); 
        if q.ndim==1 and np.issubdtype(q.dtype,np.number):cols.append(q.astype(float));names.append(k)
    if not cols:return {"group":label,"available_channels":0,"status":"UNAVAILABLE"}
    y=np.column_stack(cols);m=xm&np.all(np.isfinite(y),1);r2=linear_reconstruction_r2(y[m],x[m]);s=_summary(r2)
    return {"group":label,"available_channels":len(names),"per_channel_r2":dict(zip(names,r2.tolist())),"summary":s,"independent":bool(s["median"]<.99),"retained_sample_count":int(m.sum())}
