"""Target-blind 3-D receiver-state formation (Dev129).

This module is deliberately downstream of propagation.  It reorganizes a
frozen Dev128 packet and launch-defined bundle history; it contains no
propagator, detector raster, lensing target, or feedback path.
"""
from __future__ import annotations

from dataclasses import dataclass, field
import hashlib
import json
from pathlib import Path
from typing import Mapping

import numpy as np

EPS = np.finfo(np.float64).eps
SCALES = (1, 2, 4, 8, 16, 32)
FAMILIES = tuple(f"C{i}" for i in range(10))
VALIDITY_POLICY = "VALID_NEIGHBOR_ONLY"


@dataclass
class ReceiverState3D:
    """Canonical grouped received-light state; every channel is per ray."""

    ray_identity: dict = field(default_factory=dict)
    arrival_state: dict = field(default_factory=dict)
    launch_correspondence: dict = field(default_factory=dict)
    trajectory_state: dict = field(default_factory=dict)
    local_receiver_relations: dict = field(default_factory=dict)
    bundle_relations: dict = field(default_factory=dict)
    channel_bank: dict = field(default_factory=dict)
    validity_masks: dict = field(default_factory=dict)
    metadata: dict = field(default_factory=dict)

    def family(self, name):
        return self.channel_bank[name]

    @property
    def ray_count(self):
        return len(self.ray_identity["ray_index"])


def _arrays(receipt):
    if isinstance(receipt, (str, Path)):
        with np.load(receipt) as z:
            return {k: z[k] for k in z.files}
    if hasattr(receipt, "endpoint"):
        out = {f"endpoint_{k}": v for k, v in receipt.endpoint.items()}
        out.update({f"path_{k}": v for k, v in receipt.path_summary.items()})
        for native, values in receipt.native_path_summary.items():
            out.update({f"native_{native}_{k}": v for k, v in values.items()})
        return out
    return {str(k): np.asarray(v) for k, v in receipt.items()}


def _vector(a, key, n=None):
    x = np.asarray(a[key], dtype=np.float64)
    if x.ndim != 2 or x.shape[1] != 3 or (n is not None and len(x) != n):
        raise ValueError(f"{key} must have shape (ray_count, 3)")
    return x


def _scalar(a, key, n):
    x = np.asarray(a[key])
    if x.shape != (n,):
        raise ValueError(f"{key} must have shape ({n},)")
    return x.copy()


def _invariants(cuu, cuv, cvv):
    trace = cuu + cvv
    disc = np.sqrt(np.maximum(0.0, (cuu-cvv)**2 + 4*cuv*cuv))
    e1, e2 = .5*(trace+disc), .5*(trace-disc)
    return {"trace": trace, "determinant": cuu*cvv-cuv*cuv,
            "eigenvalue_1": e1, "eigenvalue_2": e2,
            "anisotropy": np.divide(e1-e2, e1+e2+EPS),
            "orientation": .5*np.arctan2(2*cuv, cuu-cvv)}


def local_arrival_relations(launch_uv, receive_position, final_direction, side,
                            scales=SCALES):
    """Launch-topology neighborhood statistics with complete windows only.

    A scale is a Chebyshev launch-cell radius.  Boundary windows are undefined
    and represented by NaN plus an explicit mask, never by fabricated rays.
    """
    uv = np.asarray(launch_uv, float).reshape(side, side, 2)
    pos = np.asarray(receive_position, float).reshape(side, side, 3)
    direction = np.asarray(final_direction, float).reshape(side, side, 3)
    if side*side != np.asarray(receive_position).shape[0]:
        raise ValueError("side does not match ray count")
    result, masks = {}, {}
    # Integral image gives O(scales * rays * fields), independent of window area.
    def box_mean(x, radius):
        x=np.asarray(x,float); pad=np.pad(x,((1,0),(1,0)),constant_values=0).cumsum(0).cumsum(1)
        out=np.full((side,side),np.nan); lo=radius; hi=side-radius
        if lo >= hi: return out
        y0=np.arange(lo,hi)-radius; y1=np.arange(lo,hi)+radius+1
        x0=y0; x1=y1
        sums=pad[y1[:,None],x1[None,:]]-pad[y0[:,None],x1[None,:]]-pad[y1[:,None],x0[None,:]]+pad[y0[:,None],x0[None,:]]
        out[lo:hi,lo:hi]=sums/float((2*radius+1)**2)
        return out
    for scale in scales:
        s=int(scale); valid=np.zeros((side,side),bool)
        if 2*s < side: valid[s:side-s,s:side-s]=True
        masks[f"C7_s{s}"]=valid.ravel()
        mu=[box_mean(pos[...,q],s) for q in range(3)]
        md=[box_mean(direction[...,q],s) for q in range(3)]
        pcuu=box_mean(pos[...,0]**2,s)-mu[0]**2
        pcuv=box_mean(pos[...,0]*pos[...,1],s)-mu[0]*mu[1]
        pcvv=box_mean(pos[...,1]**2,s)-mu[1]**2
        dcuu=box_mean(direction[...,0]**2,s)-md[0]**2
        dcuv=box_mean(direction[...,0]*direction[...,1],s)-md[0]*md[1]
        dcvv=box_mean(direction[...,1]**2,s)-md[1]**2
        pinv=_invariants(pcuu,pcuv,pcvv); dinv=_invariants(dcuu,dcuv,dcvv)
        launch_cuu=box_mean(uv[...,0]**2,s)-box_mean(uv[...,0],s)**2
        launch_cvv=box_mean(uv[...,1]**2,s)-box_mean(uv[...,1],s)**2
        area0=np.sqrt(np.maximum(0,launch_cuu*launch_cvv))
        area=np.sqrt(np.maximum(0,pinv["determinant"]))
        prefix=f"s{s}_"
        values={
            "mean_neighbor_displacement_u":mu[0]-pos[...,0],
            "mean_neighbor_displacement_v":mu[1]-pos[...,1],
            "position_cov_uu":pcuu,"position_cov_uv":pcuv,"position_cov_vv":pcvv,
            **{f"position_{k}":v for k,v in pinv.items()},
            "direction_cov_uu":dcuu,"direction_cov_uv":dcuv,"direction_cov_vv":dcvv,
            **{f"direction_{k}":v for k,v in dinv.items()},
            "mean_pair_separation":box_mean(np.hypot(pos[...,0]-mu[0],pos[...,1]-mu[1]),s),
            "separation_variance":pinv["trace"],
            "local_received_area_ratio":np.divide(area,area0,out=np.full_like(area,np.nan),where=area0>EPS),
            "local_anisotropy":pinv["anisotropy"],
            "mean_receive_w":mu[2],
            "var_receive_w":box_mean(pos[...,2]**2,s)-mu[2]**2,
            "cov_u_w":box_mean(pos[...,0]*pos[...,2],s)-mu[0]*mu[2],
            "cov_v_w":box_mean(pos[...,1]*pos[...,2],s)-mu[1]*mu[2],
        }
        result.update({prefix+k: v.ravel() for k,v in values.items()})
    return result, masks


def bundle_relation_channels(bundle, ray_count=None):
    """Preserve final components and depth summaries from Dev128 bundle data."""
    if bundle is None: return {}
    b=_arrays(bundle); out={}
    def flat(name, value):
        a=np.asarray(value); out[name]=a.reshape(-1).astype(np.float64,copy=False)
    J=np.asarray(b["jacobian"],float); H=np.asarray(b["hessian"],float)
    for i in range(2):
        for j in range(2):
            q=J[...,i,j];flat(f"final_J{i+1}{j+1}",q[-1]);flat(f"mean_J{i+1}{j+1}",q.mean(0));flat(f"rms_J{i+1}{j+1}",np.sqrt(np.mean(q*q,axis=0)))
    labels=("uu","uv","vv")
    for i in range(2):
        for j,label in enumerate(labels):
            q=H[...,i,j];flat(f"final_H{i+1}_{label}",q[-1]);flat(f"mean_H{i+1}_{label}",q.mean(0));flat(f"rms_H{i+1}_{label}",np.sqrt(np.mean(q*q,axis=0)))
    for src,dst in (("area_ratio","area_ratio"),("anisotropy","anisotropy"),("orientation","orientation"),
                    ("second_order_norm","second_order_norm")):
        q=np.asarray(b[src],float);flat(f"final_{dst}",q[-1]);flat(f"minimum_{dst}",np.nanmin(q,axis=0));flat(f"maximum_{dst}",np.nanmax(q,axis=0))
    if "orientation" in b: flat("area_orientation_change",np.asarray(b["orientation"])[-1]-np.asarray(b["orientation"])[0])
    if ray_count is not None and any(len(v)!=ray_count for v in out.values()):
        raise ValueError("bundle grid does not match receiver population")
    return out


def build_receiver_state(receipt, *, side=None, launch_grid_index=None, validity=None,
                         termination_status="provided_final_state", bundle=None,
                         scales=SCALES):
    """Build all C0--C9 families from a frozen receipt and optional bundle."""
    a=_arrays(receipt); launch=_vector(a,"endpoint_launch_position"); n=len(launch)
    receive=_vector(a,"endpoint_receive_position",n); initial=_vector(a,"endpoint_initial_direction",n)
    final=_vector(a,"endpoint_final_direction",n); side=int(round(np.sqrt(n))) if side is None else int(side)
    if side*side != n: raise ValueError("receiver relations require the frozen square launch topology")
    valid=np.ones(n,bool) if validity is None else np.asarray(validity,bool)
    if valid.shape != (n,): raise ValueError("validity must have one value per ray")
    norm=np.linalg.norm(final,axis=1); final_unit=np.divide(final,norm[:,None],out=np.zeros_like(final),where=norm[:,None]>0)
    delta=receive-launch
    path={k[5:]:np.asarray(v).copy() for k,v in a.items() if k.startswith("path_")}
    native={k[7:]:np.asarray(v).copy() for k,v in a.items() if k.startswith("native_")}
    c0={f"receive_{q}":receive[:,i] for i,q in enumerate("uvw")}
    c1={f"final_dir_{q}":final[:,i] for i,q in enumerate("uvw")}
    if not np.allclose(norm[valid],1.0,rtol=0,atol=32*EPS): c1.update({f"normalized_dir_{q}":final_unit[:,i] for i,q in enumerate("uvw")})
    c2={f"delta_{q}":delta[:,i] for i,q in enumerate("uvw")}
    c2.update(transverse_displacement_norm=np.linalg.norm(delta[:,:2],axis=1),full_displacement_norm=np.linalg.norm(delta,axis=1))
    def take(*names): return {name:path[name] for name in names if name in path}
    c3=take("path_length","straight_line_distance","path_excess")
    if {"path_length","straight_line_distance"} <= path.keys():
        c3["path_excess_fraction"]=(path["path_length"]-path["straight_line_distance"])/(path["straight_line_distance"]+EPS)
    c4=take("net_direction_change","total_direction_change","maximum_local_direction_change","path_location_of_maximum_direction_change")
    if {"net_direction_change","total_direction_change"} <= path.keys(): c4["direction_history_ratio"]=path["total_direction_change"]/(path["net_direction_change"]+EPS)
    curvature_names=("path_curvature_integral","path_curvature_squared_integral","curvature_mean","curvature_rms","curvature_max","curvature_variance","curvature_path_centroid")
    c5={k.replace("path_curvature","curvature"):path[k] for k in curvature_names if k in path}
    # C6 also carries packet fields not assigned to a more specific family;
    # this is the lossless receiver contract, not a feature-selection step.
    assigned=set(c3)|set(c4)|set(curvature_names)
    c6=dict(native)
    c6.update({k:v for k,v in path.items() if k not in assigned})
    c6.update({f"initial_dir_{q}":initial[:,i] for i,q in enumerate("uvw")})
    launch_uv=launch[:,:2]; c7,masks=local_arrival_relations(launch_uv,receive,final_unit,side,scales)
    c8=bundle_relation_channels(bundle,n) if bundle is not None else {}
    # Frozen direct cross-family set; unavailable structural inputs are omitted.
    c9={"displacement_direction_alignment":np.sum(delta*final_unit,axis=1)/(np.linalg.norm(delta,axis=1)+EPS)}
    kmag=c5.get("curvature_integral",c5.get("curvature_rms",np.zeros(n)))
    c9["displacement_curvature_alignment"]=c2["transverse_displacement_norm"]*kmag
    c9["direction_history_curvature_magnitude"]=c4.get("total_direction_change",np.zeros(n))*kmag
    c9["path_excess_curvature_magnitude"]=c3.get("path_excess",np.zeros(n))*kmag
    if "final_anisotropy" in c8: c9["bundle_anisotropy_curvature_magnitude"]=c8["final_anisotropy"]*kmag
    if "final_area_ratio" in c8:
        dispersion=next((v for k,v in c7.items() if k.endswith("direction_trace")),np.zeros(n))
        c9["bundle_area_ratio_direction_dispersion"]=c8["final_area_ratio"]*dispersion
    bank={f"C{i}":c for i,c in enumerate((c0,c1,c2,c3,c4,c5,c6,c7,c8,c9))}
    base_masks={f"C{i}":valid.copy() for i in range(10)};base_masks.update(masks)
    identity={"ray_index":np.arange(n,dtype=np.int64),
              "launch_grid_index":np.arange(n,dtype=np.int64) if launch_grid_index is None else np.asarray(launch_grid_index),
              "launch_u":launch[:,0],"launch_v":launch[:,1],"validity":valid,
              "termination_status":np.full(n,termination_status,dtype=object)}
    return ReceiverState3D(identity,{"receive_position":receive,"final_direction":final,
        "normalized_final_direction":final_unit,"global_launch_position":launch,
        "global_receive_position":receive,"global_initial_direction":initial,"global_final_direction":final},
        {"delta":delta},path,c7,c8,bank,base_masks,
        {"coordinate_contract":"native receiver basis (u,v,w) = frozen global Cartesian components",
         "direction_normalization":"COPIED_UNIT" if len(c1)==3 else "DERIVED_NORMALIZED",
         "validity_policy":VALIDITY_POLICY,"scales":list(scales),"target_access":False})


def channel_manifest(state):
    rows=[];counter=0
    for family in FAMILIES:
        for name,value in sorted(state.channel_bank.get(family,{}).items()):
            scale=None
            if family=="C7" and name.startswith("s"): scale=int(name.split("_",1)[0][1:])
            source="Dev128 trajectory packet" if family in ("C0","C1","C3","C4","C5","C6") else ("Dev128 bundle history" if family=="C8" else "receiver relation")
            kind="PRIMARY" if family in ("C0","C1","C3","C4","C5","C6","C8") and not name.startswith("normalized_") else "DERIVED"
            rows.append({"channel_id":f"DEV129_{counter:04d}","family":family,"name":name,
                "description":name.replace("_"," "),"source_fields":[source],"formula":"direct copy" if kind=="PRIMARY" else "deterministic receiver relation",
                "dtype":str(np.asarray(value).dtype),"units":"native or dimensionless as named","scale":scale,
                "classification":kind,"validity_requirements":VALIDITY_POLICY if scale else "valid ray",
                "target_access":False});counter+=1
    return sorted(rows,key=lambda r:(r["family"],-1 if r["scale"] is None else r["scale"],r["channel_id"]))


def manifest_sha256(manifest):
    payload=json.dumps(manifest,sort_keys=True,separators=(",",":"),allow_nan=False).encode()
    return hashlib.sha256(payload).hexdigest()


def receiver_matrix(state, families=FAMILIES, *, remove_aliases=True, shared_valid=True):
    """Return raw canonical matrix, names, mask, and constant/alias audit."""
    columns=[];names=[];mask=np.asarray(state.ray_identity["validity"],bool).copy()
    for family in families:
        for name,value in sorted(state.channel_bank.get(family,{}).items()):
            x=np.asarray(value,float);columns.append(x);names.append(f"{family}:{name}");mask &= np.isfinite(x)
    X=np.column_stack(columns) if columns else np.empty((state.ray_count,0))
    if not shared_valid: mask=np.asarray(state.ray_identity["validity"],bool)
    constants=[];aliases=[];keep=[];seen={}
    for j,name in enumerate(names):
        q=X[mask,j];constant=len(q)==0 or np.ptp(q)==0
        if constant: constants.append(name);continue
        digest=hashlib.sha256(np.ascontiguousarray(q).tobytes()).hexdigest()
        if remove_aliases and digest in seen and np.array_equal(q,X[mask,seen[digest]]): aliases.append({"alias":name,"primary":names[seen[digest]]});continue
        seen[digest]=j;keep.append(j)
    return X[:,keep], [names[j] for j in keep], mask, {"constant_channels":constants,"alias_channels":aliases,"retained_sample_count":int(mask.sum())}
