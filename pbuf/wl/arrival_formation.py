"""Physical receiver-plane to continuous 2-D arrival events (Dev130).

This module is target blind and downstream only.  It neither propagates rays
nor rasterizes events.  The native basis is the frozen Dev124/Dev129 global
Cartesian receiver basis; its rows are (e_u, e_v, normal).
"""
from __future__ import annotations

from dataclasses import dataclass, field
import hashlib
import json
from typing import Mapping

import numpy as np

SCALES = (1, 2, 4, 8, 16, 32)
FAMILIES = ("A0", "A1", "A2", "A3", "A4", "A5", "A6")
EPS_PARALLEL = 64 * np.finfo(np.float64).eps
EPS_T = 128 * np.finfo(np.float64).eps
BASIS_TOL = 256 * np.finfo(np.float64).eps
SURFACE_TOL = 1e-12
NEAR_SURFACE_TOL = 1e-8


@dataclass(frozen=True)
class ReceiverPlane:
    origin: np.ndarray
    e_u: np.ndarray
    e_v: np.ndarray
    normal: np.ndarray
    handedness: str = "RIGHT_HANDED_EU_CROSS_EV_EQUALS_NORMAL"
    source_module: str = "pbuf.labs.foundation.los_consistent_ray_geometry001._propagate_g3d"
    construction_method: str = "native Cartesian basis; nominal fixed-step termination w=step*(steps-1)"

    def __post_init__(self):
        for name in ("origin", "e_u", "e_v", "normal"):
            x = np.asarray(getattr(self, name), dtype=np.float64)
            if x.shape != (3,) or not np.all(np.isfinite(x)):
                raise ValueError(f"{name} must be a finite 3-vector")
            object.__setattr__(self, name, x)
        basis = np.stack((self.e_u, self.e_v, self.normal))
        gram = basis @ basis.T
        if not np.allclose(gram, np.eye(3), rtol=0, atol=BASIS_TOL):
            raise ValueError("receiver basis is not orthonormal")
        handed = np.dot(np.cross(self.e_u, self.e_v), self.normal)
        if self.handedness == "RIGHT_HANDED_EU_CROSS_EV_EQUALS_NORMAL" and not np.isclose(handed, 1, rtol=0, atol=BASIS_TOL):
            raise ValueError("receiver basis handedness mismatch")

    def manifest(self):
        b = np.stack((self.e_u, self.e_v, self.normal)); gram = b @ b.T
        return {"origin": self.origin.tolist(), "e_u": self.e_u.tolist(),
                "e_v": self.e_v.tolist(), "normal": self.normal.tolist(),
                "handedness": self.handedness, "source_module": self.source_module,
                "construction_method": self.construction_method,
                "orthogonality_errors": {"eu_ev": float(abs(gram[0,1])), "eu_n": float(abs(gram[0,2])), "ev_n": float(abs(gram[1,2]))},
                "normalization_errors": {"e_u": float(abs(np.linalg.norm(b[0])-1)), "e_v": float(abs(np.linalg.norm(b[1])-1)), "normal": float(abs(np.linalg.norm(b[2])-1))},
                "target_access": False, "hst_pixel_access": False}


def native_receiver_plane(step=.03, steps=160):
    return ReceiverPlane(np.array([0., 0., step*(steps-1)]), np.array([1.,0.,0.]),
                         np.array([0.,1.,0.]), np.array([0.,0.,1.]))


@dataclass
class ArrivalEvent2D:
    event_geometry: dict = field(default_factory=dict)
    local_relations: dict = field(default_factory=dict)
    validity_masks: dict = field(default_factory=dict)
    receiver_reference: dict = field(default_factory=dict)
    metadata: dict = field(default_factory=dict)

    @property
    def ray_count(self): return len(self.event_geometry["ray_index"])


def _vectors(x, name):
    a=np.asarray(x,dtype=np.float64)
    if a.ndim != 2 or a.shape[1] != 3: raise ValueError(f"{name} must have shape (N,3)")
    return a


def intersect_rays(position, direction, plane, *, eps_parallel=EPS_PARALLEL, eps_t=EPS_T):
    """Vectorized oriented ray-plane intersection without direction reversal."""
    x=_vectors(position,"position"); d=_vectors(direction,"direction")
    if x.shape != d.shape: raise ValueError("position and direction shapes differ")
    finite=np.all(np.isfinite(x),1)&np.all(np.isfinite(d),1)
    norms=np.linalg.norm(d,axis=1); valid_direction=finite&(norms>eps_parallel)
    dh=np.divide(d,norms[:,None],out=np.zeros_like(d),where=norms[:,None]>eps_parallel)
    denominator=dh@plane.normal; parallel=valid_direction&(np.abs(denominator)<eps_parallel)
    usable=valid_direction&~parallel
    t=np.full(len(x),np.nan); t[usable]=((plane.origin-x[usable])@plane.normal)/denominator[usable]
    p=np.full_like(x,np.nan); p[usable]=x[usable]+t[usable,None]*dh[usable]
    status=np.full(len(x),"NONFINITE",dtype="U24")
    status[finite&~valid_direction]="INVALID_DIRECTION";status[parallel]="PARALLEL"
    status[usable&(t>eps_t)]="FORWARD_INTERSECTION"
    status[usable&(np.abs(t)<=eps_t)]="ON_SURFACE"
    status[usable&(t < -eps_t)]="BACKWARD_INTERSECTION"
    return p,t,dh,denominator,status


def _invariants(cuu,cuv,cvv):
    tr=cuu+cvv; disc=np.sqrt(np.maximum(0,(cuu-cvv)**2+4*cuv*cuv)); l1=.5*(tr+disc);l2=.5*(tr-disc)
    return {"trace":tr,"determinant":cuu*cvv-cuv*cuv,"eigenvalue_1":l1,"eigenvalue_2":l2,
            "anisotropy":np.divide(l1-l2,l1+l2,out=np.zeros_like(tr),where=np.abs(l1+l2)>np.finfo(float).eps),
            "orientation":.5*np.arctan2(2*cuv,cuu-cvv)}


def arrival_relations(launch_uv, arrival_uv, direction_uv, side, scales=SCALES):
    launch=np.asarray(launch_uv,float).reshape(side,side,2); pos=np.asarray(arrival_uv,float).reshape(side,side,2);dire=np.asarray(direction_uv,float).reshape(side,side,2)
    def mean(q,s):
        pad=np.pad(q,((1,0),(1,0)),constant_values=0).cumsum(0).cumsum(1);out=np.full((side,side),np.nan);lo=s;hi=side-s
        if lo>=hi:return out
        a=np.arange(lo,hi)-s;b=np.arange(lo,hi)+s+1
        out[lo:hi,lo:hi]=(pad[b[:,None],b[None,:]]-pad[a[:,None],b[None,:]]-pad[b[:,None],a[None,:]]+pad[a[:,None],a[None,:]])/float((2*s+1)**2)
        return out
    out={};masks={}
    for s0 in scales:
        s=int(s0); valid=np.zeros((side,side),bool)
        if 2*s<side:valid[s:-s,s:-s]=True
        masks[f"A34_s{s}"]=valid.ravel(); mp=[mean(pos[...,i],s) for i in range(2)];md=[mean(dire[...,i],s) for i in range(2)]
        pc=(mean(pos[...,0]**2,s)-mp[0]**2,mean(pos[...,0]*pos[...,1],s)-mp[0]*mp[1],mean(pos[...,1]**2,s)-mp[1]**2)
        dc=(mean(dire[...,0]**2,s)-md[0]**2,mean(dire[...,0]*dire[...,1],s)-md[0]*md[1],mean(dire[...,1]**2,s)-md[1]**2)
        cross=(mean(pos[...,0]*dire[...,0],s)-mp[0]*md[0],mean(pos[...,0]*dire[...,1],s)-mp[0]*md[1],mean(pos[...,1]*dire[...,0],s)-mp[1]*md[0],mean(pos[...,1]*dire[...,1],s)-mp[1]*md[1])
        lc0=mean(launch[...,0]**2,s)-mean(launch[...,0],s)**2;lc1=mean(launch[...,1]**2,s)-mean(launch[...,1],s)**2
        pi=_invariants(*pc);di=_invariants(*dc); area=np.sqrt(np.maximum(0,pi["determinant"]));area0=np.sqrt(np.maximum(0,lc0*lc1))
        vals={"cov_u_du":cross[0],"cov_u_dv":cross[1],"cov_v_du":cross[2],"cov_v_dv":cross[3],
              "arrival_cov_uu":pc[0],"arrival_cov_uv":pc[1],"arrival_cov_vv":pc[2],"arrival_anisotropy":pi["anisotropy"],"arrival_orientation":pi["orientation"],
              "arrival_trace":pi["trace"],"arrival_determinant":pi["determinant"],"arrival_eigenvalue_1":pi["eigenvalue_1"],"arrival_eigenvalue_2":pi["eigenvalue_2"],
              "arrival_area_ratio":np.divide(area,area0,out=np.full_like(area,np.nan),where=area0>np.finfo(float).eps),
              "direction_cov_uu":dc[0],"direction_cov_uv":dc[1],"direction_cov_vv":dc[2],"direction_anisotropy":di["anisotropy"],"direction_orientation":di["orientation"],
              "direction_trace":di["trace"],"direction_determinant":di["determinant"]}
        out.update({f"s{s}_{k}":v.ravel() for k,v in vals.items()})
    return out,masks


def form_arrival_events(position, direction, plane, *, ray_index=None, receiver_row_index=None,
                        launch_uv=None, launch_grid_index=None, side=None, scales=SCALES):
    p,t,dh,dn,status=intersect_rays(position,direction,plane);n=len(t);delta=p-plane.origin
    u=delta@plane.e_u;v=delta@plane.e_v;w=delta@plane.normal;du=dh@plane.e_u;dv=dh@plane.e_v;mu=np.abs(dn)
    validity=np.full(n,"NONFINITE",dtype="U20");validity[status=="INVALID_DIRECTION"]="INVALID_DIRECTION";validity[status=="PARALLEL"]="PARALLEL";validity[status=="BACKWARD_INTERSECTION"]="BACKWARD_ONLY";validity[status=="FORWARD_INTERSECTION"]="VALID_FORWARD";validity[status=="ON_SURFACE"]="VALID_ON_SURFACE"
    idx=np.arange(n,dtype=np.int64) if ray_index is None else np.asarray(ray_index,dtype=np.int64);ridx=idx.copy() if receiver_row_index is None else np.asarray(receiver_row_index,dtype=np.int64)
    geometry={"ray_index":idx,"receiver_row_index":ridx,"arrival_u":u,"arrival_v":v,"intersection_t":t,
      "arrival_dir_u":du,"arrival_dir_v":dv,"arrival_dir_n":dn,"receiver_incidence_cosine":mu,
      "receiver_incidence_angle":np.arccos(np.clip(mu,0,1)),"forward_distance":t,"surface_residual":w,
      "validity":validity,"intersection_status":status}
    ref={"receiver_row_index":ridx,"ray_index":idx}
    if launch_uv is not None:
        luv=np.asarray(launch_uv,float);ref.update(launch_u=luv[:,0],launch_v=luv[:,1],launch_grid_index=np.arange(n) if launch_grid_index is None else np.asarray(launch_grid_index))
    local={};masks={"primary":np.isin(validity,["VALID_FORWARD","VALID_ON_SURFACE"])}
    if side is not None and launch_uv is not None:local,m=arrival_relations(launch_uv,np.column_stack((u,v)),np.column_stack((du,dv)),int(side),scales);masks.update(m)
    return ArrivalEvent2D(geometry,local,masks,ref,{"receiver_attachment":"foreign key to immutable Dev129 row","target_access":False,"hst_pixel_access":False})


def boundary_audit(position, plane):
    w=(_vectors(position,"position")-plane.origin)@plane.normal;aw=np.abs(w)
    if np.max(aw)<=SURFACE_TOL:cls="ON_RECEIVER_PLANE"
    elif np.max(aw)<=NEAR_SURFACE_TOL:cls="NEAR_RECEIVER_PLANE"
    elif np.all(w<0):cls="BEFORE_RECEIVER_PLANE"
    elif np.all(w>0):cls="AFTER_RECEIVER_PLANE"
    else:cls="UNKNOWN"
    return {"classification":cls,"normal_residual":{"minimum":float(w.min()),"median":float(np.median(w)),"maximum":float(w.max()),"median_abs":float(np.median(aw)),"maximum_abs":float(aw.max())},"thresholds":{"on_surface":SURFACE_TOL,"near_surface":NEAR_SURFACE_TOL}}


def channel_manifest(scales=SCALES):
    specs={"A0":["arrival_u","arrival_v","intersection_t"],"A1":["arrival_dir_u","arrival_dir_v","arrival_dir_n"],"A2":["receiver_incidence_cosine","receiver_incidence_angle","forward_distance","surface_residual"]}
    rows=[];counter=0
    for family,names in specs.items():
        for name in names:
            rows.append({"channel_id":f"DEV130_{counter:04d}","family":family,"name":name,"formula":"physical ray-plane intersection or basis projection","source":"Dev129 endpoint plus frozen receiver plane","scale":None,"units":"native length, radians, or dimensionless as named","classification":"PRIMARY" if family in ("A0","A1") else "DERIVED","validity_requirements":"valid intersection","target_access":False});counter+=1
    for s in scales:
        for family,names in (("A3",["cov_u_du","cov_u_dv","cov_v_du","cov_v_dv"]),("A4",["arrival_cov_uu","arrival_cov_uv","arrival_cov_vv","arrival_anisotropy","arrival_orientation","arrival_area_ratio","arrival_trace","arrival_determinant","arrival_eigenvalue_1","arrival_eigenvalue_2","direction_cov_uu","direction_cov_uv","direction_cov_vv","direction_anisotropy","direction_orientation","direction_trace","direction_determinant"])):
            for name in names:
                rows.append({"channel_id":f"DEV130_{counter:04d}","family":family,"name":name,"formula":f"launch-topology complete-window covariance, radius {s}","source":"ArrivalEvent2D launch-neighborhood relation","scale":s,"units":"native or dimensionless as named","classification":"DERIVED","validity_requirements":"complete finite launch neighborhood","target_access":False});counter+=1
    for family,name in (("A5","dev129_trajectory_attachment"),("A6","dev129_receiver_attachment")):
        rows.append({"channel_id":f"DEV130_{counter:04d}","family":family,"name":name,"formula":"foreign-key reference by receiver_row_index","source":"immutable Dev129 artifacts","scale":None,"units":"reference","classification":"REFERENCE","validity_requirements":"matching row identity","target_access":False});counter+=1
    return rows


def canonical_sha(value):
    return hashlib.sha256(json.dumps(value,sort_keys=True,separators=(",",":"),allow_nan=False).encode()).hexdigest()
