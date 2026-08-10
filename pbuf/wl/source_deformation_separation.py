"""Target-blind source-coherence / deformation separation for received rays.

This module deliberately has no benchmark, lens-registration, or target imports.
All scales, weights, and thresholds are deterministic functions of received state.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
from scipy import ndimage

from .observer_volume3d import construct_volume


RAW_CHANNELS = ("occupancy", "mean_delta_u", "mean_delta_v", "mean_delta_w",
                "mean_dir_u", "mean_dir_v", "mean_dir_w")


def robust_z(a, valid):
    x = np.asarray(a, float); q = x[valid & np.isfinite(x)]
    med = np.median(q) if q.size else 0.; mad = np.median(np.abs(q-med)) if q.size else 1.
    scale = 1.4826*mad
    if not np.isfinite(scale) or scale <= np.finfo(float).eps:
        scale = np.std(q) if q.size else 1.
    return np.nan_to_num((x-med)/max(scale, np.finfo(float).eps))


def neighbor_mean(a, valid):
    a=np.nan_to_num(np.asarray(a,float)); v=valid.astype(float)
    kernel=np.zeros((3,3,3)); kernel[1,1,:]=1; kernel[1,:,1]=1; kernel[:,1,1]=1; kernel[1,1,1]=0
    n=ndimage.convolve(v,kernel,mode="constant"); s=ndimage.convolve(a*v,kernel,mode="constant")
    return np.divide(s,n,out=np.zeros_like(s),where=n>0), n


def _variance(volume, name):
    mean=np.nan_to_num(volume[f"mean_{name}"]); occ=volume["occupancy"]
    second=np.divide(volume[f"sum_squared_{name}"],occ,out=np.zeros_like(occ),where=occ>0)
    return np.maximum(second-mean*mean,0)


def build_feature_bank(volume):
    """Construct the frozen primitive bank and complementary continuous fields."""
    valid=np.asarray(volume["occupancy"])>0
    raw=[]
    for name in RAW_CHANNELS:
        a=np.asarray(volume[name],float); raw.append(robust_z(np.nan_to_num(a),valid))
    agreement=[]; residual=[]; connected=None
    for a in raw:
        mean,n=neighbor_mean(a,valid); residual.append(np.abs(a-mean)); agreement.append(np.exp(-np.abs(a-mean)))
        connected=n/6.
    du,dv,dw=(raw[i] for i in (1,2,3)); ru,rv,rw=(residual[i] for i in (1,2,3))
    dru,drv,drw=(raw[i] for i in (4,5,6)); rdu,rdv,rdw=(residual[i] for i in (4,5,6))
    dvar=np.stack([_variance(volume,f"delta_{a}") for a in "uvw"],-1)
    qvar=np.stack([_variance(volume,f"dir_{a}") for a in "uvw"],-1)
    delta_cov=np.linalg.norm(dvar,axis=-1); direction_cov=np.linalg.norm(qvar,axis=-1)
    anis=np.std(dvar,axis=-1)/(np.mean(dvar,axis=-1)+1e-12)
    cross=[]
    for a in "uvw":
        for b in "uvw":
            cross.append(np.nan_to_num(volume[f"mean_delta_{a}_dir_{b}"])-np.nan_to_num(volume[f"mean_delta_{a}"])*np.nan_to_num(volume[f"mean_dir_{b}"]))
    cross_norm=np.linalg.norm(np.stack(cross,-1),axis=-1)
    depth_var=_variance(volume,"received_depth")
    covdiag=np.stack((dvar,qvar),-1).reshape(valid.shape+(6,))
    p=covdiag/(covdiag.sum(-1,keepdims=True)+1e-12); effective_rank=np.exp(-np.sum(np.where(p>0,p*np.log(p+1e-30),0),-1))
    spatial_cont=np.exp(-np.sqrt(ru*ru+rv*rv+rw*rw)); direction_cont=np.exp(-np.sqrt(rdu*rdu+rdv*rdv+rdw*rdw))
    occ_cont=agreement[0]; neighbor_agreement=np.mean(agreement,axis=0)
    neighbor_asym=np.std(np.stack(residual,axis=-1),axis=-1)
    occupancy_residual=residual[0]; transverse_residual=np.hypot(ru,rv); direction_residual=np.sqrt(rdu*rdu+rdv*rdv+rdw*rdw)
    depth_mean,_=neighbor_mean(raw[3],valid); depth_residual=np.abs(raw[3]-depth_mean)
    source={"occupancy_continuity":occ_cont,"direction_continuity":direction_cont,
            "spatial_continuity":spatial_cont,"neighbor_feature_agreement":neighbor_agreement,
            "local_connectedness":connected}
    deform={"occupancy_residual":occupancy_residual,"transverse_displacement_residual":transverse_residual,
            "direction_residual":direction_residual,"cross_tensor_norm":cross_norm,
            "anisotropy":anis,"depth_residual":depth_residual}
    source_score=np.mean([np.maximum(robust_z(v,valid),0) for v in source.values()],axis=0)
    deformation_score=np.mean([np.abs(robust_z(v,valid)) for v in deform.values()],axis=0)
    source_score[~valid]=0; deformation_score[~valid]=0
    primitives={"occupancy":volume["occupancy"],"occupancy_residual":occupancy_residual,
        "mean_delta_u":volume["mean_delta_u"],"mean_delta_v":volume["mean_delta_v"],"mean_delta_w":volume["mean_delta_w"],
        "mean_dir_u":volume["mean_dir_u"],"mean_dir_v":volume["mean_dir_v"],"mean_dir_w":volume["mean_dir_w"],
        "delta_covariance":delta_cov,"direction_covariance":direction_cov,"cross_covariance":cross_norm,
        "anisotropy_magnitude":anis,"depth_variance":depth_var,"local_effective_rank":effective_rank,
        "neighbor_continuity":neighbor_agreement,"neighbor_asymmetry":neighbor_asym}
    return source,deform,primitives,source_score,deformation_score


def latent_decomposition(primitives, valid, components=8):
    names=tuple(primitives); X=np.column_stack([robust_z(primitives[n],valid)[valid] for n in names])
    X-=X.mean(0); _,s,vt=np.linalg.svd(X,full_matrices=False); k=min(components,len(s))
    scores=X@vt[:k].T; fields=np.zeros(valid.shape+(k,)); fields[valid]=scores
    explained=s*s/(np.sum(s*s)+1e-30)
    return {"latent_fields":np.moveaxis(fields,-1,0),"singular_values":s,"svd_components":vt[:k],
            "pca_components":vt[:k],"pca_explained_variance_ratio":explained[:k],
            "feature_names":np.array(names),"methods":np.array(["centered_SVD","PCA_via_centered_SVD"])}


def project(field, occupancy):
    occ=np.asarray(occupancy,float); den=occ.sum(2)
    return {"sum":np.sum(field,2),"occupancy_weighted_mean":np.divide(np.sum(field*occ,2),den,out=np.zeros_like(den),where=den>0),
            "rms":np.sqrt(np.mean(field*field,2))}


def decompose_rays(rays, shape=(64,64,64), uv_bounds=(-8.,8.), with_latent=True):
    volume,metadata=construct_volume(rays,shape=shape,uv_bounds=uv_bounds)
    source,deform,primitives,S,D=build_feature_bank(volume); latent=latent_decomposition(primitives,volume["occupancy"]>0) if with_latent else {}
    return {"volume":volume,"metadata":metadata,"source":source,"deformation":deform,"primitives":primitives,
            "S":S,"D":D,"latent":latent,"S2":project(S,volume["occupancy"]),"D2":project(D,volume["occupancy"])}


def structural_sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()
