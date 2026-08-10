"""Geometric invariants of a frozen two-dimensional transport map."""
from __future__ import annotations
import numpy as np


def geometry_from_derivatives(first, second=None):
    """Return area, rotation, stretch and spin-2 shape without naming them shear."""
    j11=1+np.asarray(first["d_u_delta_u"],float);j12=np.asarray(first["d_v_delta_u"],float)
    j21=np.asarray(first["d_u_delta_v"],float);j22=1+np.asarray(first["d_v_delta_v"],float)
    J=np.stack((j11,j12,j21,j22),axis=-1).reshape(j11.shape+(2,2))
    area=np.linalg.det(J)
    u,sv,vh=np.linalg.svd(J);R=u@vh
    orientation=np.arctan2(R[...,1,0],R[...,0,0])
    C=J@np.swapaxes(J,-1,-2);tr=C[...,0,0]+C[...,1,1]+1e-30
    q1=(C[...,0,0]-C[...,1,1])/tr;q2=2*C[...,0,1]/tr
    anisotropy=np.hypot(q1,q2);axis=.5*np.arctan2(q2,q1)
    curvature=np.zeros_like(area)
    if second is not None:
        keys=[k for k in second if k.startswith("d_")]
        curvature=np.sqrt(sum(np.asarray(second[k],float)**2 for k in keys))
    return {"transport_area_ratio":area,"local_area_change":area-1,
            "local_orientation":orientation,"orientation_change":orientation,
            "local_anisotropy":anisotropy,"local_curvature":curvature,
            "axis_ratio":np.divide(sv[...,1],sv[...,0],out=np.ones_like(area),where=sv[...,0]>0),
            "principal_axis":axis,"spin2_shape_q1":q1,"spin2_shape_q2":q2,
            "jacobian":J}


def orientation_spin2(theta_received,theta_launch=0.):
    d=np.asarray(theta_received)-np.asarray(theta_launch)
    return np.cos(2*d),np.sin(2*d)


def geometry_feature_matrix(g):
    keys=("transport_area_ratio","local_orientation","local_anisotropy","local_curvature",
          "axis_ratio","spin2_shape_q1","spin2_shape_q2")
    return np.column_stack([np.asarray(g[k]).ravel() for k in keys])
