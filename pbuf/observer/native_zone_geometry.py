"""Unfitted classifications for frozen DEV206 N6-zone geometry."""
from __future__ import annotations
import numpy as np
from .native_n6_zone import zone_tensor, symmetric_antisymmetric, axial_dual

def zone_geometry(directed: np.ndarray) -> dict[str,np.ndarray]:
    g=zone_tensor(directed); s,a=symmetric_antisymmetric(g)
    return {'tensor':g,'symmetric':s,'antisymmetric':a,'axial':axial_dual(a)}

def exact_transversality(v: np.ndarray, direction: np.ndarray) -> str:
    n=np.asarray(direction,float); n/=np.linalg.norm(n)
    longitudinal=np.sum(np.asarray(v)*n,axis=-1)
    if np.array_equal(longitudinal,np.zeros_like(longitudinal)): return 'PURE_TRANSVERSE'
    return 'MIXED'

def handedness(polar: np.ndarray, axial: np.ndarray, direction: np.ndarray) -> str:
    n=np.asarray(direction,float); n/=np.linalg.norm(n)
    sign=np.sum(np.cross(polar,axial)*n,axis=-1)
    if np.all(sign>0) or np.all(sign<0): return 'UNIDIRECTIONAL'
    if np.all(sign==0): return 'CANCELING'
    return 'BIDIRECTIONAL'
