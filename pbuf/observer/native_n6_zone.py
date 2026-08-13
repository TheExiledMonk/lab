"""Fixed-topology, read-only N6-zone operators for DEV206."""
from __future__ import annotations
import numpy as np

N6_ORIENTATIONS = np.array(((1,0,0),(-1,0,0),(0,1,0),(0,-1,0),(0,0,1),(0,0,-1)), dtype=float)

def zone_tensor(directed: np.ndarray) -> np.ndarray:
    """G_ij=sum_alpha ehat_alpha,i delta-r_alpha,j / 2.

    The fixed 1/2 is the exact opposite-pair normalization, not a fitted
    coefficient.  It is identically DEV203's directional_tensor.
    """
    return np.einsum('ai,...aj->...ij', N6_ORIENTATIONS, np.asarray(directed, float)) / 2

def zone_balance(directed: np.ndarray) -> np.ndarray:
    return np.sum(np.asarray(directed, float), axis=-2)

def opposite_pair_components(directed: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    x=np.asarray(directed, float)
    return (x[...,0::2,:]+x[...,1::2,:])/2, (x[...,0::2,:]-x[...,1::2,:])/2

def symmetric_antisymmetric(tensor: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    x=np.asarray(tensor, float)
    return (x+np.swapaxes(x,-1,-2))/2, (x-np.swapaxes(x,-1,-2))/2

def axial_dual(antisymmetric: np.ndarray) -> np.ndarray:
    a=np.asarray(antisymmetric,float)
    return np.stack((a[...,2,1],a[...,0,2],a[...,1,0]),axis=-1)
