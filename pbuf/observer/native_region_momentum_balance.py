"""Exact N6 bond-cut accounting for a fixed node region."""
from __future__ import annotations
import numpy as np
from pbuf.excitation.native_vector_pair_dynamics import pair_forces

def bond_cut(displacement: np.ndarray, mask: np.ndarray) -> np.ndarray:
    fp=pair_forces(displacement); total=np.zeros(3)
    for axis in range(3):
        inside=mask; other=np.roll(mask,-1,axis=axis); cut=inside & ~other
        total += fp[...,axis,:][cut].sum(axis=0)
        reverse=~inside & other
        total -= fp[...,axis,:][reverse].sum(axis=0)
    return total
