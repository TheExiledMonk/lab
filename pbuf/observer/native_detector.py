"""DEV187 exact receipt-surface detector accumulation.

This is a native transfer-response state, not an astronomical image or a
spin-2 observer.  Cells are existing DEV168 receipt-cell identities.
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np

@dataclass(frozen=True)
class NativeDetectorState:
    """Additive C100 response accumulated by exact native receipt cells."""
    detector_cell_ids: np.ndarray
    detector_coordinates: np.ndarray
    receipt_count: np.ndarray
    weight_sum: np.ndarray
    momentum_sum: np.ndarray
    flux_sum: np.ndarray
    displacement_sum: np.ndarray
    w02_sum: np.ndarray
    representation: str = 'NativeDetectorState/v1'

def accumulate(received_positions, native_cell_ids, weights, momentum, flux,
               displacement, w02, *, basis: np.ndarray | None = None) -> NativeDetectorState:
    """Exact, order-invariant additive receipt aggregation; no pixelization."""
    cell=np.asarray(native_cell_ids, np.int64); ids, inverse=np.unique(cell, return_inverse=True); n=len(ids)
    pos=np.asarray(received_positions,float); w=np.asarray(weights,float)
    if basis is None:  # DEV168 receipt face: native x normal, y/z detector axes.
        q=pos[:,1:3]
    else:
        q=(pos @ np.asarray(basis,float).T)
    count=np.bincount(inverse,minlength=n).astype(np.int64)
    def total(x):
        x=np.asarray(x,float)
        if x.ndim==1:return np.bincount(inverse,weights=x,minlength=n)
        return np.stack([np.bincount(inverse,weights=x[:,i],minlength=n) for i in range(x.shape[1])],axis=1)
    # Exact receipt-cell coordinate: all BOND_FLUX events in a cell share y/z.
    coord=np.stack([np.bincount(inverse,weights=q[:,i],minlength=n)/count for i in range(2)],axis=1)
    return NativeDetectorState(ids,coord,count,total(w),total(momentum),total(flux),total(displacement),total(w02))

def weighted_tensor(points, measure):
    """Central second moment of a positive native receipt measure; no spin-2."""
    p=np.asarray(points,float); m=np.asarray(measure,float); total=float(m.sum())
    if total <= 0:return np.full((2,2),np.nan),np.full(2,np.nan)
    centre=(p*m[:,None]).sum(0)/total; d=p-centre
    return (d.T*m)@d/total,centre
