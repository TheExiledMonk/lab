"""Native, threshold-free reductions for the DEV197 cross-event audit."""
from __future__ import annotations
import numpy as np


def exact_region(shape, x_index=2):
    """The frozen DEV182 7x7 B launch support, on its fixed launch face."""
    out=np.zeros(shape, dtype=bool); out[x_index,2:9,2:9]=True
    return out


def shell_index(shape, origin=(2,5,5)):
    grid=np.indices(shape)
    d=[np.minimum((grid[i]-origin[i])%shape[i], (origin[i]-grid[i])%shape[i]) for i in range(3)]
    return np.maximum.reduce(d)


def ratios(effect, baseline):
    """Pointwise magnitude ratios, with undefined (never regularized) zero bases."""
    e=np.linalg.norm(effect,axis=-1) if effect.ndim and effect.shape[-1]==3 else np.abs(effect)
    b=np.linalg.norm(baseline,axis=-1) if baseline.ndim and baseline.shape[-1]==3 else np.abs(baseline)
    defined=b != 0
    out=np.full(b.shape,np.nan); np.divide(e,b,out=out,where=defined)
    return out, defined


def parallel_transverse(effect, baseline):
    """Exact projection of a vector effect onto its nonzero native baseline."""
    norm=np.linalg.norm(baseline,axis=-1)
    defined=norm != 0
    unit=np.zeros_like(baseline); np.divide(baseline,norm[...,None],out=unit,where=defined[...,None])
    parallel=np.sum(effect*unit,axis=-1)
    transverse=effect-parallel[...,None]*unit
    return parallel, transverse, defined


def summary(values, mask=None):
    a=np.asarray(values)
    if mask is not None:
        a=a[:,mask] if a.ndim > mask.ndim else a[mask]
    finite=a[np.isfinite(a)]
    return {'defined_count':int(finite.size),'max':float(np.max(finite)) if finite.size else None,
            'l2':float(np.linalg.norm(finite)) if finite.size else None}
