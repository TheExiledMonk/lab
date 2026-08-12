"""Read-only native field reductions used by DEV198."""
from __future__ import annotations
import numpy as np


def region_norm(field, mask):
    """L2, max, and mean magnitude; no tail/support threshold is involved."""
    a=np.asarray(field)[mask]
    m=np.linalg.norm(a,axis=-1) if a.ndim and a.shape[-1]==3 else np.abs(a)
    return {'l2':float(np.linalg.norm(m)), 'max':float(np.max(m)),
            'mean':float(np.mean(m)), 'raw_count':int(m.size)}


def force_orientation(residual_force, fresh_force, mask):
    """Cellwise native cosine, undefined where either exact norm is zero."""
    a=np.asarray(residual_force)[mask]; b=np.asarray(fresh_force)[mask]
    an=np.linalg.norm(a,axis=-1); bn=np.linalg.norm(b,axis=-1); defined=(an!=0)&(bn!=0)
    c=np.full(an.shape,np.nan); np.divide(np.sum(a*b,axis=-1),an*bn,out=c,where=defined)
    return c,defined


def ratio(effect, baseline):
    """Unregularized magnitude ratio, retaining undefined exact zero bases."""
    e=np.linalg.norm(effect,axis=-1) if effect.shape[-1:]==(3,) else np.abs(effect)
    b=np.linalg.norm(baseline,axis=-1) if baseline.shape[-1:]==(3,) else np.abs(baseline)
    out=np.full(b.shape,np.nan); np.divide(e,b,out=out,where=b!=0)
    return out,b!=0
