"""Target-blind forward native-depth fingerprint tools."""
from __future__ import annotations
import numpy as np

FINGERPRINT_IDS=tuple(f"F{i:02d}" for i in range(1,21))

def monotonicity(values,depths,minimum_fraction=.8):
    q=np.asarray(values,float);z=np.asarray(depths,float)
    if len(q)!=len(z) or len(q)<2 or np.any(~np.isfinite(q)):return {"classification":"UNSTABLE","monotonic_fraction":0.}
    d=np.diff(q)/np.diff(z); pos=np.mean(d>=0);neg=np.mean(d<=0); frac=float(max(pos,neg))
    cls="STRICTLY_MONOTONIC" if frac==1 else "MOSTLY_MONOTONIC" if frac>=minimum_fraction else "NON_MONOTONIC"
    return {"classification":cls,"monotonic_fraction":frac,"derivative":d.tolist()}

def build_fingerprint_bank(depths,feature_rows,control_ids,split_fn):
    x=np.asarray(feature_rows,float);z=np.asarray(depths,float)
    if x.ndim!=2 or len(x)!=len(z) or len(control_ids)!=len(z):raise ValueError("fingerprint row mismatch")
    split=np.array([split_fn(str(i)) for i in control_ids]); return {"depths":z,"features":x,"control_ids":np.asarray(control_ids),"split":split}

def nearest_fingerprint(observation,bank,*,training_only=True):
    mask=bank["split"]=="TRAIN" if training_only else np.ones(len(bank["depths"]),bool)
    if not np.any(mask):raise ValueError("no training fingerprints")
    x=bank["features"][mask];o=np.asarray(observation,float); scale=np.std(x,axis=0);scale[scale==0]=1
    j=int(np.argmin(np.sum(((x-o)/scale)**2,axis=1))); return float(bank["depths"][mask][j])
