"""Dimensionless spatial wavelength-ratio histories and stopping contracts."""
from __future__ import annotations
import numpy as np

def redshift_history_from_log_shift(log_shift):
    return np.expm1(np.asarray(log_shift,float))

def stopping_candidates(path, redshift, target, *, atol=1e-12):
    s=np.asarray(path,float); z=np.asarray(redshift,float)
    if s.shape != z.shape or s.ndim != 1: raise ValueError("matching 1-D arrays required")
    y=z-float(target); roots=[]
    for i in range(len(s)-1):
        if abs(y[i])<=atol: roots.append(float(s[i]))
        if y[i]*y[i+1] < 0:
            roots.append(float(s[i]+(s[i+1]-s[i])*(-y[i])/(y[i+1]-y[i])))
    if len(s) and abs(y[-1])<=atol: roots.append(float(s[-1]))
    return tuple(dict.fromkeys(round(x,14) for x in roots))

def spatial_redshift_stop(path, redshift, target, *, mechanism, scale_free):
    roots=stopping_candidates(path,redshift,target)
    return {"contract":"PBUF_SPATIAL_REDSHIFT_STOP_V1","candidate_mechanism":mechanism,
      "target_redshift":float(target),"native_stop_depth":roots[0] if len(roots)==1 else None,
      "stop_candidates":list(roots),"absolute_length_required":not scale_free,"time_required":False,
      "scale_free":bool(scale_free),"ambiguity":"NO_SPATIAL_REDSHIFT_STOP" if not roots else
      "MULTIPLE_SPATIAL_REDSHIFT_STOPS" if len(roots)>1 else "UNIQUE","path_provenance":"caller supplied"}

def multipath_comparison(stops):
    vals=[x for x in stops if x is not None]
    return {"path_count":len(stops),"common_stop":bool(vals and np.allclose(vals,vals[0])),
            "spread":float(np.ptp(vals)) if vals else None}
