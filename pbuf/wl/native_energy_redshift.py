"""Energy-ratio redshift histories and native-path stop inversion."""
from __future__ import annotations
import numpy as np


def redshift_from_energy_ratio(energy_ratio):
    r=np.asarray(energy_ratio,float)
    if np.any(~np.isfinite(r)) or np.any(r<=0): raise ValueError("energy ratio must be finite and positive")
    return 1.0/r-1.0


def redshift_from_mode_energy_proxy(mode_energy_ratio, *, proxy_established: bool):
    if not proxy_established:
        raise ValueError("mode-energy redshift is gated on an established proxy")
    return redshift_from_energy_ratio(mode_energy_ratio)


def energy_ratio_from_redshift(redshift):
    z=np.asarray(redshift,float)
    if np.any(~np.isfinite(z)) or np.any(z<=-1): raise ValueError("redshift must be finite and greater than -1")
    return 1.0/(1.0+z)


def stopping_candidates(path, redshift, target, *, atol=1e-12):
    s=np.asarray(path,float); z=np.asarray(redshift,float)
    if s.ndim != 1 or z.shape != s.shape or len(s)<2: raise ValueError("matching 1-D path and redshift required")
    y=z-float(target); roots=[]
    for i in range(len(s)-1):
        if abs(y[i])<=atol: roots.append(float(s[i]))
        if y[i]*y[i+1]<0: roots.append(float(s[i]+np.diff(s)[i]*(-y[i])/(y[i+1]-y[i])))
    if abs(y[-1])<=atol: roots.append(float(s[-1]))
    return tuple(dict.fromkeys(round(x,14) for x in roots))


def energy_redshift_stop(path, energy_ratio, target_redshift, *, mechanism=None, scale_free=False):
    roots=stopping_candidates(path,redshift_from_energy_ratio(energy_ratio),target_redshift)
    classification="NO_STOP" if not roots else "UNIQUE_STOP" if len(roots)==1 else "MULTIPLE_STOPS"
    return {"contract":"PBUF_ZERO_MASS_ENERGY_REDSHIFT_STOP_V1","target_redshift":float(target_redshift),
            "stop_candidates":list(roots),"native_stop_depth":roots[0] if len(roots)==1 else None,
            "classification":classification,"ambiguity":"NO_ENERGY_REDSHIFT_STOP" if not roots else classification,
            "candidate_mechanism":mechanism,"scale_free":bool(scale_free),"time_required":False}


def multipath_comparison(source_energy_ratios):
    vals=np.asarray(source_energy_ratios,float)
    if vals.size==0: return {"path_count":0,"compatible":False,"branch_spread":None}
    return {"path_count":int(vals.size),"compatible":bool(np.allclose(vals,vals[0],rtol=1e-8,atol=1e-10)),
            "branch_spread":float(np.ptp(vals))}
