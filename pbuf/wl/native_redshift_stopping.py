"""Dimensionless shift-history crossing and footprint helpers for Dev140."""
from __future__ import annotations
import numpy as np


def stopping_depths(path, shift, target, *, atol=1e-12):
    s=np.asarray(path,float); z=np.asarray(shift,float)
    if s.ndim != 1 or z.shape != s.shape or len(s)<2 or np.any(np.diff(s)<=0): raise ValueError("ordered 1-D history required")
    roots=[]
    for i in range(len(s)-1):
        a=z[i]-target; b=z[i+1]-target
        if abs(a)<=atol: roots.append(float(s[i]))
        if a*b < 0:
            roots.append(float(s[i]+(s[i+1]-s[i])*(-a)/(b-a)))
    if abs(z[-1]-target)<=atol: roots.append(float(s[-1]))
    roots=list(dict.fromkeys(round(x,14) for x in roots))
    status="NO_REDSHIFT_STOP_SOLUTION" if not roots else "MULTIPLE_REDSHIFT_STOP_CANDIDATES" if len(roots)>1 else "UNIQUE_REDSHIFT_STOP_CANDIDATE"
    return {"status":status,"stop_candidates":roots}


def monotonicity(shift):
    d=np.diff(np.asarray(shift,float))
    if np.all(d>0) or np.all(d<0): return "STRICT_MONOTONIC"
    if np.mean(d>=0)>=.9 or np.mean(d<=0)>=.9: return "MOSTLY_MONOTONIC"
    return "NON_MONOTONIC"


def footprint(points):
    p=np.asarray(points,float); c=p.mean(0); x=p-c; cov=np.cov(x.T,bias=True); eig=np.sort(np.maximum(np.linalg.eigvalsh(cov),0))[::-1]
    major,minor=2*np.sqrt(eig)
    return {"centroid":c.tolist(),"RMS_radius":float(np.sqrt(np.mean(np.sum(x*x,axis=1)))),
            "major_axis":float(major),"minor_axis":float(minor),"axis_ratio":float(minor/major) if major else 0,
            "area_equivalent_radius":float(np.sqrt(major*minor)),"component_count":1}


def multipath_consistency(depths):
    x=np.asarray(depths,float); m=float(np.median(x)); cv=float(np.std(x)/abs(np.mean(x))) if np.mean(x) else float("inf")
    return {"median_stopping_depth":m,"branch_to_branch_CV":cv,"common_depth_support":bool(cv<=.1)}
