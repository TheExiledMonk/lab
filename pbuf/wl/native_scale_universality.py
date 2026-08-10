"""Deterministic universality, convergence, agreement and clustering tools."""
from __future__ import annotations
import math
import numpy as np

def stability(values, controls=None):
    x=np.asarray(values,dtype=float)
    if not len(x) or np.any(~np.isfinite(x)) or np.any(x<=0):
        return {"classification":"NOT_TESTABLE","cv":None,"log_slope":None,"values":x.tolist()}
    cv=float(np.std(x)/np.mean(x))
    if controls is None: controls=np.arange(1,len(x)+1,dtype=float)
    slope=float(np.polyfit(np.log(np.asarray(controls,dtype=float)),np.log(x),1)[0]) if len(x)>1 else 0.
    cls="STABLE" if cv<=.02 else "WEAKLY_DRIFTING" if cv<=.10 else "FAIL"
    return {"classification":cls,"cv":cv,"log_slope":slope,"values":x.tolist()}

def resolution_stability(values,resolutions=(32,48,64,96,128)):
    out=stability(values,resolutions)
    if out["log_slope"] is not None and abs(abs(out["log_slope"])-1)<.1:
      out.update(classification="RESOLUTION_DEPENDENT",rejection_reason="REJECT_GRID_RESOLUTION_ARTIFACT")
    return out

def pairwise_agreement(candidates):
    finite=[c for c in candidates if c.get("L0_m_per_native") is not None]
    ids=[c["candidate_id"] for c in finite]
    matrix=[]
    for a in finite:
      matrix.append([abs(math.log(a["L0_m_per_native"]/b["L0_m_per_native"])) for b in finite])
    return {"candidate_ids":ids,"delta_abs_log_matrix":matrix}

def cluster_candidates(candidates,tolerance=.10):
    finite=sorted((c for c in candidates if c.get("L0_m_per_native") and c.get("independent_of_target",True) and c.get("independent_of_lcdm",True)),key=lambda x:x["L0_m_per_native"])
    if not finite:return {"outcome":"NO_FINITE_CANDIDATES","clusters":[]}
    groups=[]
    for c in finite:
      if not groups or abs(math.log(c["L0_m_per_native"]/groups[-1][-1]["L0_m_per_native"]))>tolerance: groups.append([c])
      else: groups[-1].append(c)
    outcome="ONE_DOMINANT_CLUSTER" if len(groups)==1 else "MULTIPLE_INCOMPATIBLE_CLUSTERS"
    return {"outcome":outcome,"clusters":[{"candidate_ids":[x["candidate_id"] for x in g],"independence_classes":sorted(set(x["independence_class"] for x in g)),"geometric_mean_m_per_native":float(np.exp(np.mean(np.log([x["L0_m_per_native"] for x in g]))))} for g in groups]}
