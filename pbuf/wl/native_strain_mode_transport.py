"""Passive histories and controls for already-produced strain packets."""
from __future__ import annotations
import numpy as np


def scalar_history(path_position_native, values, *, packet_uid="packet-0"):
    s, v = np.asarray(path_position_native, float), np.asarray(values, float)
    if s.ndim != 1 or v.shape != s.shape or np.any(~np.isfinite(v)): raise ValueError("matching finite 1-D arrays required")
    return {"packet_uid": packet_uid, "path_position_native": s, "propagation_progression_index": np.arange(len(s)),
            "candidate_scalar_state": v}


def conservation_classification(values, tolerance=0.05):
    v = np.asarray(values, float)
    if v.size == 0 or np.any(~np.isfinite(v)): return {"classification":"NOT_WELL_DEFINED", "cv":None}
    mean = float(np.mean(v)); cv = float(np.std(v) / abs(mean)) if mean else (0.0 if np.all(v == 0) else float("inf"))
    return {"classification":"CONSERVED" if cv <= tolerance else "NUMERICALLY_DRIFTING", "cv":cv}


def direction_scalar_control(directions, scalar_values, tolerance=0.05):
    d, q = np.asarray(directions, float), np.asarray(scalar_values, float)
    if d.ndim != 2 or len(d) != len(q): raise ValueError("one direction per scalar required")
    norms=np.linalg.norm(d,axis=1)
    if np.any(norms == 0): raise ValueError("directions must be nonzero")
    turns=np.arccos(np.clip(np.sum((d/norms[:,None])[:-1]*(d/norms[:,None])[1:],axis=1),-1,1))
    c=conservation_classification(q,tolerance)
    return {**c,"direction_changed":bool(np.any(turns>1e-10)),"direction_independent":c["classification"]=="CONSERVED"}


def entry_exit_control(before, inside, after, tolerance=0.05):
    b,i,a=map(float,(before,inside,after)); scale=max(abs(b),1e-15)
    if abs(a-b)/scale <= tolerance: kind="TEMPORARY_EXCHANGE" if abs(i-b)/scale>tolerance else "ZERO_NET_CHANGE"
    else: kind="PERMANENT_MODE_CHANGE"
    return {"before":b,"inside":i,"after":a,"classification":kind,"relative_closure_error":abs(a-b)/scale}


def multi_region_control(initial, after_first, after_second, reverse_order_final=None, tolerance=0.05):
    vals=np.asarray([initial,after_first,after_second],float)
    result={"values":vals.tolist(),"classification":"MULTIPLICATIVE" if np.all(vals>0) else "NONLINEAR"}
    if reverse_order_final is not None:
        result["order_dependent"] = bool(abs(float(reverse_order_final)-vals[-1]) > tolerance*max(abs(vals[-1]),1e-15))
    return result


def forward_reverse_control(initial, final, tolerance=0.05):
    err=abs(float(final)-float(initial))/max(abs(float(initial)),1e-15)
    return {"closure_error":err,"conservative":bool(err<=tolerance)}
