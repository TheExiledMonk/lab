"""Dev141 spatial wave-state and constitutive candidate audit.

Synthetic laws supported here test mathematics only.  They do not establish a
PBUF wave law without independently proven wave-state and coefficient provenance.
"""
from __future__ import annotations
import numpy as np

WAVE_CLASSIFICATIONS = ("PHYSICAL_NATIVE_SPATIAL_WAVE_STATE", "DIMENSIONLESS_MODE_STRUCTURE",
 "NUMERICAL_ONLY", "DERIVED_FROM_EXTERNAL_PHYSICAL_INPUT", "NOT_WAVE_STATE", "UNUSABLE")

_CANDIDATES = (
 ("W01","local medium loading","LOCAL_LOADING","c_state"),
 ("W02","accumulated medium response","ACCUMULATED_RESPONSE","u"),
 ("W03","response gradient","GRADIENT","grad_u"),
 ("W04","fast channel","FAST_SLOW","delta_u_fast"),
 ("W05","slow channel","FAST_SLOW","delta_u_slow"),
 ("W06","frozen combined transfer","FAST_SLOW","0.03*du_fast+0.003*du_slow"),
 ("W07","strain","STRAIN","epsilon"), ("W08","strain gradient","STRAIN","grad_epsilon"),
 ("W09","trajectory curvature","TRAJECTORY","kappa"),
 ("W10","curvature integral","TRAJECTORY","integral(kappa ds)"),
 ("W11","path excess","PATH","delta_s"),
 ("W12","homogeneous medium evolution","HOMOGENEOUS_MEDIUM","homogeneous_state"))

def spatial_wave_inventory():
    return [
      {"name":"transport mode/iteration indices","classification":"NUMERICAL_ONLY","usable":False},
      {"name":"trajectory spatial path","classification":"NOT_WAVE_STATE","usable":False},
      {"name":"local variable named phase in D07","classification":"NUMERICAL_ONLY","usable":False},
      {"name":"physical k_n or lambda_n","classification":"UNUSABLE","usable":False,"reason":"absent"}]

def candidate_registry():
    rows=[]
    for cid,name,group,driver in _CANDIDATES:
        inverse_length = cid in ("W03","W08","W09")
        rows.append({"candidate_id":cid,"name":name,"dependency_class":group,"driver":driver,
          "wave_quantity_required":"lambda_n or k_n", "native_dimensions":"spatial/native",
          "coefficient_dimensions":"none only if driver is proven inverse native length" if inverse_length else "required and unresolved",
          "L0_cancels":inverse_length,"absolute_length_required":not inverse_length,
          "dimensionless_redshift_possible":inverse_length,
          "status":"MISSING_NATIVE_STATE", "provenance_established":False,
          "reversible":None,"free_coefficient":not inverse_length})
    return rows

def accumulate_log_wavelength(path, q, *, orientation="forward"):
    s=np.asarray(path,float); v=np.asarray(q,float)
    if s.ndim != 1 or v.shape != s.shape or len(s)<2: raise ValueError("matching 1-D path and q required")
    sign=1.0 if orientation=="forward" else -1.0 if orientation=="reverse" else None
    if sign is None: raise ValueError("orientation must be forward or reverse")
    increments=.5*(v[1:]+v[:-1])*np.diff(s)*sign
    return np.r_[0.,np.cumsum(increments)]

def evolve_wavelength(path, q, wavelength0=1., *, orientation="forward", reversible=True):
    if wavelength0 <= 0: raise ValueError("wavelength must be positive")
    if orientation=="reverse" and not reversible:
        raise ValueError("NON_REVERSIBLE_CANDIDATE")
    return wavelength0*np.exp(accumulate_log_wavelength(path,q,orientation=orientation))

def scale_cancellation(q_native, ds_native, L0):
    if L0 <= 0: raise ValueError("L0 must be positive")
    return {"native":float(q_native*ds_native),"physical":float((q_native/L0)*(L0*ds_native)),
            "exact_cancellation":bool(q_native*ds_native==(q_native/L0)*(L0*ds_native))}
