"""Dev142 inventory and mathematical controls for zero-mass energy transport.

Candidate propagation is available only for auditing a caller-supplied Q_E.  It
does not assert that any Q_E is present in PBUF.
"""
from __future__ import annotations
import numpy as np

STATUSES = ("ESTABLISHED", "DERIVABLE", "RELATION_ONLY", "MISSING_NATIVE_STATE",
            "MISSING_CONSTITUTIVE_LAW", "MISSING_NORMALIZATION", "NON_IDENTIFIABLE",
            "CIRCULAR", "REDUNDANT", "NONUNIVERSAL", "MULTIVALUED", "NOT_APPLICABLE")


def energy_state_inventory():
    return [
      {"name":"bounded-strain medium energy density W(epsilon)","equation_origin":"bounded-strain constitutive law","dimensional_role":"medium energy density","native_or_physical":"physical after K normalization","scope":"local","conservation":"not a propagated-mode invariant","static_or_dynamic":"static medium state","source_dependent":True,"propagation_dependent":False,"status":"RELATION_ONLY"},
      {"name":"fast/slow response transfer","equation_origin":"frozen 0.03 delta_u_fast + 0.003 delta_u_slow","dimensional_role":"response increment","native_or_physical":"native","scope":"local","conservation":"not defined","static_or_dynamic":"transport update","source_dependent":True,"propagation_dependent":True,"status":"MISSING_CONSTITUTIVE_LAW"},
      {"name":"zero-mass mode energy","equation_origin":None,"dimensional_role":"energy","native_or_physical":None,"scope":None,"conservation":"unknown","static_or_dynamic":None,"source_dependent":None,"propagation_dependent":None,"status":"MISSING_NATIVE_STATE"},
      {"name":"action/Hamiltonian/complex amplitude","equation_origin":None,"dimensional_role":"absent","native_or_physical":None,"scope":None,"conservation":"unknown","static_or_dynamic":None,"source_dependent":None,"propagation_dependent":None,"status":"MISSING_NATIVE_STATE"}]


def momentum_state_inventory():
    return [
      {"name":"unit propagation direction","component":"momentum direction only","available":True,"status":"ESTABLISHED"},
      {"name":"scalar momentum magnitude","component":"magnitude","available":False,"status":"MISSING_NATIVE_STATE"},
      {"name":"momentum vector p n_hat","component":"vector","available":False,"status":"MISSING_NATIVE_STATE"}]


def direction_magnitude_audit(directions, magnitude=1.0):
    n=np.asarray(directions,float)
    if n.ndim != 2 or n.shape[1] not in (2,3) or len(n)<2: raise ValueError("directions must be Nx2 or Nx3")
    norms=np.linalg.norm(n,axis=1)
    if np.any(norms==0): raise ValueError("directions must be nonzero")
    unit=n/norms[:,None]; p=unit*float(magnitude)
    return {"MOMENTUM_DIRECTION_AVAILABLE":True,
            "MOMENTUM_DIRECTION_CHANGED":bool(np.max(np.linalg.norm(unit-unit[0],axis=1))>1e-12),
            "MOMENTUM_MAGNITUDE_CHANGED":bool(np.ptp(np.linalg.norm(p,axis=1))>1e-12),
            "momentum_magnitude_physically_defined":False}


def accumulate_log_energy(path, q_energy, *, orientation="forward"):
    s=np.asarray(path,float); q=np.asarray(q_energy,float)
    if s.ndim != 1 or q.shape != s.shape or len(s)<2: raise ValueError("matching 1-D path and Q_E required")
    sign={"forward":1.0,"reverse":-1.0}.get(orientation)
    if sign is None: raise ValueError("orientation must be forward or reverse")
    return np.r_[0.0,np.cumsum(.5*(q[1:]+q[:-1])*np.diff(s)*sign)]


def propagate_energy_ratio(path, q_energy, *, orientation="forward", reversible=True):
    if orientation == "reverse" and not reversible: raise ValueError("NON_REVERSIBLE_CANDIDATE")
    return np.exp(accumulate_log_energy(path,q_energy,orientation=orientation))


def scale_cancellation(q_native, ds_native, L0):
    if L0 <= 0: raise ValueError("L0 must be positive")
    native=float(q_native*ds_native); physical=float((q_native/L0)*(L0*ds_native))
    return {"native_log_increment":native,"physical_log_increment":physical,
            "exact_cancellation":bool(np.isclose(native,physical,rtol=0,atol=1e-15))}


def candidate_registry():
    qnames=["existing zero-mass energy state","existing zero-mass momentum state","trajectory direction as momentum direction","momentum magnitude from transport state","energy from momentum and c","momentum from energy and c","wave number from p/hbar","wave number from E/(hbar c)","wavelength from hc/E","redshift from energy ratio","redshift from momentum ratio","redshift from k ratio","local medium loading -> energy evolution","accumulated response -> energy evolution","response gradient -> energy evolution","fast transfer -> energy evolution","slow transfer -> energy evolution","combined fast/slow transfer -> energy evolution","strain -> energy evolution","strain gradient -> energy evolution","trajectory curvature -> energy evolution","curvature integral -> energy evolution","path excess -> energy evolution","bundle deformation -> energy evolution","entry/exit energy transfer","energy conservation in uniform medium","forward/reverse energy closure","scale-free energy-ratio transport","redshift stopping from energy evolution","source reconstruction at energy stop","multipath common-source consistency","energy-depth vs geometry-depth convergence","inferred spatial wave-state closure","corrected quantum-to-PBUF bridge classification","missing constitutive quantity identification"]
    established={3,26,34,35}; relation={5,6,7,8,9,10,11,12}; rows=[]
    for i,name in enumerate(qnames,1):
        status="ESTABLISHED" if i in established else "RELATION_ONLY" if i in relation else "MISSING_NATIVE_STATE" if i in (1,2,4,33) else "MISSING_CONSTITUTIVE_LAW"
        rows.append({"candidate_id":f"Q{i:02d}","name":name,"status":status})
    return rows


def energy_driver_registry():
    names=["response gradient","strain gradient","trajectory curvature","c_state","accumulated response u","delta_u_fast","delta_u_slow","frozen fast+slow transfer","strain epsilon","path excess","bundle area change","bundle divergence","local Hessian structure","response-curvature cross term","medium energy-density difference","bounded-strain energy difference"]
    return [{"candidate_id":f"E{i:02d}","name":n,
             "dimensionally_inverse_length":i in (1,2,3),
             "physical_magnitude_coupling_established":False,
             "status":"MISSING_CONSTITUTIVE_LAW"} for i,n in enumerate(names,1)]
