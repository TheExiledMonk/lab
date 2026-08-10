"""Passive packet and spatial-progression diagnostics for Dev146.

Packet construction is synthetic; no function advances a packet dynamically or
accepts beta as an input.
"""
from __future__ import annotations
import numpy as np

CENTER_NAMES=("excitation centroid","excitation-energy-like centroid","absolute-state centroid",
 "packet-support centroid","peak position","conserved-norm centroid")
INTERNAL_NAMES=("sign alternation","transverse orientation change","local excitation cycling",
 "packet deformation without centroid motion","link-state exchange","fast/slow exchange",
 "loading/excitation exchange","bounded local state cycling")
MECHANISM_NAMES=("excitation centroid progression","excitation packet group progression",
 "internal-vs-translational partition","excitation/loading exchange","conserved total excitation norm",
 "conserved local pair norm","link excitation transfer","standing+traveling state decomposition",
 "transverse excitation transport","fast/slow excitation decomposition","loading-modulated packet progression",
 "bounded-strain-modulated excitation","spatial periodicity change","localized mode plus translating envelope",
 "packet dispersion from loading","internal cycling with reduced translation","native quadratic progression norm",
 "native nonquadratic progression norm","excitation definition found, speed still missing","current state insufficient")

def gaussian_packet(points, center=0.0, width=1.0, magnitude=1.0):
    x=np.asarray(points,dtype=float); width=float(width); magnitude=float(magnitude)
    if width<=0 or not np.isfinite(width) or not np.isfinite(magnitude): raise ValueError("finite magnitude and positive width required")
    return magnitude*np.exp(-.5*((x-float(center))/width)**2)

def center_diagnostics(points,state):
    x,q=np.broadcast_arrays(np.asarray(points,dtype=float),np.asarray(state,dtype=float)); a=np.abs(q)
    if a.sum()==0: raise ValueError("centroid undefined for zero packet")
    support=a > a.max()*1e-6
    return {"C01":float(np.sum(x*q)/np.sum(q)) if not np.isclose(np.sum(q),0) else None,
            "C02":None,"C03":float(np.sum(x*a)/np.sum(a)),"C04":float(np.mean(x[support])),
            "C05":float(x[np.argmax(a)]),"C06":None,
            "C02_status":"ENERGY_LIKE_WEIGHT_UNESTABLISHED","C06_status":"CONSERVED_NORM_UNESTABLISHED"}

def progression_ratio(unloaded_center_delta, loaded_center_delta):
    """Measure beta only from independently supplied native histories."""
    maximum=float(unloaded_center_delta); translation=float(loaded_center_delta)
    if not np.isfinite(maximum) or not np.isfinite(translation) or maximum<=0: raise ValueError("valid measured progressions required")
    beta=translation/maximum
    if beta<0 or beta>1: raise ValueError("measured progression violates causal audit range")
    return beta

def mechanism_audit():
    rows=[]
    for i,name in enumerate(MECHANISM_NAMES,1):
        status="MEASUREMENT_DEFINITION_ONLY" if i in (1,2,13) else "MISSING_DYNAMIC_EXCITATION_LAW" if i<19 else "ESTABLISHED"
        rows.append({"id":f"D{i:02d}","name":name,"attempted":True,"status":status,"beta_emerged":False})
    return rows

def spatial_dispersion_audit():
    return {"native_relation":"D(k_n,X,L)=0","independent_mode_quantities_available":False,
            "spatial_wavelength_available":"SYNTHETIC_INPUT_ONLY","omega_created":False,
            "loaded_dispersion_established":False,"status":"NOT_DERIVABLE_FROM_STATIC_STATE"}

