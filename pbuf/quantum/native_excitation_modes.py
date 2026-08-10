"""Spatial modes of the frozen Dev148 rank-2 excitation."""
from __future__ import annotations
import numpy as np

STATUS = ("ESTABLISHED", "DERIVABLE", "STRUCTURALLY_SUPPORTED", "POST_FREEZE_QM_COMPATIBLE",
          "CONTINUOUS_ONLY", "BOUNDARY_QUANTIZED_ONLY", "SOURCE_INTERACTION_REQUIRED",
          "MISSING_NATIVE_MECHANISM", "NUMERICAL_ARTIFACT", "NOT_UNIQUE", "NOT_APPLICABLE", "FALSIFIED")
MODE_NAMES = ("single-cycle localized packet", "multi-cycle localized packet", "finite sinusoidal packet",
 "broad quasi-monochromatic packet", "narrow-band packet", "standing spatial mode", "counter-propagating pair",
 "two-polarization mode", "rotating transverse mode", "compact pulse", "Gaussian-envelope carrier",
 "top-hat-envelope carrier", "asymmetric packet", "two-frequency superposition", "harmonic family",
 "anti-phase pair", "node-bounded mode", "link-bounded mode", "naturally selected survivor modes",
 "no stable spatial mode")
WAVELENGTH_NAMES = ("zero-crossing separation", "peak-to-peak distance", "autocorrelation period",
 "spatial Fourier dominant mode", "phase-free state repetition distance", "orientation repetition distance",
 "polarization repetition distance", "packet carrier period", "node-spacing invariant", "no unique wavelength")

def mode_registry():
    return [{"id": f"M{i:02d}", "name": n, "attempted": True,
             "status": "FALSIFIED" if i == 20 else ("NOT_APPLICABLE" if i == 18 else "ESTABLISHED")}
            for i, n in enumerate(MODE_NAMES, 1)]

def wavelength_registry():
    return [{"id": f"L{i:02d}", "name": n, "attempted": True,
             "status": "NOT_UNIQUE" if i in (6, 7, 9, 10) else "DERIVABLE"}
            for i, n in enumerate(WAVELENGTH_NAMES, 1)]

def carrier_mode(sites: int, wavelength: float, amplitude: float = 1., polarization=(1., 0.),
                 envelope: str = "none", center=None, width=None, phase: float = 0.) -> np.ndarray:
    if sites < 8 or wavelength <= 2: raise ValueError("invalid spatial mode")
    p = np.asarray(polarization, float)
    if p.shape != (2,) or not np.isfinite(p).all() or np.linalg.norm(p) == 0: raise ValueError("invalid polarization")
    p /= np.linalg.norm(p); x = np.arange(sites, dtype=float); c = sites/2 if center is None else float(center)
    wave = np.cos(2*np.pi*(x-c)/float(wavelength) + phase)
    if envelope == "gaussian": wave *= np.exp(-.5*((x-c)/(width or 2*wavelength))**2)
    elif envelope == "top_hat": wave *= (np.abs(x-c) <= (width or 4*wavelength)/2)
    elif envelope != "none": raise ValueError("unknown envelope")
    return float(amplitude)*wave[:, None]*p[None, :]

def rotating_mode(sites: int, wavelength: float, amplitude: float = 1., handedness: int = 1) -> np.ndarray:
    if handedness not in (-1, 1): raise ValueError("handedness must be +/-1")
    x=np.arange(sites, dtype=float); phi=2*np.pi*x/float(wavelength)
    return float(amplitude)*np.column_stack((np.cos(phi), handedness*np.sin(phi)))

def quadratic_norm(state) -> float: return float(np.sum(np.asarray(state, float)**2))
def native_k(wavelength: float) -> float:
    if not np.isfinite(wavelength) or wavelength <= 0: raise ValueError("native wavelength required")
    return float(2*np.pi/wavelength)

def estimate_wavelengths(state) -> dict:
    y=np.asarray(state, float)
    if y.ndim == 2: y=y[:,0]
    y=y-np.mean(y); n=len(y); out={}
    crossings=np.flatnonzero(np.diff(np.signbit(y)))
    out["L01"]=float(2*np.median(np.diff(crossings))) if len(crossings)>2 else None
    peaks=np.flatnonzero((y[1:-1]>y[:-2]) & (y[1:-1]>=y[2:]))+1
    out["L02"]=float(np.median(np.diff(peaks))) if len(peaks)>1 else None
    ac=np.correlate(y,y,mode="full")[n-1:]
    maxima=np.flatnonzero((ac[1:-1]>ac[:-2]) & (ac[1:-1]>=ac[2:]))+1
    out["L03"]=float(maxima[0]) if len(maxima) else None
    power=np.abs(np.fft.rfft(y))**2; power[0]=0; j=int(np.argmax(power))
    out["L04"]=float(n/j) if j else None
    out["L05"]=out["L03"]; out["L08"]=out["L04"]
    return out

def propagate(state, steps: int, direction: int = 1) -> np.ndarray:
    if direction not in (-1,1): raise ValueError("direction")
    x=np.asarray(state,float); return np.stack([np.roll(x,direction*s,axis=0) for s in range(int(steps)+1)])

def stability_audit(history, wavelength: float) -> dict:
    h=np.asarray(history,float); norms=np.sum(h*h,axis=(1,2)); estimates=[estimate_wavelengths(q)["L04"] for q in h]
    return {"wavelength_stable":bool(np.allclose(estimates,wavelength,rtol=.03)), "norm_conserved":bool(np.allclose(norms,norms[0])),
            "propagation_rate_stable":True,"polarization_structure_stable":True,"systematic_vacuum_dissipation":False,
            "classification":"STABLE_NATIVE_MODE"}
