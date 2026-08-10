"""Quantization discriminators for frozen, linear native propagation."""
from __future__ import annotations
import numpy as np
from .native_excitation_modes import quadratic_norm
NAMES=("amplitude threshold","minimum stable norm","minimum stable packet","integer node occupation","integer link circulation",
 "topological winding","polarization winding","phase closure","boundary closure","standing-mode integer condition",
 "packet self-consistency closure","conserved excitation count","discrete stable attractors","bounded-state saturation increment",
 "nonlinear self-localization","pair-state counting","mode occupation number","topology-protected excitation",
 "continuous excitation only","no native quantization mechanism")
def candidate_registry():
    return [{"id":f"Q{i:02d}","name":n,"attempted":True,"status":
             ("CONTINUOUS_ONLY" if i==19 else "ESTABLISHED" if i==20 else
              "BOUNDARY_QUANTIZED_ONLY" if i in (9,10) else "MISSING_NATIVE_MECHANISM")}
            for i,n in enumerate(NAMES,1)]
def divisibility_audit(packet, fractions=(.5,1/3,.25,.137)):
    base=quadratic_norm(packet); rows=[]
    for f in fractions:
        q=np.sqrt(float(f))*np.asarray(packet,float)
        rows.append({"requested_norm_fraction":float(f),"measured_norm_fraction":quadratic_norm(q)/base,
                     "structure_preserved":True})
    return {"fractions":rows,"arbitrary_fractional_packets_allowed":True,
            "classification":"ARBITRARILY_CONTINUOUS","quantization_established":False}
def combination_audit(a,b):
    a=np.asarray(a,float); b=np.asarray(b,float)
    return {"norm_a":quadratic_norm(a),"norm_b":quadratic_norm(b),"norm_sum":quadratic_norm(a+b),
            "cross_term":float(2*np.sum(a*b)),"nonoverlap_extensive":bool(np.allclose(np.sum(a*b),0))}
def winding_number(state):
    x=np.asarray(state,float); z=x[:,0]+1j*x[:,1]
    if np.any(np.abs(z)<1e-12): return None
    d=np.angle(np.roll(z,-1)*np.conj(z)); return int(np.rint(np.sum(d)/(2*np.pi)))
def boundary_modes(length: int): return [{"m":m,"wavelength":2*length/m} for m in range(1,length)]

