"""Native excitation state definitions independent of rays and effective EM fields."""
from __future__ import annotations
from dataclasses import dataclass, field
import numpy as np

NAMES=("signed node scalar","signed link scalar","node vector","oriented link vector",
 "two-component node state","two-component link state","transverse two-component state",
 "antisymmetric pair state","coupled scalar pair","coupled vector pair","local orientation + magnitude",
 "link orientation + magnitude","excitation/response conjugate pair","displacement-like + circulation-like pair",
 "strain-like dynamic tensor perturbation","antisymmetric tensor-like excitation","three-axis pair state",
 "six-component native mode state","minimal state inferred from EM rank","no tested state sufficient")

def state_registry():
    rows=[]
    for i,name in enumerate(NAMES,1):
        if i in (1,2,9): status="FAILS_TWO_TRANSVERSE_DOF"
        elif i in (3,4,10,17,18): status="OVERCOMPLETE"
        elif i in (5,6,7,19): status="STRUCTURALLY_SUPPORTED"
        elif i in (8,11,12,13,14,15,16): status="DERIVABLE"
        else: status="NOT_APPLICABLE"
        rows.append({"id":f"X{i:02d}","name":name,"status":status,"attempted":True})
    return rows

def localized_packet(sites=128, center=None, width=7., amplitude=1., polarization=(1.,0.)):
    """Localized signed rank-2 packet on native spatial support."""
    if sites < 8 or width <= 0: raise ValueError("invalid packet support")
    center=(sites//4 if center is None else float(center)); x=np.arange(sites,dtype=float)
    p=np.asarray(polarization,float)
    if p.shape != (2,) or not np.isfinite(p).all() or np.linalg.norm(p)==0:
        raise ValueError("polarization must be a finite nonzero two-vector")
    p=p/np.linalg.norm(p)
    return float(amplitude)*np.exp(-.5*((x-center)/float(width))**2)[:,None]*p[None,:]

@dataclass
class NativeExcitationState:
    values: np.ndarray
    progression_step: int = 0
    history: list[np.ndarray] = field(default_factory=list)
    def __post_init__(self):
        self.values=np.asarray(self.values,dtype=np.float64)
        if self.values.ndim != 2 or self.values.shape[1] != 2 or not np.isfinite(self.values).all():
            raise ValueError("values must have shape (sites, 2) and be finite")
        self.history=[self.values.copy()] if not self.history else [np.asarray(x,float).copy() for x in self.history]
    @property
    def rank(self): return 2
    @property
    def location(self): return "NODE_STATE"

