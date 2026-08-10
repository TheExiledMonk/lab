"""Dev155: executable N6 operators for the established two-mode excitation.

Operators here are audited as mathematical candidates.  Registration does not
promote an operator to a physical evolution law.
"""
from __future__ import annotations
from dataclasses import dataclass, field
import numpy as np

N6_OFFSETS=((1,0,0),(-1,0,0),(0,1,0),(0,-1,0),(0,0,1),(0,0,-1))

@dataclass
class NativeExcitationN6State:
    values: np.ndarray
    progression_step: int=0
    history: list[np.ndarray]=field(default_factory=list)
    def __post_init__(self):
        self.values=np.asarray(self.values,dtype=np.float64)
        if self.values.ndim != 4 or self.values.shape[-1] != 2 or not np.isfinite(self.values).all():
            raise ValueError("values must be finite with shape (Nx, Ny, Nz, 2)")
        if min(self.values.shape[:3]) < 3: raise ValueError("each N6 axis needs at least three sites")
        self.history=[self.values.copy()] if not self.history else [np.asarray(x,float).copy() for x in self.history]
    @property
    def rank(self): return 2
    @property
    def topology(self): return "N6_3D_PERIODIC"

def shift(x,offset):
    """Periodic pull from the neighbor at offset; an exact permutation."""
    out=np.asarray(x,float)
    for axis,amount in enumerate(offset):
        if amount: out=np.roll(out,amount,axis=axis)
    return out

def neighbor_stack(x): return np.stack([shift(x,o) for o in N6_OFFSETS],axis=0)
def neighbor_sum(x): return np.sum(neighbor_stack(x),axis=0)
def neighbor_mean(x): return neighbor_sum(x)/6.0
def vector_laplacian(x): return neighbor_sum(x)-6.0*np.asarray(x,float)
def central_gradient(x):
    return np.stack([(shift(x,tuple(1 if a==axis else 0 for a in range(3)))-shift(x,tuple(-1 if a==axis else 0 for a in range(3))))/2.0 for axis in range(3)],axis=0)
def forward_differences(x): return neighbor_stack(x)-np.asarray(x,float)[None,...]

def transverse_rotation(x,quarter_turns=1):
    """Structural two-mode rotation, not an N6 allocation rule."""
    q=int(quarter_turns)%4
    rotation=np.linalg.matrix_power(np.array([[0.,-1.],[1.,0.]]),q)
    return np.asarray(x,float)@rotation.T

def execute_operator(operator_id,x,direction=(1,0,0)):
    """Execute every Dev148 spatial-operator candidate without promotion."""
    a=np.asarray(x,float)
    if operator_id=="O01": return forward_differences(a)
    if operator_id=="O02": return central_gradient(a)
    if operator_id=="O03": return neighbor_sum(a)-6*a
    if operator_id=="O04": return {"status":"STATE_PLACEMENT_MISMATCH","reason":"oriented circulation requires link/edge values; X is node based"}
    if operator_id=="O05": return {"status":"STATE_RANK_MISMATCH","reason":"3D curl requires a spatial 3-vector or defined link cochain; X has two internal components"}
    if operator_id=="O06": return np.stack([transverse_rotation(shift(a,o)) for o in N6_OFFSETS],axis=0)
    if operator_id=="O07": return np.stack([shift(a,o)-shift(a,tuple(-v for v in o)) for o in N6_OFFSETS[::2]],axis=0)
    if operator_id=="O08": return vector_laplacian(a)
    if operator_id=="O09": return {"status":"OPERATOR_UNDERDEFINED","reason":"no native map from spatial divergence/circulation outputs back to the two internal modes"}
    if operator_id=="O10": return neighbor_stack(a)
    raise ValueError(operator_id)

def propagate_directional(state,steps,direction):
    """N6 direction-permutation control; not selected as the unique law."""
    if tuple(direction) not in N6_OFFSETS: raise ValueError("direction must be one N6 offset")
    for _ in range(int(steps)):
        state.values=shift(state.values,tuple(direction)); state.progression_step+=1; state.history.append(state.values.copy())
    return state

def gaussian_packet(shape=(24,24,24),center=None,width=2.5,amplitude=1.,polarization=(1.,0.)):
    shape=tuple(map(int,shape)); center=np.array([(n-1)/2 for n in shape]) if center is None else np.asarray(center,float)
    grids=np.meshgrid(*[np.arange(n) for n in shape],indexing="ij"); r2=sum((g-center[i])**2 for i,g in enumerate(grids))
    p=np.asarray(polarization,float); p/=np.linalg.norm(p)
    return float(amplitude)*np.exp(-.5*r2/float(width)**2)[...,None]*p

def quadratic_norm(x): return float(np.sum(np.asarray(x,float)**2))
def density(x): return np.sum(np.asarray(x,float)**2,axis=-1)
def centroid(x):
    q=density(x); total=q.sum(); grids=np.meshgrid(*[np.arange(n) for n in q.shape],indexing="ij")
    return np.array([np.sum(g*q)/total for g in grids]) if total else np.full(3,np.nan)

def operator_registry():
    names=("nearest-neighbor differences","central differences","divergence-like node sum","oriented link circulation","plaquette curl-like operator","transverse neighbor rotation","antisymmetric pair exchange","six-neighbor vector Laplacian","coupled divergence/circulation","existing medium N6 topology")
    return [{"id":f"O{i:02d}","name":n,"executed":True} for i,n in enumerate(names,1)]
