"""Synthetic, scale-free source controls for Dev138.

Truth packages deliberately have a different type from reconstruction-facing
observations.  In particular, :class:`BlindReceivedState` has no truth field.
"""
from __future__ import annotations
from dataclasses import dataclass
import hashlib
import numpy as np

MORPHOLOGIES=("compact_isotropic_blob","elliptical_source","asymmetric_single_component",
 "two_component_source","ring_like_source","irregular_connected_source","elongated_source","sparse_multi_knot_source")
SOURCE_SIZES=(.25,.5,1.,1.5,2.)
DEPTH_OFFSETS=(.25,.5,.75,1.,1.5,2.,3.,4.,6.,8.)
LENS_FAMILIES=("symmetric_single","compact_strong","diffuse","asymmetric","two_source","offset_two_source","weak_response","strong_unsaturated")

@dataclass(frozen=True)
class BlindReceivedState:
    control_id: str
    positions: np.ndarray
    directions: np.ndarray
    event_uids: np.ndarray
    bundle_ids: np.ndarray|None=None
    trajectory_history: tuple=()

@dataclass(frozen=True)
class SourceTruth:
    control_id: str
    depth_native: float
    positions_native: np.ndarray
    lens_depth_native: float
    lens_radius_native: float

def deterministic_split(control_id: str)->str:
    """Fixed 60/40 split without depending on Python's salted hash."""
    return "TRAIN" if int(hashlib.sha256(control_id.encode()).hexdigest()[:8],16)%10<6 else "VALIDATION"

def source_cloud(morphology="compact_isotropic_blob", radius=1., n=64):
    """Return a deterministic generic event cloud; no galaxy templates."""
    t=2*np.pi*(np.arange(n)+.5)/n
    r=radius*np.sqrt((np.arange(n)+.5)/n)
    x=np.c_[r*np.cos(t),r*np.sin(t)]
    if morphology=="elliptical_source": x[:,1]*=.55
    elif morphology=="asymmetric_single_component": x[:,0]+=0.25*radius*(x[:,1]/radius)**2
    elif morphology=="two_component_source": x[:,0]+=np.where(np.arange(n)%2, .55,-.55)*radius
    elif morphology=="ring_like_source": x=radius*np.c_[np.cos(t),np.sin(t)]
    elif morphology=="irregular_connected_source": x*=1+.2*np.sin(3*t)[:,None]
    elif morphology=="elongated_source": x[:,1]*=.25
    elif morphology=="sparse_multi_knot_source": x=np.round(2*x/radius)*radius/2
    return x

def blind_package(truth: SourceTruth, forward):
    """Apply a frozen caller-supplied forward operator and erase truth."""
    p,d=forward(truth.positions_native,truth.depth_native)
    n=len(p)
    return BlindReceivedState(truth.control_id,np.asarray(p,float),np.asarray(d,float),
       np.array([f"{truth.control_id}:{i}" for i in range(n)]),np.arange(n,dtype=int))

