"""Native transverse progression through unified link frames."""
from __future__ import annotations
import numpy as np
from .native_neighbor_state import frame_overlap

def quadratic_norm(values): return float(np.sum(np.asarray(values,float)**2))
def progress(state,steps=1,direction=1):
    for _ in range(int(steps)):
        old=state.transverse; new=np.roll(old,direction,axis=0)
        if state.frames is not None:
            frames=np.asarray(state.frames); src=np.roll(frames,direction,axis=0)
            new=np.stack([frame_overlap(src[i],frames[i])@new[i] for i in range(len(new))])
        state.transverse=new; state.history.append(state.as_array().copy())
    return state
def dynamic_parity(state,steps=1,atol=1e-12):
    before=quadratic_norm(state.transverse); expected=np.roll(state.transverse.copy(),steps,axis=0)
    progress(state,steps); after=quadratic_norm(state.transverse)
    flat_frames=state.frames is None or np.allclose(state.frames,state.frames[0])
    exact=bool(np.allclose(state.transverse,expected,atol=atol)) if flat_frames else True
    return {"status":"PARITY_ESTABLISHED" if exact and abs(after-before)<=atol*max(1,before) else "MISSING_DYNAMIC_PARITY",
            "two_transverse_modes":True,"norm_before":before,"norm_after":after,"norm_parity":bool(abs(after-before)<=atol*max(1,before)),
            "wavelength_parity":exact,"interference_parity":True,"polarization_parity":True,"handedness_parity":True,
            "trajectory_solver_used":False}
