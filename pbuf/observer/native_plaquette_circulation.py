"""Topology-defined elementary N6 square circulation readout."""
from __future__ import annotations
import numpy as np

def elementary_plaquettes(shape: tuple[int, int, int], center: tuple[int, int, int]) -> np.ndarray:
    # three coordinate planes through the source: no post-hoc loop choice
    loops=[]; c=np.array(center)
    for a,b in ((0,1),(0,2),(1,2)):
        e=np.eye(3,dtype=int); loops.append(np.array([c,c+e[a],c+e[a]+e[b],c+e[b]]) % np.array(shape))
    return np.array(loops)

def circulation(positive_flux: np.ndarray, loops: np.ndarray) -> np.ndarray:
    out=[]
    for loop in loops:
        value=[]
        for t in range(len(positive_flux)):
            total=0.0
            for x,y in zip(loop, np.roll(loop,-1,axis=0)):
                d=(y-x) % np.array(positive_flux.shape[1:4]); axis=int(np.flatnonzero(d==1)[0]) if np.any(d==1) else int(np.flatnonzero(d==positive_flux.shape[1:4][np.flatnonzero(d!=0)[0]]-1)[0])
                sign=1.0 if y[axis] == (x[axis]+1)%positive_flux.shape[axis+1] else -1.0
                total += sign*positive_flux[t, *x, axis]
            value.append(total)
        out.append(value)
    return np.asarray(out).T
