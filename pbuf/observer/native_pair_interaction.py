"""Native, coefficient-free two packet diagnostics for DEV207."""
from __future__ import annotations
import numpy as np
from pbuf.excitation.native_vector_pair_dynamics import VectorPairState, net_force
from pbuf.observer.sequential_event_independence import inject, support_mask
from pbuf.observer.local_state_cross_event import four_state_cross_term

def reflected_x(packet: np.ndarray) -> np.ndarray:
    """Exact x reflection of a polar node-vector packet (not a sign choice)."""
    q=np.flip(np.asarray(packet), axis=0).copy(); q[...,0]*=-1
    return q

def orientation_packets(u: np.ndarray, p: np.ndarray) -> dict[str, tuple[np.ndarray,np.ndarray]]:
    return {'SAME': (u.copy(),p.copy()), 'REVERSED': (reflected_x(u),reflected_x(p))}

def evolve(state: VectorPairState, dt: float, steps: int, external, step) -> dict:
    rows={k:[] for k in ('displacement','momentum','force')}
    s=state
    for n in range(steps+1):
        rows['displacement'].append(s.displacement.copy()); rows['momentum'].append(s.momentum.copy()); rows['force'].append(net_force(s.displacement))
        if n<steps: s=step(s,dt,external)
    return {k:np.asarray(v) for k,v in rows.items()}

def centroid(mask: np.ndarray, displacement: np.ndarray) -> np.ndarray:
    # Fixed canonical bookkeeping support; no response-selected threshold.
    x=np.indices(mask.shape).transpose(1,2,3,0).astype(float)
    w=np.linalg.norm(displacement,axis=-1)*mask
    return (x*w[...,None]).sum(axis=(0,1,2))/w.sum() if w.sum() else x[mask].mean(axis=0)

def support_momentum(momentum: np.ndarray, mask: np.ndarray) -> np.ndarray:
    return momentum[:,mask,:].sum(axis=1)

def torque(momentum_cross: np.ndarray, displacement: np.ndarray, mask: np.ndarray) -> np.ndarray:
    grid=np.indices(mask.shape).transpose(1,2,3,0).astype(float)
    out=[]
    for p,u in zip(momentum_cross,displacement):
        c=centroid(mask,u); out.append(np.cross(grid[mask]-c,p[mask]).sum(axis=0))
    return np.asarray(out)

def four_state_trajectory(z0, za, zb, zab):
    """Exact state inclusion-exclusion and bond-force cross residual."""
    xcross=zab['displacement']-za['displacement']-zb['displacement']+z0['displacement']
    pcross=zab['momentum']-za['momentum']-zb['momentum']+z0['momentum']
    fcross=zab['force']-za['force']-zb['force']+z0['force']
    bonds=[four_state_cross_term(z0['displacement'][t],za['displacement'][t],zb['displacement'][t],zab['displacement'][t]) for t in range(len(xcross))]
    return {'displacement':xcross,'momentum':pcross,'force':fcross,
            'bond_force':np.asarray([b['force_cross'] for b in bonds]),
            'bond_strain':np.asarray([b['strain_cross'] for b in bonds]),
            'constitutive':np.asarray([b['constitutive_cross'] for b in bonds]),
            'geometric':np.asarray([b['geometric_cross'] for b in bonds])}

def canonical_support(u,p): return support_mask(u,p)
