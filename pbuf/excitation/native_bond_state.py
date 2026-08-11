"""Dev156 oriented N6 bond representation and relational diagnostics."""
from __future__ import annotations
import numpy as np

POSITIVE_AXES = ((1, 0, 0), (0, 1, 0), (0, 0, 1))
N6_OFFSETS = ((1,0,0),(-1,0,0),(0,1,0),(0,-1,0),(0,0,1),(0,0,-1))


def neighbor(q: np.ndarray, offset: tuple[int, int, int]) -> np.ndarray:
    out = np.asarray(q, dtype=float)
    for axis, amount in enumerate(offset):
        if amount:
            out = np.roll(out, -amount, axis=axis)
    return out


def relational_differences(q: np.ndarray) -> np.ndarray:
    """Six directed d[a,offset]=q[a+offset]-q[a]."""
    q = np.asarray(q, dtype=float)
    return np.stack([neighbor(q, o) - q for o in N6_OFFSETS], axis=-1)


def positive_gradient(q: np.ndarray) -> np.ndarray:
    return np.stack([neighbor(q, o) - q for o in POSITIVE_AXES], axis=-1)


def gradient_adjoint(bonds: np.ndarray) -> np.ndarray:
    """Adjoint of positive_gradient for the periodic lattice."""
    b = np.asarray(bonds, dtype=float)
    if b.ndim != 4 or b.shape[-1] != 3:
        raise ValueError("bonds must have shape (Nx, Ny, Nz, 3)")
    return sum(np.roll(b[..., axis], 1, axis=axis) - b[..., axis] for axis in range(3))


def relational_imbalance(q: np.ndarray) -> np.ndarray:
    return relational_differences(q).sum(axis=-1)


def axis_antisymmetry(q: np.ndarray) -> np.ndarray:
    d = relational_differences(q)
    return np.stack((d[...,0]-d[...,1], d[...,2]-d[...,3], d[...,4]-d[...,5]), axis=-1)


def full_directed_bonds(positive_bonds: np.ndarray) -> np.ndarray:
    """Materialize both orientations; negative links obey tau_ab=-tau_ba."""
    b = np.asarray(positive_bonds, dtype=float)
    return np.stack((b[...,0], -np.roll(b[...,0],1,0), b[...,1],
                     -np.roll(b[...,1],1,1), b[...,2],-np.roll(b[...,2],1,2)), axis=-1)


def antisymmetry_error(positive_bonds: np.ndarray) -> float:
    d = full_directed_bonds(positive_bonds)
    errors = []
    for plus, minus, axis in ((0,1,0),(2,3,1),(4,5,2)):
        errors.append(np.max(np.abs(d[...,plus] + np.roll(d[...,minus],-1,axis=axis))))
    return float(max(errors))
