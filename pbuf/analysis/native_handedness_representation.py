"""Coefficient-free representations of a 3D antisymmetric relational tensor."""
from __future__ import annotations

import numpy as np


def axial_dual(a: np.ndarray) -> np.ndarray:
    """Return omega=(A_yz,-A_xz,A_xy), with epsilon_xyz=+1."""
    a = np.asarray(a)
    return np.stack((a[..., 1, 2], -a[..., 0, 2], a[..., 0, 1]), axis=-1)


def frobenius_neighbor_relation(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Canonical Euclidean A:B contraction; no normalization or threshold."""
    return np.einsum("...ij,...ij->...", a, b)


def axial_neighbor_relation(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """The equivalent axial-vector relation under the fixed orientation."""
    return np.einsum("...i,...i->...", axial_dual(a), axial_dual(b))
