"""Read-only DEV205 candidate construction from frozen DEV203--204 arrays."""
from __future__ import annotations

import numpy as np


def transverse_relative_motion(delta_relation: np.ndarray, direction: np.ndarray) -> np.ndarray:
    """Return the already-defined transverse part of every directed N6 relation."""
    n = np.asarray(direction, dtype=float)
    n = n / np.linalg.norm(n)
    return np.asarray(delta_relation, dtype=float) - np.sum(delta_relation * n, axis=-1, keepdims=True) * n


def axial_dual(antisymmetric: np.ndarray) -> np.ndarray:
    """Dual of a 3x3 antisymmetric native directional tensor.

    This is only the fixed orientation identity Q_i=1/2 eps_ijk A_jk;
    the explicit components avoid introducing any fitted convention.
    """
    a = np.asarray(antisymmetric, dtype=float)
    if a.shape[-2:] != (3, 3):
        raise ValueError("axial dual requires trailing 3x3 tensor axes")
    return np.stack((a[..., 2, 1], a[..., 0, 2], a[..., 1, 0]), axis=-1)
