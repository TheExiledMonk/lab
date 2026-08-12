"""Homogeneous N6 scale-deformation diagnostics for DEV208."""
from __future__ import annotations

import numpy as np

N6_DIRECTIONS = ("+x", "-x", "+y", "-y", "+z", "-z")


def homogeneous_relations(scale: float) -> np.ndarray:
    """Six isotropic relation vectors after r_ab -> scale r_ab."""
    if not np.isfinite(scale) or scale <= 0:
        raise ValueError("scale must be finite and positive")
    return scale * np.array(((1, 0, 0), (-1, 0, 0), (0, 1, 0),
                             (0, -1, 0), (0, 0, 1), (0, 0, -1)), dtype=float)


def extension(scale: float) -> float:
    """Extension of a unit-reference native relation."""
    if not np.isfinite(scale):
        raise ValueError("scale must be finite")
    return float(scale - 1.0)


def admissible(scale: float) -> bool:
    return abs(extension(scale)) < 1.0
