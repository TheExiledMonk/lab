"""Map relative distributed source content onto the stationary N6 constraint."""
from __future__ import annotations
import numpy as np


def distributed_source_imposed_excursion(source: np.ndarray) -> np.ndarray:
    """Linear extension of Dev159's local contact plus zero-mode projection."""
    forcing = np.asarray(source, float).copy()
    if forcing.ndim != 3 or not np.isfinite(forcing).all():
        raise ValueError("source must be a finite three-dimensional array")
    forcing -= forcing.mean()
    return forcing
