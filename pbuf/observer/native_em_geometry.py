"""Coefficient-free discrete geometry diagnostics used by DEV205."""
from __future__ import annotations

import numpy as np


def project(vector: np.ndarray, direction: np.ndarray) -> np.ndarray:
    n = np.asarray(direction, dtype=float); n /= np.linalg.norm(n)
    return np.sum(vector * n, axis=-1)


def pair_geometry(polar: np.ndarray, axial: np.ndarray) -> dict[str, np.ndarray]:
    """Evaluate every directed-bond polar value against its node axial value."""
    q = np.expand_dims(np.asarray(axial, float), axis=-2)
    p = np.asarray(polar, float)
    return {"dot": np.sum(p * q, axis=-1), "cross": np.cross(p, q)}


def n6_directional_difference(vector: np.ndarray) -> np.ndarray:
    """Exact central directional difference on the existing periodic N6 lattice."""
    v = np.asarray(vector, float)
    return np.stack([(np.roll(v, -1, axis=i) - np.roll(v, 1, axis=i)) / 2 for i in range(3)], axis=-2)


def native_source_diagnostic(vector: np.ndarray) -> np.ndarray:
    """Native source/sink readout: diagonal contraction of N6 directional differences."""
    d = n6_directional_difference(vector)
    return np.trace(d, axis1=-2, axis2=-1)
