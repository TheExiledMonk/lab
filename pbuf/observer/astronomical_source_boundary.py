"""DEV189 guardrails for the astronomy-to-native source boundary.

This module deliberately supplies no astronomical source reconstruction.  It
only implements exact, coefficient-free checks used by the DEV189 audit.
"""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np


@dataclass(frozen=True)
class SimilarityBridge:
    """A pre-fixed coordinate convention: isotropic scale, rotation, origin."""
    scale: float
    rotation: np.ndarray
    origin: np.ndarray

    def __post_init__(self) -> None:
        r = np.asarray(self.rotation, dtype=np.float64)
        o = np.asarray(self.origin, dtype=np.float64)
        if self.scale <= 0 or r.shape != (2, 2) or o.shape != (2,):
            raise ValueError("invalid similarity bridge dimensions")
        if not np.allclose(r.T @ r, np.eye(2), rtol=0, atol=1e-12):
            raise ValueError("rotation must be orthogonal")
        object.__setattr__(self, "rotation", r)
        object.__setattr__(self, "origin", o)

    @property
    def matrix(self) -> np.ndarray:
        return self.scale * self.rotation

    def transform(self, points: np.ndarray) -> np.ndarray:
        return np.asarray(points, dtype=np.float64) @ self.matrix.T + self.origin

    def injects_anisotropy(self) -> bool:
        return not np.allclose(self.matrix.T @ self.matrix,
                               self.scale ** 2 * np.eye(2), rtol=0, atol=1e-12)


def normalized_second_moment(points: np.ndarray, weights: np.ndarray) -> np.ndarray:
    """Dimensionless central second moment for synthetic coordinate controls."""
    q, w = np.asarray(points, float), np.asarray(weights, float)
    centroid = (q * w[:, None]).sum(axis=0) / w.sum()
    cov = ((q - centroid).T * w) @ (q - centroid) / w.sum()
    return cov / np.trace(cov)
