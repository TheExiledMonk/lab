"""Read-only relation classifiers; this module never assigns a phase variable."""
from __future__ import annotations
import numpy as np

def normalized_dot(a: np.ndarray, b: np.ndarray) -> float:
    d = float(np.linalg.norm(a) * np.linalg.norm(b))
    return float(np.vdot(a, b).real / d) if d else float("nan")

def neighbor_correlations(features: np.ndarray, pairs: np.ndarray) -> np.ndarray:
    return np.array([normalized_dot(features[:, i] - features[:, i].mean(0),
                                    features[:, j] - features[:, j].mean(0)) for i, j in pairs])
