"""Deterministic DEV205 classifications; no EM coefficient or fit is used."""
from __future__ import annotations

import numpy as np


def exact_transversality(values: np.ndarray, direction: np.ndarray) -> str:
    x = np.sum(np.asarray(values) * np.asarray(direction), axis=-1)
    return "EXACT" if np.array_equal(x, np.zeros_like(x)) else "MIXED"


def orthogonality(dot: np.ndarray) -> str:
    return "EXACT" if np.array_equal(dot, np.zeros_like(dot)) else "MIXED"


def handedness(cross: np.ndarray, direction: np.ndarray) -> str:
    signed = np.sum(cross * np.asarray(direction), axis=-1)
    if np.array_equal(cross, np.zeros_like(cross)):
        return "NOT_DEFINED"
    nonzero = signed[signed != 0]
    if nonzero.size == 0:
        return "MIXED"
    if np.all(nonzero > 0): return "ALIGNED_WITH_PROPAGATION"
    if np.all(nonzero < 0): return "ANTI_ALIGNED"
    return "BIDIRECTIONAL"


def source_structure(source: np.ndarray) -> str:
    return "SOURCE_FREE" if np.array_equal(source, np.zeros_like(source)) else "DISTRIBUTED"
