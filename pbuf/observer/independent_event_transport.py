"""DEV194 event-ensemble algebra, intentionally separate from state evolution."""
from __future__ import annotations
import numpy as np

def expected_accumulation(weight_kernel: np.ndarray, event_weights: np.ndarray) -> np.ndarray:
    """Finite expected receipt accumulation ``K lambda``; no packet superposition."""
    k, lam = np.asarray(weight_kernel, float), np.asarray(event_weights, float)
    if k.ndim != 2 or lam.shape != (k.shape[1],) or np.any(lam < 0):
        raise ValueError("event weights must be nonnegative and match launch columns")
    return k @ lam

def repeated_response_accumulation(columns: np.ndarray, counts: np.ndarray) -> np.ndarray:
    """Exact detector-side sum of independently supplied response columns."""
    return expected_accumulation(columns, counts)
