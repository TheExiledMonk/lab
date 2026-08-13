"""Exact DEV207 orientation and contrast helpers for DEV209."""
from __future__ import annotations

import numpy as np


def reflected_x(packet: np.ndarray) -> np.ndarray:
    """The DEV207 exact x-lattice reflection for a polar node vector."""
    q = np.flip(np.asarray(packet), axis=0).copy()
    q[..., 0] *= -1.0
    return q


def orientation_packets(displacement: np.ndarray, momentum: np.ndarray) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    return {"SAME": (np.asarray(displacement).copy(), np.asarray(momentum).copy()),
            "REVERSED": (reflected_x(displacement), reflected_x(momentum))}


def reflected_state(value: np.ndarray) -> np.ndarray:
    """Reflect a node-vector history; scalar/bond axes are retained unchanged."""
    q = np.flip(np.asarray(value), axis=-4).copy()
    if q.shape[-1] == 3:
        q[..., 0] *= -1.0
    return q
