"""Native DEV204 relational/stress tensor helpers."""
from __future__ import annotations

import numpy as np

from pbuf.observer.relational_wave_content import native_symmetric_antisymmetric


def force_response_tensors(delta_force: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Apply the unchanged DEV203 opposite-bond tensor construction to dF."""
    return native_symmetric_antisymmetric(delta_force)


def l2(value: np.ndarray) -> float:
    return float(np.linalg.norm(value))
