"""Read-only ordered-N6 relational-pattern mismatch observers for DEV223."""
from __future__ import annotations

import numpy as np

from pbuf.observer.native_n6_field import n6_field


N6_ORDER = ("+x", "-x", "+y", "-y", "+z", "-z")


def n6_signature(displacement: np.ndarray) -> np.ndarray:
    """The canonical signed strain signature, retaining directed order."""
    return n6_field(displacement, np.zeros_like(displacement))["strain"]


def neighbor_mismatch(signature: np.ndarray) -> np.ndarray:
    """D_ab=S_b-S_a for every directed, fixed-topology N6 edge."""
    result = np.empty(signature.shape + (6,), dtype=signature.dtype)
    for direction, (axis, shift) in enumerate(((0, -1), (0, 1), (1, -1), (1, 1), (2, -1), (2, 1))):
        result[..., direction, :] = np.roll(signature, shift, axis=axis) - signature
    return result


def longitudinal_component_profile(mismatch: np.ndarray) -> np.ndarray:
    """Unweighted componentwise sum over each canonical x plane."""
    return mismatch.sum(axis=(1, 2))
