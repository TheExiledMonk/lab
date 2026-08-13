"""Predeclared fixed-plane receiver geometry for DEV209."""
from __future__ import annotations

import numpy as np


def receiver_planes(shape: tuple[int, int, int], transverse_mask: np.ndarray) -> dict[str, dict[str, np.ndarray]]:
    """Fixed x planes 4, 5, 6 and their exact x-reflected counterparts."""
    if shape[0] != 11:
        raise ValueError("DEV209 reuses the fixed canonical 11-node DEV182 lattice")
    yz = np.asarray(transverse_mask, dtype=bool)
    if yz.shape != shape[1:]:
        raise ValueError("transverse mask shape must match receiver plane")
    result = {}
    for label, x in zip(("R1", "R2", "R3"), (4, 5, 6)):
        same = np.zeros(shape, dtype=bool); same[x] = yz
        reversed_mask = np.zeros(shape, dtype=bool); reversed_mask[shape[0] - 1 - x] = yz
        result[label] = {"SAME": same, "REVERSED": reversed_mask,
                         "same_plane_x": x, "reversed_plane_x": shape[0] - 1 - x}
    return result
