"""DEV214 torque diagnostic using fixed initial provenance masks."""
from __future__ import annotations
import numpy as np


def support_torque(node_force: np.ndarray, mask: np.ndarray, center: np.ndarray) -> np.ndarray:
    grid = np.indices(mask.shape).transpose(1, 2, 3, 0).astype(float)
    return np.sum(np.cross(grid[mask] - center, np.asarray(node_force)[mask]), axis=0)
