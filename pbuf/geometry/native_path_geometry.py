"""Topological-path controls which do not presume a loaded embedding."""
from __future__ import annotations

import numpy as np


def straight_path(start, axis: int, steps: int) -> np.ndarray:
    start = np.asarray(start, dtype=float)
    if start.shape != (3,) or axis not in (0, 1, 2) or steps < 1:
        raise ValueError("start is length 3, axis is 0..2, and steps is positive")
    points = np.repeat(start[None, :], steps + 1, axis=0)
    points[:, axis] += np.arange(steps + 1)
    return points


def path_diagnostics(points: np.ndarray) -> dict:
    points = np.asarray(points, dtype=float)
    bonds = np.diff(points, axis=0)
    lengths = np.linalg.norm(bonds, axis=1)
    tangents = bonds / lengths[:, None]
    dots = np.sum(tangents[:-1] * tangents[1:], axis=1)
    turns = np.arccos(np.clip(dots, -1.0, 1.0))
    return {"point_count": len(points), "bond_lengths": lengths.tolist(),
            "turning_radians": turns.tolist(),
            "maximum_turning_radians": float(turns.max(initial=0.0))}
