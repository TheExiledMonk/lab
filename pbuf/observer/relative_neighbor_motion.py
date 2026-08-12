"""Read-only DEV203 differences of the canonical directed N6 state."""
from __future__ import annotations
import numpy as np

from pbuf.observer.native_n6_field import n6_field


def relative_neighbor_motion(u0: np.ndarray, p0: np.ndarray,
                             u1: np.ndarray, p1: np.ndarray) -> dict[str, np.ndarray]:
    """Materialize exact one-update changes; this does not evolve a state."""
    a, b = n6_field(u0, p0), n6_field(u1, p1)
    return {
        "relation": a["relation"], "next_relation": b["relation"],
        "delta_relation": b["relation"] - a["relation"],
        "strain": a["strain"], "next_strain": b["strain"],
        "delta_strain": b["strain"] - a["strain"],
        "unit": a["unit"], "next_unit": b["unit"],
        "delta_unit": b["unit"] - a["unit"],
        "force": a["force"], "next_force": b["force"],
        "delta_force": b["force"] - a["force"],
        "momentum_difference": directed_momentum_difference(p0),
    }


def directed_momentum_difference(momentum: np.ndarray) -> np.ndarray:
    """p_b-p_a in canonical (+x,-x,+y,-y,+z,-z) order."""
    p = np.asarray(momentum, dtype=float)
    parts = []
    for axis in range(3):
        plus = np.roll(p, -1, axis=axis) - p
        minus = np.roll(p, 1, axis=axis) - p
        parts.extend((plus, minus))
    return np.stack(parts, axis=-2)


def motion_state(u0: np.ndarray, p0: np.ndarray, u1: np.ndarray, p1: np.ndarray) -> dict[str, np.ndarray]:
    """Return the DEV203 W_ab descriptor with correctly directed momentum."""
    out = relative_neighbor_motion(u0, p0, u1, p1)
    out["momentum_difference"] = directed_momentum_difference(p0)
    return out
