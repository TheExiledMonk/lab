"""DEV216 fixed-region N6 bond-cut pair-force observer.

The functions here are observer-only: they materialize no new mechanics and
select no region from trajectory data.
"""
from __future__ import annotations
import numpy as np
from pbuf.observer.native_region_momentum_balance import bond_cut


def n6_ball(shape: tuple[int, int, int], center: np.ndarray, radius: int = 2) -> np.ndarray:
    """Exact periodic N6 graph ball used by DEV215 (radius two)."""
    points = np.indices(shape).reshape(3, -1).T
    shape_a = np.asarray(shape)
    delta = np.abs((points - np.asarray(center) + shape_a // 2) % shape_a - shape_a // 2)
    mask = np.zeros(shape, dtype=bool)
    mask[tuple(points[np.sum(delta, axis=1) <= radius].T)] = True
    return mask


def interaction_cut(u0, ua, ub, uab, mask: np.ndarray) -> np.ndarray:
    """Exact four-state observer residual at one fixed boundary."""
    return bond_cut(uab, mask) - bond_cut(ua, mask) - bond_cut(ub, mask) + bond_cut(u0, mask)


def bond_force_cut(force: np.ndarray, mask: np.ndarray) -> np.ndarray:
    """Boundary sum for an already-materialized positive-bond force array."""
    total = np.zeros(3)
    for axis in range(3):
        inside, other = mask, np.roll(mask, -1, axis=axis)
        total += np.asarray(force)[..., axis, :][inside & ~other].sum(axis=0)
        total -= np.asarray(force)[..., axis, :][~inside & other].sum(axis=0)
    return total


def cut_reciprocity(displacement: np.ndarray, mask: np.ndarray) -> float:
    """The complementary cut must be the negative of the given cut."""
    return float(np.max(np.abs(bond_cut(displacement, mask) + bond_cut(displacement, ~mask))))
