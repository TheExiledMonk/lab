"""Read-only DEV195 local force-balance diagnostics for DEV167 states."""
from __future__ import annotations

import numpy as np


def lattice_distance(shape: tuple[int, int, int], center: tuple[int, int, int]) -> np.ndarray:
    """Exact graph distance on the periodic N6 lattice."""
    axes = [np.minimum((np.arange(n) - c) % n, (c - np.arange(n)) % n) for n, c in zip(shape, center)]
    return axes[0][:, None, None] + axes[1][None, :, None] + axes[2][None, None, :]


def shell_partition(shape, center):
    d = lattice_distance(tuple(shape), tuple(center))
    return d, tuple(d == r for r in range(int(d.max()) + 1))


def region_l2(field: np.ndarray, mask: np.ndarray) -> float:
    return float(np.sqrt(np.sum(np.asarray(field)[mask] ** 2)))


def directed_outward_flux(power_flux: np.ndarray, distance: np.ndarray) -> np.ndarray:
    """Native bond-power flux, signed outward by graph-distance increase.

    Bonds that do not change distance (or cross a periodic distance tie) are
    retained as zero rather than being assigned an invented radial direction.
    """
    out = np.zeros_like(power_flux)
    for axis in range(3):
        delta = np.roll(distance, -1, axis=axis) - distance
        out[..., axis] = np.where(delta > 0, power_flux[..., axis],
                                  np.where(delta < 0, -power_flux[..., axis], 0.0))
    return out
