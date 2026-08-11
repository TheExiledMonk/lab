"""Threshold-independent spatial-support diagnostics for Dev157."""
from __future__ import annotations

import numpy as np


def periodic_distance_grids(shape: tuple[int, int, int], center: tuple[int, int, int]):
    axes = []
    for n, c in zip(shape, center):
        x = np.arange(n) - c
        axes.append((x + n // 2) % n - n // 2)
    return np.meshgrid(*axes, indexing="ij")


def support_metrics(density: np.ndarray, center: tuple[int, int, int] | None = None) -> dict:
    """RMS radius and participation scale for an explicitly nonnegative density."""
    rho = np.asarray(density, float)
    if rho.ndim != 3 or np.any(rho < -1e-15) or not np.isfinite(rho).all():
        raise ValueError("density must be a finite nonnegative three-dimensional array")
    rho = np.maximum(rho, 0.0); total = float(rho.sum())
    if total == 0.0:
        return {"rms_radius": None, "participation_volume": 0.0,
                "participation_linear_scale": 0.0}
    center = center or tuple(n // 2 for n in rho.shape)
    grids = periodic_distance_grids(rho.shape, center)
    radius2 = sum(x.astype(float) ** 2 for x in grids)
    p = rho / total
    volume = float(1.0 / np.sum(p*p))
    return {"rms_radius": float(np.sqrt(np.sum(radius2 * p))),
            "participation_volume": volume,
            "participation_linear_scale": float(volume ** (1.0/3.0))}


def signed_periodic_correlation(field: np.ndarray) -> np.ndarray:
    """Circular autocorrelation, normalized so C(0)=1 when nonzero."""
    q = np.asarray(field, float)
    transformed = np.fft.fftn(q)
    c = np.fft.ifftn(np.abs(transformed) ** 2).real
    return c / c[(0, 0, 0)] if c[(0, 0, 0)] else c


def radial_correlation(field: np.ndarray) -> dict:
    c = signed_periodic_correlation(field)
    grids = periodic_distance_grids(c.shape, (0, 0, 0))
    radius = np.rint(np.sqrt(sum(x.astype(float) ** 2 for x in grids))).astype(int)
    count = np.bincount(radius.ravel())
    values = np.bincount(radius.ravel(), weights=c.ravel()) / np.maximum(count, 1)
    zeros = np.flatnonzero(values <= 0)
    below = np.flatnonzero(values <= np.exp(-1))
    return {"radius": list(map(int, range(len(values)))), "correlation": values.tolist(),
            "first_zero_crossing": int(zeros[0]) if len(zeros) else None,
            "first_one_over_e_crossing": int(below[0]) if len(below) else None,
            "integral_scale_positive_lobes": float(np.sum(np.maximum(values, 0.0))),
            "measure":"signed q autocorrelation; no fabricated complex phase"}
