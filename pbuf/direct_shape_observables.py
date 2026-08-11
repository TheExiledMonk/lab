"""Coefficient-free DEV176 spin-2 observables from frozen receipt state.

This module is deliberately upstream of the observer and contains no
astrometry fitting, shear calibration, propagation, or observational logic.
"""
from __future__ import annotations

import numpy as np


def project_to_screen(vectors, u_axis, v_axis, *, native_xyz=True):
    """Project native vectors onto the frozen observer screen.

    DEV168 stores positions as ``(x, y, z)`` whereas the frozen observer basis
    is expressed in its historical ``(y, z, x)`` ordering.
    """
    a = np.asarray(vectors, dtype=float)
    if a.shape[-1] != 3:
        raise ValueError("vectors must end in three native components")
    if native_xyz:
        a = a[..., [1, 2, 0]]
    return np.stack((a @ np.asarray(u_axis, float), a @ np.asarray(v_axis, float)), axis=-1)


def weighted_second_moment_tensor(points, weights=None):
    """Return a symmetric central 2D moment tensor, centroid and support."""
    p = np.asarray(points, dtype=float)
    if p.size == 0:
        p = p.reshape(0, 2)
    if p.ndim != 2 or p.shape[1] != 2:
        raise ValueError("points must have shape (N,2)")
    w = np.ones(len(p), dtype=float) if weights is None else np.asarray(weights, dtype=float)
    good = np.isfinite(p).all(1) & np.isfinite(w) & (w >= 0)
    p, w = p[good], w[good]
    total = float(w.sum())
    if not len(p) or total <= 0:
        return np.full((2, 2), np.nan), np.full(2, np.nan), 0
    centre = (p * w[:, None]).sum(0) / total
    delta = p - centre
    return (delta.T * w) @ delta / total, centre, len(p)


def quadrupole_tensor(vectors, weights=None):
    """Return the raw 2D directional tensor (no mean subtraction)."""
    v = np.asarray(vectors, dtype=float)
    if v.ndim != 2 or v.shape[1] != 2:
        raise ValueError("vectors must have shape (N,2)")
    w = np.ones(len(v), dtype=float) if weights is None else np.asarray(weights, dtype=float)
    good = np.isfinite(v).all(1) & np.isfinite(w) & (w >= 0)
    v, w = v[good], w[good]
    total = float(w.sum())
    if not len(v) or total <= 0:
        return np.full((2, 2), np.nan), 0
    return (v.T * w) @ v / total, len(v)


def spin2_from_tensor(tensor):
    """Trace-normalized ``(e1,e2)``; undefined trace is returned as NaN."""
    t = np.asarray(tensor, dtype=float)
    if t.shape != (2, 2):
        raise ValueError("tensor must have shape (2,2)")
    trace = float(np.trace(t))
    if not np.isfinite(trace) or abs(trace) <= np.finfo(float).eps:
        return np.array([np.nan, np.nan])
    return np.array([(t[0, 0] - t[1, 1]) / trace, 2.0 * t[0, 1] / trace])


def local_deformation_tensor(source_points, received_points):
    """Least-squares source-to-receipt Jacobian and symmetric traceless part."""
    s, r = np.asarray(source_points, float), np.asarray(received_points, float)
    if s.shape != r.shape or s.ndim != 2 or s.shape[1] != 2:
        raise ValueError("source and received points must both have shape (N,2)")
    good = np.isfinite(s).all(1) & np.isfinite(r).all(1)
    s, r = s[good], r[good]
    if len(s) < 3:
        return None
    ds, dr = s - s.mean(0), r - r.mean(0)
    if np.linalg.matrix_rank(ds) < 2:
        return None
    # dr = ds @ J.T, retaining the conventional output-by-input layout.
    j_t, _, rank, _ = np.linalg.lstsq(ds, dr, rcond=None)
    if rank < 2:
        return None
    j = j_t.T
    symmetric = 0.5 * (j + j.T)
    stf = symmetric - np.eye(2) * np.trace(symmetric) / 2.0
    return j, stf


def spin2_rotate(e, phi):
    """Frozen passive-screen convention used by DEV176 G5."""
    e = np.asarray(e, float)
    c, s = np.cos(2 * phi), np.sin(2 * phi)
    return np.array([e[0] * c + e[1] * s, -e[0] * s + e[1] * c])
