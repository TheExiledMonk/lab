"""Coefficient-free DEV190 scalar conditional-transfer geometry.

This module intentionally operates on the positive scalar receipt measure only.
It does not create a norm on the heterogeneous DEV188 multichannel stack and it
does not define an astronomical or GR lensing observable.
"""
from __future__ import annotations

import numpy as np


def conditional_centroids_and_covariances(weight: np.ndarray, coordinates: np.ndarray):
    """Per-launch conditional centroids/covariances; undefined columns are NaN."""
    weight, q = np.asarray(weight, float), np.asarray(coordinates, float)
    eta = weight.sum(axis=0)
    n = weight.shape[1]
    centroids = np.full((n, 2), np.nan)
    covariance = np.full((n, 2, 2), np.nan)
    valid = eta > 0
    p = weight[:, valid] / eta[valid]
    centroids[valid] = p.T @ q
    for column in np.flatnonzero(valid):
        d = q - centroids[column]
        covariance[column] = (d.T * (weight[:, column] / eta[column])) @ d
    return eta, centroids, covariance, ~valid


def periodic_centroid_jacobian(centroids: np.ndarray, side: int = 11) -> np.ndarray:
    """The unique symmetric nearest-neighbour central difference on Z_side²."""
    y = np.asarray(centroids, float).reshape(side, side, 2)
    # Input-coordinate order is (native y, native z), matching launch storage.
    dy = (np.roll(y, -1, axis=0) - np.roll(y, 1, axis=0)) / 2.0
    dz = (np.roll(y, -1, axis=1) - np.roll(y, 1, axis=1)) / 2.0
    return np.stack((dy, dz), axis=-1).reshape(side * side, 2, 2)


def polar_fields(jacobian: np.ndarray):
    """Right polar fields, deliberately local 2x2 coordinate geometry only."""
    j = np.asarray(jacobian, float)
    n = len(j)
    u = np.full((n, 2, 2), np.nan); rotation = np.full((n, 2, 2), np.nan)
    eigen = np.full((n, 2), np.nan); anisotropy = np.full(n, np.nan)
    angle = np.full(n, np.nan); det = np.full(n, np.nan)
    for i, a in enumerate(j):
        if not np.isfinite(a).all():
            continue
        det[i] = np.linalg.det(a)
        left, s, right_t = np.linalg.svd(a)
        rotation[i] = left @ right_t
        u[i] = right_t.T @ np.diag(s) @ right_t
        eigen[i] = np.sort(s)[::-1]
        anisotropy[i] = (eigen[i, 0] - eigen[i, 1]) / (eigen[i].sum()) if eigen[i].sum() else np.nan
        angle[i] = np.arctan2(rotation[i, 1, 0], rotation[i, 0, 0])
    return u, rotation, eigen, anisotropy, det, angle


def affine_centroids(a: np.ndarray, side: int = 11) -> np.ndarray:
    """Synthetic affine fixture on exact periodic labels; central differences are exact."""
    grid = np.indices((side, side)).transpose(1, 2, 0).astype(float)
    return (grid @ np.asarray(a, float).T).reshape(side * side, 2)
