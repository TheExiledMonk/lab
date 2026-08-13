"""Structural transverse bases and unthresholded field moments for DEV202."""
from __future__ import annotations
import numpy as np


def transverse_basis(direction: np.ndarray) -> np.ndarray:
    n = np.asarray(direction, float); n = n / np.linalg.norm(n)
    seed = np.eye(3)[np.argmin(np.abs(n))]
    e1 = seed - n * np.dot(seed, n); e1 /= np.linalg.norm(e1)
    return np.stack((e1, np.cross(n, e1)))


def weighted_periodic_centroid(weights: np.ndarray) -> np.ndarray:
    """Circular centroid, defined from all signed-absolute field weights."""
    w = np.asarray(weights, float)
    shape = w.shape
    result = []
    for axis, size in enumerate(shape):
        marginal = w.sum(axis=tuple(i for i in range(3) if i != axis))
        angle = 2*np.pi*np.arange(size)/size
        z = np.sum(marginal*np.exp(1j*angle))
        result.append(float((np.angle(z) % (2*np.pi))*size/(2*np.pi)) if z else float("nan"))
    return np.asarray(result)
