"""Coefficient-free DEV203 relational pattern decompositions."""
from __future__ import annotations
import numpy as np


def decompose(vector: np.ndarray, direction: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    n = np.asarray(direction, float); n = n / np.linalg.norm(n)
    parallel = np.sum(vector * n, axis=-1, keepdims=True) * n
    return parallel, vector - parallel


def opposite_components(directed: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Opposite-pair sum/difference, with axis order x,y,z, retaining vectors."""
    x = np.asarray(directed, float)
    return (x[..., 0::2, :] + x[..., 1::2, :]) / 2, (x[..., 0::2, :] - x[..., 1::2, :]) / 2


def directional_tensor(directed: np.ndarray) -> np.ndarray:
    """Native opposite-bond directional matrix M_ij=(dr_+i,j-dr_-i,j)/2."""
    return (np.asarray(directed)[..., 0::2, :] - np.asarray(directed)[..., 1::2, :]) / 2


def native_symmetric_antisymmetric(directed: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    m = directional_tensor(directed)
    return (m + np.swapaxes(m, -1, -2)) / 2, (m - np.swapaxes(m, -1, -2)) / 2


def normalized_dot(a: np.ndarray, b: np.ndarray) -> float:
    a, b = np.ravel(a), np.ravel(b)
    den = np.linalg.norm(a) * np.linalg.norm(b)
    return float(np.dot(a, b) / den) if den else float("nan")
