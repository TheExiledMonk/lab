"""Coefficient-free, single-object DEV203 longitudinal geometry observers."""
from __future__ import annotations

import numpy as np

from pbuf.observer.native_n6_field import n6_field


def signed_axis_coordinates(n: int, center: int) -> np.ndarray:
    """Exact periodic lattice coordinates relative to an already-fixed center."""
    return (np.arange(n) - center + n // 2) % n - n // 2


def strain_magnitude_density(displacement: np.ndarray) -> np.ndarray:
    """q_epsilon(x)=sum_b |epsilon_xb|, the canonical scalar N6 geometry."""
    return np.abs(n6_field(displacement, np.zeros_like(displacement))["strain"]).sum(axis=-1)


def geometry_moments(q: np.ndarray, center: tuple[int, int, int]) -> dict[str, np.ndarray | float]:
    """Raw longitudinal and second spatial moments without fitted centering."""
    coords = [signed_axis_coordinates(q.shape[i], center[i]) for i in range(3)]
    grid = np.stack(np.meshgrid(*coords, indexing="ij"), axis=-1)
    s = grid[..., 0]
    plus, minus = s > 0, s < 0
    return {"s": s, "q_plus": float(q[plus].sum()), "q_minus": float(q[minus].sum()),
            "end_asymmetry": float(q[plus].sum() - q[minus].sum()),
            "first_moment": float((s * q).sum()),
            "shape_tensor": np.einsum("abc,abci,abcj->ij", q, grid, grid)}


def profile(q: np.ndarray, center_x: int) -> tuple[np.ndarray, np.ndarray]:
    s = signed_axis_coordinates(q.shape[0], center_x)
    return s, q.sum(axis=(1, 2))


def transverse_moments(q: np.ndarray, center: tuple[int, int, int]) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    y = signed_axis_coordinates(q.shape[1], center[1])[None, :, None]
    z = signed_axis_coordinates(q.shape[2], center[2])[None, None, :]
    return ((q * y * y).sum(axis=(1, 2)), (q * z * z).sum(axis=(1, 2)),
            (q * y * z).sum(axis=(1, 2)))
