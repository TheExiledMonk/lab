"""Dev156 experimental scalar state on the native Cartesian N6 lattice."""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np


@dataclass
class NativeRelationalState:
    """Node state and optional retained change; neither is an EM quantity."""

    local_state: np.ndarray
    retained_change: np.ndarray | None = None
    progression_step: int = 0

    def __post_init__(self) -> None:
        self.local_state = np.asarray(self.local_state, dtype=np.float64)
        if self.local_state.ndim != 3 or min(self.local_state.shape) < 3:
            raise ValueError("local_state must be a 3D N6 lattice")
        if not np.isfinite(self.local_state).all():
            raise ValueError("local_state must be finite")
        if self.retained_change is not None:
            self.retained_change = np.asarray(self.retained_change, dtype=np.float64)
            if self.retained_change.shape != self.local_state.shape:
                raise ValueError("retained_change must match local_state")
            if not np.isfinite(self.retained_change).all():
                raise ValueError("retained_change must be finite")


def centered_coordinates(shape: tuple[int, int, int]) -> tuple[np.ndarray, ...]:
    return tuple(np.arange(n, dtype=float) - n // 2 for n in shape)


def perturbation(kind: str, shape=(17, 17, 17), amplitude=1.0) -> np.ndarray:
    """Target-blind P01--P05 initial node states."""
    shape = tuple(map(int, shape))
    q = np.zeros(shape, dtype=np.float64)
    c = tuple(n // 2 for n in shape)
    if kind == "P01":
        q[c] = amplitude
    elif kind == "P02":
        q[c] = amplitude
        q[(c[0] - 1, c[1], c[2])] = -amplitude
    elif kind == "P03":
        q[(c[0] + 1, c[1], c[2])] = amplitude
        q[(c[0] - 1, c[1], c[2])] = -amplitude
    elif kind == "P04":
        axes = centered_coordinates(shape)
        x, y, z = np.meshgrid(*axes, indexing="ij")
        # Width is fixed by the smallest nontrivial compact lattice stencil.
        q = amplitude * np.exp(-(x*x + y*y + z*z) / 8.0)
    elif kind == "P05":
        q[(c[0] + 1, c[1], c[2])] = amplitude
        q[(c[0] - 1, c[1], c[2])] = -amplitude
    else:
        raise ValueError(f"unknown perturbation {kind}")
    return q
