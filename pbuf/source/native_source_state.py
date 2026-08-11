"""Minimal lattice-native source state used by the Dev159 audit."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence


@dataclass(frozen=True)
class NativeSourceState:
    """A one-cell source; position and amplitude are native dimensionless data."""

    position: tuple[int, int, int]
    amplitude: float = 1.0
    geometry: str = "ONE_CELL"

    def __post_init__(self) -> None:
        if len(self.position) != 3:
            raise ValueError("position must contain three lattice indices")
        if self.geometry != "ONE_CELL":
            raise ValueError("Dev159 promotes only the minimal ONE_CELL geometry")

    def wrapped(self, shape: Sequence[int]) -> "NativeSourceState":
        if len(shape) != 3 or any(int(n) < 2 for n in shape):
            raise ValueError("shape must be a three-dimensional lattice")
        return NativeSourceState(tuple(int(x) % int(n) for x, n in zip(self.position, shape)),
                                 float(self.amplitude), self.geometry)
