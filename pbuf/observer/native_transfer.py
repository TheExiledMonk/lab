"""DEV188 finite native incident-distribution transfer operators.

These objects operate only on the frozen discrete DEV183 launch IDs.  They do
not define an astronomical image, a continuous source plane, or an intensity
calibration.
"""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np


@dataclass(frozen=True)
class NativeIncidentDistribution:
    """Nonnegative content indexed explicitly by frozen native launch IDs."""
    launch_ids: tuple[str, ...]
    values: np.ndarray
    normalization: str = "unnormalized_native_content"
    representation: str = "NativeIncidentDistribution/v1"

    def __post_init__(self) -> None:
        values = np.asarray(self.values, dtype=np.float64)
        if values.ndim != 1 or len(values) != len(self.launch_ids):
            raise ValueError("values must be one coefficient for every explicit launch ID")
        if np.any(values < 0):
            raise ValueError("native incident coefficients must be nonnegative")
        object.__setattr__(self, "values", values)

    @property
    def total_content(self) -> float:
        return float(self.values.sum())


@dataclass(frozen=True)
class NativeTransferOperator:
    """One realization's exact DEV184-to-DEV187 receipt-cell response."""
    realization_id: str
    launch_ids: tuple[str, ...]
    detector_cell_ids: np.ndarray
    detector_coordinates: np.ndarray
    weight_kernel: np.ndarray
    additive_channel_kernels: dict[str, np.ndarray]
    source_domain_representation: str = "NativeIncidentDistribution/v1"
    detector_representation: str = "NativeDetectorState/v1"
    representation: str = "NativeTransferOperator/v1"

    def __post_init__(self) -> None:
        weight = np.asarray(self.weight_kernel, dtype=np.float64)
        if weight.shape != (len(self.detector_cell_ids), len(self.launch_ids)):
            raise ValueError("weight kernel shape is not aligned with explicit IDs")
        if np.any(weight < 0):
            raise ValueError("positive-crossing weight kernel must be nonnegative")
        object.__setattr__(self, "weight_kernel", weight)

    def pushforward(self, source: NativeIncidentDistribution | np.ndarray) -> np.ndarray:
        """Return the scalar native positive-receipt measure K s."""
        if isinstance(source, NativeIncidentDistribution):
            if source.launch_ids != self.launch_ids:
                raise ValueError("source launch IDs do not exactly match operator launch IDs")
            values = source.values
        else:
            values = np.asarray(source, dtype=np.float64)
            if values.shape != (len(self.launch_ids),):
                raise ValueError("source array must have one value per explicit launch ID")
            if np.any(values < 0):
                raise ValueError("native incident coefficients must be nonnegative")
        return self.weight_kernel @ values

    def conditional_weight_kernel(self) -> tuple[np.ndarray, np.ndarray]:
        """Return P(c|launch) and a zero-throughput mask (undefined columns are NaN)."""
        throughput = self.weight_kernel.sum(axis=0)
        conditional = np.full_like(self.weight_kernel, np.nan)
        positive = throughput > 0
        conditional[:, positive] = self.weight_kernel[:, positive] / throughput[positive]
        return conditional, ~positive


def weighted_detector_geometry(coordinates: np.ndarray, measure: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Centroid and central second moment, undefined for zero total measure."""
    q = np.asarray(coordinates, dtype=np.float64)
    d = np.asarray(measure, dtype=np.float64)
    total = float(d.sum())
    if total <= 0:
        return np.full(2, np.nan), np.full((2, 2), np.nan)
    centroid = (d[:, None] * q).sum(axis=0) / total
    delta = q - centroid
    return centroid, (delta.T * d) @ delta / total
