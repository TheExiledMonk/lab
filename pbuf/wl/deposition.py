"""Target-blind detector deposition rules for observer stability audits.

Coordinates are never modified.  Every weighted rule normalizes a valid
in-domain ray at detector edges, so its (possibly signed) value is conserved.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Protocol

import numpy as np


class DepositionMethod(Protocol):
    name: str

    def deposit(
        self,
        u: np.ndarray,
        v: np.ndarray,
        values: np.ndarray | None,
        *,
        bins: int,
        extent: float,
    ) -> np.ndarray: ...


def _inputs(u, v, values, bins: int, extent: float):
    if not isinstance(bins, (int, np.integer)) or bins <= 0:
        raise ValueError("bins must be a positive integer (square detectors only)")
    if not np.isfinite(extent) or extent <= 0:
        raise ValueError("extent must be positive and finite")
    x = np.asarray(u, dtype=np.float64)
    y = np.asarray(v, dtype=np.float64)
    if x.ndim != 1 or y.ndim != 1 or x.shape != y.shape:
        raise ValueError("u and v must be equal-length one-dimensional arrays")
    w = np.ones_like(x) if values is None else np.asarray(values, dtype=np.float64)
    if w.shape != x.shape:
        raise ValueError("values must have the same shape as u and v")
    valid = np.isfinite(x) & np.isfinite(y) & np.isfinite(w)
    valid &= (x >= -extent) & (x <= extent) & (y >= -extent) & (y <= extent)
    return x, y, w, valid, 2.0 * float(extent) / bins


def _accumulate(rows, cols, weights, bins):
    out = np.zeros((bins, bins), dtype=np.float64)
    np.add.at(out, (rows, cols), weights)
    return out


@dataclass(frozen=True)
class HardBinCurrent:
    name: str = "hard_bin_current"

    def deposit(self, u, v, values=None, *, bins, extent):
        x, y, w, valid, _ = _inputs(u, v, values, bins, extent)
        edges = np.linspace(-extent, extent, bins + 1, dtype=np.float64)
        return np.histogram2d(y[valid], x[valid], bins=(edges, edges), weights=w[valid])[0]


@dataclass(frozen=True)
class HardBinHalfOpen:
    name: str = "hard_bin_half_open"

    def deposit(self, u, v, values=None, *, bins, extent):
        x, y, w, valid, width = _inputs(u, v, values, bins, extent)
        # Explicit [edge_i, edge_i+1), with the final upper endpoint included.
        col = np.floor((x[valid] + extent) / width).astype(np.int64)
        row = np.floor((y[valid] + extent) / width).astype(np.int64)
        col = np.minimum(col, bins - 1)
        row = np.minimum(row, bins - 1)
        return _accumulate(row, col, w[valid], bins)


@dataclass(frozen=True)
class NearestCenter:
    name: str = "nearest_center"

    def deposit(self, u, v, values=None, *, bins, extent):
        x, y, w, valid, width = _inputs(u, v, values, bins, extent)
        # ceil(q - 1/2) makes exact equidistant ties choose the lower index.
        qx, qy = (x[valid] + extent) / width, (y[valid] + extent) / width
        col = np.clip(np.ceil(qx - 1.0).astype(np.int64), 0, bins - 1)
        row = np.clip(np.ceil(qy - 1.0).astype(np.int64), 0, bins - 1)
        return _accumulate(row, col, w[valid], bins)


class _WeightedMethod:
    def _axis(self, coordinate, bins, extent, width):
        raise NotImplementedError

    def deposit(self, u, v, values=None, *, bins, extent):
        x, y, w, valid, width = _inputs(u, v, values, bins, extent)
        x, y, w = x[valid], y[valid], w[valid]
        out = np.zeros((bins, bins), dtype=np.float64)
        x_indices, x_weights = self._axis(x, bins, extent, width)
        y_indices, y_weights = self._axis(y, bins, extent, width)
        norm = np.zeros(x.size, dtype=np.float64)
        for ix, wx in zip(x_indices, x_weights):
            for iy, wy in zip(y_indices, y_weights):
                ok = (ix >= 0) & (ix < bins) & (iy >= 0) & (iy < bins)
                norm[ok] += wx[ok] * wy[ok]
        if np.any(norm <= 0):
            raise RuntimeError("valid detector ray has zero deposition support")
        for ix, wx in zip(x_indices, x_weights):
            for iy, wy in zip(y_indices, y_weights):
                ok = (ix >= 0) & (ix < bins) & (iy >= 0) & (iy < bins)
                np.add.at(out, (iy[ok], ix[ok]), w[ok] * wx[ok] * wy[ok] / norm[ok])
        return out


@dataclass(frozen=True)
class BilinearCIC(_WeightedMethod):
    name: str = "bilinear_cic"

    def _axis(self, coordinate, bins, extent, width):
        q = (coordinate + extent) / width - 0.5
        lo = np.floor(q).astype(np.int64)
        f = q - lo
        return (lo, lo + 1), (1.0 - f, f)


def _tsc_weight(distance):
    a = np.abs(distance)
    return np.where(a <= 0.5, 0.75 - a * a,
                    np.where(a <= 1.5, 0.5 * (1.5 - a) ** 2, 0.0))


@dataclass(frozen=True)
class TSC3x3(_WeightedMethod):
    name: str = "tsc_3x3"

    def _axis(self, coordinate, bins, extent, width):
        q = (coordinate + extent) / width - 0.5
        center = np.floor(q + 0.5).astype(np.int64)
        indices = (center - 1, center, center + 1)
        return indices, tuple(_tsc_weight(q - index) for index in indices)


@dataclass(frozen=True)
class GaussianSigmaHalfCell(_WeightedMethod):
    name: str = "gaussian_sigma_half_cell"

    def _axis(self, coordinate, bins, extent, width):
        q = (coordinate + extent) / width - 0.5
        center = np.floor(q + 0.5).astype(np.int64)
        # sigma=.5 cell and radius=3 sigma admits offsets center +/- 2,
        # filtering exact distances outside 1.5 cell widths.
        indices = tuple(center + offset for offset in (-2, -1, 0, 1, 2))
        weights = tuple(np.where(np.abs(q - index) <= 1.5,
                                 np.exp(-0.5 * ((q - index) / 0.5) ** 2), 0.0)
                        for index in indices)
        return indices, weights


METHODS: tuple[DepositionMethod, ...] = (
    HardBinCurrent(), HardBinHalfOpen(), NearestCenter(), BilinearCIC(),
    TSC3x3(), GaussianSigmaHalfCell(),
)
METHOD_BY_NAME = {method.name: method for method in METHODS}


def get_deposition_method(method: str | DepositionMethod | None) -> DepositionMethod:
    if method is None:
        return METHOD_BY_NAME["hard_bin_current"]
    if isinstance(method, str):
        try:
            return METHOD_BY_NAME[method]
        except KeyError as exc:
            raise ValueError(f"unknown deposition method: {method}") from exc
    if not hasattr(method, "name") or not callable(getattr(method, "deposit", None)):
        raise TypeError("deposition method must implement name and deposit")
    return method
