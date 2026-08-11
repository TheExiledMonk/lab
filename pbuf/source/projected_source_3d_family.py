"""Projection-equivalent diagnostic LOS realizations of a 2D source field.

Array order is ``(z, y, x)``.  The profiles are dimensionless fixtures in
native cells; none is an estimate of physical depth.
"""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np


@dataclass(frozen=True)
class SourceRealization:
    name: str
    family: str
    source: np.ndarray
    depth_cells: tuple[int, ...]
    physical_truth_claimed: bool = False


def project(source: np.ndarray) -> np.ndarray:
    source = np.asarray(source, float)
    if source.ndim != 3:
        raise ValueError("source must have (z,y,x) shape")
    return source.sum(axis=0)


def _normalized(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, float)
    if np.any(values < 0) or not values.sum() > 0:
        raise ValueError("LOS weights must be nonnegative and nonzero")
    return values / values.sum()


def _separable(image: np.ndarray, weights: np.ndarray) -> np.ndarray:
    return _normalized(weights)[:, None, None] * image[None, :, :]


def diagnostic_family(image: np.ndarray, nz: int = 17) -> list[SourceRealization]:
    """Return six distinct, exactly projection-equivalent ambiguity fixtures."""
    image = np.asarray(image, float)
    if image.ndim != 2 or np.any(image < 0) or not image.sum() > 0 or nz < 9 or nz % 2 == 0:
        raise ValueError("image must be nonnegative 2D and nz an odd integer >= 9")
    c = nz // 2
    def one(indices, values=None):
        w = np.zeros(nz); w[np.asarray(indices)] = 1 if values is None else values
        return w
    rows = [
        SourceRealization("Z01_THIN_CENTRAL", "THIN", _separable(image, one([c])), (0,)),
        SourceRealization("Z02_UNIFORM_DEPTH_2", "FINITE_SYMMETRIC", _separable(image, one(range(c-2,c+3))), tuple(range(-2,3))),
        SourceRealization("Z02_UNIFORM_DEPTH_4", "FINITE_SYMMETRIC", _separable(image, one(range(c-4,c+5))), tuple(range(-4,5))),
        SourceRealization("Z03_SYMMETRIC_COMPACT", "SYMMETRIC_COMPACT", _separable(image, one(range(c-3,c+4), [1,2,3,4,3,2,1])), tuple(range(-3,4))),
        SourceRealization("Z04_DOUBLE_LAYER", "DOUBLE_LAYER", _separable(image, one([c-4,c+4])), (-4,4)),
        SourceRealization("Z05_ASYMMETRIC", "ASYMMETRIC", _separable(image, one([c,c+2,c+4], [1,2,5])), (0,2,4)),
    ]
    yy, xx = np.indices(image.shape)
    # A fixed image-coordinate partition: no lens output or target feature enters it.
    layer = np.where(((xx // 8 + yy // 8) % 2) == 0, c-4, c+3)
    spatial = np.zeros((nz,) + image.shape)
    np.put_along_axis(spatial, layer[None, ...], image[None, ...], axis=0)
    rows.append(SourceRealization("Z06_SPATIALLY_VARIABLE", "SPATIALLY_VARIABLE_DEPTH", spatial, (-4,3)))
    return rows


def projection_error(realization: SourceRealization, image: np.ndarray) -> float:
    return float(np.max(np.abs(project(realization.source) - np.asarray(image, float))))
