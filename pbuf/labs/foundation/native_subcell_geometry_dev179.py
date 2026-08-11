"""DEV179 audit boundary for native source coordinates.

This module deliberately does *not* construct a sub-cell loading kernel.  The
frozen DEV167 operation has a one-cell integer contact centre; accepting an
off-node coordinate would require a new source-to-medium coupling law.
"""
from __future__ import annotations

import numpy as np

from pbuf.excitation.native_vector_pair_dynamics import source_contact_force


class SubcellSourceRepresentationNotDerived(ValueError):
    """Raised when a coordinate lacks a frozen native loading meaning."""


def native_node_coordinate(position: np.ndarray | tuple[float, float, float]) -> tuple[int, int, int]:
    """Validate the only source coordinate currently defined by DEV167/168."""
    p = np.asarray(position, dtype=np.float64)
    if p.shape != (3,) or not np.isfinite(p).all():
        raise ValueError("source position must be one finite 3-vector")
    integer = np.rint(p)
    if not np.array_equal(p, integer):
        raise SubcellSourceRepresentationNotDerived(
            "DEV167 defines source contact at an integer native node only; "
            "no source-to-node coupling or cell-ownership law is frozen"
        )
    return tuple(int(x) for x in integer)


def node_contact_loading(shape: tuple[int, int, int], position, magnitude: float) -> np.ndarray:
    """The exact existing source-loading operation, with no new semantics."""
    return source_contact_force(shape, native_node_coordinate(position), magnitude)
