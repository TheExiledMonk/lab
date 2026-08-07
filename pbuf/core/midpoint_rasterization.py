"""M10 — Midpoint Rasterization (interface field).

CORRECTION-001: The corrected rasterizer is now the dedicated
M10 implementation, not a thin re-export. It uses ``[:-1]`` for the
valid source slice (the last valid pair adjacent to the upper
boundary is no longer omitted).
"""
from __future__ import annotations
import numpy as np

from .conventions import EPS_FLOAT, EPS_ZERO, N6_POSITIVE_DIRECTIONS
from .pair_transfer import (
    rasterize_interface_field,
    rasterize_interface_field_reference,
    interface_pair_count_audit,
    expected_interface_pair_count,
    consumed_interface_pair_count,
)

__all__ = [
    "rasterize_interface_field",
    "rasterize_interface_field_reference",
    "interface_pair_count_audit",
    "expected_interface_pair_count",
    "consumed_interface_pair_count",
]
