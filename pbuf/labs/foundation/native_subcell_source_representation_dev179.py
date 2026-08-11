"""Isolated DEV179 source-representation audit helpers.

The result is intentionally negative: current code supports a discrete source
contact, not a continuous coordinate-to-loading map.  Keeping that result in
code prevents accidental interpolation from being mistaken for physics.
"""
from __future__ import annotations

from .native_subcell_geometry_dev179 import (
    SubcellSourceRepresentationNotDerived,
    native_node_coordinate,
    node_contact_loading,
)

__all__ = ["SubcellSourceRepresentationNotDerived", "native_node_coordinate", "node_contact_loading"]
