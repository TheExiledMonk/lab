"""M01 — Conventions Registry.

Owns ALL shared conventions for the verified numerical core. No other
module may define its own axis, component, stencil, or tolerance
convention.

Correction pass FOUNDATION-001-CORRECTION-001
---------------------------------------------
* Removed EPS_HASH (hashing is exact, no numerical epsilon needed).
* Split single EPS_* constants into purpose-specific tolerances:
  - EPS_EXACT_COMPARISON   — for exact-zero / exact-identity checks.
  - EPS_VARIANCE_UNDEFINED — for "variance undefined" guard.
  - EPS_NORM_RELATIVE      — for relative-norm tolerance.
  - EPS_TRANSFORM          — for transform round-trips.
  - EPS_CLOSURE            — for global closure checks.
* Added backward-compatible aliases (EPS_FLOAT, EPS_ZERO) for legacy
  imports. The new purpose-specific tolerances are authoritative.
"""
from __future__ import annotations
import numpy as np
from types import MappingProxyType

__all__ = [
    "ARRAY_AXIS_ORDER", "VECTOR_COMPONENT_ORDER", "TENSOR_COMPONENT_ORDER",
    "AXIS_TO_ARRAY_INDEX", "COMPONENT_TO_INDEX",
    "EPS_FLOAT", "EPS_ZERO",
    "EPS_EXACT_COMPARISON", "EPS_VARIANCE_UNDEFINED",
    "EPS_NORM_RELATIVE", "EPS_TRANSFORM", "EPS_CLOSURE",
    "CONVENTIONS_VERSION",
    "RC_TRANSFORMS", "RC_MATRICES_FWD", "RC_MATRICES_INV",
    "N6_DIRECTIONS", "N6_POSITIVE_DIRECTIONS", "STENCIL_OFFSETS",
    "AXIS_OFFSETS", "ARRAY_ORDER",
]

CONVENTIONS_VERSION = "1.1.0-correction001"

# Array axes are (z, y, x) — axis 0 = z, axis 1 = y, axis 2 = x.
ARRAY_AXIS_ORDER = ("z", "y", "x")

# Vector components are stored in (x, y, z) order — index 0 = x, 1 = y, 2 = z.
VECTOR_COMPONENT_ORDER = ("x", "y", "z")

# Symmetric tensor components, in the canonical 6-vector order:
#   (xx, xy, xz, yy, yz, zz)
TENSOR_COMPONENT_ORDER = ("xx", "xy", "xz", "yy", "yz", "zz")

# Mapping from axis label (string) to integer array axis.
AXIS_TO_ARRAY_INDEX = MappingProxyType({
    "z": 0,
    "y": 1,
    "x": 2,
})

# Mapping from vector component label to its integer index in (x, y, z).
COMPONENT_TO_INDEX = MappingProxyType({
    "x": 0,
    "y": 1,
    "z": 2,
})

# ----------------------------------------------------------------------
# Purpose-specific numerical tolerances (CORRECTION-001).
# ----------------------------------------------------------------------
EPS_EXACT_COMPARISON = 0.0     # for exact-zero / exact-identity checks
EPS_VARIANCE_UNDEFINED = 1e-15  # variance floor for "undefined correlation"
EPS_NORM_RELATIVE = 1e-14      # relative tolerance for norm comparisons
EPS_TRANSFORM = 1e-14          # tolerance for transform round-trips
EPS_CLOSURE = 1e-12            # global closure tolerance

# Backward-compatible aliases.
EPS_FLOAT = EPS_NORM_RELATIVE
EPS_ZERO = EPS_EXACT_COMPARISON


# ----------------------------------------------------------------------
# Coordinate transforms — RC0..RC6.
# ----------------------------------------------------------------------
RC_TRANSFORMS = ("RC0", "RC1", "RC2", "RC3", "RC4", "RC5", "RC6")


def _build_RC_matrices() -> dict:
    """Construct the seven RC matrices exactly as the corrected lab does."""
    Q = {
        "RC0": np.eye(3, dtype=np.float64),
        # Coordinate swaps (permutations with det = -1):
        "RC1": np.array([[0, 1, 0], [1, 0, 0], [0, 0, 1]], dtype=np.float64),
        "RC2": np.array([[0, 0, 1], [0, 1, 0], [1, 0, 0]], dtype=np.float64),
        "RC3": np.array([[1, 0, 0], [0, 0, 1], [0, 1, 0]], dtype=np.float64),
        # 90-degree rotations (det = +1):
        # RC4: +90° about x  →  (x, y, z) → (x, -z, y)
        "RC4": np.array([[1, 0, 0], [0, 0, -1], [0, 1, 0]], dtype=np.float64),
        # RC5: +90° about y  →  (x, y, z) → (z, y, -x)
        "RC5": np.array([[0, 0, 1], [0, 1, 0], [-1, 0, 0]], dtype=np.float64),
        # RC6: +90° about z  →  (x, y, z) → (-y, x, z)
        "RC6": np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]], dtype=np.float64),
    }
    return Q


RC_MATRICES_FWD = MappingProxyType(_build_RC_matrices())
RC_MATRICES_INV = MappingProxyType({
    k: v.T for k, v in _build_RC_matrices().items()
})

# ----------------------------------------------------------------------
# Static closed-form expected mapping table (CORRECTION-001, §16.1).
# This table is independent of any production helper or registry
# and is used by M02/M03/M04 independent reference validation.
# ----------------------------------------------------------------------
# Each entry: new (x, y, z) coordinate expressed in terms of OLD (x, y, z).
# The "permutation" gives the array-axis remapping in (z, y, x) array order
# AFTER the permutation step. "flip_array_axis" gives the array axis
# (in (z, y, x) order, AFTER the permutation step) that np.flip is
# applied to. This matches the production implementation convention.
EXPECTED_AXIS_MAPPING = MappingProxyType({
    "RC0": {
        "x_new": "x",
        "y_new": "y",
        "z_new": "z",
        "det": +1, "is_rotation": False,
        "permutation": (0, 1, 2),
        "flip_array_axis": (),
    },
    "RC1": {
        "x_new": "y",
        "y_new": "x",
        "z_new": "z",
        "det": -1, "is_rotation": False,
        "permutation": (0, 2, 1),
        "flip_array_axis": (),
    },
    "RC2": {
        "x_new": "z",
        "y_new": "y",
        "z_new": "x",
        "det": -1, "is_rotation": False,
        "permutation": (2, 1, 0),
        "flip_array_axis": (),
    },
    "RC3": {
        "x_new": "x",
        "y_new": "z",
        "z_new": "y",
        "det": -1, "is_rotation": False,
        "permutation": (1, 0, 2),
        "flip_array_axis": (),
    },
    # 90° rotations. flip_array_axis is the array axis (after perm)
    # where np.flip is applied.
    "RC4": {
        "x_new": "x",
        "y_new": "-z",
        "z_new": "y",
        "det": +1, "is_rotation": True,
        "permutation": (1, 0, 2),
        "flip_array_axis": (1,),
    },
    "RC5": {
        "x_new": "z",
        "y_new": "y",
        "z_new": "-x",
        "det": +1, "is_rotation": True,
        "permutation": (2, 1, 0),
        "flip_array_axis": (0,),
    },
    "RC6": {
        "x_new": "-y",
        "y_new": "x",
        "z_new": "z",
        "det": +1, "is_rotation": True,
        "permutation": (0, 2, 1),
        "flip_array_axis": (2,),
    },
})

# ----------------------------------------------------------------------
# Six N6 neighbour directions expressed as integer triples (dx, dy, dz)
# in (x, y, z) component order.
# ----------------------------------------------------------------------
N6_DIRECTIONS = MappingProxyType({
    "xp": np.array([+1, 0, 0], dtype=np.int64),
    "xm": np.array([-1, 0, 0], dtype=np.int64),
    "yp": np.array([0, +1, 0], dtype=np.int64),
    "ym": np.array([0, -1, 0], dtype=np.int64),
    "zp": np.array([0, 0, +1], dtype=np.int64),
    "zm": np.array([0, 0, -1], dtype=np.int64),
})

N6_POSITIVE_DIRECTIONS = ("xp", "yp", "zp")

# Mapping from N6 direction label to the array-axis along which the
# neighbour lies (axis 0 = z, axis 1 = y, axis 2 = x).
STENCIL_OFFSETS = MappingProxyType({
    "xp": ("x", +1), "xm": ("x", -1),
    "yp": ("y", +1), "ym": ("y", -1),
    "zp": ("z", +1), "zm": ("z", -1),
})

# Offsets expressed in ARRAY axes (z, y, x) for compatibility with the
# storage layout. Built from STENCIL_OFFSETS by remapping the axis index.
AXIS_OFFSETS = MappingProxyType({
    label: (AXIS_TO_ARRAY_INDEX[axis_name], sign)
    for label, (axis_name, sign) in STENCIL_OFFSETS.items()
})

# Convenience alias to express array-axis order.
ARRAY_ORDER = ARRAY_AXIS_ORDER


def get_coordinate_matrix(transform_id: str, inverse: bool = False) -> np.ndarray:
    """Return the 3x3 orthogonal matrix for ``transform_id``."""
    if transform_id not in RC_MATRICES_FWD:
        raise ValueError(f"unknown transform_id: {transform_id!r}")
    if inverse:
        return RC_MATRICES_INV[transform_id].copy()
    return RC_MATRICES_FWD[transform_id].copy()


def validate_transform_id(transform_id: str) -> None:
    """Raise ValueError if ``transform_id`` is not a registered RC name."""
    if transform_id not in RC_TRANSFORMS:
        raise ValueError(f"unknown transform_id: {transform_id!r}")


# ----------------------------------------------------------------------
# Pair symmetrisation (PS) lanes (CORRECTION-001, §8).
# PS1-A — raw single-endpoint directional diagnostic.
#           v_ij = P_i n̂_ij  (NOT antisymmetrised).
# PS1   — antisymmetrised source-local: 0.5 (v_i - v_j).
# PS1-B — midpoint antisymmetrised: (v_i - v_j)/2 with project of v_j
#           through partner's projector at j.
# PS2   — midpoint-symmetrised projector: 0.5 (P_i + P_j) n̂.
#
# All four lanes are DECLARED and DISTINCT.  PS1-B is mathematically
# identical to PS2 for unit-magnitude (PM1) outputs (because
# normalising a linear combination of two unit vectors is not equal
# to normalising each of them separately) — the underlying unscaled
# R_ij BEFORE magnitude normalisation differs.
# ----------------------------------------------------------------------
PS_LANES = ("PS1-A", "PS1", "PS1-B", "PS2")

# Magnitude formulations.
PM_LANES = ("PM1", "PM2")


if __name__ == "__main__":
    # Self-check: all 7 RC matrices are orthogonal and det = ±1.
    for rc in RC_TRANSFORMS:
        Q = RC_MATRICES_FWD[rc]
        assert np.allclose(Q @ Q.T, np.eye(3)), f"{rc} not orthogonal"
        det = np.linalg.det(Q)
        assert abs(abs(det) - 1.0) < 1e-12, f"{rc} det != ±1"
        # Cross-check against static EXPECTED_AXIS_MAPPING.
        m = EXPECTED_AXIS_MAPPING[rc]
        assert m["det"] == (1 if det > 0 else -1), f"{rc} det mismatch"
    # N6 directions are unit and antiparallel partners.
    assert np.allclose(N6_DIRECTIONS["xp"], -N6_DIRECTIONS["xm"])
    assert np.allclose(N6_DIRECTIONS["yp"], -N6_DIRECTIONS["ym"])
    assert np.allclose(N6_DIRECTIONS["zp"], -N6_DIRECTIONS["zm"])
    print("M01 conventions: self-check passed")
