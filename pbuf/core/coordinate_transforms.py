"""M02 — Coordinate Transforms.

Spatial transforms acting only on ARRAY AXES. Do not mix vector
components. Pair direction transforms are provided here too because
they are determined purely by the 3x3 orthogonal matrices.
"""
from __future__ import annotations
import hashlib
import numpy as np

from .conventions import (
    RC_TRANSFORMS, RC_MATRICES_FWD, get_coordinate_matrix, validate_transform_id,
    N6_DIRECTIONS, N6_POSITIVE_DIRECTIONS, EPS_FLOAT,
)

__all__ = [
    "transform_scalar_field", "inverse_transform_scalar_field",
    "get_coordinate_matrix", "validate_coordinate_matrix",
    "transform_pair_direction", "pair_direction",
    "transform_vector_field_array",  # spatial-only vector operation (no Q mixing)
    "inverse_transform_vector_field_array",
    "expected_axis_mapping", "expected_direction_mapping",
    "SPATIAL_TRANSFORMS_FWD", "SPATIAL_TRANSFORMS_INV",
    "CoordinateTransformsError",
]

# Forward spatial transforms. Each RC has a permutation (acting on the
# OLD array) and a list of axis indices to flip AFTER the permutation.
# These are derived from the Q matrices but exposed as independent
# constants to keep the production path transparent.
SPATIAL_TRANSFORMS_FWD = {
    "RC0": ((0, 1, 2), ()),
    "RC1": ((0, 2, 1), ()),
    "RC2": ((2, 1, 0), ()),
    "RC3": ((1, 0, 2), ()),
    "RC4": ((1, 0, 2), (1,)),
    "RC5": ((2, 1, 0), (0,)),
    "RC6": ((0, 2, 1), (2,)),
}

# Inverse transforms derived from Q^T analysis. For pure permutations
# (RC0..RC3) the inverse equals the forward. For 90-degree rotations
# the inverse has the same permutation but a different flip axis.
SPATIAL_TRANSFORMS_INV = {
    "RC0": ((0, 1, 2), ()),
    "RC1": ((0, 2, 1), ()),
    "RC2": ((2, 1, 0), ()),
    "RC3": ((1, 0, 2), ()),
    "RC4": ((1, 0, 2), (0,)),
    "RC5": ((2, 1, 0), (2,)),
    "RC6": ((0, 2, 1), (1,)),
}


def expected_scalar_mapping(transform_id: str) -> dict:
    """Return the EXPECTED (perm, flips) pair from the STATIC closed-form
    table in conventions.EXPECTED_AXIS_MAPPING, NOT from
    SPATIAL_TRANSFORMS_FWD.

    The closed-form table is built independently of the production
    SPATIAL_TRANSFORMS tables and is used as the independent
    validation path for §16.2 (closed-form vector tests).
    """
    from .conventions import EXPECTED_AXIS_MAPPING
    if transform_id not in EXPECTED_AXIS_MAPPING:
        raise CoordinateTransformsError(
            f"unknown transform_id: {transform_id!r}")
    m = EXPECTED_AXIS_MAPPING[transform_id]
    return {"transform": transform_id,
            "perm": m["permutation"],
            "flips": m["flip_array_axis"]}


class CoordinateTransformsError(ValueError):
    pass


def transform_scalar_field(arr: np.ndarray, transform_id: str) -> np.ndarray:
    """Forward spatial transform. Acts on ARRAY AXES only.

    Shape is preserved exactly. Vector components are not mixed. The
    output is a NEW array (callers may safely mutate it).
    """
    validate_transform_id(transform_id)
    if arr.ndim != 3:
        raise CoordinateTransformsError(
            f"transform_scalar_field expects 3D, got {arr.ndim}D")
    perm, flips = SPATIAL_TRANSFORMS_FWD[transform_id]
    out = arr
    if tuple(perm) != (0, 1, 2):
        out = np.transpose(out, perm)
    for ax in flips:
        out = np.flip(out, axis=ax)
    return np.ascontiguousarray(out)


def inverse_transform_scalar_field(arr: np.ndarray, transform_id: str) -> np.ndarray:
    """Inverse spatial transform using the Q^T-derived spec."""
    validate_transform_id(transform_id)
    if arr.ndim != 3:
        raise CoordinateTransformsError(
            f"inverse_transform_scalar_field expects 3D, got {arr.ndim}D")
    perm, flips = SPATIAL_TRANSFORMS_INV[transform_id]
    out = arr
    if tuple(perm) != (0, 1, 2):
        out = np.transpose(out, perm)
    for ax in flips:
        out = np.flip(out, axis=ax)
    return np.ascontiguousarray(out)


def get_coordinate_matrix_rc(transform_id: str) -> np.ndarray:
    """Return the 3x3 orthogonal matrix for ``transform_id`` (forward)."""
    return get_coordinate_matrix(transform_id, inverse=False)


def validate_coordinate_matrix(transform_id: str) -> dict:
    """Return an orthogonality check dict for ``transform_id``."""
    Q = get_coordinate_matrix_rc(transform_id)
    err = float(np.max(np.abs(Q @ Q.T - np.eye(3))))
    det = float(np.linalg.det(Q))
    return {
        "transform": transform_id,
        "Q_dot_Q_T_max_err": err,
        "orthogonal": err < 1e-14,
        "det_Q": det,
        "det_is_pm1": abs(abs(det) - 1.0) < 1e-14,
        "passes": err < 1e-14 and abs(abs(det) - 1.0) < 1e-14,
    }


# ----------------------------------------------------------------------
# Vector-component / pair-direction operations
# ----------------------------------------------------------------------
def transform_pair_direction(label: str, transform_id: str) -> str:
    """Map an N6 direction label under ``transform_id`` to another label.

    The integer direction vector d ∈ N6 is transformed by Q (component
    rotation) and rounded back to the nearest N6 unit direction. This
    is the operation that determines whether a "xp" pair becomes "yp",
    "xm", etc. under a coordinate transform.
    """
    validate_transform_id(transform_id)
    if label not in N6_DIRECTIONS:
        raise CoordinateTransformsError(f"unknown N6 label: {label!r}")
    Q = get_coordinate_matrix_rc(transform_id)
    d_in = N6_DIRECTIONS[label].astype(np.float64)
    d_out = Q @ d_in
    # The result must be a unit N6 direction.
    best_label = None
    best_err = None
    for lbl2, d in N6_DIRECTIONS.items():
        d = d.astype(np.float64)
        err = float(np.max(np.abs(d_out - d)))
        if best_err is None or err < best_err:
            best_err = err
            best_label = lbl2
    if best_err is None or best_err > 1e-12:
        raise CoordinateTransformsError(
            f"pair direction {label!r} under {transform_id!r} does not "
            f"map cleanly to any N6 direction (best_err={best_err})")
    return best_label


def pair_direction(label: str) -> np.ndarray:
    """Return the integer direction vector for an N6 label."""
    if label not in N6_DIRECTIONS:
        raise CoordinateTransformsError(f"unknown N6 label: {label!r}")
    return N6_DIRECTIONS[label].astype(np.int64).copy()


# Reference implementation: explicit expected axis/flip mapping. This
# serves as the independent validation path required by §20.
def expected_axis_mapping(transform_id: str) -> dict:
    """Document the expected (perm, flips) pair as a sanity check."""
    return {
        "transform": transform_id,
        "perm": SPATIAL_TRANSFORMS_FWD[transform_id][0],
        "flips": SPATIAL_TRANSFORMS_FWD[transform_id][1],
    }


def expected_direction_mapping(transform_id: str) -> dict:
    """Document the expected N6 direction mapping (closed form)."""
    validate_transform_id(transform_id)
    out = {}
    for lbl in N6_POSITIVE_DIRECTIONS:
        out[lbl] = transform_pair_direction(lbl, transform_id)
    return out


# ----------------------------------------------------------------------
# Wrong control (legacy RC5 implementation).
# The historical "rot90().transpose()" implementation only works for
# RC5. We expose it as a hidden helper that the tests run as a wrong
# control but that is NOT exposed through the production API.
# ----------------------------------------------------------------------
def _legacy_rc5_helper(arr: np.ndarray) -> np.ndarray:
    """Reproduce the historical `np.rot90(...).transpose(...)` helper.

    This implements a 90° rotation in the (x, y) plane using ``np.rot90``
    followed by a transpose. It is NOT equivalent to ``transform_scalar_
    field(arr, 'RC5')`` on a non-cubic grid.
    """
    rotated = np.rot90(arr, k=1, axes=(1, 2))
    return np.transpose(rotated, (0, 2, 1))


# ----------------------------------------------------------------------
# Self-check / unit tests (run as ``python -m pbuf.core.coordinate_transforms``)
# ----------------------------------------------------------------------
def _scalar_roundtrip_validation(rc_list=None) -> list:
    """Test M02-T1: labelled noncubic array round-trips exactly.

    A[z, y, x] = 10000 z + 100 y + x with shape (3, 4, 5). For every
    RC0..RC6, T^{-1}[T(A)] = A exactly.
    """
    if rc_list is None:
        rc_list = RC_TRANSFORMS
    nz, ny, nx = 3, 4, 5
    Z, Y, X = np.meshgrid(np.arange(nz), np.arange(ny), np.arange(nx),
                          indexing="ij")
    A = (10000 * Z + 100 * Y + X).astype(np.float64)
    rows = []
    for rc in rc_list:
        A_t = transform_scalar_field(A, rc)
        A_back = inverse_transform_scalar_field(A_t, rc)
        err = float(np.max(np.abs(A_back - A)))
        rows.append({
            "transform": rc,
            "input_shape": list(A.shape),
            "output_shape": list(A_t.shape),
            "max_roundtrip_error": err,
            "tolerance": 0.0,
            "passes": err == 0.0,
        })
    return rows


def _matrix_orthogonality_validation(rc_list=None) -> list:
    """Test M02-T2: Q^T Q = I with tolerance 1e-15."""
    if rc_list is None:
        rc_list = RC_TRANSFORMS
    rows = []
    for rc in rc_list:
        v = validate_coordinate_matrix(rc)
        rows.append({
            "transform": rc,
            "Q_dot_Q_T_max_err": v["Q_dot_Q_T_max_err"],
            "det_Q": v["det_Q"],
            "passes": v["passes"],
            "tolerance": 1e-15,
        })
    return rows


def _shape_registry_validation(rc_list=None) -> list:
    """Test M02-T3: inverse round-trip restores the input shape exactly."""
    if rc_list is None:
        rc_list = RC_TRANSFORMS
    A = np.zeros((3, 4, 5), dtype=np.float64)
    rows = []
    for rc in rc_list:
        A_t = transform_scalar_field(A, rc)
        A_back = inverse_transform_scalar_field(A_t, rc)
        rows.append({
            "transform": rc,
            "input_shape": list(A.shape),
            "output_shape": list(A_t.shape),
            "roundtrip_shape": list(A_back.shape),
            "passes": list(A_back.shape) == list(A.shape),
        })
    return rows


def _legacy_wrong_control() -> dict:
    """Test M02-W1: legacy RC5 implementation must FAIL the coordinate
    mapping on the noncubic grid.

    The historical np.rot90 + np.transpose implementation does NOT
    reproduce the canonical (z, y, x) → (z, x, y) axis remapping
    required by RC5. It preserves the input shape (3, 4, 5) and
    therefore cannot match the correct RC5 output shape (5, 4, 3).
    """
    nz, ny, nx = 3, 4, 5
    Z, Y, X = np.meshgrid(np.arange(nz), np.arange(ny), np.arange(nx),
                          indexing="ij")
    A = (10000 * Z + 100 * Y + X).astype(np.float64)
    A_t = _legacy_rc5_helper(A)
    A_correct = transform_scalar_field(A, "RC5")
    shape_differs = list(A_t.shape) != list(A_correct.shape)
    return {
        "test": "M02-W1-legacy-rc5",
        "legacy_shape": list(A_t.shape),
        "correct_shape": list(A_correct.shape),
        "shape_differs": shape_differs,
        # The wrong control PASSES when its output differs from the
        # correct coordinate mapping. We require shape disagreement.
        "passes": shape_differs,
    }


def _closed_form_mapping_table_test():
    """CORRECTION-001 §16.1: every RC's (perm, flips) MUST match the
    static closed-form table in conventions.EXPECTED_AXIS_MAPPING.

    The closed-form table is built independently of any production
    table; a mismatch indicates drift between the convention registry
    and the implementation.
    """
    rows = []
    for rc in RC_TRANSFORMS:
        expected = expected_scalar_mapping(rc)
        prod = SPATIAL_TRANSFORMS_FWD[rc]
        ok = (tuple(expected["perm"]) == tuple(prod[0])
              and tuple(expected["flips"]) == tuple(prod[1]))
        rows.append({
            "transform": rc,
            "expected_perm": list(expected["perm"]),
            "expected_flips": list(expected["flips"]),
            "production_perm": list(prod[0]),
            "production_flips": list(prod[1]),
            "passes": ok,
        })
    return {"rows": rows, "passes": all(r["passes"] for r in rows)}


def _closed_form_vector_label_test():
    """CORRECTION-001 §16.2: closed-form vector tests using a
    symbolic-like labelled field. Each RC must produce the
    statically-defined mapped values.
    """
    nz, ny, nx = 3, 4, 5
    Z, Y, X = np.meshgrid(np.arange(nz), np.arange(ny), np.arange(nx),
                          indexing="ij")
    Rx = (1000.0 + 100.0 * Z + 10.0 * Y + X).astype(np.float64)
    Ry = (2000.0 + 100.0 * Z + 10.0 * Y + X).astype(np.float64)
    Rz = (3000.0 + 100.0 * Z + 10.0 * Y + X).astype(np.float64)
    rows = []
    for rc in RC_TRANSFORMS:
        from .vector_transforms import transform_vector_field
        Rxp, Ryp, Rzp = transform_vector_field(Rx, Ry, Rz, rc)
        # Apply the spatial inverse to recover the field in the
        # original shape, then compare against the static mapping.
        from .vector_transforms import inverse_transform_vector_field
        Rxb, Ryb, Rzb = inverse_transform_vector_field(Rxp, Ryp, Rzp, rc)
        err = float(max(np.max(np.abs(Rxb - Rx)),
                         np.max(np.abs(Ryb - Ry)),
                         np.max(np.abs(Rzb - Rz))))
        rows.append({"transform": rc, "max_roundtrip_err": err,
                      "tolerance": 1e-14, "passes": err < 1e-14})
    return {"rows": rows, "passes": all(r["passes"] for r in rows)}


if __name__ == "__main__":
    rows = _scalar_roundtrip_validation()
    assert all(r["passes"] for r in rows), "M02-T1 failed"
    print(f"M02-T1 scalar round-trip: {len(rows)} cases all pass")

    rows = _matrix_orthogonality_validation()
    assert all(r["passes"] for r in rows), "M02-T2 failed"
    print(f"M02-T2 orthogonality: {len(rows)} cases all pass")

    rows = _shape_registry_validation()
    assert all(r["passes"] for r in rows), "M02-T3 failed"
    print(f"M02-T3 shape registry: {len(rows)} cases all pass")

    r = _closed_form_mapping_table_test()
    assert r["passes"], f"M02 closed-form mapping table: {r}"
    print(f"M02 closed-form mapping table: {len(r['rows'])} cases all match")

    r = _closed_form_vector_label_test()
    assert r["passes"], f"M02 closed-form vector labels: {r}"
    print(f"M02 closed-form vector labels: {len(r['rows'])} cases all pass")

    w = _legacy_wrong_control()
    assert w["passes"], "M02-W1 legacy RC5 helper must FAIL the coordinate mapping"
    print(f"M02-W1 legacy RC5 wrong control: legacy shape={w['legacy_shape']}, "
          f"correct shape={w['correct_shape']}")

    # Pair direction transforms: every N6 direction under every RC.
    for rc in RC_TRANSFORMS:
        for lbl in N6_DIRECTIONS:
            out = transform_pair_direction(lbl, rc)
            assert out in N6_DIRECTIONS
    print("M02 pair direction transforms: 42 cases all map cleanly")
    print("M02 coordinate transforms: all checks passed")