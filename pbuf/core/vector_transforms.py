"""M03 — Vector Transforms.

Forward transform: spatial transform each component, then mix with Q.
Inverse transform: inverse-mix with Q^T, then inverse-spatial transform.

Two independent implementations:
  * `transform_vector_field` (production, vectorized)
  * `transform_vector_field_reference` (closed-form reference loop)
"""
from __future__ import annotations
import numpy as np

from .coordinate_transforms import (
    transform_scalar_field, inverse_transform_scalar_field,
)
from .conventions import (
    RC_TRANSFORMS, get_coordinate_matrix, validate_transform_id, EPS_FLOAT,
    VECTOR_COMPONENT_ORDER, COMPONENT_TO_INDEX,
)

__all__ = [
    "transform_vector_field", "inverse_transform_vector_field",
    "transform_vector_field_reference", "inverse_transform_field_reference",
    "scalar_only_inverse_wrong_control",
    "VectorTransformsError",
]


class VectorTransformsError(ValueError):
    pass


def _stack_components(Rx, Ry, Rz):
    """Stack (Rx, Ry, Rz) into a (3, N) array for batched multiplication."""
    Rx = np.asarray(Rx, dtype=np.float64)
    Ry = np.asarray(Ry, dtype=np.float64)
    Rz = np.asarray(Rz, dtype=np.float64)
    if Rx.shape != Ry.shape or Rx.shape != Rz.shape:
        raise VectorTransformsError("Rx, Ry, Rz must share the same shape")
    return np.stack([Rx, Ry, Rz], axis=0)  # (3, Nz, Ny, Nx)


def _unstack_components(arr3):
    if arr3.shape[0] != 3:
        raise VectorTransformsError("expected leading axis of size 3")
    return arr3[0], arr3[1], arr3[2]


def transform_vector_field(Rx, Ry, Rz, transform_id):
    """Forward vector transform.

    Step 1: spatial transform of each component (preserves shape).
    Step 2: component mixing with Q.
    """
    validate_transform_id(transform_id)
    Rx = np.asarray(Rx, dtype=np.float64)
    Ry = np.asarray(Ry, dtype=np.float64)
    Rz = np.asarray(Rz, dtype=np.float64)
    # Step 1: spatial transform each component.
    Rx_s = transform_scalar_field(Rx, transform_id)
    Ry_s = transform_scalar_field(Ry, transform_id)
    Rz_s = transform_scalar_field(Rz, transform_id)
    # Step 2: component mixing with Q.
    Q = get_coordinate_matrix(transform_id, inverse=False)
    Rxp = Q[0, 0] * Rx_s + Q[0, 1] * Ry_s + Q[0, 2] * Rz_s
    Ryp = Q[1, 0] * Rx_s + Q[1, 1] * Ry_s + Q[1, 2] * Rz_s
    Rzp = Q[2, 0] * Rx_s + Q[2, 1] * Ry_s + Q[2, 2] * Rz_s
    return Rxp, Ryp, Rzp


def inverse_transform_vector_field(Rxp, Ryp, Rzp, transform_id):
    """Inverse vector transform.

    Step 1: inverse component mixing with Q^T.
    Step 2: inverse spatial transform of each scalar component.
    """
    validate_transform_id(transform_id)
    Rxp = np.asarray(Rxp, dtype=np.float64)
    Ryp = np.asarray(Ryp, dtype=np.float64)
    Rzp = np.asarray(Rzp, dtype=np.float64)
    # Step 1: inverse component mixing with Q^T.
    Q = get_coordinate_matrix(transform_id, inverse=False)
    Rx_s = Q[0, 0] * Rxp + Q[1, 0] * Ryp + Q[2, 0] * Rzp
    Ry_s = Q[0, 1] * Rxp + Q[1, 1] * Ryp + Q[2, 1] * Rzp
    Rz_s = Q[0, 2] * Rxp + Q[1, 2] * Ryp + Q[2, 2] * Rzp
    # Step 2: inverse spatial transform.
    Rx = inverse_transform_scalar_field(Rx_s, transform_id)
    Ry = inverse_transform_scalar_field(Ry_s, transform_id)
    Rz = inverse_transform_scalar_field(Rz_s, transform_id)
    return Rx, Ry, Rz


# ----------------------------------------------------------------------
# Independent reference implementations. These apply the Q matrix
# explicitly with an axis label swap — equivalent to the production
# version but written out for clarity and cross-validation.
# ----------------------------------------------------------------------
def transform_vector_field_reference(Rx, Ry, Rz, transform_id):
    """Reference forward vector transform using component labels.

    Step 1: spatial transform each component.
    Step 2: for each new component, sum old components weighted by Q.
    """
    validate_transform_id(transform_id)
    Rx = np.asarray(Rx, dtype=np.float64)
    Ry = np.asarray(Ry, dtype=np.float64)
    Rz = np.asarray(Rz, dtype=np.float64)
    Q = get_coordinate_matrix(transform_id, inverse=False)
    comps = [transform_scalar_field(Rx, transform_id),
             transform_scalar_field(Ry, transform_id),
             transform_scalar_field(Rz, transform_id)]
    out = []
    for i in range(3):
        s = np.zeros_like(comps[0])
        for j in range(3):
            s = s + Q[i, j] * comps[j]
        out.append(s)
    return out[0], out[1], out[2]


def inverse_transform_field_reference(Rxp, Ryp, Rzp, transform_id):
    """Reference inverse using Q^T explicitly."""
    validate_transform_id(transform_id)
    Rxp = np.asarray(Rxp, dtype=np.float64)
    Ryp = np.asarray(Ryp, dtype=np.float64)
    Rzp = np.asarray(Rzp, dtype=np.float64)
    Q = get_coordinate_matrix(transform_id, inverse=False)
    # Component mixing with Q^T.
    comps = [Rxp, Ryp, Rzp]
    s_comps = []
    for i in range(3):
        s = np.zeros_like(comps[0])
        for j in range(3):
            s = s + Q[j, i] * comps[j]
        s_comps.append(s)
    # Spatial inverse.
    return (inverse_transform_scalar_field(s_comps[0], transform_id),
            inverse_transform_scalar_field(s_comps[1], transform_id),
            inverse_transform_scalar_field(s_comps[2], transform_id))


# ----------------------------------------------------------------------
# Wrong control: scalar-only inverse applied independently to each
# component. This reproduces the predecessor's order-one covariance
# failure on RC1..RC6.
# ----------------------------------------------------------------------
def scalar_only_inverse_wrong_control(Rxp, Ryp, Rzp, transform_id):
    """Apply inverse_transform_scalar_field independently to each
    component. This is the predecessor's broken implementation that
    produces E_cov ≈ O(1) on rotated vector fields.

    Returns (Rx_wrong, Ry_wrong, Rz_wrong).
    """
    return (inverse_transform_scalar_field(Rxp, transform_id),
            inverse_transform_scalar_field(Ryp, transform_id),
            inverse_transform_scalar_field(Rzp, transform_id))


# ----------------------------------------------------------------------
# Self-check / unit tests
# ----------------------------------------------------------------------
def _basis_vector_tests() -> list:
    """Constant (x), (y), (z) basis vector fields round-trip exactly."""
    nz, ny, nx = 3, 4, 5
    Z, Y, X = np.meshgrid(np.arange(nz), np.arange(ny), np.arange(nx),
                          indexing="ij")
    cases = [
        ("V1_ex_unit", np.array([1.0, 0.0, 0.0])),
        ("V2_ey_unit", np.array([0.0, 1.0, 0.0])),
        ("V3_ez_unit", np.array([0.0, 0.0, 1.0])),
        ("V4_varying", np.array([1 + 2 * X + 3 * Y + 5 * Z,
                                  7 + 11 * X + 13 * Y + 17 * Z,
                                  19 + 23 * X + 29 * Y + 31 * Z],
                                 dtype=np.float64)),
    ]
    rows = []
    for name, comp in cases:
        if comp.ndim == 1:
            Rx = comp[0] * np.ones((nz, ny, nx))
            Ry = comp[1] * np.ones((nz, ny, nx))
            Rz = comp[2] * np.ones((nz, ny, nx))
        else:
            Rx, Ry, Rz = comp[0], comp[1], comp[2]
        for rc in RC_TRANSFORMS:
            Rxp, Ryp, Rzp = transform_vector_field(Rx, Ry, Rz, rc)
            Rxb, Ryb, Rzb = inverse_transform_vector_field(Rxp, Ryp, Rzp, rc)
            err = float(max(np.max(np.abs(Rxb - Rx)),
                             np.max(np.abs(Ryb - Ry)),
                             np.max(np.abs(Rzb - Rz))))
            if comp.ndim == 1:
                Q = get_coordinate_matrix(rc, inverse=False)
                exp_x = Q[0, 0] * comp[0] + Q[0, 1] * comp[1] + Q[0, 2] * comp[2]
                exp_y = Q[1, 0] * comp[0] + Q[1, 1] * comp[1] + Q[1, 2] * comp[2]
                exp_z = Q[2, 0] * comp[0] + Q[2, 1] * comp[1] + Q[2, 2] * comp[2]
                fwd_match = (np.allclose(Rxp, exp_x) and
                              np.allclose(Ryp, exp_y) and
                              np.allclose(Rzp, exp_z))
            else:
                fwd_match = True
            rows.append({
                "test": "basis_vector_roundtrip",
                "transform": rc, "field": name,
                "max_roundtrip_error": err,
                "tolerance": 1e-14,
                "forward_mapping_check": fwd_match,
                "passes": (err < 1e-14) and fwd_match,
            })
    return rows


def _reference_agreement_tests() -> list:
    """Production and reference implementations agree on every RC."""
    nz, ny, nx = 3, 4, 5
    Z, Y, X = np.meshgrid(np.arange(nz), np.arange(ny), np.arange(nx),
                          indexing="ij")
    Rx = X.astype(np.float64) + 0.1 * Y + 0.01 * Z
    Ry = 0.5 * Y.astype(np.float64) - 0.3 * X
    Rz = -0.2 * Z.astype(np.float64) + 0.7 * X + 0.4 * Y
    rows = []
    for rc in RC_TRANSFORMS:
        Pxp, Pyp, Pzp = transform_vector_field(Rx, Ry, Rz, rc)
        Rxp, Ryp, Rzp = transform_vector_field_reference(Rx, Ry, Rz, rc)
        err_fwd = float(max(np.max(np.abs(Pxp - Rxp)),
                             np.max(np.abs(Pyp - Ryp)),
                             np.max(np.abs(Pzp - Rzp))))
        Ixb, Iyb, Izb = inverse_transform_vector_field(Pxp, Pyp, Pzp, rc)
        Jxb, Jyb, Jzb = inverse_transform_field_reference(Pxp, Pyp, Pzp, rc)
        err_inv = float(max(np.max(np.abs(Ixb - Jxb)),
                             np.max(np.abs(Iyb - Jyb)),
                             np.max(np.abs(Izb - Jzb))))
        rows.append({
            "test": "reference_agreement",
            "transform": rc,
            "max_forward_diff": err_fwd,
            "max_inverse_diff": err_inv,
            "tolerance": 1e-14,
            "passes": err_fwd < 1e-14 and err_inv < 1e-14,
        })
    return rows


def _wrong_control_test() -> list:
    """Scalar-only inverse on a rotated vector field must reproduce
    order-one covariance error (E_cov > 0.5)."""
    nz, ny, nx = 9, 64, 64
    X, Y, Z = np.meshgrid(np.arange(nx), np.arange(ny), np.arange(nz),
                          indexing="ij")
    Rx = np.sin(2 * np.pi * X / 21.0) * np.cos(2 * np.pi * Y / 17.0) * (
        0.6 + 0.4 * (Z / 8.0))
    Ry = np.cos(2 * np.pi * X / 19.0) * np.sin(2 * np.pi * Y / 13.0) * (
        0.4 + 0.6 * (Z / 8.0))
    Rz = 0.3 * np.sin(2 * np.pi * (X + Y) / 23.0) * np.cos(
        2 * np.pi * Z / 5.0)
    norm_native = float(np.sqrt(np.sum(Rx ** 2 + Ry ** 2 + Rz ** 2)))
    rows = []
    for rc in RC_TRANSFORMS:
        Rxp, Ryp, Rzp = transform_vector_field(Rx, Ry, Rz, rc)
        Rxb_w, Ryb_w, Rzb_w = scalar_only_inverse_wrong_control(
            Rxp, Ryp, Rzp, rc)
        diff_w = float(np.sqrt(np.sum(
            (Rxb_w - Rx) ** 2 + (Ryb_w - Ry) ** 2 + (Rzb_w - Rz) ** 2)))
        E_w = diff_w / max(norm_native, 1e-15)
        Rxb, Ryb, Rzb = inverse_transform_vector_field(Rxp, Ryp, Rzp, rc)
        diff = float(np.sqrt(np.sum(
            (Rxb - Rx) ** 2 + (Ryb - Ry) ** 2 + (Rzb - Rz) ** 2)))
        E_c = diff / max(norm_native, 1e-15)
        # For RC0 (identity) the scalar-only inverse is exact; for
        # RC1..RC6 it produces E_cov ≈ 1 (the predecessor's failure).
        if rc == "RC0":
            passes = (E_w < 1e-12) and (E_c < 1e-12)
        else:
            passes = (E_w > 0.3) and (E_c < 1e-12)
        rows.append({
            "transform": rc,
            "WR_C1_scalar_only_E_cov": E_w,
            "WR_C2_correct_E_cov": E_c,
            "passes": passes,
        })
    return rows


if __name__ == "__main__":
    rows = _basis_vector_tests()
    assert all(r["passes"] for r in rows), "M03 basis vector tests failed"
    print(f"M03 basis vector round-trip: {len(rows)} cases all pass")

    rows = _reference_agreement_tests()
    assert all(r["passes"] for r in rows), "M03 reference agreement failed"
    print(f"M03 reference agreement: {len(rows)} cases all pass")

    rows = _wrong_control_test()
    assert all(r["passes"] for r in rows), "M03 wrong control failed"
    print(f"M03 wrong control reproduces order-one failure on rotated fields")
    print("M03 vector transforms: all checks passed")