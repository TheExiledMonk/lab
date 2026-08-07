"""M03 — Vector Transforms.

Scientific contract
-------------------
Array axes are ``(z, y, x)`` while vector components are ``(x, y, z)``.
A passive forward transform therefore has two distinct operations:

1. spatially remap each scalar component field;
2. mix vector components according to the physical coordinate transform.

The production implementation uses the canonical Q matrices from M01.
The independent validation path in this module deliberately does NOT use Q
for component expectations.  It uses a frozen, closed-form component map for
RC0..RC6 and only relies on the already-verified M02 scalar spatial transform.
That separation is important: a shared error in Q must not automatically pass
M03's secondary validation.
"""
from __future__ import annotations

import numpy as np

from .coordinate_transforms import (
    transform_scalar_field,
    inverse_transform_scalar_field,
)
from .conventions import (
    RC_TRANSFORMS,
    get_coordinate_matrix,
    validate_transform_id,
    EPS_TRANSFORM,
)

__all__ = [
    "transform_vector_field",
    "inverse_transform_vector_field",
    "transform_vector_field_reference",
    "inverse_transform_field_reference",
    "transform_vector_field_closed_form",
    "inverse_transform_vector_field_closed_form",
    "scalar_only_inverse_wrong_control",
    "VectorTransformsError",
]


class VectorTransformsError(ValueError):
    pass


# ----------------------------------------------------------------------
# Independent closed-form component registry.
# ----------------------------------------------------------------------
# Each new component is expressed as (source_component, sign) after the
# scalar spatial transform has been applied.  This table is intentionally
# literal and MUST NOT be generated from the production Q matrices.
_COMPONENT_MAP_FWD = {
    "RC0": ((0, +1.0), (1, +1.0), (2, +1.0)),
    "RC1": ((1, +1.0), (0, +1.0), (2, +1.0)),
    "RC2": ((2, +1.0), (1, +1.0), (0, +1.0)),
    "RC3": ((0, +1.0), (2, +1.0), (1, +1.0)),
    # +90 deg about x: (Rx, Ry, Rz) -> (Rx, -Rz, Ry)
    "RC4": ((0, +1.0), (2, -1.0), (1, +1.0)),
    # +90 deg about y: (Rx, Ry, Rz) -> (Rz, Ry, -Rx)
    "RC5": ((2, +1.0), (1, +1.0), (0, -1.0)),
    # +90 deg about z: (Rx, Ry, Rz) -> (-Ry, Rx, Rz)
    "RC6": ((1, -1.0), (0, +1.0), (2, +1.0)),
}

# Explicit inverse component maps.  These are written independently rather
# than generated from the forward table so an accidental forward-table edit
# does not silently rewrite the inverse expectation.
_COMPONENT_MAP_INV = {
    "RC0": ((0, +1.0), (1, +1.0), (2, +1.0)),
    "RC1": ((1, +1.0), (0, +1.0), (2, +1.0)),
    "RC2": ((2, +1.0), (1, +1.0), (0, +1.0)),
    "RC3": ((0, +1.0), (2, +1.0), (1, +1.0)),
    # inverse of RC4: (Rx', Ry', Rz') -> (Rx', Rz', -Ry')
    "RC4": ((0, +1.0), (2, +1.0), (1, -1.0)),
    # inverse of RC5: (Rx', Ry', Rz') -> (-Rz', Ry', Rx')
    "RC5": ((2, -1.0), (1, +1.0), (0, +1.0)),
    # inverse of RC6: (Rx', Ry', Rz') -> (Ry', -Rx', Rz')
    "RC6": ((1, +1.0), (0, -1.0), (2, +1.0)),
}


def _validate_component_shapes(Rx, Ry, Rz):
    Rx = np.asarray(Rx, dtype=np.float64)
    Ry = np.asarray(Ry, dtype=np.float64)
    Rz = np.asarray(Rz, dtype=np.float64)
    if Rx.shape != Ry.shape or Rx.shape != Rz.shape:
        raise VectorTransformsError("Rx, Ry, Rz must share the same shape")
    if Rx.ndim != 3:
        raise VectorTransformsError("vector components must be 3-D fields")
    return Rx, Ry, Rz


# ----------------------------------------------------------------------
# Production implementation.
# ----------------------------------------------------------------------
def transform_vector_field(Rx, Ry, Rz, transform_id):
    """Forward vector transform: scalar spatial map followed by Q mixing."""
    validate_transform_id(transform_id)
    Rx, Ry, Rz = _validate_component_shapes(Rx, Ry, Rz)

    Rx_s = transform_scalar_field(Rx, transform_id)
    Ry_s = transform_scalar_field(Ry, transform_id)
    Rz_s = transform_scalar_field(Rz, transform_id)

    Q = get_coordinate_matrix(transform_id, inverse=False)
    Rxp = Q[0, 0] * Rx_s + Q[0, 1] * Ry_s + Q[0, 2] * Rz_s
    Ryp = Q[1, 0] * Rx_s + Q[1, 1] * Ry_s + Q[1, 2] * Rz_s
    Rzp = Q[2, 0] * Rx_s + Q[2, 1] * Ry_s + Q[2, 2] * Rz_s
    return Rxp, Ryp, Rzp


def inverse_transform_vector_field(Rxp, Ryp, Rzp, transform_id):
    """Inverse vector transform: Q^T mixing followed by inverse spatial map."""
    validate_transform_id(transform_id)
    Rxp, Ryp, Rzp = _validate_component_shapes(Rxp, Ryp, Rzp)

    Q = get_coordinate_matrix(transform_id, inverse=False)
    Rx_s = Q[0, 0] * Rxp + Q[1, 0] * Ryp + Q[2, 0] * Rzp
    Ry_s = Q[0, 1] * Rxp + Q[1, 1] * Ryp + Q[2, 1] * Rzp
    Rz_s = Q[0, 2] * Rxp + Q[1, 2] * Ryp + Q[2, 2] * Rzp

    return (
        inverse_transform_scalar_field(Rx_s, transform_id),
        inverse_transform_scalar_field(Ry_s, transform_id),
        inverse_transform_scalar_field(Rz_s, transform_id),
    )


# ----------------------------------------------------------------------
# Independent closed-form validation implementation.
# ----------------------------------------------------------------------
def transform_vector_field_closed_form(Rx, Ry, Rz, transform_id):
    """Independent forward expectation with no production-Q dependency.

    M02 is responsible for the spatial remapping.  M03 independently checks
    only the component transformation using the literal RC component table
    above.  This prevents a wrong Q matrix from validating itself.
    """
    validate_transform_id(transform_id)
    Rx, Ry, Rz = _validate_component_shapes(Rx, Ry, Rz)

    spatial = (
        transform_scalar_field(Rx, transform_id),
        transform_scalar_field(Ry, transform_id),
        transform_scalar_field(Rz, transform_id),
    )
    mapping = _COMPONENT_MAP_FWD[transform_id]
    return tuple(sign * spatial[src] for src, sign in mapping)


def inverse_transform_vector_field_closed_form(Rxp, Ryp, Rzp, transform_id):
    """Independent inverse expectation with no production-Q dependency."""
    validate_transform_id(transform_id)
    Rxp, Ryp, Rzp = _validate_component_shapes(Rxp, Ryp, Rzp)
    transformed = (Rxp, Ryp, Rzp)
    mapping = _COMPONENT_MAP_INV[transform_id]
    spatial_native = tuple(sign * transformed[src] for src, sign in mapping)
    return tuple(
        inverse_transform_scalar_field(comp, transform_id)
        for comp in spatial_native
    )


# Backward-compatible names.  They now point to the genuinely independent
# closed-form route rather than another Q-based implementation.
def transform_vector_field_reference(Rx, Ry, Rz, transform_id):
    return transform_vector_field_closed_form(Rx, Ry, Rz, transform_id)


def inverse_transform_field_reference(Rxp, Ryp, Rzp, transform_id):
    return inverse_transform_vector_field_closed_form(Rxp, Ryp, Rzp, transform_id)


# ----------------------------------------------------------------------
# Wrong control: predecessor's scalar-only inverse.
# ----------------------------------------------------------------------
def scalar_only_inverse_wrong_control(Rxp, Ryp, Rzp, transform_id):
    """Broken predecessor control: spatial inverse without component unmixing."""
    return (
        inverse_transform_scalar_field(Rxp, transform_id),
        inverse_transform_scalar_field(Ryp, transform_id),
        inverse_transform_scalar_field(Rzp, transform_id),
    )


# ----------------------------------------------------------------------
# Self-checks.
# ----------------------------------------------------------------------
def _max_component_error(a, b):
    return float(max(np.max(np.abs(x - y)) for x, y in zip(a, b)))


def _basis_vector_tests() -> list:
    """Constant basis vectors and a nonsymmetric varying field."""
    nz, ny, nx = 3, 4, 5
    Z, Y, X = np.meshgrid(
        np.arange(nz), np.arange(ny), np.arange(nx), indexing="ij"
    )

    cases = [
        ("V1_ex_unit", (
            np.ones((nz, ny, nx)),
            np.zeros((nz, ny, nx)),
            np.zeros((nz, ny, nx)),
        )),
        ("V2_ey_unit", (
            np.zeros((nz, ny, nx)),
            np.ones((nz, ny, nx)),
            np.zeros((nz, ny, nx)),
        )),
        ("V3_ez_unit", (
            np.zeros((nz, ny, nx)),
            np.zeros((nz, ny, nx)),
            np.ones((nz, ny, nx)),
        )),
        ("V4_varying", (
            1000.0 + 100.0 * Z + 10.0 * Y + X,
            2000.0 + 100.0 * Z + 10.0 * Y + X,
            3000.0 + 100.0 * Z + 10.0 * Y + X,
        )),
    ]

    rows = []
    for name, field in cases:
        for rc in RC_TRANSFORMS:
            prod = transform_vector_field(*field, rc)
            closed = transform_vector_field_closed_form(*field, rc)
            forward_err = _max_component_error(prod, closed)

            back = inverse_transform_vector_field(*prod, rc)
            roundtrip_err = _max_component_error(back, field)

            closed_back = inverse_transform_vector_field_closed_form(*prod, rc)
            inverse_ref_err = _max_component_error(closed_back, field)

            rows.append({
                "test": "closed_form_component_mapping",
                "transform": rc,
                "field": name,
                "max_forward_diff": forward_err,
                "max_roundtrip_error": roundtrip_err,
                "max_closed_form_inverse_error": inverse_ref_err,
                "tolerance": EPS_TRANSFORM,
                "passes": (
                    forward_err <= EPS_TRANSFORM
                    and roundtrip_err <= EPS_TRANSFORM
                    and inverse_ref_err <= EPS_TRANSFORM
                ),
            })
    return rows


def _symbolic_component_mapping_tests() -> list:
    """Directly verify the expected component identities for every RC.

    The three native component fields are deliberately separated by large
    offsets (1000, 2000, 3000).  A component swap or sign error therefore
    cannot hide behind similar numerical values.
    """
    nz, ny, nx = 3, 4, 5
    Z, Y, X = np.meshgrid(
        np.arange(nz), np.arange(ny), np.arange(nx), indexing="ij"
    )
    native = (
        1000.0 + 100.0 * Z + 10.0 * Y + X,
        2000.0 + 100.0 * Z + 10.0 * Y + X,
        3000.0 + 100.0 * Z + 10.0 * Y + X,
    )

    rows = []
    for rc in RC_TRANSFORMS:
        prod = transform_vector_field(*native, rc)
        spatial = tuple(transform_scalar_field(c, rc) for c in native)
        mapping = _COMPONENT_MAP_FWD[rc]
        expected = tuple(sign * spatial[src] for src, sign in mapping)
        err = _max_component_error(prod, expected)
        rows.append({
            "test": "static_symbolic_component_mapping",
            "transform": rc,
            "max_error": err,
            "tolerance": EPS_TRANSFORM,
            "passes": err <= EPS_TRANSFORM,
        })
    return rows


def _wrong_control_test() -> list:
    """Scalar-only inverse must reproduce the predecessor's large error."""
    nz, ny, nx = 9, 16, 17
    Z, Y, X = np.meshgrid(
        np.arange(nz), np.arange(ny), np.arange(nx), indexing="ij"
    )
    Rx = np.sin(2 * np.pi * X / 11.0) * np.cos(2 * np.pi * Y / 7.0) * (0.6 + 0.4 * Z / 8.0)
    Ry = np.cos(2 * np.pi * X / 13.0) * np.sin(2 * np.pi * Y / 9.0) * (0.4 + 0.6 * Z / 8.0)
    Rz = 0.3 * np.sin(2 * np.pi * (X + Y) / 15.0) * np.cos(2 * np.pi * Z / 5.0)
    native = (Rx, Ry, Rz)
    norm_native = float(np.sqrt(np.sum(Rx ** 2 + Ry ** 2 + Rz ** 2)))

    rows = []
    for rc in RC_TRANSFORMS:
        transformed = transform_vector_field(*native, rc)
        wrong = scalar_only_inverse_wrong_control(*transformed, rc)
        correct = inverse_transform_vector_field(*transformed, rc)

        wrong_norm = float(np.sqrt(sum(np.sum((a - b) ** 2) for a, b in zip(wrong, native))))
        correct_norm = float(np.sqrt(sum(np.sum((a - b) ** 2) for a, b in zip(correct, native))))
        E_wrong = wrong_norm / max(norm_native, 1e-15)
        E_correct = correct_norm / max(norm_native, 1e-15)

        if rc == "RC0":
            passes = E_wrong < 1e-12 and E_correct < 1e-12
        else:
            passes = E_wrong > 0.3 and E_correct < 1e-12

        rows.append({
            "test": "scalar_only_inverse_wrong_control",
            "transform": rc,
            "WR_C1_scalar_only_E_cov": E_wrong,
            "WR_C2_correct_E_cov": E_correct,
            "passes": passes,
        })
    return rows


if __name__ == "__main__":
    rows = _basis_vector_tests()
    assert all(r["passes"] for r in rows), f"M03 closed-form tests failed: {rows}"
    print(f"M03 closed-form vector validation: {len(rows)} cases pass")

    rows = _symbolic_component_mapping_tests()
    assert all(r["passes"] for r in rows), f"M03 symbolic mapping failed: {rows}"
    print(f"M03 static symbolic component mapping: {len(rows)} RC cases pass")

    rows = _wrong_control_test()
    assert all(r["passes"] for r in rows), f"M03 wrong control failed: {rows}"
    print("M03 wrong control reproduces scalar-only inverse failure")
    print("M03 vector transforms: all checks passed")
