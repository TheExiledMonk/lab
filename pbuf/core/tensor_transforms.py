"""M04 — Tensor Transforms.

Symmetric tensor transform with the canonical identity:
    P' = Q P Q^T

Correction pass FOUNDATION-001-CORRECTION-001
---------------------------------------------
* Inverse round-trip is enforced for EVERY RC0..RC6 transform with no
  waivers (the previous test silently passed rotations despite
  boundary-cell issues).
* The reference implementation uses an explicit per-voxel matrix
  loop with a Python-level matrix product ``Q @ P_native @ Q.T``,
  NOT the same ``einsum`` as production.
* Projector identity test for P_T = I - ê_L ê_L^T under transforms.
"""
from __future__ import annotations
import numpy as np

from .coordinate_transforms import (
    transform_scalar_field, inverse_transform_scalar_field,
)
from .conventions import (
    RC_TRANSFORMS, get_coordinate_matrix, validate_transform_id,
    TENSOR_COMPONENT_ORDER,
)

__all__ = [
    "transform_symmetric_tensor_field",
    "inverse_transform_symmetric_tensor_field",
    "transform_symmetric_tensor_field_reference",
    "inverse_transform_symmetric_tensor_field_reference",
    "TensorTransformsError",
]


class TensorTransformsError(ValueError):
    pass


def _stack_symmetric(Pxx, Pxy, Pxz, Pyy, Pyz, Pzz):
    arr = np.array([
        [np.asarray(Pxx, dtype=np.float64),
         np.asarray(Pxy, dtype=np.float64),
         np.asarray(Pxz, dtype=np.float64)],
        [np.asarray(Pxy, dtype=np.float64),
         np.asarray(Pyy, dtype=np.float64),
         np.asarray(Pyz, dtype=np.float64)],
        [np.asarray(Pxz, dtype=np.float64),
         np.asarray(Pyz, dtype=np.float64),
         np.asarray(Pzz, dtype=np.float64)],
    ], dtype=np.float64)
    return arr


def transform_symmetric_tensor_field(Pxx, Pxy, Pxz, Pyy, Pyz, Pzz,
                                        transform_id):
    """Forward transform of a 3x3 symmetric tensor field.

    Step 1: spatial transform of each component.
    Step 2: P' = Q P Q^T (per-voxel matmul).
    """
    validate_transform_id(transform_id)
    Pxx_s = transform_scalar_field(Pxx, transform_id)
    Pxy_s = transform_scalar_field(Pxy, transform_id)
    Pxz_s = transform_scalar_field(Pxz, transform_id)
    Pyy_s = transform_scalar_field(Pyy, transform_id)
    Pyz_s = transform_scalar_field(Pyz, transform_id)
    Pzz_s = transform_scalar_field(Pzz, transform_id)
    Q = get_coordinate_matrix(transform_id, inverse=False)
    P = _stack_symmetric(Pxx_s, Pxy_s, Pxz_s, Pyy_s, Pyz_s, Pzz_s)
    P_p = np.einsum('ia,abklm,jb->ijklm', Q, P, Q)
    return (P_p[0, 0], P_p[0, 1], P_p[0, 2],
            P_p[1, 1], P_p[1, 2], P_p[2, 2])


def inverse_transform_symmetric_tensor_field(Pxx_p, Pxy_p, Pxz_p,
                                                Pyy_p, Pyz_p, Pzz_p,
                                                transform_id):
    """Inverse: P = Q^T P_p Q."""
    validate_transform_id(transform_id)
    P_p = _stack_symmetric(Pxx_p, Pxy_p, Pxz_p, Pyy_p, Pyz_p, Pzz_p)
    Q = get_coordinate_matrix(transform_id, inverse=False)
    # CORRECTED (CORRECTION-001): einsum must use 'bj' for the second Q
    # so the matrix multiplication reads Q[b, j] (matching Q @ (Q.T @ P_p)),
    # not Q[j, b].
    P = np.einsum('ia,abklm,bj->ijklm', Q.T, P_p, Q)
    Pxx_s = P[0, 0]; Pxy_s = P[0, 1]; Pxz_s = P[0, 2]
    Pyy_s = P[1, 1]; Pyz_s = P[1, 2]; Pzz_s = P[2, 2]
    Pxx = inverse_transform_scalar_field(Pxx_s, transform_id)
    Pxy = inverse_transform_scalar_field(Pxy_s, transform_id)
    Pxz = inverse_transform_scalar_field(Pxz_s, transform_id)
    Pyy = inverse_transform_scalar_field(Pyy_s, transform_id)
    Pyz = inverse_transform_scalar_field(Pyz_s, transform_id)
    Pzz = inverse_transform_scalar_field(Pzz_s, transform_id)
    return Pxx, Pxy, Pxz, Pyy, Pyz, Pzz


# ----------------------------------------------------------------------
# Independent reference — explicit per-voxel matrix loop.
# This does NOT use einsum; it builds a Python-level 3x3 matrix and
# performs the multiplication with ``@``.
# ----------------------------------------------------------------------
def transform_symmetric_tensor_field_reference(Pxx, Pxy, Pxz, Pyy, Pyz, Pzz,
                                                  transform_id):
    """Reference forward transform using explicit index-mixing loop
    and a Python-level ``Q @ P_native @ Q.T`` per voxel.

    The reference does NOT use einsum. The algebra is identical to
    production but the code path is genuinely independent.
    """
    validate_transform_id(transform_id)
    Q = get_coordinate_matrix(transform_id, inverse=False)
    Pxx_s = transform_scalar_field(Pxx, transform_id)
    Pxy_s = transform_scalar_field(Pxy, transform_id)
    Pxz_s = transform_scalar_field(Pxz, transform_id)
    Pyy_s = transform_scalar_field(Pyy, transform_id)
    Pyz_s = transform_scalar_field(Pyz, transform_id)
    Pzz_s = transform_scalar_field(Pzz, transform_id)
    nz, ny, nx = Pxx_s.shape
    Pxx_p = np.zeros_like(Pxx_s)
    Pxy_p = np.zeros_like(Pxx_s)
    Pxz_p = np.zeros_like(Pxx_s)
    Pyy_p = np.zeros_like(Pxx_s)
    Pyz_p = np.zeros_like(Pxx_s)
    Pzz_p = np.zeros_like(Pxx_s)
    for iz in range(nz):
        for iy in range(ny):
            for ix in range(nx):
                P_native = np.array([
                    [Pxx_s[iz, iy, ix], Pxy_s[iz, iy, ix], Pxz_s[iz, iy, ix]],
                    [Pxy_s[iz, iy, ix], Pyy_s[iz, iy, ix], Pyz_s[iz, iy, ix]],
                    [Pxz_s[iz, iy, ix], Pyz_s[iz, iy, ix], Pzz_s[iz, iy, ix]],
                ], dtype=np.float64)
                P_out = Q @ P_native @ Q.T
                Pxx_p[iz, iy, ix] = P_out[0, 0]
                Pxy_p[iz, iy, ix] = P_out[0, 1]
                Pxz_p[iz, iy, ix] = P_out[0, 2]
                Pyy_p[iz, iy, ix] = P_out[1, 1]
                Pyz_p[iz, iy, ix] = P_out[1, 2]
                Pzz_p[iz, iy, ix] = P_out[2, 2]
    return Pxx_p, Pxy_p, Pxz_p, Pyy_p, Pyz_p, Pzz_p


def inverse_transform_symmetric_tensor_field_reference(Pxx_p, Pxy_p, Pxz_p,
                                                          Pyy_p, Pyz_p, Pzz_p,
                                                          transform_id):
    """Reference inverse using explicit per-voxel matrix loop and Q^T."""
    validate_transform_id(transform_id)
    Q = get_coordinate_matrix(transform_id, inverse=False)
    nz, ny, nx = Pxx_p.shape
    Pxx_s = np.zeros_like(Pxx_p)
    Pxy_s = np.zeros_like(Pxx_p)
    Pxz_s = np.zeros_like(Pxx_p)
    Pyy_s = np.zeros_like(Pxx_p)
    Pyz_s = np.zeros_like(Pxx_p)
    Pzz_s = np.zeros_like(Pxx_p)
    for iz in range(nz):
        for iy in range(ny):
            for ix in range(nx):
                P_native = np.array([
                    [Pxx_p[iz, iy, ix], Pxy_p[iz, iy, ix], Pxz_p[iz, iy, ix]],
                    [Pxy_p[iz, iy, ix], Pyy_p[iz, iy, ix], Pyz_p[iz, iy, ix]],
                    [Pxz_p[iz, iy, ix], Pyz_p[iz, iy, ix], Pzz_p[iz, iy, ix]],
                ], dtype=np.float64)
                P_out = Q.T @ P_native @ Q
                Pxx_s[iz, iy, ix] = P_out[0, 0]
                Pxy_s[iz, iy, ix] = P_out[0, 1]
                Pxz_s[iz, iy, ix] = P_out[0, 2]
                Pyy_s[iz, iy, ix] = P_out[1, 1]
                Pyz_s[iz, iy, ix] = P_out[1, 2]
                Pzz_s[iz, iy, ix] = P_out[2, 2]
    Pxx = inverse_transform_scalar_field(Pxx_s, transform_id)
    Pxy = inverse_transform_scalar_field(Pxy_s, transform_id)
    Pxz = inverse_transform_scalar_field(Pxz_s, transform_id)
    Pyy = inverse_transform_scalar_field(Pyy_s, transform_id)
    Pyz = inverse_transform_scalar_field(Pyz_s, transform_id)
    Pzz = inverse_transform_scalar_field(Pzz_s, transform_id)
    return Pxx, Pxy, Pxz, Pyy, Pyz, Pzz


# ----------------------------------------------------------------------
# Tests (CORRECTION-001 §9.4)
# ----------------------------------------------------------------------
def _tensor_roundtrip_validation():
    nz, ny, nx = 3, 4, 5
    Z, Y, X = np.meshgrid(np.arange(nz), np.arange(ny), np.arange(nx),
                          indexing="ij")
    rows = []

    def tensor_id_components():
        return (np.ones_like(X, dtype=np.float64),
                np.zeros_like(X), np.zeros_like(X),
                np.ones_like(X), np.zeros_like(X),
                np.ones_like(X))

    def diagonal_aniso_components():
        return (2.0 + 0.1 * X, np.zeros_like(X), np.zeros_like(X),
                1.0 + 0.2 * Y, np.zeros_like(X),
                3.0 + 0.05 * Z)

    def varying_components():
        return (1 + 0.1 * X + 0.2 * Y + 0.3 * Z,
                0.05 * (X * Y) - 0.1 * (Y * Z),
                0.03 * (X * Z),
                1 + 0.4 * Y - 0.1 * Z,
                0.07 * X + 0.05 * Y,
                1 + 0.6 * Z)

    def transverse_projector_components():
        eL_x = np.full_like(X, 0.6, dtype=np.float64)
        eL_y = np.full_like(X, 0.8, dtype=np.float64)
        eL_z = np.zeros_like(X, dtype=np.float64)
        norm = np.sqrt(eL_x ** 2 + eL_y ** 2 + eL_z ** 2)
        eL_x /= norm; eL_y /= norm; eL_z /= norm
        return (1.0 - eL_x * eL_x,
                -eL_x * eL_y,
                -eL_x * eL_z,
                1.0 - eL_y * eL_y,
                -eL_y * eL_z,
                1.0 - eL_z * eL_z)

    def rank_one_components():
        # vv^T with v = (1, 2, 3) (constant).
        v1 = 1.0 + 0.01 * X + 0.02 * Y
        v2 = 2.0 + 0.03 * Y
        v3 = 3.0 + 0.04 * Z
        return (v1 * v1, v1 * v2, v1 * v3,
                v2 * v2, v2 * v3, v3 * v3)

    cases = [
        ("identity", tensor_id_components()),
        ("diag_aniso", diagonal_aniso_components()),
        ("varying", varying_components()),
        ("transverse_projector", transverse_projector_components()),
        ("rank_one", rank_one_components()),
    ]
    for name, comps in cases:
        Pxx, Pxy, Pxz, Pyy, Pyz, Pzz = comps
        for rc in RC_TRANSFORMS:
            Pxxp, Pxyp, Pxzp, Pyyp, Pyzp, Pzzp = transform_symmetric_tensor_field(
                Pxx, Pxy, Pxz, Pyy, Pyz, Pzz, rc)
            Rxxp, Rxyp, Rxzp, Ryyp, Ryzp, Rzzp = transform_symmetric_tensor_field_reference(
                Pxx, Pxy, Pxz, Pyy, Pyz, Pzz, rc)
            # Forward agreement (production vs reference).
            err_ref = float(max(
                np.max(np.abs(Pxxp - Rxxp)),
                np.max(np.abs(Pxyp - Rxyp)),
                np.max(np.abs(Pxzp - Rxzp)),
                np.max(np.abs(Pyyp - Ryyp)),
                np.max(np.abs(Pyzp - Ryzp)),
                np.max(np.abs(Pzzp - Rzzp))))
            # Inverse round-trip — NO WAIVERS for any RC.
            Pxxb, Pxyb, Pxzb, Pyyb, Pyzb, Pzzb = inverse_transform_symmetric_tensor_field(
                Pxxp, Pxyp, Pxzp, Pyyp, Pyzp, Pzzp, rc)
            err_back = float(max(
                np.max(np.abs(Pxxb - Pxx)),
                np.max(np.abs(Pxyb - Pxy)),
                np.max(np.abs(Pxzb - Pxz)),
                np.max(np.abs(Pyyb - Pyy)),
                np.max(np.abs(Pyzb - Pyz)),
                np.max(np.abs(Pzzb - Pzz))))
            rows.append({
                "test": "tensor_roundtrip",
                "case": name, "transform": rc,
                "max_reference_diff": err_ref,
                "max_roundtrip_error": err_back,
                "tolerance": 1e-14,
                "is_rotation": rc in ("RC4", "RC5", "RC6"),
                # CORRECTION-001 §9.6: NO exceptions for rotations.
                "passes": (err_ref < 1e-14 and err_back < 1e-14),
            })
    return rows


def _projector_identity_test():
    """P_T = I - ê_L ê_L^T transforms as Q P_T Q^T = I - (Q ê_L)(Q ê_L)^T.

    The transformed longitudinal unit vector must be Q ê_L exactly.
    """
    nz, ny, nx = 3, 4, 5
    Z, Y, X = np.meshgrid(np.arange(nz), np.arange(ny), np.arange(nx),
                          indexing="ij")
    eL_x = np.full_like(X, 0.6, dtype=np.float64)
    eL_y = np.full_like(X, 0.8, dtype=np.float64)
    eL_z = np.zeros_like(X, dtype=np.float64)
    norm = np.sqrt(eL_x ** 2 + eL_y ** 2 + eL_z ** 2)
    eL_x /= norm; eL_y /= norm; eL_z /= norm
    # Build P_T components.
    Pxx = 1.0 - eL_x * eL_x
    Pxy = -eL_x * eL_y
    Pxz = -eL_x * eL_z
    Pyy = 1.0 - eL_y * eL_y
    Pyz = -eL_y * eL_z
    Pzz = 1.0 - eL_z * eL_z
    rows = []
    for rc in RC_TRANSFORMS:
        Q = get_coordinate_matrix(rc, inverse=False)
        # Reference longitudinal vector: Q ê_L.
        eL_x_t = Q[0, 0] * eL_x + Q[0, 1] * eL_y + Q[0, 2] * eL_z
        eL_y_t = Q[1, 0] * eL_x + Q[1, 1] * eL_y + Q[1, 2] * eL_z
        eL_z_t = Q[2, 0] * eL_x + Q[2, 1] * eL_y + Q[2, 2] * eL_z
        # Transformed projector: full FWD then INV spatial inverse
        # so we compare against the static reference shape.
        Pxx_p, Pxy_p, Pxz_p, Pyy_p, Pyz_p, Pzz_p = transform_symmetric_tensor_field(
            Pxx, Pxy, Pxz, Pyy, Pyz, Pzz, rc)
        Pxx_p = inverse_transform_scalar_field(Pxx_p, rc)
        Pxy_p = inverse_transform_scalar_field(Pxy_p, rc)
        Pxz_p = inverse_transform_scalar_field(Pxz_p, rc)
        Pyy_p = inverse_transform_scalar_field(Pyy_p, rc)
        Pyz_p = inverse_transform_scalar_field(Pyz_p, rc)
        Pzz_p = inverse_transform_scalar_field(Pzz_p, rc)
        # Reference: I - (Q eL)(Q eL)^T.
        Pxx_ref = 1.0 - eL_x_t * eL_x_t
        Pxy_ref = -eL_x_t * eL_y_t
        Pxz_ref = -eL_x_t * eL_z_t
        Pyy_ref = 1.0 - eL_y_t * eL_y_t
        Pyz_ref = -eL_y_t * eL_z_t
        Pzz_ref = 1.0 - eL_z_t * eL_z_t
        err = float(max(
            np.max(np.abs(Pxx_p - Pxx_ref)),
            np.max(np.abs(Pxy_p - Pxy_ref)),
            np.max(np.abs(Pxz_p - Pxz_ref)),
            np.max(np.abs(Pyy_p - Pyy_ref)),
            np.max(np.abs(Pyz_p - Pyz_ref)),
            np.max(np.abs(Pzz_p - Pzz_ref))))
        rows.append({
            "test": "projector_identity",
            "transform": rc,
            "max_err": err,
            "tolerance": 1e-14,
            "passes": err < 1e-14,
        })
    return rows


if __name__ == "__main__":
    rows = _tensor_roundtrip_validation()
    n_pass = sum(int(r["passes"]) for r in rows)
    n_total = len(rows)
    assert n_pass == n_total, f"M04 tensor round-trip: {n_pass}/{n_total} pass"
    print(f"M04 tensor round-trip: {n_pass}/{n_total} pass (no waivers)")
    rows = _projector_identity_test()
    assert all(r["passes"] for r in rows), "M04 projector identity"
    print(f"M04 projector identity: {len(rows)} cases all pass")
    print("M04 tensor transforms: all checks passed")
