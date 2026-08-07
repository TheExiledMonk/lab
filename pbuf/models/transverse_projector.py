"""M07 — Transverse Projector.

Build the per-voxel transverse projector P_T = I - ê_L ê_L^T from a
longitudinal unit-vector field.
"""
from __future__ import annotations
import numpy as np

from ..core.conventions import EPS_FLOAT, EPS_ZERO

__all__ = [
    "build_longitudinal_direction", "build_transverse_projector",
    "validate_transverse_projector", "project_pair_direction",
    "TransverseProjectorError",
]


class TransverseProjectorError(ValueError):
    pass


def build_longitudinal_direction(scalar_field, spacing=(1.0, 1.0, 1.0)):
    """Compute the per-voxel longitudinal unit vector ê_L = ∇f / |∇f|.

    Parameters
    ----------
    scalar_field : ndarray of shape (nz, ny, nx)
    spacing : (dz, dy, dx)

    Returns (eL_x, eL_y, eL_z, valid_mask, grad_magnitude).
    """
    dz, dy, dx = spacing
    gz, gy, gx = np.gradient(scalar_field, dz, dy, dx, edge_order=1)
    g_mag = np.sqrt(gx ** 2 + gy ** 2 + gz ** 2)
    valid = g_mag > 1e-12
    safe = np.where(valid, g_mag, 1.0)
    eL_x = np.where(valid, gx / safe, 0.0)
    eL_y = np.where(valid, gy / safe, 0.0)
    eL_z = np.where(valid, gz / safe, 0.0)
    return eL_x, eL_y, eL_z, valid, g_mag


def build_transverse_projector(eL_x, eL_y, eL_z):
    """P_T = I - ê_L ê_L^T as a 3x3 tensor field.

    Returns (Pxx, Pxy, Pxz, Pyy, Pyz, Pzz).
    """
    return (
        1.0 - eL_x * eL_x,
        -eL_x * eL_y,
        -eL_x * eL_z,
        1.0 - eL_y * eL_y,
        -eL_y * eL_z,
        1.0 - eL_z * eL_z,
    )


def validate_transverse_projector(Pxx, Pxy, Pxz, Pyy, Pyz, Pzz,
                                    eL_x, eL_y, eL_z):
    """Check P_T^2 = P_T, P_T^T = P_T (symmetric by construction), and
    P_T ê_L = 0."""
    # Idempotence: P² = P.
    Pxx2 = Pxx * Pxx + Pxy * Pxy + Pxz * Pxz
    Pxy2 = Pxx * Pxy + Pxy * Pyy + Pxz * Pyz
    Pxz2 = Pxx * Pxz + Pxy * Pyz + Pxz * Pzz
    Pyy2 = Pxy * Pxy + Pyy * Pyy + Pyz * Pyz
    Pyz2 = Pxy * Pxz + Pyy * Pyz + Pyz * Pzz
    Pzz2 = Pxz * Pxz + Pyz * Pyz + Pzz * Pzz
    err_idem = float(max(
        np.max(np.abs(Pxx2 - Pxx)),
        np.max(np.abs(Pxy2 - Pxy)),
        np.max(np.abs(Pxz2 - Pxz)),
        np.max(np.abs(Pyy2 - Pyy)),
        np.max(np.abs(Pyz2 - Pyz)),
        np.max(np.abs(Pzz2 - Pzz)),
    ))
    # P ê_L = 0.
    PeL_x = Pxx * eL_x + Pxy * eL_y + Pxz * eL_z
    PeL_y = Pxy * eL_x + Pyy * eL_y + Pyz * eL_z
    PeL_z = Pxz * eL_x + Pyz * eL_y + Pzz * eL_z
    err_long = float(max(
        np.max(np.abs(PeL_x)),
        np.max(np.abs(PeL_y)),
        np.max(np.abs(PeL_z)),
    ))
    return {
        "err_idempotence": err_idem,
        "err_longitudinal": err_long,
        "passes_idempotence": err_idem < 1e-14,
        "passes_longitudinal": err_long < 1e-14,
        "passes": err_idem < 1e-14 and err_long < 1e-14,
    }


def project_pair_direction(projector, direction_xyz):
    """Compute P_T · n where n is a single direction (dx, dy, dz).

    Returns the projected vector (vx, vy, vz). The caller may apply
    PM1 (unit-magnitude) or PM2 (raw) normalisation.
    """
    Pxx, Pxy, Pxz, Pyy, Pyz, Pzz = projector
    dx, dy, dz = direction_xyz
    vx = Pxx * dx + Pxy * dy + Pxz * dz
    vy = Pxy * dx + Pyy * dy + Pyz * dz
    vz = Pxz * dx + Pyz * dy + Pzz * dz
    return vx, vy, vz


# ----------------------------------------------------------------------
# Self-check
# ----------------------------------------------------------------------
def _uniform_longitudinal_test():
    """For uniform ê_L = (1, 0, 0), P_T is well-defined and idempotent.
    Projecting (1, 0, 0) must give zero."""
    nz, ny, nx = 4, 5, 6
    eL_x = np.ones((nz, ny, nx))
    eL_y = np.zeros((nz, ny, nx))
    eL_z = np.zeros((nz, ny, nx))
    Pxx, Pxy, Pxz, Pyy, Pyz, Pzz = build_transverse_projector(
        eL_x, eL_y, eL_z)
    v = validate_transverse_projector(Pxx, Pxy, Pxz, Pyy, Pyz, Pzz,
                                        eL_x, eL_y, eL_z)
    vx, vy, vz = project_pair_direction(
        (Pxx, Pxy, Pxz, Pyy, Pyz, Pzz), (1, 0, 0))
    return {"validation": v, "projection_max": float(max(
        np.max(np.abs(vx)), np.max(np.abs(vy)), np.max(np.abs(vz)))),
        "passes": v["passes"] and float(max(np.max(np.abs(vx)),
                                              np.max(np.abs(vy)),
                                              np.max(np.abs(vz)))) < 1e-14}


def _varying_longitudinal_test():
    """For a varying ê_L, P_T is idempotent and ê_L is in its null-space."""
    nz, ny, nx = 5, 6, 7
    Z, Y, X = np.meshgrid(np.arange(nz), np.arange(ny), np.arange(nx),
                          indexing="ij")
    scalar = (X.astype(np.float64) + 2 * Y.astype(np.float64) +
              3 * Z.astype(np.float64))
    eL_x, eL_y, eL_z, valid, g_mag = build_longitudinal_direction(scalar)
    proj = build_transverse_projector(eL_x, eL_y, eL_z)
    v = validate_transverse_projector(*proj, eL_x, eL_y, eL_z)
    # Project the local longitudinal direction at every voxel.
    proj_field = project_pair_direction(proj, (eL_x, eL_y, eL_z))
    max_proj = float(max(np.max(np.abs(proj_field[0])),
                          np.max(np.abs(proj_field[1])),
                          np.max(np.abs(proj_field[2]))))
    n_valid = int(valid.sum())
    n_zero_grad = int((~valid).sum())
    # Project a perpendicular direction at every voxel.
    perp_x = -eL_y
    perp_y = eL_x
    perp_z = np.zeros_like(eL_x)
    perp_proj = project_pair_direction(proj, (perp_x, perp_y, perp_z))
    # Re-projected magnitude must equal the input magnitude on valid
    # voxels (P_T^2 = P_T).
    input_mag = np.sqrt(perp_x ** 2 + perp_y ** 2 + perp_z ** 2)
    output_mag = np.sqrt(perp_proj[0] ** 2 + perp_proj[1] ** 2 +
                          perp_proj[2] ** 2)
    diff = output_mag - input_mag
    max_idem = float(np.max(np.abs(diff[valid]))) if n_valid > 0 else 0.0
    return {
        "validation": v,
        "longitudinal_projection_max": max_proj,
        "perpendicular_idempotence_max": max_idem,
        "n_valid": n_valid,
        "n_zero_grad": n_zero_grad,
        "passes": v["passes"] and max_proj < 1e-13 and max_idem < 1e-13,
    }


if __name__ == "__main__":
    r = _uniform_longitudinal_test()
    assert r["passes"], f"uniform longitudinal failed: {r}"
    print(f"M07 uniform longitudinal: idempotence+null-space OK; "
          f"projection_max={r['projection_max']:.3e}")
    r = _varying_longitudinal_test()
    assert r["passes"], f"varying longitudinal failed: {r}"
    print(f"M07 varying longitudinal: {r['n_valid']} valid voxels, "
          f"idem_max={r['perpendicular_idempotence_max']:.3e}, "
          f"long_max={r['longitudinal_projection_max']:.3e}")
    print("M07 transverse projector: all checks passed")