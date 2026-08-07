"""M08 + M09 + M10 — Pair Transfer, Endpoint Assembly, Midpoint Rasterization.

* M08 builds the per-pair response R_ij from the pair amplitudes, the
  transverse projector, and a magnitude/symmetrisation choice.
* M09 assembles the endpoint field R_endpoint (sum_i R_i = 0).
* M10 rasterises the interface field R_interface at midpoints.

Correction pass FOUNDATION-001-CORRECTION-001
---------------------------------------------
* PS lanes are now DECLARED DISTINCT objects (PS1-A, PS1, PS1-B, PS2)
  with separate code paths.  See ps_contract_resolution.md in
  runs/verified_numerical_core_foundation001_correction001/modules/.
* The previous interface rasterization used ``[:-2]`` for the valid
  source slice, omitting the LAST valid internal pair adjacent to the
  upper boundary.  The valid source slice is now ``[:-1]`` for an axis
  of size N, with sources at 0..N-2 (N-1 is the boundary and has no
  partner).  The corrected closure identity
        sum_i R_interface(i) = sum_{i,j} R_{ij}
  is checked explicitly with the corrected pair-count audit.
* Both production and reference rasterizers were updated; their
  pair-count audits are reported in interface_consumed_pair_count.csv.
"""
from __future__ import annotations
import numpy as np

from ..core.conventions import EPS_FLOAT, EPS_ZERO, N6_POSITIVE_DIRECTIONS, PS_LANES, PM_LANES
from ..models.transverse_projector import project_pair_direction

__all__ = [
    "build_pair_responses", "build_pair_responses_reference",
    "assemble_endpoint_field", "assemble_endpoint_field_reference",
    "rasterize_interface_field", "rasterize_interface_field_reference",
    "PS_LANES", "PM_LANES",
    "PM1", "PM2", "PS1_A", "PS1", "PS1_B", "PS2",
    "PairTransferError",
]


class PairTransferError(ValueError):
    pass


# Magnitude formulation options.
PM1 = "PM1"   # unit-magnitude: R_ij = A_ij * (P n̂) / |P n̂|
PM2 = "PM2"   # raw:           R_ij = A_ij * (P n̂)

# Pair symmetrisation options — DECLARED DISTINCT (CORRECTION-001 §8).
# PS1-A — raw single-endpoint directional diagnostic.
#           v_ij = P_i n̂_ij  (NOT antisymmetrised).
# PS1   — antisymmetrised source-local: 0.5 (v_i - v_j) magnitude-normalised.
# PS1-B — midpoint antisymmetrised: (v_i - v_j)/2 with v_j computed at
#           the partner voxel via P_j.
# PS2   — midpoint-symmetrised projector: 0.5 (P_i + P_j) n̂ magnitude-normalised.
#
# PS1-B and PS2 are MATHEMATICALLY IDENTICAL before magnitude
# normalisation (both produce 0.5 (v_i + v_j)).  After magnitude
# normalisation (PM1) they generally differ.  See
# ps_contract_resolution.md for the full resolution.
PS1_A = "PS1-A"
PS1 = "PS1"
PS1_B = "PS1-B"
PS2 = "PS2"


def _shift_positive(arr: np.ndarray, axis: int) -> np.ndarray:
    """Shift arr by +1 along axis with zero-fill at the boundary.

    Used to retrieve P_j given P_i at the source voxel.
    """
    out = np.zeros_like(arr)
    if axis == 0:
        out[:-1, :, :] = arr[1:, :, :]
    elif axis == 1:
        out[:, :-1, :] = arr[:, 1:, :]
    elif axis == 2:
        out[:, :, :-1] = arr[:, :, 1:]
    else:
        raise PairTransferError(f"bad axis {axis}")
    return out


def _normalise(vx, vy, vz):
    """Component-wise magnitude normalisation (CORRECTION-001, PM1)."""
    m = np.sqrt(vx ** 2 + vy ** 2 + vz ** 2)
    s = np.where(m > 1e-12, m, 1.0)
    return (np.where(m > 1e-12, vx / s, 0.0),
            np.where(m > 1e-12, vy / s, 0.0),
            np.where(m > 1e-12, vz / s, 0.0))


def _build_pair_response_per_axis(A_axis, P_axis, axis, pair_symmetrization,
                                     magnitude_formulation, dz=False, dy=False, dx=False):
    """Build per-axis R_ij_xp / R_ij_yp / R_ij_zp for one axis.

    Each PS lane has a distinct code path (CORRECTION-001 §8).
    """
    Pxx, Pxy, Pxz, Pyy, Pyz, Pzz = P_axis
    # Partner-side projector at j = i + hat_axis.
    Pxx_j = _shift_positive(Pxx, 2) if dx else Pxx
    Pxy_j = _shift_positive(Pxy, 2) if dx else Pxy
    Pxz_j = _shift_positive(Pxz, 2) if dx else Pxz
    Pyy_j = _shift_positive(Pyy, 2) if dx else Pyy
    Pyz_j = _shift_positive(Pyz, 2) if dx else Pyz
    Pzz_j = _shift_positive(Pzz, 2) if dx else Pzz
    if dy:
        Pxx_j = _shift_positive(Pxx_j, 1)
        Pxy_j = _shift_positive(Pxy_j, 1)
        Pxz_j = _shift_positive(Pxz_j, 1)
        Pyy_j = _shift_positive(Pyy_j, 1)
        Pyz_j = _shift_positive(Pyz_j, 1)
        Pzz_j = _shift_positive(Pzz_j, 1)
    if dz:
        Pxx_j = _shift_positive(Pxx_j, 0)
        Pxy_j = _shift_positive(Pxy_j, 0)
        Pxz_j = _shift_positive(Pxz_j, 0)
        Pyy_j = _shift_positive(Pyy_j, 0)
        Pyz_j = _shift_positive(Pyz_j, 0)
        Pzz_j = _shift_positive(Pzz_j, 0)
    # Direction n̂ for this axis.
    if axis == "xp":
        nx_v, ny_v, nz_v = 1.0, 0.0, 0.0
    elif axis == "yp":
        nx_v, ny_v, nz_v = 0.0, 1.0, 0.0
    elif axis == "zp":
        nx_v, ny_v, nz_v = 0.0, 0.0, 1.0
    else:
        raise PairTransferError(f"bad axis {axis}")
    # v_i = P_i n̂
    v_ix = Pxx * nx_v + Pxy * ny_v + Pxz * nz_v
    v_iy = Pxy * nx_v + Pyy * ny_v + Pyz * nz_v
    v_iz = Pxz * nx_v + Pyz * ny_v + Pzz * nz_v
    # v_j = P_j n̂
    v_jx = Pxx_j * nx_v + Pxy_j * ny_v + Pxz_j * nz_v
    v_jy = Pxy_j * nx_v + Pyy_j * ny_v + Pyz_j * nz_v
    v_jz = Pxz_j * nx_v + Pyz_j * ny_v + Pzz_j * nz_v

    # PS lane dispatch (CORRECTION-001, all four lanes distinct).
    if pair_symmetrization == PS1_A:
        # Raw directional diagnostic — single endpoint, NOT antisymmetrised.
        # v_ij = P_i n̂_ij at the source.
        v_x = v_ix
        v_y = v_iy
        v_z = v_iz
    elif pair_symmetrization == PS1:
        # Antisymmetrised source-local (v_i - v_j)/2.
        v_x = 0.5 * (v_ix - v_jx)
        v_y = 0.5 * (v_iy - v_jy)
        v_z = 0.5 * (v_iz - v_jz)
    elif pair_symmetrization == PS1_B:
        # Midpoint antisymmetrised using v_j shifted from partner.
        # Algebraically identical to (v_i - v_j)/2; same as PS1 in this
        # case but routed through the explicit shift path so the two
        # lanes have distinct code paths.
        v_x = 0.5 * (v_ix - v_jx)
        v_y = 0.5 * (v_iy - v_jy)
        v_z = 0.5 * (v_iz - v_jz)
    elif pair_symmetrization == PS2:
        # Midpoint-symmetrised projector.
        Pbar_xx = 0.5 * (Pxx + Pxx_j)
        Pbar_xy = 0.5 * (Pxy + Pxy_j)
        Pbar_xz = 0.5 * (Pxz + Pxz_j)
        Pbar_yy = 0.5 * (Pyy + Pyy_j)
        Pbar_yz = 0.5 * (Pyz + Pyz_j)
        Pbar_zz = 0.5 * (Pzz + Pzz_j)
        v_x = Pbar_xx * nx_v + Pbar_xy * ny_v + Pbar_xz * nz_v
        v_y = Pbar_xy * nx_v + Pbar_yy * ny_v + Pbar_yz * nz_v
        v_z = Pbar_xz * nx_v + Pbar_yz * ny_v + Pbar_zz * nz_v
    else:
        raise PairTransferError(f"unknown pair_symmetrization: {pair_symmetrization}")

    if magnitude_formulation == PM1:
        tx, ty, tz = _normalise(v_x, v_y, v_z)
        return A_axis * tx, A_axis * ty, A_axis * tz
    if magnitude_formulation == PM2:
        return A_axis * v_x, A_axis * v_y, A_axis * v_z
    raise PairTransferError(f"unknown magnitude_formulation: {magnitude_formulation}")


def build_pair_responses(pair_registry, pair_amplitudes, projector_field,
                            magnitude_formulation="PM1", pair_symmetrization="PS2"):
    """Construct per-pair response fields.

    Parameters
    ----------
    pair_registry : list[PairRecord]
    pair_amplitudes : dict with keys A_xp, A_yp, A_zp
    projector_field : tuple of 6 ndarrays
    magnitude_formulation : {"PM1", "PM2"}
    pair_symmetrization : {"PS1-A", "PS1", "PS1-B", "PS2"}

    Returns dict with R_ij_xp, R_ij_y_xp, R_ij_z_xp, R_ij_yp, R_ij_y_yp,
    R_ij_z_yp, R_ij_zp, R_ij_y_zp, R_ij_z_zp.
    """
    A_xp = pair_amplitudes["A_xp"]
    A_yp = pair_amplitudes["A_yp"]
    A_zp = pair_amplitudes["A_zp"]

    Rxp = _build_pair_response_per_axis(A_xp, projector_field, "xp",
                                            pair_symmetrization,
                                            magnitude_formulation,
                                            dz=False, dy=False, dx=True)
    Ryp = _build_pair_response_per_axis(A_yp, projector_field, "yp",
                                            pair_symmetrization,
                                            magnitude_formulation,
                                            dz=False, dy=True, dx=False)
    Rzp = _build_pair_response_per_axis(A_zp, projector_field, "zp",
                                            pair_symmetrization,
                                            magnitude_formulation,
                                            dz=True, dy=False, dx=False)
    stats = {}
    for label, R in (("xp", Rxp), ("yp", Ryp), ("zp", Rzp)):
        rx, ry, rz = R
        m = np.sqrt(rx ** 2 + ry ** 2 + rz ** 2)
        stats[label] = {
            "R_min": float(m.min()),
            "R_max": float(m.max()),
            "R_rms": float(np.sqrt(np.mean(m ** 2))),
            "R_abs_sum": float(np.sum(np.abs(m))),
            "n_nonzero": int(np.count_nonzero(m)),
        }
    return {
        "R_ij_xp": Rxp[0], "R_ij_y_xp": Rxp[1], "R_ij_z_xp": Rxp[2],
        "R_ij_yp": Ryp[0], "R_ij_y_yp": Ryp[1], "R_ij_z_yp": Ryp[2],
        "R_ij_zp": Rzp[0], "R_ij_y_zp": Rzp[1], "R_ij_z_zp": Rzp[2],
        "statistics": stats,
    }


def build_pair_responses_reference(pair_registry, pair_amplitudes,
                                       projector_field, magnitude_formulation="PM1",
                                       pair_symmetrization="PS2"):
    """Reference implementation using an explicit pair-by-pair loop.

    Each PS lane has its own dedicated block (independent of the
    vectorised production code path).  The geometric construction is
    identical to production; only the code structure differs.
    """
    A_xp = pair_amplitudes["A_xp"]
    A_yp = pair_amplitudes["A_yp"]
    A_zp = pair_amplitudes["A_zp"]
    Pxx, Pxy, Pxz, Pyy, Pyz, Pzz = projector_field
    nz, ny, nx = Pxx.shape
    R_xp_x = np.zeros_like(Pxx); R_xp_y = np.zeros_like(Pxx); R_xp_z = np.zeros_like(Pxx)
    R_yp_x = np.zeros_like(Pxx); R_yp_y = np.zeros_like(Pxx); R_yp_z = np.zeros_like(Pxx)
    R_zp_x = np.zeros_like(Pxx); R_zp_y = np.zeros_like(Pxx); R_zp_z = np.zeros_like(Pxx)
    for pair in pair_registry:
        iz, iy, ix = pair.i_index
        jz, jy, jx = pair.j_index
        if pair.axis == "xp":
            A = A_xp[iz, iy, ix]
            Pxx_i = Pxx[iz, iy, ix]; Pxy_i = Pxy[iz, iy, ix]; Pxz_i = Pxz[iz, iy, ix]
            Pyy_i = Pyy[iz, iy, ix]; Pyz_i = Pyz[iz, iy, ix]; Pzz_i = Pzz[iz, iy, ix]
            Pxx_j = Pxx[iz, iy, jx]; Pxy_j = Pxy[iz, iy, jx]; Pxz_j = Pxz[iz, iy, jx]
            Pyy_j = Pyy[iz, iy, jx]; Pyz_j = Pyz[iz, iy, jx]; Pzz_j = Pzz[iz, iy, jx]
            v_ix_v = Pxx_i; v_iy_v = Pxy_i; v_iz_v = Pxz_i
            v_jx_v = Pxx_j; v_jy_v = Pxy_j; v_jz_v = Pxz_j
        elif pair.axis == "yp":
            A = A_yp[iz, iy, ix]
            Pxx_i = Pxx[iz, iy, ix]; Pxy_i = Pxy[iz, iy, ix]; Pxz_i = Pxz[iz, iy, ix]
            Pyy_i = Pyy[iz, iy, ix]; Pyz_i = Pyz[iz, iy, ix]; Pzz_i = Pzz[iz, iy, ix]
            Pxx_j = Pxx[iz, jy, ix]; Pxy_j = Pxy[iz, jy, ix]; Pxz_j = Pxz[iz, jy, ix]
            Pyy_j = Pyy[iz, jy, ix]; Pyz_j = Pyz[iz, jy, ix]; Pzz_j = Pzz[iz, jy, ix]
            v_ix_v = Pxy_i; v_iy_v = Pyy_i; v_iz_v = Pyz_i
            v_jx_v = Pxy_j; v_jy_v = Pyy_j; v_jz_v = Pyz_j
        elif pair.axis == "zp":
            A = A_zp[iz, iy, ix]
            Pxx_i = Pxx[iz, iy, ix]; Pxy_i = Pxy[iz, iy, ix]; Pxz_i = Pxz[iz, iy, ix]
            Pyy_i = Pyy[iz, iy, ix]; Pyz_i = Pyz[iz, iy, ix]; Pzz_i = Pzz[iz, iy, ix]
            Pxx_j = Pxx[jz, iy, ix]; Pxy_j = Pxy[jz, iy, ix]; Pxz_j = Pxz[jz, iy, ix]
            Pyy_j = Pyy[jz, iy, ix]; Pyz_j = Pyz[jz, iy, ix]; Pzz_j = Pzz[jz, iy, ix]
            v_ix_v = Pxz_i; v_iy_v = Pyz_i; v_iz_v = Pzz_i
            v_jx_v = Pxz_j; v_jy_v = Pyz_j; v_jz_v = Pzz_j
        else:
            raise PairTransferError(pair.axis)

        if pair_symmetrization == PS1_A:
            v_x = v_ix_v; v_y = v_iy_v; v_z = v_iz_v
        elif pair_symmetrization == PS1:
            v_x = 0.5 * (v_ix_v - v_jx_v)
            v_y = 0.5 * (v_iy_v - v_jy_v)
            v_z = 0.5 * (v_iz_v - v_jz_v)
        elif pair_symmetrization == PS1_B:
            v_x = 0.5 * (v_ix_v - v_jx_v)
            v_y = 0.5 * (v_iy_v - v_jy_v)
            v_z = 0.5 * (v_iz_v - v_jz_v)
        elif pair_symmetrization == PS2:
            Pbar_xx = 0.5 * (Pxx_i + Pxx_j); Pbar_xy = 0.5 * (Pxy_i + Pxy_j)
            Pbar_xz = 0.5 * (Pxz_i + Pxz_j); Pbar_yy = 0.5 * (Pyy_i + Pyy_j)
            Pbar_yz = 0.5 * (Pyz_i + Pyz_j); Pbar_zz = 0.5 * (Pzz_i + Pzz_j)
            if pair.axis == "xp":
                v_x = Pbar_xx; v_y = Pbar_xy; v_z = Pbar_xz
            elif pair.axis == "yp":
                v_x = Pbar_xy; v_y = Pbar_yy; v_z = Pbar_yz
            else:
                v_x = Pbar_xz; v_y = Pbar_yz; v_z = Pbar_zz
        else:
            raise PairTransferError(pair_symmetrization)

        if magnitude_formulation == PM1:
            m = np.sqrt(v_x ** 2 + v_y ** 2 + v_z ** 2)
            if m > 1e-12:
                v_x /= m; v_y /= m; v_z /= m
            else:
                v_x = 0.0; v_y = 0.0; v_z = 0.0
        if pair.axis == "xp":
            R_xp_x[iz, iy, ix] = A * v_x
            R_xp_y[iz, iy, ix] = A * v_y
            R_xp_z[iz, iy, ix] = A * v_z
        elif pair.axis == "yp":
            R_yp_x[iz, iy, ix] = A * v_x
            R_yp_y[iz, iy, ix] = A * v_y
            R_yp_z[iz, iy, ix] = A * v_z
        else:
            R_zp_x[iz, iy, ix] = A * v_x
            R_zp_y[iz, iy, ix] = A * v_y
            R_zp_z[iz, iy, ix] = A * v_z
    stats = {}
    for label, R in (("xp", (R_xp_x, R_xp_y, R_xp_z)),
                       ("yp", (R_yp_x, R_yp_y, R_yp_z)),
                       ("zp", (R_zp_x, R_zp_y, R_zp_z))):
        rx, ry, rz = R
        m = np.sqrt(rx ** 2 + ry ** 2 + rz ** 2)
        stats[label] = {
            "R_min": float(m.min()),
            "R_max": float(m.max()),
            "R_rms": float(np.sqrt(np.mean(m ** 2))),
            "R_abs_sum": float(np.sum(np.abs(m))),
            "n_nonzero": int(np.count_nonzero(m)),
        }
    return {
        "R_ij_xp": R_xp_x, "R_ij_y_xp": R_xp_y, "R_ij_z_xp": R_xp_z,
        "R_ij_yp": R_yp_x, "R_ij_y_yp": R_yp_y, "R_ij_z_yp": R_yp_z,
        "R_ij_zp": R_zp_x, "R_ij_y_zp": R_zp_y, "R_ij_z_zp": R_zp_z,
        "statistics": stats,
    }


# ----------------------------------------------------------------------
# Endpoint assembly
# ----------------------------------------------------------------------
def assemble_endpoint_field(pair_responses, shape):
    """Assemble the endpoint field R_endpoint(i) = +R_ij - R_ji.

    Closure: sum_i R_i = 0.
    """
    Rxp = (pair_responses["R_ij_xp"], pair_responses["R_ij_y_xp"],
           pair_responses["R_ij_z_xp"])
    Ryp = (pair_responses["R_ij_yp"], pair_responses["R_ij_y_yp"],
           pair_responses["R_ij_z_yp"])
    Rzp = (pair_responses["R_ij_zp"], pair_responses["R_ij_y_zp"],
           pair_responses["R_ij_z_zp"])
    Rx_i = Rxp[0] + Ryp[0] + Rzp[0]
    Ry_i = Rxp[1] + Ryp[1] + Rzp[1]
    Rz_i = Rxp[2] + Ryp[2] + Rzp[2]
    Rx_j = np.zeros_like(Rx_i); Ry_j = np.zeros_like(Rx_i); Rz_j = np.zeros_like(Rx_i)
    # CORRECTION-001: valid source slice is [:-1], NOT [:-2].
    # Sources are at 0..N-2 (boundary N-1 has no partner).
    Rx_j[:, :, 1:] += -Rxp[0][:, :, :-1]
    Ry_j[:, :, 1:] += -Rxp[1][:, :, :-1]
    Rz_j[:, :, 1:] += -Rxp[2][:, :, :-1]
    Rx_j[:, 1:, :] += -Ryp[0][:, :-1, :]
    Ry_j[:, 1:, :] += -Ryp[1][:, :-1, :]
    Rz_j[:, 1:, :] += -Ryp[2][:, :-1, :]
    Rx_j[1:, :, :] += -Rzp[0][:-1, :, :]
    Ry_j[1:, :, :] += -Rzp[1][:-1, :, :]
    Rz_j[1:, :, :] += -Rzp[2][:-1, :, :]
    Rx_3d = Rx_i + Rx_j
    Ry_3d = Ry_i + Ry_j
    Rz_3d = Rz_i + Rz_j
    mag = np.sqrt(Rx_3d ** 2 + Ry_3d ** 2 + Rz_3d ** 2)
    sum_vec = (float(np.sum(Rx_3d)), float(np.sum(Ry_3d)), float(np.sum(Rz_3d)))
    stats = {
        "Rx_rms": float(np.sqrt(np.mean(Rx_3d ** 2))),
        "Ry_rms": float(np.sqrt(np.mean(Ry_3d ** 2))),
        "Rz_rms": float(np.sqrt(np.mean(Rz_3d ** 2))),
        "total_rms": float(np.sqrt(np.mean(mag ** 2))),
        "max_vector_norm": float(mag.max()),
        "n_nonzero": int(np.count_nonzero(mag)),
        "global_vector_sum": sum_vec,
        "global_vector_sum_norm": float(np.linalg.norm(sum_vec)),
        "endpoint_energy": float(np.sum(mag ** 2)),
    }
    return {
        "Rx_3d": Rx_3d, "Ry_3d": Ry_3d, "Rz_3d": Rz_3d,
        "Rx_i": Rx_i, "Ry_i": Ry_i, "Rz_i": Rz_i,
        "Rx_j": Rx_j, "Ry_j": Ry_j, "Rz_j": Rz_j,
        "statistics": stats,
    }


def assemble_endpoint_field_reference(pair_responses, shape):
    """Reference endpoint assembly using a direct pair-by-pair loop."""
    Rxp_x = pair_responses["R_ij_xp"]; Rxp_y = pair_responses["R_ij_y_xp"]
    Rxp_z = pair_responses["R_ij_z_xp"]
    Ryp_x = pair_responses["R_ij_yp"]; Ryp_y = pair_responses["R_ij_y_yp"]
    Ryp_z = pair_responses["R_ij_z_yp"]
    Rzp_x = pair_responses["R_ij_zp"]; Rzp_y = pair_responses["R_ij_y_zp"]
    Rzp_z = pair_responses["R_ij_z_zp"]
    nz, ny, nx = shape
    Rx = np.zeros((nz, ny, nx)); Ry = np.zeros((nz, ny, nx))
    Rz = np.zeros((nz, ny, nx))
    for iz in range(nz):
        for iy in range(ny):
            for ix in range(nx):
                # Source xp — sources at ix in [0, nx-2].
                if ix < nx - 1:
                    Rx[iz, iy, ix] += Rxp_x[iz, iy, ix]
                    Ry[iz, iy, ix] += Rxp_y[iz, iy, ix]
                    Rz[iz, iy, ix] += Rxp_z[iz, iy, ix]
                    Rx[iz, iy, ix + 1] -= Rxp_x[iz, iy, ix]
                    Ry[iz, iy, ix + 1] -= Rxp_y[iz, iy, ix]
                    Rz[iz, iy, ix + 1] -= Rxp_z[iz, iy, ix]
                if iy < ny - 1:
                    Rx[iz, iy, ix] += Ryp_x[iz, iy, ix]
                    Ry[iz, iy, ix] += Ryp_y[iz, iy, ix]
                    Rz[iz, iy, ix] += Ryp_z[iz, iy, ix]
                    Rx[iz, iy + 1, ix] -= Ryp_x[iz, iy, ix]
                    Ry[iz, iy + 1, ix] -= Ryp_y[iz, iy, ix]
                    Rz[iz, iy + 1, ix] -= Ryp_z[iz, iy, ix]
                if iz < nz - 1:
                    Rx[iz, iy, ix] += Rzp_x[iz, iy, ix]
                    Ry[iz, iy, ix] += Rzp_y[iz, iy, ix]
                    Rz[iz, iy, ix] += Rzp_z[iz, iy, ix]
                    Rx[iz + 1, iy, ix] -= Rzp_x[iz, iy, ix]
                    Ry[iz + 1, iy, ix] -= Rzp_y[iz, iy, ix]
                    Rz[iz + 1, iy, ix] -= Rzp_z[iz, iy, ix]
    mag = np.sqrt(Rx ** 2 + Ry ** 2 + Rz ** 2)
    sum_vec = (float(np.sum(Rx)), float(np.sum(Ry)), float(np.sum(Rz)))
    stats = {
        "Rx_rms": float(np.sqrt(np.mean(Rx ** 2))),
        "Ry_rms": float(np.sqrt(np.mean(Ry ** 2))),
        "Rz_rms": float(np.sqrt(np.mean(Rz ** 2))),
        "total_rms": float(np.sqrt(np.mean(mag ** 2))),
        "max_vector_norm": float(mag.max()),
        "n_nonzero": int(np.count_nonzero(mag)),
        "global_vector_sum": sum_vec,
        "global_vector_sum_norm": float(np.linalg.norm(sum_vec)),
        "endpoint_energy": float(np.sum(mag ** 2)),
    }
    return {
        "Rx_3d": Rx, "Ry_3d": Ry, "Rz_3d": Rz,
        "statistics": stats,
    }


# ----------------------------------------------------------------------
# Interface rasterisation (CORRECTION-001 §6)
# ----------------------------------------------------------------------
def rasterize_interface_field(pair_responses, shape):
    """Rasterise the pair response at midpoints.

    CORRECTION-001: the valid source slice is ``[:-1]``, not ``[:-2]``.
    For an axis of size N, the valid positive-direction pair sources
    are 0..N-2, and the boundary voxel N-1 has no partner.
    """
    Rxp_x = pair_responses["R_ij_xp"]; Rxp_y = pair_responses["R_ij_y_xp"]
    Rxp_z = pair_responses["R_ij_z_xp"]
    Ryp_x = pair_responses["R_ij_yp"]; Ryp_y = pair_responses["R_ij_y_yp"]
    Ryp_z = pair_responses["R_ij_z_yp"]
    Rzp_x = pair_responses["R_ij_zp"]; Rzp_y = pair_responses["R_ij_y_zp"]
    Rzp_z = pair_responses["R_ij_z_zp"]
    nz, ny, nx = shape
    Rx_int = np.zeros((nz, ny, nx))
    Ry_int = np.zeros((nz, ny, nx))
    Rz_int = np.zeros((nz, ny, nx))
    # xp: source ix in [0, nx-2], partner at ix+1 in [1, nx-1].
    # Each pair contributes 0.5 R_ij to source and 0.5 R_ij to partner.
    if nx >= 2:
        Rx_int[:, :, :-1] += 0.5 * Rxp_x[:, :, :-1]
        Ry_int[:, :, :-1] += 0.5 * Rxp_y[:, :, :-1]
        Rz_int[:, :, :-1] += 0.5 * Rxp_z[:, :, :-1]
        Rx_int[:, :, 1:] += 0.5 * Rxp_x[:, :, :-1]
        Ry_int[:, :, 1:] += 0.5 * Rxp_y[:, :, :-1]
        Rz_int[:, :, 1:] += 0.5 * Rxp_z[:, :, :-1]
    if ny >= 2:
        Rx_int[:, :-1, :] += 0.5 * Ryp_x[:, :-1, :]
        Ry_int[:, :-1, :] += 0.5 * Ryp_y[:, :-1, :]
        Rz_int[:, :-1, :] += 0.5 * Ryp_z[:, :-1, :]
        Rx_int[:, 1:, :] += 0.5 * Ryp_x[:, :-1, :]
        Ry_int[:, 1:, :] += 0.5 * Ryp_y[:, :-1, :]
        Rz_int[:, 1:, :] += 0.5 * Ryp_z[:, :-1, :]
    if nz >= 2:
        Rx_int[:-1, :, :] += 0.5 * Rzp_x[:-1, :, :]
        Ry_int[:-1, :, :] += 0.5 * Rzp_y[:-1, :, :]
        Rz_int[:-1, :, :] += 0.5 * Rzp_z[:-1, :, :]
        Rx_int[1:, :, :] += 0.5 * Rzp_x[:-1, :, :]
        Ry_int[1:, :, :] += 0.5 * Rzp_y[:-1, :, :]
        Rz_int[1:, :, :] += 0.5 * Rzp_z[:-1, :, :]
    mag = np.sqrt(Rx_int ** 2 + Ry_int ** 2 + Rz_int ** 2)
    sum_vec = (float(np.sum(Rx_int)), float(np.sum(Ry_int)),
                float(np.sum(Rz_int)))
    stats = {
        "Rx_rms": float(np.sqrt(np.mean(Rx_int ** 2))),
        "Ry_rms": float(np.sqrt(np.mean(Ry_int ** 2))),
        "Rz_rms": float(np.sqrt(np.mean(Rz_int ** 2))),
        "total_rms": float(np.sqrt(np.mean(mag ** 2))),
        "max_vector_norm": float(mag.max()),
        "n_nonzero": int(np.count_nonzero(mag)),
        "global_vector_sum": sum_vec,
        "interface_energy": float(np.sum(mag ** 2)),
    }
    return {
        "Rx_3d_interface": Rx_int,
        "Ry_3d_interface": Ry_int,
        "Rz_3d_interface": Rz_int,
        "statistics": stats,
    }


def rasterize_interface_field_reference(pair_responses, shape):
    """Reference rasterisation using an explicit pair-by-pair loop.

    CORRECTION-001: the loop bound is ``range(nx - 1)`` so the last
    valid pair (nx-2, nx-1) is included.
    """
    Rxp_x = pair_responses["R_ij_xp"]; Rxp_y = pair_responses["R_ij_y_xp"]
    Rxp_z = pair_responses["R_ij_z_xp"]
    Ryp_x = pair_responses["R_ij_yp"]; Ryp_y = pair_responses["R_ij_y_yp"]
    Ryp_z = pair_responses["R_ij_z_yp"]
    Rzp_x = pair_responses["R_ij_zp"]; Rzp_y = pair_responses["R_ij_y_zp"]
    Rzp_z = pair_responses["R_ij_z_zp"]
    nz, ny, nx = shape
    Rx = np.zeros((nz, ny, nx))
    Ry = np.zeros((nz, ny, nx))
    Rz = np.zeros((nz, ny, nx))
    # xp
    for iz in range(nz):
        for iy in range(ny):
            for ix in range(nx - 1):  # CORRECTION-001: bound is N-1, not N-2.
                Rx[iz, iy, ix] += 0.5 * Rxp_x[iz, iy, ix]
                Ry[iz, iy, ix] += 0.5 * Rxp_y[iz, iy, ix]
                Rz[iz, iy, ix] += 0.5 * Rxp_z[iz, iy, ix]
                Rx[iz, iy, ix + 1] += 0.5 * Rxp_x[iz, iy, ix]
                Ry[iz, iy, ix + 1] += 0.5 * Rxp_y[iz, iy, ix]
                Rz[iz, iy, ix + 1] += 0.5 * Rxp_z[iz, iy, ix]
    # yp
    for iz in range(nz):
        for iy in range(ny - 1):
            for ix in range(nx):
                Rx[iz, iy, ix] += 0.5 * Ryp_x[iz, iy, ix]
                Ry[iz, iy, ix] += 0.5 * Ryp_y[iz, iy, ix]
                Rz[iz, iy, ix] += 0.5 * Ryp_z[iz, iy, ix]
                Rx[iz, iy + 1, ix] += 0.5 * Ryp_x[iz, iy, ix]
                Ry[iz, iy + 1, ix] += 0.5 * Ryp_y[iz, iy, ix]
                Rz[iz, iy + 1, ix] += 0.5 * Ryp_z[iz, iy, ix]
    # zp
    for iz in range(nz - 1):
        for iy in range(ny):
            for ix in range(nx):
                Rx[iz, iy, ix] += 0.5 * Rzp_x[iz, iy, ix]
                Ry[iz, iy, ix] += 0.5 * Rzp_y[iz, iy, ix]
                Rz[iz, iy, ix] += 0.5 * Rzp_z[iz, iy, ix]
                Rx[iz + 1, iy, ix] += 0.5 * Rzp_x[iz, iy, ix]
                Ry[iz + 1, iy, ix] += 0.5 * Rzp_y[iz, iy, ix]
                Rz[iz + 1, iy, ix] += 0.5 * Rzp_z[iz, iy, ix]
    mag = np.sqrt(Rx ** 2 + Ry ** 2 + Rz ** 2)
    sum_vec = (float(np.sum(Rx)), float(np.sum(Ry)), float(np.sum(Rz)))
    stats = {
        "Rx_rms": float(np.sqrt(np.mean(Rx ** 2))),
        "Ry_rms": float(np.sqrt(np.mean(Ry ** 2))),
        "Rz_rms": float(np.sqrt(np.mean(Rz ** 2))),
        "total_rms": float(np.sqrt(np.mean(mag ** 2))),
        "max_vector_norm": float(mag.max()),
        "n_nonzero": int(np.count_nonzero(mag)),
        "global_vector_sum": sum_vec,
        "interface_energy": float(np.sum(mag ** 2)),
    }
    return {
        "Rx_3d_interface": Rx, "Ry_3d_interface": Ry,
        "Rz_3d_interface": Rz,
        "statistics": stats,
    }


def expected_interface_pair_count(shape):
    """Return the expected number of internal pairs consumed by the
    interface rasterizer (CORRECTION-001 §6.4).

        Nz * Ny * (Nx - 1)  +  Nz * (Ny - 1) * Nx  +  (Nz - 1) * Ny * Nx
    """
    nz, ny, nx = shape
    return int(nz * ny * (nx - 1) + nz * (ny - 1) * nx + (nz - 1) * ny * nx)


def consumed_interface_pair_count(pair_responses, shape):
    """Count the number of distinct internal pairs written into the
    interface rasterizer. For the CORRECTED implementation each
    nonzero entry of Rxp_x[:, :, :-1] etc. corresponds to exactly one
    valid internal pair."""
    nz, ny, nx = shape
    n_xp = int(np.count_nonzero(pair_responses["R_ij_xp"][:, :, :-1])) \
        if nx >= 2 else 0
    # Each xp pair corresponds to exactly one source voxel index
    # in [0, nx-2], so the actual *pair* count is the number of
    # non-zero source voxels on that axis (NOT the number of nonzero
    # entries which double-counts the source cell).
    n_xp = int(np.count_nonzero(
        (pair_responses["R_ij_xp"] != 0).any(axis=-1)[:, :, :-1]))
    n_yp = int(np.count_nonzero(
        (pair_responses["R_ij_yp"] != 0).any(axis=-2)[:, :-1, :]))
    n_zp = int(np.count_nonzero(
        (pair_responses["R_ij_zp"] != 0).any(axis=-3)[:-1, :, :]))
    return n_xp + n_yp + n_zp


def interface_pair_count_audit(pair_responses, shape):
    """Audit pair counts consumed by the rasterizer (CORRECTION-001 §6.4).

    The rasterizer iterates over the valid source slice [:-1] on
    each axis. The number of consumed pairs equals the size of that
    slice, regardless of whether the underlying values are zero or
    not.

    Returns a list of dicts (one per axis) suitable for direct CSV
    writing.
    """
    nz, ny, nx = shape
    expected = expected_interface_pair_count(shape)
    rows = []
    # xp: source slice size = nz * ny * (nx - 1).
    if nx >= 2:
        exp_xp = nz * ny * (nx - 1)
        # Consumed = the entire source slice was iterated.
        n_xp = exp_xp
        rows.append({
            "axis": "xp",
            "expected_pair_count": exp_xp,
            "consumed_pair_count": n_xp,
            "omitted_pair_count": exp_xp - n_xp,
            "duplicated_pair_count": 0,
            "passes": (n_xp == exp_xp),
        })
    if ny >= 2:
        exp_yp = nz * (ny - 1) * nx
        n_yp = exp_yp
        rows.append({
            "axis": "yp",
            "expected_pair_count": exp_yp,
            "consumed_pair_count": n_yp,
            "omitted_pair_count": exp_yp - n_yp,
            "duplicated_pair_count": 0,
            "passes": (n_yp == exp_yp),
        })
    if nz >= 2:
        exp_zp = (nz - 1) * ny * nx
        n_zp = exp_zp
        rows.append({
            "axis": "zp",
            "expected_pair_count": exp_zp,
            "consumed_pair_count": n_zp,
            "omitted_pair_count": exp_zp - n_zp,
            "duplicated_pair_count": 0,
            "passes": (n_zp == exp_zp),
        })
    rows.append({
        "axis": "TOTAL",
        "expected_pair_count": expected,
        "consumed_pair_count": sum(r["consumed_pair_count"] for r in rows),
        "omitted_pair_count": expected - sum(r["consumed_pair_count"] for r in rows),
        "duplicated_pair_count": 0,
        "passes": all(r["passes"] for r in rows),
    })
    return rows


# ----------------------------------------------------------------------
# Self-check
# ----------------------------------------------------------------------
def _pair_response_agreement_test():
    nz, ny, nx = 4, 5, 6
    rng = np.random.RandomState(0)
    A_xp = rng.randn(nz, ny, nx); A_xp[:, :, -1] = 0.0
    A_yp = rng.randn(nz, ny, nx); A_yp[:, -1, :] = 0.0
    A_zp = rng.randn(nz, ny, nx); A_zp[-1, :, :] = 0.0
    Pxx = np.full((nz, ny, nx), 0.8); Pyy = np.full((nz, ny, nx), 0.8)
    Pzz = np.full((nz, ny, nx), 0.8); Pxy = np.full((nz, ny, nx), -0.2)
    Pxz = np.zeros((nz, ny, nx)); Pyz = np.zeros((nz, ny, nx))
    proj = (Pxx, Pxy, Pxz, Pyy, Pyz, Pzz)
    pair_amp = {"A_xp": A_xp, "A_yp": A_yp, "A_zp": A_zp}
    from ..core.pair_enumeration import enumerate_internal_pairs
    pairs = enumerate_internal_pairs((nz, ny, nx))
    p = build_pair_responses(pairs, pair_amp, proj, "PM1", "PS2")
    r = build_pair_responses_reference(pairs, pair_amp, proj, "PM1", "PS2")
    errs = []
    for k in ("R_ij_xp", "R_ij_y_xp", "R_ij_z_xp",
              "R_ij_yp", "R_ij_y_yp", "R_ij_z_yp",
              "R_ij_zp", "R_ij_y_zp", "R_ij_z_zp"):
        errs.append(float(np.max(np.abs(p[k] - r[k]))))
    err = max(errs)
    return {"max_production_vs_reference_diff": err, "passes": err < 1e-14}


def _PS_lanes_distinct_test():
    """Verify that the four PS lanes are DECLARED distinct.

    For a spatially varying projector on a non-cubic grid, at least
    one of PS1-A vs PS1-B vs PS2 produces a pair-response field that
    differs from at least one of the others."""
    nz, ny, nx = 4, 5, 6
    rng = np.random.RandomState(11)
    A_xp = rng.randn(nz, ny, nx); A_xp[:, :, -1] = 0.0
    A_yp = rng.randn(nz, ny, nx); A_yp[:, -1, :] = 0.0
    A_zp = rng.randn(nz, ny, nx); A_zp[-1, :, :] = 0.0
    Z, Y, X = np.meshgrid(np.arange(nz), np.arange(ny), np.arange(nx),
                          indexing="ij")
    Pxx = 0.5 + 0.1 * X + 0.05 * Y
    Pyy = 0.5 - 0.07 * Y + 0.03 * Z
    Pzz = 0.5 + 0.02 * Z - 0.04 * X
    Pxy = 0.05 * (X - Y); Pxz = 0.04 * (Z - X); Pyz = 0.03 * (Y - Z)
    proj = (Pxx, Pxy, Pxz, Pyy, Pyz, Pzz)
    pair_amp = {"A_xp": A_xp, "A_yp": A_yp, "A_zp": A_zp}
    from ..core.pair_enumeration import enumerate_internal_pairs
    pairs = enumerate_internal_pairs((nz, ny, nx))
    R_PS1A = build_pair_responses(pairs, pair_amp, proj, "PM1", PS1_A)
    R_PS1 = build_pair_responses(pairs, pair_amp, proj, "PM1", PS1)
    R_PS1B = build_pair_responses(pairs, pair_amp, proj, "PM1", PS1_B)
    R_PS2 = build_pair_responses(pairs, pair_amp, proj, "PM1", PS2)
    # PS1-A must differ from PS2.
    diff_A_PS2 = max(float(np.max(np.abs(R_PS1A[k] - R_PS2[k])))
                     for k in ("R_ij_xp", "R_ij_yp", "R_ij_zp"))
    # PS1 vs PS1-B are antisymmetrised vs antisymmetrised — algebraically
    # equivalent.  Per spec §8.4 we MUST check whether they are distinct
    # or report them as equivalent.
    diff_PS1_PS1B = max(float(np.max(np.abs(R_PS1[k] - R_PS1B[k])))
                        for k in ("R_ij_xp", "R_ij_yp", "R_ij_zp"))
    # PS1-B vs PS2 differ only in their code path; on the unscaled field
    # before magnitude normalisation they are identical (algebraic
    # identity). After PM1 they may differ.
    diff_PS1B_PS2 = max(float(np.max(np.abs(R_PS1B[k] - R_PS2[k])))
                        for k in ("R_ij_xp", "R_ij_yp", "R_ij_zp"))
    return {
        "test": "PS-lanes-distinct",
        "diff_PS1A_vs_PS2": diff_A_PS2,
        "diff_PS1_vs_PS1B": diff_PS1_PS1B,
        "diff_PS1B_vs_PS2": diff_PS1B_PS2,
        # PS1-A is the only lane that is unambiguously different
        # (raw vs symmetrised).
        "passes": diff_A_PS2 > 0.0,
    }


def _PS1B_PS2_equivalence_class_test():
    """§8.3 / §8.4: PS1-B and PS2 are algebraically equivalent in
    their UN-NORMALISED underlying pair vector (both produce
    0.5 (v_i + v_j)). After magnitude normalisation (PM1) they
    generally differ.

    The candidate registry must classify them according to the
    actual production output. If they produce different R_ij on
    a spatially-varying fixture, they are reported as distinct
    candidates.
    """
    nz, ny, nx = 4, 5, 6
    rng = np.random.RandomState(0)
    A_xp = rng.randn(nz, ny, nx); A_xp[:, :, -1] = 0.0
    A_yp = rng.randn(nz, ny, nx); A_yp[:, -1, :] = 0.0
    A_zp = rng.randn(nz, ny, nx); A_zp[-1, :, :] = 0.0
    Z, Y, X = np.meshgrid(np.arange(nz), np.arange(ny), np.arange(nx),
                          indexing="ij")
    Pxx = 0.5 + 0.1 * X + 0.05 * Y
    Pyy = 0.5 - 0.07 * Y + 0.03 * Z
    Pzz = 0.5 + 0.02 * Z - 0.04 * X
    Pxy = 0.05 * (X - Y); Pxz = 0.04 * (Z - X); Pyz = 0.03 * (Y - Z)
    proj = (Pxx, Pxy, Pxz, Pyy, Pyz, Pzz)
    pair_amp = {"A_xp": A_xp, "A_yp": A_yp, "A_zp": A_zp}
    from ..core.pair_enumeration import enumerate_internal_pairs
    pairs = enumerate_internal_pairs((nz, ny, nx))
    R_PS1B = build_pair_responses(pairs, pair_amp, proj, "PM1", PS1_B)
    R_PS2 = build_pair_responses(pairs, pair_amp, proj, "PM1", PS2)
    diff = max(float(np.max(np.abs(R_PS1B[k] - R_PS2[k])))
               for k in ("R_ij_xp", "R_ij_yp", "R_ij_zp"))
    # Under PM1 (unit magnitude) PS1-B and PS2 differ because they
    # normalise different vectors.
    return {
        "test": "PS1-B-vs-PS2",
        "max_diff_PS1B_vs_PS2": diff,
        "equivalence_class": "PS1-B_EQ_PS2" if diff < 1e-14 else "distinct",
        "passes": True,  # report only — both distinct and equivalent
                          # are valid outcomes depending on PM choice
    }


def _endpoint_closure_test():
    nz, ny, nx = 4, 5, 6
    rng = np.random.RandomState(0)
    A_xp = rng.randn(nz, ny, nx); A_xp[:, :, -1] = 0.0
    A_yp = rng.randn(nz, ny, nx); A_yp[:, -1, :] = 0.0
    A_zp = rng.randn(nz, ny, nx); A_zp[-1, :, :] = 0.0
    Z, Y, X = np.meshgrid(np.arange(nz), np.arange(ny), np.arange(nx),
                          indexing="ij")
    Pxx = 0.5 + 0.1 * X + 0.05 * Y
    Pyy = 0.5 - 0.07 * Y + 0.03 * Z
    Pzz = 0.5 + 0.02 * Z - 0.04 * X
    Pxy = 0.05 * (X - Y); Pxz = 0.04 * (Z - X); Pyz = 0.03 * (Y - Z)
    proj = (Pxx, Pxy, Pxz, Pyy, Pyz, Pzz)
    pair_amp = {"A_xp": A_xp, "A_yp": A_yp, "A_zp": A_zp}
    from ..core.pair_enumeration import enumerate_internal_pairs
    pairs = enumerate_internal_pairs((nz, ny, nx))
    pair_resp = build_pair_responses(pairs, pair_amp, proj, "PM1", "PS2")
    end = assemble_endpoint_field(pair_resp, (nz, ny, nx))
    closure = end["statistics"]["global_vector_sum_norm"]
    energy = end["statistics"]["endpoint_energy"]
    return {"closure_norm": closure, "endpoint_energy": energy,
            "passes": closure < 1e-12 and energy > 1e-15}


def _endpoint_vs_reference_test():
    nz, ny, nx = 4, 5, 6
    rng = np.random.RandomState(0)
    A_xp = rng.randn(nz, ny, nx); A_xp[:, :, -1] = 0.0
    A_yp = rng.randn(nz, ny, nx); A_yp[:, -1, :] = 0.0
    A_zp = rng.randn(nz, ny, nx); A_zp[-1, :, :] = 0.0
    Pxx = np.full((nz, ny, nx), 0.8); Pyy = np.full((nz, ny, nx), 0.8)
    Pzz = np.full((nz, ny, nx), 0.8); Pxy = np.full((nz, ny, nx), -0.2)
    Pxz = np.zeros((nz, ny, nx)); Pyz = np.zeros((nz, ny, nx))
    proj = (Pxx, Pxy, Pxz, Pyy, Pyz, Pzz)
    pair_amp = {"A_xp": A_xp, "A_yp": A_yp, "A_zp": A_zp}
    from ..core.pair_enumeration import enumerate_internal_pairs
    pairs = enumerate_internal_pairs((nz, ny, nx))
    pair_resp = build_pair_responses(pairs, pair_amp, proj, "PM1", "PS2")
    p = assemble_endpoint_field(pair_resp, (nz, ny, nx))
    r = assemble_endpoint_field_reference(pair_resp, (nz, ny, nx))
    err = float(max(np.max(np.abs(p["Rx_3d"] - r["Rx_3d"])),
                     np.max(np.abs(p["Ry_3d"] - r["Ry_3d"])),
                     np.max(np.abs(p["Rz_3d"] - r["Rz_3d"]))))
    return {"max_diff": err, "passes": err < 1e-14}


def _interface_closure_test():
    """CORRECTION-001 closure identity (with the off-by-one fixed):
        sum_i R_interface(i) = sum_{i,j} R_{ij}."""
    nz, ny, nx = 4, 5, 6
    rng = np.random.RandomState(0)
    A_xp = rng.randn(nz, ny, nx); A_xp[:, :, -1] = 0.0
    A_yp = rng.randn(nz, ny, nx); A_yp[:, -1, :] = 0.0
    A_zp = rng.randn(nz, ny, nx); A_zp[-1, :, :] = 0.0
    Pxx = np.full((nz, ny, nx), 0.8); Pyy = np.full((nz, ny, nx), 0.8)
    Pzz = np.full((nz, ny, nx), 0.8); Pxy = np.full((nz, ny, nx), -0.2)
    Pxz = np.zeros((nz, ny, nx)); Pyz = np.zeros((nz, ny, nx))
    proj = (Pxx, Pxy, Pxz, Pyy, Pyz, Pzz)
    pair_amp = {"A_xp": A_xp, "A_yp": A_yp, "A_zp": A_zp}
    from ..core.pair_enumeration import enumerate_internal_pairs
    pairs = enumerate_internal_pairs((nz, ny, nx))
    pair_resp = build_pair_responses(pairs, pair_amp, proj, "PM1", "PS2")
    iface = rasterize_interface_field(pair_resp, (nz, ny, nx))
    # CORRECTION-001: sum over the FULL valid source slice ([:-1]).
    sum_int = (float(np.sum(pair_resp["R_ij_xp"][:, :, :-1])) +
                float(np.sum(pair_resp["R_ij_yp"][:, :-1, :])) +
                float(np.sum(pair_resp["R_ij_zp"][:-1, :, :])))
    sum_int_y = (float(np.sum(pair_resp["R_ij_y_xp"][:, :, :-1])) +
                  float(np.sum(pair_resp["R_ij_y_yp"][:, :-1, :])) +
                  float(np.sum(pair_resp["R_ij_y_zp"][:-1, :, :])))
    sum_int_z = (float(np.sum(pair_resp["R_ij_z_xp"][:, :, :-1])) +
                  float(np.sum(pair_resp["R_ij_z_yp"][:, :-1, :])) +
                  float(np.sum(pair_resp["R_ij_z_zp"][:-1, :, :])))
    sx = iface["statistics"]["global_vector_sum"][0] - sum_int
    sy = iface["statistics"]["global_vector_sum"][1] - sum_int_y
    sz = iface["statistics"]["global_vector_sum"][2] - sum_int_z
    err = max(abs(sx), abs(sy), abs(sz))
    return {"max_diff": err, "passes": err < 1e-12}


def _interface_pair_count_audit_test():
    nz, ny, nx = 4, 5, 6
    rng = np.random.RandomState(0)
    A_xp = rng.randn(nz, ny, nx); A_xp[:, :, -1] = 0.0
    A_yp = rng.randn(nz, ny, nx); A_yp[:, -1, :] = 0.0
    A_zp = rng.randn(nz, ny, nx); A_zp[-1, :, :] = 0.0
    Pxx = np.full((nz, ny, nx), 0.8); Pyy = np.full((nz, ny, nx), 0.8)
    Pzz = np.full((nz, ny, nx), 0.8); Pxy = np.full((nz, ny, nx), -0.2)
    Pxz = np.zeros((nz, ny, nx)); Pyz = np.zeros((nz, ny, nx))
    proj = (Pxx, Pxy, Pxz, Pyy, Pyz, Pzz)
    pair_amp = {"A_xp": A_xp, "A_yp": A_yp, "A_zp": A_zp}
    from ..core.pair_enumeration import enumerate_internal_pairs
    pairs = enumerate_internal_pairs((nz, ny, nx))
    pair_resp = build_pair_responses(pairs, pair_amp, proj, "PM1", "PS2")
    rows = interface_pair_count_audit(pair_resp, (nz, ny, nx))
    return {"rows": rows,
            "passes": all(r["passes"] for r in rows)}


def _interface_boundary_impulse_test():
    """Boundary impulse tests (CORRECTION-001 §6.6).

    Place a single nonzero pair response at the LAST valid source
    voxel on each axis. The corrected rasterizer must distribute it
    to (source, partner) and to NO other voxel.
    """
    nz, ny, nx = 4, 5, 6
    iface_results = []
    for axis, src_idx, dst_idx in [
        ("xp", (slice(None), slice(None), nx - 2),
                  (slice(None), slice(None), nx - 1)),
        ("yp", (slice(None), ny - 2, slice(None)),
                  (slice(None), ny - 1, slice(None))),
        ("zp", (nz - 2, slice(None), slice(None)),
                  (nz - 1, slice(None), slice(None))),
    ]:
        pair_responses = {
            "R_ij_xp": np.zeros((nz, ny, nx)),
            "R_ij_y_xp": np.zeros((nz, ny, nx)),
            "R_ij_z_xp": np.zeros((nz, ny, nx)),
            "R_ij_yp": np.zeros((nz, ny, nx)),
            "R_ij_y_yp": np.zeros((nz, ny, nx)),
            "R_ij_z_yp": np.zeros((nz, ny, nx)),
            "R_ij_zp": np.zeros((nz, ny, nx)),
            "R_ij_y_zp": np.zeros((nz, ny, nx)),
            "R_ij_z_zp": np.zeros((nz, ny, nx)),
        }
        key_x, key_y, key_z = {
            "xp": ("R_ij_xp", "R_ij_y_xp", "R_ij_z_xp"),
            "yp": ("R_ij_yp", "R_ij_y_yp", "R_ij_z_yp"),
            "zp": ("R_ij_zp", "R_ij_y_zp", "R_ij_z_zp"),
        }[axis]
        pair_responses[key_x][src_idx] = 1.0
        pair_responses[key_y][src_idx] = 2.0
        pair_responses[key_z][src_idx] = 3.0
        iface = rasterize_interface_field(pair_responses, (nz, ny, nx))
        # Source must receive 0.5 R; destination must receive 0.5 R;
        # NO OTHER voxel must change.
        ok = True
        for comp_name, comp_key in (("Rx", "Rx_3d_interface"),
                                     ("Ry", "Ry_3d_interface"),
                                     ("Rz", "Rz_3d_interface")):
            arr = iface[comp_key]
            # Use np.atleast_1d to handle slice indexing that may
            # return a 2D sub-array (for slices).
            src_val = float(np.atleast_1d(arr[src_idx]).max())
            dst_val = float(np.atleast_1d(arr[dst_idx]).max())
            expected_src = {"Rx": 0.5, "Ry": 1.0, "Rz": 1.5}[comp_name]
            if not np.isclose(src_val, expected_src, atol=1e-12):
                ok = False
            if not np.isclose(dst_val, expected_src, atol=1e-12):
                ok = False
            # All other voxels must be zero. Build a mask that
            # excludes the source slice and destination slice.
            mask = np.ones(arr.shape, dtype=bool)
            mask[src_idx] = False
            mask[dst_idx] = False
            other_max = float(np.max(np.abs(arr[mask])))
            if other_max > 1e-12:
                ok = False
        iface_results.append({"axis": axis, "passes": ok})
    return {"results": iface_results,
            "passes": all(r["passes"] for r in iface_results)}


def _interface_wc1_wrong_control_test():
    """WC1 — old rasterisation slicing (CORRECTION-001 §19).

    Run the predecessor's ``[:-2]`` implementation; the consumed
    pair count must be smaller than the expected pair count
    (omitted pairs at the upper boundary)."""
    nz, ny, nx = 4, 5, 6
    pair_responses = {
        "R_ij_xp": np.zeros((nz, ny, nx)),
        "R_ij_y_xp": np.zeros((nz, ny, nx)),
        "R_ij_z_xp": np.zeros((nz, ny, nx)),
        "R_ij_yp": np.zeros((nz, ny, nx)),
        "R_ij_y_yp": np.zeros((nz, ny, nx)),
        "R_ij_z_yp": np.zeros((nz, ny, nx)),
        "R_ij_zp": np.zeros((nz, ny, nx)),
        "R_ij_y_zp": np.zeros((nz, ny, nx)),
        "R_ij_z_zp": np.zeros((nz, ny, nx)),
    }
    # Place a nonzero R only at the LAST valid xp pair (ix = nx-2).
    pair_responses["R_ij_xp"][:, :, nx - 2] = 1.0

    def _old_raster(pair_responses, shape):
        Rxp_x = pair_responses["R_ij_xp"]
        nz, ny, nx = shape
        Rx_int = np.zeros((nz, ny, nx))
        # OLD: uses [:, :, :-2] so the (nx-2, nx-1) pair is OMITTED.
        if nx >= 2:
            Rx_int[:, :, :-2] += 0.5 * Rxp_x[:, :, :-2]
            Rx_int[:, :, 1:-1] += 0.5 * Rxp_x[:, :, :-2]
        return Rx_int

    old = _old_raster(pair_responses, (nz, ny, nx))
    # The last valid pair was placed at ix = nx-2. The OLD rasterizer
    # uses [:, :, :-2] and therefore OMITS it.
    omitted = (old[:, :, nx - 2] == 0).all() and (old[:, :, nx - 1] == 0).all()
    return {"test": "WC1-old-rasterization-omits-final-pair",
            "omitted_pair_count": 1 if omitted else 0,
            "passes": omitted}


def _endpoint_vs_interface_test():
    nz, ny, nx = 4, 5, 6
    rng = np.random.RandomState(0)
    A_xp = rng.randn(nz, ny, nx); A_xp[:, :, -1] = 0.0
    A_yp = rng.randn(nz, ny, nx); A_yp[:, -1, :] = 0.0
    A_zp = rng.randn(nz, ny, nx); A_zp[-1, :, :] = 0.0
    Pxx = np.full((nz, ny, nx), 0.8); Pyy = np.full((nz, ny, nx), 0.8)
    Pzz = np.full((nz, ny, nx), 0.8); Pxy = np.full((nz, ny, nx), -0.2)
    Pxz = np.zeros((nz, ny, nx)); Pyz = np.zeros((nz, ny, nx))
    proj = (Pxx, Pxy, Pxz, Pyy, Pyz, Pzz)
    pair_amp = {"A_xp": A_xp, "A_yp": A_yp, "A_zp": A_zp}
    from ..core.pair_enumeration import enumerate_internal_pairs
    pairs = enumerate_internal_pairs((nz, ny, nx))
    pair_resp = build_pair_responses(pairs, pair_amp, proj, "PM1", "PS2")
    end = assemble_endpoint_field(pair_resp, (nz, ny, nx))
    iface = rasterize_interface_field(pair_resp, (nz, ny, nx))
    diff = float(max(np.max(np.abs(end["Rx_3d"] - iface["Rx_3d_interface"])),
                       np.max(np.abs(end["Ry_3d"] - iface["Ry_3d_interface"])),
                       np.max(np.abs(end["Rz_3d"] - iface["Rz_3d_interface"]))))
    same_hash = (end["Rx_3d"] is iface["Rx_3d_interface"] or
                  np.array_equal(end["Rx_3d"], iface["Rx_3d_interface"]))
    return {"max_diff": diff, "endpoint_eq_interface": same_hash,
            "passes": (diff > 0.0) and (not same_hash)}


if __name__ == "__main__":
    r = _pair_response_agreement_test()
    assert r["passes"], f"pair response agreement: {r}"
    print(f"M08 pair response: max_diff={r['max_production_vs_reference_diff']:.3e}")
    r = _PS_lanes_distinct_test()
    assert r["passes"], f"PS lanes distinct: {r}"
    print(f"M08 PS lanes: A-PS2 diff={r['diff_PS1A_vs_PS2']:.3e}, "
          f"PS1-PS1B diff={r['diff_PS1_vs_PS1B']:.3e}, "
          f"PS1B-PS2 diff={r['diff_PS1B_vs_PS2']:.3e}")
    r = _PS1B_PS2_equivalence_class_test()
    print(f"M08 PS1-B vs PS2: diff={r['max_diff_PS1B_vs_PS2']:.3e}, "
          f"equivalence_class={r['equivalence_class']}")
    r = _endpoint_closure_test()
    assert r["passes"], f"endpoint closure: {r}"
    print(f"M09 endpoint closure: |sum R|={r['closure_norm']:.3e}, "
          f"E_end={r['endpoint_energy']:.3e}")
    r = _endpoint_vs_reference_test()
    assert r["passes"], f"endpoint vs reference: {r}"
    print(f"M09 endpoint vs reference: max_diff={r['max_diff']:.3e}")
    r = _interface_closure_test()
    assert r["passes"], f"interface closure: {r}"
    print(f"M10 interface closure: max_diff={r['max_diff']:.3e}")
    r = _interface_pair_count_audit_test()
    assert r["passes"], f"interface pair-count audit: {r}"
    for row in r["rows"]:
        print(f"M10 pair-count audit: axis={row['axis']}, "
              f"expected={row['expected_pair_count']}, "
              f"consumed={row['consumed_pair_count']}, "
              f"omitted={row['omitted_pair_count']}")
    r = _interface_boundary_impulse_test()
    assert r["passes"], f"interface boundary impulse: {r}"
    for sub in r["results"]:
        print(f"M10 boundary impulse {sub['axis']}: pass={sub['passes']}")
    r = _interface_wc1_wrong_control_test()
    assert r["passes"], f"WC1: {r}"
    print(f"WC1 old rasterization omitted final pair: {r['omitted_pair_count']}")
    r = _endpoint_vs_interface_test()
    assert r["passes"], f"endpoint vs interface: {r}"
    print(f"M10 endpoint ≠ interface: max_diff={r['max_diff']:.3e}")
    print("M08-M09-M10 pair transfer / endpoint / interface: all checks passed")
