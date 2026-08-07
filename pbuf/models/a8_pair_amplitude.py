"""M06 — Frozen A8 Pair Amplitude.

Decomposes the frozen A8/T1 voxel update into per-pair signed
amplitudes stored on the three positive N6 directions only.
"""
from __future__ import annotations
import numpy as np

from .a8_state import (
    A8_INIT_DT, A8_INIT_OMEGA, A8_INIT_K, A8_INIT_FAST_TIMESCALE,
    A8_INIT_SLOW_TIMESCALE, neighbours6_face_reflective_3d,
)
from ..core.conventions import EPS_FLOAT, N6_POSITIVE_DIRECTIONS

__all__ = [
    "PairAmplitudeRecord", "compute_a8_pair_amplitudes",
    "compute_a8_pair_amplitudes_reference",
    "longitudinal_axis_from_scalar",
    "A8PairAmplitudeError",
]


class A8PairAmplitudeError(ValueError):
    pass


class PairAmplitudeRecord:
    """Per-pair signed amplitude record."""

    __slots__ = ("pair_id", "i_index", "j_index", "axis",
                 "source_endpoint", "destination_endpoint",
                 "direction_xyz", "slow_contrib", "fast_contrib",
                 "coupling_contrib", "amplitude", "timestep")

    def __init__(self, pair_id, i_index, j_index, axis,
                 source_endpoint, destination_endpoint, direction_xyz,
                 slow_contrib, fast_contrib, coupling_contrib,
                 amplitude, timestep):
        self.pair_id = int(pair_id)
        self.i_index = tuple(int(x) for x in i_index)
        self.j_index = tuple(int(x) for x in j_index)
        self.axis = str(axis)
        self.source_endpoint = tuple(int(x) for x in source_endpoint)
        self.destination_endpoint = tuple(int(x) for x in destination_endpoint)
        self.direction_xyz = tuple(int(x) for x in direction_xyz)
        self.slow_contrib = float(slow_contrib)
        self.fast_contrib = float(fast_contrib)
        self.coupling_contrib = float(coupling_contrib)
        self.amplitude = float(amplitude)
        self.timestep = int(timestep)


def compute_a8_pair_amplitudes(u_slow, u_fast, c_state, pair_registry,
                                  frozen_parameters=None):
    """Compute per-pair signed amplitudes A_ij.

    Parameters
    ----------
    u_slow, u_fast, c_state : ndarray of identical shape (nz, ny, nx)
        Frozen A8 state.
    pair_registry : list[PairRecord]
        Enumerated unordered pairs from M05.
    frozen_parameters : dict, optional
        Override the frozen T1 coefficients. Default uses a8_state values.

    Returns
    -------
    dict with keys:
      * A_xp, A_yp, A_zp : ndarray of shape (nz, ny, nx) with the
        per-direction signed amplitudes. Boundary cells (where there
        is no valid neighbour) are zeroed.
      * records : list[PairAmplitudeRecord]
      * statistics : dict of summary statistics
    """
    if frozen_parameters is None:
        frozen_parameters = {
            "coef_fast": A8_INIT_DT * A8_INIT_OMEGA * A8_INIT_K,
            "coef_slow": A8_INIT_DT * A8_INIT_SLOW_TIMESCALE,
            "timestep": A8_INIT_DT * A8_INIT_SLOW_TIMESCALE * 160,  # placeholder
        }

    coef_fast = float(frozen_parameters["coef_fast"])
    coef_slow = float(frozen_parameters["coef_slow"])

    nz, ny, nx = u_fast.shape
    p_fast = np.pad(u_fast, ((1, 1), (1, 1), (1, 1)), mode="reflect")
    p_slow = np.pad(u_slow, ((1, 1), (1, 1), (1, 1)), mode="reflect")
    fast_xp = p_fast[1:-1, 1:-1, 2:]
    fast_yp = p_fast[1:-1, 2:, 1:-1]
    fast_zp = p_fast[2:, 1:-1, 1:-1]
    slow_xp = p_slow[1:-1, 1:-1, 2:]
    slow_yp = p_slow[1:-1, 2:, 1:-1]
    slow_zp = p_slow[2:, 1:-1, 1:-1]

    A_xp = coef_fast * (fast_xp - u_fast) + coef_slow * (slow_xp - u_slow)
    A_yp = coef_fast * (fast_yp - u_fast) + coef_slow * (slow_yp - u_slow)
    A_zp = coef_fast * (fast_zp - u_fast) + coef_slow * (slow_zp - u_slow)

    # Boundary cells where there is no valid neighbour must be zero.
    A_xp[:, :, -1] = 0.0
    A_yp[:, -1, :] = 0.0
    A_zp[-1, :, :] = 0.0

    # Build per-pair records by sampling each axis-aligned A_ij field
    # at the source voxel of every pair.
    records = []
    pid_axis = {"xp": A_xp, "yp": A_yp, "zp": A_zp}
    for pair in pair_registry:
        iz, iy, ix = pair.i_index
        amp = float(pid_axis[pair.axis][iz, iy, ix])
        records.append(PairAmplitudeRecord(
            pair_id=pair.pair_id,
            i_index=pair.i_index, j_index=pair.j_index,
            axis=pair.axis,
            source_endpoint=pair.i_index,
            destination_endpoint=pair.j_index,
            direction_xyz=pair.direction_xyz,
            slow_contrib=coef_slow * (0.0),  # derived: see below
            fast_contrib=coef_fast * (0.0),
            coupling_contrib=0.0,
            amplitude=amp,
            timestep=160,
        ))

    # Refine slow/fast/coupling decomposition by sampling the
    # neighbour value explicitly. This is the structural T1 update
    # decomposition; the "amplitude" is the combined signed transfer.
    for rec in records:
        iz, iy, ix = rec.i_index
        jz, jy, jx = rec.j_index
        axis = rec.axis
        if axis == "xp":
            f_j = u_fast[iz, iy, jx - ix]; f_i = u_fast[iz, iy, ix]
            s_j = u_slow[iz, iy, jx - ix]; s_i = u_slow[iz, iy, ix]
        elif axis == "yp":
            f_j = u_fast[iz, jy, ix]; f_i = u_fast[iz, iy, ix]
            s_j = u_slow[iz, jy, ix]; s_i = u_slow[iz, iy, ix]
        elif axis == "zp":
            f_j = u_fast[jz, iy, ix]; f_i = u_fast[iz, iy, ix]
            s_j = u_slow[jz, iy, ix]; s_i = u_slow[iz, iy, ix]
        else:
            raise A8PairAmplitudeError(f"unknown axis {axis}")
        rec.fast_contrib = coef_fast * (f_j - f_i)
        rec.slow_contrib = coef_slow * (s_j - s_i)
        rec.coupling_contrib = rec.fast_contrib + rec.slow_contrib - rec.amplitude

    amps = np.array([r.amplitude for r in records], dtype=np.float64)
    nonzero_mask = amps != 0.0
    n_nonzero = int(nonzero_mask.sum())
    stats = {
        "n_pairs": len(records),
        "n_nonzero": n_nonzero,
        "A_min": float(amps.min()) if len(amps) else 0.0,
        "A_max": float(amps.max()) if len(amps) else 0.0,
        "A_mean": float(amps.mean()) if len(amps) else 0.0,
        "A_rms": float(np.sqrt(np.mean(amps ** 2))) if len(amps) else 0.0,
        "A_abs_sum": float(np.sum(np.abs(amps))),
        "sum_abs_nonzero": float(np.sum(np.abs(amps[nonzero_mask])))
            if n_nonzero else 0.0,
    }
    return {"A_xp": A_xp, "A_yp": A_yp, "A_zp": A_zp,
            "records": records, "statistics": stats}


def compute_a8_pair_amplitudes_reference(u_slow, u_fast, c_state,
                                            pair_registry,
                                            frozen_parameters=None):
    """Reference implementation using an explicit pair-by-pair loop.

    Used for cross-validation against the production vectorized version.
    """
    if frozen_parameters is None:
        frozen_parameters = {
            "coef_fast": A8_INIT_DT * A8_INIT_OMEGA * A8_INIT_K,
            "coef_slow": A8_INIT_DT * A8_INIT_SLOW_TIMESCALE,
        }
    coef_fast = float(frozen_parameters["coef_fast"])
    coef_slow = float(frozen_parameters["coef_slow"])

    nz, ny, nx = u_fast.shape
    A_xp = np.zeros_like(u_fast)
    A_yp = np.zeros_like(u_fast)
    A_zp = np.zeros_like(u_fast)
    records = []
    pid = 0
    for iz in range(nz):
        for iy in range(ny):
            for ix in range(nx - 1):
                pid += 1
                f_i = u_fast[iz, iy, ix]; f_j = u_fast[iz, iy, ix + 1]
                s_i = u_slow[iz, iy, ix]; s_j = u_slow[iz, iy, ix + 1]
                amp = coef_fast * (f_j - f_i) + coef_slow * (s_j - s_i)
                A_xp[iz, iy, ix] = amp
                records.append(PairAmplitudeRecord(
                    pair_id=pid, i_index=(iz, iy, ix),
                    j_index=(iz, iy, ix + 1), axis="xp",
                    source_endpoint=(iz, iy, ix),
                    destination_endpoint=(iz, iy, ix + 1),
                    direction_xyz=(+1, 0, 0),
                    slow_contrib=coef_slow * (s_j - s_i),
                    fast_contrib=coef_fast * (f_j - f_i),
                    coupling_contrib=0.0,
                    amplitude=amp, timestep=160,
                ))
    for iz in range(nz):
        for iy in range(ny - 1):
            for ix in range(nx):
                pid += 1
                f_i = u_fast[iz, iy, ix]; f_j = u_fast[iz, iy + 1, ix]
                s_i = u_slow[iz, iy, ix]; s_j = u_slow[iz, iy + 1, ix]
                amp = coef_fast * (f_j - f_i) + coef_slow * (s_j - s_i)
                A_yp[iz, iy, ix] = amp
                records.append(PairAmplitudeRecord(
                    pair_id=pid, i_index=(iz, iy, ix),
                    j_index=(iz, iy + 1, ix), axis="yp",
                    source_endpoint=(iz, iy, ix),
                    destination_endpoint=(iz, iy + 1, ix),
                    direction_xyz=(0, +1, 0),
                    slow_contrib=coef_slow * (s_j - s_i),
                    fast_contrib=coef_fast * (f_j - f_i),
                    coupling_contrib=0.0,
                    amplitude=amp, timestep=160,
                ))
    for iz in range(nz - 1):
        for iy in range(ny):
            for ix in range(nx):
                pid += 1
                f_i = u_fast[iz, iy, ix]; f_j = u_fast[iz + 1, iy, ix]
                s_i = u_slow[iz, iy, ix]; s_j = u_slow[iz + 1, iy, ix]
                amp = coef_fast * (f_j - f_i) + coef_slow * (s_j - s_i)
                A_zp[iz, iy, ix] = amp
                records.append(PairAmplitudeRecord(
                    pair_id=pid, i_index=(iz, iy, ix),
                    j_index=(iz + 1, iy, ix), axis="zp",
                    source_endpoint=(iz, iy, ix),
                    destination_endpoint=(iz + 1, iy, ix),
                    direction_xyz=(0, 0, +1),
                    slow_contrib=coef_slow * (s_j - s_i),
                    fast_contrib=coef_fast * (f_j - f_i),
                    coupling_contrib=0.0,
                    amplitude=amp, timestep=160,
                ))
    amps = np.array([r.amplitude for r in records], dtype=np.float64)
    nonzero_mask = amps != 0.0
    stats = {
        "n_pairs": len(records),
        "n_nonzero": int(nonzero_mask.sum()),
        "A_min": float(amps.min()) if len(amps) else 0.0,
        "A_max": float(amps.max()) if len(amps) else 0.0,
        "A_mean": float(amps.mean()) if len(amps) else 0.0,
        "A_rms": float(np.sqrt(np.mean(amps ** 2))) if len(amps) else 0.0,
        "A_abs_sum": float(np.sum(np.abs(amps))),
    }
    return {"A_xp": A_xp, "A_yp": A_yp, "A_zp": A_zp,
            "records": records, "statistics": stats}


def longitudinal_axis_from_scalar(scalar, spacing=(1.0, 1.0, 1.0)):
    """Compute the longitudinal unit-vector field ê_L from a scalar.

    Returns (eL_x, eL_y, eL_z, valid_mask).
    """
    dz, dy, dx = spacing
    gz, gy, gx = np.gradient(scalar, dz, dy, dx, edge_order=1)
    g_mag = np.sqrt(gx ** 2 + gy ** 2 + gz ** 2)
    valid = g_mag > 1e-12
    safe = np.where(valid, g_mag, 1.0)
    eL_x = np.where(valid, gx / safe, 0.0)
    eL_y = np.where(valid, gy / safe, 0.0)
    eL_z = np.where(valid, gz / safe, 0.0)
    return eL_x, eL_y, eL_z, valid


# ----------------------------------------------------------------------
# Self-check
# ----------------------------------------------------------------------
def _antisymmetry_view_test():
    """A_ji = -A_ij must hold as a derived endpoint view of one stored
    amplitude. We construct the partner view explicitly from the T1
    update at the partner voxel and confirm the negation.

    Only the POSITIVE-direction amplitudes are stored (xp, yp, zp).
    The negative-direction amplitudes are DERIVED: A_xm[i] = -A_xp[j]
    where j = i - axis_x. This test verifies the derivation.
    """
    nz, ny, nx = 5, 6, 7
    u_slow = np.linspace(0, 1, nz * ny * nx).reshape(nz, ny, nx)
    u_fast = u_slow + 0.1
    c_state = 0.5 * (u_slow + u_fast)
    from ..core.pair_enumeration import enumerate_internal_pairs
    pairs = enumerate_internal_pairs((nz, ny, nx))
    res = compute_a8_pair_amplitudes(u_slow, u_fast, c_state, pairs)
    n = res["statistics"]["n_nonzero"]
    # Build a map (i, j) -> A_ij from the records.
    amp_map = {(r.i_index, r.j_index): r.amplitude for r in res["records"]}

    # The stored field holds A_ij at the SOURCE voxel. The partner
    # view A_ji (from j toward i) is derived from the same T1 update
    # at the partner voxel using the SAME amplitude coefficients.
    coef_fast = A8_INIT_DT * A8_INIT_OMEGA * A8_INIT_K
    coef_slow = A8_INIT_DT * A8_INIT_SLOW_TIMESCALE
    err = 0.0
    for (i, j), a_ij in amp_map.items():
        # Compute the partner view A_ji from scratch (at voxel j
        # toward voxel i).
        iz, iy, ix = j
        jz, jy, jx = i
        f_j = u_fast[iz, iy, ix]; f_i = u_fast[jz, jy, jx]
        s_j = u_slow[iz, iy, ix]; s_i = u_slow[jz, jy, jx]
        a_ji = coef_fast * (f_i - f_j) + coef_slow * (s_i - s_j)
        err = max(err, abs(a_ij + a_ji))
    return {"antisymmetry_max_error": err, "passes": err == 0.0,
            "n_nonzero": n}


def _production_vs_reference_test():
    nz, ny, nx = 4, 5, 6
    rng = np.random.RandomState(42)
    u_slow = rng.randn(nz, ny, nx)
    u_fast = u_slow + 0.1 * rng.randn(nz, ny, nx)
    c_state = 0.5 * (u_slow + u_fast)
    from ..core.pair_enumeration import enumerate_internal_pairs
    pairs = enumerate_internal_pairs((nz, ny, nx))
    p = compute_a8_pair_amplitudes(u_slow, u_fast, c_state, pairs)
    r = compute_a8_pair_amplitudes_reference(u_slow, u_fast, c_state, pairs)
    err_xp = float(np.max(np.abs(p["A_xp"] - r["A_xp"])))
    err_yp = float(np.max(np.abs(p["A_yp"] - r["A_yp"])))
    err_zp = float(np.max(np.abs(p["A_zp"] - r["A_zp"])))
    err = max(err_xp, err_yp, err_zp)
    return {"max_production_vs_reference_diff": err,
            "n_nonzero_prod": p["statistics"]["n_nonzero"],
            "n_nonzero_ref": r["statistics"]["n_nonzero"],
            "passes": err < 1e-14}


if __name__ == "__main__":
    r = _antisymmetry_view_test()
    assert r["passes"], f"antisymmetry view failed: {r}"
    print(f"M06 antisymmetry view: err={r['antisymmetry_max_error']:.3e}")
    r = _production_vs_reference_test()
    assert r["passes"], f"production vs reference failed: {r}"
    print(f"M06 production vs reference: max_diff={r['max_production_vs_reference_diff']:.3e}")
    print("M06 A8 pair amplitude: all checks passed")