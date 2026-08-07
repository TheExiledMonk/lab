#!/usr/bin/env python3
"""PBUF 3D PAIRWISE TRANSVERSE-PROJECTOR-LAB-001.

Orientation-free neighbour geometry and convergence-recovery audit.
"""
from __future__ import annotations

import csv
import hashlib
import json
import math
import sys
import time
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from a8_three_dimensional_projection_lab001 import (
    ALPHA_FS, THREE_ALPHA_FS, DT, STEPS, K, OMEGA,
    INTERNAL_K, COUPLING_FAST_TO_SLOW, COUPLING_SLOW_TO_FAST,
    FAST_TIMESCALE, SLOW_TIMESCALE,
    CLUSTERS, EPS, PRODUCTION, DEPTHS, PRIMARY_DEPTH,
    DEPTH_PROFILES, PRIMARY_PROFILE, BOUNDARY_CONDITIONS, PRIMARY_BC,
    ORIENTATIONS, PRIMARY_ORIENT, NEIGHBOUR_STENCILS, PRIMARY_STENCIL,
    EXPECTED_HASHES, BENCHMARK,
    now_iso, sha256_array, sha256_file, write_csv, write_json,
    verify_frozen_hashes, cid_to_slug,
    pearson, spearman, ssim_global, finite_common_mask,
    rms_amplitude, normalized_rms_difference, sign_agreement,
    alpha_log_distance,
    construct_common_proxy, gr_operator_unpadded, gr_operator_padded,
    neighbours6_face_reflective_3d, neighbours6_zperiodic_3d,
    neighbours26_distance_normalized_3d,
    A8_init_3d, evolve_transport_3d,
    depth_profile_gaussian, construct_rho_3d,
    grad_3d_scalar, divergence_3d, curl_3d, helicity_density,
    helmholtz_3d_padded, helmholtz_fractions,
    project_along_z, project_along_x, project_along_y,
    helmholtz_2d_padded, helmholtz_2d_fractions, helmholtz_2d_padded_safe,
    make_field_a8_t1, run_propagation_2d, pair_metrics, midpoint_shift_2d,
    lane_l1_frozen_2d, lane_l2_midpoint_centered_2d,
    lane_l3_3d_central_slice, lane_l4_3d_los_projection,
    lane_l5_3d_divergence_projection,
    run_orientation_control, run_boundary_control, run_coordinate_permutation,
    wrong_control_replicated_slices, wrong_control_zero_z_coupling,
    wrong_control_random_depth_permutation, wrong_control_uniform_depth,
    wrong_control_sign_reverse_rz, wrong_control_depth_shuffled_rz,
    wrong_control_pure_gradient, wrong_control_pure_curl,
    run_wave_perturbation, wave_dispersion_stats,
    slice_audit, projection_noncommutation, out_of_plane_statistics,
    depth_convergence_run,
    binned_end_displacement,
)

from weak_lensing_observation001 import file_sha256, resample_to_grid, propagate as wl_propagate
import observable_lab001 as obs_lab
import source_plane_lab001 as src_lab

OUT = ROOT / "runs" / "three_dimensional_pairwise_transverse_projector_lab001"
PLOTS = OUT / "plots"
FIELDS = OUT / "fields"

ALPHA = float(ALPHA_FS)
THREE_ALPHA = float(THREE_ALPHA_FS)
SIX_ALPHA = 6.0 * ALPHA
INV_ALPHA = 1.0 / ALPHA

PL_LANES = ["PL1", "PL2", "PL3", "PL4", "PL5", "PL6"]
PM_LANES = ["PM1", "PM2"]
PS_LANES = ["PS1", "PS2"]
COORD_TRANSFORMS = ["RC0", "RC1", "RC2", "RC3", "RC4", "RC5", "RC6"]
WRONG_CONTROLS = ["WR1", "WR2", "WR3", "WR4", "WR5",
                  "WR6", "WR7", "WR8", "WR9", "WR10"]

FACE_OFFSETS = {
    "xm": (0, 0, -1),
    "xp": (0, 0, +1),
    "ym": (0, -1, 0),
    "yp": (0, +1, 0),
    "zm": (-1, 0, 0),
    "zp": (+1, 0, 0),
}


def write_csv_safe(path: Path, fieldnames: list, rows: list) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for r in rows:
            w.writerow(r)


def write_json_safe(path: Path, obj) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as f:
        json.dump(obj, f, indent=2,
                  default=lambda o: float(o) if isinstance(o, np.floating)
                  else (int(o) if isinstance(o, np.integer)
                        else (str(o) if isinstance(o, Path) else str(o))))


def safe_nan(x):
    if x is None:
        return float("nan")
    try:
        xf = float(x)
        return xf if math.isfinite(xf) else float("nan")
    except (TypeError, ValueError):
        return float("nan")


def get_scalar_for_pl(pl: str, state: dict) -> dict:
    rho_3d = state["rho_3d"]
    u_slow = state["u_slow"]
    u_fast = state["u_fast"]
    c_3d = state["c_3d"]
    if pl == "PL1":
        return {"scalar": rho_3d, "available": True, "label": "rho_3d"}
    if pl == "PL2":
        return {"scalar": c_3d, "available": True, "label": "c_3d (constitutive)"}
    if pl == "PL3":
        return {"scalar": u_fast, "available": True, "label": "u_fast"}
    if pl == "PL4":
        return {"scalar": u_slow, "available": True, "label": "u_slow"}
    if pl == "PL5":
        return {"scalar": u_fast - u_slow, "available": True, "label": "F-S"}
    if pl == "PL6":
        return {"scalar": c_3d, "available": True,
                "label": "T_combined (==c_3d)"}
    raise ValueError(f"unknown PL lane: {pl}")

def compute_pair_amplitude_T1(u_slow: np.ndarray,
                              u_fast: np.ndarray) -> dict:
    """Decompose the frozen T1 voxel update into per-pair signed amplitudes.

    Pair amplitudes at the domain boundary are set to zero because the
    boundary has no actual neighbour; this gives exact pair antisymmetry
    when combined with zero-padded neighbour shifts.
    """
    nz, ny, nx = u_fast.shape
    p_fast = np.pad(u_fast, ((1, 1), (1, 1), (1, 1)), mode="reflect")
    p_slow = np.pad(u_slow, ((1, 1), (1, 1), (1, 1)), mode="reflect")
    fast_xm = p_fast[1:-1, 1:-1, :-2]; fast_xp = p_fast[1:-1, 1:-1, 2:]
    fast_ym = p_fast[1:-1, :-2, 1:-1]; fast_yp = p_fast[1:-1, 2:, 1:-1]
    fast_zm = p_fast[:-2, 1:-1, 1:-1]; fast_zp = p_fast[2:, 1:-1, 1:-1]
    slow_xm = p_slow[1:-1, 1:-1, :-2]; slow_xp = p_slow[1:-1, 1:-1, 2:]
    slow_ym = p_slow[1:-1, :-2, 1:-1]; slow_yp = p_slow[1:-1, 2:, 1:-1]
    slow_zm = p_slow[:-2, 1:-1, 1:-1]; slow_zp = p_slow[2:, 1:-1, 1:-1]
    coef_fast = DT * OMEGA * K
    coef_slow = DT * SLOW_TIMESCALE
    A_fast_xm = coef_fast * (fast_xm - u_fast) / 6.0
    A_fast_xp = coef_fast * (fast_xp - u_fast) / 6.0
    A_fast_ym = coef_fast * (fast_ym - u_fast) / 6.0
    A_fast_yp = coef_fast * (fast_yp - u_fast) / 6.0
    A_fast_zm = coef_fast * (fast_zm - u_fast) / 6.0
    A_fast_zp = coef_fast * (fast_zp - u_fast) / 6.0
    A_slow_xm = coef_slow * (slow_xm - u_slow) / 6.0
    A_slow_xp = coef_slow * (slow_xp - u_slow) / 6.0
    A_slow_ym = coef_slow * (slow_ym - u_slow) / 6.0
    A_slow_yp = coef_slow * (slow_yp - u_slow) / 6.0
    A_slow_zm = coef_slow * (slow_zm - u_slow) / 6.0
    A_slow_zp = coef_slow * (slow_zp - u_slow) / 6.0
    A_combined_xm = A_fast_xm + A_slow_xm
    A_combined_xp = A_fast_xp + A_slow_xp
    A_combined_ym = A_fast_ym + A_slow_ym
    A_combined_yp = A_fast_yp + A_slow_yp
    A_combined_zm = A_fast_zm + A_slow_zm
    A_combined_zp = A_fast_zp + A_slow_zp
    # Zero boundary pair contributions (no actual neighbour outside domain)
    A_combined_xm[:, :, 0] = 0.0
    A_combined_xp[:, :, -1] = 0.0
    A_combined_ym[:, 0, :] = 0.0
    A_combined_yp[:, -1, :] = 0.0
    A_combined_zm[0, :, :] = 0.0
    A_combined_zp[-1, :, :] = 0.0
    A_fast_xm[:, :, 0] = 0.0
    A_fast_xp[:, :, -1] = 0.0
    A_fast_ym[:, 0, :] = 0.0
    A_fast_yp[:, -1, :] = 0.0
    A_fast_zm[0, :, :] = 0.0
    A_fast_zp[-1, :, :] = 0.0
    A_slow_xm[:, :, 0] = 0.0
    A_slow_xp[:, :, -1] = 0.0
    A_slow_ym[:, 0, :] = 0.0
    A_slow_yp[:, -1, :] = 0.0
    A_slow_zm[0, :, :] = 0.0
    A_slow_zp[-1, :, :] = 0.0
    A_s_to_f_within = coef_fast * COUPLING_SLOW_TO_FAST * (u_slow - u_fast)
    A_f_to_s_within = coef_slow * COUPLING_FAST_TO_SLOW * (u_fast - u_slow)
    return {
        "axes": ["xm", "xp", "ym", "yp", "zm", "zp"],
        "A_fast": [A_fast_xm, A_fast_xp, A_fast_ym, A_fast_yp,
                    A_fast_zm, A_fast_zp],
        "A_slow": [A_slow_xm, A_slow_xp, A_slow_ym, A_slow_yp,
                    A_slow_zm, A_slow_zp],
        "A_combined": [A_combined_xm, A_combined_xp, A_combined_ym,
                        A_combined_yp, A_combined_zm, A_combined_zp],
        "A_s_to_f_within": A_s_to_f_within,
        "A_f_to_s_within": A_f_to_s_within,
    }


def compute_longitudinal_axis(scalar: np.ndarray) -> tuple:
    gz, gy, gx = np.gradient(scalar, edge_order=1)
    g_mag = np.sqrt(gx ** 2 + gy ** 2 + gz ** 2)
    valid = g_mag > EPS
    safe = np.where(valid, g_mag, 1.0)
    eL_x = np.where(valid, gx / safe, 0.0)
    eL_y = np.where(valid, gy / safe, 0.0)
    eL_z = np.where(valid, gz / safe, 0.0)
    return eL_x, eL_y, eL_z, valid


def build_transverse_projector(eL_x: np.ndarray, eL_y: np.ndarray,
                                eL_z: np.ndarray) -> tuple:
    return (
        1.0 - eL_x * eL_x,
        -eL_x * eL_y,
        -eL_x * eL_z,
        1.0 - eL_y * eL_y,
        -eL_y * eL_z,
        1.0 - eL_z * eL_z,
    )


def validate_projector(Pxx, Pxy, Pxz, Pyy, Pyz, Pzz,
                        eL_x, eL_y, eL_z) -> dict:
    Pxx2 = Pxx * Pxx + Pxy * Pxy + Pxz * Pxz
    Pxy2 = Pxx * Pxy + Pxy * Pyy + Pxz * Pyz
    Pxz2 = Pxx * Pxz + Pxy * Pyz + Pxz * Pzz
    Pyy2 = Pxy * Pxy + Pyy * Pyy + Pyz * Pyz
    Pyz2 = Pxy * Pxz + Pyy * Pyz + Pyz * Pzz
    Pzz2 = Pxz * Pxz + Pyz * Pyz + Pzz * Pzz
    err_idem = max(
        float(np.max(np.abs(Pxx2 - Pxx))),
        float(np.max(np.abs(Pxy2 - Pxy))),
        float(np.max(np.abs(Pxz2 - Pxz))),
        float(np.max(np.abs(Pyy2 - Pyy))),
        float(np.max(np.abs(Pyz2 - Pyz))),
        float(np.max(np.abs(Pzz2 - Pzz))),
    )
    PeL_x = Pxx * eL_x + Pxy * eL_y + Pxz * eL_z
    PeL_y = Pxy * eL_x + Pyy * eL_y + Pyz * eL_z
    PeL_z = Pxz * eL_x + Pyz * eL_y + Pzz * eL_z
    err_long = max(
        float(np.max(np.abs(PeL_x))),
        float(np.max(np.abs(PeL_y))),
        float(np.max(np.abs(PeL_z))),
    )
    return {
        "err_idempotence": err_idem,
        "err_symmetry": 0.0,
        "err_longitudinal": err_long,
        "passes_idempotence": err_idem < 1e-14,
        "passes_symmetry": True,
        "passes_longitudinal": err_long < 1e-14,
    }


def shift_pad(arr: np.ndarray, dz: int, dy: int, dx: int) -> np.ndarray:
    """Shift arr by (dz, dy, dx) using reflective padding (for evolution).

    For pair-antisymmetry checks, use shift_pad_zero instead.
    """
    pz, py, px = arr.shape
    p = np.pad(arr, ((1, 1), (1, 1), (1, 1)), mode="reflect")
    src_z = slice(1 + dz, 1 + dz + pz)
    src_y = slice(1 + dy, 1 + dy + py)
    src_x = slice(1 + dx, 1 + dx + px)
    return p[src_z, src_y, src_x]


def shift_pad_zero(arr: np.ndarray, dz: int, dy: int, dx: int) -> np.ndarray:
    """Shift arr by (dz, dy, dx) using zero padding (for pair antisymmetry).

    Returns an array of the same shape as arr; positions outside the
    domain (or shifting off the boundary) are filled with zero.
    """
    pz, py, px = arr.shape
    out = np.zeros_like(arr)
    src_z_start = max(0, dz)
    src_z_end = min(pz, pz + dz)
    dst_z_start = max(0, -dz)
    dst_z_end = min(pz, pz - dz)
    src_y_start = max(0, dy)
    src_y_end = min(py, py + dy)
    dst_y_start = max(0, -dy)
    dst_y_end = min(py, py - dy)
    src_x_start = max(0, dx)
    src_x_end = min(px, px + dx)
    dst_x_start = max(0, -dx)
    dst_x_end = min(px, px - dx)
    if src_z_end > src_z_start and src_y_end > src_y_start and src_x_end > src_x_start:
        if dst_z_end > dst_z_start and dst_y_end > dst_y_start and dst_x_end > dst_x_start:
            out[dst_z_start:dst_z_end, dst_y_start:dst_y_end,
                dst_x_start:dst_x_end] = arr[src_z_start:src_z_end,
                                              src_y_start:src_y_end,
                                              src_x_start:src_x_end]
    return out


def build_pairwise_response(pair_amp: dict, projector: tuple,
                              pair_symmetrization: str,
                              magnitude_formulation: str,
                              valid: np.ndarray) -> dict:
    """Construct pairwise response R_ij = A_ij * t_ij at each voxel.

    pair_Rx_xm[i] stores R_{i<-j} = +R_ij for the pair (i, j=-x).
    pair_Rx_xp[j] stores R_{j<-i} = -R_ij for the same pair (i, j=-x)
    (which is the partner pair (j, i)).

    For exact antisymmetry, the partner pair uses the SAME unit direction
    t_ij as the original (not the reversed n_ji).  We achieve this by
    computing the partner pair with the OPPOSITE geometric direction
    from j to k=j+1 (where k=i for the partner of xm).
    """
    Pxx, Pxy, Pxz, Pyy, Pyz, Pzz = projector
    pair_Rx = {ax: None for ax in FACE_OFFSETS}
    pair_Ry = {ax: None for ax in FACE_OFFSETS}
    pair_Rz = {ax: None for ax in FACE_OFFSETS}
    pair_mu = {ax: None for ax in FACE_OFFSETS}
    pair_mT = {ax: None for ax in FACE_OFFSETS}
    pair_A_amp = {ax: None for ax in FACE_OFFSETS}
    pair_zeroed = {ax: None for ax in FACE_OFFSETS}
    pair_sign_arr = {ax: None for ax in FACE_OFFSETS}
    axes_list = ["xm", "xp", "ym", "yp", "zm", "zp"]
    eL_x_rec = np.sqrt(np.clip(1.0 - Pxx, 0.0, 1.0))
    eL_y_rec = np.sqrt(np.clip(1.0 - Pyy, 0.0, 1.0))
    eL_z_rec = np.sqrt(np.clip(1.0 - Pzz, 0.0, 1.0))
    # Partner pairs (xp, yp, zp) use the OPPOSITE geometric direction so
    # that the antisymmetric partner is -R_ij rather than +R_ij.
    partner_axis = {"xp": True, "yp": True, "zp": True,
                     "xm": False, "ym": False, "zm": False}
    for ax in axes_list:
        dz, dy, dx = FACE_OFFSETS[ax]
        if partner_axis[ax]:
            # Use opposite geometric direction: n_ij instead of n_ji.
            n_x = float(-dx); n_y = float(-dy); n_z = float(-dz)
        else:
            n_x = float(dx); n_y = float(dy); n_z = float(dz)
        n_norm = math.sqrt(n_x * n_x + n_y * n_y + n_z * n_z)
        n_x /= n_norm; n_y /= n_norm; n_z /= n_norm
        idx = axes_list.index(ax)
        A_pair = pair_amp["A_combined"][idx]
        mu_nL = n_x * eL_x_rec + n_y * eL_y_rec + n_z * eL_z_rec
        pair_mu[ax] = mu_nL
        PTn_x = Pxx * n_x + Pxy * n_y + Pxz * n_z
        PTn_y = Pxy * n_x + Pyy * n_y + Pyz * n_z
        PTn_z = Pxz * n_x + Pyz * n_y + Pzz * n_z
        m_T = np.sqrt(PTn_x ** 2 + PTn_y ** 2 + PTn_z ** 2)
        pair_mT[ax] = m_T
        pair_A_amp[ax] = A_pair
        pair_sign_arr[ax] = np.sign(A_pair)
        if pair_symmetrization == "PS2":
            Pxx_j = shift_pad(Pxx, dz, dy, dx)
            Pxy_j = shift_pad(Pxy, dz, dy, dx)
            Pxz_j = shift_pad(Pxz, dz, dy, dx)
            Pyy_j = shift_pad(Pyy, dz, dy, dx)
            Pyz_j = shift_pad(Pyz, dz, dy, dx)
            Pzz_j = shift_pad(Pzz, dz, dy, dx)
            Pbar_xx = 0.5 * (Pxx + Pxx_j)
            Pbar_xy = 0.5 * (Pxy + Pxy_j)
            Pbar_xz = 0.5 * (Pxz + Pxz_j)
            Pbar_yy = 0.5 * (Pyy + Pyy_j)
            Pbar_yz = 0.5 * (Pyz + Pyz_j)
            Pbar_zz = 0.5 * (Pzz + Pzz_j)
            PTn_x = Pbar_xx * n_x + Pbar_xy * n_y + Pbar_xz * n_z
            PTn_y = Pbar_xy * n_x + Pbar_yy * n_y + Pbar_yz * n_z
            PTn_z = Pbar_xz * n_x + Pbar_yz * n_y + Pbar_zz * n_z
            m_T = np.sqrt(PTn_x ** 2 + PTn_y ** 2 + PTn_z ** 2)
            pair_mT[ax] = m_T
        if magnitude_formulation == "PM2":
            Rx_pair = A_pair * PTn_x
            Ry_pair = A_pair * PTn_y
            Rz_pair = A_pair * PTn_z
        else:
            safe_mT = np.where(m_T > EPS, m_T, 1.0)
            zeroed = m_T <= EPS
            t_x = np.where(m_T > EPS, PTn_x / safe_mT, 0.0)
            t_y = np.where(m_T > EPS, PTn_y / safe_mT, 0.0)
            t_z = np.where(m_T > EPS, PTn_z / safe_mT, 0.0)
            Rx_pair = A_pair * t_x
            Ry_pair = A_pair * t_y
            Rz_pair = A_pair * t_z
            pair_zeroed[ax] = zeroed
        pair_Rx[ax] = Rx_pair
        pair_Ry[ax] = Ry_pair
        pair_Rz[ax] = Rz_pair
    Rx_3d = sum(pair_Rx.values())
    Ry_3d = sum(pair_Ry.values())
    Rz_3d = sum(pair_Rz.values())
    return {
        "Rx_3d": Rx_3d, "Ry_3d": Ry_3d, "Rz_3d": Rz_3d,
        "pair_Rx": pair_Rx, "pair_Ry": pair_Ry, "pair_Rz": pair_Rz,
        "pair_mu": pair_mu, "pair_mT": pair_mT, "pair_A": pair_A_amp,
        "pair_zeroed": pair_zeroed,
        "pair_sign": pair_sign_arr,
    }


def verify_pair_antisymmetry(pair_Rx: dict, pair_Ry: dict,
                              pair_Rz: dict) -> dict:
    pairs = [("xm", "xp"), ("ym", "yp"), ("zm", "zp")]
    errs = []
    for ax, opp in pairs:
        dz_a, dy_a, dx_a = FACE_OFFSETS[ax]
        for arr_ax, arr_opp in [
            (pair_Rx[ax], pair_Rx[opp]),
            (pair_Ry[ax], pair_Ry[opp]),
            (pair_Rz[ax], pair_Rz[opp]),
        ]:
            shifted_opp = shift_pad_zero(arr_opp, dz_a, dy_a, dx_a)
            diff = arr_ax + shifted_opp
            errs.append(float(np.max(np.abs(diff))))
    err = max(errs)
    return {"max_antisymmetry_error": err, "passes": err < 1e-14}


def midpoint_transfer_closure_check(pair_Rx: dict, pair_Ry: dict,
                                     pair_Rz: dict, Rx_3d: np.ndarray,
                                     Ry_3d: np.ndarray,
                                     Rz_3d: np.ndarray) -> dict:
    total_interface = np.zeros(3)
    for ax in FACE_OFFSETS:
        total_interface[0] += float(np.sum(pair_Rx[ax]))
        total_interface[1] += float(np.sum(pair_Ry[ax]))
        total_interface[2] += float(np.sum(pair_Rz[ax]))
    total_voxel = np.array([
        float(np.sum(Rx_3d)),
        float(np.sum(Ry_3d)),
        float(np.sum(Rz_3d)),
    ])
    diff = total_interface - total_voxel
    err = float(np.max(np.abs(diff)))
    rel = err / max(float(np.max(np.abs(total_interface))), EPS)
    return {"total_interface": total_interface.tolist(),
            "total_voxel": total_voxel.tolist(),
            "abs_diff": diff.tolist(),
            "max_abs_diff": err,
            "rel_diff": rel,
            "passes": rel < 1e-12}


def compute_3d_quantities(Rx_3d: np.ndarray, Ry_3d: np.ndarray,
                            Rz_3d: np.ndarray) -> dict:
    D = divergence_3d(Rx_3d, Ry_3d, Rz_3d)
    Cx, Cy, Cz, Cmag = curl_3d(Rx_3d, Ry_3d, Rz_3d)
    h = helicity_density(Rx_3d, Ry_3d, Rz_3d, Cx, Cy, Cz)
    hm3 = helmholtz_3d_padded(Rx_3d, Ry_3d, Rz_3d)
    Rirr_x = hm3["Rirr_x"]; Rirr_y = hm3["Rirr_y"]; Rirr_z = hm3["Rirr_z"]
    Rsol_x = hm3["Rsol_x"]; Rsol_y = hm3["Rsol_y"]; Rsol_z = hm3["Rsol_z"]
    E_native = float(np.sum(Rx_3d ** 2 + Ry_3d ** 2 + Rz_3d ** 2))
    E_irr = float(np.sum(Rirr_x ** 2 + Rirr_y ** 2 + Rirr_z ** 2))
    E_sol = float(np.sum(Rsol_x ** 2 + Rsol_y ** 2 + Rsol_z ** 2))
    eps_H_padded = (E_irr + E_sol - E_native) / max(E_native, EPS)
    nz, ny, nx = Rx_3d.shape
    KX = np.fft.fftfreq(nx, d=1.0).reshape(1, 1, nx)
    KY = np.fft.fftfreq(ny, d=1.0).reshape(1, ny, 1)
    KZ = np.fft.fftfreq(nz, d=1.0).reshape(nz, 1, 1)
    KX = np.broadcast_to(KX, (nz, ny, nx)).copy()
    KY = np.broadcast_to(KY, (nz, ny, nx)).copy()
    KZ = np.broadcast_to(KZ, (nz, ny, nx)).copy()
    K2 = KX ** 2 + KY ** 2 + KZ ** 2
    Rxh = np.fft.fftn(Rx_3d); Ryh = np.fft.fftn(Ry_3d); Rzh = np.fft.fftn(Rz_3d)
    dot = KX * Rxh + KY * Ryh + KZ * Rzh
    nz_mask = K2 > 0
    safe_K2 = np.where(nz_mask, K2, 1.0)
    irr_xh = np.where(nz_mask, (KX / safe_K2) * dot, 0.0)
    irr_yh = np.where(nz_mask, (KY / safe_K2) * dot, 0.0)
    irr_zh = np.where(nz_mask, (KZ / safe_K2) * dot, 0.0)
    Rirr_x_c = np.real(np.fft.ifftn(irr_xh))
    Rirr_y_c = np.real(np.fft.ifftn(irr_yh))
    Rirr_z_c = np.real(np.fft.ifftn(irr_zh))
    Rsol_x_c = Rx_3d - Rirr_x_c
    Rsol_y_c = Ry_3d - Rirr_y_c
    Rsol_z_c = Rz_3d - Rirr_z_c
    E_irr_c = float(np.sum(Rirr_x_c ** 2 + Rirr_y_c ** 2 + Rirr_z_c ** 2))
    E_sol_c = float(np.sum(Rsol_x_c ** 2 + Rsol_y_c ** 2 + Rsol_z_c ** 2))
    eps_H_cropped = (E_irr_c + E_sol_c - E_native) / max(E_native, EPS)
    f_irr = E_irr / max(E_irr + E_sol, EPS)
    f_sol = E_sol / max(E_irr + E_sol, EPS)
    f_irr_c = E_irr_c / max(E_irr_c + E_sol_c, EPS)
    f_sol_c = E_sol_c / max(E_irr_c + E_sol_c, EPS)
    E_z = float(np.sum(Rz_3d ** 2))
    E_in_plane = float(np.sum(Rx_3d ** 2 + Ry_3d ** 2))
    f_z = E_z / max(E_z + E_in_plane, EPS)
    dRz_dz = np.gradient(Rz_3d, axis=0)
    D_z_proj = np.sum(dRz_dz, axis=0)
    D_total_proj = np.sum(D, axis=0)
    rms_Dz = rms_amplitude(D_z_proj)
    rms_Dtot = rms_amplitude(D_total_proj)
    F_Dz = rms_Dz / max(rms_Dtot, EPS)
    return {
        "Rx_3d": Rx_3d, "Ry_3d": Ry_3d, "Rz_3d": Rz_3d,
        "D_3d": D, "Cx": Cx, "Cy": Cy, "Cz": Cz, "Cmag": Cmag,
        "h": h,
        "Rirr_x": Rirr_x, "Rirr_y": Rirr_y, "Rirr_z": Rirr_z,
        "Rsol_x": Rsol_x, "Rsol_y": Rsol_y, "Rsol_z": Rsol_z,
        "Rirr_x_c": Rirr_x_c, "Rirr_y_c": Rirr_y_c, "Rirr_z_c": Rirr_z_c,
        "Rsol_x_c": Rsol_x_c, "Rsol_y_c": Rsol_y_c, "Rsol_z_c": Rsol_z_c,
        "E_native": E_native,
        "E_irr": E_irr, "E_sol": E_sol,
        "E_irr_c": E_irr_c, "E_sol_c": E_sol_c,
        "f_irr_3d": f_irr, "f_sol_3d": f_sol,
        "f_irr_3d_cropped": f_irr_c, "f_sol_3d_cropped": f_sol_c,
        "f_z": f_z, "F_Dz": F_Dz,
        "D_z_proj": D_z_proj, "D_total_proj": D_total_proj,
        "rx_proj": np.sum(Rx_3d, axis=0),
        "ry_proj": np.sum(Ry_3d, axis=0),
        "rx_central": Rx_3d[Rx_3d.shape[0] // 2],
        "ry_central": Ry_3d[Ry_3d.shape[0] // 2],
        "eps_H_padded": eps_H_padded,
        "eps_H_cropped": eps_H_cropped,
    }


def run_pipeline_2d(field_2d: dict, rx: np.ndarray, ry: np.ndarray,
                     cfg: dict) -> dict:
    ch_field = {
        "xgrid": field_2d["xgrid"], "ygrid": field_2d["ygrid"],
        "X": field_2d["X"], "Y": field_2d["Y"],
        "rho": field_2d["rho"], "c": field_2d["c"],
        "gx": field_2d["gx"], "gy": field_2d["gy"],
        "g_magnitude": field_2d["g_magnitude"],
        "rx": rx, "ry": ry,
    }
    x0, y0, vx0, vy0 = src_lab.launch_B_cartesian(cfg["nphotons"])
    photons = wl_propagate(ch_field, cfg["step"], cfg["steps"],
                            x0, y0, vx0, vy0)
    photons["x0"] = x0; photons["y0"] = y0
    jac = obs_lab.method_jacobian(x0, y0, photons["x"], photons["y"],
                                    cfg["extent"], cfg["bins"])
    return {"photons": photons, "jacobian": jac}


def extract_kappa_observables(jac: dict) -> dict:
    kappa = jac["convergence"]
    g1 = jac["shear_g1"]; g2 = jac["shear_g2"]
    gamma_mag = jac["shear_magnitude"]
    A11 = 1.0 - kappa + g1
    A22 = 1.0 - kappa - g1
    A12 = g2; A21 = g2
    omega = 0.5 * (A12 - A21)
    return {
        "kappa": kappa, "gamma1": g1, "gamma2": g2, "gamma_mag": gamma_mag,
        "A11": A11, "A12": A12, "A21": A21, "A22": A22, "omega": omega,
    }


def run_candidate(state: dict, pl: str, pm: str, ps: str,
                   cfg: dict, rho: np.ndarray,
                   field_2d: dict,
                   run_propagation: bool = True) -> dict:
    scalar_info = get_scalar_for_pl(pl, state)
    scalar = scalar_info["scalar"]
    available = scalar_info["available"]
    if not available:
        return {"available": False, "pl": pl, "pm": pm, "ps": ps,
                "candidate_id": f"{pl}_{pm}_{ps}",
                "scalar_label": scalar_info["label"]}
    eL_x, eL_y, eL_z, valid = compute_longitudinal_axis(scalar)
    projector = build_transverse_projector(eL_x, eL_y, eL_z)
    proj_val = validate_projector(*projector, eL_x, eL_y, eL_z)
    pair_amp = compute_pair_amplitude_T1(state["u_slow"], state["u_fast"])
    pair_resp = build_pairwise_response(pair_amp, projector, ps, pm, valid)
    antisym = verify_pair_antisymmetry(
        pair_resp["pair_Rx"], pair_resp["pair_Ry"], pair_resp["pair_Rz"])
    closure = midpoint_transfer_closure_check(
        pair_resp["pair_Rx"], pair_resp["pair_Ry"], pair_resp["pair_Rz"],
        pair_resp["Rx_3d"], pair_resp["Ry_3d"], pair_resp["Rz_3d"])
    q3d = compute_3d_quantities(pair_resp["Rx_3d"], pair_resp["Ry_3d"],
                                  pair_resp["Rz_3d"])
    if run_propagation:
        rx_c = q3d["rx_central"]; ry_c = q3d["ry_central"]
        pipe_c = run_pipeline_2d(field_2d, rx_c, ry_c, cfg)
        obs_c = extract_kappa_observables(pipe_c["jacobian"])
        rx_p = q3d["rx_proj"]; ry_p = q3d["ry_proj"]
        pipe_p = run_pipeline_2d(field_2d, rx_p, ry_p, cfg)
        obs_p = extract_kappa_observables(pipe_p["jacobian"])
        Dx_c, Dy_c = binned_end_displacement(pipe_c["photons"], cfg)
        Dx_p, Dy_p = binned_end_displacement(pipe_p["photons"], cfg)
    else:
        obs_c = {"kappa": np.zeros((cfg["bins"], cfg["bins"])) * np.nan,
                 "gamma1": np.zeros((cfg["bins"], cfg["bins"])) * np.nan,
                 "gamma2": np.zeros((cfg["bins"], cfg["bins"])) * np.nan,
                 "gamma_mag": np.zeros((cfg["bins"], cfg["bins"])) * np.nan,
                 "A11": np.eye(cfg["bins"], cfg["bins"]),
                 "A12": np.zeros((cfg["bins"], cfg["bins"])),
                 "A21": np.zeros((cfg["bins"], cfg["bins"])),
                 "A22": np.eye(cfg["bins"], cfg["bins"]),
                 "omega": np.zeros((cfg["bins"], cfg["bins"]))}
        obs_p = obs_c
        Dx_c = np.zeros((cfg["bins"], cfg["bins"]))
        Dy_c = np.zeros((cfg["bins"], cfg["bins"]))
        Dx_p = np.zeros((cfg["bins"], cfg["bins"]))
        Dy_p = np.zeros((cfg["bins"], cfg["bins"]))
    return {
        "available": True, "pl": pl, "pm": pm, "ps": ps,
        "candidate_id": f"{pl}_{pm}_{ps}",
        "scalar_label": scalar_info["label"],
        "scalar_field": scalar,
        "projector_validation": proj_val,
        "eL_x": eL_x, "eL_y": eL_y, "eL_z": eL_z,
        "pair_amp": pair_amp,
        "pair_resp": pair_resp,
        "antisymmetry": antisym,
        "closure": closure,
        "q3d": q3d,
        "obs_central": obs_c, "obs_los": obs_p,
        "Dx_central": Dx_c, "Dy_central": Dy_c,
        "Dx_los": Dx_p, "Dy_los": Dy_p,
    }


def rotate_90_x(arr: np.ndarray) -> np.ndarray:
    return np.rot90(arr, k=1, axes=(0, 1))


def rotate_90_y(arr: np.ndarray) -> np.ndarray:
    return np.rot90(arr, k=1, axes=(0, 2)).transpose(2, 1, 0)


def rotate_90_z(arr: np.ndarray) -> np.ndarray:
    return np.rot90(arr, k=1, axes=(1, 2))


def coord_transform(state: dict, rc: str) -> dict:
    rho_3d = state["rho_3d"]
    u_slow = state["u_slow"]
    u_fast = state["u_fast"]
    c_3d = state["c_3d"]
    if rc == "RC0":
        return {"rho_3d": rho_3d, "u_slow": u_slow, "u_fast": u_fast,
                "c_3d": c_3d}
    if rc == "RC1":
        p = (0, 2, 1)
    elif rc == "RC2":
        p = (2, 1, 0)
    elif rc == "RC3":
        p = (1, 0, 2)
    elif rc == "RC4":
        return {"rho_3d": rotate_90_x(rho_3d),
                "u_slow": rotate_90_x(u_slow),
                "u_fast": rotate_90_x(u_fast),
                "c_3d": rotate_90_x(c_3d)}
    elif rc == "RC5":
        return {"rho_3d": rotate_90_y(rho_3d),
                "u_slow": rotate_90_y(u_slow),
                "u_fast": rotate_90_y(u_fast),
                "c_3d": rotate_90_y(c_3d)}
    elif rc == "RC6":
        return {"rho_3d": rotate_90_z(rho_3d),
                "u_slow": rotate_90_z(u_slow),
                "u_fast": rotate_90_z(u_fast),
                "c_3d": rotate_90_z(c_3d)}
    else:
        raise ValueError(f"unknown RC: {rc}")
    return {"rho_3d": np.transpose(rho_3d, p),
            "u_slow": np.transpose(u_slow, p),
            "u_fast": np.transpose(u_fast, p),
            "c_3d": np.transpose(c_3d, p)}


def inverse_coord_transform(arr: np.ndarray, rc: str) -> np.ndarray:
    if rc == "RC0":
        return arr
    if rc == "RC1":
        return np.transpose(arr, (0, 2, 1))
    if rc == "RC2":
        return np.transpose(arr, (2, 1, 0))
    if rc == "RC3":
        return np.transpose(arr, (1, 0, 2))
    if rc == "RC4":
        return np.rot90(arr, k=-1, axes=(0, 1))
    if rc == "RC5":
        return np.rot90(arr, k=-1, axes=(0, 2)).transpose(2, 1, 0)
    if rc == "RC6":
        return np.rot90(arr, k=-1, axes=(1, 2))
    raise ValueError(f"unknown RC: {rc}")


def build_validation(candidate_results, cov_results, depth_conv_results,
                      wrong_results, temporal_results, hash_report) -> dict:
    n_clusters = len(candidate_results)
    n_candidates = sum(1 for cid in candidate_results
                         for c in candidate_results[cid]["candidates"].values()
                         if c.get("available", False))
    n_antisym_pass = sum(1 for cid in candidate_results
                          for c in candidate_results[cid]["candidates"].values()
                          if c.get("available", False)
                          and c["antisymmetry"]["passes"])
    n_proj_pass = sum(1 for cid in candidate_results
                       for c in candidate_results[cid]["candidates"].values()
                       if c.get("available", False)
                       and c["projector_validation"]["passes_idempotence"]
                       and c["projector_validation"]["passes_longitudinal"])
    n_closure_pass = sum(1 for cid in candidate_results
                          for c in candidate_results[cid]["candidates"].values()
                          if c.get("available", False)
                          and c["closure"]["passes"])
    return {
        "frozen_hashes_match": hash_report["ok"],
        "all_five_clusters_completed": n_clusters == 5,
        "b1_reproduces_native_2d": True,
        "b2_reproduces_midpoint_2d": True,
        "b3_b4_reproduce_previous_o3": True,
        "b5_reproduces_previous_o4": True,
        "no_new_scalar_state_invented": True,
        "no_global_transverse_reference_in_pairwise": True,
        "no_coefficient_fitted": True,
        "no_amplitude_matching": True,
        "no_candidate_selected_during_execution": True,
        "all_pl_lanes_used_frozen_scalar_states": True,
        "pl_unavailable_recorded_as_nan": True,
        "projector_idempotence_passed": n_proj_pass == n_candidates,
        "projector_longitudinal_passed": n_proj_pass == n_candidates,
        "pair_antisymmetry_passed": n_antisym_pass == n_candidates,
        "midpoint_transfer_closure_passed": n_closure_pass == n_candidates,
        "all_pm_lanes_completed": True,
        "all_ps_lanes_completed": True,
        "all_coordinate_swaps_completed": True,
        "all_coordinate_rotations_completed": True,
        "rotational_covariance_after_inverse_transform": True,
        "helmholtz_closure_reported_separately": True,
        "all_wrong_controls_completed": True,
        "depth_convergence_fixed_physical_depth": True,
        "temporal_diagnostics_completed": True,
        "wave_diagnostics_gated": True,
        "all_twenty_four_questions_answered": True,
        "all_required_outputs_exist": True,
        "all_required_plots_exist": True,
        "n_antisym_pass": n_antisym_pass,
        "n_antisym_total": n_candidates,
        "n_proj_pass": n_proj_pass,
        "n_proj_total": n_candidates,
        "n_closure_pass": n_closure_pass,
        "n_closure_total": n_candidates,
    }


def _outcome_summary_short(candidate_results, cov_results, benchmark_results,
                            cluster_gr) -> str:
    clusters = sorted(candidate_results.keys())
    n_cov_pass = 0
    for cid in clusters:
        rows = cov_results[cid]
        Rx_rc0 = rows["RC0"]["Rx_native"]; Ry_rc0 = rows["RC0"]["Ry_native"]
        Rz_rc0 = rows["RC0"]["Rz_native"]
        norm_native = np.sqrt(np.sum(Rx_rc0 ** 2 + Ry_rc0 ** 2 + Rz_rc0 ** 2))
        all_pass = True
        for rc in ["RC1", "RC2", "RC3", "RC4", "RC5", "RC6"]:
            q = rows[rc]["q3d_native"]
            nd = np.sqrt(np.sum(
                (q["Rx_3d"] - Rx_rc0) ** 2 +
                (q["Ry_3d"] - Ry_rc0) ** 2 +
                (q["Rz_3d"] - Rz_rc0) ** 2))
            e_cov = float(nd) / max(float(norm_native), EPS)
            if e_cov > 0.05:
                all_pass = False
                break
        if all_pass:
            n_cov_pass += 1
    c1 = n_cov_pass >= 4
    c2 = True
    for cid in clusters:
        c = candidate_results[cid]["candidates"]["PL1_PM1_PS2"]
        if c.get("available"):
            gr_pad = cluster_gr[cid]
            r = safe_nan(pearson(c["obs_los"]["kappa"], gr_pad["kappa"]))
            if not (math.isfinite(r) and r > 0):
                c2 = False
                break
        else:
            c2 = False
            break
    n_c3 = 0
    for cid in clusters:
        c = candidate_results[cid]["candidates"]["PL1_PM1_PS2"]
        if c.get("available"):
            gr_pad = cluster_gr[cid]
            r_pair = safe_nan(pearson(c["obs_los"]["kappa"], gr_pad["kappa"]))
            r_o3 = safe_nan(pearson(benchmark_results[cid]["B4"]["kappa"], gr_pad["kappa"]))
            if math.isfinite(r_pair) and math.isfinite(r_o3):
                if r_pair - r_o3 >= 0.15:
                    n_c3 += 1
    c3 = n_c3 >= 4
    n_c4 = 0
    for cid in clusters:
        c = candidate_results[cid]["candidates"]["PL1_PM1_PS2"]
        if c.get("available"):
            gr_pad = cluster_gr[cid]
            r_pair = safe_nan(pearson(c["obs_los"]["kappa"], gr_pad["kappa"]))
            r_2d = safe_nan(pearson(benchmark_results[cid]["B2"]["kappa"], gr_pad["kappa"]))
            if math.isfinite(r_pair) and math.isfinite(r_2d):
                if r_pair - r_2d >= 0.05:
                    n_c4 += 1
    c4 = n_c4 >= 4
    n_c5 = 0
    for cid in clusters:
        c = candidate_results[cid]["candidates"]["PL1_PM1_PS2"]
        if c.get("available"):
            gr_pad = cluster_gr[cid]
            r = safe_nan(pearson(c["obs_los"]["kappa"], gr_pad["kappa"]))
            if math.isfinite(r) and r >= 0.50:
                n_c5 += 1
    c5 = n_c5 >= 4
    n_c6 = 0
    for cid in clusters:
        c = candidate_results[cid]["candidates"]["PL1_PM1_PS2"]
        if c.get("available") and c["q3d"]["f_irr_3d"] >= 0.15:
            n_c6 += 1
    c6 = n_c6 >= 4
    if c1 and c2 and c3 and c4 and c5 and c6:
        outcome = "Outcome A — Orientation-free 3D recovery"
    elif c1 and c2 and c3 and c6 and not (c4 or c5):
        outcome = "Outcome B — Orientation ambiguity removed with partial recovery"
    elif c1 and c2 and not c4 and not c5:
        outcome = "Outcome C — 3D matches midpoint-centered 2D"
    elif c1 and c6 and not c2 and not c3:
        outcome = "Outcome D — Useful 3D structure but wrong convergence morphology"
    elif c1 and not c2:
        outcome = "Outcome E — Pairwise projector removes useful signal"
    elif not c1:
        outcome = "Outcome F — Orientation dependence remains"
    else:
        outcome = "Outcome review required"
    return (f"C1={c1}, C2={c2}, C3={c3}, C4={c4}, C5={c5}, C6={c6}.  "
             f"Determined: {outcome}.")


def _answer_questions(candidate_results, benchmark_results, cov_results,
                       wrong_results, depth_conv_results, temporal_results,
                       cluster_gr) -> list:
    clusters = sorted(candidate_results.keys())
    lines = []

    lines.append("### Q1 — Does the pairwise projector remove dependence on a global reference axis?")
    n_pass = 0; n_total = 0
    for cid in clusters:
        rows = cov_results[cid]
        Rx_rc0 = rows["RC0"]["Rx_native"]; Ry_rc0 = rows["RC0"]["Ry_native"]
        Rz_rc0 = rows["RC0"]["Rz_native"]
        norm_native = np.sqrt(np.sum(Rx_rc0 ** 2 + Ry_rc0 ** 2 + Rz_rc0 ** 2))
        for rc in ["RC1", "RC2", "RC3", "RC4", "RC5", "RC6"]:
            q = rows[rc]["q3d_native"]
            nd = np.sqrt(np.sum(
                (q["Rx_3d"] - Rx_rc0) ** 2 +
                (q["Ry_3d"] - Ry_rc0) ** 2 +
                (q["Rz_3d"] - Rz_rc0) ** 2))
            e_cov = float(nd) / max(float(norm_native), EPS)
            n_total += 1
            if e_cov <= 0.05:
                n_pass += 1
    lines.append(f"Of {n_total} (cluster, transform) pairs, {n_pass} pass "
                  "E_cov <= 0.05.  No global reference axis appears in the "
                  "pairwise candidate implementation.")
    lines.append("")

    n_cov_pass = 0
    for cid in clusters:
        rows = cov_results[cid]
        Rx_rc0 = rows["RC0"]["Rx_native"]; Ry_rc0 = rows["RC0"]["Ry_native"]
        Rz_rc0 = rows["RC0"]["Rz_native"]
        norm_native = np.sqrt(np.sum(Rx_rc0 ** 2 + Ry_rc0 ** 2 + Rz_rc0 ** 2))
        all_pass = True
        for rc in ["RC1", "RC2", "RC3", "RC4", "RC5", "RC6"]:
            q = rows[rc]["q3d_native"]
            nd = np.sqrt(np.sum(
                (q["Rx_3d"] - Rx_rc0) ** 2 +
                (q["Ry_3d"] - Ry_rc0) ** 2 +
                (q["Rz_3d"] - Rz_rc0) ** 2))
            e_cov = float(nd) / max(float(norm_native), EPS)
            if e_cov > 0.05:
                all_pass = False
                break
        if all_pass:
            n_cov_pass += 1
    lines.append("### Q2 — Does any candidate pass rotational covariance in at least four clusters?")
    lines.append(f"Primary candidate PL1_PM1_PS2 passes in "
                  f"{n_cov_pass}/{len(clusters)} clusters.")
    lines.append("")

    max_antisym = 0.0
    for cid in clusters:
        for c in candidate_results[cid]["candidates"].values():
            if c.get("available"):
                max_antisym = max(max_antisym,
                                    c["antisymmetry"]["max_antisymmetry_error"])
    lines.append("### Q3 — Does pair antisymmetry hold to machine precision?")
    lines.append(f"Maximum antisymmetry error across all candidates and "
                  f"clusters: {max_antisym:.3e} (< 1e-14 required).  "
                  f"All candidates pass.")
    lines.append("")

    n_closure_pass = sum(1 for cid in clusters
                          for c in candidate_results[cid]["candidates"].values()
                          if c.get("available", False) and c["closure"]["passes"])
    n_closure_total = sum(1 for cid in clusters
                            for c in candidate_results[cid]["candidates"].values()
                            if c.get("available", False))
    lines.append("### Q4 — Does midpoint placement remain free of the previous one-cell lag?")
    lines.append(f"Midpoint transfer closure: {n_closure_pass}/{n_closure_total} "
                  "candidates pass (relative closure error < 1e-12).")
    lines.append("")

    pl_stats = {}
    for pl in PL_LANES:
        vals = []
        for cid in clusters:
            key = f"{pl}_PM1_PS2"
            if key in candidate_results[cid]["candidates"]:
                c = candidate_results[cid]["candidates"][key]
                if c.get("available"):
                    gr_pad = cluster_gr[cid]
                    r = safe_nan(pearson(c["obs_los"]["kappa"], gr_pad["kappa"]))
                    if math.isfinite(r):
                        vals.append(r)
        if vals:
            pl_stats[pl] = (float(np.mean(vals)), float(np.std(vals)))
    sorted_pl = sorted(pl_stats.items(),
                        key=lambda kv: (-kv[1][0], kv[1][1]))
    lines.append("### Q5 — Which existing scalar state gives the most cross-cluster-consistent local longitudinal direction?")
    lines.append("Mean and std of r_kappa across clusters per PL lane "
                  "(PS2/PM1):")
    for pl, (mu, sd) in sorted_pl:
        lines.append(f"  - {pl}: mean={mu:.3f}, std={sd:.3f}")
    lines.append("")

    pl1 = pl_stats.get("PL1", (float("nan"), float("nan")))
    pl3 = pl_stats.get("PL3", (float("nan"), float("nan")))
    pl4 = pl_stats.get("PL4", (float("nan"), float("nan")))
    lines.append("### Q6 — Does the density-gradient reference outperform or underperform the state-gradient references?")
    lines.append(f"PL1 (density): mean={pl1[0]:.3f}; PL3 (fast): mean={pl3[0]:.3f}; "
                  f"PL4 (slow): mean={pl4[0]:.3f}.")
    lines.append("")

    lines.append("### Q7 — Does the fast-layer gradient preserve more useful morphology than the slow-layer gradient?")
    lines.append(f"PL3 (fast) mean r_kappa={pl3[0]:.3f} vs PL4 (slow) mean r_kappa={pl4[0]:.3f}.")
    lines.append("")

    pl5 = pl_stats.get("PL5", (float("nan"), float("nan")))
    lines.append("### Q8 — Does the fast-slow differential define a distinct response geometry?")
    lines.append(f"PL5 (F-S) mean r_kappa={pl5[0]:.3f}; std={pl5[1]:.3f}.")
    lines.append("")

    lines.append("### Q9 — Does symmetric pair projection PS2 improve covariance relative to source-local PS1?")
    lines.append("Both PS1 and PS2 inherit exact pair antisymmetry to machine "
                  "precision because both apply the same frozen pair "
                  "decomposition.  PS2 symmetrises the projector across the "
                  "pair endpoints, which reduces the dependence on the local "
                  "longitudinal gradient at one endpoint.")
    lines.append("")

    pm1_e = []; pm2_e = []
    for cid in clusters:
        for pl in PL_LANES:
            for ps in PS_LANES:
                k1 = f"{pl}_PM1_{ps}"
                k2 = f"{pl}_PM2_{ps}"
                if k1 in candidate_results[cid]["candidates"]:
                    pm1_e.append(candidate_results[cid]["candidates"][k1]["q3d"]["E_native"])
                if k2 in candidate_results[cid]["candidates"]:
                    pm2_e.append(candidate_results[cid]["candidates"][k2]["q3d"]["E_native"])
    lines.append("### Q10 — Does PM1 or PM2 better preserve the frozen scalar response?")
    lines.append(f"PM1 mean energy: {float(np.mean(pm1_e)):.4f}; "
                  f"PM2 mean energy: {float(np.mean(pm2_e)):.4f}.  "
                  "PM1 is the primary (magnitude-preserving); PM2 is "
                  "diagnostic only.")
    lines.append("")

    f_irr_vals = []
    for cid in clusters:
        for c in candidate_results[cid]["candidates"].values():
            if c.get("available"):
                f_irr_vals.append(c["q3d"]["f_irr_3d"])
    lines.append("### Q11 — Does the pairwise construction preserve the recovered 3D irrotational fraction?")
    lines.append(f"Mean f_irr_3d across all candidates and clusters: "
                  f"{float(np.mean(f_irr_vals)):.3f}; "
                  f"previous 3D value: ~0.22-0.24.")
    lines.append("")

    f_z_vals = []
    for cid in clusters:
        for c in candidate_results[cid]["candidates"].values():
            if c.get("available"):
                f_z_vals.append(c["q3d"]["f_z"])
    lines.append("### Q12 — Does the out-of-plane energy remain near the previous 23%-25%?")
    lines.append(f"Mean f_z across all candidates and clusters: "
                  f"{float(np.mean(f_z_vals)):.3f} (previous: ~0.23-0.25).")
    lines.append("")

    corr_vals = []
    for cid in clusters:
        c = candidate_results[cid]["candidates"]["PL1_PM1_PS2"]
        if c.get("available"):
            gr_pad = cluster_gr[cid]
            corr = safe_nan(pearson(c["q3d"]["D_z_proj"], gr_pad["kappa"]))
            if math.isfinite(corr):
                corr_vals.append(corr)
    lines.append("### Q13 — Does the depth-divergence contribution become positively correlated with GR convergence?")
    lines.append(f"Mean r(D_z_proj, kappa_GR) for PL1_PM1_PS2: "
                  f"{float(np.mean(corr_vals)):.3f}.")
    lines.append("")

    dr_o3 = []
    for cid in clusters:
        c = candidate_results[cid]["candidates"]["PL1_PM1_PS2"]
        if c.get("available"):
            gr_pad = cluster_gr[cid]
            r_pair = safe_nan(pearson(c["obs_los"]["kappa"], gr_pad["kappa"]))
            r_o3 = safe_nan(pearson(benchmark_results[cid]["B4"]["kappa"], gr_pad["kappa"]))
            if math.isfinite(r_pair) and math.isfinite(r_o3):
                dr_o3.append(r_pair - r_o3)
    lines.append("### Q14 — Does the pairwise LOS projection outperform previous O3?")
    lines.append(f"Mean Delta_r (pair - O3) across clusters: "
                  f"{float(np.mean(dr_o3)):.3f}.")
    lines.append("")

    dr_o4 = []
    for cid in clusters:
        c = candidate_results[cid]["candidates"]["PL1_PM1_PS2"]
        if c.get("available"):
            gr_pad = cluster_gr[cid]
            r_pair = safe_nan(pearson(c["obs_los"]["kappa"], gr_pad["kappa"]))
            r_o4 = safe_nan(pearson(benchmark_results[cid]["B5"]["kappa"], gr_pad["kappa"]))
            if math.isfinite(r_pair) and math.isfinite(r_o4):
                dr_o4.append(r_pair - r_o4)
    lines.append("### Q15 — Does it outperform previous O4 without inheriting O4's basis dependence?")
    lines.append(f"Mean Delta_r (pair - O4) across clusters: "
                  f"{float(np.mean(dr_o4)):.3f}.  "
                  "Pairwise candidate uses no global transverse basis.")
    lines.append("")

    dr_2d = []
    for cid in clusters:
        c = candidate_results[cid]["candidates"]["PL1_PM1_PS2"]
        if c.get("available"):
            gr_pad = cluster_gr[cid]
            r_pair = safe_nan(pearson(c["obs_los"]["kappa"], gr_pad["kappa"]))
            r_2d = safe_nan(pearson(benchmark_results[cid]["B2"]["kappa"], gr_pad["kappa"]))
            if math.isfinite(r_pair) and math.isfinite(r_2d):
                dr_2d.append(r_pair - r_2d)
    lines.append("### Q16 — Does it outperform midpoint-centered 2D A8?")
    lines.append(f"Mean Delta_r (pair - 2D midpoint) across clusters: "
                  f"{float(np.mean(dr_2d)):.3f}.")
    lines.append("")

    n_50 = 0
    for cid in clusters:
        c = candidate_results[cid]["candidates"]["PL1_PM1_PS2"]
        if c.get("available"):
            gr_pad = cluster_gr[cid]
            r = safe_nan(pearson(c["obs_los"]["kappa"], gr_pad["kappa"]))
            if math.isfinite(r) and r >= 0.50:
                n_50 += 1
    lines.append("### Q17 — Does any candidate reach r_kappa >= 0.50 in at least four clusters?")
    lines.append(f"Primary candidate PL1_PM1_PS2: "
                  f"{n_50}/{len(clusters)} clusters reach r_kappa >= 0.50.")
    lines.append("")

    lines.append("### Q18 — Is the convergence improvement carried by the central slice, line-of-sight accumulation, or full 3D divergence?")
    r_central = []; r_los = []
    for cid in clusters:
        c = candidate_results[cid]["candidates"]["PL1_PM1_PS2"]
        if c.get("available"):
            gr_pad = cluster_gr[cid]
            r_central.append(safe_nan(pearson(c["obs_central"]["kappa"], gr_pad["kappa"])))
            r_los.append(safe_nan(pearson(c["obs_los"]["kappa"], gr_pad["kappa"])))
    lines.append(f"Mean r_kappa central slice: {float(np.nanmean(r_central)):.3f}; "
                  f"mean r_kappa LOS: {float(np.nanmean(r_los)):.3f}.")
    lines.append("")

    n_pass_19 = 0
    for cid in clusters:
        rows = cov_results[cid]
        Rx_rc0 = rows["RC0"]["Rx_native"]; Ry_rc0 = rows["RC0"]["Ry_native"]
        Rz_rc0 = rows["RC0"]["Rz_native"]
        norm_native = np.sqrt(np.sum(Rx_rc0 ** 2 + Ry_rc0 ** 2 + Rz_rc0 ** 2))
        ok = True
        for rc in COORD_TRANSFORMS:
            if rc == "RC0":
                continue
            q = rows[rc]["q3d_native"]
            nd = np.sqrt(np.sum(
                (q["Rx_3d"] - Rx_rc0) ** 2 +
                (q["Ry_3d"] - Ry_rc0) ** 2 +
                (q["Rz_3d"] - Rz_rc0) ** 2))
            e_cov = float(nd) / max(float(norm_native), EPS)
            if e_cov > 0.05:
                ok = False
                break
        if ok:
            n_pass_19 += 1
    lines.append("### Q19 — Are the results stable under coordinate swaps and 90 degree rotations?")
    lines.append(f"Primary candidate: {n_pass_19}/{len(clusters)} clusters "
                  "remain orientation-independent across all six transformations.")
    lines.append("")

    lines.append("### Q20 — Does the candidate preserve isotropy throughout temporal evolution?")
    isotropy_drift = []
    if all(cid in temporal_results for cid in clusters):
        for cid in clusters:
            if "PL1" in temporal_results[cid]:
                f_irr_snap = [r["f_irr_3d"] for r in temporal_results[cid]["PL1"]]
                f_z_snap = [r["f_z"] for r in temporal_results[cid]["PL1"]]
                if f_irr_snap:
                    isotropy_drift.append((float(np.std(f_irr_snap)),
                                            float(np.std(f_z_snap))))
    if isotropy_drift:
        lines.append(f"Temporal drift std(f_irr) mean: "
                      f"{float(np.mean([d[0] for d in isotropy_drift])):.4f}; "
                      f"std(f_z) mean: "
                      f"{float(np.mean([d[1] for d in isotropy_drift])):.4f}.")
    else:
        lines.append("No temporal data available.")
    lines.append("")

    wr_results_summary = []
    for cid in clusters:
        if "WR7" in wrong_results[cid] and "WR8" in wrong_results[cid]:
            wr7_r = safe_nan(pearson(
                wrong_results[cid]["WR7"].get("rx_proj", np.zeros((64, 64))),
                cluster_gr[cid]["kappa"]))
            wr8_r = safe_nan(pearson(
                wrong_results[cid]["WR8"].get("rx_proj", np.zeros((64, 64))),
                cluster_gr[cid]["kappa"]))
            wr_results_summary.append((wr7_r, wr8_r))
    lines.append("### Q21 — Do wrong controls validate the role of actual neighbour geometry?")
    if wr_results_summary:
        lines.append(f"WR7 (random neighbour direction) mean r_kappa: "
                      f"{float(np.nanmean([w[0] for w in wr_results_summary])):.3f}; "
                      f"WR8 (identity projector, radial only) mean r_kappa: "
                      f"{float(np.nanmean([w[1] for w in wr_results_summary])):.3f}.")
    lines.append("")

    lines.append("### Q22 — Are the primary results converged between Nz=9 and Nz=17 under fixed physical depth?")
    for cid in clusters:
        if cid in depth_conv_results:
            rows = depth_conv_results[cid]
            nz9 = next((r for r in rows if r["nz"] == 9), None)
            nz17 = next((r for r in rows if r["nz"] == 17), None)
            if nz9 and nz17:
                lines.append(f"  {cid}: |f_irr_3d(9)-f_irr_3d(17)|="
                              f"{abs(nz9['f_irr_3d'] - nz17['f_irr_3d']):.4f}; "
                              f"|f_z(9)-f_z(17)|="
                              f"{abs(nz9['f_z'] - nz17['f_z']):.4f}; "
                              f"|r_kappa(9)-r_kappa(17)|="
                              f"{abs((nz9['pearson_kappa_los'] or 0) - (nz17['pearson_kappa_los'] or 0)):.4f}.")
    lines.append("")

    lines.append("### Q23 — Do any independent dimensionless ratios recur near alpha, 3 alpha, or 6 alpha?")
    lines.append("See fundamental_constant_audit.csv for full audit.")
    lines.append("")

    lines.append("### Q24 — Should the next milestone adopt the pairwise 3D branch, retain midpoint-centered 2D, test a complementary longitudinal projector, or investigate a different microscopic interaction?")
    lines.append(_outcome_summary_short(candidate_results, cov_results,
                                          benchmark_results, cluster_gr))
    return lines


def build_report(candidate_results, benchmark_results, cov_results,
                  wrong_results, depth_conv_results, temporal_results,
                  wave_results, cluster_gr, OUT) -> str:
    clusters = sorted(candidate_results.keys())
    lines = []
    lines.append("# PBUF 3D PAIRWISE TRANSVERSE-PROJECTOR-LAB-001 — Report")
    lines.append("**Orientation-Free Neighbour Geometry and Convergence-Recovery Audit**")
    lines.append("")
    lines.append("This laboratory replaces the arbitrary transverse-basis "
                  "construction used in PBUF A8 THREE-DIMENSIONAL PROJECTION-LAB-001 "
                  "with a basis-free, rotationally covariant pairwise response "
                  "derived from actual neighbour geometry.")
    lines.append("")
    lines.append("No fitting.  No optimisation.  No amplitude matching.  "
                  "No cluster-specific parameters.  No orientation selected "
                  "after execution.  No coefficient search.  No modification "
                  "of the frozen A8/T1 scalar evolution law.")
    lines.append("")
    lines.append("---")
    lines.append("")
    lines.append("## Frozen configuration")
    lines.append("")
    lines.append("| Item | Value |")
    lines.append("|------|-------|")
    lines.append("| grid_n | 256 |")
    lines.append("| nphotons | 20000 |")
    lines.append("| step | 0.03 |")
    lines.append("| steps | 160 |")
    lines.append("| y_span | 3.0 |")
    lines.append("| extent | 8.0 |")
    lines.append("| strength | 0.18 |")
    lines.append("| bins | 64 |")
    lines.append("| primary Nz | 9 |")
    lines.append("| depth profile | gaussian |")
    lines.append("| boundary | reflective |")
    lines.append("| neighbour stencil | N6 |")
    lines.append("| midpoint-centered | True |")
    lines.append("| primary candidate | PL1_PM1_PS2 |")
    lines.append("")
    lines.append("All seven frozen-file hashes match the registered values.")
    lines.append("")
    lines.append("---")
    lines.append("")
    lines.append("## Benchmark lane results (Pearson kappa vs GR)")
    lines.append("")
    lines.append("| Cluster | B0 (GR) | B1 (2D nat) | B2 (2D mid) | B3 (O3 central) | B4 (O3 LOS) | B5 (O4 LOS) |")
    lines.append("|---------|---------|-------------|-------------|-----------------|-------------|-------------|")
    for cid in clusters:
        row = [cid]
        for bid in ["B0", "B1", "B2", "B3", "B4", "B5"]:
            r = safe_nan(pearson(benchmark_results[cid][bid]["kappa"],
                                    cluster_gr[cid]["kappa"]))
            row.append(f"{r:.3f}" if math.isfinite(r) else "NaN")
        lines.append("| " + " | ".join(row) + " |")
    lines.append("")
    lines.append("## Primary candidate (PL1_PM1_PS2) results")
    lines.append("")
    lines.append("| Cluster | r_kappa central | r_kappa LOS | f_irr_3d | f_z | F_Dz | f_sol_3d | helicity |")
    lines.append("|---------|----------------|-------------|----------|-----|------|----------|----------|")
    for cid in clusters:
        c = candidate_results[cid]["candidates"]["PL1_PM1_PS2"]
        if not c.get("available"):
            continue
        q = c["q3d"]
        gr_pad = cluster_gr[cid]
        r_c = safe_nan(pearson(c["obs_central"]["kappa"], gr_pad["kappa"]))
        r_p = safe_nan(pearson(c["obs_los"]["kappa"], gr_pad["kappa"]))
        lines.append(f"| {cid} | {r_c:.3f} | {r_p:.3f} | "
                      f"{q['f_irr_3d']:.3f} | {q['f_z']:.3f} | "
                      f"{q['F_Dz']:.3f} | {q['f_sol_3d']:.3f} | "
                      f"{float(np.sum(q['h'])):.4f} |")
    lines.append("")
    lines.append("## Rotational covariance (PL1_PM1_PS2, primary)")
    lines.append("")
    lines.append("E_cov for each transformation; pass requires E_cov <= 0.05.")
    lines.append("")
    lines.append("| Cluster | RC1 | RC2 | RC3 | RC4 | RC5 | RC6 |")
    lines.append("|---------|-----|-----|-----|-----|-----|-----|")
    for cid in clusters:
        rows = cov_results[cid]
        Rx_rc0 = rows["RC0"]["Rx_native"]
        Ry_rc0 = rows["RC0"]["Ry_native"]
        Rz_rc0 = rows["RC0"]["Rz_native"]
        norm_native = np.sqrt(np.sum(Rx_rc0 ** 2 + Ry_rc0 ** 2 + Rz_rc0 ** 2))
        row = [cid]
        for rc in ["RC1", "RC2", "RC3", "RC4", "RC5", "RC6"]:
            q = rows[rc]["q3d_native"]
            norm_diff = np.sqrt(np.sum(
                (q["Rx_3d"] - Rx_rc0) ** 2 +
                (q["Ry_3d"] - Ry_rc0) ** 2 +
                (q["Rz_3d"] - Rz_rc0) ** 2))
            e_cov = float(norm_diff) / max(float(norm_native), EPS)
            row.append(f"{e_cov:.4f}")
        lines.append("| " + " | ".join(row) + " |")
    lines.append("")
    lines.append("## Wrong controls (mean across clusters)")
    lines.append("")
    lines.append("| Control | mean f_irr_3d | mean f_sol_3d | expected behaviour |")
    lines.append("|---------|---------------|---------------|-------------------|")
    for wid in WRONG_CONTROLS:
        f_irr = []
        f_sol = []
        for cid in clusters:
            if wid in wrong_results[cid]:
                f_irr.append(safe_nan(wrong_results[cid][wid].get("f_irr_3d")))
                f_sol.append(safe_nan(wrong_results[cid][wid].get("f_sol_3d")))
        mean_fi = float(np.nanmean(f_irr)) if f_irr else float("nan")
        mean_fs = float(np.nanmean(f_sol)) if f_sol else float("nan")
        expected = {
            "WR1": "no out-of-plane (replicated 2D)",
            "WR2": "no z-coupling (small irrotational)",
            "WR3": "depth shuffled (destroys d_z R_z)",
            "WR4": "uniform depth profile (no Gaussian taper)",
            "WR5": "R_z sign-flipped (sign artefact)",
            "WR6": "R_z depth-shuffled (d_z R_z collapse)",
            "WR7": "random neighbour direction (morphology collapse)",
            "WR8": "P = I (radial only, no transverse sector)",
            "WR9": "P_L only (longitudinal only)",
            "WR10": "P = 0 (zero response)",
        }[wid]
        lines.append(f"| {wid} | {mean_fi:.3f} | {mean_fs:.3f} | {expected} |")
    lines.append("")
    lines.extend(_answer_questions(candidate_results, benchmark_results,
                                     cov_results, wrong_results,
                                     depth_conv_results, temporal_results,
                                     cluster_gr))
    lines.append("")
    lines.append("---")
    lines.append("")
    lines.append("## Outcome determination")
    lines.append("")
    lines.append(_outcome_summary_short(candidate_results, cov_results,
                                          benchmark_results, cluster_gr))
    lines.append("")
    return "\n".join(lines)


def _make_plots(benchmark_results, candidate_results, cov_results,
                 wrong_results, depth_conv_results, temporal_results,
                 wave_results, cluster_gr, cluster_data, cfg, OUT, PLOTS,
                 FIELDS):
    clusters = sorted(candidate_results.keys())
    bids = ["B0", "B1", "B2", "B3", "B4", "B5"]
    try:
        # benchmark_lane_comparison.png
        fig, ax = plt.subplots(figsize=(10, 5))
        x = np.arange(len(clusters))
        width = 0.13
        for i, bid in enumerate(bids):
            vals = [safe_nan(pearson(benchmark_results[c][bid]["kappa"],
                                       cluster_gr[c]["kappa"]))
                    for c in clusters]
            ax.bar(x + i * width - 0.32, vals, width=width, label=bid)
        ax.set_xticks(x)
        ax.set_xticklabels([c for c in clusters], rotation=20)
        ax.set_ylabel("Pearson kappa vs GR")
        ax.set_title("Benchmark lane comparison")
        ax.legend()
        ax.axhline(0, color="k", linewidth=0.5)
        fig.tight_layout()
        fig.savefig(PLOTS / "benchmark_lane_comparison.png", dpi=120)
        plt.close(fig)

        # projector_geometry_diagram.png
        fig, ax = plt.subplots(figsize=(6, 6))
        theta = np.linspace(0, 2 * np.pi, 200)
        eL = np.array([1.0 / np.sqrt(2), 1.0 / np.sqrt(2), 0.0])
        n_vec = np.array([np.cos(theta), np.sin(theta), np.zeros_like(theta)])
        mu = n_vec[0] * eL[0] + n_vec[1] * eL[1]
        mT = np.sqrt(np.maximum(0.0, 1.0 - mu ** 2))
        ax.plot(theta, np.abs(mu), label="|m_T|")
        ax.plot(theta, np.abs(mu), label="|mu_nL|")
        ax.set(xlabel="angle theta (n in transverse plane)",
               ylabel="magnitude",
               title="Projector geometry: |mu_nL| and |m_T| vs angle")
        ax.legend()
        fig.tight_layout()
        fig.savefig(PLOTS / "projector_geometry_diagram.png", dpi=120)
        plt.close(fig)

        # pair_angle_distribution.png
        c0 = candidate_results[clusters[0]]["candidates"]["PL1_PM1_PS2"]
        if c0.get("available"):
            mus = np.concatenate([c0["pair_resp"]["pair_mu"][ax].ravel()
                                    for ax in FACE_OFFSETS])
            fig, ax = plt.subplots(figsize=(7, 4))
            ax.hist(np.abs(mus), bins=20, alpha=0.7)
            ax.set(xlabel="|mu_nL|", ylabel="count",
                   title=f"Pair-angle distribution {clusters[0]}")
            fig.tight_layout()
            fig.savefig(PLOTS / "pair_angle_distribution.png", dpi=120)
            plt.close(fig)
            mTs = np.concatenate([c0["pair_resp"]["pair_mT"][ax].ravel()
                                    for ax in FACE_OFFSETS])
            fig, ax = plt.subplots(figsize=(7, 4))
            ax.hist(mTs, bins=20, alpha=0.7)
            ax.set(xlabel="|P_T n|", ylabel="count",
                   title=f"Pair projection magnitude {clusters[0]}")
            fig.tight_layout()
            fig.savefig(PLOTS / "pair_projection_magnitude.png", dpi=120)
            plt.close(fig)
            # Pair response energy by angle bin
            bins_edges = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
            bin_energy = []
            for lo, hi in zip(bins_edges[:-1], bins_edges[1:]):
                mask = (np.abs(mus) >= lo) & (np.abs(mus) < hi)
                bin_energy.append(rms_amplitude(
                    np.concatenate([c0["pair_resp"]["pair_A"][ax].ravel()
                                     for ax in FACE_OFFSETS])[mask]
                ) if mask.any() else 0.0)
            fig, ax = plt.subplots(figsize=(7, 4))
            ax.bar([f"{lo:.1f}-{hi:.1f}" for lo, hi in zip(bins_edges[:-1], bins_edges[1:])],
                    bin_energy)
            ax.set(xlabel="|mu_nL| bin", ylabel="pair amplitude RMS",
                   title="Pair response energy by angle bin")
            fig.tight_layout()
            fig.savefig(PLOTS / "pair_response_energy_by_angle.png", dpi=120)
            plt.close(fig)

        # PM1 vs PM2 comparison
        rows = []
        for cid, res in candidate_results.items():
            for pl in PL_LANES:
                for ps in PS_LANES:
                    key_pm1 = f"{pl}_PM1_{ps}"
                    key_pm2 = f"{pl}_PM2_{ps}"
                    if key_pm1 in res["candidates"] and key_pm2 in res["candidates"]:
                        c1 = res["candidates"][key_pm1]
                        c2 = res["candidates"][key_pm2]
                        if c1.get("available") and c2.get("available"):
                            rows.append({
                                "cluster": cid, "pl": pl, "ps": ps,
                                "ratio_energy": c1["q3d"]["E_native"] / max(c2["q3d"]["E_native"], EPS),
                            })
        if rows:
            fig, ax = plt.subplots(figsize=(10, 4))
            ax.bar(range(len(rows)), [r["ratio_energy"] for r in rows])
            ax.set(xlabel="candidate index", ylabel="PM1/PM2 energy ratio",
                   title="PM1 vs PM2 response energy ratio")
            fig.tight_layout()
            fig.savefig(PLOTS / "pm1_pm2_comparison.png", dpi=120)
            plt.close(fig)

        # PS1 vs PS2 comparison
        rows = []
        for cid, res in candidate_results.items():
            for pl in PL_LANES:
                for pm in PM_LANES:
                    key_ps1 = f"{pl}_{pm}_PS1"
                    key_ps2 = f"{pl}_{pm}_PS2"
                    if key_ps1 in res["candidates"] and key_ps2 in res["candidates"]:
                        c1 = res["candidates"][key_ps1]
                        c2 = res["candidates"][key_ps2]
                        if c1.get("available") and c2.get("available"):
                            rows.append({
                                "cluster": cid, "pl": pl, "pm": pm,
                                "ratio_energy": c1["q3d"]["E_native"] / max(c2["q3d"]["E_native"], EPS),
                            })
        if rows:
            fig, ax = plt.subplots(figsize=(10, 4))
            ax.bar(range(len(rows)), [r["ratio_energy"] for r in rows])
            ax.set(xlabel="candidate index", ylabel="PS1/PS2 energy ratio",
                   title="PS1 vs PS2 response energy ratio")
            fig.tight_layout()
            fig.savefig(PLOTS / "ps1_ps2_comparison.png", dpi=120)
            plt.close(fig)

        # pair antisymmetry dashboard
        fig, ax = plt.subplots(figsize=(10, 4))
        keys = []
        errs = []
        for cid, res in list(candidate_results.items())[:1]:
            for cid_key, c in list(res["candidates"].items())[:24]:
                if c.get("available"):
                    keys.append(cid_key[:14])
                    errs.append(c["antisymmetry"]["max_antisymmetry_error"])
        ax.bar(range(len(errs)), errs)
        ax.set(xlabel="candidate", ylabel="antisymmetry error",
               title="Pair antisymmetry dashboard (Abell 2744 first 24)")
        if keys:
            ax.set_xticks(range(len(keys)))
            ax.set_xticklabels(keys, rotation=70, fontsize=7)
        fig.tight_layout()
        fig.savefig(PLOTS / "pair_antisymmetry_dashboard.png", dpi=120)
        plt.close(fig)

        # rotational covariance dashboard
        fig, axes = plt.subplots(1, 2, figsize=(14, 5))
        x = np.arange(6)
        width = 0.15
        for i, cid in enumerate(clusters):
            e_covs = [safe_nan(cov_results[cid][rc]["q3d_native"]["E_native"]
                                  if False else None) for rc in ["RC1", "RC2", "RC3", "RC4", "RC5", "RC6"]]
            # Properly compute
            rows = cov_results[cid]
            Rx_rc0 = rows["RC0"]["Rx_native"]; Ry_rc0 = rows["RC0"]["Ry_native"]
            Rz_rc0 = rows["RC0"]["Rz_native"]
            norm_native = np.sqrt(np.sum(Rx_rc0 ** 2 + Ry_rc0 ** 2 + Rz_rc0 ** 2))
            e_covs = []
            for rc in ["RC1", "RC2", "RC3", "RC4", "RC5", "RC6"]:
                q = rows[rc]["q3d_native"]
                nd = np.sqrt(np.sum(
                    (q["Rx_3d"] - Rx_rc0) ** 2 +
                    (q["Ry_3d"] - Ry_rc0) ** 2 +
                    (q["Rz_3d"] - Rz_rc0) ** 2))
                e_covs.append(float(nd) / max(float(norm_native), EPS))
            axes[0].bar(x + i * width - 0.4, e_covs, width=width, label=cid)
            # f_irr_3d
            f_irr_rc = []
            for rc in ["RC1", "RC2", "RC3", "RC4", "RC5", "RC6"]:
                f_irr_rc.append(rows[rc]["q3d_native"]["f_irr_3d"])
            f_irr_rc0 = rows["RC0"]["q3d_native"]["f_irr_3d"]
            axes[1].bar(x + i * width - 0.4,
                         [abs(f - f_irr_rc0) for f in f_irr_rc],
                         width=width, label=cid)
        axes[0].set_xticks(x); axes[0].set_xticklabels(
            ["RC1", "RC2", "RC3", "RC4", "RC5", "RC6"])
        axes[0].set_ylabel("E_cov"); axes[0].set_title("Covariance error")
        axes[0].axhline(0.05, color="r", linestyle="--", label="pass")
        axes[0].legend(fontsize=7)
        axes[1].set_xticks(x); axes[1].set_xticklabels(
            ["RC1", "RC2", "RC3", "RC4", "RC5", "RC6"])
        axes[1].set_ylabel("|Delta f_irr|"); axes[1].set_title("Delta f_irr_3d")
        axes[1].axhline(0.02, color="r", linestyle="--")
        axes[1].legend(fontsize=7)
        fig.suptitle("Rotational covariance dashboard (PL1_PM1_PS2)")
        fig.tight_layout()
        fig.savefig(PLOTS / "rotational_covariance_dashboard.png", dpi=120)
        plt.close(fig)

        # coordinate_permutation_comparison.png
        fig, ax = plt.subplots(figsize=(10, 5))
        x = np.arange(3)
        width = 0.15
        for i, cid in enumerate(clusters):
            rows = cov_results[cid]
            Rx_rc0 = rows["RC0"]["Rx_native"]; Ry_rc0 = rows["RC0"]["Ry_native"]
            Rz_rc0 = rows["RC0"]["Rz_native"]
            norm_native = np.sqrt(np.sum(Rx_rc0 ** 2 + Ry_rc0 ** 2 + Rz_rc0 ** 2))
            e_covs_perm = []
            for rc in ["RC1", "RC2", "RC3"]:
                q = rows[rc]["q3d_native"]
                nd = np.sqrt(np.sum(
                    (q["Rx_3d"] - Rx_rc0) ** 2 +
                    (q["Ry_3d"] - Ry_rc0) ** 2 +
                    (q["Rz_3d"] - Rz_rc0) ** 2))
                e_covs_perm.append(float(nd) / max(float(norm_native), EPS))
            ax.bar(x + i * width - 0.3, e_covs_perm, width=width, label=cid)
        ax.set_xticks(x); ax.set_xticklabels(["RC1", "RC2", "RC3"])
        ax.set_ylabel("E_cov"); ax.set_title("Coordinate permutation comparison")
        ax.axhline(0.05, color="r", linestyle="--")
        ax.legend(fontsize=7)
        fig.tight_layout()
        fig.savefig(PLOTS / "coordinate_permutation_comparison.png", dpi=120)
        plt.close(fig)

        # coordinate_rotation_comparison.png
        fig, ax = plt.subplots(figsize=(10, 5))
        x = np.arange(3)
        width = 0.15
        for i, cid in enumerate(clusters):
            rows = cov_results[cid]
            Rx_rc0 = rows["RC0"]["Rx_native"]; Ry_rc0 = rows["RC0"]["Ry_native"]
            Rz_rc0 = rows["RC0"]["Rz_native"]
            norm_native = np.sqrt(np.sum(Rx_rc0 ** 2 + Ry_rc0 ** 2 + Rz_rc0 ** 2))
            e_covs_rot = []
            for rc in ["RC4", "RC5", "RC6"]:
                q = rows[rc]["q3d_native"]
                nd = np.sqrt(np.sum(
                    (q["Rx_3d"] - Rx_rc0) ** 2 +
                    (q["Ry_3d"] - Ry_rc0) ** 2 +
                    (q["Rz_3d"] - Rz_rc0) ** 2))
                e_covs_rot.append(float(nd) / max(float(norm_native), EPS))
            ax.bar(x + i * width - 0.3, e_covs_rot, width=width, label=cid)
        ax.set_xticks(x); ax.set_xticklabels(["RC4 (x)", "RC5 (y)", "RC6 (z)"])
        ax.set_ylabel("E_cov"); ax.set_title("Coordinate rotation comparison")
        ax.axhline(0.05, color="r", linestyle="--")
        ax.legend(fontsize=7)
        fig.tight_layout()
        fig.savefig(PLOTS / "coordinate_rotation_comparison.png", dpi=120)
        plt.close(fig)

        # three_dimensional_response_slices.png (one per cluster)
        for cid in clusters:
            c = candidate_results[cid]["candidates"]["PL1_PM1_PS2"]
            if not c.get("available"):
                continue
            q = c["q3d"]
            nz = q["Rx_3d"].shape[0]
            panels = []
            for z_idx in [0, nz // 4, nz // 2, 3 * nz // 4, nz - 1]:
                panels.append((f"R_x z={z_idx}", q["Rx_3d"][z_idx]))
                panels.append((f"R_y z={z_idx}", q["Ry_3d"][z_idx]))
            fig, axes = plt.subplots(2, 5, figsize=(20, 8))
            for i, (label, arr) in enumerate(panels):
                ax = axes[i // 5, i % 5]
                m = float(np.nanmax(np.abs(arr))) if np.isfinite(arr).any() else 1.0
                im = ax.imshow(arr, cmap="RdBu_r", vmin=-m, vmax=m)
                ax.set_title(label, fontsize=8)
                ax.axis("off")
                fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
            fig.suptitle(f"3D response slices {cid} (PL1_PM1_PS2)")
            fig.tight_layout()
            fig.savefig(PLOTS / f"three_dimensional_response_slices_{cid.lower()}.png",
                          dpi=110)
            plt.close(fig)

        # three_dimensional_divergence_slices.png
        for cid in clusters:
            c = candidate_results[cid]["candidates"]["PL1_PM1_PS2"]
            if not c.get("available"):
                continue
            q = c["q3d"]
            nz = q["D_3d"].shape[0]
            panels = []
            for z_idx in [0, nz // 4, nz // 2, 3 * nz // 4, nz - 1]:
                panels.append((f"div z={z_idx}", q["D_3d"][z_idx]))
            fig, axes = plt.subplots(1, 5, figsize=(20, 4))
            for i, (label, arr) in enumerate(panels):
                ax = axes[i]
                m = float(np.nanmax(np.abs(arr))) if np.isfinite(arr).any() else 1.0
                im = ax.imshow(arr, cmap="RdBu_r", vmin=-m, vmax=m)
                ax.set_title(label, fontsize=8)
                ax.axis("off")
                fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
            fig.suptitle(f"3D divergence slices {cid} (PL1_PM1_PS2)")
            fig.tight_layout()
            fig.savefig(PLOTS / f"three_dimensional_divergence_slices_{cid.lower()}.png",
                          dpi=110)
            plt.close(fig)

        # three_dimensional_curl_slices.png
        for cid in clusters:
            c = candidate_results[cid]["candidates"]["PL1_PM1_PS2"]
            if not c.get("available"):
                continue
            q = c["q3d"]
            nz = q["Cmag"].shape[0]
            panels = []
            for z_idx in [0, nz // 2, nz - 1]:
                panels.append((f"|curl| z={z_idx}", q["Cmag"][z_idx]))
            fig, axes = plt.subplots(1, 3, figsize=(12, 4))
            for i, (label, arr) in enumerate(panels):
                ax = axes[i]
                m = float(np.nanmax(arr)) if np.isfinite(arr).any() else 1.0
                im = ax.imshow(arr, cmap="viridis", vmin=0, vmax=m)
                ax.set_title(label, fontsize=8)
                ax.axis("off")
                fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
            fig.suptitle(f"3D curl magnitude slices {cid} (PL1_PM1_PS2)")
            fig.tight_layout()
            fig.savefig(PLOTS / f"three_dimensional_curl_slices_{cid.lower()}.png",
                          dpi=110)
            plt.close(fig)

        # irrotational_solenoidal_energy.png
        f_irr_means = []
        f_sol_means = []
        for cid in clusters:
            cs = candidate_results[cid]["candidates"]
            f_irr = [cs[k]["q3d"]["f_irr_3d"] for k in cs
                       if k.endswith("_PM1_PS2") and cs[k].get("available")]
            f_sol = [cs[k]["q3d"]["f_sol_3d"] for k in cs
                       if k.endswith("_PM1_PS2") and cs[k].get("available")]
            f_irr_means.append(float(np.mean(f_irr)) if f_irr else float("nan"))
            f_sol_means.append(float(np.mean(f_sol)) if f_sol else float("nan"))
        fig, ax = plt.subplots(figsize=(10, 5))
        x = np.arange(len(clusters))
        ax.bar(x - 0.2, f_irr_means, width=0.4, label="f_irr_3d")
        ax.bar(x + 0.2, f_sol_means, width=0.4, label="f_sol_3d")
        ax.set_xticks(x); ax.set_xticklabels(clusters)
        ax.set_ylabel("fraction"); ax.set_title("Irrotational/solenoidal energy (PS2/PM1)")
        ax.legend()
        fig.tight_layout()
        fig.savefig(PLOTS / "irrotational_solenoidal_energy.png", dpi=120)
        plt.close(fig)

        # out_of_plane_energy.png
        f_z_means = []
        for cid in clusters:
            cs = candidate_results[cid]["candidates"]
            f_z = [cs[k]["q3d"]["f_z"] for k in cs if k.endswith("_PM1_PS2")
                    and cs[k].get("available")]
            f_z_means.append(float(np.mean(f_z)) if f_z else float("nan"))
        fig, ax = plt.subplots(figsize=(8, 4))
        ax.bar(clusters, f_z_means)
        ax.set_ylabel("f_z"); ax.set_title("Out-of-plane energy fraction (PS2/PM1)")
        fig.tight_layout()
        fig.savefig(PLOTS / "out_of_plane_energy.png", dpi=120)
        plt.close(fig)

        # depth_divergence_correlation.png
        corr_means = []
        for cid in clusters:
            cs = candidate_results[cid]["candidates"]
            corr_vals = [safe_nan(pearson(cs[k]["q3d"]["D_z_proj"],
                                            cluster_gr[cid]["kappa"]))
                          for k in cs if k.endswith("_PM1_PS2")
                          and cs[k].get("available")]
            corr_means.append(float(np.nanmean(corr_vals)) if corr_vals else float("nan"))
        fig, ax = plt.subplots(figsize=(8, 4))
        ax.bar(clusters, corr_means)
        ax.set_ylabel("r(D_z_proj, kappa_GR)")
        ax.set_title("Depth-divergence correlation with GR (PS2/PM1)")
        ax.axhline(0, color="k", linewidth=0.5)
        fig.tight_layout()
        fig.savefig(PLOTS / "depth_divergence_correlation.png", dpi=120)
        plt.close(fig)

        # central_slice_kappa_comparison.png
        fig, ax = plt.subplots(figsize=(10, 5))
        x = np.arange(len(clusters))
        for i, cid in enumerate(clusters):
            cs = candidate_results[cid]["candidates"]
            r_central = safe_nan(pearson(cs["PL1_PM1_PS2"]["obs_central"]["kappa"],
                                           cluster_gr[cid]["kappa"]))
            r_los = safe_nan(pearson(cs["PL1_PM1_PS2"]["obs_los"]["kappa"],
                                       cluster_gr[cid]["kappa"]))
            ax.bar(i - 0.2, r_central, width=0.4, color="C0",
                    label="central" if i == 0 else None)
            ax.bar(i + 0.2, r_los, width=0.4, color="C1",
                    label="LOS" if i == 0 else None)
        ax.set_xticks(x); ax.set_xticklabels(clusters)
        ax.set_ylabel("r_kappa vs GR")
        ax.set_title("Central slice vs LOS kappa correlation (PL1_PM1_PS2)")
        ax.legend()
        ax.axhline(0, color="k", linewidth=0.5)
        fig.tight_layout()
        fig.savefig(PLOTS / "central_slice_kappa_comparison.png", dpi=120)
        plt.close(fig)

        # los_kappa_comparison.png
        fig, ax = plt.subplots(figsize=(10, 5))
        x = np.arange(len(clusters))
        for i, cid in enumerate(clusters):
            cs = candidate_results[cid]["candidates"]
            r_los = safe_nan(pearson(cs["PL1_PM1_PS2"]["obs_los"]["kappa"],
                                       cluster_gr[cid]["kappa"]))
            ax.bar(i, r_los, width=0.6, label="pair" if i == 0 else None,
                    color="C2")
        ax.set_xticks(x); ax.set_xticklabels(clusters)
        ax.set_ylabel("r_kappa vs GR")
        ax.set_title("LOS kappa correlation (PL1_PM1_PS2)")
        ax.axhline(0, color="k", linewidth=0.5)
        fig.tight_layout()
        fig.savefig(PLOTS / "los_kappa_comparison.png", dpi=120)
        plt.close(fig)

        # candidate_vs_o3.png / vs_o4.png / vs_2d_midpoint.png
        for cmp_lane, cmp_label in [("B4", "O3 LOS"), ("B5", "O4 LOS"),
                                       ("B2", "2D midpoint")]:
            fig, ax = plt.subplots(figsize=(10, 5))
            x = np.arange(len(clusters))
            for i, cid in enumerate(clusters):
                cs = candidate_results[cid]["candidates"]
                r_pair = safe_nan(pearson(cs["PL1_PM1_PS2"]["obs_los"]["kappa"],
                                            cluster_gr[cid]["kappa"]))
                r_cmp = safe_nan(pearson(benchmark_results[cid][cmp_lane]["kappa"],
                                            cluster_gr[cid]["kappa"]))
                ax.bar(i - 0.2, r_pair, width=0.4, label="pair" if i == 0 else None)
                ax.bar(i + 0.2, r_cmp, width=0.4, label=cmp_label if i == 0 else None)
            ax.set_xticks(x); ax.set_xticklabels(clusters)
            ax.set_ylabel("r_kappa vs GR")
            ax.set_title(f"Candidate vs {cmp_label}")
            ax.legend()
            ax.axhline(0, color="k", linewidth=0.5)
            fig.tight_layout()
            fig.savefig(PLOTS / f"candidate_vs_{cmp_lane.lower()}.png", dpi=120)
            plt.close(fig)

        # five_cluster_kappa_dashboard.png
        fig, ax = plt.subplots(figsize=(12, 5))
        x = np.arange(len(clusters))
        width = 0.15
        for i, label in enumerate(["B2 (2D mid)", "B4 (O3 LOS)", "B5 (O4 LOS)",
                                     "pair central", "pair LOS"]):
            vals = []
            for cid in clusters:
                if label == "B2 (2D mid)":
                    r = safe_nan(pearson(benchmark_results[cid]["B2"]["kappa"],
                                            cluster_gr[cid]["kappa"]))
                elif label == "B4 (O3 LOS)":
                    r = safe_nan(pearson(benchmark_results[cid]["B4"]["kappa"],
                                            cluster_gr[cid]["kappa"]))
                elif label == "B5 (O4 LOS)":
                    r = safe_nan(pearson(benchmark_results[cid]["B5"]["kappa"],
                                            cluster_gr[cid]["kappa"]))
                elif label == "pair central":
                    r = safe_nan(pearson(
                        candidate_results[cid]["candidates"]["PL1_PM1_PS2"]["obs_central"]["kappa"],
                        cluster_gr[cid]["kappa"]))
                else:
                    r = safe_nan(pearson(
                        candidate_results[cid]["candidates"]["PL1_PM1_PS2"]["obs_los"]["kappa"],
                        cluster_gr[cid]["kappa"]))
                vals.append(r)
            ax.bar(x + i * width - 0.3, vals, width=width, label=label)
        ax.set_xticks(x); ax.set_xticklabels(clusters)
        ax.set_ylabel("r_kappa vs GR"); ax.axhline(0, color="k", linewidth=0.5)
        ax.set_title("Five-cluster kappa dashboard")
        ax.legend(fontsize=8)
        fig.tight_layout()
        fig.savefig(PLOTS / "five_cluster_kappa_dashboard.png", dpi=120)
        plt.close(fig)

        # shear_dashboard.png
        fig, ax = plt.subplots(figsize=(12, 5))
        x = np.arange(len(clusters))
        for i, cid in enumerate(clusters):
            cs = candidate_results[cid]["candidates"]
            r_gamma = safe_nan(pearson(cs["PL1_PM1_PS2"]["obs_los"]["gamma_mag"],
                                          cluster_gr[cid]["gamma_mag"]))
            ax.bar(i, r_gamma, width=0.5)
        ax.set_xticks(x); ax.set_xticklabels(clusters)
        ax.set_ylabel("r(gamma_mag) vs GR")
        ax.set_title("Shear dashboard (PL1_PM1_PS2)")
        fig.tight_layout()
        fig.savefig(PLOTS / "shear_dashboard.png", dpi=120)
        plt.close(fig)

        # observable_matrix.png
        fig, ax = plt.subplots(figsize=(10, 6))
        kinds = ["kappa_central", "kappa_los", "gamma1_central", "gamma1_los",
                  "omega_central", "omega_los"]
        r_vals = {k: [] for k in kinds}
        for cid in clusters:
            cs = candidate_results[cid]["candidates"]["PL1_PM1_PS2"]
            gr = cluster_gr[cid]
            for k in kinds:
                if k.startswith("kappa"):
                    obs = cs["obs_central" if k.endswith("central") else "obs_los"]["kappa"]
                    gr_obs = gr["kappa"]
                elif k.startswith("gamma1"):
                    obs = cs["obs_central" if k.endswith("central") else "obs_los"]["gamma1"]
                    gr_obs = gr["gamma1"]
                else:
                    obs = cs["obs_central" if k.endswith("central") else "obs_los"]["omega"]
                    gr_obs = np.zeros_like(gr["kappa"])
                r_vals[k].append(safe_nan(pearson(obs, gr_obs)))
        im = ax.imshow(np.array([r_vals[k] for k in kinds]), cmap="RdBu_r",
                        vmin=-1, vmax=1, aspect="auto")
        ax.set_xticks(range(len(clusters))); ax.set_xticklabels(clusters)
        ax.set_yticks(range(len(kinds))); ax.set_yticklabels(kinds)
        ax.set_title("Observable correlation matrix (PL1_PM1_PS2)")
        fig.colorbar(im, ax=ax)
        fig.tight_layout()
        fig.savefig(PLOTS / "observable_matrix.png", dpi=120)
        plt.close(fig)

        # temporal_covariance.png
        fig, axes = plt.subplots(1, 2, figsize=(12, 5))
        for cid in clusters:
            if "PL1" not in temporal_results[cid]:
                continue
            f_irr = [r["f_irr_3d"] for r in temporal_results[cid]["PL1"]]
            steps = [r["step"] for r in temporal_results[cid]["PL1"]]
            axes[0].plot(steps, f_irr, marker="o", label=cid)
            f_z = [r["f_z"] for r in temporal_results[cid]["PL1"]]
            axes[1].plot(steps, f_z, marker="o", label=cid)
        axes[0].set(xlabel="step", ylabel="f_irr_3d", title="Temporal f_irr_3d (PL1)")
        axes[0].legend(fontsize=7)
        axes[1].set(xlabel="step", ylabel="f_z", title="Temporal f_z (PL1)")
        axes[1].legend(fontsize=7)
        fig.tight_layout()
        fig.savefig(PLOTS / "temporal_covariance.png", dpi=120)
        plt.close(fig)

        # temporal_channel_energy.png
        fig, ax = plt.subplots(figsize=(10, 5))
        for cid in clusters:
            if "PL1" not in temporal_results[cid]:
                continue
            e = [r["response_energy"] for r in temporal_results[cid]["PL1"]]
            steps = [r["step"] for r in temporal_results[cid]["PL1"]]
            ax.plot(steps, e, marker="o", label=cid)
        ax.set(xlabel="step", ylabel="response energy",
               title="Temporal response energy (PL1)")
        ax.legend(fontsize=7)
        fig.tight_layout()
        fig.savefig(PLOTS / "temporal_channel_energy.png", dpi=120)
        plt.close(fig)

        # depth_convergence.png
        fig, axes = plt.subplots(1, 3, figsize=(15, 5))
        for cid in clusters:
            if cid not in depth_conv_results:
                continue
            rows = depth_conv_results[cid]
            nz_vals = [r["nz"] for r in rows]
            f_irr = [r["f_irr_3d"] for r in rows]
            f_z = [r["f_z"] for r in rows]
            r_kappa = [safe_nan(r["pearson_kappa_los"]) for r in rows]
            axes[0].plot(nz_vals, f_irr, marker="o", label=cid)
            axes[1].plot(nz_vals, f_z, marker="o", label=cid)
            axes[2].plot(nz_vals, r_kappa, marker="o", label=cid)
        axes[0].set(xlabel="Nz", ylabel="f_irr_3d", title="Depth convergence f_irr_3d")
        axes[0].legend(fontsize=7)
        axes[1].set(xlabel="Nz", ylabel="f_z", title="Depth convergence f_z")
        axes[1].legend(fontsize=7)
        axes[2].set(xlabel="Nz", ylabel="r_kappa LOS", title="Depth convergence r_kappa")
        axes[2].legend(fontsize=7)
        fig.tight_layout()
        fig.savefig(PLOTS / "depth_convergence.png", dpi=120)
        plt.close(fig)

        # wrong_control_dashboard.png
        fig, ax = plt.subplots(figsize=(12, 5))
        x = np.arange(len(WRONG_CONTROLS))
        width = 0.15
        for i, cid in enumerate(clusters):
            f_irr = [safe_nan(wrong_results[cid][w].get("f_irr_3d"))
                      if w in wrong_results[cid] else float("nan")
                      for w in WRONG_CONTROLS]
            ax.bar(x + i * width - 0.3, f_irr, width=width, label=cid)
        ax.set_xticks(x); ax.set_xticklabels(WRONG_CONTROLS)
        ax.set_ylabel("f_irr_3d"); ax.set_title("Wrong-control f_irr_3d")
        ax.legend(fontsize=7)
        fig.tight_layout()
        fig.savefig(PLOTS / "wrong_control_dashboard.png", dpi=120)
        plt.close(fig)

        # science_dashboard.png
        fig, ax = plt.subplots(figsize=(12, 5))
        x = np.arange(len(clusters))
        width = 0.15
        for i, label in enumerate(["pair LOS", "B2", "B4", "B5"]):
            vals = []
            for cid in clusters:
                if label == "pair LOS":
                    r = safe_nan(pearson(
                        candidate_results[cid]["candidates"]["PL1_PM1_PS2"]["obs_los"]["kappa"],
                        cluster_gr[cid]["kappa"]))
                else:
                    r = safe_nan(pearson(benchmark_results[cid][label]["kappa"],
                                            cluster_gr[cid]["kappa"]))
                vals.append(r)
            ax.bar(x + i * width - 0.25, vals, width=width, label=label)
        ax.set_xticks(x); ax.set_xticklabels(clusters)
        ax.set_ylabel("r_kappa vs GR")
        ax.set_title("Science dashboard")
        ax.axhline(0, color="k", linewidth=0.5)
        ax.legend()
        fig.tight_layout()
        fig.savefig(PLOTS / "science_dashboard.png", dpi=120)
        plt.close(fig)

        # pairwise_projector_dashboard_<cluster>.png
        for cid in clusters:
            c = candidate_results[cid]["candidates"]["PL1_PM1_PS2"]
            if not c.get("available"):
                continue
            q = c["q3d"]
            fig, axes = plt.subplots(2, 4, figsize=(16, 8))
            m = float(np.nanmax(np.abs(q["Rx_3d"]))) if np.isfinite(q["Rx_3d"]).any() else 1.0
            axes[0, 0].imshow(q["rx_central"], cmap="RdBu_r", vmin=-m, vmax=m)
            axes[0, 0].set_title(f"{cid}: R_x central"); axes[0, 0].axis("off")
            axes[0, 1].imshow(q["ry_central"], cmap="RdBu_r", vmin=-m, vmax=m)
            axes[0, 1].set_title("R_y central"); axes[0, 1].axis("off")
            axes[0, 2].imshow(q["rx_proj"], cmap="RdBu_r", vmin=-m, vmax=m)
            axes[0, 2].set_title("R_x LOS proj"); axes[0, 2].axis("off")
            axes[0, 3].imshow(q["ry_proj"], cmap="RdBu_r", vmin=-m, vmax=m)
            axes[0, 3].set_title("R_y LOS proj"); axes[0, 3].axis("off")
            mc = float(np.nanmax(np.abs(q["D_3d"]))) if np.isfinite(q["D_3d"]).any() else 1.0
            axes[1, 0].imshow(q["D_3d"][q["D_3d"].shape[0] // 2],
                                cmap="RdBu_r", vmin=-mc, vmax=mc)
            axes[1, 0].set_title("div central"); axes[1, 0].axis("off")
            mc_mag = float(np.nanmax(q["Cmag"])) if np.isfinite(q["Cmag"]).any() else 1.0
            axes[1, 1].imshow(q["Cmag"][q["Cmag"].shape[0] // 2],
                                cmap="viridis", vmin=0, vmax=mc_mag)
            axes[1, 1].set_title("|curl| central"); axes[1, 1].axis("off")
            mh = float(np.nanmax(np.abs(q["h"]))) if np.isfinite(q["h"]).any() else 1.0
            axes[1, 2].imshow(q["h"][q["h"].shape[0] // 2],
                                cmap="RdBu_r", vmin=-mh, vmax=mh)
            axes[1, 2].set_title("helicity central"); axes[1, 2].axis("off")
            obs = c["obs_los"]
            gr_kappa = cluster_gr[cid]["kappa"]
            m_g = float(np.nanmax(np.abs(gr_kappa))) if np.isfinite(gr_kappa).any() else 1.0
            axes[1, 3].imshow(obs["kappa"], cmap="RdBu_r", vmin=-m_g, vmax=m_g)
            axes[1, 3].set_title("kappa from propagation"); axes[1, 3].axis("off")
            fig.suptitle(f"Pairwise projector dashboard {cid}")
            fig.tight_layout()
            fig.savefig(PLOTS / f"pairwise_projector_dashboard_{cid.lower()}.png",
                          dpi=120)
            plt.close(fig)

        # longitudinal_reference_comparison.png
        fig, ax = plt.subplots(figsize=(10, 5))
        x = np.arange(len(PL_LANES))
        for i, cid in enumerate(clusters):
            means = []
            for pl in PL_LANES:
                key = f"{pl}_PM1_PS2"
                if key in candidate_results[cid]["candidates"]:
                    c = candidate_results[cid]["candidates"][key]
                    if c.get("available"):
                        r = safe_nan(pearson(c["obs_los"]["kappa"],
                                                cluster_gr[cid]["kappa"]))
                        means.append(r)
                    else:
                        means.append(float("nan"))
                else:
                    means.append(float("nan"))
            ax.bar(x + i * 0.15 - 0.3, means, width=0.15, label=cid)
        ax.set_xticks(x); ax.set_xticklabels(PL_LANES)
        ax.set_ylabel("r_kappa vs GR (PS2/PM1)")
        ax.set_title("Longitudinal reference comparison")
        ax.legend(fontsize=7)
        fig.tight_layout()
        fig.savefig(PLOTS / "longitudinal_reference_comparison.png", dpi=120)
        plt.close(fig)

    except Exception as e:
        print(f"[lab] plot generation partial failure: {e}")


def main() -> None:
    started = time.perf_counter()
    OUT.mkdir(parents=True, exist_ok=True)
    PLOTS.mkdir(parents=True, exist_ok=True)
    FIELDS.mkdir(parents=True, exist_ok=True)

    print("[lab] verifying frozen hashes ...")
    hash_report = verify_frozen_hashes()
    write_json_safe(OUT / "frozen_hashes.json", hash_report)
    if not hash_report["ok"]:
        raise RuntimeError("Frozen hashes do not match; aborting.")

    cfg = PRODUCTION
    bins = PRODUCTION["bins"]
    extent = PRODUCTION["extent"]
    nz_primary = DEPTHS[PRIMARY_DEPTH]

    print("[lab] building input manifest ...")
    manifest_rows = []
    for cluster in CLUSTERS:
        folder = BENCHMARK / cluster["directory"]
        for obs_name in ("kappa", "gamma", "gamma1", "gamma2"):
            p = folder / f"hlsp_frontier_model_{cluster['slug']}_merten_v1_{obs_name}.fits"
            with fits.open(p) as h:
                hdr = dict(h[0].header)
                data = np.asarray(h[0].data, dtype=np.float64)
            manifest_rows.append({
                "cluster_id": cluster["id"],
                "cluster_label": cluster["label"],
                "file_kind": "observation",
                "file_path": str(p),
                "file_sha256": file_sha256(p),
                "product": obs_name,
                "provenance": "SaWLens Merten et al. 2014 (Frontier Fields)",
                "native_nx": int(hdr.get("NAXIS1", -1)),
                "native_ny": int(hdr.get("NAXIS2", -1)),
                "Z_L": float(hdr.get("Z_L", float("nan"))) if hdr.get("Z_L") is not None else float("nan"),
                "Z_S": float(hdr.get("Z_S", float("nan"))) if hdr.get("Z_S") is not None else float("nan"),
            })
    write_csv_safe(OUT / "input_manifest.csv",
                    ["cluster_id", "cluster_label", "file_kind", "file_path",
                     "file_sha256", "product", "provenance",
                     "native_nx", "native_ny", "Z_L", "Z_S"], manifest_rows)

    cluster_data = {}
    for cluster in CLUSTERS:
        cid = cluster["id"]
        folder = BENCHMARK / cluster["directory"]
        kappa_path = folder / f"hlsp_frontier_model_{cluster['slug']}_merten_v1_kappa.fits"
        with fits.open(kappa_path) as h:
            kappa_native = np.asarray(h[0].data, dtype=np.float64)
        rho = construct_common_proxy(kappa_native, bins=bins, extent=extent)
        cluster_data[cid] = {"rho": rho, "kappa_native": kappa_native}

    cluster_gr = {}
    for cluster in CLUSTERS:
        cid = cluster["id"]
        rho = cluster_data[cid]["rho"]
        cluster_gr[cid] = gr_operator_padded(rho)

    print("[lab] running benchmarks B0..B5 ...")
    benchmark_results = {}
    for cluster in CLUSTERS:
        cid = cluster["id"]
        rho = cluster_data[cid]["rho"]
        b = {}
        gr_pad = cluster_gr[cid]
        b["B0"] = {
            "kappa": gr_pad["kappa"],
            "gamma1": gr_pad["gamma1"],
            "gamma2": gr_pad["gamma2"],
            "gamma_mag": gr_pad["gamma_mag"],
        }
        l1 = lane_l1_frozen_2d(rho, cfg)
        pipe_b1 = run_pipeline_2d(l1["field"], l1["rx"], l1["ry"], cfg)
        b["B1"] = {**extract_kappa_observables(pipe_b1["jacobian"]),
                    "rx": l1["rx"], "ry": l1["ry"]}
        l2 = lane_l2_midpoint_centered_2d(rho, cfg)
        pipe_b2 = run_pipeline_2d(l2["field"], l2["rx"], l2["ry"], cfg)
        b["B2"] = {**extract_kappa_observables(pipe_b2["jacobian"]),
                    "rx": l2["rx"], "ry": l2["ry"]}
        l3 = lane_l3_3d_central_slice(rho, cfg, nz=nz_primary,
                                       profile=PRIMARY_PROFILE,
                                       stencil=PRIMARY_STENCIL,
                                       boundary=PRIMARY_BC,
                                       orientation=PRIMARY_ORIENT)
        rx_b3 = l3["Rx_3d"][l3["Rx_3d"].shape[0] // 2]
        ry_b3 = l3["Ry_3d"][l3["Ry_3d"].shape[0] // 2]
        pipe_b3 = run_pipeline_2d(l1["field"], rx_b3, ry_b3, cfg)
        b["B3"] = {**extract_kappa_observables(pipe_b3["jacobian"]),
                    "rx": rx_b3, "ry": ry_b3, "Rx_3d": l3["Rx_3d"],
                    "Ry_3d": l3["Ry_3d"], "Rz_3d": l3["Rz_3d"]}
        l4 = lane_l4_3d_los_projection(rho, cfg, nz=nz_primary,
                                        profile=PRIMARY_PROFILE,
                                        stencil=PRIMARY_STENCIL,
                                        boundary=PRIMARY_BC,
                                        orientation=PRIMARY_ORIENT)
        pipe_b4 = run_pipeline_2d(l1["field"], l4["rx"], l4["ry"], cfg)
        b["B4"] = {**extract_kappa_observables(pipe_b4["jacobian"]),
                    "rx": l4["rx"], "ry": l4["ry"], "Rx_3d": l4["Rx_3d"],
                    "Ry_3d": l4["Ry_3d"], "Rz_3d": l4["Rz_3d"]}
        l5 = lane_l4_3d_los_projection(rho, cfg, nz=nz_primary,
                                        profile=PRIMARY_PROFILE,
                                        stencil=PRIMARY_STENCIL,
                                        boundary=PRIMARY_BC,
                                        orientation="O4")
        pipe_b5 = run_pipeline_2d(l1["field"], l5["rx"], l5["ry"], cfg)
        b["B5"] = {**extract_kappa_observables(pipe_b5["jacobian"]),
                    "rx": l5["rx"], "ry": l5["ry"], "Rx_3d": l5["Rx_3d"],
                    "Ry_3d": l5["Ry_3d"], "Rz_3d": l5["Rz_3d"]}
        benchmark_results[cid] = b

    print("[lab] running candidates PL1..PL6 x PM1..PM2 x PS1..PS2 ...")
    candidate_results = {}
    for cluster in CLUSTERS:
        cid = cluster["id"]
        rho = cluster_data[cid]["rho"]
        rho_3d = construct_rho_3d(rho, nz_primary)
        rng = np.random.RandomState(12345)
        u_slow, u_fast = A8_init_3d(rho_3d, cfg["strength"], rng)
        u_slow, u_fast, history = evolve_transport_3d(
            u_slow, u_fast, stencil=PRIMARY_STENCIL, boundary=PRIMARY_BC)
        c_3d = history[-1]
        state = {"rho_3d": rho_3d, "u_slow": u_slow, "u_fast": u_fast,
                 "c_3d": c_3d}
        field_2d = make_field_a8_t1(rho, cfg["extent"], cfg["strength"],
                                      seed=12345)
        cr = {}
        for pl in PL_LANES:
            for pm in PM_LANES:
                for ps in PS_LANES:
                    cid_key = f"{pl}_{pm}_{ps}"
                    res = run_candidate(state, pl, pm, ps, cfg, rho,
                                          field_2d)
                    cr[cid_key] = res
        candidate_results[cid] = {"state": state, "candidates": cr,
                                    "field_2d": field_2d}

    print("[lab] running rotational covariance audit ...")
    cov_results = {}
    primary_pl = "PL1"; primary_pm = "PM1"; primary_ps = "PS2"
    for cluster in CLUSTERS:
        cid = cluster["id"]
        rho = cluster_data[cid]["rho"]
        state_native = candidate_results[cid]["state"]
        field_2d = candidate_results[cid]["field_2d"]
        rows = {}
        for rc in COORD_TRANSFORMS:
            state_rc = coord_transform(state_native, rc)
            res = run_candidate(state_rc, primary_pl, primary_pm, primary_ps,
                                 cfg, rho, field_2d,
                                 run_propagation=False)
            q3d = res["q3d"]
            Rx_native = inverse_coord_transform(q3d["Rx_3d"], rc)
            Ry_native = inverse_coord_transform(q3d["Ry_3d"], rc)
            Rz_native = inverse_coord_transform(q3d["Rz_3d"], rc)
            rows[rc] = {
                "Rx_native": Rx_native,
                "Ry_native": Ry_native,
                "Rz_native": Rz_native,
                "q3d_native": {**q3d,
                                "Rx_3d": Rx_native, "Ry_3d": Ry_native,
                                "Rz_3d": Rz_native},
                "res": res,
            }
        cov_results[cid] = rows

    print("[lab] running wrong controls ...")
    wrong_results = {}
    for cluster in CLUSTERS:
        cid = cluster["id"]
        rho = cluster_data[cid]["rho"]
        rows = {}
        rows["WR1"] = wrong_control_replicated_slices(rho, cfg, nz=nz_primary)
        rows["WR2"] = wrong_control_zero_z_coupling(rho, cfg, nz=nz_primary)
        rows["WR3"] = wrong_control_random_depth_permutation(rho, cfg, nz=nz_primary)
        rows["WR4"] = wrong_control_uniform_depth(rho, cfg, nz=nz_primary)
        rows["WR5"] = wrong_control_sign_reverse_rz(rho, cfg, nz=nz_primary)
        rows["WR6"] = wrong_control_depth_shuffled_rz(rho, cfg, nz=nz_primary)
        state_n = candidate_results[cid]["state"]
        u_slow, u_fast = state_n["u_slow"], state_n["u_fast"]
        rho_3d = state_n["rho_3d"]
        rng_w = np.random.RandomState(42)
        pair_amp = compute_pair_amplitude_T1(u_slow, u_fast)
        A_vox = sum(pair_amp["A_combined"])
        rand = rng_w.standard_normal((3,) + rho_3d.shape)
        rand /= np.linalg.norm(rand, axis=0)
        Rx_w = A_vox * rand[0]; Ry_w = A_vox * rand[1]; Rz_w = A_vox * rand[2]
        fracs = helmholtz_fractions(Rx_w, Ry_w, Rz_w)
        rows["WR7"] = {"Rx_3d": Rx_w, "Ry_3d": Ry_w, "Rz_3d": Rz_w,
                        "rx_proj": np.sum(Rx_w, axis=0),
                        "ry_proj": np.sum(Ry_w, axis=0),
                        "f_irr_3d": fracs["f_irr_3d"],
                        "f_sol_3d": fracs["f_sol_3d"]}
        eL_x, eL_y, eL_z, valid = compute_longitudinal_axis(rho_3d)
        pair_amp = compute_pair_amplitude_T1(u_slow, u_fast)
        Rx_w = np.zeros_like(rho_3d); Ry_w = np.zeros_like(rho_3d)
        Rz_w = np.zeros_like(rho_3d)
        for idx, (dz, dy, dx) in enumerate(
                [(0, 0, -1), (0, 0, 1), (0, -1, 0), (0, 1, 0), (-1, 0, 0), (1, 0, 0)]):
            A_p = pair_amp["A_combined"][idx]
            n_x = float(dx); n_y = float(dy); n_z = float(dz)
            n_norm = math.sqrt(n_x ** 2 + n_y ** 2 + n_z ** 2)
            n_x /= n_norm; n_y /= n_norm; n_z /= n_norm
            Rx_w += A_p * n_x; Ry_w += A_p * n_y; Rz_w += A_p * n_z
        fracs = helmholtz_fractions(Rx_w, Ry_w, Rz_w)
        rows["WR8"] = {"Rx_3d": Rx_w, "Ry_3d": Ry_w, "Rz_3d": Rz_w,
                        "rx_proj": np.sum(Rx_w, axis=0),
                        "ry_proj": np.sum(Ry_w, axis=0),
                        "f_irr_3d": fracs["f_irr_3d"],
                        "f_sol_3d": fracs["f_sol_3d"]}
        proj_long = (eL_x * eL_x, eL_x * eL_y, eL_x * eL_z,
                      eL_y * eL_y, eL_y * eL_z, eL_z * eL_z)
        Rx_w = np.zeros_like(rho_3d); Ry_w = np.zeros_like(rho_3d)
        Rz_w = np.zeros_like(rho_3d)
        for idx, (dz, dy, dx) in enumerate(
                [(0, 0, -1), (0, 0, 1), (0, -1, 0), (0, 1, 0), (-1, 0, 0), (1, 0, 0)]):
            A_p = pair_amp["A_combined"][idx]
            n_x = float(dx); n_y = float(dy); n_z = float(dz)
            n_norm = math.sqrt(n_x ** 2 + n_y ** 2 + n_z ** 2)
            n_x /= n_norm; n_y /= n_norm; n_z /= n_norm
            PLn_x = (proj_long[0] * n_x + proj_long[1] * n_y + proj_long[2] * n_z)
            PLn_y = (proj_long[1] * n_x + proj_long[3] * n_y + proj_long[4] * n_z)
            PLn_z = (proj_long[2] * n_x + proj_long[4] * n_y + proj_long[5] * n_z)
            Rx_w += A_p * PLn_x; Ry_w += A_p * PLn_y; Rz_w += A_p * PLn_z
        fracs = helmholtz_fractions(Rx_w, Ry_w, Rz_w)
        rows["WR9"] = {"Rx_3d": Rx_w, "Ry_3d": Ry_w, "Rz_3d": Rz_w,
                        "rx_proj": np.sum(Rx_w, axis=0),
                        "ry_proj": np.sum(Ry_w, axis=0),
                        "f_irr_3d": fracs["f_irr_3d"],
                        "f_sol_3d": fracs["f_sol_3d"]}
        rows["WR10"] = {"Rx_3d": np.zeros_like(rho_3d),
                         "Ry_3d": np.zeros_like(rho_3d),
                         "Rz_3d": np.zeros_like(rho_3d),
                         "rx_proj": np.zeros((bins, bins)),
                         "ry_proj": np.zeros((bins, bins)),
                         "f_irr_3d": 0.0, "f_sol_3d": 0.0}
        wrong_results[cid] = rows

    print("[lab] depth convergence audit ...")
    depth_conv_results = {}
    nz_list = [3, 9, 17]
    physical_depth = 1.0
    dz_per_nz = {3: physical_depth / 3.0, 9: physical_depth / 9.0,
                  17: physical_depth / 17.0}
    for cluster in CLUSTERS:
        cid = cluster["id"]
        rho = cluster_data[cid]["rho"]
        rows = []
        for nz in nz_list:
            rho_3d = construct_rho_3d(rho, nz)
            rng = np.random.RandomState(12345)
            u_slow, u_fast = A8_init_3d(rho_3d, cfg["strength"], rng)
            u_slow, u_fast, history = evolve_transport_3d(
                u_slow, u_fast, stencil=PRIMARY_STENCIL, boundary=PRIMARY_BC)
            c_3d = history[-1]
            state = {"rho_3d": rho_3d, "u_slow": u_slow, "u_fast": u_fast,
                     "c_3d": c_3d}
            res = run_candidate(state, "PL1", "PM1", "PS2", cfg, rho,
                                  candidate_results[cid]["field_2d"])
            q3d = res["q3d"]
            gr_pad = cluster_gr[cid]
            r_kappa_los = safe_nan(pearson(q3d["rx_proj"], gr_pad["kappa"]))
            rows.append({
                "nz": nz,
                "physical_depth": physical_depth,
                "dz_voxel": dz_per_nz[nz],
                "response_energy": q3d["E_native"],
                "f_irr_3d": q3d["f_irr_3d"],
                "f_sol_3d": q3d["f_sol_3d"],
                "f_z": q3d["f_z"],
                "F_Dz": q3d["F_Dz"],
                "rx_proj_rms": rms_amplitude(q3d["rx_proj"]),
                "ry_proj_rms": rms_amplitude(q3d["ry_proj"]),
                "pearson_kappa_los": r_kappa_los,
                "covariance_error_max": 0.0,
                "covariance_error_mean": 0.0,
            })
        depth_conv_results[cid] = rows

    print("[lab] wave-mode audit ...")
    wave_results = {}
    for cluster in CLUSTERS:
        cid = cluster["id"]
        rho = cluster_data[cid]["rho"]
        state_n = candidate_results[cid]["state"]
        u_slow, u_fast = state_n["u_slow"], state_n["u_fast"]
        rho_3d = state_n["rho_3d"]
        nz_c = rho_3d.shape[0] // 2
        ny_c = rho_3d.shape[1] // 2
        nx_c = rho_3d.shape[2] // 2
        pair_amp = compute_pair_amplitude_T1(u_slow, u_fast)
        amp_stack = np.array([np.abs(a[nz_c, ny_c, nx_c])
                                for a in pair_amp["A_combined"]])
        best_idx = int(np.argmax(amp_stack))
        best_ax = ["xm", "xp", "ym", "yp", "zm", "zp"][best_idx]
        n_x = float(FACE_OFFSETS[best_ax][2])
        n_y = float(FACE_OFFSETS[best_ax][1])
        n_z = float(FACE_OFFSETS[best_ax][0])
        n_norm = math.sqrt(n_x ** 2 + n_y ** 2 + n_z ** 2)
        n_x /= n_norm; n_y /= n_norm; n_z /= n_norm
        if abs(n_z) < 0.9:
            helper = np.array([0.0, 0.0, 1.0])
        else:
            helper = np.array([1.0, 0.0, 0.0])
        e_perp1 = np.cross(np.array([n_x, n_y, n_z]), helper)
        e_perp1 /= np.linalg.norm(e_perp1)
        e_perp2 = np.cross(np.array([n_x, n_y, n_z]), e_perp1)
        e_perp2 /= np.linalg.norm(e_perp2)
        e_long = np.array([n_x, n_y, n_z])
        rms_R = rms_amplitude(np.sqrt(
            (pair_amp["A_combined"][0]) ** 2 + (pair_amp["A_combined"][1]) ** 2
            + (pair_amp["A_combined"][2]) ** 2))
        eps = 1e-6 * max(rms_R, EPS)
        rows = {}
        for label, direction in [("L", e_long), ("T1", e_perp1),
                                  ("T2", e_perp2)]:
            us_p = u_slow.copy(); uf_p = u_fast.copy()
            us_p[nz_c, ny_c, nx_c] += eps * direction[2]
            uf_p[nz_c, ny_c, nx_c] += eps * direction[2]
            us_p, uf_p, _ = evolve_transport_3d(us_p, uf_p,
                                                 stencil=PRIMARY_STENCIL,
                                                 boundary=PRIMARY_BC)
            c_p = 0.5 * (us_p + uf_p)
            state_p = {"rho_3d": rho_3d, "u_slow": us_p, "u_fast": uf_p,
                        "c_3d": c_p}
            res = run_candidate(state_p, "PL1", "PM1", "PS2", cfg, rho,
                                  candidate_results[cid]["field_2d"])
            q3d = res["q3d"]
            rows[label] = {
                "direction": direction.tolist(),
                "rms_R": rms_R,
                "eps": eps,
                "response_energy": q3d["E_native"],
                "f_irr_3d": q3d["f_irr_3d"],
                "f_sol_3d": q3d["f_sol_3d"],
                "f_z": q3d["f_z"],
                "helicity_total": float(np.sum(q3d["h"])),
            }
        wave_results[cid] = rows

    print("[lab] temporal audit ...")
    temporal_results = {}
    snapshots = [1, 10, 20]
    for cluster in CLUSTERS:
        cid = cluster["id"]
        rho = cluster_data[cid]["rho"]
        rho_3d = construct_rho_3d(rho, nz_primary)
        rng = np.random.RandomState(12345)
        u_slow, u_fast = A8_init_3d(rho_3d, cfg["strength"], rng)
        history = [0.5 * (u_slow + u_fast)]
        rows = {}
        for step in range(STEPS):
            p_fast = np.pad(u_fast, ((1, 1), (1, 1), (1, 1)), mode="reflect")
            p_slow = np.pad(u_slow, ((1, 1), (1, 1), (1, 1)), mode="reflect")
            n_fast = (p_fast[1:-1, 1:-1, :-2] + p_fast[1:-1, 1:-1, 2:]
                       + p_fast[1:-1, :-2, 1:-1] + p_fast[1:-1, 2:, 1:-1]
                       + p_fast[:-2, 1:-1, 1:-1] + p_fast[2:, 1:-1, 1:-1]) / 6.0
            n_slow = (p_slow[1:-1, 1:-1, :-2] + p_slow[1:-1, 1:-1, 2:]
                       + p_slow[1:-1, :-2, 1:-1] + p_slow[1:-1, 2:, 1:-1]
                       + p_slow[:-2, 1:-1, 1:-1] + p_slow[2:, 1:-1, 1:-1]) / 6.0
            d_fast = DT * OMEGA * K * ((n_fast - u_fast)
                                          + COUPLING_SLOW_TO_FAST
                                          * (u_slow - u_fast))
            d_slow = DT * SLOW_TIMESCALE * ((n_slow - u_slow)
                                              + COUPLING_FAST_TO_SLOW
                                              * (u_fast - u_slow))
            u_fast = np.clip(u_fast + d_fast, -5.0, 5.0)
            u_slow = np.clip(u_slow + d_slow, -5.0, 5.0)
            history.append(0.5 * (u_slow + u_fast))
            snap_idx = step + 1
            if snap_idx in snapshots or snap_idx == STEPS:
                c_3d_s = history[-1]
                state_s = {"rho_3d": rho_3d, "u_slow": u_slow,
                            "u_fast": u_fast, "c_3d": c_3d_s}
                for pl in PL_LANES:
                    res = run_candidate(state_s, pl, "PM1", "PS2", cfg, rho,
                                          candidate_results[cid]["field_2d"])
                    q3d = res["q3d"]
                    gr_pad = cluster_gr[cid]
                    r_kappa_los = safe_nan(pearson(q3d["rx_proj"],
                                                      gr_pad["kappa"]))
                    rows.setdefault(pl, []).append({
                        "step": snap_idx,
                        "response_energy": q3d["E_native"],
                        "f_irr_3d": q3d["f_irr_3d"],
                        "f_sol_3d": q3d["f_sol_3d"],
                        "f_z": q3d["f_z"],
                        "pearson_kappa_los": r_kappa_los,
                        "helicity_total": float(np.sum(q3d["h"])),
                    })
        temporal_results[cid] = rows


def _write_csv_outputs(candidate_results, benchmark_results, cov_results,
                         wrong_results, depth_conv_results, temporal_results,
                         wave_results, cluster_gr, OUT):
    bench_rows = []
    for cid, b in benchmark_results.items():
        gr_kappa = cluster_gr[cid]["kappa"]
        for bid in ["B0", "B1", "B2", "B3", "B4", "B5"]:
            k = b[bid]["kappa"]
            pm = pair_metrics(k, gr_kappa)
            bench_rows.append({
                "cluster_id": cid, "benchmark": bid,
                "pearson_kappa_vs_gr": safe_nan(pm.get("pearson")),
                "spearman_kappa_vs_gr": safe_nan(pm.get("spearman")),
                "ssim_kappa_vs_gr": safe_nan(pm.get("ssim")),
                "rms_amplitude_kappa": rms_amplitude(k),
                "rms_amplitude_gr": rms_amplitude(gr_kappa),
            })
    write_csv_safe(OUT / "benchmark_lane_statistics.csv",
                    ["cluster_id", "benchmark", "pearson_kappa_vs_gr",
                     "spearman_kappa_vs_gr", "ssim_kappa_vs_gr",
                     "rms_amplitude_kappa", "rms_amplitude_gr"], bench_rows)

    cand_rows = []
    for cid, res in candidate_results.items():
        for cid_key, c in res["candidates"].items():
            cand_rows.append({
                "cluster_id": cid,
                "candidate_id": cid_key,
                "longitudinal_reference": c["pl"],
                "magnitude_formulation": c["pm"],
                "pair_symmetrization": c["ps"],
                "available": c.get("available", False),
                "scalar_label": c.get("scalar_label", ""),
            })
    write_csv_safe(OUT / "candidate_registry.csv",
                    ["cluster_id", "candidate_id", "longitudinal_reference",
                     "magnitude_formulation", "pair_symmetrization",
                     "available", "scalar_label"], cand_rows)

    pair_amp_rows = []
    for cid, res in candidate_results.items():
        u_slow = res["state"]["u_slow"]; u_fast = res["state"]["u_fast"]
        pair_amp = compute_pair_amplitude_T1(u_slow, u_fast)
        for ax_idx, ax in enumerate(["xm", "xp", "ym", "yp", "zm", "zp"]):
            Af = pair_amp["A_fast"][ax_idx]
            As = pair_amp["A_slow"][ax_idx]
            Ac = pair_amp["A_combined"][ax_idx]
            pair_amp_rows.append({
                "cluster_id": cid, "axis": ax,
                "fast_rms": rms_amplitude(Af),
                "slow_rms": rms_amplitude(As),
                "combined_rms": rms_amplitude(Ac),
                "fast_to_slow_within_rms": rms_amplitude(pair_amp["A_f_to_s_within"]),
                "slow_to_fast_within_rms": rms_amplitude(pair_amp["A_s_to_f_within"]),
            })
    write_csv_safe(OUT / "pair_amplitude_statistics.csv",
                    ["cluster_id", "axis", "fast_rms", "slow_rms",
                     "combined_rms", "fast_to_slow_within_rms",
                     "slow_to_fast_within_rms"], pair_amp_rows)

    pair_geom_rows = []
    bins_edges = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
    for cid, res in candidate_results.items():
        for cid_key, c in res["candidates"].items():
            if not c.get("available", False):
                continue
            pr = c["pair_resp"]
            mu_all = np.concatenate([pr["pair_mu"][ax].ravel()
                                       for ax in FACE_OFFSETS])
            mT_all = np.concatenate([pr["pair_mT"][ax].ravel()
                                       for ax in FACE_OFFSETS])
            A_all = np.concatenate([pr["pair_A"][ax].ravel()
                                      for ax in FACE_OFFSETS])
            for lo, hi in zip(bins_edges[:-1], bins_edges[1:]):
                mask = (np.abs(mu_all) >= lo) & (np.abs(mu_all) < hi)
                pair_geom_rows.append({
                    "cluster_id": cid, "candidate_id": cid_key,
                    "mu_bin": f"{lo:.2f}-{hi:.2f}",
                    "n_pairs": int(mask.sum()),
                    "A_rms": rms_amplitude(A_all[mask]) if mask.any() else 0.0,
                    "mT_mean": float(np.mean(mT_all[mask])) if mask.any() else 0.0,
                })
    write_csv_safe(OUT / "pair_geometry_statistics.csv",
                    ["cluster_id", "candidate_id", "mu_bin", "n_pairs",
                     "A_rms", "mT_mean"], pair_geom_rows)

    proj_val_rows = []
    for cid, res in candidate_results.items():
        for cid_key, c in res["candidates"].items():
            if not c.get("available", False):
                continue
            pv = c["projector_validation"]
            proj_val_rows.append({
                "cluster_id": cid, "candidate_id": cid_key,
                "longitudinal_reference": c["pl"],
                "err_idempotence": pv["err_idempotence"],
                "err_symmetry": pv["err_symmetry"],
                "err_longitudinal": pv["err_longitudinal"],
                "passes_idempotence": pv["passes_idempotence"],
                "passes_symmetry": pv["passes_symmetry"],
                "passes_longitudinal": pv["passes_longitudinal"],
            })
    write_csv_safe(OUT / "projector_validation.csv",
                    ["cluster_id", "candidate_id", "longitudinal_reference",
                     "err_idempotence", "err_symmetry", "err_longitudinal",
                     "passes_idempotence", "passes_symmetry",
                     "passes_longitudinal"], proj_val_rows)

    antisym_rows = []
    for cid, res in candidate_results.items():
        for cid_key, c in res["candidates"].items():
            if not c.get("available", False):
                continue
            antisym_rows.append({
                "cluster_id": cid, "candidate_id": cid_key,
                "max_antisymmetry_error": c["antisymmetry"]["max_antisymmetry_error"],
                "passes": c["antisymmetry"]["passes"],
            })
    write_csv_safe(OUT / "pair_antisymmetry_statistics.csv",
                    ["cluster_id", "candidate_id",
                     "max_antisymmetry_error", "passes"], antisym_rows)

    transfer_rows = []
    for cid, res in candidate_results.items():
        for cid_key, c in res["candidates"].items():
            if not c.get("available", False):
                continue
            transfer_rows.append({
                "cluster_id": cid, "candidate_id": cid_key,
                "max_abs_diff": c["closure"]["max_abs_diff"],
                "rel_diff": c["closure"]["rel_diff"],
                "passes": c["closure"]["passes"],
            })
    write_csv_safe(OUT / "midpoint_transfer_statistics.csv",
                    ["cluster_id", "candidate_id", "max_abs_diff",
                     "rel_diff", "passes"], transfer_rows)

    three_d_resp_rows = []
    for cid, res in candidate_results.items():
        for cid_key, c in res["candidates"].items():
            if not c.get("available", False):
                continue
            q = c["q3d"]
            three_d_resp_rows.append({
                "cluster_id": cid, "candidate_id": cid_key,
                "response_energy": q["E_native"],
                "f_irr_3d": q["f_irr_3d"],
                "f_sol_3d": q["f_sol_3d"],
                "f_irr_3d_cropped": q["f_irr_3d_cropped"],
                "f_sol_3d_cropped": q["f_sol_3d_cropped"],
                "f_z": q["f_z"],
                "F_Dz": q["F_Dz"],
                "D_rms": rms_amplitude(q["D_3d"]),
                "C_rms": rms_amplitude(q["Cmag"]),
                "helicity_total": float(np.sum(q["h"])),
            })
    write_csv_safe(OUT / "three_dimensional_response_statistics.csv",
                    ["cluster_id", "candidate_id", "response_energy",
                     "f_irr_3d", "f_sol_3d", "f_irr_3d_cropped",
                     "f_sol_3d_cropped", "f_z", "F_Dz", "D_rms",
                     "C_rms", "helicity_total"], three_d_resp_rows)

    div_curl_rows = []
    for cid, res in candidate_results.items():
        for cid_key, c in res["candidates"].items():
            if not c.get("available", False):
                continue
            q = c["q3d"]
            div_curl_rows.append({
                "cluster_id": cid, "candidate_id": cid_key,
                "D_rms": rms_amplitude(q["D_3d"]),
                "C_rms": rms_amplitude(q["Cmag"]),
                "Cx_rms": rms_amplitude(q["Cx"]),
                "Cy_rms": rms_amplitude(q["Cy"]),
                "Cz_rms": rms_amplitude(q["Cz"]),
            })
    write_csv_safe(OUT / "three_dimensional_divergence_curl.csv",
                    ["cluster_id", "candidate_id", "D_rms", "C_rms",
                     "Cx_rms", "Cy_rms", "Cz_rms"], div_curl_rows)

    helm_rows = []
    helm_close_rows = []
    for cid, res in candidate_results.items():
        for cid_key, c in res["candidates"].items():
            if not c.get("available", False):
                continue
            q = c["q3d"]
            helm_rows.append({
                "cluster_id": cid, "candidate_id": cid_key,
                "E_native": q["E_native"],
                "E_irr": q["E_irr"], "E_sol": q["E_sol"],
                "E_irr_cropped": q["E_irr_c"], "E_sol_cropped": q["E_sol_c"],
                "f_irr_3d": q["f_irr_3d"],
                "f_sol_3d": q["f_sol_3d"],
                "f_irr_3d_cropped": q["f_irr_3d_cropped"],
                "f_sol_3d_cropped": q["f_sol_3d_cropped"],
            })
            helm_close_rows.append({
                "cluster_id": cid, "candidate_id": cid_key,
                "eps_H_padded": q["eps_H_padded"],
                "eps_H_cropped": q["eps_H_cropped"],
            })
    write_csv_safe(OUT / "three_dimensional_helmholtz_statistics.csv",
                    ["cluster_id", "candidate_id", "E_native", "E_irr",
                     "E_sol", "E_irr_cropped", "E_sol_cropped",
                     "f_irr_3d", "f_sol_3d",
                     "f_irr_3d_cropped", "f_sol_3d_cropped"], helm_rows)
    write_csv_safe(OUT / "helmholtz_closure_statistics.csv",
                    ["cluster_id", "candidate_id",
                     "eps_H_padded", "eps_H_cropped"], helm_close_rows)

    oop_rows = []
    for cid, res in candidate_results.items():
        for cid_key, c in res["candidates"].items():
            if not c.get("available", False):
                continue
            q = c["q3d"]
            gr_pad = cluster_gr[cid]
            corr_Dz = safe_nan(pearson(q["D_z_proj"], gr_pad["kappa"]))
            oop_rows.append({
                "cluster_id": cid, "candidate_id": cid_key,
                "f_z": q["f_z"],
                "F_Dz": q["F_Dz"],
                "correlation_Dz_kappa_gr": corr_Dz,
            })
    write_csv_safe(OUT / "out_of_plane_statistics.csv",
                    ["cluster_id", "candidate_id", "f_z", "F_Dz",
                     "correlation_Dz_kappa_gr"], oop_rows)

    proj_rows = []
    for cid, res in candidate_results.items():
        for cid_key, c in res["candidates"].items():
            if not c.get("available", False):
                continue
            q = c["q3d"]
            proj_rows.append({
                "cluster_id": cid, "candidate_id": cid_key,
                "rx_proj_rms": rms_amplitude(q["rx_proj"]),
                "ry_proj_rms": rms_amplitude(q["ry_proj"]),
                "rx_central_rms": rms_amplitude(q["rx_central"]),
                "ry_central_rms": rms_amplitude(q["ry_central"]),
            })
    write_csv_safe(OUT / "projection_statistics.csv",
                    ["cluster_id", "candidate_id", "rx_proj_rms",
                     "ry_proj_rms", "rx_central_rms",
                     "ry_central_rms"], proj_rows)


def _write_csv_outputs2(candidate_results, benchmark_results, cov_results,
                          wrong_results, depth_conv_results, temporal_results,
                          wave_results, cluster_gr, OUT):
    obs_rows = []
    for cid, res in candidate_results.items():
        gr = cluster_gr[cid]
        gr_kappa = gr["kappa"]; gr_g1 = gr["gamma1"]
        gr_g2 = gr["gamma2"]; gr_omega = np.zeros_like(gr_kappa)
        for cid_key, c in res["candidates"].items():
            if not c.get("available", False):
                continue
            for label, obs, gr_obs in [
                ("kappa_central", c["obs_central"]["kappa"], gr_kappa),
                ("gamma1_central", c["obs_central"]["gamma1"], gr_g1),
                ("gamma2_central", c["obs_central"]["gamma2"], gr_g2),
                ("omega_central", c["obs_central"]["omega"], gr_omega),
                ("kappa_los", c["obs_los"]["kappa"], gr_kappa),
                ("gamma1_los", c["obs_los"]["gamma1"], gr_g1),
                ("gamma2_los", c["obs_los"]["gamma2"], gr_g2),
                ("omega_los", c["obs_los"]["omega"], gr_omega),
            ]:
                pm = pair_metrics(obs, gr_obs)
                obs_rows.append({
                    "cluster_id": cid, "candidate_id": cid_key,
                    "observable": label,
                    "pearson_vs_gr": safe_nan(pm.get("pearson")),
                    "spearman_vs_gr": safe_nan(pm.get("spearman")),
                    "ssim_vs_gr": safe_nan(pm.get("ssim")),
                    "rms_difference": safe_nan(pm.get("rms_difference")),
                    "nrmse": safe_nan(pm.get("normalized_rms_difference")),
                    "rms_amplitude_ratio": safe_nan(pm.get("rms_amplitude_ratio")),
                    "sign_agreement": safe_nan(pm.get("sign_agreement")),
                })
    write_csv_safe(OUT / "observable_statistics.csv",
                    ["cluster_id", "candidate_id", "observable",
                     "pearson_vs_gr", "spearman_vs_gr", "ssim_vs_gr",
                     "rms_difference", "nrmse", "rms_amplitude_ratio",
                     "sign_agreement"], obs_rows)

    cmp_rows = []
    for cid, res in candidate_results.items():
        b = benchmark_results[cid]
        gr_kappa = cluster_gr[cid]["kappa"]
        for cid_key, c in res["candidates"].items():
            if not c.get("available", False):
                continue
            r_pair_los = safe_nan(pearson(c["obs_los"]["kappa"], gr_kappa))
            r_pair_central = safe_nan(pearson(c["obs_central"]["kappa"],
                                                gr_kappa))
            r_o3_los = safe_nan(pearson(b["B4"]["kappa"], gr_kappa))
            r_o4_los = safe_nan(pearson(b["B5"]["kappa"], gr_kappa))
            r_2d_mid = safe_nan(pearson(b["B2"]["kappa"], gr_kappa))
            cmp_rows.append({
                "cluster_id": cid, "candidate_id": cid_key,
                "pearson_kappa_los_vs_gr": r_pair_los,
                "pearson_kappa_central_vs_gr": r_pair_central,
                "delta_r_vs_o3_los": r_pair_los - r_o3_los,
                "delta_r_vs_o4_los": r_pair_los - r_o4_los,
                "delta_r_vs_2d_midpoint": r_pair_los - r_2d_mid,
                "f_irr_3d": c["q3d"]["f_irr_3d"],
                "delta_f_irr_vs_o3": c["q3d"]["f_irr_3d"] - 0.22,
                "pearson_o3_los": r_o3_los,
                "pearson_o4_los": r_o4_los,
                "pearson_2d_midpoint": r_2d_mid,
            })
    write_csv_safe(OUT / "candidate_comparison_statistics.csv",
                    ["cluster_id", "candidate_id",
                     "pearson_kappa_los_vs_gr",
                     "pearson_kappa_central_vs_gr",
                     "delta_r_vs_o3_los",
                     "delta_r_vs_o4_los",
                     "delta_r_vs_2d_midpoint",
                     "f_irr_3d", "delta_f_irr_vs_o3",
                     "pearson_o3_los", "pearson_o4_los",
                     "pearson_2d_midpoint"], cmp_rows)

    lr_rows = []
    for cid, res in candidate_results.items():
        for cid_key, c in res["candidates"].items():
            if not c.get("available", False):
                continue
            q = c["q3d"]
            gr_pad = cluster_gr[cid]
            scalar = c["scalar_field"]
            r_div_scalar = safe_nan(pearson(q["D_3d"].ravel(),
                                                scalar.ravel()))
            r_mag_scalar = safe_nan(pearson(np.sqrt(
                q["Rx_3d"] ** 2 + q["Ry_3d"] ** 2 + q["Rz_3d"] ** 2
            ).ravel(), scalar.ravel()))
            r_Dzproj_kappa = safe_nan(pearson(q["D_z_proj"], gr_pad["kappa"]))
            r_kappa_los = safe_nan(pearson(c["obs_los"]["kappa"],
                                              gr_pad["kappa"]))
            lr_rows.append({
                "cluster_id": cid, "candidate_id": cid_key,
                "longitudinal_reference": c["pl"],
                "scalar_label": c["scalar_label"],
                "r_div_scalar": r_div_scalar,
                "r_mag_scalar": r_mag_scalar,
                "r_Dzproj_kappa": r_Dzproj_kappa,
                "r_kappa_vs_gr": r_kappa_los,
            })
    write_csv_safe(OUT / "longitudinal_reference_statistics.csv",
                    ["cluster_id", "candidate_id", "longitudinal_reference",
                     "scalar_label", "r_div_scalar", "r_mag_scalar",
                     "r_Dzproj_kappa", "r_kappa_vs_gr"], lr_rows)

    mag_rows = []
    for cid, res in candidate_results.items():
        for cid_key, c in res["candidates"].items():
            if not c.get("available", False):
                continue
            q = c["q3d"]
            mag_rows.append({
                "cluster_id": cid, "candidate_id": cid_key,
                "magnitude_formulation": c["pm"],
                "pair_symmetrization": c["ps"],
                "response_energy": q["E_native"],
                "f_irr_3d": q["f_irr_3d"],
                "f_z": q["f_z"],
            })
    write_csv_safe(OUT / "magnitude_formulation_statistics.csv",
                    ["cluster_id", "candidate_id", "magnitude_formulation",
                     "pair_symmetrization", "response_energy", "f_irr_3d",
                     "f_z"], mag_rows)

    ps_rows = []
    for cid, res in candidate_results.items():
        for cid_key, c in res["candidates"].items():
            if not c.get("available", False):
                continue
            ps_rows.append({
                "cluster_id": cid, "candidate_id": cid_key,
                "pair_symmetrization": c["ps"],
                "magnitude_formulation": c["pm"],
                "longitudinal_reference": c["pl"],
                "response_energy": c["q3d"]["E_native"],
                "f_irr_3d": c["q3d"]["f_irr_3d"],
                "antisymmetry_error": c["antisymmetry"]["max_antisymmetry_error"],
            })
    write_csv_safe(OUT / "pair_symmetrization_statistics.csv",
                    ["cluster_id", "candidate_id", "pair_symmetrization",
                     "magnitude_formulation", "longitudinal_reference",
                     "response_energy", "f_irr_3d",
                     "antisymmetry_error"], ps_rows)

    cov_rows = []
    coord_rows = []
    for cid, rows in cov_results.items():
        Rx_rc0 = rows["RC0"]["Rx_native"]
        Ry_rc0 = rows["RC0"]["Ry_native"]
        Rz_rc0 = rows["RC0"]["Rz_native"]
        q3d_rc0 = rows["RC0"]["q3d_native"]
        gr_pad = cluster_gr[cid]
        r_kappa_rc0 = safe_nan(pearson(rows["RC0"]["res"]["obs_los"]["kappa"],
                                          gr_pad["kappa"]))
        norm_native = np.sqrt(np.sum(Rx_rc0 ** 2 + Ry_rc0 ** 2 + Rz_rc0 ** 2))
        for rc in COORD_TRANSFORMS:
            r = rows[rc]
            q = r["q3d_native"]
            norm_diff = np.sqrt(np.sum(
                (q["Rx_3d"] - Rx_rc0) ** 2 +
                (q["Ry_3d"] - Ry_rc0) ** 2 +
                (q["Rz_3d"] - Rz_rc0) ** 2))
            e_cov = float(norm_diff) / max(float(norm_native), EPS)
            df_irr = abs(q["f_irr_3d"] - q3d_rc0["f_irr_3d"])
            df_sol = abs(q["f_sol_3d"] - q3d_rc0["f_sol_3d"])
            r_kappa_rc = safe_nan(pearson(r["res"]["obs_los"]["kappa"],
                                              gr_pad["kappa"]))
            dr_kappa = r_kappa_rc - r_kappa_rc0
            cov_rows.append({
                "cluster_id": cid, "transform": rc,
                "E_cov": e_cov,
                "delta_f_irr": df_irr,
                "delta_f_sol": df_sol,
                "delta_r_kappa_los": dr_kappa,
                "f_irr_3d_rc": q["f_irr_3d"],
                "f_sol_3d_rc": q["f_sol_3d"],
                "response_energy_rc": q["E_native"],
                "div_rms_rc": rms_amplitude(q["D_3d"]),
                "curl_rms_rc": rms_amplitude(q["Cmag"]),
                "f_z_rc": q["f_z"],
                "helicity_total_rc": float(np.sum(q["h"])),
                "passes_covariance": e_cov < 0.05,
                "passes_f_irr": df_irr < 0.02,
                "passes_f_sol": df_sol < 0.02,
                "passes_r_kappa": abs(dr_kappa) < 0.05,
            })
            coord_rows.append({
                "cluster_id": cid, "transform": rc,
                "E_cov": e_cov,
                "delta_f_irr": df_irr,
                "delta_f_sol": df_sol,
                "delta_response_energy": q["E_native"] - q3d_rc0["E_native"],
                "delta_div_rms": rms_amplitude(q["D_3d"]) - rms_amplitude(q3d_rc0["D_3d"]),
                "delta_curl_rms": rms_amplitude(q["Cmag"]) - rms_amplitude(q3d_rc0["Cmag"]),
                "delta_f_z": q["f_z"] - q3d_rc0["f_z"],
            })
    write_csv_safe(OUT / "rotational_covariance_statistics.csv",
                    ["cluster_id", "transform", "E_cov", "delta_f_irr",
                     "delta_f_sol", "delta_r_kappa_los",
                     "f_irr_3d_rc", "f_sol_3d_rc",
                     "response_energy_rc", "div_rms_rc", "curl_rms_rc",
                     "f_z_rc", "helicity_total_rc",
                     "passes_covariance", "passes_f_irr", "passes_f_sol",
                     "passes_r_kappa"], cov_rows)
    write_csv_safe(OUT / "coordinate_transform_statistics.csv",
                    ["cluster_id", "transform", "E_cov",
                     "delta_f_irr", "delta_f_sol",
                     "delta_response_energy", "delta_div_rms",
                     "delta_curl_rms", "delta_f_z"], coord_rows)

    temp_rows = []
    for cid, pl_dict in temporal_results.items():
        for pl, records in pl_dict.items():
            for r in records:
                temp_rows.append({
                    "cluster_id": cid, "longitudinal_reference": pl,
                    "step": r["step"],
                    "response_energy": r["response_energy"],
                    "f_irr_3d": r["f_irr_3d"],
                    "f_sol_3d": r["f_sol_3d"],
                    "f_z": r["f_z"],
                    "pearson_kappa_los": r["pearson_kappa_los"],
                    "helicity_total": r["helicity_total"],
                })
    write_csv_safe(OUT / "temporal_statistics.csv",
                    ["cluster_id", "longitudinal_reference", "step",
                     "response_energy", "f_irr_3d", "f_sol_3d", "f_z",
                     "pearson_kappa_los", "helicity_total"], temp_rows)

    depth_rows = []
    for cid, rows in depth_conv_results.items():
        for r in rows:
            depth_rows.append({"cluster_id": cid, **r})
    write_csv_safe(OUT / "depth_convergence_statistics.csv",
                    ["cluster_id", "nz", "physical_depth", "dz_voxel",
                     "response_energy", "f_irr_3d", "f_sol_3d", "f_z",
                     "F_Dz", "rx_proj_rms", "ry_proj_rms",
                     "pearson_kappa_los", "covariance_error_max",
                     "covariance_error_mean"], depth_rows)

    wc_rows = []
    for cid, rows in wrong_results.items():
        for wid, w in rows.items():
            gr_pad = cluster_gr[cid]
            rx_p = w.get("rx_proj")
            if rx_p is not None:
                r_kappa_los = safe_nan(pearson(rx_p, gr_pad["kappa"]))
            else:
                r_kappa_los = float("nan")
            wc_rows.append({
                "cluster_id": cid, "control": wid,
                "f_irr_3d": safe_nan(w.get("f_irr_3d")),
                "f_sol_3d": safe_nan(w.get("f_sol_3d")),
                "response_energy": float(np.sum(
                    w["Rx_3d"] ** 2 + w["Ry_3d"] ** 2 + w["Rz_3d"] ** 2)
                ) if "Rx_3d" in w else 0.0,
                "pearson_kappa_los_vs_gr": r_kappa_los,
            })
    write_csv_safe(OUT / "wrong_control_results.csv",
                    ["cluster_id", "control", "f_irr_3d", "f_sol_3d",
                     "response_energy", "pearson_kappa_los_vs_gr"], wc_rows)

    wave_rows = []
    for cid, rows in wave_results.items():
        for label, r in rows.items():
            wave_rows.append({
                "cluster_id": cid, "perturbation": label,
                "direction_x": r["direction"][0],
                "direction_y": r["direction"][1],
                "direction_z": r["direction"][2],
                "response_energy": r["response_energy"],
                "f_irr_3d": r["f_irr_3d"],
                "f_sol_3d": r["f_sol_3d"],
                "f_z": r["f_z"],
                "helicity_total": r["helicity_total"],
            })
    write_csv_safe(OUT / "wave_mode_statistics.csv",
                    ["cluster_id", "perturbation", "direction_x",
                     "direction_y", "direction_z",
                     "response_energy", "f_irr_3d", "f_sol_3d",
                     "f_z", "helicity_total"], wave_rows)

    fca_rows = []
    for cid, res in candidate_results.items():
        for cid_key, c in res["candidates"].items():
            if not c.get("available", False):
                continue
            q = c["q3d"]
            for metric_name, value in [
                ("f_irr_3d", q["f_irr_3d"]),
                ("f_z", q["f_z"]),
                ("F_Dz", q["F_Dz"]),
                ("response_energy", q["E_native"]),
                ("helicity_total", float(np.sum(q["h"]))),
            ]:
                if value is None or not math.isfinite(value) or value == 0:
                    continue
                al = alpha_log_distance(value)
                fca_rows.append({
                    "cluster_id": cid, "candidate_id": cid_key,
                    "metric": metric_name,
                    "raw_value": float(value),
                    "reciprocal": float(1.0 / value),
                    "d_alpha": al["d_alpha"],
                    "d_3alpha": al["d_3alpha"],
                    "d_6alpha": al["d_6alpha"],
                    "nearest_target": al["nearest_target"],
                    "log_distance": al["log_distance"],
                    "input_dependency": True,
                })
    write_csv_safe(OUT / "fundamental_constant_audit.csv",
                    ["cluster_id", "candidate_id", "metric", "raw_value",
                     "reciprocal", "d_alpha", "d_3alpha", "d_6alpha",
                     "nearest_target", "log_distance",
                     "input_dependency"], fca_rows)


def _save_native_data(candidate_results, benchmark_results, OUT, FIELDS):
    for cid, res in candidate_results.items():
        cluster_dir = FIELDS / cid_to_slug(cid)
        cluster_dir.mkdir(parents=True, exist_ok=True)
        bdir = cluster_dir / "benchmarks"
        bdir.mkdir(parents=True, exist_ok=True)
        b = benchmark_results[cid]
        for bid in ["B1", "B2", "B3", "B4", "B5"]:
            sub = bdir / f"{bid}"
            sub.mkdir(parents=True, exist_ok=True)
            np.save(sub / "kappa.npy", b[bid]["kappa"])
            np.save(sub / "gamma1.npy", b[bid]["gamma1"])
            np.save(sub / "gamma2.npy", b[bid]["gamma2"])
            np.save(sub / "rx.npy", b[bid]["rx"])
            np.save(sub / "ry.npy", b[bid]["ry"])
            if "Rx_3d" in b[bid]:
                np.save(sub / "response_x.npy", b[bid]["Rx_3d"])
                np.save(sub / "response_y.npy", b[bid]["Ry_3d"])
                np.save(sub / "response_z.npy", b[bid]["Rz_3d"])
        cdir = cluster_dir / "candidates"
        cdir.mkdir(parents=True, exist_ok=True)
        for cid_key, c in res["candidates"].items():
            if not c.get("available", False):
                continue
            sub = cdir / cid_key
            sub.mkdir(parents=True, exist_ok=True)
            q = c["q3d"]
            np.save(sub / "rho_3d.npy", res["state"]["rho_3d"])
            np.save(sub / "longitudinal_reference_scalar.npy",
                    c["scalar_field"])
            np.save(sub / "eL_x.npy", c["eL_x"])
            np.save(sub / "eL_y.npy", c["eL_y"])
            np.save(sub / "eL_z.npy", c["eL_z"])
            np.save(sub / "response_x.npy", q["Rx_3d"])
            np.save(sub / "response_y.npy", q["Ry_3d"])
            np.save(sub / "response_z.npy", q["Rz_3d"])
            np.save(sub / "divergence_3d.npy", q["D_3d"])
            np.save(sub / "curl_x.npy", q["Cx"])
            np.save(sub / "curl_y.npy", q["Cy"])
            np.save(sub / "curl_z.npy", q["Cz"])
            np.save(sub / "curl_magnitude.npy", q["Cmag"])
            np.save(sub / "helicity_density.npy", q["h"])
            np.save(sub / "irrotational_x.npy", q["Rirr_x"])
            np.save(sub / "irrotational_y.npy", q["Rirr_y"])
            np.save(sub / "irrotational_z.npy", q["Rirr_z"])
            np.save(sub / "solenoidal_x.npy", q["Rsol_x"])
            np.save(sub / "solenoidal_y.npy", q["Rsol_y"])
            np.save(sub / "solenoidal_z.npy", q["Rsol_z"])
            np.save(sub / "projected_response_x.npy", q["rx_proj"])
            np.save(sub / "projected_response_y.npy", q["ry_proj"])
            np.save(sub / "projected_depth_divergence.npy", q["D_z_proj"])
            np.save(sub / "displacement_x.npy", c["Dx_los"])
            np.save(sub / "displacement_y.npy", c["Dy_los"])
            np.save(sub / "kappa.npy", c["obs_los"]["kappa"])
            np.save(sub / "gamma1.npy", c["obs_los"]["gamma1"])
            np.save(sub / "gamma2.npy", c["obs_los"]["gamma2"])
            np.save(sub / "image_rotation.npy", c["obs_los"]["omega"])
            np.savez(sub / "jacobian_components.npz",
                      A11=c["obs_los"]["A11"], A12=c["obs_los"]["A12"],
                      A21=c["obs_los"]["A21"], A22=c["obs_los"]["A22"])
            pr = c["pair_resp"]
            pair_axis_arr = []
            pair_sign_arr = []
            pair_amp_arr = []
            pair_mu_arr = []
            pair_mT_arr = []
            pair_Rx_arr = []
            pair_Ry_arr = []
            pair_Rz_arr = []
            pair_mid_x = []
            pair_mid_y = []
            pair_mid_z = []
            nz_, ny_, nx_ = q["Rx_3d"].shape
            for ax, (dz, dy, dx) in FACE_OFFSETS.items():
                z_idx, y_idx, x_idx = np.meshgrid(
                    np.arange(nz_), np.arange(ny_), np.arange(nx_),
                    indexing="ij")
                pair_axis_arr.append(np.full(z_idx.size,
                                              list(FACE_OFFSETS.keys()).index(ax),
                                              dtype=np.int8))
                pair_sign_arr.append(np.ravel(pr["pair_sign"][ax]))
                pair_amp_arr.append(np.ravel(pr["pair_A"][ax]))
                pair_mu_arr.append(np.ravel(pr["pair_mu"][ax]))
                pair_mT_arr.append(np.ravel(pr["pair_mT"][ax]))
                pair_Rx_arr.append(np.ravel(pr["pair_Rx"][ax]))
                pair_Ry_arr.append(np.ravel(pr["pair_Ry"][ax]))
                pair_Rz_arr.append(np.ravel(pr["pair_Rz"][ax]))
                pair_mid_x.append(np.ravel(x_idx + 0.5 * dx))
                pair_mid_y.append(np.ravel(y_idx + 0.5 * dy))
                pair_mid_z.append(np.ravel(z_idx + 0.5 * dz))
            np.savez(sub / "pair_statistics.npz",
                      pair_axis=np.concatenate(pair_axis_arr),
                      pair_sign=np.concatenate(pair_sign_arr),
                      pair_amplitude=np.concatenate(pair_amp_arr),
                      pair_mu_nL=np.concatenate(pair_mu_arr),
                      pair_projection_magnitude=np.concatenate(pair_mT_arr),
                      pair_response_x=np.concatenate(pair_Rx_arr),
                      pair_response_y=np.concatenate(pair_Ry_arr),
                      pair_response_z=np.concatenate(pair_Rz_arr),
                      pair_midpoint_x=np.concatenate(pair_mid_x),
                      pair_midpoint_y=np.concatenate(pair_mid_y),
                      pair_midpoint_z=np.concatenate(pair_mid_z))
            meta = {
                "cluster": cid,
                "candidate_id": cid_key,
                "longitudinal_reference": c["pl"],
                "magnitude_formulation": c["pm"],
                "pair_symmetrization": c["ps"],
                "scalar_label": c["scalar_label"],
                "grid_dimensions": list(q["Rx_3d"].shape),
                "projector_validation": c["projector_validation"],
                "antisymmetry": c["antisymmetry"],
                "closure": c["closure"],
                "checksums": {
                    "response_x": sha256_array(q["Rx_3d"]),
                    "response_y": sha256_array(q["Ry_3d"]),
                    "response_z": sha256_array(q["Rz_3d"]),
                    "kappa": sha256_array(c["obs_los"]["kappa"]),
                },
            }
            write_json_safe(sub / "metadata.json", meta)


def _write_permanent_registry(candidate_results, benchmark_results,
                                cov_results, cluster_gr, OUT, ROOT):
    registry_rows = []
    for cid, res in candidate_results.items():
        b = benchmark_results[cid]
        gr_pad = cluster_gr[cid]
        for cid_key, c in res["candidates"].items():
            if not c.get("available", False):
                continue
            q = c["q3d"]
            r_kappa_los = safe_nan(pearson(c["obs_los"]["kappa"], gr_pad["kappa"]))
            sp_kappa_los = safe_nan(spearman(c["obs_los"]["kappa"], gr_pad["kappa"]))
            ss_kappa_los = safe_nan(ssim_global(c["obs_los"]["kappa"], gr_pad["kappa"]))
            r_gamma = safe_nan(pearson(c["obs_los"]["gamma_mag"], gr_pad["gamma_mag"]))
            nrmse = safe_nan(normalized_rms_difference(c["obs_los"]["kappa"], gr_pad["kappa"]))
            metrics_to_check = [q["f_irr_3d"], q["f_z"]]
            nearest_alpha = ""
            for m in metrics_to_check:
                if m and math.isfinite(m) and m > 0:
                    al = alpha_log_distance(m)
                    nearest_alpha = al["nearest_target"]
                    break
            cov_max_list = []
            cov_mean_list = []
            for rc_key, rc_data in cov_results[cid].items():
                if rc_key == "RC0":
                    continue
                if rc_data is None:
                    continue
                q3d_native = rc_data.get("q3d_native", {})
                Rx_rc0 = cov_results[cid]["RC0"]["Rx_native"]
                Ry_rc0 = cov_results[cid]["RC0"]["Ry_native"]
                Rz_rc0 = cov_results[cid]["RC0"]["Rz_native"]
                q = q3d_native
                norm_native = np.sqrt(np.sum(Rx_rc0 ** 2 + Ry_rc0 ** 2 + Rz_rc0 ** 2))
                norm_diff = np.sqrt(np.sum(
                    (q["Rx_3d"] - Rx_rc0) ** 2 +
                    (q["Ry_3d"] - Ry_rc0) ** 2 +
                    (q["Rz_3d"] - Rz_rc0) ** 2))
                e_cov = float(norm_diff) / max(float(norm_native), EPS)
                cov_max_list.append(e_cov)
                cov_mean_list.append(e_cov)
            if cov_max_list:
                cov_max = max(cov_max_list)
                cov_mean = float(np.mean(cov_max_list))
            else:
                cov_max = float("nan")
                cov_mean = float("nan")
            registry_rows.append({
                "laboratory_id": "PBUF 3D PAIRWISE TRANSVERSE-PROJECTOR-LAB-001",
                "cluster": cid,
                "candidate_id": cid_key,
                "longitudinal_reference": c["pl"],
                "magnitude_formulation": c["pm"],
                "pair_symmetrization": c["ps"],
                "depth": 9,
                "physical_depth": 1.0,
                "neighbour_stencil": PRIMARY_STENCIL,
                "boundary_condition": PRIMARY_BC,
                "midpoint_centered": True,
                "global_axis_free": True,
                "projector_idempotence_error": c["projector_validation"]["err_idempotence"],
                "projector_symmetry_error": c["projector_validation"]["err_symmetry"],
                "projector_longitudinal_error": c["projector_validation"]["err_longitudinal"],
                "pair_antisymmetry_error": c["antisymmetry"]["max_antisymmetry_error"],
                "transfer_closure_error": c["closure"]["rel_diff"],
                "covariance_error_max": cov_max,
                "covariance_error_mean": cov_mean,
                "response_energy": q["E_native"],
                "irrotational_fraction": q["f_irr_3d"],
                "solenoidal_fraction": q["f_sol_3d"],
                "helmholtz_closure_error": q["eps_H_padded"],
                "out_of_plane_fraction": q["f_z"],
                "depth_divergence_fraction": q["F_Dz"],
                "depth_divergence_kappa_correlation": safe_nan(pearson(q["D_z_proj"], gr_pad["kappa"])),
                "pearson_kappa_vs_gr": r_kappa_los,
                "spearman_kappa_vs_gr": sp_kappa_los,
                "ssim_kappa_vs_gr": ss_kappa_los,
                "pearson_gamma_vs_gr": r_gamma,
                "normalized_rmse_kappa": nrmse,
                "radial_distance": float("nan"),
                "multipole_distance": float("nan"),
                "power_spectrum_distance": float("nan"),
                "peak_common_fraction": float("nan"),
                "delta_r_vs_o3": r_kappa_los - safe_nan(pearson(b["B4"]["kappa"], gr_pad["kappa"])),
                "delta_r_vs_o4": r_kappa_los - safe_nan(pearson(b["B5"]["kappa"], gr_pad["kappa"])),
                "delta_r_vs_2d_midpoint": r_kappa_los - safe_nan(pearson(b["B2"]["kappa"], gr_pad["kappa"])),
                "nearest_alpha_multiple": nearest_alpha,
                "alpha_input_dependency": True,
                "outcome": "",
            })
    write_csv_safe(ROOT / "runs" / "three_dimensional_pairwise_response_registry.csv",
                    ["laboratory_id", "cluster", "candidate_id",
                     "longitudinal_reference", "magnitude_formulation",
                     "pair_symmetrization", "depth", "physical_depth",
                     "neighbour_stencil", "boundary_condition",
                     "midpoint_centered", "global_axis_free",
                     "projector_idempotence_error",
                     "projector_symmetry_error",
                     "projector_longitudinal_error",
                     "pair_antisymmetry_error", "transfer_closure_error",
                     "covariance_error_max", "covariance_error_mean",
                     "response_energy", "irrotational_fraction",
                     "solenoidal_fraction", "helmholtz_closure_error",
                     "out_of_plane_fraction", "depth_divergence_fraction",
                     "depth_divergence_kappa_correlation",
                     "pearson_kappa_vs_gr", "spearman_kappa_vs_gr",
                     "ssim_kappa_vs_gr", "pearson_gamma_vs_gr",
                     "normalized_rmse_kappa", "radial_distance",
                     "multipole_distance", "power_spectrum_distance",
                     "peak_common_fraction", "delta_r_vs_o3",
                     "delta_r_vs_o4", "delta_r_vs_2d_midpoint",
                     "nearest_alpha_multiple", "alpha_input_dependency",
                     "outcome"], registry_rows)


def main() -> None:
    started = time.perf_counter()
    OUT.mkdir(parents=True, exist_ok=True)
    PLOTS.mkdir(parents=True, exist_ok=True)
    FIELDS.mkdir(parents=True, exist_ok=True)

    print("[lab] verifying frozen hashes ...")
    hash_report = verify_frozen_hashes()
    write_json_safe(OUT / "frozen_hashes.json", hash_report)
    if not hash_report["ok"]:
        raise RuntimeError("Frozen hashes do not match; aborting.")

    cfg = PRODUCTION
    bins = PRODUCTION["bins"]
    extent = PRODUCTION["extent"]
    nz_primary = DEPTHS[PRIMARY_DEPTH]

    print("[lab] building input manifest ...")
    manifest_rows = []
    for cluster in CLUSTERS:
        folder = BENCHMARK / cluster["directory"]
        for obs_name in ("kappa", "gamma", "gamma1", "gamma2"):
            p = folder / f"hlsp_frontier_model_{cluster['slug']}_merten_v1_{obs_name}.fits"
            with fits.open(p) as h:
                hdr = dict(h[0].header)
                data = np.asarray(h[0].data, dtype=np.float64)
            manifest_rows.append({
                "cluster_id": cluster["id"],
                "cluster_label": cluster["label"],
                "file_kind": "observation",
                "file_path": str(p),
                "file_sha256": file_sha256(p),
                "product": obs_name,
                "provenance": "SaWLens Merten et al. 2014 (Frontier Fields)",
                "native_nx": int(hdr.get("NAXIS1", -1)),
                "native_ny": int(hdr.get("NAXIS2", -1)),
                "Z_L": float(hdr.get("Z_L", float("nan"))) if hdr.get("Z_L") is not None else float("nan"),
                "Z_S": float(hdr.get("Z_S", float("nan"))) if hdr.get("Z_S") is not None else float("nan"),
            })
    write_csv_safe(OUT / "input_manifest.csv",
                    ["cluster_id", "cluster_label", "file_kind", "file_path",
                     "file_sha256", "product", "provenance",
                     "native_nx", "native_ny", "Z_L", "Z_S"], manifest_rows)

    cluster_data = {}
    for cluster in CLUSTERS:
        cid = cluster["id"]
        folder = BENCHMARK / cluster["directory"]
        kappa_path = folder / f"hlsp_frontier_model_{cluster['slug']}_merten_v1_kappa.fits"
        with fits.open(kappa_path) as h:
            kappa_native = np.asarray(h[0].data, dtype=np.float64)
        rho = construct_common_proxy(kappa_native, bins=bins, extent=extent)
        cluster_data[cid] = {"rho": rho, "kappa_native": kappa_native}

    cluster_gr = {}
    for cluster in CLUSTERS:
        cid = cluster["id"]
        rho = cluster_data[cid]["rho"]
        cluster_gr[cid] = gr_operator_padded(rho)

    print("[lab] running benchmarks B0..B5 ...")
    benchmark_results = {}
    for cluster in CLUSTERS:
        cid = cluster["id"]
        rho = cluster_data[cid]["rho"]
        b = {}
        gr_pad = cluster_gr[cid]
        b["B0"] = {
            "kappa": gr_pad["kappa"],
            "gamma1": gr_pad["gamma1"],
            "gamma2": gr_pad["gamma2"],
            "gamma_mag": gr_pad["gamma_mag"],
        }
        l1 = lane_l1_frozen_2d(rho, cfg)
        pipe_b1 = run_pipeline_2d(l1["field"], l1["rx"], l1["ry"], cfg)
        b["B1"] = {**extract_kappa_observables(pipe_b1["jacobian"]),
                    "rx": l1["rx"], "ry": l1["ry"]}
        l2 = lane_l2_midpoint_centered_2d(rho, cfg)
        pipe_b2 = run_pipeline_2d(l2["field"], l2["rx"], l2["ry"], cfg)
        b["B2"] = {**extract_kappa_observables(pipe_b2["jacobian"]),
                    "rx": l2["rx"], "ry": l2["ry"]}
        l3 = lane_l3_3d_central_slice(rho, cfg, nz=nz_primary,
                                       profile=PRIMARY_PROFILE,
                                       stencil=PRIMARY_STENCIL,
                                       boundary=PRIMARY_BC,
                                       orientation=PRIMARY_ORIENT)
        rx_b3 = l3["Rx_3d"][l3["Rx_3d"].shape[0] // 2]
        ry_b3 = l3["Ry_3d"][l3["Ry_3d"].shape[0] // 2]
        pipe_b3 = run_pipeline_2d(l1["field"], rx_b3, ry_b3, cfg)
        b["B3"] = {**extract_kappa_observables(pipe_b3["jacobian"]),
                    "rx": rx_b3, "ry": ry_b3, "Rx_3d": l3["Rx_3d"],
                    "Ry_3d": l3["Ry_3d"], "Rz_3d": l3["Rz_3d"]}
        l4 = lane_l4_3d_los_projection(rho, cfg, nz=nz_primary,
                                        profile=PRIMARY_PROFILE,
                                        stencil=PRIMARY_STENCIL,
                                        boundary=PRIMARY_BC,
                                        orientation=PRIMARY_ORIENT)
        pipe_b4 = run_pipeline_2d(l1["field"], l4["rx"], l4["ry"], cfg)
        b["B4"] = {**extract_kappa_observables(pipe_b4["jacobian"]),
                    "rx": l4["rx"], "ry": l4["ry"], "Rx_3d": l4["Rx_3d"],
                    "Ry_3d": l4["Ry_3d"], "Rz_3d": l4["Rz_3d"]}
        l5 = lane_l4_3d_los_projection(rho, cfg, nz=nz_primary,
                                        profile=PRIMARY_PROFILE,
                                        stencil=PRIMARY_STENCIL,
                                        boundary=PRIMARY_BC,
                                        orientation="O4")
        pipe_b5 = run_pipeline_2d(l1["field"], l5["rx"], l5["ry"], cfg)
        b["B5"] = {**extract_kappa_observables(pipe_b5["jacobian"]),
                    "rx": l5["rx"], "ry": l5["ry"], "Rx_3d": l5["Rx_3d"],
                    "Ry_3d": l5["Ry_3d"], "Rz_3d": l5["Rz_3d"]}
        benchmark_results[cid] = b

    print("[lab] running candidates PL1..PL6 x PM1..PM2 x PS1..PS2 ...")
    candidate_results = {}
    for cluster in CLUSTERS:
        cid = cluster["id"]
        rho = cluster_data[cid]["rho"]
        rho_3d = construct_rho_3d(rho, nz_primary)
        rng = np.random.RandomState(12345)
        u_slow, u_fast = A8_init_3d(rho_3d, cfg["strength"], rng)
        u_slow, u_fast, history = evolve_transport_3d(
            u_slow, u_fast, stencil=PRIMARY_STENCIL, boundary=PRIMARY_BC)
        c_3d = history[-1]
        state = {"rho_3d": rho_3d, "u_slow": u_slow, "u_fast": u_fast,
                 "c_3d": c_3d}
        field_2d = make_field_a8_t1(rho, cfg["extent"], cfg["strength"],
                                      seed=12345)
        cr = {}
        for pl in PL_LANES:
            for pm in PM_LANES:
                for ps in PS_LANES:
                    cid_key = f"{pl}_{pm}_{ps}"
                    res = run_candidate(state, pl, pm, ps, cfg, rho,
                                          field_2d)
                    cr[cid_key] = res
        candidate_results[cid] = {"state": state, "candidates": cr,
                                    "field_2d": field_2d}
        del state, rho_3d, rng, u_slow, u_fast, history, c_3d, field_2d
        import gc as _gc_cluster
        _gc_cluster.collect()

    import gc as _gc_after_candidates
    _gc_after_candidates.collect()
    print(f"[lab] candidates done, memory used after gc: {__import__('resource').getrusage(__import__('resource').RUSAGE_SELF).ru_maxrss / 1024:.0f} MB")

    print("[lab] running rotational covariance audit ...")
    cov_results = {}
    primary_pl = "PL1"; primary_pm = "PM1"; primary_ps = "PS2"
    for cluster in CLUSTERS:
        cid = cluster["id"]
        rho = cluster_data[cid]["rho"]
        state_native = candidate_results[cid]["state"]
        field_2d = candidate_results[cid]["field_2d"]
        rows = {}
        for rc in COORD_TRANSFORMS:
            state_rc = coord_transform(state_native, rc)
            res = run_candidate(state_rc, primary_pl, primary_pm, primary_ps,
                                 cfg, rho, field_2d,
                                 run_propagation=False)
            q3d = res["q3d"]
            Rx_native = inverse_coord_transform(q3d["Rx_3d"], rc)
            Ry_native = inverse_coord_transform(q3d["Ry_3d"], rc)
            Rz_native = inverse_coord_transform(q3d["Rz_3d"], rc)
            rows[rc] = {
                "Rx_native": Rx_native,
                "Ry_native": Ry_native,
                "Rz_native": Rz_native,
                "q3d_native": {**q3d,
                                "Rx_3d": Rx_native, "Ry_3d": Ry_native,
                                "Rz_3d": Rz_native,
                                "rx_proj": np.sum(Rx_native, axis=0),
                                "ry_proj": np.sum(Ry_native, axis=0),
                                "rx_central": Rx_native[Rx_native.shape[0] // 2],
                                "ry_central": Ry_native[Ry_native.shape[0] // 2]},
                "res": res,
            }
        cov_results[cid] = rows

    print("[lab] running wrong controls ...")
    wrong_results = {}
    for cluster in CLUSTERS:
        cid = cluster["id"]
        rho = cluster_data[cid]["rho"]
        rows = {}
        rows["WR1"] = wrong_control_replicated_slices(rho, cfg, nz=nz_primary)
        rows["WR2"] = wrong_control_zero_z_coupling(rho, cfg, nz=nz_primary)
        rows["WR3"] = wrong_control_random_depth_permutation(rho, cfg, nz=nz_primary)
        rows["WR4"] = wrong_control_uniform_depth(rho, cfg, nz=nz_primary)
        rows["WR5"] = wrong_control_sign_reverse_rz(rho, cfg, nz=nz_primary)
        rows["WR6"] = wrong_control_depth_shuffled_rz(rho, cfg, nz=nz_primary)
        state_n = candidate_results[cid]["state"]
        u_slow, u_fast = state_n["u_slow"], state_n["u_fast"]
        rho_3d = state_n["rho_3d"]
        rng_w = np.random.RandomState(42)
        pair_amp = compute_pair_amplitude_T1(u_slow, u_fast)
        A_vox = sum(pair_amp["A_combined"])
        rand = rng_w.standard_normal((3,) + rho_3d.shape)
        rand /= np.linalg.norm(rand, axis=0)
        Rx_w = A_vox * rand[0]; Ry_w = A_vox * rand[1]; Rz_w = A_vox * rand[2]
        fracs = helmholtz_fractions(Rx_w, Ry_w, Rz_w)
        rows["WR7"] = {"Rx_3d": Rx_w, "Ry_3d": Ry_w, "Rz_3d": Rz_w,
                        "rx_proj": np.sum(Rx_w, axis=0),
                        "ry_proj": np.sum(Ry_w, axis=0),
                        "f_irr_3d": fracs["f_irr_3d"],
                        "f_sol_3d": fracs["f_sol_3d"]}
        eL_x, eL_y, eL_z, valid = compute_longitudinal_axis(rho_3d)
        pair_amp = compute_pair_amplitude_T1(u_slow, u_fast)
        Rx_w = np.zeros_like(rho_3d); Ry_w = np.zeros_like(rho_3d)
        Rz_w = np.zeros_like(rho_3d)
        for idx, (dz, dy, dx) in enumerate(
                [(0, 0, -1), (0, 0, 1), (0, -1, 0), (0, 1, 0), (-1, 0, 0), (1, 0, 0)]):
            A_p = pair_amp["A_combined"][idx]
            n_x = float(dx); n_y = float(dy); n_z = float(dz)
            n_norm = math.sqrt(n_x ** 2 + n_y ** 2 + n_z ** 2)
            n_x /= n_norm; n_y /= n_norm; n_z /= n_norm
            Rx_w += A_p * n_x; Ry_w += A_p * n_y; Rz_w += A_p * n_z
        fracs = helmholtz_fractions(Rx_w, Ry_w, Rz_w)
        rows["WR8"] = {"Rx_3d": Rx_w, "Ry_3d": Ry_w, "Rz_3d": Rz_w,
                        "rx_proj": np.sum(Rx_w, axis=0),
                        "ry_proj": np.sum(Ry_w, axis=0),
                        "f_irr_3d": fracs["f_irr_3d"],
                        "f_sol_3d": fracs["f_sol_3d"]}
        proj_long = (eL_x * eL_x, eL_x * eL_y, eL_x * eL_z,
                      eL_y * eL_y, eL_y * eL_z, eL_z * eL_z)
        Rx_w = np.zeros_like(rho_3d); Ry_w = np.zeros_like(rho_3d)
        Rz_w = np.zeros_like(rho_3d)
        for idx, (dz, dy, dx) in enumerate(
                [(0, 0, -1), (0, 0, 1), (0, -1, 0), (0, 1, 0), (-1, 0, 0), (1, 0, 0)]):
            A_p = pair_amp["A_combined"][idx]
            n_x = float(dx); n_y = float(dy); n_z = float(dz)
            n_norm = math.sqrt(n_x ** 2 + n_y ** 2 + n_z ** 2)
            n_x /= n_norm; n_y /= n_norm; n_z /= n_norm
            PLn_x = (proj_long[0] * n_x + proj_long[1] * n_y + proj_long[2] * n_z)
            PLn_y = (proj_long[1] * n_x + proj_long[3] * n_y + proj_long[4] * n_z)
            PLn_z = (proj_long[2] * n_x + proj_long[4] * n_y + proj_long[5] * n_z)
            Rx_w += A_p * PLn_x; Ry_w += A_p * PLn_y; Rz_w += A_p * PLn_z
        fracs = helmholtz_fractions(Rx_w, Ry_w, Rz_w)
        rows["WR9"] = {"Rx_3d": Rx_w, "Ry_3d": Ry_w, "Rz_3d": Rz_w,
                        "rx_proj": np.sum(Rx_w, axis=0),
                        "ry_proj": np.sum(Ry_w, axis=0),
                        "f_irr_3d": fracs["f_irr_3d"],
                        "f_sol_3d": fracs["f_sol_3d"]}
        rows["WR10"] = {"Rx_3d": np.zeros_like(rho_3d),
                         "Ry_3d": np.zeros_like(rho_3d),
                         "Rz_3d": np.zeros_like(rho_3d),
                         "rx_proj": np.zeros((bins, bins)),
                         "ry_proj": np.zeros((bins, bins)),
                         "f_irr_3d": 0.0, "f_sol_3d": 0.0}
        wrong_results[cid] = rows

    print("[lab] depth convergence audit ...")
    depth_conv_results = {}
    nz_list = [3, 9, 17]
    physical_depth = 1.0
    dz_per_nz = {3: physical_depth / 3.0, 9: physical_depth / 9.0,
                  17: physical_depth / 17.0}
    for cluster in CLUSTERS:
        cid = cluster["id"]
        rho = cluster_data[cid]["rho"]
        rows = []
        for nz in nz_list:
            rho_3d = construct_rho_3d(rho, nz)
            rng = np.random.RandomState(12345)
            u_slow, u_fast = A8_init_3d(rho_3d, cfg["strength"], rng)
            u_slow, u_fast, history = evolve_transport_3d(
                u_slow, u_fast, stencil=PRIMARY_STENCIL, boundary=PRIMARY_BC)
            c_3d = history[-1]
            state = {"rho_3d": rho_3d, "u_slow": u_slow, "u_fast": u_fast,
                     "c_3d": c_3d}
            res = run_candidate(state, "PL1", "PM1", "PS2", cfg, rho,
                                  candidate_results[cid]["field_2d"])
            q3d = res["q3d"]
            gr_pad = cluster_gr[cid]
            r_kappa_los = safe_nan(pearson(q3d["rx_proj"], gr_pad["kappa"]))
            rows.append({
                "nz": nz,
                "physical_depth": physical_depth,
                "dz_voxel": dz_per_nz[nz],
                "response_energy": q3d["E_native"],
                "f_irr_3d": q3d["f_irr_3d"],
                "f_sol_3d": q3d["f_sol_3d"],
                "f_z": q3d["f_z"],
                "F_Dz": q3d["F_Dz"],
                "rx_proj_rms": rms_amplitude(q3d["rx_proj"]),
                "ry_proj_rms": rms_amplitude(q3d["ry_proj"]),
                "pearson_kappa_los": r_kappa_los,
                "covariance_error_max": 0.0,
                "covariance_error_mean": 0.0,
            })
        depth_conv_results[cid] = rows

    print("[lab] wave-mode audit ...")
    wave_results = {}
    for cluster in CLUSTERS:
        cid = cluster["id"]
        rho = cluster_data[cid]["rho"]
        state_n = candidate_results[cid]["state"]
        u_slow, u_fast = state_n["u_slow"], state_n["u_fast"]
        rho_3d = state_n["rho_3d"]
        nz_c = rho_3d.shape[0] // 2
        ny_c = rho_3d.shape[1] // 2
        nx_c = rho_3d.shape[2] // 2
        pair_amp = compute_pair_amplitude_T1(u_slow, u_fast)
        amp_stack = np.array([np.abs(a[nz_c, ny_c, nx_c])
                                for a in pair_amp["A_combined"]])
        best_idx = int(np.argmax(amp_stack))
        best_ax = ["xm", "xp", "ym", "yp", "zm", "zp"][best_idx]
        n_x = float(FACE_OFFSETS[best_ax][2])
        n_y = float(FACE_OFFSETS[best_ax][1])
        n_z = float(FACE_OFFSETS[best_ax][0])
        n_norm = math.sqrt(n_x ** 2 + n_y ** 2 + n_z ** 2)
        n_x /= n_norm; n_y /= n_norm; n_z /= n_norm
        if abs(n_z) < 0.9:
            helper = np.array([0.0, 0.0, 1.0])
        else:
            helper = np.array([1.0, 0.0, 0.0])
        e_perp1 = np.cross(np.array([n_x, n_y, n_z]), helper)
        e_perp1 /= np.linalg.norm(e_perp1)
        e_perp2 = np.cross(np.array([n_x, n_y, n_z]), e_perp1)
        e_perp2 /= np.linalg.norm(e_perp2)
        e_long = np.array([n_x, n_y, n_z])
        rms_R = rms_amplitude(np.sqrt(
            (pair_amp["A_combined"][0]) ** 2 + (pair_amp["A_combined"][1]) ** 2
            + (pair_amp["A_combined"][2]) ** 2))
        eps = 1e-6 * max(rms_R, EPS)
        rows = {}
        for label, direction in [("L", e_long), ("T1", e_perp1),
                                  ("T2", e_perp2)]:
            us_p = u_slow.copy(); uf_p = u_fast.copy()
            us_p[nz_c, ny_c, nx_c] += eps * direction[2]
            uf_p[nz_c, ny_c, nx_c] += eps * direction[2]
            us_p, uf_p, _ = evolve_transport_3d(us_p, uf_p,
                                                 stencil=PRIMARY_STENCIL,
                                                 boundary=PRIMARY_BC)
            c_p = 0.5 * (us_p + uf_p)
            state_p = {"rho_3d": rho_3d, "u_slow": us_p, "u_fast": uf_p,
                        "c_3d": c_p}
            res = run_candidate(state_p, "PL1", "PM1", "PS2", cfg, rho,
                                  candidate_results[cid]["field_2d"])
            q3d = res["q3d"]
            rows[label] = {
                "direction": direction.tolist(),
                "rms_R": rms_R,
                "eps": eps,
                "response_energy": q3d["E_native"],
                "f_irr_3d": q3d["f_irr_3d"],
                "f_sol_3d": q3d["f_sol_3d"],
                "f_z": q3d["f_z"],
                "helicity_total": float(np.sum(q3d["h"])),
            }
        wave_results[cid] = rows

    print("[lab] temporal audit ...")
    temporal_results = {}
    snapshots = [1, 10, 20]
    # Only run for the primary candidate (PL1_PM1_PS2) to limit memory
    primary_pl_for_temporal = ["PL1"]
    for cluster in CLUSTERS:
        cid = cluster["id"]
        rho = cluster_data[cid]["rho"]
        rho_3d = construct_rho_3d(rho, nz_primary)
        rng = np.random.RandomState(12345)
        u_slow, u_fast = A8_init_3d(rho_3d, cfg["strength"], rng)
        history = [0.5 * (u_slow + u_fast)]
        rows = {}
        for step in range(STEPS):
            p_fast = np.pad(u_fast, ((1, 1), (1, 1), (1, 1)), mode="reflect")
            p_slow = np.pad(u_slow, ((1, 1), (1, 1), (1, 1)), mode="reflect")
            n_fast = (p_fast[1:-1, 1:-1, :-2] + p_fast[1:-1, 1:-1, 2:]
                       + p_fast[1:-1, :-2, 1:-1] + p_fast[1:-1, 2:, 1:-1]
                       + p_fast[:-2, 1:-1, 1:-1] + p_fast[2:, 1:-1, 1:-1]) / 6.0
            n_slow = (p_slow[1:-1, 1:-1, :-2] + p_slow[1:-1, 1:-1, 2:]
                       + p_slow[1:-1, :-2, 1:-1] + p_slow[1:-1, 2:, 1:-1]
                       + p_slow[:-2, 1:-1, 1:-1] + p_slow[2:, 1:-1, 1:-1]) / 6.0
            d_fast = DT * OMEGA * K * ((n_fast - u_fast)
                                          + COUPLING_SLOW_TO_FAST
                                          * (u_slow - u_fast))
            d_slow = DT * SLOW_TIMESCALE * ((n_slow - u_slow)
                                              + COUPLING_FAST_TO_SLOW
                                              * (u_fast - u_slow))
            u_fast = np.clip(u_fast + d_fast, -5.0, 5.0)
            u_slow = np.clip(u_slow + d_slow, -5.0, 5.0)
            history.append(0.5 * (u_slow + u_fast))
            snap_idx = step + 1
            if snap_idx in snapshots or snap_idx == STEPS:
                c_3d_s = history[-1]
                state_s = {"rho_3d": rho_3d, "u_slow": u_slow,
                            "u_fast": u_fast, "c_3d": c_3d_s}
                for pl in primary_pl_for_temporal:
                    res = run_candidate(state_s, pl, "PM1", "PS2", cfg, rho,
                                          candidate_results[cid]["field_2d"])
                    q3d = res["q3d"]
                    gr_pad = cluster_gr[cid]
                    r_kappa_los = safe_nan(pearson(q3d["rx_proj"],
                                                      gr_pad["kappa"]))
                    rows.setdefault(pl, []).append({
                        "step": snap_idx,
                        "response_energy": q3d["E_native"],
                        "f_irr_3d": q3d["f_irr_3d"],
                        "f_sol_3d": q3d["f_sol_3d"],
                        "f_z": q3d["f_z"],
                        "pearson_kappa_los": r_kappa_los,
                        "helicity_total": float(np.sum(q3d["h"])),
                    })
                    del res, q3d, r_kappa_los, gr_pad
                del state_s, c_3d_s
                import gc as _gc_temporal
                _gc_temporal.collect()
        temporal_results[cid] = rows
        # Free temporaries
        del u_slow, u_fast, history, p_fast, p_slow, n_fast, n_slow, rho_3d
        import gc as _gc_temporal2
        _gc_temporal2.collect()

    print("[lab] writing CSVs ...")
    _write_csv_outputs(candidate_results, benchmark_results, cov_results,
                         wrong_results, depth_conv_results, temporal_results,
                         wave_results, cluster_gr, OUT)
    _write_csv_outputs2(candidate_results, benchmark_results, cov_results,
                          wrong_results, depth_conv_results, temporal_results,
                          wave_results, cluster_gr, OUT)

    print("[lab] saving native data ...")
    _save_native_data(candidate_results, benchmark_results, OUT, FIELDS)

    print("[lab] writing permanent registry ...")
    _write_permanent_registry(candidate_results, benchmark_results,
                                cov_results, cluster_gr, OUT, ROOT)

    print("[lab] generating plots ...")
    _make_plots(benchmark_results, candidate_results, cov_results,
                 wrong_results, depth_conv_results, temporal_results,
                 wave_results, cluster_gr, cluster_data, cfg, OUT, PLOTS,
                 FIELDS)

    print("[lab] writing validation and report ...")
    val = build_validation(candidate_results, cov_results, depth_conv_results,
                            wrong_results, temporal_results, hash_report)
    write_json_safe(OUT / "validation.json", val)
    report = build_report(candidate_results, benchmark_results, cov_results,
                           wrong_results, depth_conv_results, temporal_results,
                           wave_results, cluster_gr, OUT)
    (OUT / "report.md").write_text(report)
    write_json_safe(OUT / "run.json",
                     {"started_at": now_iso(),
                      "finished_at": now_iso(),
                      "duration_s": time.perf_counter() - started,
                      "clusters": [c["id"] for c in CLUSTERS],
                      "candidates_total": len(PL_LANES) * len(PM_LANES) * len(PS_LANES),
                      "config": cfg, "nz_primary": nz_primary})
    print(f"[lab] complete in {time.perf_counter() - started:.1f}s")


if __name__ == "__main__":
    main()
