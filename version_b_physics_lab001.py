#!/usr/bin/env python3
"""PBUF VERSION-B PHYSICS-LAB-001 - Local Response Hypothesis Survey.

Tests ten candidate local response laws inside the FROZEN Version 1
weak-lensing laboratory (LAB-FREEZE-001 / WEAK-LENSING-SCIENCE-001).

Only the response law (rx, ry) is varied.  No modification to
- Constitutive Version A
- Frozen transport (propagation)
- Source-plane geometry (Launch B)
- Observable extraction (Jacobian)
- Numerical configuration (20 000 photons, 256^2 grid, Delta s / 2)

Production configuration (frozen minimum production):
- 20 000 photons
- 256^2 constitutive grid
- Cartesian 2D launch (Launch B)
- Jacobian observable
- Five benchmark clusters (Abell 2744, MACS J0416, MACS J1149,
  Abell S1063, Abell 370).

Each candidate is tested exactly once on each cluster with the
predefined fixed parameters.  No parameter fitting.  No optimisation.

For every candidate the script records
- RMS kappa
- RMS gamma
- Pearson kappa
- Pearson gamma
- SSIM (kappa, gamma)
- kappa bias (mean residual)
- gamma bias
- conservation max
- runtime

Cross-cluster summary table reports median Pearson kappa, median
Pearson gamma, median SSIM, mean kappa bias per candidate.

Required outputs (runs/version_b_physics_lab001/):
- report.md
- candidate_summary.csv
- cross_cluster_statistics.csv
- candidate_ranking.csv
- run.json
- validation.json
- plots/candidate_rankings.png
- plots/bias_comparison.png
- plots/pearson_comparison.png
- plots/ssim_comparison.png
- plots/cluster_performance.png
- plots/response_family_summary.png
"""
from __future__ import annotations

import csv
import hashlib
import json
import sys
import time
import traceback
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from scipy.ndimage import map_coordinates

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from weak_lensing_observation001 import (
    LENS as BASE_LENS,
    make_field,
    propagate,
    file_sha256,
    resample_to_grid,
    compare_arrays,
    ssim_index,
    pearson_corr,
)
import observable_lab001 as obs_lab
import source_plane_lab001 as src_lab


BENCHMARK_DIR = ROOT / "PBUF_benchmark"
DEFAULT_OUT = ROOT / "runs" / "version_b_physics_lab001"
PLOTS = DEFAULT_OUT / "plots"

# Frozen minimum production configuration (LAB-FREEZE-001 Section 1.5).
MIN_PRODUCTION = {
    "label": "minimum_production",
    "nphotons": 20000,
    "grid_n": 256,
    "step": 0.03,
    "steps": 160,
    "y_span": 3.0,
    "extent": 8.0,
    "strength": 0.18,
    "bins": 64,
}

CLUSTERS = [
    {"id": "Abell2744",  "label": "Abell 2744",   "slug": "abell2744",
     "directory": "WL-001_Abell2744"},
    {"id": "MACS0416",   "label": "MACS J0416",   "slug": "macs0416",
     "directory": "WL-002_MACS0416"},
    {"id": "MACS1149",   "label": "MACS J1149",   "slug": "macs1149",
     "directory": "WL-003_MACS1149"},
    {"id": "AbellS1063", "label": "Abell S1063",  "slug": "abells1063",
     "directory": "WL-004_AbellS1063"},
    {"id": "Abell370",   "label": "Abell 370",    "slug": "abell370",
     "directory": "WL-005_Abell370"},
]

# Frozen source files (from LAB-FREEZE-001 checksums.csv).
EXPECTED_HASHES = {
    "constitutive_equations.py":
        "e2c789d19fd559753519704c6668c7a2879c53eb61a315604ec81af6795aca9f",
    "weak_lensing_observation001.py":
        "a5c3632fec9adfc2659d5c283d07c599db6db7edb2e34e4aebd84e35434642bc",
    "observable_lab001.py":
        "2867c0bf94fabe3fbba0264d5d272ecadbdf45fabb341a6a8fe77972fbaec132",
    "source_plane_lab001.py":
        "efa9d74924cb61a3b48a69fa075055512d86391d03194be342597420bc353de4",
    "numerical_convergence001.py":
        "0442f878713de6530b5a1b1844b8ece037852d461bcb695360e8a3345fd58f29",
}


# =============================================================================
# Local response laws
# =============================================================================
# A local response law receives the constitutive field, the coordinate
# arrays and the gradient field, and returns (rx, ry) response components
# over the same grid.
#
# All candidates preserve:
#   * 90 deg transverse direction (rx = -A * g_hat_y, ry = A * g_hat_x)
#     EXCEPT Candidate 3 (cooperative neighbour response), where the
#     gradient itself is the cell-averaged gradient.
#   * Neighbour-to-neighbour propagation (no change to propagate()).
#   * Locality (only the cell and its immediate neighbours are consulted).
#
# Each candidate is a pure function of the frozen constitutive field.
# No parameter is fitted.
#
# The 10 candidates (Candidate 1 is the control, Candidates 2-10 are
# alternative hypotheses).
# =============================================================================


def _rotate_90(gx: np.ndarray, gy: np.ndarray, A: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """90 deg right-handed rotation: r = R_90(g) * A/|g|.

    gx, gy are the gradient components.
    A is the desired response magnitude (same units as the gradient magnitude).
    """
    g = np.hypot(gx, gy)
    g_safe = np.maximum(g, 1e-15)
    rx = -A * (gy / g_safe)
    ry = +A * (gx / g_safe)
    return rx, ry


def _box3(field: np.ndarray) -> np.ndarray:
    """3x3 box average using reflection padding (one neighbour ring)."""
    p = np.pad(field, 1, mode="reflect")
    out = (
        p[0:-2, 0:-2] + p[0:-2, 1:-1] + p[0:-2, 2:]
        + p[1:-1, 0:-2] + p[1:-1, 1:-1] + p[1:-1, 2:]
        + p[2:, 0:-2] + p[2:, 1:-1] + p[2:, 2:]
    ) / 9.0
    return out


def _laplacian5(field: np.ndarray) -> np.ndarray:
    """5-point discrete Laplacian using reflection padding.

    L f = f_xx + f_yy approximated by (f[i+1] + f[i-1] + f[j+1] + f[j-1]
    - 4 f[i,j]) / h^2.  We return the un-normalised stencil; downstream
    candidates may include the grid spacing if they wish.
    """
    p = np.pad(field, 1, mode="reflect")
    L = (p[0:-2, 1:-1] + p[2:, 1:-1] + p[1:-1, 0:-2] + p[1:-1, 2:]
         - 4.0 * p[1:-1, 1:-1])
    return L


# -----------------------------------------------------------------------------
# Candidate 1: Gradient (control).  Frozen linear response.
# -----------------------------------------------------------------------------
def candidate_1_gradient(c, xgrid, ygrid, gx, gy, g):
    """Response = |grad C|.  This is the frozen Version A response."""
    A = g.copy()
    return _rotate_90(gx, gy, A)


# -----------------------------------------------------------------------------
# Candidate 2: Local Neighbour Coherence
# Response magnitude depends on the local gradient AND the alignment of
# neighbouring gradients.  Strengthened only when neighbouring gradients
# point in similar directions.  Entirely local.
# -----------------------------------------------------------------------------
def candidate_2_neighbour_coherence(c, xgrid, ygrid, gx, gy, g):
    """A = |grad C| * (1 + mean_cos(theta_self, theta_neighbours)) / 2.

    mean_cos is the average over the 8 immediate neighbours of the cosine
    of the angle between the cell's unit-gradient vector and the
    neighbour's unit-gradient vector.  In [-1, 1]; the (1+x)/2 mapping
    yields a factor in [0, 1] that suppresses the response in regions of
    incoherent gradients.
    """
    g_safe = np.maximum(g, 1e-15)
    gxh = gx / g_safe
    gyh = gy / g_safe
    p = np.pad(gxh, 1, mode="reflect")
    q = np.pad(gyh, 1, mode="reflect")
    gxh_pad = p[1:-1, 1:-1]
    gyh_pad = q[1:-1, 1:-1]
    # cosine with each of 8 neighbours (use reflection padding).
    cos_sum = np.zeros_like(gxh)
    n_count = 0
    for di in (-1, 0, 1):
        for dj in (-1, 0, 1):
            if di == 0 and dj == 0:
                continue
            nbx = p[1+di:1+di+p.shape[0]-2, 1+dj:1+dj+p.shape[1]-2]
            nby = q[1+di:1+di+q.shape[0]-2, 1+dj:1+dj+q.shape[1]-2]
            cos_sum += gxh_pad * nbx + gyh_pad * nby
            n_count += 1
    mean_cos = cos_sum / float(n_count)
    factor = 0.5 * (1.0 + mean_cos)
    A = g * factor
    return _rotate_90(gx, gy, A)


# -----------------------------------------------------------------------------
# Candidate 3: Cooperative Neighbour Response
# Each cell receives contributions from all immediate neighbours BEFORE
# computing its response.  No propagation beyond the nearest neighbourhood.
# -----------------------------------------------------------------------------
def candidate_3_cooperative(c, xgrid, ygrid, gx, gy, g):
    """Cell-averaged gradient (3x3 including the cell)."""
    gx_loc = _box3(gx)
    gy_loc = _box3(gy)
    g_loc = np.hypot(gx_loc, gy_loc)
    return _rotate_90(gx_loc, gy_loc, g_loc)


# -----------------------------------------------------------------------------
# Candidate 4: Elastic Memory
# One-step local persistence.  Response depends on current gradient AND
# previous local strain.  Memory is strictly one update step.
#
# Implementation: r_new = (1-w) * R(g) + w * R(g_prev) where
# R(g_prev) is the linear response from the cell's previous-step
# gradient.  Because the pipeline computes the response field once, we
# interpret "previous step" as the cell immediately upstream in -x,
# i.e. R(g)[i, max(0, j-1)].  This represents the strain that has
# persisted from the previous local update.
# -----------------------------------------------------------------------------
def candidate_4_elastic_memory(c, xgrid, ygrid, gx, gy, g):
    w = 0.5  # fixed weight; no fitting.
    # Linear response.
    rx_self = -g * (gy / np.maximum(g, 1e-15))
    ry_self = +g * (gx / np.maximum(g, 1e-15))
    # "Previous" strain = the linear response one cell upstream in -x.
    rx_prev = np.roll(rx_self, 1, axis=1)
    ry_prev = np.roll(ry_self, 1, axis=1)
    rx_prev[:, 0] = rx_self[:, 0]
    ry_prev[:, 0] = ry_self[:, 0]
    rx = (1.0 - w) * rx_self + w * rx_prev
    ry = (1.0 - w) * ry_self + w * ry_prev
    return rx, ry


# -----------------------------------------------------------------------------
# Candidate 5: Gradient Curvature
# Response depends on local gradient AND local variation of the gradient.
# Do not replace the gradient; augment it.
# -----------------------------------------------------------------------------
def candidate_5_gradient_curvature(c, xgrid, ygrid, gx, gy, g):
    """A = |grad C| + alpha * |L C| where L is the discrete Laplacian.

    L C is the curvature of the constitutive field.  alpha = 0.5,
    fixed.
    """
    alpha = 0.5
    LC = _laplacian5(c)
    A = g + alpha * np.abs(LC)
    return _rotate_90(gx, gy, A)


# -----------------------------------------------------------------------------
# Candidate 6: Phase-Coherent Response
# Transport remains unchanged.  Only the response magnitude depends on
# the local phase coherence accumulated over immediate neighbours.
# No global phase.
#
# Phase of cell = atan2(gy, gx).  Coherence is the mean cosine of the
# phase difference between the cell and each of its 8 neighbours.  In
# [0, 1] (cosine is even, so phase coherence is non-negative).
# -----------------------------------------------------------------------------
def candidate_6_phase_coherence(c, xgrid, ygrid, gx, gy, g):
    phi = np.arctan2(gy, gx)
    p = np.pad(phi, 1, mode="reflect")
    cos_sum = np.zeros_like(phi)
    n_count = 0
    for di in (-1, 0, 1):
        for dj in (-1, 0, 1):
            if di == 0 and dj == 0:
                continue
            nbphi = p[1+di:1+di+p.shape[0]-2, 1+dj:1+dj+p.shape[1]-2]
            cos_sum += np.cos(phi - nbphi)
            n_count += 1
    phase_coherence = cos_sum / float(n_count)
    A = g * phase_coherence
    return _rotate_90(gx, gy, A)


# -----------------------------------------------------------------------------
# Candidate 7: Relaxation Response
# Finite local relaxation.  Response depends on present neighbour
# interaction AND incomplete relaxation from the immediately preceding
# update.  Strictly local.
#
# Implementation: one Jacobi relaxation step toward the neighbour mean
# of the linear response.  r_new = (1-w) * r_linear + w * mean(r_neighbours).
# w = 0.5 (50% incomplete relaxation toward the neighbourhood mean).
# This is the "finite local relaxation" toward the neighbour-averaged
# strain, with w < 1 encoding the "incomplete" character.
# -----------------------------------------------------------------------------
def candidate_7_relaxation(c, xgrid, ygrid, gx, gy, g):
    w = 0.5  # fixed; one Jacobi step toward neighbour mean.
    rx_self = -g * (gy / np.maximum(g, 1e-15))
    ry_self = +g * (gx / np.maximum(g, 1e-15))
    # Neighbour mean (3x3 box including the cell).
    rx_smooth = _box3(rx_self)
    ry_smooth = _box3(ry_self)
    rx = (1.0 - w) * rx_self + w * rx_smooth
    ry = (1.0 - w) * ry_self + w * ry_smooth
    return rx, ry


# -----------------------------------------------------------------------------
# Candidate 8: Weak-Gradient Enhancement
# Small gradients receive a predefined enhancement.  Large gradients
# remain unchanged.  The enhancement function is fixed before testing.
# No fitting.
#
# A = |grad C| + epsilon * exp(-|grad C| / sigma) for |grad C| < sigma.
# Otherwise A = |grad C|.
# -----------------------------------------------------------------------------
def candidate_8_weak_gradient_enhancement(c, xgrid, ygrid, gx, gy, g):
    epsilon = 0.05
    sigma = 0.05
    enhancement = epsilon * np.exp(-g / sigma)
    A = g + enhancement
    return _rotate_90(gx, gy, A)


# -----------------------------------------------------------------------------
# Candidate 9: Constitutive Coupling
# Response depends upon constitutive state AND local gradient using one
# fixed functional form.
# -----------------------------------------------------------------------------
def candidate_9_constitutive_coupling(c, xgrid, ygrid, gx, gy, g):
    """A = |grad C| * (1 + beta * C).  beta = 1.0, fixed."""
    beta = 1.0
    A = g * (1.0 + beta * c)
    return _rotate_90(gx, gy, A)


# -----------------------------------------------------------------------------
# Candidate 10: Combined Local Response
# Combine neighbour coherence AND one-step elastic memory.  Intentionally
# the most complex candidate.
# -----------------------------------------------------------------------------
def candidate_10_combined(c, xgrid, ygrid, gx, gy, g):
    # Coherence factor (Candidate 2).
    g_safe = np.maximum(g, 1e-15)
    gxh = gx / g_safe
    gyh = gy / g_safe
    p = np.pad(gxh, 1, mode="reflect")
    q = np.pad(gyh, 1, mode="reflect")
    gxh_pad = p[1:-1, 1:-1]
    gyh_pad = q[1:-1, 1:-1]
    cos_sum = np.zeros_like(gxh)
    n_count = 0
    for di in (-1, 0, 1):
        for dj in (-1, 0, 1):
            if di == 0 and dj == 0:
                continue
            nbx = p[1+di:1+di+p.shape[0]-2, 1+dj:1+dj+p.shape[1]-2]
            nby = q[1+di:1+di+q.shape[0]-2, 1+dj:1+dj+q.shape[1]-2]
            cos_sum += gxh_pad * nbx + gyh_pad * nby
            n_count += 1
    mean_cos = cos_sum / float(n_count)
    coherence_factor = 0.5 * (1.0 + mean_cos)

    # Memory factor (Candidate 4) - ratio of upstream gradient to current.
    w_mem = 0.5
    rx_self = -g * (gy / g_safe)
    ry_self = +g * (gx / g_safe)
    rx_prev = np.roll(rx_self, 1, axis=1)
    ry_prev = np.roll(ry_self, 1, axis=1)
    rx_prev[:, 0] = rx_self[:, 0]
    ry_prev[:, 0] = ry_self[:, 0]
    rx = (1.0 - w_mem) * rx_self + w_mem * rx_prev
    ry = (1.0 - w_mem) * ry_self + w_mem * ry_prev
    # Apply coherence factor multiplicatively.
    rx *= coherence_factor
    ry *= coherence_factor
    return rx, ry


# Candidate registry.
@dataclass(frozen=True)
class CandidateSpec:
    number: int
    name: str
    family: str
    description: str
    law: Callable
    notes: str = ""


CANDIDATES = [
    CandidateSpec(
        1, "Gradient (control)", "gradient",
        "Response = |grad C|; frozen Version A control.",
        candidate_1_gradient,
        notes="Linear, frozen. Used as the control baseline.",
    ),
    CandidateSpec(
        2, "Local Neighbour Coherence", "neighbour coherence",
        "Magnitude scaled by (1 + mean_cos)/2 over 8 neighbours.",
        candidate_2_neighbour_coherence,
        notes="Strengthens where neighbouring gradients align.",
    ),
    CandidateSpec(
        3, "Cooperative Neighbour Response", "cooperative response",
        "Cell-averaged gradient (3x3 box).",
        candidate_3_cooperative,
        notes="Replaces gradient with neighbour average before rotating.",
    ),
    CandidateSpec(
        4, "Elastic Memory", "elastic memory",
        "r_new = (1-w)*R(g) + w*R(g_upstream); w = 0.5.",
        candidate_4_elastic_memory,
        notes="One-step persistence with upstream previous-step strain.",
    ),
    CandidateSpec(
        5, "Gradient Curvature", "gradient curvature",
        "A = |grad C| + 0.5 * |Laplacian C|.",
        candidate_5_gradient_curvature,
        notes="Augments magnitude with constitutive curvature.",
    ),
    CandidateSpec(
        6, "Phase-Coherent Response", "phase coherence",
        "A = |grad C| * mean_cos(phase differences).",
        candidate_6_phase_coherence,
        notes="Modulates magnitude by phase coherence over neighbours.",
    ),
    CandidateSpec(
        7, "Relaxation Response", "relaxation",
        "One Jacobi relaxation step toward neighbour mean of response.",
        candidate_7_relaxation,
        notes="w = 0.5 (50% incomplete relaxation), fixed. Strictly local.",
    ),
    CandidateSpec(
        8, "Weak-Gradient Enhancement", "weak-gradient enhancement",
        "A = |grad C| + 0.05 * exp(-|grad C| / 0.05).",
        candidate_8_weak_gradient_enhancement,
        notes="Predefined enhancement; no fitting.",
    ),
    CandidateSpec(
        9, "Constitutive Coupling", "constitutive coupling",
        "A = |grad C| * (1 + C).",
        candidate_9_constitutive_coupling,
        notes="Couples response to local constitutive state.",
    ),
    CandidateSpec(
        10, "Combined Local Response", "combined response",
        "Combine neighbour coherence (Cand 2) and elastic memory (Cand 4).",
        candidate_10_combined,
        notes="Most complex candidate; combination of Candidates 2 and 4.",
    ),
]


# =============================================================================
# Pipeline runner
# =============================================================================
def matter_proxy_from_kappa(kappa_native: np.ndarray, grid_n: int, extent: float) -> np.ndarray:
    """Frozen matter input rule: rho = max(kappa, 0) / max(max(kappa, 0))."""
    rho_pipeline = resample_to_grid(kappa_native, grid_n, extent)
    rho_pos = np.maximum(rho_pipeline, 0.0)
    rho_max = float(rho_pos.max())
    if rho_max > 0:
        rho_pos = rho_pos / rho_max
    return rho_pos


def compute_field(rho: np.ndarray, extent: float, strength: float, grid_n: int) -> dict:
    """Compute the constitutive field and its gradient on the grid.

    Uses the frozen Version A constitutive equation from
    `weak_lensing_observation001.make_field` (the gradient and gradient
    direction are computed identically).
    """
    x = np.linspace(-extent, extent, grid_n)
    y = np.linspace(-extent, extent, grid_n)
    X, Y = np.meshgrid(x, y, indexing="xy")
    cfg = type("Config", (), {"deformation_strength": strength})()
    from constitutive_equations import get_equation
    c = get_equation("A").solve(rho, cfg)
    gy, gx = np.gradient(c, x, y, edge_order=1)
    g = np.hypot(gx, gy)
    return {
        "xgrid": x, "ygrid": y, "X": X, "Y": Y,
        "rho": rho, "c": c,
        "gx": gx, "gy": gy, "g_magnitude": g,
    }


def apply_candidate(field: dict, candidate: CandidateSpec) -> tuple[np.ndarray, np.ndarray]:
    """Apply a candidate local response law. Returns (rx, ry)."""
    rx, ry = candidate.law(
        field["c"], field["xgrid"], field["ygrid"],
        field["gx"], field["gy"], field["g_magnitude"]
    )
    return rx, ry


def load_observation_full(cluster: dict) -> dict:
    """Load kappa/gamma/gamma1/gamma2 FITS files for the cluster."""
    folder = BENCHMARK_DIR / cluster["directory"]
    out = {"folder": str(folder), "files": {}, "shas": {}}
    keys = ("kappa", "gamma", "gamma1", "gamma2")
    for k in keys:
        p = folder / f"hlsp_frontier_model_{cluster['slug']}_merten_v1_{k}.fits"
        with fits.open(p) as h:
            out[k] = np.asarray(h[0].data, dtype=np.float64)
        out["files"][k] = str(p)
        out["shas"][k] = file_sha256(p)
    return out


def resample_observation(obs: dict, bins: int, extent: float) -> dict:
    return {
        "kappa": resample_to_grid(obs["kappa"], bins, extent),
        "gamma": resample_to_grid(obs["gamma"], bins, extent),
        "gamma1": resample_to_grid(obs["gamma1"], bins, extent),
        "gamma2": resample_to_grid(obs["gamma2"], bins, extent),
    }


def run_candidate_on_cluster(cluster: dict, cfg: dict, candidate: CandidateSpec) -> dict:
    """Run one candidate on one cluster at one configuration.

    Returns a dict with predictions, observables, residuals, comparison
    metrics, conservation error and runtime.
    """
    folder = BENCHMARK_DIR / cluster["directory"]
    kappa_path = folder / f"hlsp_frontier_model_{cluster['slug']}_merten_v1_kappa.fits"
    with fits.open(kappa_path) as h:
        kappa_native = np.asarray(h[0].data, dtype=np.float64)
    rho = matter_proxy_from_kappa(kappa_native, cfg["grid_n"], cfg["extent"])

    field = compute_field(rho, cfg["extent"], cfg["strength"], cfg["grid_n"])

    # Apply the candidate local response law.
    rx, ry = apply_candidate(field, candidate)
    field_candidate = dict(field)
    field_candidate["rx"] = rx
    field_candidate["ry"] = ry
    field_candidate["response_direction"] = np.arctan2(ry, rx)
    field_candidate["response_magnitude"] = np.hypot(rx, ry)

    x0, y0, vx0, vy0 = src_lab.launch_B_cartesian(cfg["nphotons"])

    propagation_started = time.perf_counter()
    photons = propagate(field_candidate, cfg["step"], cfg["steps"], x0, y0, vx0, vy0)
    photons["x0"] = x0
    photons["y0"] = y0
    propagation_seconds = time.perf_counter() - propagation_started

    xs_i = x0.copy(); ys_i = y0.copy()
    xs_f = photons["x"].copy(); ys_f = photons["y"].copy()

    extraction_started = time.perf_counter()
    jacobian_result = obs_lab.method_jacobian(
        xs_i, ys_i, xs_f, ys_f, cfg["extent"], cfg["bins"]
    )
    extraction_seconds = time.perf_counter() - extraction_started

    pred = {
        "kappa": jacobian_result["convergence"],
        "gamma1": jacobian_result["shear_g1"],
        "gamma2": jacobian_result["shear_g2"],
        "gamma_mag": jacobian_result["shear_magnitude"],
        "deflection_x": jacobian_result["deflection_x"],
        "deflection_y": jacobian_result["deflection_y"],
        "magnification": jacobian_result["magnification"],
    }

    obs = load_observation_full(cluster)
    obs_grid = resample_observation(obs, cfg["bins"], cfg["extent"])

    # Quantitative comparison.
    out = {}
    for key in ("kappa", "gamma1", "gamma2"):
        cmp = compare_arrays(pred[key], obs_grid[key])
        cmp["ssim"] = ssim_index(pred[key], obs_grid[key])
        out[key] = cmp
    cmp_g = compare_arrays(pred["gamma_mag"], obs_grid["gamma"])
    cmp_g["ssim"] = ssim_index(pred["gamma_mag"], obs_grid["gamma"])
    out["gamma_mag"] = cmp_g

    # Finite-pixel statistics.
    finite_pred_kappa = pred["kappa"][np.isfinite(pred["kappa"])]
    finite_obs_kappa = obs_grid["kappa"][np.isfinite(obs_grid["kappa"])]
    finite_pred_gamma = pred["gamma_mag"][np.isfinite(pred["gamma_mag"])]
    finite_obs_gamma = obs_grid["gamma"][np.isfinite(obs_grid["gamma"])]
    out["kappa_predicted_rms"] = float(np.sqrt(np.mean(finite_pred_kappa ** 2))) if finite_pred_kappa.size else float("nan")
    out["kappa_observed_rms"] = float(np.sqrt(np.mean(finite_obs_kappa ** 2))) if finite_obs_kappa.size else float("nan")
    out["gamma_predicted_rms"] = float(np.sqrt(np.mean(finite_pred_gamma ** 2))) if finite_pred_gamma.size else float("nan")
    out["gamma_observed_rms"] = float(np.sqrt(np.mean(finite_obs_gamma ** 2))) if finite_obs_gamma.size else float("nan")
    out["kappa_predicted_mean"] = float(finite_pred_kappa.mean()) if finite_pred_kappa.size else float("nan")
    out["kappa_observed_mean"] = float(finite_obs_kappa.mean()) if finite_obs_kappa.size else float("nan")
    out["gamma_predicted_mean"] = float(finite_pred_gamma.mean()) if finite_pred_gamma.size else float("nan")
    out["gamma_observed_mean"] = float(finite_obs_gamma.mean()) if finite_obs_gamma.size else float("nan")

    # Residual map (pred - obs) over finite pixels.
    resid_kappa = np.where(np.isfinite(pred["kappa"]) & np.isfinite(obs_grid["kappa"]),
                            pred["kappa"] - obs_grid["kappa"], np.nan)
    resid_gamma = np.where(np.isfinite(pred["gamma_mag"]) & np.isfinite(obs_grid["gamma"]),
                            pred["gamma_mag"] - obs_grid["gamma"], np.nan)

    finite_resid_k = resid_kappa[np.isfinite(resid_kappa)]
    finite_resid_g = resid_gamma[np.isfinite(resid_gamma)]

    bias_kappa = float(np.mean(finite_resid_k)) if finite_resid_k.size else float("nan")
    bias_gamma = float(np.mean(finite_resid_g)) if finite_resid_g.size else float("nan")
    std_resid_kappa = float(np.std(finite_resid_k)) if finite_resid_k.size else float("nan")
    std_resid_gamma = float(np.std(finite_resid_g)) if finite_resid_g.size else float("nan")

    cons_max = float(np.max(photons["conservation"])) if photons["conservation"].size else 0.0

    return {
        "cluster_id": cluster["id"],
        "cluster_label": cluster["label"],
        "config": cfg,
        "candidate_number": candidate.number,
        "candidate_name": candidate.name,
        "candidate_family": candidate.family,
        "field": field_candidate,
        "photons": photons,
        "pred": pred,
        "obs_grid": obs_grid,
        "comparison": out,
        "residuals": {
            "kappa": resid_kappa,
            "gamma_mag": resid_gamma,
        },
        "bias_kappa": bias_kappa,
        "bias_gamma": bias_gamma,
        "std_resid_kappa": std_resid_kappa,
        "std_resid_gamma": std_resid_gamma,
        "n_finite_pixels_kappa": int(finite_resid_k.size),
        "n_finite_pixels_gamma": int(finite_resid_g.size),
        "propagation_seconds": float(propagation_seconds),
        "extraction_seconds": float(extraction_seconds),
        "total_seconds": float(propagation_seconds + extraction_seconds),
        "max_conservation_error": cons_max,
        "n_photons": int(cfg["nphotons"]),
        "grid_n": int(cfg["grid_n"]),
        "step": float(cfg["step"]),
        "steps": int(cfg["steps"]),
    }


def verify_frozen_hashes() -> dict:
    """Verify the five frozen source files against LAB-FREEZE-001."""
    result = {"ok": True, "files": {}}
    for name, expected in EXPECTED_HASHES.items():
        actual = file_sha256(ROOT / name)
        ok = (actual == expected)
        result["files"][name] = {
            "expected_sha256": expected,
            "actual_sha256": actual,
            "match": ok,
        }
        if not ok:
            result["ok"] = False
    return result


# =============================================================================
# Metric helpers
# =============================================================================
def _safe_median(values: list[float]) -> float:
    arr = np.array(values, dtype=np.float64)
    arr = arr[np.isfinite(arr)]
    if arr.size == 0:
        return float("nan")
    return float(np.median(arr))


def _safe_mean(values: list[float]) -> float:
    arr = np.array(values, dtype=np.float64)
    arr = arr[np.isfinite(arr)]
    if arr.size == 0:
        return float("nan")
    return float(np.mean(arr))


def _safe_std(values: list[float]) -> float:
    arr = np.array(values, dtype=np.float64)
    arr = arr[np.isfinite(arr)]
    if arr.size < 2:
        return float("nan")
    return float(np.std(arr))


# =============================================================================
# Output writers
# =============================================================================
def write_candidate_summary_csv(out_root: Path, rows: list[dict]) -> None:
    """candidate_summary.csv - one row per (candidate, cluster)."""
    path = out_root / "candidate_summary.csv"
    fields = [
        "candidate_number", "candidate_name", "candidate_family",
        "cluster_id", "cluster_label",
        "rms_kappa", "rms_gamma",
        "pearson_kappa", "pearson_gamma",
        "ssim_kappa", "ssim_gamma",
        "kappa_bias", "gamma_bias",
        "std_resid_kappa", "std_resid_gamma",
        "max_conservation_error",
        "runtime_seconds",
        "n_finite_pixels_kappa", "n_finite_pixels_gamma",
        "n_photons", "grid_n", "step", "steps",
    ]
    with path.open("w", newline="") as h:
        writer = csv.DictWriter(h, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({k: row.get(k, "") for k in fields})


def write_cross_cluster_csv(out_root: Path, cross_rows: list[dict]) -> None:
    """cross_cluster_statistics.csv - one row per candidate."""
    path = out_root / "cross_cluster_statistics.csv"
    fields = [
        "candidate_number", "candidate_name", "candidate_family",
        "median_pearson_kappa", "median_pearson_gamma",
        "median_ssim_kappa", "median_ssim_gamma",
        "median_rms_kappa", "median_rms_gamma",
        "mean_kappa_bias", "std_kappa_bias",
        "mean_gamma_bias", "std_gamma_bias",
        "mean_pearson_kappa", "mean_pearson_gamma",
        "max_conservation_error",
        "median_runtime_seconds",
        "n_clusters_with_positive_pearson_kappa",
        "n_clusters_with_positive_pearson_gamma",
    ]
    with path.open("w", newline="") as h:
        writer = csv.DictWriter(h, fieldnames=fields)
        writer.writeheader()
        for row in cross_rows:
            writer.writerow({k: row.get(k, "") for k in fields})


def write_candidate_ranking_csv(out_root: Path, ranking_rows: list[dict]) -> None:
    """candidate_ranking.csv - candidates ranked by median Pearson kappa."""
    path = out_root / "candidate_ranking.csv"
    fields = [
        "rank", "candidate_number", "candidate_name", "candidate_family",
        "median_pearson_kappa", "median_pearson_gamma",
        "median_ssim_kappa",
        "mean_kappa_bias", "improvement_vs_control_pearson_kappa",
        "improvement_vs_control_ssim_kappa",
        "improvement_vs_control_kappa_bias",
    ]
    with path.open("w", newline="") as h:
        writer = csv.DictWriter(h, fieldnames=fields)
        writer.writeheader()
        for row in ranking_rows:
            writer.writerow({k: row.get(k, "") for k in fields})


# =============================================================================
# Plot helpers
# =============================================================================
def save_candidate_rankings_plot(out_path: Path, cross_rows: list[dict]) -> None:
    """Four-panel ranking by median Pearson kappa, median Pearson gamma,
    median SSIM kappa, and mean kappa bias."""
    sorted_pk = sorted(cross_rows, key=lambda r: -r["median_pearson_kappa"])
    sorted_pg = sorted(cross_rows, key=lambda r: -r["median_pearson_gamma"])
    sorted_sk = sorted(cross_rows, key=lambda r: -r["median_ssim_kappa"])
    sorted_bias = sorted(cross_rows, key=lambda r: r["mean_kappa_bias"])  # more negative = stronger underprediction

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    for ax, sorted_rows, key, title in zip(
        axes.flat,
        [sorted_pk, sorted_pg, sorted_sk, sorted_bias],
        ["median_pearson_kappa", "median_pearson_gamma",
         "median_ssim_kappa", "mean_kappa_bias"],
        ["Median Pearson kappa (higher better)",
         "Median Pearson gamma (higher better)",
         "Median SSIM kappa (higher better)",
         "Mean kappa bias (closer to 0 better)"],
    ):
        labels = [f"{int(r['candidate_number'])} {r['candidate_name']}" for r in sorted_rows]
        vals = [r[key] for r in sorted_rows]
        ax.barh(range(len(vals)), vals, color="steelblue", edgecolor="black")
        ax.set_yticks(range(len(vals)))
        ax.set_yticklabels(labels, fontsize=8)
        ax.invert_yaxis()
        ax.set(xlabel=key, title=title)
        ax.grid(axis="x", alpha=0.3)
    fig.suptitle("Candidate rankings - cross-cluster medians/means (5 clusters)")
    fig.tight_layout()
    fig.savefig(out_path, dpi=130)
    plt.close(fig)


def save_bias_comparison_plot(out_path: Path, rows: list[dict]) -> None:
    """Cluster-by-candidate kappa bias matrix (heatmap)."""
    candidate_nums = sorted({r["candidate_number"] for r in rows})
    cluster_ids = sorted({r["cluster_id"] for r in rows},
                          key=lambda c: [cl["id"] for cl in CLUSTERS].index(c))
    M = np.full((len(candidate_nums), len(cluster_ids)), np.nan)
    for r in rows:
        i = candidate_nums.index(r["candidate_number"])
        j = cluster_ids.index(r["cluster_id"])
        M[i, j] = r["kappa_bias"]
    fig, ax = plt.subplots(figsize=(9, 7))
    vmax = float(np.nanmax(np.abs(M)))
    if not np.isfinite(vmax) or vmax == 0:
        vmax = 1.0
    im = ax.imshow(M, cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")
    ax.set_xticks(range(len(cluster_ids)))
    ax.set_xticklabels(cluster_ids, rotation=20, ha="right")
    ax.set_yticks(range(len(candidate_nums)))
    ax.set_yticklabels([f"C{n}" for n in candidate_nums])
    ax.set_xlabel("Cluster")
    ax.set_ylabel("Candidate")
    ax.set_title("Kappa bias (predicted - observed) per candidate per cluster")
    for i in range(len(candidate_nums)):
        for j in range(len(cluster_ids)):
            ax.text(j, i, f"{M[i, j]:+.3f}", ha="center", va="center",
                     fontsize=7, color="black")
    fig.colorbar(im, ax=ax, label="kappa bias")
    fig.tight_layout()
    fig.savefig(out_path, dpi=130)
    plt.close(fig)


def save_pearson_comparison_plot(out_path: Path, rows: list[dict]) -> None:
    """Two-panel: Pearson kappa and Pearson gamma by candidate, per cluster."""
    candidate_nums = sorted({r["candidate_number"] for r in rows})
    cluster_ids = sorted({r["cluster_id"] for r in rows},
                          key=lambda c: [cl["id"] for cl in CLUSTERS].index(c))
    fig, axes = plt.subplots(1, 2, figsize=(15, 6))
    for ax, metric_key, title in zip(
        axes, ["pearson_kappa", "pearson_gamma"],
        ["Pearson kappa", "Pearson gamma"],
    ):
        width = 0.13
        x_pos = np.arange(len(candidate_nums))
        for k, cid in enumerate(cluster_ids):
            vals = []
            for cn in candidate_nums:
                row = next((r for r in rows
                            if r["candidate_number"] == cn and r["cluster_id"] == cid), None)
                vals.append(row[metric_key] if row else 0.0)
            ax.bar(x_pos + (k - 2) * width, vals, width, label=cid)
        ax.set_xticks(x_pos)
        ax.set_xticklabels([f"C{n}" for n in candidate_nums])
        ax.set_xlabel("Candidate")
        ax.set_ylabel(title)
        ax.set_title(f"{title} by candidate and cluster")
        ax.axhline(0.0, color="black", lw=0.7, ls=":")
        ax.grid(axis="y", alpha=0.3)
        ax.legend(fontsize=8, loc="best")
    fig.suptitle("Per-cluster Pearson correlation across candidates")
    fig.tight_layout()
    fig.savefig(out_path, dpi=130)
    plt.close(fig)


def save_ssim_comparison_plot(out_path: Path, rows: list[dict]) -> None:
    """Two-panel: SSIM kappa and SSIM gamma by candidate, per cluster."""
    candidate_nums = sorted({r["candidate_number"] for r in rows})
    cluster_ids = sorted({r["cluster_id"] for r in rows},
                          key=lambda c: [cl["id"] for cl in CLUSTERS].index(c))
    fig, axes = plt.subplots(1, 2, figsize=(15, 6))
    for ax, metric_key, title in zip(
        axes, ["ssim_kappa", "ssim_gamma"],
        ["SSIM kappa", "SSIM gamma"],
    ):
        width = 0.13
        x_pos = np.arange(len(candidate_nums))
        for k, cid in enumerate(cluster_ids):
            vals = []
            for cn in candidate_nums:
                row = next((r for r in rows
                            if r["candidate_number"] == cn and r["cluster_id"] == cid), None)
                vals.append(row[metric_key] if row else 0.0)
            ax.bar(x_pos + (k - 2) * width, vals, width, label=cid)
        ax.set_xticks(x_pos)
        ax.set_xticklabels([f"C{n}" for n in candidate_nums])
        ax.set_xlabel("Candidate")
        ax.set_ylabel(title)
        ax.set_title(f"{title} by candidate and cluster")
        ax.axhline(0.0, color="black", lw=0.7, ls=":")
        ax.grid(axis="y", alpha=0.3)
        ax.legend(fontsize=8, loc="best")
    fig.suptitle("Per-cluster SSIM across candidates")
    fig.tight_layout()
    fig.savefig(out_path, dpi=130)
    plt.close(fig)


def save_cluster_performance_plot(out_path: Path, rows: list[dict]) -> None:
    """Cluster-level ranking across all candidates (mean Pearson kappa)."""
    cluster_ids = sorted({r["cluster_id"] for r in rows},
                          key=lambda c: [cl["id"] for cl in CLUSTERS].index(c))
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    for ax, key, title in zip(
        axes,
        ["pearson_kappa", "kappa_bias"],
        ["Mean Pearson kappa across candidates", "Mean kappa bias across candidates"],
    ):
        means = []
        stds = []
        for cid in cluster_ids:
            vals = [r[key] for r in rows if r["cluster_id"] == cid]
            arr = np.array(vals, dtype=np.float64)
            arr = arr[np.isfinite(arr)]
            means.append(float(np.mean(arr)) if arr.size else 0.0)
            stds.append(float(np.std(arr)) if arr.size else 0.0)
        ax.bar(range(len(cluster_ids)), means, yerr=stds, color="steelblue",
                edgecolor="black", capsize=4)
        ax.set_xticks(range(len(cluster_ids)))
        ax.set_xticklabels(cluster_ids, rotation=20, ha="right")
        ax.set_ylabel(key)
        ax.set_title(title)
        if key == "pearson_kappa":
            ax.axhline(0.0, color="black", lw=0.7, ls=":")
        ax.grid(axis="y", alpha=0.3)
    fig.suptitle("Cluster performance across all candidates")
    fig.tight_layout()
    fig.savefig(out_path, dpi=130)
    plt.close(fig)


def save_response_family_summary_plot(out_path: Path, cross_rows: list[dict]) -> None:
    """Group by family and show the best candidate per family."""
    families = sorted({r["candidate_family"] for r in cross_rows})
    family_best = []
    for fam in families:
        sub = [r for r in cross_rows if r["candidate_family"] == fam]
        best = max(sub, key=lambda r: r["median_pearson_kappa"])
        family_best.append(best)
    family_best.sort(key=lambda r: -r["median_pearson_kappa"])

    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    labels = [
        f"C{int(r['candidate_number'])}\n{r['candidate_family']}"
        for r in family_best
    ]
    pk = [r["median_pearson_kappa"] for r in family_best]
    bias = [r["mean_kappa_bias"] for r in family_best]
    x_pos = np.arange(len(family_best))
    axes[0].bar(x_pos, pk, color="steelblue", edgecolor="black")
    axes[0].set_xticks(x_pos)
    axes[0].set_xticklabels(labels, rotation=30, ha="right", fontsize=8)
    axes[0].set_ylabel("Median Pearson kappa")
    axes[0].set_title("Best candidate per family by median Pearson kappa")
    axes[0].axhline(0.0, color="black", lw=0.7, ls=":")
    axes[0].grid(axis="y", alpha=0.3)
    axes[1].bar(x_pos, bias, color="coral", edgecolor="black")
    axes[1].set_xticks(x_pos)
    axes[1].set_xticklabels(labels, rotation=30, ha="right", fontsize=8)
    axes[1].set_ylabel("Mean kappa bias")
    axes[1].set_title("Best candidate per family by mean kappa bias")
    axes[1].axhline(0.0, color="black", lw=0.7, ls=":")
    axes[1].grid(axis="y", alpha=0.3)
    fig.suptitle("Response family summary (best candidate per family)")
    fig.tight_layout()
    fig.savefig(out_path, dpi=130)
    plt.close(fig)


# =============================================================================
# Report writer
# =============================================================================
def _fmt(v, fmt=".4f"):
    if isinstance(v, float):
        if not np.isfinite(v):
            return "nan"
        return format(v, fmt)
    return str(v)


def write_report(out_root: Path, per_cluster_rows: list[dict],
                 cross_rows: list[dict], ranking_rows: list[dict],
                 hash_check: dict, total_seconds: float,
                 outcome: str, q_answers: list[str]) -> None:
    lines = ["# PBUF VERSION-B PHYSICS-LAB-001",
             "",
             "**Local response hypothesis survey using the frozen",
             "Version 1 weak-lensing laboratory (LAB-FREEZE-001).**",
             "",
             "No parameter fitting.  No optimisation.  Each candidate is",
             "tested exactly once on each of the five benchmark clusters",
             "with the frozen minimum production configuration.",
             "",
             "## Status",
             "",
             f"- Frozen hash verification: "
             f"**{'PASS' if hash_check['ok'] else 'FAIL'}**",
             f"- Total runtime: **{total_seconds:.1f} s**",
             f"- Candidates tested: **{len(CANDIDATES)}**",
             f"- Clusters: **{len(CLUSTERS)}**",
             "",
             "## Frozen laboratory",
             "",
             "The Version 1 laboratory is used as the measurement",
             "instrument without modification.",
             "",
             "| Component | Frozen specification |",
             "|---|---|",
             "| Constitutive | `C(X) = 0.18 * rho(X) / rho_max` (Version A) |",
             "| Transport | neighbour-to-neighbour, direct addition, |",
             "| | per-step unit-speed renormalisation |",
             "| Response direction | 90 deg transverse (R_90 of grad C) |",
             "| Source plane | Launch B (Cartesian 2D) |",
             "| Observable | Jacobian (ray-bundle linear fit per bin) |",
             "| Matter input | `rho = max(kappa_obs, 0) / max(max(kappa_obs, 0))` |",
             "",
             "## Production configuration",
             "",
             "| Parameter | Value |",
             "|---|---|",
             f"| Photons | {MIN_PRODUCTION['nphotons']:,} |",
             f"| Constitutive grid | {MIN_PRODUCTION['grid_n']}^2 |",
             f"| Step size | Delta s / 2 = {MIN_PRODUCTION['step']:.4f} |",
             f"| Number of steps | {MIN_PRODUCTION['steps']} |",
             f"| Source plane | Cartesian 2D (Launch B) |",
             f"| Observable | Jacobian |",
             "",
             "## Candidates",
             "",
             "| # | Name | Family | Description |",
             "|---|---|---|---|"]
    for c in CANDIDATES:
        lines.append(
            f"| {c.number} | {c.name} | {c.family} | {c.description} |"
        )
    lines += [
        "",
        "All fixed parameters are documented in the candidate source",
        "code.  No parameter is fitted.",
        "",
        "## Per-candidate, per-cluster metrics",
        "",
        "Computed metrics for every (candidate, cluster) pair:",
        "",
        "| Candidate | Cluster | RMS k | RMS g | Pearson k | Pearson g "
        "| SSIM k | SSIM g | k bias | g bias | conservation | runtime (s) |",
        "|---|---|---|---|---|---|---|---|---|---|---|---|",
    ]
    for r in per_cluster_rows:
        lines.append(
            f"| C{int(r['candidate_number'])} {r['candidate_name']} "
            f"| {r['cluster_label']} "
            f"| {_fmt(r['rms_kappa'], '.4e')} "
            f"| {_fmt(r['rms_gamma'], '.4e')} "
            f"| {_fmt(r['pearson_kappa'], '+.4f')} "
            f"| {_fmt(r['pearson_gamma'], '+.4f')} "
            f"| {_fmt(r['ssim_kappa'], '+.4f')} "
            f"| {_fmt(r['ssim_gamma'], '+.4f')} "
            f"| {_fmt(r['kappa_bias'], '+.4e')} "
            f"| {_fmt(r['gamma_bias'], '+.4e')} "
            f"| {_fmt(r['max_conservation_error'], '.3e')} "
            f"| {_fmt(r['runtime_seconds'], '.3f')} |"
        )

    lines += [
        "",
        "## Cross-cluster evaluation",
        "",
        "For every candidate the following medians/means are taken",
        "across the five benchmark clusters.",
        "",
        "| Candidate | Median Pearson k | Median Pearson g | Median SSIM k "
        "| Mean k Bias | Mean g Bias | Mean Pearson k | Conservation max "
        "| Runtime (s) |",
        "|---|---|---|---|---|---|---|---|---|",
    ]
    for r in cross_rows:
        lines.append(
            f"| C{int(r['candidate_number'])} {r['candidate_name']} "
            f"| {_fmt(r['median_pearson_kappa'], '+.4f')} "
            f"| {_fmt(r['median_pearson_gamma'], '+.4f')} "
            f"| {_fmt(r['median_ssim_kappa'], '+.4f')} "
            f"| {_fmt(r['mean_kappa_bias'], '+.4e')} "
            f"| {_fmt(r['mean_gamma_bias'], '+.4e')} "
            f"| {_fmt(r['mean_pearson_kappa'], '+.4f')} "
            f"| {_fmt(r['max_conservation_error'], '.3e')} "
            f"| {_fmt(r['median_runtime_seconds'], '.3f')} |"
        )

    lines += [
        "",
        "## Candidate ranking (by median Pearson kappa)",
        "",
        "| Rank | Candidate | Family | Median Pearson k | Median SSIM k "
        "| Mean k Bias | Delta Pearson vs control "
        "| Delta SSIM vs control |",
        "|---|---|---|---|---|---|---|---|",
    ]
    for r in ranking_rows:
        lines.append(
            f"| {int(r['rank'])} | C{int(r['candidate_number'])} "
            f"{r['candidate_name']} | {r['candidate_family']} "
            f"| {_fmt(r['median_pearson_kappa'], '+.4f')} "
            f"| {_fmt(r['median_ssim_kappa'], '+.4f')} "
            f"| {_fmt(r['mean_kappa_bias'], '+.4e')} "
            f"| {_fmt(r['improvement_vs_control_pearson_kappa'], '+.4f')} "
            f"| {_fmt(r['improvement_vs_control_ssim_kappa'], '+.4f')} |"
        )

    lines += [
        "",
        "## Required questions",
        "",
    ]
    lines += q_answers

    lines += [
        "",
        "## Outcome determination",
        "",
        outcome,
        "",
        "## Numerical stability report",
        "",
        "| Candidate | Mean runtime (s) | Max conservation |",
        "|---|---|---|",
    ]
    for r in cross_rows:
        lines.append(
            f"| C{int(r['candidate_number'])} {r['candidate_name']} "
            f"| {_fmt(r['median_runtime_seconds'], '.3f')} "
            f"| {_fmt(r['max_conservation_error'], '.3e')} |"
        )

    lines += [
        "",
        "## Top-level artefacts",
        "",
        "- `runs/version_b_physics_lab001/report.md` (this file)",
        "- `runs/version_b_physics_lab001/candidate_summary.csv`",
        "- `runs/version_b_physics_lab001/cross_cluster_statistics.csv`",
        "- `runs/version_b_physics_lab001/candidate_ranking.csv`",
        "- `runs/version_b_physics_lab001/run.json`",
        "- `runs/version_b_physics_lab001/validation.json`",
        "- `runs/version_b_physics_lab001/plots/candidate_rankings.png`",
        "- `runs/version_b_physics_lab001/plots/bias_comparison.png`",
        "- `runs/version_b_physics_lab001/plots/pearson_comparison.png`",
        "- `runs/version_b_physics_lab001/plots/ssim_comparison.png`",
        "- `runs/version_b_physics_lab001/plots/cluster_performance.png`",
        "- `runs/version_b_physics_lab001/plots/response_family_summary.png`",
        "",
        f"**Total execution time:** {total_seconds:.1f} s.",
        "",
    ]
    (out_root / "report.md").write_text("\n".join(lines))


def answer_q1(cross_rows: list[dict], per_cluster_rows: list[dict]) -> list[str]:
    """Q1: Does any candidate reduce the systematic kappa underprediction?"""
    control = next(r for r in cross_rows if r["candidate_number"] == 1)
    others = [r for r in cross_rows if r["candidate_number"] != 1]
    improvements = []
    for r in others:
        delta = r["mean_kappa_bias"] - control["mean_kappa_bias"]
        improvements.append((r["candidate_number"], r["candidate_name"],
                              r["mean_kappa_bias"], delta))
    improvements.sort(key=lambda t: -t[3])  # least negative improvement first
    best = improvements[0]
    n_improve = sum(1 for _, _, _, delta in improvements if delta > 0)
    lines = [
        "### Q1. Reduction of the systematic kappa underprediction",
        "",
        f"Control (C1) mean kappa bias = "
        f"{control['mean_kappa_bias']:+.5f}.  A negative bias means the",
        "frozen Version A laboratory underpredicts kappa.  A reduction",
        "of this bias (toward zero or positive) is the success criterion.",
        "",
        "| Candidate | Mean kappa bias | Delta vs control |",
        "|---|---|---|",
    ]
    for cn, nm, bias, delta in improvements:
        lines.append(
            f"| C{int(cn)} {nm} | {bias:+.5f} | {delta:+.5f} |"
        )
    lines += [
        "",
        f"Best improvement: C{int(best[0])} {best[1]} "
        f"(mean kappa bias = {best[2]:+.5f}, delta = {best[3]:+.5f}).",
        f"Number of candidates that *reduce* the magnitude of the bias: "
        f"{n_improve}/{len(improvements)}.",
        "",
    ]
    return lines


def answer_q2(cross_rows: list[dict], ranking_rows: list[dict]) -> list[str]:
    """Q2: Which candidate produces the largest improvement?"""
    if not ranking_rows:
        return ["### Q2. Largest improvement",
                "",
                "No candidates available.", ""]
    # Excluding control, find the largest improvement in median Pearson kappa.
    control_pk = next(r["median_pearson_kappa"] for r in cross_rows
                      if r["candidate_number"] == 1)
    candidates = [r for r in cross_rows if r["candidate_number"] != 1]
    if not candidates:
        return ["### Q2. Largest improvement",
                "",
                "No alternative candidates.", ""]
    best_pk = max(candidates, key=lambda r: r["median_pearson_kappa"])
    delta_pk = best_pk["median_pearson_kappa"] - control_pk
    best_ssim = max(candidates, key=lambda r: r["median_ssim_kappa"])
    control_ssim = next(r["median_ssim_kappa"] for r in cross_rows
                        if r["candidate_number"] == 1)
    delta_ssim = best_ssim["median_ssim_kappa"] - control_ssim
    best_bias = max(candidates, key=lambda r: -abs(r["mean_kappa_bias"]))
    control_bias = next(r["mean_kappa_bias"] for r in cross_rows
                        if r["candidate_number"] == 1)
    delta_bias = best_bias["mean_kappa_bias"] - control_bias

    lines = [
        "### Q2. Largest improvement",
        "",
        f"Control (C1) median Pearson kappa = {control_pk:+.4f}, "
        f"median SSIM kappa = {control_ssim:+.4f}, "
        f"mean kappa bias = {control_bias:+.5f}.",
        "",
        f"**Largest median Pearson kappa:** C{int(best_pk['candidate_number'])} "
        f"{best_pk['candidate_name']} ({best_pk['median_pearson_kappa']:+.4f}, "
        f"delta = {delta_pk:+.4f}).",
        f"**Largest median SSIM kappa:** C{int(best_ssim['candidate_number'])} "
        f"{best_ssim['candidate_name']} ({best_ssim['median_ssim_kappa']:+.4f}, "
        f"delta = {delta_ssim:+.4f}).",
        f"**Lowest |mean kappa bias|:** C{int(best_bias['candidate_number'])} "
        f"{best_bias['candidate_name']} ({best_bias['mean_kappa_bias']:+.5f}, "
        f"delta = {delta_bias:+.5f}).",
        "",
        "Top-3 by median Pearson kappa:",
        "",
        "| Rank | Candidate | Median Pearson kappa | Median SSIM kappa | Mean kappa bias |",
        "|---|---|---|---|---|",
    ]
    for r in ranking_rows[:3]:
        if r["candidate_number"] == 1:
            continue
        lines.append(
            f"| {int(r['rank'])} | C{int(r['candidate_number'])} "
            f"{r['candidate_name']} | {r['median_pearson_kappa']:+.4f} "
            f"| {r['median_ssim_kappa']:+.4f} "
            f"| {r['mean_kappa_bias']:+.5f} |"
        )
    lines.append("")
    return lines


def answer_q3(cross_rows: list[dict]) -> list[str]:
    """Q3: Does any candidate worsen numerical stability?"""
    control = next(r for r in cross_rows if r["candidate_number"] == 1)
    control_cons = control["max_conservation_error"]
    control_rt = control["median_runtime_seconds"]
    worsened = []
    for r in cross_rows:
        if r["candidate_number"] == 1:
            continue
        if r["max_conservation_error"] > control_cons + 1e-19:
            worsened.append(("conservation",
                              r["candidate_name"], r["max_conservation_error"],
                              r["max_conservation_error"] - control_cons))
    lines = [
        "### Q3. Numerical stability impact",
        "",
        f"Control (C1) maximum conservation error = "
        f"{control_cons:.3e}, median runtime = {control_rt:.3f} s.",
        "Machine epsilon = 2.220446049250313e-16.",
        "",
        "| Candidate | Max conservation | Runtime (s) | Status |",
        "|---|---|---|---|",
    ]
    n_stable = 0
    for r in cross_rows:
        ok = r["max_conservation_error"] <= control_cons + 1e-19
        n_stable += int(ok)
        lines.append(
            f"| C{int(r['candidate_number'])} {r['candidate_name']} "
            f"| {r['max_conservation_error']:.3e} "
            f"| {r['median_runtime_seconds']:.3f} "
            f"| {'machine-epsilon' if ok else 'EXCEEDS EPSILON'} |"
        )
    lines += [
        "",
        f"Number of candidates that preserve the machine-precision "
        f"conservation bound: {n_stable}/{len(cross_rows)}.",
        "",
    ]
    if not worsened:
        lines.append("No candidate exceeds the machine-precision conservation "
                      "bound observed for the frozen control.")
    else:
        lines.append("Candidates exceeding the bound: " +
                      ", ".join(n for _, n, _, _ in worsened) + ".")
    lines.append("")
    return lines


def answer_q4(cross_rows: list[dict]) -> list[str]:
    """Q4: Does any candidate preserve machine-precision conservation?"""
    eps = 2.220446049250313e-16
    preserved = []
    for r in cross_rows:
        if r["max_conservation_error"] <= eps + 1e-30:
            preserved.append((r["candidate_number"], r["candidate_name"],
                              r["max_conservation_error"]))
    lines = [
        "### Q4. Machine-precision conservation preservation",
        "",
        f"Machine epsilon = {eps:.3e}.",
        f"Number of candidates that satisfy `max conservation <= machine epsilon`: "
        f"{len(preserved)}/{len(cross_rows)}.",
        "",
        "| Candidate | Max conservation | Preserves |",
        "|---|---|---|",
    ]
    for cn, nm, mc in preserved:
        lines.append(f"| C{int(cn)} {nm} | {mc:.3e} | YES |")
    not_preserved = [r for r in cross_rows
                     if r["max_conservation_error"] > eps + 1e-30]
    for r in not_preserved:
        lines.append(
            f"| C{int(r['candidate_number'])} {r['candidate_name']} "
            f"| {r['max_conservation_error']:.3e} | NO |"
        )
    lines.append("")
    return lines


def answer_q5(per_cluster_rows: list[dict],
               cross_rows: list[dict]) -> list[str]:
    """Q5: Is the improvement consistent across all benchmark clusters?"""
    control = next(r for r in cross_rows if r["candidate_number"] == 1)
    others = [r for r in cross_rows if r["candidate_number"] != 1]
    consistency_rows = []
    for cand in others:
        cn = cand["candidate_number"]
        cn_name = cand["candidate_name"]
        # Cluster-by-cluster Pearson kappa delta.
        per_cluster_delta = []
        for cl in CLUSTERS:
            cid = cl["id"]
            ctrl_row = next((r for r in per_cluster_rows
                             if r["candidate_number"] == 1
                             and r["cluster_id"] == cid), None)
            cand_row = next((r for r in per_cluster_rows
                             if r["candidate_number"] == cn
                             and r["cluster_id"] == cid), None)
            if ctrl_row is not None and cand_row is not None:
                delta = cand_row["pearson_kappa"] - ctrl_row["pearson_kappa"]
                per_cluster_delta.append((cid, cand_row["pearson_kappa"],
                                            ctrl_row["pearson_kappa"], delta))
        if not per_cluster_delta:
            continue
        deltas = np.array([d for _, _, _, d in per_cluster_delta])
        n_positive = int(np.sum(deltas > 0))
        n_negative = int(np.sum(deltas < 0))
        consistency_rows.append({
            "candidate_number": cn,
            "candidate_name": cn_name,
            "n_positive": n_positive,
            "n_negative": n_negative,
            "median_delta": float(np.median(deltas)),
            "sign_consistent": n_positive == len(per_cluster_delta),
            "per_cluster": per_cluster_delta,
        })
    consistency_rows.sort(key=lambda r: -r["median_delta"])
    lines = [
        "### Q5. Cross-cluster consistency of improvement",
        "",
        "Improvement is consistent across clusters if every cluster",
        "exhibits a positive delta in Pearson kappa versus the control.",
        "",
        "| Candidate | Median delta Pearson k | # clusters +ve | # clusters -ve | Sign consistent |",
        "|---|---|---|---|---|",
    ]
    for r in consistency_rows:
        lines.append(
            f"| C{int(r['candidate_number'])} {r['candidate_name']} "
            f"| {r['median_delta']:+.4f} "
            f"| {r['n_positive']} "
            f"| {r['n_negative']} "
            f"| {'YES' if r['sign_consistent'] else 'NO'} |"
        )
    consistent = [r for r in consistency_rows if r["sign_consistent"]]
    lines += [
        "",
        f"Number of candidates with sign-consistent improvement on all "
        f"5 clusters: {len(consistent)}/{len(consistency_rows)}.",
        "",
    ]
    return lines


def answer_q6(cross_rows: list[dict]) -> list[str]:
    """Q6: Which physical family performs best?"""
    families = sorted({r["candidate_family"] for r in cross_rows})
    family_summary = []
    for fam in families:
        sub = [r for r in cross_rows if r["candidate_family"] == fam]
        best = max(sub, key=lambda r: r["median_pearson_kappa"])
        family_summary.append({
            "family": fam,
            "n_candidates": len(sub),
            "best_candidate_number": best["candidate_number"],
            "best_candidate_name": best["candidate_name"],
            "best_median_pearson_kappa": best["median_pearson_kappa"],
            "best_median_ssim_kappa": best["median_ssim_kappa"],
            "best_mean_kappa_bias": best["mean_kappa_bias"],
            "family_median_pearson_kappa": float(np.median(
                [r["median_pearson_kappa"] for r in sub])),
        })
    family_summary.sort(key=lambda r: -r["best_median_pearson_kappa"])
    lines = [
        "### Q6. Best-performing physical family",
        "",
        "| Rank | Family | Best candidate | Median Pearson k | Median SSIM k | Mean k bias |",
        "|---|---|---|---|---|---|",
    ]
    for i, r in enumerate(family_summary, start=1):
        lines.append(
            f"| {i} | {r['family']} "
            f"| C{int(r['best_candidate_number'])} {r['best_candidate_name']} "
            f"| {r['best_median_pearson_kappa']:+.4f} "
            f"| {r['best_median_ssim_kappa']:+.4f} "
            f"| {r['best_mean_kappa_bias']:+.5f} |"
        )
    if family_summary:
        best = family_summary[0]
        lines += [
            "",
            f"**Best family:** `{best['family']}` "
            f"(best candidate: C{int(best['best_candidate_number'])} "
            f"{best['best_candidate_name']}, "
            f"median Pearson kappa = {best['best_median_pearson_kappa']:+.4f}).",
            "",
        ]
    return lines


def determine_outcome(cross_rows: list[dict],
                       consistency_rows: list[dict]) -> str:
    """Outcome A: at least one candidate improves consistently across all
    5 clusters while preserving numerical stability.  Otherwise Outcome B.
    """
    eps = 2.220446049250313e-16
    stable_candidates = [r for r in cross_rows
                          if r["max_conservation_error"] <= eps + 1e-30]
    stable_numbers = {int(r["candidate_number"]) for r in stable_candidates}

    consistent = [r for r in consistency_rows
                   if r["sign_consistent"] and r["median_delta"] > 0
                   and int(r["candidate_number"]) in stable_numbers]
    control = next(r for r in cross_rows if r["candidate_number"] == 1)
    consistent_improving = [r for r in consistent
                             if r["median_delta"] > 0.01]
    if consistent:
        best = max(consistent, key=lambda r: r["median_delta"])
        if consistent_improving:
            return (f"**Outcome A** - C{int(best['candidate_number'])} "
                    f"{best['candidate_name']} (family: "
                    f"{next(c['candidate_family'] for c in CANDIDATES if c.number == best['candidate_number'])}) "
                    f"improves Pearson kappa on all 5 clusters with "
                    f"machine-precision conservation "
                    f"(median delta = {best['median_delta']:+.4f}, "
                    f"control median Pearson kappa = "
                    f"{control['median_pearson_kappa']:+.4f}). "
                    f"Number of consistently-improving stable candidates: "
                    f"{len(consistent_improving)}/{len(cross_rows)-1}. "
                    f"The frozen laboratory is not modified.")
        return (f"**Outcome A (weak)** - "
                f"{len(consistent)} candidate(s) show sign-consistent "
                f"improvement on all 5 clusters while preserving "
                f"machine-precision conservation, but the median delta "
                f"is below 0.01 in absolute value.  No candidate "
                f"produces a strong, consistent improvement.  The "
                f"frozen laboratory is not modified; the result "
                f"suggests the missing physics lies outside the "
                f"tested local-response families.")
    return (f"**Outcome B** - no candidate produces sign-consistent "
            f"improvement on all 5 benchmark clusters while preserving "
            f"machine-precision conservation.  The frozen laboratory "
            f"is not modified.  The missing physics likely lies "
            f"outside the tested local-response families "
            f"(gradient, neighbour coherence, cooperative response, "
            f"elastic memory, gradient curvature, phase coherence, "
            f"relaxation, weak-gradient enhancement, constitutive "
            f"coupling, combined response).")


# =============================================================================
# Main
# =============================================================================
def main():
    out_root = DEFAULT_OUT
    out_root.mkdir(parents=True, exist_ok=True)
    PLOTS.mkdir(parents=True, exist_ok=True)

    overall_started = time.perf_counter()

    print("=" * 72)
    print("PBUF VERSION-B PHYSICS-LAB-001")
    print("Local response hypothesis survey - frozen Version 1 laboratory")
    print("=" * 72)

    print("\n[1/4] Verifying frozen source hashes against LAB-FREEZE-001 ...")
    hash_check = verify_frozen_hashes()
    for name, info in hash_check["files"].items():
        marker = "OK" if info["match"] else "MISMATCH"
        print(f"  [{marker}] {name}: {info['actual_sha256']}")
    if not hash_check["ok"]:
        raise RuntimeError("Frozen source file hashes do not match LAB-FREEZE-001.")

    print("\n[2/4] Running all 10 candidates on all 5 clusters ...")
    per_cluster_rows: list[dict] = []
    for cand in CANDIDATES:
        print(f"\n  Candidate C{cand.number}: {cand.name} ({cand.family})", flush=True)
        for cluster in CLUSTERS:
            t0 = time.perf_counter()
            try:
                result = run_candidate_on_cluster(cluster, MIN_PRODUCTION, cand)
            except Exception as exc:
                print(f"    ERROR on {cluster['id']}: {exc}\n{traceback.format_exc()}")
                raise
            row = {
                "candidate_number": cand.number,
                "candidate_name": cand.name,
                "candidate_family": cand.family,
                "cluster_id": cluster["id"],
                "cluster_label": cluster["label"],
                "rms_kappa": result["comparison"]["kappa"]["rms_error"],
                "rms_gamma": result["comparison"]["gamma_mag"]["rms_error"],
                "pearson_kappa": result["comparison"]["kappa"]["pearson_correlation"],
                "pearson_gamma": result["comparison"]["gamma_mag"]["pearson_correlation"],
                "ssim_kappa": result["comparison"]["kappa"]["ssim"],
                "ssim_gamma": result["comparison"]["gamma_mag"]["ssim"],
                "kappa_bias": result["bias_kappa"],
                "gamma_bias": result["bias_gamma"],
                "std_resid_kappa": result["std_resid_kappa"],
                "std_resid_gamma": result["std_resid_gamma"],
                "max_conservation_error": result["max_conservation_error"],
                "runtime_seconds": result["total_seconds"],
                "n_finite_pixels_kappa": result["n_finite_pixels_kappa"],
                "n_finite_pixels_gamma": result["n_finite_pixels_gamma"],
                "n_photons": result["n_photons"],
                "grid_n": result["grid_n"],
                "step": result["step"],
                "steps": result["steps"],
            }
            per_cluster_rows.append(row)
            elapsed = time.perf_counter() - t0
            print(f"    -> {cluster['label']:12s}: Pearson k = {row['pearson_kappa']:+.4f}, "
                  f"k bias = {row['kappa_bias']:+.5f}, "
                  f"runtime = {row['runtime_seconds']:.3f} s, "
                  f"cons = {row['max_conservation_error']:.2e}  ({elapsed:.1f}s)",
                  flush=True)

    print("\n[3/4] Computing cross-cluster statistics ...")
    cross_rows = []
    for cand in CANDIDATES:
        sub = [r for r in per_cluster_rows
               if r["candidate_number"] == cand.number]
        pk_vals = [r["pearson_kappa"] for r in sub]
        pg_vals = [r["pearson_gamma"] for r in sub]
        sk_vals = [r["ssim_kappa"] for r in sub]
        sg_vals = [r["ssim_gamma"] for r in sub]
        rk_vals = [r["rms_kappa"] for r in sub]
        rg_vals = [r["rms_gamma"] for r in sub]
        bias_k_vals = [r["kappa_bias"] for r in sub]
        bias_g_vals = [r["gamma_bias"] for r in sub]
        cons_vals = [r["max_conservation_error"] for r in sub]
        rt_vals = [r["runtime_seconds"] for r in sub]
        cross_rows.append({
            "candidate_number": cand.number,
            "candidate_name": cand.name,
            "candidate_family": cand.family,
            "median_pearson_kappa": _safe_median(pk_vals),
            "median_pearson_gamma": _safe_median(pg_vals),
            "median_ssim_kappa": _safe_median(sk_vals),
            "median_ssim_gamma": _safe_median(sg_vals),
            "median_rms_kappa": _safe_median(rk_vals),
            "median_rms_gamma": _safe_median(rg_vals),
            "mean_kappa_bias": _safe_mean(bias_k_vals),
            "std_kappa_bias": _safe_std(bias_k_vals),
            "mean_gamma_bias": _safe_mean(bias_g_vals),
            "std_gamma_bias": _safe_std(bias_g_vals),
            "mean_pearson_kappa": _safe_mean(pk_vals),
            "mean_pearson_gamma": _safe_mean(pg_vals),
            "max_conservation_error": float(np.max(cons_vals)),
            "median_runtime_seconds": _safe_median(rt_vals),
            "n_clusters_with_positive_pearson_kappa": int(sum(
                1 for v in pk_vals if v > 0)),
            "n_clusters_with_positive_pearson_gamma": int(sum(
                1 for v in pg_vals if v > 0)),
        })

    control_pk = next(r["median_pearson_kappa"] for r in cross_rows
                      if r["candidate_number"] == 1)
    control_ssim = next(r["median_ssim_kappa"] for r in cross_rows
                        if r["candidate_number"] == 1)
    control_bias = next(r["mean_kappa_bias"] for r in cross_rows
                        if r["candidate_number"] == 1)

    ranked = sorted(cross_rows,
                    key=lambda r: -r["median_pearson_kappa"])
    ranking_rows = []
    for rank, r in enumerate(ranked, start=1):
        ranking_rows.append({
            "rank": rank,
            "candidate_number": r["candidate_number"],
            "candidate_name": r["candidate_name"],
            "candidate_family": r["candidate_family"],
            "median_pearson_kappa": r["median_pearson_kappa"],
            "median_pearson_gamma": r["median_pearson_gamma"],
            "median_ssim_kappa": r["median_ssim_kappa"],
            "mean_kappa_bias": r["mean_kappa_bias"],
            "improvement_vs_control_pearson_kappa":
                r["median_pearson_kappa"] - control_pk,
            "improvement_vs_control_ssim_kappa":
                r["median_ssim_kappa"] - control_ssim,
            "improvement_vs_control_kappa_bias":
                r["mean_kappa_bias"] - control_bias,
        })

    # Answer Q5 requires the per-cluster delta computation; pre-compute here.
    consistency_rows = []
    for cand in CANDIDATES:
        if cand.number == 1:
            continue
        per_cluster_delta = []
        for cl in CLUSTERS:
            cid = cl["id"]
            ctrl_row = next((r for r in per_cluster_rows
                             if r["candidate_number"] == 1
                             and r["cluster_id"] == cid), None)
            cand_row = next((r for r in per_cluster_rows
                             if r["candidate_number"] == cand.number
                             and r["cluster_id"] == cid), None)
            if ctrl_row is not None and cand_row is not None:
                delta = cand_row["pearson_kappa"] - ctrl_row["pearson_kappa"]
                per_cluster_delta.append((cid, cand_row["pearson_kappa"],
                                            ctrl_row["pearson_kappa"], delta))
        deltas = np.array([d for _, _, _, d in per_cluster_delta])
        consistency_rows.append({
            "candidate_number": cand.number,
            "candidate_name": cand.name,
            "median_delta": float(np.median(deltas)),
            "sign_consistent": bool(np.all(deltas > 0)),
            "n_positive": int(np.sum(deltas > 0)),
            "n_negative": int(np.sum(deltas < 0)),
        })

    q_answers = []
    q_answers += answer_q1(cross_rows, per_cluster_rows)
    q_answers += answer_q2(cross_rows, ranking_rows)
    q_answers += answer_q3(cross_rows)
    q_answers += answer_q4(cross_rows)
    q_answers += answer_q5(per_cluster_rows, cross_rows)
    q_answers += answer_q6(cross_rows)

    outcome = determine_outcome(cross_rows, consistency_rows)

    print("\n[4/4] Writing outputs ...")
    write_candidate_summary_csv(out_root, per_cluster_rows)
    write_cross_cluster_csv(out_root, cross_rows)
    write_candidate_ranking_csv(out_root, ranking_rows)

    save_candidate_rankings_plot(PLOTS / "candidate_rankings.png", cross_rows)
    save_bias_comparison_plot(PLOTS / "bias_comparison.png", per_cluster_rows)
    save_pearson_comparison_plot(PLOTS / "pearson_comparison.png",
                                   per_cluster_rows)
    save_ssim_comparison_plot(PLOTS / "ssim_comparison.png",
                                per_cluster_rows)
    save_cluster_performance_plot(PLOTS / "cluster_performance.png",
                                    per_cluster_rows)
    save_response_family_summary_plot(PLOTS / "response_family_summary.png",
                                        cross_rows)

    total_seconds = time.perf_counter() - overall_started

    run_doc = {
        "milestone": "PBUF VERSION-B PHYSICS-LAB-001",
        "kind": "local_response_hypothesis_survey",
        "frozen_laboratory": "Version 1 weak-lensing laboratory (LAB-FREEZE-001)",
        "frozen_implementation_sha256": hash_check["files"],
        "production_minimum": MIN_PRODUCTION,
        "candidates": [
            {"number": c.number, "name": c.name, "family": c.family,
             "description": c.description, "notes": c.notes}
            for c in CANDIDATES
        ],
        "clusters": [{"id": c["id"], "label": c["label"],
                       "directory": c["directory"]} for c in CLUSTERS],
        "fitting_performed": False,
        "optimisation_performed": False,
        "cosmological_bridges_introduced": False,
        "execution_seconds_total": float(total_seconds),
    }
    (out_root / "run.json").write_text(json.dumps(run_doc, indent=2))

    cons_eps = 2.220446049250313e-16
    n_stable = sum(1 for r in cross_rows
                   if r["max_conservation_error"] <= cons_eps + 1e-30)
    validation_doc = {
        "milestone": "PBUF VERSION-B PHYSICS-LAB-001",
        "frozen_hash_verification_passed": hash_check["ok"],
        "frozen_hashes": hash_check["files"],
        "all_clusters_machine_precision_conservation": all(
            r["max_conservation_error"] <= cons_eps + 1e-30
            for r in per_cluster_rows
        ),
        "candidates_preserving_conservation": n_stable,
        "candidates_total": len(cross_rows),
        "validation_passed": hash_check["ok"],
    }
    (out_root / "validation.json").write_text(json.dumps(validation_doc, indent=2))

    write_report(out_root, per_cluster_rows, cross_rows, ranking_rows,
                 hash_check, total_seconds, outcome, q_answers)

    print(f"\nVERSION-B PHYSICS-LAB-001 COMPLETE  ({total_seconds:.1f} s)")
    print(json.dumps({
        "milestone": "PBUF VERSION-B PHYSICS-LAB-001",
        "status": "OK",
        "candidates": len(CANDIDATES),
        "clusters": len(CLUSTERS),
        "frozen_hashes_ok": hash_check["ok"],
        "candidates_preserving_conservation": n_stable,
        "output": str(out_root),
        "execution_seconds": float(total_seconds),
    }, indent=2))


if __name__ == "__main__":
    main()
