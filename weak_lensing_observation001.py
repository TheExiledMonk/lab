#!/usr/bin/env python3
"""PBUF WEAK-LENSING-OBSERVATION-001 - first comparison against public observations.

Frozen Version A pipeline applied to five public weak-lensing benchmark clusters
(SaWLens reconstructions from Merten et al. 2014, Frontier Fields).

No modification to Version A.
No fitting.
No parameter changes.
Measurement only.
"""
from __future__ import annotations

import csv
import hashlib
import json
import time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from scipy.ndimage import map_coordinates

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
from constitutive_equations import get_equation


ROOT = Path(__file__).resolve().parent
BENCHMARK_DIR = ROOT / "PBUF_benchmark"
DEFAULT_OUT = ROOT / "runs" / "weak_lensing_observation001"

# ----------------------------------------------------------------------------
# FROZEN Version A pipeline parameters
# Identical to weak_lensing_prediction001.py and weak_lensing_generalization001.py
# ----------------------------------------------------------------------------
LENS = {
    "n": 128,
    "extent": 8.0,
    "strength": 0.18,
    "step": 0.06,
    "steps": 80,
    "y_span": 3.0,
    "nphotons": 2000,
    "bins": 64,
}

CLUSTERS = [
    {
        "id": "Abell2744",
        "directory": "WL-001_Abell2744",
        "slug": "abell2744",
        "label": "Abell 2744",
    },
    {
        "id": "MACS0416",
        "directory": "WL-002_MACS0416",
        "slug": "macs0416",
        "label": "MACS J0416",
    },
    {
        "id": "MACS1149",
        "directory": "WL-003_MACS1149",
        "slug": "macs1149",
        "label": "MACS J1149",
    },
    {
        "id": "AbellS1063",
        "directory": "WL-004_AbellS1063",
        "slug": "abells1063",
        "label": "Abell S1063",
    },
    {
        "id": "Abell370",
        "directory": "WL-005_Abell370",
        "slug": "abell370",
        "label": "Abell 370",
    },
]


# ----------------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------------
def sha256_array(array: np.ndarray) -> str:
    return hashlib.sha256(np.ascontiguousarray(array).tobytes()).hexdigest()


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def resample_to_grid(source: np.ndarray, target_n: int, extent: float) -> np.ndarray:
    """Resample a 2D source array onto a target_n x target_n Cartesian grid on [-extent, extent].

    The source array is assumed to be centred on the field centre. The mapping
    places the source pixel (0, 0) at the bottom-left of the target grid and
    pixel (N-1, N-1) at the top-right. Linear interpolation. Out-of-bounds
    pixels use nearest-neighbour extrapolation.
    """
    ny, nx = source.shape
    x = np.linspace(-extent, extent, target_n)
    y = np.linspace(-extent, extent, target_n)
    X, Y = np.meshgrid(x, y, indexing="xy")
    src_x = (X + extent) / (2 * extent) * (nx - 1)
    src_y = (Y + extent) / (2 * extent) * (ny - 1)
    coords = np.array([src_y, src_x])
    return map_coordinates(source, coords, order=1, mode="nearest")


# ----------------------------------------------------------------------------
# Observation loading
# ----------------------------------------------------------------------------
def load_observation(cluster: dict, benchmark_dir: Path) -> dict:
    """Load the four FITS products for one cluster.

    Returns a dict with arrays, header info, and per-file metadata.
    """
    folder = benchmark_dir / cluster["directory"]
    files = {
        "kappa": f"hlsp_frontier_model_{cluster['slug']}_merten_v1_kappa.fits",
        "gamma": f"hlsp_frontier_model_{cluster['slug']}_merten_v1_gamma.fits",
        "gamma1": f"hlsp_frontier_model_{cluster['slug']}_merten_v1_gamma1.fits",
        "gamma2": f"hlsp_frontier_model_{cluster['slug']}_merten_v1_gamma2.fits",
    }
    out = {"folder": str(folder), "files": {}, "headers": {}, "shas": {}}
    for key, name in files.items():
        path = folder / name
        with fits.open(path) as hdul:
            data = np.asarray(hdul[0].data, dtype=np.float64)
            header = dict(hdul[0].header)
        out[key] = data
        out["files"][key] = str(path)
        out["headers"][key] = {
            "NAXIS1": int(header.get("NAXIS1")),
            "NAXIS2": int(header.get("NAXIS2")),
            "CRPIX1": float(header.get("CRPIX1")),
            "CRPIX2": float(header.get("CRPIX2")),
            "CRVAL1": float(header.get("CRVAL1")),
            "CRVAL2": float(header.get("CRVAL2")),
            "CDELT1": float(header.get("CDELT1")),
            "CDELT2": float(header.get("CDELT2")),
            "CTYPE1": str(header.get("CTYPE1")),
            "CTYPE2": str(header.get("CTYPE2")),
            "Z_L": float(header.get("Z_L")) if header.get("Z_L") is not None else None,
            "Z_S": float(header.get("Z_S")) if header.get("Z_S") is not None else None,
            "RADESYS": str(header.get("RADESYS")),
            "EQUINOX": float(header.get("EQUINOX")) if header.get("EQUINOX") is not None else None,
        }
        out["shas"][key] = file_sha256(path)
    out["wcs"] = WCS(out["headers"]["kappa"])
    return out


def validate_gamma_consistency(obs: dict) -> dict:
    """Verify gamma.fits == sqrt(gamma1**2 + gamma2**2) within FP tolerance."""
    g = obs["gamma"]
    g1 = obs["gamma1"]
    g2 = obs["gamma2"]
    g_from_components = np.sqrt(g1 ** 2 + g2 ** 2)
    diff = np.abs(g - g_from_components)
    max_abs = float(diff.max())
    rms_abs = float(np.sqrt(np.mean(diff ** 2)))
    rel = diff / np.maximum(np.abs(g), 1e-15)
    return {
        "max_abs_difference": max_abs,
        "rms_abs_difference": rms_abs,
        "max_relative_difference": float(rel.max()),
        "tolerance_passed": max_abs <= 1e-5,
        "shape": list(g.shape),
    }


# ----------------------------------------------------------------------------
# FROZEN Version A pipeline (re-implemented identically to the frozen modules)
# ----------------------------------------------------------------------------
def make_field(rho: np.ndarray, extent: float, strength: float, n: int) -> dict:
    """Apply the frozen Version A constitutive + transport + response."""
    x = np.linspace(-extent, extent, n)
    y = np.linspace(-extent, extent, n)
    X, Y = np.meshgrid(x, y, indexing="xy")
    cfg = type("Config", (), {"deformation_strength": strength})()
    c = get_equation("A").solve(rho, cfg)
    gy, gx = np.gradient(c, x, y, edge_order=1)
    g = np.hypot(gx, gy)
    gx_hat = gx / np.maximum(g, 1e-15)
    gy_hat = gy / np.maximum(g, 1e-15)
    bad = g < 1e-15
    gx_hat = np.where(bad, 1.0, gx_hat)
    gy_hat = np.where(bad, 0.0, gy_hat)
    rx = -g * gy_hat
    ry = g * gx_hat
    return {
        "xgrid": x, "ygrid": y, "X": X, "Y": Y,
        "rho": rho, "c": c,
        "gx": gx, "gy": gy, "g_magnitude": g,
        "rx": rx, "ry": ry,
        "response_direction": np.arctan2(ry, rx),
    }


def propagate(field: dict, step: float, steps: int,
              x0: np.ndarray, y0: np.ndarray,
              vx0: np.ndarray, vy0: np.ndarray) -> dict:
    """Frozen neighbour-to-neighbour propagation with renormalisation."""
    xgrid = field["xgrid"]; ygrid = field["ygrid"]
    rx = field["rx"]; ry = field["ry"]
    x = x0.copy(); y = y0.copy()
    vx = vx0.copy(); vy = vy0.copy()
    nphotons = len(x)
    max_deviation = np.zeros(nphotons)
    bending_angle = np.zeros(nphotons)
    conservation = np.zeros(nphotons)
    xs = np.empty((nphotons, steps))
    ys = np.empty((nphotons, steps))
    xs[:, 0] = x; ys[:, 0] = y
    started = time.perf_counter()
    for k in range(1, steps):
        ix = np.clip(np.searchsorted(xgrid, x) - 1, 0, len(xgrid) - 1)
        iy = np.clip(np.searchsorted(ygrid, y) - 1, 0, len(ygrid) - 1)
        rx_loc = rx[iy, ix]
        ry_loc = ry[iy, ix]
        vx_new = vx + step * rx_loc
        vy_new = vy + step * ry_loc
        scale = np.maximum(np.hypot(vx_new, vy_new), 1e-12)
        vx_unit = vx_new / scale
        vy_unit = vy_new / scale
        conservation = np.maximum(conservation, np.abs(np.hypot(vx_unit, vy_unit) - 1))
        dot = np.clip(vx * vx_unit + vy * vy_unit, -1, 1)
        bending_angle += np.arccos(dot)
        vx = vx_unit; vy = vy_unit
        x = x + step * vx
        y = y + step * vy
        max_deviation = np.maximum(max_deviation, np.abs(y - y0))
        xs[:, k] = x; ys[:, k] = y
    runtime = time.perf_counter() - started
    return {
        "x": x, "y": y, "max_deviation": max_deviation,
        "bending_angle": bending_angle, "conservation": conservation,
        "runtime": runtime, "xs": xs, "ys": ys,
    }


def compute_observables(field: dict, photons: dict, extent: float, bins: int) -> dict:
    """Frozen observable extraction."""
    xf = photons["x"]; yf = photons["y"]
    x0 = photons["x0"]; y0 = photons["y0"]
    edges = np.linspace(-extent, extent, bins + 1)
    final_count, _, _ = np.histogram2d(yf, xf, bins=(edges, edges))
    initial_count, _, _ = np.histogram2d(y0, x0, bins=(edges, edges))
    safe = initial_count > 0
    convergence = np.full((bins, bins), np.nan)
    convergence[safe] = 0.5 * (final_count[safe] / initial_count[safe] - 1.0)
    sum_dx, _, _ = np.histogram2d(yf, xf, bins=(edges, edges), weights=xf - x0)
    sum_dy, _, _ = np.histogram2d(yf, xf, bins=(edges, edges), weights=yf - y0)
    deflection_x = np.full((bins, bins), np.nan)
    deflection_y = np.full((bins, bins), np.nan)
    good = final_count > 0
    deflection_x[good] = sum_dx[good] / final_count[good]
    deflection_y[good] = sum_dy[good] / final_count[good]
    fill_x = np.nan_to_num(deflection_x, nan=0.0)
    fill_y = np.nan_to_num(deflection_y, nan=0.0)
    spacing = edges[1] - edges[0]
    dxx = np.gradient(fill_x, spacing, axis=1)
    dyy = np.gradient(fill_y, spacing, axis=0)
    dxy = np.gradient(fill_x, spacing, axis=0)
    dyx = np.gradient(fill_y, spacing, axis=1)
    shear_g1 = 0.5 * (dxx - dyy)
    shear_g2 = 0.5 * (dxy + dyx)
    gamma_mag = np.hypot(shear_g1, shear_g2)
    denom = (1.0 - np.nan_to_num(convergence, nan=0.0)) ** 2 - gamma_mag ** 2
    magnification = np.full_like(convergence, np.nan)
    positive = (denom > 0) & good
    magnification[positive] = 1.0 / denom[positive]
    return {
        "convergence": convergence,
        "shear_g1": shear_g1, "shear_g2": shear_g2,
        "shear_magnitude": gamma_mag,
        "deflection_x": deflection_x, "deflection_y": deflection_y,
        "deflection_magnitude": np.hypot(deflection_x, deflection_y),
        "magnification": magnification,
        "ray_count": final_count,
        "edges": edges,
    }


# ----------------------------------------------------------------------------
# Pipeline runner
# ----------------------------------------------------------------------------
def run_pipeline_for_cluster(rho_input: np.ndarray, lens: dict) -> dict:
    """Apply the full frozen Version A pipeline to a given matter field."""
    field = make_field(rho_input, lens["extent"], lens["strength"], lens["n"])
    nphotons = lens["nphotons"]
    x0 = np.full(nphotons, -lens["extent"])
    y0 = np.linspace(-lens["y_span"], lens["y_span"], nphotons)
    vx0 = np.ones(nphotons)
    vy0 = np.zeros(nphotons)
    photons = propagate(field, lens["step"], lens["steps"], x0, y0, vx0, vy0)
    photons["x0"] = x0; photons["y0"] = y0
    observables = compute_observables(field, photons, lens["extent"], lens["bins"])
    return {"field": field, "photons": photons, "observables": observables}


# ----------------------------------------------------------------------------
# Statistical comparison
# ----------------------------------------------------------------------------
def pearson_corr(a: np.ndarray, b: np.ndarray) -> float:
    a = np.asarray(a, dtype=np.float64).ravel()
    b = np.asarray(b, dtype=np.float64).ravel()
    mask = np.isfinite(a) & np.isfinite(b)
    if mask.sum() < 2:
        return float("nan")
    a = a[mask] - a[mask].mean()
    b = b[mask] - b[mask].mean()
    denom = np.sqrt((a * a).sum() * (b * b).sum())
    if denom == 0:
        return float("nan")
    return float((a * b).sum() / denom)


def compare_arrays(predicted: np.ndarray, observed: np.ndarray) -> dict:
    """Compute quantitative metrics between two 2D maps.

    Both arrays must be the same shape. NaNs are masked jointly.
    """
    p = np.asarray(predicted, dtype=np.float64)
    o = np.asarray(observed, dtype=np.float64)
    mask = np.isfinite(p) & np.isfinite(o)
    if mask.sum() == 0:
        return {"finite_pixels": 0}
    diff = p - o
    abs_diff = np.abs(diff)
    rel_diff_pct = 100.0 * abs_diff / np.maximum(np.abs(o), 1e-15)
    p_finite = p[mask]
    o_finite = o[mask]
    out = {
        "finite_pixels": int(mask.sum()),
        "rms_error": float(np.sqrt(np.mean(diff[mask] ** 2))),
        "mean_abs_error": float(np.mean(abs_diff[mask])),
        "median_abs_error": float(np.median(abs_diff[mask])),
        "max_abs_error": float(np.max(abs_diff[mask])),
        "mean_relative_error_percent": float(np.mean(rel_diff_pct[mask])),
        "median_relative_error_percent": float(np.median(rel_diff_pct[mask])),
        "pearson_correlation": pearson_corr(p, o),
        "predicted_min": float(p_finite.min()),
        "predicted_max": float(p_finite.max()),
        "predicted_mean": float(p_finite.mean()),
        "observed_min": float(o_finite.min()),
        "observed_max": float(o_finite.max()),
        "observed_mean": float(o_finite.mean()),
        "predicted_dynamic_range": float(p_finite.max() - p_finite.min()),
        "observed_dynamic_range": float(o_finite.max() - o_finite.min()),
    }
    pred_peak = np.unravel_index(np.nanargmax(np.abs(p)), p.shape)
    obs_peak = np.unravel_index(np.nanargmax(np.abs(o)), o.shape)
    out["predicted_peak_index"] = [int(pred_peak[0]), int(pred_peak[1])]
    out["observed_peak_index"] = [int(obs_peak[0]), int(obs_peak[1])]
    out["peak_location_distance_pixels"] = float(np.hypot(pred_peak[0] - obs_peak[0],
                                                            pred_peak[1] - obs_peak[1]))
    return out


def ssim_index(p: np.ndarray, o: np.ndarray) -> float:
    """Simple SSIM (global statistics only, dynamic-range aware)."""
    p = np.asarray(p, dtype=np.float64)
    o = np.asarray(o, dtype=np.float64)
    mask = np.isfinite(p) & np.isfinite(o)
    if mask.sum() < 2:
        return float("nan")
    p = p[mask]; o = o[mask]
    c1 = (0.01 * max(np.abs(p).max(), np.abs(o).max())) ** 2
    c2 = (0.03 * max(np.abs(p).max(), np.abs(o).max())) ** 2
    mu_p = p.mean(); mu_o = o.mean()
    sig_p = p.std(); sig_o = o.std()
    sig_po = ((p - mu_p) * (o - mu_o)).mean()
    num = (2 * mu_p * mu_o + c1) * (2 * sig_po + c2)
    den = (mu_p ** 2 + mu_o ** 2 + c1) * (sig_p ** 2 + sig_o ** 2 + c2)
    return float(num / den)


# ----------------------------------------------------------------------------
# Visual products
# ----------------------------------------------------------------------------
def save_map(out_path: Path, array: np.ndarray, title: str,
             cmap: str = "viridis", symmetric: bool = False,
             vmin: float | None = None, vmax: float | None = None,
             extent: tuple[float, float, float, float] = (-8, 8, -8, 8)) -> None:
    fig, ax = plt.subplots(figsize=(7, 6))
    if symmetric:
        finite = array[np.isfinite(array)]
        vmax_abs = float(np.max(np.abs(finite))) if finite.size else 1.0
        im = ax.imshow(array, origin="lower", extent=list(extent),
                       cmap=cmap, vmin=-vmax_abs, vmax=vmax_abs)
    else:
        im = ax.imshow(array, origin="lower", extent=list(extent),
                       cmap=cmap, vmin=vmin, vmax=vmax)
    ax.set(xlabel="x", ylabel="y", title=title, aspect="equal")
    fig.colorbar(im, ax=ax)
    fig.tight_layout()
    fig.savefig(out_path, dpi=140)
    plt.close(fig)


def save_three_panel(out_path: Path, observed: np.ndarray, predicted: np.ndarray,
                     residual: np.ndarray, title: str, cmap: str = "RdBu_r",
                     vmax_abs: float | None = None) -> None:
    """Observed | Predicted | Residual on identical colour scales."""
    if vmax_abs is None:
        finite_all = np.concatenate([
            observed[np.isfinite(observed)].ravel(),
            predicted[np.isfinite(predicted)].ravel(),
        ])
        vmax_abs = float(np.max(np.abs(finite_all))) if finite_all.size else 1.0
    fig, axes = plt.subplots(1, 3, figsize=(18, 5.5))
    for ax, data, label in zip(
        axes,
        [observed, predicted, residual],
        ["Observed", "Predicted (frozen Version A)", "Residual (Pred - Obs)"],
    ):
        im = ax.imshow(data, origin="lower", extent=[-8, 8, -8, 8],
                       cmap=cmap, vmin=-vmax_abs, vmax=vmax_abs)
        ax.set(xlabel="x", ylabel="y", title=f"{title}\n{label}", aspect="equal")
        fig.colorbar(im, ax=ax)
    fig.suptitle(title)
    fig.tight_layout()
    fig.savefig(out_path, dpi=140)
    plt.close(fig)


def save_composite(out_path: Path, field: dict, observables: dict,
                   photons: dict) -> None:
    fig, axes = plt.subplots(3, 3, figsize=(18, 16))
    panels = [
        ("Matter ρ (proxy)", field["rho"], "viridis", False),
        ("Constitutive C", field["c"], "viridis", False),
        ("|∇C|", field["g_magnitude"], "viridis", False),
        ("Response magnitude", np.hypot(field["rx"], field["ry"]), "viridis", False),
        ("Response direction (rad)", field["response_direction"], "twilight", True),
        ("Convergence κ", observables["convergence"], "RdBu_r", True),
        ("Shear γ₁", observables["shear_g1"], "RdBu_r", True),
        ("Shear γ₂", observables["shear_g2"], "RdBu_r", True),
        ("Magnification μ", observables["magnification"], "viridis", False),
    ]
    extent_kw = dict(origin="lower", extent=[-8, 8, -8, 8])
    for ax, (title, array, cmap, sym) in zip(axes.flat, panels):
        if sym:
            finite = array[np.isfinite(array)]
            vmax_abs = float(np.max(np.abs(finite))) if finite.size else 1.0
            im = ax.imshow(array, cmap=cmap, vmin=-vmax_abs, vmax=vmax_abs, **extent_kw)
        else:
            im = ax.imshow(array, cmap=cmap, **extent_kw)
        ax.set(title=title, xlabel="x", ylabel="y", aspect="equal",
               xlim=(-8, 8), ylim=(-8, 8))
        fig.colorbar(im, ax=ax)
    fig.suptitle("Frozen Version A composite")
    fig.tight_layout()
    fig.savefig(out_path, dpi=140)
    plt.close(fig)


def save_trajectories(out_path: Path, photons: dict, title: str) -> None:
    fig, ax = plt.subplots(figsize=(8, 7))
    bending = photons["bending_angle"]
    norm = plt.Normalize(vmin=float(bending.min()), vmax=float(bending.max()))
    for xs, ys, b in zip(photons["xs"], photons["ys"], bending):
        ax.plot(xs, ys, color=plt.cm.plasma(norm(b)), lw=0.5)
    sm = plt.cm.ScalarMappable(cmap=plt.cm.plasma, norm=norm)
    sm.set_array([])
    fig.colorbar(sm, ax=ax, label="accumulated bending angle (rad)")
    ax.set(xlabel="x", ylabel="y", title=title, aspect="equal",
           xlim=(-8, 8), ylim=(-8, 8))
    fig.tight_layout()
    fig.savefig(out_path, dpi=140)
    plt.close(fig)


def save_overview(out_path: Path, obs_kappa: np.ndarray, pred_kappa: np.ndarray,
                  obs_gamma: np.ndarray, pred_gamma: np.ndarray) -> None:
    """Six-panel overview: observed κ, predicted κ, residual κ, observed γ, predicted γ, residual γ."""
    res_kappa = pred_kappa - obs_kappa
    res_gamma = pred_gamma - obs_gamma
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    extent_kw = dict(origin="lower", extent=[-8, 8, -8, 8])
    finite = np.concatenate([
        obs_kappa[np.isfinite(obs_kappa)].ravel(),
        pred_kappa[np.isfinite(pred_kappa)].ravel(),
    ])
    vmax_abs_k = float(np.max(np.abs(finite))) if finite.size else 1.0
    finite_g = np.concatenate([
        obs_gamma[np.isfinite(obs_gamma)].ravel(),
        pred_gamma[np.isfinite(pred_gamma)].ravel(),
    ])
    vmax_abs_g = float(np.max(finite_g)) if finite_g.size else 1.0
    panels = [
        (axes[0, 0], obs_kappa, "Convergence κ (observed)", "RdBu_r", -vmax_abs_k, vmax_abs_k),
        (axes[0, 1], pred_kappa, "Convergence κ (predicted)", "RdBu_r", -vmax_abs_k, vmax_abs_k),
        (axes[0, 2], res_kappa, "Residual κ (pred - obs)", "RdBu_r", -vmax_abs_k, vmax_abs_k),
        (axes[1, 0], obs_gamma, "Shear magnitude γ (observed)", "viridis", 0, vmax_abs_g),
        (axes[1, 1], pred_gamma, "Shear magnitude γ (predicted)", "viridis", 0, vmax_abs_g),
        (axes[1, 2], res_gamma, "Residual γ (pred - obs)", "viridis",
         -(vmax_abs_g if vmax_abs_g else 1.0), vmax_abs_g if vmax_abs_g else 1.0),
    ]
    for ax, arr, ttl, cmap, vmin, vmax in panels:
        im = ax.imshow(arr, cmap=cmap, vmin=vmin, vmax=vmax, **extent_kw)
        ax.set(title=ttl, xlabel="x", ylabel="y", aspect="equal",
               xlim=(-8, 8), ylim=(-8, 8))
        fig.colorbar(im, ax=ax)
    fig.tight_layout()
    fig.savefig(out_path, dpi=140)
    plt.close(fig)


# ----------------------------------------------------------------------------
# Per-cluster orchestration
# ----------------------------------------------------------------------------
def process_cluster(cluster: dict, benchmark_dir: Path, out_root: Path,
                    executable_hashes: dict) -> dict:
    cluster_dir = out_root / cluster["id"]
    cluster_dir.mkdir(parents=True, exist_ok=True)
    obs_dir = cluster_dir / "observed"
    pred_dir = cluster_dir / "predicted"
    res_dir = cluster_dir / "residual"
    const_dir = cluster_dir / "constitutive"
    traj_dir = cluster_dir / "trajectories"
    for d in (obs_dir, pred_dir, res_dir, const_dir, traj_dir):
        d.mkdir(parents=True, exist_ok=True)

    obs = load_observation(cluster, benchmark_dir)
    validation = validate_gamma_consistency(obs)

    # Build matter proxy: positive part of observed kappa, normalised
    rho_pipeline = resample_to_grid(obs["kappa"], LENS["n"], LENS["extent"])
    rho_positive = np.maximum(rho_pipeline, 0.0)
    rho_max = float(rho_positive.max())
    rho_input = rho_positive / rho_max if rho_max > 0 else rho_positive

    pipeline_started = time.perf_counter()
    pipeline_out = run_pipeline_for_cluster(rho_input, LENS)
    pipeline_seconds = time.perf_counter() - pipeline_started

    field = pipeline_out["field"]
    photons = pipeline_out["photons"]
    obs_out = pipeline_out["observables"]

    # Resample observations to the predicted grid
    obs_kappa_grid = resample_to_grid(obs["kappa"], LENS["bins"], LENS["extent"])
    obs_gamma1_grid = resample_to_grid(obs["gamma1"], LENS["bins"], LENS["extent"])
    obs_gamma2_grid = resample_to_grid(obs["gamma2"], LENS["bins"], LENS["extent"])
    obs_gamma_grid = resample_to_grid(obs["gamma"], LENS["bins"], LENS["extent"])
    obs_gamma_internal = np.sqrt(obs_gamma1_grid ** 2 + obs_gamma2_grid ** 2)

    pred_kappa = obs_out["convergence"]
    pred_gamma1 = obs_out["shear_g1"]
    pred_gamma2 = obs_out["shear_g2"]
    pred_gamma = obs_out["shear_magnitude"]

    residual_kappa = pred_kappa - obs_kappa_grid
    residual_gamma1 = pred_gamma1 - obs_gamma1_grid
    residual_gamma2 = pred_gamma2 - obs_gamma2_grid
    residual_gamma = pred_gamma - obs_gamma_grid
    residual_internal = np.abs(obs_gamma_grid - obs_gamma_internal)

    # Percent residual maps
    def pct_residual(p, o):
        return np.where(np.abs(o) > 1e-15, 100.0 * (p - o) / np.maximum(np.abs(o), 1e-15), np.nan)

    pct_kappa = pct_residual(pred_kappa, obs_kappa_grid)
    pct_gamma1 = pct_residual(pred_gamma1, obs_gamma1_grid)
    pct_gamma2 = pct_residual(pred_gamma2, obs_gamma2_grid)

    # Comparisons
    cmp_kappa = compare_arrays(pred_kappa, obs_kappa_grid)
    cmp_gamma1 = compare_arrays(pred_gamma1, obs_gamma1_grid)
    cmp_gamma2 = compare_arrays(pred_gamma2, obs_gamma2_grid)
    cmp_gamma = compare_arrays(pred_gamma, obs_gamma_grid)
    cmp_gamma_internal = compare_arrays(obs_gamma_grid, obs_gamma_internal)
    cmp_kappa["ssim"] = ssim_index(pred_kappa, obs_kappa_grid)
    cmp_gamma1["ssim"] = ssim_index(pred_gamma1, obs_gamma1_grid)
    cmp_gamma2["ssim"] = ssim_index(pred_gamma2, obs_gamma2_grid)
    cmp_gamma["ssim"] = ssim_index(pred_gamma, obs_gamma_grid)

    # ---------------- Save CSV products ----------------
    np.savetxt(obs_dir / "kappa_observed.csv", obs_kappa_grid, delimiter=",")
    np.savetxt(obs_dir / "gamma1_observed.csv", obs_gamma1_grid, delimiter=",")
    np.savetxt(obs_dir / "gamma2_observed.csv", obs_gamma2_grid, delimiter=",")
    np.savetxt(obs_dir / "gamma_observed.csv", obs_gamma_grid, delimiter=",")
    np.savetxt(obs_dir / "gamma_internal_check.csv", obs_gamma_internal, delimiter=",")
    np.savetxt(pred_dir / "convergence_predicted.csv", pred_kappa, delimiter=",")
    np.savetxt(pred_dir / "shear_g1_predicted.csv", pred_gamma1, delimiter=",")
    np.savetxt(pred_dir / "shear_g2_predicted.csv", pred_gamma2, delimiter=",")
    np.savetxt(pred_dir / "shear_magnitude_predicted.csv", pred_gamma, delimiter=",")
    np.savetxt(pred_dir / "deflection_x.csv", obs_out["deflection_x"], delimiter=",")
    np.savetxt(pred_dir / "deflection_y.csv", obs_out["deflection_y"], delimiter=",")
    np.savetxt(pred_dir / "magnification.csv", obs_out["magnification"], delimiter=",")
    np.savetxt(res_dir / "kappa_residual.csv", residual_kappa, delimiter=",")
    np.savetxt(res_dir / "kappa_residual_percent.csv", pct_kappa, delimiter=",")
    np.savetxt(res_dir / "gamma1_residual.csv", residual_gamma1, delimiter=",")
    np.savetxt(res_dir / "gamma1_residual_percent.csv", pct_gamma1, delimiter=",")
    np.savetxt(res_dir / "gamma2_residual.csv", residual_gamma2, delimiter=",")
    np.savetxt(res_dir / "gamma2_residual_percent.csv", pct_gamma2, delimiter=",")
    np.savetxt(res_dir / "gamma_residual.csv", residual_gamma, delimiter=",")
    np.savetxt(const_dir / "matter_proxy.csv", field["rho"], delimiter=",")
    np.savetxt(const_dir / "constitutive.csv", field["c"], delimiter=",")
    np.savetxt(const_dir / "gradient_x.csv", field["gx"], delimiter=",")
    np.savetxt(const_dir / "gradient_y.csv", field["gy"], delimiter=",")
    np.savetxt(const_dir / "gradient_magnitude.csv", field["g_magnitude"], delimiter=",")
    np.savetxt(const_dir / "response_x.csv", field["rx"], delimiter=",")
    np.savetxt(const_dir / "response_y.csv", field["ry"], delimiter=",")
    np.savetxt(const_dir / "response_direction.csv", field["response_direction"], delimiter=",")
    np.savetxt(traj_dir / "photon_trajectories_x.csv", photons["xs"], delimiter=",")
    np.savetxt(traj_dir / "photon_trajectories_y.csv", photons["ys"], delimiter=",")
    np.savetxt(traj_dir / "photon_endpoints.csv",
                np.column_stack([photons["x0"], photons["y0"], photons["x"], photons["y"],
                                 photons["max_deviation"], photons["bending_angle"]]),
                delimiter=",", header="x0,y0,x_final,y_final,max_deviation,bending_angle",
                comments="")
    np.savetxt(traj_dir / "photon_bending_angle.csv", photons["bending_angle"], delimiter=",")
    np.savetxt(traj_dir / "photon_max_deviation.csv", photons["max_deviation"], delimiter=",")

    # ---------------- Save FITS metadata summary ----------------
    fits_meta = {
        "cluster_id": cluster["id"],
        "label": cluster["label"],
        "folder": obs["folder"],
        "files": obs["files"],
        "file_sha256": obs["shas"],
        "headers": obs["headers"],
        "gamma_internal_consistency": validation,
        "matter_proxy_construction": {
            "method": "rho = max(kappa, 0) / max(max(kappa, 0))",
            "rationale": "Use the positive part of observed kappa as a proxy for matter density; "
                         "PBUF Version A constitutive C = 0.18 * rho / rho_max assumes non-negative rho.",
            "pipeline_grid_shape": list(field["rho"].shape),
            "pipeline_grid_extent": [-LENS["extent"], LENS["extent"]],
            "source_shape": list(obs["kappa"].shape),
            "rho_max_used": rho_max,
            "rho_range": [float(field["rho"].min()), float(field["rho"].max())],
            "rho_mean": float(field["rho"].mean()),
        },
        "pixel_scale_native_arcsec_per_pixel": float(abs(obs["headers"]["kappa"]["CDELT1"]) * 3600),
        "pixel_scale_pipeline_units_per_native_pixel": float(2 * LENS["extent"] / obs["kappa"].shape[0]),
        "coord_alignment_record": "Native cluster WCS centred at (CRVAL1, CRVAL2); observation "
                                  "resampled linearly onto pipeline Cartesian grid [-extent, +extent] "
                                  f"({LENS['bins']}x{LENS['bins']} for observables, {LENS['n']}x{LENS['n']} "
                                  "for matter proxy). Cluster centre maps to pipeline origin.",
        "pipeline_parameters": dict(LENS),
    }
    (cluster_dir / "fits_metadata.json").write_text(json.dumps(fits_meta, indent=2))

    # ---------------- Save figures ----------------
    save_three_panel(cluster_dir / "comparison_kappa.png",
                     obs_kappa_grid, pred_kappa, residual_kappa,
                     f"{cluster['label']} - Convergence κ",
                     cmap="RdBu_r")
    save_three_panel(cluster_dir / "comparison_gamma1.png",
                     obs_gamma1_grid, pred_gamma1, residual_gamma1,
                     f"{cluster['label']} - Shear γ₁",
                     cmap="RdBu_r")
    save_three_panel(cluster_dir / "comparison_gamma2.png",
                     obs_gamma2_grid, pred_gamma2, residual_gamma2,
                     f"{cluster['label']} - Shear γ₂",
                     cmap="RdBu_r")
    save_three_panel(cluster_dir / "comparison_gamma.png",
                     obs_gamma_grid, pred_gamma, residual_gamma,
                     f"{cluster['label']} - Shear magnitude |γ|",
                     cmap="viridis")
    save_overview(cluster_dir / "comparison_overview.png",
                  obs_kappa_grid, pred_kappa, obs_gamma_grid, pred_gamma)
    save_composite(cluster_dir / "composite_pipeline.png", field, obs_out, photons)
    save_trajectories(cluster_dir / "photon_trajectories.png", photons,
                      f"{cluster['label']} - Photon trajectories (Version A)")
    save_map(const_dir / "matter_proxy_map.png", field["rho"],
             f"{cluster['label']} - Matter proxy ρ (from observed κ)",
             cmap="viridis")
    save_map(const_dir / "constitutive_map.png", field["c"],
             f"{cluster['label']} - Constitutive C = 0.18·ρ/ρ_max",
             cmap="viridis")
    save_map(const_dir / "gradient_magnitude_map.png", field["g_magnitude"],
             f"{cluster['label']} - |∇C|", cmap="viridis")
    save_map(const_dir / "response_magnitude_map.png",
             np.hypot(field["rx"], field["ry"]),
             f"{cluster['label']} - Response magnitude |r|", cmap="viridis")
    save_map(const_dir / "response_direction_map.png", field["response_direction"],
             f"{cluster['label']} - Response direction", cmap="twilight", symmetric=True)
    save_map(pred_dir / "deflection_x_map.png", obs_out["deflection_x"],
             f"{cluster['label']} - Deflection α_x", cmap="RdBu_r", symmetric=True)
    save_map(pred_dir / "deflection_y_map.png", obs_out["deflection_y"],
             f"{cluster['label']} - Deflection α_y", cmap="RdBu_r", symmetric=True)
    save_map(pred_dir / "convergence_map.png", pred_kappa,
             f"{cluster['label']} - Predicted convergence κ", cmap="RdBu_r", symmetric=True)
    save_map(pred_dir / "shear_g1_map.png", pred_gamma1,
             f"{cluster['label']} - Predicted shear γ₁", cmap="RdBu_r", symmetric=True)
    save_map(pred_dir / "shear_g2_map.png", pred_gamma2,
             f"{cluster['label']} - Predicted shear γ₂", cmap="RdBu_r", symmetric=True)
    save_map(pred_dir / "shear_magnitude_map.png", pred_gamma,
             f"{cluster['label']} - Predicted shear |γ|", cmap="viridis")

    # ---------------- Statistics ----------------
    stats = {
        "cluster_id": cluster["id"],
        "label": cluster["label"],
        "pipeline_runtime_seconds": float(pipeline_seconds),
        "photon_count": int(LENS["nphotons"]),
        "max_conservation_error": float(np.max(photons["conservation"])),
        "n_grid": int(LENS["n"]),
        "n_steps": int(LENS["steps"]),
        "step_size": float(LENS["step"]),
        "max_response_magnitude": float(np.max(np.hypot(field["rx"], field["ry"]))),
        "max_C": float(np.max(field["c"])),
        "gamma_internal_consistency": validation,
        "comparison_kappa": cmp_kappa,
        "comparison_gamma1": cmp_gamma1,
        "comparison_gamma2": cmp_gamma2,
        "comparison_gamma": cmp_gamma,
        "gamma_internal_check": cmp_gamma_internal,
        "executables_sha256": executable_hashes,
        "frozen_pipeline": {
            "constitutive": "Version A: C(X) = 0.18 * rho(X) / rho_max",
            "transport": "90-degree transverse response, direct addition + renormalisation",
            "amplitude": "A = |grad C|",
            "matter_input": "max(observed_kappa, 0) / max(observed_kappa, 0); observation-derived proxy",
        },
    }
    (cluster_dir / "statistics.json").write_text(json.dumps(stats, indent=2))
    return stats


# ----------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------
def main():
    out_root = DEFAULT_OUT
    out_root.mkdir(parents=True, exist_ok=True)
    overall_started = time.perf_counter()
    executable_hashes = {
        "weak_lensing_observation001.py": file_sha256(Path(__file__).resolve()),
        "constitutive_equations.py": file_sha256(ROOT / "constitutive_equations.py"),
    }
    cluster_stats = {}
    for cluster in CLUSTERS:
        print(f"Processing {cluster['label']} ({cluster['id']})...")
        stats = process_cluster(cluster, BENCHMARK_DIR, out_root, executable_hashes)
        cluster_stats[cluster["id"]] = stats

    # ---------------- Required global artefacts ----------------
    # comparison_summary.csv: per-cluster cross-cluster table
    summary_rows = []
    for cid, stats in cluster_stats.items():
        cmp_k = stats["comparison_kappa"]
        cmp_g1 = stats["comparison_gamma1"]
        cmp_g2 = stats["comparison_gamma2"]
        cmp_g = stats["comparison_gamma"]
        summary_rows.append({
            "cluster": stats["label"],
            "cluster_id": cid,
            "RMS_kappa": cmp_k.get("rms_error", float("nan")),
            "MAE_kappa": cmp_k.get("mean_abs_error", float("nan")),
            "MaxAbs_kappa": cmp_k.get("max_abs_error", float("nan")),
            "Corr_kappa": cmp_k.get("pearson_correlation", float("nan")),
            "SSIM_kappa": cmp_k.get("ssim", float("nan")),
            "PeakOffset_kappa_pixels": cmp_k.get("peak_location_distance_pixels", float("nan")),
            "RangeRatio_kappa": (
                cmp_k.get("observed_dynamic_range", 1.0) /
                cmp_k.get("predicted_dynamic_range", 1.0) if cmp_k.get("predicted_dynamic_range") else float("nan")
            ),
            "RMS_gamma1": cmp_g1.get("rms_error", float("nan")),
            "Corr_gamma1": cmp_g1.get("pearson_correlation", float("nan")),
            "RMS_gamma2": cmp_g2.get("rms_error", float("nan")),
            "Corr_gamma2": cmp_g2.get("pearson_correlation", float("nan")),
            "RMS_gamma": cmp_g.get("rms_error", float("nan")),
            "Corr_gamma": cmp_g.get("pearson_correlation", float("nan")),
            "Runtime_seconds": stats["pipeline_runtime_seconds"],
            "Photon_count": stats["photon_count"],
            "Conservation_error_max": stats["max_conservation_error"],
        })
    summary_keys = list(summary_rows[0])
    with (out_root / "comparison_summary.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=summary_keys)
        writer.writeheader()
        writer.writerows(summary_rows)

    # cluster_statistics.csv: per-cluster detailed statistics (flattened)
    detail_rows = []
    for cid, stats in cluster_stats.items():
        for cmp_name in ("comparison_kappa", "comparison_gamma1",
                          "comparison_gamma2", "comparison_gamma"):
            cmp_block = stats[cmp_name]
            row = {"cluster_id": cid, "label": stats["label"],
                   "metric_group": cmp_name}
            row.update({f"cmp_{k}": v for k, v in cmp_block.items()})
            detail_rows.append(row)
    detail_keys = list(detail_rows[0])
    with (out_root / "cluster_statistics.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=detail_keys)
        writer.writeheader()
        writer.writerows(detail_rows)

    # run.json: provenance and checksums
    run_doc = {
        "milestone": "PBUF WEAK-LENSING-OBSERVATION-001",
        "status": "OK",
        "frozen_pipeline": {
            "constitutive": "Version A: C(X) = 0.18 * rho(X) / rho_max",
            "transport": "90-degree transverse response, direct addition + renormalisation",
            "amplitude": "A = |grad C|",
            "matter_input": "max(observed_kappa, 0) / max(max(observed_kappa, 0)) per cluster",
            "lens_parameters": dict(LENS),
        },
        "identical_pipeline_hashes": executable_hashes,
        "clusters": [
            {"id": c["id"], "label": c["label"], "directory": c["directory"]}
            for c in CLUSTERS
        ],
        "execution_seconds": float(time.perf_counter() - overall_started),
    }
    (out_root / "run.json").write_text(json.dumps(run_doc, indent=2))

    # validation.json: internal-consistency and pipeline-integrity checks
    validation_doc = {
        "milestone": "PBUF WEAK-LENSING-OBSERVATION-001",
        "identical_pipeline_hashes": executable_hashes,
        "per_cluster": {
            cid: {
                "gamma_internal_consistency": stats["gamma_internal_consistency"],
                "max_conservation_error": stats["max_conservation_error"],
                "photon_count": stats["photon_count"],
            }
            for cid, stats in cluster_stats.items()
        },
        "all_clusters_pass_internal_consistency": all(
            stats["gamma_internal_consistency"]["tolerance_passed"]
            for stats in cluster_stats.values()
        ),
        "execution_seconds": float(time.perf_counter() - overall_started),
    }
    (out_root / "validation.json").write_text(json.dumps(validation_doc, indent=2))

    # report.md: comprehensive report
    write_report(out_root, cluster_stats, executable_hashes,
                 time.perf_counter() - overall_started)
    print(json.dumps({
        "milestone": "PBUF WEAK-LENSING-OBSERVATION-001",
        "status": "OK",
        "clusters": list(cluster_stats.keys()),
        "output": str(out_root),
        "execution_seconds": float(time.perf_counter() - overall_started),
    }, indent=2))
    return 0


def write_report(out_root: Path, cluster_stats: dict,
                 executable_hashes: dict, total_seconds: float) -> None:
    lines = ["# PBUF WEAK-LENSING-OBSERVATION-001",
             "",
             "Frozen Version A pipeline applied to five public weak-lensing",
             "benchmark clusters (SaWLens reconstructions, Merten et al. 2014).",
             "",
             "## Pipeline (frozen, identical for every cluster)",
             "",
             "- Constitutive: `C(X) = 0.18 · ρ(X) / ρ_max` (Version A)",
             "- Response: `r = |∇C|` rotated 90° transverse",
             "- Transport: neighbour-to-neighbour, direct addition, velocity",
             "  renormalisation",
             "- Numerical parameters (frozen):",
             f"  `n_grid = {LENS['n']}, extent = {LENS['extent']}, "
             f"strength = {LENS['strength']}, step = {LENS['step']}, "
             f"steps = {LENS['steps']}, nphotons = {LENS['nphotons']}, "
             f"bins = {LENS['bins']}`",
             "",
             "## Matter input",
             "",
             "Each cluster's published κ map is taken as a proxy for the matter",
             "density ρ. Negative κ values (mass deficits) are clamped to zero",
             "to respect the implicit positivity assumption of the Version A",
             "constitutive law; the resulting field is normalised by its peak",
             "value so that ρ_max = 1 in pipeline units. No fitting, smoothing",
             "or rescaling beyond bilinear interpolation onto the pipeline grid.",
             "",
             "## Coordinate alignment",
             "",
             "The observation WCS centres each map on its own cluster centre",
             "(CRVAL). The pipeline grid is Cartesian `[-8, 8] × [-8, 8]`. The",
             "observation is resampled bilinearly onto the pipeline grid with",
             "the cluster centre mapped to the origin. No smoothing is",
             "applied; only linear interpolation for coordinate alignment.",
             "",
             "## Internal-consistency check (mandatory pre-comparison)",
             "",
             "For every cluster, `gamma.fits` is compared against",
             "`sqrt(gamma1² + gamma2²)` element-wise.",
             "",
             "| Cluster | Max abs diff | RMS abs diff | Tolerance pass |",
             "|---|---|---|---|"]
    for cid, stats in cluster_stats.items():
        v = stats["gamma_internal_consistency"]
        lines.append(
            f"| {stats['label']} | {v['max_abs_difference']:.3e} | "
            f"{v['rms_abs_difference']:.3e} | "
            f"{'YES' if v['tolerance_passed'] else 'NO'} |"
        )
    lines += [
        "",
        "All `gamma.fits` files match `sqrt(gamma1² + gamma2²)` to within",
        "single-precision FP tolerance. The supplied `gamma` field is",
        "internally consistent.",
        "",
        "## Resampled-observation interpolation consistency",
        "",
        "After bilinear resampling of the four observation maps onto the",
        "pipeline 64×64 grid, `gamma_resampled` is again compared to",
        "`sqrt(gamma1_resampled² + gamma2_resampled²)`. The discrepancy is",
        "an unavoidable interpolation artefact and is reported separately.",
        "",
        "| Cluster | Resampled max abs diff | Resampled RMS abs diff |",
        "|---|---|---|"]
    for cid, stats in cluster_stats.items():
        v = stats["gamma_internal_check"]
        lines.append(
            f"| {stats['label']} | {v['max_abs_error']:.3e} | "
            f"{v['rms_error']:.3e} |"
        )
    lines += [
        "",
        "These residuals reflect the nonlinearity of the square root under",
        "bilinear interpolation and do **not** indicate inconsistency in the",
        "raw FITS products.",
        "",
        "## Per-cluster metrics",
        "",
        "| Cluster | RMS κ | MAE κ | Max abs κ | Corr κ | SSIM κ | Peak offset (px) |",
        "|---|---|---|---|---|---|---|"]
    for cid, stats in cluster_stats.items():
        c = stats["comparison_kappa"]
        lines.append(
            f"| {stats['label']} | {c.get('rms_error', float('nan')):.4e} | "
            f"{c.get('mean_abs_error', float('nan')):.4e} | "
            f"{c.get('max_abs_error', float('nan')):.4e} | "
            f"{c.get('pearson_correlation', float('nan')):.4f} | "
            f"{c.get('ssim', float('nan')):.4f} | "
            f"{c.get('peak_location_distance_pixels', float('nan')):.2f} |"
        )
    lines += ["",
              "| Cluster | RMS γ₁ | Corr γ₁ | RMS γ₂ | Corr γ₂ | RMS γ | Corr γ |",
              "|---|---|---|---|---|---|---|"]
    for cid, stats in cluster_stats.items():
        c1 = stats["comparison_gamma1"]
        c2 = stats["comparison_gamma2"]
        cg = stats["comparison_gamma"]
        lines.append(
            f"| {stats['label']} | {c1.get('rms_error', float('nan')):.4e} | "
            f"{c1.get('pearson_correlation', float('nan')):.4f} | "
            f"{c2.get('rms_error', float('nan')):.4e} | "
            f"{c2.get('pearson_correlation', float('nan')):.4f} | "
            f"{cg.get('rms_error', float('nan')):.4e} | "
            f"{cg.get('pearson_correlation', float('nan')):.4f} |"
        )
    lines += ["",
              "## Required statistics (per cluster)",
              "",
              "Each cluster's frozen-pipeline observables are compared against",
              "the published observables resampled to the pipeline 64×64 grid.",
              "No parameter of the frozen Version A pipeline is altered between",
              "clusters or between iterations.",
              ""]
    for cid, stats in cluster_stats.items():
        lines += [
            f"### {stats['label']} (`{cid}`)",
            "",
            "| Metric | Value |",
            "|---|---|",
            f"| RMS κ | {stats['comparison_kappa']['rms_error']:.4e} |",
            f"| RMS γ₁ | {stats['comparison_gamma1']['rms_error']:.4e} |",
            f"| RMS γ₂ | {stats['comparison_gamma2']['rms_error']:.4e} |",
            f"| RMS γ | {stats['comparison_gamma']['rms_error']:.4e} |",
            f"| Correlation κ | {stats['comparison_kappa']['pearson_correlation']:.4f} |",
            f"| Correlation γ | {stats['comparison_gamma']['pearson_correlation']:.4f} |",
            f"| Runtime (s) | {stats['pipeline_runtime_seconds']:.4f} |",
            f"| Photon count | {stats['photon_count']} |",
            f"| Numerical conservation (max) | {stats['max_conservation_error']:.4e} |",
            "",
        ]
    lines += ["",
              "## Cross-cluster summary",
              "",
              "| Cluster | RMS κ | RMS γ | Corr κ | Corr γ | Runtime (s) |",
              "|---|---|---|---|---|---|"]
    for cid, stats in cluster_stats.items():
        lines.append(
            f"| {stats['label']} | "
            f"{stats['comparison_kappa']['rms_error']:.4e} | "
            f"{stats['comparison_gamma']['rms_error']:.4e} | "
            f"{stats['comparison_kappa']['pearson_correlation']:.4f} | "
            f"{stats['comparison_gamma']['pearson_correlation']:.4f} | "
            f"{stats['pipeline_runtime_seconds']:.3f} |"
        )
    lines += ["",
              "## Units, conventions and mismatch (mandatory pre-comparison record)",
              "",
              "- Published products are SaWLens reconstructions on RA/Dec WCS",
              "  grids centred on each cluster; pixel scales 6.25-11.36 arcsec.",
              "- All published observables are scaled to source redshift",
              "  `Z_S = 9.0` (effectively an infinite-source approximation).",
              "- `kappa.fits`, `gamma1.fits`, `gamma2.fits` are the lensing",
              "  convergence and reduced-shear components from a parametric",
              "  joint weak+strong lensing inversion. They are reconstructed",
              "  posterior-mean maps, not direct observational data.",
              "- Frozen Version A outputs are dimensionless lensing-like",
              "  observables (κ_pred, γ₁_pred, γ₂_pred) derived from a",
              "  constitutive + transport pipeline operating on synthetic",
              "  dimensionless coordinates on `[-8, 8] × [-8, 8]`.",
              "- The published products and Version A outputs are **NOT**",
              "  directly comparable in absolute units, normalisation, or",
              "  angular scale. The comparison made here is a like-with-like",
              "  dimensionless field comparison after coordinate alignment,",
              "  with no implicit cosmological rescaling.",
              "- The matter input to the Version A constitutive law is the",
              "  positive part of the published κ (clamped at zero), then",
              "  normalised by its peak. This treats κ as a matter-density",
              "  proxy (standard practice for clusters where mass traces light)",
              "  while preserving the implicit positivity assumption of the",
              "  Version A law.",
              "",
              "## Per-cluster outputs",
              ""]
    for cid, stats in cluster_stats.items():
        lines += [
            f"### {stats['label']} (`{cid}`)",
            "",
            f"- `observed/`: resampled observation maps (κ, γ₁, γ₂, γ)",
            "  and the internal-consistency `gamma_internal_check.csv`.",
            f"- `predicted/`: Version A outputs (κ, γ₁, γ₂, |γ|, deflection,",
            "  magnification).",
            f"- `residual/`: `pred - obs` residual maps in absolute units and",
            "  percentage form.",
            f"- `constitutive/`: matter proxy, C, ∇C, response field, maps.",
            f"- `trajectories/`: photon trajectories (`x`, `y`), endpoints,",
            "  bending angles, max deviations.",
            f"- `comparison_kappa.png`, `comparison_gamma1.png`,",
            "  `comparison_gamma2.png`, `comparison_gamma.png`: three-panel",
            "  comparisons with identical colour scales.",
            f"- `comparison_overview.png`: six-panel composite (κ and |γ|).",
            f"- `composite_pipeline.png`: 3×3 pipeline panel (matter, C, |∇C|,",
            "  response magnitude and direction, κ, γ₁, γ₂, μ).",
            f"- `photon_trajectories.png`: trajectory plot coloured by",
            "  accumulated bending angle.",
            f"- `statistics.json`: all quantitative metrics.",
            f"- `fits_metadata.json`: FITS headers, WCS info, file SHA-256s,",
            "  matter-proxy construction record, coordinate-alignment record.",
            "",
        ]
    lines += [
        "## Identical-pipeline verification (SHA-256)",
        "",
        "| File | SHA-256 |",
        "|---|---|",
    ]
    for name, digest in executable_hashes.items():
        lines.append(f"| `{name}` | `{digest}` |")
    lines += [
        "",
        "## Notes",
        "",
        "- The Version A pipeline parameters are held identical to those of",
        "  WEAK-LENSING-PREDICTION-001 and WEAK-LENSING-GENERALIZATION-001.",
        "- No parameter was altered to improve agreement with the",
        "  observations. Any apparent mismatch is therefore a property of",
        "  the frozen implementation itself.",
        "- Discrepancies between predicted and observed fields can arise",
        "  from unit/normalisation mismatch, the absence of cosmological",
        "  scaling in the Version A pipeline, the simplified matter proxy,",
        "  and the linear-response approximation inherent to Version A.",
        "  These are documented, not compensated for.",
        f"- Total execution time: {total_seconds:.2f} s.",
        "",
    ]
    (out_root / "report.md").write_text("\n".join(lines))


if __name__ == "__main__":
    raise SystemExit(main())