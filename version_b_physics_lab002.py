#!/usr/bin/env python3
"""PBUF VERSION-B PHYSICS-LAB-002 - Response Family Decomposition.

Decomposes C10 (Combined Local Response) from PHYSICS-LAB-001 into the
two physical mechanisms it contains:

    Mechanism A : Neighbour Coherence  (multiplicative factor)
    Mechanism B : Elastic Memory       (one-step persistence mix)

Tests four candidate local response laws inside the FROZEN Version 1
weak-lensing laboratory (LAB-FREEZE-001 / WEAK-LENSING-SCIENCE-001).

Only the response law (rx, ry) is varied.  No modification to
- Constitutive Version A
- Frozen transport (propagation)
- Source-plane geometry (Launch B)
- Observable extraction (Jacobian)
- Numerical configuration (20 000 photons, 256^2 grid, Delta s / 2)

Phase 1 - Decompose C10

| Candidate | Mechanism A | Mechanism B |
|-----------|-------------|-------------|
| Control   | no          | no          |
| C10-A     | yes         | no          |
| C10-B     | no          | yes         |
| C10-C     | yes         | yes         |

Phase 2 - Interaction matrix
Pair AB = Coherence + Memory (the only pair contained in C10).

Phase 3 - Leave-one-out
With N=2 mechanisms: N=2 variants (one per removed component).

For every variant compute
- Median Pearson kappa
- Median Pearson gamma
- RMS kappa
- RMS gamma
- SSIM (kappa, gamma)
- kappa bias
- gamma bias
- Conservation
- Runtime

Contribution analysis reports
- Individual contribution
- Pairwise contribution
- Combined contribution
- Interaction effect

Reports whether the C10 improvement is additive, synergistic, redundant.

Required outputs (runs/version_b_physics_lab002/):
- report.md
- component_contributions.csv
- interaction_matrix.csv
- leave_one_out.csv
- family_ranking.csv
- run.json
- validation.json
- plots/component_importance.png
- plots/interaction_heatmap.png
- plots/leave_one_out.png
- plots/synergy_matrix.png
- plots/family_contributions.png
"""
from __future__ import annotations

import csv
import json
import sys
import time
import traceback
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits

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
DEFAULT_OUT = ROOT / "runs" / "version_b_physics_lab002"
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
# Decomposition of C10
# =============================================================================
# C10 = coherence_factor * ((1-w) * R_90(g) + w * R_90(g_prev))
#
# Mechanism A : Neighbour Coherence
#   coherence_factor = 0.5 * (1 + mean_cos(theta_self, theta_neighbours))
#   In [0, 1]; suppresses response in regions of incoherent gradients.
#
# Mechanism B : Elastic Memory
#   r_mem = (1-w) * r_self + w * r_prev
#   w = 0.5; r_prev is the linear response one cell upstream in -x.
#
# Pure gradient response (Control) is just R_90(g).
#
# The four candidates below apply every subset of {A, B} to the linear
# baseline, with no modification to transport, constitutive, source, or
# observable extraction.
# =============================================================================


def _coherence_factor(gx: np.ndarray, gy: np.ndarray) -> np.ndarray:
    """Mechanism A : Neighbour coherence factor in [0, 1]."""
    g_safe = np.maximum(np.hypot(gx, gy), 1e-15)
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
    return 0.5 * (1.0 + mean_cos)


def _linear_response(gx: np.ndarray, gy: np.ndarray,
                     g: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Baseline linear response R_90(g).  Equal to (rx_self, ry_self)
    in C10."""
    g_safe = np.maximum(g, 1e-15)
    rx = -g * (gy / g_safe)
    ry = +g * (gx / g_safe)
    return rx, ry


def _memory_mix(rx_self: np.ndarray, ry_self: np.ndarray,
                w: float = 0.5) -> tuple[np.ndarray, np.ndarray]:
    """Mechanism B : one-step persistence, w = 0.5 (fixed)."""
    rx_prev = np.roll(rx_self, 1, axis=1)
    ry_prev = np.roll(ry_self, 1, axis=1)
    rx_prev[:, 0] = rx_self[:, 0]
    ry_prev[:, 0] = ry_self[:, 0]
    rx = (1.0 - w) * rx_self + w * rx_prev
    ry = (1.0 - w) * ry_self + w * ry_prev
    return rx, ry


# Candidate 1: Control (pure gradient).  No coherence, no memory.
def response_control(c, xgrid, ygrid, gx, gy, g):
    rx, ry = _linear_response(gx, gy, g)
    return rx, ry


# Candidate C10-A: Neighbour Coherence only.
# Magnitude scaled by (1 + mean_cos)/2 over the 8 immediate neighbours.
def response_c10_a(c, xgrid, ygrid, gx, gy, g):
    rx, ry = _linear_response(gx, gy, g)
    factor = _coherence_factor(gx, gy)
    rx *= factor
    ry *= factor
    return rx, ry


# Candidate C10-B: Elastic Memory only.
# r_new = (1-w) * R(g) + w * R(g_prev), w = 0.5.
def response_c10_b(c, xgrid, ygrid, gx, gy, g):
    rx, ry = _linear_response(gx, gy, g)
    rx, ry = _memory_mix(rx, ry, w=0.5)
    return rx, ry


# Candidate C10-C: Combined Local Response (original C10).
# factor * ((1-w) * R(g) + w * R(g_prev)).
def response_c10_c(c, xgrid, ygrid, gx, gy, g):
    rx, ry = _linear_response(gx, gy, g)
    rx, ry = _memory_mix(rx, ry, w=0.5)
    factor = _coherence_factor(gx, gy)
    rx *= factor
    ry *= factor
    return rx, ry


@dataclass(frozen=True)
class CandidateSpec:
    code: str           # 'CONTROL', 'C10-A', 'C10-B', 'C10-C'
    name: str
    family: str
    description: str
    coherence: bool     # mechanism A
    memory: bool        # mechanism B
    law: Callable
    notes: str = ""


CANDIDATES = [
    CandidateSpec(
        "CONTROL", "Gradient (control)", "gradient",
        "Response = |grad C|; frozen Version A control.",
        False, False, response_control,
        notes="Pure gradient baseline. No C10 mechanism enabled.",
    ),
    CandidateSpec(
        "C10-A", "Neighbour Coherence only", "neighbour coherence",
        "Magnitude scaled by (1 + mean_cos)/2 over 8 neighbours.",
        True, False, response_c10_a,
        notes="Mechanism A only (no memory).",
    ),
    CandidateSpec(
        "C10-B", "Elastic Memory only", "elastic memory",
        "r_new = (1-w)*R(g) + w*R(g_prev); w = 0.5.",
        False, True, response_c10_b,
        notes="Mechanism B only (no coherence).",
    ),
    CandidateSpec(
        "C10-C", "Combined (original C10)", "combined response",
        "Coherence factor times elastic memory mix (original C10).",
        True, True, response_c10_c,
        notes="Both mechanisms enabled; identical to C10 from LAB-001.",
    ),
]


# =============================================================================
# Pipeline helpers (identical to LAB-001, no frozen-component modifications)
# =============================================================================
def matter_proxy_from_kappa(kappa_native: np.ndarray, grid_n: int,
                             extent: float) -> np.ndarray:
    """Frozen matter input rule: rho = max(kappa, 0) / max(max(kappa, 0))."""
    rho_pipeline = resample_to_grid(kappa_native, grid_n, extent)
    rho_pos = np.maximum(rho_pipeline, 0.0)
    rho_max = float(rho_pos.max())
    if rho_max > 0:
        rho_pos = rho_pos / rho_max
    return rho_pos


def compute_field(rho: np.ndarray, extent: float, strength: float,
                   grid_n: int) -> dict:
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


def apply_candidate(field: dict, candidate: CandidateSpec
                     ) -> tuple[np.ndarray, np.ndarray]:
    rx, ry = candidate.law(
        field["c"], field["xgrid"], field["ygrid"],
        field["gx"], field["gy"], field["g_magnitude"]
    )
    return rx, ry


def load_observation_full(cluster: dict) -> dict:
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


def run_candidate_on_cluster(cluster: dict, cfg: dict,
                              candidate: CandidateSpec) -> dict:
    folder = BENCHMARK_DIR / cluster["directory"]
    kappa_path = folder / f"hlsp_frontier_model_{cluster['slug']}_merten_v1_kappa.fits"
    with fits.open(kappa_path) as h:
        kappa_native = np.asarray(h[0].data, dtype=np.float64)
    rho = matter_proxy_from_kappa(kappa_native, cfg["grid_n"], cfg["extent"])

    field = compute_field(rho, cfg["extent"], cfg["strength"], cfg["grid_n"])

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

    out = {}
    for key in ("kappa", "gamma1", "gamma2"):
        cmp = compare_arrays(pred[key], obs_grid[key])
        cmp["ssim"] = ssim_index(pred[key], obs_grid[key])
        out[key] = cmp
    cmp_g = compare_arrays(pred["gamma_mag"], obs_grid["gamma"])
    cmp_g["ssim"] = ssim_index(pred["gamma_mag"], obs_grid["gamma"])
    out["gamma_mag"] = cmp_g

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
        "candidate_code": candidate.code,
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


def _safe_median(values):
    arr = np.array(list(values), dtype=np.float64)
    arr = arr[np.isfinite(arr)]
    return float(np.median(arr)) if arr.size else float("nan")


def _safe_mean(values):
    arr = np.array(list(values), dtype=np.float64)
    arr = arr[np.isfinite(arr)]
    return float(np.mean(arr)) if arr.size else float("nan")


def _safe_std(values):
    arr = np.array(list(values), dtype=np.float64)
    arr = arr[np.isfinite(arr)]
    return float(np.std(arr)) if arr.size > 1 else float("nan")


# =============================================================================
# Output writers
# =============================================================================
def per_cluster_row(result: dict) -> dict:
    return {
        "candidate_code": result["candidate_code"],
        "candidate_name": result["candidate_name"],
        "candidate_family": result["candidate_family"],
        "cluster_id": result["cluster_id"],
        "cluster_label": result["cluster_label"],
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


def write_component_contributions_csv(out_root: Path,
                                       per_cluster_rows: list[dict]) -> None:
    """One row per (mechanism, cluster)."""
    path = out_root / "component_contributions.csv"
    fields = [
        "candidate_code", "candidate_name", "candidate_family",
        "cluster_id", "cluster_label",
        "rms_kappa", "rms_gamma",
        "pearson_kappa", "pearson_gamma",
        "ssim_kappa", "ssim_gamma",
        "kappa_bias", "gamma_bias",
        "std_resid_kappa", "std_resid_gamma",
        "max_conservation_error",
        "runtime_seconds",
        "n_photons", "grid_n", "step", "steps",
    ]
    with path.open("w", newline="") as h:
        writer = csv.DictWriter(h, fieldnames=fields)
        writer.writeheader()
        for row in per_cluster_rows:
            writer.writerow({k: row.get(k, "") for k in fields})


def write_interaction_matrix_csv(out_root: Path,
                                   cross_rows: list[dict],
                                   control_row: dict) -> None:
    """All pairwise combinations of mechanisms in C10.

    For N=2 there is one pair (AB) plus the four single-component
    reference rows used to compute interaction effects.
    """
    path = out_root / "interaction_matrix.csv"
    fields = [
        "combination", "n_mechanisms",
        "mechanism_a", "mechanism_b",
        "candidate_code", "candidate_name",
        "median_pearson_kappa",
        "median_pearson_gamma",
        "median_ssim_kappa",
        "median_ssim_gamma",
        "median_rms_kappa",
        "median_rms_gamma",
        "mean_kappa_bias",
        "mean_gamma_bias",
        "max_conservation_error",
        "median_runtime_seconds",
        "delta_pearson_kappa_vs_control",
        "delta_ssim_kappa_vs_control",
        "delta_kappa_bias_vs_control",
        "interaction_pearson_kappa",
        "interaction_ssim_kappa",
        "interaction_kappa_bias",
    ]
    by_code = {r["candidate_code"]: r for r in cross_rows}
    with path.open("w", newline="") as h:
        writer = csv.DictWriter(h, fieldnames=fields)
        writer.writeheader()
        # Baseline (control).
        ctrl = control_row
        for code in ("CONTROL", "C10-A", "C10-B", "C10-C"):
            r = by_code[code]
            row = {
                "combination": code,
                "n_mechanisms": int(
                    (1 if code in ("C10-A", "C10-C") else 0)
                    + (1 if code in ("C10-B", "C10-C") else 0)
                ),
                "mechanism_a": "yes" if code in ("C10-A", "C10-C") else "no",
                "mechanism_b": "yes" if code in ("C10-B", "C10-C") else "no",
                "candidate_code": code,
                "candidate_name": r["candidate_name"],
                "median_pearson_kappa": r["median_pearson_kappa"],
                "median_pearson_gamma": r["median_pearson_gamma"],
                "median_ssim_kappa": r["median_ssim_kappa"],
                "median_ssim_gamma": r["median_ssim_gamma"],
                "median_rms_kappa": r["median_rms_kappa"],
                "median_rms_gamma": r["median_rms_gamma"],
                "mean_kappa_bias": r["mean_kappa_bias"],
                "mean_gamma_bias": r["mean_gamma_bias"],
                "max_conservation_error": r["max_conservation_error"],
                "median_runtime_seconds": r["median_runtime_seconds"],
                "delta_pearson_kappa_vs_control":
                    r["median_pearson_kappa"] - ctrl["median_pearson_kappa"],
                "delta_ssim_kappa_vs_control":
                    r["median_ssim_kappa"] - ctrl["median_ssim_kappa"],
                "delta_kappa_bias_vs_control":
                    r["mean_kappa_bias"] - ctrl["mean_kappa_bias"],
                "interaction_pearson_kappa": "",
                "interaction_ssim_kappa": "",
                "interaction_kappa_bias": "",
            }
            writer.writerow(row)
        # Pairwise combinations (only AB exists in C10).
        a = by_code["C10-A"]
        b = by_code["C10-B"]
        ab = by_code["C10-C"]
        # interaction = combined - control - [(A-control) + (B-control)]
        delta_a_pk = a["median_pearson_kappa"] - ctrl["median_pearson_kappa"]
        delta_b_pk = b["median_pearson_kappa"] - ctrl["median_pearson_kappa"]
        delta_ab_pk = ab["median_pearson_kappa"] - ctrl["median_pearson_kappa"]
        interaction_pk = delta_ab_pk - (delta_a_pk + delta_b_pk)
        delta_a_sk = a["median_ssim_kappa"] - ctrl["median_ssim_kappa"]
        delta_b_sk = b["median_ssim_kappa"] - ctrl["median_ssim_kappa"]
        delta_ab_sk = ab["median_ssim_kappa"] - ctrl["median_ssim_kappa"]
        interaction_sk = delta_ab_sk - (delta_a_sk + delta_b_sk)
        delta_a_bk = a["mean_kappa_bias"] - ctrl["mean_kappa_bias"]
        delta_b_bk = b["mean_kappa_bias"] - ctrl["mean_kappa_bias"]
        delta_ab_bk = ab["mean_kappa_bias"] - ctrl["mean_kappa_bias"]
        interaction_bk = delta_ab_bk - (delta_a_bk + delta_b_bk)
        ab_row = {
            "combination": "AB",
            "n_mechanisms": 2,
            "mechanism_a": "yes",
            "mechanism_b": "yes",
            "candidate_code": "C10-C",
            "candidate_name": ab["candidate_name"],
            "median_pearson_kappa": ab["median_pearson_kappa"],
            "median_pearson_gamma": ab["median_pearson_gamma"],
            "median_ssim_kappa": ab["median_ssim_kappa"],
            "median_ssim_gamma": ab["median_ssim_gamma"],
            "median_rms_kappa": ab["median_rms_kappa"],
            "median_rms_gamma": ab["median_rms_gamma"],
            "mean_kappa_bias": ab["mean_kappa_bias"],
            "mean_gamma_bias": ab["mean_gamma_bias"],
            "max_conservation_error": ab["max_conservation_error"],
            "median_runtime_seconds": ab["median_runtime_seconds"],
            "delta_pearson_kappa_vs_control": delta_ab_pk,
            "delta_ssim_kappa_vs_control": delta_ab_sk,
            "delta_kappa_bias_vs_control": delta_ab_bk,
            "interaction_pearson_kappa": interaction_pk,
            "interaction_ssim_kappa": interaction_sk,
            "interaction_kappa_bias": interaction_bk,
        }
        writer.writerow(ab_row)


def write_leave_one_out_csv(out_root: Path,
                              cross_rows: list[dict],
                              control_row: dict) -> None:
    """Leave-one-out: each variant removes exactly one mechanism."""
    path = out_root / "leave_one_out.csv"
    fields = [
        "removed_mechanism", "remaining_combination", "candidate_code",
        "median_pearson_kappa", "median_pearson_gamma",
        "median_ssim_kappa", "median_ssim_gamma",
        "mean_kappa_bias", "mean_gamma_bias",
        "delta_pearson_kappa_vs_control",
        "delta_ssim_kappa_vs_control",
        "delta_kappa_bias_vs_control",
        "delta_pearson_kappa_vs_full_c10",
        "delta_ssim_kappa_vs_full_c10",
        "delta_kappa_bias_vs_full_c10",
        "lost_pearson_kappa",
        "lost_ssim_kappa",
        "lost_kappa_bias",
    ]
    by_code = {r["candidate_code"]: r for r in cross_rows}
    ctrl = control_row
    full = by_code["C10-C"]
    full_pk = full["median_pearson_kappa"]
    full_sk = full["median_ssim_kappa"]
    full_bk = full["mean_kappa_bias"]
    with path.open("w", newline="") as h:
        writer = csv.DictWriter(h, fieldnames=fields)
        writer.writeheader()
        # Remove Coherence (A) -> remaining = Memory only = C10-B
        rem = by_code["C10-B"]
        lost_pk = full_pk - rem["median_pearson_kappa"]
        lost_sk = full_sk - rem["median_ssim_kappa"]
        lost_bk = full_bk - rem["mean_kappa_bias"]
        row = {
            "removed_mechanism": "Coherence (A)",
            "remaining_combination": "Memory only",
            "candidate_code": "C10-B",
            "median_pearson_kappa": rem["median_pearson_kappa"],
            "median_pearson_gamma": rem["median_pearson_gamma"],
            "median_ssim_kappa": rem["median_ssim_kappa"],
            "median_ssim_gamma": rem["median_ssim_gamma"],
            "mean_kappa_bias": rem["mean_kappa_bias"],
            "mean_gamma_bias": rem["mean_gamma_bias"],
            "delta_pearson_kappa_vs_control":
                rem["median_pearson_kappa"] - ctrl["median_pearson_kappa"],
            "delta_ssim_kappa_vs_control":
                rem["median_ssim_kappa"] - ctrl["median_ssim_kappa"],
            "delta_kappa_bias_vs_control":
                rem["mean_kappa_bias"] - ctrl["mean_kappa_bias"],
            "delta_pearson_kappa_vs_full_c10":
                rem["median_pearson_kappa"] - full_pk,
            "delta_ssim_kappa_vs_full_c10":
                rem["median_ssim_kappa"] - full_sk,
            "delta_kappa_bias_vs_full_c10":
                rem["mean_kappa_bias"] - full_bk,
            "lost_pearson_kappa": lost_pk,
            "lost_ssim_kappa": lost_sk,
            "lost_kappa_bias": lost_bk,
        }
        writer.writerow(row)
        # Remove Memory (B) -> remaining = Coherence only = C10-A
        rem = by_code["C10-A"]
        lost_pk = full_pk - rem["median_pearson_kappa"]
        lost_sk = full_sk - rem["median_ssim_kappa"]
        lost_bk = full_bk - rem["mean_kappa_bias"]
        row = {
            "removed_mechanism": "Elastic Memory (B)",
            "remaining_combination": "Coherence only",
            "candidate_code": "C10-A",
            "median_pearson_kappa": rem["median_pearson_kappa"],
            "median_pearson_gamma": rem["median_pearson_gamma"],
            "median_ssim_kappa": rem["median_ssim_kappa"],
            "median_ssim_gamma": rem["median_ssim_gamma"],
            "mean_kappa_bias": rem["mean_kappa_bias"],
            "mean_gamma_bias": rem["mean_gamma_bias"],
            "delta_pearson_kappa_vs_control":
                rem["median_pearson_kappa"] - ctrl["median_pearson_kappa"],
            "delta_ssim_kappa_vs_control":
                rem["median_ssim_kappa"] - ctrl["median_ssim_kappa"],
            "delta_kappa_bias_vs_control":
                rem["mean_kappa_bias"] - ctrl["mean_kappa_bias"],
            "delta_pearson_kappa_vs_full_c10":
                rem["median_pearson_kappa"] - full_pk,
            "delta_ssim_kappa_vs_full_c10":
                rem["median_ssim_kappa"] - full_sk,
            "delta_kappa_bias_vs_full_c10":
                rem["mean_kappa_bias"] - full_bk,
            "lost_pearson_kappa": lost_pk,
            "lost_ssim_kappa": lost_sk,
            "lost_kappa_bias": lost_bk,
        }
        writer.writerow(row)


def write_family_ranking_csv(out_root: Path, ranking_rows: list[dict]) -> None:
    """Family ranking by individual contribution (median Pearson kappa
    improvement over control)."""
    path = out_root / "family_ranking.csv"
    fields = [
        "rank", "candidate_code", "candidate_name", "candidate_family",
        "median_pearson_kappa", "median_pearson_gamma",
        "median_ssim_kappa", "median_ssim_gamma",
        "mean_kappa_bias", "mean_gamma_bias",
        "alone_contribution_pearson_kappa",
        "alone_contribution_kappa_bias",
        "combined_contribution_pearson_kappa",
        "lost_when_removed_pearson_kappa",
        "is_redundant",
    ]
    with path.open("w", newline="") as h:
        writer = csv.DictWriter(h, fieldnames=fields)
        writer.writeheader()
        for r in ranking_rows:
            writer.writerow({k: r.get(k, "") for k in fields})


# =============================================================================
# Plot helpers
# =============================================================================
def save_component_importance_plot(out_path: Path,
                                     cross_rows: list[dict],
                                     control_row: dict) -> None:
    """Bar chart of individual mechanism contributions (median Pearson
    kappa improvement over control)."""
    codes = ["C10-A", "C10-B", "C10-C"]
    labels = [r["candidate_name"] for r in
              [next(x for x in cross_rows if x["candidate_code"] == c)
               for c in codes]]
    pk = [r["median_pearson_kappa"]
          for r in [next(x for x in cross_rows if x["candidate_code"] == c)
                    for c in codes]]
    delta_pk = [v - control_row["median_pearson_kappa"] for v in pk]
    sk = [r["median_ssim_kappa"]
          for r in [next(x for x in cross_rows if x["candidate_code"] == c)
                    for c in codes]]
    delta_sk = [v - control_row["median_ssim_kappa"] for v in sk]
    bias = [r["mean_kappa_bias"]
            for r in [next(x for x in cross_rows if x["candidate_code"] == c)
                      for c in codes]]
    delta_bias = [v - control_row["mean_kappa_bias"] for v in bias]

    fig, axes = plt.subplots(1, 3, figsize=(15, 6))
    x_pos = np.arange(len(codes))
    colours = ["steelblue", "darkorange", "seagreen"]
    axes[0].bar(x_pos, delta_pk, color=colours, edgecolor="black")
    axes[0].set_xticks(x_pos)
    axes[0].set_xticklabels(codes, fontsize=9)
    axes[0].set_ylabel("Median Pearson kappa improvement vs control")
    axes[0].set_title("Component importance (Pearson kappa)", pad=15)
    axes[0].axhline(0.0, color="black", lw=0.7, ls=":")
    axes[0].grid(axis="y", alpha=0.3)
    ymin0, ymax0 = axes[0].get_ylim()
    yr0 = max(abs(ymax0), abs(ymin0), 1e-12)
    for i, v in enumerate(delta_pk):
        offset = 0.03 * yr0
        if v >= 0:
            label_y = v + offset
            va = "bottom"
        else:
            label_y = v - offset
            va = "top"
        axes[0].text(i, label_y, f"{v:+.4f}", ha="center",
                     va=va, fontsize=9)
    axes[1].bar(x_pos, delta_sk, color=colours, edgecolor="black")
    axes[1].set_xticks(x_pos)
    axes[1].set_xticklabels(codes, fontsize=9)
    axes[1].set_ylabel("Median SSIM kappa improvement vs control")
    axes[1].set_title("Component importance (SSIM kappa)", pad=15)
    axes[1].axhline(0.0, color="black", lw=0.7, ls=":")
    axes[1].grid(axis="y", alpha=0.3)
    ymin1, ymax1 = axes[1].get_ylim()
    yr1 = max(abs(ymax1), abs(ymin1), 1e-12)
    for i, v in enumerate(delta_sk):
        offset = 0.03 * yr1
        if v >= 0:
            label_y = v + offset
            va = "bottom"
        else:
            label_y = v - offset
            va = "top"
        axes[1].text(i, label_y, f"{v:+.4f}", ha="center",
                     va=va, fontsize=9)
    axes[2].bar(x_pos, delta_bias, color=colours, edgecolor="black")
    axes[2].set_xticks(x_pos)
    axes[2].set_xticklabels(codes, fontsize=9)
    axes[2].set_ylabel("Mean kappa bias delta vs control")
    axes[2].set_title("Component importance (kappa bias)", pad=15)
    axes[2].axhline(0.0, color="black", lw=0.7, ls=":")
    axes[2].grid(axis="y", alpha=0.3)
    ymin2, ymax2 = axes[2].get_ylim()
    yr2 = max(abs(ymax2), abs(ymin2), 1e-12)
    for i, v in enumerate(delta_bias):
        offset = 0.03 * yr2
        if v >= 0:
            label_y = v + offset
            va = "bottom"
        else:
            label_y = v - offset
            va = "top"
        axes[2].text(i, label_y, f"{v:+.5f}", ha="center",
                     va=va, fontsize=9)
    fig.subplots_adjust(top=0.85, bottom=0.15, wspace=0.3)
    fig.text(0.5, 0.95,
              "C10 component importance (decomposition)",
              ha="center", fontsize=12)
    fig.savefig(out_path, dpi=130)
    plt.close(fig)


def save_interaction_heatmap_plot(out_path: Path,
                                    cross_rows: list[dict],
                                    control_row: dict) -> None:
    """Heatmap of pairwise interaction effects, one per metric."""
    by_code = {r["candidate_code"]: r for r in cross_rows}
    ctrl = control_row
    a = by_code["C10-A"]
    b = by_code["C10-B"]
    ab = by_code["C10-C"]

    def deltas(metric: str):
        da = a[metric] - ctrl[metric]
        db = b[metric] - ctrl[metric]
        dab = ab[metric] - ctrl[metric]
        return da, db, dab

    delta_a_pk, delta_b_pk, delta_ab_pk = deltas("median_pearson_kappa")
    delta_a_sk, delta_b_sk, delta_ab_sk = deltas("median_ssim_kappa")
    delta_a_bk, delta_b_bk, delta_ab_bk = deltas("mean_kappa_bias")

    matrices = [
        (np.array([
            [0.0, delta_a_pk, delta_ab_pk],
            [delta_a_pk, 0.0, delta_b_pk],
            [delta_ab_pk, delta_b_pk, 0.0],
        ]), "Median Pearson kappa delta vs control"),
        (np.array([
            [0.0, delta_a_sk, delta_ab_sk],
            [delta_a_sk, 0.0, delta_b_sk],
            [delta_ab_sk, delta_b_sk, 0.0],
        ]), "Median SSIM kappa delta vs control"),
        (np.array([
            [0.0, delta_a_bk, delta_ab_bk],
            [delta_a_bk, 0.0, delta_b_bk],
            [delta_ab_bk, delta_b_bk, 0.0],
        ]), "Mean kappa bias delta vs control"),
    ]
    interactions = (
        delta_ab_pk - (delta_a_pk + delta_b_pk),
        delta_ab_sk - (delta_a_sk + delta_b_sk),
        delta_ab_bk - (delta_a_bk + delta_b_bk),
    )

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    labels = ["Control\n(none)", "Coherence\n(A)", "Coherence+Memory\n(AB)"]
    for ax, (M, title), interaction in zip(axes, matrices, interactions):
        vmax = float(np.nanmax(np.abs(M)))
        if not np.isfinite(vmax) or vmax == 0:
            vmax = 1e-6
        im = ax.imshow(M, cmap="RdBu_r",
                        vmin=-vmax, vmax=vmax, aspect="auto")
        ax.set_xticks(range(3))
        ax.set_yticks(range(3))
        ax.set_xticklabels(labels, fontsize=8)
        ax.set_yticklabels(labels, fontsize=8)
        for i in range(3):
            for j in range(3):
                ax.text(j, i, f"{M[i, j]:+.4f}", ha="center", va="center",
                         fontsize=9, color="black")
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        ax.set_title(f"{title}\nInteraction = {interaction:+.5f}",
                      fontsize=10)
    fig.suptitle(
        "C10 interaction matrix - pairwise improvements vs control\n"
        "(off-diagonal: marginal pair, diagonal: combined)",
        fontsize=10
    )
    fig.tight_layout()
    fig.savefig(out_path, dpi=130)
    plt.close(fig)


def save_leave_one_out_plot(out_path: Path, cross_rows: list[dict],
                              control_row: dict) -> None:
    """Leave-one-out: median Pearson kappa and lost performance."""
    by_code = {r["candidate_code"]: r for r in cross_rows}
    ctrl = control_row
    full = by_code["C10-C"]
    loo_a = by_code["C10-B"]  # coherence removed
    loo_b = by_code["C10-A"]  # memory removed

    labels = ["Control", "Coherence\nonly (C10-A)",
              "Memory\nonly (C10-B)", "Combined\n(C10-C)"]
    pk = [ctrl["median_pearson_kappa"],
          loo_b["median_pearson_kappa"],
          loo_a["median_pearson_kappa"],
          full["median_pearson_kappa"]]
    sk = [ctrl["median_ssim_kappa"],
          loo_b["median_ssim_kappa"],
          loo_a["median_ssim_kappa"],
          full["median_ssim_kappa"]]
    bias = [ctrl["mean_kappa_bias"],
            loo_b["mean_kappa_bias"],
            loo_a["mean_kappa_bias"],
            full["mean_kappa_bias"]]

    fig, axes = plt.subplots(1, 3, figsize=(15, 6))
    x_pos = np.arange(len(labels))
    axes[0].bar(x_pos, pk, color=["gray", "steelblue", "darkorange", "seagreen"],
                 edgecolor="black")
    axes[0].set_xticks(x_pos)
    axes[0].set_xticklabels(labels, fontsize=8)
    axes[0].set_ylabel("Median Pearson kappa")
    axes[0].set_title("Median Pearson kappa")
    axes[0].grid(axis="y", alpha=0.3)
    for i, v in enumerate(pk):
        axes[0].text(i, v + 0.0005, f"{v:+.4f}", ha="center", fontsize=8)
    axes[1].bar(x_pos, sk, color=["gray", "steelblue", "darkorange", "seagreen"],
                 edgecolor="black")
    axes[1].set_xticks(x_pos)
    axes[1].set_xticklabels(labels, fontsize=8)
    axes[1].set_ylabel("Median SSIM kappa")
    axes[1].set_title("Median SSIM kappa")
    axes[1].grid(axis="y", alpha=0.3)
    for i, v in enumerate(sk):
        axes[1].text(i, v + 0.0005, f"{v:+.4f}", ha="center", fontsize=8)
    axes[2].bar(x_pos, bias, color=["gray", "steelblue", "darkorange", "seagreen"],
                 edgecolor="black")
    axes[2].set_xticks(x_pos)
    axes[2].set_xticklabels(labels, fontsize=8)
    axes[2].set_ylabel("Mean kappa bias")
    axes[2].set_title("Mean kappa bias (closer to 0 better)")
    axes[2].grid(axis="y", alpha=0.3)
    for i, v in enumerate(bias):
        axes[2].text(i, v + 0.0005, f"{v:+.5f}", ha="center", fontsize=8)

    fig.suptitle(
        "Leave-one-out analysis\n"
        f"Lost Pearson kappa when removing Coherence = "
        f"{full['median_pearson_kappa'] - loo_a['median_pearson_kappa']:+.5f}; "
        f"when removing Memory = "
        f"{full['median_pearson_kappa'] - loo_b['median_pearson_kappa']:+.5f}"
    )
    fig.tight_layout()
    fig.savefig(out_path, dpi=130)
    plt.close(fig)


def save_synergy_matrix_plot(out_path: Path, cross_rows: list[dict],
                               control_row: dict) -> None:
    """Synergy matrix: shows how the combined improvement differs from
    the sum of individual improvements, and shows the magnitude of
    synergy per mechanism."""
    by_code = {r["candidate_code"]: r for r in cross_rows}
    ctrl = control_row
    a = by_code["C10-A"]
    b = by_code["C10-B"]
    ab = by_code["C10-C"]

    d_pk_a = a["median_pearson_kappa"] - ctrl["median_pearson_kappa"]
    d_pk_b = b["median_pearson_kappa"] - ctrl["median_pearson_kappa"]
    d_pk_ab = ab["median_pearson_kappa"] - ctrl["median_pearson_kappa"]
    d_sk_a = a["median_ssim_kappa"] - ctrl["median_ssim_kappa"]
    d_sk_b = b["median_ssim_kappa"] - ctrl["median_ssim_kappa"]
    d_sk_ab = ab["median_ssim_kappa"] - ctrl["median_ssim_kappa"]
    d_bk_a = a["mean_kappa_bias"] - ctrl["mean_kappa_bias"]
    d_bk_b = b["mean_kappa_bias"] - ctrl["mean_kappa_bias"]
    d_bk_ab = ab["mean_kappa_bias"] - ctrl["mean_kappa_bias"]

    M = np.array([
        [d_pk_a, d_pk_b, d_pk_ab],
        [d_sk_a, d_sk_b, d_sk_ab],
        [d_bk_a, d_bk_b, d_bk_ab],
    ])

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    metric_labels = ["Pearson kappa", "SSIM kappa", "Kappa bias"]
    bar_labels = ["A alone", "B alone", "A+B combined"]
    colours = ["steelblue", "darkorange", "seagreen"]
    x_pos = np.arange(3)
    for i, (ax, vals, label) in enumerate(zip(axes, M, metric_labels)):
        ax.bar(x_pos, vals, color=colours, edgecolor="black")
        ax.set_xticks(x_pos)
        ax.set_xticklabels(bar_labels, fontsize=8)
        ax.set_ylabel(label)
        ax.axhline(0.0, color="black", lw=0.7, ls=":")
        ax.grid(axis="y", alpha=0.3)
        ymin, ymax = ax.get_ylim()
        yrange = max(abs(ymax), abs(ymin), 1e-12)
        for j, v in enumerate(vals):
            offset = 0.025 * yrange
            if v >= 0:
                label_y = v + offset
                va = "bottom"
            else:
                label_y = v - offset
                va = "top"
            ax.text(j, label_y, f"{v:+.5f}", ha="center",
                     va=va, fontsize=8)
        interaction = vals[2] - (vals[0] + vals[1])
        ax.set_title(f"{label}  (interaction={interaction:+.5f})",
                      fontsize=10, pad=15)

    fig.subplots_adjust(top=0.86, bottom=0.15, wspace=0.3)
    fig.text(0.5, 0.95,
              "C10 synergy matrix (improvement over control per metric)",
              ha="center", fontsize=12)
    fig.savefig(out_path, dpi=130)
    plt.close(fig)


def save_family_contributions_plot(out_path: Path,
                                     per_cluster_rows: list[dict]) -> None:
    """Family-level contributions to C10 (per-cluster Pearson kappa
    deltas vs control)."""
    cluster_ids = sorted({r["cluster_id"] for r in per_cluster_rows},
                          key=lambda c: [cl["id"] for cl in CLUSTERS].index(c))
    codes = ["C10-A", "C10-B", "C10-C"]
    colours = {"C10-A": "steelblue", "C10-B": "darkorange", "C10-C": "seagreen"}
    titles = {
        "C10-A": "Coherence only (C10-A)",
        "C10-B": "Memory only (C10-B)",
        "C10-C": "Combined (C10-C)",
    }

    fig, axes = plt.subplots(1, 3, figsize=(15, 6))
    for ax, code in zip(axes, codes):
        per_cluster_delta = []
        for cid in cluster_ids:
            ctrl_row = next((r for r in per_cluster_rows
                              if r["candidate_code"] == "CONTROL"
                              and r["cluster_id"] == cid), None)
            cand_row = next((r for r in per_cluster_rows
                              if r["candidate_code"] == code
                              and r["cluster_id"] == cid), None)
            if ctrl_row is not None and cand_row is not None:
                per_cluster_delta.append(
                    cand_row["pearson_kappa"] - ctrl_row["pearson_kappa"]
                )
            else:
                per_cluster_delta.append(0.0)
        x_pos = np.arange(len(cluster_ids))
        ax.bar(x_pos, per_cluster_delta, color=colours[code], edgecolor="black")
        ax.set_xticks(x_pos)
        ax.set_xticklabels(cluster_ids, rotation=20, ha="right", fontsize=8)
        ax.set_ylabel("Pearson kappa delta vs control")
        ax.set_title(titles[code], fontsize=10, pad=15)
        ax.axhline(0.0, color="black", lw=0.7, ls=":")
        ax.grid(axis="y", alpha=0.3)
        ymin, ymax = ax.get_ylim()
        yr = max(abs(ymax), abs(ymin), 1e-12)
        for j, v in enumerate(per_cluster_delta):
            offset = 0.05 * yr
            if v >= 0:
                label_y = v + offset
                va = "bottom"
            else:
                label_y = v - offset
                va = "top"
            ax.text(j, label_y, f"{v:+.4f}", ha="center",
                     va=va, fontsize=8)
    fig.subplots_adjust(top=0.86, bottom=0.20, wspace=0.3)
    fig.text(0.5, 0.95,
              "Family contributions to C10 (per-cluster Pearson kappa delta)",
              ha="center", fontsize=12)
    fig.savefig(out_path, dpi=130)
    plt.close(fig)


# =============================================================================
# Report helpers
# =============================================================================
def _fmt(v, fmt=".4f"):
    if isinstance(v, float):
        if not np.isfinite(v):
            return "nan"
        return format(v, fmt)
    return str(v)


def answer_q1(cross_rows: list[dict], control_row: dict) -> list[str]:
    """Q1: Which individual mechanism contributes the largest
    improvement?

    Individual mechanisms are the C10-A (Coherence) and C10-B (Memory)
    variants.  The combined C10-C is reported for context but excluded
    from the individual ranking.
    """
    candidates = [r for r in cross_rows
                  if r["candidate_code"] not in ("CONTROL", "C10-C")]
    contributions = []
    for r in candidates:
        d_pk = r["median_pearson_kappa"] - control_row["median_pearson_kappa"]
        d_sk = r["median_ssim_kappa"] - control_row["median_ssim_kappa"]
        d_bk = r["mean_kappa_bias"] - control_row["mean_kappa_bias"]
        contributions.append((r["candidate_code"], r["candidate_name"],
                              d_pk, d_sk, d_bk))
    contributions.sort(key=lambda t: -t[2])

    combined = next(r for r in cross_rows if r["candidate_code"] == "C10-C")
    combined_d_pk = (combined["median_pearson_kappa"]
                     - control_row["median_pearson_kappa"])
    combined_d_sk = (combined["median_ssim_kappa"]
                     - control_row["median_ssim_kappa"])
    combined_d_bk = (combined["mean_kappa_bias"]
                     - control_row["mean_kappa_bias"])

    lines = [
        "### Q1. Largest individual contribution",
        "",
        f"Control median Pearson kappa = "
        f"{control_row['median_pearson_kappa']:+.5f}, "
        f"median SSIM kappa = {control_row['median_ssim_kappa']:+.5f}, "
        f"mean kappa bias = {control_row['mean_kappa_bias']:+.5f}.",
        "",
        "Individual mechanism contributions (excluding the combined):",
        "",
        "| Mechanism | Median Pearson k delta | Median SSIM k delta | Mean k bias delta |",
        "|---|---|---|---|",
    ]
    for code, name, dpk, dsk, dbk in contributions:
        lines.append(
            f"| {code} {name} | {dpk:+.5f} | {dsk:+.5f} | {dbk:+.5f} |"
        )
    lines += [
        "",
        f"For reference, the combined C10-C delta = "
        f"{combined_d_pk:+.5f} (Pearson k), "
        f"{combined_d_sk:+.5f} (SSIM k), "
        f"{combined_d_bk:+.5f} (k bias).",
        "",
    ]
    best = contributions[0]
    lines += [
        f"**Largest individual contribution:** `{best[0]}` "
        f"({best[1]}) with median Pearson kappa delta = {best[2]:+.5f}.",
        "",
    ]
    return lines


def answer_q2(cross_rows: list[dict], control_row: dict) -> list[str]:
    """Q2: Can C10's improvement be explained entirely by one mechanism?"""
    full = next(r for r in cross_rows if r["candidate_code"] == "C10-C")
    a = next(r for r in cross_rows if r["candidate_code"] == "C10-A")
    b = next(r for r in cross_rows if r["candidate_code"] == "C10-B")
    ctrl = control_row
    full_d_pk = full["median_pearson_kappa"] - ctrl["median_pearson_kappa"]
    a_d_pk = a["median_pearson_kappa"] - ctrl["median_pearson_kappa"]
    b_d_pk = b["median_pearson_kappa"] - ctrl["median_pearson_kappa"]
    full_d_sk = full["median_ssim_kappa"] - ctrl["median_ssim_kappa"]
    a_d_sk = a["median_ssim_kappa"] - ctrl["median_ssim_kappa"]
    b_d_sk = b["median_ssim_kappa"] - ctrl["median_ssim_kappa"]
    full_d_bk = full["mean_kappa_bias"] - ctrl["mean_kappa_bias"]
    a_d_bk = a["mean_kappa_bias"] - ctrl["mean_kappa_bias"]
    b_d_bk = b["mean_kappa_bias"] - ctrl["mean_kappa_bias"]

    pk_a_share = abs(a_d_pk) / max(abs(full_d_pk), 1e-15) if full_d_pk != 0 else float("nan")
    pk_b_share = abs(b_d_pk) / max(abs(full_d_pk), 1e-15) if full_d_pk != 0 else float("nan")

    best_alone = max(
        ((abs(a_d_pk), "Coherence (A)", a_d_pk),
         (abs(b_d_pk), "Memory (B)", b_d_pk))
    )
    explain_share = best_alone[0] / max(abs(full_d_pk), 1e-15)
    explained = explain_share >= 0.8

    lines = [
        "### Q2. Single-mechanism explanation of C10",
        "",
        "| Mechanism | Alone delta Pearson k | Combined delta Pearson k | Share |",
        "|---|---|---|---|",
        f"| Coherence (A) | {a_d_pk:+.5f} | {full_d_pk:+.5f} | {pk_a_share*100:.1f}% |",
        f"| Memory (B)    | {b_d_pk:+.5f} | {full_d_pk:+.5f} | {pk_b_share*100:.1f}% |",
        "",
        f"Best single-mechanism explanation share = {explain_share*100:.1f}% "
        f"(mechanism `{best_alone[1]}`).",
        f"C10 improvement explained entirely by one mechanism: "
        f"**{'YES' if explained else 'NO'}**.",
        "",
        "Additional metrics:",
        "",
        "| Mechanism | SSIM k delta alone | SSIM k delta combined | Bias delta alone | Bias delta combined |",
        "|---|---|---|---|---|",
        f"| Coherence (A) | {a_d_sk:+.5f} | {full_d_sk:+.5f} | {a_d_bk:+.5f} | {full_d_bk:+.5f} |",
        f"| Memory (B)    | {b_d_sk:+.5f} | {full_d_sk:+.5f} | {b_d_bk:+.5f} | {full_d_bk:+.5f} |",
        "",
    ]
    return lines


def answer_q3(cross_rows: list[dict], control_row: dict) -> list[str]:
    """Q3: Do multiple mechanisms reinforce one another?"""
    full = next(r for r in cross_rows if r["candidate_code"] == "C10-C")
    a = next(r for r in cross_rows if r["candidate_code"] == "C10-A")
    b = next(r for r in cross_rows if r["candidate_code"] == "C10-B")
    ctrl = control_row
    a_d = a["median_pearson_kappa"] - ctrl["median_pearson_kappa"]
    b_d = b["median_pearson_kappa"] - ctrl["median_pearson_kappa"]
    ab_d = full["median_pearson_kappa"] - ctrl["median_pearson_kappa"]
    reinforce = ab_d > max(a_d, b_d) + 1e-6
    lines = [
        "### Q3. Mutual reinforcement",
        "",
        "Does the combined improvement exceed the best single-mechanism",
        "improvement?",
        "",
        f"- Coherence alone delta: {a_d:+.5f}",
        f"- Memory alone delta:    {b_d:+.5f}",
        f"- Combined delta:        {ab_d:+.5f}",
        f"- Best single delta:     {max(a_d, b_d):+.5f}",
        "",
        f"Reinforcement present (combined > best alone): "
        f"**{'YES' if reinforce else 'NO'}**.",
        "",
    ]
    return lines


def answer_q4(cross_rows: list[dict], control_row: dict) -> list[str]:
    """Q4: Is any mechanism redundant?"""
    full = next(r for r in cross_rows if r["candidate_code"] == "C10-C")
    a = next(r for r in cross_rows if r["candidate_code"] == "C10-A")
    b = next(r for r in cross_rows if r["candidate_code"] == "C10-B")
    ctrl = control_row
    ab_d = full["median_pearson_kappa"] - ctrl["median_pearson_kappa"]
    a_d = a["median_pearson_kappa"] - ctrl["median_pearson_kappa"]
    b_d = b["median_pearson_kappa"] - ctrl["median_pearson_kappa"]

    threshold = 0.005 * max(abs(ab_d), 1e-12)
    a_redundant = abs(ab_d - a_d) < threshold
    b_redundant = abs(ab_d - b_d) < threshold

    lines = [
        "### Q4. Redundancy check",
        "",
        f"Combined delta Pearson kappa = {ab_d:+.5f}.",
        f"Redundancy threshold = {threshold*100:.2f}% of |combined delta|.",
        "",
        "| Mechanism | Combined - alone | Redundant? |",
        "|---|---|---|",
        f"| Coherence (A) | {ab_d - a_d:+.5f} | {'YES' if a_redundant else 'NO'} |",
        f"| Memory (B)    | {ab_d - b_d:+.5f} | {'YES' if b_redundant else 'NO'} |",
        "",
    ]
    return lines


def answer_q5(per_cluster_rows: list[dict], cross_rows: list[dict],
                control_row: dict) -> list[str]:
    """Q5: Is C2 the only universally sign-consistent improvement?"""
    consistency = []
    for cand in CANDIDATES:
        if cand.code == "CONTROL":
            continue
        deltas = []
        for cl in CLUSTERS:
            cid = cl["id"]
            ctrl_row = next((r for r in per_cluster_rows
                              if r["candidate_code"] == "CONTROL"
                              and r["cluster_id"] == cid), None)
            cand_row = next((r for r in per_cluster_rows
                              if r["candidate_code"] == cand.code
                              and r["cluster_id"] == cid), None)
            if ctrl_row is not None and cand_row is not None:
                deltas.append(cand_row["pearson_kappa"]
                               - ctrl_row["pearson_kappa"])
        if not deltas:
            continue
        consistency.append({
            "code": cand.code,
            "name": cand.name,
            "n_positive": int(np.sum(np.array(deltas) > 0)),
            "n_negative": int(np.sum(np.array(deltas) < 0)),
            "median_delta": float(np.median(deltas)),
            "sign_consistent": bool(np.all(np.array(deltas) > 0)),
            "deltas": deltas,
        })
    consistency.sort(key=lambda r: -r["median_delta"])

    lines = [
        "### Q5. Universally sign-consistent improvement",
        "",
        "| Variant | Median delta Pearson k | # clusters +ve | # clusters -ve | Sign consistent |",
        "|---|---|---|---|---|",
    ]
    for r in consistency:
        lines.append(
            f"| {r['code']} {r['name']} "
            f"| {r['median_delta']:+.5f} "
            f"| {r['n_positive']} "
            f"| {r['n_negative']} "
            f"| {'YES' if r['sign_consistent'] else 'NO'} |"
        )
    consistent = [r for r in consistency if r["sign_consistent"]]
    lines += [
        "",
        f"Number of variants with sign-consistent improvement on all 5 "
        f"clusters: {len(consistent)}/{len(consistency)}.",
        "",
    ]
    # Cross-check with the original C2 from LAB-001.
    c10_a = next(r for r in consistency if r["code"] == "C10-A")
    lines += [
        "Cross-check with VERSION-B PHYSICS-LAB-001 (Candidate 2 = C10-A):",
        f"  C2 / C10-A sign-consistent: "
        f"**{'YES' if c10_a['sign_consistent'] else 'NO'}** "
        f"({c10_a['n_positive']} +ve, {c10_a['n_negative']} -ve, "
        f"median delta = {c10_a['median_delta']:+.5f}).",
        "",
    ]
    return lines


def answer_q6(cross_rows: list[dict], control_row: dict) -> list[str]:
    """Q6: Is the interaction additive, nonlinear, or negligible?"""
    full = next(r for r in cross_rows if r["candidate_code"] == "C10-C")
    a = next(r for r in cross_rows if r["candidate_code"] == "C10-A")
    b = next(r for r in cross_rows if r["candidate_code"] == "C10-B")
    ctrl = control_row
    a_d = a["median_pearson_kappa"] - ctrl["median_pearson_kappa"]
    b_d = b["median_pearson_kappa"] - ctrl["median_pearson_kappa"]
    ab_d = full["median_pearson_kappa"] - ctrl["median_pearson_kappa"]
    interaction_pk = ab_d - (a_d + b_d)

    a_d_sk = a["median_ssim_kappa"] - ctrl["median_ssim_kappa"]
    b_d_sk = b["median_ssim_kappa"] - ctrl["median_ssim_kappa"]
    ab_d_sk = full["median_ssim_kappa"] - ctrl["median_ssim_kappa"]
    interaction_sk = ab_d_sk - (a_d_sk + b_d_sk)

    a_d_bk = a["mean_kappa_bias"] - ctrl["mean_kappa_bias"]
    b_d_bk = b["mean_kappa_bias"] - ctrl["mean_kappa_bias"]
    ab_d_bk = full["mean_kappa_bias"] - ctrl["mean_kappa_bias"]
    interaction_bk = ab_d_bk - (a_d_bk + b_d_bk)

    magnitude = max(abs(interaction_pk),
                     abs(interaction_sk),
                     abs(interaction_bk),
                     1e-12)
    ref_magnitude = max(abs(ab_d), 1e-12)
    ratio = magnitude / ref_magnitude

    if ratio < 0.10:
        verdict = "NEGLIGIBLE"
    elif interaction_pk > 0 and interaction_sk > 0 and interaction_bk > 0:
        verdict = "SYNERGISTIC"
    elif interaction_pk < 0 and interaction_sk < 0 and interaction_bk < 0:
        verdict = "ANTAGONISTIC / REDUNDANT"
    else:
        verdict = "NONLINEAR / MIXED"

    lines = [
        "### Q6. Nature of the interaction",
        "",
        "| Metric | A alone | B alone | A+B combined | Interaction |",
        "|---|---|---|---|---|",
        f"| Pearson kappa delta | {a_d:+.5f} | {b_d:+.5f} | {ab_d:+.5f} | {interaction_pk:+.5f} |",
        f"| SSIM kappa delta    | {a_d_sk:+.5f} | {b_d_sk:+.5f} | {ab_d_sk:+.5f} | {interaction_sk:+.5f} |",
        f"| Kappa bias delta    | {a_d_bk:+.5f} | {b_d_bk:+.5f} | {ab_d_bk:+.5f} | {interaction_bk:+.5f} |",
        "",
        f"Interaction magnitude relative to combined improvement: "
        f"{ratio*100:.1f}%.",
        "",
        f"Verdict: **{verdict}**.",
        "",
    ]
    return lines


def determine_outcome(cross_rows: list[dict],
                       per_cluster_rows: list[dict]) -> str:
    """Determine Outcome A / B / C from the decomposition results."""
    full = next(r for r in cross_rows if r["candidate_code"] == "C10-C")
    a = next(r for r in cross_rows if r["candidate_code"] == "C10-A")
    b = next(r for r in cross_rows if r["candidate_code"] == "C10-B")
    ctrl = next(r for r in cross_rows if r["candidate_code"] == "CONTROL")
    eps = 2.220446049250313e-16

    full_d_pk = full["median_pearson_kappa"] - ctrl["median_pearson_kappa"]
    a_d_pk = a["median_pearson_kappa"] - ctrl["median_pearson_kappa"]
    b_d_pk = b["median_pearson_kappa"] - ctrl["median_pearson_kappa"]
    interaction_pk = full_d_pk - (a_d_pk + b_d_pk)
    ab_alone_share = max(abs(a_d_pk), abs(b_d_pk)) / max(abs(full_d_pk), 1e-15)
    rel_interaction = abs(interaction_pk) / max(abs(full_d_pk), 1e-15)

    if ab_alone_share >= 0.7 and rel_interaction < 0.30:
        return (
            f"**Outcome A** - one physical mechanism is responsible for "
            f"the majority of C10's improvement. Best single-mechanism "
            f"share = {ab_alone_share*100:.1f}% of |combined delta|, "
            f"relative interaction = {rel_interaction*100:.1f}%. "
            f"Future work should focus on the dominant mechanism."
        )
    if rel_interaction < 0.30:
        return (
            f"**Outcome B** - several mechanisms contribute independently. "
            f"Best single-mechanism share = {ab_alone_share*100:.1f}%, "
            f"relative interaction = {rel_interaction*100:.1f}%. "
            f"Future work should investigate the physical relationship "
            f"between the contributing mechanisms."
        )
    return (
        f"**Outcome C** - the apparent improvement cannot be attributed "
        f"to any individual mechanism. Best single-mechanism share = "
        f"{ab_alone_share*100:.1f}%, relative interaction = "
        f"{rel_interaction*100:.1f}%. The combined behaviour is emergent "
        f"and requires further investigation."
    )


# =============================================================================
# Main
# =============================================================================
def main():
    out_root = DEFAULT_OUT
    out_root.mkdir(parents=True, exist_ok=True)
    PLOTS.mkdir(parents=True, exist_ok=True)

    overall_started = time.perf_counter()

    print("=" * 72)
    print("PBUF VERSION-B PHYSICS-LAB-002")
    print("C10 response family decomposition - frozen Version 1 laboratory")
    print("=" * 72)

    print("\n[1/5] Verifying frozen source hashes against LAB-FREEZE-001 ...")
    hash_check = verify_frozen_hashes()
    for name, info in hash_check["files"].items():
        marker = "OK" if info["match"] else "MISMATCH"
        print(f"  [{marker}] {name}: {info['actual_sha256']}")
    if not hash_check["ok"]:
        raise RuntimeError("Frozen source file hashes do not match LAB-FREEZE-001.")

    print("\n[2/5] Running 4 candidates on all 5 clusters ...")
    per_cluster_rows: list[dict] = []
    for cand in CANDIDATES:
        mech_str = "A" if cand.coherence else "-"
        mech_str += "B" if cand.memory else "-"
        print(f"\n  {cand.code} ({mech_str}): {cand.name} ({cand.family})", flush=True)
        for cluster in CLUSTERS:
            t0 = time.perf_counter()
            try:
                result = run_candidate_on_cluster(cluster, MIN_PRODUCTION, cand)
            except Exception as exc:
                print(f"    ERROR on {cluster['id']}: {exc}\n{traceback.format_exc()}")
                raise
            row = per_cluster_row(result)
            per_cluster_rows.append(row)
            elapsed = time.perf_counter() - t0
            print(f"    -> {cluster['label']:12s}: Pearson k = {row['pearson_kappa']:+.4f}, "
                  f"k bias = {row['kappa_bias']:+.5f}, "
                  f"runtime = {row['runtime_seconds']:.3f} s, "
                  f"cons = {row['max_conservation_error']:.2e}  ({elapsed:.1f}s)",
                  flush=True)

    print("\n[3/5] Computing cross-cluster statistics ...")
    cross_rows = []
    for cand in CANDIDATES:
        sub = [r for r in per_cluster_rows
               if r["candidate_code"] == cand.code]
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
            "candidate_code": cand.code,
            "candidate_name": cand.name,
            "candidate_family": cand.family,
            "mechanism_a": cand.coherence,
            "mechanism_b": cand.memory,
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

    by_code = {r["candidate_code"]: r for r in cross_rows}
    control_row = by_code["CONTROL"]

    # Family ranking: rank by alone contribution (delta vs control).
    alone_contributions = []
    for r in cross_rows:
        if r["candidate_code"] == "CONTROL":
            continue
        full = by_code["C10-C"]
        alone_contributions.append({
            "candidate_code": r["candidate_code"],
            "candidate_name": r["candidate_name"],
            "candidate_family": r["candidate_family"],
            "median_pearson_kappa": r["median_pearson_kappa"],
            "median_pearson_gamma": r["median_pearson_gamma"],
            "median_ssim_kappa": r["median_ssim_kappa"],
            "median_ssim_gamma": r["median_ssim_gamma"],
            "mean_kappa_bias": r["mean_kappa_bias"],
            "mean_gamma_bias": r["mean_gamma_bias"],
            "alone_contribution_pearson_kappa":
                r["median_pearson_kappa"] - control_row["median_pearson_kappa"],
            "alone_contribution_kappa_bias":
                r["mean_kappa_bias"] - control_row["mean_kappa_bias"],
            "combined_contribution_pearson_kappa":
                full["median_pearson_kappa"] - control_row["median_pearson_kappa"],
            "lost_when_removed_pearson_kappa":
                full["median_pearson_kappa"] - r["median_pearson_kappa"],
            "is_redundant": bool(
                abs(full["median_pearson_kappa"] - r["median_pearson_kappa"])
                < 0.005 * max(abs(full["median_pearson_kappa"]
                                    - control_row["median_pearson_kappa"]), 1e-12)
            ),
        })
    alone_contributions.sort(
        key=lambda r: -r["alone_contribution_pearson_kappa"]
    )
    for i, r in enumerate(alone_contributions, start=1):
        r["rank"] = i
    # Re-order fields.
    ranking_rows = [{
        "rank": r["rank"],
        "candidate_code": r["candidate_code"],
        "candidate_name": r["candidate_name"],
        "candidate_family": r["candidate_family"],
        "median_pearson_kappa": r["median_pearson_kappa"],
        "median_pearson_gamma": r["median_pearson_gamma"],
        "median_ssim_kappa": r["median_ssim_kappa"],
        "median_ssim_gamma": r["median_ssim_gamma"],
        "mean_kappa_bias": r["mean_kappa_bias"],
        "mean_gamma_bias": r["mean_gamma_bias"],
        "alone_contribution_pearson_kappa":
            r["alone_contribution_pearson_kappa"],
        "alone_contribution_kappa_bias":
            r["alone_contribution_kappa_bias"],
        "combined_contribution_pearson_kappa":
            r["combined_contribution_pearson_kappa"],
        "lost_when_removed_pearson_kappa":
            r["lost_when_removed_pearson_kappa"],
        "is_redundant": r["is_redundant"],
    } for r in alone_contributions]

    print("\n[4/5] Answering required questions ...")
    q_answers = []
    q_answers += answer_q1(cross_rows, control_row)
    q_answers += answer_q2(cross_rows, control_row)
    q_answers += answer_q3(cross_rows, control_row)
    q_answers += answer_q4(cross_rows, control_row)
    q_answers += answer_q5(per_cluster_rows, cross_rows, control_row)
    q_answers += answer_q6(cross_rows, control_row)

    outcome = determine_outcome(cross_rows, per_cluster_rows)

    print("\n[5/5] Writing outputs ...")
    write_component_contributions_csv(out_root, per_cluster_rows)
    write_interaction_matrix_csv(out_root, cross_rows, control_row)
    write_leave_one_out_csv(out_root, cross_rows, control_row)
    write_family_ranking_csv(out_root, ranking_rows)

    # Feed per-cluster rows to the family_contributions plot via module attr.
    save_component_importance_plot(PLOTS / "component_importance.png",
                                     cross_rows, control_row)
    save_interaction_heatmap_plot(PLOTS / "interaction_heatmap.png",
                                   cross_rows, control_row)
    save_leave_one_out_plot(PLOTS / "leave_one_out.png",
                             cross_rows, control_row)
    save_synergy_matrix_plot(PLOTS / "synergy_matrix.png",
                              cross_rows, control_row)
    save_family_contributions_plot(PLOTS / "family_contributions.png",
                                    per_cluster_rows)

    # Report
    lines = [
        "# PBUF VERSION-B PHYSICS-LAB-002",
        "",
        "**Response family decomposition of C10 inside the frozen",
        "Version 1 weak-lensing laboratory (LAB-FREEZE-001).**",
        "",
        "C10 (Combined Local Response) from PHYSICS-LAB-001 is decomposed",
        "into the two physical mechanisms it contains: Neighbour Coherence",
        "(A) and Elastic Memory (B). Each mechanism and their combination",
        "are evaluated inside the frozen laboratory without modifying any",
        "other component.",
        "",
        "## Status",
        "",
        f"- Frozen hash verification: **{'PASS' if hash_check['ok'] else 'FAIL'}**",
        f"- Variants tested: **{len(CANDIDATES)}** (Control + 3 mechanism subsets of {{A, B}})",
        f"- Clusters: **{len(CLUSTERS)}**",
        "",
        "## Frozen laboratory",
        "",
        "The Version 1 laboratory is used as the measurement instrument",
        "without modification.  All frozen source files are verified by",
        "SHA-256 against LAB-FREEZE-001.",
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
        "## Decomposed mechanisms in C10",
        "",
        "C10 in PHYSICS-LAB-001 was implemented as",
        "",
        "    r = coherence_factor * ((1-w) * R_90(g) + w * R_90(g_prev))",
        "",
        "with `w = 0.5` and",
        "",
        "    coherence_factor = 0.5 * (1 + mean_cos(theta_self, theta_8nn))",
        "",
        "Therefore C10 contains exactly two physical mechanisms:",
        "",
        "| Code | Mechanism | Implementation |",
        "|---|---|---|",
        "| A | Neighbour Coherence | Multiplicative factor over 8 nearest neighbours |",
        "| B | Elastic Memory      | One-step persistence mix (w = 0.5) |",
        "",
        "## Variants",
        "",
        "Every subset of {A, B} is realised as a local response law.  No",
        "new physics is introduced.",
        "",
        "| Code | Mechanism A | Mechanism B | Description |",
        "|---|---|---|---|",
    ]
    for c in CANDIDATES:
        lines.append(
            f"| {c.code} | "
            f"{'yes' if c.coherence else 'no'} | "
            f"{'yes' if c.memory else 'no'} | "
            f"{c.description} |"
        )
    lines += [
        "",
        "All fixed parameters are documented in the candidate source",
        "code.  No parameter is fitted.",
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
        "## Per-variant, per-cluster metrics",
        "",
        "Computed metrics for every (variant, cluster) pair:",
        "",
        "| Variant | Cluster | RMS k | RMS g | Pearson k | Pearson g "
        "| SSIM k | SSIM g | k bias | g bias | conservation | runtime (s) |",
        "|---|---|---|---|---|---|---|---|---|---|---|---|",
    ]
    for r in per_cluster_rows:
        lines.append(
            f"| {r['candidate_code']} {r['candidate_name']} "
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
        "For every variant the following medians/means are taken across",
        "the five benchmark clusters.",
        "",
        "| Variant | Median Pearson k | Median Pearson g | Median SSIM k "
        "| Mean k Bias | Mean g Bias | Mean Pearson k | Conservation max "
        "| Runtime (s) |",
        "|---|---|---|---|---|---|---|---|---|",
    ]
    for r in cross_rows:
        lines.append(
            f"| {r['candidate_code']} {r['candidate_name']} "
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
        "## Family ranking by individual contribution",
        "",
        "Variants ranked by `median Pearson kappa` minus the control",
        "(individual contribution).  Combined contribution and",
        "lost-when-removed are reported alongside.",
        "",
        "| Rank | Code | Name | Family | Alone delta Pearson k "
        "| Combined delta | Lost when removed | Redundant |",
        "|---|---|---|---|---|---|---|---|",
    ]
    for r in ranking_rows:
        lines.append(
            f"| {int(r['rank'])} | {r['candidate_code']} "
            f"| {r['candidate_name']} | {r['candidate_family']} "
            f"| {_fmt(r['alone_contribution_pearson_kappa'], '+.5f')} "
            f"| {_fmt(r['combined_contribution_pearson_kappa'], '+.5f')} "
            f"| {_fmt(r['lost_when_removed_pearson_kappa'], '+.5f')} "
            f"| {'YES' if r['is_redundant'] else 'NO'} |"
        )

    lines += [
        "",
        "## Contribution analysis",
        "",
        "Decomposition of C10 improvement (vs frozen Version A control):",
        "",
        "| Source | Delta Pearson k | Delta SSIM k | Delta kappa bias |",
        "|---|---|---|---|",
    ]
    full = by_code["C10-C"]
    a = by_code["C10-A"]
    b = by_code["C10-B"]
    ctrl = control_row
    a_d_pk = a["median_pearson_kappa"] - ctrl["median_pearson_kappa"]
    b_d_pk = b["median_pearson_kappa"] - ctrl["median_pearson_kappa"]
    ab_d_pk = full["median_pearson_kappa"] - ctrl["median_pearson_kappa"]
    interaction_pk = ab_d_pk - (a_d_pk + b_d_pk)
    a_d_sk = a["median_ssim_kappa"] - ctrl["median_ssim_kappa"]
    b_d_sk = b["median_ssim_kappa"] - ctrl["median_ssim_kappa"]
    ab_d_sk = full["median_ssim_kappa"] - ctrl["median_ssim_kappa"]
    interaction_sk = ab_d_sk - (a_d_sk + b_d_sk)
    a_d_bk = a["mean_kappa_bias"] - ctrl["mean_kappa_bias"]
    b_d_bk = b["mean_kappa_bias"] - ctrl["mean_kappa_bias"]
    ab_d_bk = full["mean_kappa_bias"] - ctrl["mean_kappa_bias"]
    interaction_bk = ab_d_bk - (a_d_bk + b_d_bk)
    lines += [
        f"| A alone (Coherence) | {a_d_pk:+.5f} | {a_d_sk:+.5f} | {a_d_bk:+.5f} |",
        f"| B alone (Memory)    | {b_d_pk:+.5f} | {b_d_sk:+.5f} | {b_d_bk:+.5f} |",
        f"| Sum A + B           | {a_d_pk + b_d_pk:+.5f} | {a_d_sk + b_d_sk:+.5f} | {a_d_bk + b_d_bk:+.5f} |",
        f"| Combined (C10-C)    | {ab_d_pk:+.5f} | {ab_d_sk:+.5f} | {ab_d_bk:+.5f} |",
        f"| Interaction         | {interaction_pk:+.5f} | {interaction_sk:+.5f} | {interaction_bk:+.5f} |",
        "",
    ]

    lines += [
        "## Required questions",
        "",
    ] + q_answers + [
        "",
        "## Outcome determination",
        "",
        outcome,
        "",
        "## Numerical stability report",
        "",
        "| Variant | Median runtime (s) | Max conservation |",
        "|---|---|---|",
    ]
    for r in cross_rows:
        lines.append(
            f"| {r['candidate_code']} {r['candidate_name']} "
            f"| {_fmt(r['median_runtime_seconds'], '.3f')} "
            f"| {_fmt(r['max_conservation_error'], '.3e')} |"
        )

    lines += [
        "",
        "## Top-level artefacts",
        "",
        "- `runs/version_b_physics_lab002/report.md` (this file)",
        "- `runs/version_b_physics_lab002/component_contributions.csv`",
        "- `runs/version_b_physics_lab002/interaction_matrix.csv`",
        "- `runs/version_b_physics_lab002/leave_one_out.csv`",
        "- `runs/version_b_physics_lab002/family_ranking.csv`",
        "- `runs/version_b_physics_lab002/run.json`",
        "- `runs/version_b_physics_lab002/validation.json`",
        "- `runs/version_b_physics_lab002/plots/component_importance.png`",
        "- `runs/version_b_physics_lab002/plots/interaction_heatmap.png`",
        "- `runs/version_b_physics_lab002/plots/leave_one_out.png`",
        "- `runs/version_b_physics_lab002/plots/synergy_matrix.png`",
        "- `runs/version_b_physics_lab002/plots/family_contributions.png`",
        "",
    ]
    total_seconds = time.perf_counter() - overall_started
    lines += [f"**Total execution time:** {total_seconds:.1f} s.", ""]

    (out_root / "report.md").write_text("\n".join(lines))

    run_doc = {
        "milestone": "PBUF VERSION-B PHYSICS-LAB-002",
        "kind": "response_family_decomposition",
        "frozen_laboratory": "Version 1 weak-lensing laboratory (LAB-FREEZE-001)",
        "frozen_implementation_sha256": hash_check["files"],
        "production_minimum": MIN_PRODUCTION,
        "target_candidate": "C10 (Combined Local Response) from PHYSICS-LAB-001",
        "mechanisms_in_c10": [
            {"code": "A", "name": "Neighbour Coherence",
             "implementation": "multiplicative factor over 8 nearest neighbours"},
            {"code": "B", "name": "Elastic Memory",
             "implementation": "one-step persistence mix, w = 0.5"},
        ],
        "variants": [
            {"code": c.code, "name": c.name, "family": c.family,
             "description": c.description,
             "mechanism_a": c.coherence, "mechanism_b": c.memory,
             "notes": c.notes}
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
        "milestone": "PBUF VERSION-B PHYSICS-LAB-002",
        "frozen_hash_verification_passed": hash_check["ok"],
        "frozen_hashes": hash_check["files"],
        "all_clusters_machine_precision_conservation": all(
            r["max_conservation_error"] <= cons_eps + 1e-30
            for r in per_cluster_rows
        ),
        "variants_preserving_conservation": n_stable,
        "variants_total": len(cross_rows),
        "validation_passed": hash_check["ok"],
    }
    (out_root / "validation.json").write_text(json.dumps(validation_doc, indent=2))

    print(f"\nVERSION-B PHYSICS-LAB-002 COMPLETE  ({total_seconds:.1f} s)")
    print(json.dumps({
        "milestone": "PBUF VERSION-B PHYSICS-LAB-002",
        "status": "OK",
        "variants": len(CANDIDATES),
        "clusters": len(CLUSTERS),
        "frozen_hashes_ok": hash_check["ok"],
        "variants_preserving_conservation": n_stable,
        "output": str(out_root),
        "execution_seconds": float(total_seconds),
    }, indent=2))


if __name__ == "__main__":
    main()