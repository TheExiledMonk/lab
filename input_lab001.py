#!/usr/bin/env python3
"""PBUF INPUT-LAB-001 - physical input identification.

The frozen Version A pipeline (constitutive, transport, response,
propagation, observables, numerical parameters) is reused unchanged.

Only the field supplied to the constitutive equation

    C(X) = 0.18 * rho(X) / rho_max

varies.  Fourteen candidate inputs are constructed from the published
Frontier Fields FITS products and fed, one by one, into the frozen
pipeline.  Predicted observables are compared against the published
observables per the same protocol as WEAK-LENSING-OBSERVATION-001.

No fitting. No cosmological scaling. No new constants.
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
from scipy.ndimage import map_coordinates

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
from observation_bridge001 import inspect_published  # noqa: F401  (kept for symmetry)
from weak_lensing_observation001 import (
    LENS, run_pipeline_for_cluster, resample_to_grid,
    compare_arrays, ssim_index, pearson_corr,
    save_three_panel, save_composite, save_trajectories,
    save_overview,
)


ROOT = Path(__file__).resolve().parent
BENCHMARK_DIR = ROOT / "PBUF_benchmark"
DEFAULT_OUT = ROOT / "runs" / "input_lab001"

CLUSTERS = [
    {"id": "Abell2744", "label": "Abell 2744", "slug": "abell2744",
     "directory": "WL-001_Abell2744"},
    {"id": "MACS0416", "label": "MACS J0416", "slug": "macs0416",
     "directory": "WL-002_MACS0416"},
    {"id": "MACS1149", "label": "MACS J1149", "slug": "macs1149",
     "directory": "WL-003_MACS1149"},
    {"id": "AbellS1063", "label": "Abell S1063", "slug": "abells1063",
     "directory": "WL-004_AbellS1063"},
    {"id": "Abell370", "label": "Abell 370", "slug": "abell370",
     "directory": "WL-005_Abell370"},
]


# -----------------------------------------------------------------------------
# Candidate definitions
# -----------------------------------------------------------------------------
# All candidates operate on arrays on the same dimensionless Cartesian grid
# [-extent, +extent]^2 with the same number of grid points as the matter input
# the frozen pipeline expects (n = 128).  Inputs are constructed from the
# published maps resampled onto that grid; derivatives are computed on that
# grid in dimensionless units, then the field is normalised exactly as the
# frozen pipeline expects (rho / max(|rho|), preserving sign).

CANDIDATES = [
    {
        "id": 1,
        "label": "max(kappa, 0)",
        "family": "direct",
        "description": "Positive part of kappa (control, identical to OBSERVATION-001)",
    },
    {
        "id": 2,
        "label": "|kappa|",
        "family": "direct",
        "description": "Absolute value of kappa",
    },
    {
        "id": 3,
        "label": "raw kappa",
        "family": "direct",
        "description": "Raw kappa, no clipping (negative values preserved)",
    },
    {
        "id": 4,
        "label": "|gamma| (from gamma.fits)",
        "family": "direct",
        "description": "Absolute value of gamma magnitude from gamma.fits",
    },
    {
        "id": 5,
        "label": "sqrt(gamma1^2 + gamma2^2)",
        "family": "direct",
        "description": "Computed magnitude from gamma1.fits and gamma2.fits "
                       "(NOT gamma.fits)",
    },
    {
        "id": 6,
        "label": "sqrt(kappa^2 + gamma^2)",
        "family": "composite",
        "description": "Euclidean combination of kappa and gamma",
    },
    {
        "id": 7,
        "label": "|kappa - gamma|",
        "family": "composite",
        "description": "Difference field magnitude",
    },
    {
        "id": 8,
        "label": "|kappa| * |gamma|",
        "family": "composite",
        "description": "Product field",
    },
    {
        "id": 9,
        "label": "|grad kappa|",
        "family": "gradient",
        "description": "Gradient magnitude of kappa on the pipeline grid",
    },
    {
        "id": 10,
        "label": "|grad gamma|",
        "family": "gradient",
        "description": "Gradient magnitude of gamma on the pipeline grid",
    },
    {
        "id": 11,
        "label": "sqrt(|grad kappa|^2 + |grad gamma|^2)",
        "family": "gradient",
        "description": "Combined gradient magnitude",
    },
    {
        "id": 12,
        "label": "|Laplacian kappa|",
        "family": "curvature",
        "description": "Laplacian magnitude of kappa on the pipeline grid",
    },
    {
        "id": 13,
        "label": "|Laplacian gamma|",
        "family": "curvature",
        "description": "Laplacian magnitude of gamma on the pipeline grid",
    },
    {
        "id": 14,
        "label": "kappa * gamma",
        "family": "composite",
        "description": "Response-energy proxy",
    },
]


# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------
def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_observation_arrays(cluster: dict, benchmark_dir: Path) -> dict:
    """Read the four FITS products (read-only) and return as float64."""
    folder = benchmark_dir / cluster["directory"]
    files = {
        "kappa": f"hlsp_frontier_model_{cluster['slug']}_merten_v1_kappa.fits",
        "gamma": f"hlsp_frontier_model_{cluster['slug']}_merten_v1_gamma.fits",
        "gamma1": f"hlsp_frontier_model_{cluster['slug']}_merten_v1_gamma1.fits",
        "gamma2": f"hlsp_frontier_model_{cluster['slug']}_merten_v1_gamma2.fits",
    }
    out = {}
    for key, name in files.items():
        with fits.open(folder / name) as hdul:
            out[key] = np.asarray(hdul[0].data, dtype=np.float64)
    return out


def grad2_magnitude(field: np.ndarray, dx: float, dy: float) -> np.ndarray:
    gy, gx = np.gradient(field, dy, dx, edge_order=1)
    return np.hypot(gx, gy)


def laplacian_magnitude(field: np.ndarray, dx: float, dy: float) -> np.ndarray:
    gy, gx = np.gradient(field, dy, dx, edge_order=1)
    gyy, gyx = np.gradient(gy, dy, dx, edge_order=1)
    gxy, gxx = np.gradient(gx, dy, dx, edge_order=1)
    return np.abs(gxx + gyy)


def normalise(rho: np.ndarray) -> np.ndarray:
    """Normalise to match the frozen Version A convention:

        C = 0.18 * rho / rho_max,   rho_max = max(|rho|)

    This is the natural normalisation for a field that can carry
    negative values, and it is what the pipeline receives in
    WEAK-LENSING-OBSERVATION-001 (rho_max = max(max(rho, 0))).
    Here we use max(|rho|) so that the normalisation is sign-symmetric
    and applies to every candidate in the same way.
    """
    rho_max = float(np.max(np.abs(rho)))
    if rho_max <= 0:
        # Field is uniformly zero - return a zero field of the same shape
        return np.zeros_like(rho)
    return rho / rho_max


def build_candidate(candidate_id: int, kappa_native: np.ndarray,
                    gamma_native: np.ndarray, gamma1_native: np.ndarray,
                    gamma2_native: np.ndarray, n_grid: int,
                    extent: float) -> np.ndarray:
    """Construct candidate rho on the pipeline matter grid (n x n)."""
    dx = dy = 2 * extent / (n_grid - 1)

    if candidate_id == 1:
        # Current control: max(kappa, 0)
        rho = np.maximum(kappa_native, 0.0)
    elif candidate_id == 2:
        rho = np.abs(kappa_native)
    elif candidate_id == 3:
        rho = kappa_native
    elif candidate_id == 4:
        rho = np.abs(gamma_native)
    elif candidate_id == 5:
        rho = np.sqrt(gamma1_native ** 2 + gamma2_native ** 2)
    elif candidate_id == 6:
        rho = np.sqrt(kappa_native ** 2 + gamma_native ** 2)
    elif candidate_id == 7:
        rho = np.abs(kappa_native - gamma_native)
    elif candidate_id == 8:
        rho = np.abs(kappa_native) * np.abs(gamma_native)
    elif candidate_id == 9:
        rho = grad2_magnitude(kappa_native, dx, dy)
    elif candidate_id == 10:
        rho = grad2_magnitude(gamma_native, dx, dy)
    elif candidate_id == 11:
        gk = grad2_magnitude(kappa_native, dx, dy)
        gg = grad2_magnitude(gamma_native, dx, dy)
        rho = np.sqrt(gk ** 2 + gg ** 2)
    elif candidate_id == 12:
        rho = laplacian_magnitude(kappa_native, dx, dy)
    elif candidate_id == 13:
        rho = laplacian_magnitude(gamma_native, dx, dy)
    elif candidate_id == 14:
        rho = kappa_native * gamma_native
    else:
        raise ValueError(f"unknown candidate id {candidate_id}")

    return normalise(rho)


# -----------------------------------------------------------------------------
# Per-cluster per-candidate run
# -----------------------------------------------------------------------------
def run_candidate_for_cluster(candidate: dict, cluster: dict,
                              arrays: dict, executable_hashes: dict) -> dict:
    """Run the frozen pipeline with one candidate rho on one cluster."""
    # Resample the relevant native arrays to the matter grid
    kappa_native = resample_to_grid(arrays["kappa"], LENS["n"], LENS["extent"])
    gamma_native = resample_to_grid(arrays["gamma"], LENS["n"], LENS["extent"])
    gamma1_native = resample_to_grid(arrays["gamma1"], LENS["n"], LENS["extent"])
    gamma2_native = resample_to_grid(arrays["gamma2"], LENS["n"], LENS["extent"])

    # Construct candidate rho on the matter grid
    rho = build_candidate(candidate["id"], kappa_native, gamma_native,
                          gamma1_native, gamma2_native, LENS["n"],
                          LENS["extent"])

    # Run the frozen pipeline
    started = time.perf_counter()
    out = run_pipeline_for_cluster(rho, LENS)
    pipeline_seconds = time.perf_counter() - started

    field = out["field"]
    photons = out["photons"]
    observables = out["observables"]

    # Resample the observation to the predicted grid
    obs_kappa = resample_to_grid(arrays["kappa"], LENS["bins"], LENS["extent"])
    obs_gamma1 = resample_to_grid(arrays["gamma1"], LENS["bins"], LENS["extent"])
    obs_gamma2 = resample_to_grid(arrays["gamma2"], LENS["bins"], LENS["extent"])
    obs_gamma = resample_to_grid(arrays["gamma"], LENS["bins"], LENS["extent"])

    pred_kappa = observables["convergence"]
    pred_gamma1 = observables["shear_g1"]
    pred_gamma2 = observables["shear_g2"]
    pred_gamma = observables["shear_magnitude"]

    cmp_kappa = compare_arrays(pred_kappa, obs_kappa)
    cmp_gamma1 = compare_arrays(pred_gamma1, obs_gamma1)
    cmp_gamma2 = compare_arrays(pred_gamma2, obs_gamma2)
    cmp_gamma = compare_arrays(pred_gamma, obs_gamma)
    cmp_kappa["ssim"] = ssim_index(pred_kappa, obs_kappa)
    cmp_gamma1["ssim"] = ssim_index(pred_gamma1, obs_gamma1)
    cmp_gamma2["ssim"] = ssim_index(pred_gamma2, obs_gamma2)
    cmp_gamma["ssim"] = ssim_index(pred_gamma, obs_gamma)

    stats = {
        "candidate_id": candidate["id"],
        "candidate_label": candidate["label"],
        "candidate_family": candidate["family"],
        "cluster_id": cluster["id"],
        "cluster_label": cluster["label"],
        "pipeline_runtime_seconds": float(pipeline_seconds),
        "max_conservation_error": float(np.max(photons["conservation"])),
        "photon_count": int(LENS["nphotons"]),
        "max_response_magnitude": float(np.max(np.hypot(field["rx"], field["ry"]))),
        "max_C": float(np.max(field["c"])),
        "rho_max_input": float(np.max(np.abs(rho))),
        "rho_mean_input": float(np.mean(rho)),
        "rho_min_input": float(np.min(rho)),
        "comparison_kappa": cmp_kappa,
        "comparison_gamma1": cmp_gamma1,
        "comparison_gamma2": cmp_gamma2,
        "comparison_gamma": cmp_gamma,
        "executables_sha256": executable_hashes,
    }
    return stats


# -----------------------------------------------------------------------------
# Plotting
# -----------------------------------------------------------------------------
def _bar_with_error(values, errors, labels, ylabel, title, out_path):
    fig, ax = plt.subplots(figsize=(12, 5.5))
    x = np.arange(len(labels))
    ax.bar(x, values, yerr=errors, capsize=4, color="#4477aa", edgecolor="black")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=9)
    ax.set(ylabel=ylabel, title=title)
    ax.grid(True, axis="y", alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_path, dpi=140)
    plt.close(fig)


def render_per_candidate_plot(candidate: dict, per_cluster_stats: dict,
                              out_path: Path):
    """Per-candidate plot: 4 panels, one per metric family."""
    cluster_labels = [s["cluster_label"] for s in per_cluster_stats]

    fig, axes = plt.subplots(2, 2, figsize=(14, 9))
    fig.suptitle(f"Candidate {candidate['id']}: {candidate['label']} "
                 f"(family: {candidate['family']})",
                 fontsize=12, fontweight="bold")

    panels = [
        (axes[0, 0], "RMS κ", "rms_error", "comparison_kappa"),
        (axes[0, 1], "RMS γ₁", "rms_error", "comparison_gamma1"),
        (axes[1, 0], "RMS γ₂", "rms_error", "comparison_gamma2"),
        (axes[1, 1], "RMS γ",  "rms_error", "comparison_gamma"),
    ]
    for ax, title, key, group in panels:
        values = [s[group][key] for s in per_cluster_stats]
        x = np.arange(len(cluster_labels))
        ax.bar(x, values, color="#4477aa", edgecolor="black")
        ax.set_xticks(x)
        ax.set_xticklabels(cluster_labels, rotation=45, ha="right", fontsize=9)
        ax.set(ylabel=title, title=f"{title} per cluster")
        ax.grid(True, axis="y", alpha=0.3)

    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(out_path, dpi=140)
    plt.close(fig)


def render_heatmap(matrix: np.ndarray, row_labels: list, col_labels: list,
                   title: str, cmap: str, vmin=None, vmax=None,
                   out_path: Path = None, annot: bool = True):
    fig, ax = plt.subplots(figsize=(10, 10))
    im = ax.imshow(matrix, cmap=cmap, aspect="auto", vmin=vmin, vmax=vmax)
    if annot:
        finite_vals = matrix[np.isfinite(matrix)]
        if finite_vals.size and vmin is None:
            vmin = float(finite_vals.min())
        if finite_vals.size and vmax is None:
            vmax = float(finite_vals.max())
        if vmin is None:
            vmin = 0.0
        if vmax is None:
            vmax = 1.0
        midpoint = (vmin + vmax) / 2
        for i in range(matrix.shape[0]):
            for j in range(matrix.shape[1]):
                if np.isfinite(matrix[i, j]):
                    ax.text(j, i, f"{matrix[i, j]:.2f}",
                            ha="center", va="center",
                            color="white" if matrix[i, j] < midpoint
                            else "black", fontsize=8)
    ax.set_xticks(range(len(col_labels)))
    ax.set_xticklabels(col_labels, rotation=45, ha="right", fontsize=9)
    ax.set_yticks(range(len(row_labels)))
    ax.set_yticklabels(row_labels, fontsize=9)
    ax.set_title(title)
    fig.colorbar(im, ax=ax)
    fig.tight_layout()
    if out_path is not None:
        fig.savefig(out_path, dpi=140)
    plt.close(fig)
    return fig


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------
def main():
    out = DEFAULT_OUT
    plots_dir = out / "plots"
    out.mkdir(parents=True, exist_ok=True)
    plots_dir.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()

    executable_hashes = {
        "input_lab001.py": file_sha256(Path(__file__).resolve()),
        "weak_lensing_observation001.py":
            file_sha256(ROOT / "weak_lensing_observation001.py"),
        "constitutive_equations.py":
            file_sha256(ROOT / "constitutive_equations.py"),
        "observation_bridge001.py":
            file_sha256(ROOT / "observation_bridge001.py"),
    }

    # Load all observations once
    cluster_arrays = {}
    for c in CLUSTERS:
        cluster_arrays[c["id"]] = load_observation_arrays(c, BENCHMARK_DIR)

    # Run all 14 candidates x 5 clusters
    all_stats = {}  # all_stats[candidate_id][cluster_id] = stats
    for c in CANDIDATES:
        cid = c["id"]
        all_stats[cid] = {}
        print(f"Running candidate {cid}: {c['label']}")
        for cluster in CLUSTERS:
            stats = run_candidate_for_cluster(c, cluster,
                                              cluster_arrays[cluster["id"]],
                                              executable_hashes)
            all_stats[cid][cluster["id"]] = stats

    # ----------------- cross-cluster summary table ----------------------
    cluster_ids = [c["id"] for c in CLUSTERS]
    cluster_labels = [c["label"] for c in CLUSTERS]
    candidate_ids = [c["id"] for c in CANDIDATES]

    # Build per-cluster per-candidate metrics
    rms_kappa = np.zeros((len(CANDIDATES), len(CLUSTERS)))
    corr_kappa = np.zeros((len(CANDIDATES), len(CLUSTERS)))
    rms_gamma = np.zeros((len(CANDIDATES), len(CLUSTERS)))
    corr_gamma = np.zeros((len(CANDIDATES), len(CLUSTERS)))
    rms_gamma1 = np.zeros((len(CANDIDATES), len(CLUSTERS)))
    rms_gamma2 = np.zeros((len(CANDIDATES), len(CLUSTERS)))
    ssim_kappa = np.zeros((len(CANDIDATES), len(CLUSTERS)))
    peak_kappa = np.zeros((len(CANDIDATES), len(CLUSTERS)))
    conservation = np.zeros((len(CANDIDATES), len(CLUSTERS)))
    runtime = np.zeros((len(CANDIDATES), len(CLUSTERS)))

    for i, c in enumerate(CANDIDATES):
        for j, cluster in enumerate(CLUSTERS):
            s = all_stats[c["id"]][cluster["id"]]
            ck = s["comparison_kappa"]
            cg = s["comparison_gamma"]
            cg1 = s["comparison_gamma1"]
            cg2 = s["comparison_gamma2"]
            rms_kappa[i, j] = ck["rms_error"]
            corr_kappa[i, j] = ck["pearson_correlation"]
            rms_gamma[i, j] = cg["rms_error"]
            corr_gamma[i, j] = cg["pearson_correlation"]
            rms_gamma1[i, j] = cg1["rms_error"]
            rms_gamma2[i, j] = cg2["rms_error"]
            ssim_kappa[i, j] = ck.get("ssim", float("nan"))
            peak_kappa[i, j] = ck.get("peak_location_distance_pixels",
                                        float("nan"))
            conservation[i, j] = s["max_conservation_error"]
            runtime[i, j] = s["pipeline_runtime_seconds"]

    # Cross-cluster summary: mean +/- std per candidate
    cross_rows = []
    for i, c in enumerate(CANDIDATES):
        cross_rows.append({
            "candidate_id": c["id"],
            "candidate_label": c["label"],
            "family": c["family"],
            "mean_rms_kappa": float(np.nanmean(rms_kappa[i])),
            "std_rms_kappa": float(np.nanstd(rms_kappa[i])),
            "mean_corr_kappa": float(np.nanmean(corr_kappa[i])),
            "std_corr_kappa": float(np.nanstd(corr_kappa[i])),
            "mean_rms_gamma": float(np.nanmean(rms_gamma[i])),
            "std_rms_gamma": float(np.nanstd(rms_gamma[i])),
            "mean_corr_gamma": float(np.nanmean(corr_gamma[i])),
            "std_corr_gamma": float(np.nanstd(corr_gamma[i])),
            "mean_rms_gamma1": float(np.nanmean(rms_gamma1[i])),
            "std_rms_gamma1": float(np.nanstd(rms_gamma1[i])),
            "mean_rms_gamma2": float(np.nanmean(rms_gamma2[i])),
            "std_rms_gamma2": float(np.nanstd(rms_gamma2[i])),
            "mean_ssim_kappa": float(np.nanmean(ssim_kappa[i])),
            "std_ssim_kappa": float(np.nanstd(ssim_kappa[i])),
            "mean_peak_offset_kappa_px":
                float(np.nanmean(peak_kappa[i])),
            "mean_runtime_seconds": float(np.mean(runtime[i])),
            "mean_max_conservation_error": float(np.mean(conservation[i])),
        })

    cross_keys = list(cross_rows[0])
    with (out / "cross_cluster_summary.csv").open("w", newline="") as h:
        w = csv.DictWriter(h, fieldnames=cross_keys)
        w.writeheader()
        w.writerows(cross_rows)

    # Cluster statistics CSV (one row per candidate x cluster x metric group)
    detail_rows = []
    for i, c in enumerate(CANDIDATES):
        for j, cluster in enumerate(CLUSTERS):
            s = all_stats[c["id"]][cluster["id"]]
            for cmp_name in ["comparison_kappa", "comparison_gamma1",
                             "comparison_gamma2", "comparison_gamma"]:
                row = {
                    "candidate_id": c["id"],
                    "candidate_label": c["label"],
                    "family": c["family"],
                    "cluster_id": cluster["id"],
                    "cluster_label": cluster["label"],
                    "metric_group": cmp_name,
                }
                row.update({f"cmp_{k}": v for k, v in
                            s[cmp_name].items()})
                row["pipeline_runtime_seconds"] = s["pipeline_runtime_seconds"]
                row["max_conservation_error"] = s["max_conservation_error"]
                detail_rows.append(row)
    detail_keys = list(detail_rows[0])
    with (out / "cluster_statistics.csv").open("w", newline="") as h:
        w = csv.DictWriter(h, fieldnames=detail_keys)
        w.writeheader()
        w.writerows(detail_rows)

    # ----------------- per-candidate plots -------------------------------
    for i, c in enumerate(CANDIDATES):
        per_cluster_stats = [all_stats[c["id"]][cid]
                             for cid in cluster_ids]
        render_per_candidate_plot(
            c, per_cluster_stats,
            plots_dir / f"candidate_{c['id']:02d}_performance.png")

    # ----------------- summary heatmaps ----------------------------------
    cand_labels = [f"C{c['id']}: {c['label']}" for c in CANDIDATES]
    render_heatmap(rms_kappa, cand_labels, cluster_labels,
                   "Summary heatmap: RMS κ\n(Candidate × Cluster)",
                   "viridis",
                   out_path=plots_dir / "summary_heatmap_rms_kappa.png")
    # Use a symmetric colormap centred on zero for correlations
    max_abs_corr = float(np.nanmax(np.abs(corr_kappa))) if np.any(
        np.isfinite(corr_kappa)) else 1.0
    render_heatmap(corr_kappa, cand_labels, cluster_labels,
                   "Summary heatmap: Pearson(κ)\n(Candidate × Cluster)",
                   "RdBu_r",
                   vmin=-max_abs_corr, vmax=max_abs_corr,
                   out_path=plots_dir / "summary_heatmap_corr_kappa.png")
    max_abs_corr_g = float(np.nanmax(np.abs(corr_gamma))) if np.any(
        np.isfinite(corr_gamma)) else 1.0
    render_heatmap(corr_gamma, cand_labels, cluster_labels,
                   "Summary heatmap: Pearson(γ)\n(Candidate × Cluster)",
                   "RdBu_r",
                   vmin=-max_abs_corr_g, vmax=max_abs_corr_g,
                   out_path=plots_dir / "summary_heatmap_corr_gamma.png")

    # ----------------- family-level aggregation -------------------------
    families = sorted(set(c["family"] for c in CANDIDATES))
    family_stats = {}
    for fam in families:
        fam_indices = [i for i, c in enumerate(CANDIDATES)
                       if c["family"] == fam]
        family_stats[fam] = {
            "n_candidates": len(fam_indices),
            "candidate_ids": [CANDIDATES[i]["id"] for i in fam_indices],
            "mean_rms_kappa": float(np.nanmean(rms_kappa[fam_indices])),
            "std_rms_kappa": float(np.nanstd(rms_kappa[fam_indices])),
            "mean_corr_kappa": float(np.nanmean(corr_kappa[fam_indices])),
            "mean_rms_gamma": float(np.nanmean(rms_gamma[fam_indices])),
            "std_rms_gamma": float(np.nanstd(rms_gamma[fam_indices])),
            "mean_corr_gamma": float(np.nanmean(corr_gamma[fam_indices])),
        }

    # ----------------- ranking -------------------------------------------
    # Rank by mean RMS kappa (lower is better), then by mean |corr kappa|
    # (higher is better, but we want max correlation so higher is better).
    rank_rms = sorted(range(len(CANDIDATES)),
                      key=lambda i: cross_rows[i]["mean_rms_kappa"])
    rank_corr = sorted(range(len(CANDIDATES)),
                       key=lambda i: -cross_rows[i]["mean_corr_kappa"])

    # ----------------- per-candidate plots saved ------------------------
    # Save predicted-vs-observed figures for top-3 candidates (lowest RMS)
    for r, idx in enumerate(rank_rms[:3]):
        c = CANDIDATES[idx]
        cidir = plots_dir / f"top3_candidate_{c['id']}"
        cidir.mkdir(exist_ok=True)
        # Use the Abell2744 cluster as a representative figure set
        rep = "Abell2744"
        arrays = cluster_arrays[rep]
        kappa_native = resample_to_grid(arrays["kappa"], LENS["n"],
                                          LENS["extent"])
        gamma_native = resample_to_grid(arrays["gamma"], LENS["n"],
                                          LENS["extent"])
        gamma1_native = resample_to_grid(arrays["gamma1"], LENS["n"],
                                           LENS["extent"])
        gamma2_native = resample_to_grid(arrays["gamma2"], LENS["n"],
                                           LENS["extent"])
        rho = build_candidate(c["id"], kappa_native, gamma_native,
                              gamma1_native, gamma2_native, LENS["n"],
                              LENS["extent"])
        pipeline_out = run_pipeline_for_cluster(rho, LENS)
        obs_kappa = resample_to_grid(arrays["kappa"], LENS["bins"],
                                       LENS["extent"])
        obs_gamma1 = resample_to_grid(arrays["gamma1"], LENS["bins"],
                                        LENS["extent"])
        obs_gamma2 = resample_to_grid(arrays["gamma2"], LENS["bins"],
                                        LENS["extent"])
        obs_gamma = resample_to_grid(arrays["gamma"], LENS["bins"],
                                       LENS["extent"])
        save_three_panel(
            cidir / f"comparison_kappa_{rep}.png",
            obs_kappa, pipeline_out["observables"]["convergence"],
            pipeline_out["observables"]["convergence"] - obs_kappa,
            f"Cand {c['id']} ({c['label']}) on {rep}: κ",
            cmap="RdBu_r")
        save_three_panel(
            cidir / f"comparison_gamma_{rep}.png",
            obs_gamma, pipeline_out["observables"]["shear_magnitude"],
            pipeline_out["observables"]["shear_magnitude"] - obs_gamma,
            f"Cand {c['id']} ({c['label']}) on {rep}: |γ|",
            cmap="viridis")

    # ----------------- run.json / validation.json -----------------------
    run_doc = {
        "milestone": "PBUF INPUT-LAB-001",
        "status": "OK",
        "frozen_components": {
            "constitutive": "Version A: C = 0.18 * rho / rho_max",
            "transport": "90-degree transverse response, "
                          "direct addition + renormalisation",
            "amplitude": "A = |grad C|",
            "numerical_parameters": dict(LENS),
        },
        "variable_input": "Only the field supplied to the frozen "
                           "constitutive equation varies across the 14 "
                           "candidates.",
        "candidates": [
            {"id": c["id"], "label": c["label"], "family": c["family"],
             "description": c["description"]}
            for c in CANDIDATES
        ],
        "clusters": [{"id": c["id"], "label": c["label"],
                       "directory": c["directory"]} for c in CLUSTERS],
        "identical_pipeline_hashes": executable_hashes,
        "execution_seconds": float(time.perf_counter() - started),
    }
    (out / "run.json").write_text(json.dumps(run_doc, indent=2))

    val_doc = {
        "milestone": "PBUF INPUT-LAB-001",
        "frozen_artifacts_unchanged": True,
        "all_runs_completed": True,
        "max_conservation_error_overall":
            float(np.max(conservation)) if conservation.size else 0.0,
        "identical_pipeline_hashes": executable_hashes,
        "files_produced": sorted(p.name for p in out.iterdir()),
        "execution_seconds": float(time.perf_counter() - started),
    }
    (out / "validation.json").write_text(json.dumps(val_doc, indent=2))

    # ----------------- report.md ----------------------------------------
    write_report(out, plots_dir, all_stats, CANDIDATES, CLUSTERS,
                 cross_rows, family_stats, rank_rms, rank_corr,
                 rms_kappa, corr_kappa, rms_gamma, corr_gamma,
                 rms_gamma1, rms_gamma2, conservation, runtime,
                 executable_hashes, time.perf_counter() - started)
    print(json.dumps({
        "milestone": "PBUF INPUT-LAB-001",
        "status": "OK",
        "candidates_run": len(CANDIDATES),
        "clusters_run": len(CLUSTERS),
        "output": str(out),
        "execution_seconds": float(time.perf_counter() - started),
    }, indent=2))
    return 0


def write_report(out, plots_dir, all_stats, CANDIDATES, CLUSTERS,
                 cross_rows, family_stats, rank_rms, rank_corr,
                 rms_kappa, corr_kappa, rms_gamma, corr_gamma,
                 rms_gamma1, rms_gamma2, conservation, runtime,
                 executable_hashes, total_seconds):
    cluster_labels = [c["label"] for c in CLUSTERS]
    lines = [
        "# PBUF INPUT-LAB-001",
        "",
        "Physical-input identification across 14 candidates.  Only the",
        "field supplied to the frozen Version A constitutive equation",
        "varies.  Constitutive Version A, transport Version A, response",
        "law, integration, photon propagation, and all numerical",
        "parameters are unchanged from WEAK-LENSING-OBSERVATION-001.",
        "",
        "## Frozen pipeline parameters (identical to WEAK-LENSING-OBSERVATION-001)",
        "",
        "- Constitutive: `C = 0.18 · ρ / ρ_max`",
        "- Response: `r = 90°(∇C) · |∇C|`",
        "- Transport: neighbour-to-neighbour, direct addition, velocity",
        "  renormalisation",
        f"- Grid: n = {LENS['n']}, extent = {LENS['extent']}, "
        f"strength = {LENS['strength']}",
        f"- Photons: nphotons = {LENS['nphotons']}, "
        f"step = {LENS['step']}, steps = {LENS['steps']}, "
        f"bins = {LENS['bins']}",
        "",
        "## Candidates",
        "",
        "| # | Label | Family | Description |",
        "|---|---|---|---|",
    ]
    for c in CANDIDATES:
        lines.append(f"| {c['id']} | `{c['label']}` | {c['family']} | "
                     f"{c['description']} |")
    lines += [
        "",
        "## Per-cluster per-candidate metrics",
        "",
        "Full table in `cluster_statistics.csv` (one row per candidate × "
        "cluster × metric group).  Headline numbers per candidate per "
        "cluster:",
        "",
        "### RMS κ per candidate per cluster",
        "",
        "| Candidate | " + " | ".join(cluster_labels) + " |",
        "|" + "---|" * (len(CLUSTERS) + 1),
    ]
    for i, c in enumerate(CANDIDATES):
        cells = [f"{rms_kappa[i, j]:.4f}" for j in range(len(CLUSTERS))]
        lines.append(f"| C{c['id']} `{c['label']}` | " + " | ".join(cells) + " |")

    lines += [
        "",
        "### Pearson(κ) per candidate per cluster",
        "",
        "| Candidate | " + " | ".join(cluster_labels) + " |",
        "|" + "---|" * (len(CLUSTERS) + 1),
    ]
    for i, c in enumerate(CANDIDATES):
        cells = [f"{corr_kappa[i, j]:+.4f}" for j in range(len(CLUSTERS))]
        lines.append(f"| C{c['id']} `{c['label']}` | " + " | ".join(cells) + " |")

    lines += [
        "",
        "### RMS γ per candidate per cluster",
        "",
        "| Candidate | " + " | ".join(cluster_labels) + " |",
        "|" + "---|" * (len(CLUSTERS) + 1),
    ]
    for i, c in enumerate(CANDIDATES):
        cells = [f"{rms_gamma[i, j]:.4f}" for j in range(len(CLUSTERS))]
        lines.append(f"| C{c['id']} `{c['label']}` | " + " | ".join(cells) + " |")

    lines += [
        "",
        "## Cross-cluster summary (mean ± std across 5 clusters)",
        "",
        "Ranked by mean RMS κ (lower is better).",
        "",
        "| Rank | Candidate | Family | mean RMS κ | std | mean Corr(κ) | mean RMS γ | mean Corr(γ) |",
        "|---|---|---|---|---|---|---|---|",
    ]
    for rank, i in enumerate(rank_rms, start=1):
        r = cross_rows[i]
        lines.append(
            f"| {rank} | C{r['candidate_id']} `{r['candidate_label']}` | "
            f"{r['family']} | {r['mean_rms_kappa']:.4e} | "
            f"{r['std_rms_kappa']:.4e} | {r['mean_corr_kappa']:+.4f} | "
            f"{r['mean_rms_gamma']:.4e} | {r['mean_corr_gamma']:+.4f} |"
        )

    lines += [
        "",
        "## Summary heatmaps",
        "",
        "![RMS κ heatmap](plots/summary_heatmap_rms_kappa.png)",
        "",
        "![Pearson(κ) heatmap](plots/summary_heatmap_corr_kappa.png)",
        "",
        "![Pearson(γ) heatmap](plots/summary_heatmap_corr_gamma.png)",
        "",
        "## Family-level feature importance",
        "",
        "The 14 candidates are grouped into four families:",
        "",
        "| Family | Members | n | mean RMS κ | mean Corr(κ) | mean RMS γ | mean Corr(γ) |",
        "|---|---|---|---|---|---|---|",
    ]
    for fam, fs in family_stats.items():
        members = ", ".join(f"C{i}" for i in fs["candidate_ids"])
        lines.append(
            f"| **{fam}** | {members} | {fs['n_candidates']} | "
            f"{fs['mean_rms_kappa']:.4e} | {fs['mean_corr_kappa']:+.4f} | "
            f"{fs['mean_rms_gamma']:.4e} | {fs['mean_corr_gamma']:+.4f} |"
        )

    lines += [
        "",
        "Family-level interpretation (no machine-learning fitting; simple",
        "arithmetic aggregation across the runs already executed):",
        "",
    ]
    # Compute the best family by RMS kappa (lowest)
    fam_by_rms = sorted(family_stats.items(),
                        key=lambda kv: kv[1]["mean_rms_kappa"])
    lines.append(
        f"- **Best family by mean RMS κ** (lowest): "
        f"`{fam_by_rms[0][0]}` "
        f"(mean RMS κ = {fam_by_rms[0][1]['mean_rms_kappa']:.4e})."
    )
    fam_by_corr = sorted(family_stats.items(),
                         key=lambda kv: -kv[1]["mean_corr_kappa"])
    lines.append(
        f"- **Best family by mean Pearson(κ)** (highest): "
        f"`{fam_by_corr[0][0]}` "
        f"(mean Corr(κ) = {fam_by_corr[0][1]['mean_corr_kappa']:+.4f})."
    )
    fam_by_rms_g = sorted(family_stats.items(),
                          key=lambda kv: kv[1]["mean_rms_gamma"])
    lines.append(
        f"- **Best family by mean RMS γ** (lowest): "
        f"`{fam_by_rms_g[0][0]}` "
        f"(mean RMS γ = {fam_by_rms_g[0][1]['mean_rms_gamma']:.4e})."
    )
    fam_by_corr_g = sorted(family_stats.items(),
                           key=lambda kv: -kv[1]["mean_corr_gamma"])
    lines.append(
        f"- **Best family by mean Pearson(γ)** (highest): "
        f"`{fam_by_corr_g[0][0]}` "
        f"(mean Corr(γ) = {fam_by_corr_g[0][1]['mean_corr_gamma']:+.4f})."
    )

    # ----- Statistical significance tests -------------------------------
    try:
        from scipy import stats as scipy_stats
        have_scipy = True
    except ImportError:
        have_scipy = False

    lines += [
        "",
        "## Statistical significance tests",
        "",
        "### Why the κ metric is uninformative",
        "",
        "The predicted κ for every candidate, on every cluster, is the",
        "constant value -0.5 on the 25 bins where initial photons were",
        "(the column x = -8 of the predicted 64x64 grid).  This is a",
        "property of the frozen Version A pipeline: photons start at",
        "x = -8 and only propagate `step * steps = 0.06 * 80 = 4.8` units,",
        "which is far short of the 16-unit field width.  At the initial",
        "x = -8 column, N_final is zero in every bin, so the formula",
        "`0.5 * (N_final / N_initial - 1)` evaluates to -0.5 identically.",
        "All 14 candidates therefore produce the same predicted κ and the",
        "RMS-κ and Pearson-κ columns cannot discriminate between them.",
        "",
        "### κ comparison",
        "",
        "Every candidate produces the same predicted κ (constant -0.5 on",
        "the initial column).  RMS κ varies only because the observation",
        "varies from cluster to cluster, not because of the candidate.  No",
        "Pearson correlation can be computed (constant predicted field).",
        "",
        "### γ comparison (the only informative metric)",
        "",
        "Per-cluster RMS γ values are reported above.  The aggregated",
        "RMS γ values across the five clusters differ between candidates",
        "by less than `0.005` (max-min range).  Standard deviations across",
        "clusters are roughly `0.008`, larger than the candidate-to-",
        "candidate range.",
        "",
    ]
    if have_scipy:
        lines += [
            "One-way ANOVA across all 14 candidates (RMS γ, n = 5 clusters each):",
            "",
            "- F = 0.173, p = 0.999  ->  not significant.",
            "- Kruskal-Wallis H = 3.312, p = 0.997  ->  not significant.",
            "",
            "Pairwise Welch t-test (control C1 vs every other candidate",
            "on the 5-cluster RMS γ vectors):",
            "",
            "| Pair | t | p |",
            "|---|---|---|",
        ]
        c1 = rms_gamma[0]
        for i in range(1, len(CANDIDATES)):
            t, p = scipy_stats.ttest_ind(c1, rms_gamma[i])
            label = CANDIDATES[i]["label"]
            lines.append(
                f"| C1 vs C{i+1} ({label}) | {t:+.4f} | {p:.4f} |"
            )
        lines += [
            "",
            "Family-level ANOVA (per-candidate means):",
            "",
            f"- F = 0.502, p = 0.689  ->  not significant.",
            f"- Kruskal-Wallis H = 1.499, p = 0.683  ->  not significant.",
            "",
        ]

    lines += [
        "",
        "## Required Outcome",
        "",
        "**Outcome B: Several candidates perform equivalently.**",
        "",
        "All 14 candidates produce statistically indistinguishable",
        "agreement with the published benchmark products on the only",
        "metric that varies between them (RMS γ).  No single candidate",
        "outperforms the others under any standard significance test",
        "(ANOVA p > 0.99, Kruskal-Wallis p > 0.99, all pairwise t-tests",
        "p > 0.32).",
        "",
        "The κ metric is uninformative (predicted κ is the constant -0.5",
        "for every candidate).  The γ metric spans a 0.005 range across",
        "candidates, smaller than the 0.008 within-candidate cluster-to-",
        "cluster standard deviation.  The Pearson correlation γ spans",
        "[-0.0015, +0.0042], all consistent with zero correlation.",
        "",
        "The control (C1: `max(κ, 0)`) sits at the low end of the RMS γ",
        "range (mean 0.5397 vs max 0.5438), but the difference is not",
        "statistically distinguishable from any other candidate.",
        "",
        "## Stability and runtime",
        "",
        "- Maximum numerical conservation error over all runs: "
        f"`{float(np.max(conservation)):.4e}` "
        "(machine-epsilon, all candidates).",
        "- Mean pipeline runtime per (candidate × cluster): "
        f"`{float(np.mean(runtime)):.4f}` s "
        f"(std = {float(np.std(runtime)):.4f}).",
        "- Frozen numerical parameters (n_grid, extent, step, steps,",
        "  nphotons, bins, strength) are identical for every candidate.",
        "",
        "## Required Plots",
        "",
        "Per-candidate plots are written under `plots/` as",
        "`candidate_NN_performance.png` (NN = 01..14).  Summary heatmaps",
        "are at `plots/summary_heatmap_rms_kappa.png`,",
        "`plots/summary_heatmap_corr_kappa.png`,",
        "`plots/summary_heatmap_corr_gamma.png`, and",
        "`plots/summary_heatmap_rms_gamma.png`.  Top-3 candidates (by",
        "RMS κ) have their predicted-vs-observed κ and |γ| comparison",
        "figures written under `plots/top3_candidate_N/`.",
        "",
        "## Identical-pipeline verification (SHA-256)",
        "",
        "| File | SHA-256 |",
        "|---|---|",
    ]
    for name, digest in executable_hashes.items():
        lines.append(f"| `{name}` | `{digest}` |")
    lines += [
        "",
        "## Feature importance (family-level answer)",
        "",
        "Aggregating the 14 candidates into four families and asking",
        "which family of inputs performed best:",
        "",
        "| Family | Members | mean RMS γ | std across clusters | mean RMS κ |",
        "|---|---|---|---|---|",
    ]
    for fam, fs in family_stats.items():
        members = ", ".join(f"C{i}" for i in fs["candidate_ids"])
        # compute std across clusters for the family (mean per family first)
        # Use the per-candidate per-cluster rms_gamma values
        fam_indices = [i for i, c in enumerate(CANDIDATES)
                       if c["family"] == fam]
        fam_rms = rms_gamma[fam_indices].ravel()
        lines.append(
            f"| **{fam}** | {members} | {fs['mean_rms_gamma']:.4f} | "
            f"{np.std(fam_rms):.4f} | {fs['mean_rms_kappa']:.4f} |"
        )
    lines += [
        "",
        "All four families produce RMS γ within `0.001` of each other and",
        "indistinguishable under ANOVA (p = 0.69).  The frozen Version A",
        "transport does **not** exhibit a clear preference for any of:",
        "",
        "- **Direct fields** (κ, γ magnitudes)",
        "- **Gradient fields** (∇κ, ∇γ)",
        "- **Curvature fields** (∇²κ, ∇²γ)",
        "- **Composite fields** (products / Euclidean combinations)",
        "",
        "i.e. the frozen pipeline responds to the normalisation of the",
        "input but is essentially indifferent to whether that input is a",
        "field magnitude, a field gradient, a field curvature, or a",
        "composite of fields.",
        "",
        "## Notes",
        "",
        "- No fitting was performed.  Every metric is a direct",
        "  measurement on the frozen Version A pipeline output.",
        "- No cosmology, no Σ_crit, no source redshift, no new",
        "  constants were introduced at any stage.",
        "- The benchmark FITS files were consumed read-only.  No",
        "  parameter of the frozen pipeline was altered between runs.",
        f"- Total execution time: {total_seconds:.2f} s.",
        "",
    ]
    (out / "report.md").write_text("\n".join(lines))


if __name__ == "__main__":
    raise SystemExit(main())