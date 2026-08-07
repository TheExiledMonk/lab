#!/usr/bin/env python3
"""PBUF FOUNDATION — M10 INFORMATION-LOSS AUDIT 001.

Five-cluster diagnostic audit of the already-requalified physical chain:

    M10 interface field
      -> native-z LOS
      -> ray propagation
      -> binned ray displacement
      -> ray-bundle Jacobian kappa

The purpose is to locate where morphology may be lost between the coordinate-
safe M10 physical field and the sparse Jacobian observable.  This lab does NOT
change physics, tune parameters, fit amplitudes, choose candidates, or change
source-plane/ray settings.

The exact primary-candidate construction and frozen ray/Jacobian path are reused
from primary_candidate_science_rerun_m10_001.  Correlations are diagnostics,
never pass/fail gates.
"""
from __future__ import annotations

import csv
import json
import math
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))

from a8_three_dimensional_projection_lab001 import CLUSTERS
from pbuf.labs.foundation import primary_candidate_science_rerun_m10_001 as PRIMARY

LAB_ID = "PBUF-FOUNDATION-M10-INFORMATION-LOSS-AUDIT-001"
OUT = ROOT / "runs" / "m10_information_loss_audit001"
EXPECTED_CLUSTERS = 5


def _git(*args: str) -> str:
    return subprocess.check_output(["git", *args], cwd=str(ROOT), text=True).strip()


def _repo_state() -> dict:
    return {
        "repository": "TheExiledMonk/lab",
        "branch": _git("rev-parse", "--abbrev-ref", "HEAD"),
        "head_sha": _git("rev-parse", "HEAD"),
        "tracked_changes": _git("diff", "--name-only"),
        "staged_changes": _git("diff", "--name-only", "--cached"),
    }


def _json_default(obj):
    if isinstance(obj, np.floating):
        return float(obj)
    if isinstance(obj, np.integer):
        return int(obj)
    if isinstance(obj, np.bool_):
        return bool(obj)
    if isinstance(obj, Path):
        return str(obj)
    if isinstance(obj, tuple):
        return list(obj)
    return str(obj)


def _write_json(path: Path, obj) -> None:
    path.write_text(json.dumps(obj, indent=2, default=_json_default))


def _write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        path.write_text("")
        return
    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)


def _finite_pair(a, b, extra_mask=None):
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    if a.shape != b.shape:
        raise ValueError(f"shape mismatch {a.shape} != {b.shape}")
    mask = np.isfinite(a) & np.isfinite(b)
    if extra_mask is not None:
        m = np.asarray(extra_mask, dtype=bool)
        if m.shape != a.shape:
            raise ValueError(f"mask shape mismatch {m.shape} != {a.shape}")
        mask &= m
    return a[mask], b[mask], mask


def _pearson(a, b, extra_mask=None) -> float:
    x, y, _ = _finite_pair(a, b, extra_mask)
    if x.size < 2:
        return float("nan")
    x = x - x.mean()
    y = y - y.mean()
    den = math.sqrt(float(np.sum(x*x) * np.sum(y*y)))
    if den <= 0:
        return float("nan")
    return float(np.sum(x*y) / den)


def _average_ranks(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=np.float64)
    order = np.argsort(x, kind="mergesort")
    ranks = np.empty(x.size, dtype=np.float64)
    i = 0
    while i < x.size:
        j = i + 1
        while j < x.size and x[order[j]] == x[order[i]]:
            j += 1
        # zero-based average rank is sufficient for correlation
        r = 0.5 * (i + j - 1)
        ranks[order[i:j]] = r
        i = j
    return ranks


def _spearman(a, b, extra_mask=None) -> float:
    x, y, _ = _finite_pair(a, b, extra_mask)
    if x.size < 2:
        return float("nan")
    return _pearson(_average_ranks(x), _average_ranks(y))


def _corr(prefix: str, a, b, mask=None) -> dict:
    _, _, common = _finite_pair(a, b, mask)
    return {
        f"{prefix}_finite_count": int(np.count_nonzero(common)),
        f"{prefix}_pearson": _pearson(a, b, mask),
        f"{prefix}_spearman": _spearman(a, b, mask),
    }


def _los_divergence_half(Rx: np.ndarray, Ry: np.ndarray, extent: float) -> np.ndarray:
    Rx = np.asarray(Rx, dtype=np.float64)
    Ry = np.asarray(Ry, dtype=np.float64)
    if Rx.shape != Ry.shape or Rx.ndim != 2:
        raise ValueError("LOS components must be same-shape 2D arrays")
    ny, nx = Rx.shape
    dx = (2.0 * extent) / max(nx - 1, 1)
    dy = (2.0 * extent) / max(ny - 1, 1)
    dRx_dx = np.gradient(Rx, dx, axis=1, edge_order=1)
    dRy_dy = np.gradient(Ry, dy, axis=0, edge_order=1)
    return 0.5 * (dRx_dx + dRy_dy)


def _binned_ray_displacement(photons: dict, bins: int, extent: float):
    x0 = np.asarray(photons["x0"], dtype=np.float64)
    y0 = np.asarray(photons["y0"], dtype=np.float64)
    xf = np.asarray(photons["x"], dtype=np.float64)
    yf = np.asarray(photons["y"], dtype=np.float64)
    edges = np.linspace(-extent, extent, bins + 1)

    counts, _, _ = np.histogram2d(y0, x0, bins=(edges, edges))
    sum_dx, _, _ = np.histogram2d(y0, x0, bins=(edges, edges), weights=xf-x0)
    sum_dy, _, _ = np.histogram2d(y0, x0, bins=(edges, edges), weights=yf-y0)

    dx_map = np.full((bins, bins), np.nan, dtype=np.float64)
    dy_map = np.full((bins, bins), np.nan, dtype=np.float64)
    good = counts > 0
    dx_map[good] = sum_dx[good] / counts[good]
    dy_map[good] = sum_dy[good] / counts[good]
    return counts, dx_map, dy_map, edges


def _masked_divergence_half(dx_map: np.ndarray, dy_map: np.ndarray, spacing: float):
    """Central-difference half-divergence only where all four neighbours exist."""
    dx_map = np.asarray(dx_map, dtype=np.float64)
    dy_map = np.asarray(dy_map, dtype=np.float64)
    if dx_map.shape != dy_map.shape or dx_map.ndim != 2:
        raise ValueError("deflection components must be same-shape 2D arrays")
    out = np.full(dx_map.shape, np.nan, dtype=np.float64)
    if dx_map.shape[0] < 3 or dx_map.shape[1] < 3:
        return out

    centre = (
        np.isfinite(dx_map[1:-1, 1:-1]) & np.isfinite(dy_map[1:-1, 1:-1]) &
        np.isfinite(dx_map[1:-1, :-2]) & np.isfinite(dx_map[1:-1, 2:]) &
        np.isfinite(dy_map[:-2, 1:-1]) & np.isfinite(dy_map[2:, 1:-1])
    )
    ddx = (dx_map[1:-1, 2:] - dx_map[1:-1, :-2]) / (2.0 * spacing)
    ddy = (dy_map[2:, 1:-1] - dy_map[:-2, 1:-1]) / (2.0 * spacing)
    val = 0.5 * (ddx + ddy)
    block = out[1:-1, 1:-1]
    block[centre] = val[centre]
    return out


def _coverage_bbox(mask: np.ndarray, extent: float) -> dict:
    mask = np.asarray(mask, dtype=bool)
    yy, xx = np.where(mask)
    if xx.size == 0:
        return {
            "coverage_row_min": None, "coverage_row_max": None,
            "coverage_col_min": None, "coverage_col_max": None,
            "coverage_x_min": None, "coverage_x_max": None,
            "coverage_y_min": None, "coverage_y_max": None,
        }
    ny, nx = mask.shape
    xcentres = -extent + (np.arange(nx) + 0.5) * (2.0*extent/nx)
    ycentres = -extent + (np.arange(ny) + 0.5) * (2.0*extent/ny)
    return {
        "coverage_row_min": int(yy.min()), "coverage_row_max": int(yy.max()),
        "coverage_col_min": int(xx.min()), "coverage_col_max": int(xx.max()),
        "coverage_x_min": float(xcentres[xx.min()]),
        "coverage_x_max": float(xcentres[xx.max()]),
        "coverage_y_min": float(ycentres[yy.min()]),
        "coverage_y_max": float(ycentres[yy.max()]),
    }


def _core_overlap(observed: np.ndarray, support_mask: np.ndarray) -> dict:
    observed = np.asarray(observed, dtype=np.float64)
    support = np.asarray(support_mask, dtype=bool)
    positive = observed[np.isfinite(observed) & (observed > 0)]
    out = {}
    for label, q in (("top25", 0.75), ("top10", 0.90)):
        if positive.size == 0:
            out[f"{label}_threshold_positive_kappa"] = float("nan")
            out[f"{label}_core_pixel_count"] = 0
            out[f"{label}_core_coverage_fraction"] = float("nan")
            out[f"support_fraction_inside_{label}_core"] = float("nan")
            continue
        threshold = float(np.quantile(positive, q))
        core = np.isfinite(observed) & (observed >= threshold) & (observed > 0)
        overlap = core & support
        ncore = int(np.count_nonzero(core))
        nsupport = int(np.count_nonzero(support))
        out[f"{label}_threshold_positive_kappa"] = threshold
        out[f"{label}_core_pixel_count"] = ncore
        out[f"{label}_core_coverage_fraction"] = (
            float(np.count_nonzero(overlap) / ncore) if ncore else float("nan")
        )
        out[f"support_fraction_inside_{label}_core"] = (
            float(np.count_nonzero(overlap) / nsupport) if nsupport else float("nan")
        )
    return out


def _run_cluster(cluster: dict) -> dict:
    cid = cluster["id"]
    real = PRIMARY._load_cluster(cluster)
    initial = PRIMARY._initial_state(real["rho3"])
    state = PRIMARY._evolve(initial)
    candidate = PRIMARY._candidate(state)
    vector = PRIMARY._interface_vector(candidate)
    ray = PRIMARY._ray_and_observable(cid, vector, real)

    observed = np.asarray(ray["observed_kappa"], dtype=np.float64)
    kappa = np.asarray(ray["observable"]["kappa"], dtype=np.float64)
    Rx = np.asarray(ray["los"]["comp_1"], dtype=np.float64)
    Ry = np.asarray(ray["los"]["comp_2"], dtype=np.float64)
    los_div = _los_divergence_half(Rx, Ry, PRIMARY.CFG["extent"])
    los_mag = np.hypot(Rx, Ry)

    counts, ray_dx, ray_dy, edges = _binned_ray_displacement(
        ray["photons"], PRIMARY.CFG["bins"], PRIMARY.CFG["extent"]
    )
    spacing = float(edges[1] - edges[0])
    ray_div = _masked_divergence_half(ray_dx, ray_dy, spacing)

    source_support = counts > 0
    jac_support = counts >= 6
    jac_finite = np.isfinite(kappa)
    support_matches_jacobian = bool(np.array_equal(jac_support, jac_finite))

    metrics = {
        "cluster_id": cid,
        "cluster_label": cluster["label"],
        "candidate_id": PRIMARY.CANDIDATE_ID,
        "physical_source_representation": PRIMARY.PHYSICAL_SOURCE,
        "source_launch": "launch_B_cartesian_frozen",
        "bins": int(PRIMARY.CFG["bins"]),
        "total_pixel_count": int(kappa.size),
        "source_support_pixel_count": int(np.count_nonzero(source_support)),
        "source_support_fraction": float(np.mean(source_support)),
        "jacobian_ge6_support_pixel_count": int(np.count_nonzero(jac_support)),
        "jacobian_ge6_support_fraction": float(np.mean(jac_support)),
        "jacobian_finite_pixel_count": int(np.count_nonzero(jac_finite)),
        "jacobian_finite_fraction": float(np.mean(jac_finite)),
        "jacobian_support_mask_equals_ge6_source_mask": support_matches_jacobian,
        "ray_divergence_finite_pixel_count": int(np.count_nonzero(np.isfinite(ray_div))),
        "source_count_min_nonzero": float(np.min(counts[counts > 0])) if np.any(counts > 0) else float("nan"),
        "source_count_max": float(np.max(counts)),
        "source_count_mean_nonzero": float(np.mean(counts[counts > 0])) if np.any(counts > 0) else float("nan"),
        **_coverage_bbox(jac_support, PRIMARY.CFG["extent"]),
        **_core_overlap(observed, jac_support),
    }

    # Morphology at the M10 -> LOS checkpoint.  Report both signed and sign-
    # reversed divergence explicitly; neither is selected as "best".
    metrics.update(_corr("los_div_vs_observed_full", los_div, observed))
    metrics.update(_corr("neg_los_div_vs_observed_full", -los_div, observed))
    metrics.update(_corr("los_mag_vs_observed_full", los_mag, observed))
    metrics.update(_corr("los_div_vs_observed_jacmask", los_div, observed, jac_finite))
    metrics.update(_corr("neg_los_div_vs_observed_jacmask", -los_div, observed, jac_finite))
    metrics.update(_corr("los_mag_vs_observed_jacmask", los_mag, observed, jac_finite))

    # Ray-displacement checkpoint and final Jacobian observable.
    metrics.update(_corr("ray_div_vs_observed", ray_div, observed))
    metrics.update(_corr("neg_ray_div_vs_observed", -ray_div, observed))
    metrics.update(_corr("jacobian_kappa_vs_observed", kappa, observed))

    # Stage-to-stage retention diagnostics.
    metrics.update(_corr("los_div_vs_ray_div", los_div, ray_div))
    metrics.update(_corr("los_div_vs_jacobian_kappa", los_div, kappa))
    metrics.update(_corr("ray_div_vs_jacobian_kappa", ray_div, kappa))

    # Observed field statistics on full map vs the exact Jacobian footprint.
    obs_full = observed[np.isfinite(observed)]
    obs_masked = observed[np.isfinite(observed) & jac_finite]
    metrics.update({
        "observed_full_finite_count": int(obs_full.size),
        "observed_full_mean": float(np.mean(obs_full)) if obs_full.size else float("nan"),
        "observed_full_rms": PRIMARY._rms(obs_full) if obs_full.size else float("nan"),
        "observed_jacmask_finite_count": int(obs_masked.size),
        "observed_jacmask_mean": float(np.mean(obs_masked)) if obs_masked.size else float("nan"),
        "observed_jacmask_rms": PRIMARY._rms(obs_masked) if obs_masked.size else float("nan"),
    })

    cdir = OUT / "clusters" / cid
    cdir.mkdir(parents=True, exist_ok=True)
    _write_json(cdir / "information_loss_metrics.json", metrics)
    np.savez_compressed(
        cdir / "checkpoint_fields.npz",
        observed_kappa=observed,
        los_rx=Rx,
        los_ry=Ry,
        los_divergence_half=los_div,
        los_magnitude=los_mag,
        source_counts=counts,
        source_support_mask=source_support,
        jacobian_ge6_support_mask=jac_support,
        jacobian_finite_mask=jac_finite,
        ray_deflection_x=ray_dx,
        ray_deflection_y=ray_dy,
        ray_deflection_divergence_half=ray_div,
        jacobian_kappa=kappa,
    )
    return metrics


def main() -> int:
    started = time.perf_counter()
    OUT.mkdir(parents=True, exist_ok=True)
    repo = _repo_state()
    _write_json(OUT / "repository_state.json", repo)

    if repo["branch"] != "main" or repo["tracked_changes"] or repo["staged_changes"]:
        validation = {
            "lab_id": LAB_ID,
            "outcome": "REPOSITORY_GATE_FAILURE",
            "head_sha": repo["head_sha"],
            "science_change_authorized": False,
        }
        _write_json(OUT / "validation.json", validation)
        print(json.dumps(validation, indent=2))
        return 2

    rows = []
    failures = []
    for cluster in CLUSTERS:
        try:
            print(f"[audit] {cluster['id']}")
            rows.append(_run_cluster(cluster))
        except Exception as exc:
            failures.append({
                "cluster_id": cluster.get("id"),
                "exception_type": type(exc).__name__,
                "message": str(exc),
            })
            # Fail closed: surface the runtime problem rather than silently
            # continuing with a partial scientific sample.
            _write_json(OUT / "cluster_failures.json", failures)
            raise

    _write_csv(OUT / "information_loss_summary.csv", rows)
    _write_json(OUT / "cluster_failures.json", failures)

    coverage_same = bool(rows) and len({r["jacobian_finite_pixel_count"] for r in rows}) == 1
    all_support_matches = bool(rows) and all(
        r["jacobian_support_mask_equals_ge6_source_mask"] for r in rows
    )
    complete = len(rows) == EXPECTED_CLUSTERS and not failures

    validation = {
        "lab_id": LAB_ID,
        "outcome": (
            "Outcome A — M10-TO-JACOBIAN INFORMATION-LOSS AUDIT COMPLETE"
            if complete else "Outcome D — INFORMATION-LOSS AUDIT INCOMPLETE"
        ),
        "head_sha": repo["head_sha"],
        "candidate_id": PRIMARY.CANDIDATE_ID,
        "physical_source_representation": PRIMARY.PHYSICAL_SOURCE,
        "cluster_count_expected": EXPECTED_CLUSTERS,
        "cluster_count_completed": len(rows),
        "all_jacobian_masks_match_ge6_source_support": all_support_matches,
        "jacobian_finite_count_identical_across_clusters": coverage_same,
        "jacobian_finite_pixel_count": rows[0]["jacobian_finite_pixel_count"] if coverage_same and rows else None,
        "jacobian_finite_fraction": rows[0]["jacobian_finite_fraction"] if coverage_same and rows else None,
        "science_interpretation_required": True,
        "physics_change_authorized": False,
        "candidate_change_authorized": False,
        "source_plane_change_authorized": False,
        "duration_seconds": time.perf_counter() - started,
    }
    _write_json(OUT / "validation.json", validation)
    _write_json(OUT / "run.json", {
        "lab_id": LAB_ID,
        "head_sha": repo["head_sha"],
        "candidate_id": PRIMARY.CANDIDATE_ID,
        "physical_source_representation": PRIMARY.PHYSICAL_SOURCE,
        "fixed_primary_config": PRIMARY.CFG,
        "duration_seconds": validation["duration_seconds"],
    })

    report = [
        f"# {LAB_ID}", "",
        f"**Head:** `{repo['head_sha']}`", "",
        f"**Outcome:** {validation['outcome']}", "",
        "This is a diagnostic audit. Correlations are measurements, not gates.", "",
        "| Cluster | Jacobian finite | fraction | mask==source>=6 | LOS div r full | LOS div r jacmask | ray div r | Jacobian kappa r | LOS->ray r | LOS->kappa r | ray->kappa r |",
        "|---|---:|---:|---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for r in rows:
        report.append(
            f"| {r['cluster_id']} | {r['jacobian_finite_pixel_count']} | "
            f"{r['jacobian_finite_fraction']:.6f} | "
            f"{r['jacobian_support_mask_equals_ge6_source_mask']} | "
            f"{r['los_div_vs_observed_full_pearson']:.6f} | "
            f"{r['los_div_vs_observed_jacmask_pearson']:.6f} | "
            f"{r['ray_div_vs_observed_pearson']:.6f} | "
            f"{r['jacobian_kappa_vs_observed_pearson']:.6f} | "
            f"{r['los_div_vs_ray_div_pearson']:.6f} | "
            f"{r['los_div_vs_jacobian_kappa_pearson']:.6f} | "
            f"{r['ray_div_vs_jacobian_kappa_pearson']:.6f} |"
        )
    report += [
        "", "No physics/source-plane/candidate change is authorized by this lab.",
    ]
    (OUT / "report.md").write_text("\n".join(report) + "\n")

    print(json.dumps(validation, indent=2))
    return 0 if complete else 1


if __name__ == "__main__":
    raise SystemExit(main())
