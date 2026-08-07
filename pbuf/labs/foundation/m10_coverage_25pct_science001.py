#!/usr/bin/env python3
"""PBUF FOUNDATION — M10 COVERAGE 25 PERCENT SCIENCE 001.

Paired five-cluster coverage experiment after M10 information-loss audit 001.

The physical candidate and propagation law are frozen.  For each real cluster
we construct the PL1_PM1_PS2 M10 interface field once, then run two observation
lanes through the same M14/M15/photon/Jacobian path:

  C07 control  : historical Launch-B geometry (~7.03% of 64x64 bins)
  C25 expanded : left-edge anchored 8 x 8 source rectangle, targeting
                 32 x 32 = 1024 Jacobian bins (25% of 64x64)

C25 preserves approximately the C07 photons-per-supported-bin density.  The
only intended experimental change is source-plane spatial coverage plus the
photon count required to retain that density.  No physics coefficient,
candidate, M10 field, propagation step, propagation length, Jacobian rule,
observable definition, fitting, optimisation, or amplitude scaling is changed.
"""
from __future__ import annotations

import csv
import hashlib
import json
import math
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
from astropy.io import fits

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))

from a8_three_dimensional_projection_lab001 import (
    CLUSTERS, PRODUCTION, construct_common_proxy, construct_rho_3d,
)
from weak_lensing_observation001 import propagate as wl_propagate, resample_to_grid
import observable_lab001 as obs_lab
import source_plane_lab001 as src_lab

from pbuf.core import pair_enumeration as M05
from pbuf.core import pair_transfer as M08
from pbuf.core import midpoint_rasterization as M10
from pbuf.core import los_projection as M14
from pbuf.core import ray_interface as M15
from pbuf.core import observable_extraction as M16
from pbuf.models import a8_state as M06_state
from pbuf.models import a8_pair_amplitude as M06
from pbuf.models import transverse_projector as M07

LAB_ID = "PBUF-FOUNDATION-M10-COVERAGE-25PCT-SCIENCE-001"
OUT = ROOT / "runs" / "m10_coverage_25pct_science001"
BENCHMARK = ROOT / "PBUF_benchmark"

CANDIDATE_ID = "PL1_PM1_PS2"
PHYSICAL_SOURCE = "M10_interface_field"
ENDPOINT_ROLE = "conservation_bookkeeping_only"
NZ = 9
PROFILE = "gaussian"
STRENGTH = 0.18
SEED = 12345
CFG = dict(PRODUCTION)

OBS_BINS = int(CFG["bins"])
TOTAL_PIXELS = OBS_BINS * OBS_BINS
CONTROL_NPHOTONS = int(CFG["nphotons"])
CONTROL_EXPECTED_PIXELS = 288
TARGET_PIXELS = 1024
TARGET_FRACTION = TARGET_PIXELS / TOTAL_PIXELS

# Historical source rectangle is 3 x 6 = 18 square units.  The expanded
# rectangle is 8 x 8 = 64 square units.  Scale photon count by area so the
# average ray density per supported observable bin stays approximately fixed.
CONTROL_SOURCE_AREA = 18.0
EXPANDED_X_MIN = -float(CFG["extent"])
EXPANDED_X_MAX = 0.0
EXPANDED_Y_MIN = -4.0
EXPANDED_Y_MAX = 4.0
EXPANDED_SOURCE_AREA = (EXPANDED_X_MAX - EXPANDED_X_MIN) * (EXPANDED_Y_MAX - EXPANDED_Y_MIN)
EXPANDED_SIDE = int(round(math.sqrt(CONTROL_NPHOTONS * EXPANDED_SOURCE_AREA / CONTROL_SOURCE_AREA)))
EXPANDED_NPHOTONS = EXPANDED_SIDE * EXPANDED_SIDE


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
    if isinstance(obj, np.floating): return float(obj)
    if isinstance(obj, np.integer): return int(obj)
    if isinstance(obj, np.bool_): return bool(obj)
    if isinstance(obj, Path): return str(obj)
    if isinstance(obj, tuple): return list(obj)
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
        w.writeheader(); w.writerows(rows)


def _sha_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _sha_arr(a) -> str:
    return hashlib.sha256(np.ascontiguousarray(np.asarray(a, dtype=np.float64)).tobytes()).hexdigest()


def _vec_sha(v) -> str:
    h = hashlib.sha256()
    for a in v:
        h.update(np.ascontiguousarray(np.asarray(a, dtype=np.float64)).tobytes())
    return h.hexdigest()


def _rms(a) -> float:
    x = np.asarray(a, dtype=np.float64)
    return float(np.sqrt(np.mean(x*x))) if x.size else float("nan")


def _corr(a, b) -> tuple[float, float, int]:
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    mask = np.isfinite(a) & np.isfinite(b)
    n = int(mask.sum())
    if n < 2:
        return float("nan"), float("nan"), n
    return (
        float(M16.pearson(a[mask], b[mask])),
        float(M16.spearman(a[mask], b[mask])),
        n,
    )


def _load_cluster(cluster: dict) -> dict:
    path = BENCHMARK / cluster["directory"] / f"hlsp_frontier_model_{cluster['slug']}_merten_v1_kappa.fits"
    if not path.exists():
        raise FileNotFoundError(path)
    with fits.open(path) as hdul:
        kappa = np.asarray(hdul[0].data, dtype=np.float64)
    rho2 = construct_common_proxy(kappa, bins=OBS_BINS, extent=CFG["extent"])
    rho3 = construct_rho_3d(rho2, NZ, profile=PROFILE)
    return {
        "path": path,
        "kappa": kappa,
        "rho2": rho2,
        "rho3": rho3,
        "observed_kappa": resample_to_grid(kappa, OBS_BINS, CFG["extent"]),
    }


def _initial_state(rho3: np.ndarray) -> dict:
    rng = np.random.RandomState(SEED)
    eq = STRENGTH * rho3
    noise = M06_state.A8_INIT_INJECTION_NOISE * STRENGTH * rng.randn(*rho3.shape)
    return {"rho_3d": rho3.copy(), "u_slow0": eq.copy(), "u_fast0": eq + noise}


def _evolve(initial: dict) -> dict:
    us, uf, history = M06_state.evolve_a8_transport_3d(
        initial["u_slow0"].copy(), initial["u_fast0"].copy(),
        stencil="N6", boundary="reflective",
    )
    return {"rho_3d": initial["rho_3d"].copy(), "u_slow": us, "u_fast": uf, "c_state": history[-1]}


def _candidate(state: dict) -> dict:
    shape = tuple(state["c_state"].shape)
    pairs = M05.enumerate_internal_pairs(shape)
    ex, ey, ez, valid, gmag = M07.build_longitudinal_direction(state["c_state"])
    projector = M07.build_transverse_projector(ex, ey, ez)
    amp = M06.compute_a8_pair_amplitudes(state["u_slow"], state["u_fast"], state["c_state"], pairs)
    response = M08.build_pair_responses(pairs, amp, projector, magnitude_formulation="PM1", pair_symmetrization="PS2")
    endpoint = M08.assemble_endpoint_field(response, shape)
    interface = M10.rasterize_interface_field(response, shape)
    return {
        "shape": shape, "pairs": pairs, "valid_count": int(np.count_nonzero(valid)),
        "gradient_rms": _rms(gmag), "endpoint": endpoint, "interface": interface,
    }


def _interface_vector(candidate: dict):
    x = candidate["interface"]
    return x["Rx_3d_interface"], x["Ry_3d_interface"], x["Rz_3d_interface"]


def _launch_expanded_25pct():
    """Deterministic 266x266 Cartesian launch over an 8x8 left-edge rectangle."""
    x_edges = np.linspace(EXPANDED_X_MIN, EXPANDED_X_MAX, EXPANDED_SIDE + 1)
    y_edges = np.linspace(EXPANDED_Y_MIN, EXPANDED_Y_MAX, EXPANDED_SIDE + 1)
    xc = 0.5 * (x_edges[:-1] + x_edges[1:])
    yc = 0.5 * (y_edges[:-1] + y_edges[1:])
    X, Y = np.meshgrid(xc, yc, indexing="xy")
    x0 = X.ravel(); y0 = Y.ravel()
    vx0 = np.ones_like(x0); vy0 = np.zeros_like(x0)
    return x0, y0, vx0, vy0


def _source_counts(x0, y0) -> tuple[np.ndarray, np.ndarray]:
    edges = np.linspace(-CFG["extent"], CFG["extent"], OBS_BINS + 1)
    counts, _, _ = np.histogram2d(y0, x0, bins=(edges, edges))
    return counts, edges


def _coverage_metrics(counts: np.ndarray, edges: np.ndarray) -> dict:
    support = counts >= 6
    rr, cc = np.where(support)
    if rr.size:
        row_min, row_max = int(rr.min()), int(rr.max())
        col_min, col_max = int(cc.min()), int(cc.max())
        x_min, x_max = float(edges[col_min]), float(edges[col_max+1])
        y_min, y_max = float(edges[row_min]), float(edges[row_max+1])
    else:
        row_min=row_max=col_min=col_max=None
        x_min=x_max=y_min=y_max=None
    nz = counts[counts > 0]
    return {
        "support_pixel_count": int(support.sum()),
        "support_fraction": float(support.mean()),
        "source_count_min_nonzero": float(nz.min()) if nz.size else 0.0,
        "source_count_max": float(nz.max()) if nz.size else 0.0,
        "source_count_mean_nonzero": float(nz.mean()) if nz.size else 0.0,
        "coverage_row_min": row_min, "coverage_row_max": row_max,
        "coverage_col_min": col_min, "coverage_col_max": col_max,
        "coverage_x_min": x_min, "coverage_x_max": x_max,
        "coverage_y_min": y_min, "coverage_y_max": y_max,
    }


def _run_lane(cluster_id: str, lane: str, vector, observed_kappa: np.ndarray) -> dict:
    los = M14.project_vector_to_image_plane(*vector, los_axis="z")
    Rx, Ry = los["comp_1"], los["comp_2"]
    if lane == "C07":
        x0, y0, vx0, vy0 = src_lab.launch_B_cartesian(CONTROL_NPHOTONS)
    elif lane == "C25":
        x0, y0, vx0, vy0 = _launch_expanded_25pct()
    else:
        raise ValueError(lane)

    counts, edges = _source_counts(x0, y0)
    coverage = _coverage_metrics(counts, edges)
    metadata = {
        "candidate_id": CANDIDATE_ID,
        "cluster_id": cluster_id,
        "transform_id": "RC0",
        "role": "los",
        "physical_source_representation": PHYSICAL_SOURCE,
        "coverage_lane": lane,
        "source_artifact_ids": [f"{cluster_id}_same_M10_interface_field"],
    }
    artifact = M15.prepare_ray_input(Rx, Ry, metadata, require_nontrivial=True)
    grid = np.linspace(-CFG["extent"], CFG["extent"], Rx.shape[0])
    field = {"xgrid": grid, "ygrid": grid, "rx": Rx, "ry": Ry}
    photons = wl_propagate(field, CFG["step"], CFG["steps"], x0, y0, vx0, vy0)
    photons["x0"] = x0; photons["y0"] = y0
    jac = obs_lab.method_jacobian(x0, y0, photons["x"], photons["y"], CFG["extent"], OBS_BINS)
    observable = M16.package_lensing_observables(
        jac["convergence"], jac["shear_g1"], jac["shear_g2"], reference_kappa=observed_kappa,
    )

    kappa = np.asarray(observable["kappa"], dtype=np.float64)
    finite = np.isfinite(kappa)
    displacement = np.hypot(photons["x"]-x0, photons["y"]-y0)
    traj = hashlib.sha256()
    for key in ("xs", "ys", "x", "y", "conservation"):
        traj.update(np.ascontiguousarray(np.asarray(photons[key], dtype=np.float64)).tobytes())

    # Full M10 LOS magnitude is unchanged between lanes; correlation on the
    # lane's Jacobian footprint shows how representative that sampled region is.
    los_mag = np.hypot(Rx, Ry)
    p_mag, s_mag, n_mag = _corr(los_mag[finite], observed_kappa[finite])
    metrics = {
        "lane": lane,
        "n_photons": int(len(x0)),
        **coverage,
        "jacobian_finite_pixel_count": int(finite.sum()),
        "jacobian_finite_fraction": float(finite.mean()),
        "jacobian_mask_equals_ge6_support": bool(np.array_equal(finite, counts >= 6)),
        "ray_classification": artifact.statistics["ray_classification"],
        "mean_endpoint_displacement": float(np.mean(displacement)),
        "max_endpoint_displacement": float(np.max(displacement)),
        "conservation_max": float(np.max(photons["conservation"])),
        "kappa_variance_finite": float(np.var(kappa[finite])) if finite.sum() >= 2 else float("nan"),
        "kappa_rms_finite": _rms(kappa[finite]) if finite.any() else float("nan"),
        "pearson_vs_observed": observable.get("pearson_vs_reference", float("nan")),
        "spearman_vs_observed": observable.get("spearman_vs_reference", float("nan")),
        "los_mag_vs_observed_on_jacmask_pearson": p_mag,
        "los_mag_vs_observed_on_jacmask_spearman": s_mag,
        "los_mag_vs_observed_on_jacmask_count": n_mag,
        "trajectory_sha256": traj.hexdigest(),
    }
    return {
        "metrics": metrics, "kappa": kappa, "finite": finite,
        "source_counts": counts, "los_Rx": Rx, "los_Ry": Ry,
    }


def _core_overlap(observed: np.ndarray, mask: np.ndarray) -> dict:
    pos = observed[np.isfinite(observed) & (observed > 0)]
    out = {}
    for q, label in ((75.0, "top25"), (90.0, "top10")):
        thr = float(np.percentile(pos, q)) if pos.size else float("nan")
        core = np.isfinite(observed) & (observed >= thr)
        ncore = int(core.sum())
        out[f"{label}_threshold_positive_kappa"] = thr
        out[f"{label}_core_pixel_count"] = ncore
        out[f"{label}_core_coverage_fraction"] = float((core & mask).sum()/ncore) if ncore else float("nan")
        out[f"support_fraction_inside_{label}_core"] = float((core & mask).sum()/max(int(mask.sum()),1))
    return out


def _run_cluster(cluster: dict) -> dict:
    cid = cluster["id"]
    real = _load_cluster(cluster)
    state = _evolve(_initial_state(real["rho3"]))
    candidate = _candidate(state)
    vector = _interface_vector(candidate)
    expected_pairs = int(M08.expected_interface_pair_count(candidate["shape"]))
    consumed_pairs = int(candidate["interface"]["statistics"]["consumed_pair_count_total"])

    c07 = _run_lane(cid, "C07", vector, real["observed_kappa"])
    c25 = _run_lane(cid, "C25", vector, real["observed_kappa"])
    common = c07["finite"] & c25["finite"]
    p_cross, s_cross, n_cross = _corr(c07["kappa"][common], c25["kappa"][common])

    m07, m25 = c07["metrics"], c25["metrics"]
    row = {
        "cluster_id": cid,
        "cluster_label": cluster["label"],
        "candidate_id": CANDIDATE_ID,
        "physical_source_representation": PHYSICAL_SOURCE,
        "m10_interface_sha256": _vec_sha(vector),
        "n_pairs": int(len(candidate["pairs"])),
        "interface_expected_pair_count": expected_pairs,
        "interface_consumed_pair_count": consumed_pairs,
        "interface_pair_count_ok": bool(expected_pairs == consumed_pairs == len(candidate["pairs"])),
        "interface_energy": float(candidate["interface"]["statistics"]["interface_energy"]),
        "endpoint_closure_bookkeeping": float(candidate["endpoint"]["statistics"]["global_vector_sum_norm"]),
        "C07_n_photons": m07["n_photons"],
        "C07_jacobian_finite_pixel_count": m07["jacobian_finite_pixel_count"],
        "C07_jacobian_finite_fraction": m07["jacobian_finite_fraction"],
        "C07_pearson_vs_observed": m07["pearson_vs_observed"],
        "C07_spearman_vs_observed": m07["spearman_vs_observed"],
        "C07_source_count_mean_nonzero": m07["source_count_mean_nonzero"],
        "C25_n_photons": m25["n_photons"],
        "C25_jacobian_finite_pixel_count": m25["jacobian_finite_pixel_count"],
        "C25_jacobian_finite_fraction": m25["jacobian_finite_fraction"],
        "C25_coverage_target_1024_met": bool(m25["jacobian_finite_pixel_count"] == TARGET_PIXELS),
        "C25_jacobian_mask_equals_ge6_support": m25["jacobian_mask_equals_ge6_support"],
        "C25_source_count_min_nonzero": m25["source_count_min_nonzero"],
        "C25_source_count_max": m25["source_count_max"],
        "C25_source_count_mean_nonzero": m25["source_count_mean_nonzero"],
        "C25_coverage_row_min": m25["coverage_row_min"],
        "C25_coverage_row_max": m25["coverage_row_max"],
        "C25_coverage_col_min": m25["coverage_col_min"],
        "C25_coverage_col_max": m25["coverage_col_max"],
        "C25_coverage_x_min": m25["coverage_x_min"],
        "C25_coverage_x_max": m25["coverage_x_max"],
        "C25_coverage_y_min": m25["coverage_y_min"],
        "C25_coverage_y_max": m25["coverage_y_max"],
        "C25_pearson_vs_observed": m25["pearson_vs_observed"],
        "C25_spearman_vs_observed": m25["spearman_vs_observed"],
        "delta_pearson_C25_minus_C07": float(m25["pearson_vs_observed"] - m07["pearson_vs_observed"]),
        "delta_spearman_C25_minus_C07": float(m25["spearman_vs_observed"] - m07["spearman_vs_observed"]),
        "C25_los_mag_vs_observed_on_jacmask_pearson": m25["los_mag_vs_observed_on_jacmask_pearson"],
        "C25_los_mag_vs_observed_on_jacmask_spearman": m25["los_mag_vs_observed_on_jacmask_spearman"],
        "C07_vs_C25_kappa_commonmask_count": n_cross,
        "C07_vs_C25_kappa_commonmask_pearson": p_cross,
        "C07_vs_C25_kappa_commonmask_spearman": s_cross,
        "C25_ray_classification": m25["ray_classification"],
        "C25_mean_endpoint_displacement": m25["mean_endpoint_displacement"],
        "C25_max_endpoint_displacement": m25["max_endpoint_displacement"],
        "C25_conservation_max": m25["conservation_max"],
        "C25_kappa_variance_finite": m25["kappa_variance_finite"],
        "C25_kappa_rms_finite": m25["kappa_rms_finite"],
        **{f"C25_{k}": v for k,v in _core_overlap(real["observed_kappa"], c25["finite"]).items()},
    }

    cluster_dir = OUT / "clusters" / cid
    cluster_dir.mkdir(parents=True, exist_ok=True)
    _write_json(cluster_dir / "input_provenance.json", {
        "cluster_id": cid,
        "fits_path": str(real["path"].relative_to(ROOT)),
        "fits_sha256": _sha_file(real["path"]),
        "rho2_sha256": _sha_arr(real["rho2"]),
        "rho3_sha256": _sha_arr(real["rho3"]),
        "m10_interface_sha256": row["m10_interface_sha256"],
    })
    _write_json(cluster_dir / "control_07_metrics.json", m07)
    _write_json(cluster_dir / "expanded_25_metrics.json", m25)
    _write_json(cluster_dir / "paired_comparison.json", row)
    np.savez_compressed(
        cluster_dir / "paired_observables.npz",
        observed_kappa=real["observed_kappa"],
        C07_kappa=c07["kappa"], C07_mask=c07["finite"], C07_source_counts=c07["source_counts"],
        C25_kappa=c25["kappa"], C25_mask=c25["finite"], C25_source_counts=c25["source_counts"],
        los_Rx=c25["los_Rx"], los_Ry=c25["los_Ry"],
    )
    integrity = bool(
        row["interface_pair_count_ok"]
        and m07["jacobian_finite_pixel_count"] == CONTROL_EXPECTED_PIXELS
        and m07["jacobian_mask_equals_ge6_support"]
        and m25["jacobian_mask_equals_ge6_support"]
        and m25["jacobian_finite_pixel_count"] == TARGET_PIXELS
        and m25["ray_classification"] in ("structured_small", "structured_normal", "constant_nonzero")
        and np.isfinite(m25["conservation_max"])
    )
    row["integrity_pass"] = integrity
    return row


def main() -> int:
    started = time.perf_counter()
    OUT.mkdir(parents=True, exist_ok=True)
    repo = _repo_state(); _write_json(OUT / "repository_state.json", repo)
    if repo["branch"] != "main" or repo["tracked_changes"] or repo["staged_changes"]:
        v = {"lab_id": LAB_ID, "outcome": "REPOSITORY_GATE_FAILURE", "head_sha": repo["head_sha"]}
        _write_json(OUT / "validation.json", v); print(json.dumps(v, indent=2)); return 2

    _write_json(OUT / "run_config.json", {
        "lab_id": LAB_ID,
        "candidate_id": CANDIDATE_ID,
        "physical_source_representation": PHYSICAL_SOURCE,
        "endpoint_role": ENDPOINT_ROLE,
        "cluster_ids": [c["id"] for c in CLUSTERS],
        "nz": NZ, "profile": PROFILE, "strength": STRENGTH, "seed": SEED,
        "control_nphotons_requested": CONTROL_NPHOTONS,
        "control_expected_jacobian_pixels": CONTROL_EXPECTED_PIXELS,
        "expanded_source_rectangle": [EXPANDED_X_MIN, EXPANDED_X_MAX, EXPANDED_Y_MIN, EXPANDED_Y_MAX],
        "expanded_source_area": EXPANDED_SOURCE_AREA,
        "expanded_grid_side": EXPANDED_SIDE,
        "expanded_nphotons": EXPANDED_NPHOTONS,
        "expanded_target_pixels": TARGET_PIXELS,
        "expanded_target_fraction": TARGET_FRACTION,
        "density_rule": "scale photon count by source-rectangle area relative to C07",
        "propagation_step": CFG["step"], "propagation_steps": CFG["steps"],
    })

    rows=[]; failures=[]
    for cluster in CLUSTERS:
        print(f"[{cluster['id']}] same M10 field -> C07 control + C25 expanded coverage")
        try:
            rows.append(_run_cluster(cluster))
        except Exception as exc:
            failures.append({"cluster_id": cluster["id"], "error": repr(exc)})
            raise

    _write_csv(OUT / "coverage_25pct_summary.csv", rows)
    _write_json(OUT / "cluster_failures.json", failures)
    all_integrity = bool(len(rows)==len(CLUSTERS) and all(r["integrity_pass"] for r in rows))
    all_target = bool(rows and all(r["C25_coverage_target_1024_met"] for r in rows))
    p07 = [r["C07_pearson_vs_observed"] for r in rows if np.isfinite(r["C07_pearson_vs_observed"])]
    p25 = [r["C25_pearson_vs_observed"] for r in rows if np.isfinite(r["C25_pearson_vs_observed"])]
    s07 = [r["C07_spearman_vs_observed"] for r in rows if np.isfinite(r["C07_spearman_vs_observed"])]
    s25 = [r["C25_spearman_vs_observed"] for r in rows if np.isfinite(r["C25_spearman_vs_observed"])]

    outcome = (
        "Outcome A — 25-PERCENT COVERAGE TARGET ACHIEVED; PAIRED SCIENCE COMPLETE"
        if all_integrity and all_target else
        "Outcome B — 25-PERCENT COVERAGE TARGET/INTEGRITY FAILURE"
    )
    validation = {
        "lab_id": LAB_ID,
        "outcome": outcome,
        "head_sha": repo["head_sha"],
        "candidate_id": CANDIDATE_ID,
        "physical_source_representation": PHYSICAL_SOURCE,
        "cluster_count_expected": len(CLUSTERS),
        "cluster_count_completed": len(rows),
        "all_cluster_integrity_pass": all_integrity,
        "all_clusters_C25_target_1024_pixels": all_target,
        "C07_expected_pixel_count": CONTROL_EXPECTED_PIXELS,
        "C25_target_pixel_count": TARGET_PIXELS,
        "C25_target_fraction": TARGET_FRACTION,
        "C25_nphotons": EXPANDED_NPHOTONS,
        "C07_pearson_mean_finite": float(np.mean(p07)) if p07 else float("nan"),
        "C25_pearson_mean_finite": float(np.mean(p25)) if p25 else float("nan"),
        "C07_spearman_mean_finite": float(np.mean(s07)) if s07 else float("nan"),
        "C25_spearman_mean_finite": float(np.mean(s25)) if s25 else float("nan"),
        "mean_delta_pearson_C25_minus_C07": float(np.mean([r["delta_pearson_C25_minus_C07"] for r in rows])) if rows else float("nan"),
        "mean_delta_spearman_C25_minus_C07": float(np.mean([r["delta_spearman_C25_minus_C07"] for r in rows])) if rows else float("nan"),
        "science_interpretation_required": True,
        "physics_change_authorized": False,
        "candidate_change_authorized": False,
        "next_coverage_change_authorized": False,
        "duration_seconds": time.perf_counter()-started,
    }
    _write_json(OUT / "validation.json", validation)
    _write_json(OUT / "run.json", {"lab_id": LAB_ID, "head_sha": repo["head_sha"], "duration_seconds": validation["duration_seconds"]})

    report=[f"# {LAB_ID}","",f"**Outcome:** {outcome}","",f"C25 target: {TARGET_PIXELS}/{TOTAL_PIXELS} = {TARGET_FRACTION:.4f}",f"C25 photons: {EXPANDED_NPHOTONS} ({EXPANDED_SIDE}x{EXPANDED_SIDE})","", "| cluster | C07 px | C25 px | C07 r | C25 r | delta r | C07 rho | C25 rho | delta rho | integrity |", "|---|---:|---:|---:|---:|---:|---:|---:|---:|---|"]
    for r in rows:
        report.append(f"| {r['cluster_id']} | {r['C07_jacobian_finite_pixel_count']} | {r['C25_jacobian_finite_pixel_count']} | {r['C07_pearson_vs_observed']:.4f} | {r['C25_pearson_vs_observed']:.4f} | {r['delta_pearson_C25_minus_C07']:.4f} | {r['C07_spearman_vs_observed']:.4f} | {r['C25_spearman_vs_observed']:.4f} | {r['delta_spearman_C25_minus_C07']:.4f} | {r['integrity_pass']} |")
    report += ["", "Correlation is a science measurement, not an execution acceptance gate. No further coverage or physics change is authorized by this lab."]
    (OUT / "report.md").write_text("\n".join(report)+"\n")
    print(json.dumps(validation, indent=2))
    return 0 if all_integrity and all_target else 1


if __name__ == "__main__":
    raise SystemExit(main())
