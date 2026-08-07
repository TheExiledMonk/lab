#!/usr/bin/env python3
"""PBUF FOUNDATION — JACOBIAN OBSERVABLE CONSISTENCY AUDIT 001.

Diagnostic-only five-cluster audit after the C07/C25 coverage experiment.

The requalified PL1_PM1_PS2 / M10 field and the C25 source geometry are held
fixed.  This lab isolates the final observable mapping by reconstructing the
per-bin ray Jacobian independently from the same propagated rays and testing

    J = I + D
    det(J) = 1 + tr(D) + det(D)                    (exact in 2D)
    kappa_production = 1 - det(J)
                     = -tr(D) - det(D)             (exact identity)

It also reports the weak/linear trace lanes -tr(D), -0.5 tr(D), and their sign
reversals as diagnostics only.  No sign or formula is selected after seeing
correlations.  No physics/source-plane/coverage/candidate implementation is
changed by this lab.
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
sys.path.insert(0, str(Path(__file__).resolve().parent))

import m10_coverage_25pct_science001 as COV
from weak_lensing_observation001 import propagate as wl_propagate
import observable_lab001 as obs_lab
from pbuf.core import los_projection as M14
from pbuf.core import ray_interface as M15
from pbuf.core import observable_extraction as M16

LAB_ID = "PBUF-FOUNDATION-JACOBIAN-OBSERVABLE-CONSISTENCY-001"
OUT = ROOT / "runs" / "jacobian_observable_consistency001"
C25_NPHOTONS = int(COV.EXPANDED_NPHOTONS)
BINS = int(COV.OBS_BINS)
EXTENT = float(COV.CFG["extent"])
IDENTITY_TOL = 1e-10
LINEAR_RETENTION_THRESHOLD = 0.99


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


def _json_default(o):
    if isinstance(o, np.floating): return float(o)
    if isinstance(o, np.integer): return int(o)
    if isinstance(o, np.bool_): return bool(o)
    if isinstance(o, Path): return str(o)
    if isinstance(o, tuple): return list(o)
    return str(o)


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


def _rms(a) -> float:
    x = np.asarray(a, dtype=np.float64)
    x = x[np.isfinite(x)]
    if not x.size:
        return float("nan")
    return float(np.sqrt(np.mean(x*x)))


def _corr(a, b) -> tuple[float, float, int]:
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    if a.shape != b.shape:
        raise RuntimeError(f"correlation shape mismatch {a.shape} != {b.shape}")
    mask = np.isfinite(a) & np.isfinite(b)
    return (
        float(M16.safe_pearson(a, b)),
        float(M16.safe_spearman(a, b)),
        int(mask.sum()),
    )


def _relative_error(ref, test) -> float:
    a = np.asarray(ref, dtype=np.float64)
    b = np.asarray(test, dtype=np.float64)
    mask = np.isfinite(a) & np.isfinite(b)
    if not mask.any():
        return float("inf")
    den = float(np.linalg.norm(a[mask]))
    num = float(np.linalg.norm((b-a)[mask]))
    return num / max(den, 1e-15)


def _independent_jacobian_fit(x0, y0, xf, yf) -> dict:
    """Reconstruct per-bin J without calling observable_lab.method_jacobian."""
    edges = np.linspace(-EXTENT, EXTENT, BINS + 1)
    shape = (BINS, BINS)
    J00 = np.full(shape, np.nan)
    J01 = np.full(shape, np.nan)
    J10 = np.full(shape, np.nan)
    J11 = np.full(shape, np.nan)
    counts = np.zeros(shape, dtype=np.int64)

    for i in range(BINS):
        for j in range(BINS):
            mask = (
                (x0 >= edges[j]) & (x0 < edges[j+1]) &
                (y0 >= edges[i]) & (y0 < edges[i+1])
            )
            n = int(mask.sum())
            counts[i, j] = n
            if n < 6:
                continue
            xi = x0[mask]; yi = y0[mask]
            xo = xf[mask]; yo = yf[mask]
            A = np.column_stack([xi-xi.mean(), yi-yi.mean()])
            bx = xo-xo.mean(); by = yo-yo.mean()
            try:
                cx, *_ = np.linalg.lstsq(A, bx, rcond=None)
                cy, *_ = np.linalg.lstsq(A, by, rcond=None)
            except np.linalg.LinAlgError:
                continue
            J00[i,j] = cx[0]; J01[i,j] = cx[1]
            J10[i,j] = cy[0]; J11[i,j] = cy[1]

    D00 = J00 - 1.0
    D01 = J01
    D10 = J10
    D11 = J11 - 1.0
    trace_D = D00 + D11
    det_D = D00*D11 - D01*D10
    det_J = J00*J11 - J01*J10

    kappa_det = 1.0 - det_J
    kappa_exact_identity = -trace_D - det_D
    kappa_linear_trace = -trace_D
    kappa_standard_half_trace = -0.5 * trace_D
    kappa_positive_trace = trace_D
    kappa_positive_half_trace = 0.5 * trace_D

    return {
        "counts": counts,
        "J00": J00, "J01": J01, "J10": J10, "J11": J11,
        "D00": D00, "D01": D01, "D10": D10, "D11": D11,
        "trace_D": trace_D, "det_D": det_D, "det_J": det_J,
        "kappa_det": kappa_det,
        "kappa_exact_identity": kappa_exact_identity,
        "kappa_linear_trace": kappa_linear_trace,
        "kappa_standard_half_trace": kappa_standard_half_trace,
        "kappa_positive_trace": kappa_positive_trace,
        "kappa_positive_half_trace": kappa_positive_half_trace,
    }


def _affine_fixture() -> dict:
    """Analytic identity fixture independent of cluster/ray propagation."""
    cases = [
        ("identity", np.array([[1.0,0.0],[0.0,1.0]])),
        ("weak_diag", np.array([[1.01,0.0],[0.0,0.98]])),
        ("weak_shear", np.array([[1.005,0.02],[-0.01,0.995]])),
        ("moderate", np.array([[1.08,0.03],[0.02,0.94]])),
    ]
    rows = []
    passed = True
    for name, J in cases:
        D = J - np.eye(2)
        lhs = 1.0 - float(np.linalg.det(J))
        rhs = -float(np.trace(D)) - float(np.linalg.det(D))
        err = abs(lhs-rhs)
        rows.append({"fixture":name,"kappa_det":lhs,"kappa_identity":rhs,"abs_error":err})
        passed &= err <= 1e-14
    return {"rows": rows, "passes": bool(passed)}


def _run_cluster(cluster: dict) -> dict:
    cid = cluster["id"]
    real = COV._load_cluster(cluster)
    state = COV._evolve(COV._initial_state(real["rho3"]))
    candidate = COV._candidate(state)
    vector = COV._interface_vector(candidate)
    los = M14.project_vector_to_image_plane(*vector, los_axis="z")
    Rx = np.asarray(los["comp_1"], dtype=np.float64)
    Ry = np.asarray(los["comp_2"], dtype=np.float64)

    x0, y0, vx0, vy0 = COV._launch_expanded_25pct()
    if len(x0) != C25_NPHOTONS:
        raise RuntimeError(f"{cid}: C25 launch count {len(x0)} != {C25_NPHOTONS}")

    metadata = {
        "candidate_id": COV.CANDIDATE_ID,
        "cluster_id": cid,
        "transform_id": "RC0",
        "role": "los",
        "physical_source_representation": COV.PHYSICAL_SOURCE,
        "coverage_lane": "C25",
        "audit": LAB_ID,
        "source_artifact_ids": [f"{cid}_same_M10_interface_field"],
    }
    artifact = M15.prepare_ray_input(Rx, Ry, metadata, require_nontrivial=True)
    grid = np.linspace(-EXTENT, EXTENT, Rx.shape[0])
    field = {"xgrid": grid, "ygrid": grid, "rx": Rx, "ry": Ry}
    photons = wl_propagate(field, COV.CFG["step"], COV.CFG["steps"], x0, y0, vx0, vy0)

    prod = obs_lab.method_jacobian(x0, y0, photons["x"], photons["y"], EXTENT, BINS)
    k_prod = np.asarray(prod["convergence"], dtype=np.float64)
    fit = _independent_jacobian_fit(x0, y0, photons["x"], photons["y"])
    observed = np.asarray(real["observed_kappa"], dtype=np.float64)

    finite_prod = np.isfinite(k_prod)
    finite_fit = np.isfinite(fit["kappa_det"])
    mask_match = bool(np.array_equal(finite_prod, finite_fit))
    exact_identity_error = _relative_error(fit["kappa_det"], fit["kappa_exact_identity"])
    production_reconstruction_error = _relative_error(k_prod, fit["kappa_det"])

    p_prod_fit, s_prod_fit, n_prod_fit = _corr(k_prod, fit["kappa_det"])
    p_prod_linear, s_prod_linear, n_prod_linear = _corr(k_prod, fit["kappa_linear_trace"])
    p_prod_half, s_prod_half, n_prod_half = _corr(k_prod, fit["kappa_standard_half_trace"])

    p_obs_prod, s_obs_prod, n_obs_prod = _corr(k_prod, observed)
    p_obs_exact, s_obs_exact, n_obs_exact = _corr(fit["kappa_exact_identity"], observed)
    p_obs_linear, s_obs_linear, n_obs_linear = _corr(fit["kappa_linear_trace"], observed)
    p_obs_half, s_obs_half, n_obs_half = _corr(fit["kappa_standard_half_trace"], observed)
    p_obs_pos, s_obs_pos, n_obs_pos = _corr(fit["kappa_positive_trace"], observed)
    p_obs_poshalf, s_obs_poshalf, n_obs_poshalf = _corr(fit["kappa_positive_half_trace"], observed)

    trace_rms = _rms(fit["trace_D"])
    detD_rms = _rms(fit["det_D"])
    nonlinear_ratio = detD_rms / max(trace_rms, 1e-15)
    los_mag = np.hypot(Rx, Ry)
    p_los_obs, s_los_obs, n_los_obs = _corr(los_mag, observed)
    p_los_prod, s_los_prod, n_los_prod = _corr(los_mag, k_prod)

    stats = {
        "cluster_id": cid,
        "cluster_label": cluster["label"],
        "candidate_id": COV.CANDIDATE_ID,
        "physical_source_representation": COV.PHYSICAL_SOURCE,
        "coverage_lane": "C25",
        "n_photons": int(len(x0)),
        "finite_production_kappa_count": int(finite_prod.sum()),
        "finite_independent_fit_count": int(finite_fit.sum()),
        "finite_masks_match": mask_match,
        "production_vs_independent_relative_error": production_reconstruction_error,
        "exact_det_identity_relative_error": exact_identity_error,
        "production_vs_independent_pearson": p_prod_fit,
        "production_vs_independent_spearman": s_prod_fit,
        "production_vs_independent_count": n_prod_fit,
        "production_vs_linear_trace_pearson": p_prod_linear,
        "production_vs_linear_trace_spearman": s_prod_linear,
        "production_vs_linear_trace_count": n_prod_linear,
        "production_vs_half_trace_pearson": p_prod_half,
        "production_vs_half_trace_spearman": s_prod_half,
        "production_vs_half_trace_count": n_prod_half,
        "trace_D_rms": trace_rms,
        "det_D_rms": detD_rms,
        "nonlinear_detD_to_traceD_rms_ratio": nonlinear_ratio,
        "production_kappa_rms": _rms(k_prod),
        "linear_trace_kappa_rms": _rms(fit["kappa_linear_trace"]),
        "determinant_correction_rms": detD_rms,
        "observed_vs_production_pearson": p_obs_prod,
        "observed_vs_production_spearman": s_obs_prod,
        "observed_vs_production_count": n_obs_prod,
        "observed_vs_exact_identity_pearson": p_obs_exact,
        "observed_vs_exact_identity_spearman": s_obs_exact,
        "observed_vs_exact_identity_count": n_obs_exact,
        "observed_vs_linear_minus_trace_pearson": p_obs_linear,
        "observed_vs_linear_minus_trace_spearman": s_obs_linear,
        "observed_vs_linear_minus_trace_count": n_obs_linear,
        "observed_vs_standard_minus_half_trace_pearson": p_obs_half,
        "observed_vs_standard_minus_half_trace_spearman": s_obs_half,
        "observed_vs_standard_minus_half_trace_count": n_obs_half,
        "observed_vs_positive_trace_pearson": p_obs_pos,
        "observed_vs_positive_trace_spearman": s_obs_pos,
        "observed_vs_positive_trace_count": n_obs_pos,
        "observed_vs_positive_half_trace_pearson": p_obs_poshalf,
        "observed_vs_positive_half_trace_spearman": s_obs_poshalf,
        "observed_vs_positive_half_trace_count": n_obs_poshalf,
        "los_magnitude_vs_observed_pearson": p_los_obs,
        "los_magnitude_vs_observed_spearman": s_los_obs,
        "los_magnitude_vs_observed_count": n_los_obs,
        "los_magnitude_vs_production_kappa_pearson": p_los_prod,
        "los_magnitude_vs_production_kappa_spearman": s_los_prod,
        "los_magnitude_vs_production_kappa_count": n_los_prod,
        "ray_classification": artifact.statistics["ray_classification"],
        "identity_gate_pass": bool(
            mask_match
            and production_reconstruction_error <= IDENTITY_TOL
            and exact_identity_error <= IDENTITY_TOL
        ),
        "linear_trace_retains_production_morphology": bool(
            np.isfinite(p_prod_linear) and p_prod_linear >= LINEAR_RETENTION_THRESHOLD
        ),
    }

    cdir = OUT / "clusters" / cid
    cdir.mkdir(parents=True, exist_ok=True)
    _write_json(cdir / "jacobian_consistency_metrics.json", stats)
    np.savez_compressed(
        cdir / "jacobian_checkpoint_fields.npz",
        production_kappa=k_prod,
        observed_kappa=observed,
        los_Rx=Rx,
        los_Ry=Ry,
        **{k:v for k,v in fit.items() if isinstance(v, np.ndarray)},
    )
    return stats


def main() -> int:
    started = time.perf_counter()
    OUT.mkdir(parents=True, exist_ok=True)
    repo = _repo_state()
    _write_json(OUT / "repository_state.json", repo)
    if repo["branch"] != "main" or repo["tracked_changes"] or repo["staged_changes"]:
        v = {"lab_id":LAB_ID,"outcome":"REPOSITORY_GATE_FAILURE","head_sha":repo["head_sha"]}
        _write_json(OUT / "validation.json", v)
        print(json.dumps(v, indent=2))
        return 2

    fixture = _affine_fixture()
    _write_csv(OUT / "analytic_identity_fixtures.csv", fixture["rows"])
    if not fixture["passes"]:
        raise RuntimeError("analytic determinant identity fixture failed")

    rows = []
    failures = []
    for cluster in COV.CLUSTERS:
        print(f"[{cluster['id']}] C25 rays -> independent J -> determinant/trace audit")
        try:
            rows.append(_run_cluster(cluster))
        except Exception as exc:
            failures.append({"cluster_id":cluster["id"],"error_type":type(exc).__name__,"message":str(exc)})
            _write_json(OUT / "cluster_failures.json", failures)
            raise

    _write_csv(OUT / "jacobian_consistency_summary.csv", rows)
    _write_json(OUT / "cluster_failures.json", failures)

    identity_all = all(r["identity_gate_pass"] for r in rows)
    linear_all = all(r["linear_trace_retains_production_morphology"] for r in rows)
    if not identity_all:
        outcome = "Outcome C — PRODUCTION JACOBIAN IMPLEMENTATION INCONSISTENCY"
    elif linear_all:
        outcome = "Outcome A — JACOBIAN DETERMINANT MAPPING INTERNALLY CONSISTENT; MORPHOLOGY LOSS PRECEDES DETERMINANT EXTRACTION"
    else:
        outcome = "Outcome B — JACOBIAN DETERMINANT NONLINEARITY MATERIALLY CHANGES RAY-DERIVATIVE MORPHOLOGY"

    def mean_finite(key):
        vals = np.asarray([r[key] for r in rows], dtype=np.float64)
        vals = vals[np.isfinite(vals)]
        return float(np.mean(vals)) if vals.size else float("nan")

    validation = {
        "lab_id": LAB_ID,
        "outcome": outcome,
        "head_sha": repo["head_sha"],
        "candidate_id": COV.CANDIDATE_ID,
        "physical_source_representation": COV.PHYSICAL_SOURCE,
        "coverage_lane": "C25",
        "cluster_count_expected": len(COV.CLUSTERS),
        "cluster_count_completed": len(rows),
        "analytic_identity_fixture_pass": bool(fixture["passes"]),
        "all_cluster_identity_gates_pass": bool(identity_all),
        "all_clusters_linear_trace_retains_production_morphology": bool(linear_all),
        "identity_tolerance": IDENTITY_TOL,
        "linear_retention_threshold": LINEAR_RETENTION_THRESHOLD,
        "mean_production_vs_independent_relative_error": mean_finite("production_vs_independent_relative_error"),
        "mean_exact_det_identity_relative_error": mean_finite("exact_det_identity_relative_error"),
        "mean_production_vs_linear_trace_pearson": mean_finite("production_vs_linear_trace_pearson"),
        "mean_nonlinear_detD_to_traceD_rms_ratio": mean_finite("nonlinear_detD_to_traceD_rms_ratio"),
        "mean_observed_vs_production_pearson": mean_finite("observed_vs_production_pearson"),
        "mean_observed_vs_linear_minus_trace_pearson": mean_finite("observed_vs_linear_minus_trace_pearson"),
        "mean_los_magnitude_vs_observed_pearson": mean_finite("los_magnitude_vs_observed_pearson"),
        "physics_change_authorized": False,
        "candidate_change_authorized": False,
        "jacobian_change_authorized": False,
        "next_experiment_authorized": False,
        "science_interpretation_required": True,
        "duration_seconds": time.perf_counter()-started,
    }
    _write_json(OUT / "validation.json", validation)
    _write_json(OUT / "run.json", {
        "lab_id":LAB_ID,
        "head_sha":repo["head_sha"],
        "coverage_lane":"C25",
        "n_photons":C25_NPHOTONS,
        "bins":BINS,
        "duration_seconds":validation["duration_seconds"],
    })

    report = [
        f"# {LAB_ID}", "", f"**Outcome:** {outcome}", "",
        "Exact 2D identity audited: `1-det(J) = -tr(J-I)-det(J-I)`.",
        "Trace/sign lanes are diagnostic only; no formula change is authorized by this lab.",
    ]
    (OUT / "report.md").write_text("\n".join(report)+"\n")
    print(json.dumps(validation, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
