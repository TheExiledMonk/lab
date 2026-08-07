#!/usr/bin/env python3
"""PBUF FOUNDATION — G3D DISPLACEMENT TOPOLOGY 001.

Diagnostic-only follow-up to LOS-consistent ray geometry 001.

The PBUF candidate, M10 interface field, C25 source geometry, photon count,
normalized propagation law, step size/count, and G3D LOS-consistent geometry
are frozen.  No conventional gravitational law is introduced into the PBUF
pipeline.

This lab treats the public observed-kappa map only as an external morphology
benchmark.  It decomposes the PBUF-produced transverse displacement-gradient
tensor into purely geometric 2D modes:

  divergence-like       = d(ax)/dx + d(ay)/dy
  curl-like             = d(ay)/dx - d(ax)/dy
  symmetric-plus-like   = d(ax)/dx - d(ay)/dy
  symmetric-cross-like  = d(ax)/dy + d(ay)/dx

where a=(ax,ay) is the G3D transverse displacement field.  Signed and magnitude
lanes are reported without selecting a preferred mode or sign.

The decomposition obeys the exact 2x2 Frobenius identity

  ||D||_F^2 = 1/2 * (div^2 + curl^2 + plus^2 + cross^2)

for D = grad(a).  The lab verifies that identity and independently verifies
that -divergence matches the prior G3D linear-deformation diagnostic.  These
are numerical/geometric consistency gates only, not gravitational assumptions.
"""
from __future__ import annotations

import csv
import json
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))

import pbuf.labs.foundation.m10_coverage_25pct_science001 as BASE
import pbuf.labs.foundation.los_consistent_ray_geometry001 as GEO
from pbuf.core import los_projection as M14
from pbuf.core import observable_extraction as M16

LAB_ID = "PBUF-FOUNDATION-G3D-DISPLACEMENT-TOPOLOGY-001"
OUT = ROOT / "runs" / "g3d_displacement_topology001"
CHECKPOINTS = GEO.CHECKPOINTS
EXPECTED_SUPPORT = GEO.EXPECTED_SUPPORT
UNIT_SPEED_TOL = GEO.UNIT_SPEED_TOL
FIRST_STEP_TOL = GEO.FIRST_STEP_TOL
TOPOLOGY_IDENTITY_TOL = 1e-12
LINEAR_CONSISTENCY_TOL = 1e-12


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
    return str(o)


def _write_json(path: Path, obj) -> None:
    path.write_text(json.dumps(obj, indent=2, default=_json_default))


def _write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        path.write_text("")
        return
    keys: list[str] = []
    for row in rows:
        for key in row:
            if key not in keys:
                keys.append(key)
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)


def _corr(a, b):
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    if a.shape != b.shape:
        raise RuntimeError(f"correlation shape mismatch {a.shape} vs {b.shape}")
    mask = np.isfinite(a) & np.isfinite(b)
    n = int(mask.sum())
    if n < 2:
        return float("nan"), float("nan"), n
    return (
        float(M16.safe_pearson(a[mask], b[mask])),
        float(M16.safe_spearman(a[mask], b[mask])),
        n,
    )


def _rms(a) -> float:
    x = np.asarray(a, dtype=np.float64)
    return float(np.sqrt(np.mean(x*x))) if x.size else float("nan")


def _safe_rel_rms(diff, reference) -> float:
    d = np.asarray(diff, dtype=np.float64)
    r = np.asarray(reference, dtype=np.float64)
    mask = np.isfinite(d) & np.isfinite(r)
    if not np.any(mask):
        return float("nan")
    denom = _rms(r[mask])
    num = _rms(d[mask])
    if denom <= 1e-30:
        return 0.0 if num <= 1e-30 else float("inf")
    return float(num / denom)


def _gradient_tensor_maps(x0, y0, xf, yf, groups):
    """Fit D=grad(final-initial transverse displacement) in each source bin."""
    shape = (BASE.OBS_BINS, BASE.OBS_BINS)
    names = ("Dxx", "Dxy", "Dyx", "Dyy")
    out = {name: np.full(shape, np.nan, dtype=np.float64) for name in names}

    dx_all = np.asarray(xf, dtype=np.float64) - np.asarray(x0, dtype=np.float64)
    dy_all = np.asarray(yf, dtype=np.float64) - np.asarray(y0, dtype=np.float64)

    for q, idx in groups.items():
        r, c = divmod(q, BASE.OBS_BINS)
        xi = np.asarray(x0[idx], dtype=np.float64)
        yi = np.asarray(y0[idx], dtype=np.float64)
        A = np.column_stack([xi - xi.mean(), yi - yi.mean()])
        bx = dx_all[idx] - np.mean(dx_all[idx])
        by = dy_all[idx] - np.mean(dy_all[idx])
        try:
            gx, *_ = np.linalg.lstsq(A, bx, rcond=None)
            gy, *_ = np.linalg.lstsq(A, by, rcond=None)
        except np.linalg.LinAlgError:
            continue
        out["Dxx"][r, c] = float(gx[0])
        out["Dxy"][r, c] = float(gx[1])
        out["Dyx"][r, c] = float(gy[0])
        out["Dyy"][r, c] = float(gy[1])

    Dxx, Dxy, Dyx, Dyy = out["Dxx"], out["Dxy"], out["Dyx"], out["Dyy"]
    divergence = Dxx + Dyy
    curl = Dyx - Dxy
    plus = Dxx - Dyy
    cross = Dxy + Dyx
    shear_like_mag = np.hypot(plus, cross)
    frobenius = np.sqrt(Dxx*Dxx + Dxy*Dxy + Dyx*Dyx + Dyy*Dyy)

    # Orthogonal 2x2 decomposition energies per pixel.
    isotropic_energy = 0.5 * divergence * divergence
    rotational_energy = 0.5 * curl * curl
    symmetric_traceless_energy = 0.5 * (plus*plus + cross*cross)
    total_energy = frobenius * frobenius
    identity_rhs = isotropic_energy + rotational_energy + symmetric_traceless_energy

    out.update({
        "divergence_like": divergence,
        "negative_divergence_like": -divergence,
        "abs_divergence_like": np.abs(divergence),
        "curl_like": curl,
        "abs_curl_like": np.abs(curl),
        "symmetric_plus_like": plus,
        "abs_symmetric_plus_like": np.abs(plus),
        "symmetric_cross_like": cross,
        "abs_symmetric_cross_like": np.abs(cross),
        "symmetric_traceless_mag": shear_like_mag,
        "gradient_frobenius_mag": frobenius,
        "isotropic_mode_energy": isotropic_energy,
        "rotational_mode_energy": rotational_energy,
        "symmetric_traceless_mode_energy": symmetric_traceless_energy,
        "total_gradient_energy": total_energy,
        "topology_identity_rhs": identity_rhs,
    })
    return out


def _displacement_maps(x0, y0, xf, yf, groups):
    dx = np.asarray(xf, dtype=np.float64) - np.asarray(x0, dtype=np.float64)
    dy = np.asarray(yf, dtype=np.float64) - np.asarray(y0, dtype=np.float64)
    return {
        "displacement_x": GEO._mean_map(dx, groups),
        "displacement_y": GEO._mean_map(dy, groups),
        "abs_displacement_x": GEO._mean_map(np.abs(dx), groups),
        "abs_displacement_y": GEO._mean_map(np.abs(dy), groups),
        "displacement_mag": GEO._mean_map(np.hypot(dx, dy), groups),
    }


def _mode_correlations(row: dict, fields: dict, observed, los_mag) -> None:
    compare_names = (
        "displacement_x",
        "displacement_y",
        "abs_displacement_x",
        "abs_displacement_y",
        "displacement_mag",
        "divergence_like",
        "negative_divergence_like",
        "abs_divergence_like",
        "curl_like",
        "abs_curl_like",
        "symmetric_plus_like",
        "abs_symmetric_plus_like",
        "symmetric_cross_like",
        "abs_symmetric_cross_like",
        "symmetric_traceless_mag",
        "gradient_frobenius_mag",
    )
    for name in compare_names:
        arr = fields[name]
        p, s, n = _corr(arr, observed)
        p2, s2, n2 = _corr(arr, los_mag)
        finite = arr[np.isfinite(arr)]
        row[f"{name}_vs_observed_pearson"] = p
        row[f"{name}_vs_observed_spearman"] = s
        row[f"{name}_vs_observed_count"] = n
        row[f"{name}_vs_los_mag_pearson"] = p2
        row[f"{name}_vs_los_mag_spearman"] = s2
        row[f"{name}_vs_los_mag_count"] = n2
        row[f"{name}_rms"] = _rms(finite) if finite.size else float("nan")


def _energy_fractions(fields: dict) -> dict:
    mask = np.isfinite(fields["total_gradient_energy"])
    total = float(np.sum(fields["total_gradient_energy"][mask])) if np.any(mask) else 0.0
    iso = float(np.sum(fields["isotropic_mode_energy"][mask])) if np.any(mask) else 0.0
    rot = float(np.sum(fields["rotational_mode_energy"][mask])) if np.any(mask) else 0.0
    sym = float(np.sum(fields["symmetric_traceless_mode_energy"][mask])) if np.any(mask) else 0.0
    if total <= 1e-30:
        return {
            "isotropic_mode_energy_fraction": float("nan"),
            "rotational_mode_energy_fraction": float("nan"),
            "symmetric_traceless_mode_energy_fraction": float("nan"),
        }
    return {
        "isotropic_mode_energy_fraction": iso / total,
        "rotational_mode_energy_fraction": rot / total,
        "symmetric_traceless_mode_energy_fraction": sym / total,
    }


def _checkpoint_topology(cid, k, snap, x0, y0, groups, observed, los_mag):
    disp = _displacement_maps(x0, y0, snap["x"], snap["y"], groups)
    topo = _gradient_tensor_maps(x0, y0, snap["x"], snap["y"], groups)
    fields = {**disp, **topo}

    row = {
        "cluster_id": cid,
        "step_index": int(k),
        "propagation_distance": float(k * BASE.CFG["step"]),
    }
    _mode_correlations(row, fields, observed, los_mag)
    row.update(_energy_fractions(fields))

    lhs = fields["total_gradient_energy"]
    rhs = fields["topology_identity_rhs"]
    row["topology_identity_relative_rms_error"] = _safe_rel_rms(lhs-rhs, lhs)
    return row, fields


def _run_cluster(cluster):
    cid = cluster["id"]
    real = BASE._load_cluster(cluster)
    state = BASE._evolve(BASE._initial_state(real["rho3"]))
    candidate = BASE._candidate(state)
    vector = BASE._interface_vector(candidate)
    los = M14.project_vector_to_image_plane(*vector, los_axis="z")
    Rx = np.asarray(los["comp_1"], dtype=np.float64)
    Ry = np.asarray(los["comp_2"], dtype=np.float64)
    los_mag = np.hypot(Rx, Ry)
    observed = np.asarray(real["observed_kappa"], dtype=np.float64)
    grid = np.linspace(-BASE.CFG["extent"], BASE.CFG["extent"], Rx.shape[0])
    field = {"xgrid": grid, "ygrid": grid, "rx": Rx, "ry": Ry}

    x0, y0, _, _ = BASE._launch_expanded_25pct()
    groups = GEO._source_groups(x0, y0)
    if len(groups) != EXPECTED_SUPPORT:
        raise RuntimeError(f"{cid}: expected {EXPECTED_SUPPORT} source bins, got {len(groups)}")

    checkpoints, g3d = GEO._propagate_g3d(
        field, BASE.CFG["step"], BASE.CFG["steps"], x0, y0
    )
    if g3d["max_unit_speed_error"] > UNIT_SPEED_TOL:
        raise RuntimeError(f"{cid}: G3D unit-speed gate failed: {g3d['max_unit_speed_error']}")

    first = GEO._first_step_geometry(field, x0, y0, checkpoints[1], observed, los_mag)
    if not first["first_step_exact_pass"]:
        raise RuntimeError(f"{cid}: G3D first-step exact geometry gate failed")

    rows = []
    all_fields = {}
    for k in CHECKPOINTS:
        row, fields = _checkpoint_topology(
            cid, k, checkpoints[k], x0, y0, groups, observed, los_mag
        )
        if np.isfinite(row["topology_identity_relative_rms_error"]) and row["topology_identity_relative_rms_error"] > TOPOLOGY_IDENTITY_TOL:
            raise RuntimeError(
                f"{cid}: topology identity failed at step {k}: "
                f"{row['topology_identity_relative_rms_error']}"
            )
        rows.append(row)
        all_fields[f"step_{k}"] = fields

    final_fields = all_fields[f"step_{CHECKPOINTS[-1]}"]
    prior_linear = GEO._linear_kappa_map(x0, y0, g3d["x"], g3d["y"], groups)
    neg_div = final_fields["negative_divergence_like"]
    linear_rel = _safe_rel_rms(neg_div-prior_linear, prior_linear)
    if linear_rel > LINEAR_CONSISTENCY_TOL:
        raise RuntimeError(f"{cid}: negative-divergence/prior-linear consistency failed: {linear_rel}")

    final_row = rows[-1]
    summary = {
        "cluster_id": cid,
        "candidate_id": BASE.CANDIDATE_ID,
        "physical_source_representation": BASE.PHYSICAL_SOURCE,
        "geometry_lane": "G3D_LOS_consistent_diagnostic",
        "benchmark_role": "external_morphology_comparison_only",
        "n_photons": int(len(x0)),
        "source_supported_bins": int(len(groups)),
        "checkpoint_count": len(rows),
        "g3d_unit_speed_max_error": float(g3d["max_unit_speed_error"]),
        "g3d_unit_speed_pass": bool(g3d["max_unit_speed_error"] <= UNIT_SPEED_TOL),
        "first_step_exact_max_vector_error": first["first_step_exact_max_vector_error"],
        "first_step_exact_pass": first["first_step_exact_pass"],
        "los_mag_vs_observed_pearson": _corr(los_mag, observed)[0],
        "los_mag_vs_observed_spearman": _corr(los_mag, observed)[1],
        "final_topology_identity_relative_rms_error": final_row["topology_identity_relative_rms_error"],
        "negative_divergence_vs_prior_linear_relative_rms_error": linear_rel,
        "negative_divergence_vs_prior_linear_pearson": _corr(neg_div, prior_linear)[0],
        "negative_divergence_vs_prior_linear_spearman": _corr(neg_div, prior_linear)[1],
        "final_displacement_mag_vs_observed_pearson": final_row["displacement_mag_vs_observed_pearson"],
        "final_displacement_mag_vs_observed_spearman": final_row["displacement_mag_vs_observed_spearman"],
        "final_divergence_like_vs_observed_pearson": final_row["divergence_like_vs_observed_pearson"],
        "final_divergence_like_vs_observed_spearman": final_row["divergence_like_vs_observed_spearman"],
        "final_negative_divergence_like_vs_observed_pearson": final_row["negative_divergence_like_vs_observed_pearson"],
        "final_negative_divergence_like_vs_observed_spearman": final_row["negative_divergence_like_vs_observed_spearman"],
        "final_abs_divergence_like_vs_observed_pearson": final_row["abs_divergence_like_vs_observed_pearson"],
        "final_abs_divergence_like_vs_observed_spearman": final_row["abs_divergence_like_vs_observed_spearman"],
        "final_curl_like_vs_observed_pearson": final_row["curl_like_vs_observed_pearson"],
        "final_curl_like_vs_observed_spearman": final_row["curl_like_vs_observed_spearman"],
        "final_abs_curl_like_vs_observed_pearson": final_row["abs_curl_like_vs_observed_pearson"],
        "final_abs_curl_like_vs_observed_spearman": final_row["abs_curl_like_vs_observed_spearman"],
        "final_symmetric_plus_like_vs_observed_pearson": final_row["symmetric_plus_like_vs_observed_pearson"],
        "final_symmetric_cross_like_vs_observed_pearson": final_row["symmetric_cross_like_vs_observed_pearson"],
        "final_symmetric_traceless_mag_vs_observed_pearson": final_row["symmetric_traceless_mag_vs_observed_pearson"],
        "final_symmetric_traceless_mag_vs_observed_spearman": final_row["symmetric_traceless_mag_vs_observed_spearman"],
        "final_gradient_frobenius_mag_vs_observed_pearson": final_row["gradient_frobenius_mag_vs_observed_pearson"],
        "final_gradient_frobenius_mag_vs_observed_spearman": final_row["gradient_frobenius_mag_vs_observed_spearman"],
        "final_isotropic_mode_energy_fraction": final_row["isotropic_mode_energy_fraction"],
        "final_rotational_mode_energy_fraction": final_row["rotational_mode_energy_fraction"],
        "final_symmetric_traceless_mode_energy_fraction": final_row["symmetric_traceless_mode_energy_fraction"],
    }
    return summary, rows, all_fields, {
        "los_mag": los_mag,
        "observed_benchmark": observed,
        "g3d_x_final": g3d["x"],
        "g3d_y_final": g3d["y"],
        "g3d_z_final": g3d["z"],
        "prior_linear_diagnostic": prior_linear,
    }


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
        }
        _write_json(OUT / "validation.json", validation)
        print(json.dumps(validation, indent=2))
        return 2

    summaries = []
    checkpoint_rows = []
    failures = []
    for cluster in BASE.CLUSTERS:
        cid = cluster["id"]
        print(f"[{cid}] G3D transverse-displacement topology decomposition")
        try:
            summary, rows, fields, final_arrays = _run_cluster(cluster)
            summaries.append(summary)
            checkpoint_rows.extend(rows)
            cdir = OUT / "clusters" / cid
            cdir.mkdir(parents=True, exist_ok=True)
            _write_json(cdir / "topology_summary.json", summary)
            _write_csv(cdir / "topology_checkpoints.csv", rows)
            npz = {}
            for step_name, fd in fields.items():
                for name, arr in fd.items():
                    npz[f"{step_name}__{name}"] = arr
            npz.update(final_arrays)
            np.savez_compressed(cdir / "topology_checkpoint_fields.npz", **npz)
        except Exception as exc:
            failures.append({"cluster_id": cid, "error": repr(exc)})
            _write_json(OUT / "cluster_failures.json", failures)
            raise

    _write_csv(OUT / "topology_summary.csv", summaries)
    _write_csv(OUT / "topology_checkpoint_summary.csv", checkpoint_rows)
    _write_json(OUT / "cluster_failures.json", failures)

    validation = {
        "lab_id": LAB_ID,
        "outcome": "Outcome A — G3D DISPLACEMENT TOPOLOGY AUDIT COMPLETE",
        "head_sha": repo["head_sha"],
        "candidate_id": BASE.CANDIDATE_ID,
        "physical_source_representation": BASE.PHYSICAL_SOURCE,
        "geometry_lane": "G3D_LOS_consistent_diagnostic",
        "benchmark_role": "external_morphology_comparison_only",
        "cluster_count_expected": len(BASE.CLUSTERS),
        "cluster_count_completed": len(summaries),
        "checkpoint_steps": list(CHECKPOINTS),
        "topology_identity_tolerance": TOPOLOGY_IDENTITY_TOL,
        "linear_consistency_tolerance": LINEAR_CONSISTENCY_TOL,
        "all_cluster_g3d_unit_speed_pass": bool(all(r["g3d_unit_speed_pass"] for r in summaries)),
        "all_cluster_first_step_exact_pass": bool(all(r["first_step_exact_pass"] for r in summaries)),
        "all_cluster_topology_identity_pass": bool(all(r["final_topology_identity_relative_rms_error"] <= TOPOLOGY_IDENTITY_TOL for r in summaries)),
        "all_cluster_negative_divergence_linear_consistency_pass": bool(all(r["negative_divergence_vs_prior_linear_relative_rms_error"] <= LINEAR_CONSISTENCY_TOL for r in summaries)),
        "mean_los_mag_vs_observed_pearson": float(np.mean([r["los_mag_vs_observed_pearson"] for r in summaries])),
        "mean_final_displacement_mag_vs_observed_pearson": float(np.mean([r["final_displacement_mag_vs_observed_pearson"] for r in summaries])),
        "mean_final_negative_divergence_vs_observed_pearson": float(np.mean([r["final_negative_divergence_like_vs_observed_pearson"] for r in summaries])),
        "mean_final_abs_divergence_vs_observed_pearson": float(np.mean([r["final_abs_divergence_like_vs_observed_pearson"] for r in summaries])),
        "mean_final_abs_curl_vs_observed_pearson": float(np.mean([r["final_abs_curl_like_vs_observed_pearson"] for r in summaries])),
        "mean_final_symmetric_traceless_mag_vs_observed_pearson": float(np.mean([r["final_symmetric_traceless_mag_vs_observed_pearson"] for r in summaries])),
        "mean_final_gradient_frobenius_mag_vs_observed_pearson": float(np.mean([r["final_gradient_frobenius_mag_vs_observed_pearson"] for r in summaries])),
        "mean_final_isotropic_mode_energy_fraction": float(np.mean([r["final_isotropic_mode_energy_fraction"] for r in summaries])),
        "mean_final_rotational_mode_energy_fraction": float(np.mean([r["final_rotational_mode_energy_fraction"] for r in summaries])),
        "mean_final_symmetric_traceless_mode_energy_fraction": float(np.mean([r["final_symmetric_traceless_mode_energy_fraction"] for r in summaries])),
        "physics_change_authorized": False,
        "candidate_change_authorized": False,
        "production_geometry_change_authorized": False,
        "topology_mode_selection_authorized": False,
        "next_experiment_authorized": False,
        "science_interpretation_required": True,
        "duration_seconds": float(time.perf_counter() - started),
    }
    _write_json(OUT / "validation.json", validation)
    _write_json(OUT / "run.json", {
        "lab_id": LAB_ID,
        "head_sha": repo["head_sha"],
        "candidate_id": BASE.CANDIDATE_ID,
        "physical_source_representation": BASE.PHYSICAL_SOURCE,
        "geometry_lane": "G3D_LOS_consistent_diagnostic",
        "benchmark_role": "external_morphology_comparison_only",
        "checkpoint_steps": list(CHECKPOINTS),
        "note": "No conventional gravitational law is injected; observed kappa is an external morphology benchmark only.",
    })
    print(json.dumps(validation, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
