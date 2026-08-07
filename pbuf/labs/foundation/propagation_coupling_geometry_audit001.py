#!/usr/bin/env python3
"""PBUF FOUNDATION — PROPAGATION COUPLING GEOMETRY AUDIT 001.

Diagnostic-only five-cluster C25 audit following propagation morphology
localization. The frozen PL1_PM1_PS2/M10 field, C25 launch, propagation rule,
Jacobian rule, and all physical coefficients remain unchanged.

This lab asks which part of the M10 image-plane response is actually coupled
into unit-speed photon motion. At selected checkpoints it decomposes the local
response relative to the instantaneous photon direction v:

    R_parallel = (R dot v) v
    R_perp     = R - R_parallel

and records signed/image-plane components, parallel/perpendicular magnitudes,
actual velocity changes, longitudinal/transverse displacements, and the local
linear deformation kappa = -tr(J-I).

For the first update, where v0=(1,0), it also verifies the exact normalized
update and the small-step relation dy ~= h^2 Ry. No propagation or physics
change is authorized by this lab.
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

from weak_lensing_observation001 import propagate as production_propagate
from pbuf.core import los_projection as M14
from pbuf.core import observable_extraction as M16
import pbuf.labs.foundation.m10_coverage_25pct_science001 as BASE

LAB_ID = "PBUF-FOUNDATION-PROPAGATION-COUPLING-GEOMETRY-AUDIT-001"
OUT = ROOT / "runs" / "propagation_coupling_geometry_audit001"
CHECKPOINTS = (0, 1, 5, 10, 20, 40, 80, 120, 159)
PARITY_TOL = 1e-12
FIRST_STEP_EXACT_TOL = 1e-14


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
        w.writeheader(); w.writerows(rows)


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


def _rms(a):
    x = np.asarray(a, dtype=np.float64)
    return float(np.sqrt(np.mean(x*x))) if x.size else float("nan")


def _sample_field(field, x, y):
    xgrid = field["xgrid"]; ygrid = field["ygrid"]
    ix = np.clip(np.searchsorted(xgrid, x) - 1, 0, len(xgrid) - 1)
    iy = np.clip(np.searchsorted(ygrid, y) - 1, 0, len(ygrid) - 1)
    return field["rx"][iy, ix], field["ry"][iy, ix]


def _source_groups(x0, y0):
    edges = np.linspace(-BASE.CFG["extent"], BASE.CFG["extent"], BASE.OBS_BINS + 1)
    col = np.searchsorted(edges, x0, side="right") - 1
    row = np.searchsorted(edges, y0, side="right") - 1
    valid = (row >= 0) & (row < BASE.OBS_BINS) & (col >= 0) & (col < BASE.OBS_BINS)
    flat = row * BASE.OBS_BINS + col
    groups = {}
    for q in np.unique(flat[valid]):
        idx = np.where(valid & (flat == q))[0]
        if idx.size >= 6:
            groups[int(q)] = idx
    return groups


def _mean_map(values, groups):
    arr = np.full((BASE.OBS_BINS, BASE.OBS_BINS), np.nan)
    values = np.asarray(values, dtype=np.float64)
    for q, idx in groups.items():
        r, c = divmod(q, BASE.OBS_BINS)
        arr[r, c] = float(np.mean(values[idx]))
    return arr


def _linear_kappa_map(x0, y0, x, y, groups):
    out = np.full((BASE.OBS_BINS, BASE.OBS_BINS), np.nan)
    for q, idx in groups.items():
        r, c = divmod(q, BASE.OBS_BINS)
        xi = x0[idx]; yi = y0[idx]; xf = x[idx]; yf = y[idx]
        A = np.column_stack([xi - xi.mean(), yi - yi.mean()])
        try:
            Jx, *_ = np.linalg.lstsq(A, xf - xf.mean(), rcond=None)
            Jy, *_ = np.linalg.lstsq(A, yf - yf.mean(), rcond=None)
        except np.linalg.LinAlgError:
            continue
        J = np.array([[Jx[0], Jx[1]], [Jy[0], Jy[1]]], dtype=np.float64)
        out[r, c] = -float(np.trace(J - np.eye(2)))
    return out


def _decompose(rx, ry, vx, vy):
    dot = rx*vx + ry*vy
    parx = dot*vx; pary = dot*vy
    perpx = rx - parx; perpy = ry - pary
    return {
        "dot": dot,
        "parallel_x": parx,
        "parallel_y": pary,
        "parallel_mag": np.hypot(parx, pary),
        "perp_x": perpx,
        "perp_y": perpy,
        "perp_mag": np.hypot(perpx, perpy),
    }


def _checkpoint_propagate(field, step, steps, x0, y0, vx0, vy0):
    wanted = set(int(k) for k in CHECKPOINTS)
    if max(wanted) >= steps:
        raise RuntimeError(f"checkpoint >= steps: max={max(wanted)} steps={steps}")
    x = x0.copy(); y = y0.copy(); vx = vx0.copy(); vy = vy0.copy()
    last_dvx = np.zeros_like(vx); last_dvy = np.zeros_like(vy)
    out = {}

    def capture(k):
        rx, ry = _sample_field(field, x, y)
        dec = _decompose(rx, ry, vx, vy)
        out[int(k)] = {
            "x": x.copy(), "y": y.copy(), "vx": vx.copy(), "vy": vy.copy(),
            "rx": rx.copy(), "ry": ry.copy(),
            "last_dvx": last_dvx.copy(), "last_dvy": last_dvy.copy(),
            **{key: np.asarray(val).copy() for key, val in dec.items()},
        }

    if 0 in wanted: capture(0)
    for k in range(1, steps):
        rx, ry = _sample_field(field, x, y)
        vx_old = vx.copy(); vy_old = vy.copy()
        vx_new = vx + step*rx
        vy_new = vy + step*ry
        scale = np.maximum(np.hypot(vx_new, vy_new), 1e-12)
        vx = vx_new/scale; vy = vy_new/scale
        last_dvx = vx - vx_old; last_dvy = vy - vy_old
        x = x + step*vx; y = y + step*vy
        if k in wanted: capture(k)
    return out, (x, y)


def _field_metrics(prefix, field, observed, los_mag, row):
    p, s, n = _corr(field, observed)
    row[f"{prefix}_vs_observed_pearson"] = p
    row[f"{prefix}_vs_observed_spearman"] = s
    row[f"{prefix}_vs_observed_count"] = n
    p, s, n = _corr(field, los_mag)
    row[f"{prefix}_vs_los_mag_pearson"] = p
    row[f"{prefix}_vs_los_mag_spearman"] = s
    row[f"{prefix}_vs_los_mag_count"] = n
    finite = np.asarray(field)[np.isfinite(field)]
    row[f"{prefix}_rms"] = _rms(finite) if finite.size else float("nan")


def _checkpoint_metrics(cid, k, snap, x0, y0, vx0, vy0, groups, los_mag, observed):
    response_mag = np.hypot(snap["rx"], snap["ry"])
    accumulated_dvx = snap["vx"] - vx0
    accumulated_dvy = snap["vy"] - vy0
    longitudinal_displacement = snap["x"] - x0 - float(k*BASE.CFG["step"])
    transverse_displacement = snap["y"] - y0

    raw = {
        "response_mag": response_mag,
        "response_x": snap["rx"],
        "response_y": snap["ry"],
        "parallel_mag": snap["parallel_mag"],
        "perp_mag": snap["perp_mag"],
        "perp_x": snap["perp_x"],
        "perp_y": snap["perp_y"],
        "last_delta_vx": snap["last_dvx"],
        "last_delta_vy": snap["last_dvy"],
        "last_delta_v_mag": np.hypot(snap["last_dvx"], snap["last_dvy"]),
        "accumulated_delta_vx": accumulated_dvx,
        "accumulated_delta_vy": accumulated_dvy,
        "accumulated_delta_v_mag": np.hypot(accumulated_dvx, accumulated_dvy),
        "longitudinal_displacement_residual": longitudinal_displacement,
        "transverse_displacement": transverse_displacement,
        "abs_transverse_displacement": np.abs(transverse_displacement),
    }
    maps = {name: _mean_map(value, groups) for name, value in raw.items()}
    maps["linear_kappa"] = _linear_kappa_map(x0, y0, snap["x"], snap["y"], groups)

    row = {
        "cluster_id": cid,
        "step_index": int(k),
        "propagation_distance": float(k*BASE.CFG["step"]),
    }
    for name, field in maps.items():
        _field_metrics(name, field, observed, los_mag, row)
    return row, maps


def _first_step_audit(field, checkpoints, x0, y0, groups, observed, los_mag):
    h = float(BASE.CFG["step"])
    rx0, ry0 = _sample_field(field, x0, y0)
    denom = np.maximum(np.hypot(1.0 + h*rx0, h*ry0), 1e-12)
    vy_exact = h*ry0/denom
    dy_exact = h*vy_exact
    dy_small = h*h*ry0
    actual_dy = checkpoints[1]["y"] - y0
    exact_err = actual_dy - dy_exact
    small_err = actual_dy - dy_small

    actual_map = _mean_map(actual_dy, groups)
    exact_map = _mean_map(dy_exact, groups)
    small_map = _mean_map(dy_small, groups)
    ry_map = _mean_map(ry0, groups)
    rx_map = _mean_map(rx0, groups)
    abs_ry_map = _mean_map(np.abs(ry0), groups)
    abs_rx_map = _mean_map(np.abs(rx0), groups)

    p_ae, s_ae, n_ae = _corr(actual_map, exact_map)
    p_as, s_as, n_as = _corr(actual_map, small_map)
    out = {
        "first_step_exact_max_abs_error": float(np.max(np.abs(exact_err))),
        "first_step_exact_rms_error": _rms(exact_err),
        "first_step_small_approx_max_abs_error": float(np.max(np.abs(small_err))),
        "first_step_small_approx_rms_error": _rms(small_err),
        "actual_dy_vs_exact_pearson": p_ae,
        "actual_dy_vs_exact_spearman": s_ae,
        "actual_dy_vs_exact_count": n_ae,
        "actual_dy_vs_small_pearson": p_as,
        "actual_dy_vs_small_spearman": s_as,
        "actual_dy_vs_small_count": n_as,
        "first_step_exact_pass": bool(np.max(np.abs(exact_err)) <= FIRST_STEP_EXACT_TOL),
    }
    for name, fmap in (
        ("initial_Rx", rx_map), ("initial_Ry", ry_map),
        ("initial_abs_Rx", abs_rx_map), ("initial_abs_Ry", abs_ry_map),
        ("first_step_actual_dy", actual_map),
        ("first_step_exact_dy", exact_map),
        ("first_step_small_dy", small_map),
    ):
        p, s, n = _corr(fmap, observed)
        out[f"{name}_vs_observed_pearson"] = p
        out[f"{name}_vs_observed_spearman"] = s
        out[f"{name}_vs_observed_count"] = n
        p, s, n = _corr(fmap, los_mag)
        out[f"{name}_vs_los_mag_pearson"] = p
        out[f"{name}_vs_los_mag_spearman"] = s
        out[f"{name}_vs_los_mag_count"] = n
    return out, {
        "initial_Rx": rx_map, "initial_Ry": ry_map,
        "initial_abs_Rx": abs_rx_map, "initial_abs_Ry": abs_ry_map,
        "first_step_actual_dy": actual_map,
        "first_step_exact_dy": exact_map,
        "first_step_small_dy": small_map,
    }


def _run_cluster(cluster):
    cid = cluster["id"]
    real = BASE._load_cluster(cluster)
    state = BASE._evolve(BASE._initial_state(real["rho3"]))
    cand = BASE._candidate(state)
    vector = BASE._interface_vector(cand)
    los = M14.project_vector_to_image_plane(*vector, los_axis="z")
    Rx = np.asarray(los["comp_1"], dtype=np.float64)
    Ry = np.asarray(los["comp_2"], dtype=np.float64)
    los_mag = np.hypot(Rx, Ry)
    field = {
        "xgrid": np.linspace(-BASE.CFG["extent"], BASE.CFG["extent"], Rx.shape[0]),
        "ygrid": np.linspace(-BASE.CFG["extent"], BASE.CFG["extent"], Ry.shape[0]),
        "rx": Rx, "ry": Ry,
    }
    x0, y0, vx0, vy0 = BASE._launch_expanded_25pct()
    groups = _source_groups(x0, y0)
    checkpoints, final_local = _checkpoint_propagate(
        field, BASE.CFG["step"], BASE.CFG["steps"], x0, y0, vx0, vy0
    )
    prod = production_propagate(
        field, BASE.CFG["step"], BASE.CFG["steps"], x0, y0, vx0, vy0
    )
    dx = np.asarray(final_local[0]) - np.asarray(prod["x"])
    dy = np.asarray(final_local[1]) - np.asarray(prod["y"])
    parity_rms = float(np.sqrt(np.mean(dx*dx + dy*dy)))
    parity_max = float(np.max(np.hypot(dx, dy)))
    if parity_max > PARITY_TOL:
        raise RuntimeError(f"{cid}: production parity failure max={parity_max}")

    rows = []
    checkpoint_fields = {}
    for k in CHECKPOINTS:
        row, maps = _checkpoint_metrics(
            cid, k, checkpoints[k], x0, y0, vx0, vy0, groups,
            los_mag, real["observed_kappa"],
        )
        rows.append(row)
        for name, arr in maps.items():
            checkpoint_fields[f"step_{k}_{name}"] = arr

    first, first_fields = _first_step_audit(
        field, checkpoints, x0, y0, groups, real["observed_kappa"], los_mag
    )
    if not first["first_step_exact_pass"]:
        raise RuntimeError(
            f"{cid}: first-step normalized-update identity failure "
            f"max={first['first_step_exact_max_abs_error']}"
        )
    checkpoint_fields.update(first_fields)
    plos, slos, _ = _corr(los_mag, real["observed_kappa"])

    summary = {
        "cluster_id": cid,
        "n_photons": int(len(x0)),
        "source_supported_bins": int(len(groups)),
        "checkpoint_count": len(rows),
        "production_parity_rms": parity_rms,
        "production_parity_max": parity_max,
        "production_parity_pass": True,
        "first_step_exact_pass": first["first_step_exact_pass"],
        "first_step_exact_max_abs_error": first["first_step_exact_max_abs_error"],
        "first_step_small_approx_rms_error": first["first_step_small_approx_rms_error"],
        "los_mag_vs_observed_pearson": plos,
        "los_mag_vs_observed_spearman": slos,
        "initial_Rx_vs_observed_pearson": first["initial_Rx_vs_observed_pearson"],
        "initial_Ry_vs_observed_pearson": first["initial_Ry_vs_observed_pearson"],
        "initial_abs_Rx_vs_observed_pearson": first["initial_abs_Rx_vs_observed_pearson"],
        "initial_abs_Ry_vs_observed_pearson": first["initial_abs_Ry_vs_observed_pearson"],
        "first_step_actual_dy_vs_observed_pearson": first["first_step_actual_dy_vs_observed_pearson"],
        "first_step_actual_dy_vs_small_pearson": first["actual_dy_vs_small_pearson"],
        "final_perp_mag_vs_observed_pearson": rows[-1]["perp_mag_vs_observed_pearson"],
        "final_accumulated_delta_vy_vs_observed_pearson": rows[-1]["accumulated_delta_vy_vs_observed_pearson"],
        "final_transverse_displacement_vs_observed_pearson": rows[-1]["transverse_displacement_vs_observed_pearson"],
        "final_linear_kappa_vs_observed_pearson": rows[-1]["linear_kappa_vs_observed_pearson"],
    }
    return {
        "summary": summary,
        "rows": rows,
        "first_step": first,
        "fields": checkpoint_fields,
    }


def main() -> int:
    started = time.perf_counter()
    OUT.mkdir(parents=True, exist_ok=True)
    repo = _repo_state()
    _write_json(OUT / "repository_state.json", repo)
    if repo["branch"] != "main" or repo["tracked_changes"] or repo["staged_changes"]:
        validation = {"lab_id": LAB_ID, "outcome": "REPOSITORY_GATE_FAILURE", "head_sha": repo["head_sha"]}
        _write_json(OUT / "validation.json", validation)
        print(json.dumps(validation, indent=2)); return 2

    summaries = []
    all_rows = []
    first_rows = []
    failures = []
    for cluster in BASE.CLUSTERS:
        cid = cluster["id"]
        print(f"[{cid}] C25 response/photon coupling geometry audit")
        try:
            result = _run_cluster(cluster)
            summaries.append(result["summary"])
            all_rows.extend(result["rows"])
            first_rows.append({"cluster_id": cid, **result["first_step"]})
            cdir = OUT / "clusters" / cid
            cdir.mkdir(parents=True, exist_ok=True)
            _write_json(cdir / "coupling_summary.json", result["summary"])
            _write_json(cdir / "first_step_geometry.json", result["first_step"])
            _write_csv(cdir / "coupling_checkpoints.csv", result["rows"])
            np.savez_compressed(cdir / "coupling_checkpoint_fields.npz", **result["fields"])
        except Exception as exc:
            failures.append({"cluster_id": cid, "error": repr(exc)})
            raise

    _write_csv(OUT / "cluster_summary.csv", summaries)
    _write_csv(OUT / "coupling_checkpoint_summary.csv", all_rows)
    _write_csv(OUT / "first_step_geometry_summary.csv", first_rows)
    _write_json(OUT / "cluster_failures.json", failures)

    all_parity = bool(all(r["production_parity_pass"] for r in summaries))
    all_first = bool(all(r["first_step_exact_pass"] for r in summaries))
    validation = {
        "lab_id": LAB_ID,
        "outcome": "Outcome A — RESPONSE-PHOTON COUPLING GEOMETRY AUDIT COMPLETE" if all_parity and all_first else "Outcome D — COUPLING GEOMETRY INTEGRITY FAILURE",
        "candidate_id": BASE.CANDIDATE_ID,
        "physical_source_representation": BASE.PHYSICAL_SOURCE,
        "coverage_lane": "C25",
        "cluster_count_expected": len(BASE.CLUSTERS),
        "cluster_count_completed": len(summaries),
        "checkpoint_steps": list(CHECKPOINTS),
        "production_parity_tolerance": PARITY_TOL,
        "first_step_exact_tolerance": FIRST_STEP_EXACT_TOL,
        "all_cluster_production_parity_pass": all_parity,
        "all_cluster_first_step_exact_pass": all_first,
        "physics_change_authorized": False,
        "candidate_change_authorized": False,
        "propagation_change_authorized": False,
        "next_experiment_authorized": False,
        "science_interpretation_required": True,
        "duration_seconds": float(time.perf_counter() - started),
    }
    _write_json(OUT / "validation.json", validation)
    _write_json(OUT / "run.json", {"validation": validation, "repository_state": repo})
    print(json.dumps(validation, indent=2))
    return 0 if all_parity and all_first else 1


if __name__ == "__main__":
    raise SystemExit(main())
