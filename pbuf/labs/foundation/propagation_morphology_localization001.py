#!/usr/bin/env python3
"""PBUF FOUNDATION — PROPAGATION MORPHOLOGY LOCALIZATION 001.

Diagnostic-only five-cluster C25 lab.  The reviewed PL1_PM1_PS2 / M10 field,
C25 source geometry, propagation parameters, and Jacobian implementation are
left unchanged.  This lab duplicates the frozen propagation loop only to expose
selected internal checkpoints and verifies that its final ray positions match
production propagation before any checkpoint result is accepted.

At selected propagation steps it measures, on the fixed initial-source 64x64
bins:
  * instantaneous sampled LOS-response magnitude seen by the rays;
  * accumulated velocity-change magnitude;
  * accumulated displacement magnitude;
  * local linear deformation kappa = -tr(J-I) from the current ray map;
  * morphology correlations against observed kappa and against the fixed LOS
    magnitude field.

No physics/source/Jacobian change is authorized by this lab.
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

LAB_ID = "PBUF-FOUNDATION-PROPAGATION-MORPHOLOGY-LOCALIZATION-001"
OUT = ROOT / "runs" / "propagation_morphology_localization001"
CHECKPOINTS = (0, 1, 5, 10, 20, 40, 80, 120, 159)
PARITY_TOL = 1e-12


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
        for k in row:
            if k not in keys: keys.append(k)
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
    return float(M16.safe_pearson(a[mask], b[mask])), float(M16.safe_spearman(a[mask], b[mask])), n


def _rms(a):
    x = np.asarray(a, dtype=np.float64)
    return float(np.sqrt(np.mean(x*x))) if x.size else float("nan")


def _sample_field(field, x, y):
    xgrid = field["xgrid"]; ygrid = field["ygrid"]
    ix = np.clip(np.searchsorted(xgrid, x) - 1, 0, len(xgrid) - 1)
    iy = np.clip(np.searchsorted(ygrid, y) - 1, 0, len(ygrid) - 1)
    return field["rx"][iy, ix], field["ry"][iy, ix]


def _checkpoint_propagate(field, step, steps, x0, y0, vx0, vy0):
    wanted = set(int(k) for k in CHECKPOINTS)
    if max(wanted) >= steps:
        raise RuntimeError(f"checkpoint >= steps: max={max(wanted)} steps={steps}")
    x = x0.copy(); y = y0.copy(); vx = vx0.copy(); vy = vy0.copy()
    out = {}

    def capture(k):
        rx_loc, ry_loc = _sample_field(field, x, y)
        out[int(k)] = {
            "x": x.copy(), "y": y.copy(), "vx": vx.copy(), "vy": vy.copy(),
            "rx_sample": rx_loc.copy(), "ry_sample": ry_loc.copy(),
        }

    if 0 in wanted: capture(0)
    for k in range(1, steps):
        rx_loc, ry_loc = _sample_field(field, x, y)
        vx_new = vx + step * rx_loc
        vy_new = vy + step * ry_loc
        scale = np.maximum(np.hypot(vx_new, vy_new), 1e-12)
        vx = vx_new / scale
        vy = vy_new / scale
        x = x + step * vx
        y = y + step * vy
        if k in wanted: capture(k)
    return out, (x, y)


def _source_bins(x0, y0):
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
        D = J - np.eye(2)
        out[r, c] = -float(np.trace(D))
    return out


def _metrics_for_checkpoint(cluster_id, k, snap, x0, y0, vx0, vy0, groups, los_mag, observed):
    sampled_mag = np.hypot(snap["rx_sample"], snap["ry_sample"])
    vel_delta_mag = np.hypot(snap["vx"] - vx0, snap["vy"] - vy0)
    displacement_mag = np.hypot(snap["x"] - x0, snap["y"] - y0)
    sampled_map = _mean_map(sampled_mag, groups)
    vel_map = _mean_map(vel_delta_mag, groups)
    disp_map = _mean_map(displacement_mag, groups)
    linear_kappa = _linear_kappa_map(x0, y0, snap["x"], snap["y"], groups)

    row = {"cluster_id": cluster_id, "step_index": int(k), "propagation_distance": float(k * BASE.CFG["step"])}
    for name, field in (
        ("sampled_response_mag", sampled_map),
        ("velocity_change_mag", vel_map),
        ("displacement_mag", disp_map),
        ("linear_kappa", linear_kappa),
    ):
        p, s, n = _corr(field, observed)
        row[f"{name}_vs_observed_pearson"] = p
        row[f"{name}_vs_observed_spearman"] = s
        row[f"{name}_vs_observed_count"] = n
        p2, s2, n2 = _corr(field, los_mag)
        row[f"{name}_vs_los_mag_pearson"] = p2
        row[f"{name}_vs_los_mag_spearman"] = s2
        row[f"{name}_vs_los_mag_count"] = n2
        finite = field[np.isfinite(field)]
        row[f"{name}_rms"] = _rms(finite) if finite.size else float("nan")

    return row, {
        "sampled_response_mag": sampled_map,
        "velocity_change_mag": vel_map,
        "displacement_mag": disp_map,
        "linear_kappa": linear_kappa,
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
    groups = _source_bins(x0, y0)
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
    parity_pass = bool(parity_max <= PARITY_TOL)
    if not parity_pass:
        raise RuntimeError(f"{cid}: checkpoint propagator parity failed max={parity_max}")

    rows = []
    fields = {}
    for k in CHECKPOINTS:
        row, f = _metrics_for_checkpoint(
            cid, k, checkpoints[k], x0, y0, vx0, vy0, groups, los_mag, real["observed_kappa"]
        )
        rows.append(row)
        fields[f"step_{k}"] = f

    finite_linear = [r for r in rows if np.isfinite(r["linear_kappa_vs_observed_pearson"])]
    return {
        "rows": rows,
        "fields": fields,
        "summary": {
            "cluster_id": cid,
            "n_photons": int(len(x0)),
            "checkpoint_count": len(rows),
            "production_parity_rms": parity_rms,
            "production_parity_max": parity_max,
            "production_parity_pass": parity_pass,
            "source_supported_bins": int(len(groups)),
            "los_mag_vs_observed_pearson": _corr(los_mag, real["observed_kappa"])[0],
            "los_mag_vs_observed_spearman": _corr(los_mag, real["observed_kappa"])[1],
            "first_finite_linear_step": int(finite_linear[0]["step_index"]) if finite_linear else None,
            "final_linear_vs_observed_pearson": rows[-1]["linear_kappa_vs_observed_pearson"],
            "final_linear_vs_observed_spearman": rows[-1]["linear_kappa_vs_observed_spearman"],
            "final_displacement_mag_vs_observed_pearson": rows[-1]["displacement_mag_vs_observed_pearson"],
            "final_velocity_change_mag_vs_observed_pearson": rows[-1]["velocity_change_mag_vs_observed_pearson"],
        },
    }


def main() -> int:
    started = time.perf_counter()
    OUT.mkdir(parents=True, exist_ok=True)
    repo = _repo_state()
    _write_json(OUT / "repository_state.json", repo)
    if repo["branch"] != "main" or repo["tracked_changes"] or repo["staged_changes"]:
        v = {"lab_id": LAB_ID, "outcome": "REPOSITORY_GATE_FAILURE", "head_sha": repo["head_sha"]}
        _write_json(OUT / "validation.json", v); print(json.dumps(v, indent=2)); return 2

    all_rows = []
    summaries = []
    failures = []
    for cluster in BASE.CLUSTERS:
        cid = cluster["id"]
        print(f"[{cid}] C25 propagation checkpoint localization")
        try:
            result = _run_cluster(cluster)
            all_rows.extend(result["rows"])
            summaries.append(result["summary"])
            cdir = OUT / "clusters" / cid
            cdir.mkdir(parents=True, exist_ok=True)
            _write_json(cdir / "propagation_summary.json", result["summary"])
            _write_csv(cdir / "propagation_checkpoints.csv", result["rows"])
            payload = {}
            for sk, bundle in result["fields"].items():
                for name, arr in bundle.items():
                    payload[f"{sk}_{name}"] = arr
            np.savez_compressed(cdir / "propagation_checkpoint_fields.npz", **payload)
        except Exception as exc:
            failures.append({"cluster_id": cid, "error": repr(exc)})
            raise

    _write_csv(OUT / "propagation_checkpoint_summary.csv", all_rows)
    _write_csv(OUT / "cluster_summary.csv", summaries)
    _write_json(OUT / "cluster_failures.json", failures)

    all_parity = all(bool(x["production_parity_pass"]) for x in summaries) and len(summaries) == len(BASE.CLUSTERS)
    validation = {
        "lab_id": LAB_ID,
        "outcome": "Outcome A — PROPAGATION MORPHOLOGY LOCALIZATION COMPLETE" if all_parity else "Outcome B — PROPAGATION CHECKPOINT PARITY FAILURE",
        "head_sha": repo["head_sha"],
        "candidate_id": BASE.CANDIDATE_ID,
        "physical_source_representation": BASE.PHYSICAL_SOURCE,
        "coverage_lane": "C25",
        "cluster_count_expected": len(BASE.CLUSTERS),
        "cluster_count_completed": len(summaries),
        "checkpoint_steps": list(CHECKPOINTS),
        "production_parity_tolerance": PARITY_TOL,
        "all_cluster_production_parity_pass": all_parity,
        "physics_change_authorized": False,
        "candidate_change_authorized": False,
        "propagation_change_authorized": False,
        "next_experiment_authorized": False,
        "science_interpretation_required": True,
        "duration_seconds": time.perf_counter() - started,
    }
    _write_json(OUT / "validation.json", validation)
    _write_json(OUT / "run.json", {
        "lab_id": LAB_ID, "head_sha": repo["head_sha"], "checkpoint_steps": list(CHECKPOINTS),
        "config": BASE.CFG, "duration_seconds": validation["duration_seconds"],
    })
    print(json.dumps(validation, indent=2))
    return 0 if all_parity else 1


if __name__ == "__main__":
    raise SystemExit(main())
