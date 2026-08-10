#!/usr/bin/env python3
"""PBUF FOUNDATION — LOS-CONSISTENT RAY GEOMETRY 001.

Paired diagnostic after propagation-coupling geometry audit 001.

The PL1_PM1_PS2 / M10 physical field, C25 source plane, photon surface density,
step size, step count, and Jacobian observable are frozen.  Two geometry lanes
are compared without fitting or selection:

  G2D — current production geometry.  The image-plane x coordinate is also the
        forward propagation coordinate and photons launch with v=(1,0).

  G3D — LOS-consistent diagnostic geometry.  The two image-plane coordinates
        remain transverse coordinates (x,y), the photon launches along an
        independent LOS coordinate z with v=(0,0,1), and the M10 LOS response
        enters as R=(Rx,Ry,0).

G3D uses the same normalized unit-velocity update law as production, generalized
only by the explicit LOS dimension.  It is a diagnostic lane, not a replacement
for production.  Correlation with observed kappa is measurement only and is
never an execution gate.
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
import observable_lab001 as obs_lab
import pbuf.labs.foundation.m10_coverage_25pct_science001 as BASE
from pbuf.core import los_projection as M14
from pbuf.core import observable_extraction as M16

LAB_ID = "PBUF-FOUNDATION-LOS-CONSISTENT-RAY-GEOMETRY-001"
OUT = ROOT / "runs" / "los_consistent_ray_geometry001"
CHECKPOINTS = (0, 1, 5, 10, 20, 40, 80, 120, 159)
FIRST_STEP_TOL = 1e-14
UNIT_SPEED_TOL = 1e-12
EXPECTED_SUPPORT = 1024


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


def _sample(field, x, y):
    xgrid = field["xgrid"]
    ygrid = field["ygrid"]
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
    out = np.full((BASE.OBS_BINS, BASE.OBS_BINS), np.nan)
    values = np.asarray(values, dtype=np.float64)
    for q, idx in groups.items():
        r, c = divmod(q, BASE.OBS_BINS)
        out[r, c] = float(np.mean(values[idx]))
    return out


def _linear_kappa_map(x0, y0, xf, yf, groups):
    out = np.full((BASE.OBS_BINS, BASE.OBS_BINS), np.nan)
    for q, idx in groups.items():
        r, c = divmod(q, BASE.OBS_BINS)
        xi = x0[idx]; yi = y0[idx]
        xo = xf[idx]; yo = yf[idx]
        A = np.column_stack([xi - xi.mean(), yi - yi.mean()])
        Jx, *_ = np.linalg.lstsq(A, xo - xo.mean(), rcond=None)
        Jy, *_ = np.linalg.lstsq(A, yo - yo.mean(), rcond=None)
        J = np.array([[Jx[0], Jx[1]], [Jy[0], Jy[1]]], dtype=np.float64)
        out[r, c] = -float(np.trace(J - np.eye(2)))
    return out


def _observer_array(value):
    """Return an isolated read-only observer payload.

    Observer callbacks are outside the frozen numerical path.  Copies prevent
    even a deliberately hostile callback from mutating propagator state.
    """
    payload = np.asarray(value, dtype=np.float64).copy()
    payload.flags.writeable = False
    return payload


def _propagate_g3d(field, step, steps, x0, y0, step_observer=None):
    """Same normalized response update as production with explicit LOS z."""
    wanted = set(CHECKPOINTS)
    if max(wanted) >= steps:
        raise RuntimeError(f"checkpoint >= steps: {max(wanted)} >= {steps}")

    x = np.asarray(x0, dtype=np.float64).copy()
    y = np.asarray(y0, dtype=np.float64).copy()
    z = np.zeros_like(x)
    vx = np.zeros_like(x)
    vy = np.zeros_like(x)
    vz = np.ones_like(x)
    checkpoints = {}
    max_unit_error = 0.0

    if step_observer is not None:
        launch_rx, launch_ry = _sample(field, x, y)
        step_observer.observe_launch(
            ray_index=np.arange(x.size, dtype=np.int64),
            position=_observer_array(np.column_stack((x, y, z))),
            direction=_observer_array(np.column_stack((vx, vy, vz))),
            launch_coordinates=_observer_array(np.column_stack((x0, y0))),
            native_state={"rx_sample": _observer_array(launch_rx),
                          "ry_sample": _observer_array(launch_ry)},
        )

    def capture(k):
        rx, ry = _sample(field, x, y)
        checkpoints[int(k)] = {
            "x": x.copy(), "y": y.copy(), "z": z.copy(),
            "vx": vx.copy(), "vy": vy.copy(), "vz": vz.copy(),
            "rx_sample": rx.copy(), "ry_sample": ry.copy(),
        }

    if 0 in wanted:
        capture(0)

    for k in range(1, steps):
        rx, ry = _sample(field, x, y)
        vx_raw = vx + step * rx
        vy_raw = vy + step * ry
        vz_raw = vz
        scale = np.maximum(np.sqrt(vx_raw*vx_raw + vy_raw*vy_raw + vz_raw*vz_raw), 1e-12)
        vx = vx_raw / scale
        vy = vy_raw / scale
        vz = vz_raw / scale
        unit_error = np.max(np.abs(np.sqrt(vx*vx + vy*vy + vz*vz) - 1.0))
        max_unit_error = max(max_unit_error, float(unit_error))
        x = x + step * vx
        y = y + step * vy
        z = z + step * vz
        if step_observer is not None:
            # POST_STEP convention.  rx/ry are the native values actually used
            # by this update; no duplicate interpolation is performed.
            step_observer.observe_step(
                ray_index=np.arange(x.size, dtype=np.int64), step_index=int(k),
                position=_observer_array(np.column_stack((x, y, z))),
                direction=_observer_array(np.column_stack((vx, vy, vz))),
                native_state={"rx_sample": _observer_array(rx),
                              "ry_sample": _observer_array(ry)},
                ds=float(step),
            )
        if k in wanted:
            capture(k)

    if step_observer is not None:
        step_observer.observe_termination(
            ray_index=np.arange(x.size, dtype=np.int64),
            termination_status="completed", final_step_index=int(steps - 1))

    return checkpoints, {
        "x": x, "y": y, "z": z,
        "vx": vx, "vy": vy, "vz": vz,
        "max_unit_speed_error": max_unit_error,
    }


def _first_step_geometry(field, x0, y0, snap1, observed, los_mag):
    h = float(BASE.CFG["step"])
    rx0, ry0 = _sample(field, x0, y0)
    norm = np.sqrt(1.0 + (h*rx0)**2 + (h*ry0)**2)
    dx_exact = h*h*rx0 / norm
    dy_exact = h*h*ry0 / norm
    dz_exact = h / norm
    dx_actual = snap1["x"] - x0
    dy_actual = snap1["y"] - y0
    dz_actual = snap1["z"]
    err = np.sqrt((dx_actual-dx_exact)**2 + (dy_actual-dy_exact)**2 + (dz_actual-dz_exact)**2)

    groups = _source_groups(x0, y0)
    trans_actual = np.hypot(dx_actual, dy_actual)
    trans_map = _mean_map(trans_actual, groups)
    p_obs, s_obs, n_obs = _corr(trans_map, observed)
    p_los, s_los, n_los = _corr(trans_map, los_mag)

    return {
        "first_step_exact_max_vector_error": float(np.max(err)),
        "first_step_exact_rms_vector_error": _rms(err),
        "first_step_exact_pass": bool(float(np.max(err)) <= FIRST_STEP_TOL),
        "first_step_transverse_mag_vs_observed_pearson": p_obs,
        "first_step_transverse_mag_vs_observed_spearman": s_obs,
        "first_step_transverse_mag_vs_observed_count": n_obs,
        "first_step_transverse_mag_vs_los_mag_pearson": p_los,
        "first_step_transverse_mag_vs_los_mag_spearman": s_los,
        "first_step_transverse_mag_vs_los_mag_count": n_los,
        "first_step_dx_vs_Rx_pearson": _corr(dx_actual, rx0)[0],
        "first_step_dy_vs_Ry_pearson": _corr(dy_actual, ry0)[0],
    }


def _checkpoint_metrics(cid, k, snap, x0, y0, groups, observed, los_mag):
    sampled_mag = np.hypot(snap["rx_sample"], snap["ry_sample"])
    trans_v_mag = np.hypot(snap["vx"], snap["vy"])
    trans_disp_mag = np.hypot(snap["x"]-x0, snap["y"]-y0)
    sampled_map = _mean_map(sampled_mag, groups)
    vel_map = _mean_map(trans_v_mag, groups)
    disp_map = _mean_map(trans_disp_mag, groups)
    linear = _linear_kappa_map(x0, y0, snap["x"], snap["y"], groups)

    row = {
        "cluster_id": cid,
        "step_index": int(k),
        "propagation_distance": float(k * BASE.CFG["step"]),
    }
    for name, arr in (
        ("sampled_response_mag", sampled_map),
        ("transverse_velocity_mag", vel_map),
        ("transverse_displacement_mag", disp_map),
        ("linear_kappa", linear),
    ):
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
    return row, {
        "sampled_response_mag": sampled_map,
        "transverse_velocity_mag": vel_map,
        "transverse_displacement_mag": disp_map,
        "linear_kappa": linear,
    }


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
    observed = real["observed_kappa"]
    grid = np.linspace(-BASE.CFG["extent"], BASE.CFG["extent"], Rx.shape[0])
    field = {"xgrid": grid, "ygrid": grid, "rx": Rx, "ry": Ry}

    x0, y0, _, _ = BASE._launch_expanded_25pct()
    groups = _source_groups(x0, y0)
    if len(groups) != EXPECTED_SUPPORT:
        raise RuntimeError(f"{cid}: expected {EXPECTED_SUPPORT} source bins, got {len(groups)}")

    # G2D: exact production geometry, frozen.
    vx0_2d = np.ones_like(x0)
    vy0_2d = np.zeros_like(y0)
    g2d_ph = production_propagate(
        field, BASE.CFG["step"], BASE.CFG["steps"], x0, y0, vx0_2d, vy0_2d
    )
    g2d_jac = obs_lab.method_jacobian(
        x0, y0, g2d_ph["x"], g2d_ph["y"], BASE.CFG["extent"], BASE.OBS_BINS
    )
    g2d_kappa = np.asarray(g2d_jac["convergence"], dtype=np.float64)

    # G3D: same law, explicit LOS direction and two transverse coordinates.
    checkpoints, g3d = _propagate_g3d(
        field, BASE.CFG["step"], BASE.CFG["steps"], x0, y0
    )
    if g3d["max_unit_speed_error"] > UNIT_SPEED_TOL:
        raise RuntimeError(f"{cid}: G3D unit-speed gate failed: {g3d['max_unit_speed_error']}")
    first = _first_step_geometry(field, x0, y0, checkpoints[1], observed, los_mag)
    if not first["first_step_exact_pass"]:
        raise RuntimeError(f"{cid}: G3D first-step exact geometry gate failed")

    g3d_jac = obs_lab.method_jacobian(
        x0, y0, g3d["x"], g3d["y"], BASE.CFG["extent"], BASE.OBS_BINS
    )
    g3d_kappa = np.asarray(g3d_jac["convergence"], dtype=np.float64)
    mask2 = np.isfinite(g2d_kappa)
    mask3 = np.isfinite(g3d_kappa)
    if int(mask3.sum()) != EXPECTED_SUPPORT:
        raise RuntimeError(f"{cid}: G3D finite kappa count {int(mask3.sum())} != {EXPECTED_SUPPORT}")

    checkpoint_rows = []
    checkpoint_fields = {}
    for k in CHECKPOINTS:
        row, fields = _checkpoint_metrics(cid, k, checkpoints[k], x0, y0, groups, observed, los_mag)
        checkpoint_rows.append(row)
        checkpoint_fields[f"step_{k}"] = fields

    g3d_linear_final = checkpoint_fields[f"step_{CHECKPOINTS[-1]}"]["linear_kappa"]
    trans_disp = np.hypot(g3d["x"]-x0, g3d["y"]-y0)
    trans_disp_map = _mean_map(trans_disp, groups)
    trans_v_map = _mean_map(np.hypot(g3d["vx"], g3d["vy"]), groups)

    common = mask2 & mask3
    summary = {
        "cluster_id": cid,
        "candidate_id": BASE.CANDIDATE_ID,
        "physical_source_representation": BASE.PHYSICAL_SOURCE,
        "n_photons": int(len(x0)),
        "source_supported_bins": int(len(groups)),
        "g2d_finite_kappa_count": int(mask2.sum()),
        "g3d_finite_kappa_count": int(mask3.sum()),
        "g3d_unit_speed_max_error": float(g3d["max_unit_speed_error"]),
        "g3d_unit_speed_pass": bool(g3d["max_unit_speed_error"] <= UNIT_SPEED_TOL),
        **first,
        "los_mag_vs_observed_pearson": _corr(los_mag, observed)[0],
        "los_mag_vs_observed_spearman": _corr(los_mag, observed)[1],
        "g2d_kappa_vs_observed_pearson": _corr(g2d_kappa, observed)[0],
        "g2d_kappa_vs_observed_spearman": _corr(g2d_kappa, observed)[1],
        "g3d_kappa_vs_observed_pearson": _corr(g3d_kappa, observed)[0],
        "g3d_kappa_vs_observed_spearman": _corr(g3d_kappa, observed)[1],
        "g3d_linear_kappa_vs_observed_pearson": _corr(g3d_linear_final, observed)[0],
        "g3d_linear_kappa_vs_observed_spearman": _corr(g3d_linear_final, observed)[1],
        "g3d_exact_vs_linear_kappa_pearson": _corr(g3d_kappa, g3d_linear_final)[0],
        "g3d_exact_vs_linear_kappa_spearman": _corr(g3d_kappa, g3d_linear_final)[1],
        "g2d_vs_g3d_kappa_commonmask_count": int(common.sum()),
        "g2d_vs_g3d_kappa_commonmask_pearson": _corr(g2d_kappa[common], g3d_kappa[common])[0],
        "g2d_vs_g3d_kappa_commonmask_spearman": _corr(g2d_kappa[common], g3d_kappa[common])[1],
        "g3d_transverse_displacement_mag_vs_observed_pearson": _corr(trans_disp_map, observed)[0],
        "g3d_transverse_displacement_mag_vs_observed_spearman": _corr(trans_disp_map, observed)[1],
        "g3d_transverse_velocity_mag_vs_observed_pearson": _corr(trans_v_map, observed)[0],
        "g3d_transverse_velocity_mag_vs_observed_spearman": _corr(trans_v_map, observed)[1],
        "g2d_kappa_rms": _rms(g2d_kappa[mask2]),
        "g3d_kappa_rms": _rms(g3d_kappa[mask3]),
        "g3d_final_mean_z": float(np.mean(g3d["z"])),
        "g3d_final_z_rms_spread": float(np.std(g3d["z"])),
    }
    return summary, checkpoint_rows, checkpoint_fields, {
        "g2d_kappa": g2d_kappa,
        "g3d_kappa": g3d_kappa,
        "los_mag": los_mag,
        "observed_kappa": observed,
        "g3d_x_final": g3d["x"],
        "g3d_y_final": g3d["y"],
        "g3d_z_final": g3d["z"],
    }


def main() -> int:
    started = time.perf_counter()
    OUT.mkdir(parents=True, exist_ok=True)
    repo = _repo_state()
    _write_json(OUT / "repository_state.json", repo)
    if repo["branch"] != "main" or repo["tracked_changes"] or repo["staged_changes"]:
        validation = {"lab_id": LAB_ID, "outcome": "REPOSITORY_GATE_FAILURE", "head_sha": repo["head_sha"]}
        _write_json(OUT / "validation.json", validation)
        print(json.dumps(validation, indent=2))
        return 2

    summaries = []
    all_checkpoints = []
    failures = []
    for cluster in BASE.CLUSTERS:
        cid = cluster["id"]
        print(f"[{cid}] paired G2D production vs G3D LOS-consistent geometry")
        try:
            summary, rows, fields, final_fields = _run_cluster(cluster)
            summaries.append(summary)
            all_checkpoints.extend(rows)
            cdir = OUT / "clusters" / cid
            cdir.mkdir(parents=True, exist_ok=True)
            _write_json(cdir / "geometry_summary.json", summary)
            _write_csv(cdir / "g3d_checkpoints.csv", rows)
            npz = {}
            for step_name, fd in fields.items():
                for name, arr in fd.items():
                    npz[f"{step_name}__{name}"] = arr
            npz.update(final_fields)
            np.savez_compressed(cdir / "geometry_checkpoint_fields.npz", **npz)
        except Exception as exc:
            failures.append({"cluster_id": cid, "error": repr(exc)})
            _write_json(OUT / "cluster_failures.json", failures)
            raise

    _write_csv(OUT / "geometry_summary.csv", summaries)
    _write_csv(OUT / "g3d_checkpoint_summary.csv", all_checkpoints)
    _write_json(OUT / "cluster_failures.json", failures)

    validation = {
        "lab_id": LAB_ID,
        "outcome": "Outcome A — PAIRED 2D/3D LOS RAY GEOMETRY AUDIT COMPLETE",
        "head_sha": repo["head_sha"],
        "candidate_id": BASE.CANDIDATE_ID,
        "physical_source_representation": BASE.PHYSICAL_SOURCE,
        "coverage_lane": "C25",
        "cluster_count_expected": len(BASE.CLUSTERS),
        "cluster_count_completed": len(summaries),
        "checkpoint_steps": list(CHECKPOINTS),
        "first_step_exact_tolerance": FIRST_STEP_TOL,
        "unit_speed_tolerance": UNIT_SPEED_TOL,
        "all_cluster_first_step_exact_pass": bool(all(r["first_step_exact_pass"] for r in summaries)),
        "all_cluster_g3d_unit_speed_pass": bool(all(r["g3d_unit_speed_pass"] for r in summaries)),
        "all_cluster_g3d_support_1024": bool(all(r["g3d_finite_kappa_count"] == EXPECTED_SUPPORT for r in summaries)),
        "mean_los_mag_vs_observed_pearson": float(np.mean([r["los_mag_vs_observed_pearson"] for r in summaries])),
        "mean_g2d_kappa_vs_observed_pearson": float(np.mean([r["g2d_kappa_vs_observed_pearson"] for r in summaries])),
        "mean_g3d_kappa_vs_observed_pearson": float(np.mean([r["g3d_kappa_vs_observed_pearson"] for r in summaries])),
        "mean_g2d_kappa_vs_observed_spearman": float(np.mean([r["g2d_kappa_vs_observed_spearman"] for r in summaries])),
        "mean_g3d_kappa_vs_observed_spearman": float(np.mean([r["g3d_kappa_vs_observed_spearman"] for r in summaries])),
        "mean_first_step_transverse_mag_vs_los_mag_pearson": float(np.mean([r["first_step_transverse_mag_vs_los_mag_pearson"] for r in summaries])),
        "mean_first_step_transverse_mag_vs_observed_pearson": float(np.mean([r["first_step_transverse_mag_vs_observed_pearson"] for r in summaries])),
        "physics_change_authorized": False,
        "candidate_change_authorized": False,
        "production_geometry_change_authorized": False,
        "next_experiment_authorized": False,
        "science_interpretation_required": True,
        "duration_seconds": float(time.perf_counter() - started),
    }
    _write_json(OUT / "validation.json", validation)
    _write_json(OUT / "run.json", {"lab_id": LAB_ID, "validation": validation})
    print(json.dumps(validation, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
