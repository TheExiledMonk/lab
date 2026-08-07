#!/usr/bin/env python3
"""PBUF TRANSPORT-CHIRALITY-001 - local resolution of the transverse orientation."""
from __future__ import annotations

import csv
import hashlib
import json
import math
import time
from dataclasses import dataclass, asdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from constitutive_equations import get_equation


ROOT = Path(__file__).resolve().parent
DEFAULT_OUT = ROOT / "runs" / "transport_chirality001"
BASE = {
    "extent": 8.0,
    "mass_x": -0.65,
    "mass_y": 0.0,
    "mass_sigma": 0.75,
    "strength": 0.18,
    "step": 0.06,
    "steps": 80,
    "y_span": 3.0,
    "n": 128,
}
CANDIDATES = [
    ("Cand 1", "Global +90 (control)",        "global_plus"),
    ("Cand 2", "Global -90",                  "global_minus"),
    ("Cand 3", "Local geometry selection",    "local_laplacian"),
    ("Cand 4", "Symmetric dual (diagnostic)", "symmetric_dual"),
    ("Cand 5", "Local centre-seeking rule",   "local_centre"),
]


def sha(array: np.ndarray) -> str:
    return hashlib.sha256(np.ascontiguousarray(array).tobytes()).hexdigest()


def make_field(n: int, extent: float, center: np.ndarray, mass_center: np.ndarray,
               mass_sigma: float, strength: float, dtype):
    x = np.linspace(center[0] - extent, center[0] + extent, n, dtype=dtype)
    y = np.linspace(center[1] - extent, center[1] + extent, n, dtype=dtype)
    X, Y = np.meshgrid(x, y, indexing="xy")
    rho = np.exp(-((X - mass_center[0]) ** 2 + (Y - mass_center[1]) ** 2) /
                 (2 * mass_sigma ** 2)).astype(dtype)
    cfg = type("Config", (), {"deformation_strength": strength})()
    c = get_equation("A").solve(rho, cfg)
    gy, gx = np.gradient(c, x, y, edge_order=1)
    spacing = np.asarray(x[1] - x[0], dtype=dtype)
    gyy, gyx = np.gradient(gy, spacing, spacing, edge_order=1)
    gxy, gxx = np.gradient(gx, spacing, spacing, edge_order=1)
    laplacian = gxx + gyy
    g = np.hypot(gx, gy)
    gx_hat = gx / np.maximum(g, np.asarray(1e-15, dtype=dtype))
    gy_hat = gy / np.maximum(g, np.asarray(1e-15, dtype=dtype))
    bad = g < np.asarray(1e-15, dtype=dtype)
    gx_hat = np.where(bad, np.asarray(1.0, dtype=dtype), gx_hat)
    gy_hat = np.where(bad, np.asarray(0.0, dtype=dtype), gy_hat)
    rx_plus = -g * gy_hat
    ry_plus = g * gx_hat
    rx_minus = g * gy_hat
    ry_minus = -g * gx_hat
    return {"xgrid": x, "ygrid": y, "rho": rho, "c": c, "gx": gx, "gy": gy,
            "laplacian": laplacian, "g_magnitude": g, "rx_plus": rx_plus,
            "ry_plus": ry_plus, "rx_minus": rx_minus, "ry_minus": ry_minus,
            "mass_center": np.asarray(mass_center, dtype=dtype)}


def direct_addition(vx, vy, rx, ry, step, dtype):
    vx_new = vx + np.asarray(step, dtype=dtype) * rx
    vy_new = vy + np.asarray(step, dtype=dtype) * ry
    scale = np.maximum(np.hypot(vx_new, vy_new), np.asarray(1e-12, dtype=dtype))
    return vx_new / scale, vy_new / scale, scale


def select_sign_plus(rx_plus, ry_plus, rx_minus, ry_minus, laplacian, mass_center, x, y, dtype):
    return rx_plus, ry_plus


def select_sign_minus(rx_plus, ry_plus, rx_minus, ry_minus, laplacian, mass_center, x, y, dtype):
    return rx_minus, ry_minus


def select_sign_local_laplacian(rx_plus, ry_plus, rx_minus, ry_minus, laplacian, mass_center, x, y, dtype):
    return np.where(laplacian >= np.asarray(0.0, dtype=dtype), rx_plus, rx_minus), \
           np.where(laplacian >= np.asarray(0.0, dtype=dtype), ry_plus, ry_minus)


def select_sign_centre(rx_plus, ry_plus, rx_minus, ry_minus, laplacian, mass_center, x, y, dtype):
    mc_x = np.asarray(mass_center[0], dtype=dtype)
    mc_y = np.asarray(mass_center[1], dtype=dtype)
    dx = mc_x - x
    dy = mc_y - y
    dot_plus = rx_plus * dx + ry_plus * dy
    return np.where(dot_plus >= np.asarray(0.0, dtype=dtype), rx_plus, rx_minus), \
           np.where(dot_plus >= np.asarray(0.0, dtype=dtype), ry_plus, ry_minus)


SELECTORS = {
    "global_plus": select_sign_plus,
    "global_minus": select_sign_minus,
    "local_laplacian": select_sign_local_laplacian,
    "local_centre": select_sign_centre,
}


def propagate(field, selector, step, steps, x0, y0, vx0, vy0, record):
    xgrid = field["xgrid"]
    ygrid = field["ygrid"]
    rx_plus = field["rx_plus"]
    ry_plus = field["ry_plus"]
    rx_minus = field["rx_minus"]
    ry_minus = field["ry_minus"]
    laplacian = field["laplacian"]
    mass_center = field["mass_center"]
    dtype = np.result_type(xgrid.dtype, ygrid.dtype, rx_plus.dtype,
                           x0.dtype, y0.dtype, vx0.dtype)
    x = np.asarray(x0, dtype=dtype).copy()
    y = np.asarray(y0, dtype=dtype).copy()
    vx = np.asarray(vx0, dtype=dtype).copy()
    vy = np.asarray(vy0, dtype=dtype).copy()
    nphotons = len(x)
    max_deviation = np.zeros(nphotons, dtype=dtype)
    bending_angle = np.zeros(nphotons, dtype=dtype)
    conservation = np.zeros(nphotons, dtype=dtype)
    xs = np.empty((nphotons, steps), dtype=dtype) if record else None
    ys = np.empty((nphotons, steps), dtype=dtype) if record else None
    if record:
        xs[:, 0] = x; ys[:, 0] = y
    started = time.perf_counter()
    for k in range(1, steps):
        ix = np.clip(np.searchsorted(xgrid, x) - 1, 0, len(xgrid) - 1)
        iy = np.clip(np.searchsorted(ygrid, y) - 1, 0, len(ygrid) - 1)
        rx_p = rx_plus[iy, ix]
        ry_p = ry_plus[iy, ix]
        rx_m = rx_minus[iy, ix]
        ry_m = ry_minus[iy, ix]
        lap = laplacian[iy, ix]
        rx, ry = selector(rx_p, ry_p, rx_m, ry_m, lap, mass_center, x, y, dtype)
        vx_new, vy_new, scale = direct_addition(vx, vy, rx, ry, step, dtype)
        conservation = np.maximum(conservation, np.abs(np.hypot(vx_new, vy_new) - 1))
        dot = np.clip(vx * vx_new + vy * vy_new, -1, 1)
        bending_angle += np.arccos(dot)
        vx = vx_new
        vy = vy_new
        x = x + np.asarray(step, dtype=dtype) * vx
        y = y + np.asarray(step, dtype=dtype) * vy
        max_deviation = np.maximum(max_deviation, np.abs(y - np.asarray(y0, dtype=dtype)))
        if record:
            xs[:, k] = x; ys[:, k] = y
    runtime = time.perf_counter() - started
    return {"x": x, "y": y, "max_deviation": max_deviation,
            "bending_angle": bending_angle, "conservation": conservation,
            "runtime": runtime, "xs": xs, "ys": ys}


def propagate_dual(field, step, steps, x0, y0, vx0, vy0, record):
    plus = propagate(field, SELECTORS["global_plus"], step, steps, x0, y0, vx0, vy0, record)
    minus = propagate(field, SELECTORS["global_minus"], step, steps, x0, y0, vx0, vy0, record)
    dtype = plus["x"].dtype
    combined_x = (plus["x"] + minus["x"]) / np.asarray(2.0, dtype=dtype)
    combined_y = (plus["y"] + minus["y"]) / np.asarray(2.0, dtype=dtype)
    combined_max = np.maximum(plus["max_deviation"], minus["max_deviation"])
    combined_bend = np.maximum(plus["bending_angle"], minus["bending_angle"])
    combined_cons = np.maximum(plus["conservation"], minus["conservation"])
    runtime = plus["runtime"] + minus["runtime"]
    combined_xs = None
    combined_ys = None
    if record:
        combined_xs = (plus["xs"] + minus["xs"]) / np.asarray(2.0, dtype=dtype)
        combined_ys = (plus["ys"] + minus["ys"]) / np.asarray(2.0, dtype=dtype)
    return {"x": combined_x, "y": combined_y, "max_deviation": combined_max,
            "bending_angle": combined_bend, "conservation": combined_cons,
            "runtime": runtime, "xs": combined_xs, "ys": combined_ys,
            "plus": plus, "minus": minus}


def run_case(candidate_id, nphotons=9, record=True, n=BASE["n"],
             extent=BASE["extent"], launch_x=None, mass_center=None,
             ray_transform=None, translation=None, dtype=np.float64):
    if mass_center is None:
        mass_center = np.array([BASE["mass_x"], BASE["mass_y"]], dtype=dtype)
    else:
        mass_center = np.asarray(mass_center, dtype=dtype)
    if ray_transform is None:
        ray_transform = np.eye(2, dtype=dtype)
    if translation is None:
        translation = np.zeros(2, dtype=dtype)
    if launch_x is None:
        launch_x = -extent
    field = make_field(n, extent, np.zeros(2, dtype=dtype), mass_center,
                       BASE["mass_sigma"], BASE["strength"], dtype)
    x0 = np.full(nphotons, launch_x, dtype=dtype)
    y0 = np.linspace(-BASE["y_span"], BASE["y_span"], nphotons, dtype=dtype)
    vx0 = np.ones(nphotons, dtype=dtype)
    vy0 = np.zeros(nphotons, dtype=dtype)
    initial = np.column_stack((x0, y0)) @ ray_transform.T + translation
    velocity = np.column_stack((vx0, vy0)) @ ray_transform.T
    if candidate_id == "symmetric_dual":
        result = propagate_dual(field, BASE["step"], BASE["steps"],
                                initial[:, 0], initial[:, 1],
                                velocity[:, 0], velocity[:, 1], record)
    else:
        selector = SELECTORS[candidate_id]
        result = propagate(field, selector, BASE["step"], BASE["steps"],
                           initial[:, 0], initial[:, 1],
                           velocity[:, 0], velocity[:, 1], record)
    result["candidate_id"] = candidate_id
    result.update({"field": field, "x0": x0, "y0": y0, "initial": initial,
                   "velocity": velocity, "transform": ray_transform,
                   "translation": translation, "candidate_id": candidate_id,
                   "dtype": np.dtype(dtype).name})
    return result


def inverse_transform_paths(result):
    transform = np.asarray(result["transform"])
    translation = np.asarray(result["translation"])
    points = np.stack((result["xs"], result["ys"]), axis=-1)
    points = (points - translation) @ transform
    return points[..., 0], points[..., 1]


def summary(result):
    return {
        "bend_max": float(np.max(result["max_deviation"])),
        "bend_mean": float(np.mean(result["max_deviation"])),
        "bending_angle_max": float(np.max(result["bending_angle"])),
        "conservation": float(np.max(result["conservation"])),
        "finite": bool(np.isfinite(result["x"]).all() and np.isfinite(result["y"]).all()),
        "runtime": float(result["runtime"]),
    }


def max_path_difference(a, b):
    return float(max(np.max(np.abs(a["xs"] - b["xs"])), np.max(np.abs(a["ys"] - b["ys"]))))


def symmetry_pass(delta, bend_max):
    relative = delta / max(bend_max, 1e-15)
    return delta <= 1e-7 or relative <= 0.15


def validate_candidate(candidate_id, out_dir):
    base = run_case(candidate_id, record=True)
    repeat = run_case(candidate_id, record=True)
    repeat_delta = max_path_difference(base, repeat)
    mirror = run_case(candidate_id, record=True, ray_transform=np.diag([1.0, -1.0]))
    mx, my = inverse_transform_paths(mirror)
    mirror_delta = float(max(np.max(np.abs(base["xs"] - mx)), np.max(np.abs(base["ys"] - my))))
    rotation = run_case(candidate_id, record=True, ray_transform=matrix_rotate(np.pi / 2))
    rx, ry = inverse_transform_paths(rotation)
    rotation_delta = float(max(np.max(np.abs(base["xs"] - rx)), np.max(np.abs(base["ys"] - ry))))
    translation = run_case(candidate_id, record=True, translation=np.array([0.5, 0.4]),
                           mass_center=np.array([0.5, 0.4]) + np.array([BASE["mass_x"], BASE["mass_y"]]))
    tx, ty = inverse_transform_paths(translation)
    translation_delta = float(max(np.max(np.abs(base["xs"] - tx)), np.max(np.abs(base["ys"] - ty))))
    base_s = summary(base)
    mirror_s = summary(mirror)
    bend_diff = abs(base_s["bend_max"] - mirror_s["bend_max"])
    return {
        "candidate_id": candidate_id,
        "base": base,
        "repeat": repeat,
        "mirror": mirror,
        "tests": {
            "repeatability": {
                "delta": repeat_delta,
                "passed": repeat_delta <= 1e-14,
            },
            "translation": {
                "delta": translation_delta,
                "relative": translation_delta / max(base_s["bend_max"], 1e-15),
                "passed": symmetry_pass(translation_delta, base_s["bend_max"]),
            },
            "rotation": {
                "delta": rotation_delta,
                "relative": rotation_delta / max(base_s["bend_max"], 1e-15),
                "passed": symmetry_pass(rotation_delta, base_s["bend_max"]),
            },
            "mirror": {
                "delta": mirror_delta,
                "relative": mirror_delta / max(base_s["bend_max"], 1e-15),
                "passed": symmetry_pass(mirror_delta, base_s["bend_max"]),
            },
        },
        "measurements": {
            "bend_max": base_s["bend_max"],
            "bending_angle_max": base_s["bending_angle_max"],
            "conservation": base_s["conservation"],
            "mirror_bend_max": mirror_s["bend_max"],
            "bending_diff": bend_diff,
            "runtime": base_s["runtime"],
            "finite": base_s["finite"],
        },
    }


def matrix_rotate(theta):
    c, s = math.cos(theta), math.sin(theta)
    return np.array([[c, -s], [s, c]], dtype=np.float64)


def write_csv(path, rows):
    if not rows:
        return
    keys = []
    for row in rows:
        for key in row:
            if key not in keys:
                keys.append(key)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)


def overlay_plot(out_dir, candidate_results):
    fig, axes = plt.subplots(len(candidate_results), 3, figsize=(15, 4 * len(candidate_results)),
                             squeeze=False)
    for i, (label, _desc, candidate_id) in enumerate(CANDIDATES):
        result = candidate_results[candidate_id]
        base = result["base"]
        mirror = result["mirror"]
        mx, my = inverse_transform_paths(mirror)
        axes[i, 0].set_title(f"{label}: baseline")
        for xs, ys in zip(base["xs"], base["ys"]):
            axes[i, 0].plot(xs, ys, lw=0.8)
        axes[i, 0].set(xlabel="x", ylabel="y", aspect="equal")
        axes[i, 0].set_xlim(-8, 8); axes[i, 0].set_ylim(-8, 8)
        axes[i, 1].set_title(f"{label}: mirrored")
        for xs, ys in zip(mx, my):
            axes[i, 1].plot(xs, ys, lw=0.8)
        axes[i, 1].set(xlabel="x", ylabel="y", aspect="equal")
        axes[i, 1].set_xlim(-8, 8); axes[i, 1].set_ylim(-8, 8)
        diff = np.hypot(base["xs"] - mx, base["ys"] - my)
        im = axes[i, 2].imshow(diff, origin="lower", aspect="auto",
                                extent=[0, BASE["steps"], -BASE["y_span"], BASE["y_span"]])
        axes[i, 2].set_title(f"{label}: overlay diff (max={diff.max():.2e})")
        axes[i, 2].set(xlabel="step", ylabel="y0")
        fig.colorbar(im, ax=axes[i, 2])
    fig.tight_layout()
    fig.savefig(out_dir / "trajectory_overlay.png", dpi=130)
    plt.close(fig)


def main():
    parser_args = __import__("argparse").ArgumentParser()
    parser_args.add_argument("--output", type=Path, default=DEFAULT_OUT)
    args = parser_args.parse_args()
    out_dir = args.output
    out_dir.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()
    candidate_results = {}
    comparison_rows = []
    symmetry_rows = []
    for label, desc, candidate_id in CANDIDATES:
        print(f"Running {label} ({desc}, {candidate_id})...")
        result = validate_candidate(candidate_id, out_dir)
        candidate_results[candidate_id] = result
        m = result["measurements"]
        t = result["tests"]
        comparison_rows.append({
            "candidate": label,
            "candidate_id": candidate_id,
            "bend_max": m["bend_max"],
            "bending_angle_max": m["bending_angle_max"],
            "conservation": m["conservation"],
            "mirror_bend_max": m["mirror_bend_max"],
            "bending_diff": m["bending_diff"],
            "runtime": m["runtime"],
            "mirror_delta": t["mirror"]["delta"],
            "mirror_relative": t["mirror"]["relative"],
            "mirror_passed": t["mirror"]["passed"],
            "translation_delta": t["translation"]["delta"],
            "rotation_delta": t["rotation"]["delta"],
            "repeatability_delta": t["repeatability"]["delta"],
        })
        for test_name, test_data in t.items():
            symmetry_rows.append({
                "candidate": label,
                "candidate_id": candidate_id,
                "test": test_name,
                "delta": test_data["delta"],
                "passed": test_data["passed"],
            })
    write_csv(out_dir / "comparison_table.csv", comparison_rows)
    write_csv(out_dir / "symmetry_summary.csv", symmetry_rows)
    overlay_plot(out_dir, candidate_results)
    summary_payload = {
        "milestone": "PBUF TRANSPORT-CHIRALITY-001",
        "frozen_conditions": {
            "constitutive": "u = 0.18 rho/rho_max (Version A)",
            "transport_magnitude": "|grad C|",
            "response_angle": "90 degrees",
            "update_rule": "direct addition + renormalisation (transport_lab007.upd_direct_addition)",
            "timestep": BASE["step"],
            "normalisation": "max A preserved",
            "lens": "Lens-001 dataset",
        },
        "candidates": [
            {"id": c[1], "label": c[0], "rule": c[1]} for c in CANDIDATES
        ],
        "comparison": comparison_rows,
        "symmetry": symmetry_rows,
        "execution_seconds": time.perf_counter() - started,
    }
    (out_dir / "validation.json").write_text(json.dumps(summary_payload, indent=2))
    (out_dir / "run.json").write_text(json.dumps({
        "milestone": "PBUF TRANSPORT-CHIRALITY-001",
        "status": "OK",
        "execution_seconds": time.perf_counter() - started,
    }, indent=2))
    mirror_pass_per_candidate = {row["candidate_id"]: row["mirror_passed"] for row in comparison_rows}
    outcome_a = any(passed for cid, passed in mirror_pass_per_candidate.items() if cid != "symmetric_dual")
    print(json.dumps({
        "milestone": "PBUF TRANSPORT-CHIRALITY-001",
        "mirror_pass_per_candidate": mirror_pass_per_candidate,
        "outcome_a_local_restores_mirror": outcome_a,
        "output": str(out_dir),
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
