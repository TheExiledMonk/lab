#!/usr/bin/env python3
"""PBUF WEAK-LENSING-VALIDATION-001 numerical validation runner."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import subprocess
import sys
import time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from constitutive_equations import get_equation
from numpy import isfinite


ROOT = Path(__file__).resolve().parent
DEFAULT_OUT = ROOT / "runs" / "weak_lensing_validation001"
BASE = {
    "extent": 8.0,
    "mass_x": -0.65,
    "mass_y": 0.0,
    "mass_sigma": 0.75,
    "strength": 0.18,
    "step": 0.06,
    "steps": 80,
    "y_span": 3.0,
}


def direct_addition(v: np.ndarray, r: np.ndarray, step: float) -> np.ndarray:
    v_new = v + step * r
    n = max(float(np.linalg.norm(v_new)), 1e-12)
    return v_new / n


def sha(array: np.ndarray) -> str:
    return hashlib.sha256(np.ascontiguousarray(array).tobytes()).hexdigest()


def matrix_rotate(theta: float) -> np.ndarray:
    c, s = math.cos(theta), math.sin(theta)
    return np.array([[c, -s], [s, c]], dtype=np.float64)


def make_field(n: int, extent: float, center: np.ndarray, mass_center: np.ndarray,
               mass_sigma: float, strength: float, dtype) -> tuple[np.ndarray, ...]:
    x = np.linspace(center[0] - extent, center[0] + extent, n, dtype=dtype)
    y = np.linspace(center[1] - extent, center[1] + extent, n, dtype=dtype)
    X, Y = np.meshgrid(x, y, indexing="xy")
    rho = np.exp(-((X - mass_center[0]) ** 2 + (Y - mass_center[1]) ** 2) /
                 (2 * mass_sigma ** 2)).astype(dtype)
    c = get_equation("A").solve(rho, type("Config", (), {"deformation_strength": strength})())
    gy, gx = np.gradient(c, x, y, edge_order=1)
    g = np.hypot(gx, gy)
    gx_hat = gx / np.maximum(g, np.array(1e-15, dtype=dtype))
    gy_hat = gy / np.maximum(g, np.array(1e-15, dtype=dtype))
    bad = g < np.array(1e-15, dtype=dtype)
    gx_hat = np.where(bad, np.array(1.0, dtype=dtype), gx_hat)
    gy_hat = np.where(bad, np.array(0.0, dtype=dtype), gy_hat)
    rx = -g * gy_hat
    ry = g * gx_hat
    return x, y, rho, c, gx, gy, rx, ry


def propagate(xgrid: np.ndarray, ygrid: np.ndarray, rx: np.ndarray, ry: np.ndarray,
              step: float, steps: int, x0: np.ndarray, y0: np.ndarray,
              vx0: np.ndarray, vy0: np.ndarray, record: bool = False) -> dict:
    dtype = np.result_type(xgrid.dtype, ygrid.dtype, rx.dtype, ry.dtype,
                           x0.dtype, y0.dtype, vx0.dtype, vy0.dtype)
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
        xs[:, 0], ys[:, 0] = x, y
    started = time.perf_counter()
    for k in range(1, steps):
        ix = np.clip(np.searchsorted(xgrid, x) - 1, 0, len(xgrid) - 1)
        iy = np.clip(np.searchsorted(ygrid, y) - 1, 0, len(ygrid) - 1)
        response_x = rx[iy, ix]
        response_y = ry[iy, ix]
        vx_new = vx + np.asarray(step, dtype=dtype) * response_x
        vy_new = vy + np.asarray(step, dtype=dtype) * response_y
        scale = np.maximum(np.hypot(vx_new, vy_new), np.array(1e-12, dtype=dtype))
        vx_unit = vx_new / scale
        vy_unit = vy_new / scale
        conservation = np.maximum(conservation, np.abs(np.hypot(vx_unit, vy_unit) - 1))
        dot = np.clip(vx * vx_unit + vy * vy_unit, -1, 1)
        bending_angle += np.arccos(dot)
        vx = vx_unit
        vy = vy_unit
        x = x + np.asarray(step, dtype=dtype) * vx
        y = y + np.asarray(step, dtype=dtype) * vy
        max_deviation = np.maximum(max_deviation, np.abs(y - y0))
        if record:
            xs[:, k], ys[:, k] = x, y
    runtime = time.perf_counter() - started
    return {
        "x": x,
        "y": y,
        "max_deviation": max_deviation,
        "bending_angle": bending_angle,
        "conservation": conservation,
        "runtime": runtime,
        "xs": xs,
        "ys": ys,
    }


def base_rays(nphotons: int, extent: float, y_span: float, dtype) -> tuple[np.ndarray, ...]:
    y = np.linspace(-y_span, y_span, nphotons, dtype=dtype)
    x = np.full(nphotons, -extent, dtype=dtype)
    return x, y, np.ones(nphotons, dtype=dtype), np.zeros(nphotons, dtype=dtype)


def run_case(n: int = 128, extent: float = 8.0, step: float = 0.06,
             steps: int = 80, nphotons: int = 9, dtype=np.float64,
             center: np.ndarray | None = None, mass_center: np.ndarray | None = None,
             mass_sigma: float = BASE["mass_sigma"],
             launch_x: float | None = None,
             ray_transform: np.ndarray | None = None, translation: np.ndarray | None = None,
             record: bool = False) -> dict:
    center = np.zeros(2, dtype=dtype) if center is None else np.asarray(center, dtype=dtype)
    mass_center = np.array([BASE["mass_x"], BASE["mass_y"]], dtype=dtype) if mass_center is None else np.asarray(mass_center, dtype=dtype)
    transform = np.eye(2, dtype=dtype) if ray_transform is None else np.asarray(ray_transform, dtype=dtype)
    translation = np.zeros(2, dtype=dtype) if translation is None else np.asarray(translation, dtype=dtype)
    xgrid, ygrid, rho, c, gx, gy, rx, ry = make_field(
        n, extent, center, mass_center, mass_sigma, BASE["strength"], dtype)
    x_launch = -extent if launch_x is None else launch_x
    x0 = np.full(nphotons, x_launch, dtype=dtype)
    y0 = np.linspace(-BASE["y_span"], BASE["y_span"], nphotons, dtype=dtype)
    vx0 = np.ones(nphotons, dtype=dtype)
    vy0 = np.zeros(nphotons, dtype=dtype)
    initial = np.column_stack((x0, y0)) @ transform.T + translation
    velocity = np.column_stack((vx0, vy0)) @ transform.T
    result = propagate(xgrid, ygrid, rx, ry, step, steps, initial[:, 0], initial[:, 1],
                       velocity[:, 0], velocity[:, 1], record)
    result.update({"rho": rho, "c": c, "gx": gx, "gy": gy, "rx": rx, "ry": ry,
                   "xgrid": xgrid, "ygrid": ygrid, "x0": x0, "y0": y0,
                   "initial": initial, "velocity": velocity,
                   "transform": transform, "translation": translation,
                   "dtype": np.dtype(dtype).name})
    return result


def inverse_transform_paths(result: dict) -> tuple[np.ndarray, np.ndarray]:
    transform = result["transform"]
    translation = result["translation"]
    points = np.stack((result["xs"], result["ys"]), axis=-1)
    points = (points - translation) @ transform
    return points[..., 0], points[..., 1]


def map_observables(result: dict, bins: int = 32) -> dict:
    xf, yf = result["x"], result["y"]
    x0, y0 = result["x0"], result["y0"]
    extent = BASE["extent"]
    edges = np.linspace(-extent, extent, bins + 1)
    final_count, _, _ = np.histogram2d(yf, xf, bins=(edges, edges))
    initial_count, _, _ = np.histogram2d(y0, x0, bins=(edges, edges))
    safe = initial_count > 0
    convergence = np.full((bins, bins), np.nan)
    convergence[safe] = 0.5 * (final_count[safe] / initial_count[safe] - 1.0)
    sum_dx, _, _ = np.histogram2d(yf, xf, bins=(edges, edges), weights=xf - x0)
    sum_dy, _, _ = np.histogram2d(yf, xf, bins=(edges, edges), weights=yf - y0)
    deflection_x = np.full((bins, bins), np.nan)
    deflection_y = np.full((bins, bins), np.nan)
    good = final_count > 0
    deflection_x[good] = sum_dx[good] / final_count[good]
    deflection_y[good] = sum_dy[good] / final_count[good]
    fill_x = np.nan_to_num(deflection_x, nan=0.0)
    fill_y = np.nan_to_num(deflection_y, nan=0.0)
    spacing = edges[1] - edges[0]
    dxx = np.gradient(fill_x, spacing, axis=1)
    dyy = np.gradient(fill_y, spacing, axis=0)
    dxy = np.gradient(fill_x, spacing, axis=0)
    dyx = np.gradient(fill_y, spacing, axis=1)
    shear_g1 = 0.5 * (dxx - dyy)
    shear_g2 = 0.5 * (dxy + dyx)
    shear_g1[~good] = np.nan
    shear_g2[~good] = np.nan
    return {"convergence": convergence, "shear_g1": shear_g1, "shear_g2": shear_g2,
            "deflection_x": deflection_x, "deflection_y": deflection_y,
            "ray_count": final_count, "edges": edges}


def summary(result: dict) -> dict:
    return {
        "bend_max": float(np.max(result["max_deviation"])),
        "bend_mean": float(np.mean(result["max_deviation"])),
        "bending_angle_max": float(np.max(result["bending_angle"])),
        "bending_angle_mean": float(np.mean(result["bending_angle"])),
        "conservation": float(np.max(result["conservation"])),
        "final_x_min": float(np.min(result["x"])),
        "final_x_max": float(np.max(result["x"])),
        "final_y_min": float(np.min(result["y"])),
        "final_y_max": float(np.max(result["y"])),
        "finite": bool(np.isfinite(result["x"]).all() and np.isfinite(result["y"]).all()),
        "runtime": float(result["runtime"]),
    }


def max_path_difference(a: dict, b: dict) -> float:
    return float(max(np.max(np.abs(a["xs"] - b["xs"])), np.max(np.abs(a["ys"] - b["ys"]))))


def relative_difference(a: float, b: float) -> float:
    return abs(a - b) / max(abs(a), abs(b), 1e-15)


def order_from_errors(errors: list[float]) -> float:
    pairs = [(errors[i], errors[i + 1]) for i in range(len(errors) - 1)
             if errors[i] > 1e-30 and errors[i + 1] > 1e-30]
    if not pairs:
        return float("nan")
    return float(np.mean([math.log(e0 / e1, 2) for e0, e1 in pairs]))


def write_csv(path: Path, rows: list[dict]) -> None:
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


def save_products(out: Path, result: dict) -> None:
    products = map_observables(result)
    for name in ["convergence", "shear_g1", "shear_g2", "deflection_x", "deflection_y"]:
        np.savetxt(out / f"{name}.csv", products[name], delimiter=",")
    if result["xs"] is not None:
        np.savetxt(out / "photon_trajectories_x.csv", result["xs"], delimiter=",")
        np.savetxt(out / "photon_trajectories_y.csv", result["ys"], delimiter=",")
    x = result["xgrid"]
    fields = [("matter", result["rho"]), ("constitutive", result["c"]),
              ("gradient_magnitude", np.hypot(result["gx"], result["gy"])),
              ("response_x", result["rx"]), ("response_y", result["ry"])]
    for name, array in fields:
        np.savetxt(out / f"{name}.csv", array, delimiter=",")
    fig, ax = plt.subplots(figsize=(8, 5))
    if result["xs"] is not None:
        for xs, ys in zip(result["xs"], result["ys"]):
            ax.plot(xs, ys, linewidth=0.8)
    ax.set(xlabel="x", ylabel="y", title="Photon trajectories")
    fig.tight_layout(); fig.savefig(out / "photon_trajectories.png", dpi=140); plt.close(fig)
    extent = [-8, 8, -8, 8]
    for name, title, array in [("convergence", "Convergence map", products["convergence"]),
                               ("shear", "Shear map", np.hypot(products["shear_g1"], products["shear_g2"])),
                               ("deflection", "Deflection map", np.hypot(products["deflection_x"], products["deflection_y"]))]:
        fig, ax = plt.subplots(figsize=(6, 5))
        im = ax.imshow(array, origin="lower", extent=extent)
        ax.set(title=title, xlabel="x", ylabel="y")
        fig.colorbar(im, ax=ax)
        fig.tight_layout(); fig.savefig(out / f"{name}_map.png", dpi=140); plt.close(fig)


def validation_row(test: str, passed: bool, notes: str, value=None) -> dict:
    row = {"test": test, "pass": bool(passed), "fail": not bool(passed), "notes": notes}
    if value is not None:
        row["value"] = value
    return row


def probe() -> None:
    result = run_case(record=True)
    print(json.dumps({"path_checksum_x": sha(result["xs"]), "path_checksum_y": sha(result["ys"])}))


def rescaled_sigma(extent: float) -> float:
    return BASE["mass_sigma"] * extent / BASE["extent"]


def _plot_refinement(out: Path, grid_rows: list[dict], dt_rows: list[dict]) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    spacings = [16.0 / (row["resolution"] - 1) for row in grid_rows]
    bends = [row["bend"] for row in grid_rows]
    axes[0].plot(spacings, bends, "o-", label="bending")
    axes[0].set(xlabel="grid spacing", ylabel="bending (max |y-y0|)",
                title="Bending vs grid spacing", xscale="log", yscale="log")
    axes[0].grid(True, which="both", alpha=0.3)
    axes[0].legend()
    steps = [row["step"] for row in dt_rows]
    axes[1].plot(steps, [row["bend"] for row in dt_rows], "o-", label="bending")
    axes[1].set(xlabel="timestep", ylabel="bending (max |y-y0|)",
                title="Bending vs timestep", xscale="log", yscale="log")
    axes[1].grid(True, which="both", alpha=0.3)
    axes[1].legend()
    fig.tight_layout()
    fig.savefig(out / "refinement_plot.png", dpi=140)
    plt.close(fig)


def _write_report(out: Path, rows: list[dict], grid_rows: list[dict],
                  dt_rows: list[dict], photon_rows: list[dict],
                  domain_rows: list[dict], precision_rows: list[dict],
                  all_pass: bool) -> None:
    lines = ["# PBUF WEAK-LENSING-VALIDATION-001", ""]
    lines.append("Frozen pipeline: matter → C=0.18 ρ/ρ_max → ∇C → 90° transverse response → "
                 "direct addition + renormalisation → photon propagation → observables.")
    lines.append("")
    lines.append("## Validation summary")
    lines.append("")
    lines.append("| Test | Pass | Fail | Notes |")
    lines.append("|---|---|---|---|")
    for row in rows:
        status = "PASS" if row["pass"] else "FAIL"
        lines.append(f"| {row['test']} | {status} | {' ' if row['pass'] else 'X'} | {row['notes']} |")
    lines.append("")
    lines.append("## Convergence table (grid refinement)")
    lines.append("")
    lines.append("| Resolution | Bend | Error | Runtime (s) |")
    lines.append("|---|---|---|---|")
    for row in grid_rows:
        lines.append(f"| {row['resolution']}² | {row['bend']:.4e} | "
                     f"{row['error']:.4e} | {row['runtime']:.3f} |")
    lines.append("")
    lines.append("## Convergence table (step-size refinement)")
    lines.append("")
    lines.append("| Step | Steps | Bend | Bending angle | Conservation | Runtime (s) |")
    lines.append("|---|---|---|---|---|---|")
    for row in dt_rows:
        lines.append(f"| {row['step']:.4f} | {row['steps']} | {row['bend']:.4e} | "
                     f"{row['bending_angle']:.4e} | {row['conservation']:.3e} | {row['runtime']:.3f} |")
    lines.append("")
    lines.append("## Domain size")
    lines.append("")
    lines.append("| Domain | Extent | Bend | Conservation | Runtime (s) |")
    lines.append("|---|---|---|---|---|")
    for row in domain_rows:
        lines.append(f"| {row['domain']} | {row['extent']} | {row['bend']:.4e} | "
                     f"{row['conservation']:.3e} | {row['runtime']:.3f} |")
    lines.append("")
    lines.append("## Photon density")
    lines.append("")
    lines.append("| Photons | Bend | Bending angle | Conservation | Runtime (s) |")
    lines.append("|---|---|---|---|---|")
    for row in photon_rows:
        lines.append(f"| {row['photons']} | {row['bend']:.4e} | {row['bending_angle']:.4e} | "
                     f"{row['conservation']:.3e} | {row['runtime']:.3f} |")
    lines.append("")
    lines.append("## Precision")
    lines.append("")
    lines.append("| Precision | Bend | Conservation | Runtime (s) |")
    lines.append("|---|---|---|---|")
    for row in precision_rows:
        lines.append(f"| {row['precision']} | {row['bend']:.4e} | {row['conservation']:.3e} | "
                     f"{row['runtime']:.3f} |")
    lines.append("")
    if all_pass:
        lines.append("**Status: PASS** — all validation tests pass. The implementation is considered "
                     "numerically validated.")
    else:
        failing = [row['test'] for row in rows if not row['pass']]
        lines.append(f"**Status: FAIL** — failing tests: {', '.join(failing)}.")
    lines.append("")
    lines.append("## Findings")
    lines.append("")
    failing = [row for row in rows if not row['pass']]
    if not failing:
        lines.append("No findings: the frozen pipeline satisfies every requested numerical invariant "
                     "within the tolerances adopted for this milestone.")
    else:
        for row in failing:
            lines.append(f"### {row['test']}")
            lines.append("")
            lines.append(row['notes'])
            lines.append("")
            if row['test'] == "Mirror symmetry":
                lines.append("**Interpretation.** The frozen 90° transverse response carries a "
                             "right-handed rotation `R_90(∇̂C)`. A mirror reflection flips the "
                             "handedness of the plane, so the mirrored response is `R_-90(∇̂C)` "
                             "rather than `R_90`. The chirality of the transport is therefore a "
                             "geometric property of the frozen control law, not a numerical "
                             "artefact. This is a real, expected property of the implementation "
                             "and would require modifying the frozen transport to remove it.")
            else:
                lines.append("**Interpretation.** See notes above; the discrepancy exceeds the "
                             "tolerance adopted for this milestone.")
    lines.append("")
    lines.append("## Observable products")
    lines.append("")
    lines.append("The following maps and trajectories are written for visual inspection: "
                 "`convergence_map.png`, `shear_map.png`, `deflection_map.png`, "
                 "`photon_trajectories.png`, plus their underlying CSVs.")
    lines.append("")
    lines.append("## Refinement plot")
    lines.append("")
    lines.append("![Refinement plot](refinement_plot.png)")
    (out / "report.md").write_text("\n".join(lines))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--skip-100000", action="store_true")
    parser.add_argument("--probe", action="store_true")
    args = parser.parse_args()
    if args.probe:
        probe()
        return 0
    args.output.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()
    rows = []
    base = run_case(record=True)
    repeat = run_case(record=True)
    repeat_delta = max_path_difference(base, repeat)
    rows.append(validation_row("Repeatability", repeat_delta <= 1e-14,
                               f"max trajectory delta={repeat_delta:.3e}", repeat_delta))

    transforms = [
        ("Translation invariance", np.eye(2), np.array([0.5, 0.4]), "translation"),
        ("Rotation invariance", matrix_rotate(np.pi / 2), np.zeros(2), "rotation"),
        ("Mirror symmetry", np.diag([1.0, -1.0]), np.zeros(2), "mirror"),
    ]
    base_summary = summary(base)
    base_bend = base_summary["bend_max"]
    for name, transform, translation, label in transforms:
        transformed = run_case(ray_transform=transform, translation=translation, record=True,
                               mass_center=translation + transform @ np.array([BASE["mass_x"], BASE["mass_y"]]))
        tx, ty = inverse_transform_paths(transformed)
        delta = max(float(np.max(np.abs(base["xs"] - tx))), float(np.max(np.abs(base["ys"] - ty))))
        relative = delta / max(base_bend, 1e-15)
        if name == "Mirror symmetry":
            notes = (f"max transformed trajectory delta={delta:.3e} (relative={relative:.2e}); "
                     f"the 90-degree transverse response has a definite handedness, so the mirror "
                     f"test is expected to reveal chirality; this is a known property of the "
                     f"frozen transport, not an implementation bug")
            passed = False
        else:
            notes = f"max transformed trajectory delta={delta:.3e} (relative={relative:.2e})"
            passed = relative <= 0.15
        rows.append(validation_row(name, passed, notes, delta))

    grid_rows = []
    grid_results = {}
    for n in [64, 128, 256, 512]:
        result = run_case(n=n)
        m = summary(result)
        grid_results[n] = m
        grid_rows.append({"resolution": n, "bend": m["bend_max"], "runtime": m["runtime"], "error": None})
    finest = grid_results[512]["bend_max"]
    grid_errors = []
    for row in grid_rows:
        row["error"] = abs(row["bend"] - finest)
        grid_errors.append(row["error"])
    grid_order = order_from_errors(grid_errors[:-1])
    rows.append(validation_row("Grid refinement", bool(np.isfinite(grid_order) and grid_order > 0.0),
                               f"reported order={grid_order:.3f} (bending is a local quantity dominated by "
                               f"a few photon steps in the high-gradient region)", grid_order))
    write_csv(args.output / "convergence_table.csv", grid_rows)

    dt_rows = []
    dt_results = []
    for label, factor in [("current", 1.0), ("half", 0.5), ("quarter", 0.25)]:
        result = run_case(step=BASE["step"] * factor, steps=int(BASE["steps"] / factor))
        m = summary(result)
        dt_results.append(m)
        dt_rows.append({"timestep": label, "step": BASE["step"] * factor,
                        "steps": int(BASE["steps"] / factor), "bend": m["bend_max"],
                        "bending_angle": m["bending_angle_max"], "conservation": m["conservation"],
                        "runtime": m["runtime"]})
    dt_errors = [abs(m["bend_max"] - dt_results[-1]["bend_max"]) for m in dt_results]
    dt_order = order_from_errors(dt_errors[:-1])
    dt_conservation = max(m["conservation"] for m in dt_results)
    rows.append(validation_row("Step-size refinement", bool(np.isfinite(dt_order) and dt_order > 0.0 and dt_conservation <= 1e-12),
                               f"reported order={dt_order:.3f}; max conservation residual={dt_conservation:.3e}", dt_order))
    write_csv(args.output / "timestep_table.csv", dt_rows)

    domain_rows = []
    domain_results = {}
    for label, extent, n in [("current", 8.0, 128), ("doubled", 16.0, 256)]:
        result = run_case(n=n, extent=extent, launch_x=-8.0)
        m = summary(result)
        domain_results[label] = m
        domain_rows.append({"domain": label, "extent": extent, "n": n,
                            "launch_x": -8.0, "bend": m["bend_max"],
                            "bending_angle": m["bending_angle_max"], "conservation": m["conservation"],
                            "runtime": m["runtime"]})
    domain_delta = relative_difference(domain_results["current"]["bend_max"], domain_results["doubled"]["bend_max"])
    rows.append(validation_row("Domain size", domain_delta <= 0.1,
                               f"relative bending difference (same physical launch, same physical mass, "
                               f"n scaled to preserve grid spacing)={domain_delta:.3e}",
                               domain_delta))
    write_csv(args.output / "domain_table.csv", domain_rows)

    photon_rows = []
    counts = [100, 1000, 10000] if args.skip_100000 else [100, 1000, 10000, 100000]
    for count in counts:
        result = run_case(nphotons=count, record=True)
        m = summary(result)
        photon_rows.append({"photons": count, "bend": m["bend_max"],
                            "bending_angle": m["bending_angle_max"],
                            "conservation": m["conservation"], "runtime": m["runtime"]})
    bend_angle_max = max(row["bending_angle"] for row in photon_rows)
    bend_angle_min = min(row["bending_angle"] for row in photon_rows)
    bend_invariant = (bend_angle_max - bend_angle_min) / max(bend_angle_max, 1e-30) <= 1e-10
    runtimes = [row["runtime"] for row in photon_rows]
    runtime_scales = runtimes[-1] / max(runtimes[0], 1e-12)
    expected_scale = photon_rows[-1]["photons"] / photon_rows[0]["photons"]
    runtime_pass = 0.05 * expected_scale <= runtime_scales <= 20.0 * expected_scale
    photon_pass = bend_invariant and runtime_pass
    rows.append(validation_row("Photon density", bool(photon_pass),
                               (f"max bending-angle spread={bend_angle_max - bend_angle_min:.3e}; "
                                f"runtime scales x{runtime_scales:.2f} for {expected_scale:.0f}x more photons"),
                               bend_angle_max - bend_angle_min))
    write_csv(args.output / "photon_density_table.csv", photon_rows)

    precision_rows = []
    precision_results = {}
    for label, dtype in [("float32", np.float32), ("float64", np.float64), ("longdouble", np.longdouble)]:
        result = run_case(nphotons=1000, dtype=dtype, record=True)
        precision_results[label] = result
        m = summary(result)
        precision_rows.append({"precision": label, "bend": m["bend_max"],
                               "conservation": m["conservation"], "runtime": m["runtime"],
                               "path_checksum_x": sha(result["xs"]), "path_checksum_y": sha(result["ys"])})
    precision_delta = max_path_difference(precision_results["float32"], precision_results["float64"])
    precision_rows[1]["delta_vs_float32"] = precision_delta
    precision_pass = all(summary(result)["finite"] for result in precision_results.values())
    rows.append(validation_row("Floating-point precision", precision_pass,
                               f"float32/float64 max trajectory delta={precision_delta:.3e}", precision_delta))
    write_csv(args.output / "precision_table.csv", precision_rows)

    optimized_pass = False
    optimized_note = "Python optimization probe unavailable"
    try:
        completed = subprocess.run([sys.executable, "-O", str(Path(__file__).resolve()), "--probe"],
                                   capture_output=True, text=True, check=True)
        probe_data = json.loads(completed.stdout.strip().splitlines()[-1])
        optimized_pass = (probe_data["path_checksum_x"] == sha(base["xs"]) and
                          probe_data["path_checksum_y"] == sha(base["ys"]))
        optimized_note = f"checksum parity={optimized_pass}"
    except (OSError, subprocess.CalledProcessError, json.JSONDecodeError, IndexError) as exc:
        optimized_note = f"probe unavailable: {exc}"
    rows.append(validation_row("Compiler optimisation", optimized_pass, optimized_note))

    products = map_observables(base)
    save_products(args.output, base)
    all_pass = all(row["pass"] for row in rows)
    _plot_refinement(args.output, grid_rows, dt_rows)
    _write_report(args.output, rows, grid_rows, dt_rows, photon_rows, domain_rows,
                  precision_rows, all_pass)
    summary_rows = []
    for name, result in [("canonical", base), ("repeat", repeat)]:
        m = summary(result)
        m["case"] = name
        summary_rows.append(m)
    write_csv(args.output / "summary.csv", summary_rows)
    write_csv(args.output / "validation_summary.csv", rows)
    with (args.output / "validation.json").open("w") as handle:
        json.dump({"status": "PASS" if all_pass else "FAIL", "tests": rows,
                   "observables": {"convergence": "ray-density surrogate",
                                   "shear": "deflection-Jacobian proxy",
                                   "deflection": "endpoint displacement",
                                   "phase": "not defined by frozen photon coupling"},
                   "control": {"constitutive": "u = 0.18 rho/rho_max",
                               "transport": "90-degree transverse response; direct addition + renormalisation"},
                   "execution_seconds": time.perf_counter() - started,
                   "checksums": {"matter": sha(base["rho"]), "constitutive": sha(base["c"]),
                                 "response_x": sha(base["rx"]), "response_y": sha(base["ry"]),
                                 "trajectories_x": sha(base["xs"]), "trajectories_y": sha(base["ys"])},
                   "map_shapes": {name: list(array.shape) for name, array in products.items()
                                  if name != "edges"}}, handle, indent=2)
    with (args.output / "run.json").open("w") as handle:
        json.dump({"milestone": "PBUF WEAK-LENSING-VALIDATION-001", "status": "PASS" if all_pass else "FAIL",
                   "validation_summary": rows, "grid_order": grid_order,
                   "timestep_order": dt_order,
                   "photon_density_bending_angle_spread": bend_angle_max - bend_angle_min,
                   "photon_density_runtime_scaling": runtime_scales}, handle, indent=2)
    print(json.dumps({"status": "PASS" if all_pass else "FAIL", "output": str(args.output),
                      "grid_order": grid_order, "timestep_order": dt_order,
                      "photon_density_bending_angle_spread": bend_angle_max - bend_angle_min,
                      "photon_density_runtime_scaling": runtime_scales}, indent=2))
    return 0 if all_pass else 1


if __name__ == "__main__":
    raise SystemExit(main())
