#!/usr/bin/env python3
"""PBUF WEAK-LENSING-PREDICTION-001 - first complete end-to-end prediction of Version A."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import time
from dataclasses import asdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from constitutive_equations import get_equation


ROOT = Path(__file__).resolve().parent
DEFAULT_OUT = ROOT / "runs" / "weak_lensing_prediction001"
LENS = {
    "extent": 8.0,
    "mass_x": -0.65,
    "mass_y": 0.0,
    "mass_sigma": 0.75,
    "mass_amplitude": 1.0,
    "n": 128,
    "strength": 0.18,
    "step": 0.06,
    "steps": 80,
    "y_span": 3.0,
    "nphotons": 2000,
}


def sha(array: np.ndarray) -> str:
    return hashlib.sha256(np.ascontiguousarray(array).tobytes()).hexdigest()


def make_field(n, extent, mass_center, mass_sigma, strength, dtype=np.float64):
    x = np.linspace(-extent, extent, n, dtype=dtype)
    y = np.linspace(-extent, extent, n, dtype=dtype)
    X, Y = np.meshgrid(x, y, indexing="xy")
    rho = (mass_sigma
           * np.exp(-((X - mass_center[0]) ** 2 + (Y - mass_center[1]) ** 2) /
                    (2 * mass_sigma ** 2))).astype(dtype)
    cfg = type("Config", (), {"deformation_strength": strength})()
    c = get_equation("A").solve(rho, cfg)
    gy, gx = np.gradient(c, x, y, edge_order=1)
    g = np.hypot(gx, gy)
    gx_hat = gx / np.maximum(g, np.asarray(1e-15, dtype=dtype))
    gy_hat = gy / np.maximum(g, np.asarray(1e-15, dtype=dtype))
    bad = g < np.asarray(1e-15, dtype=dtype)
    gx_hat = np.where(bad, np.asarray(1.0, dtype=dtype), gx_hat)
    gy_hat = np.where(bad, np.asarray(0.0, dtype=dtype), gy_hat)
    rx = -g * gy_hat
    ry = g * gx_hat
    return {"xgrid": x, "ygrid": y, "X": X, "Y": Y, "rho": rho, "c": c,
            "gx": gx, "gy": gy, "g_magnitude": g, "rx": rx, "ry": ry,
            "response_direction": np.arctan2(ry, rx)}


def propagate(xgrid, ygrid, rx, ry, step, steps, x0, y0, vx0, vy0, record=True):
    dtype = np.result_type(xgrid.dtype, ygrid.dtype, rx.dtype,
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
        rx_loc = rx[iy, ix]
        ry_loc = ry[iy, ix]
        vx_new = vx + np.asarray(step, dtype=dtype) * rx_loc
        vy_new = vy + np.asarray(step, dtype=dtype) * ry_loc
        scale = np.maximum(np.hypot(vx_new, vy_new), np.asarray(1e-12, dtype=dtype))
        vx_unit = vx_new / scale
        vy_unit = vy_new / scale
        conservation = np.maximum(conservation, np.abs(np.hypot(vx_unit, vy_unit) - 1))
        dot = np.clip(vx * vx_unit + vy * vy_unit, -1, 1)
        bending_angle += np.arccos(dot)
        vx = vx_unit
        vy = vy_unit
        x = x + np.asarray(step, dtype=dtype) * vx
        y = y + np.asarray(step, dtype=dtype) * vy
        max_deviation = np.maximum(max_deviation, np.abs(y - np.asarray(y0, dtype=dtype)))
        if record:
            xs[:, k] = x; ys[:, k] = y
    runtime = time.perf_counter() - started
    return {"x": x, "y": y, "max_deviation": max_deviation,
            "bending_angle": bending_angle, "conservation": conservation,
            "runtime": runtime, "xs": xs, "ys": ys}


def compute_observables(field, photons, bins=64):
    xf = photons["x"]; yf = photons["y"]
    x0 = photons["x0"]; y0 = photons["y0"]
    extent = LENS["extent"]
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
    gamma_mag = np.hypot(shear_g1, shear_g2)
    denom = (1.0 - np.nan_to_num(convergence, nan=0.0)) ** 2 - gamma_mag ** 2
    magnification = np.full_like(convergence, np.nan)
    positive = (denom > 0) & good
    magnification[positive] = 1.0 / denom[positive]
    return {"convergence": convergence, "shear_g1": shear_g1, "shear_g2": shear_g2,
            "shear_magnitude": gamma_mag, "deflection_x": deflection_x,
            "deflection_y": deflection_y, "deflection_magnitude": np.hypot(deflection_x, deflection_y),
            "magnification": magnification, "ray_count": final_count, "edges": edges}


def save_field_csv(out, name, array):
    np.savetxt(out / f"{name}.csv", array, delimiter=",")


def make_overlay(field, photons, out):
    fig, ax = plt.subplots(figsize=(8, 7))
    x = field["xgrid"]
    extent = [x[0], x[-1], x[0], x[-1]]
    cs = ax.contourf(field["X"], field["Y"], field["c"], levels=20, alpha=0.4, cmap="viridis")
    fig.colorbar(cs, ax=ax, label="C")
    bending = photons["bending_angle"]
    norm = plt.Normalize(vmin=float(bending.min()), vmax=float(bending.max()))
    cmap = plt.cm.plasma
    for i in range(0, len(photons["xs"]), max(1, len(photons["xs"]) // 200)):
        xs = photons["xs"][i]
        ys = photons["ys"][i]
        points = np.array([xs, ys]).T.reshape(-1, 1, 2)
        segments = np.concatenate([points[:-1], points[1:]], axis=1)
        lc = __import__("matplotlib").collections.LineCollection(
            segments, cmap=cmap, norm=norm, linewidth=0.6)
        lc.set_array(np.linspace(bending[i], bending[i], len(segments)))
        ax.add_collection(lc)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    fig.colorbar(sm, ax=ax, label="accumulated bending angle (rad)")
    ax.set(xlabel="x", ylabel="y", title="Photon trajectories over constitutive field",
           xlim=(extent[0], extent[1]), ylim=(extent[2], extent[3]), aspect="equal")
    fig.tight_layout()
    fig.savefig(out / "photon_trajectories_overlay.png", dpi=140)
    plt.close(fig)


def make_trajectory_color(field, photons, out):
    fig, ax = plt.subplots(figsize=(8, 7))
    bending = photons["bending_angle"]
    norm = plt.Normalize(vmin=float(bending.min()), vmax=float(bending.max()))
    for xs, ys, b in zip(photons["xs"], photons["ys"], bending):
        ax.plot(xs, ys, color=plt.cm.plasma(norm(b)), lw=0.5)
    sm = plt.cm.ScalarMappable(cmap=plt.cm.plasma, norm=norm)
    sm.set_array([])
    fig.colorbar(sm, ax=ax, label="accumulated bending angle (rad)")
    ax.set(xlabel="x", ylabel="y", title="Photon trajectories coloured by accumulated deflection",
           xlim=(-8, 8), ylim=(-8, 8), aspect="equal")
    fig.tight_layout()
    fig.savefig(out / "photon_trajectories.png", dpi=140)
    plt.close(fig)


def make_map(out, name, array, title, cmap="viridis", symmetric=False, vmin=None, vmax=None):
    fig, ax = plt.subplots(figsize=(7, 6))
    if symmetric:
        vmax_abs = float(np.nanmax(np.abs(array)))
        im = ax.imshow(array, origin="lower", extent=[-8, 8, -8, 8],
                       cmap=cmap, vmin=-vmax_abs, vmax=vmax_abs)
    else:
        im = ax.imshow(array, origin="lower", extent=[-8, 8, -8, 8],
                       cmap=cmap, vmin=vmin, vmax=vmax)
    ax.set(xlabel="x", ylabel="y", title=title, aspect="equal")
    fig.colorbar(im, ax=ax)
    fig.tight_layout()
    fig.savefig(out / f"{name}.png", dpi=140)
    plt.close(fig)


def make_deflection_quiver(out, def_x, def_y):
    fig, ax = plt.subplots(figsize=(7, 6))
    x = np.linspace(-8, 8, def_x.shape[1])
    y = np.linspace(-8, 8, def_x.shape[0])
    X, Y = np.meshgrid(x, y, indexing="xy")
    stride = max(1, def_x.shape[0] // 16)
    ax.quiver(X[::stride, ::stride], Y[::stride, ::stride],
              def_x[::stride, ::stride], def_y[::stride, ::stride],
              np.hypot(def_x[::stride, ::stride], def_y[::stride, ::stride]),
              cmap="plasma", pivot="middle")
    ax.set(xlabel="x", ylabel="y", title="Deflection vectors (α_x, α_y)", aspect="equal",
           xlim=(-8, 8), ylim=(-8, 8))
    fig.tight_layout()
    fig.savefig(out / "deflection_vectors.png", dpi=140)
    plt.close(fig)


def make_composite(out, field, observables, photons):
    fig, axes = plt.subplots(3, 3, figsize=(18, 16))
    extent_kw = dict(origin="lower", extent=[-8, 8, -8, 8])
    fields = [
        ("Matter ρ", field["rho"], "viridis", False),
        ("Constitutive C", field["c"], "viridis", False),
        ("|∇C|", field["g_magnitude"], "viridis", False),
        ("Response magnitude", np.hypot(field["rx"], field["ry"]), "viridis", False),
        ("Response direction (rad)", field["response_direction"], "twilight", True),
        ("Convergence κ", observables["convergence"], "RdBu_r", True),
        ("Shear γ₁", observables["shear_g1"], "RdBu_r", True),
        ("Shear γ₂", observables["shear_g2"], "RdBu_r", True),
        ("Magnification μ", observables["magnification"], "viridis", False),
    ]
    for ax, (title, array, cmap, sym) in zip(axes.flat, fields):
        if sym:
            vmax_abs = float(np.nanmax(np.abs(array))) if np.any(np.isfinite(array)) else 1.0
            im = ax.imshow(array, cmap=cmap, vmin=-vmax_abs, vmax=vmax_abs, **extent_kw)
        else:
            im = ax.imshow(array, cmap=cmap, **extent_kw)
        ax.set(title=title, xlabel="x", ylabel="y", aspect="equal",
               xlim=(-8, 8), ylim=(-8, 8))
        fig.colorbar(im, ax=ax)
    fig.tight_layout()
    fig.savefig(out / "composite_observables.png", dpi=140)
    plt.close(fig)


def compute_statistics(field, photons, observables):
    deflection = observables["deflection_magnitude"]
    convergence = observables["convergence"]
    shear = observables["shear_magnitude"]
    bending = photons["bending_angle"]
    return {
        "max_bending": float(np.max(photons["max_deviation"])),
        "mean_bending": float(np.mean(photons["max_deviation"])),
        "rms_bending": float(np.sqrt(np.mean(photons["max_deviation"] ** 2))),
        "max_bending_angle": float(np.max(bending)),
        "mean_bending_angle": float(np.mean(bending)),
        "max_deflection": float(np.nanmax(deflection)) if np.any(np.isfinite(deflection)) else 0.0,
        "max_convergence": float(np.nanmax(np.abs(convergence))) if np.any(np.isfinite(convergence)) else 0.0,
        "mean_convergence": float(np.nanmean(convergence)) if np.any(np.isfinite(convergence)) else 0.0,
        "max_shear": float(np.nanmax(shear)) if np.any(np.isfinite(shear)) else 0.0,
        "mean_shear": float(np.nanmean(shear)) if np.any(np.isfinite(shear)) else 0.0,
        "max_conservation": float(np.max(photons["conservation"])),
        "photon_count": int(len(photons["x"])),
        "propagation_runtime_seconds": float(photons["runtime"]),
        "n_grid": LENS["n"],
        "n_steps": LENS["steps"],
        "step_size": LENS["step"],
        "max_response_magnitude": float(np.max(np.hypot(field["rx"], field["ry"]))),
        "max_C": float(np.max(field["c"])),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()
    out = args.output
    out.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()
    field = make_field(LENS["n"], LENS["extent"],
                       np.array([LENS["mass_x"], LENS["mass_y"]]),
                       LENS["mass_sigma"], LENS["strength"])
    x0 = np.full(LENS["nphotons"], -LENS["extent"])
    y0 = np.linspace(-LENS["y_span"], LENS["y_span"], LENS["nphotons"])
    vx0 = np.ones(LENS["nphotons"])
    vy0 = np.zeros(LENS["nphotons"])
    photons = propagate(field["xgrid"], field["ygrid"], field["rx"], field["ry"],
                        LENS["step"], LENS["steps"], x0, y0, vx0, vy0, record=True)
    photons["x0"] = x0
    photons["y0"] = y0
    observables = compute_observables(field, photons, bins=64)
    statistics = compute_statistics(field, photons, observables)
    save_field_csv(out, "matter", field["rho"])
    save_field_csv(out, "constitutive", field["c"])
    save_field_csv(out, "gradient_magnitude", field["g_magnitude"])
    save_field_csv(out, "gradient_x", field["gx"])
    save_field_csv(out, "gradient_y", field["gy"])
    save_field_csv(out, "response_x", field["rx"])
    save_field_csv(out, "response_y", field["ry"])
    save_field_csv(out, "response_direction", field["response_direction"])
    save_field_csv(out, "convergence", observables["convergence"])
    save_field_csv(out, "shear_g1", observables["shear_g1"])
    save_field_csv(out, "shear_g2", observables["shear_g2"])
    save_field_csv(out, "shear_magnitude", observables["shear_magnitude"])
    save_field_csv(out, "deflection_x", observables["deflection_x"])
    save_field_csv(out, "deflection_y", observables["deflection_y"])
    save_field_csv(out, "deflection_magnitude", observables["deflection_magnitude"])
    save_field_csv(out, "magnification", observables["magnification"])
    save_field_csv(out, "photon_trajectories_x", photons["xs"])
    save_field_csv(out, "photon_trajectories_y", photons["ys"])
    save_field_csv(out, "photon_max_deviation", photons["max_deviation"])
    save_field_csv(out, "photon_bending_angle", photons["bending_angle"])
    np.savetxt(out / "photon_endpoints.csv",
               np.column_stack([photons["x0"], photons["y0"], photons["x"], photons["y"],
                                photons["max_deviation"], photons["bending_angle"]]),
               delimiter=",", header="x0,y0,x_final,y_final,max_deviation,bending_angle",
               comments="")
    make_map(out, "matter_map", field["rho"], "Matter density ρ", "viridis")
    make_map(out, "constitutive_map", field["c"], "Constitutive field C", "viridis")
    make_map(out, "gradient_map", field["g_magnitude"], "|∇C|", "viridis")
    make_map(out, "response_magnitude_map", np.hypot(field["rx"], field["ry"]),
             "Response magnitude |r|", "viridis")
    make_map(out, "response_direction_map", field["response_direction"],
             "Response direction (rad)", "twilight", symmetric=True)
    make_map(out, "convergence_map", observables["convergence"],
             "Convergence κ", "RdBu_r", symmetric=True)
    make_map(out, "shear_g1_map", observables["shear_g1"],
             "Shear γ₁", "RdBu_r", symmetric=True)
    make_map(out, "shear_g2_map", observables["shear_g2"],
             "Shear γ₂", "RdBu_r", symmetric=True)
    make_map(out, "shear_magnitude_map", observables["shear_magnitude"],
             "|γ| = √(γ₁² + γ₂²)", "viridis")
    make_map(out, "deflection_magnitude_map", observables["deflection_magnitude"],
             "Deflection magnitude |α|", "viridis")
    make_map(out, "magnification_map", observables["magnification"],
             "Magnification μ = 1/((1-κ)² - |γ|²)", "viridis")
    make_deflection_quiver(out, observables["deflection_x"], observables["deflection_y"])
    make_trajectory_color(field, photons, out)
    make_overlay(field, photons, out)
    make_composite(out, field, observables, photons)
    stats_with_checksums = dict(statistics)
    stats_with_checksums["checksums"] = {
        "matter": sha(field["rho"]),
        "constitutive": sha(field["c"]),
        "gradient_magnitude": sha(field["g_magnitude"]),
        "response_x": sha(field["rx"]),
        "response_y": sha(field["ry"]),
        "photon_trajectories_x": sha(photons["xs"]),
        "photon_trajectories_y": sha(photons["ys"]),
        "convergence": sha(observables["convergence"]),
        "shear_g1": sha(observables["shear_g1"]),
        "shear_g2": sha(observables["shear_g2"]),
        "magnification": sha(observables["magnification"]),
    }
    stats_with_checksums["frozen_pipeline"] = {
        "constitutive": "Version A: C(X) = 0.18 * rho(X) / rho_max",
        "transport": "90-degree transverse response, direct addition + renormalisation",
        "amplitude": "A = |grad C|",
        "lens": "Lens-001 frozen dataset",
    }
    stats_with_checksums["execution_seconds"] = time.perf_counter() - started
    (out / "prediction_statistics.json").write_text(json.dumps(stats_with_checksums, indent=2))
    (out / "run.json").write_text(json.dumps({
        "milestone": "PBUF WEAK-LENSING-PREDICTION-001",
        "status": "OK",
        "execution_seconds": time.perf_counter() - started,
    }, indent=2))
    with (out / "summary.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["metric", "value"])
        for key, value in statistics.items():
            writer.writerow([key, value])
    report_lines = [
        "# PBUF WEAK-LENSING-PREDICTION-001",
        "",
        "Frozen pipeline: matter → C = 0.18 · ρ/ρ_max → ∇C → |∇C| response → 90° transverse "
        "response → direct addition + renormalisation → photon propagation → observables.",
        "",
        "## Integrated statistics",
        "",
        "| Metric | Value |",
        "|---|---|",
    ]
    for key, value in statistics.items():
        if isinstance(value, float):
            report_lines.append(f"| {key} | {value:.6e} |")
        else:
            report_lines.append(f"| {key} | {value} |")
    report_lines += [
        "",
        "## Products",
        "",
        "- `matter_map.png`, `constitutive_map.png`, `gradient_map.png`",
        "- `response_magnitude_map.png`, `response_direction_map.png`",
        "- `photon_trajectories.png`, `photon_trajectories_overlay.png`",
        "- `convergence_map.png`, `shear_g1_map.png`, `shear_g2_map.png`, `shear_magnitude_map.png`",
        "- `deflection_magnitude_map.png`, `deflection_vectors.png`",
        "- `magnification_map.png`",
        "- `composite_observables.png` (3×3 panel)",
        "- Underlying CSVs for every field and map",
        "",
        "## Notes",
        "",
        "All products were generated by the frozen Version A pipeline. No comparison with "
        "ΛCDM, GR, or observations is made at this milestone.",
    ]
    (out / "report.md").write_text("\n".join(report_lines))
    print(json.dumps({"milestone": "PBUF WEAK-LENSING-PREDICTION-001",
                      "status": "OK", "output": str(out),
                      "photon_count": statistics["photon_count"],
                      "runtime_seconds": statistics["propagation_runtime_seconds"]}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
