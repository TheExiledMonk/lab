#!/usr/bin/env python3
"""PBUF WEAK-LENSING-GENERALIZATION-001 - multi-dataset prediction validation."""
from __future__ import annotations

import csv
import hashlib
import json
import time
from dataclasses import dataclass
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from constitutive_equations import get_equation


ROOT = Path(__file__).resolve().parent
DEFAULT_OUT = ROOT / "runs" / "weak_lensing_generalization001"
LENS = {
    "n": 128,
    "extent": 8.0,
    "strength": 0.18,
    "step": 0.06,
    "steps": 80,
    "y_span": 3.0,
    "nphotons": 2000,
    "bins": 64,
}
DATASETS = [
    {
        "id": "isolated",
        "label": "Isolated cluster",
        "rho_fn": lambda X, Y: gaussian(X, Y, (-0.65, 0.0), (0.75, 0.75), 1.0),
    },
    {
        "id": "binary",
        "label": "Binary cluster",
        "rho_fn": lambda X, Y: (gaussian(X, Y, (-2.0, 0.0), (0.7, 0.7), 0.6)
                               + gaussian(X, Y, (2.0, 0.5), (0.7, 0.7), 0.6)),
    },
    {
        "id": "elongated",
        "label": "Elongated cluster",
        "rho_fn": lambda X, Y: gaussian(X, Y, (0.0, 0.0), (2.0, 0.4), 1.0),
    },
    {
        "id": "asymmetric",
        "label": "Asymmetric cluster",
        "rho_fn": lambda X, Y: (gaussian(X, Y, (0.0, 0.0), (0.7, 0.7), 1.0)
                               + gaussian(X, Y, (1.6, 0.6), (0.5, 0.5), 0.4)),
    },
    {
        "id": "sparse",
        "label": "Sparse field",
        "rho_fn": lambda X, Y: gaussian(X, Y, (0.0, 0.0), (2.5, 2.5), 0.8),
    },
    {
        "id": "dense",
        "label": "Dense field",
        "rho_fn": lambda X, Y: gaussian(X, Y, (0.0, 0.0), (0.35, 0.35), 1.5),
    },
]


def gaussian(X, Y, center, sigma, amplitude):
    return amplitude * np.exp(-(((X - center[0]) / sigma[0]) ** 2 +
                                ((Y - center[1]) / sigma[1]) ** 2) / 2.0)


def sha(array: np.ndarray) -> str:
    return hashlib.sha256(np.ascontiguousarray(array).tobytes()).hexdigest()


def file_sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def make_field(rho: np.ndarray, extent: float, strength: float, n: int):
    x = np.linspace(-extent, extent, n)
    y = np.linspace(-extent, extent, n)
    X, Y = np.meshgrid(x, y, indexing="xy")
    cfg = type("Config", (), {"deformation_strength": strength})()
    c = get_equation("A").solve(rho, cfg)
    gy, gx = np.gradient(c, x, y, edge_order=1)
    g = np.hypot(gx, gy)
    gx_hat = gx / np.maximum(g, 1e-15)
    gy_hat = gy / np.maximum(g, 1e-15)
    bad = g < 1e-15
    gx_hat = np.where(bad, 1.0, gx_hat)
    gy_hat = np.where(bad, 0.0, gy_hat)
    rx = -g * gy_hat
    ry = g * gx_hat
    return {"xgrid": x, "ygrid": y, "X": X, "Y": Y, "rho": rho, "c": c,
            "gx": gx, "gy": gy, "g_magnitude": g, "rx": rx, "ry": ry,
            "response_direction": np.arctan2(ry, rx)}


def propagate(xgrid, ygrid, rx, ry, step, steps, x0, y0, vx0, vy0, record=True):
    x = x0.copy(); y = y0.copy()
    vx = vx0.copy(); vy = vy0.copy()
    nphotons = len(x)
    max_deviation = np.zeros(nphotons)
    bending_angle = np.zeros(nphotons)
    conservation = np.zeros(nphotons)
    xs = np.empty((nphotons, steps)) if record else None
    ys = np.empty((nphotons, steps)) if record else None
    if record:
        xs[:, 0] = x; ys[:, 0] = y
    started = time.perf_counter()
    for k in range(1, steps):
        ix = np.clip(np.searchsorted(xgrid, x) - 1, 0, len(xgrid) - 1)
        iy = np.clip(np.searchsorted(ygrid, y) - 1, 0, len(ygrid) - 1)
        rx_loc = rx[iy, ix]
        ry_loc = ry[iy, ix]
        vx_new = vx + step * rx_loc
        vy_new = vy + step * ry_loc
        scale = np.maximum(np.hypot(vx_new, vy_new), 1e-12)
        vx_unit = vx_new / scale
        vy_unit = vy_new / scale
        conservation = np.maximum(conservation, np.abs(np.hypot(vx_unit, vy_unit) - 1))
        dot = np.clip(vx * vx_unit + vy * vy_unit, -1, 1)
        bending_angle += np.arccos(dot)
        vx = vx_unit
        vy = vy_unit
        x = x + step * vx
        y = y + step * vy
        max_deviation = np.maximum(max_deviation, np.abs(y - y0))
        if record:
            xs[:, k] = x; ys[:, k] = y
    runtime = time.perf_counter() - started
    return {"x": x, "y": y, "max_deviation": max_deviation,
            "bending_angle": bending_angle, "conservation": conservation,
            "runtime": runtime, "xs": xs, "ys": ys}


def compute_observables(field, photons, bins):
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


def make_map(out, name, array, title, cmap="viridis", symmetric=False):
    fig, ax = plt.subplots(figsize=(7, 6))
    if symmetric:
        vmax_abs = float(np.nanmax(np.abs(array))) if np.any(np.isfinite(array)) else 1.0
        im = ax.imshow(array, origin="lower", extent=[-8, 8, -8, 8],
                       cmap=cmap, vmin=-vmax_abs, vmax=vmax_abs)
    else:
        im = ax.imshow(array, origin="lower", extent=[-8, 8, -8, 8], cmap=cmap)
    ax.set(xlabel="x", ylabel="y", title=title, aspect="equal",
           xlim=(-8, 8), ylim=(-8, 8))
    fig.colorbar(im, ax=ax)
    fig.tight_layout()
    fig.savefig(out / f"{name}.png", dpi=140)
    plt.close(fig)


def make_trajectory_plot(out, photons):
    fig, ax = plt.subplots(figsize=(8, 7))
    bending = photons["bending_angle"]
    norm = plt.Normalize(vmin=float(bending.min()), vmax=float(bending.max()))
    for xs, ys, b in zip(photons["xs"], photons["ys"], bending):
        ax.plot(xs, ys, color=plt.cm.plasma(norm(b)), lw=0.5)
    sm = plt.cm.ScalarMappable(cmap=plt.cm.plasma, norm=norm)
    sm.set_array([])
    fig.colorbar(sm, ax=ax, label="accumulated bending angle (rad)")
    ax.set(xlabel="x", ylabel="y",
           title="Photon trajectories (coloured by accumulated deflection)",
           xlim=(-8, 8), ylim=(-8, 8), aspect="equal")
    fig.tight_layout()
    fig.savefig(out / "photon_trajectories.png", dpi=140)
    plt.close(fig)


def make_composite(out, field, observables):
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
    bending = photons["max_deviation"]
    finite_outputs = bool(np.isfinite(photons["x"]).all()
                          and np.isfinite(photons["y"]).all()
                          and np.isfinite(observables["convergence"][~np.isnan(observables["convergence"])]).all() if np.any(~np.isnan(observables["convergence"])) else True)
    return {
        "photon_count": int(len(photons["x"])),
        "propagation_runtime_seconds": float(photons["runtime"]),
        "max_bending": float(np.max(bending)),
        "rms_bending": float(np.sqrt(np.mean(bending ** 2))),
        "max_bending_angle": float(np.max(photons["bending_angle"])),
        "max_convergence": float(np.nanmax(np.abs(convergence))) if np.any(np.isfinite(convergence)) else 0.0,
        "rms_convergence": float(np.sqrt(np.nanmean(convergence ** 2))) if np.any(np.isfinite(convergence)) else 0.0,
        "max_shear": float(np.nanmax(shear)) if np.any(np.isfinite(shear)) else 0.0,
        "rms_shear": float(np.sqrt(np.nanmean(shear ** 2))) if np.any(np.isfinite(shear)) else 0.0,
        "max_deflection": float(np.nanmax(deflection)) if np.any(np.isfinite(deflection)) else 0.0,
        "max_conservation": float(np.max(photons["conservation"])),
        "finite_outputs": finite_outputs,
        "max_C": float(np.max(field["c"])),
        "max_response_magnitude": float(np.max(np.hypot(field["rx"], field["ry"]))),
    }


def run_dataset(dataset, out_root, executable_hashes):
    dataset_dir = out_root / dataset["id"]
    dataset_dir.mkdir(parents=True, exist_ok=True)
    x = np.linspace(-LENS["extent"], LENS["extent"], LENS["n"])
    y = np.linspace(-LENS["extent"], LENS["extent"], LENS["n"])
    X, Y = np.meshgrid(x, y, indexing="xy")
    rho = dataset["rho_fn"](X, Y)
    field = make_field(rho, LENS["extent"], LENS["strength"], LENS["n"])
    x0 = np.full(LENS["nphotons"], -LENS["extent"])
    y0 = np.linspace(-LENS["y_span"], LENS["y_span"], LENS["nphotons"])
    vx0 = np.ones(LENS["nphotons"])
    vy0 = np.zeros(LENS["nphotons"])
    photons = propagate(field["xgrid"], field["ygrid"], field["rx"], field["ry"],
                        LENS["step"], LENS["steps"], x0, y0, vx0, vy0, record=True)
    photons["x0"] = x0
    photons["y0"] = y0
    observables = compute_observables(field, photons, bins=LENS["bins"])
    statistics = compute_statistics(field, photons, observables)
    np.savetxt(dataset_dir / "matter.csv", field["rho"], delimiter=",")
    np.savetxt(dataset_dir / "constitutive.csv", field["c"], delimiter=",")
    np.savetxt(dataset_dir / "gradient_magnitude.csv", field["g_magnitude"], delimiter=",")
    np.savetxt(dataset_dir / "gradient_x.csv", field["gx"], delimiter=",")
    np.savetxt(dataset_dir / "gradient_y.csv", field["gy"], delimiter=",")
    np.savetxt(dataset_dir / "response_x.csv", field["rx"], delimiter=",")
    np.savetxt(dataset_dir / "response_y.csv", field["ry"], delimiter=",")
    np.savetxt(dataset_dir / "convergence.csv", observables["convergence"], delimiter=",")
    np.savetxt(dataset_dir / "shear_g1.csv", observables["shear_g1"], delimiter=",")
    np.savetxt(dataset_dir / "shear_g2.csv", observables["shear_g2"], delimiter=",")
    np.savetxt(dataset_dir / "shear_magnitude.csv", observables["shear_magnitude"], delimiter=",")
    np.savetxt(dataset_dir / "deflection_x.csv", observables["deflection_x"], delimiter=",")
    np.savetxt(dataset_dir / "deflection_y.csv", observables["deflection_y"], delimiter=",")
    np.savetxt(dataset_dir / "deflection_magnitude.csv", observables["deflection_magnitude"], delimiter=",")
    np.savetxt(dataset_dir / "magnification.csv", observables["magnification"], delimiter=",")
    np.savetxt(dataset_dir / "photon_trajectories_x.csv", photons["xs"], delimiter=",")
    np.savetxt(dataset_dir / "photon_trajectories_y.csv", photons["ys"], delimiter=",")
    np.savetxt(dataset_dir / "photon_max_deviation.csv", photons["max_deviation"], delimiter=",")
    np.savetxt(dataset_dir / "photon_bending_angle.csv", photons["bending_angle"], delimiter=",")
    make_map(dataset_dir, "matter_map", field["rho"], f"{dataset['label']}: matter ρ")
    make_map(dataset_dir, "constitutive_map", field["c"], f"{dataset['label']}: C")
    make_map(dataset_dir, "gradient_map", field["g_magnitude"], f"{dataset['label']}: |∇C|")
    make_map(dataset_dir, "response_magnitude_map", np.hypot(field["rx"], field["ry"]),
             f"{dataset['label']}: response magnitude")
    make_map(dataset_dir, "response_direction_map", field["response_direction"],
             f"{dataset['label']}: response direction (rad)", symmetric=True)
    make_map(dataset_dir, "convergence_map", observables["convergence"],
             f"{dataset['label']}: κ", symmetric=True)
    make_map(dataset_dir, "shear_g1_map", observables["shear_g1"],
             f"{dataset['label']}: γ₁", symmetric=True)
    make_map(dataset_dir, "shear_g2_map", observables["shear_g2"],
             f"{dataset['label']}: γ₂", symmetric=True)
    make_map(dataset_dir, "shear_magnitude_map", observables["shear_magnitude"],
             f"{dataset['label']}: |γ|")
    make_map(dataset_dir, "deflection_magnitude_map", observables["deflection_magnitude"],
             f"{dataset['label']}: |α|")
    make_map(dataset_dir, "magnification_map", observables["magnification"],
             f"{dataset['label']}: μ")
    make_trajectory_plot(dataset_dir, photons)
    make_composite(dataset_dir, field, observables)
    statistics["checksums"] = {
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
    statistics["identical_pipeline"] = executable_hashes
    (dataset_dir / "dataset_statistics.json").write_text(json.dumps(statistics, indent=2))
    return statistics


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()
    out_root = args.output
    out_root.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()
    executable_hashes = {
        "weak_lensing_generalization001.py": file_sha(Path(__file__).resolve()),
        "constitutive_equations.py": file_sha(ROOT / "constitutive_equations.py"),
    }
    cross_rows = []
    per_dataset_statistics = {}
    for dataset in DATASETS:
        print(f"Running {dataset['label']} ({dataset['id']})...")
        statistics = run_dataset(dataset, out_root, executable_hashes)
        per_dataset_statistics[dataset["id"]] = statistics
        cross_rows.append({
            "dataset": dataset["label"],
            "dataset_id": dataset["id"],
            "runtime": statistics["propagation_runtime_seconds"],
            "max_bending": statistics["max_bending"],
            "rms_bending": statistics["rms_bending"],
            "max_convergence": statistics["max_convergence"],
            "rms_convergence": statistics["rms_convergence"],
            "max_shear": statistics["max_shear"],
            "rms_shear": statistics["rms_shear"],
            "max_conservation": statistics["max_conservation"],
            "finite_outputs": statistics["finite_outputs"],
            "photon_count": statistics["photon_count"],
        })
    with (out_root / "cross_dataset_summary.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(cross_rows[0]))
        writer.writeheader()
        writer.writerows(cross_rows)
    summary_lines = [
        "# PBUF WEAK-LENSING-GENERALIZATION-001",
        "",
        "Frozen pipeline: matter → C = 0.18 · ρ/ρ_max → ∇C → |∇C| response → "
        "90° transverse response → direct addition + renormalisation → photon propagation "
        "→ observables.",
        "",
        "## Identical pipeline check (SHA-256)",
        "",
        "| File | SHA-256 |",
        "|---|---|",
    ]
    for name, digest in executable_hashes.items():
        summary_lines.append(f"| {name} | `{digest}` |")
    summary_lines += [
        "",
        "## Cross-dataset summary",
        "",
        "| Dataset | Runtime (s) | Max bend | RMS bend | Max κ | RMS κ | Max γ | RMS γ | Cons. |",
        "|---|---|---|---|---|---|---|---|---|",
    ]
    for row in cross_rows:
        summary_lines.append(
            f"| {row['dataset']} | {row['runtime']:.3f} | {row['max_bending']:.3e} | "
            f"{row['rms_bending']:.3e} | {row['max_convergence']:.3e} | {row['rms_convergence']:.3e} | "
            f"{row['max_shear']:.3e} | {row['rms_shear']:.3e} | {row['max_conservation']:.2e} |"
        )
    summary_lines += [
        "",
        "## Per-dataset products",
        "",
    ]
    for dataset in DATASETS:
        ds_dir = out_root / dataset["id"]
        summary_lines.append(f"### {dataset['label']} (`{dataset['id']}`)")
        summary_lines.append("")
        summary_lines.append(
            f"- CSVs: `matter.csv`, `constitutive.csv`, `gradient_*.csv`, `response_*.csv`, "
            f"`convergence.csv`, `shear_g1.csv`, `shear_g2.csv`, `deflection_*.csv`, "
            f"`magnification.csv`, `photon_*.csv`"
        )
        summary_lines.append(
            f"- Maps: `matter_map.png`, `constitutive_map.png`, `gradient_map.png`, "
            f"`response_magnitude_map.png`, `response_direction_map.png`, "
            f"`convergence_map.png`, `shear_g1_map.png`, `shear_g2_map.png`, "
            f"`shear_magnitude_map.png`, `deflection_magnitude_map.png`, `magnification_map.png`"
        )
        summary_lines.append(
            f"- Trajectories: `photon_trajectories.png`; composite: `composite_observables.png`"
        )
        summary_lines.append(f"- Statistics: `dataset_statistics.json`")
        summary_lines.append("")
    summary_lines += [
        "## Notes",
        "",
        "All datasets were processed by the identical frozen Version A pipeline. No "
        "parameter was altered between datasets. No comparison with ΛCDM, GR, or "
        "observations was performed at this milestone.",
    ]
    (out_root / "report.md").write_text("\n".join(summary_lines))
    (out_root / "run.json").write_text(json.dumps({
        "milestone": "PBUF WEAK-LENSING-GENERALIZATION-001",
        "status": "OK",
        "datasets": [d["id"] for d in DATASETS],
        "identical_pipeline": executable_hashes,
        "execution_seconds": time.perf_counter() - started,
    }, indent=2))
    (out_root / "validation.json").write_text(json.dumps({
        "milestone": "PBUF WEAK-LENSING-GENERALIZATION-001",
        "frozen_pipeline": {
            "constitutive": "Version A: C(X) = 0.18 * rho(X) / rho_max",
            "transport": "90-degree transverse response, direct addition + renormalisation",
            "amplitude": "A = |grad C|",
        },
        "datasets": [
            {"id": d["id"], "label": d["label"], "statistics": per_dataset_statistics[d["id"]]}
            for d in DATASETS
        ],
        "identical_pipeline_hashes": executable_hashes,
        "execution_seconds": time.perf_counter() - started,
    }, indent=2))
    all_pass = all(row["finite_outputs"] and row["max_conservation"] <= 1e-12
                    for row in cross_rows)
    print(json.dumps({
        "milestone": "PBUF WEAK-LENSING-GENERALIZATION-001",
        "status": "PASS" if all_pass else "FAIL",
        "datasets_run": len(cross_rows),
        "output": str(out_root),
        "identical_pipeline_sha256": executable_hashes,
    }, indent=2))
    return 0 if all_pass else 1


if __name__ == "__main__":
    raise SystemExit(main())
