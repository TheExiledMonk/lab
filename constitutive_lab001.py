#!/usr/bin/env python3
"""PBUF CONSTITUTIVE-LAB-001 — Constitutive-field generation rules.

Frozen transport (from TRANSPORT-LAB-008):
  - neighbour-to-neighbour transport
  - 90° transverse response
  - linear amplitude law A = |∇C|
  - direct addition + renormalisation update
  - identical kernel, integration, normalisation, timestep

Variable: the rule that turns matter(x) into C(x).
"""
from __future__ import annotations

import argparse
import csv
import json
from dataclasses import dataclass, asdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from transport_lab001 import FrozenInputs, load_inputs
from transport_lab007 import (
    upd_direct_addition, propagate_with_update,
)


ROOT = Path(__file__).resolve().parent
DEFAULT_OUT = ROOT / "runs" / "constitutive_lab001"


# ----------------------------------------------------------------------------
# Candidate generators: matter -> C (un-normalised).
# ----------------------------------------------------------------------------

def cand_existing_control(matter, spacing):
    """Candidate 1: existing deformation field (already computed).

    Loaded directly from runs/wl001/deformation.csv.
    """
    return None  # handled specially in main()


def cand_local_linear(matter, spacing):
    """Candidate 2: purely local, no kernel.
    C(x) = strength * matter(x). After peak-normalisation this is
    proportional to matter, i.e. mathematically the same shape as Candidate 1.
    """
    return matter.copy()


def _build_kernel(n: int, extent: float, kernel_fn, max_radius=None):
    spacing = 2 * extent / (n - 1)
    i, j = np.meshgrid(np.arange(n), np.arange(n), indexing="ij")
    i0 = j0 = n // 2
    r = np.sqrt((i - i0) ** 2 + (j - j0) ** 2) * spacing
    k = kernel_fn(r)
    if max_radius is not None:
        k = np.where(r <= max_radius, k, 0.0)
    # Zero DC bias via subtracting the mean (so convolution doesn't
    # create a uniform background). We will renormalise anyway.
    return k


def _convolve(matter, kernel):
    return np.real(np.fft.ifft2(np.fft.fft2(matter) * np.fft.fft2(kernel)))


def cand_finite_range(matter, spacing, radius: float = 1.5):
    """Candidate 3: finite-range boxcar kernel.

    Average matter within a radius (in grid units). Outside the radius,
    C = 0.
    """
    n = matter.shape[0]
    extent = 8.0
    k = _build_kernel(n, extent, lambda r: np.ones_like(r), max_radius=radius)
    return _convolve(matter, k)


def cand_exponential(matter, spacing, length: float = 1.0):
    """Candidate 4: exponential decay kernel."""
    n = matter.shape[0]
    extent = 8.0
    k = _build_kernel(n, extent, lambda r: np.exp(-r / length))
    return _convolve(matter, k)


def cand_gaussian(matter, spacing, sigma: float = 1.0):
    """Candidate 5: Gaussian kernel."""
    n = matter.shape[0]
    extent = 8.0
    k = _build_kernel(n, extent, lambda r: np.exp(-r ** 2 / (2 * sigma ** 2)))
    return _convolve(matter, k)


def cand_inverse_distance(matter, spacing, length: float = 1.0):
    """Candidate 6: inverse-distance kernel."""
    n = matter.shape[0]
    extent = 8.0
    def kfn(r):
        out = np.where(r > 0, 1.0 / (1.0 + r / length), 1.0 / length)
        return out
    k = _build_kernel(n, extent, kfn)
    return _convolve(matter, k)


def cand_compact_support(matter, spacing, support: float = 2.0):
    """Candidate 7: compact-support (Wendland C2) kernel.
    K(r) = (1 - r/R)^4 * (4 r/R + 1) for r < R, else 0.
    """
    n = matter.shape[0]
    extent = 8.0
    R = support
    def kfn(r):
        rho = np.where(R > 0, r / R, 0.0)
        out = np.where(rho < 1.0, (1.0 - rho) ** 4 * (4.0 * rho + 1.0), 0.0)
        return out
    k = _build_kernel(n, extent, kfn, max_radius=R)
    return _convolve(matter, k)


# ----------------------------------------------------------------------------
# Driver.
# ----------------------------------------------------------------------------

@dataclass
class M:
    candidate: str
    label: str
    C_max: float
    C_mean: float
    gradC_max: float
    gradC_mean: float
    bend_max: float
    bend_mean: float
    conservation: float
    speed_drift_pre_max: float
    direction_drift_mean: float
    position_error: float
    stable: bool
    runtime: float


def _normalize_to_target(C, target_max):
    if C.max() <= 0:
        return np.zeros_like(C)
    return C * (target_max / C.max())


def stats(C, gradC):
    return {
        "max_C": float(C.max()),
        "mean_C": float(C.mean()),
        "max_gradC": float(np.hypot(gradC[0], gradC[1]).max()),
        "mean_gradC": float(np.hypot(gradC[0], gradC[1]).mean()),
    }


def compute_gradC(C, fi):
    n = fi.n
    extent = fi.extent
    spacing = 2 * extent / (n - 1)
    gy, gx = np.gradient(C, spacing, spacing)
    return gx, gy


def run_candidate(name, label, C, fi, matter, spacing) -> M:
    import time
    gx, gy = compute_gradC(C, fi)
    A_field = np.hypot(gx, gy)
    g_norm = np.maximum(A_field, 1e-15)
    gx_hat = gx / g_norm
    gy_hat = gy / g_norm
    bad = A_field < 1e-15
    gx_hat = np.where(bad, 1.0, gx_hat)
    gy_hat = np.where(bad, 0.0, gy_hat)

    started = time.time()
    paths, diag, runtime = propagate_with_update(
        fi, A_field, gx_hat, gy_hat, upd_direct_addition)
    actual_runtime = time.time() - started

    devs = np.array([float(np.max(np.abs(y - y[0]))) for x, y in paths])
    pre = np.concatenate(diag["speed_drift_pre"])
    dirs = np.concatenate(diag["direction_drift"])
    errs = np.concatenate(diag["pos_error"])
    finite = all(np.isfinite(np.concatenate([p[0], p[1]])).all() for p in paths)
    stable = bool(finite and pre.max() < 1.0)
    s = stats(C, (gx, gy))
    return M(
        candidate=name, label=label,
        C_max=s["max_C"], C_mean=s["mean_C"],
        gradC_max=s["max_gradC"], gradC_mean=s["mean_gradC"],
        bend_max=float(devs.max()), bend_mean=float(devs.mean()),
        conservation=float(pre.max()),
        speed_drift_pre_max=float(pre.max()),
        direction_drift_mean=float(dirs.mean()),
        position_error=float(errs.sum()),
        stable=stable, runtime=runtime,
    )


CANDIDATES = [
    ("Cand 1", "Existing control (linear local)", cand_existing_control),
    ("Cand 2", "Local linear (no kernel)", cand_local_linear),
    ("Cand 3", "Finite-range boxcar (R = 1.5)", lambda m, s: cand_finite_range(m, s, 1.5)),
    ("Cand 4", "Exponential kernel (L = 1.0)", lambda m, s: cand_exponential(m, s, 1.0)),
    ("Cand 5", "Gaussian kernel (σ = 1.0)", lambda m, s: cand_gaussian(m, s, 1.0)),
    ("Cand 6", "Inverse-distance (L = 1.0)", lambda m, s: cand_inverse_distance(m, s, 1.0)),
    ("Cand 7", "Compact-support (R = 2.0)", lambda m, s: cand_compact_support(m, s, 2.0)),
]


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--output", type=Path, default=DEFAULT_OUT)
    a = p.parse_args()
    a.output.mkdir(parents=True, exist_ok=True)
    fi = load_inputs()
    matter = fi.matter
    spacing = 2 * fi.extent / (fi.n - 1)
    control_C = fi.deformation
    control_max = float(control_C.max())
    print(f"Control max C = {control_max:.4e}")

    measurements = []
    C_fields = {}

    for name, label, fn in CANDIDATES:
        if fn is cand_existing_control:
            C = control_C.copy()
        else:
            C = fn(matter, spacing)
            C = _normalize_to_target(C, control_max)
        C_fields[name] = C
        m = run_candidate(name, label, C, fi, matter, spacing)
        measurements.append(m)
        print(f"  {name}: C_max={m.C_max:.4e}  gradC_max={m.gradC_max:.4e}  "
              f"bend={m.bend_max:.4e}  cons={m.conservation:.4e}")

    rows = [asdict(m) for m in measurements]
    keys = list(rows[0].keys())
    with (a.output / "measurements.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        for r in rows:
            w.writerow(r)
    (a.output / "measurements.json").write_text(json.dumps(rows, indent=2))

    # Validation: every candidate has max(C) == control max
    validation = []
    all_pass = True
    for m in measurements:
        ok = abs(m.C_max - control_max) / control_max < 1e-10
        validation.append((m.candidate, m.label, m.C_max, m.C_mean, ok))
        if not ok:
            all_pass = False
    with (a.output / "validation.csv").open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["Candidate", "Label", "Max C", "Mean C", "Pass/Fail"])
        for v in validation:
            w.writerow([v[0], v[1], f"{v[2]:.4e}", f"{v[3]:.4e}",
                        "PASS" if v[4] else "FAIL"])

    if not all_pass:
        print("\nExperiment invalid: candidate constitutive fields not normalised.")
        return 1

    # Constitutive statistics table
    with (a.output / "constitutive_statistics.csv").open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["Candidate", "Max C", "Mean C", "Max |∇C|", "Mean |∇C|"])
        for m in measurements:
            w.writerow([m.candidate,
                        f"{m.C_max:.4e}", f"{m.C_mean:.4e}",
                        f"{m.gradC_max:.4e}", f"{m.gradC_mean:.4e}"])

    # Weak-lensing performance table
    with (a.output / "weak_lensing_performance.csv").open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["Candidate", "Bend", "Conservation", "Runtime (s)"])
        for m in measurements:
            w.writerow([m.candidate,
                        f"{m.bend_max:.4e}",
                        f"{m.conservation:.4e}",
                        f"{m.runtime:.3f}"])

    # Relative comparison
    control = measurements[0]
    rel = []
    for m in measurements:
        bend_pct = (m.bend_max - control.bend_max) / control.bend_max * 100
        cons_pct = (m.conservation - control.conservation) / control.conservation * 100 if control.conservation > 0 else 0
        grad_max_pct = (m.gradC_max - control.gradC_max) / control.gradC_max * 100
        grad_mean_pct = (m.gradC_mean - control.gradC_mean) / control.gradC_mean * 100
        rel.append({"candidate": m.candidate, "label": m.label,
                    "bend_pct": bend_pct, "cons_pct": cons_pct,
                    "grad_max_pct": grad_max_pct, "grad_mean_pct": grad_mean_pct})
    with (a.output / "relative_vs_control.csv").open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["Candidate", "Label", "Bend Δ%", "Conservation Δ%",
                    "|∇C| max Δ%", "|∇C| mean Δ%"])
        for r in rel:
            w.writerow([r["candidate"], r["label"],
                        f"{r['bend_pct']:+.4f}%",
                        f"{r['cons_pct']:+.4f}%",
                        f"{r['grad_max_pct']:+.4f}%",
                        f"{r['grad_mean_pct']:+.4f}%"])

    # Indistinguishable group (5% on bend_max)
    indistinguishable = []
    for m in measurements:
        if abs(m.bend_max - control.bend_max) / control.bend_max < 0.05:
            indistinguishable.append(m.candidate)
    (a.output / "indistinguishable.json").write_text(json.dumps(
        {"tolerance": "5% on bend_max vs control",
         "candidates_within_tolerance": indistinguishable}, indent=2))

    # Visualisations
    cmap = "viridis"
    fig, axes = plt.subplots(2, 4, figsize=(16, 8))
    extent_xy = [-fi.extent, fi.extent, -fi.extent, fi.extent]
    for ax, name in zip(axes[0], [c[0] for c in CANDIDATES[:4]]):
        C = C_fields[name]
        im = ax.imshow(C, origin="lower", extent=extent_xy, cmap=cmap)
        ax.set_title(f"{name}\nC (constitutive)")
        fig.colorbar(im, ax=ax, shrink=.7)
    for ax, name in zip(axes[1], [c[0] for c in CANDIDATES[:4]]):
        C = C_fields[name]
        gx, gy = compute_gradC(C, fi)
        im = ax.imshow(np.hypot(gx, gy), origin="lower", extent=extent_xy, cmap=cmap)
        ax.set_title(f"{name}\n|∇C|")
        fig.colorbar(im, ax=ax, shrink=.7)
    fig.tight_layout()
    fig.savefig(a.output / "constitutive_fields_1to4.png", dpi=130)
    plt.close(fig)

    fig, axes = plt.subplots(2, 3, figsize=(12, 8))
    for ax, name in zip(axes[0], [c[0] for c in CANDIDATES[4:]]):
        C = C_fields[name]
        im = ax.imshow(C, origin="lower", extent=extent_xy, cmap=cmap)
        ax.set_title(f"{name}\nC (constitutive)")
        fig.colorbar(im, ax=ax, shrink=.7)
    for ax, name in zip(axes[1], [c[0] for c in CANDIDATES[4:]]):
        C = C_fields[name]
        gx, gy = compute_gradC(C, fi)
        im = ax.imshow(np.hypot(gx, gy), origin="lower", extent=extent_xy, cmap=cmap)
        ax.set_title(f"{name}\n|∇C|")
        fig.colorbar(im, ax=ax, shrink=.7)
    fig.tight_layout()
    fig.savefig(a.output / "constitutive_fields_5to7.png", dpi=130)
    plt.close(fig)

    print(f"\n=== Indistinguishable from control (5% on bend_max) ===")
    for c in indistinguishable:
        print(f"  {c}")
    print(f"\nArtefacts written to {a.output}/")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())