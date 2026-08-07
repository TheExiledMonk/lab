#!/usr/bin/env python3
"""PBUF TRANSPORT-LAB-004 — Continuous angular sweep of the local response.

Parameterize:
    r(θ) = cos(θ) * ĝ + sin(θ) * ĝ_perp
where ĝ = ∇C / |∇C| and ĝ_perp is the 90° rotation of ĝ.
The response magnitude is |∇C| (matches TRANSPORT-LAB-001/2/3 scaling).

θ = 0°   : pure gradient-following
θ = 90°  : pure transverse
θ = 180° : anti-parallel to gradient

Same kernel, same inputs, same normalization as TRANSPORT-LAB-003.
Only the response direction varies.
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
from transport_lab003 import propagate_diag


ROOT = Path(__file__).resolve().parent
DEFAULT_OUT = ROOT / "runs" / "transport_lab004"


COARSE_ANGLES = [0, 15, 30, 45, 60, 75, 90, 105, 120, 135, 150, 165, 180]
FINE_RANGE = (80, 100)  # inclusive


def make_response(theta_deg: float):
    theta = np.radians(theta_deg)
    cos_t = float(np.cos(theta))
    sin_t = float(np.sin(theta))

    def fn(ix, iy, fi: FrozenInputs, state):
        gx = float(fi.gradient_x[iy, ix])
        gy = float(fi.gradient_y[iy, ix])
        g_mag = float(np.hypot(gx, gy))
        if g_mag < 1e-15:
            return 0.0, 0.0, "subtract"
        gx_h = gx / g_mag
        gy_h = gy / g_mag
        gpx_h = -gy_h
        gpy_h = gx_h
        rx = cos_t * gx_h + sin_t * gpx_h
        ry = cos_t * gy_h + sin_t * gpy_h
        return rx * g_mag, ry * g_mag, "subtract"

    return fn


@dataclass
class M:
    angle: float
    bend_max: float
    bend_mean: float
    conservation: float
    speed_drift_pre_max: float
    speed_drift_pre_mean: float
    direction_drift_mean: float
    position_error: float
    stable: bool
    runtime: float


def measure_one(theta_deg, fi):
    paths, diag, runtime = propagate_diag(fi, make_response(theta_deg))
    devs = np.array([float(np.max(np.abs(y - y[0]))) for x, y in paths])
    pre = np.concatenate(diag["speed_drift_pre"])
    dirs = np.concatenate(diag["direction_drift"])
    errs = np.concatenate(diag["pos_error"])
    finite = all(np.isfinite(np.concatenate([p[0], p[1]])).all() for p in paths)
    stable = bool(finite and pre.max() < 1.0)
    return M(
        angle=theta_deg,
        bend_max=float(devs.max()),
        bend_mean=float(devs.mean()),
        conservation=float(pre.max()),
        speed_drift_pre_max=float(pre.max()),
        speed_drift_pre_mean=float(pre.mean()),
        direction_drift_mean=float(dirs.mean()),
        position_error=float(errs.sum()),
        stable=stable,
        runtime=runtime,
    )


def run_sweep(fi, angles):
    out = []
    for a in angles:
        out.append(measure_one(a, fi))
    return out


def detect_optima(coarse):
    """Numerical detection of the maximum-bending and best-conservation angles.

    The coarse sweep covers [0, 180]°. Maximum bending may occur twice
    (near 90° and near -90° = 270° ≡ -90°). We report the angle in [0, 180]
    where the value is highest (or, if there's a unique maximum, that angle).
    """
    bends = np.array([m.bend_max for m in coarse])
    conservs = np.array([m.conservation for m in coarse])
    angles = np.array([m.angle for m in coarse])
    idx_max_bend = int(np.argmax(bends))
    idx_min_cons = int(np.argmin(conservs))
    return angles[idx_max_bend], angles[idx_min_cons]


def make_plots(coarse, fine, out: Path):
    angles_c = np.array([m.angle for m in coarse])
    bends_c = np.array([m.bend_max for m in coarse])
    conservs_c = np.array([m.conservation for m in coarse])
    pre_c = np.array([m.speed_drift_pre_max for m in coarse])
    dirs_c = np.array([m.direction_drift_mean for m in coarse])

    fig, axes = plt.subplots(2, 2, figsize=(12, 9))
    axes[0, 0].plot(angles_c, bends_c, "o-", color="C0")
    axes[0, 0].set_xlabel("θ (deg)")
    axes[0, 0].set_ylabel("bend_max")
    axes[0, 0].set_title("Angle vs bending")
    axes[0, 0].grid(True, alpha=0.3)
    axes[0, 0].axvline(90, color="grey", ls="--", alpha=0.5)

    axes[0, 1].plot(angles_c, conservs_c, "s-", color="C1")
    axes[0, 1].set_xlabel("θ (deg)")
    axes[0, 1].set_ylabel("conservation residual (pre-renorm max)")
    axes[0, 1].set_title("Angle vs conservation")
    axes[0, 1].grid(True, alpha=0.3)
    axes[0, 1].axvline(90, color="grey", ls="--", alpha=0.5)

    axes[1, 0].plot(angles_c, pre_c, "^-", color="C2")
    axes[1, 0].set_xlabel("θ (deg)")
    axes[1, 0].set_ylabel("speed drift pre-renorm (max)")
    axes[1, 0].set_title("Angle vs speed drift")
    axes[1, 0].grid(True, alpha=0.3)
    axes[1, 0].axvline(90, color="grey", ls="--", alpha=0.5)

    axes[1, 1].plot(angles_c, dirs_c, "d-", color="C3")
    axes[1, 1].set_xlabel("θ (deg)")
    axes[1, 1].set_ylabel("direction drift (mean)")
    axes[1, 1].set_title("Angle vs direction drift")
    axes[1, 1].grid(True, alpha=0.3)
    axes[1, 1].axvline(90, color="grey", ls="--", alpha=0.5)

    fig.tight_layout()
    fig.savefig(out / "angle_sweep_coarse.png", dpi=140)
    plt.close(fig)

    if fine:
        angles_f = np.array([m.angle for m in fine])
        bends_f = np.array([m.bend_max for m in fine])
        conservs_f = np.array([m.conservation for m in fine])
        pre_f = np.array([m.speed_drift_pre_max for m in fine])
        dirs_f = np.array([m.direction_drift_mean for m in fine])

        fig2, axes2 = plt.subplots(2, 2, figsize=(12, 9))
        axes2[0, 0].plot(angles_f, bends_f, "o-", color="C0")
        axes2[0, 0].set_xlabel("θ (deg)")
        axes2[0, 0].set_ylabel("bend_max")
        axes2[0, 0].set_title(f"Fine sweep {FINE_RANGE[0]}°–{FINE_RANGE[1]}°: bending")
        axes2[0, 0].grid(True, alpha=0.3)
        axes2[0, 0].axvline(90, color="grey", ls="--", alpha=0.5)

        axes2[0, 1].plot(angles_f, conservs_f, "s-", color="C1")
        axes2[0, 1].set_xlabel("θ (deg)")
        axes2[0, 1].set_ylabel("conservation residual")
        axes2[0, 1].set_title(f"Fine sweep: conservation")
        axes2[0, 1].grid(True, alpha=0.3)
        axes2[0, 1].axvline(90, color="grey", ls="--", alpha=0.5)

        axes2[1, 0].plot(angles_f, pre_f, "^-", color="C2")
        axes2[1, 0].set_xlabel("θ (deg)")
        axes2[1, 0].set_ylabel("speed drift pre-renorm (max)")
        axes2[1, 0].set_title(f"Fine sweep: speed drift")
        axes2[1, 0].grid(True, alpha=0.3)
        axes2[1, 0].axvline(90, color="grey", ls="--", alpha=0.5)

        axes2[1, 1].plot(angles_f, dirs_f, "d-", color="C3")
        axes2[1, 1].set_xlabel("θ (deg)")
        axes2[1, 1].set_ylabel("direction drift (mean)")
        axes2[1, 1].set_title(f"Fine sweep: direction drift")
        axes2[1, 1].grid(True, alpha=0.3)
        axes2[1, 1].axvline(90, color="grey", ls="--", alpha=0.5)

        fig2.tight_layout()
        fig2.savefig(out / "angle_sweep_fine.png", dpi=140)
        plt.close(fig2)


def write_table(out: Path, coarse, fine):
    with (out / "sweep_table.csv").open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["Angle (deg)", "Bend", "Conservation",
                    "Speed Drift", "Stable"])
        for m in coarse:
            w.writerow([f"{m.angle:.1f}",
                        f"{m.bend_max:.6e}",
                        f"{m.conservation:.6e}",
                        f"{m.speed_drift_pre_max:.6e}",
                        "yes" if m.stable else "no"])
        if fine:
            w.writerow(["--- fine sweep ---"] * 5)
            for m in fine:
                w.writerow([f"{m.angle:.1f}",
                            f"{m.bend_max:.6e}",
                            f"{m.conservation:.6e}",
                            f"{m.speed_drift_pre_max:.6e}",
                            "yes" if m.stable else "no"])


def detect_plateau(coarse, key: str, tol: float = 0.05):
    """Detect whether the optimum is a broad plateau.

    Returns the (lo, hi) interval in degrees where the metric is within
    `tol` (fraction) of the best value, or None if no plateau detected.
    """
    vals = np.array([getattr(m, key) for m in coarse])
    angles = np.array([m.angle for m in coarse])
    if key in ("bend_max", "direction_drift_mean", "position_error"):
        # Larger is better? For bending, we report where the metric is at its
        # maximum. The plateau is where the metric stays within tol of the max.
        best = float(vals.max())
        threshold = best * (1 - tol)
        mask = vals >= threshold
    else:
        # Smaller is better (conservation, speed drift)
        best = float(vals.min())
        threshold = best * (1 + tol)
        mask = vals <= threshold
    if not mask.any():
        return None
    return float(angles[mask].min()), float(angles[mask].max())


def combined_optimum(coarse):
    """Combined optimum: rank by normalized bending + normalized conservation.

    Each metric is normalised to [0, 1] across the sweep (1 = best).
    The combined score is the sum of the two normalised scores.
    """
    bends = np.array([m.bend_max for m in coarse])
    conservs = np.array([m.conservation for m in coarse])
    b_max = bends.max()
    c_min = conservs.min()
    b_norm = bends / b_max  # 1 at max bend
    c_norm = c_min / np.where(conservs > 0, conservs, 1e-30)  # 1 at min conservation
    score = b_norm + c_norm
    idx = int(np.argmax(score))
    return coarse[idx].angle, score[idx], b_norm[idx], c_norm[idx]


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--output", type=Path, default=DEFAULT_OUT)
    p.add_argument("--no-fine", action="store_true",
                   help="skip the fine sweep near 90 degrees")
    a = p.parse_args()
    a.output.mkdir(parents=True, exist_ok=True)
    fi = load_inputs()

    print(f"Coarse sweep: {COARSE_ANGLES}")
    coarse = run_sweep(fi, COARSE_ANGLES)
    for m in coarse:
        print(f"  θ={m.angle:6.1f}°  bend={m.bend_max:.4e}  "
              f"cons={m.conservation:.4e}  speed={m.speed_drift_pre_max:.4e}  "
              f"stable={m.stable}")

    fine = []
    if not a.no_fine:
        fine_angles = list(range(FINE_RANGE[0], FINE_RANGE[1] + 1))
        print(f"\nFine sweep: {fine_angles}")
        fine = run_sweep(fi, fine_angles)
        for m in fine:
            print(f"  θ={m.angle:6.1f}°  bend={m.bend_max:.4e}  "
                  f"cons={m.conservation:.4e}  speed={m.speed_drift_pre_max:.4e}  "
                  f"stable={m.stable}")

    # Persist
    rows = [asdict(m) for m in coarse]
    if fine:
        rows.extend([asdict(m) for m in fine])
    keys = list(rows[0].keys())
    with (a.output / "measurements.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        for r in rows:
            w.writerow(r)
    (a.output / "measurements.json").write_text(json.dumps(rows, indent=2))

    # Table
    write_table(a.output, coarse, fine)

    # Plots
    make_plots(coarse, fine, a.output)

    # Optima
    max_bend_angle, min_cons_angle = detect_optima(coarse)
    bend_plateau = detect_plateau(coarse, "bend_max", tol=0.05)
    cons_plateau = detect_plateau(coarse, "conservation", tol=0.05)
    combined_angle, combined_score, b_norm, c_norm = combined_optimum(coarse)

    summary = {
        "coarse_max_bending_angle_deg": float(max_bend_angle),
        "coarse_min_conservation_angle_deg": float(min_cons_angle),
        "coarse_combined_optimum_angle_deg": float(combined_angle),
        "coarse_combined_score": float(combined_score),
        "coarse_combined_bending_normalised": float(b_norm),
        "coarse_combined_conservation_normalised": float(c_norm),
        "bending_plateau_within_5pct": (
            [float(bend_plateau[0]), float(bend_plateau[1])]
            if bend_plateau else None
        ),
        "conservation_plateau_within_5pct": (
            [float(cons_plateau[0]), float(cons_plateau[1])]
            if cons_plateau else None
        ),
        "fine_sweep": (
            {"performed": True,
             "range_deg": [FINE_RANGE[0], FINE_RANGE[1]],
             "argmax_bend_deg": float(max(m.angle for m in fine
                                          if m.bend_max ==
                                          max(fm.bend_max for fm in fine))),
             "argmin_conservation_deg": float(min(m.angle for m in fine
                                                  if m.conservation ==
                                                  min(fm.conservation for fm in fine)))}
            if fine else {"performed": False}
        ),
    }
    (a.output / "optima.json").write_text(json.dumps(summary, indent=2))

    print("\n=== Optima (coarse) ===")
    print(f"  Max bending angle:   θ = {max_bend_angle}°")
    print(f"  Min conservation:    θ = {min_cons_angle}°")
    print(f"  Combined optimum:    θ = {combined_angle}°")
    if bend_plateau:
        print(f"  Bending plateau:     {bend_plateau[0]}°–{bend_plateau[1]}°  (within 5%)")
    if cons_plateau:
        print(f"  Conservation plateau: {cons_plateau[0]}°–{cons_plateau[1]}°  (within 5%)")

    print(f"\nArtefacts written to {a.output}/")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
