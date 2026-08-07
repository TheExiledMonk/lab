#!/usr/bin/env python3
"""PBUF TRANSPORT-LAB-006 — Pure direction ablation (magnitude-normalised).

Critical fix vs LAB-005: every candidate produces a response with
exactly the same magnitude A = |∇C| (the control's magnitude field).
Only the response direction varies.

For each candidate direction field g_i:

    g_hat_i = g_i / |g_i|              (unit direction; safe at |g|≈0)
    r_i     = A * R_90(g_hat_i)         (A = |∇C| everywhere)

R_90(x, y) = (-y, x).

Each candidate thus has |r_i| = A = |∇C| at every cell — verified
numerically before any propagation runs.
"""
from __future__ import annotations

import argparse
import csv
import json
from dataclasses import dataclass, asdict
from pathlib import Path

import numpy as np

from transport_lab001 import FrozenInputs, load_inputs
from transport_lab003 import propagate_diag


ROOT = Path(__file__).resolve().parent
DEFAULT_OUT = ROOT / "runs" / "transport_lab006"


# ----------------------------------------------------------------------------
# Candidate direction fields (unit vectors, |g|=1 wherever |g|>0).
# ----------------------------------------------------------------------------

def cand_gradC(fi: FrozenInputs):
    gx, gy = fi.gradient_x, fi.gradient_y
    return _unit(gx, gy)


def cand_stress_grad(fi: FrozenInputs):
    n = fi.n
    extent = fi.extent
    spacing = 2 * extent / (n - 1)
    P = np.hypot(fi.gradient_x, fi.gradient_y)
    gx_P = np.zeros_like(P)
    gy_P = np.zeros_like(P)
    gx_P[:, 1:-1] = (P[:, 2:] - P[:, :-2]) / (2 * spacing)
    gy_P[1:-1, :] = (P[2:, :] - P[:-2, :]) / (2 * spacing)
    eps = 1e-3
    gx = gx_P / (fi.deformation + eps)
    gy = gy_P / (fi.deformation + eps)
    return _unit(gx, gy)


def cand_energy_grad(fi: FrozenInputs):
    gx = fi.deformation * fi.gradient_x
    gy = fi.deformation * fi.gradient_y
    return _unit(gx, gy)


def cand_traction(fi: FrozenInputs):
    n = fi.n
    extent = fi.extent
    x = np.linspace(-extent, extent, n)
    X, Y = np.meshgrid(x, x, indexing="xy")
    rx = X - fi.mass_x
    ry = Y - fi.mass_y
    r = np.hypot(rx, ry)
    r_safe = np.where(r > 1e-12, r, 1e-12)
    Nx = np.where(r > 1e-12, rx / r_safe, 0.0)
    Ny = np.where(r > 1e-12, ry / r_safe, 0.0)
    return Nx, Ny  # already unit magnitude


def cand_force_density(fi: FrozenInputs):
    n = fi.n
    extent = fi.extent
    spacing = 2 * extent / (n - 1)
    C = fi.deformation
    lap = np.zeros_like(C)
    lap[:, 1:-1] = (C[:, 2:] - 2 * C[:, 1:-1] + C[:, :-2]) / spacing**2
    lap[1:-1, :] += (C[2:, :] - 2 * C[1:-1, :] + C[:-2, :]) / spacing**2
    gx = np.zeros_like(lap)
    gy = np.zeros_like(lap)
    gx[:, 1:-1] = (lap[:, 2:] - lap[:, :-2]) / (2 * spacing)
    gy[1:-1, :] = (lap[2:, :] - lap[:-2, :]) / (2 * spacing)
    return _unit(gx, gy)


def _principal_eigenvector_unit(fi: FrozenInputs, F: np.ndarray, G: np.ndarray):
    n = fi.n
    extent = fi.extent
    spacing = 2 * extent / (n - 1)
    gx_F = np.zeros_like(F)
    gy_F = np.zeros_like(F)
    gx_F[:, 1:-1] = (F[:, 2:] - F[:, :-2]) / (2 * spacing)
    gy_F[1:-1, :] = (F[2:, :] - F[:-2, :]) / (2 * spacing)

    Hxx = np.zeros_like(F)
    Hyy = np.zeros_like(F)
    Hxy = np.zeros_like(F)
    Hxx[:, 1:-1] = (G[:, 2:] - 2 * G[:, 1:-1] + G[:, :-2]) / spacing**2
    Hyy[1:-1, :] = (G[2:, :] - 2 * G[1:-1, :] + G[:-2, :]) / spacing**2
    Hxy[1:-1, 1:-1] = (G[2:, 2:] - G[2:, :-2] - G[:-2, 2:] + G[:-2, :-2]) / (4 * spacing**2)

    alpha = 1.0
    Txx = gx_F * gx_F + alpha * Hxx
    Tyy = gy_F * gy_F + alpha * Hyy
    Txy = gx_F * gy_F + alpha * Hxy
    trace = Txx + Tyy
    det = Txx * Tyy - Txy * Txy
    disc = np.sqrt(np.maximum(trace * trace - 4 * det, 0.0))
    lam_max = (trace + disc) / 2.0
    vx = Txy
    vy = lam_max - Txx
    return _unit(vx, vy)


def cand_principal_strain(fi: FrozenInputs):
    return _principal_eigenvector_unit(fi, fi.deformation, fi.deformation)


def cand_principal_stress(fi: FrozenInputs):
    W = 0.5 * fi.deformation ** 2
    return _principal_eigenvector_unit(fi, W, W)


def _unit(gx, gy):
    mag = np.hypot(gx, gy)
    safe = mag > 1e-15
    ux = np.where(safe, gx / np.where(safe, mag, 1.0), 1.0)
    uy = np.where(safe, gy / np.where(safe, mag, 1.0), 0.0)
    return ux, uy


# ----------------------------------------------------------------------------
# Build response = A * R_90(g_hat), with A = |∇C| (the control magnitude).
# ----------------------------------------------------------------------------

def make_response_magnitude_normalised(gx_hat, gy_hat, A_field):
    def fn(ix, iy, fi: FrozenInputs, state):
        A = float(A_field[iy, ix])
        ux = float(gx_hat[iy, ix])
        uy = float(gy_hat[iy, ix])
        return -A * uy, A * ux, "subtract"
    return fn


# ----------------------------------------------------------------------------
# Driver.
# ----------------------------------------------------------------------------

@dataclass
class M:
    candidate: str
    reference_label: str
    bend_max: float
    bend_mean: float
    conservation: float
    speed_drift_pre_max: float
    speed_drift_pre_mean: float
    direction_drift_mean: float
    position_error: float
    stable: bool
    runtime: float


def validate_magnitude(gx_hat, gy_hat, A_field, label: str) -> tuple[bool, float, float]:
    """Confirm |r| = A everywhere (to machine precision).

    r = A * R_90(g_hat), |R_90(g_hat)| = |g_hat| = 1 wherever |g_hat| is
    well-defined. So |r| = A wherever the unit direction is valid.
    """
    mag_r = A_field * np.hypot(gx_hat, gy_hat)
    safe = A_field > 0
    diff = np.abs(mag_r[safe] - A_field[safe])
    mean_diff = float(diff.mean()) if diff.size else 0.0
    max_diff = float(diff.max()) if diff.size else 0.0
    return max_diff < 1e-12, mean_diff, max_diff


def measure(cand, label, gx_hat, gy_hat, A_field, fi) -> M:
    response = make_response_magnitude_normalised(gx_hat, gy_hat, A_field)
    paths, diag, runtime = propagate_diag(fi, response)
    devs = np.array([float(np.max(np.abs(y - y[0]))) for x, y in paths])
    pre = np.concatenate(diag["speed_drift_pre"])
    dirs = np.concatenate(diag["direction_drift"])
    errs = np.concatenate(diag["pos_error"])
    finite = all(np.isfinite(np.concatenate([p[0], p[1]])).all() for p in paths)
    stable = bool(finite and pre.max() < 1.0)
    return M(
        candidate=cand, reference_label=label,
        bend_max=float(devs.max()), bend_mean=float(devs.mean()),
        conservation=float(pre.max()),
        speed_drift_pre_max=float(pre.max()),
        speed_drift_pre_mean=float(pre.mean()),
        direction_drift_mean=float(dirs.mean()),
        position_error=float(errs.sum()),
        stable=stable, runtime=runtime,
    )


CANDIDATES = [
    ("Candidate 1", "Constitutive gradient ∇C",                    cand_gradC),
    ("Candidate 2", "Stress gradient ∇|∇C|/(C+ε)",                 cand_stress_grad),
    ("Candidate 3", "Energy gradient C·∇C",                        cand_energy_grad),
    ("Candidate 4", "Traction (∇C·N̂)N̂",                            cand_traction),
    ("Candidate 5", "Force density ∇(∇²C)",                        cand_force_density),
    ("Candidate 6", "Principal strain direction",                   cand_principal_strain),
    ("Candidate 7", "Principal stress direction",                   cand_principal_stress),
]


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--output", type=Path, default=DEFAULT_OUT)
    a = p.parse_args()
    a.output.mkdir(parents=True, exist_ok=True)
    fi = load_inputs()

    A_field = np.hypot(fi.gradient_x, fi.gradient_y)

    # ---- Validation pass -----------------------------------------------
    print("=== Magnitude validation ===")
    validation_rows = []
    all_pass = True
    for cand, label, fn in CANDIDATES:
        gx_hat, gy_hat = fn(fi)
        passed, mean_diff, max_diff = validate_magnitude(gx_hat, gy_hat, A_field, label)
        validation_rows.append({
            "candidate": cand, "label": label,
            "mean_|r_minus_A|": mean_diff,
            "max_|r_minus_A|": max_diff,
            "passed": passed,
        })
        status = "PASS" if passed else "FAIL"
        print(f"  {cand}: max|r-A| = {max_diff:.3e}  [{status}]")
        if not passed:
            all_pass = False

    with (a.output / "validation.csv").open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["Candidate", "Label", "Mean |r-A|", "Max |r-A|", "Pass/Fail"])
        for v in validation_rows:
            w.writerow([v["candidate"], v["label"],
                        f"{v['mean_|r_minus_A|']:.3e}",
                        f"{v['max_|r_minus_A|']:.3e}",
                        "PASS" if v["passed"] else "FAIL"])

    if not all_pass:
        print("\nExperiment invalid: response magnitude differs between candidates.")
        (a.output / "status.json").write_text(json.dumps(
            {"status": "INVALID", "reason": "magnitude mismatch"}, indent=2))
        return 1

    # ---- Propagation ---------------------------------------------------
    print("\n=== Propagation ===")
    measurements = []
    for cand, label, fn in CANDIDATES:
        gx_hat, gy_hat = fn(fi)
        m = measure(cand, label, gx_hat, gy_hat, A_field, fi)
        measurements.append(m)
        print(f"  {cand}: bend={m.bend_max:.4e}  cons={m.conservation:.4e}  "
              f"stable={m.stable}  runtime={m.runtime:.3f}s")

    rows = [asdict(m) for m in measurements]
    keys = list(rows[0].keys())
    with (a.output / "measurements.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        for r in rows:
            w.writerow(r)
    (a.output / "measurements.json").write_text(json.dumps(rows, indent=2))

    # Performance table
    with (a.output / "performance_table.csv").open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["Candidate", "Bend", "Conservation", "Speed Drift", "Stable"])
        for m in measurements:
            w.writerow([m.candidate,
                        f"{m.bend_max:.4e}",
                        f"{m.conservation:.4e}",
                        f"{m.speed_drift_pre_max:.4e}",
                        "yes" if m.stable else "no"])

    # Relative comparison vs control
    control = measurements[0]
    rel = []
    for m in measurements:
        bend_pct = (m.bend_max - control.bend_max) / control.bend_max * 100
        cons_pct = (m.conservation - control.conservation) / control.conservation * 100
        speed_pct = (m.speed_drift_pre_max - control.speed_drift_pre_max) / control.speed_drift_pre_max * 100
        rel.append({
            "candidate": m.candidate,
            "reference_label": m.reference_label,
            "bend_pct_vs_control": bend_pct,
            "cons_pct_vs_control": cons_pct,
            "speed_pct_vs_control": speed_pct,
        })
    with (a.output / "relative_vs_control.csv").open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["Candidate", "Reference Direction",
                    "Bend % vs control", "Conservation % vs control",
                    "Speed drift % vs control"])
        for r in rel:
            w.writerow([r["candidate"], r["reference_label"],
                        f"{r['bend_pct_vs_control']:+.4f}%",
                        f"{r['cons_pct_vs_control']:+.4f}%",
                        f"{r['speed_pct_vs_control']:+.4f}%"])

    # Indistinguishable group detection (within tolerance of control)
    bend_tol = 0.05  # 5% of control
    cons_tol = 0.05
    indistinguishable = []
    for m in measurements:
        bend_within = abs(m.bend_max - control.bend_max) / control.bend_max < bend_tol
        cons_within = abs(m.conservation - control.conservation) / control.conservation < cons_tol
        if bend_within and cons_within:
            indistinguishable.append(m.candidate)

    (a.output / "indistinguishable.json").write_text(json.dumps(
        {"tolerance": "5% on both bend_max and conservation",
         "candidates_within_tolerance": indistinguishable}, indent=2))

    # Output status
    status = {
        "validation": "all candidates pass (|r|=A within machine precision)",
        "control_candidate": "Candidate 1 (∇C)",
        "indistinguishable_from_control": indistinguishable,
        "max_bend_candidate": max(measurements, key=lambda m: m.bend_max).candidate,
        "min_conservation_candidate": min(measurements, key=lambda m: m.conservation).candidate,
    }
    (a.output / "status.json").write_text(json.dumps(status, indent=2))

    print(f"\n=== Indistinguishable from control (5% on bend AND cons) ===")
    for c in indistinguishable:
        print(f"  {c}")
    print(f"\nMax bend:   {status['max_bend_candidate']}")
    print(f"Min conservation: {status['min_conservation_candidate']}")
    print(f"\nArtefacts written to {a.output}/")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
