#!/usr/bin/env python3
"""PBUF TRANSPORT-LAB-008 — Local response amplitude law ablation.

Frozen: response direction (90° transverse), update rule (direct addition
with renormalisation), kernel, integration, timestep, normalisation, Lens-001.

Variable: the scalar amplitude law A(|∇C|). All candidates are normalised
so max(A) = max_control = max(|∇C|).
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
DEFAULT_OUT = ROOT / "runs" / "transport_lab008"


# ----------------------------------------------------------------------------
# Candidate amplitude laws.
# ----------------------------------------------------------------------------

def cand_linear(g, max_g):
    """Control: A = |∇C|. max(A) = max_g by construction."""
    return g.copy()


def cand_sqrt(g, max_g):
    """A = sqrt(|∇C|), normalised to max = max_g."""
    A = np.sqrt(np.maximum(g, 0.0))
    return _normalize(A, max_g)


def cand_quadratic(g, max_g):
    """A = |∇C|², normalised."""
    A = g ** 2
    return _normalize(A, max_g)


def cand_log(g, max_g):
    """A = log(1 + |∇C|), normalised."""
    A = np.log1p(np.maximum(g, 0.0))
    return _normalize(A, max_g)


def cand_saturating(g, max_g, beta=None):
    """A = |∇C| / (1 + β|∇C|), normalised.

    Default β = 5/max_g so that the function is well into saturation by
    the time |∇C| reaches its maximum, then rescaled to peak = max_g.
    """
    if beta is None:
        beta = 5.0 / max_g
    A = g / (1.0 + beta * g)
    return _normalize(A, max_g)


def cand_threshold(g, max_g, threshold_frac: float):
    """A = 0 below threshold, A = |∇C| above. max(A) = max_g by construction."""
    threshold = threshold_frac * max_g
    return np.where(g < threshold, 0.0, g)


def cand_piecewise(g, max_g, g_break_frac: float = 0.5):
    """A = |∇C| for g < g_break, A = g_break + 0.3·(g - g_break) above, normalised."""
    g_break = g_break_frac * max_g
    A = np.where(g <= g_break, g, g_break + 0.3 * (g - g_break))
    return _normalize(A, max_g)


def _normalize(A, max_g):
    """Rescale so max(A) = max_g."""
    A_max = float(A.max())
    if A_max < 1e-30:
        return np.zeros_like(A)
    return A * (max_g / A_max)


# ----------------------------------------------------------------------------
# Driver.
# ----------------------------------------------------------------------------

@dataclass
class M:
    candidate: str
    label: str
    bend_max: float
    bend_mean: float
    conservation: float
    speed_drift_pre_max: float
    speed_drift_pre_mean: float
    direction_drift_mean: float
    position_error: float
    stable: bool
    runtime: float
    A_max: float
    A_mean: float


def run_candidate(name, label, A_field, fi, gx_hat, gy_hat) -> M:
    import time
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
    return M(
        candidate=name, label=label,
        bend_max=float(devs.max()), bend_mean=float(devs.mean()),
        conservation=float(pre.max()),
        speed_drift_pre_max=float(pre.max()),
        speed_drift_pre_mean=float(pre.mean()),
        direction_drift_mean=float(dirs.mean()),
        position_error=float(errs.sum()),
        stable=stable, runtime=runtime,
        A_max=float(A_field.max()), A_mean=float(A_field.mean()),
    )


CANDIDATES_BASIC = [
    ("Cand 1", "Linear (control)", lambda g, mg: cand_linear(g, mg)),
    ("Cand 2", "Square root (normalised)", lambda g, mg: cand_sqrt(g, mg)),
    ("Cand 3", "Quadratic (normalised)", lambda g, mg: cand_quadratic(g, mg)),
    ("Cand 4", "Logarithmic (normalised)", lambda g, mg: cand_log(g, mg)),
    ("Cand 5", "Saturating |∇C|/(1+β|∇C|)", lambda g, mg: cand_saturating(g, mg)),
    ("Cand 7", "Piecewise linear (break at 50%)", lambda g, mg: cand_piecewise(g, mg)),
]


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--output", type=Path, default=DEFAULT_OUT)
    a = p.parse_args()
    a.output.mkdir(parents=True, exist_ok=True)
    fi = load_inputs()

    g = np.hypot(fi.gradient_x, fi.gradient_y)
    max_g = float(g.max())
    gx_hat = fi.gradient_x / np.maximum(g, 1e-15)
    gy_hat = fi.gradient_y / np.maximum(g, 1e-15)
    bad = g < 1e-15
    gx_hat = np.where(bad, 1.0, gx_hat)
    gy_hat = np.where(bad, 0.0, gy_hat)

    print(f"Control max|∇C| = {max_g:.4e}")
    print(f"Control mean|∇C| = {g.mean():.4e}")

    # ---- Validation ----------------------------------------------------
    print("\n=== Validation (every candidate must have max(A) = max_g) ===")
    A_fields = {}
    validation_rows = []
    all_pass = True

    for name, label, fn in CANDIDATES_BASIC:
        A = fn(g, max_g)
        A_fields[name] = A
        ok = abs(float(A.max()) - max_g) < 1e-12 * max(max_g, 1e-30)
        validation_rows.append((name, label, float(A.max()), float(A.mean()), ok))
        print(f"  {name}: max(A)={A.max():.4e}  mean(A)={A.mean():.4e}  "
              f"[{'PASS' if ok else 'FAIL'}]")
        if not ok:
            all_pass = False

    # Threshold candidates — natural max already equals max_g.
    for thr in [0.05, 0.10, 0.20]:
        name = f"Cand 6.{int(thr*100):02d}"
        label = f"Threshold at {int(thr*100)}% of max|∇C|"
        A = cand_threshold(g, max_g, thr)
        A_fields[name] = A
        ok = abs(float(A.max()) - max_g) < 1e-12 * max(max_g, 1e-30)
        validation_rows.append((name, label, float(A.max()), float(A.mean()), ok))
        print(f"  {name}: max(A)={A.max():.4e}  mean(A)={A.mean():.4e}  "
              f"[{'PASS' if ok else 'FAIL'}]")
        if not ok:
            all_pass = False

    with (a.output / "validation.csv").open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["Candidate", "Label", "Max A", "Mean A", "Pass/Fail"])
        for row in validation_rows:
            w.writerow([row[0], row[1], f"{row[2]:.4e}",
                        f"{row[3]:.4e}", "PASS" if row[4] else "FAIL"])

    if not all_pass:
        print("\nExperiment invalid: amplitude laws not normalised.")
        return 1

    # ---- Propagation --------------------------------------------------
    print("\n=== Propagation ===")
    measurements = []
    for name, label, fn in CANDIDATES_BASIC:
        A = A_fields[name]
        m = run_candidate(name, label, A, fi, gx_hat, gy_hat)
        measurements.append(m)
        print(f"  {name}: bend={m.bend_max:.4e}  cons={m.conservation:.4e}  "
              f"stable={m.stable}  runtime={m.runtime:.3f}s")
    for thr in [0.05, 0.10, 0.20]:
        name = f"Cand 6.{int(thr*100):02d}"
        label = f"Threshold at {int(thr*100)}% of max|∇C|"
        A = A_fields[name]
        m = run_candidate(name, label, A, fi, gx_hat, gy_hat)
        measurements.append(m)
        print(f"  {name}: bend={m.bend_max:.4e}  cons={m.conservation:.4e}  "
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
        w.writerow(["Candidate", "Bend", "Conservation", "Runtime (s)"])
        for m in measurements:
            w.writerow([m.candidate, f"{m.bend_max:.4e}",
                        f"{m.conservation:.4e}",
                        f"{m.runtime:.3f}"])

    # Relative comparison vs control
    control = measurements[0]
    rel = []
    for m in measurements:
        bend_pct = (m.bend_max - control.bend_max) / control.bend_max * 100
        cons_pct = (m.conservation - control.conservation) / control.conservation * 100 if control.conservation > 0 else 0.0
        runtime_pct = (m.runtime - control.runtime) / control.runtime * 100
        rel.append({"candidate": m.candidate,
                    "label": m.label,
                    "bend_pct_vs_control": bend_pct,
                    "cons_pct_vs_control": cons_pct,
                    "runtime_pct_vs_control": runtime_pct})
    with (a.output / "relative_vs_control.csv").open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["Candidate", "Label", "Bend Δ%",
                    "Conservation Δ%", "Runtime Δ%"])
        for r in rel:
            w.writerow([r["candidate"], r["label"],
                        f"{r['bend_pct_vs_control']:+.4f}%",
                        f"{r['cons_pct_vs_control']:+.4f}%",
                        f"{r['runtime_pct_vs_control']:+.4f}%"])

    # Plot A vs |∇C| for every candidate
    g_sorted = np.sort(g.ravel())
    fig, ax = plt.subplots(figsize=(8, 5))
    for name, _, fn in CANDIDATES_BASIC:
        A = A_fields[name]
        A_sorted = np.sort(A.ravel())
        # Subsample to make plot manageable
        idx = np.linspace(0, len(g_sorted) - 1, 500).astype(int)
        ax.plot(g_sorted[idx], A_sorted[idx], label=name, lw=1.2)
    for thr in [0.05, 0.10, 0.20]:
        name = f"Cand 6.{int(thr*100):02d}"
        A = A_fields[name]
        A_sorted = np.sort(A.ravel())
        idx = np.linspace(0, len(g_sorted) - 1, 500).astype(int)
        ax.plot(g_sorted[idx], A_sorted[idx], "--", label=name, lw=1.0)
    ax.set_xlabel("|∇C|")
    ax.set_ylabel("A")
    ax.set_title("Response amplitude A vs |∇C| (normalised to peak)")
    ax.legend(loc="best", fontsize=8)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(a.output / "response_curves.png", dpi=140)
    plt.close(fig)

    # Indistinguishable group (5% tolerance on bend_max)
    indistinguishable = []
    for m in measurements:
        if abs(m.bend_max - control.bend_max) / control.bend_max < 0.05:
            indistinguishable.append(m.candidate)
    (a.output / "indistinguishable.json").write_text(json.dumps(
        {"tolerance": "5% on bend_max vs control",
         "candidates_within_tolerance": indistinguishable}, indent=2))

    print(f"\n=== Indistinguishable from control (5% on bend_max) ===")
    for c in indistinguishable:
        print(f"  {c}")
    print(f"\nArtefacts written to {a.output}/")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())