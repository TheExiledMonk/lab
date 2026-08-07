#!/usr/bin/env python3
"""PBUF TRANSPORT-LAB-007 — Local neighbour update mechanism ablation.

Frozen: Lens-001, constitutive solution, neighbour graph, kernel envelope,
response magnitude, response angle (90°), integration, timestep, normalisation.

Variable: how the neighbour updates its state after receiving the response.

The response is the magnitude-normalised 90° transverse response from
LAB-006: r = |∇C| · R_90(∇̂C). Same response fed into all 7 update rules.
"""
from __future__ import annotations

import argparse
import csv
import json
import time
from dataclasses import dataclass, asdict
from pathlib import Path

import numpy as np

from transport_lab001 import FrozenInputs, load_inputs


ROOT = Path(__file__).resolve().parent
DEFAULT_OUT = ROOT / "runs" / "transport_lab007"


# ----------------------------------------------------------------------------
# Update rules. Each takes (v, r, step) and returns v_new.
# ----------------------------------------------------------------------------

def upd_direct_addition(v, r, step):
    """Exp 1: v_new = v + step * r, then normalise."""
    v_new = v + step * r
    n = max(np.linalg.norm(v_new), 1e-12)
    return v_new / n


def upd_pure_rotation(v, r, step):
    """Exp 2: v_new = R(step * |r|) * v."""
    angle = step * float(np.linalg.norm(r))
    c, s = np.cos(angle), np.sin(angle)
    return np.array([v[0] * c - v[1] * s, v[0] * s + v[1] * c])


def upd_projection(v, r, step):
    """Exp 3: v_new = v + step * (r - (r . v_hat) v_hat)."""
    n = max(np.linalg.norm(v), 1e-12)
    v_hat = v / n
    r_perp = r - float(np.dot(r, v_hat)) * v_hat
    v_new = v + step * r_perp
    n = max(np.linalg.norm(v_new), 1e-12)
    return v_new / n


def upd_shear(v, r, step):
    """Exp 4: incremental shear. Decompose v into v_para and v_perp
    relative to r; scale v_perp by (1 + step * |r|)."""
    rn = max(np.linalg.norm(r), 1e-12)
    r_hat = r / rn
    r_perp_hat = np.array([-r_hat[1], r_hat[0]])
    v_para = float(np.dot(v, r_hat)) * r_hat
    v_perp = v - v_para
    shear = step * rn
    v_new = v_para + v_perp * (1.0 + shear)
    n = max(np.linalg.norm(v_new), 1e-12)
    return v_new / n


def upd_momentum(v, r, step, mass: float = 0.5):
    """Exp 5: local momentum transfer. v_new = (m * v + step * r) / m."""
    p_new = mass * v + step * r
    v_new = p_new / mass
    n = max(np.linalg.norm(v_new), 1e-12)
    return v_new / n


def upd_phase(v, r, step):
    """Exp 6: phase transfer. v_new = R(step * arg(r)) * v."""
    rn = float(np.linalg.norm(r))
    if rn < 1e-15:
        return v.copy()
    r_angle = float(np.arctan2(r[1], r[0]))
    angle = step * r_angle
    c, s = np.cos(angle), np.sin(angle)
    return np.array([v[0] * c - v[1] * s, v[0] * s + v[1] * c])


def upd_mixed(v, r, step):
    """Exp 7: simultaneous rotation + translation + shear, equal weight."""
    # Rotation (1/3 weight)
    rot_angle = step * float(np.linalg.norm(r)) / 3.0
    c, s = np.cos(rot_angle), np.sin(rot_angle)
    v_rot = np.array([v[0] * c - v[1] * s, v[0] * s + v[1] * c])
    # Translation (1/3 weight)
    v_trans = v_rot + step * r / 3.0
    # Shear (1/3 weight): project r onto v_trans_perp
    n = max(np.linalg.norm(v_trans), 1e-12)
    v_hat = v_trans / n
    r_perp = r / 3.0 - float(np.dot(r / 3.0, v_hat)) * v_hat
    v_new = v_trans + step * r_perp
    n = max(np.linalg.norm(v_new), 1e-12)
    return v_new / n


UPDATES = [
    ("Exp 1", "Direct vector addition (baseline)", upd_direct_addition),
    ("Exp 2", "Pure rotation",                     upd_pure_rotation),
    ("Exp 3", "Projection update (r ⊥ v only)",    upd_projection),
    ("Exp 4", "Incremental shear",                  upd_shear),
    ("Exp 5", "Local momentum transfer (m=0.5)",    upd_momentum),
    ("Exp 6", "Phase transfer",                     upd_phase),
    ("Exp 7", "Mixed (rotation + translation + shear)", upd_mixed),
]


# ----------------------------------------------------------------------------
# Propagation kernel: pluggable update rule.
# ----------------------------------------------------------------------------

def propagate_with_update(fi: FrozenInputs, A_field, gx_hat, gy_hat,
                          update_fn, n_photons: int = 9, y_span: float = 3.0):
    n = fi.n
    x = np.linspace(-fi.extent, fi.extent, n)
    step = fi.photon_step_size
    n_steps = fi.photon_steps

    paths = []
    speed_drift_pre = []
    speed_drift_post = []
    direction_drift = []
    pos_error = []
    started = time.time()

    for y0 in np.linspace(-y_span, y_span, n_photons):
        xs = np.full(n_steps, -fi.extent)
        ys = np.full(n_steps, y0)
        vx = np.ones(n_steps)
        vy = np.zeros(n_steps)
        per_speed_pre = [0.0]
        per_speed_post = [0.0]
        per_dir = [0.0]
        per_err = [0.0]

        for k in range(1, n_steps):
            ix = int(np.clip(np.searchsorted(x, xs[k - 1]) - 1, 0, n - 1))
            iy = int(np.clip(np.searchsorted(x, ys[k - 1]) - 1, 0, n - 1))
            A = float(A_field[iy, ix])
            ux = float(gx_hat[iy, ix])
            uy = float(gy_hat[iy, ix])
            rx, ry = -A * uy, A * ux  # R_90 of (ux, uy), scaled by A
            v = np.array([vx[k - 1], vy[k - 1]])
            r = np.array([rx, ry])
            v_new = update_fn(v, r, step)
            v_new = np.asarray(v_new, dtype=float)
            n_v = max(np.linalg.norm(v_new), 1e-12)
            pre_norm = n_v
            per_speed_pre.append(abs(pre_norm - 1.0))
            v_new_unit = v_new / n_v
            post_norm = float(np.linalg.norm(v_new_unit))
            per_speed_post.append(abs(post_norm - 1.0))
            vx[k], vy[k] = v_new_unit[0], v_new_unit[1]
            xs[k] = xs[k - 1] + step * vx[k]
            ys[k] = ys[k - 1] + step * vy[k]

            v_old = np.array([vx[k - 1], vy[k - 1]])
            v_n = np.array([vx[k], vy[k]])
            cos_t = float(np.clip(
                v_old @ v_n / (np.linalg.norm(v_old) * np.linalg.norm(v_n) + 1e-15),
                -1.0, 1.0))
            per_dir.append(float(np.arccos(cos_t)))

            x_straight = -fi.extent + step * k
            y_straight = y0
            per_err.append(float(np.hypot(xs[k] - x_straight, ys[k] - y_straight)))

        paths.append((xs.copy(), ys.copy()))
        speed_drift_pre.append(np.array(per_speed_pre))
        speed_drift_post.append(np.array(per_speed_post))
        direction_drift.append(np.array(per_dir))
        pos_error.append(np.array(per_err))

    runtime = time.time() - started
    return paths, {
        "speed_drift_pre": speed_drift_pre,
        "speed_drift_post": speed_drift_post,
        "direction_drift": direction_drift,
        "pos_error": pos_error,
    }, runtime


# ----------------------------------------------------------------------------
# Measurement helpers
# ----------------------------------------------------------------------------

@dataclass
class M:
    rule: str
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
    behaviour: str


def behaviour_class(paths) -> str:
    """Classify behaviour as convergent/divergent, oscillatory/monotonic.

    Convergent: max position error is finite and below 1.0
    Divergent: any path diverges (|y| > extent)
    Oscillatory: photons change y-direction at least once
    Monotonic: photons move in one direction (no sign change in vy)
    """
    converged = True
    oscillatory = False
    for xs, ys in paths:
        if np.max(np.abs(ys)) > 16.0 or not np.isfinite(ys).all():
            converged = False
        # Detect sign changes in vy by differencing
        vy = np.diff(ys)
        if np.any(vy[:-1] * vy[1:] < 0):
            oscillatory = True
    s = "stable" if converged else "unstable"
    s += " / oscillatory" if oscillatory else " / monotonic"
    if not converged:
        s = "unstable / divergent"
    return s


def measure(label, update_fn, fi, A_field, gx_hat, gy_hat) -> M:
    paths, diag, runtime = propagate_with_update(fi, A_field, gx_hat, gy_hat, update_fn)
    devs = np.array([float(np.max(np.abs(y - y[0]))) for x, y in paths])
    pre = np.concatenate(diag["speed_drift_pre"])
    post = np.concatenate(diag["speed_drift_post"])
    dirs = np.concatenate(diag["direction_drift"])
    errs = np.concatenate(diag["pos_error"])
    finite = all(np.isfinite(np.concatenate([p[0], p[1]])).all() for p in paths)
    stable = bool(finite and pre.max() < 1.0)
    return M(
        rule=label[0], label=label[1],
        bend_max=float(devs.max()),
        bend_mean=float(devs.mean()),
        conservation=float(pre.max()),
        speed_drift_pre_max=float(pre.max()),
        speed_drift_pre_mean=float(pre.mean()),
        direction_drift_mean=float(dirs.mean()),
        position_error=float(errs.sum()),
        stable=stable,
        runtime=runtime,
        behaviour=behaviour_class(paths),
    )


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--output", type=Path, default=DEFAULT_OUT)
    a = p.parse_args()
    a.output.mkdir(parents=True, exist_ok=True)
    fi = load_inputs()

    # Use the magnitude-normalised 90° transverse response (LAB-006 control).
    gx_hat = fi.gradient_x / np.maximum(np.hypot(fi.gradient_x, fi.gradient_y), 1e-15)
    gy_hat = fi.gradient_y / np.maximum(np.hypot(fi.gradient_x, fi.gradient_y), 1e-15)
    A_field = np.hypot(fi.gradient_x, fi.gradient_y)
    # Where |g| is below the safety threshold, fall back to (1, 0).
    bad = np.hypot(fi.gradient_x, fi.gradient_y) < 1e-15
    gx_hat = np.where(bad, 1.0, gx_hat)
    gy_hat = np.where(bad, 0.0, gy_hat)

    measurements = []
    for rule_id, label, fn in UPDATES:
        m = measure((rule_id, label), fn, fi, A_field, gx_hat, gy_hat)
        measurements.append(m)
        print(f"  {rule_id}: bend={m.bend_max:.4e}  cons={m.conservation:.4e}  "
              f"stable={m.stable}  behaviour={m.behaviour}  runtime={m.runtime:.3f}s")

    rows = [asdict(m) for m in measurements]
    keys = list(rows[0].keys())
    with (a.output / "measurements.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        for r in rows:
            w.writerow(r)
    (a.output / "measurements.json").write_text(json.dumps(rows, indent=2))

    # Ranking table
    score = []
    for m in measurements:
        bends = np.array([x.bend_max for x in measurements])
        conservs = np.array([x.conservation for x in measurements])
        b_max = bends.max()
        c_min = conservs.min()
        c_eff = m.conservation if m.conservation > 0 else 1e-30
        s = (m.bend_max / b_max) + (c_min / c_eff)
        score.append(s)
    order = np.argsort(score)[::-1]
    with (a.output / "ranking.csv").open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["Rank", "Rule", "Update Mechanism", "Bend",
                    "Conservation", "Stable", "Runtime (s)", "Score"])
        for rank_idx, idx in enumerate(order):
            m = measurements[idx]
            w.writerow([rank_idx + 1, m.rule, m.label,
                        f"{m.bend_max:.4e}", f"{m.conservation:.4e}",
                        "yes" if m.stable else "no",
                        f"{m.runtime:.3f}", f"{score[idx]:.3f}"])

    # Behaviour summary
    with (a.output / "behaviour_summary.csv").open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["Rule", "Mechanism", "Stability", "Behaviour"])
        for m in measurements:
            parts = m.behaviour.split(" / ")
            stability = parts[0]
            behaviour = parts[1] if len(parts) > 1 else "—"
            w.writerow([m.rule, m.label, stability, behaviour])

    print(f"\n=== Ranking ===")
    for rank_idx, idx in enumerate(order):
        m = measurements[idx]
        print(f"  #{rank_idx + 1}: {m.rule} ({m.label})  "
              f"bend={m.bend_max:.3e}  cons={m.conservation:.3e}  score={score[idx]:.3f}")
    print(f"\nArtefacts written to {a.output}/")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
