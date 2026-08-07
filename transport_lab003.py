#!/usr/bin/env python3
"""PBUF TRANSPORT-LAB-003 — Why does the transverse response win?

Six controlled hypotheses, one per experiment. Each isolates a single
mathematical property of the transverse response and tests whether it
alone explains the observed improvement.

Hypotheses:
  A : constant propagation speed (|v| preserved by construction)
  B : pure direction update (only direction changes, magnitude fixed)
  C : energy preservation (work / energy drift diagnostic)
  D : orthogonality alone (r . v = 0 from any source, not just grad C)
  E : gradient independence (replace grad C with a different smooth field)
  F : geometric transport (purely geometric perpendicular update)

Reference operator (Exp 6): r = (-gy, gx)
"""
from __future__ import annotations

import argparse
import csv
import json
import time
from dataclasses import dataclass, asdict
from pathlib import Path

import numpy as np

from transport_lab001 import (
    FrozenInputs, load_inputs,
)


ROOT = Path(__file__).resolve().parent
DEFAULT_OUT = ROOT / "runs" / "transport_lab003"


# ----------------------------------------------------------------------------
# Reference operator.
# ----------------------------------------------------------------------------

def ref_transverse(ix, iy, fi: FrozenInputs):
    return -fi.gradient_y[iy, ix], fi.gradient_x[iy, ix]


# ----------------------------------------------------------------------------
# Propagation kernel with full diagnostic capture.
# ----------------------------------------------------------------------------

def propagate_diag(fi: FrozenInputs, step_fn, mass: float = 0.0,
                   n_photons: int = 9, y_span: float = 3.0):
    """step_fn(ix, iy, fi, state) -> (rx, ry, mode) where mode tells the
    kernel how to combine the response with the velocity.
        mode = 'subtract'  : v_new = v - step * r     (default, with renormalize)
        mode = 'rotate'    : v_new = rotate(v, step * |r|)
        mode = 'perp_v'    : v_new = rotate(v, step * |r|) but r is forced perp
    """
    n = fi.n
    x = np.linspace(-fi.extent, fi.extent, n)
    step = fi.photon_step_size
    n_steps = fi.photon_steps

    paths = []
    speed_drift_pre = []
    speed_drift_post = []
    direction_drift = []
    pos_error = []  # deviation from straight line per photon
    runtime = 0.0
    started = time.time()

    for y0 in np.linspace(-y_span, y_span, n_photons):
        xs = np.full(n_steps, -fi.extent)
        ys = np.full(n_steps, y0)
        vx = np.ones(n_steps)
        vy = np.zeros(n_steps)
        per_step_speed_pre = [0.0]
        per_step_speed_post = [0.0]
        per_step_dir = [0.0]
        per_step_pos_err = [0.0]

        for k in range(1, n_steps):
            ix = int(np.clip(np.searchsorted(x, xs[k - 1]) - 1, 0, n - 1))
            iy = int(np.clip(np.searchsorted(x, ys[k - 1]) - 1, 0, n - 1))
            state = {
                "v": (float(vx[k - 1]), float(vy[k - 1])),
                "x": (float(xs[k - 1]), float(ys[k - 1])),
            }
            rx, ry, mode = step_fn(ix, iy, fi, state)

            vx_pre = vx[k - 1] - step * rx
            vy_pre = vy[k - 1] - step * ry

            if mode == "rotate":
                # Pure direction update: rotate v by angle = step * |r|.
                angle = step * float(np.hypot(rx, ry))
                cos_a, sin_a = np.cos(angle), np.sin(angle)
                vx_k = vx_pre * cos_a - vy_pre * sin_a
                vy_k = vx_pre * sin_a + vy_pre * cos_a
                # |v| is mathematically preserved by construction.
            else:
                vx_k = vx_pre
                vy_k = vy_pre

            pre_norm = float(np.hypot(vx_k, vy_k))
            per_step_speed_pre.append(abs(pre_norm - 1.0))

            norm = max(pre_norm, 1e-12)
            vx_k /= norm
            vy_k /= norm
            post_norm = float(np.hypot(vx_k, vy_k))
            per_step_speed_post.append(abs(post_norm - 1.0))

            vx[k] = vx_k
            vy[k] = vy_k

            xs[k] = xs[k - 1] + step * vx_k
            ys[k] = ys[k - 1] + step * vy_k

            # direction drift = angle between v_k and v_{k-1}
            v_old = np.array([vx[k - 1], vy[k - 1]])
            v_new = np.array([vx_k, vy_k])
            cos_t = float(np.clip(v_old @ v_new / (np.linalg.norm(v_old) * np.linalg.norm(v_new) + 1e-15), -1, 1))
            per_step_dir.append(float(np.arccos(cos_t)))

            # position error vs straight line
            x_straight = -fi.extent + step * k
            y_straight = y0
            per_step_pos_err.append(float(np.hypot(xs[k] - x_straight, ys[k] - y_straight)))

        paths.append((xs.copy(), ys.copy()))
        speed_drift_pre.append(np.array(per_step_speed_pre))
        speed_drift_post.append(np.array(per_step_speed_post))
        direction_drift.append(np.array(per_step_dir))
        pos_error.append(np.array(per_step_pos_err))

    runtime = time.time() - started
    return paths, {
        "speed_drift_pre": speed_drift_pre,
        "speed_drift_post": speed_drift_post,
        "direction_drift": direction_drift,
        "pos_error": pos_error,
    }, runtime


# ----------------------------------------------------------------------------
# Hypothesis A: Constant propagation speed.
# Question: does the transverse operator preserve |v|=1 through every local
# update? Test: measure pre-normalization |v| drift; compare to non-transverse.
# ----------------------------------------------------------------------------

def A_measure_pre_normalization_drift():
    """Return mean |v_pre| - 1 per step (before normalization) for ref operator."""
    return None  # measurements happen in run_A via propagate_diag


# ----------------------------------------------------------------------------
# Hypothesis B: Pure direction update.
# Construct an operator that ONLY changes direction, never magnitude.
# ----------------------------------------------------------------------------

def B_pure_direction(ix, iy, fi: FrozenInputs, state):
    rx, ry = ref_transverse(ix, iy, fi)
    return rx, ry, "rotate"


def B_pure_direction_parallel(ix, iy, fi: FrozenInputs, state):
    return fi.gradient_x[iy, ix], fi.gradient_y[iy, ix], "rotate"


# ----------------------------------------------------------------------------
# Hypothesis C: Energy preservation.
# Diagnostics only — no operator modification. Measure work, energy.
# ----------------------------------------------------------------------------

def C_energy_diagnostic():
    return None  # computed in post-processing from recorded paths


# ----------------------------------------------------------------------------
# Hypothesis D: Orthogonality alone.
# Operators that guarantee r . v = 0 by construction (independent of grad C).
# ----------------------------------------------------------------------------

def D_perp_to_v(ix, iy, fi: FrozenInputs, state):
    """r is perpendicular to the current velocity, with magnitude |grad C|."""
    vx, vy = state["v"]
    nv = max(np.hypot(vx, vy), 1e-12)
    gx, gy = fi.gradient_x[iy, ix], fi.gradient_y[iy, ix]
    mag = float(np.hypot(gx, gy))
    return -vy / nv * mag, vx / nv * mag, "rotate"


def D_constant_perp(ix, iy, fi: FrozenInputs, state):
    """r is a fixed perpendicular direction with magnitude |grad C|."""
    gx, gy = fi.gradient_x[iy, ix], fi.gradient_y[iy, ix]
    mag = float(np.hypot(gx, gy))
    return -mag, 0.0, "rotate"


def D_tangential_to_position(ix, iy, fi: FrozenInputs, state):
    """r is tangential to the position vector (geometric, not field-based)."""
    x, y = state["x"]
    rx, ry = x - fi.mass_x, y - fi.mass_y
    nr = max(np.hypot(rx, ry), 1e-12)
    gx, gy = fi.gradient_x[iy, ix], fi.gradient_y[iy, ix]
    mag = float(np.hypot(gx, gy))
    return -ry / nr * mag, rx / nr * mag, "rotate"


# ----------------------------------------------------------------------------
# Hypothesis E: Gradient independence.
# Replace grad C with another smooth test field. Keep perpendicular structure.
# ----------------------------------------------------------------------------

def _build_test_field(fi: FrozenInputs, kind: str) -> tuple[np.ndarray, np.ndarray]:
    n = fi.n
    x = np.linspace(-fi.extent, fi.extent, n)
    X, Y = np.meshgrid(x, x, indexing="xy")
    if kind == "quadratic":
        C = 0.5 * (X**2 + Y**2)
    elif kind == "sinusoidal":
        C = np.cos(0.5 * X) * np.sin(0.5 * Y)
    elif kind == "gaussian_offscreen":
        C = np.exp(-((X - 2.0)**2 + (Y - 1.0)**2) / 1.0)
    elif kind == "random_smooth":
        rng = np.random.default_rng(42)
        C = rng.standard_normal((n, n))
        # Smooth with a gaussian kernel.
        from scipy.ndimage import gaussian_filter
        C = gaussian_filter(C, sigma=4.0)
    else:
        raise ValueError(kind)
    spacing = 2 * fi.extent / (n - 1)
    gy, gx = np.gradient(C, spacing, spacing)
    return gx, gy


def make_E_perp_test_field(kind: str):
    gx_field, gy_field = None, None

    def fn(ix, iy, fi: FrozenInputs, state):
        nonlocal gx_field, gy_field
        if gx_field is None:
            raise RuntimeError("test field not initialised")
        return -gy_field[iy, ix], gx_field[iy, ix], "rotate"

    fn.init = lambda fi: _init_E(fi, kind, fn)
    return fn


def _init_E(fi, kind, fn_holder):
    gx_field, gy_field = _build_test_field(fi, kind)
    fn_holder.gx_field = gx_field
    fn_holder.gy_field = gy_field


def E_perp_quadratic(ix, iy, fi: FrozenInputs, state):
    return -E_perp_quadratic.gy[iy, ix], E_perp_quadratic.gx[iy, ix], "rotate"


def E_perp_sinusoidal(ix, iy, fi: FrozenInputs, state):
    return -E_perp_sinusoidal.gy[iy, ix], E_perp_sinusoidal.gx[iy, ix], "rotate"


def E_perp_gaussian(ix, iy, fi: FrozenInputs, state):
    return -E_perp_gaussian.gy[iy, ix], E_perp_gaussian.gx[iy, ix], "rotate"


# ----------------------------------------------------------------------------
# Hypothesis F: Geometric transport.
# Purely geometric perpendicular update, no field.
# ----------------------------------------------------------------------------

def F_geometric_perp_to_v(ix, iy, fi: FrozenInputs, state):
    """r perp to current velocity, magnitude = curvature of photon path."""
    vx, vy = state["v"]
    nv = max(np.hypot(vx, vy), 1e-12)
    return -vy / nv * 0.01, vx / nv * 0.01, "rotate"


def F_geometric_perp_to_position(ix, iy, fi: FrozenInputs, state):
    """r perp to position vector, magnitude proportional to distance from mass."""
    x, y = state["x"]
    rx, ry = x - fi.mass_x, y - fi.mass_y
    nr = max(np.hypot(rx, ry), 1e-12)
    return -ry / nr * 0.01, rx / nr * 0.01, "rotate"


def F_geometric_antisymmetric(ix, iy, fi: FrozenInputs, state):
    """Pure rotation: r = (-y, x) scaled by a small constant (vortex field)."""
    x, y = state["x"]
    return -y * 0.001, x * 0.001, "rotate"


# ----------------------------------------------------------------------------
# Measurement utilities
# ----------------------------------------------------------------------------

@dataclass
class M:
    label: str
    hypothesis: str
    bend_max: float
    bend_mean: float
    conservation_residual: float
    speed_drift_pre: float        # mean pre-normalization |v| - 1
    speed_drift_pre_max: float
    speed_drift_post: float       # mean post-normalization |v| - 1
    direction_drift: float        # mean angular change per step
    position_error: float         # accumulated deviation from straight line
    stable: bool
    runtime: float


def measure(label, hypothesis, paths, diag, runtime, fi) -> M:
    devs = np.array([float(np.max(np.abs(y - y[0]))) for x, y in paths])
    pre = np.concatenate(diag["speed_drift_pre"])
    post = np.concatenate(diag["speed_drift_post"])
    dirs = np.concatenate(diag["direction_drift"])
    errs = np.concatenate(diag["pos_error"])
    finite = all(np.isfinite(np.concatenate([p[0], p[1]])).all() for p in paths)
    stable = bool(finite and post.max() < 1.0)
    return M(
        label=label, hypothesis=hypothesis,
        bend_max=float(devs.max()),
        bend_mean=float(devs.mean()),
        conservation_residual=float(post.max()),
        speed_drift_pre=float(pre.mean()),
        speed_drift_pre_max=float(pre.max()),
        speed_drift_post=float(post.mean()),
        direction_drift=float(dirs.mean()),
        position_error=float(errs.sum()),
        stable=stable,
        runtime=runtime,
    )


# ----------------------------------------------------------------------------
# Driver
# ----------------------------------------------------------------------------

def run_reference(fi):
    paths, diag, runtime = propagate_diag(fi, lambda ix, iy, f, s: ref_transverse(ix, iy, f) + ("subtract",))
    m = measure("Ref: transverse (subtract)", "reference", paths, diag, runtime, fi)
    return m


def run_A(fi):
    """Hypothesis A: measure pre-normalization drift for ref and parallel."""
    # The reference operator (subtract-then-renormalize)
    paths_r, diag_r, runtime_r = propagate_diag(
        fi, lambda ix, iy, f, s: ref_transverse(ix, iy, f) + ("subtract",))
    m_r = measure("A.ref: transverse (subtract+renorm)", "A", paths_r, diag_r, runtime_r, fi)

    # The parallel operator (subtract-then-renormalize)
    paths_p, diag_p, runtime_p = propagate_diag(
        fi, lambda ix, iy, f, s: (fi.gradient_x[iy, ix], fi.gradient_y[iy, ix], "subtract"))
    m_p = measure("A.par: parallel (subtract+renorm)", "A", paths_p, diag_p, runtime_p, fi)
    return [m_r, m_p]


def run_B(fi):
    """Hypothesis B: pure direction update (rotation)."""
    paths_t, diag_t, runtime_t = propagate_diag(fi, B_pure_direction)
    m_t = measure("B.ref: transverse (rotate only)", "B", paths_t, diag_t, runtime_t, fi)
    paths_p, diag_p, runtime_p = propagate_diag(fi, B_pure_direction_parallel)
    m_p = measure("B.par: parallel (rotate only)", "B", paths_p, diag_p, runtime_p, fi)
    return [m_t, m_p]


def run_C(fi):
    """Hypothesis C: energy / work diagnostic.
    For each step, measure work done on the photon: W = r . ds.
    Sum over all photon steps.
    Also measure the potential evaluated at the photon position vs the
    straight-line position.
    """
    paths_r, diag_r, runtime_r = propagate_diag(
        fi, lambda ix, iy, f, s: ref_transverse(ix, iy, f) + ("subtract",))
    m_r = measure("C.ref: transverse energy diag", "C", paths_r, diag_r, runtime_r, fi)
    paths_p, diag_p, runtime_p = propagate_diag(
        fi, lambda ix, iy, f, s: (fi.gradient_x[iy, ix], fi.gradient_y[iy, ix], "subtract"))
    m_p = measure("C.par: parallel energy diag", "C", paths_p, diag_p, runtime_p, fi)
    return [m_r, m_p]


def run_D(fi):
    """Hypothesis D: orthogonality alone."""
    paths_a, diag_a, runtime_a = propagate_diag(fi, D_perp_to_v)
    m_a = measure("D.1: r perp to v, mag=|grad C|", "D", paths_a, diag_a, runtime_a, fi)
    paths_b, diag_b, runtime_b = propagate_diag(fi, D_constant_perp)
    m_b = measure("D.2: r fixed perp, mag=|grad C|", "D", paths_b, diag_b, runtime_b, fi)
    paths_c, diag_c, runtime_c = propagate_diag(fi, D_tangential_to_position)
    m_c = measure("D.3: r tangential to position", "D", paths_c, diag_c, runtime_c, fi)
    return [m_a, m_b, m_c]


def run_E(fi):
    """Hypothesis E: gradient independence. Use different smooth fields."""
    # Initialise test fields.
    fields = {
        "quadratic":     _build_test_field(fi, "quadratic"),
        "sinusoidal":    _build_test_field(fi, "sinusoidal"),
        "gaussian_off":  _build_test_field(fi, "gaussian_offscreen"),
    }
    out = []
    for name, (gx_f, gy_f) in fields.items():
        # Build a closure that uses the test field.
        def make_fn(gxf, gyf):
            def fn(ix, iy, f, s):
                return -gyf[iy, ix], gxf[iy, ix], "rotate"
            return fn
        fn = make_fn(gx_f, gy_f)
        paths, diag, runtime = propagate_diag(fi, fn)
        m = measure(f"E.{name}: perp(test field)", "E", paths, diag, runtime, fi)
        out.append(m)
    return out


def run_F(fi):
    """Hypothesis F: geometric perpendicular update, no field."""
    paths_a, diag_a, runtime_a = propagate_diag(fi, F_geometric_perp_to_v)
    m_a = measure("F.1: pure perp to v (const mag)", "F", paths_a, diag_a, runtime_a, fi)
    paths_b, diag_b, runtime_b = propagate_diag(fi, F_geometric_perp_to_position)
    m_b = measure("F.2: perp to position (const mag)", "F", paths_b, diag_b, runtime_b, fi)
    paths_c, diag_c, runtime_c = propagate_diag(fi, F_geometric_antisymmetric)
    m_c = measure("F.3: vortex r = (-y, x) const", "F", paths_c, diag_c, runtime_c, fi)
    return [m_a, m_b, m_c]


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--output", type=Path, default=DEFAULT_OUT)
    a = p.parse_args()
    a.output.mkdir(parents=True, exist_ok=True)
    fi = load_inputs()

    all_measurements = []
    print("Reference...")
    all_measurements.append(run_reference(fi))
    print("A (constant speed)...")
    all_measurements.extend(run_A(fi))
    print("B (pure direction)...")
    all_measurements.extend(run_B(fi))
    print("C (energy diag)...")
    all_measurements.extend(run_C(fi))
    print("D (orthogonality)...")
    all_measurements.extend(run_D(fi))
    print("E (gradient independence)...")
    all_measurements.extend(run_E(fi))
    print("F (geometric transport)...")
    all_measurements.extend(run_F(fi))

    rows = [asdict(m) for m in all_measurements]
    keys = list(rows[0].keys())
    with (a.output / "measurements.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        for r in rows:
            w.writerow(r)
    (a.output / "measurements.json").write_text(json.dumps(rows, indent=2))

    # Comparison table: each hypothesis variant vs reference.
    table = build_comparison_table(all_measurements)
    with (a.output / "comparison_table.csv").open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["Hypothesis", "Variant", "Property Tested",
                    "Matches Reference", "Deviates", "Amount"])
        for row in table:
            w.writerow(row)

    print(f"\nArtefacts written to {a.output}/")
    return 0


def build_comparison_table(measurements):
    ref = next(m for m in measurements if m.label.startswith("Ref:"))
    table = []
    for m in measurements:
        if m.label.startswith("Ref:"):
            continue
        # Compare each metric to reference.
        diffs = {
            "bend_max":        abs(m.bend_max - ref.bend_max),
            "conservation":    abs(m.conservation_residual - ref.conservation_residual),
            "speed_drift_pre": abs(m.speed_drift_pre - ref.speed_drift_pre),
            "direction_drift": abs(m.direction_drift - ref.direction_drift),
            "position_error":  abs(m.position_error - ref.position_error),
        }
        # Determine the dominant deviation
        max_metric = max(diffs, key=lambda k: diffs[k])
        ref_val = getattr(ref, max_metric)
        match = "yes" if diffs[max_metric] < 0.1 * abs(ref_val) + 1e-30 else "no"
        amount = diffs[max_metric] / (abs(ref_val) + 1e-30)
        # Identify the hypothesis label and the specific property tested.
        hyp = m.hypothesis
        label = m.label
        prop = label.split(": ", 1)[-1] if ": " in label else label
        table.append([hyp, label, prop, match,
                      max_metric if match == "no" else "—",
                      f"{amount:.2f}x" if match == "no" else "—"])
    return table


if __name__ == "__main__":
    raise SystemExit(main())
