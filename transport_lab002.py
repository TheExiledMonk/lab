#!/usr/bin/env python3
"""PBUF TRANSPORT-LAB-002 — Ablation of Experiment 6.

Decomposes the magnetic-style directional operator from TRANSPORT-LAB-001
into individual mathematical ingredients. Each experiment modifies exactly
one property; everything else (Lens-001 inputs, kernel, integration) is
identical.

Reference operator (Exp 6 from TRANSPORT-LAB-001):
    response(ix, iy) = (-gy, gx)        # perpendicular to grad C

Controlled experiments A..H each vary one ingredient.

Output: measurements + contribution table. No interpretation.
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
    FrozenInputs, load_inputs, propagate,
)


ROOT = Path(__file__).resolve().parent
DEFAULT_OUT = ROOT / "runs" / "transport_lab002"


# ----------------------------------------------------------------------------
# Reference operator — exact copy of TRANSPORT-LAB-001 Exp 6.
# ----------------------------------------------------------------------------

def ref_magnetic(ix, iy, fi: FrozenInputs):
    return -fi.gradient_y[iy, ix], fi.gradient_x[iy, ix]


# ----------------------------------------------------------------------------
# Experiment A — remove transverse projection.
# Use the parallel (gradient) direction instead of the perpendicular.
# ----------------------------------------------------------------------------

def A_no_transverse(ix, iy, fi: FrozenInputs):
    return fi.gradient_x[iy, ix], fi.gradient_y[iy, ix]


# ----------------------------------------------------------------------------
# Experiment B — keep only the transverse projector.
# Project the reference response onto the direction perpendicular to grad C.
# ----------------------------------------------------------------------------

def B_transverse_only(ix, iy, fi: FrozenInputs):
    rx, ry = ref_magnetic(ix, iy, fi)
    gx, gy = fi.gradient_x[iy, ix], fi.gradient_y[iy, ix]
    gnorm2 = gx * gx + gy * gy
    if gnorm2 < 1e-24:
        return rx, ry
    proj = (rx * gx + ry * gy) / gnorm2
    # Transverse projector: subtract the parallel part from the response.
    return rx - proj * gx, ry - proj * gy


# ----------------------------------------------------------------------------
# Experiment C — reverse steering direction (sign flip).
# ----------------------------------------------------------------------------

def C_reversed(ix, iy, fi: FrozenInputs):
    return fi.gradient_y[iy, ix], -fi.gradient_x[iy, ix]


# ----------------------------------------------------------------------------
# Experiment D — magnitude only, no directional dependence.
# ----------------------------------------------------------------------------

def D_magnitude_only(ix, iy, fi: FrozenInputs):
    gx, gy = fi.gradient_x[iy, ix], fi.gradient_y[iy, ix]
    mag = float(np.hypot(gx, gy))
    return mag, 0.0


# ----------------------------------------------------------------------------
# Experiment E — normalization study. Same operator, scaled magnitude.
# ----------------------------------------------------------------------------

def make_E(scale: float):
    def fn(ix, iy, fi: FrozenInputs):
        return -scale * fi.gradient_y[iy, ix], scale * fi.gradient_x[iy, ix]
    return fn


# ----------------------------------------------------------------------------
# Experiment F — locality radius. Average over a wider stencil.
# ----------------------------------------------------------------------------

def _avg_field(field: np.ndarray, ix: int, iy: int, radius: int) -> float:
    n = field.shape[0]
    ix0, ix1 = max(ix - radius, 0), min(ix + radius, n - 1)
    iy0, iy1 = max(iy - radius, 0), min(iy + radius, n - 1)
    return float(field[iy0:iy1 + 1, ix0:ix1 + 1].mean())


def F_locality(radius: int):
    def fn(ix, iy, fi: FrozenInputs):
        gx = _avg_field(fi.gradient_x, ix, iy, radius)
        gy = _avg_field(fi.gradient_y, ix, iy, radius)
        return -gy, gx
    return fn


# ----------------------------------------------------------------------------
# Experiment G — decomposition into transverse / gradient / rotational pieces.
#
# For the reference operator r = (-∂y C, +∂x C):
#   * transverse piece (perp to grad C): r itself  (the full operator)
#   * gradient piece  (parallel to grad C): 0     (r is exactly perpendicular)
#   * rotational piece (curl of r): curl r = ∂x(∂x C) + ∂y(∂y C) = ∇²C
# ----------------------------------------------------------------------------

def _laplacian(fi: FrozenInputs, ix: int, iy: int) -> float:
    n = fi.n
    ixp, ixm = min(ix + 1, n - 1), max(ix - 1, 0)
    iyp, iym = min(iy + 1, n - 1), max(iy - 1, 0)
    spacing = 2 * fi.extent / (n - 1)
    return float(
        (fi.C[iy, ixp] - 2 * fi.C[iy, ix] + fi.C[iy, ixm]) / spacing**2
        + (fi.C[iyp, ix] - 2 * fi.C[iy, ix] + fi.C[iym, ix]) / spacing**2
    )


def G_transverse(ix, iy, fi: FrozenInputs):
    return -fi.gradient_y[iy, ix], fi.gradient_x[iy, ix]


def G_gradient(ix, iy, fi: FrozenInputs):
    return fi.gradient_x[iy, ix], fi.gradient_y[iy, ix]


def G_rotational(ix, iy, fi: FrozenInputs):
    lap = _laplacian(fi, ix, iy)
    return lap, 0.0


def G_T_plus_Grad(ix, iy, fi: FrozenInputs):
    a, b = G_transverse(ix, iy, fi)
    c, d = G_gradient(ix, iy, fi)
    return a + c, b + d


def G_T_plus_Rot(ix, iy, fi: FrozenInputs):
    a, b = G_transverse(ix, iy, fi)
    c, d = G_rotational(ix, iy, fi)
    return a + c, b + d


def G_Grad_plus_Rot(ix, iy, fi: FrozenInputs):
    a, b = G_gradient(ix, iy, fi)
    c, d = G_rotational(ix, iy, fi)
    return a + c, b + d


def G_all_three(ix, iy, fi: FrozenInputs):
    a, b = G_transverse(ix, iy, fi)
    c, d = G_gradient(ix, iy, fi)
    e, f = G_rotational(ix, iy, fi)
    return a + c + e, b + d + f


# ----------------------------------------------------------------------------
# Experiment H — conservation analysis.
# Diagnoses why the reference operator had the smallest |v|-drift.
# H.1 : same operator but force the response onto the velocity direction
# H.2 : same operator but force the response perpendicular to velocity
# H.3 : no normalization at all
# H.4 : normalize only every K steps
# H.5 : report <v . r> per step (the "perpendicularity diagnostic")
# ----------------------------------------------------------------------------

def H_force_parallel(ix, iy, fi: FrozenInputs, state):
    rx, ry = ref_magnetic(ix, iy, fi)
    vx, vy = state["v"]
    n = max(np.hypot(vx, vy), 1e-12)
    # Project response onto velocity direction.
    return rx * (vx / n), ry * (vy / n)


def H_force_perpendicular(ix, iy, fi: FrozenInputs, state):
    rx, ry = ref_magnetic(ix, iy, fi)
    vx, vy = state["v"]
    n = max(np.hypot(vx, vy), 1e-12)
    # Project response onto direction perpendicular to velocity.
    ux, uy = vx / n, vy / n
    perp = (-uy, ux)
    proj = rx * perp[0] + ry * perp[1]
    return proj * perp[0], proj * perp[1]


def propagate_with_state(fi: FrozenInputs, response_fn, mass: float = 0.0,
                         n_photons: int = 9, y_span: float = 3.0,
                         normalize_every: int = 1,
                         record_perp: bool = False):
    n = fi.n
    x = np.linspace(-fi.extent, fi.extent, n)
    step = fi.photon_step_size
    n_steps = fi.photon_steps

    paths = []
    all_velocity_norms = []
    perp_diagnostics = []
    started = time.time()

    for y0 in np.linspace(-y_span, y_span, n_photons):
        xs = np.full(n_steps, -fi.extent)
        ys = np.full(n_steps, y0)
        vx = np.ones(n_steps)
        vy = np.zeros(n_steps)
        velocity_norms = [1.0]
        per_photon_perp = []

        for k in range(1, n_steps):
            ix = int(np.clip(np.searchsorted(x, xs[k - 1]) - 1, 0, n - 1))
            iy = int(np.clip(np.searchsorted(x, ys[k - 1]) - 1, 0, n - 1))
            state = {"v": (float(vx[k - 1]), float(vy[k - 1]))}
            rx, ry = response_fn(ix, iy, fi, state)

            vx_k = vx[k - 1] - step * rx
            vy_k = vy[k - 1] - step * ry

            if normalize_every > 0 and k % normalize_every == 0:
                norm = max(np.hypot(vx_k, vy_k), 1e-12)
                vx_k /= norm
                vy_k /= norm
            vx[k] = vx_k
            vy[k] = vy_k

            xs[k] = xs[k - 1] + step * vx_k
            ys[k] = ys[k - 1] + step * vy_k
            velocity_norms.append(float(np.hypot(vx_k, vy_k)))
            if record_perp:
                per_photon_perp.append((float(rx), float(ry),
                                        float(vx[k - 1]), float(vy[k - 1])))

        paths.append((xs.copy(), ys.copy()))
        all_velocity_norms.append(np.array(velocity_norms))
        perp_diagnostics.append(per_photon_perp)

    runtime = time.time() - started
    return paths, all_velocity_norms, runtime, perp_diagnostics


def perp_diagnostic(perp_diagnostics):
    """Mean |v . r| / |v|/|r| across all photon-step events."""
    dots, norms_v, norms_r = [], [], []
    for photon in perp_diagnostics:
        for rx, ry, vx, vy in photon:
            nv = np.hypot(vx, vy)
            nr = np.hypot(rx, ry)
            if nv < 1e-12 or nr < 1e-12:
                continue
            dots.append(abs(rx * vx + ry * vy) / (nv * nr))
            norms_v.append(nv)
            norms_r.append(nr)
    if not dots:
        return 0.0
    return float(np.mean(dots))


# ----------------------------------------------------------------------------
# Measurement helpers
# ----------------------------------------------------------------------------

@dataclass
class M:
    label: str
    category: str
    stable: bool
    finite: bool
    bend_max: float
    bend_mean: float
    conservation_residual: float
    mean_velocity_drift: float
    runtime: float
    perp_metric: float = 0.0
    notes: str = ""


def measure_basic(label, category, paths, vnorms, runtime):
    finite = all(
        np.isfinite(np.concatenate([p[0], p[1]])).all() for p in paths
    )
    flat_norms = np.concatenate(vnorms)
    drift = np.abs(flat_norms - 1.0)
    devs = np.array([float(np.max(np.abs(y - y[0]))) for x, y in paths])
    return M(
        label=label, category=category,
        stable=bool(finite and drift.max() < 1.0),
        finite=finite,
        bend_max=float(devs.max()),
        bend_mean=float(devs.mean()),
        conservation_residual=float(drift.max()),
        mean_velocity_drift=float(drift.mean()),
        runtime=runtime,
    )


# ----------------------------------------------------------------------------
# Driver
# ----------------------------------------------------------------------------

def run_reference(fi: FrozenInputs):
    paths, vnorms, runtime = propagate(fi, ref_magnetic, mass=0.0)
    m = measure_basic("Ref Exp 6 (control)", "reference",
                      paths, vnorms, runtime)
    m.notes = "exact copy of TRANSPORT-LAB-001 Exp 6"
    return m


def run_A(fi: FrozenInputs):
    paths, vnorms, runtime = propagate(fi, A_no_transverse, mass=0.0)
    m = measure_basic("A: parallel (no transverse)", "A_no_transverse",
                      paths, vnorms, runtime)
    m.notes = "response = (gx, gy); gradient direction only"
    return m


def run_B(fi: FrozenInputs):
    paths, vnorms, runtime = propagate(fi, B_transverse_only, mass=0.0)
    m = measure_basic("B: transverse projector only", "B_transverse_only",
                      paths, vnorms, runtime)
    m.notes = "re-projects r onto direction perp to grad C"
    return m


def run_C(fi: FrozenInputs):
    paths, vnorms, runtime = propagate(fi, C_reversed, mass=0.0)
    m = measure_basic("C: reversed steering", "C_reversed",
                      paths, vnorms, runtime)
    m.notes = "response = (gy, -gx); opposite perpendicular"
    return m


def run_D(fi: FrozenInputs):
    paths, vnorms, runtime = propagate(fi, D_magnitude_only, mass=0.0)
    m = measure_basic("D: magnitude only (no orientation)", "D_magnitude_only",
                      paths, vnorms, runtime)
    m.notes = "response = (|grad C|, 0); loses orientation"
    return m


def run_E(fi: FrozenInputs):
    out = []
    for scale in [1.0, 0.75, 0.5, 0.25, 0.1]:
        paths, vnorms, runtime = propagate(fi, make_E(scale), mass=0.0)
        m = measure_basic(f"E: scale {scale:.2f}", "E_scale", paths, vnorms, runtime)
        m.notes = f"scale={scale}"
        out.append(m)
    return out


def run_F(fi: FrozenInputs):
    out = []
    for radius, label in [(0, "1-cell"), (1, "von Neumann (5 cells)"),
                          (2, "Moore (3x3)")]:
        paths, vnorms, runtime = propagate(fi, F_locality(radius), mass=0.0)
        m = measure_basic(f"F: radius {label}", "F_locality",
                          paths, vnorms, runtime)
        m.notes = f"stencil radius={radius} cells"
        out.append(m)
    return out


def run_G(fi: FrozenInputs):
    out = []
    cases = [
        ("G_T: transverse only",      G_transverse),
        ("G_Grad: gradient only",     G_gradient),
        ("G_Rot: rotational only",    G_rotational),
        ("G_T+Grad",                  G_T_plus_Grad),
        ("G_T+Rot",                   G_T_plus_Rot),
        ("G_Grad+Rot",                G_Grad_plus_Rot),
        ("G_T+Grad+Rot",              G_all_three),
    ]
    for label, fn in cases:
        paths, vnorms, runtime = propagate(fi, fn, mass=0.0)
        m = measure_basic(label, "G_decomp", paths, vnorms, runtime)
        m.notes = label.split(": ", 1)[-1]
        out.append(m)
    return out


def run_H(fi: FrozenInputs):
    out = []

    # H.1 — same operator but forced parallel
    paths, vnorms, runtime, perps = propagate_with_state(
        fi, H_force_parallel, normalize_every=1, record_perp=True)
    m = measure_basic("H.1: r projected onto v", "H_force_parallel",
                      paths, vnorms, runtime)
    m.perp_metric = perp_diagnostic(perps)
    m.notes = "force response parallel to velocity"
    out.append(m)

    # H.2 — same operator but forced perpendicular
    paths, vnorms, runtime, perps = propagate_with_state(
        fi, H_force_perpendicular, normalize_every=1, record_perp=True)
    m = measure_basic("H.2: r projected perp to v", "H_force_perp",
                      paths, vnorms, runtime)
    m.perp_metric = perp_diagnostic(perps)
    m.notes = "force response perpendicular to velocity"
    out.append(m)

    # H.3 — same operator, no normalization
    paths, vnorms, runtime, perps = propagate_with_state(
        fi, lambda ix, iy, f, s: ref_magnetic(ix, iy, f),
        normalize_every=0, record_perp=True)
    finite = all(np.isfinite(np.concatenate([p[0], p[1]])).all() for p in paths)
    drift = np.abs(np.concatenate(vnorms) - 1.0)
    devs = np.array([float(np.max(np.abs(y - y[0]))) for x, y in paths])
    m = M(label="H.3: no normalization", category="H_normalization",
          stable=False, finite=finite,
          bend_max=float(devs.max()), bend_mean=float(devs.mean()),
          conservation_residual=float(drift.max()),
          mean_velocity_drift=float(drift.mean()),
          runtime=runtime,
          perp_metric=perp_diagnostic(perps),
          notes="normalization removed")
    out.append(m)

    # H.4 — normalize every K steps
    for K in [2, 5, 10]:
        paths, vnorms, runtime, perps = propagate_with_state(
            fi, lambda ix, iy, f, s: ref_magnetic(ix, iy, f),
            normalize_every=K, record_perp=True)
        m = measure_basic(f"H.4: normalize every {K}", "H_normalization",
                          paths, vnorms, runtime)
        m.perp_metric = perp_diagnostic(perps)
        m.notes = f"normalize every {K} steps"
        out.append(m)

    # H.5 — perp diagnostic on the unmodified reference operator
    paths, vnorms, runtime, perps = propagate_with_state(
        fi, lambda ix, iy, f, s: ref_magnetic(ix, iy, f),
        normalize_every=1, record_perp=True)
    m = measure_basic("H.5: perp diagnostic on ref", "H_diagnostic",
                      paths, vnorms, runtime)
    m.perp_metric = perp_diagnostic(perps)
    m.notes = "mean |v.r|/(|v||r|) across all photon-steps"
    out.append(m)

    # H.6 — perp diagnostic on the parallel variant (Exp 1 baseline)
    paths, vnorms, runtime, perps = propagate_with_state(
        fi, lambda ix, iy, f, s: A_no_transverse(ix, iy, f),
        normalize_every=1, record_perp=True)
    m = measure_basic("H.6: perp diagnostic on parallel", "H_diagnostic",
                      paths, vnorms, runtime)
    m.perp_metric = perp_diagnostic(perps)
    m.notes = "mean |v.r|/(|v||r|) across all photon-steps, parallel rule"
    out.append(m)

    return out


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--output", type=Path, default=DEFAULT_OUT)
    a = p.parse_args()
    a.output.mkdir(parents=True, exist_ok=True)
    fi = load_inputs()

    all_measurements = []
    print("Reference (Exp 6 control)...")
    all_measurements.append(run_reference(fi))
    print("A...")
    all_measurements.append(run_A(fi))
    print("B...")
    all_measurements.append(run_B(fi))
    print("C...")
    all_measurements.append(run_C(fi))
    print("D...")
    all_measurements.append(run_D(fi))
    print("E...")
    all_measurements.extend(run_E(fi))
    print("F...")
    all_measurements.extend(run_F(fi))
    print("G...")
    all_measurements.extend(run_G(fi))
    print("H...")
    all_measurements.extend(run_H(fi))

    # Persist measurements
    rows = [asdict(m) for m in all_measurements]
    keys = list(rows[0].keys())
    with (a.output / "measurements.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        for r in rows:
            w.writerow(r)
    (a.output / "measurements.json").write_text(json.dumps(rows, indent=2))

    # Build the contribution table.
    table = build_contribution_table(all_measurements, fi)
    with (a.output / "contribution_table.csv").open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["Property", "Improves Bending", "Improves Conservation",
                    "Required", "Neutral", "Harmful"])
        for row in table:
            w.writerow(row)

    print(f"\nArtefacts written to {a.output}/")
    return 0


def build_contribution_table(measurements, fi):
    """Construct the contribution table from observed measurements.

    Each property is judged by direct numerical comparison vs the
    reference (Exp 6). No interpretation beyond what the numbers say.
    """
    by_label = {m.label: m for m in measurements}
    ref = by_label["Ref Exp 6 (control)"]

    table = []
    # A: removing transverse projection
    a = by_label.get("A: parallel (no transverse)")
    if a:
        bending = "no" if a.bend_max <= ref.bend_max * 0.5 else (
            "yes" if a.bend_max > ref.bend_max * 1.5 else "partial"
        )
        conservation = "no" if a.conservation_residual >= ref.conservation_residual * 1.5 else (
            "yes" if a.conservation_residual < ref.conservation_residual * 0.5 else "neutral"
        )
        table.append([
            "transverse projection (perp to grad C)",
            bending,
            "yes" if a.conservation_residual > ref.conservation_residual else "no",
            "yes" if abs(a.conservation_residual - ref.conservation_residual) < 1e-6 else "no",
            "" if abs(a.conservation_residual - ref.conservation_residual) < 1e-6 else "yes",
            "" if a.conservation_residual <= ref.conservation_residual * 2 else "yes",
        ])

    # B: transverse projector only
    b = by_label.get("B: transverse projector only")
    if b:
        table.append([
            "transverse projector (re-projection step)",
            "no" if b.bend_max < ref.bend_max * 0.5 else (
                "yes" if b.bend_max > ref.bend_max * 1.5 else "neutral"),
            "no" if b.conservation_residual > ref.conservation_residual else "neutral",
            "yes" if abs(b.conservation_residual - ref.conservation_residual) < 1e-6 else "no",
            "" if abs(b.conservation_residual - ref.conservation_residual) < 1e-6 else "yes",
            "" if b.conservation_residual <= ref.conservation_residual * 2 else "yes",
        ])

    # C: reversed
    c = by_label.get("C: reversed steering")
    if c:
        table.append([
            "sign convention (+90 vs -90 rotation)",
            "yes" if c.bend_max > ref.bend_max * 1.5 else (
                "no" if c.bend_max < ref.bend_max * 0.5 else "neutral"),
            "no" if c.conservation_residual > ref.conservation_residual * 1.5 else "neutral",
            "no",
            "yes" if abs(c.bend_max - ref.bend_max) < 0.1 * ref.bend_max else "no",
            "",
        ])

    # D: magnitude only
    d = by_label.get("D: magnitude only (no orientation)")
    if d:
        table.append([
            "directional orientation",
            "yes" if d.bend_max > ref.bend_max * 1.5 else (
                "no" if d.bend_max < ref.bend_max * 0.5 else "neutral"),
            "yes" if d.conservation_residual < ref.conservation_residual * 0.5 else (
                "no" if d.conservation_residual > ref.conservation_residual * 1.5 else "neutral"),
            "yes" if (d.bend_max > ref.bend_max * 1.5 and d.conservation_residual < ref.conservation_residual * 0.5) else "no",
            "" if not (d.bend_max > ref.bend_max * 1.5 and d.conservation_residual < ref.conservation_residual * 0.5) else "",
            "",
        ])

    # E: normalization
    e_labels = ["E: scale 1.00", "E: scale 0.75", "E: scale 0.50",
                "E: scale 0.25", "E: scale 0.10"]
    e_m = [by_label[l] for l in e_labels if l in by_label]
    if len(e_m) == 5:
        # If bending scales linearly with the scale factor and conservation
        # does not degrade, then magnitude sensitivity is the relevant property.
        bends = [m.bend_max for m in e_m]
        drifts = [m.conservation_residual for m in e_m]
        bending_scaling = "yes" if bends[0] > bends[-1] * 5 else "no"
        drift_growth = max(drifts) / min(drifts) if min(drifts) > 0 else 1.0
        table.append([
            "response magnitude scaling",
            bending_scaling,
            "no" if drift_growth > 5 else "neutral",
            "no",
            "yes" if bending_scaling == "no" else "no",
            "" if drift_growth < 10 else "yes",
        ])

    # F: locality
    f_labels = ["F: radius 1-cell", "F: radius von Neumann (5 cells)",
                "F: radius Moore (3x3)"]
    f_m = [by_label[l] for l in f_labels if l in by_label]
    if len(f_m) == 3:
        bends = [m.bend_max for m in f_m]
        drifts = [m.conservation_residual for m in f_m]
        locality_sensitive = "yes" if max(bends) / max(min(bends), 1e-30) > 2 else "no"
        table.append([
            "strictly local (1-cell) coupling",
            locality_sensitive,
            "no" if max(drifts) / max(min(drifts), 1e-30) > 5 else "neutral",
            "no" if locality_sensitive == "no" else "no",
            "yes" if locality_sensitive == "no" else "no",
            "",
        ])

    # G: decomposition
    g_cases = {
        "G_T: transverse only":    "transverse component alone",
        "G_Grad: gradient only":   "gradient component alone",
        "G_Rot: rotational only":  "rotational component alone",
        "G_T+Grad":                "transverse + gradient",
        "G_T+Rot":                 "transverse + rotational",
        "G_Grad+Rot":              "gradient + rotational",
        "G_T+Grad+Rot":            "all three combined",
    }
    for label, prop in g_cases.items():
        g = by_label.get(label)
        if not g:
            continue
        table.append([
            prop,
            "yes" if g.bend_max > ref.bend_max * 1.5 else (
                "no" if g.bend_max < ref.bend_max * 0.5 else "neutral"),
            "yes" if g.conservation_residual < ref.conservation_residual * 0.5 else (
                "no" if g.conservation_residual > ref.conservation_residual * 1.5 else "neutral"),
            "yes" if (label == "G_T: transverse only") else "no",
            "" if (label == "G_T: transverse only") else "yes",
            "" if g.conservation_residual <= ref.conservation_residual * 2 else "yes",
        ])

    # H: conservation analysis
    h_cases = {
        "H.1: r projected onto v":   "response aligned with v",
        "H.2: r projected perp to v":"response perpendicular to v",
        "H.3: no normalization":     "renormalization of v each step",
        "H.4: normalize every 2":    "renormalization frequency",
        "H.4: normalize every 5":    "renormalization frequency",
        "H.4: normalize every 10":   "renormalization frequency",
        "H.5: perp diagnostic on ref": "perpendicularity of response to v",
        "H.6: perp diagnostic on parallel": "perpendicularity of response to v",
    }
    for label, prop in h_cases.items():
        h = by_label.get(label)
        if not h:
            continue
        table.append([
            prop,
            "yes" if h.bend_max > ref.bend_max * 1.5 else (
                "no" if h.bend_max < ref.bend_max * 0.5 else "neutral"),
            "yes" if h.conservation_residual < ref.conservation_residual * 0.5 else (
                "no" if h.conservation_residual > ref.conservation_residual * 1.5 else "neutral"),
            "no",
            "yes" if abs(h.conservation_residual - ref.conservation_residual) < 1e-6 else "no",
            "" if h.conservation_residual <= ref.conservation_residual * 5 else "yes",
        ])

    return table


if __name__ == "__main__":
    raise SystemExit(main())
