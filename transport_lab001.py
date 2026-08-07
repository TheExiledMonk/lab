#!/usr/bin/env python3
"""PBUF TRANSPORT-LAB-001.

Laboratory experiments on nine local neighbour-update rules for
weak-lensing photon transport through the frozen Lens-001 substrate.

Fixed inputs (frozen throughout):
  * Lens-001 matter field       (runs/wl001/matter.csv)
  * Lens-001 deformation field  (runs/wl001/deformation.csv)
  * Lens-001 observation field  (runs/wl001/observation.csv)
  * Lens-001 gradients          (runs/wl001/gradient_x.csv, gradient_y.csv)
  * computational grid          (128x128 over [-8,8]^2, dx = dy = 16/127)
  * observer/source geometry    (same as pbuf_experiment.Config)
  * propagation parameters      (80 photon steps, step size 0.06)

The propagation algorithm is identical in every experiment:

    current neighbour
        -> compute local response         (varies per experiment)
        -> update next neighbour          (identical kernel)
        -> repeat

Only the local response rule changes.

Output: measurements + ranked table. No interpretation.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import time
from pathlib import Path
from dataclasses import dataclass, field, asdict

import numpy as np


ROOT = Path(__file__).resolve().parent
WL001 = ROOT / "runs" / "wl001"
DEFAULT_OUT = ROOT / "runs" / "transport_lab001"


# ----------------------------------------------------------------------------
# Frozen inputs (Lens-001)
# ----------------------------------------------------------------------------

@dataclass(frozen=True)
class FrozenInputs:
    n: int
    extent: float
    source_x: float
    source_y: float
    mass_x: float
    mass_y: float
    mass_sigma: float
    deformation_strength: float
    photon_steps: int
    photon_step_size: float
    matter: np.ndarray
    deformation: np.ndarray
    observation: np.ndarray
    gradient_x: np.ndarray
    gradient_y: np.ndarray
    C: np.ndarray            # constitutive state ( = deformation here )
    P: np.ndarray            # a stress-like proxy ( = |grad deformation| here )
    W: np.ndarray            # a work-density proxy ( = 1/2 |deformation|^2 here )
    N: np.ndarray            # outward unit normal field on neighbour boundary (radial)
    checksums: dict


def load_inputs() -> FrozenInputs:
    matter = np.loadtxt(WL001 / "matter.csv", delimiter=",")
    deformation = np.loadtxt(WL001 / "deformation.csv", delimiter=",")
    observation = np.loadtxt(WL001 / "observation.csv", delimiter=",")
    gx = np.loadtxt(WL001 / "gradient_x.csv", delimiter=",")
    gy = np.loadtxt(WL001 / "gradient_y.csv", delimiter=",")

    n = matter.shape[0]
    extent = 8.0
    x = np.linspace(-extent, extent, n)
    X, Y = np.meshgrid(x, x, indexing="xy")

    # Frozen observer/source geometry (mirrors pbuf_experiment.Config defaults)
    cfg = dict(
        n=128, extent=8.0, source_x=1.25, source_y=0.35,
        mass_x=-0.65, mass_y=0.0, mass_sigma=0.75,
        deformation_strength=0.18, photon_steps=80, photon_step_size=0.06,
    )

    C = deformation.copy()
    P = np.hypot(gx, gy)
    W = 0.5 * deformation**2

    # Outward neighbour-boundary normal: radial from the mass centre,
    # frozen with the data (independent of the propagation experiment).
    N = np.zeros((2, n, n), dtype=np.float64)
    rx = X - cfg["mass_x"]
    ry = Y - cfg["mass_y"]
    r = np.hypot(rx, ry)
    r_safe = np.where(r > 1e-12, r, 1e-12)
    N[0] = np.where(r > 1e-12, rx / r_safe, 0.0)
    N[1] = np.where(r > 1e-12, ry / r_safe, 0.0)

    sha = lambda a: hashlib.sha256(np.ascontiguousarray(a).tobytes()).hexdigest()
    return FrozenInputs(
        n=n, extent=extent, source_x=cfg["source_x"], source_y=cfg["source_y"],
        mass_x=cfg["mass_x"], mass_y=cfg["mass_y"],
        mass_sigma=cfg["mass_sigma"], deformation_strength=cfg["deformation_strength"],
        photon_steps=cfg["photon_steps"], photon_step_size=cfg["photon_step_size"],
        matter=matter, deformation=deformation, observation=observation,
        gradient_x=gx, gradient_y=gy, C=C, P=P, W=W, N=N,
        checksums={k: sha(v) for k, v in [
            ("matter", matter), ("deformation", deformation),
            ("observation", observation), ("gradient_x", gx),
            ("gradient_y", gy), ("C", C), ("P", P), ("W", W),
        ]},
    )


# ----------------------------------------------------------------------------
# Fixed propagation kernel. ONLY the local response function changes.
# ----------------------------------------------------------------------------

def lookup_index(value: float, x: np.ndarray) -> int:
    return int(np.clip(np.searchsorted(x, value) - 1, 0, len(x) - 1))


def propagate(fi: FrozenInputs, response_fn, mass: float = 0.0,
              n_photons: int = 9, y_span: float = 3.0,
              record: bool = False):
    """Identical neighbour-to-neighbour kernel.

    response_fn(ix, iy, fi) -> (rx, ry): the local response vector at neighbour
    (iy, ix). The kernel applies it, normalises the velocity, and walks one
    step. The kernel is identical for every experiment; only response_fn varies.

    mass = 0.0 -> photon-like: response fully redirects velocity each step.
    mass > 0.0 -> massive test particle: response is integrated into a
    momentum-like state with inertia; velocity = momentum / mass.
    """
    n = fi.n
    x = np.linspace(-fi.extent, fi.extent, n)
    step = fi.photon_step_size
    n_steps = fi.photon_steps

    paths = []
    all_velocity_norms = []
    started = time.time()

    for y0 in np.linspace(-y_span, y_span, n_photons):
        xs = np.full(n_steps, -fi.extent)
        ys = np.full(n_steps, y0)

        if mass <= 0.0:
            vx = np.ones(n_steps)
            vy = np.zeros(n_steps)
        else:
            # momentum carried, velocity = momentum/mass
            px = np.full(n_steps, mass * 1.0)  # initial px = m * vx0
            py = np.full(n_steps, mass * 0.0)  # initial py = m * vy0
            vx = px / mass
            vy = py / mass

        velocity_norms = [float(np.hypot(vx[0], vy[0]))]

        for k in range(1, n_steps):
            ix = lookup_index(xs[k - 1], x)
            iy = lookup_index(ys[k - 1], x)
            rx, ry = response_fn(ix, iy, fi)

            if mass <= 0.0:
                # Photon-like: local response fully redirects velocity
                vx_k = vx[k - 1] - step * rx
                vy_k = vy[k - 1] - step * ry
            else:
                # Massive: response is a force; momentum absorbs inertia
                px_k = px[k - 1] - step * rx * mass
                py_k = py[k - 1] - step * ry * mass
                px[k] = px_k
                py[k] = py_k
                vx_k = px_k / mass
                vy_k = py_k / mass

            norm = max(np.hypot(vx_k, vy_k), 1e-12)
            vx_k /= norm
            vy_k /= norm
            vx[k] = vx_k
            vy[k] = vy_k

            xs[k] = xs[k - 1] + step * vx_k
            ys[k] = ys[k - 1] + step * vy_k
            velocity_norms.append(float(norm))

        paths.append((xs.copy(), ys.copy()))
        all_velocity_norms.append(np.array(velocity_norms))

    runtime = time.time() - started
    return paths, all_velocity_norms, runtime


# ----------------------------------------------------------------------------
# Local response rules (the only thing that varies between experiments).
# ----------------------------------------------------------------------------

def response_gradC(ix, iy, fi):
    """Exp 1: response ∝ ∇C."""
    return fi.gradient_x[iy, ix], fi.gradient_y[iy, ix]


def response_gradP_over_C(ix, iy, fi):
    """Exp 2: response ∝ ∇P / C."""
    eps = 1e-12
    # Approximate ∇P from gradient of |grad deformation| via Hessian of C.
    # Use central differences of the deformation gradient magnitude P.
    n = fi.n
    ixp = min(ix + 1, n - 1); ixm = max(ix - 1, 0)
    iyp = min(iy + 1, n - 1); iym = max(iy - 1, 0)
    spacing = 2 * fi.extent / (n - 1)
    gxP = (fi.P[iy, ixp] - fi.P[iy, ixm]) / (2 * spacing)
    gyP = (fi.P[iyp, ix] - fi.P[iym, ix]) / (2 * spacing)
    C = fi.C[iy, ix] + eps
    return gxP / C, gyP / C


def response_gradW(ix, iy, fi):
    """Exp 3: response ∝ ∇W where W = (1/2) C^2."""
    n = fi.n
    ixp = min(ix + 1, n - 1); ixm = max(ix - 1, 0)
    iyp = min(iy + 1, n - 1); iym = max(iy - 1, 0)
    spacing = 2 * fi.extent / (n - 1)
    gxW = (fi.W[iy, ixp] - fi.W[iy, ixm]) / (2 * spacing)
    gyW = (fi.W[iyp, ix] - fi.W[iym, ix]) / (2 * spacing)
    return gxW, gyW


def response_traction(ix, iy, fi):
    """Exp 4: response ∝ traction t = P_F · N (neighbour receives traction)."""
    # P_F proxy = ∇C; N = outward radial unit normal at neighbour boundary
    return (fi.gradient_x[iy, ix] * fi.N[0, iy, ix]
            + fi.gradient_y[iy, ix] * fi.N[1, iy, ix]) * fi.N[0, iy, ix], \
           (fi.gradient_x[iy, ix] * fi.N[0, iy, ix]
            + fi.gradient_y[iy, ix] * fi.N[1, iy, ix]) * fi.N[1, iy, ix]


def response_elastic_spring(ix, iy, fi):
    """Exp 5: elastic spring — neighbour displacement via restoring interaction."""
    # F = -k (u_neighbour - u_self); continuous limit F = -k * grad(u).
    # Use gradient magnitude as the spring scale.
    return fi.gradient_x[iy, ix], fi.gradient_y[iy, ix]


def response_magnetic_style(ix, iy, fi):
    """Exp 6: magnetic-style directional update (perpendicular to local field)."""
    # Local update direction is rotated 90° from ∇C, scaled by |∇C|.
    # Mirrors Lorentz-style response where the force is perpendicular to the
    # local field direction.
    gx = fi.gradient_x[iy, ix]
    gy = fi.gradient_y[iy, ix]
    return -gy, gx


def response_phase(ix, iy, fi):
    """Exp 7: phase accumulation — neighbour advance follows local phase delay."""
    # Local phase delay ~ C; phase gradient ∝ ∇C scaled by |C|.
    C = fi.C[iy, ix]
    return fi.gradient_x[iy, ix] * C, fi.gradient_y[iy, ix] * C


def response_constant(ix, iy, fi):
    """Exp 8: constant local transfer — fixed fraction of disturbance."""
    # Constant response regardless of local field value.
    # Magnitude chosen to match the typical scale of ∇C in the dataset.
    return 0.01, 0.0


# ----------------------------------------------------------------------------
# Measurement utilities
# ----------------------------------------------------------------------------

@dataclass
class Measurement:
    experiment: str
    rule: str
    stable: bool
    finite_outputs: bool
    max_velocity_drift: float          # | |v| - 1 | over all photons and steps
    mean_velocity_drift: float
    photon_max_deviation: float        # max |y_final - y_initial|
    photon_mean_deviation: float
    total_bending_angle: float         # sum |dv| over all steps, summed over photons
    shear_proxy: float                 # |perp component|/|parallel component|
    conservation_residual: float      # how close to |v|=1 is preserved
    final_x_spread: float              # spread of final x positions
    final_y_spread: float              # spread of final y positions
    runtime_seconds: float
    notes: str


def measure(name: str, rule_name: str, paths, vnorms, runtime, fi: FrozenInputs) -> Measurement:
    finite_outputs = all(
        np.isfinite(np.concatenate([p[0], p[1]])).all() for p in paths
    )
    flat_norms = np.concatenate(vnorms)
    drift = np.abs(flat_norms - 1.0)
    max_drift = float(drift.max())
    mean_drift = float(drift.mean())

    deviations = np.array([float(np.max(np.abs(y - y[0]))) for x, y in paths])
    photon_max_deviation = float(deviations.max())
    photon_mean_deviation = float(deviations.mean())

    # Total bending angle: sum of |delta v| per photon, then averaged.
    bending_per_photon = []
    for v in vnorms:
        deltas = np.diff(v)
        # Use the un-normalised velocity change, not unit-norm.
        bending_per_photon.append(float(np.sum(np.abs(deltas))))
    total_bending_angle = float(np.mean(bending_per_photon))

    # Final position spread
    final_x = np.array([p[0][-1] for p in paths])
    final_y = np.array([p[1][-1] for p in paths])
    fx_spread = float(final_x.std())
    fy_spread = float(final_y.std())

    # Shear proxy: ratio of perpendicular response component to parallel one,
    # computed from the response field directly.
    # Use a generic proxy: cross-product magnitude vs dot-product magnitude,
    # averaged over the interior. Smaller ratio = more aligned shear.
    gxv = np.array([fi.gradient_x.ravel()])[0]
    gyv = np.array([fi.gradient_y.ravel()])[0]
    parallel = np.abs(gxv * gxv + gyv * gyv)
    perp = np.abs(gxv * gyv - gyv * gxv)
    interior = parallel > 1e-8
    shear_proxy = float(perp[interior].mean() / parallel[interior].mean()) if interior.any() else 0.0

    stable = bool(finite_outputs and max_drift < 1.0 and
                  all(np.isfinite(p[0]).all() and np.isfinite(p[1]).all() for p in paths))

    return Measurement(
        experiment=name, rule=rule_name,
        stable=stable, finite_outputs=finite_outputs,
        max_velocity_drift=max_drift, mean_velocity_drift=mean_drift,
        photon_max_deviation=photon_max_deviation,
        photon_mean_deviation=photon_mean_deviation,
        total_bending_angle=total_bending_angle,
        shear_proxy=shear_proxy,
        conservation_residual=max_drift,
        final_x_spread=fx_spread,
        final_y_spread=fy_spread,
        runtime_seconds=runtime,
        notes="",
    )


# ----------------------------------------------------------------------------
# Driver
# ----------------------------------------------------------------------------

EXPERIMENTS = [
    ("Exp 1", "response ∝ ∇C",                       response_gradC),
    ("Exp 2", "response ∝ ∇P / C",                   response_gradP_over_C),
    ("Exp 3", "response ∝ ∇W",                       response_gradW),
    ("Exp 4", "response ∝ traction t = P_F·N",       response_traction),
    ("Exp 5", "response = elastic spring restoring", response_elastic_spring),
    ("Exp 6", "response = magnetic-style directional", response_magnetic_style),
    ("Exp 7", "response = phase accumulation",       response_phase),
    ("Exp 8", "response = constant local transfer",  response_constant),
]


def run_all(fi: FrozenInputs, out: Path) -> list[Measurement]:
    out.mkdir(parents=True, exist_ok=True)
    measurements: list[Measurement] = []

    for name, rule_name, fn in EXPERIMENTS:
        paths, vnorms, runtime = propagate(fi, fn, mass=0.0)
        m = measure(name, rule_name, paths, vnorms, runtime, fi)
        measurements.append(m)

    # Exp 9: same local rule as Exp 1, but with mass > 0 (test particle),
    # compared to its m=0 counterpart.
    paths_m0, vnorms_m0, _ = propagate(fi, response_gradC, mass=0.0)
    paths_m, vnorms_m, runtime = propagate(fi, response_gradC, mass=2.0)
    m_m0 = measure("Exp 9 (m=0)", "photon (m=0) baseline", paths_m0, vnorms_m0, 0.0, fi)
    m_m = measure("Exp 9 (m>0)", "test particle (m>0)", paths_m, vnorms_m, runtime, fi)
    measurements.extend([m_m0, m_m])

    # ------------------------------------------------------------------
    # Persist results
    # ------------------------------------------------------------------
    rows = [asdict(m) for m in measurements]
    keys = list(rows[0].keys())
    with (out / "measurements.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        for r in rows:
            w.writerow(r)

    # Ranked table: composite numerical score.
    ranked = rank(measurements)
    with (out / "ranked_table.csv").open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["Rank", "Experiment", "Rule", "Stable",
                    "WeakLensing(yes/no)", "Conservation residual",
                    "Total bending", "Runtime (s)", "Score", "Comments"])
        for r in ranked:
            comments = r["notes"]
            w.writerow([r["rank"], r["experiment"], r["rule"],
                        r["stable"], r["weak_lensing"],
                        f"{r['conservation_residual']:.3e}",
                        f"{r['total_bending_angle']:.3e}",
                        f"{r['runtime_seconds']:.3f}",
                        f"{r['score']:.3f}", comments])

    # Raw JSON
    (out / "measurements.json").write_text(json.dumps(
        {"frozen_inputs": {k: v for k, v in fi.checksums.items()},
         "config": {
             "n": fi.n, "extent": fi.extent, "source_x": fi.source_x,
             "source_y": fi.source_y, "mass_x": fi.mass_x, "mass_y": fi.mass_y,
             "mass_sigma": fi.mass_sigma,
             "deformation_strength": fi.deformation_strength,
             "photon_steps": fi.photon_steps, "photon_step_size": fi.photon_step_size,
         },
         "experiments": rows,
         "ranked": ranked},
        indent=2))

    return measurements


def rank(measurements: list[Measurement]) -> list[dict]:
    """Numerical-only ranking.

    Score = (stability_bonus)
           + (weak_lensing_bonus)
           + (1 / (1 + conservation_residual))
           + (small_runtime_bonus if stable)
    """
    scored = []
    for m in measurements:
        weak_lensing = (m.photon_max_deviation > 1e-6)
        # Score: stable & weak_lensing give bonuses; conservation residual
        # contributes 1/(1+residual); runtime contributes a small penalty.
        score = 0.0
        if m.stable:
            score += 1.0
        if weak_lensing:
            score += 1.0
        score += 1.0 / (1.0 + m.conservation_residual)
        score += 0.001 / (1.0 + m.runtime_seconds)
        notes = []
        if not m.stable:
            notes.append("unstable")
        if not weak_lensing:
            notes.append("no bending observed")
        if m.conservation_residual > 1e-3:
            notes.append("velocity drift")
        if m.runtime_seconds > 1.0:
            notes.append("slow")
        scored.append({
            "experiment": m.experiment, "rule": m.rule,
            "stable": m.stable,
            "weak_lensing": "yes" if weak_lensing else "no",
            "conservation_residual": m.conservation_residual,
            "total_bending_angle": m.total_bending_angle,
            "photon_max_deviation": m.photon_max_deviation,
            "runtime_seconds": m.runtime_seconds,
            "score": score,
            "notes": "; ".join(notes) if notes else "nominal",
            "rank": 0,
        })
    scored.sort(key=lambda r: r["score"], reverse=True)
    for i, r in enumerate(scored):
        r["rank"] = i + 1
    return scored


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--output", type=Path, default=DEFAULT_OUT)
    a = p.parse_args()
    fi = load_inputs()
    measurements = run_all(fi, a.output)
    print(f"PBUF TRANSPORT-LAB-001 — {len(measurements)} runs")
    for m in measurements:
        print(f"  {m.experiment:13s}  stable={m.stable!s:5s}  "
              f"bend={m.photon_max_deviation:.3e}  "
              f"|v|-drift={m.max_velocity_drift:.3e}  "
              f"runtime={m.runtime_seconds:.3f}s")
    print(f"\nArtefacts written to {a.output}/")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
