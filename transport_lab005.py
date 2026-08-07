#!/usr/bin/env python3
"""PBUF TRANSPORT-LAB-005 — Which constitutive direction governs the response?

Frozen: neighbour-to-neighbour transport, Lens-001, kernel, 90° response.
Variable: the reference direction g to which the response is constrained
to be perpendicular. The response is r = R_90(g) = (-g_y, g_x).

Candidates:
  1. Constitutive gradient  g = ∇C                 (control)
  2. Stress gradient        g = ∇(|∇C|) / (C + ε)
  3. Energy gradient        g = ∇(½ C²) = C · ∇C
  4. Traction direction     g = (∇C · N̂) N̂
  5. Force density          g = ∇(∇² C)
  6. Principal strain       g = principal eigenvector of ∇⊗∇C
  7. Principal stress       g = principal eigenvector of ∇⊗∇(½ C²)
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
from transport_lab003 import propagate_diag


ROOT = Path(__file__).resolve().parent
DEFAULT_OUT = ROOT / "runs" / "transport_lab005"


# ----------------------------------------------------------------------------
# Per-candidate reference vector fields.
# ----------------------------------------------------------------------------

def candidate_constitutive_gradient(fi: FrozenInputs):
    return fi.gradient_x.copy(), fi.gradient_y.copy()


def candidate_stress_gradient(fi: FrozenInputs):
    """g = ∇(|∇C|) / (C + ε).

    P := |∇C| (gradient magnitude). Then g = (∇P) / C, regularised.
    """
    n = fi.n
    extent = fi.extent
    spacing = 2 * extent / (n - 1)
    P = np.hypot(fi.gradient_x, fi.gradient_y)
    gx_P = np.zeros_like(P)
    gy_P = np.zeros_like(P)
    gx_P[:, 1:-1] = (P[:, 2:] - P[:, :-2]) / (2 * spacing)
    gy_P[1:-1, :] = (P[2:, :] - P[:-2, :]) / (2 * spacing)
    eps = 1e-3
    C_reg = fi.deformation + eps
    return gx_P / C_reg, gy_P / C_reg


def candidate_energy_gradient(fi: FrozenInputs):
    """g = ∇(½ C²) = C · ∇C."""
    return fi.deformation * fi.gradient_x, fi.deformation * fi.gradient_y


def candidate_traction(fi: FrozenInputs):
    """g = (∇C · N̂) N̂. Use the radial N̂ from transport_lab001."""
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
    # t = (∇C · N̂) N̂
    cos_t = fi.gradient_x * Nx + fi.gradient_y * Ny
    return cos_t * Nx, cos_t * Ny


def candidate_force_density(fi: FrozenInputs):
    """g = ∇(∇² C)."""
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
    return gx, gy


def _principal_eigenvector_field(fi: FrozenInputs, F: np.ndarray,
                                 G: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """For a 2x2 symmetric tensor T_ij = F_i F_j + G H_ij (Hessian of a scalar),
    return the principal eigenvector at each grid point.
    Here we use: T = (∇F)(∇F)^T + α·Hessian(G). F is a scalar that we
    differentiate, G is the scalar whose Hessian we add.
    """
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

    # Tensor: T = (∇F)(∇F)^T + α·Hessian(G). Use α = 1.
    alpha = 1.0
    Txx = gx_F * gx_F + alpha * Hxx
    Tyy = gy_F * gy_F + alpha * Hyy
    Txy = gx_F * gy_F + alpha * Hxy

    # Principal eigenvector of [[Txx, Txy], [Txy, Tyy]].
    # Eigenvalues: λ = (trace ± sqrt(trace² - 4 det)) / 2.
    # Eigenvector for λ_max: (Txy, λ_max - Txx) or similar.
    trace = Txx + Tyy
    det = Txx * Tyy - Txy * Txy
    disc = np.sqrt(np.maximum(trace * trace - 4 * det, 0.0))
    lam_max = (trace + disc) / 2.0
    # v = (Txy, λ_max - Txx); if zero, fall back to (1, 0).
    vx = Txy
    vy = lam_max - Txx
    norm = np.hypot(vx, vy)
    safe = norm > 1e-15
    vx_safe = np.where(safe, vx / np.where(safe, norm, 1.0), 1.0)
    vy_safe = np.where(safe, vy / np.where(safe, norm, 1.0), 0.0)
    # Sign consistency: align with the gradient direction.
    align = np.sign(vx_safe * gx_F + vy_safe * gy_F)
    align = np.where(align == 0, 1.0, align)
    vx_safe *= align
    vy_safe *= align
    return vx_safe, vy_safe


def candidate_principal_strain(fi: FrozenInputs):
    """g = principal eigenvector of (∇C)(∇C)^T + ∇⊗∇C."""
    return _principal_eigenvector_field(fi, fi.deformation, fi.deformation)


def candidate_principal_stress(fi: FrozenInputs):
    """g = principal eigenvector of (∇W)(∇W)^T + ∇⊗∇W, where W = ½C²."""
    W = 0.5 * fi.deformation ** 2
    return _principal_eigenvector_field(fi, W, W)


# ----------------------------------------------------------------------------
# Build the response = R_90(g) for a given candidate.
# ----------------------------------------------------------------------------

def make_response_from_field(gx_field, gy_field):
    def fn(ix, iy, fi: FrozenInputs, state):
        return -gy_field[iy, ix], gx_field[iy, ix], "subtract"
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
    reference_magnitude_mean: float
    reference_magnitude_max: float
    stable: bool
    runtime: float


def measure(candidate, label, gx_f, gy_f, fi) -> M:
    response = make_response_from_field(gx_f, gy_f)
    paths, diag, runtime = propagate_diag(fi, response)
    devs = np.array([float(np.max(np.abs(y - y[0]))) for x, y in paths])
    pre = np.concatenate(diag["speed_drift_pre"])
    dirs = np.concatenate(diag["direction_drift"])
    errs = np.concatenate(diag["pos_error"])
    finite = all(np.isfinite(np.concatenate([p[0], p[1]])).all() for p in paths)
    stable = bool(finite and pre.max() < 1.0)
    ref_mag = np.hypot(gx_f, gy_f)
    return M(
        candidate=candidate, reference_label=label,
        bend_max=float(devs.max()), bend_mean=float(devs.mean()),
        conservation=float(pre.max()),
        speed_drift_pre_max=float(pre.max()),
        speed_drift_pre_mean=float(pre.mean()),
        direction_drift_mean=float(dirs.mean()),
        position_error=float(errs.sum()),
        reference_magnitude_mean=float(ref_mag.mean()),
        reference_magnitude_max=float(ref_mag.max()),
        stable=stable, runtime=runtime,
    )


CANDIDATES = [
    ("Candidate 1", "Constitutive gradient  ∇C",                candidate_constitutive_gradient),
    ("Candidate 2", "Stress gradient  ∇|∇C| / (C + ε)",         candidate_stress_gradient),
    ("Candidate 3", "Energy gradient  ∇(½C²) = C·∇C",           candidate_energy_gradient),
    ("Candidate 4", "Traction  (∇C·N̂)N̂",                       candidate_traction),
    ("Candidate 5", "Force density  ∇(∇²C)",                    candidate_force_density),
    ("Candidate 6", "Principal strain direction",                candidate_principal_strain),
    ("Candidate 7", "Principal stress direction",                candidate_principal_stress),
]


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--output", type=Path, default=DEFAULT_OUT)
    a = p.parse_args()
    a.output.mkdir(parents=True, exist_ok=True)
    fi = load_inputs()

    measurements = []
    for cand, label, fn in CANDIDATES:
        gx, gy = fn(fi)
        print(f"{cand}: {label}  |g|max={np.hypot(gx, gy).max():.3e}")
        m = measure(cand, label, gx, gy, fi)
        measurements.append(m)
        print(f"  bend_max={m.bend_max:.4e}  cons={m.conservation:.4e}  "
              f"stable={m.stable}  runtime={m.runtime:.3f}s")

    rows = [asdict(m) for m in measurements]
    keys = list(rows[0].keys())
    with (a.output / "measurements.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        for r in rows:
            w.writerow(r)
    (a.output / "measurements.json").write_text(json.dumps(rows, indent=2))

    # Ranking
    control = measurements[0]
    score = []
    for m in measurements:
        # Higher bending is better, lower conservation is better. Normalise.
        # Score = bending / max_bending + min_cons / conservation.
        bends = np.array([x.bend_max for x in measurements])
        conservs = np.array([x.conservation for x in measurements])
        b_max = bends.max()
        c_min = conservs.min()
        c_eff = m.conservation if m.conservation > 0 else 1e-30
        s = (m.bend_max / b_max) + (c_min / c_eff)
        score.append(s)
    order = np.argsort(score)[::-1]
    rank_table = []
    for rank_idx, idx in enumerate(order):
        m = measurements[idx]
        rank_table.append({
            "rank": rank_idx + 1,
            "candidate": m.candidate,
            "reference_label": m.reference_label,
            "bend_max": m.bend_max,
            "conservation": m.conservation,
            "speed_drift_pre_max": m.speed_drift_pre_max,
            "stable": m.stable,
            "score": float(score[idx]),
        })

    with (a.output / "ranking.csv").open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["Rank", "Candidate", "Reference Direction", "Bend",
                    "Conservation", "Speed Drift", "Stable", "Score"])
        for r in rank_table:
            w.writerow([r["rank"], r["candidate"], r["reference_label"],
                        f"{r['bend_max']:.4e}", f"{r['conservation']:.4e}",
                        f"{r['speed_drift_pre_max']:.4e}",
                        "yes" if r["stable"] else "no",
                        f"{r['score']:.3f}"])

    # Relative comparison vs control
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
                        f"{r['bend_pct_vs_control']:+.2f}%",
                        f"{r['cons_pct_vs_control']:+.2f}%",
                        f"{r['speed_pct_vs_control']:+.2f}%"])

    # Correlation analysis
    bends = np.array([m.bend_max for m in measurements])
    conservs = np.array([m.conservation for m in measurements])
    speeds = np.array([m.speed_drift_pre_max for m in measurements])
    cor = {}
    if bends.std() > 0 and conservs.std() > 0:
        cor["bend_vs_conservation_pearson"] = float(np.corrcoef(bends, conservs)[0, 1])
    if bends.std() > 0 and speeds.std() > 0:
        cor["bend_vs_speed_pearson"] = float(np.corrcoef(bends, speeds)[0, 1])
    if conservs.std() > 0 and speeds.std() > 0:
        cor["conservation_vs_speed_pearson"] = float(np.corrcoef(conservs, speeds)[0, 1])
    (a.output / "correlations.json").write_text(json.dumps(cor, indent=2))

    print(f"\n=== Ranking ===")
    for r in rank_table:
        print(f"  #{r['rank']}: {r['candidate']} ({r['reference_label']})  "
              f"bend={r['bend_max']:.3e}  cons={r['conservation']:.3e}  score={r['score']:.3f}")
    print(f"\n=== Correlations ===")
    for k, v in cor.items():
        print(f"  {k}: {v:+.4f}")
    print(f"\nArtefacts written to {a.output}/")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
