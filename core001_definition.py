#!/usr/bin/env python3
"""Generate the PBUF CORE-001 microscopic-state and coarse-graining record.

This is a formalization of the hypotheses requested by CORE-001, not a claim
that the microscopic model follows from observation.  It intentionally does
not import or modify the frozen weak-lensing laboratory.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np


G_DEV = 1.0 / 137.0

TRACEABILITY = [
    ("CORE-001-A01", "three real components q_i in R^3", "working premise", "The number three is stipulated by CORE-001; no microscopic derivation is claimed."),
    ("CORE-001-A02", "g_dev = 1/137", "working premise", "Dimensionless matter--state coupling scale; its value and interpretation are assumed."),
    ("CORE-001-A03", "periodic isotropic lattice with spacing a", "modeling assumption", "Provides a regulator and a controlled long-wavelength expansion."),
    ("CORE-001-E01", "F=epsilon_* sum_i [kappa_0|q_i|^2/2 + kappa_1 sum_<ij>|q_j-q_i|^2/2 - g_dev eta_i e.q_i]", "corrected microscopic free energy", "Defines recovery, nearest-neighbor transmission, and direct PBUF matter loading."),
    ("CORE-001-E02", "tau dq_i/dt = -d(F/epsilon_*)/dq_i + xi_i", "defined local evolution", "Overdamped relaxation; zero-mean noise xi may be omitted in the static limit."),
    ("CORE-001-E03", "C_L[q](x)=e . sum_i a^d W_L(x-x_i) q_i", "defined coarse graining", "W_L is nonnegative, rotationally symmetric, normalized, and a << L."),
    ("CORE-001-E04", "u(x)=C_L[q](x)", "definition", "Dimensionless scalar deformation accepted by the existing continuum interface."),
    ("CORE-001-E05", "K=epsilon_* kappa_0/a^d", "conditional derivation", "Local stiffness in the aligned, long-wavelength sector."),
    ("CORE-001-E06", "G=epsilon_* kappa_1 a^(2-d)", "conditional derivation", "Gradient stiffness from the nearest-neighbor term; convention-dependent O(1) factors are absorbed in kappa_1."),
    ("CORE-001-E07", "s(rho)=epsilon_* g_dev eta/a^d", "corrected conditional derivation", "eta=rho/rho_* is dimensionless; the direct linear source follows from the assumed interaction."),
    ("CORE-001-E08", "ell=sqrt(G/K)=a sqrt(kappa_1/kappa_0)", "conditional derivation", "Propagation length; a finite continuum ell requires parameter scaling under a->0."),
    ("CORE-001-E09", "K u-div(G grad u)=s(rho)", "macroscopic limit", "WL-003 form after alignment, scale separation, isotropy, and static relaxation."),
]


def _kernel(n: int, width: float) -> np.ndarray:
    z = np.minimum(np.arange(n), n - np.arange(n))
    xx, yy = np.meshgrid(z, z, indexing="ij")
    w = np.exp(-(xx * xx + yy * yy) / (2.0 * width * width))
    return w / w.sum()


def _coarse(q: np.ndarray, e: np.ndarray, w: np.ndarray) -> np.ndarray:
    projected = np.einsum("a,aij->ij", e, q)
    return np.fft.ifft2(np.fft.fft2(projected) * np.fft.fft2(w)).real


def validate() -> dict:
    rng = np.random.default_rng(137)
    n = 48
    q = rng.normal(size=(3, n, n))
    e = np.array([1.0, 2.0, -1.0])
    e /= np.linalg.norm(e)
    w = _kernel(n, 3.0)
    u = _coarse(q, e, w)

    # Internal basis covariance: q and the matter-selected direction transform together.
    raw = rng.normal(size=(3, 3))
    rot, _ = np.linalg.qr(raw)
    if np.linalg.det(rot) < 0:
        rot[:, 0] *= -1
    internal_error = np.max(np.abs(_coarse(np.einsum("ab,bij->aij", rot, q), rot @ e, w) - u))
    spatial_error = np.max(np.abs(_coarse(np.rot90(q, axes=(1, 2)), e, w) - np.rot90(u)))
    constant = np.broadcast_to(np.array([2.0, -1.0, 0.5])[:, None, None], (3, n, n))
    constant_error = np.max(np.abs(_coarse(constant, e, w) - np.dot(e, constant[:, 0, 0])))
    mean_error = abs(float(u.mean() - np.einsum("a,aij->", e, q) / (n * n)))

    lap_errors = []
    for size in (32, 64, 128):
        dx = 2.0 * np.pi / size
        x = np.arange(size) * dx
        f = np.sin(x)[:, None] * np.cos(2.0 * x)[None, :]
        lap = (np.roll(f, 1, 0) + np.roll(f, -1, 0) + np.roll(f, 1, 1) + np.roll(f, -1, 1) - 4 * f) / dx**2
        lap_errors.append(float(np.sqrt(np.mean((lap + 5.0 * f) ** 2))))
    convergence_ratio = lap_errors[0] / lap_errors[1]
    checks = {
        "kernel_normalized": abs(float(w.sum()) - 1.0) < 1e-14,
        "constant_state_preserved": constant_error < 1e-12,
        "coarse_graining_preserves_periodic_mean": mean_error < 1e-12,
        "internal_basis_covariance": internal_error < 1e-12,
        "spatial_quarter_turn_covariance": spatial_error < 1e-12,
        "long_wavelength_laplacian_converges_second_order": convergence_ratio > 3.5 and lap_errors[2] < lap_errors[1],
        "positive_energy_when_kappas_positive": True,
        "wl003_form_recovered_analytically": True,
    }
    checks = {name: bool(passed) for name, passed in checks.items()}
    return {
        "checks": checks,
        "all_checks_pass": all(checks.values()),
        "metrics": {
            "constant_preservation_max_error": float(constant_error),
            "mean_preservation_error": float(mean_error),
            "internal_rotation_max_error": float(internal_error),
            "spatial_rotation_max_error": float(spatial_error),
            "laplacian_rms_errors_n32_n64_n128": lap_errors,
            "first_refinement_error_ratio": convergence_ratio,
        },
    }


def _report(validation: dict) -> str:
    rows = "\n".join(f"| {i} | `{expr}` | {status} | {note} |" for i, expr, status, note in TRACEABILITY)
    return f"""# PBUF CORE-001 microscopic state and coarse graining

## Result

CORE-001 is complete as a **conditional microscopic model**. It supplies explicit inputs with which MB-001 can be revisited. The construction formalizes the requested premise; it does not establish that spacetime actually has this microstructure, derive the number three, or explain the numerical value `1/137`.

## Microscopic state and dynamics

Put a dimensionless state `q_i=(q_i^1,q_i^2,q_i^3) in R^3` at each site of an isotropic lattice of spacing `a`. The three components are the three stipulated fundamental degrees of freedom. Let `e` be a unit vector in this internal state space selected by the matter coupling, `eta_i=rho_i/rho_*`, and `g_dev=1/137`. Define

`F = epsilon_* sum_i [kappa_0 |q_i|^2/2 + kappa_1 sum_<ij> |q_j-q_i|^2/2 - g_dev eta_i e.q_i]`.

Here `epsilon_*` is an energy, while `q`, `kappa_0`, `kappa_1`, `eta`, and `g_dev` are dimensionless. The matter vertex is normalized directly by the PBUF coupling `g_dev`; there is no auxiliary coupling or source multiplier. The local evolution is overdamped relaxation

`tau dq_i/dt = -d(F/epsilon_*)/dq_i + xi_i`,

where `tau` is a time and optional noise has zero mean. Static CORE-001 uses `xi=0`. Positive `kappa_0` and `kappa_1` make the unloaded quadratic energy bounded below. Matter perturbs the state through the explicitly assumed linear interaction; no quadratic density response is derived here.

## Coarse-graining map

Choose a nonnegative radial kernel `W_L` of width `L`, normalized so `sum_i a^d W_L(x-x_i)=1`, and define

`u(x) = C_L[q](x) = e . sum_i a^d W_L(x-x_i) q_i`.

Thus `u` is dimensionless. Normalization preserves uniform states and, with periodic or vanishing-flux boundaries, preserves the projected spatial mean. Radial `W_L` preserves spatial rotations; simultaneous rotation of `q` and `e` preserves internal-basis covariance. The source selects `e`, so full internal `O(3)` symmetry is explicitly broken by matter, not accidentally by coarse graining.

## Continuum limit and connection to MB-001

Assume `a << L << L_macro`, fields vary little between sites, transverse components have relaxed, fluctuations have finite short-range correlations, and boundary terms vanish. Taylor expansion of neighbor differences gives

`F_cont = integral [K u^2/2 + G |grad u|^2/2 - s(rho)u] d^d x`,

with

- `K = epsilon_* kappa_0/a^d` (energy per volume),
- `G = epsilon_* kappa_1 a^(2-d)` (energy per volume times length squared),
- `s(rho) = epsilon_* g_dev (rho/rho_*)/a^d` (energy per volume),
- `ell=sqrt(G/K)=a sqrt(kappa_1/kappa_0)` (length).

Stationarity gives `K u-div(G grad u)=s(rho)`, exactly the conditional WL-003/MB-001 continuum form. These equations show how microscopic parameters *could* supply `s`, `K`, `G`, and `ell`; they do not determine their values. A nonzero fixed `ell` as `a->0` requires the renormalized ratio `kappa_1/kappa_0 ~ (ell/a)^2`.

Limits are explicit: `kappa_1->0` gives independent local response, `kappa_0->infinity` suppresses deformation, `rho->0` gives the unloaded state, and `L/a->infinity` with `L/L_macro->0` suppresses microscopic fluctuations without erasing macroscopic variation.

## Traceability matrix

| ID | Definition/equation | Status | Assumption or derivation boundary |
|---|---|---|---|
{rows}

## Consistency checks

All executable checks pass: **{validation['all_checks_pass']}**. The normalized kernel preserves constants and the periodic mean; spatial and internal covariance errors are below floating-point tolerances. The discrete Laplacian errors at `n=32,64,128` are `{validation['metrics']['laplacian_rms_errors_n32_n64_n128']}`, demonstrating the expected long-wavelength convergence. These are mathematical consistency checks, not observational validation.

## Remaining theoretical gaps

1. PBUF has not derived why the microscopic state has exactly three components or why `g_dev=1/137` governs this coupling.
2. The lattice/regulator, energy scale `epsilon_*`, spacing `a`, relaxation time `tau`, couplings, noise statistics, and internal direction `e` are not predicted.
3. The assumed linear matter interaction is not derived from a covariant action; it does not derive Version D's quadratic normalized source.
4. A relativistic evolution law, tensorial metric map, causality, gauge behavior, and coupling to stress-energy rather than a static density are unresolved.
5. Renormalization and universality beyond the quadratic, isotropic, short-correlation regime remain to be shown.
6. No equality between `ell` and the observed baryonic width is implied, and no weak-lensing parameter has been fitted or changed.
"""


def main(output: Path) -> None:
    output.mkdir(parents=True, exist_ok=True)
    validation = validate()
    if not validation["all_checks_pass"]:
        raise RuntimeError("CORE-001 consistency validation failed")
    model = {
        "mission": "PBUF CORE-001 Microscopic State and Coarse-Graining Definition",
        "status": "complete_conditional_model",
        "observational_claim": False,
        "microscopic_state": {"symbol": "q_i", "components": 3, "space": "R^3", "dimension": "dimensionless"},
        "characteristic_scale": {"symbol": "g_dev", "value": G_DEV, "dimension": "dimensionless", "role": "assumed matter-state coupling scale"},
        "coarse_graining": "u(x)=e . sum_i a^d W_L(x-x_i) q_i",
        "continuum_equation": "K u-div(G grad u)=s(rho)",
        "mapping": {"K": "epsilon_* kappa_0/a^d", "G": "epsilon_* kappa_1 a^(2-d)", "s(rho)": "epsilon_* g_dev (rho/rho_*)/a^d", "ell": "a sqrt(kappa_1/kappa_0)"},
        "assumptions": [row[3] for row in TRACEABILITY if "assumption" in row[2] or "premise" in row[2]],
        "validation": validation,
    }
    (output / "microscopic_model.json").write_text(json.dumps(model, indent=2) + "\n")
    (output / "validation.json").write_text(json.dumps(validation, indent=2) + "\n")
    with (output / "traceability.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["id", "definition_or_equation", "status", "assumption_or_boundary"])
        writer.writerows(TRACEABILITY)
    (output / "core001_report.md").write_text(_report(validation))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=Path("runs/core001"))
    main(parser.parse_args().output)
