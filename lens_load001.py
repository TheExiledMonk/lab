#!/usr/bin/env python3
"""PBUF LENS-LOAD-001: local single-lens reconstruction readiness gate.

The frozen corpus does not contain a selected medium-to-metric map.  This
program therefore never turns lensing pixels into a native load by assumption.
It writes the exact conditional mechanical inverse and, when a placement is
provided, evaluates the frozen weak-field Candidate-A load reconstruction.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent
DEFAULT_OUT = ROOT / "runs/lens_load001"
FROZEN = {
    "CONSTITUTIVE-CONSTRUCTION-001": "runs/constitutive_construction001/constitutive_construction_report.md",
    "WEAK-LENSING-LOCALITY-001": "runs/weak_lensing_locality001/weak_lensing_locality_report.md",
    "LOCAL-STATE-001": "runs/local_state001/local_state_report.md",
    "INVERSE-SOURCE-001": "runs/inverse_source001/inverse_source_report.md",
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_field(path: Path) -> np.ndarray:
    a = np.loadtxt(path, delimiter=",")
    if a.ndim != 2 or min(a.shape) < 3 or not np.isfinite(a).all():
        raise ValueError(f"{path}: expected a finite 2-D field of at least 3x3")
    return a


def derivatives(a: np.ndarray, dx: float, dy: float) -> tuple[np.ndarray, np.ndarray]:
    ay, ax = np.gradient(a, dy, dx, edge_order=2)
    return ax, ay


def conditional_scalar_load(u: np.ndarray, dx: float, dy: float, modulus: float) -> np.ndarray:
    """Scalar anti-plane restriction of -Div(P_F), for verification only.

    In this restriction E_xz=u_x/2, E_yz=u_y/2 and Candidate A gives a
    shear response proportional to mu0.  It is not an optical map and it is
    used only when u has independently been supplied.
    """
    load = np.zeros_like(u)
    load[1:-1, 1:-1] = modulus * (
        (2*u[1:-1, 1:-1]-u[1:-1, :-2]-u[1:-1, 2:]) / dx**2
        + (2*u[1:-1, 1:-1]-u[:-2, 1:-1]-u[2:, 1:-1]) / dy**2
    )
    return load


def forward_scalar(load: np.ndarray, dx: float, dy: float, modulus: float) -> np.ndarray:
    """Solve -mu Laplacian(u)=b with zero Dirichlet boundary by Jacobi CG."""
    ny, nx = load.shape
    iy, ix = ny - 2, nx - 2
    n = iy * ix
    if n <= 0:
        raise ValueError("grid has no interior")

    def apply(v: np.ndarray) -> np.ndarray:
        z = np.zeros((ny, nx)); z[1:-1, 1:-1] = v.reshape(iy, ix)
        az = modulus * ((2*z[1:-1, 1:-1]-z[1:-1, :-2]-z[1:-1, 2:]) / dx**2
                        + (2*z[1:-1, 1:-1]-z[:-2, 1:-1]-z[2:, 1:-1]) / dy**2)
        return az.ravel()

    rhs = load[1:-1, 1:-1].ravel()
    x = np.zeros(n); r = rhs.copy(); p = r.copy(); rr = float(r @ r)
    tol = max(1e-24, rr * 1e-24)
    for _ in range(max(100, 4*n)):
        ap = apply(p); denom = float(p @ ap)
        if denom <= 0 or rr <= tol:
            break
        alpha = rr / denom; x += alpha*p; r -= alpha*ap
        rr_new = float(r @ r)
        if rr_new <= tol:
            rr = rr_new; break
        p = r + (rr_new/rr)*p; rr = rr_new
    u = np.zeros((ny, nx)); u[1:-1, 1:-1] = x.reshape(iy, ix)
    return u


def write_csv(path: Path, a: np.ndarray) -> None:
    np.savetxt(path, a, delimiter=",")


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--output", type=Path, default=DEFAULT_OUT)
    p.add_argument("--matter", type=Path, default=ROOT / "runs/wl001/matter.csv")
    p.add_argument("--lensing", type=Path, default=ROOT / "runs/wl001/observation.csv")
    p.add_argument("--placement", type=Path, help="independently reconstructed scalar placement; not inferred from lensing")
    p.add_argument("--extent", type=float, default=8.0)
    p.add_argument("--mu0", type=float, default=1.0, help="frozen modulus in the supplied unit system")
    a = p.parse_args()
    if a.extent <= 0 or a.mu0 <= 0:
        p.error("extent and mu0 must be positive")
    missing = [ROOT / v for v in FROZEN.values() if not (ROOT / v).is_file()]
    if missing:
        raise FileNotFoundError(missing)
    matter, lensing = read_field(a.matter), read_field(a.lensing)
    if matter.shape != lensing.shape:
        raise ValueError("matter and lensing grids must have the same shape")
    out = a.output; out.mkdir(parents=True, exist_ok=True)
    ny, nx = matter.shape; dx = 2*a.extent/(nx-1); dy = 2*a.extent/(ny-1)

    # The archived optical paths and both supplied maps occupy this square.
    # A zero-placement boundary removes rigid modes, but is only a declared
    # finite isolated-boundary approximation, never an exact infinity claim.
    domain = {
        "lens_id": "Lens 001", "dimension": 2,
        "bounds": {"x": [-a.extent, a.extent], "y": [-a.extent, a.extent]},
        "grid": {"nx": nx, "ny": ny, "dx": dx, "dy": dy},
        "contains": ["supplied baryonic map", "supplied lensing map", "archived optical-path window"],
        "boundary": "zero placement on all four sides (declared isolated finite-truncation approximation)",
        "justification": "the data footprint is the smallest evidenced region; full Dirichlet data make the local elastic restriction coercive and remove rigid modes",
        "limitation": "the frozen theory supplies no falloff rate, so boundary adequacy needs a padding/convergence study with real data",
    }
    (out / "computational_domain.json").write_text(json.dumps(domain, indent=2)+"\n")
    write_csv(out / "baryonic_mass.csv", matter)
    write_csv(out / "observed_lensing.csv", lensing)

    status = {
        "milestone": "PBUF LENS-LOAD-001", "lens_id": "Lens 001",
        "end_to_end_reconstruction": "BLOCKED_BY_FROZEN_IDENTIFIABILITY",
        "load_dataset_status": "NOT_RECONSTRUCTED_FROM_LENSING",
        "missing_required_inputs": [
            "selected frozen-compatible medium-to-metric map G with normalization and support",
            "measurement/ray operator mapping G to the reported weak-lensing samples",
            "observational covariance/uncertainties (none supplied)",
            "induced elastic boundary work or an independently justified boundary model",
        ],
        "reason": "weak lensing does not determine placement y; without y, Pi_req=A_Omega(y) cannot be evaluated",
        "no_physical_interpretation": True,
        "no_new_law_or_metric": True,
        "inputs": {"matter": str(a.matter), "lensing": str(a.lensing)},
        "checksums": {"matter": sha256(a.matter), "lensing": sha256(a.lensing)},
        "frozen_sources": FROZEN,
    }

    rows = [["quantity", "result", "qualification"],
            ["spatial_distribution", "not identifiable", "placement-to-optics map unselected"],
            ["magnitude", "not identifiable", "metric-map normalization unselected"],
            ["symmetry", "not identifiable", "a regularizer cannot create information"],
            ["localization", "not identifiable", "metric-map support unselected"],
            ["uncertainty", "unbounded structural component", "no covariance supplied; optical null space remains"],
            ["mass_correlation", "not computable", "there is no reconstructed load field"]]
    with (out / "load_characterization.csv").open("w", newline="") as f:
        csv.writer(f).writerows(rows)

    if a.placement:
        u = read_field(a.placement)
        if u.shape != matter.shape:
            raise ValueError("placement grid must match matter grid")
        if not (np.allclose(u[[0, -1], :], 0) and np.allclose(u[:, [0, -1]], 0)):
            raise ValueError("placement violates the declared zero Dirichlet boundary")
        load = conditional_scalar_load(u, dx, dy, a.mu0)
        # Boundary finite differences are not a bulk load reconstruction.
        load[[0, -1], :] = np.nan; load[:, [0, -1]] = np.nan
        solve_load = np.nan_to_num(load)
        recovered = forward_scalar(solve_load, dx, dy, a.mu0)
        residual = recovered-u
        write_csv(out / "conditional_native_load.csv", load)
        write_csv(out / "conditional_forward_placement.csv", recovered)
        write_csv(out / "conditional_forward_residual.csv", residual)
        interior = residual[1:-1, 1:-1]
        status["conditional_mechanical_reconstruction"] = {
            "status": "PASS" if np.sqrt(np.mean(interior**2)) < 1e-7 else "FAIL",
            "source": str(a.placement), "not_from_lensing": True,
            "rmse": float(np.sqrt(np.mean(interior**2))),
            "max_abs": float(np.max(np.abs(interior))),
            "formula": "b_req=-Div_0 P_F[y] (scalar anti-plane Candidate-A restriction)",
        }
        fig, axes = plt.subplots(1, 3, figsize=(12, 4))
        for ax, field, title in zip(axes, (u, load, residual), ("Supplied placement", "Conditional load", "Forward residual")):
            im = ax.imshow(field, origin="lower", extent=[-a.extent, a.extent]*2)
            ax.set_title(title); fig.colorbar(im, ax=ax, shrink=.8)
        fig.tight_layout(); fig.savefig(out / "conditional_load_maps.png", dpi=150); plt.close(fig)
    (out / "reconstruction_status.json").write_text(json.dumps(status, indent=2)+"\n")
    return 0 if status.get("conditional_mechanical_reconstruction", {}).get("status") != "FAIL" else 1


if __name__ == "__main__":
    raise SystemExit(main())
