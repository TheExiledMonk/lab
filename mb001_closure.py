#!/usr/bin/env python3
"""MB-001 micro--macro closure audit and frozen-laboratory validation.

The supplied PBUF sources name microscopic quantities but provide no equations
or dimensions that map them to the continuum source and moduli.  MB-001 must
therefore report Outcome C rather than turn Version D's empirical choices into
purported derivations.  The unchanged candidate is rerun only as a frozen-
laboratory reproducibility control.
"""
from __future__ import annotations

import argparse
import csv
import json
from dataclasses import fields
from pathlib import Path

import numpy as np

from constitutive_equations import get_equation
from pbuf_experiment import Config, grid, propagate, run


EQUATIONS = [
    {
        "equation_id": "MB-001-E01",
        "equation": "K u - div(G grad u) = s(rho)",
        "status": "conditional continuum balance; not a closed PBUF law",
        "originating_pbuf_equations": "none supplied; WL-003 conditional balance",
        "physical_interpretation": "local recovery plus divergence of elastic flux balances matter loading",
        "units": "[K]=E L^-d u^-2; [G]=[K] L^2; [s]=[K]u (u dimensionless in WL)",
        "assumptions": "static scalar isotropic medium; differentiable field; positive K and G",
        "limiting_cases": "G->0 gives local response Ku=s; s->0 gives homogeneous recovery",
        "closure_gap": "PBUF supplies no map from microscopic variables to s, K, or G",
    },
    {
        "equation_id": "MB-001-E02",
        "equation": "J = -G grad u; K u + div J = s(rho)",
        "status": "conservation form of MB-001-E01",
        "originating_pbuf_equations": "algebraic rewriting of MB-001-E01",
        "physical_interpretation": "J is deformation-transmission flux",
        "units": "[J]=[G]u/L; [div J]=[K]u",
        "assumptions": "same as MB-001-E01",
        "limiting_cases": "uniform u has J=0; G->0 removes propagation",
        "closure_gap": "the microscopic transported quantity and coefficient G are undefined",
    },
    {
        "equation_id": "MB-001-E03",
        "equation": "ell^2 = G/K",
        "status": "derived algebraically, conditional on MB-001-E01",
        "originating_pbuf_equations": "MB-001-E01 with constant positive K and G",
        "physical_interpretation": "recovery/propagation balance length",
        "units": "[G]/[K]=L^2, hence [ell]=L",
        "assumptions": "K>0 and G>0 are constant scalar coefficients",
        "limiting_cases": "G->0 gives ell->0; K->infinity gives ell->0",
        "closure_gap": "neither K nor G is supplied by PBUF, so ell has no predicted value",
    },
    {
        "equation_id": "MB-001-E04",
        "equation": "s(rho)/K = u0 (rho/rho_max)^2",
        "status": "empirical Version-D identification; not derived",
        "originating_pbuf_equations": "WL-002 Version D only",
        "physical_interpretation": "normalized quadratic matter loading",
        "units": "both sides have units of u; rho/rho_max is dimensionless",
        "assumptions": "rho_max>0; quadratic exponent; normalization and u0 are inherited",
        "limiting_cases": "rho=0 gives zero loading; rho=rho_max gives u0",
        "closure_gap": "microscopic matter-medium interaction/action or response law is absent",
    },
    {
        "equation_id": "MB-001-E05",
        "equation": "ell = sigma_rho",
        "status": "empirical Version-D identification; not derived",
        "originating_pbuf_equations": "WL-002 Version D only",
        "physical_interpretation": "uses observed lens width as propagation length",
        "units": "both quantities are lengths",
        "assumptions": "medium correlation length equals source geometry",
        "limiting_cases": "point-source limit forces ell->0, exposing the unproved identification",
        "closure_gap": "PBUF gives no correlation-length or dispersion relation",
    },
]


def _load(path: Path) -> np.ndarray:
    return np.loadtxt(path, delimiter=",")


def _metric(old: np.ndarray, new: np.ndarray) -> dict:
    delta = new - old
    old_flat, new_flat = old.ravel(), new.ravel()
    if np.std(old_flat) == 0 and np.std(new_flat) == 0:
        correlation = 1.0 if np.array_equal(old, new) else None
    else:
        correlation = float(np.corrcoef(old_flat, new_flat)[0, 1])
    return {
        "max_abs_delta": float(np.max(np.abs(delta))),
        "rmse_delta": float(np.sqrt(np.mean(delta * delta))),
        "topology_correlation": correlation,
        "checksum_equal": bool(np.array_equal(old, new)),
    }


def _config(meta: dict, requested_n: int | None) -> Config:
    values = meta["config"].copy()
    if requested_n is not None:
        values["n"] = requested_n
    allowed = {item.name for item in fields(Config)}
    return Config(**{key: value for key, value in values.items() if key in allowed})


def compare(archive: Path, rerun: Path, old_meta: dict, new_meta: dict, config: Config) -> dict:
    names = (
        "deformation", "gradient_x", "gradient_y", "pbuf", "residual",
        "observation_minus_pbuf", "pbuf_minus_lcdm",
    )
    arrays = {name: _metric(_load(archive / f"{name}.csv"), _load(rerun / f"{name}.csv")) for name in names}
    X, Y = grid(config)
    old_paths = propagate(_load(archive / "deformation.csv"), X, Y, config)
    new_paths = propagate(_load(rerun / "deformation.csv"), X, Y, config)
    old_xy = np.stack([np.stack(pair) for pair in old_paths])
    new_xy = np.stack([np.stack(pair) for pair in new_paths])
    return {
        "archive": str(archive),
        "same_grid_as_archive": config.n == old_meta["config"]["n"],
        "equation_id_unchanged": old_meta["equation_id"] == new_meta["equation_id"],
        "field_and_residual_comparison": arrays,
        "photon_trajectories": _metric(old_xy, new_xy),
        "diagnostics": {
            key: {"archived": old_meta["diagnostics"][key], "rerun": new_meta["diagnostics"][key],
                  "change": new_meta["diagnostics"][key] - old_meta["diagnostics"][key]}
            for key in ("deformation_mean", "deformation_std", "gradient_rms", "gradient_max", "photon_max_deviation")
        },
        "rmse": {
            key: {"archived": old_meta["summary"][key], "rerun": new_meta["summary"][key],
                  "change": new_meta["summary"][key] - old_meta["summary"][key]}
            for key in ("rmse_pbuf", "rmse_baryonic_gr", "rmse_lcdm")
        },
        "validation": new_meta["validation"],
        "all_validation_gates_pass": all(new_meta["validation"].values()),
    }


def report(result: dict) -> str:
    validation = result["validation"]
    comparisons = validation["field_and_residual_comparison"]
    rows = "\n".join(
        f"| {e['equation_id']} | `{e['equation']}` | {e['status']} | {e['originating_pbuf_equations']} | {e['closure_gap']} |"
        for e in EQUATIONS
    )
    evidence = "\n".join(
        f"| {name} | {item['max_abs_delta']:.3g} | {item['rmse_delta']:.3g} | {item['topology_correlation']} |"
        for name, item in comparisons.items()
    )
    return f"""# PBUF MB-001 micro--macro closure report

## Result: Outcome C

Existing supplied PBUF theory is insufficient to close the frozen macroscopic law. No supplied document defines a microscopic deformation variable, a quantitative matter--medium interaction, a free energy, a correlation function, or a dispersion relation from which `s(rho)`, `K`, or `G` can be calculated. The named V11 quantity `alpha_T(a)`, together with `epsilon_0,T(a)`, `k_max(a)`, and thermal rigidity, is explicitly only a candidate input in the Physics Starting Pack and is supplied there without equations, definitions, numerical values, or dimensions. This reference does not identify `alpha_T(a)` with any microscopic matter--state vertex.

Accordingly, no constitutive code was changed. Version D remains empirical, and its unchanged rerun is a reproducibility control—not validation of a new closure.

## What can be established conditionally

For the scalar laboratory variable `u`, static local balance may be written as MB-001-E01. Defining the flux in MB-001-E02 gives conservation form. If `K` and `G` are positive constants, division by `K` gives `(1 - ell^2 Laplacian)u = s/K` and MB-001-E03 follows. These statements determine the form and dimensions of a possible closure but do not derive its coefficients from PBUF.

## Research-task disposition

1. **Microscopic quantity:** absent. The laboratory's dimensionless scalar `u` is only a proxy; no supplied PBUF equation maps a microstate, strain, occupancy, or metric perturbation to it.
2. **Matter loading:** absent. “Mass loads spacetime” is qualitative; no action, conjugate force, susceptibility, or response function derives `s(rho)`. MB-001-E04 remains empirical.
3. **Effective stiffness:** absent. Thermal rigidity is named but never quantitatively related to `K`, and its dimensions are not supplied.
4. **Propagation coefficient:** absent. No gradient-energy coefficient, microscopic coupling, correlation function, or dispersion relation derives `G`.
5. **Propagation length:** only the conditional ratio `ell=sqrt(G/K)` emerges. Its value does not; MB-001-E05 remains empirical.
6. **Conservation form:** MB-001-E02 supplies the conditional continuum form.

## Equation-to-PBUF traceability matrix

| ID | Equation | Status | Origin | Unclosed relationship |
|---|---|---|---|---|
{rows}

Full physical interpretations, dimensional justifications, assumptions, and limiting cases are in `equation_traceability.csv` and `closure_equations.json`.

## Frozen-laboratory validation against archived Version D

The rerun used the archived configuration: `{validation['same_grid_as_archive']}`. Equation ID unchanged: `{validation['equation_id_unchanged']}`. All validation gates pass: `{validation['all_validation_gates_pass']}`.

| Artifact | Maximum absolute change | Delta RMSE | Topology correlation |
|---|---:|---:|---:|
{evidence}
| photon trajectories | {validation['photon_trajectories']['max_abs_delta']:.3g} | {validation['photon_trajectories']['rmse_delta']:.3g} | {validation['photon_trajectories']['topology_correlation']} |

PBUF RMSE was `{validation['rmse']['rmse_pbuf']['archived']:.12g}` archived and `{validation['rmse']['rmse_pbuf']['rerun']:.12g}` on rerun (change `{validation['rmse']['rmse_pbuf']['change']:.3g}`). Detailed deformation, gradient, trajectory, topology, residual, RMSE, stability, and gate comparisons are in `validation.json`.

## Exact missing physical law and next milestone

PBUF must supply a dimensionally explicit microscopic energy/action or response functional `F[microstate, rho]` together with a coarse-graining definition `u=C[microstate]`. Its long-wavelength expansion must independently predict the matter-conjugate source `s(rho)`, local curvature `K`, and gradient coefficient `G`. This would make `ell=sqrt(G/K)` a prediction and test—rather than assume—the Version-D relations `s/K=u0(rho/rho_max)^2` and `ell=sigma_rho`.

Until that law is supplied, retain Version D only as the leading empirical candidate and do not update the constitutive law.
"""


def main(output: Path, archive: Path, requested_n: int | None) -> None:
    meta_path = archive / "run.json"
    if not meta_path.exists():
        raise FileNotFoundError(f"archived Version D not found at {archive}")
    old_meta = json.loads(meta_path.read_text())
    config = _config(old_meta, requested_n)
    if config.n != old_meta["config"]["n"]:
        raise ValueError("comparison requires the archived grid size; omit --n or use the archived value")
    output.mkdir(parents=True, exist_ok=True)
    rerun = output / "version_d_reproducibility"
    new_meta = run(rerun, config, get_equation("D"))
    validation = compare(archive, rerun, old_meta, new_meta, config)
    result = {
        "mission": "PBUF MB-001 Micro-Macro Closure Development",
        "outcome": "C",
        "constitutive_law_updated": False,
        "reason": "No supplied PBUF equation quantitatively maps microscopic quantities to u, s(rho), K, or G.",
        "conditional_closure": "K u - div(G grad u) = s(rho); ell^2=G/K",
        "missing_physical_law": "A microscopic energy/action or response plus coarse graining that predicts s(rho), K, and G.",
        "empirical_version_d_relations": ["s(rho)/K=u0(rho/rho_max)^2", "ell=sigma_rho"],
        "equations": EQUATIONS,
        "validation": validation,
        "recommendation": "Define and dimensionally validate the microscopic free energy/response and coarse-graining map before changing Version D.",
    }
    (output / "closure_equations.json").write_text(json.dumps(result, indent=2) + "\n")
    with (output / "equation_traceability.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(EQUATIONS[0]))
        writer.writeheader()
        writer.writerows(EQUATIONS)
    (output / "validation.json").write_text(json.dumps(validation, indent=2) + "\n")
    (output / "micro_macro_closure_report.md").write_text(report(result))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=Path("runs/mb001"))
    parser.add_argument("--archive", type=Path, default=Path("runs/wl003/version_d_reproducibility"))
    parser.add_argument("--n", type=int, default=None, help="must match the archived grid size")
    args = parser.parse_args()
    main(args.output, args.archive, args.n)
