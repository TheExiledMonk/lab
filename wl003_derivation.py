#!/usr/bin/env python3
"""WL-003 physical-derivation audit and frozen-laboratory validation.

The supplied PBUF documents do not contain a quantitative micro--macro
closure.  This program therefore does not silently promote Version D to a
derived law.  It records the conditional conservation-law derivation, names
the missing closure, reruns the unchanged leading candidate, and compares the
rerun with its archived WL-002 result.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np

from constitutive_equations import get_equation
from pbuf_experiment import Config, run


SYMBOLS = [
    ("u", "scalar spacetime-deformation proxy", "dimensionless", "field passed to the frozen photon-gradient interface", "WL-001 Physics Starting Pack: local spacetime deformation field", "scalar isotropic reduction"),
    ("rho", "identified baryonic matter density/map", "input-map units", "matter loading of the medium", "WL-001 Physics Starting Pack: mass loads spacetime", "the supplied 2-D map is an adequate loading variable"),
    ("q=rho/rho_max", "normalized matter loading", "dimensionless", "removes the input map's amplitude unit", "Version A/WL-001 implemented normalization", "rho_max > 0; zero input handled separately"),
    ("u0", "deformation-strength scale", "dimensionless", "peak local response scale", "fixed WL-001 PBUF constant (deformation_strength)", "its numerical value is inherited, not derived in WL-003"),
    ("s(rho)", "matter-produced generalized source", "K (energy-density per u)", "work conjugate to u", "PBUF matter--medium interaction concept", "no quantitative PBUF source law was supplied"),
    ("K", "local recovery stiffness", "energy density per u^2", "penalizes local deformation", "PBUF thermal-rigidity concept", "positive and spatially constant for Version D"),
    ("G", "gradient stiffness", "K times length^2", "penalizes spatial variation and transmits deformation", "existing elasticity concept requested by WL-003", "positive, isotropic and spatially constant"),
    ("ell=sqrt(G/K)", "propagation length", "length", "balance scale of gradient and local recovery", "follows algebraically from the elasticity balance", "G and K must be supplied by a PBUF closure"),
    ("sigma_rho", "observed baryonic Gaussian width", "length", "Version D's current proxy for ell", "WL-002 Version D", "observed geometry scale; equality ell=sigma_rho is not derived"),
    ("nabla", "spatial gradient", "inverse length", "measures spatial variation", "standard continuum conservation notation", "flat 2-D frozen laboratory coordinates"),
]


def _load_csv(path: Path) -> np.ndarray:
    return np.loadtxt(path, delimiter=",")


def _comparison(archive: Path, rerun: Path, new_meta: dict) -> dict:
    old_meta = json.loads((archive / "run.json").read_text())
    field_metrics = {}
    for name in ("deformation", "gradient_x", "gradient_y"):
        old, new = _load_csv(archive / f"{name}.csv"), _load_csv(rerun / f"{name}.csv")
        delta = new - old
        field_metrics[name] = {
            "max_abs_delta": float(np.max(np.abs(delta))),
            "rmse_delta": float(np.sqrt(np.mean(delta * delta))),
            "topology_correlation": float(np.corrcoef(old.ravel(), new.ravel())[0, 1]),
        }
    old_d, new_d = old_meta["diagnostics"], new_meta["diagnostics"]
    return {
        "archived_equation_id": old_meta["equation_id"],
        "rerun_equation_id": new_meta["equation_id"],
        "rmse_pbuf_archived": old_meta["summary"]["rmse_pbuf"],
        "rmse_pbuf_rerun": new_meta["summary"]["rmse_pbuf"],
        "rmse_pbuf_change": new_meta["summary"]["rmse_pbuf"] - old_meta["summary"]["rmse_pbuf"],
        "field_comparison": field_metrics,
        "photon_max_deviation_archived": old_d["photon_max_deviation"],
        "photon_max_deviation_rerun": new_d["photon_max_deviation"],
        "photon_max_deviation_change": new_d["photon_max_deviation"] - old_d["photon_max_deviation"],
        "finite_outputs": new_d["finite_outputs"],
        "all_validation_gates_pass": all(new_meta["validation"].values()),
    }


def _report(comparison: dict) -> str:
    f = comparison["field_comparison"]
    rows = "\n".join(f"| {s} | {d} | {u} | {m} | {o} | {a} |" for s, d, u, m, o, a in SYMBOLS)
    return f"""# PBUF WL-003 physical derivation report

## Result: Outcome B

Version D can be placed in conservation form, but it cannot be fully derived from the supplied PBUF theory. The exact missing theoretical component is a **quantitative micro--macro constitutive closure** mapping existing PBUF quantities (including thermal rigidity, `alpha_T`, `epsilon_0,T`, and `k_max`) to the matter source `s(rho)` and continuum coefficients `K` and `G`.

No replacement equation has therefore been introduced. Version D remains explicitly empirical pending that closure.

## Conditional derivation

Let `u` be the scalar deformation proxy accepted by the frozen laboratory. The least continuum free-energy density containing local recovery, elastic transmission, and matter work is

`psi(u, grad u; rho) = (K/2) u^2 + (G/2) |grad u|^2 - s(rho) u`.

Stationarity of its spatial integral gives

`K u - div(G grad u) = s(rho)`.

This is conservation form: the elastic flux is `J = -G grad u`, and local equilibrium is `K u + div J = s`. Under the additional Version-D assumptions that `K` and `G` are positive constants,

`(1 - ell^2 Laplacian) u = s(rho)/K`, where `ell^2 = G/K`.

Version D is recovered only by separately postulating

`s(rho)/K = u0 (rho/rho_max)^2` and `ell = sigma_rho`.

The documents qualitatively motivate loading, propagation, and rigidity. They do not derive the quadratic exponent, give `K` or `G` in terms of PBUF quantities, or prove `sqrt(G/K) = sigma_rho`. Consequently the displayed steps are a conditional continuum derivation, not a completed PBUF derivation.

## Symbol traceability

| Symbol | Definition | Units | Physical meaning | Origin | Assumption/status |
|---|---|---|---|---|---|
{rows}

## Research-question answers

1. The executable deforms a dimensionless scalar proxy `u`; the supplied theory does not establish its microscopic identity or a tensorial parent field.
2. Matter produces a generalized source `s(rho)`, but PBUF supplies no quantitative source law; Version D postulates a normalized quadratic law.
3. Deformation propagates when spatial variation costs elastic energy (`G > 0`), producing the flux `-G grad u`.
4. A Helmholtz operator follows from balancing that gradient flux against local recovery `K u` in static equilibrium.
5. The propagation length is `sqrt(G/K)`. It is derived from continuum coefficients, but neither coefficient nor its equality to the observed mass width is derived from supplied PBUF quantities.
6. Thermal rigidity is a plausible origin for `K`, but the needed units and mapping equation are absent.
7. Yes: `K u - div(G grad u) = s(rho)`.

## Frozen-laboratory reproducibility validation

Because no derived replacement exists, the correct non-tuning check is a rerun of unchanged Version D against the archived Version D:

- PBUF RMSE: `{comparison['rmse_pbuf_archived']:.12g}` -> `{comparison['rmse_pbuf_rerun']:.12g}` (change `{comparison['rmse_pbuf_change']:.3g}`).
- Deformation topology correlation: `{f['deformation']['topology_correlation']:.12g}`; max absolute field change `{f['deformation']['max_abs_delta']:.3g}`.
- Gradient topology correlations: x `{f['gradient_x']['topology_correlation']:.12g}`, y `{f['gradient_y']['topology_correlation']:.12g}`.
- Photon maximum deviation: `{comparison['photon_max_deviation_archived']:.12g}` -> `{comparison['photon_max_deviation_rerun']:.12g}` (change `{comparison['photon_max_deviation_change']:.3g}`).
- Stability: finite outputs = `{comparison['finite_outputs']}`; all validation gates pass = `{comparison['all_validation_gates_pass']}`.

This is a reproducibility comparison, not evidence that the empirical closures are physically derived.

## Recommendation

Retain Version D as the leading empirical candidate. Do not label it derived or replace it until PBUF provides one dimensionally explicit closure giving `s(rho)`, `K`, and `G` from its existing microscopic/thermal quantities. That single closure would test both disputed identifications: `s/K = u0 q^2` and `G/K = sigma_rho^2`.
"""


def main(output: Path, archive: Path, n: int) -> None:
    if not (archive / "run.json").exists():
        raise FileNotFoundError(f"archived Version D not found at {archive}")
    output.mkdir(parents=True, exist_ok=True)
    rerun = output / "version_d_reproducibility"
    meta = run(rerun, Config(n=n), get_equation("D"))
    comparison = _comparison(archive, rerun, meta)
    result = {
        "outcome": "B",
        "equation_replaced": False,
        "recommendation": "Version D remains empirical pending a quantitative PBUF micro--macro closure.",
        "missing_relationship": "A dimensionally explicit map from PBUF microscopic/thermal quantities to s(rho), K, and G.",
        "conditional_equation": "K u - div(G grad u) = s(rho); ell^2 = G/K",
        "version_d_empirical_identifications": ["s(rho)/K = u0 (rho/rho_max)^2", "G/K = sigma_rho^2"],
        "comparison": comparison,
    }
    (output / "derivation.json").write_text(json.dumps(result, indent=2) + "\n")
    with (output / "symbol_traceability.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["symbol", "definition", "units", "physical_meaning", "origin", "assumption_status"])
        writer.writerows(SYMBOLS)
    (output / "derivation_report.md").write_text(_report(comparison))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=Path("runs/wl003"))
    parser.add_argument("--archive", type=Path, default=Path("runs/wl002/version_d"))
    parser.add_argument("--n", type=int, default=128)
    args = parser.parse_args()
    main(args.output, args.archive, args.n)
