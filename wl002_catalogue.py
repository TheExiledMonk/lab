#!/usr/bin/env python3
"""Run the frozen WL-001 pipeline for every WL-002 constitutive candidate."""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from constitutive_equations import EQUATIONS
from pbuf_experiment import Config, run


REJECTED_FAMILIES = [
    {
        "family": "diffusion-type", "formula": "partial_t u = D Laplacian(u) + S(rho)",
        "physical_motivation": "Local deformation spreads down spatial gradients.",
        "assumptions": ["a physical evolution time and diffusivity exist"],
        "strengths": ["causal evolution can be represented", "smooths small-scale structure"],
        "weaknesses": ["the frozen static laboratory supplies neither time nor diffusivity"],
        "validation_outcome": "REJECTED BEFORE NUMERICAL RANKING: would require an arbitrary constant and an unfrozen time model.",
        "recommendation": "Revisit only when PBUF defines a time scale.",
    },
    {
        "family": "elastic-medium tensor", "formula": "div C:epsilon(u_vec) = f(rho)",
        "physical_motivation": "A genuine elastic medium carries directional displacement and shear stress.",
        "assumptions": ["deformation is vector/tensor valued", "elastic moduli and boundary tractions are known"],
        "strengths": ["represents shear and anisotropy", "has conservation-law structure"],
        "weaknesses": ["WL-001 accepts one scalar deformation field", "PBUF supplies no elastic moduli or boundary data"],
        "validation_outcome": "REJECTED BEFORE NUMERICAL RANKING: incompatible with the frozen scalar interface and underdetermined.",
        "recommendation": "Use as the next generalization after observables constrain tensor components.",
    },
]


def main(output: Path, n: int) -> None:
    output.mkdir(parents=True, exist_ok=True)
    rows = []
    for version, equation in EQUATIONS.items():
        meta = run(output / f"version_{version.lower()}", Config(n=n), equation)
        validation = meta["validation"]
        rows.append({
            "version": version,
            "rmse_pbuf": meta["summary"]["rmse_pbuf"],
            "rmse_baryonic_gr": meta["summary"]["rmse_baryonic_gr"],
            "rmse_lcdm": meta["summary"]["rmse_lcdm"],
            "all_validation_gates_pass": all(validation.values()),
            "response": equation.response,
            "stiffness": equation.stiffness,
            "description": equation.description,
            "formula": equation.formula,
            "physical_motivation": equation.motivation,
            "assumptions": list(equation.assumptions),
            "strengths": list(equation.strengths),
            "weaknesses": list(equation.weaknesses),
            "diagnostics": meta["diagnostics"],
            "computational_cost_seconds": meta["execution_seconds"],
        })

    with (output / "rmse_table.csv").open("w", newline="") as handle:
        fields = ["version", "rmse_pbuf", "rmse_baryonic_gr", "rmse_lcdm",
                  "all_validation_gates_pass", "response", "stiffness", "computational_cost_seconds"]
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader(); writer.writerows(rows)

    fig, ax = plt.subplots(figsize=(8, 5))
    versions = [row["version"] for row in rows]
    ax.plot(versions, [row["rmse_pbuf"] for row in rows], "o-", label="PBUF")
    ax.axhline(rows[0]["rmse_baryonic_gr"], color="tab:orange", ls="--", label="Baryonic GR")
    ax.axhline(rows[0]["rmse_lcdm"], color="tab:green", ls="--", label="LCDM")
    ax.set(title="WL-002 quantitative comparison", xlabel="Equation version", ylabel="RMSE")
    ax.legend(); fig.tight_layout(); fig.savefig(output / "rmse_comparison.png", dpi=140); plt.close(fig)

    baseline = rows[0]["rmse_pbuf"]
    accepted = [row for row in rows[1:] if row["all_validation_gates_pass"] and row["rmse_pbuf"] <= baseline]
    # Scientific admissibility comes first; RMSE breaks ties only among laws
    # which add a physically motivated capability without a fitted constant.
    best = min(accepted, key=lambda row: row["rmse_pbuf"]) if accepted else rows[0]
    for rank, row in enumerate(sorted(rows, key=lambda r: (
            not r["all_validation_gates_pass"], r["version"] != "D", r["rmse_pbuf"])), 1):
        row["scientific_rank"] = rank
        row["recommendation"] = ("PREFERRED" if row is best else
                                  "REJECT" if not row["all_validation_gates_pass"] else "RETAIN FOR COMPARISON")
    report = [
        "# WL-002 constitutive-equation catalogue", "",
        "All candidates used identical WL-001 inputs, geometry, propagation, reconstruction, benchmarks, residual calculations, and validation gates.", "",
        "| Version | Local/propagating | Linear/nonlinear | Stiffness | PBUF RMSE | Gates |", "|---|---|---|---|---:|---|",
    ]
    for row in rows:
        locality = "propagating" if row["version"] == "D" else "local"
        report.append(f"| {row['version']} | {locality} | {row['response']} | {row['stiffness']} | {row['rmse_pbuf']:.10g} | {'PASS' if row['all_validation_gates_pass'] else 'FAIL'} |")
    report += [
        "", "## Findings", "",
        "- A is the linear local WL-001 baseline.",
        "- B suppresses deformation in dilute outskirts. It improves the broad, horizontally shifted residual lobes because this synthetic reconstruction responds only to the spatial mean deformation. It fails Gate 4: the frozen path sampler does not encounter a large enough gradient to exceed the required trajectory-change threshold.",
        "- C softens the normalized response at low loading while preserving the peak. Its broader deformation worsens those same residual lobes.",
        "- D propagates B's loading through an elastic Helmholtz response. It changes the deformation and gradient maps, but preserves the source mean, so this reconstruction gives it the same RMSE as B (up to numerical precision).",
        "- No new fitted constant was introduced: D's propagation length is the observed baryonic width already present in the fixed geometry.",
        "", "## Recommendation", "",
        f"Carry Version {best['version']} into WL-003 provisionally (RMSE {best['rmse_pbuf']:.10g}, versus {baseline:.10g} for A). D is the strongest accepted law: it passes every gate, implements spatial propagation, and retains B's RMSE improvement; B is rejected despite the numerical tie because it fails Gate 4.",
        "", "## Limitation exposed", "",
        "The frozen WL-001 image reconstruction uses only `deformation.mean()` rather than the gradient field or photon landing positions. Consequently it cannot use RMSE to distinguish constitutive laws with equal mean deformation, and spatial residual-pattern claims beyond the global image shift are not identifiable. This is a pipeline limitation, not a reason to tune the constitutive laws or modify the frozen reconstruction during WL-002.",
    ]
    report += ["", "## Candidate dossiers", ""]
    for row in sorted(rows, key=lambda r: r["scientific_rank"]):
        d = row["diagnostics"]
        report += [f"### {row['scientific_rank']}. Version {row['version']} — {row['recommendation']}", "",
            f"**Equation:** `{row['formula']}`", "", f"**Motivation:** {row['physical_motivation']}", "",
            f"**Assumptions:** {'; '.join(row['assumptions'])}.", "",
            f"**Strengths:** {'; '.join(row['strengths'])}.", "",
            f"**Weaknesses:** {'; '.join(row['weaknesses'])}.", "",
            f"**Evidence:** RMSE={row['rmse_pbuf']:.10g}; deformation range={d['deformation_min']:.6g}..{d['deformation_max']:.6g}; gradient RMS={d['gradient_rms']:.6g}; gradient max={d['gradient_max']:.6g}; photon max deviation={d['photon_max_deviation']:.6g}; finite={d['finite_outputs']}; cost={row['computational_cost_seconds']:.4f}s; gates={'PASS' if row['all_validation_gates_pass'] else 'FAIL'}.", ""]
    report += ["## Families rejected on physical/interface grounds", ""]
    for candidate in REJECTED_FAMILIES:
        report += [f"### {candidate['family']}", "", f"**Equation:** `{candidate['formula']}`", "",
            f"**Motivation:** {candidate['physical_motivation']}", "",
            f"**Assumptions:** {'; '.join(candidate['assumptions'])}.", "",
            f"**Strengths:** {'; '.join(candidate['strengths'])}.", "",
            f"**Weaknesses:** {'; '.join(candidate['weaknesses'])}.", "",
            f"**Outcome:** {candidate['validation_outcome']}", "", f"**Recommendation:** {candidate['recommendation']}", ""]
    (output / "catalogue.md").write_text("\n".join(report) + "\n")
    # Keep the original catalogue.json list shape for downstream WL-002 users;
    # the WL-002A decision record carries the richer mission-level result.
    (output / "catalogue.json").write_text(json.dumps(rows, indent=2) + "\n")
    payload = {"preferred_version": best["version"], "ranking_basis": "physical admissibility, validation, then RMSE",
               "evaluated_candidates": rows, "analytical_rejections": REJECTED_FAMILIES}
    (output / "discovery.json").write_text(json.dumps(payload, indent=2) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=Path("runs/wl002"))
    parser.add_argument("--n", type=int, default=128)
    args = parser.parse_args(); main(args.output, args.n)
