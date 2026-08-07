#!/usr/bin/env python3
"""Produce the theory-only PBUF FND-002 assumption audit.

The audit does not alter CORE-001 or the frozen weak-lensing laboratory.  It
records which pieces of the CORE construction are consequences of a small set
of continuum requirements and which are merely representatives of a
universality class.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


LABELS = {"DERIVED", "SUPPORTED", "OPTIONAL", "NEW POSTULATE"}

AUDIT = [
    {
        "id": "FND-002-A01", "assumption": "Exactly three microscopic degrees of freedom", "label": "NEW POSTULATE",
        "explicit_in_pbuf": "Yes: CORE-001 calls it an existing PBUF working premise, while explicitly saying it is not derived.",
        "derivable": "No; neither the supplied equations nor symmetry data select dim(q)=3.",
        "weaken": "Yes: require only a microscopic state with at least one matter-coupled projection.",
        "remove": "From the scalar closure, yes; from the stated three-degree PBUF ontology, no.",
        "observable": "Extra components matter only if they have distinct sources, modes, or couplings; the frozen laboratory observes one projection.",
        "changed": "CORE-001 internal O(3) interpretation changes, but its scalar long-wave equation need not.",
    },
    {
        "id": "FND-002-A02", "assumption": "Characteristic coupling scale g_dev=1/137", "label": "NEW POSTULATE",
        "explicit_in_pbuf": "Yes: CORE-001 identifies 1/137 as an existing characteristic-scale premise, not a derivation.",
        "derivable": "No. Corrected CORE-001 places g_dev directly in the normalized matter vertex; its numerical value remains a PBUF premise rather than a consequence of the other assumptions.",
        "weaken": "No within the stated PBUF ontology: replacing g_dev by an effective coefficient would discard the fundamental-coupling premise.",
        "remove": "No without removing the stated PBUF matter coupling; the zero-coupling limit removes matter loading.",
        "observable": "The normalized microscopic source and conditional coarse source scale directly with g_dev; downstream observability still requires response and access maps.",
        "changed": "The former source-normalization degeneracy is withdrawn; lack of a derivation of the numerical value remains.",
    },
    {
        "id": "FND-002-A03", "assumption": "Isotropic lattice/regulator", "label": "OPTIONAL",
        "explicit_in_pbuf": "No earlier supplied material requires a lattice; CORE-001 introduced it as a regulator.",
        "derivable": "No specific regulator is derivable. Isotropic leading-order propagation follows conditionally from spatial rotational symmetry.",
        "weaken": "Yes: any statistically homogeneous, isotropic local substrate with a valid long-wave limit suffices.",
        "remove": "Yes: a continuum microfield or isotropic graph Laplacian gives the same effective operator.",
        "observable": "Finite-spacing anisotropy and dispersion can appear beyond the long-wave regime.",
        "changed": "Only regulator-level corrections change if the same isotropic continuum limit is retained.",
    },
    {
        "id": "FND-002-A04", "assumption": "Linear matter-state interaction", "label": "OPTIONAL",
        "explicit_in_pbuf": "No quantitative matter-state law exists in the supplied pre-CORE material.",
        "derivable": "Only as the first Taylor term of an analytic interaction about a reference state, assuming a nonzero linear coefficient.",
        "weaken": "Yes: require a differentiable local source whose leading perturbative term is linear.",
        "remove": "Yes, but a purely nonlinear source changes the weak-loading closure.",
        "observable": "It fixes superposition and leading response to density; nonlinear terms create amplitude-dependent response.",
        "changed": "The exact Helmholtz form and density scaling generally cease to hold outside the linearized regime.",
    },
    {
        "id": "FND-002-A05", "assumption": "Quadratic local recovery energy", "label": "DERIVED",
        "explicit_in_pbuf": "CORE-001 defines it; earlier material asks for stiffness/recovery but does not fix the potential.",
        "derivable": "Yes as the universal leading nonconstant term of a smooth stable energy expanded about a nondegenerate equilibrium.",
        "weaken": "Yes: assume only smoothness, stability, and positive Hessian at equilibrium.",
        "remove": "The exact quadratic model can be removed; some stabilizing recovery is required for finite static susceptibility.",
        "observable": "Leading response is linear; higher powers produce nonlinear saturation away from equilibrium.",
        "changed": "Without positive local curvature the unloaded state may be unstable or the Helmholtz mass term vanishes.",
    },
    {
        "id": "FND-002-A06", "assumption": "Nearest-neighbour transmission", "label": "OPTIONAL",
        "explicit_in_pbuf": "No; it is a CORE-001 discretization choice.",
        "derivable": "No. A gradient-squared term is the leading local isotropic spatial correction, but many microscopic couplings generate it.",
        "weaken": "Yes: finite-range or sufficiently decaying, symmetric couplings with a finite second moment.",
        "remove": "Yes: graph, nonlocal, or direct continuum operators may share the same small-wavenumber expansion.",
        "observable": "Higher-order dispersion and short-scale propagation depend on the coupling stencil.",
        "changed": "Long-wave G is renormalized; nonlocal tails or anisotropy can change the continuum operator itself.",
    },
    {
        "id": "FND-002-A07", "assumption": "Overdamped relaxation", "label": "OPTIONAL",
        "explicit_in_pbuf": "No supplied theory selects first-order time evolution; CORE-001 chose it for static relaxation.",
        "derivable": "No. The static Euler-Lagrange equation is independent of whether relaxation is overdamped, inertial, or otherwise equilibrating.",
        "weaken": "Yes: require dynamics that converge to stationary extrema of the effective energy.",
        "remove": "Yes for the static milestone and frozen laboratory.",
        "observable": "Transient response, mode spectrum, causality, and damping differ.",
        "changed": "Nothing static if the same equilibrium is reached; time-dependent predictions change completely.",
    },
    {
        "id": "FND-002-A08", "assumption": "Choice of coarse-graining kernel", "label": "OPTIONAL",
        "explicit_in_pbuf": "Coarse graining is required by CORE-001, but no supplied PBUF source fixes a kernel.",
        "derivable": "No unique kernel. Normalization and symmetry constraints follow from preserving constants and rotations.",
        "weaken": "Yes: any normalized, localized operator with vanishing first moment and controlled higher moments.",
        "remove": "A distinct kernel can be removed in a direct continuum formulation; a scale-separation map remains conceptually necessary for discrete models.",
        "observable": "Kernel moments control smoothing and finite-resolution corrections, not the leading constant/slow-field limit.",
        "changed": "A nonnormalized or biased kernel fails constant preservation; anisotropic kernels imprint directional artifacts.",
    },
    {
        "id": "FND-002-A09", "assumption": "Scalar continuum field u(x)", "label": "SUPPORTED",
        "explicit_in_pbuf": "Yes at the existing continuum/frozen-laboratory interface; earlier discovery material still lists scalar versus tensor character as open microscopically.",
        "derivable": "No fundamental scalar ontology follows. A scalar effective field follows conditionally when the laboratory couples to one invariant projection and other modes decouple.",
        "weaken": "Yes: u may be the scalar observable sector of vector or tensor internal states.",
        "remove": "Not without changing the frozen laboratory interface; it can be emergent rather than fundamental.",
        "observable": "Additional unsuppressed tensor/vector sectors would generate responses absent from a single scalar field.",
        "changed": "The current constitutive and gradient interfaces require revision if extra continuum modes couple appreciably.",
    },
]

DERIVATIONS = [
    ("FND-002-D01", "Stable smooth local energy at q=0", "V(q)=V(0)+q^T H q/2+O(|q|^3), H>0", "Quadratic recovery at leading order", "Fails at a degenerate/critical or nonsmooth equilibrium"),
    ("FND-002-D02", "Analytic local matter interaction", "I(eta,q)=I(0,0)-eta h.q+O(eta|q|^2,eta^2)", "Linear loading at leading order", "A symmetry may force h=0"),
    ("FND-002-D03", "Local isotropic spatial response", "Gamma(k)=K+G|k|^2+O(|k|^4)", "K u-G Laplacian(u)=s at long wavelength", "Long-range kernels can be nonanalytic; anisotropy makes G a tensor"),
    ("FND-002-D04", "Normalized centered localized kernel", "C[u]=u+(mu_2/2d)Laplacian(u)+higher moments", "Kernel-independent leading slow field", "No scale separation or divergent moments"),
    ("FND-002-D05", "One sourced light mode; others massive/decoupled", "u=e.q; integrate out transverse modes", "Scalar effective sector", "Multiple sourced light modes require vector/tensor closure"),
    ("FND-002-D06", "Any equilibrating dynamics for the same F", "stationary state satisfies delta F/delta u=0", "Static equation independent of overdamped choice", "Driven/non-equilibrium states need a dynamical postulate"),
]

ALTERNATIVES = [
    ("Regular lattice", "site state + symmetric finite differences", "Yes", "Adds regulator geometry and stencil choices", "Representative, not preferred foundationally"),
    ("Irregular graph/network", "node state + weighted graph Laplacian", "Yes, if graph homogenizes isotropically", "Adds graph ensemble/weights", "Equivalent universality class"),
    ("Continuum microfield", "local effective energy with gradient expansion", "Yes, directly", "Avoids lattice and kernel postulates; assumes locality/cutoff implicitly", "Preferred minimal effective formulation"),
    ("Tensor/vector internal state", "one scalar sourced projection plus decoupled modes", "Yes conditionally", "Adds mode masses and couplings", "No benefit until extra observables require it"),
    ("Alternative coarse graining", "normalized centered convolution, block average, spectral low-pass", "Yes at leading order", "Finite-scale corrections differ", "Operator is not uniquely selectable here"),
]

MINIMAL_POSTULATES = [
    "P1: A microscopic PBUF state exists and admits a stable, approximately local, rotationally invariant long-wavelength sector.",
    "P2: Matter sources at least one dimensionless scalar observable sector u; other microscopic components are absent, massive, or decoupled at laboratory scales.",
    "P3: The effective static response has positive finite local susceptibility and a positive leading spatial stiffness.",
    "P4: PBUF's exactly-three-state and characteristic-1/137 statements are retained as ontology premises, but existing material supplies no derivation and the scalar closure cannot identify them separately.",
]


def validate() -> dict:
    ids = [row["id"] for row in AUDIT]
    labels = [row["label"] for row in AUDIT]
    required = {"explicit_in_pbuf", "derivable", "weaken", "remove", "observable", "changed"}
    checks = {
        "exactly_nine_assumptions": len(AUDIT) == 9,
        "unique_assumption_ids": len(ids) == len(set(ids)),
        "all_labels_allowed": set(labels) <= LABELS,
        "all_four_labels_used": set(labels) == LABELS,
        "all_required_questions_answered": all(required <= row.keys() and all(row[k].strip() for k in required) for row in AUDIT),
        "alternatives_complete": len(ALTERNATIVES) == 5,
        "derivation_boundaries_documented": all(len(row) == 5 and row[-1] for row in DERIVATIONS),
    }
    return {"checks": checks, "all_checks_pass": all(checks.values()), "label_counts": {label: labels.count(label) for label in sorted(LABELS)}}


def _table(headers: list[str], rows: list[list[str]]) -> str:
    return "| " + " | ".join(headers) + " |\n|" + "|".join("---" for _ in headers) + "|\n" + "\n".join("| " + " | ".join(str(v).replace("|", "\\|") for v in row) + " |" for row in rows)


def report(validation: dict) -> str:
    audit_rows = [[r["id"], r["assumption"], r["label"], r["explicit_in_pbuf"], r["derivable"], r["weaken"], r["remove"], r["observable"], r["changed"]] for r in AUDIT]
    return f"""# PBUF FND-002 — Justification of the microscopic state

## Result

All nine CORE-001 assumptions are documented and assigned exactly one required label. The audit finds one universal leading-order derivation, five optional realizations, one interface-supported statement, and two new foundational postulates. No supplied PBUF equation derives either the number three or `1/137`; their appearance in CORE-001 establishes provenance, not mathematical support.

## Assumption audit matrix

{_table(['ID', 'Assumption', 'Label', 'Explicit?', 'Derivable?', 'Can weaken?', 'Can remove?', 'Observable consequence', 'What breaks?'], audit_rows)}

`DERIVED` means derivable only under the stated regularity and stability conditions, not derivable as an exact microscopic polynomial. `SUPPORTED` means required by the frozen interface and compatible with supplied PBUF material; it is not proof of a fundamental ontology. `OPTIONAL` marks a replaceable realization. `NEW POSTULATE` marks a claimed foundational fact that the supplied theory does not derive.

## Derivation matrix

{_table(['ID', 'Premises', 'Expansion/reduction', 'Consequence', 'Boundary'], [list(r) for r in DERIVATIONS])}

These reductions define a universality class. They derive the *leading continuum form*, not the exact microscopic substrate or coefficient values.

## Alternative models

{_table(['Model', 'Construction', 'Same continuum?', 'Assumption cost', 'Finding'], [list(r) for r in ALTERNATIVES])}

## Minimal postulate set and revised model

""" + "\n".join(f"{i}. {p}" for i, p in enumerate(MINIMAL_POSTULATES, 1)) + """

The preferred revised model is therefore effective and regulator-independent:

`F_eff[u] = integral [K u^2/2 + G |grad u|^2/2 - s(rho)u + O(u^3, u grad^2 u, grad^4)] dx`, with `K>0`, `G>0`.

Its stationary leading-order equation is `K u-div(G grad u)=s(rho)`. A lattice, nearest-neighbour stencil, Gaussian kernel, and overdamped dynamics are examples that realize this equation; none is foundational. A multi-component microscopic state is permitted, but only its sourced scalar projection belongs in the current continuum interface. This revision changes no frozen weak-lensing code or parameter.

## Irreducible postulates

The irreducible content is P1–P3 for the effective theory. P4 is additionally irreducible only if the distinctive PBUF ontology—exactly three microscopic degrees and the `1/137` association—is demanded. In corrected CORE-001, `g_dev` directly normalizes matter loading and is no longer rescaling-degenerate with an auxiliary coupling. Its numerical value is still a premise rather than a result of this assumption audit, and the two transverse components have no independent observable in the frozen laboratory.

## Recommendation for FND-003

FND-003 should be a **symmetry and identifiability derivation**, not another regulator construction. It should specify the symmetry group and representation of the three microscopic degrees and define how stress-energy (not only static density) transforms and couples. Its decision tests should be: (1) derive three and `1/137` from those structures, or formally retain them as axioms; (2) determine whether extra modes decouple; and (3) derive coefficient signs and scaling without lensing fits. Dynamic, causal evolution should be deferred until the static representation and coupling are fixed.

## Completion checks

All checks pass: **{validation['all_checks_pass']}**. Label counts: `{validation['label_counts']}`. Scope remained theory-only.
"""


def main(output: Path) -> None:
    output.mkdir(parents=True, exist_ok=True)
    validation = validate()
    if not validation["all_checks_pass"]:
        raise RuntimeError("FND-002 audit validation failed")
    with (output / "assumption_audit.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(AUDIT[0]))
        writer.writeheader(); writer.writerows(AUDIT)
    with (output / "derivation_matrix.csv").open("w", newline="") as handle:
        writer = csv.writer(handle); writer.writerow(["id", "premises", "reduction", "consequence", "boundary"]); writer.writerows(DERIVATIONS)
    with (output / "alternative_models.csv").open("w", newline="") as handle:
        writer = csv.writer(handle); writer.writerow(["model", "construction", "same_continuum", "assumption_cost", "finding"]); writer.writerows(ALTERNATIVES)
    record = {"mission": "PBUF FND-002 Justification of the Microscopic State", "status": "complete_theory_audit", "scope": "theory-only", "assumption_audit": AUDIT, "derivations": DERIVATIONS, "alternatives": ALTERNATIVES, "minimal_postulates": MINIMAL_POSTULATES, "preferred_model": "regulator-independent scalar effective field", "validation": validation}
    (output / "fnd002_audit.json").write_text(json.dumps(record, indent=2) + "\n")
    (output / "validation.json").write_text(json.dumps(validation, indent=2) + "\n")
    (output / "fnd002_report.md").write_text(report(validation))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=Path("runs/fnd002"))
    main(parser.parse_args().output)
