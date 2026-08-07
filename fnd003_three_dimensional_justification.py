#!/usr/bin/env python3
"""Produce the theory-only PBUF FND-003 dimensional-state audit.

This module deliberately has no dependency on the weak-lensing experiment.  It
tests what follows from three-dimensional spatial ontology, and keeps conditional
representation-theory results separate from PBUF premises and new postulates.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


TRACEABILITY = [
    {
        "id": "FND-003-T01", "claim": "Physical space has three spatial dimensions", "classification": "working premise",
        "source": "FND-003 mission; stated PBUF ontology", "justification": "Taken as the starting ontology, not derived by this milestone.",
    },
    {
        "id": "FND-003-T02", "claim": "The microscopic state exists", "classification": "explicitly stated in PBUF",
        "source": "CORE-001 and FND-002", "justification": "A PBUF premise; existence alone fixes neither representation nor component count.",
    },
    {
        "id": "FND-003-T03", "claim": "The microscopic state has exactly three real components", "classification": "new assumption introduced after V11",
        "source": "CORE-001 working premise", "justification": "Three-dimensional space permits scalar, vector, tensor, and larger microscopic states; dimension matching is not a theorem.",
    },
    {
        "id": "FND-003-T04", "claim": "Each component is associated with one spatial direction", "classification": "new assumption introduced after V11",
        "source": "FND-003 candidate interpretation", "justification": "Requires q to be a spatial vector/covector and a choice of basis; components are basis-dependent, while q is geometric.",
    },
    {
        "id": "FND-003-T05", "claim": "A faithful linear SO(3) vector realization needs at least three real components", "classification": "mathematically derived",
        "source": "conditional representation argument", "justification": "SO(3) has no faithful real representation in dimensions 1 or 2; its defining three-dimensional representation is faithful. The premise that q must be faithful and linear is additional.",
    },
    {
        "id": "FND-003-T06", "claim": "The three-component realization is unique", "classification": "not established",
        "source": "model comparison", "justification": "Even after selecting the vector representation, dynamics, parity, locality, couplings, and field content remain open; non-vector realizations also satisfy the base ontology.",
    },
    {
        "id": "FND-003-T07", "claim": "g_dev=1/137 is an intrinsic microscopic coupling", "classification": "working premise",
        "source": "CORE-001 and FND-002", "justification": "No supplied symmetry, normalization, action, or renormalization prescription derives or operationally isolates it.",
    },
    {
        "id": "FND-003-T08", "claim": "g_dev directly normalizes the corrected CORE-001 linear source", "classification": "mathematically derived",
        "source": "corrected CORE-001 energy", "justification": "The source term is -epsilon_* g_dev eta e.q; no independent coupling multiplier or inverse-rescaling freedom remains.",
    },
    {
        "id": "FND-003-T09", "claim": "The current theory identifies 1/137 as a separately measurable parameter", "classification": "not established",
        "source": "corrected identifiability audit", "justification": "The normalized microscopic equation is sensitive directly to g_dev, but no supplied principle derives its numerical value or links it to an independently completed downstream observable.",
    },
    {
        "id": "FND-003-T10", "claim": "u=C_L[e.q] is a scalar continuum deformation", "classification": "mathematically derived",
        "source": "CORE-001 conditional map", "justification": "It is rotationally scalar only when q and the matter-selected e transform together and the kernel is isotropic and normalized.",
    },
    {
        "id": "FND-003-T11", "claim": "A direction-free nonzero linear map from a spatial vector q to a scalar u exists", "classification": "not established",
        "source": "SO(3) invariance", "justification": "No nonzero SO(3)-invariant linear functional on the vector representation exists; one must add e, use derivatives such as div(q), or use a nonlinear invariant such as |q|.",
    },
    {
        "id": "FND-003-T12", "claim": "The long-wave scalar equation K u-div(G grad u)=s(rho) follows", "classification": "mathematically derived",
        "source": "CORE-001/FND-002", "justification": "Conditional on the projection, stable analytic energy, locality, isotropy, scale separation, and decoupling of other modes; coefficients and source law are not predicted.",
    },
]

STATE_DEFINITION = {
    "geometric_state": "q(x) in T_x Sigma (or T_x^* Sigma), with dim(Sigma)=3",
    "components": "q^a=e^a(q), a=1,2,3, in a local orthonormal frame",
    "rotation_law": "q'^a=R^a_b q^b for R in SO(3)",
    "status": "minimal conditional realization, not a consequence of dimensionality alone",
    "minimality_conditions": ["linear real state", "faithful action of spatial rotations", "independent directional response", "no additional internal sectors"],
}

CONTINUUM_MAPPINGS = [
    ("Matter-selected projection", "u(x)=integral W_L(x-y) e(y).q(y) dy", "Linear and matches CORE-001", "Requires an additional transforming vector e; a fixed background e breaks isotropy"),
    ("Longitudinal scalar", "u(x)=L integral W_L(x-y) div(q(y)) dy", "Rotation scalar without an internal direction", "Introduces a derivative and length L; does not equal the CORE-001 map"),
    ("Magnitude scalar", "u(x)=integral W_L(x-y) |q(y)| dy", "Direction-free rotation scalar", "Nonlinear, nonnegative, and changes the weak-field/source expansion"),
]

POSTULATES = [
    {"id": "P1", "text": "The microscopic state carries the defining spatial-vector (or covector) representation of SO(3), rather than a scalar, tensor, spinorial, or unrelated internal representation.", "needed_for": "identifying three components with three spatial directions"},
    {"id": "P2", "text": "The rotation action is linear and faithful and all three directional responses are independent; no additional microscopic sectors are required.", "needed_for": "conditional minimality of three"},
    {"id": "P3", "text": "A covariant scalarization mechanism is selected: a matter-provided vector e, a longitudinal derivative, or a nonlinear invariant.", "needed_for": "mapping q to the scalar continuum field u"},
    {"id": "P4", "text": "The stable, local, isotropic long-wave expansion and decoupling assumptions of CORE-001/FND-002 hold.", "needed_for": "recovering K u-div(G grad u)=s(rho)"},
    {"id": "P5", "text": "If g_dev=1/137 is to be derived rather than postulated, a microscopic principle must fix its value and specify any applicable scale behavior.", "needed_for": "turning 1/137 from a direct PBUF premise into a prediction"},
]


def validate() -> dict:
    classes = {row["classification"] for row in TRACEABILITY}
    ids = [row["id"] for row in TRACEABILITY]
    checks = {
        "all_claims_classified": all(row["classification"] and row["source"] and row["justification"] for row in TRACEABILITY),
        "required_class_boundaries_present": {"explicitly stated in PBUF", "mathematically derived", "working premise", "new assumption introduced after V11", "not established"} <= classes,
        "unique_traceability_ids": len(ids) == len(set(ids)),
        "three_state_is_not_marked_unconditional_derivation": next(row for row in TRACEABILITY if row["id"] == "FND-003-T03")["classification"] != "mathematically derived",
        "one_over_137_is_not_marked_derived": next(row for row in TRACEABILITY if row["id"] == "FND-003-T07")["classification"] != "mathematically derived",
        "direct_g_dev_source_recorded": any("directly normalizes" in row["claim"] for row in TRACEABILITY),
        "scalarization_alternatives_audited": len(CONTINUUM_MAPPINGS) == 3,
        "irreducible_postulates_listed": len(POSTULATES) == 5,
        "frozen_laboratory_untouched_by_design": True,
    }
    return {"checks": checks, "all_checks_pass": all(checks.values())}


def table(headers: list[str], rows: list[list[str]]) -> str:
    return "| " + " | ".join(headers) + " |\n|" + "|".join("---" for _ in headers) + "|\n" + "\n".join("| " + " | ".join(str(value).replace("|", "\\|") for value in row) + " |" for row in rows)


def derivation_text() -> str:
    return """# FND-003 mathematical derivation and boundary

## Conditional minimality theorem

Let the microscopic state at each point be a finite-dimensional real vector space `V`, and require spatial rotations to act through a continuous, linear, faithful representation `rho: SO(3) -> GL(V)`. A one-dimensional real representation is trivial because connected `SO(3)` has no nontrivial continuous homomorphism to `R*`. A two-dimensional faithful representation would embed the three-dimensional Lie group `SO(3)` into `GL(2,R)`; its maximal compact subgroup is conjugate to `O(2)`, whose connected part is only one-dimensional, so no such embedding exists. The defining action on `R^3` is faithful. Therefore `dim(V)>=3`, and `V=R^3` realizes the minimum.

This derives **three only from the added faithful-linear-vector premises**. Three-dimensional space by itself does not require a microscopic state to transform faithfully: a scalar state is one-dimensional, a symmetric tensor has six components, and multiple fields give arbitrarily many. The PBUF ontology supplied here does not establish those added premises, so the overall milestone cannot promote the conditional theorem to an unconditional derivation.

## Scalarization obstruction

For a linear scalar map `l:R^3->R` to be rotation invariant, `l(Rq)=l(q)` for every `R`. Writing `l(q)=e.q` implies `R^T e=e` for every rotation, hence `e=0`. Thus no nonzero direction-free linear scalarization exists. CORE-001's `u=C_L[e.q]` is covariant only if `e` is additional physical data that transforms with `q`. The divergence and norm alternatives avoid a fixed direction but define different effective theories.

## Coupling identifiability

Corrected CORE-001 contains `-epsilon_* g_dev eta e.q`. The normalized microscopic source therefore depends directly on `g_dev`, and the former inverse-rescaling degeneracy is withdrawn. This makes `g_dev` identifiable inside a fully specified microscopic response calculation, but it does not derive the value `1/137`: no supplied symmetry or dynamical principle selects that number, and downstream closure and access maps remain incomplete.
"""


def report(validation: dict) -> str:
    trace_rows = [[r["id"], r["claim"], r["classification"], r["source"], r["justification"]] for r in TRACEABILITY]
    mapping_rows = [list(row) for row in CONTINUUM_MAPPINGS]
    return f"""# PBUF FND-003 — Three-dimensional microscopic state justification

## Result: Outcome C

The existing ontology is compatible with a three-component spatial-vector microscopic state, but it does not uniquely require one. Three spatial dimensions determine the dimension of tangent vectors, not the representation carried by an otherwise unspecified microscopic state. Exactly three follows only after adding that the state transforms linearly and faithfully as a spatial vector and independently resolves all directions. Those are new postulates, so Outcome A and Outcome B are not justified.

## Ontology-to-mathematics mapping

{table(['ID', 'Claim', 'Classification', 'Source', 'Reason/boundary'], trace_rows)}

## Three-dimensional state definition

The minimal conditional realization is `{STATE_DEFINITION['geometric_state']}`. In an orthonormal frame its coordinates are `q^a`, `a=1,2,3`, with `{STATE_DEFINITION['rotation_law']}`. The directions label geometric basis components, not three invariant substances. Changing frame mixes them. No extra internal degree of freedom is introduced in this realization, but excluding scalar, tensor, spinorial, or additional sectors is itself a model choice.

The conditional representation argument in `mathematical_derivation.md` proves that three is minimal under the faithful-linear-vector premises. It also records why those premises do not follow from dimensionality alone.

## Mapping to continuum deformation

{table(['Map', 'Definition', 'Advantage', 'Boundary'], mapping_rows)}

The existing CORE-001 interface selects the first map. With normalized isotropic `W_L`, a transforming unit vector `e`, scale separation, stable local response, and decoupled transverse modes, it yields the scalar effective energy and stationary equation `K u-div(G grad u)=s(rho)`. This is a conditional coarse-graining derivation, not a derivation of the vector ontology or the coefficients. A spatially fixed `e` would add preferred-direction structure absent from the base ontology.

## Treatment of 1/137

In the corrected microscopic energy, `g_dev=1/137` appears directly in the matter vertex. The previous inverse-rescaling argument was an artefact of an auxiliary modelling choice and is withdrawn. The current record therefore classifies `1/137` as a direct working premise, not a mathematically derived value. Deriving it requires a microscopic principle that fixes the coupling and its applicable scale behavior, plus a completed response/access chain for empirical identifiability.

## Remaining irreducible postulates

""" + "\n".join(f"{p['id']}. {p['text']} Needed for: {p['needed_for']}." for p in POSTULATES) + f"""

## Recommendation for FND-004

FND-004 should formulate the covariant microscopic action and representation content. It should choose and justify the SO(3) representation, derive rather than select the scalarization/source coupling, test whether a matter-provided `e` creates forbidden preferred-direction effects, and state any applicable scale behavior for `g_dev`. Its gates should require (1) a unique representation or an explicit postulate, (2) a rotation-covariant map to `u`, (3) a mode-decoupling calculation, and (4) a direct, operationally testable coupling. No empirical fitting or weak-lensing changes should occur until those gates pass.

## Completion checks

All checks pass: **{validation['all_checks_pass']}**. The work is theory-only and imports or modifies no frozen-laboratory module.
"""


def main(output: Path) -> None:
    output.mkdir(parents=True, exist_ok=True)
    validation = validate()
    if not validation["all_checks_pass"]:
        raise RuntimeError("FND-003 validation failed")
    with (output / "ontology_traceability.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(TRACEABILITY[0]))
        writer.writeheader()
        writer.writerows(TRACEABILITY)
    record = {
        "mission": "PBUF FND-003 Three-Dimensional Microscopic State Justification",
        "scope": "theory-only", "outcome": "C", "observational_claim": False,
        "three_dimensional_state": STATE_DEFINITION, "continuum_mappings": CONTINUUM_MAPPINGS,
        "coupling_1_over_137": {"status": "direct working premise", "identified_parameter": "g_dev", "intrinsic_or_derived": False},
        "remaining_postulates": POSTULATES, "recommendation": "FND-004 covariant representation, scalarization, mode-decoupling, and coupling-identifiability development",
        "validation": validation,
    }
    (output / "fnd003_analysis.json").write_text(json.dumps(record, indent=2) + "\n")
    (output / "remaining_postulates.json").write_text(json.dumps(POSTULATES, indent=2) + "\n")
    (output / "mathematical_derivation.md").write_text(derivation_text())
    (output / "validation.json").write_text(json.dumps(validation, indent=2) + "\n")
    (output / "foundation_report.md").write_text(report(validation))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=Path("runs/fnd003"))
    main(parser.parse_args().output)
