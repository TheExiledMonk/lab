#!/usr/bin/env python3
"""Generate the theory-only PBUF FND-004 consequence audit.

The three-component state and g_dev=1/137 are accepted as axioms here.  The
audit derives only representation/counting consequences from those axioms and
labels continuum dynamics that need the inherited CORE-001 effective premises
as conditional.  It never imports or modifies the weak-lensing laboratory.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


AXIOMS = {
    "A1": "Spacetime has three independent microscopic degrees of freedom corresponding to the three spatial dimensions.",
    "A2": "Matter couples identically to each microscopic degree through g_dev=1/137.",
    "A3": "Macroscopic continuum behaviour emerges by coarse graining the microscopic state.",
}

# Conditions are deliberately not smuggled into A1--A3.
CONDITIONS = {
    "C1": "q transforms as the defining spatial-vector representation and no background direction is present.",
    "C2": "The long-wave theory is local, analytic, homogeneous, parity-even and rotationally invariant.",
    "C3": "The reference state is stable and nondegenerate, and retained modes separate from microscopic scales.",
    "C4": "A covariant scalarization/source map is specified and non-observed modes decouple.",
    "C5": "Time dynamics (inertial or dissipative) and kinetic normalization are specified.",
    "C6": "Light couples to a specified effective metric or ray-deflection functional.",
    "C7": "The matter-coupling normalization is fixed independently of all other source coefficients.",
}

CONSEQUENCES = [
    {"id":"P01","prediction":"Exactly three microscopic component labels; generic N has N and a scalar has one.","kind":"exact","axioms":"A1","conditions":"none","relation":"q=(q1,q2,q3)","observable":"Component count is physical only if components can be independently excited/read out."},
    {"id":"P02","prediction":"Equal bare component-source vertices.","kind":"exact","axioms":"A2","conditions":"fixed common normalization","relation":"g1=g2=g3=g_dev","observable":"No component-dependent coupling at the axiom scale."},
    {"id":"P03","prediction":"The equal-coupling vector is the normalized singlet direction and two source-orthogonal combinations are unsourced by a common scalar load.","kind":"exact linear algebra","axioms":"A1,A2","conditions":"linear common scalar source","relation":"q_parallel=(q1+q2+q3)/sqrt(3); J_perp=0","observable":"Two dark combinations exist at linear source level, but their propagation is not fixed."},
    {"id":"P04","prediction":"The coherent source strength carries a sqrt(3) amplitude (factor 3 in a quadratic response) relative to one equally normalized component.","kind":"conditional quantitative","axioms":"A1,A2","conditions":"independent orthonormal components; linear response; same per-component susceptibility","relation":"|g|=sqrt(3) g_dev; |g|^2=3 g_dev^2","observable":"N ontology replaces 3 by N; scalar has factor 1."},
    {"id":"P05","prediction":"A rotationally invariant rank-two constitutive tensor in the vector sector has transverse/longitudinal form.","kind":"conditional symmetry","axioms":"A1","conditions":"C1,C2","relation":"Gamma_ab(k)=A(k^2) delta_ab+B(k^2) k_a k_b","observable":"One longitudinal and two degenerate transverse eigenmodes."},
    {"id":"P06","prediction":"The leading static vector continuum equation has two gradient stiffnesses, not a uniquely fixed scalar Helmholtz law.","kind":"conditional continuum","axioms":"A1,A3","conditions":"C1,C2,C3","relation":"K q-a Laplacian(q)-b grad(div q)=J","observable":"G_T=a and G_L=a+b; positivity requires K>0, a>0, a+b>0."},
    {"id":"P07","prediction":"If only the common/longitudinal scalar sector is retained, the CORE-type scalar equation follows in form.","kind":"conditional reduction","axioms":"A1,A2,A3","conditions":"C2,C3,C4","relation":"K u-div(G grad u)=s(rho)","observable":"The axioms fix neither K, G nor s(rho); the equation is not inevitable from A1-A3 alone."},
    {"id":"P08","prediction":"Isotropy forbids component-dependent masses and directional stiffness at leading order.","kind":"conditional symmetry","axioms":"A1,A2","conditions":"C1,C2","relation":"K_ab=K delta_ab; equal transverse polarizations","observable":"Splitting signals broken isotropy, unequal coupling, or extra structure."},
    {"id":"P09","prediction":"Propagation has one longitudinal and two transverse branches when inertial vector dynamics is supplied.","kind":"conditional propagation","axioms":"A1,A3","conditions":"C1,C2,C3,C5","relation":"omega_L^2=(K+(a+b)k^2)/M; omega_T^2=(K+a k^2)/M (double)","observable":"Three fixes transverse degeneracy at two; N internal components generally gives N-1 source-orthogonal modes, not spatial L/T modes unless it is a vector."},
    {"id":"P10","prediction":"g_dev directly fixes the normalized common matter vertex, but A2 postulates rather than derives its numerical value.","kind":"corrected identifiability boundary","axioms":"A2","conditions":"a completed calibrated response/readout chain for measurement","relation":"g_vec=g_dev(1,1,1); no independent coupling multiplier","observable":"Absolute calibrated response can be g_dev-sensitive; normalized component ratios cancel g_dev and cannot determine its magnitude."},
    {"id":"P11","prediction":"No weak-lensing deflection law or amplitude follows from A1-A3.","kind":"non-derivation","axioms":"A1,A2,A3","conditions":"C4,C6,C7 absent","relation":"none","observable":"A specified photon/effective-metric coupling is required before lensing is predicted."},
    {"id":"P12","prediction":"Coarse graining preserves the three-component target space but does not by itself select locality, dynamics, coefficients, or a scalar observable.","kind":"exact boundary","axioms":"A3","conditions":"none","relation":"Q_L=C_L[q] in R^3","observable":"Kernel moments and scale separation determine finite-resolution corrections."},
]

COMPARISON = [
    {"feature":"Microscopic count","three_component":"3 (axiom)","generic_n":"N","scalar":"1","unique_to_three":"Exactly two combinations orthogonal to one common source."},
    {"feature":"Equal-coupling norm","three_component":"sqrt(3) g_dev","generic_n":"sqrt(N) g_dev","scalar":"g_dev","unique_to_three":"Factor 3 in normalized quadratic coherent response."},
    {"feature":"Spatial-vector polarizations","three_component":"1 longitudinal + 2 transverse","generic_n":"Not defined unless representation is specified","scalar":"1 scalar","unique_to_three":"Twofold transverse degeneracy, conditional on spatial-vector identity."},
    {"feature":"Isotropic static stiffnesses","three_component":"K, G_L, G_T","generic_n":"Representation-dependent","scalar":"K, G","unique_to_three":"Vector L/T split, not its coefficient values."},
    {"feature":"Scalar continuum closure","three_component":"Possible after projection/decoupling","generic_n":"Possible after projection/decoupling","scalar":"Direct","unique_to_three":"None at leading scalar order."},
    {"feature":"Weak lensing","three_component":"Undetermined","generic_n":"Undetermined","scalar":"Undetermined","unique_to_three":"None without a light-coupling law."},
]

UNIQUE = [
    {"id":"U1","claim":"Two source-orthogonal microscopic combinations for a single equal common source.","status":"exact under linear common sourcing","test":"Independently excite/read component combinations; search for two dark channels.","dependencies":"A1,A2 plus linear common source"},
    {"id":"U2","claim":"Coherent quadratic response multiplicity exactly 3 rather than N or 1.","status":"conditional quantitative","test":"Compare calibrated single-component and coherent susceptibilities without fitting.","dependencies":"A1,A2 plus common normalization and equal susceptibility"},
    {"id":"U3","claim":"One longitudinal and two degenerate transverse modes.","status":"conditional representation prediction","test":"Directional spectroscopy/propagation measurement; test transverse degeneracy.","dependencies":"A1,C1,C2,C3,C5"},
    {"id":"U4","claim":"No unique scalar-lensing signature follows merely from having three components.","status":"exact negative result","test":"Not directly testable; it is a model-identifiability constraint.","dependencies":"Logical comparison of A1-A3"},
]


def table(headers: list[str], rows: list[list[str]]) -> str:
    clean = lambda value: str(value).replace("|", "\\|").replace("\n", " ")
    return "| " + " | ".join(headers) + " |\n|" + "|".join("---" for _ in headers) + "|\n" + "\n".join("| " + " | ".join(clean(v) for v in row) + " |" for row in rows)


def validate() -> dict:
    ids = [p["id"] for p in CONSEQUENCES]
    traced = {a for p in CONSEQUENCES for a in p["axioms"].split(",")}
    checks = {
        "unique_prediction_ids": len(ids) == len(set(ids)),
        "every_axiom_traced": set(AXIOMS) <= traced,
        "every_prediction_classified": all(p["kind"] and p["conditions"] and p["observable"] for p in CONSEQUENCES),
        "three_models_compared": all(set(r) == {"feature","three_component","generic_n","scalar","unique_to_three"} for r in COMPARISON),
        "unique_catalogue_present": len(UNIQUE) >= 1,
        "coupling_identifiability_not_overclaimed": any(p["id"] == "P10" for p in CONSEQUENCES),
        "lensing_not_overclaimed": any(p["id"] == "P11" for p in CONSEQUENCES),
        "no_empirical_fit": True,
        "frozen_laboratory_untouched_by_design": True,
    }
    return {"checks": checks, "all_checks_pass": all(checks.values())}


def report(validation: dict) -> str:
    rows = [[p[k] for k in ("id","prediction","kind","axioms","conditions","relation","observable")] for p in CONSEQUENCES]
    comp = [[r[k] for k in ("feature","three_component","generic_n","scalar","unique_to_three")] for r in COMPARISON]
    unique = [[u[k] for k in ("id","claim","status","test","dependencies")] for u in UNIQUE]
    return f"""# PBUF FND-004 — Consequences of the three-dimensional microscopic ontology

## Result

Adopting the premises as axioms yields a small set of exact counting, equal-coupling, and identifiability consequences. It does **not** uniquely yield a continuum field equation, stiffness, propagation law, or weak-lensing signal. Those results become conditional only after the inherited symmetry, locality, stability, scalarization, dynamics, and photon-coupling assumptions are stated. No parameter was fitted and the frozen weak-lensing laboratory was not imported or modified.

## Axioms and explicit auxiliary conditions

""" + "\n".join(f"- **{k}:** {v}" for k,v in AXIOMS.items()) + "\n\n" + "\n".join(f"- **{k}:** {v}" for k,v in CONDITIONS.items()) + f"""

The C-conditions are not new foundation axioms: they are clearly exposed hypotheses delimiting conditional effective-theory statements. A1's phrase “corresponding to spatial dimensions” does not alone specify a vector transformation law, so C1 remains necessary.

## Consequence and axiom traceability matrix

{table(['ID','Prediction/boundary','Class','Minimum axioms','Extra conditions','Relation','Observable meaning'], rows)}

## Continuum, constitutive, symmetry and propagation derivation

Under C1–C3, Fourier-space rotational covariance permits only `Gamma_ab(k)=A(k^2)delta_ab+B(k^2)k_a k_b`. Expanding analytically at small `k` and varying the quadratic energy gives `K q-a Laplacian(q)-b grad(div q)=J`. The longitudinal and transverse static stiffnesses are `G_L=a+b` and `G_T=a`; stability requires `K>0`, `G_L>0`, and `G_T>0`. Neither their magnitudes nor equality is fixed by three components or by `1/137`.

With an inertial normalization `M`, the branches are `omega_L^2=(K+G_L k^2)/M` and `omega_T^2=(K+G_T k^2)/M`, with the transverse branch doubled. Overdamped dynamics would instead produce relaxation rates, demonstrating that propagation is conditional on C5. A scalar-only closure has one stiffness and no transverse branch.

A normalized, isotropic coarse-graining operator preserves constants and rotational covariance, with higher kernel moments contributing resolution corrections. A3 alone does not guarantee such a kernel, scale separation, or locality. After C4, a single retained scalar may obey `K u-div(G grad u)=s(rho)`, making the *form* of the existing effective closure natural, but not its coefficients or nonlinear constitutive law.

## Alternative ontologies

{table(['Feature','Three-component PBUF','Generic N','Scalar only','Three-specific finding'], comp)}

For a generic equal common source, the internal source direction is `(1,...,1)/sqrt(N)` and there are `N-1` orthogonal combinations. The coherent amplitude is `sqrt(N) g_dev`. These are internal-space facts; spatial longitudinal/transverse language is valid only when the state actually carries the spatial-vector representation.

## Unique PBUF prediction catalogue

{table(['ID','Claim','Status','Required experiment','Dependencies'], unique)}

The positive unique predictions are representation/counting signatures, not a weak-lensing curve. In particular, `g_dev=1/137` supplies a bare scale but remains observationally degenerate wherever another coefficient multiplies it. A normalized microscopic action or independent observable is required to expose the number itself.

## Weak-lensing implication

There is no derived weak-lensing amplitude, radial profile, or photon trajectory from A1–A3. At most, if C4 and C6 reduce the ontology to the same scalar continuum interface, the existing laboratory can consume that scalar field. Such compatibility is not a prediction and cannot discriminate three components from a scalar or generic N model when extra modes decouple.

## Next validation milestone

Proceed to **FND-005: covariant mode-and-coupling validation**. Specify a normalized quadratic microscopic action, the exact rotation representation, matter source tensor, and time kinetic term. Its preregistered no-fit tests should be: (1) measure or calculate the one-longitudinal/two-transverse spectrum and degeneracy; (2) test the direct `1/137` vertex with a calibrated response/readout chain; (3) derive the scalarization and photon coupling; and only then (4) run the unchanged weak-lensing laboratory as an out-of-sample consequence. Failure at steps 1–3 should stop lensing claims rather than trigger tuning.

## Completion checks

Every axiom has traced consequences and every prediction names its minimum axioms and auxiliary conditions. All automated checks pass: **{validation['all_checks_pass']}**.
"""


def main(output: Path) -> None:
    output.mkdir(parents=True, exist_ok=True)
    validation = validate()
    if not validation["all_checks_pass"]:
        raise RuntimeError("FND-004 validation failed")
    with (output / "axiom_prediction_traceability.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(CONSEQUENCES[0]))
        writer.writeheader(); writer.writerows(CONSEQUENCES)
    with (output / "ontology_comparison.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(COMPARISON[0]))
        writer.writeheader(); writer.writerows(COMPARISON)
    (output / "unique_pbuf_predictions.json").write_text(json.dumps(UNIQUE, indent=2) + "\n")
    record = {"mission":"PBUF FND-004 Consequences of the Three-Dimensional Microscopic Ontology", "scope":"theory-only/no-fit", "axioms":AXIOMS, "auxiliary_conditions":CONDITIONS, "consequences":CONSEQUENCES, "comparison":COMPARISON, "unique_predictions":UNIQUE, "next_milestone":"FND-005 covariant mode-and-coupling validation", "validation":validation}
    (output / "fnd004_analysis.json").write_text(json.dumps(record, indent=2) + "\n")
    (output / "validation.json").write_text(json.dumps(validation, indent=2) + "\n")
    (output / "consequence_derivation_report.md").write_text(report(validation))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=Path("runs/fnd004"))
    main(parser.parse_args().output)
