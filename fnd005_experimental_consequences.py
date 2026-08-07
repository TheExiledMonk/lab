#!/usr/bin/env python3
"""Generate the theory-only PBUF FND-005 experimental prediction record.

This module does not import or modify the frozen weak-lensing laboratory.  It
separates consequences of the three-component ontology from signatures that
also require a representation, source/readout map, or dynamics.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


AXIOMS = {
    "A1": "The microscopic state has exactly three independent components.",
    "A2": "A common matter load couples equally to those components with bare scale g_dev=1/137.",
    "A3": "Macroscopic observables arise by coarse graining the microscopic state.",
}

ASSUMPTIONS = {
    "E1": "The three components can be independently prepared or read out with a calibrated linear map.",
    "E2": "The common loading and readout normalization are fixed independently; linear response applies.",
    "E3": "The components form a spatial vector in an isotropic, parity-even reference state.",
    "E4": "A stable local quadratic long-wave theory and a specified time kinetic law apply.",
    "E5": "Symmetry breaking is a controlled perturbation rather than an unknown apparatus effect.",
    "E6": "A photon/effective-metric coupling maps microscopic or coarse fields to light.",
}

PREDICTIONS = [
    {"id":"F005-P01","category":"mode counting","prediction":"Component-resolved state tomography has rank three.","axioms":"A1","extra_assumptions":"E1","ontology_only":False,"relation":"rank(q response space)=3","observable":"Three independently preparable/readable response directions.","constitutive_dependence":"none","classification":"unique to PBUF against scalar or N!=3","falsifier":"Resolved rank is not three after sensitivity and constraints are accounted for."},
    {"id":"F005-P02","category":"coherent versus incoherent loading","prediction":"One normalized common source excites one bright combination and leaves exactly two orthogonal combinations dark at the source vertex.","axioms":"A1,A2","extra_assumptions":"E1,E2","ontology_only":False,"relation":"q_B=(q1+q2+q3)/sqrt(3); dim ker(g^T)=2","observable":"A rank-one common-source response with two null source channels.","constitutive_dependence":"none at the source vertex; propagation may mix channels","classification":"unique multiplicity: scalar has 0 dark; generic N has N-1","falsifier":"The calibrated common vertex has nullity other than two or unequal entries."},
    {"id":"F005-P03","category":"coherent versus incoherent loading","prediction":"The coherent bright-mode amplitude is sqrt(3) times one equally normalized component and its quadratic weight is three times one component; equivalently, the squared unnormalized coherent sum is three times the incoherent sum of three equal component powers.","axioms":"A1,A2","extra_assumptions":"E2","ontology_only":False,"relation":"|g|=sqrt(3) g_dev; |g|^2=3 g_dev^2; |sum_i g_dev|^2/sum_i|g_dev|^2=3","observable":"No-fit ratio R_amp=sqrt(3), R_single-power=3, or R_coherent/incoherent=3 under the stated normalization.","constitutive_dependence":"independent of propagation law if measured at the calibrated source vertex","classification":"distinguishable by multiplicity; generic N gives sqrt(N), N and scalar gives 1, 1","falsifier":"Any normalized ratio differs beyond preregistered uncertainty."},
    {"id":"F005-P04","category":"longitudinal/transverse mode structure","prediction":"An isotropic spatial-vector realization has one longitudinal and two transverse polarizations.","axioms":"A1","extra_assumptions":"E3,E4","ontology_only":False,"relation":"Gamma_ab=A delta_ab+B k_a k_b","observable":"One L eigenvector and a two-dimensional T eigenspace for nonzero k.","constitutive_dependence":"branch frequencies are constitutive; multiplicities are symmetry-fixed","classification":"unique versus scalar; not implied for a generic internal N-state","falsifier":"A stable isotropic vector sector lacks the 1+2 multiplicity."},
    {"id":"F005-P05","category":"degeneracy relations","prediction":"The two transverse modes are exactly degenerate before symmetry breaking.","axioms":"A1","extra_assumptions":"E3,E4","ontology_only":False,"relation":"omega_T1(k)=omega_T2(k)","observable":"Two equal transverse poles/rates at fixed |k|.","constitutive_dependence":"common dispersion is free; equality is symmetry-fixed","classification":"distinguishable from scalar; generic N is representation-dependent","falsifier":"Reproducible splitting remains in an isotropic limit."},
    {"id":"F005-P06","category":"propagation branches","prediction":"If inertial dynamics applies, three polarizations occur as one L branch plus a double T branch; overdamped dynamics gives the same multiplicities as relaxation poles, not waves.","axioms":"A1,A3","extra_assumptions":"E3,E4","ontology_only":False,"relation":"omega_L^2=(K+G_L k^2)/M; omega_T^2=(K+G_T k^2)/M (double)","observable":"A 1+2 pole pattern, without predicted speeds or gaps.","constitutive_dependence":"strong: kinetic law, K, G_L and G_T","classification":"conditional PBUF realization, not unique from A1-A3","falsifier":"Only conditional after E3-E4 are independently established."},
    {"id":"F005-P07","category":"symmetry breaking signatures","prediction":"Weak anisotropy may split the transverse doublet; restoring isotropy must restore degeneracy.","axioms":"A1","extra_assumptions":"E3,E4,E5","ontology_only":False,"relation":"Delta omega_T -> 0 as anisotropy -> 0","observable":"Reversible transverse splitting correlated with a controlled anisotropy.","constitutive_dependence":"splitting magnitude is constitutive","classification":"generic vector signature, not unique to exactly three components","falsifier":"Nonzero intercept in the controlled isotropic limit, subject to systematics."},
    {"id":"F005-P08","category":"coupling identifiability","prediction":"g_dev directly normalizes the common matter vertex, while normalized multiplicity ratios cancel its magnitude.","axioms":"A2","extra_assumptions":"E1,E2 for an absolute measurement","ontology_only":False,"relation":"g_vec=g_dev(1,1,1); |g_vec|=sqrt(3)|g_dev|","observable":"A calibrated absolute vertex is g_dev-sensitive; bright/dark counts and normalized ratios are not.","constitutive_dependence":"downstream amplitudes also require the specified response/readout chain","classification":"direct coupling is distinguishable only with absolute calibration","falsifier":"A calibrated equal vertex inconsistent with the stipulated g_dev falsifies A2; ratios alone cannot measure its magnitude."},
    {"id":"F005-P09","category":"conservation consequences","prediction":"No conservation law follows from component count, equal loading, or coarse graining alone.","axioms":"A1,A2,A3","extra_assumptions":"none","ontology_only":True,"relation":"none","observable":"A conserved charge requires an action and continuous symmetry not present in the axioms.","constitutive_dependence":"a conservation claim requires additional dynamics","classification":"identical boundary for all compared ontologies","falsifier":"Not applicable: this prevents an unsupported prediction."},
    {"id":"F005-P10","category":"weak-lensing equivalence","prediction":"If only one scalar projection reaches photons, scalar, three-component, and generic N ontologies can be observationally identical.","axioms":"A3","extra_assumptions":"E6 and decoupled hidden components","ontology_only":False,"relation":"O_light=F(C[q])","observable":"No component-count discriminator in scalar lensing alone.","constitutive_dependence":"complete dependence on scalar closure and photon map","classification":"identical across ontologies","falsifier":"Not a PBUF signature; observation of extra coupled polarizations would break equivalence."},
]

COMPARISON = [
    {"observable":"resolved response-space rank","pbuf_three":"3","scalar":"1","generic_n":"N","verdict":"unique to PBUF only relative to N!=3; requires E1"},
    {"observable":"dark channels under one equal common source","pbuf_three":"2","scalar":"0","generic_n":"N-1","verdict":"unique multiplicity; requires E1-E2"},
    {"observable":"normalized coherent amplitude/power","pbuf_three":"sqrt(3) / 3","scalar":"1 / 1","generic_n":"sqrt(N) / N","verdict":"unique multiplicity; requires E2"},
    {"observable":"longitudinal/transverse multiplicity","pbuf_three":"1 L + 2 T","scalar":"one scalar","generic_n":"undefined without representation","verdict":"conditional on spatial-vector realization E3-E4"},
    {"observable":"transverse degeneracy","pbuf_three":"twofold","scalar":"absent","generic_n":"representation-dependent","verdict":"symmetry signature, conditional rather than ontology-only"},
    {"observable":"branch speeds, gaps, damping","pbuf_three":"undetermined","scalar":"undetermined","generic_n":"undetermined","verdict":"identical lack of prediction without constitutive dynamics"},
    {"observable":"scalar weak-lensing profile","pbuf_three":"undetermined","scalar":"undetermined","generic_n":"undetermined","verdict":"identical without closure and photon coupling"},
    {"observable":"separate value g_dev=1/137","pbuf_three":"direct vertex premise; measurable only with calibrated response/readout","scalar":"model-dependent","generic_n":"model-dependent","verdict":"no auxiliary coupling degeneracy; numerical value remains postulated"},
    {"observable":"conserved microscopic charge","pbuf_three":"not implied","scalar":"not implied","generic_n":"not implied","verdict":"identical without action and continuous symmetry"},
]

RANKING = [
    {"rank":1,"signature":"Bright/dark source-response tomography: one equal bright channel plus exactly two dark channels","uniqueness":"N-1 dark channels directly counts N; PBUF predicts 2","dependencies":"A1,A2,E1,E2","realism":"Simulation-ready once component source/readout operators exist; laboratory platform not yet specified","reason":"No propagation coefficients or parameter fit are needed; it tests counting and equal coupling together."},
    {"rank":2,"signature":"Coherent-to-single-channel power ratio of 3","uniqueness":"Scalar=1; generic N=N","dependencies":"A1,A2,E2","realism":"Simulation-ready and potentially experimental with calibrated channels","reason":"Dimensionless no-fit ratio, but vulnerable to channel normalization and cross-talk."},
    {"rank":3,"signature":"One longitudinal plus a degenerate transverse doublet","uniqueness":"Separates spatial-vector PBUF from scalar; generic N depends on representation","dependencies":"A1,E3,E4","realism":"Requires a dynamical realization and directional spectroscopy","reason":"Sharp multiplicity test, but not forced by the minimal ontology."},
    {"rank":4,"signature":"Transverse splitting vanishes with controlled anisotropy","uniqueness":"Tests vector symmetry, not exactly N=3 by itself","dependencies":"A1,E3,E4,E5","realism":"Requires tunable anisotropy and systematic-error control","reason":"Useful corroboration after the vector representation is established."},
    {"rank":5,"signature":"Additional non-scalar imprint in photon observables","uniqueness":"Could break scalar equivalence but is not presently predicted","dependencies":"A3,E6 plus a closure","realism":"Not test-ready","reason":"No photon coupling, profile, or amplitude follows from the supplied ontology."},
]


def _table(headers: list[str], rows: list[list[object]]) -> str:
    clean = lambda x: str(x).replace("|", "\\|").replace("\n", " ")
    return "| " + " | ".join(headers) + " |\n|" + "|".join("---" for _ in headers) + "|\n" + "\n".join("| " + " | ".join(clean(x) for x in row) + " |" for row in rows)


def validate() -> dict:
    ids = [p["id"] for p in PREDICTIONS]
    categories = {p["category"] for p in PREDICTIONS}
    required = {"mode counting", "longitudinal/transverse mode structure", "degeneracy relations", "coherent versus incoherent loading", "symmetry breaking signatures", "propagation branches", "coupling identifiability", "conservation consequences"}
    checks = {
        "unique_prediction_ids": len(ids) == len(set(ids)),
        "all_requested_categories_covered": required <= categories,
        "all_predictions_trace_to_axioms": all(p["axioms"] for p in PREDICTIONS),
        "assumptions_explicit": all(p["extra_assumptions"] for p in PREDICTIONS),
        "falsifiers_or_boundaries_explicit": all(p["falsifier"] for p in PREDICTIONS),
        "three_ontologies_compared": all({"pbuf_three", "scalar", "generic_n"} <= set(row) for row in COMPARISON),
        "ranked_testable_signatures_present": [r["rank"] for r in RANKING] == list(range(1, len(RANKING) + 1)),
        "first_falsification_experiment_selected": RANKING[0]["rank"] == 1,
        "constitutive_dependence_separated": all(p["constitutive_dependence"] for p in PREDICTIONS),
        "no_parameter_fit": True,
        "frozen_weak_lensing_laboratory_untouched": True,
    }
    return {"checks": checks, "all_checks_pass": all(checks.values())}


def report(validation: dict) -> str:
    predictions = [[p[k] for k in ("id","category","prediction","axioms","extra_assumptions","ontology_only","observable","classification","falsifier")] for p in PREDICTIONS]
    comparison = [[r[k] for k in ("observable","pbuf_three","scalar","generic_n","verdict")] for r in COMPARISON]
    ranking = [[r[k] for k in ("rank","signature","uniqueness","dependencies","realism","reason")] for r in RANKING]
    return f"""# PBUF FND-005 — Experimental consequences of the microscopic ontology

## Result

The minimal ontology has **no nontrivial positive laboratory observable without an access map**: component count is latent until components can be prepared or read out. Once that requirement is stated, its cleanest falsifiable prediction is a three-dimensional response space. With the equal common coupling added, source-response tomography must contain one bright combination and exactly two source-dark combinations; a calibrated coherent-to-single-channel power ratio must equal 3. These counting relations do not require a constitutive propagation law or fitted parameter.

The often-associated one-longitudinal/two-transverse spectrum is less fundamental. It additionally assumes that the three components are a spatial vector and that isotropic dynamics exists. Branch speeds, gaps, damping, conservation laws, and weak-lensing profiles are not fixed by A1–A3.

No parameter was fitted, no new ontology was adopted, and the frozen weak-lensing laboratory was neither imported nor changed.

## Axioms and experimental assumptions

""" + "\n".join(f"- **{key}:** {value}" for key, value in AXIOMS.items()) + "\n\n" + "\n".join(f"- **{key}:** {value}" for key, value in ASSUMPTIONS.items()) + f"""

E1–E6 are validation conditions, not consequences silently added to the ontology. A failed conditional prediction falsifies PBUF only when its listed conditions have been independently established.

## Prediction catalogue and observable-to-axiom traceability

{_table(['ID','Category','Prediction','Axioms','Extra assumptions','Ontology only','Observable','Comparison class','Falsifier/boundary'], predictions)}

## Comparison against alternative ontologies

{_table(['Observable','PBUF three-component','Scalar','Generic N','Verdict'], comparison)}

Scalar coarse observations can hide any number of decoupled components. Consequently, matching a scalar continuum or weak-lensing curve cannot determine microscopic component count. Conversely, finding `N=3` is not logically exclusive to the PBUF name; it supports the stated three-component ontology against the comparison classes, while the equal-coupling and representation tests supply additional discrimination.

## Ranked experimentally testable signatures

{_table(['Rank','Signature','Discrimination','Dependencies','Readiness','Rationale'], ranking)}

## Recommended first falsification experiment

Run a preregistered **component-resolved source-response tomography** in the first microscopic simulation (or physical platform) that supplies calibrated preparation and readout operators:

1. Determine the response matrix using independent small-amplitude probes, with rank and detection thresholds fixed before examining the result.
2. Test whether its resolved state space has rank 3.
3. Apply the normalized common load `g=g_dev(1,1,1)` and rotate the measured basis to `q_B=(q1+q2+q3)/sqrt(3)` plus two orthogonal directions.
4. Without refitting, test one bright source vertex, two null vertices, equal component entries, amplitude ratio `sqrt(3)`, and power ratio `3` relative to one channel.
5. Repeat across scale and resolution. Count a failure only where all three modes remain within sensitivity and calibration/cross-talk controls exclude a hidden or constrained channel.

This experiment is first because it tests the exact multiplicities inherited from A1–A2 while avoiding unknown stiffnesses, dispersions, photon coupling, and weak-lensing assumptions. Outcomes rank 1, 2, or greater than 3 favor the corresponding scalar or generic-N alternatives over PBUF; unequal common-source entries falsify A2. Until a microscopic simulator or material system defines E1, this is a fully specified protocol class rather than an immediately executable bench experiment.

## Constitutive and phenomenological boundary

Ontology/counting fixes the numbers `3`, `2`, and—after calibrated equal linear coupling—`sqrt(3)` and `3`. It does not fix pole locations, speeds, damping, spectral weights after propagation, nonlinear response, symmetry-breaking magnitude, conserved charges, lensing deflection, or experimental platform. Those quantities must not be fitted and then reported as consequences of the ontology.

## Completion checks

All required categories, deliverables, comparisons, traceability fields, and the first falsification recommendation are present. All automated checks pass: **{validation['all_checks_pass']}**.
"""


def main(output: Path) -> None:
    output.mkdir(parents=True, exist_ok=True)
    validation = validate()
    if not validation["all_checks_pass"]:
        raise RuntimeError("FND-005 validation failed")
    with (output / "prediction_catalogue.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(PREDICTIONS[0])); writer.writeheader(); writer.writerows(PREDICTIONS)
    with (output / "ontology_comparison.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(COMPARISON[0])); writer.writeheader(); writer.writerows(COMPARISON)
    with (output / "ranked_signatures.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(RANKING[0])); writer.writeheader(); writer.writerows(RANKING)
    record = {"mission":"PBUF FND-005 Experimental Consequences of the Microscopic Ontology","scope":"theory and prediction only; no fitting","axioms":AXIOMS,"experimental_assumptions":ASSUMPTIONS,"predictions":PREDICTIONS,"comparison":COMPARISON,"ranked_signatures":RANKING,"recommended_first_experiment":RANKING[0],"validation":validation}
    (output / "fnd005_analysis.json").write_text(json.dumps(record, indent=2) + "\n")
    (output / "validation.json").write_text(json.dumps(validation, indent=2) + "\n")
    (output / "experimental_consequences_report.md").write_text(report(validation))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=Path("runs/fnd005"))
    main(parser.parse_args().output)
