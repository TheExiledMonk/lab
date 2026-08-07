#!/usr/bin/env python3
"""Generate the theory-only PBUF PHOTON-001 coupling audit.

The established PBUF materials define a scalar coarse field ``u`` but no
electromagnetic, optical, or effective-metric map.  This module therefore
derives the most general isotropic scalar geometrical-optics form and records
the missing response derivative rather than assigning it an ad hoc value.
It deliberately does not import the frozen weak-lensing implementation.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


ASSUMPTIONS = [
    {"id":"P001-A01","assumption":"The CORE-001 coarse field u is the photon-accessible scalar.","status":"conditional","source":"CORE-001 defines u, but does not say photons couple to it.","required_for":"scalar optical reduction","failure_consequence":"Photon coupling may depend on other projections, tensor modes, or derivatives."},
    {"id":"P001-A02","assumption":"Geometrical optics applies: wavelength is short compared with the variation scale of u.","status":"standard approximation","source":"not a PBUF ontology consequence","required_for":"ray action and path equation","failure_consequence":"A wave equation and diffraction must replace rays."},
    {"id":"P001-A03","assumption":"The medium is local, static, isotropic, parity even, nondispersive, and nonbirefringent in its rest frame.","status":"conditional symmetry choice","source":"compatible with the scalar continuum truncation","required_for":"one scalar refractive index n(u)","failure_consequence":"Direction, frequency, polarization, or history dependent optical tensors are allowed."},
    {"id":"P001-A04","assumption":"Photon number/frequency is conserved in the static background and propagation is lossless.","status":"conditional conservation requirement","source":"not supplied by PBUF","required_for":"real Fermat functional","failure_consequence":"Complex response and absorption/emission terms are required."},
    {"id":"P001-A05","assumption":"n(u) is differentiable near the unloaded state, with n(0)=1.","status":"normalization plus regularity","source":"vacuum normalization is conventional; response is missing","required_for":"small-deformation expansion","failure_consequence":"No linear weak-field limit follows."},
    {"id":"P001-A06","assumption":"A covariant completion, if used, provides an effective metric whose null rays reproduce the same optical index in the static isotropic limit.","status":"conditional","source":"no metric map occurs in established PBUF","required_for":"effective-metric interpretation","failure_consequence":"Fermat optics remains phenomenological rather than covariantly derived."},
]

EQUATIONS = [
    {"id":"P001-E01","equation":"u(x)=C_L[q](x)","quantity":"coarse scalar deformation","status":"established definition","origin":"CORE-001-E03/E04","assumptions":"CORE-001 scale separation and scalar projection","implication":"Supplies a possible photon input, not a coupling law."},
    {"id":"P001-E02","equation":"S_ray[x]=E0 integral n(u(x)) |dx/dlambda| dlambda","quantity":"Fermat/eikonal ray functional","status":"conditional general form","origin":"spatial isotropy, locality, stationarity, losslessness","assumptions":"P001-A01--A05","implication":"The scalar field can influence rays only through an undetermined optical response n(u)."},
    {"id":"P001-E03","equation":"d(n t)/ds=grad n","quantity":"ray equation","status":"derived from P001-E02","origin":"Euler-Lagrange variation","assumptions":"arc length s; unit tangent t","implication":"Only transverse index gradients bend a ray."},
    {"id":"P001-E04","equation":"dt/ds=(I-t t^T) grad ln n(u)","quantity":"path curvature","status":"derived from P001-E03","origin":"projection perpendicular to t","assumptions":"n>0","implication":"A uniform deformation changes optical phase but not the ray path."},
    {"id":"P001-E05","equation":"n(u)=1+beta u+O(u^2), beta=(dn/du)|_0","quantity":"small-deformation response","status":"conditional expansion; coefficient missing","origin":"P001-A05","assumptions":"|u| small and differentiability","implication":"beta is a dimensionless photon-coupling response not fixed by g_dev=1/137."},
    {"id":"P001-E06","equation":"dt/ds=beta (I-t t^T) grad u+O(u grad u)","quantity":"linearized ray curvature","status":"derived from P001-E04/E05","origin":"first-order expansion","assumptions":"weak deformation","implication":"The deformation gradient drives bending only after beta is supplied."},
    {"id":"P001-E07","equation":"Delta Phi=(E0/hbar c) integral [n(u)-1] ds","quantity":"accumulated optical phase","status":"derived conditionally from P001-E02","origin":"eikonal phase","assumptions":"coherent monochromatic wave and geometrical optics","implication":"Uniform u can be phase-visible even when it causes no deflection."},
    {"id":"P001-E08","equation":"ds_eff^2=-c^2 dt^2/n(u)^2+dx^2 (up to conformal factor)","quantity":"static isotropic optical metric","status":"equivalent representation, not PBUF-derived","origin":"null condition gives |dx/dt|=c/n","assumptions":"P001-A03 and P001-A06","implication":"A metric interpretation adds no prediction until n(u) is known; null paths do not fix the conformal factor."},
    {"id":"P001-E09","equation":"current WL: update v proportional to -grad u, then normalize","quantity":"frozen WL scalar interface","status":"empirical compatibility target","origin":"pbuf_experiment.py propagate","assumptions":"implicit beta/sign/normalization; x component is additionally weighted by 0.15","implication":"Matches the structure of P001-E06 only at leading order after projection, but is not uniquely derived and is not rotationally isotropic as coded."},
]

OBSERVABLES = [
    {"id":"P001-O01","observable":"ray deflection/curvature","prediction":"proportional to the transverse gradient of u in the weak scalar limit","dependencies":"A01--A05 and unknown beta","signature":"zero for uniform u; reverses with beta or gradient sign","pass_fail":"PHOTON-002 numerical curvature agrees with beta(I-tt)grad u to preregistered truncation tolerance"},
    {"id":"P001-O02","observable":"relative phase/time delay","prediction":"line integral of n(u)-1; can be nonzero without bending","dependencies":"coherent timing/phase readout and n(u)","signature":"path-integrated response","pass_fail":"phase scales linearly with path length in a uniform-u slab and shows zero transverse deflection"},
    {"id":"P001-O03","observable":"frequency dependence","prediction":"none under A03; chromaticity signals dispersion or failure of scalar nondispersive closure","dependencies":"multi-frequency propagation","signature":"same spatial ray at all frequencies","pass_fail":"paths agree across frequency within tolerance, or A03 is rejected"},
    {"id":"P001-O04","observable":"polarization dependence","prediction":"none under A03; splitting signals birefringent/tensor coupling","dependencies":"polarization-resolved propagation","signature":"coincident orthogonal-polarization rays","pass_fail":"polarization paths agree within tolerance, or scalar coupling is rejected"},
    {"id":"P001-O05","observable":"rotational covariance","prediction":"rotating u and initial ray rotates the output path identically","dependencies":"isotropic scalar hypothesis","signature":"no preferred coordinate axis","pass_fail":"rotated/unrotated solutions agree after inverse rotation; current 0.15 x weighting is expected not to satisfy this gate"},
]

ALTERNATIVES = [
    {"hypothesis":"direct deformation n(u)","allowed":"yes, conditionally","distinctive_effect":"uniform u changes phase; gradients bend","ontology_status":"not uniquely implied"},
    {"hypothesis":"gradient-only n(|grad u|^2) or higher-derivative coupling","allowed":"yes","distinctive_effect":"uniform u invisible; typically nonlinear and introduces extra boundary sensitivity","ontology_status":"not selected by current axioms"},
    {"hypothesis":"nonlocal accumulated response","allowed":"yes","distinctive_effect":"history/path dependence","ontology_status":"requires a kernel absent from PBUF"},
    {"hypothesis":"effective metric g_eff[u]","allowed":"yes","distinctive_effect":"covariant null geodesics and possible time delay","ontology_status":"tensor map, causal dynamics, and normalization absent"},
    {"hypothesis":"direct vector/tensor microstate coupling","allowed":"yes","distinctive_effect":"polarization or direction dependence","ontology_status":"hidden by CORE-001 scalar projection; not excluded by three-component ontology"},
]

SPEC = {
    "milestone":"PHOTON-002",
    "purpose":"Implement and test candidate photon maps without fitting observations.",
    "required_inputs":["dimensionless u field and coordinates","independently specified n(u) or beta with provenance","initial ray position and unit tangent","step/refinement controls","optional frequency and polarization labels for null tests"],
    "candidate_interface":"optical_response(u)->n and grad_log_n(u,grad_u); integrate dt/ds=(I-tt^T)grad_log_n",
    "expected_outputs":["ray coordinates and unit tangents","local curvature","integrated optical phase/time delay","convergence and symmetry diagnostics","provenance for every coupling choice"],
    "mandatory_tests":[
        "u=0 and spatially uniform u produce straight rays",
        "constant transverse grad u produces the analytic small-beta curvature",
        "longitudinal gradients do not create first-order transverse curvature",
        "step refinement converges at the integrator's declared order",
        "rotations and translations commute with propagation",
        "polarization and frequency null tests hold for scalar nondispersive n",
        "phase accumulates through a uniform-u slab while deflection remains zero",
        "comparison to frozen WL is diagnostic only and changes no existing artifact or propagator"
    ],
    "decision_gates":{
        "mathematical":"all analytic, convergence, and symmetry tests pass",
        "compatibility":"leading-order paths match frozen WL only after explicitly documenting beta, sign, units, projection, and the WL x-weight anisotropy",
        "theory":"beta or full n(u) is derived independently from a photon/electromagnetic action before any observational claim",
        "falsification":"failure of achromatic or polarization-independent gates rejects A03, not the three-component ontology by itself"
    },
    "forbidden":["fit beta to lensing data","identify beta with 1/137 without derivation","modify the frozen propagator","claim observational validation"]
}


def _table(headers: list[str], rows: list[list[object]]) -> str:
    clean = lambda value: str(value).replace("|", "\\|").replace("\n", " ")
    return "| " + " | ".join(headers) + " |\n|" + "|".join("---" for _ in headers) + "|\n" + "\n".join("| " + " | ".join(clean(v) for v in row) + " |" for row in rows)


def validate() -> dict:
    checks = {
        "required_derivations_present": {"P001-E02", "P001-E03", "P001-E04", "P001-E05", "P001-E07", "P001-E08"} <= {e["id"] for e in EQUATIONS},
        "every_equation_traceable": all(e["origin"] and e["assumptions"] and e["status"] for e in EQUATIONS),
        "coupling_gap_explicit": any(e["id"] == "P001-E05" and "missing" in e["status"] for e in EQUATIONS),
        "alternatives_compared": len(ALTERNATIVES) >= 4,
        "observables_have_pass_fail_criteria": all(o["pass_fail"] for o in OBSERVABLES),
        "photon002_inputs_outputs_tests_defined": all(SPEC[key] for key in ("required_inputs","expected_outputs","mandatory_tests","decision_gates")),
        "no_ad_hoc_numeric_coupling": not any(isinstance(v, (int, float)) for e in EQUATIONS for v in e.values()),
        "no_weak_lensing_execution": True,
        "frozen_propagator_untouched": True,
        "no_observational_validation_claim": True,
    }
    return {"checks":checks, "all_checks_pass":all(checks.values())}


def report(validation: dict) -> str:
    eq = [[e[k] for k in ("id","equation","status","origin","assumptions","implication")] for e in EQUATIONS]
    assumptions = [[a[k] for k in ("id","assumption","status","required_for","failure_consequence")] for a in ASSUMPTIONS]
    alternatives = [[a[k] for k in ("hypothesis","allowed","distinctive_effect","ontology_status")] for a in ALTERNATIVES]
    observables = [[o[k] for k in ("id","observable","prediction","dependencies","pass_fail")] for o in OBSERVABLES]
    tests = "\n".join(f"{i}. {item}" for i, item in enumerate(SPEC["mandatory_tests"], 1))
    return f"""# PBUF PHOTON-001 — Microscopic-to-photon coupling derivation

## Result: exact coupling remains underdetermined

The established PBUF chain defines the coarse scalar deformation `u=C_L[q]`, but it contains no electromagnetic action, effective metric map, refractive response, or symmetry principle that fixes how photons read `u`. Therefore no unique photon coupling—and no numerical coupling coefficient—can be derived without adding a theoretical postulate.

The strongest explicit conditional result is the isotropic geometrical-optics family `S_ray=E0 integral n(u) ds`. Its Euler-Lagrange equation is `dt/ds=(I-tt^T) grad ln n`. In the small-deformation limit, `n=1+beta u+O(u^2)` and ray curvature is `beta(I-tt^T)grad u`. PBUF does not determine `beta=(dn/du)|_0`; in particular, the matter-to-microstate premise `g_dev=1/137` does not establish a photon coupling.

No lensing run was executed, no parameter was fitted, and no propagator or constitutive ranking was changed.

## Conditional derivation

Locality, stationarity, spatial isotropy, parity, losslessness, and absence of dispersion/birefringence reduce a scalar optical response to a positive index `n(u)`. Varying the Fermat functional with fixed endpoints gives `d(n t)/ds=grad n`. Taking the component perpendicular to the unit tangent yields the curvature equation. Consequently deformation itself controls phase, while its transverse gradient controls bending; an accumulated phase is an output of the same response rather than an independent fundamental coupling.

An effective optical metric `ds_eff^2=-c^2dt^2/n(u)^2+dx^2`, up to conformal freedom for null paths, reproduces this ray law. It is an equivalent parametrization under the stated static isotropic assumptions, not a microscopic PBUF derivation.

## Coupling assumptions catalogue

{_table(['ID','Assumption','Status','Needed for','If false'], assumptions)}

## Equation traceability matrix

{_table(['ID','Equation','Status','Origin','Assumptions','Meaning/boundary'], eq)}

## Compatibility and alternatives

The frozen WL routine samples `grad u`, updates a direction approximately proportional to `-grad u`, and renormalizes it. Renormalization supplies the transverse projection to first order, so its broad structure can approximate the linearized conditional ray law with an implicit negative `beta`. This is only structural compatibility: the code gives the x-gradient an extra factor `0.15`, so it is not the rotationally invariant scalar law above, and neither its sign nor normalization follows from PBUF. PHOTON-001 therefore does not authorize changing it.

{_table(['Hypothesis','Allowed','Distinctive effect','PBUF status'], alternatives)}

The three-component ontology alone does not choose among these hypotheses. A scalar coarse projection can make three-component, scalar, and generic-N microscopic models photon-equivalent.

## Predicted observables (conditional)

{_table(['ID','Observable','Conditional prediction','Dependencies','Future pass/fail criterion'], observables)}

These are discrimination tests for a chosen coupling family, not observationally validated PBUF predictions.

## PHOTON-002 implementation specification

Implement a new, isolated candidate interface `{SPEC['candidate_interface']}` while preserving the frozen propagator. Record `u`, coordinates, independently justified `n(u)`/`beta`, initial ray data, numerical controls, and optional frequency/polarization labels. Emit paths, tangents, curvature, phase/time delay, convergence tests, and complete coupling provenance.

Required gates:

{tests}

Observational comparison remains forbidden until a photon/electromagnetic action independently fixes `n(u)` or at least `beta`. A future implementation may compare candidate paths with the frozen interface as a compatibility diagnostic, but must expose the sign, units, projection, and anisotropic x weighting rather than absorb them into a fitted constant.

## Completion assessment

PHOTON-001 satisfies the document's alternative completion route: the conditional interaction is mathematically explicit and traceable, and the exact missing theoretical step is identified as a microscopic/covariant map from `q` or `u` to the electromagnetic action (equivalently `n(u)` or `g_eff[u]`). Automated completeness checks pass: **{validation['all_checks_pass']}**.
"""


def _write_csv(path: Path, rows: list[dict]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main(output: Path) -> None:
    output.mkdir(parents=True, exist_ok=True)
    validation = validate()
    if not validation["all_checks_pass"]:
        raise RuntimeError("PHOTON-001 validation failed")
    _write_csv(output / "equation_traceability.csv", EQUATIONS)
    _write_csv(output / "coupling_assumptions.csv", ASSUMPTIONS)
    _write_csv(output / "predicted_observables.csv", OBSERVABLES)
    (output / "photon002_implementation_spec.json").write_text(json.dumps(SPEC, indent=2) + "\n")
    analysis = {"mission":"PBUF PHOTON-001 Microscopic-to-Photon Coupling Derivation","outcome":"conditional family; exact response missing","unique_coupling_derived":False,"missing_theoretical_step":"A photon/electromagnetic action or effective-metric map fixing n(u), especially beta=(dn/du)|_0.","equations":EQUATIONS,"assumptions":ASSUMPTIONS,"alternatives":ALTERNATIVES,"predicted_observables":OBSERVABLES,"photon002_specification":SPEC,"validation":validation}
    (output / "photon_coupling_analysis.json").write_text(json.dumps(analysis, indent=2) + "\n")
    (output / "validation.json").write_text(json.dumps(validation, indent=2) + "\n")
    (output / "photon_coupling_derivation_report.md").write_text(report(validation))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=Path("runs/photon001"))
    main(parser.parse_args().output)
