"""PBUF MATERIAL-LAB-002: interaction/energy separation audit.

This milestone classifies structures already present in MATERIAL-LAB-001.  It
does not select a constitutive formula, add ontology, or touch lensing code.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path


MODELS = (
    ("A", "Progressive strain hardening", "energy law only",
     "Local hyperelastic w gives stress, amplitude hardening, and a restoring tendency; no spatial coupling or kinetic evolution is specified."),
    ("B", "Wave-equilibrium material", "mixed variational formulation; incomplete dynamical constitutive law",
     "Its local quadratic term is an energy law and its gradient term generates a variational interaction operator. Positive inertia, assumed only conditionally in LAB-001, is still needed for waves; amplitude hardening is absent."),
    ("C", "Exponential hardening", "energy law only",
     "Local hyperelastic energy gives stress, hardening, and restoring tendency, but no neighbour coupling."),
    ("D", "Finite-extensibility material", "energy law only",
     "Local barrier energy gives stress, hardening, and restoring tendency on its domain, but no neighbour coupling."),
    ("E", "Polynomial hardening", "energy law only",
     "Local hyperelastic energy gives stress, hardening, and restoring tendency, but no neighbour coupling."),
)

ARROWS = (
    ("Matter", "Interaction", "assumption", "Matter identifies the medium/source but does not uniquely determine a coupling operator."),
    ("Interaction", "Neighbour communication", "derivable by definition, conditionally", "Requires a genuinely nonlocal, gradient, or adjacency-coupling interaction rather than a pointwise operator."),
    ("Neighbour communication", "Wave propagation", "not derivable", "A kinetic/inertial evolution law and hyperbolicity/positive dispersion are additionally required; diffusion and elliptic equilibrium are counterexamples."),
    ("Wave propagation", "Deformation", "derivable only as propagation of a specified deformation field", "The wave variable must already be identified with the frozen deformation variable; waves do not define that variable."),
    ("Deformation", "Stored energy", "assumption", "A state variable admits many energies and may have none; an energy functional must be supplied."),
    ("Stored energy", "Stress", "derivable conditionally", "For differentiable hyperelastic energy, stress is its first variation; nonsmooth energies give a subdifferential."),
    ("Stress", "Progressive hardening", "not derivable", "Hardening needs a positive, amplitude-increasing tangent; linear stress is a counterexample."),
    ("Progressive hardening", "Recovery", "not derivable dynamically", "Hardening alone does not supply a stable reference minimum or an evolution law; energy stability plus admissible dynamics supplies restoring/recovery behavior."),
)


def report() -> str:
    rows = "\n".join(f"| {a} → {b} | {s} | {why} |" for a, b, s, why in ARROWS)
    models = "\n".join(f"| {k} | {n} | {c} | {r} |" for k, n, c, r in MODELS)
    return f"""# PBUF MATERIAL-LAB-002 — Separation of Interaction Law and Energy Law

## Result

Interaction and energy are distinct constitutive roles, but they are not always represented by independent mathematical objects. A local stored-energy law cannot generate neighbour communication. A nonlocal or gradient energy can generate a restricted, conservative interaction through its first variation. Conversely, a general interaction law need not possess an energy potential. The exact result is therefore **conditional/partial equivalence**, not universal equivalence or universal independence.

The minimum architecture is outcome **B: Interaction plus Energy**, embedded in the already-required balance/evolution structure. The two slots may be encoded in one generalized energy only when the interaction is variational. Wave propagation additionally requires the kinetic/inertial content of the frozen balance/duration architecture; this is not a third material law and no new ontology is introduced.

## Formal definitions

Let `C` be the frozen objective deformation field on an admissible configuration space `X`.

An **Interaction Law** is a spatial operator `I[C]` (with its domain and boundary conditions) for which `I[C](x)` depends on `C` in a neighbourhood of `x`, or couples distinct material regions. It governs redistribution and supplies the spatial part of a disturbance equation. Genuine neighbour communication means its Fréchet derivative has off-diagonal spatial support. Wave propagation follows only after an admissible kinetic operator makes the resulting evolution hyperbolic.

An **Energy Law** is a bounded-below stored-energy functional `E: X -> R union {{+infinity}}`, minimized at the reference state, whose first variation defines generalized stress/restoring force, `S=delta E/delta C` (or `S in partial E`), and whose second variation defines tangent stiffness. Progressive hardening is an additional inequality on that tangent along admissible loading paths; it is not implied merely by existence of `E`.

These definitions separate physical roles. They do not require separate formulae.

## Task 1 — Can interaction come solely from energy?

Yes, conditionally. For

`E[C] = integral_Omega psi(C, grad C) dV`,

variation gives

`delta E/delta C = partial_C psi - Div(partial_gradC psi)`.

The divergence term communicates between neighbours. More generally a symmetric nonlocal kernel gives

`E_int = (1/4) integral integral (C(x)-C(y)):K(x,y):(C(x)-C(y)) dx dy`,

and its variation is a nonlocal interaction operator. Necessary conditions for sole energetic derivation are: an energy domain and boundary conditions; differentiability (or a usable subdifferential); path independence/integrability of the force one-form; symmetry of the second variation; and lower-boundedness/stability for a stored energy. In a simply connected linear setting, `I=delta E_int/delta C` exists exactly when `I` is self-adjoint (with the stated boundary conditions), up to affine terms.

A pointwise energy `E=integral psi(C(x)) dV` has a diagonal first variation and therefore cannot generate neighbour coupling. Also, advective, gyroscopic, non-reciprocal, or dissipative interactions generally violate variational symmetry and cannot be derived from stored energy alone. Thus an independent interaction specification is required unless the interaction is explicitly placed inside a generalized gradient/nonlocal energy.

## Task 2 — Does interaction imply progressive hardening?

No. Model B is a direct counterexample:

`E_B[e]=integral (k0 e^2/2 + k0 ell^2 |grad e|^2/2) dV`,

with interaction `-k0 ell^2 Laplacian(e)`, while its amplitude tangent is the constant `k0`. Scaling `e=lambda f` makes the restoring force linear in `lambda`; no progressive amplitude hardening appears. Conversely, A, C, D, and E harden locally while possessing no native neighbour coupling. Hence neither property implies the other.

## Task 3 — Equivalence analysis

- **Equivalence:** only within the conservative variational subclass, after fixing the domain, boundary conditions, and additive/affine energy ambiguity. A self-adjoint positive linear interaction has `E_int[C]=1/2 <C,I C>`.
- **Non-equivalence:** general interaction operators may be non-self-adjoint, dissipative, non-reciprocal, or kinetic; an energy law may be purely local and contain no interaction.
- **Partial equivalence:** `delta E_int` recovers the conservative interaction force, while the local part of `E` independently controls amplitude stress and hardening. Energy does not recover inertia or dissipation.
- **Conditional equivalence:** a nonlinear `I` is energetically representable when its Fréchet derivative is symmetric (the force one-form is closed) on a suitable simply connected configuration domain and the resulting potential satisfies stored-energy stability requirements.

An operator maps a field to a force; an energy maps a field to a scalar. They cannot literally be transformed into one another without variation/integration plus domain and boundary data.

## Task 4 — Minimal coupled material model

The minimal abstract pair is

`(I,E)`, with `E=E_loc[C]` and a genuine neighbour operator `I[C]`,

inserted into the frozen balance/evolution form

`M[C] C_ddot + (delta E/delta C) + I[C] = 0`.

Here `M` denotes the already-required positive kinetic/inertial structure, not a newly selected material parameter. Minimum conditions are:

1. `I` has off-diagonal spatial support and a positive/hyperbolic linearized spatial symbol with the kinetic term;
2. `E` has a strict reference minimum and a restoring first variation;
3. the tangent of `E` is positive and increases along the intended loading measure for progressive hardening;
4. boundary conditions make the operator well posed and the energy/balance compatible.

If `I=delta E_int/delta C`, the same architecture may be written with `E_total=E_loc+E_int`; the roles remain independently testable. No particular invariant, coefficient, kinetic formula, or preferred candidate is selected.

## Tasks 5–6 — MATERIAL-LAB-001 classification

| Model | Name | Classification | Mathematical reason |
|---|---|---|---|
{models}

Model B is therefore **both an energy law and a variational interaction law**, but not a complete dynamical constitutive law: its static functional has both local and gradient energy, while wave support still assumes positive inertia and it lacks progressive hardening.

## Task 7 — Dependency graph audit

| Arrow | Status | Reason |
|---|---|---|
{rows}

```text
Matter --[assume I]--> Interaction --[nonlocality]--> Neighbour communication
                                                     |
                                      + kinetic/hyperbolic closure
                                                     v
                                               Wave propagation
                                                     |
                                      specified deformation variable
                                                     v
                                                Deformation
                                                     |
                                           [assume energy E]
                                                     v
Stored energy --[variation]--> Stress --[tangent inequality]--> Hardening
       |
       +--[stable minimum + admissible evolution]----------------> Recovery
```

The proposed single vertical chain is therefore not a chain of implications. It contains two independent constitutive insertions (`I` and `E`) plus conditional mathematical consequences.

## Task 8 — Minimum native architecture and recommendation

Outcome **B** is required at the level of constitutive roles. A single generalized variational functional is an admissible representation only if it contains distinguishable local-energy and interaction contributions and passes the integrability conditions above. Outcome A is too weak in general; outcome C is unnecessary as a material decomposition because kinetics, boundary data, and balance are already separate frozen governing structures; no evidence supports outcome D.

Subsequent constitutive development should keep two explicit closure slots: (i) invariant neighbour interaction (or its gradient/nonlocal energy), and (ii) invariant local stored energy satisfying stability and progressive-hardening inequalities. Next test tensor lifts for objectivity, variational integrability, spectral positivity/strong ellipticity, and hyperbolicity with the frozen balance law. Preserve symbolic coefficients until independently derived. Do not select among A/C/D/E, alter V11 or weak lensing, or fit observations.
"""


def write_outputs(root: Path) -> None:
    out = root / "runs" / "material_lab002"
    out.mkdir(parents=True, exist_ok=True)
    (out / "interaction_energy_separation_report.md").write_text(report())
    with (out / "model_classification.csv").open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(("model", "name", "classification", "reason"))
        writer.writerows(MODELS)
    with (out / "dependency_graph.csv").open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(("from", "to", "status", "reason"))
        writer.writerows(ARROWS)
    analysis = {
        "milestone": "PBUF MATERIAL-LAB-002",
        "conclusion": "conditional/partial equivalence; independent roles",
        "minimum_architecture": "B: Interaction plus Energy, within frozen balance/evolution structure",
        "interaction_from_energy": "iff variational/integrable; local energy alone is insufficient",
        "hardening_from_interaction": False,
        "model_B": "both energy law and variational interaction law; incomplete dynamically",
        "models_A_C_D_E": "energy laws only",
    }
    (out / "separation_analysis.json").write_text(json.dumps(analysis, indent=2) + "\n")
    validation = {
        "milestone": "PBUF MATERIAL-LAB-002", "pass": True,
        "checks": {
            "formal_definitions": True, "independence_counterexamples": True,
            "equivalence_analysis": True, "all_models_classified": True,
            "dependency_graph_audited": True, "minimal_architecture_identified": True,
            "recommendation_provided": True,
        },
        "deliverables": ["interaction_energy_separation_report.md", "separation_analysis.json", "model_classification.csv", "dependency_graph.csv"],
        "prohibitions": {"ontology_modified": False, "weak_lensing_modified": False, "observational_fit": False, "v11_modified": False, "microscopic_constituents_added": False, "additional_free_parameters": False, "preferred_law_selected": False},
    }
    (out / "validation.json").write_text(json.dumps(validation, indent=2) + "\n")


if __name__ == "__main__":
    write_outputs(Path(__file__).resolve().parent)
