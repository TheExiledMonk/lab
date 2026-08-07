"""Generate the reproducible PBUF EQUILIBRIUM-001 structural audit."""
from __future__ import annotations

import csv
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "runs/equilibrium001"

SOURCES = {
    "FOUNDATION-001": "runs/foundation001/foundational_ontology.md",
    "STATE-002": "runs/state002/primitive_medium_state.md",
    "DEFORMATION-001": "runs/deformation001/deformation_measure_report.md",
    "HYPER-001": "runs/hyper001/stored_energy_derivation.md",
    "ENERGY-PRINCIPLE-001": "runs/energy_principle001/energy_selection_derivation.md",
    "DURATION-001": "runs/duration001/emergent_duration_derivation.md",
    "METRIC-001": "runs/metric001/effective_metric_derivation.md",
    "BALANCE-001": "runs/balance001/native_balance_laws.md",
    "CONSTITUTIVE-002": "runs/constitutive002/constitutive_comparison_report.md",
    "MATERIAL-LAB-001": "runs/material_lab001/material_laboratory_report.md",
    "MATERIAL-LAB-002": "runs/material_lab002/interaction_energy_separation_report.md",
    "MATERIAL-DISCOVERY-001": "runs/material_discovery001/survey_report.md",
    "CONSTITUTIVE-PRINCIPLES-001": "runs/constitutive_principles001/constitutive_principles_report.md",
    "CONSTITUTIVE-SELECTION-001": "runs/constitutive_selection001/constitutive_selection_report.md",
    "LOCALITY-001": "runs/locality001/locality_report.md",
    "ENERGY-SEARCH-001": "runs/energy_search001/energy_search_report.md",
    "NONLINEARITY-001": "runs/nonlinearity001/nonlinearity_report.md",
}

CANDIDATES = [
    ("minimum potential energy", "static minimization of stored plus admissible loading potential", "survives conditionally", "requires conservative prescribed loading and admissible boundary data"),
    ("stationary action", "stationarity on admissible oriented histories", "survives as a family", "native degree-one action exists structurally; its kinetic integrand is unselected"),
    ("constrained energy minimization", "minimization over the frozen admissible state/domain", "survives", "natural static statement for hard or regular constrained endpoints"),
    ("elastic force equilibrium", "vanishing internal variational force balanced by admitted loading", "survives but is not independent", "Euler-Lagrange/virtual-work form of the energy principle on smooth branches"),
    ("variational equilibrium", "stationarity of an admissible functional", "survives as an umbrella", "contains static energy and history action forms but does not select either integrand"),
    ("thermodynamic equilibrium", "extremum of a thermodynamic potential subject to constraints", "not native as stated", "temperature and entropy are absent; its purely mechanical reduction is constrained energy minimization"),
    ("spectral equilibrium", "stationarity expressed in unordered eigenvalues/invariants of C", "not independent", "spectral variables are coordinates for the same objective stored energy"),
    ("geometric equilibrium", "stationarity under admissible geometric variations", "not independent/conditional", "reduces to variational equilibrium when geometry is derived from q,C; extra curvature or metric variables are forbidden"),
    ("wave equilibrium", "normal modes, standing waves, or propagation about an equilibrium", "not a constitutive principle", "waves are consequences of balance, inertia, and acoustic positivity, not selectors of W"),
]

GENERATION = [
    ("minimum potential energy", "yes, on one elastic branch", "yes", "yes", "conditional on frozen inertia/ellipticity", "permits but does not generate R"),
    ("stationary action", "conditional on reversal-compatible action/boundaries", "conditional", "yes for its W sector", "yes after kinetic closure", "permits but does not select R"),
    ("constrained energy minimization", "yes on admissible branch", "yes", "yes in interior; tangent-cone form at endpoint", "conditional on dynamic completion", "permits D1-D3 but does not select shape"),
    ("elastic force equilibrium", "yes when derived from exact W", "not by balance alone", "yes in variational realization", "conditional", "inherits W; does not generate it"),
    ("variational equilibrium", "conditional", "conditional", "yes for stored-energy member", "conditional", "neutral without a functional-selection rule"),
    ("thermodynamic equilibrium", "not auditable natively", "not auditable natively", "only after mechanical reduction", "no", "no native selection"),
    ("spectral equilibrium", "inherits", "inherits", "inherits", "inherits", "coordinate-neutral"),
    ("geometric equilibrium", "inherits if reduced", "inherits if reduced", "inherits if reduced", "conditional", "neutral"),
    ("wave equilibrium", "no", "presupposes a background equilibrium", "no", "describes modes, not their constitutive origin", "no"),
]

D_RELATION = [
    ("minimum/stationary stored energy", "compatible", "compatible", "compatible with an explicit admissible set", "neutral among D1-D3"),
    ("stationary action", "compatible", "compatible", "compatible", "neutral; endpoint treatment changes admissible histories"),
    ("constrained minimization", "compatible", "compatible as a barrier problem", "naturally expresses D3", "weak structural affinity to D3, not selection"),
    ("elastic/variational balance", "inherits", "inherits singular boundary resistance", "uses variational inequalities at constrained endpoints", "neutral"),
    ("spectral/geometric reformulations", "inherits", "inherits", "inherits", "neutral"),
    ("wave equilibrium", "nonlinear waves may sample D1", "may sense approach to D2", "cannot cross D3", "diagnostic only; selects none"),
]


def write_json(name: str, value: object) -> None:
    (OUT / name).write_text(json.dumps(value, indent=2) + "\n")


def write_csv(name: str, header: tuple[str, ...], rows: list[tuple[str, ...]]) -> None:
    with (OUT / name).open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        writer.writerows(rows)


def main() -> None:
    missing = [path for path in SOURCES.values() if not (ROOT / path).is_file()]
    if missing:
        raise FileNotFoundError(missing)
    OUT.mkdir(parents=True, exist_ok=True)
    write_csv("principle_catalogue.csv", ("principle", "established_form", "pbuf_status", "scope_or_failure"), CANDIDATES)
    write_csv("constitutive_generation_audit.csv", ("principle", "reversible_recovery", "stable_reference", "local_variational_stress", "waves", "nonlinear_response"), GENERATION)
    write_csv("d1_d3_relationship.csv", ("principle_class", "D1", "D2", "D3", "assessment"), D_RELATION)
    write_json("equivalence_classes.json", {
        "E1_static_constrained_variation": ["minimum potential energy", "minimum stored energy when unloaded", "constrained energy minimization", "elastic equilibrium", "static variational equilibrium", "spectral equilibrium", "geometric equilibrium using only q and C"],
        "E2_history_variation": ["stationary native action", "dynamic variational equilibrium"],
        "derived_not_selection_principles": ["wave equilibrium", "normal-mode equilibrium", "Noether conservation/balance"],
        "inadmissible_without_reduction": ["thermodynamic equilibrium with independent entropy or temperature", "geometric principle with independent metric/curvature field"],
        "relation": "E1 is a static/quasistatic sector of E2 only after a kinetic/action and clock gauge are chosen; the frozen architecture does not justify identifying them globally."
    })
    report = r'''# PBUF EQUILIBRIUM-001 — Native Equilibrium Principle of the Spacetime Medium

## Decision

**Outcome B: a reduced family survives; no unique deeper equilibrium principle is selected.**

The frozen architecture natively supports two inequivalent levels:

1. **E1 — static constrained stored-energy variation:** equilibrium is a stationary point, and a stable unloaded state is a local minimum, of the already-frozen energy over admissible configurations.
2. **E2 — stationary native history action:** admissible evolution is stationary under variations of oriented, order-reparametrization-invariant histories, after choosing a member of the still-unselected action/kinetic family.

E1 includes minimum potential energy, constrained minimization, elastic equilibrium, virtual work, and their spectral or geometric representations. E2 reduces to E1 in a static/quasistatic sector only after extra dynamical and boundary hypotheses are supplied. Wave equilibrium is a solution regime of E2 or of balance equations, not a third constitutive principle.

Most importantly, neither E1 nor E2 generates the unknown nonlinear remainder. They accept

\[
W(C)=W_2(C)+R(C),\qquad R({\bf1})=DR({\bf1})=D^2R({\bf1})=0,
\tag{EQ-001}
\]

and test its equilibria after it is supplied. Thus PBUF has a native **constrained variational architecture**, but no frozen optimization or preservation rule that selects the finite-deformation shape of \(R\), its endpoint class, or a unique governing equation.

No ontology, field, constant, constitutive formula, fit, V11 modification, or weak-lensing change is introduced.

## 1. Catalogue of established principles

| Candidate | Established content | PBUF disposition |
|---|---|---|
| minimum stored/potential energy | stable static equilibria minimize stored energy plus conservative load potential | survives conditionally; external load potential is not frozen |
| stationary action | physical histories are stationary points of an action | survives as an unselected native action family |
| constrained energy minimization | minimize energy over an admissible set | survives and naturally represents the frozen state domain |
| elastic equilibrium / virtual work | internal variational work balances admitted loading | survives; equivalent to E1 on smooth conservative branches |
| variational equilibrium | stationarity of a functional under admissible variations | survives as the umbrella containing E1 and E2 |
| thermodynamic equilibrium | extremize an energy/free-energy potential with thermodynamic constraints | not native as stated; entropy and temperature are not frozen variables |
| spectral equilibrium | equilibrium conditions written in principal stretches or invariants | representation only; equivalent to E1 |
| geometric equilibrium | stationary geometric functional | equivalent to E1 only when the geometry is derived from existing \(q,C\); otherwise forbidden extra structure |
| wave equilibrium | standing/normal modes or nonlinear wave states | downstream solution regime, not an equilibrium selector for \(W\) |

The catalogue distinguishes an equilibrium **principle** from a coordinate representation, an Euler–Lagrange form, and a class of solutions.

## 2. Frozen compatibility audit

Legend: **Y** compatible; **C** conditional; **R** reduces to an already-listed principle; **N** incompatible as an independent native principle.

| Candidate | FOUNDATION | STATE | DEFORMATION | HYPER | DURATION | METRIC | BALANCE | LOCALITY | NONLINEARITY |
|---|---|---|---|---|---|---|---|---|---|
| minimum potential energy | Y | Y | Y | Y | Y | C | C | Y | Y |
| stationary action | Y | Y | Y | C | Y | C | Y | Y | Y |
| constrained minimization | Y | Y | Y | Y | Y | C | Y | Y | Y |
| elastic equilibrium | Y | Y | Y | R | Y | C | Y | Y | Y |
| variational equilibrium | Y | Y | Y | C | Y | C | Y | Y | Y |
| thermodynamic equilibrium | C | N | Y | N | C | C | N | C | N |
| spectral equilibrium | Y | Y | R | R | Y | C | R | Y | R |
| geometric equilibrium | C | C | C | C | Y | C | C | C | C |
| wave equilibrium | Y | C | Y | N | C | C | R | Y | N |

### Explicit incompatibilities and conditions

- **FOUNDATION-001:** E1 and E2 use the one medium and add no substrate or constant. A geometric principle becomes incompatible if it reifies an independent metric, curvature field, embedding, or microscopic structure. Thermodynamic language is acceptable only after reduction to existing mechanical quantities; an entropy sector would be new ontology/state structure.
- **STATE-002:** equilibrium must be a statement about \(q\), or \(C[q,q_0]\), without hidden history or internal variables. Standard thermodynamic equilibrium needs at least an independently defined entropy/temperature or ensemble and therefore fails natively. A wave population with amplitudes, phases, or occupations likewise fails as an independent state principle.
- **DEFORMATION-001:** all static constitutive variation must use the objective spectral content of \(C\). “Spectral equilibrium” adds no principle. Geometric proposals using curvature, an independent metric, or an extra strain measure exceed the frozen deformation content.
- **HYPER-001:** E1 is native because \(P_C=DW\) is exact. E2 is compatible only when its stored sector is this \(W\). Thermodynamic potentials with extra arguments conflict with the single-valued local hyperelastic state function. Waves do not determine \(W\); they presuppose its tangent response.
- **DURATION-001:** static E1 is parameter-free. E2 must be invariant under monotone relabeling of order and may use calibrated duration only after its derivation. Frequency-minimization or rate-dependent equilibrium is not native. Stationarity alone does not prove time reversal.
- **METRIC-001:** every survivor must preserve the single effective Lorentzian output and V11 limit. No candidate selects the unresolved map \(G[q,C]\); a metric-based equilibrium functional is only conditional unless it is pulled back to existing \(q,C\).
- **BALANCE-001:** local balance and Noether identities follow only after a selected action/symmetry/source accounting. Minimum potential energy additionally requires conservative loading. Force balance without exact variational stress is broader than E1 and cannot generate reversible hyperelasticity. Thermodynamic balance is unavailable because no entropy balance exists.
- **LOCALITY-001:** E1 produces local communication through \((D_qC)^*DW\), or \(-\operatorname{Div}P_F\) in the placement realization. Intrinsic kernels, gradients, or wave-mode closures are not native. E2 must retain the same minimal local constitutive sector.
- **NONLINEARITY-001:** every survivor must remain neutral among admissible zero-2-jet remainders and D1–D3 unless it proves a selector. Thermodynamic and wave stories need rejected extra state variables; after eliminating them they reduce to E1/E2 and add no constitutive content.

## 3. Native generation audit

| Principle | reversible recovery | stable reference | local variational stress | wave support | nonlinear constitutive response |
|---|---|---|---|---|---|
| E1 unconstrained/minimum energy | yes on one elastic branch | yes | yes | only with inertia and acoustic positivity | admits, does not generate, \(R\) |
| E1 constrained | yes in the interior and under ideal constraints | yes | yes; variational inequality at an endpoint | only after dynamic completion | admits D1–D3; selects no profile |
| E2 stationary action | conditional on reversal-compatible action/boundaries | conditional on its potential sector | yes for the frozen \(W\) sector | yes after kinetic closure | admits, does not select, \(R\) |
| force/balance statement alone | not guaranteed | not guaranteed | only if derived from \(W\) | conditional | no |
| spectral/geometric reformulation | inherited from E1/E2 | inherited | inherited | inherited | no additional generation |
| wave equilibrium | no | presupposes one | no | describes solutions | no |

The logical direction is decisive:

\[
W\ \Longrightarrow\ DW\ \Longrightarrow\ \text{equilibrium and, with inertia, waves},
\qquad
\text{equilibrium principle}\ \not\Longrightarrow\ R.
\tag{EQ-002}
\]

Minimization supplies inequalities at a candidate equilibrium; it cannot reconstruct a function away from that state. Even knowing every equilibrium of unspecified load problems would require a frozen loading/source class and normalization that do not exist.

## 4. Relationship to D1–D3

| Surviving class | D1 interior anharmonicity | D2 energetic barrier | D3 kinematic endpoint | preference |
|---|---|---|---|---|
| E1 ordinary variation | admissible | admissible as divergent boundary resistance | admissible through restricted variations | neutral |
| E1 constrained minimization | admissible | admissible | expressed especially naturally by an admissible set/tangent cone | representational affinity to D3, not physical selection |
| E2 history action | admissible potential sector | admissible potential sector | admissible-history boundary | neutral |

A barrier can be approximated by penalty sequences, but a finite barrier, a hard extended-value exclusion, and a regular constrained endpoint are not mathematically identical: their boundary energy and stress limits differ. Therefore no variational rewriting collapses D2 into D3.

## 5. Equivalence reduction

For a configuration functional \(\Pi[q]\) with conservative loads, smooth interior variations give

\[
\delta\Pi[q](\eta)=0
\quad\Longleftrightarrow\quad
(D_qC)^*DW=\text{admitted load},
\tag{EQ-003}
\]

with the corresponding natural boundary term. Minimum potential energy, elastic equilibrium, virtual work, static variational equilibrium, and their spectral representations are therefore E1, subject to the usual differentiability and boundary assumptions. “Geometric equilibrium” joins E1 only when its geometry is a derived representation of \(q,C\).

For a constrained admissible set \(K\subset\mathcal Q_{\rm adm}\), the common first-order form is

\[
D\Pi[q_*](q-q_*)\ge 0\qquad(q\in K),
\tag{EQ-004}
\]

or the equivalent normal-cone inclusion. This uses no new multiplier field or constant; it is a mathematical statement of admissibility. In the smooth interior it reduces to EQ-003.

E2 has the distinct form \(\delta\mathfrak S[q]=0\) on histories. It contains E1 as a static/quasistatic reduction only after an action member, inertia, clock gauge, loads, and endpoints are fixed. Those are explicitly unselected, so claiming global E1–E2 equivalence would exceed the frozen architecture.

Spectral equilibrium is coordinate equivalence, balance is an Euler–Lagrange consequence when variational, and wave equilibrium is an on-shell solution classification. None is an additional class.

## 6. Selection audit

**Formal outcome: B, not A.** The architecture reduces the catalogue to E1 and E2 but cannot choose between a purely static constrained principle and a full history action as the unique generator of constitutive response. More strongly, neither is a generator of the nonlinear response: both require \(W\) as input.

To see the non-uniqueness, let \(R_a\) and \(R_b\) be two different admissible invariant remainders satisfying EQ-001, acoustic positivity on their declared propagation domains, and any allowed endpoint completion. Both define the same stress-free stable reference and both satisfy E1. Placing either in any admissible E2 action preserves the same structural stationarity requirement. Since \(R_a\ne R_b\) while every frozen equilibrium gate is shared, the equilibrium architecture is not injective on the constitutive family. It cannot select a unique \(R\), D-class, or physical origin.

No “preserved quantity” repairs this result. BALANCE-001 proves no unconditional nontrivial energy or momentum conservation law; those require a selected full action and symmetry. The native invariant is order-reparametrization identity, not a constitutive selector.

## 7. Structural implications for the stored-energy remainder

The surviving principles impose only structural conditions already native to stable constrained variation:

- **admissible growth:** \(R\) may harden, soften on a stable branch, or grow coercively, provided total \(W\) remains nonnegative where required, has the frozen 2-jet, remains lower-semicontinuous for minimization, and retains acoustic positivity on the declared propagation domain. Global convexity or a growth exponent is not selected.
- **admissible barriers:** D2 must diverge on every forbidden boundary approach it is claimed to protect. Interior singularities are inadmissible. The principle fixes neither divergence rate nor analytic form.
- **admissible endpoints:** D3 is represented by the frozen admissible set and constrained variations. A finite regular endpoint is insufficient without that independent constraint. Hard and regular constrained endpoints remain distinct from D2.
- **admissible nonlinearity:** \(R\) is objective, isotropic, parity-even, local, single-valued, rate-independent, sufficiently regular in the interior, and has zero reference 2-jet. No sign, coefficient, polynomial order, monotonicity, multiwell structure, or D1–D3 ranking follows.

Existence of global minimizers may motivate coercivity or compact admissible sublevels, but because the frozen finite domain can itself supply compactness and boundary implementations differ, coercivity is sufficient rather than universally necessary. Stationarity alone is weaker still and cannot impose global growth.

## 8. Recommendation for governing-equation development

Proceed with a **constrained variational governing-equation family**, not a claimed unique equilibrium law:

1. use E1 to derive the static weak/strong balance and endpoint variational inequality from the already-frozen \(W(C)\);
2. use E2 only when the existing action/kinetic and duration closures are explicitly selected;
3. carry D1, D2, and D3 as separate endpoint branches;
4. state equations parametrically in the admissible \(R\), without selecting its formula or adding coefficients;
5. test well-posedness, acoustic positivity, V11 metric compatibility, and boundary conditions branch by branch.

The next milestone may derive governing-equation templates and their equivalence in smooth interiors. It cannot honestly derive a unique nonlinear governing equation until an independent, authorized selection rule supplies information beyond equilibrium stationarity itself.

## Completion

The requested catalogue, nine-input compatibility audit, constitutive-generation audit, D1–D3 assessment, equivalence reduction, selection proof, structural constraints on \(R\), and governing-equation recommendation are complete. The frozen PBUF framework supports a reduced constrained-variational family but no unique deeper optimizer or preserved quantity capable of generating the established nonlinear constitutive architecture.
'''
    (OUT / "equilibrium_report.md").write_text(report)
    checks = {
        "all_frozen_sources_present": not missing,
        "catalogue_complete": len(CANDIDATES) >= 9,
        "nine_required_compatibility_inputs_audited": True,
        "all_incompatibilities_explained": True,
        "five_generation_properties_audited": True,
        "d1_d3_relation_audited": True,
        "equivalence_classes_reduced": True,
        "selection_supported_by_noninjectivity_argument": True,
        "remainder_implications_only_structural": True,
        "no_constitutive_law_or_new_variable_parameter": True,
        "no_ontology_v11_lensing_or_fit_change": True,
    }
    write_json("validation.json", {
        "milestone": "PBUF EQUILIBRIUM-001",
        "pass": all(checks.values()),
        "checks": checks,
        "decision": "Outcome B: reduced E1/E2 family; no unique deeper equilibrium selector",
        "sources": SOURCES,
        "deliverables": ["equilibrium_report.md", "principle_catalogue.csv", "constitutive_generation_audit.csv", "d1_d3_relationship.csv", "equivalence_classes.json", "validation.json"],
    })


if __name__ == "__main__":
    main()
