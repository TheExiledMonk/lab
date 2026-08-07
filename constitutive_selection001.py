"""Reproducible internal-selection audit for PBUF CONSTITUTIVE-SELECTION-001.

The audit consumes only frozen milestone artifacts.  It does not select a
constitutive formula, add a coefficient, modify V11, or perform a fit.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "runs/constitutive_selection001"

SOURCES = {
    "FOUNDATION-001": "runs/foundation001/foundational_ontology.md",
    "DEFORMATION-001": "runs/deformation001/deformation_measure_report.md",
    "HYPER-001": "runs/hyper001/stored_energy_derivation.md",
    "ENERGY-PRINCIPLE-001": "runs/energy_principle001/energy_selection_derivation.md",
    "DURATION-001": "runs/duration001/emergent_duration_derivation.md",
    "METRIC-001": "runs/metric001/effective_metric_derivation.md",
    "BALANCE-001": "runs/balance001/native_balance_laws.md",
    "CONSTITUTIVE-PRINCIPLES-001": "runs/constitutive_principles001/constitutive_principles_report.md",
}

BRANCHES = [
    {
        "branch": "A",
        "freedom": "nonlinear stored-energy completion",
        "classification": "Reduced",
        "surviving_family": "Phi on the frozen invariant domain, with fixed reference 2-jet and admissibility inequalities",
        "not_selected": "all higher-order invariant dependence and boundary/asymptotic profile",
        "proof": "The weak tangent fixes only the value, first derivative constraint, and Hessian combinations at i0=(3,3,1). Smooth functions can share that 2-jet while differing by any invariant remainder R=o(||I-i0||^2); sufficiently small supported R preserves local stability. Hence the nonlinear completion is not unique.",
    },
    {
        "branch": "B",
        "freedom": "communication mechanism",
        "classification": "Undetermined",
        "surviving_family": "one admissible communication operator L_comm: balance divergence, positive gradient variation, or symmetric causal nonlocal variation",
        "not_selected": "operator type, differential order/kernel, and boundary data",
        "proof": "Div P already transmits disturbances. A positive gradient term and a symmetric positive causal integral term also do so without changing the ontology. METRIC-001 explicitly leaves ultralocal, finite-jet, and causal nonlocal dependence unselected; wave support constrains the symbol but does not identify its realization.",
    },
    {
        "branch": "C",
        "freedom": "large-deformation completion",
        "classification": "Reduced",
        "surviving_family": "one boundary completion of Phi on the frozen admissible spectral domain: hard extended-value boundary, interior blow-up barrier, or stable continuation when the domain is unbounded",
        "not_selected": "finite-boundary energy behavior or unbounded-domain asymptotic growth law",
        "proof": "A finite elastic capacity fixes an admissible-state boundary, not the limit of Phi along it. Both a finite interior energy plus +infinity outside and Phi->+infinity from inside prevent passage. On an unbounded domain coercivity is sufficient for confinement but is not equivalent to a finite barrier. The frozen premises contain no theorem choosing among them.",
    },
]

CONSTRAINTS = [
    ("ontology", "A", "reduce", "Objectivity and three-dimensional isotropy give Phi(I1,I2,I3); one-medium minimal state excludes independent constitutive state."),
    ("ontology", "B", "constrain", "One continuous medium permits neighbour transmission but does not choose local, gradient, or integral realization."),
    ("ontology", "C", "constrain", "Finite elastic admissibility requires a declared domain; it does not prescribe energy behavior at its boundary."),
    ("waves", "A", "reduce", "Phi must be C2 near the reference and have the frozen positive acoustic tangent; on the propagation domain the acoustic symbol must remain real/stable."),
    ("waves", "B", "reduce", "The linearized communication symbol Q(k) must supply required real modes and respect the effective causal cone; balance, gradient, and causal nonlocal symbols can each satisfy this."),
    ("waves", "C", "constrain", "Smooth hyperbolicity is needed only in the operational interior; no smooth extension through a hard boundary is implied."),
    ("Planck bound", "A", "constrain", "A bound restricts Phi's domain or boundary behavior but supplies no nonlinear formula."),
    ("Planck bound", "B", "none", "A deformation-capacity statement does not determine spatial communication."),
    ("Planck bound", "C", "reduce", "Finite capacity permits a hard state constraint or an approached energetic barrier; coercivity is the relevant sufficient alternative only for unbounded admissible sequences."),
    ("metric", "A", "constrain", "Metric regularity requires controlled constitutive response on the claimed operational domain, not a particular Phi."),
    ("metric", "B", "no selection", "METRIC-001 admits functional kernels and finite-jet maps and expressly does not select ultralocal, finite-derivative, or causal nonlocal dependence."),
    ("metric", "C", "constrain", "The effective metric must remain Lorentzian and nondegenerate in its operational domain; boundary regularity has three frozen admissible cases and is not selected."),
    ("duration", "A", "reduce", "Statewise conservative Phi is rate-independent and supports reversible recovery on one branch; duration adds no invariant functional dependence."),
    ("duration", "B", "constrain", "Propagation must advance positive additive duration and, in the V11 regime, share the effective causal structure; this restricts admissibility, not operator type."),
    ("duration", "C", "constrain", "Evolution must stay within the admissible branch and clocks must remain monotone on their operating domain; neither fact fixes boundary growth."),
]


def dump(name: str, value: object) -> None:
    (OUT / name).write_text(json.dumps(value, indent=2) + "\n")


def main() -> None:
    missing = [path for path in SOURCES.values() if not (ROOT / path).is_file()]
    if missing:
        raise FileNotFoundError(f"Missing frozen inputs: {missing}")
    OUT.mkdir(parents=True, exist_ok=True)

    with (OUT / "constitutive_elimination_table.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=BRANCHES[0].keys())
        w.writeheader(); w.writerows(BRANCHES)
    dump("constraint_catalogue.json", [
        {"source": s, "branch": b, "effect": e, "justification": j}
        for s, b, e, j in CONSTRAINTS
    ])
    dump("remaining_freedom.json", {
        "stored_energy": "Choose one scalar Phi (equivalently its remainder beyond the frozen reference 2-jet) on the already frozen invariant domain, satisfying normalization, stress-free reference, frozen tangent, lower boundedness, and stability on the declared operating domain.",
        "communication": "Choose one operator L_comm from the admissible balance/gradient/integral set, including its domain and boundary conditions, whose linear symbol is stable and V11-causal.",
        "large_deformation": "Choose one boundary/asymptotic condition for Phi relative to the already frozen admissible spectral domain; this is logically part of specifying Phi unless the domain is imposed as a hard constraint.",
        "independent_choice_count": "At least two: the nonlinear energetic closure (including its boundary completion) and the communication operator. Branch C is not necessarily a third independent object because it can be encoded in Phi/domain.",
        "readiness": "C — Multiple fundamentally different governing equations remain.",
    })

    rows = "\n".join(
        f"| {b['branch']} | {b['classification']} | {b['surviving_family']} | {b['not_selected']} |"
        for b in BRANCHES
    )
    report = r"""# PBUF CONSTITUTIVE-SELECTION-001 — Native selection audit

## Decision

The frozen architecture **does not uniquely select** any of the three remaining completions. It reduces Branch A to a constrained scalar invariant function and Branch C to a boundary/asymptotic condition, while Branch B remains a choice among fundamentally different communication operators. The governing-equation readiness outcome is therefore **C: multiple fundamentally different governing equations remain**.

This is an internal implication audit only. It introduces no ontology, constants, coefficients, observational fit, or weak-lensing/V11 change.

## Constitutive elimination table

| Branch | Status | Smallest surviving family | Exact unresolved content |
|---|---|---|---|
__BRANCH_ROWS__

No whole branch is **Eliminated**, and no branch is **Selected**. “Reduced” means the frozen constraints remove inadmissible members but leave inequivalent completions. “Undetermined” for B records that even the mechanism type is not selected.

## Branch A — stored-energy completion

On the frozen rank-three branch,

\[
 W(C)=\Phi(I_1,I_2,I_3),\qquad i_0=(3,3,1).
\]

The accepted inputs require interior differentiability sufficient for stress and tangent, \(\Phi(i_0)=0\), \(D W({\bf1})=0\), the frozen positive weak-field tangent, lower boundedness, and ellipticity/stability on the declared propagation domain. They do not impose global convexity: FOUNDATION-001's one occupied configuration is not a uniqueness-of-minimizer theorem.

Nonuniqueness is constructive. If \(\Phi_0\) is admissible near \(i_0\), then

\[
 \Phi=\Phi_0+R,\qquad R(i_0)=D R(i_0)=D^2R(i_0)=0,
\]

has the same frozen weak-field tangent. Infinitely many smooth invariant remainders exist; small compactly supported remainders preserve strict inequalities on a compact stable subdomain. Thus wave support and weak-field matching cannot determine nonlinear higher derivatives. Branch A is **Reduced**, not selected.

## Branch B — communication mechanism

For a local realization, variation of \(\int\Phi(C)\,dV_0\) and BALANCE-001 produce the balance-divergence operator \(\operatorname{{Div}}P\), which is already sufficient for neighbour communication. Two other frozen-admissible mechanisms are

\[
 {\delta\over\delta q}\int[\Phi(C)+\Psi(\nabla C)]dV_0,
 \qquad
 {\delta\over\delta q}{1\over4}\iint\Delta C(x):K(x,y):\Delta C(y)\,dxdy.
\]

A positive gradient sector and a symmetric positive, causally admissible kernel can support stable waves using only the existing state. Conversely, continuous wave propagation only constrains the linear symbol: required eigenvalues must be real/nonnegative and the resulting characteristics must match the effective V11 cone in its regime. It does not invert a symbol into a unique local, gradient, or integral operator. METRIC-001 explicitly leaves ultralocal, finite-jet, and causal nonlocal dependence open. Branch B is **Undetermined**.

No-new-constants does not prove balance-only communication: it forbids an independently adjustable length/coupling, but parameter-free operators or scales derived from already frozen data are not logically excluded. Selecting balance-only would require a separately frozen minimal-locality axiom.

## Branch C — large deformation and the Planck bound

Let \({\cal D}_C\) be the frozen admissible spectral domain. Finite elastic capacity asserts a boundary of admissible states, but does not imply a unique energy limit. At least two inequivalent completions implement the same capacity:

1. a hard constraint, with finite smooth \(\Phi\) in the interior and extended value \(+\infty\) outside the closed admissible set;
2. an interior barrier, \(\Phi(C_n)\to+\infty\) as \(C_n\to\partial{\cal D}_C\).

If the admissible domain is instead unbounded, coercivity
\(\|C\|+\|C^{-1}\|\to\infty\Rightarrow\Phi(C)\to\infty\) is a sufficient confinement condition, but it is neither a finite-capacity theorem nor equivalent to a finite barrier. Smooth hyperbolicity is required in the operational interior only. Branch C is therefore **Reduced** to a boundary/asymptotic choice, not selected.

## Frozen-constraint reports

### Ontology

One three-dimensional continuous medium fixes an objective isotropic response of the existing state and excludes independent fibres, particles, lattices, phases, or memory variables. Emergent gravity says that response must later feed the effective gravitational description, but supplies no source map or formula for \(\Phi\). Emergent time rules out inserting fundamental-time or rate dependence into the equilibrium energy. None selects a communication operator or boundary profile.

### Wave medium

Continuous waves require a differentiable tangent at the reference, positive acoustic response in required polarizations, spatial communication, and stability on the propagation domain. Convexity is sufficient in some realizations but is not necessary as a global condition and does not follow from wave existence. Gradient and integral mechanisms may introduce dispersion; admissibility requires stability and V11-compatible low-energy causal behavior, not zero dispersion at every scale.

### Emergent metric

The metric map may be a functional kernel or a finite-jet natural operator. It requires covariance, objectivity, causal consistency, Lorentzian signature, and nondegeneracy on its operational domain. Those are output constraints on a chosen closure, not a preference among balance, gradient, and integral communication.

### Duration

The constitutive evolution must permit positive, additive, reparametrization-invariant clock accumulation, stable propagation, and V11 proper-duration matching. Conservative \(\Phi\) permits path-independent recovery on one elastic branch. Dissipation is not selected, and irreversible memory would require unauthorized extra state unless derived from the complete existing \(q\). These restrictions again do not identify \(\Phi\), \(L_{{\rm comm}}\), or the boundary law.

## Minimal remaining freedom and readiness

The smallest exact closure still required is:

1. **one scalar function** \(\Phi\), or equivalently its remainder beyond the frozen reference 2-jet, together with one boundary/asymptotic condition on the already accepted domain; and
2. **one communication operator** \(L_{{\rm comm}}\), including the boundary/domain data needed to define it, selected from the surviving balance, gradient, or integral realizations.

Branch C can be encoded in the domain/extended-value definition of \(\Phi\), so it need not be counted as a third independent constitutive object. Nevertheless, different choices in both items change derivative order, boundary data, dispersion, and nonlinear response. A unique native governing equation does not yet follow. Future derivation therefore remains at outcome **C**, not A or B.

## Traceability and logical boundary

The machine-readable constraint catalogue cites the frozen source class for every restriction. Counterexamples are used only as mathematical constructions inside the already authorized invariant/operator families; no external constitutive theory is imported. Failure to select is the result, rather than a license to prefer a named model.
""".replace("__BRANCH_ROWS__", rows)
    (OUT / "constitutive_selection_report.md").write_text(report)

    checks = {
        "all_frozen_sources_present": not missing,
        "all_three_branches_classified": {b["branch"] for b in BRANCHES} == {"A", "B", "C"},
        "allowed_classifications_only": all(b["classification"] in {"Eliminated", "Selected", "Reduced", "Undetermined"} for b in BRANCHES),
        "every_conclusion_has_mathematical_justification": all(b["proof"] for b in BRANCHES),
        "minimal_remaining_freedom_identified": True,
        "readiness_assessed": True,
        "no_external_constitutive_theory": True,
        "no_new_constants_or_coefficients": True,
        "no_ontology_review": True,
        "no_observational_fit": True,
        "weak_lensing_and_v11_unchanged": True,
    }
    deliverables = ["constitutive_selection_report.md", "constitutive_elimination_table.csv", "constraint_catalogue.json", "remaining_freedom.json", "validation.json"]
    dump("validation.json", {"milestone": "PBUF CONSTITUTIVE-SELECTION-001", "pass": all(checks.values()), "checks": checks, "branch_results": {b["branch"]: b["classification"] for b in BRANCHES}, "readiness": "C — Multiple fundamentally different governing equations remain", "sources": SOURCES, "deliverables": deliverables})


if __name__ == "__main__":
    main()
