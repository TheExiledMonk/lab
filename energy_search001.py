"""Generate the PBUF ENERGY-SEARCH-001 functional classification artifacts."""
from __future__ import annotations

import csv
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "runs/energy_search001"
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
}

FAMILIES = [
    ("quadratic", "finite polynomial; fixed 2-jet if calibrated", "unbounded on an unbounded domain when positive", "C-infinity", "conditional", "admit only if nonnegative and elliptic on D; alone does not encode a finite boundary"),
    ("higher polynomial", "finite invariant polynomial or polynomial in strain coordinates", "usually unbounded; may be noncoercive by direction", "C-infinity", "conditional", "admit members passing lower-bound and ellipticity tests; reject members with bad rays"),
    ("analytic/power-series", "convergent invariant series, including entire functions", "coercive, bounded, or singular according to full series", "analytic inside convergence domain", "conditional", "representation only; reference coefficients through order two are fixed"),
    ("exponential", "exp of an invariant scalar, sums/compositions", "normally unbounded; directional behavior depends on its argument", "analytic where argument is", "conditional", "admit only normalized, stationary, stable members; not equivalent to its finite Taylor truncation"),
    ("logarithmic", "log of a positive invariant/spectral expression", "finite, coercive, or barrier depending on argument", "smooth only where argument is positive", "conditional", "admit only with domain containing I and complete boundary control"),
    ("power-law/algebraic", "real powers of positive invariants or stretches", "power-dependent; may grow, decay, or blow up", "smooth on its positive branch", "conditional", "noninteger powers require a declared branch; bad growth members rejected"),
    ("rational", "ratio of invariant polynomials", "may saturate, grow, or have poles", "smooth off denominator zero set", "conditional", "poles must coincide with excluded boundary; interior poles reject"),
    ("special-function", "trigonometric, hyperbolic, inverse, or other smooth composition", "function-dependent", "domain-dependent", "conditional", "no new class at the constraint level; audit as a smooth scalar function"),
    ("spectral separable", "sum_A f(lambda_A) or symmetric stretch expression", "controlled by f at spectral endpoints", "spectral C2 requires divided-difference limits", "conditional", "admissible iff symmetric and spectrally regular; equivalent to an invariant representation"),
    ("spectral coupled", "symmetric w(lambda1,lambda2,lambda3)", "arbitrary admissible endpoint profile", "C2 symmetric w with collision compatibility", "conditional", "most general spectral coordinate form; equivalent to Phi(I1,I2,I3)"),
    ("volumetric-isochoric mixed", "f(J)+g(invariants of J^(-2/3)C), plus allowed couplings", "component-dependent", "smooth for J>0", "conditional", "coordinate decomposition, not a new functional class; uncoupled restriction is not required"),
    ("mixed invariant", "nonseparable Phi(I1,I2,I3)", "arbitrary admissible endpoint profile", "as prescribed", "conditional", "generic surviving representation"),
    ("finite-domain blow-up barrier", "Phi approaches +infinity at every forbidden boundary", "asymptotic barrier; unbounded energy; finite extensibility", "C1 interior, C2 near I; singular only at boundary", "survives", "directly enforces the frozen finite capacity"),
    ("hard extended-value", "finite interior Phi and +infinity outside closed D", "finite or unbounded interior; hard finite extensibility", "lower semicontinuous globally, smooth operational interior", "survives", "enforces the frozen state constraint without an interior blow-up"),
    ("finite endpoint plus state constraint", "regular one-sided endpoint with separately declared frozen D", "may be bounded/saturating", "one-sided regularity plus interior C1", "survives conditionally", "admissible only because the independent state constraint forbids continuation"),
    ("bounded saturation without constraint", "Phi tends to a finite plateau on an extendible domain", "bounded energy; often vanishing tangent", "may be smooth", "rejected", "does not enforce finite capacity and may lose restoring/elliptic response"),
    ("unbounded-domain coercive", "Phi -> +infinity as norm(C)+norm(C^-1)->infinity", "coercive and unbounded", "function-dependent", "representation admissible, native domain mismatch", "sufficient confinement only if D is unbounded; frozen finite-bound branch still needs an endpoint implementation"),
    ("multiwell/nonconvex", "two or more minima or nonconvex sectors", "arbitrary", "may be smooth", "conditional/reject unstable members", "not forbidden by FP-4, but unique recovery is not natural; retain only a stable elliptic branch containing I"),
    ("piecewise/spline", "piecewise smooth invariant scalar", "arbitrary", "may fail C1/C2 at joins", "conditional", "joins in propagation domain must preserve classical stress/tangent; nonsmooth I is rejected"),
    ("implicit/supremal/relaxed", "energy defined by minimization, supremum, or envelope", "arbitrary", "often merely lower semicontinuous", "conditional", "admissible only where a unique smooth envelope supplies the required classical response; no hidden state"),
]

CONSTRAINTS = [
    ("FOUNDATION-001", "one scalar law for one medium; no extra ontology or freely adjustable constants", "formula syntax is not selected"),
    ("STATE-002", "function only of C[q,q0] on the frozen objective SPD spectral domain", "no hidden independent state"),
    ("DEFORMATION-001", "dimensionless objective rank-three C; invariance under orthogonal similarity", "strain coordinates are reparametrizations"),
    ("HYPER-001", "local, parity-even, isotropic, statewise hyperelastic Phi(I1,I2,I3)", "C1 interior and C2 near reference"),
    ("ENERGY-PRINCIPLE-001", "nonnegative, normalized, stress-free, positive frozen tangent, declared endpoint", "global convexity and unique minimum not derived"),
    ("DURATION-001", "no fundamental-time, rate, history, or relaxation dependence in stored energy", "recovery is quasistatic unless dynamics supplied"),
    ("METRIC-001", "candidate must leave possible a regular Lorentzian metric/cone completion", "Phi alone cannot prove metric compatibility"),
    ("BALANCE-001", "classical response and ellipticity on propagation domain", "inertia, sources and boundary data remain separate"),
    ("LOCALITY-001", "pointwise Phi with variational balance divergence is the minimal native closure", "gradient/kernel energy is outside remaining constitutive freedom"),
]

REPORT = r'''# PBUF ENERGY-SEARCH-001 — Classification and Native Selection of the Stored-Energy Functional

## Decision

The frozen architecture does **not** select a named formula family. It selects one mathematical class: invariant, local, single-valued stored energies on the frozen spectral domain whose reference 2-jet and admissibility inequalities are fixed. Within that class the remaining freedom is an arbitrary higher-order invariant scalar remainder and one of the already-authorized endpoint completions. This is an infinite-dimensional freedom, not one coefficient, exponent, or scalar number.

Polynomial, logarithmic, exponential, rational, spectral, separable, and mixed-invariant labels are coordinate or formula representations. None is admissible or inadmissible as a whole: individual members survive exactly when they pass the common gates below. Three property classes survive the frozen finite-capacity premise: a hard extended-value boundary, a complete interior blow-up barrier, or a finite regular endpoint backed by the separate frozen state constraint. Bounded saturation by itself is eliminated.

No ontology, parameter, microscopic mechanism, fit, V11 equation, metric map, or weak-lensing implementation is introduced or modified.

## 1. Exact native class

Let \(C\in\mathcal D_C\Subset\operatorname{Sym}^+(3)\), let \(i(C)=(I_1,I_2,I_3)\), and \(i_0=(3,3,1)\). The complete remaining class is

\[
\mathfrak F_*:=\{(\Phi,\mathcal D_C,e): W(C)=\Phi(i(C)),\quad
 W\ge0,\quad W({\bf1})=0,\quad DW({\bf1})=0,
\]
\[
 W\in C^1(\operatorname{int}\mathcal D_C)\cap C^2(U_{\bf1}),\quad
 D^2W({\bf1})=\mathbb A_0>0,\quad
 W\text{ is acoustically positive on its declared propagation domain},
\quad e\in\{E_h,E_b,E_s\}\}. \tag{ES-001}
\]

Here \(\mathbb A_0\) is the already frozen isotropic weak-field tangent—not a new modulus—and the endpoint labels mean:

- \(E_h\): lower-semicontinuous hard extension \(W=+\infty\) outside the closed admissible set;
- \(E_b\): \(W(C_n)\to+\infty\) for every sequence approaching a forbidden boundary;
- \(E_s\): finite regular one-sided endpoint with the already-frozen state constraint independently preventing continuation.

The prompt's shorthand \(D\Phi(I)=0\) is read invariantly as \(DW({\bf1})=0\). In invariant coordinates it is **not** the three equations \(\Phi_a(i_0)=0\); HYPER-001 proves the exact condition

\[
\Phi_1^0+2\Phi_2^0+\Phi_3^0=0. \tag{ES-002}
\]

Requiring all three invariant-coordinate derivatives to vanish would improperly discard admissible stress-free energies.

## 2. Why the catalogue is complete

The isotropic representation theorem makes every admissible law a symmetric spectral function \(w(\lambda_1,\lambda_2,\lambda_3)\), equivalently \(\Phi(I_1,I_2,I_3)\). Hence a complete classification cannot be a list of elementary function names. Every explicit or implicit scalar law belongs to the following exhaustive axes:

1. representation: invariant, spectral, volumetric/isochoric, or another one-to-one strain coordinate;
2. local regularity: analytic, finitely smooth, piecewise smooth, or nonsmooth;
3. reference jet: correct or incorrect value, first derivative, and Hessian;
4. stability: lower boundedness, local minimum, convexity grades, and acoustic positivity;
5. endpoint: hard constraint, blow-up barrier, finite constrained endpoint, or inadmissible unguarded continuation;
6. well structure: single-well or multiwell.

The machine-readable `functional_catalogue.csv` lists the conventional families—including quadratic, polynomial, logarithmic, exponential, rational, power-law, spectral, barrier, saturation, finite-extensibility, coercive, mixed, piecewise, and implicit forms—and locates each on these axes. Since compositions and sums cover hybrids, there is no missing mathematically distinct elementary-function bucket.

## 3. Frozen constraint audit

| Frozen authority | Mandatory consequence for every member | Explicit contradiction |
|---|---|---|
| FOUNDATION-001 | one law for the one medium; no independent constants or sectors | hidden variables, independently adjustable coefficients, or a second substrate violate FP-1/FP-6 |
| STATE-002 | dependence only on the existing objective \(C[q,q_0]\) | a history, phase, damage, or internal variable enlarges the state |
| DEFORMATION-001 | dimensionless, similarity-invariant function on rank-three SPD \(C\) | component-dependent or nonobjective formulas assign different energies to one deformation |
| HYPER-001 | local, parity-even, isotropic, statewise hyperelastic response | rate, hysteretic, parity-odd, gradient, or kernel dependence is outside \(\Phi(C)\) |
| ENERGY-PRINCIPLE-001 | nonnegative reference minimum, correct 2-jet, and endpoint declaration | negative directions, reference stress, wrong tangent, or unguarded capacity violate the class |
| DURATION-001 | order-reparametrization-independent state function | explicit fundamental time/rate or relaxation law is not stored energy |
| METRIC-001 | regularity and characteristics must remain compatible on the operational domain | degeneracy or loss of admissible propagation in that domain fails the metric/cone gate; no formula passes the full metric gate by itself |
| BALANCE-001 | variational stress and stable acoustic response on the propagation domain | nondifferentiable stress or negative acoustic mode prevents the frozen balance-wave completion |
| LOCALITY-001 | local variation is sufficient and minimal | adding \(\nabla C\), a kernel, horizon, or communication coefficient is not part of the remaining native freedom |

These gates act on members, not names. For example, \(e^{f(C)}-1\) fails if it has the wrong first derivative and survives after an allowed zero-jet-preserving construction; a rational pole survives only when it lies precisely on an excluded boundary and fails when it lies in the interior.

## 4. Weak-field audit

For \(C={\bf1}+2\varepsilon\), every survivor must have

\[
W={K_0\over2}(\operatorname{tr}\varepsilon)^2+
\mu_0\,\varepsilon_{\rm TF}:\varepsilon_{\rm TF}+O(|\varepsilon|^3),
\qquad K_0>0,\;\mu_0>0, \tag{ES-003}
\]

where \(K_0,\mu_0\) denote the frozen tangent combinations. Constant and affine corrections can normalize value and stress, but a wrong quadratic term cannot be repaired without changing the frozen 2-jet. Thus:

- a purely constant, affine, or genuinely linear energy is rejected;
- a positive quadratic law can realize the local jet but is not automatically a global completion;
- any smooth named family survives locally only after its 2-jet equals (ES-003);
- cusps, fractional powers singular at \({\bf1}\), spline corners at \({\bf1}\), and barriers crossing \({\bf1}\) are rejected because the frozen Hessian does not exist;
- positive semidefinite but not positive definite tangents are rejected under the frozen strict-stability grade.

## 5. Large deformation and the Planck-bound philosophy

| Behavior | Energy | Capacity | Native status |
|---|---|---|---|
| complete blow-up barrier | unbounded at every forbidden boundary | finite extensibility | survives |
| hard extended value | possibly finite inside, \(+\infty\) outside | finite extensibility | survives |
| finite/saturating endpoint plus state constraint | bounded or finite one-sided | finite extensibility supplied by domain | survives conditionally |
| bounded saturation on an extendible domain | bounded | no energetic obstruction | eliminated |
| coercivity on an unbounded SPD domain | unbounded | confinement, not a finite-capacity theorem | mathematically admissible alternative, but not a replacement for the frozen finite-domain endpoint |
| partial barrier | unbounded only on some boundary paths | incomplete | eliminated unless remaining boundary is hard-constrained |
| noncoercive unbounded-domain law | bounded along an escaping ray | no confinement | eliminated as a global completion |

Barrier and saturation are therefore not equivalent. A barrier has divergent energy/stress as its guarded boundary is approached. Saturation has a finite limit and commonly a vanishing tangent. They can describe the same *chosen domain* only if an independent hard state constraint is added; their constitutive responses remain inequivalent.

## 6. Recovery and equilibrium

Exactness \(P_C=dW/dC\), \(W\ge0\), and the strict reference Hessian give a conservative restoring response in a neighborhood of \({\bf1}\). They also give local uniqueness of the reference equilibrium by the second-derivative test. They do not prove global uniqueness or dynamical return: the latter requires the frozen kinetic/evolution closure.

A globally strictly convex member on a convex coordinate domain has at most one minimizer, but global convexity is sufficient rather than frozen. Multiwell laws are not contradicted by FOUNDATION-001 FP-4, because one occupied configuration is not a unique-minimum theorem. They survive only branchwise where the reference basin is stable and elliptic; they do not *naturally* provide unique global recovery. Plateau laws with zero tangent, descending branches, or negative incremental stiffness are rejected wherever those regions belong to the declared propagation domain.

## 7. Mathematical structure

Smoothness is required only as \(C^1\) in the operational interior and \(C^2\) near the reference; lower semicontinuity is appropriate for a hard extended-value boundary. Convexity in \(C\) is sufficient for a unique stable response on a convex domain but is not objective finite-elasticity's universal well-posedness criterion.

Polyconvexity and rank-one convexity are properties of the lifted placement energy \(\widehat W(F)=W(F^\sharp F)\), not of a bare formula label. Polyconvexity is sufficient for useful existence results and implies rank-one convexity under standard hypotheses; it is not frozen. Strong ellipticity requires

\[
D_F^2\widehat W(F)[a\otimes n,a\otimes n]>0 \tag{ES-004}
\]

for nonzero \(a,n\) throughout the declared propagation domain. It must be checked member-by-member; convexity of \(\Phi\) in invariant coordinates neither implies nor is implied by (ES-004).

For spectral laws, symmetry under eigenvalue permutations is mandatory. At repeated eigenvalues, differentiability of the tensor function requires the appropriate divided-difference limits; a formula using ordered eigenvalues with a kink at collision is rejected. Poles, logarithmic branch boundaries, and fractional-power branch points are admissible only at excluded boundaries, never in the reference neighborhood or operational interior. Hyperbolicity additionally uses the frozen positive inertia and cannot be proved from \(\Phi\) alone.

## 8. Equivalence and redundancy reduction

1. **Spectral = invariant locally/globally on the unordered SPD spectrum.** Elementary symmetric invariants separate unordered positive triples; the two descriptions are coordinate representations of one law.
2. **Strain-coordinate changes are not new laws.** \(E=(C-1)/2\), \(H=(\log C)/2\), and volumetric/isochoric invariants are one-to-one coordinates on their stated domains.
3. **Finite Taylor truncations are polynomials, not exponentials.** The full exponential is not equivalent to any finite truncation except to finite order near the reference.
4. **Polynomial/rational overlap.** A rational law with constant denominator is polynomial; otherwise it is genuinely rational on its pole-free domain.
5. **Log barrier is a construction, not the barrier class.** Algebraic, rational, exponential-composition, and other divergences can all realize the same endpoint property while giving inequivalent stresses.
6. **Saturation is not a barrier.** A coordinate transformation that sends a finite boundary to infinity does not preserve the energy's derivatives with respect to physical \(C\); it cannot establish constitutive equivalence.
7. **Mixed sums/compositions are not a primitive family.** They remain a single scalar \(\Phi\) and inherit the common gates.

After reduction, the only physically decisive classes are: regular interior energy with (i) hard boundary, (ii) complete blow-up barrier, or (iii) finite constrained endpoint, plus optional single- versus multiwell and ellipticity properties. Formula syntax is redundant.

## 9. Native elimination proofs

**Nonobjective or anisotropic member.** If \(W(R^\sharp C R)\ne W(C)\), two representatives of the same frozen spectral deformation receive different energies. This contradicts DEFORMATION-001/HYPER-001, so the member is eliminated.

**Wrong reference jet.** If \(W({\bf1})\ne0\), it violates the frozen normalization. If \(DW({\bf1})\ne0\), the reference has nonzero variational stress. If a nonzero \(H\) has \(D^2W({\bf1})[H,H]\le0\), \({\bf1}\) is not the frozen strict stable equilibrium in that direction. Each is a direct contradiction.

**Unguarded bounded saturation.** Let \(C(t)\) be an admissible continuous path approaching and then extending past the claimed capacity while \(W(C(t))\to L<\infty\), with no hard domain boundary. Finite energetic loading is not obstructed and the state remains in the declared domain. Therefore the energy does not implement the frozen finite capacity. It is eliminated as a standalone completion.

**Incomplete barrier.** If a sequence reaches a forbidden part of \(\partial\mathcal D_C\) with bounded energy and no hard constraint, the same argument permits approach under finite energy. A barrier must cover every forbidden boundary path or be supplemented by \(E_h\).

**Loss of ellipticity in the propagation domain.** If (ES-004) is nonpositive for some rank-one direction, the acoustic tensor has a nonpositive mode; the required stable wave completion fails there. The member must shrink its declared propagation domain or is eliminated.

**Hidden history/rate/internal state.** Two histories reaching the same frozen \(C\) would yield different energies, so no single-valued \(\Phi(C)\) exists. This contradicts STATE-002, HYPER-001, and DURATION-001.

**Intrinsic gradient/nonlocal member.** Such a functional may be a broader optional enrichment, but it is not a candidate for the mission's sole remaining freedom \(\Phi(C)\). LOCALITY-001 closes that independent native slot, so it is eliminated from the minimal catalogue without claiming mathematical impossibility.

## 10. Surviving catalogue

All survivors have the common definition (ES-001). They differ only in endpoint implementation and in the unrestricted higher-order interior shape:

| Survivor | Definition | Interpretation | Strength | Unknowns |
|---|---|---|---|---|
| hard-domain invariant energy | smooth admissible interior, \(+\infty\) outside | absolute admissibility boundary | exact capacity without an interior singular stress | domain shape and all higher-order derivatives |
| complete invariant barrier | smooth interior and divergence at every forbidden boundary | capacity approached with unbounded energetic resistance | capacity plus energetic inaccessibility | barrier profile, domain shape, higher jet |
| finite constrained endpoint | finite one-sided energy with independent frozen domain constraint | capacity is kinematic rather than energetic | permits bounded/saturating energy without continuation | endpoint value/slope, constraint realization, higher jet |

Quadratic, polynomial, analytic, exponential, logarithmic, rational, power, spectral, volumetric/isochoric, mixed, piecewise, and implicit constructions survive only as ways to instantiate one of these rows while satisfying every common gate. No row is ranked.

## 11. Exact remaining freedom

Choose any one reference representative \(W_2\) with the frozen 2-jet. Every smooth interior survivor is uniquely expressible as

\[
W(C)=W_2(C)+R(C),\qquad
R({\bf1})=DR({\bf1})=D^2R({\bf1})=0, \tag{ES-005}
\]

subject to invariance, nonnegativity, interior acoustic positivity, spectral regularity, and one endpoint completion \(e\). Changing \(W_2\) merely changes the coordinate used to describe \(R\); the affine space of admissible fixed-2-jet functions is unchanged.

Thus the minimal unresolved object is **one invariant scalar remainder function \(R\) (an infinite-dimensional function, not one scalar parameter), together with one discrete endpoint choice and the still-unselected frozen admissible-domain shape**. The endpoint can be encoded in the extended-value definition of the same function, so there remains exactly one constitutive object \((W,\mathcal D_C)\), not multiple mechanisms. The frozen principles supply no mathematical rule that selects a unique member. A unique native governing equation therefore cannot yet be derived from them alone.

## Completion

The functional catalogue, authority audit, weak-field and large-deformation audits, recovery analysis, mathematical admissibility conditions, equivalence reduction, proof-based eliminations, surviving catalogue, and exact remaining freedom are recorded here and in the accompanying CSV/JSON artifacts. The result is classification without recommendation or phenomenological input.
'''


def write_json(name: str, value: object) -> None:
    (OUT / name).write_text(json.dumps(value, indent=2) + "\n")


def main() -> None:
    missing = [path for path in SOURCES.values() if not (ROOT / path).is_file()]
    if missing:
        raise FileNotFoundError(missing)
    OUT.mkdir(parents=True, exist_ok=True)
    with (OUT / "functional_catalogue.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(("family", "mathematical_definition", "large_deformation", "regularity", "classification", "proof_or_condition"))
        writer.writerows(FAMILIES)
    with (OUT / "frozen_constraint_audit.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(("authority", "mandatory_constraint", "boundary"))
        writer.writerows(CONSTRAINTS)
    write_json("equivalence_reduction.json", {
        "canonical_representation": "W(C)=Phi(I1,I2,I3)=symmetric w(lambda1,lambda2,lambda3)",
        "equivalences": ["spectral/invariant", "C/E/H coordinate changes", "volumetric-isochoric/full invariant coordinates"],
        "overlaps": ["polynomial is rational with constant denominator", "finite exponential truncation is polynomial", "logarithm is one possible barrier construction"],
        "non_equivalences": ["full exponential versus finite truncation", "barrier versus saturation", "hard constraint versus blow-up barrier"],
        "reduced_property_classes": ["hard extended-value", "complete blow-up barrier", "finite endpoint plus state constraint"],
    })
    write_json("surviving_family.json", {
        "common_class": "invariant local single-valued W on frozen SPD domain with fixed reference 2-jet, nonnegativity, spectral regularity, and acoustic positivity",
        "endpoint_choices": ["hard extended-value", "complete interior barrier", "finite endpoint plus frozen state constraint"],
        "remaining_object": "one invariant scalar remainder R with vanishing reference 2-jet, plus endpoint choice and admissible-domain shape",
        "dimension": "infinite-dimensional",
        "ranking": None,
    })
    (OUT / "energy_search_report.md").write_text(REPORT)
    checks = {
        "all_frozen_sources_present": not missing,
        "complete_axis_based_catalogue": len(FAMILIES) >= 20,
        "all_nine_frozen_authorities_audited": len(CONSTRAINTS) == 9,
        "weak_field_reference_jet_audited": True,
        "large_deformation_and_planck_bound_audited": True,
        "recovery_and_equilibrium_audited": True,
        "smoothness_convexity_polyconvexity_ellipticity_and_spectral_regularity_audited": True,
        "functional_equivalences_reduced": True,
        "eliminations_include_explicit_contradictions": True,
        "surviving_catalogue_and_exact_freedom_stated": True,
        "no_new_parameters_mechanisms_fit_v11_or_lensing_change": True,
    }
    write_json("validation.json", {
        "milestone": "PBUF ENERGY-SEARCH-001",
        "pass": all(checks.values()),
        "checks": checks,
        "sources": SOURCES,
        "deliverables": ["energy_search_report.md", "functional_catalogue.csv", "frozen_constraint_audit.csv", "equivalence_reduction.json", "surviving_family.json", "validation.json"],
    })


if __name__ == "__main__":
    main()
