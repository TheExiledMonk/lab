# PBUF CONSTITUTIVE-CONSTRUCTION-001 — Native Construction of the Minimal Stored-Energy Functional

## Decision

The smallest interior stored energy is the frozen quadratic 2-jet itself.  Put

\[
E(C):={1\over2}(C-\mathbf1),\quad t=\operatorname{tr}E,\quad
E_{\rm TF}=E-{t\over3}\mathbf1,
\]

and let (K_0>0,mu _0>0) be the already-frozen weak-field tangent combinations.  Then

\[
\boxed{Q(C)={K_0\over2}t^2+\mu _0E_{\rm TF}:E_{\rm TF}
       ={1\over2}\mathbb A_0[E,E].} \tag{CC-001}
\]

This introduces no fitted coefficient: (K_0,mu _0), equivalently (mathbb A_0), are frozen constitutive data.  There are two tied minimal completions of (Q), one using a hard extended-value boundary and one using the separately frozen state constraint.  The third endpoint class, a complete interior barrier, cannot have a uniquely simplest explicit member because neither the boundary shape nor a barrier profile is frozen.

Thus there is a uniquely minimal **interior quadratic polynomial for a specified (mathbb A_0)**, but there is no uniquely minimal complete stored-energy pair ((W,\mathcal D_C)).  The result is minimality, not uniqueness of the nonlinear constitutive law.

## Frozen requirement checklist

1. One scalar stored-energy law for the one medium; no independent energetic sector or freely adjustable constant.
2. Dependence only on the frozen objective rank-three SPD state (C[q,q_0]).
3. Dimensionless-deformation dependence invariant under (C\mapsto R^\sharp C R), equivalently (W=\Phi(I_1,I_2,I_3)).
4. Local, parity-even, isotropic, single-valued, statewise hyperelastic, and independent of fundamental time, rate, and history.
5. An admissible (\mathcal D_C) that is path-connected, permutation invariant, contains (mathbf1), and has compact closure inside the SPD cone.
6. (C^1) regularity in (\operatorname{int}\mathcal D_C), (C^2) near (mathbf1), and lower semicontinuity for an extended-valued completion.
7. (W\geq0) and (W(\mathbf1)=0).
8. (DW(\mathbf1)=0).
9. (D^2W(\mathbf1)=\mathbb A_0>0), with the frozen isotropic channels (K_0>0,mu _0>0).
10. Stress obtained by variation and acoustic positivity on the declared propagation domain.
11. Exactly one authorized endpoint implementation: hard extended value, complete interior blow-up, or finite regular endpoint plus the separately frozen state constraint.
12. Compatibility of the operational branch with a later regular V11 metric/cone completion; the energy alone cannot establish that completion.
13. No intrinsic gradient, kernel, hidden state, dissipation, fitted coefficient, new ontology, V11 modification, or weak-lensing modification.

The CSV checklist is the definitive machine-readable version.  The statements above add no interpretation to it.

## Candidate Construction A — hard quadratic

On the already-frozen admissible set, define

\[
W_A(C)=\begin{cases}
Q(C),&C\in\overline{\mathcal D_C},\\
+\infty,&C\notin\overline{\mathcal D_C}.
\end{cases} \tag{CC-002}
\]

This is the simplest self-contained capacity completion.  Its interior has two isotropic quadratic channels, no higher-order remainder, and no nonlinear scalar function.

The conjugate response and tangent are

\[
P_C^A=DW_A={K_0\over2}t\mathbf1+\mu _0E_{\rm TF},\qquad
D P_C^A[H]={K_0\over4}\operatorname{tr}H\,\mathbf1
             +{\mu _0\over2}H_{\rm TF}. \tag{CC-003}
\]

Hence the material response is exactly linear in (C).  Under the finite-deformation realization (C=F^\sharp F), the first-Piola response (2FP_C^A) is geometrically nonlinear in (F); no additional constitutive nonlinearity has been assumed.

## Candidate Construction B — constrained finite-endpoint quadratic

Define

\[
W_B(C)=Q(C),\qquad C\in\overline{\mathcal D_C}, \tag{CC-004}
\]

and admit no continuation beyond the separately frozen state domain.  It has the same interior energy, stress, tangent, and propagation set as A.  Its one-sided boundary energy and stress are finite because (\overline{\mathcal D_C}\Subset\operatorname{Sym}^+(3)).

A and B are genuinely different extended constitutive constructions: A encodes inadmissibility in the energy, while B leaves it in the already-authorized state constraint.  They are not different interior response laws.

## Candidate Construction C — complete barrier class

Let (b) be any nonnegative invariant scalar that is smooth in the interior,
finite near the reference, and is a complete boundary barrier:

\[
b(C_n)\to+\infty
\]

for every sequence (C_n\to\partial\mathcal D_C) from the interior.  Then

\[
\boxed{W_C(C)=Q(C)+{Q(C)^2\over K_0}b(C).} \tag{CC-005}
\]

The added term is nonnegative and has zero reference 2-jet because (Q=O(|E|^2)).  Since the compact forbidden boundary excludes the interior reference and (Q>0) away from it, the term diverges on every boundary approach.  Division by the frozen (K_0) only preserves energy dimension and introduces no coefficient.  The response is

\[
P_C^C=P_C^A+{2Qb\over K_0}P_C^A+{Q^2\over K_0}Db.
\]

This constructs the complete barrier family without selecting a barrier profile.  Smooth proper exhaustion functions exist on open finite-dimensional domains; an invariant member is part of the authorized barrier class.  An explicit logarithm, reciprocal distance, rational pole, or power would additionally choose a boundary gauge or divergence profile that the frozen framework never supplied.

Candidate C is admissible on the connected component containing (mathbf1) of the set where its lifted acoustic form is positive.  A proposed (b) is a full-domain candidate only if that component equals the declared propagation domain.  Nonnegativity and the zero 2-jet alone do not prove global ellipticity.

No fourth primitive construction exists at the property level: ENERGY-SEARCH-001 proved that every admissible endpoint is A-type, B-type, or C-type.  Infinitely many formula-level constructions nevertheless remain because (b), and more generally any admissible invariant remainder with zero reference 2-jet, is infinite-dimensional.  Therefore enumeration of formulas does not literally terminate; equivalence reduction terminates at these three endpoint classes.

## Six-property audit

| Property | A: hard quadratic | B: constrained quadratic | C: complete barrier |
|---|---|---|---|
| Weak-field limit | exactly (K_0t^2/2+\mu_0E_{\rm TF}:E_{\rm TF}) | same | same (+O(|E|^4)) when (b) is finite near the reference |
| Nonlinear response | affine in (C); geometric nonlinearity only after (C=F^\sharp F) | same | (P_C^A+2QbP_C^A/K_0+Q^2Db/K_0), profile unselected |
| Recovery | conservative; (mathbf1) is a strict local minimizer | same | same locally; global recovery depends on the selected remainder and evolution |
| Acoustic positivity | guaranteed at (mathbf1); positive on its declared elliptic component | same | guaranteed at (mathbf1); must be audited across the declared domain |
| Finite deformation | finite quadratic response throughout the compact interior; hard exclusion at boundary | finite quadratic response through the one-sided endpoint | barrier hardening and divergent energy; stress divergence follows when the chosen barrier derivative diverges |
| Endpoint | finite interior trace, (+\infty) outside | finite regular one-sided value, separate constraint | (W\to+\infty) on every forbidden boundary approach |

For precision, if (\widehat W(F)=W(F^\sharp F)), a rank-one increment (a\otimes n) has acoustic form

\[
D_F^2\widehat W(F)[a\otimes n]^2
=D_C^2W[H,H]+2P_C:(|a|^2n\otimes n),
\quad H=F^\sharp a\otimes n+n\otimes F^\sharp a. \tag{CC-006}
\]

At (F=\mathbf1), (P_C=0) and the frozen positive tangent makes (CC-006) positive.  Away from the reference, neither a quadratic (W(C)) nor an unspecified barrier is automatically strongly elliptic.  Defining the propagation domain as the connected positive component is exactly the frozen declared-domain gate, not a claim of global positivity.

Recovery here means conservative restoring response and local equilibrium recovery.  Dynamical return cannot be computed from stored energy without the still-separate kinetic/evolution closure.

## Comparative complexity

| Candidate | Scalar terms | Added nonlinear operations | Added constitutive assumptions | Interior regularity | Compatibility status |
|---|---:|---:|---:|---|---|
| A | 2 (or one bilinear form (\mathbb A_0[E,E]/2)) | 0 | 0 | analytic | passes; metric compatibility remains a downstream gate |
| B | 2 (same (Q)) | 0 | 0 | analytic with finite one-sided endpoint | passes conditionally on the already-frozen state constraint |
| C | 3 including (Q^2b/K_0) | at least 1 | at least 1 function/profile | frozen minimum only | conditional member-by-member |

Counting syntax depends on notation: the two isotropic channels can be displayed as two terms or contracted into one occurrence of (\mathbb A_0).  That ambiguity does not change the ordering because A and B use no operation beyond the mandatory quadratic form, while C adds a nonzero remainder.

## Frozen compatibility audit

A and B satisfy requirements 1–11 and 13 by construction.  Requirement 12 remains an open downstream compatibility gate for every stored energy, as frozen by METRIC-001.  C satisfies the reference, invariance, locality, endpoint, and nonnegativity gates by its definition; interior regularity and full propagation-domain acoustic positivity must be checked for the selected (b).  The accompanying CSV records every candidate–requirement result explicitly.

No candidate changes V11 or weak lensing.  No new variable, modulus, endpoint value, exponent, length, barrier coefficient, or observational input has been inserted.

## Minimality result

**Proposition 1 (polynomial interior minimality).**  Among polynomial functions of (E) with prescribed value, gradient, and nonzero Hessian (\mathbb A_0) at (E=0), the minimum possible degree is two, and the degree-two member is exactly (Q=\mathbb A_0[E,E]/2).

**Proof.** A polynomial of degree zero or one has zero Hessian and violates requirement 9.  Taylor's theorem makes every degree-two polynomial

\[
c+L[E]+{1\over2}B[E,E].
\]

Requirements 7–9 force (c=0), (L=0), and (B=\mathbb A_0).  Isotropy decomposes (B) into the already-frozen volumetric and trace-free channels, giving (CC-001).  Therefore no lower-degree polynomial is admissible and no different quadratic polynomial has the frozen 2-jet.  ∎

**Proposition 2 (remainder minimality).**  Any (C^2) admissible interior energy has (W=Q+R), where (R(\mathbf1)=DR(\mathbf1)=D^2R(\mathbf1)=0).  Setting (R=0) minimizes the number of nonzero higher-order constitutive terms.

This is a rigorous minimality statement for the natural partial order “delete a nonzero zero-2-jet remainder.”  It does not claim that all admissible functions possess a language-independent operation count.

**Impossibility of absolute unique minimality.**  No representation-independent total simplicity order is frozen.  Formula length changes under invariant, spectral, (E), logarithmic-strain, or named-function notation.  Moreover A and B share exactly the same (Q), regularity, response, and frozen compatibility while using two endpoint implementations already declared admissible.  The frozen axioms contain no rule preferring energetic exclusion to a state constraint.  Consequently a uniquely minimal complete construction cannot be proved.  This is an impossibility proof from a tied pair, not a uniqueness claim.

## Dominance and recommendation

A weakly dominates C by term count, nonlinear operations, assumptions, explicit regularity, and immediately auditable compatibility.  B has identical interior complexity to A.  A is the practical recommendation because its extended-value definition packages capacity and energy in one variational object; this is a bookkeeping advantage, not new physics and not a uniqueness theorem.

For governing-equation development, use Candidate A as the baseline constitutive closure:

\[
P_C={K_0\over2}\operatorname{tr}E\,\mathbf1+\mu_0E_{\rm TF}
\]

in the interior, coupled only to the already-frozen balance and future authorized kinetic/source/metric maps.  Carry B as an exactly tied endpoint alternative.  Carry C only if a later authorized principle supplies a domain-defining barrier profile; then re-audit (CC-006) over its declared propagation domain.

The recommendation selects the smallest useful baseline.  It does not claim uniqueness, alter V11, or infer a physical mechanism.
