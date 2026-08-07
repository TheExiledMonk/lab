# PBUF NATIVE-SOURCE-001 — Physical Meaning of the Native Loading Term

## 0. Decision and scope

This milestone uses only FOUNDATION-001, DEFORMATION-001, BALANCE-001,
LOCALITY-001, CONSTITUTIVE-CONSTRUCTION-001, and
WEAK-LENSING-LOCALITY-001. It does not select or derive a map from matter,
modify V11, fit data, or add ontology.

The frozen placement problem is

\[
 -\operatorname{Div}_0 P_F=b,
 \qquad P_F=2F P_C,\qquad P_C=D_CW,
 \tag{NS-001}
\]

with weak form

\[
 \int_\Omega P_F:\operatorname{Grad}_0\eta\,dV_0
 =\langle b,\eta\rangle
 +\int_{\Gamma_N}\bar t\cdot\eta\,dA_0 .
 \tag{NS-002}
\]

Accordingly,

\[
\boxed{b\text{ is the distributed external load covector conjugate to an
admissible placement variation.}}
\tag{NS-003}
\]

If this covector has a regular pointwise representative, that representative
is a **referential body-force density**. The invariant classification is the
load covector; “body force” is its regular placement representation.

The symbol is not the generic production term \(\sigma_a\) of BALANCE-001 and
is not the unrelated barrier profile denoted by \(b(C)\) in Candidate C of
CONSTITUTIVE-CONSTRUCTION-001.

## 1. Formal definition and tensor classification

Let \(\mathcal Y_{\bar y}\) be the admissible placement space and let
\(V_0\subset T_y\mathcal Y_{\bar y}\) be its zero-Dirichlet test space. The
native source is

\[
 b\in V_0^*,\qquad \eta\longmapsto\langle b,\eta\rangle ,
 \tag{NS-004}
\]

so its fundamental rank is **one in configuration space**: it is a one-form
(generalized covector), not a scalar energy density or a rank-two stress.

When \(b\) is regular enough to be represented pointwise,

\[
 \langle b,\eta\rangle=\int_\Omega b_i(X)\eta^i(X)\,dV_0,
 \tag{NS-005}
\]

it is a spatial covector-valued scalar density on the reference carrier,
equivalently a section of \(y^*T^*\mathcal S\otimes\mathrm{Vol}(\mathcal
B_0)\). With the frozen unloaded metric used to raise an index it may be
displayed as a spatial vector per reference volume. Thus it has one spatial
index; it has no material flux index. By comparison, \(P_F\) has one spatial
and one referential index, and \(P_FN\) is a surface-load covector.

The distributional definition NS-004 is more general than NS-005. It can
include admissible singular loads without pretending that every source is a
smooth field. Boundary traction remains a separate term and is not silently
included in the body density.

### Classification exclusions

| Candidate meaning | Decision | Reason |
|---|---|---|
| body force | yes, for a regular placement-density representative | it pairs with \(\eta\) over reference volume |
| body couple | no in the frozen baseline | no independent rotation variation, couple stress, or pairing with \(\nabla\eta\) is present |
| energy density | no | \(W\) is the stored-energy density; \(b\) is a first-variation load |
| stress source / stress | no | stress is \(P_C\) or \(P_F\); \(b\) balances its divergence |
| deformation source | not directly | it does not prescribe \(C\); deformation follows after solving balance for \(y\) |
| generalized covector | yes, fundamentally | NS-004 is exactly the weak-form meaning frozen by WL-006 |

## 2. Dimensions and transformations

Because \(C\) is dimensionless, \(W\), \(K_0\), and \(\mu_0\) have the same
dimension: energy per reference volume. For a length-valued placement \(y\),

\[
 [P_F]=\frac{\text{force}}{\text{reference area}}
       =\frac{\text{energy}}{\text{reference volume}},
 \qquad
 [b]=\frac{\text{force}}{\text{reference volume}}
     =\frac{\text{energy}}{\text{reference volume}\,\text{length}}.
 \tag{NS-006}
\]

In SI base units, \([b]=\mathrm{N\,m^{-3}}=\mathrm{kg\,m^{-2}s^{-2}}\).
The coordinate-free statement is
\([\langle b,\eta\rangle]=\) energy and
\([b]=\) energy divided by reference volume and by the unit of \(y\).
This last form remains valid if a later realization uses nondimensional
placement coordinates. No numerical coupling or normalization is fixed here.

Under a change of spatial frame, \(b_i\) transforms as a covector (or its
metric-raised representative as a vector). Under a relabeling of reference
coordinates, the intrinsic load form is unchanged; if the reference volume
form is kept separate, the components transform as a material scalar, while
the combined object \(b_i dV_0\) transforms as a weight-one referential volume
density. Under the internal gauge/diffeomorphism quotient, the pairing in
NS-004 must be invariant and gauge-null variations must receive zero work.
These rules make NS-002 objective and covariant. A bare invariant scalar cannot
replace \(b\), because it cannot pair with a generic vector placement
variation.

## 3. Where the source acts

The source acts **directly on placement/displacement variations** through
\(\langle b,\eta\rangle\). It acts neither on \(C\), \(P_C\), \(P_F\), nor
\(W\) as an independent constitutive argument. Its effect is mediated as

\[
b\;\xrightarrow{\text{balance BVP}}\;y
\;\xrightarrow{F=\operatorname{Grad}_0y}\;C=F^\sharp F
\;\xrightarrow{W}\;(P_C,P_F).
\tag{NS-007}
\]

Stress appears on the opposite, internal side of virtual work. In a dynamic
extension, \(b\) would share the placement equation with the still-unselected
kinetic/inertial operator; that does not change its dual-space type. A load
potential \(-\langle b,y\rangle\) may represent a prescribed dead load, but
such a potential is not required by the frozen milestones and must not be
mistaken for the definition of \(b\).

## 4. Matter and electromagnetic loading

The frozen theory supplies one placement balance and therefore one resultant
native bulk-load slot \(b\), plus boundary loading. It does **not** prove that
ordinary matter and electromagnetic energy use an identical projection rule,
nor does it authorize a new independent medium sector.

Several input channels may therefore remain distinguishable at the missing
projection stage—for example \(\Pi_{\rm m}[T_{\rm m}]\) and
\(\Pi_{\rm EM}[T_{\rm EM}]\)—provided that any bulk mechanical effects enter
the same admissible covector space and combine into the single resultant load

\[
 b=b_{\rm m}+b_{\rm EM}+\cdots .
 \tag{NS-008}
\]

Additivity itself is mandatory only when the future projection asserts
independent superposition or a partitioned accounting; nonlinear dependence on
the total stress-energy is not excluded by the frozen corpus. Separate boundary
traction or background-metric data are also admissible representations of an
excluded exterior. An independent body-couple, microforce, heat, charge, or
second displacement channel would add structure absent from the baseline and
cannot be introduced by SOURCE-PROJECTION-001 without separate authority.

FP-5 requires the completed one-metric theory to retain V11's operational
relativistic behavior, including electromagnetic propagation, but that is a
compatibility gate—not a derivation of universal equality between the matter
and electromagnetic loading maps.

## 5. Locality and mediation

For prescribed \(b|_\Omega\), the native elastic operator is local. The
constitutive response at \(X\) uses the first jet of \(y\), and balance uses
\(-\operatorname{Div}_0P_F\). No averaging, projection, kernel, or constitutive
mediation is required **between the native load and balance**.

This does not prove that the missing map \(T^{\rm matter}\mapsto b\) is local.
WL-LOCALITY-001 explicitly leaves that map open. Local, finite-jet,
distributional, projected, averaged, and nonlocal candidates are not all
equally suitable for a future *local* lensing closure, but the frozen inputs do
not select among them. Likewise the later \(G\) map may still be local or
nonlocal. Elliptic solution dependence on all regional loads and boundary data
is propagation through the local PDE, not constitutive averaging of \(b\).

## 6. Complete dependency graph

The question mark can be replaced only by a source projection whose codomain is
the placement dual. No intermediate energy, deformation, or stress object is
frozen:

```text
physical matter / electromagnetic state
        |
        v
effective stress-energy description T^matter
        |
        v  [MISSING: Pi_source, no formula or locality class frozen]
placement-load covector S_ext in (T_y Y_0)*
        |
        v  [regular representative, when it exists]
referential body-force density b_i dV_0
        |                    + boundary traction t_bar dA_0
        v
static placement balance: -Div_0 P_F = b
        |
        v
placement y -> F=Grad_0 y -> deformation C=F^sharp F
        |
        v  [MISSING: selected metric map G]
effective metric g_eff
        |
        v  [retained V11 null propagation]
weak-lensing observables
```

\(\mathcal S_{\rm ext}\) and regular \(b_i dV_0\) are not two physical
channels: the latter is a representation of the former. The graph also makes
two corrections to the proposed shorthand. Balance solves for placement before
deformation is formed, and deformation reaches weak lensing only through the
still-missing metric map \(G\).

## 7. Constraint catalogue for a future source projection

Let

\[
 \Pi_{\rm source}:T^{\rm matter}\longmapsto
 \mathcal S_{\rm ext}\in (T_y\mathcal Y_0)^*,
 \tag{NS-009}
\]

with \(b\) its regular bulk representative when one exists.

### Mandatory from the frozen theory

1. **Correct codomain and rank.** The output must be a continuous admissible
   placement-load covector; a regular bulk output has one spatial covector
   index and referential volume density weight.
2. **Virtual-work compatibility.** \(\langle\Pi[T],\eta\rangle\) must be a
   scalar with energy dimension and must enter NS-002 with the frozen sign
   convention. It may not be inserted as \(W\), \(C\), or \(P_F\).
3. **Dimensional consistency and normalization disclosure.** Every contraction,
   derivative, scale, and coefficient needed to turn stress-energy into force
   per reference volume must be identified. FP-6 forbids an unexplained new
   free coupling; no normalization is presently frozen.
4. **Covariance, objectivity, and gauge basicness.** The map must be independent
   of coordinate/frame representatives, transform as described in section 2,
   and annihilate pure gauge/null variations.
5. **Reference-carrier compatibility.** A spacetime/effective tensor input must
   be pulled back or otherwise related to \(\mathcal B_0\) by an already
   authorized realization. The projection may not silently identify spatial,
   material, and effective-spacetime indices.
6. **Balance compatibility.** The output must belong to the dual of the chosen
   test space, have enough regularity for the weak problem, and satisfy the
   solvability conditions of the boundary choice. In particular, a pure
   Neumann problem requires resultant force and moment compatibility modulo
   rigid modes.
7. **Constitutive/domain compatibility.** Loads and boundary data must admit a
   solution in the allowed placement class, with \(C\in\overline{\mathcal D_C}\);
   claims using the elliptic theory must remain on a uniformly strongly
   elliptic branch. The map must not alter \(W\), its moduli, or its domain by
   disguise.
8. **One-medium accounting.** Contributions assigned to internal partitions
   must cancel in complete accounting; excluded regional influence may appear
   as bulk load or boundary flux, but no external physical substrate or new
   ontological sector may be added.
9. **V11 compatibility.** Together with balance and the eventual \(G\), the
   weak-field output must permit the retained local Lorentz, one-metric, and
   relativistic limits. It may not modify V11.
10. **No double counting.** The same external influence may not be counted
    simultaneously as body load, boundary traction, and background metric
    unless the decomposition and subtraction are explicitly defined.
11. **Scope separation.** Static \(b\) must not be used to smuggle in an
    unfrozen kinetic law, duration calibration, dissipation, metric map, or
    observational fit.

### Required only for a claimed local weak-lensing closure

12. **Support-controlled locality.** \(b(X)\) must depend on local or explicitly
    bounded/causal source data, with derivative order or kernel support stated.
    Otherwise WL-LOCALITY-001's regional-sufficiency gap remains open.
13. **Regional restriction consistency.** Restricting global data and then
    projecting must agree with the declared regional projection up to the
    separately supplied boundary/background representation of excluded matter.
14. **Regularity/support control.** The result must lie in a stated function or
    distribution space for which NS-002 and the desired well-posedness claim
    make sense.

### Conditional on additional claims

15. **Additivity.** Required if independent sources or partitions are claimed
    to superpose; not otherwise frozen.
16. **Conservation/equivariance.** If derived from a coupled action, the map
    must satisfy the associated Noether identities and exchange accounting.
    BALANCE-001 does not provide a universal divergence-free stress-energy or
    conserved force by itself.
17. **Causality and hyperbolicity.** Required for a time-dependent projection;
    it cannot be assessed before kinetic/duration closure and cone matching.
18. **Potential integrability.** Required only if the load is claimed to arise
    from an interaction energy; a general load covector need not be exact.
19. **Universality across matter species.** Required only if a universal
    equivalence principle is separately imposed or derived. It is not frozen by
    the six authoritative inputs.

### Explicitly not generic requirements

- **Positivity:** a vector/covector body load has no coordinate-invariant
  positivity order. Only a separately declared energy or dissipation scalar
  can carry a sign condition.
- **Source-free conservation:** \(b\) is external supply to the regional
  placement balance; it need not vanish or be divergence-free. Conservation
  statements require a selected full action, symmetry, and boundary accounting.
- **Pointwise smoothness:** the frozen weak form permits a dual-space source;
  smoothness is a well-posedness choice, not the definition.
- **Linearity, isotropy, or dependence on energy density alone:** none is
  selected. A covariant map may depend on the full admissible stress-energy and
  available frozen structures.

## 8. Readiness assessment for SOURCE-PROJECTION-001

\[
\boxed{\text{Ready at the interface level; not pre-solved at the law-selection level.}}
\tag{NS-010}
\]

The next milestone now has a precise target: construct or classify maps into
\((T_y\mathcal Y_0)^*\), with a regular output measured as force per reference
volume when applicable. It must state the reference/effective index bridge,
normalization, support, regularity, species treatment, exchange accounting, and
weak-field/V11 matching.

What remains unavailable is exactly what SOURCE-PROJECTION-001 is meant to
decide: the formula, locality class, normalization, whether ordinary and
electromagnetic contributions share one rule, and whether the load is derived
from a coupled interaction action. Consequently no physical mass-energy
distribution can yet be converted uniquely into \(b\), but there is no longer
any ambiguity about the mathematical codomain or the role of the result in the
frozen native continuum.

## 9. Traceability

| Result | Frozen authority |
|---|---|
| objective dimensionless \(C=F^\sharp F\), placement realization | DEFORMATION-001 |
| \(W\), \(P_C=DW\), \(P_F=2FP_C\), dimensions of elastic response | CONSTITUTIVE-CONSTRUCTION-001; LOCALITY-001 |
| local divergence and separation of source from constitutive response | LOCALITY-001 L-001--L-002 |
| generalized balance/source accounting and conditional conservation | BALANCE-001 B-004--B-008 |
| weak form, prescribed \(b\), boundary loading, and missing \(T^{\rm matter}\mapsto b\) | WEAK-LENSING-LOCALITY-001 WL-006, WL-008, WL-012 |
| ontology, V11 gate, and no unauthorized constants | FOUNDATION-001 FP-1, FP-5, FP-6 |

