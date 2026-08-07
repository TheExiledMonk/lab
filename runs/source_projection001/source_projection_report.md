# PBUF SOURCE-PROJECTION-001 — Native Projection of Matter into Elastic Loading

## 0. Decision

**Outcome C, with the residual family characterized as Outcome B.**  The frozen
ontology does not determine a physical map

\[
 \Pi_{\rm source}:T^{\rm matter}\longrightarrow
 \mathcal S_{\rm ext}\in (T_y\mathcal Y_0)^* . \tag{SP-001}
\]

It determines the codomain, work pairing, units, transformation law, and the
place at which the result enters balance.  It does not determine the
interaction/work functional that converts matter data into that covector.
Consequently it selects none of rest mass, energy density, stress, momentum,
or the full stress-energy tensor as the unique loading property.  The minimal
missing physical ingredient is a **universal, normalized matter–medium
interaction principle** (equivalently, a gauge-basic interaction functional or
virtual-work rule).  Its first variation with respect to placement would define
the projection.  This is a coupling law, not a new field or constitutive law.

There is also a sharp tensor obstruction: on the isotropic unloaded branch,
no nonzero point-local parity-even objective map takes a symmetric rank-two
tensor to a spatial load covector without an additional vector or derivative.
Thus a nonzero rule based only on the frozen structures must either use spatial
variation (at least a first jet), use already-authorized configuration or
boundary data that supplies a direction, or remain a generalized/nonlocal
functional.  The frozen milestones choose none of these alternatives.

This report neither derives a metric map nor uses the non-frozen MATTER-001 or
GEOMETRY-001 conclusions as premises.  It changes no V11 statement.

## 1. What the frozen theory actually fixes

For the accepted placement realization,

\[
 F=\operatorname{Grad}_0y,\qquad C=F^\sharp F,\qquad
 P_C=D_CW,\qquad P_F=2FP_C,
\]

and the loaded weak problem is

\[
 \int_\Omega P_F:\operatorname{Grad}_0\eta\,dV_0
 =\langle\Pi_{\rm source}[T],\eta\rangle
 +\int_{\Gamma_N}\bar t\cdot\eta\,dA_0 . \tag{SP-002}
\]

When the source is regular,

\[
 \langle\Pi_{\rm source}[T],\eta\rangle
 =\int_\Omega b_i[T](X)\eta^i(X)\,dV_0,
 \qquad -\operatorname{Div}_0P_F=b[T]. \tag{SP-003}
\]

Therefore the native load is external virtual work conjugate to placement,
not stored energy, stress, prescribed deformation, or a modification of the
admissible domain.  Its only necessary physical description at present is the
resultant exchange of placement momentum/work between the matter accounting
and the medium accounting.  Calling that exchange “rest-mass loading” or
“energy loading” adds a principle not present in the frozen corpus.

The permitted dependency chain is exactly

```text
matter state represented by T
        |
        v  Pi_source: MISSING interaction/work rule
placement-load covector S_ext (regular representative b when available)
        |
        v  frozen balance, without alteration
placement y -> F=Grad_0 y -> C=F^sharp F
        |
        v  frozen constitutive evaluation W, P_C, P_F
```

The source does not pass through `W`, `P_C`, or `C` on its way into balance.

## 2. Loading-candidate audit

Let a stress-energy representative have components energy density
\(T_{00}\), momentum density \(T_{0i}\), and spatial stress \(T_{ij}\) in an
authorized local effective frame.  This split is representational; it is not a
new ontology.

| Proposed driver | What is frozen | Why it is not uniquely selected |
|---|---|---|
| rest mass | no independent native rest-mass density or rest congruence | does not cover radiation; a scalar cannot supply load direction or inverse length |
| total energy / energy density | admissible as part of stress-energy | still a scalar in a chosen frame; converting it to force density needs a direction and spatial derivative or a new inverse length |
| stress | has the right energy-density unit | a rank-two tensor is not a placement covector; divergence is one possible first-jet rule, not a frozen identity |
| momentum density | is frame/slicing dependent | mapping it directly to static force also requires a duration/rate rule and an index bridge |
| full stress-energy | is the most complete already-available matter descriptor | completeness of the input does not define a natural operator into the different codomain |
| deformation compatibility | constrains admissible \(y\) and \(C\) | is not matter data and cannot determine external virtual work |

Stress-energy is therefore a suitable *domain descriptor*, but the frozen
ontology does not make it a projection formula.  Which components do work on
which placement variations is exactly the absent interaction principle.

## 3. Point-local obstruction

**Lemma (algebraic isotropic obstruction).**  At the rank-three isotropic
reference state, a linear point-local parity-even objective operator

\[
 L:\operatorname{Sym}^2(V^*)\longrightarrow V^*
\]

is zero.

**Proof.**  Objectivity requires
\(L(RTR^\top)=R^{-\top}L(T)\) for every \(R\in SO(3)\).  In components a
linear map needs an invariant rank-three tensor \(A_i{}^{jk}\), symmetric in
\(j,k\).  Isotropic parity-even tensors are generated by the metric and have
even rank.  The only isotropic rank-three tensor is the parity-odd
Levi-Civita tensor, whose contraction with symmetric \(T_{jk}\) vanishes.
Hence \(A_i{}^{jk}=0\).  Equivalently,
\(\operatorname{Sym}^2(V)=\mathbf1\oplus\mathbf5\) contains no vector
representation. \(\square\)

The nonlinear conclusion at the isotropic point is the same for any rule made
only from the eigenvalues and spectral tensors of \(T\): these supply no
distinguished sign-sensitive vector.  A supplied normal, observer velocity,
material direction, acceleration, or second field would evade the lemma, but
none may be invented here.  An already-authorized nonuniform configuration or
boundary can supply a direction conditionally; that makes the rule
state/background dependent rather than uniquely matter-only.

A first derivative removes the type obstruction.  Schematically,

\[
 b_i=\nabla^j T_{ij},\qquad
 b_i=\nabla_i(\operatorname{tr}T) \tag{SP-004}
\]

both have the correct type and dimensions after an authorized pullback.  They
are witnesses of residual freedom, not proposed laws.  Neither their relative
weight, sign, applicable stress-energy contraction, nor even the decision to
use a derivative is frozen.  In addition, a conserved total stress-energy may
have vanishing covariant divergence on shell; interpreting a nonzero
divergence as exchange requires a declared coupled action and partition.

## 4. Complete dimensional chain

With \([y]=L\), \([F]=[C]=1\), and reference volume \([dV_0]=L^3\),

\[
 [W]=[P_C]=[P_F]=E L^{-3}=F L^{-2}.
\]

Therefore

\[
 [\operatorname{Div}_0P_F]=F L^{-3}=E L^{-4}=[b], \tag{SP-005}
\]

and

\[
 [b_i\eta^i dV_0]=(E L^{-4})(L)(L^3)=E, \tag{SP-006}
\]

as required by virtual work.  The full chain is

```text
stress-energy T                         [E L^-3]
  -> derivative/pullback or other missing interaction rule
native regular load b                   [E L^-4] = [F L^-3]
  -> -Div_0 P_F = b, with P_F           [E L^-3]
placement y                             [L]
  -> F = Grad_0 y, C = F^sharp F        [1]
```

This exposes the dimensional gap.  An algebraic use of mass density
\([M L^{-3}]\) needs both an acceleration scale; an algebraic use of energy
density \([E L^{-3}]\) needs an inverse length.  No such acceleration, length,
or coupling is frozen.  A derivative of stress-energy supplies \(L^{-1}\)
without a dimensional constant, but dimensional admissibility does not select
the contraction or establish its physics.  Multiplication or division by the
frozen elastic moduli cannot solve the problem uniquely: it merely forms
dimensionless ratios and still leaves a gradient/direction and functional form
undetermined.

For a spacetime tensor, a reference pullback/Jacobian and a spatial projection
are additionally necessary to turn effective-volume data and spacetime indices
into a covector-valued reference density.  DEFORMATION-001 leaves the concrete
rank-three versus rank-four realization open, so this bridge is not silently
identifiable.  The unit conversion and index conversion are distinct gaps.

No free normalization is introduced here.  Choosing coefficient `1` in one
candidate would itself choose a normalization, not prove it.

## 5. Constitutive and balance compatibility

For every admissible projection, its entire allowed mechanical action is the
right-hand side of SP-002.  Consequently:

1. \(W(C)\), \(P_C=D_CW\), and \(P_F=2FP_C\) remain unchanged;
2. \(C=F^\sharp F\), the placement class, and \(\mathcal D_C\) remain
   unchanged;
3. \(\Pi[T]\) must be continuous on the chosen test space, or be an admitted
   distribution;
4. pure Neumann data must satisfy resultant force and moment compatibility
   modulo rigid modes; and
5. a load for which no solution remains in \(\overline{\mathcal D_C}\) is not
   made admissible by altering the constitutive law.

These conditions restrict a proposed projection but do not construct one.

## 6. Universality

FP-1 fixes one medium and SP-002 fixes one resultant load slot.  Thus ordinary
matter, electromagnetism, and any other admissible source cannot create
independent displacement fields or separate media.  Their net external work
must land in the same dual space.

That is **codomain universality**, not **law universality**.  The frozen corpus
permits, without deciding among them,

\[
 \Pi_{\rm total}[T_{\rm m},T_{\rm EM}]
\]

as a joint nonlinear rule or, if independent additivity is separately
established, \(\Pi_{\rm m}[T_{\rm m}]+\Pi_{\rm EM}[T_{\rm EM}]\).  FP-1 does
not imply equal coupling, linearity, additivity, or an equivalence principle.
A universal rule for all stress-energy therefore requires the missing
universal interaction principle.  Species-specific rules are not ruled out by
the frozen milestones, provided they sum to one resultant load and do not add
new medium sectors; they are simply not derived either.

## 7. Weak-field/V11 audit

Linearizing the frozen elastic problem at the reference gives

\[
 \mathcal L_0 u=D\Pi_{T_0}[\delta T], \tag{SP-007}
\]

where \(\mathcal L_0\) is fixed by the frozen elastic tangent.  Retained V11
behavior constrains only the eventual composite response

\[
 [D G_{q_0}\,\mathcal L_0^{-1}D\Pi_{T_0}[\delta T]]_{\rm gauge}
 =[h^{\rm V11}[\delta T]]_{\rm gauge}. \tag{SP-008}
\]

Equation SP-008 is a compatibility condition, not a derivation: \(G\) is
explicitly unselected and its derivation is forbidden in this milestone.
Hence V11 cannot isolate \(D\Pi\) from the unfixed factorization.  Rescaling or
changing the source response can be compensated by a corresponding permitted
change in the unselected linear metric map, so long as the composite match is
retained.  No observational fitting or coefficient tuning can cure that
logical underdetermination.

The derived result therefore **does not naturally recover a unique V11
weak-field load**.  It preserves the V11 gate: any future projection must make
SP-008 possible, retain local Lorentz behavior, and leave V11 unchanged.

## 8. Locality classification

The only unconditional classification is

\[
 \Pi_{\rm source}[T]\in (T_y\mathcal Y_0)^*.
\]

It may have a smooth or distributional representative.  More specifically:

* a nonzero point-local algebraic matter-only rule is excluded on the
  isotropic reference branch by the lemma;
* first-jet and higher finite-jet rules are type- and dimensionally possible,
  but no derivative order or contraction is selected;
* discontinuous or concentrated matter makes distributional loads admissible
  through the frozen weak form;
* kernel/functional rules are not demanded by one-medium continuity and would
  require a support, measure, normalization, and usually a length or causal
  prescription not frozen here.

Thus the projection remains **locality-underdetermined**.  The native elastic
operator after \(b\) is supplied remains local first-gradient elasticity;
global dependence of an elliptic solution is PDE propagation, not evidence
that \(\Pi\) is a kernel.

## 9. Symmetry proof obligations

Let \(\varphi\) be a reference relabeling, \(R\) a superposed spatial frame
change, and \(\gamma\) an internal gauge transformation.  An admissible map
must satisfy the naturality conditions

\[
 \Pi[\varphi^*T]=\varphi^*\Pi[T],\qquad
 \Pi[R\cdot T]=R^{-\top}\Pi[T], \tag{SP-009}
\]

including the pulled-back reference volume density, and

\[
 \langle\Pi[\gamma\cdot T],\gamma_*\eta\rangle
 =\langle\Pi[T],\eta\rangle,
 \qquad \langle\Pi[T],\eta_{\rm gauge}\rangle=0. \tag{SP-010}
\]

Then \(P_F:\operatorname{Grad}_0\eta\), \(b\cdot\eta\), and the boundary
traction pairing transform as scalar densities, so SP-002 is covariant and
objective.  For electromagnetic matter, gauge compatibility further requires

\[
 \Pi[T(A+d\chi)]=\Pi[T(A)] \tag{SP-011}
\]

whenever the input is represented by a gauge potential; a rule expressed only
through gauge-invariant stress-energy satisfies this condition.  These are
necessary equivariance identities.  They do not select a unique natural
operator: the distinct first-jet witnesses in SP-004 obey the same tensorial
transformation rules.

Balance is preserved because the output appears only in SP-002/SP-003.
One-medium ontology is preserved because every contribution lands in the same
placement dual and internal partition exchanges must cancel.  No additional
field equation, force channel, or metric is introduced.

## 10. Uniqueness theorem

**Theorem (SOURCE-PROJECTION-001).**  Assume exactly the frozen milestones
listed in the mission.  They determine the target space, virtual-work role,
dimensions, balance placement, and covariance/objectivity/gauge constraints of
\(\Pi_{\rm source}\).  They do not determine a nonzero physical projection
from matter stress-energy to native elastic loading.  In particular:

1. no nonzero algebraic isotropic matter-only projection exists at the
   unloaded rank-three branch;
2. finite-jet constructions become possible once derivatives are admitted,
   but more than one invariant contraction has the required type and units;
3. no frozen premise chooses derivative order, normalization, support,
   additivity, matter-species universality, or the spacetime-to-reference index
   bridge; and
4. weak-field V11 compatibility constrains only a composite containing the
   independently unselected metric map and therefore cannot remove this
   freedom.

Accordingly Outcome A is false.  The admissible family is reduced as stated in
sections 5, 8, and 9 (Outcome B as a mathematical classification), while a
physical choice requires Outcome C.

**Proof.**  SP-002–SP-003 establish the frozen codomain and placement in
balance.  SP-005–SP-006 establish the unit requirement.  The isotropic lemma
excludes a nonzero point-local linear projection without extra structure.
SP-004 gives inequivalent finite-jet witnesses satisfying the same type and
unit constraints, so those constraints cannot imply uniqueness.  No cited
frozen premise orders or equates those witnesses, supplies their interaction
normalization, or imposes universal species coupling.  Finally SP-008 contains
both \(D\Pi\) and unfixed \(DG\), so the retained weak-field output cannot
identify either factor separately. \(\square\)

### Exact remaining mathematical freedom

The freedom is not an arbitrary new ontology.  It is the choice of a natural,
gauge-basic operator

\[
 \Pi:\mathcal T_{\rm admissible}\to(T_y\mathcal Y_0)^*
\]

including: the authorized index/volume bridge; differential or functional
order and support; which invariant contractions of total or partitioned
stress-energy enter; normalization and sign; linearity/additivity; regularity;
and whether the covector is exact (derived from interaction energy) or a
general external-work form.

### Minimal missing ingredient

Supply one normalized, universal virtual-work principle, preferably in the
form of an already-variable interaction functional

\[
 \mathcal A_{\rm int}[y,\text{matter}],\qquad
 \boxed{\langle\Pi[T],\eta\rangle
       =-D_y\mathcal A_{\rm int}[y,\text{matter}]\,\eta}. \tag{SP-012}
\]

The principle must use only authorized variables and must state its pullback,
locality, normalization, and matter universality.  SP-012 would automatically
give the correct covector and exchange-work interpretation; its invariances
would supply the relevant Noether compatibility.  The frozen ontology does
not provide \(\mathcal A_{\rm int}\), so writing SP-012 is the exact closure
request, not a claimed derivation.

## 11. Completion statement

The milestone is complete at Outcome C.  One continuous elastic spacetime
medium necessitates one resultant placement-loading channel, but it does not
necessitate what property of matter excites that channel or how strongly.
Stress-energy is the admissible comprehensive input description; external
virtual work is the necessary output meaning.  The physical arrow between
them remains a coupling law and cannot be inferred from ontology, elasticity,
dimensional analysis, symmetry, locality, or the retained weak-field gate
alone.

