# PBUF EQUILIBRIUM-MISMATCH-001 — Can Matter Load Spacetime by Local Equilibrium Mismatch?

## 0. Decision

**Outcome B: equilibrium mismatch is an admissible and potentially important
source mechanism, but the frozen framework neither makes it automatic nor
determines it from matter.**

The one-medium ontology permits different *derived or prescribed comparison
states* in different regions without adding a second medium.  STATE-003,
however, permits no independent preferred-state field outside the complete
state `q`.  Such a field must therefore be either:

1. fixed reference/problem data already included in the relational definition
   of deformation; or
2. a gauge-basic functional of `q` supplied by a law.

The first option describes prestrain but does not explain why matter creates
it.  The second option is precisely a matter-to-natural-state interaction law.
No frozen milestone supplies that functional, its normalization, or even the
claim that every matter-bearing configuration has a nonzero preferred-state
offset.

Moreover, a jump in preferred equilibrium does **not** necessarily produce
stress.  It does so only if the local stress-free states cannot be assembled
into one admissible continuous configuration, or if boundary/topological
constraints frustrate such an assembly.  When mismatch is incompatible, its
linearized virtual work is equivalent to a self-equilibrated generalized
native load.  At finite deformation it is a prestrain-dependent constitutive
problem, not in general a prescribed dead body load.

Thus equilibrium mismatch can *realize* the native source channel but cannot
derive `Pi_source` from matter.  The exact obstruction is the absent natural-
state selection map

\[
 \boxed{\mathcal N:\ q\ \hbox{(or an authorized matter descriptor }T[q])
       \longmapsto F_\natural[q]\ \hbox{or }q_0[q],}
 \tag{EM-001}
\]

together with its covariance, locality, compatibility class, and
normalization.  Naming matter a preferred configuration does not construct
`\mathcal N`; it restates the missing interaction in reference-state language.

No ontology, field, medium, particle, carrier, empirical coefficient, frozen
constitutive modulus, or V11 statement is changed below.

## 1. Equilibrium compatibility audit

### 1.1 One medium and local natural states

FP-1 rules out a second material substance, not spatially inhomogeneous states
or relational reference data for the same medium.  A regional comparison
`q_0^(A)` and `q_0^(B)` may therefore be used as mathematical data for two
patches of the same carrier.  This is ontologically harmless only if they are
restrictions of one admissible global comparison configuration, or are a
derived shorthand for properties of the occupied `q`.

The qualification is essential.  Under FP-4 and STATE-003, `q` is the complete
instantaneous physical state.  An independently specifiable field of local
rest states would distinguish physical states while `q` remained fixed and
would therefore be an unauthorized internal variable.  Consequently

\[
 q_0(X)\text{ may be reference/problem data},\qquad
 q_0^{\rm matter}(X)\text{ must be a declared functional of }q
 \tag{EM-002}
\]

if the latter is to be a physical matter property rather than an externally
prescribed model input.

There is a further type issue.  Frozen `q_0` denotes a complete comparison
configuration, not a pointwise value chosen independently at every `X`.
Arbitrary patching of local values need not define any admissible global
configuration.  Failure to patch is exactly incompatibility; it cannot at the
same time be treated as an already existing global stress-free `q_0`.

### 1.2 Compatibility with the frozen constitutive construction

The frozen response is one homogeneous, objective, local hyperelastic
functional

\[
 W(C),\quad W(\mathbf1)=0,\quad DW(\mathbf1)=0,
 \quad D^2W(\mathbf1)=\mathbb A_0>0.                    \tag{EM-003}
\]

A natural distortion `F_natural(X)` can be used without choosing a new formula
for `W` by evaluating the same frozen function on elastic relative deformation

\[
 F_e=F F_\natural^{-1},\qquad
 C_e=F_e^\sharp F_e,\qquad \Psi(X,F)=W(C_e).             \tag{EM-004}
\]

Here `F_natural` is notation for reference structure, not a newly authorized
physical field.  The local zero-stress wells are

\[
 \mathcal K(X)=\{R F_\natural(X):R^\sharp R=\mathbf1\}.  \tag{EM-005}
\]

Equation EM-004 preserves the formula, moduli, and stress-free point of the
frozen `W` in elastic coordinates.  It nevertheless changes the kinematic
input from the frozen single comparison map to a spatial natural-state map.
That change is admissible only conditionally under EM-002; it is not derived by
CONSTITUTIVE-CONSTRUCTION-001.  Treating matter instead as changing `W`, its
moduli, or its admissible domain would violate the frozen source/constitutive
separation and is not considered.

### 1.3 Local balance

For fixed authorized `F_natural`, variation with respect to the one placement
gives

\[
 P(X,F;F_\natural)=D_F\Psi(X,F),\qquad
 \int_\Omega P:\operatorname{Grad}_0\eta\,dV_0
 =\langle\mathcal S_{\rm ext},\eta\rangle .             \tag{EM-006}
\]

With no separately prescribed load, `S_ext=0`; in regular form
`-Div_0 P=0`, plus the applicable boundary conditions.  Spatial explicit
dependence through `F_natural` breaks material homogeneity but does not violate
local balance.  Its effect enters through the divergence of stress and
interface traction, using the already frozen one-medium load channel.

If `F_natural=N[q]` rather than fixed data, the full first variation also
contains

\[
 D_{F_\natural}\Psi\,D_q\mathcal N[q],\delta q.        \tag{EM-007}
\]

Holding `F_natural` fixed while varying the same complete `q` would omit this
term.  Therefore a state-dependent preferred configuration cannot consistently
generate loading until `N` and its variation are specified.  EM-007 is one
precise form of the missing interaction obstruction.

## 2. Compatibility-stress theorem

Let `Omega=A union B` with coherent interface `Gamma` and unit reference
normal `N`.  Let the two natural distortions be constant
`F_natural^A` and `F_natural^B`, and assume the frozen energy has exactly the
local stress-free wells EM-005 in the branch under study.

### Theorem 1 — Continuity is not sufficient to force stress

A piecewise homogeneous, continuous, stress-free placement exists across a
planar coherent interface if and only if there are admissible rotations
`R_A,R_B` and a vector `a` such that

\[
 R_BF_\natural^{B}-R_AF_\natural^{A}=a\otimes N.        \tag{EM-008}
\]

Subject to compatible exterior boundary data, EM-008 yields zero bulk stress
on both sides.  If no such rank-one connection exists, no continuous
piecewise-affine placement can occupy both local wells, so any admissible
placement must depart from at least one well (or develop a transition/singular
set).  Strict local stability then produces nonzero elastic stress wherever
that departure occurs.

**Proof.** Continuity of a piecewise-affine placement across `Gamma` is
equivalent to the Hadamard jump condition `[F]=a tensor N`.  Zero stress in the
selected stable branch requires `F_A=R_A F_natural^A` and
`F_B=R_B F_natural^B`.  Substitution gives EM-008.  Conversely EM-008 lets the
two affine maps be translated to agree on the interface, and both gradients
lie in their zero-stress wells.  If EM-008 fails, continuity and simultaneous
well membership are impossible.  Positive `A_0` converts a sufficiently small
departure from a well into nonzero stress and positive energy.  QED.

Hence the proposition “different neighboring preferred equilibria necessarily
generate stress” is false.  Equal natural states, rotation-related states,
rank-one-connected variants, suitable free boundaries, or a globally
compatible smooth natural distortion are counterexamples.  Conversely, two
different isotropic dilatations in three spatial dimensions generally fail:
`(alpha_B-alpha_A)I` has rank three rather than rank one when the dilatations
differ.

### Smooth compatibility

For a smooth natural distortion, stress-free assembly requires a placement
and rotation field satisfying

\[
 \operatorname{Grad}_0y=R(X)F_\natural(X).              \tag{EM-009}
\]

Equivalently, the natural metric `G_natural=F_natural^sharp F_natural` must be
locally realizable as a pullback metric, subject also to global topology and
boundary conditions.  In a Euclidean material realization this requires the
appropriate vanishing metric curvature; in a distortion representation it is
the corresponding curl/integrability condition, modulo rotation.  Nonzero
incompatibility is sufficient to forbid a global stress-free placement, but
its magnitude alone does not determine the realized stress: geometry,
boundary conditions, and the frozen nonlinear `W` also enter.

## 3. Mathematical derivation in the weak-field regime

Write `F=I+Grad u` and let the infinitesimal preferred strain be the symmetric
tensor `epsilon_natural`.  The frozen tangent `A_0` gives, without selecting a
new constitutive law,

\[
 \epsilon_e=\operatorname{sym}\operatorname{Grad}u-epsilon_\natural,
 \qquad \sigma=\mathbb A_0:\epsilon_e.                  \tag{EM-010}
\]

Stationarity in the absence of an independent load is

\[
 \int_\Omega \mathbb A_0:
 (\operatorname{sym}\operatorname{Grad}u-epsilon_\natural):
 \operatorname{sym}\operatorname{Grad}\eta\,dV_0=0.    \tag{EM-011}
\]

Thus

\[
 \int_\Omega (\mathbb A_0:
 \operatorname{sym}\operatorname{Grad}u):
 \operatorname{sym}\operatorname{Grad}\eta\,dV_0
 =\langle\Pi_{\rm mismatch},\eta\rangle,               \tag{EM-012}
\]

\[
 \boxed{\langle\Pi_{\rm mismatch},\eta\rangle
 =\int_\Omega P^*:\operatorname{Grad}\eta\,dV_0,
 \qquad P^*=\mathbb A_0:\epsilon_\natural.}            \tag{EM-013}
\]

For piecewise smooth `P*`, integration by parts gives the regular and singular
representatives

\[
 b_{\rm mismatch}=-\operatorname{Div}_0P^*,\qquad
 \bar t_{\rm mismatch}=P^*N,\qquad
 t_\Gamma=[P^*]N.                                      \tag{EM-014}
\]

The last term is the coherent-interface distribution with sign fixed by the
chosen interface normal.  EM-013, rather than only EM-014, is the invariant
statement and remains meaningful for discontinuous eigenstrain.

The source has no arbitrary new coefficient: its scale is set conditionally
by the magnitude of `epsilon_natural` and the already frozen elastic tangent.
But the magnitude of `epsilon_natural` as a function of matter is itself
undetermined; setting it is the missing coupling in another notation.

## 4. Load-equivalence theorem

### Theorem 2 — Conditional generalized-load equivalence

At frozen weak-field order, every prescribed preferred strain defines the
unique generalized placement load EM-013.  It is indistinguishable in the
linear displacement equation from an external native load having the same
dual action on all admissible test fields.  Its regular decomposition is
EM-014.

This equivalence is an operator rearrangement, not an identity of physical
origin.  The mismatch source has the restricted form

\[
 \Pi_{\rm mismatch}=B^*(\mathbb A_0:\epsilon_\natural),\qquad
 B\eta=\operatorname{sym}\operatorname{Grad}\eta,       \tag{EM-015}
\]

so it is a stress-polarization/eigenstress load.  In a closed body it is
self-equilibrated: it obeys the resultant force and moment compatibility
conditions when all bulk, interface, and boundary terms are included.  It
does not span arbitrary elements of `(T_y Y_0)^*` without additional boundary
or singular structure.

At finite deformation, define a homogeneous-reference stress `P_h(F)` and the
prestrained stress `P_n(F,F_natural)`.  One may formally write

\[
 \int P_h(F):\operatorname{Grad}\eta
 =\int [P_h(F)-P_n(F,F_\natural)]:\operatorname{Grad}\eta. \tag{EM-016}
\]

The right side generally depends on the unknown `F`.  It is therefore a
configuration-dependent (follower) generalized load, not a prescribed native
body density `b(X)`.  Only after linearization, or in a special constitutive
case where the polarization is independent of `F`, is it equivalent to a dead
external load.  Fundamentally, finite equilibrium mismatch belongs to the
prestrain/reference-structure class, while `b` belongs to the placement-dual
load class.  They intersect through EM-013 but are not identical classes.

## 5. Matter interpretation assessment

| Interpretation | Frozen-framework status | Reason |
|---|---|---|
| localized equilibrium defect | admissible description, not derived | it may mean incompatibility of a natural state of the one medium; no defect charge, topology, or matter identification is frozen |
| localized equilibrium offset | admissible prescribed/derived prestrain | it adds no medium, but an independent offset field would violate STATE-003 |
| localized preferred configuration | admissible only as reference data or `N[q]` | a complete local matter law must say which invariant of `q` selects it |
| all matter is equilibrium mismatch | unsupported | STATE-003 says matter is a distinction of `q`, not that every such distinction is visible in `C` or changes its stress-free state |
| matter is a separate inclusion/substance | forbidden | contradicts FP-1 and the one complete state |

The words “defect,” “offset,” and “preferred configuration” therefore do not
add ontology when used relationally.  They also do not supply physics by
themselves.  A physical interpretation adequate to replace the interaction
principle would need to establish all of the following from authorized `q`
data:

\[
 \text{matter distinction}\longmapsto F_\natural[q]
 \longmapsto C_e[q]\longmapsto P\longmapsto\Pi_{\rm mismatch}. \tag{EM-017}
\]

Only the last three conditional arrows are provided by frozen kinematics,
constitutive response, and balance.  The first arrow is absent.  Taking
`F_natural` as an independent material field would close the equations by
adding state information and is explicitly unavailable under STATE-003.

## 6. Weak-field and V11 compatibility

### 6.1 Localized deformation

A compact incompatible `epsilon_natural` gives a nonzero right side in
EM-012 and therefore generically a localized strained core or interface plus
an elastic response.  A compatible preferred strain
`epsilon_natural=sym Grad v` with compatible free boundary data can instead be
realized by `u=v` with zero stress.  Localization is therefore possible, not
necessary solely from nonuniformity.

### 6.2 Long-range response

The frozen linear operator is local elliptic elasticity.  On an unbounded
three-dimensional branch its Green tensor scales schematically as `1/r`.
Because a compact eigenstrain enters through the divergence of a compact
eigenstress, its zero-resultant far field is normally dipolar: schematically

\[
 u(x)\sim\nabla G(x):M=O(r^{-2}),\qquad
 \operatorname{Grad}u=O(r^{-3}),                         \tag{EM-018}
\]

unless boundary conditions, a noncompact source, a singular/topological
structure, or a nonzero resultant changes the leading multipole.  Thus
equilibrium mismatch naturally permits nonlocal *solution response* through a
local PDE, but it does not automatically supply an arbitrary monopolar native
body load or a universal gravitational far field.

### 6.3 V11 gate

The retained weak-field condition remains

\[
 [DG_{q_0}\,\mathcal L_0^{-1}
   D\Pi[\delta\text{matter}]]_{\rm gauge}
 =[h^{\rm V11}[\delta\text{matter}]]_{\rm gauge}.        \tag{EM-019}
\]

For mismatch,

\[
 D\Pi=B^*\mathbb A_0\,D\epsilon_\natural[\delta\text{matter}]. \tag{EM-020}
\]

Neither `D epsilon_natural` nor the effective metric map `DG` is frozen.
Consequently mismatch is structurally compatible with V11 but does not recover
V11 behavior, normalization, local Lorentz coupling, or lensing uniquely.
Its typical compact-source multipole in EM-018 further shows why long range by
itself is not a V11 recovery theorem.  No modification of V11 is required or
permitted.

## 7. Source reformulation

If a future authorized principle uniquely supplied `N`, the source could be
derived as the composite

\[
 \Pi_{\rm source}[q]
 =B^*\!\left(\mathbb A_0:epsilon_\natural[\mathcal N[q]]\right)
 \quad\text{at weak-field order},                        \tag{EM-021}
\]

with the nonlinear placement variation of EM-004 replacing EM-021 at finite
deformation.  In that conditional theory, `Pi_source` would not be an
independent primitive mapping; it would be the variational consequence of the
natural-state rule and the frozen energy.

In the actually frozen theory, EM-021 cannot be evaluated.  The exact
obstructions are:

1. **selection:** no invariant of `q` is identified as matter's natural-state
   offset;
2. **normalization:** no rule fixes the magnitude or sign of the offset from
   matter data, and choosing unit proportionality would still be a choice;
3. **tensor/index bridge:** no map fixes a rank-three natural distortion from
   an effective spacetime matter tensor;
4. **locality/support:** no finite-jet, boundary, topological, or nonlocal class
   for `N` is selected;
5. **compatibility:** no rule decides whether matter creates compatible
   remodeling (zero stress possible) or incompatible prestrain (residual
   stress forced);
6. **variation:** when the preferred state depends on the occupied `q`, the
   chain term EM-007 is undefined without `D N`;
7. **universality:** no theorem says all matter/radiation selects the same
   natural-state rule or selects a nonzero offset;
8. **V11 factorization:** the unfrozen metric map prevents the retained
   weak-field composite constraint from isolating `N`.

These are mathematical omissions, not requests for new ontology.  They show
that equilibrium mismatch relocates the missing interaction from
`matter -> b` to `matter -> preferred state`; it does not eliminate it.

## 8. Comparison with classical continuum mechanics

The comparison is mathematical only.

| Classical concept | Common mathematical structure | Difference/status in frozen PBUF |
|---|---|---|
| eigenstrain / transformation strain | `epsilon_e=epsilon-epsilon*`; incompatible `epsilon*` produces residual stress and an equivalent eigenstress load | exact weak-field analogue of EM-010–EM-015; PBUF lacks the rule assigning `epsilon*` to matter |
| thermal mismatch | spatially varying stress-free strain, often proportional to temperature; constraint or incompatibility produces stress | same operator structure, but PBUF has no authorized temperature field or expansion coefficient and none is introduced |
| growth/remodeling incompatibility | multiplicative split `F=F_e F_natural`; non-Euclidean natural metric may be unrealizable | finite-deformation analogue of EM-004 and EM-009; no biological mechanism or extra internal variable is implied |
| defects/non-Euclidean elasticity | curl/curvature or topology obstructs a global stress-free placement | useful compatibility criterion; PBUF has frozen no defect charge or microscopic defect ontology |
| residual stress | self-equilibrated stress remains with no applied dead load | the resulting PBUF stress is internal elastic response, representable as a generalized load only after a reference-side rearrangement |

The mechanism therefore belongs to an established constitutive class:
**inhomogeneous prestrain/eigenstrain (or non-Euclidean elasticity at finite
deformation)**.  Its use in a one-spacetime-medium ontology is a physical
interpretation of that class, not a mathematically distinct formulation.
What would be distinct and still missing is a PBUF-specific theorem deriving
the natural-state field from the complete medium configuration and satisfying
the V11 composite gate.

## 9. Final decision theorem

### Theorem 3 — Equilibrium mismatch cannot uniquely replace the interaction principle

Under FOUNDATION-001, STATE-003, DEFORMATION-001, HYPER-001,
ENERGY-PRINCIPLE-001, BALANCE-001, LOCALITY-001,
CONSTITUTIVE-CONSTRUCTION-001, WEAK-LENSING-LOCALITY-001,
LOCAL-STATE-001, NATIVE-SOURCE-001, SOURCE-PROJECTION-001, and
MATTER-MEDIUM-INTERACTION-001:

1. a prescribed or derived natural-state mismatch is compatible with the
   one-medium ontology and frozen local balance;
2. incompatibility or constraint, not mere regional difference, is the
   necessary condition for forced residual elastic stress;
3. its weak-field action is a native generalized load of the restricted
   eigenstress form EM-013, while finite mismatch is not generally a dead load;
4. the frozen corpus supplies no map from matter-bearing distinctions of `q`
   to the preferred state; and
5. neither the long-range multipole nor the V11 composite response is uniquely
   fixed.

Therefore Outcome A is false.  Outcome C is also too strong, because a
reference/prestrain description can be made without a new medium or a new
constitutive formula and can generate legitimate native loading.  The unique
supported decision is

\[
 \boxed{\textbf{Outcome B: equilibrium mismatch contributes but does not fully
 determine the source.}}                                      \tag{EM-022}
\]

The missing matter–medium interaction principle survives exactly as the
natural-state selection law EM-001.  `Pi_source` becomes derived only after
that law is supplied; it is not derived by the present frozen ontology alone.
