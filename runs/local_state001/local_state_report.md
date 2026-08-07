# PBUF LOCAL-STATE-001 — Local Computational State vs Global Ontological State

## 0. Decision

The symbol \(q\) has carried two different roles downstream of STATE-003:

1. **ontological role:** \(q_{\rm global}\in\mathcal Q_{\rm phys}\) is the
   complete instantaneous physical state of the one medium; and
2. **computational role:** a representative of the unknown fields actually
   needed on a selected problem domain.

Only the first role is fixed by STATE-003.  It does not follow from ontic
completeness that every functional must inspect every spatial part of its
argument.  For a bounded region \(\Omega\), the second role is represented by
the restriction

\[
 q_\Omega:=\rho_\Omega(q_{\rm global}),                         \tag{LS-001}
\]

together with admissible boundary, reference/background, gauge, source, and
observable data.  This is notation for less information about the already
authorized state, not a new state variable or a second ontology.

The frozen results prove that \(q_\Omega\) plus such data is sufficient for the
locally unique static native elastic deformation problem.  They do **not**
prove the same unconditional statement for the source projection or the
effective metric map: both remain support-undetermined.  Consequently:

\[
\boxed{\begin{aligned}
&\text{global ontology does not impose global computation;}\\
&\text{native static deformation is regionally sufficient;}\\
&\text{end-to-end local weak lensing is sufficient only conditional on}\\
&\text{support-controlled source and metric closures.}
\end{aligned}}                                                   \tag{LS-002}
\]

No statement here changes STATE-003, V11, ontology, evolution, cosmology, the
metric map, or any state-variable inventory.

## 1. Formal distinction

### 1.1 Ontological state

STATE-003 fixes

\[
 q_{\rm global}=[\widehat q]_{\mathcal G}\in
 \mathcal Q_{\rm phys}=\mathcal A/\mathcal G,                    \tag{LS-003}
\]

where \(q_{\rm global}\) separates every instantaneous physical distinction
of the complete one-medium universe.  Thus every instantaneous observable is
a gauge-invariant functional of \(q_{\rm global}\), or relationally of
\((q_{\rm global},q_{0,{\rm global}})\).  This is a statement about what
exists and what completely specifies the occupied order state.  It is not a
claim about the support of an operator, Cauchy sufficiency, or numerical input.

### 1.2 Computational state

Choose a bounded reference region \(\Omega\Subset\mathcal B_0\), a local
representative/gauge, and the induced regional gauge relation
\(\mathcal G_\Omega\).  Define

\[
 \mathcal Q_\Omega:=\rho_\Omega(\mathcal A)/\mathcal G_\Omega,
 \qquad q_\Omega=[\widehat q|_\Omega]_{\mathcal G_\Omega}.        \tag{LS-004}
\]

This definition needs one qualification.  A naive expression
\([\widehat q]_{\mathcal G}|_\Omega\) need not be representative-independent
when a global diffeomorphism moves \(\Omega\).  The restriction is therefore
defined only after \(\Omega\) is relationally identified or a regional
representative/gauge is chosen, and then quotiented by the transformations
preserving that regional identification.  Different regional representatives
of the same class contain no different physical information.

The **computational data package** for a static elastic solve is

\[
 \mathfrak C_\Omega=
 \{q_\Omega\text{ (unknown)},q_{0,\Omega},K_0,\mu_0,\mathcal D_C,
 b_\Omega,\bar y|_{\Gamma_D},\bar t|_{\Gamma_N},
 \text{gauge/rigid-mode fixing}\}.                               \tag{LS-005}
\]

The constants, domain, source, and boundary values in LS-005 are problem and
law data, not additional ontological state variables.  In the authorized
placement realization, \(q_\Omega\) is represented for this subproblem by
\(y|_\Omega\); this does not identify the complete abstract \(q\) with
placement.

## 2. Audit of the two meanings of \(q\)

The frozen texts do not contradict STATE-003, but several formulas use global
ontological notation at computational interfaces.  The following inventory
identifies every such place in the six authoritative milestones.

| Location | Use | Classification and correction |
|---|---|---|
| FOUNDATION-001 FP-4 | one complete configuration at each order state | purely ontological; no conflation |
| STATE-003 S3-001, S3-003--S3-004 and completeness theorem | complete gauge class and argument of every instantaneous observable | purely ontological; \(\mathcal O=\mathcal O[q]\) asserts factorization through the complete state, not full spatial support |
| STATE-003 S3-005 | \(q\mapsto C[q,q_0]\) | global notation for a relational observable; computationally replace by the restriction of \(C\) when its local realization is used |
| WEAK-LENSING-LOCALITY-001 information audit, WL-005, section 9 | explicitly separates \(q\) from \(q|_\Omega\) | no conflation; this is the controlling precedent for LS-001 |
| WEAK-LENSING-LOCALITY-001 \(g^{\rm eff}=G[q,C;\mathcal D]\), WL-011 | a computational map is written with global \(q\) | deliberate unresolved global/functional notation, not a necessity; it cannot yet be replaced unconditionally by \(q_\Omega\) because \(G\)'s support is not frozen |
| NATIVE-SOURCE-001 | the elastic solve uses \(b|_\Omega\), while \(T^{\rm matter}\mapsto b\) is missing | no direct \(q\)-conflation; it proves regional deformation sufficiency but leaves source support open |
| SOURCE-PROJECTION-001 | \(\Pi_{\rm source}[T]\) | no direct \(q\)-conflation; its unrestricted functional option is precisely why a local \(\Pi_\Omega\) cannot yet be asserted as frozen |
| MATTER-MEDIUM-INTERACTION-001 MMI-002, MMI-003 and completion theorem | matter is a distinction of complete \(q\); complete exchange cancels | genuinely ontological/global statements |
| MATTER-MEDIUM-INTERACTION-001 MMI-005 | \(\mathcal A_{\rm int}[y,\text{matter within }q]\) generates a load | ontological and computational roles share one symbol; the regional functional must instead state its support and use \(q_\Omega\) plus boundary/exterior data when local |
| MATTER-MEDIUM-INTERACTION-001 sections 4 and 7 | “every operand is a functional of \(q\)” is used beside a regional placement load | correct ontologically, overbroad if read as a computational requirement; replace the latter reading by factorization through regional data |
| MATTER-MEDIUM-INTERACTION-001 proposed universal work postulate | authorized functionals of \(q\) with locality/support still to be fixed | global notation is a placeholder; locality remains an explicit missing choice |

Thus the actual conflation occurs at the **notation-to-algorithm inference**:
writing \(F[q]\) because \(q\) is ontically complete has sometimes left the
impression that evaluation requires \(q_{\rm global}\).  It does not.  The
correct question is whether \(F\) factors through a restriction map.

## 3. Local restriction theorem

**Theorem 1 (regional restriction and elastic sufficiency).**  Let
\(\Omega\Subset\mathcal B_0\) be bounded and relationally identified.  Assume:

1. the authorized placement realization on \(\Omega\);
2. the frozen first-gradient energy and constitutive law;
3. prescribed \(b_\Omega\) and admissible Dirichlet/Neumann data;
4. a solution branch in a uniformly strongly elliptic subset of
   \(\operatorname{int}\mathcal D_C\); and
5. local uniqueness after gauge and rigid modes are fixed or quotiented.

Then the regional solution and all native elastic observables in \(\Omega\),

\[
 y_\Omega,\quad F_\Omega,\quad C_\Omega,
 \quad W_\Omega,\quad P_{C,\Omega},\quad P_{F,\Omega},             \tag{LS-006}
\]

factor through \(\mathfrak C_\Omega\).  No independent value of
\(q_{\rm global}|_{\mathcal B_0\setminus\overline\Omega}\) is required.

**Proof.**  The weak problem is

\[
 \int_\Omega P_F(y_\Omega):\operatorname{Grad}_0\eta\,dV_0
 =\langle b_\Omega,\eta\rangle
 +\int_{\Gamma_N}\bar t\cdot\eta\,dA_0                         \tag{LS-007}
\]

for all admissible zero-Dirichlet variations.  Every coefficient and operand
is in LS-005 or is a pointwise algebraic/first-derivative function of the
regional unknown.  Exterior field values do not occur.  Two global states
that induce the same LS-005 data therefore induce the same regional weak
problem; assumed uniqueness on the quotient makes their regional solutions
equal.  The quantities in LS-006 are local derived functionals of that
solution.  Hence they are constant on fibers of the regional-data map and
factor through \(\mathfrak C_\Omega\). \(\square\)

The theorem proves sufficiency for native static deformation, not that
\(q_\Omega\) alone is sufficient.  A boundary-value problem also needs its
law, reference, source, and admissible boundary data.  Nor does it give a
finite-radius pointwise domain of dependence: elliptic propagation couples a
point to all regional loads and the whole boundary.

## 4. Boundary sufficiency theorem

Let \(\partial\Omega=\Gamma_D\cup\Gamma_N\) up to negligible overlap.  The
minimal frozen boundary vocabulary is

\[
 y=\bar y\ \text{on }\Gamma_D,
 \qquad P_FN=\bar t\ \text{on }\Gamma_N.                         \tag{LS-008}
\]

An admissible package is any one of:

- sufficient Dirichlet data;
- mixed Dirichlet/traction data removing null modes; or
- pure traction data satisfying resultant force and moment compatibility,
  with solutions taken modulo rigid modes.

Gauge fixing or quotienting is also required.  At an ideal isolated exterior,
the permitted symbolic alternative is approach to the unloaded/background
class, \(C\to\mathbf1\), \(P_F\to0\), plus representative normalization.  No
falloff rate is frozen.  Finite traction-free, Robin, or absorbing boundaries
are not exact frozen consequences and require an explicitly declared
approximation or additional boundary law.

**Theorem 2 (exterior replacement by boundary data).**  Under Theorem 1's
assumptions, suppose two global configurations induce identical
\(q_{0,\Omega}\), \(b_\Omega\), and admissible traces LS-008.  Then their
locally unique native elastic solutions agree in \(\Omega\), modulo the
declared gauge/rigid equivalence.  Therefore all exterior influence relevant
to this elastic problem is replaceable by those induced bulk and boundary
data.

**Proof.**  The two restrictions satisfy the identical weak equation LS-007
on the identical admissible function space.  Local uniqueness on the quotient
implies equality. \(\square\)

“Minimal” here means minimal **types** of frozen boundary information, not a
claim that both \(\bar y\) and \(\bar t\) may be freely prescribed at every
boundary point.  Their numerical values are not generated by the frozen
theory.  Exterior effects may also be represented by \(b_{\rm ext}|_\Omega\)
or a declared background, but the same effect must not be double counted.

## 5. Reformulated matter–medium interaction

### 5.1 Ontological statement, unchanged

Matter is an organization or distinction of \(q_{\rm global}\), not a second
substance.  In a complete closed one-medium partition, internal exchange
contributions cancel:

\[
 \mathcal S_{\rm m\to e}[q_{\rm global}]
 +\mathcal S_{\rm e\to m}[q_{\rm global}]=0.                    \tag{LS-009}
\]

This statement genuinely concerns complete accounting and therefore retains
the complete ontological state.

### 5.2 Regional computational statement

For the elastic subsystem on \(\Omega\), the legitimate replacement is

\[
 \mathcal S_{{\rm int},\Omega}
 =\Pi_\Omega[q_\Omega;\beta_{\partial\Omega},\beta_{\rm bg}]
 \in(T_{y_\Omega}\mathcal Y_\Omega)^*,                          \tag{LS-010}
\]

where \(\beta_{\partial\Omega}\) denotes admissible exterior-induced boundary
data and \(\beta_{\rm bg}\) any separately declared background data.  Its
regular decomposition is

\[
 \langle\mathcal S_{{\rm int},\Omega},\eta\rangle
 =\int_\Omega b_\Omega\cdot\eta\,dV_0
 +\int_{\Gamma_N}\bar t\cdot\eta\,dA_0.                        \tag{LS-011}
\]

If an exact local interaction functional is later authorized, its regional
form is

\[
 \langle\Pi_\Omega,\eta\rangle
 =-D_{y_\Omega}\mathcal A_{{\rm int},\Omega}
 [y_\Omega,\text{matter}(q_\Omega);
 \beta_{\partial\Omega},\beta_{\rm bg}]\eta.
                                                                    \tag{LS-012}
\]

LS-010--LS-012 select no interaction law, normalization, differential order,
or metric map.  They state only its proper regional type.  Exactness remains
conditional, as in MATTER-MEDIUM-INTERACTION-001.

## 6. Source operator reformulation

The ontologically correct statement is that any physical source ultimately
depends on the complete state.  The frozen mathematical interface, however,
is \(\Pi_{\rm source}:T^{\rm matter}\mapsto
(T_y\mathcal Y_0)^*\), not a newly derived direct map from \(q\).  If a future
matter descriptor is established, completeness requires
\(T^{\rm matter}=T[q_{\rm global}]\), so the relevant expression is the
composite \(\Pi_{\rm source}\circ T\).  Neither writing
\(\Pi_{\rm source}[q]\) nor \(\Pi_{\rm source}[q_\Omega]\) changes that
canonical interface.  The computationally local statement is stronger and
requires a factorization theorem.  Let

\[
 R_\Omega(q_{\rm global})=
 (q_\Omega,\beta_{\partial\Omega},\beta_{\rm bg}) .              \tag{LS-013}
\]

Then a regional source operator exists exactly when

\[
 R_\Omega(q_1)=R_\Omega(q_2)
 \Longrightarrow
 (\Pi_{\rm source}\circ T)[q_1]|_\Omega=
 (\Pi_{\rm source}\circ T)[q_2]|_\Omega.                        \tag{LS-014}
\]

Equivalently, \(\Pi_{\rm source}|_\Omega\) must be constant on every fiber of
\(R_\Omega\).  In that case there is a unique induced map on the image,

\[
 \Pi_{{\rm source},\Omega}
 [T_\Omega;\beta_{\partial\Omega},\beta_{\rm bg}]
 :=\Pi_{\rm source}[T[q_{\rm global}]]|_\Omega,                 \tag{LS-015}
\]

and the right side is independent of the chosen global extension.

**Proof.**  If LS-015 exists, equal regional data give equal outputs, proving
LS-014.  Conversely, if LS-014 holds, define LS-015 using any global extension
of the regional data.  Fiber constancy makes the definition independent of
that extension. \(\square\)

Here \(T_\Omega\) denotes exactly the regional matter information on which the
selected projection depends; it is not assumed to be an independently added
state variable.  When \(T_\Omega=T[q_\Omega,q_{0,\Omega}]\) and LS-014 holds,
one may abbreviate the composite as
\(\Pi_{{\rm source},\Omega}[q_\Omega,q_{0,\Omega};\beta]\).

The frozen framework does **not** prove LS-014.  SOURCE-PROJECTION-001 permits
finite-jet, distributional, and functional/kernel rules and explicitly leaves
support unfixed.  A kernel can distinguish two global states agreeing on
\(\Omega\) and its boundary data but differing outside.  Hence it is not yet
fundamental or frozen that \(\Pi_{\rm source}=\Pi_{\rm source}[q_\Omega]\).
That notation is justified only after support-controlled locality or a bounded
kernel domain is selected.  More precisely, neither direct notation is the
fundamental frozen signature: the authorized signature is \(\Pi[T]\), while
its composition with the ontic state is regional only if LS-014 holds.
Moreover, because deformation is relational, a
source rule using relational invariants may also require \(q_{0,\Omega}\);
writing only \(q_\Omega\) must not suppress that authorized comparison datum.

## 7. Continuum comparison

| Theory | Ontological/global object | Regional computational closure | PBUF comparison |
|---|---|---|---|
| elasticity | one connected body's configuration | constitutive law, regional loads, and displacement/traction data | frozen PBUF static elasticity has the same first-gradient, elliptic bounded-domain structure |
| fluid mechanics | global fluid fields | local balance/constitutive laws plus initial and inflow/outflow/wall data | structurally analogous, but PBUF has no frozen kinetic or evolution closure |
| electromagnetism | one electromagnetic field | Maxwell equations with regional sources and constraint-compatible initial/boundary/radiation data | global field ontology does not imply global numerical input; PBUF source and metric closures are less complete |
| General Relativity | a global spacetime metric/matter solution | local covariant equations plus constraint-satisfying initial/boundary or asymptotic data and gauge | GR has a closed local source-to-metric equation; frozen PBUF retains V11 gates but has not selected its source projection or metric map |

PBUF is therefore a genuine local continuum theory **in its frozen native
static elastic sector**.  It is not yet a closed local continuum theory from
physical matter through lensing, because locality of two constitutive/interface
maps is unselected.  This is missing mathematical closure, not a conflict with
the one-medium ontology.

## 8. Weak-lensing dependency graph

```text
q_global : complete ontological state of the one medium
   |
   | restriction after regional identification / gauge handling
   v
q_Omega + q0_Omega + boundary/background data
   |
   | Pi_source,Omega                         [CONDITIONAL: support not frozen]
   v
native load b_Omega + traction data
   |
   | local static balance + frozen W         [FROZEN LOCAL]
   v
y_Omega -> F_Omega -> C_Omega -> P_F,Omega
   |
   | G_Omega                                 [CONDITIONAL: support not frozen]
   v
g_eff on an optical tube U_gamma
   |
   | V11-compatible null propagation         [LOCAL ON U_gamma]
   v
null ray / optical bundle
   |
   | source-observer endpoints, calibration,
   | and relevant background distance data
   v
weak-lensing observable
```

For an end-to-end bounded formulation, \(\Omega\) must contain the material
solve and the support required by \(\Pi\); the optical tube
\(U_\gamma\) must contain the relevant ray bundle and the support required by
\(G\).  These regions may be different and may be joined in one computational
domain.  If \(\Pi\) or \(G\) is finite-jet local, only the corresponding local
jets and regularity/traces are needed.  If either has a bounded or causal
kernel, its declared support must be included.  If unrestricted global support
is selected, complete global input is required by that selected operator—not
by STATE-003.

## 9. Complete-universe information: exact requirement audit

Complete-universe information is genuinely required only for:

1. stating which complete physical configuration is ontically occupied;
2. asserting that all matter belongs to the one medium and that no
   instantaneous physical distinction lies outside \(q\);
3. cancellation of every internal exchange in a literally complete closed
   one-medium accounting; and
4. evaluating a future operator if, and only if, its selected support is
   genuinely global, or computing an observable explicitly defined over the
   entire universe.

It is **not** intrinsically required for:

- native deformation, stress, or static balance in a bounded \(\Omega\), once
  regional source and admissible boundary/reference data are given;
- an interaction or source evaluation after it has been proven to factor
  through regional data;
- the effective metric after a local or support-controlled \(G\) is selected;
- photon propagation outside a neighborhood of the relevant null bundle; or
- distant matter whose influence is already represented without double
  counting by regional load, boundary, or background data.

Cosmological-distance lensing may require the background metric along the
source-observer path and its distance calibration.  Time-dependent lensing may
require local history/initial data and an authorized kinetic/duration closure.
Neither statement requires the complete instantaneous universe as such, and
no cosmology or evolution is derived here.

## 10. Completion theorem

**Theorem 3 (global ontology/local computation separation).**  Under exactly
the frozen milestones named in this mission:

1. \(q_{\rm global}\) remains the complete ontological state fixed by
   STATE-003;
2. its regional restriction \(q_\Omega\), properly defined modulo induced
   regional gauge, is an authorized computational restriction and not a new
   state variable;
3. \(q_\Omega\) plus LS-005 data is sufficient for every locally unique native
   static elastic prediction in \(\Omega\);
4. admissible boundary data replace the exterior for that problem;
5. a regional source operator and regional metric map exist precisely when
   their global counterparts are constant on the appropriate regional-data
   fibers; and
6. because those support properties are not frozen, complete local weak
   lensing is conditionally admissible but presently unclosed.

Therefore repeated use of \(q\) in downstream formulas is principally
ontological notation, not proof of mathematical dependence on the complete
universe.  The only unconditional global requirements are the ontological and
complete-accounting statements listed in section 9. \(\square\)

## 11. Traceability and status

| Result | Frozen authority |
|---|---|
| complete ontological state; no implication of realization or Cauchy completeness | STATE-003 S3-001--S3-004 |
| regional static elastic sufficiency and boundary vocabulary | WEAK-LENSING-LOCALITY-001 WL-005--WL-008 |
| local constitutive support but elliptic regional dependence | WEAK-LENSING-LOCALITY-001 WL-009--WL-010 |
| no frozen end-to-end locality because \(G\) may be nonlocal | WEAK-LENSING-LOCALITY-001 WL-011 and Theorem 3 |
| source codomain and unresolved locality/support | NATIVE-SOURCE-001; SOURCE-PROJECTION-001 sections 7--10 |
| one-medium internal work transfer and complete cancellation | MATTER-MEDIUM-INTERACTION-001 MMI-001--MMI-003 |
| exact interaction functional is conditional, not frozen | MATTER-MEDIUM-INTERACTION-001 MMI-005 and sections 7--9 |

**Status: complete.**  The requested separation is established without
modifying any frozen milestone.  The native elastic sector admits a bounded
computational state.  The whole weak-lensing pipeline admits the same form
exactly when the still-missing source and metric operators satisfy the stated
regional factorization/support conditions.
