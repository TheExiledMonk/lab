# PBUF GOVERNING-EQUATION-001 — Native Governing Equation Family

## 0. Decision, scope, and completeness boundary

The frozen constitutive construction determines the complete **internal elastic
operator**, but it does not determine a unique time-evolution equation.  For
Candidate A the native equation family is

\[
 \boxed{\mathcal K_\tau[q]+(D_qC[q,q_0])^*P_C(C)+N_{\mathcal K}(q)
       \ni \mathcal S_\tau,}\qquad
 \boxed{P_C={K_0\over2}\operatorname{tr}E\,{\bf1}+\mu _0E_{\rm TF}},
 \tag{GE-001}
\]

where

\[
 E={1\over2}(C-{\bf1}),\qquad C\in\mathcal K
 :=\{q:C[q,q_0]\in\overline{\mathcal D_C}\}.                 \tag{GE-002}
\]

Here \(N_{\mathcal K}\) is the variational normal cone.  It vanishes in the
interior and records the already-frozen hard admissibility boundary; it is not
a new field or constitutive response.  \(\mathcal K_\tau\) is the still-open
kinetic/inertial covector and \(\mathcal S_\tau\) the still-open source
covector, both evaluated after an emergent-duration gauge is selected.

Thus GE-001 is the strongest complete statement authorized by the frozen
architecture.  Replacing \(\mathcal K_\tau\) by a guessed mass density or
replacing \(\mathcal S_\tau\) by a guessed gravity source would not be a
derivation.  Candidate B has exactly GE-001, with the same normal-cone
admissibility, but places the constraint in the state space rather than in an
extended-valued energy.  No V11 or weak-lensing equation is changed.

## 1. Variables, domains, operators, and assumptions

### 1.1 Variables

| Kind | Symbol | Meaning and status |
|---|---|---|
| independent | \(X\in\mathcal B_0\subset\mathbb R^3\) | reference-carrier point in the authorized rank-three local placement realization |
| independent | \(\tau\in I\) | calibrated emergent duration; before calibration use an arbitrary ordered label \(s\) |
| dependent | \(y(X,\tau)\) | placement representative of the complete state \(q(\tau)\); gauge-equivalent representatives describe the same state |
| derived | \(F=\operatorname{Grad}_0y\) | deformation gradient |
| derived | \(C=F^\sharp F\) | objective dimensionless relative deformation |
| derived | \(E=(C-{\bf1})/2\), \(t=\operatorname{tr}E\), \(E_{\rm TF}=E-t{\bf1}/3\) | strain coordinates, not new state variables |
| constitutive | \(P_C=DW_A(C)\), \(P_F=2FP_C\) | stresses conjugate to \(C\) and \(F\) |
| closure slots | \(\mathcal K_\tau[y]\), \(b\), \(\bar t\) | kinetic covector, body source, boundary traction; not constitutive data |
| operationally derived | \(g^{\rm eff}=G[q,C;\mathcal D]\) | effective metric for a selected \(G\in\mathfrak G\); \(G\) remains open |

The intrinsic equations GE-001 and the weak formulation below do not require
the placement realization.  The explicit divergence form does.  A rank-four
clock/ruler realization was not selected by the frozen rank-three
constitutive construction and is therefore not silently substituted here.

### 1.2 Operators

* \(\operatorname{Grad}_0\) and \(\operatorname{Div}_0\) are the material
  gradient and row-wise material divergence associated with the fixed
  reference carrier and measure \(dV_0\).
* \(\operatorname{tr}\), \((\cdot)_{\rm TF}\), \(:\), and \(\sharp\) are,
  respectively, trace, trace-free projection, reference inner product, and
  reference adjoint.
* \(D_qC\) is the Fréchet/Gâteaux derivative of the kinematic map and
  \((D_qC)^*\) its variational adjoint.
* \(N_{\mathcal K}(q)=\{\zeta:\langle\zeta,z-q\rangle\leq0\ \forall
  z\in\mathcal K\}\) is written for a convex admissible set.  For the frozen
  domain, which need only be path-connected, the correct object is the
  corresponding variational/limiting normal cone or, equivalently, the
  tangent-cone inequality GE-009 below.  No convexity is inferred.
* \(\mathcal K_\tau\) is an operator slot, not the scalar modulus \(K_0\).
  Its differential order and normalization await kinetic closure.

### 1.3 Assumptions actually used

1. \(q_0\) is fixed and unloaded; \(C[q_0,q_0]={\bf1}\).
2. The local rank-three placement realization is used for explicit strong
   and boundary forms: \(F=\operatorname{Grad}_0y\), \(C=F^\sharp F\).
3. \(C(X,\tau)\in\overline{\mathcal D_C}\Subset\operatorname{Sym}^+(3)\).
4. In the interior, \(y\) is regular enough for the displayed derivatives and
   integration by parts.  The weak form requires only the corresponding
   Sobolev regularity and integrability.
5. Candidate A is the extended-valued quadratic CC-002, with the already
   frozen \(K_0>0\) and \(\mu _0>0\).  There is no higher-gradient,
   nonlocal, rate, history, or dissipative constitutive term.
6. Sources and boundary loading act through prescribed dual pairings.  No
   body-source formula is assumed.
7. Evolution equations using \(\tau\) are conditional on duration calibration
   and a kinetic closure.  Before that selection, only the static variational
   equation and the native degree-one history-action family are fixed.

## 2. Constitutive and strong forms

The stored energy per reference volume is

\[
 Q(C)={K_0\over2}t^2+\mu _0E_{\rm TF}:E_{\rm TF},\qquad
 W_A=Q+I_{\overline{\mathcal D_C}},                              \tag{GE-003}
\]

where \(I\) is zero on the admissible set and \(+\infty\) outside.  In the
interior,

\[
 P_C={K_0\over2}t{\bf1}+\mu _0E_{\rm TF},\qquad
 P_F=2F P_C
 =F\bigl(K_0t{\bf1}+2\mu _0E_{\rm TF}\bigr).                    \tag{GE-004}
\]

The sign-safe intrinsic elastic force is the variational derivative

\[
 \mathcal A[y]:=D_y\int_{\mathcal B_0}W_A(C[y])\,dV_0
 =(D_yC)^*P_C+N_{\mathcal K}(y).                                 \tag{GE-005}
\]

On a smooth interior placement, \((D_yC)^*P_C=-\operatorname{Div}_0P_F\).
Consequently the strong interior family is

\[
 \boxed{\mathcal K_\tau[y]-\operatorname{Div}_0\!\left[
 F\bigl(K_0t{\bf1}+2\mu _0E_{\rm TF}\bigr)\right]=b}
 \quad\hbox{in }\mathcal B_0\times I.                            \tag{GE-006}
\]

This convention defines \(b\) as the source on the right and the internal
force as the energy gradient.  If BALANCE-001/LOCALITY-001's symbolic
\(\mathcal K_{\rm frozen}+\operatorname{Div}P_F=b\) convention is retained,
its kinetic symbol is the negative of the inertial covector used in GE-006;
the variational content is identical.  Writing GE-005 prevents a sign
convention from becoming physics.

Static equilibrium is the closed constitutive-balance equation

\[
 -\operatorname{Div}_0P_F=b                                      \tag{GE-007}
\]

in the admissible interior.  A classical evolution is not closed until
\(\mathcal K_\tau\) is selected.  For example only, if a future authorized
closure establishes \(\mathcal K_\tau[y]=\rho_0\partial_\tau^2y\), GE-006
becomes elastodynamics.  The symbol \(\rho_0\) in that conditional statement
is closure data, not a coefficient derived or introduced by this milestone.

## 3. Weak variational formulation

Let \(V\) be an admissible placement space, let \(V_D\) impose prescribed
Dirichlet trace on \(\Gamma_D\), and let \(V_0\) be its zero-trace variation
space.  For an interior solution, find \(y(\tau)\in V_D\) such that for every
\(\eta\in V_0\),

\[
 \boxed{\langle\mathcal K_\tau[y],\eta\rangle
 +\int_{\mathcal B_0}P_F(y):\operatorname{Grad}_0\eta\,dV_0
 =\langle b,\eta\rangle
 +\int_{\Gamma_N}\bar t\cdot\eta\,dA_0.}                        \tag{GE-008}
\]

It follows exactly by varying
\(\mathcal E[y]=\int W_A(C[y])dV_0\), using
\(\delta C=F^\sharp\operatorname{Grad}_0\eta+
(\operatorname{Grad}_0\eta)^\sharp F\), hence
\(P_C:\delta C=P_F:\operatorname{Grad}_0\eta\), and integrating by
parts.  The right side is external virtual work.  No source potential is
required unless conservative loading is separately declared.

At a hard endpoint the correct statement is the variational inequality: find
\(y\in V_D\cap\mathcal K\) such that, for every admissible
\(z\in V_D\cap\mathcal K\),

\[
 \boxed{\langle\mathcal K_\tau[y]-\mathcal S_\tau,z-y\rangle
 +\int_{\mathcal B_0}P_F(y):\operatorname{Grad}_0(z-y)\,dV_0\geq0.} \tag{GE-009}
\]

For Candidate A this is the subdifferential Euler condition of the
extended-valued energy.  For Candidate B it is the first-order condition for
the separately constrained state space.  Thus their smooth interior and
endpoint trajectories coincide when the same admissible set and kinetic/source
closures are used, although their extended constitutive bookkeeping differs.

If an authorized clock-gauge action
\(S_\tau=\int(T-W_A)+S_{\rm ext}\,d\tau\) is later selected, GE-008 is its
Euler equation with \(\mathcal K_\tau\) equal to the kinetic variation.  Before
that selection, GE-008 with \(\mathcal K_\tau=0\) is the fully derived E1
static weak form; the dynamic term is a closure slot inherited from E2.

## 4. Balance form

Define a generalized material momentum density \(p_\tau\) only after a kinetic
closure supplies it.  If that closure has the local balance representation
\(\mathcal K_\tau=\partial_\tau p_\tau\), GE-006 reads

\[
 \boxed{\partial_\tau p_\tau+\operatorname{Div}_0J_p=\sigma_p,\qquad
 J_p=-P_F,\quad\sigma_p=b.}                                      \tag{GE-010}
\]

Equivalently, with the conventional momentum flux \(P_F\), one writes
\(\partial_\tau p_\tau-\operatorname{Div}_0P_F=b\).  The classification is:

| Item | Status |
|---|---|
| \(p_\tau\) | potentially conserved density, but undefined until kinetic closure |
| \(P_F=2FP_C\) | completely fixed constitutive stress/flux contribution |
| \(b\) | external/body source, not fixed by constitutive construction |
| \(\bar t=P_FN\) | boundary flux/loading |
| \(N_{\mathcal K}\) | ideal constraint reaction; it enforces admissibility and does no work in admissible tangent directions |

With \(b=0\) and zero net traction, integrated momentum is conserved only if
the selected full action actually has the required translation symmetry.
BALANCE-001 proves that this is not unconditional.

The only unconditional energy statement remains the storage chain rule

\[
 \partial_\tau W_A=P_C:\partial_\tau C                   \tag{GE-011}
\]

on a smooth elastic branch.  A total-energy conservation law additionally
needs kinetic energy, power/source pairing, clock-translation symmetry, and
zero boundary supply.  None is constitutive.  In an effective representation,
any selected balance is \(\nabla^{\rm eff}_\mu J^\mu=\Sigma\); neither current
nor source is generated by \(W_A\) alone.

## 5. Boundary and endpoint conditions

Let \(\partial\mathcal B_0=\overline{\Gamma_D}\cup
\overline{\Gamma_N}\) with disjoint interiors.

1. **Dirichlet (essential):** \(y=\bar y\) on \(\Gamma_D\), with
   \(C[\bar y]\) admitting an extension in \(\overline{\mathcal D_C}\).
   Variations vanish there.
2. **Neumann (natural):** \(P_FN=\bar t\) on \(\Gamma_N\).  The traction-free
   case is \(P_FN=0\).  These and only these spatial natural data follow from
   the local first-gradient energy.
3. **Mixed:** Dirichlet conditions on \(\Gamma_D\) and Neumann conditions on
   \(\Gamma_N\).  A Robin law on the same boundary point would require a
   boundary energy or boundary constitutive law and is not derived.
4. **History endpoints:** fixed \(q(\tau_-)\), \(q(\tau_+)\) make history
   variations vanish.  Free temporal endpoints would yield the natural
   transversality condition supplied by a selected kinetic action; because
   that action is open, no explicit momentum endpoint value is frozen.
5. **Elastic-domain endpoint:** at
   \(C\in\partial\mathcal D_C\), admissible variations lie in the tangent cone
   and GE-009 holds, equivalently a normal reaction appears in GE-001.  No
   prescribed limiting stretch, multiplier field, rebound rule, or boundary
   stress law follows.

Items 1–4 are identical for Candidates A and B.  Item 5 is encoded by
\(\partial W_A\) for A and by the separate state constraint for B; it is the
same variational inequality but not the same extended-energy definition.

## 6. Initial-value structure

A specified evolution problem requires the following distinct data.

| Class | Required data | Frozen/open status |
|---|---|---|
| constitutive | \(K_0,\mu_0\), \(W_A\), \(\mathcal D_C\), \(q_0\), and the placement kinematics | frozen (the domain is frozen data, not newly selected here) |
| kinetic | operator \(\mathcal K_\tau\), its normalization, constraints, and duration gauge | open closure |
| initial | \(q(\tau_0)=q_i\); plus as many independent kinetic data as the selected temporal order requires (e.g. velocity for a nondegenerate second-order closure) | problem data, closure-dependent |
| boundary | admissible \(\bar y\), \(\bar t\), boundary partition, and compatibility at \(\tau_0\) | problem data |
| source | body/source covector, boundary supply, and any native-to-effective source projection | open/problem data |
| metric | one \(G\in\mathfrak G\), clock/ruler calibration, and cone match | open closure |

Initial data must satisfy gauge and kinetic constraints, the Dirichlet trace,
\(C[q_i,q_0]\in\overline{\mathcal D_C}\), and any compatibility conditions
at the intersection of initial and spatial boundaries.  At a hard state-space
boundary, the initial tangent must lie in the admissible tangent cone.  It is
not possible to state a universal number of initial derivatives before the
temporal order of \(\mathcal K_\tau\) is chosen.

## 7. Reference linearization and waves

Let \(y(X,\tau)=X+u(X,\tau)\).  To first order,

\[
 F={\bf1}+\nabla u,\qquad E=\varepsilon=\operatorname{sym}\nabla u,
\qquad
 \delta P_F=K_0\operatorname{tr}\varepsilon\,{\bf1}
                 +2\mu_0\varepsilon_{\rm TF}.                    \tag{GE-012}
\]

Thus the frozen weak-field constitutive response is recovered exactly:

\[
 \boxed{\sigma_0=2\mu_0\varepsilon+
 \left(K_0-{2\mu_0\over3}\right)(\nabla\!\cdot u){\bf1}.}       \tag{GE-013}
\]

Writing \(\mathcal K_0\) and \(f\) for the linearized kinetic and source
closures gives the native linear family

\[
 \boxed{\mathcal K_0[u]-\mu_0\Delta u-
 \left(K_0+{\mu_0\over3}\right)\nabla(\nabla\!\cdot u)=f.}       \tag{GE-014}
\]

For the conditional local inertial closure
\(\mathcal K_0=\rho_0\partial_\tau^2\), Helmholtz decomposition gives

\[
 \rho_0\partial_\tau^2u_T-\mu_0\Delta u_T=f_T,\qquad
 \rho_0\partial_\tau^2u_L-\left(K_0+{4\mu_0\over3}\right)
 \Delta u_L=f_L,                                                 \tag{GE-015}
\]

and the acoustic speeds would be
\(c_T^2=\mu_0/\rho_0\) and
\(c_L^2=(K_0+4\mu_0/3)/\rho_0\).  GE-015 demonstrates the local wave
equations promised by LOCALITY-001 without adopting \(\rho_0\) as a frozen
constant.  Positivity of \(K_0,\mu_0\), together with positive inertia, makes
both sectors hyperbolic.

V11 operational compatibility is precisely the pair of gates

\[
 [D_qG|_{q_0}\,\delta q]_{\rm gauge}=[h^{\rm V11}]_{\rm gauge},
 \qquad
 \operatorname{Char}(\mathcal K_0+\mathcal A_0)
 =\{k:g^{\rm eff\,\mu\nu}_0k_\mu k_\nu=0\}                       \tag{GE-016}
\]

for the universal V11 signal sector.  GE-012–GE-014 recover the required
symmetric weak perturbation and frozen elastic tangent.  Actual equality of
the characteristic cone, polarization content, and normalization cannot be
proved until \(\mathcal K_0\) and \(G\) are selected.  This conditional result
is the complete honest demonstration of V11 compatibility; asserting more
would resolve two explicitly open closures by assumption.

## 8. Mathematical classification

| Sector | Classification |
|---|---|
| constitutive map \(C\mapsto P_C\) | affine linear in \(C\) in the interior; isotropic and hyperelastic |
| placement operator | second order in space and quasilinear/geometrically nonlinear through \(C=F^\sharp F\) and \(P_F=2FP_C\); no higher spatial order |
| static smooth interior | nonlinear elliptic system wherever the lifted Legendre–Hadamard/acoustic form CC-006 is positive |
| hard endpoint | nonlinear variational inequality/differential inclusion with a state constraint |
| dynamic smooth interior | hyperbolic only for a positive, nondegenerate, second-order inertial closure; other kinetic choices can change temporal order and PDE type |
| linearized static | strongly elliptic Navier system because \(K_0>0,\mu_0>0\) |
| linearized inertial | longitudinal/transverse hyperbolic wave system conditional on positive inertia |
| mixed formulations | displacement–stress or displacement–constraint formulations are equivalent representations; an endpoint reaction may be represented in the dual normal cone, but no new physical multiplier field is required |
| variational character | E1 is the Euler/subdifferential equation of \(\mathcal E\); E2 is variational only after a kinetic action is selected |

Finite-deformation strong ellipticity is guaranteed at the reference and only
on the declared connected positive acoustic component thereafter.  The
quadratic dependence on \(C\) does not make the placement problem globally
linear or globally elliptic.

## 9. Well-posedness requirements

### 9.1 Frozen contributions

* \(K_0>0\) and \(\mu_0>0\) give positive reference tangent, local coercivity
  modulo gauge/rigid modes, and reference strong ellipticity.
* \(W_A\geq0\), lower semicontinuity of its extended completion, and the
  compactly contained SPD admissible closure provide the frozen energy/domain
  controls.
* Objectivity creates gauge/rigid null directions; uniqueness can only be
  stated after essential boundary conditions or quotienting removes them.
* Local first-gradient structure supplies the correct weak space and traction
  pairing without extra boundary data.

### 9.2 Additional mathematical/problem requirements

Existence of static weak solutions requires a nonempty admissible class,
compatible boundary data, sequential weak compactness/coercivity, lower
semicontinuity in the chosen placement topology, and a load continuous in the
dual topology.  The frozen compact state domain does not by itself prove weak
closure of gradients or quasiconvexity of \(W_A(F^\sharp F)\); these must be
proved on the selected realization/operating domain.

Uniqueness requires strict monotonicity/convexity on the admissible quotient or
a suitable local inverse theorem.  Reference positivity supplies local, not
global, uniqueness.  Multiple finite-deformation equilibria are not excluded.

Continuous dependence requires uniform coercivity/strong ellipticity,
Lipschitz or adequate continuous dependence of the operator on state, stable
boundary conditions, and source control.  Approaching loss of ellipticity or a
nonsmooth hard boundary can invalidate classical estimates and requires
variational-inequality stability theory.

Evolution existence and uniqueness additionally require a selected kinetic
operator with positive energy, a well-defined constrained phase space,
regularity compatible with the elastic operator, and source regularity.
Admissible propagation requires real characteristic roots, finite propagation
speed, constraint preservation, acoustic positivity throughout the operating
domain, and agreement of the universal signal cone with the selected
\(g^{\rm eff}\).  All but reference acoustic positivity are closure- or
solution-domain requirements, not further constitutive uncertainties.

## 10. Remaining intentionally open closures

| Open component | What must eventually be supplied | Why it is not constitutive uncertainty |
|---|---|---|
| kinetic/inertial closure | \(\mathcal K_\tau\), temporal order, momentum map, normalization, possible constraints | DYNAMICS-001 explicitly leaves the action/kinetic family open; \(W_A\) is rate-free |
| duration gauge | calibration from ordered histories to \(\tau\) | DURATION-001 separates emergent duration from the arbitrary order label |
| source specification | \(b\), boundary work, balance production, native-to-effective source projection | BALANCE-001 separates sources from stress and derives no universal source density |
| effective metric map | one \(G\in\mathfrak G\), its regularity and linear normalization | METRIC-001 defines compatibility gates, not a unique map |
| concrete global realization | function spaces, gauge fixing, reference/event correspondence, domain regularity | STATE/DEFORMATION fix the abstract state and \(C\), not all analytic realization data |
| loading/IBVP data | initial state/tangent, boundary partition and values | problem data, not material response |
| global propagation domain | proof of CC-006 positivity along the intended branch | Candidate A fixes the formula but reference positivity does not imply global strong ellipticity |

There is no remaining higher-order constitutive function, fitted coefficient,
communication operator, barrier profile, microscopic variable, or weak-lensing
term in the canonical baseline.  Candidate A has closed those constitutive
choices.  The table lists orthogonal kinetic, source, metric, realization, and
problem-data slots already separated by the frozen architecture.

## 11. Complete term traceability audit

| Equation term/object | Originating frozen milestone(s) | Derivation route |
|---|---|---|
| \(q,q_0,\mathcal Q_{\rm adm}\) | FOUNDATION-001; STATE-002 | one complete medium state, unloaded reference, admissible quotient |
| \(X,\mathcal B_0,dV_0\) | HYPER-001; DYNAMICS-001 | authorized rank-three local reference realization |
| \(F=\operatorname{Grad}_0y\) | DEFORMATION-001; LOCALITY-001 | placement realization |
| \(C=F^\sharp F\), \(C(q_0)={\bf1}\) | DEFORMATION-001 | objective relative-deformation endomorphism |
| \(E,t,E_{\rm TF}\) | DEFORMATION-001; CONSTITUTIVE-CONSTRUCTION-001 | equivalent strain coordinates and isotropic split |
| \(\mathcal D_C\), hard endpoint | HYPER-001; ENERGY-PRINCIPLE-001; NONLINEARITY-001; ENERGY-SEARCH-001 | frozen finite-capacity admissible domain and authorized endpoint class |
| \(K_0,\mu_0>0\), \(Q\) | HYPER-001; ENERGY-PRINCIPLE-001; CONSTITUTIVE-PRINCIPLES-001; CONSTITUTIVE-CONSTRUCTION-001 | frozen isotropic reference tangent and selected minimal quadratic 2-jet |
| \(P_C=DW_A\) | HYPER-001; BALANCE-001; CONSTITUTIVE-CONSTRUCTION-001 | hyperelastic derivative/storage chain rule |
| \(P_F=2FP_C\) | LOCALITY-001; CONSTITUTIVE-CONSTRUCTION-001 | chain rule through \(C=F^\sharp F\) |
| \(-\operatorname{Div}_0P_F\) | LOCALITY-001; BALANCE-001 | adjoint of the local kinematic derivative and localization |
| weak internal-work integral | ENERGY-PRINCIPLE-001; EQUILIBRIUM-001; LOCALITY-001 | E1 first variation and integration by parts |
| \(N_{\mathcal K}\)/GE-009 | STATE-002; ENERGY-SEARCH-001; EQUILIBRIUM-001; CONSTITUTIVE-CONSTRUCTION-001 | admissible tangent cone and hard extended-value/subdifferential endpoint |
| \(\mathcal K_\tau\) slot | DYNAMICS-001; DURATION-001; BALANCE-001 | E2/clock-gauge kinetic closure explicitly left unselected |
| \(b,\bar t,\mathcal S_\tau\) slots | BALANCE-001; EQUILIBRIUM-001; LOCALITY-001 | source and external virtual-work classifications |
| \(\partial_\tau p+\operatorname{Div}J=\sigma\) | DURATION-001; BALANCE-001 | calibrated local balance template |
| Dirichlet/Neumann/mixed data | EQUILIBRIUM-001; LOCALITY-001 | admissible variations and the first-gradient boundary term |
| linear stress and Navier operator | DEFORMATION-001; HYPER-001; BALANCE-001; LOCALITY-001; CONSTITUTIVE-CONSTRUCTION-001 | reference expansion and frozen tangent |
| conditional longitudinal/shear waves | LOCALITY-001; CONSTITUTIVE-PRINCIPLES-001 | positive inertia plus acoustic tensor; kinetic normalization remains open |
| \(g^{\rm eff}=G[q,C;\mathcal D]\), cone and weak map | DURATION-001; METRIC-001 | operational clock/ruler and V11 compatibility family |
| variational/static vs history family | ENERGY-PRINCIPLE-001; DYNAMICS-001; EQUILIBRIUM-001 | E1 is fixed; E2 awaits action/kinetic selection |

Dependency graph:

```text
FOUNDATION-001 + STATE-002
  -> q, q0, admissible histories and gauge quotient
  -> DEFORMATION-001
       -> F = Grad0 y -> C = F^sharp F -> E, t, E_TF
       -> HYPER-001 + ENERGY-PRINCIPLE-001
            -> local objective W(C), PC = DW, stable reference tangent
            -> CONSTITUTIVE-002 / MATERIAL-DISCOVERY-001
            -> CONSTITUTIVE-PRINCIPLES-001 / SELECTION-001
            -> LOCALITY-001 + ENERGY-SEARCH-001 + NONLINEARITY-001
            -> EQUILIBRIUM-001
            -> CONSTITUTIVE-CONSTRUCTION-001
                 -> WA = Q + hard domain, K0, mu0
                 -> PC -> PF = 2 F PC -> -Div0 PF
                 -> static strong/weak equation + natural traction
                 -> endpoint normal cone / variational inequality

DYNAMICS-001 + DURATION-001 -----> K_tau closure slot ----\
BALANCE-001 ----------------------> source/balance slots ---+-> GE-001 family
METRIC-001 -----------------------> G and cone gate --------/

reference linearization
  -> frozen weak-field stress
  -> Navier operator
  -> waves only after positive kinetic closure
  -> V11 only after kinetic-cone and metric-map gates pass
```

Every displayed term therefore has a frozen origin or is visibly labeled as a
previously recognized non-constitutive closure slot.  No microscopic ontology,
new constitutive function, coefficient, phenomenological fit, weak-lensing
change, or V11 modification has entered the family.

## 12. Completion statement

The canonical interior elastic equation, weak and balance forms, all natural
spatial and admissible endpoint conditions, closure-dependent initial-value
structure, exact reference linearization, mathematical classification,
well-posedness requirements, remaining closure inventory, and term-level
dependency graph are complete.  The absence of a unique fully numerical
evolution PDE is a derived closure boundary, not an incompleteness of the
constitutive construction.
