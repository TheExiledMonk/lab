# PBUF DYNAMICS-001 — Derivation of the Native Evolution Principle

## 0. Decision and boundary

FOUNDATION-001, DEFORMATION-001, HYPER-001, ENERGY-PRINCIPLE-001, and
STATE-002 are fixed inputs. This milestone introduces no field, coupling,
constant, constitutive choice, or change to V11.

The native variational object is an oriented, unparameterized admissible curve
in the physical configuration space. If a representative parameter (s) is
used, its action density must be positively homogeneous of degree one in the
tangent. Thus the smallest action family is

\[
 \boxed{\mathfrak S[q]=\int_S \mathscr L(q(s),\dot q(s);q_0)\,ds,
 \qquad \mathscr L(q,a v;q_0)=a\mathscr L(q,v;q_0)\quad(a>0).}       \tag{D-001}
\]

Here spatial integration, when a local representative exists, is already part
of \(\mathscr L\). Equation (D-001) is an action-space definition, not a field
equation. The accepted inputs do not select a particular \(\mathscr L\), a
kinetic metric, a clock, or an inertia normalization.

There is also a sharp obstruction: the conventional expression
\(\int(T-W)ds\), with quadratic \(T\), is not invariant under arbitrary
monotone relabeling of emergent order. A native \(T-W\) action therefore cannot
be derived at this milestone. It becomes available only after an already
identified clock/duration structure is supplied and a parameter gauge is fixed;
adding a lapse, clock field, Jacobi energy, or new scale here is forbidden.

## 1. Native action space

Let \({\cal Q}_{\rm adm}\subset{\cal Q}_{\rm phys}\) be (S2-004), and let
\(S=[s_-,s_+]\) be an oriented interval. The minimum differentiable history
space needed by a first-order local action is

\[
 \widetilde{\cal H}^{,1}:=
 \{q:S\to{\cal Q}_{\rm adm}\mid q\text{ is piecewise }C^1\},
 \qquad
 {\cal H}^{,1}_{\rm phys}:=widetilde{\cal H}^{,1}/\operatorname{Diff}_+(S).
                                                                    \tag{D-002}
\]

This is the differentiable refinement of STATE-002's ordered continuous
history space. It does not promote \(s\) to time. Constant histories and
isolated zero tangents require a continuous extension of the integrand to the
zero section; a regular degree-one integrand is generally non-smooth there.

For each \(q\), physical velocities are tangent classes

\[
 v=[\dot{\widehat q}]\in T_q{\cal Q}_{\rm adm},                       \tag{D-003}
\]

where vertical gauge directions are quotiented out. Under \(s'=f(s)\),
\(v'=dq/ds'=v/f'(s)\). The ray \([v]_+=\{a v:a>0\}\), together with the
orientation of the history, is parameter independent. Velocity is therefore
an admissible mathematical tangent, not an observable rate until a clock is
derived.

### Local representative

If \(q\) has a local realization on the rank-three reference carrier
\({\cal B}_0\), locality permits

\[
 \mathscr L(q,v;q_0)=\int_{{\cal B}_0}
 \ell\!\left(j^r q(x),j^r v(x),C[q,q_0](x),W(C(x))\right)dV_0 .       \tag{D-004}
\]

Only arguments already obtained from \(q\), \(v\), and fixed \(q_0\) are
admissible. No independent metric, lapse, connection, matter variable, memory
variable, or higher-gradient energy is licensed. The minimum spatial order
\(r\) is realization-dependent: the accepted corpus fixes no component form of
\(q\) and hence no jet order for \(C[q,q_0]\). “Local” means dependence on a
finite jet at the same material point; ultralocal dependence is allowed.

The reference measure \(dV_0\) is the covariant material measure already used
by HYPER-001. If a concrete realization cannot supply it without extra
structure, (D-001) remains the intrinsic statement and that realization is not
action-ready.

### Admissibility conditions

An admissible \(\mathscr L:T{\cal Q}_{\rm adm}\to\mathbb R\) must:

1. descend to the STATE-002 gauge quotient and be scalar under reference
   relabelings;
2. be invariant under the accepted internal objectivity action, with tensor
   arguments transformed covariantly rather than componentwise fixed;
3. satisfy the positive degree-one condition in (D-001), exactly equivalent to
   orientation-preserving reparametrization invariance;
4. be local in the sense of (D-004), unless nonlocality is separately derived;
5. be continuous on its domain and differentiable in \(q,v\) along admissible
   nonzero tangents sufficiently to define a first variation;
6. have a continuous zero-section extension if zero-velocity segments are
   admitted; and
7. use \(W=\Phi(I_1,I_2,I_3)\) only as the accepted state scalar, without
   selecting \(\Phi\) or adding an energy term.

For a classical first variation, \(\mathscr L\) is Gateaux differentiable in
admissible directions and its derivatives are integrable. A strong classical
form would additionally need realization-specific spatial regularity and
integration-by-parts hypotheses; these are requirements for FIELD-001, not
deductions here. Twice differentiability is needed only for a second variation
or local well-posedness analysis.

## 2. Minimal kinetic structure

The minimum evolution object is the tangent bundle of the physical quotient,
\(T{\cal Q}_{\rm adm}\), restricted to admissible tangent cones at a hard
elastic boundary. No vector space coordinates for \(q\) are assumed.

The smallest kinetic family is the family \({\cal K}_1\) of objective,
gauge-basic, local, differentiable fiber functions

\[
 {\cal K}_1:=\{K:T{\cal Q}_{\rm adm}\to\mathbb R_{≥0}mid
 K(q,av)=aK(q,v),\ a>0\}.                                           \tag{D-005}
\]

Positivity or strict positivity away from gauge/zero directions is an
admissibility/stability condition, not a numerical choice. A Finsler norm is a
regular member of this family. A quadratic form
\(G_q(v,v)\) is not itself native because it is degree two; its square root is
degree one, but no \(G\) is selected by the accepted milestones. Degenerate
members are allowed when covariance produces constraint directions.

For any chosen admissible integrand, generalized momentum is the cotangent
functional

\[
 \boxed{p:=D_v\mathscr L(q,v)\in T_q^*{\cal Q}_{\rm adm},\qquad
 \langle p,\delta q\rangle:=D_v\mathscr L(q,v)[\delta q].}           \tag{D-006}
\]

This duality is the canonical pairing. It requires no metric and introduces no
new field: \(p\) is derived from \((q,v)\) after an action is chosen. Euler's
homogeneous-function identity gives the unavoidable structural relation

\[
 \langle p,v\rangle=\mathscr L(q,v),                                \tag{D-007}
\]

where differentiable. Consequently the canonical Hamiltonian
\(\langle p,v\rangle-\mathscr L\) vanishes identically on the image of the
Legendre map. This is a reparametrization identity, not an evolution equation.
The Legendre map cannot be fully regular in the radial velocity direction, so
a constrained variational formulation is mathematically unavoidable.

## 3. Energy decomposition audit

Let \(T(q,av)=a^kT(q,v)\). Under \(s'=f(s)\), invariance of
\(\int(T-W)ds\) for arbitrary positive \(f'\) requires both terms in the
integrand to be degree one in velocity. But \(W(q)\) is degree zero and a
standard quadratic kinetic energy has \(k=2\). Hence

\[
 \boxed{\text{native reparametrization invariance}\quad\not\Rightarrow\quad
 \mathscr L=T-W;\quad\text{indeed it excludes that form without duration
 structure}.}                                                       \tag{D-008}
\]

The smallest admissible native family is therefore

\[
 {\cal L}_{\rm nat}:=\{\mathscr L(q,v;C,W)\mid
 \mathscr L(q,av;C,W)=a\mathscr L(q,v;C,W),
 \text{ conditions 1--7 above hold}\}.                              \tag{D-009}
\]

This permits, but does not require, a symbolic split
\(\mathscr L=K+U\) only when every summand is separately degree one; a bare
\(-W\) is inadmissible. Multiplying \(W\) by a degree-one clock/rate functional
would restore covariance, but the accepted corpus supplies no such functional.
Likewise, a Jacobi square-root construction needs a kinetic metric and a fixed
energy level. Neither is authorized. These are closure alternatives for a
later milestone, not members derived here.

After a relational clock is independently derived and the gauge is fixed to
its duration \(\tau\), the gauge-fixed conditional family may take

\[
 S_\tau[q]=\int d\tau\,[T(q,dq/d\tau)-W(C[q,q_0])],                 \tag{D-010}
\]

with an objective local kinetic functional \(T\). Equation (D-010) records the
condition under which \(T-W\) becomes meaningful; it is not the native action
and does not select \(T\).

## 4. Variational framework

An admissible variation is a differentiable two-parameter family
\(q_\epsilon(s)\in{\cal Q}_{\rm adm}\), with
\(\eta=\partial_\epsilon q_\epsilon|_0\) a tangent section along \(q\). It must:

- respect the gauge quotient and any elastic-domain tangent cone;
- keep fixed endpoints, \(\eta(s_-)=\eta(s_+)=0\), or obey separately declared
  natural boundary conditions;
- preserve fixed boundary data on \(\partial{\cal B}_0\), unless corresponding
  boundary work is part of the authorized action; and
- not vary the fixed comparison prescription \(q_0\) unless a previously
  accepted reference family explicitly makes it dependent on the history.

Stationarity means only

\[
 \boxed{\delta\mathfrak S[q;\eta]
 :=\left.{d\over d\epsilon}\mathfrak S[q_\epsilon]\right|_{\epsilon=0}=0
 \quad\text{for every admissible }\eta.}                            \tag{D-011}
\]

Existence of (D-011) requires differentiability under the integrals and an
integrable dominating bound (or an equivalent functional-analytic theorem).
Deriving a local strong statement would additionally require a specified
realization, function spaces, boundary regularity, and legitimate integration
by parts. No Euler–Lagrange or field equation is written here.

Gauge-related variations lie in the null directions of the physical first
variation. Reparametrization invariance makes the tangential history variation
dependent, consistent with (D-007). At a hard boundary, stationarity is taken
over the admissible tangent cone and may be one-sided; barrier and asymptotic
endpoint classes instead restrict the domain or limiting behavior.

## 5. Symmetry audit

Noether conclusions are conditional on differentiability, an invariant
measure, admissible boundary behavior, and a symmetry acting on the chosen
action—not merely on \(W\).

| Accepted structure | Variational consequence | Conservation-law status |
|---|---|---|
| internal objectivity / rigid frame action | first variation annihilates symmetry directions; objective momentum map | angular-momentum-type law only if the realization makes this a continuous global action |
| material relabeling / spatial covariance | gauge identities among variations | Noether identities or constraint relations, not automatically a new charge |
| locality | localizable first variation and boundary flux structure once a realization is fixed | enables local balance form; does not itself create a conserved quantity |
| order relabeling | degree-one action, (D-007), degenerate Legendre map | parameter-Hamiltonian constraint; no physical energy conservation follows |
| shift of an independently derived clock | action independence of that clock | energy conservation only after such a clock and symmetry are established |
| homogeneity of the medium/reference | spatial translation symmetry | momentum-type conservation only if that symmetry is actually imposed on the full action |

Isotropy/objectivity alone does not imply spatial homogeneity. Emergent time
specifically prevents assuming ordinary time-translation symmetry. Local
Lorentz invariance in the V11 effective regime is a downstream matching gate;
until the one-metric map exists, no native Lorentz Noether current can be
claimed.

## 6. Compatibility and closure audit

| Input | Compatibility result |
|---|---|
| FOUNDATION-001 FP-1 | one action is assigned to histories of the one medium; no substrate is added |
| FP-2 | the action selects medium histories without postulating a fundamental gravitational interaction |
| FP-3 | histories are quotient by \(\operatorname{Diff}_+(S)\); no external clock or duration is smuggled in |
| FP-4 | each point on a history is one complete \(q\); history and tangent spaces remain representational |
| FP-5 / V11 | covariance and the effective local Lorentz/GR limits are mandatory matching tests, not derived here |
| FP-6 | no inertia scale, lapse, Jacobi energy, coefficient, or coupling is introduced |
| STATE-002 | \(q\in{\cal Q}_{\rm adm}\), quotient structure, reference state, and ordered histories are preserved |
| DEFORMATION-001 | only the accepted rank-three objective \(C[q,q_0]\) is used |
| HYPER-001 | \(W=\Phi(I_1,I_2,I_3)\) remains state-local, rate-free, and constitutively unselected |
| ENERGY-PRINCIPLE-001 | normalization, stability, domain, and endpoint restrictions are inherited unchanged |

The ingredients still required before FIELD-001 are: (i) a concrete realization
and function spaces for \(q\); (ii) a choice of \(\mathscr L\in{\cal L}_{\rm nat}\)
or a derived relational clock supporting (D-010); (iii) an inertial/kinetic
structure and normalization derived without a new free constant; (iv) boundary
and endpoint data; (v) the one-metric clock/ruler map and source projection;
(vi) matter action/coupling if matter is included; and (vii) proof of causal
well-posedness and recovery of V11's operational Lorentz and GR limits.

## 7. Dependency graph and traceability

\[
 \boxed{q\ \longrightarrow\ (q,[v]_+,C[q,q_0],W(C))
 \longrightarrow\ \mathfrak S[q]\ \longrightarrow\ \delta\mathfrak S[q].}
                                                                    \tag{D-012}
\]

| Equation | Content | Premises | Downstream use |
|---|---|---|---|
| D-001 | general native action and degree-one condition | FP-3; S2-012--013 | all variational dynamics |
| D-002 | physical action/history space | STATE-002 | admissible curves and variations |
| D-003 | generalized velocity | STATE-002 quotient | kinetic domain |
| D-004 | local representative | HYPER-001 measure/locality | realization gate |
| D-005 | kinetic family | D-001; objectivity | kinetic closure gate |
| D-006 | momentum and pairing | differentiable action | canonical data |
| D-007 | homogeneous identity | D-001, D-006 | constraint structure |
| D-008 | \(T-W\) no-go | FP-3; rate-free \(W\) | excludes hidden absolute time |
| D-009 | smallest native family | D-001--D-008 | FIELD-001 action choice |
| D-010 | conditional clock-gauge family | future clock principle | not yet native |
| D-011 | stationarity requirement | D-002, D-009 | FIELD-001 variation |
| D-012 | required dependency chain | all accepted inputs | milestone summary |

## 8. Completion statement

The native action space, evolution tangent, momentum pairing, admissible
degree-one action/kinetic families, variational conditions, and conditional
symmetry consequences are established. The result deliberately leaves the
action non-unique: selecting a kinetic geometry or physical clock is not a
mathematical consequence of the accepted inputs. The exact obstruction to a
native \(T-W\) decomposition is recorded, and every remaining FIELD-001 gate is
explicit. No field equation has been derived.
