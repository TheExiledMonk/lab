# PBUF METRIC-001 — Derivation of the Effective Relativistic Metric

## 0. Result and scope

The effective metric is not the physical medium and is not an additional
fundamental field.  It is the single Lorentzian, operational constitutive
response that packages the clock, ruler, propagation-cone, and synchronization
content of a medium state.  Its precise mathematical type is

\[
 \boxed{\;G:{\cal A}\longrightarrow
 \Gamma\!\left(\operatorname{Lor}(T^*M)\right),\qquad
 g^{\rm eff}=G[q,C[q,q_0];\mathcal D]\;}                         \tag{M-001}
\]

where \({\cal A}\) is the admissible set of accepted state/comparison data,
\(\operatorname{Lor}(T^*M)\subset\operatorname{Sym}^2T^*M\) is the open bundle
of nondegenerate tensors of signature \((-+++ )\), and \(\mathcal D\) denotes
the already accepted duration structure of DURATION-001.  The semicolon records
that \(C\) is derived from \(q\), not independent data.  Because duration is
itself a functional of \(q\) and physical processes, the requested shorthand
\(G[q,C]\) is valid; \(\mathcal D\) is shown only to expose the clock
identification that a rank-three \(C\) cannot supply.

The accepted inputs determine the admissible class \(\mathfrak G\), not a unique
member of it.  This is a derivation of the metric-map family and its exact
constraints.  Choosing a particular nonlinear response, derivative order, or
normalization beyond the V11 matching conditions would introduce unauthorized
constitutive information.

## 1. Mathematical role

For each state, \(g^{\rm eff}\) is simultaneously:

1. an **operational geometry**, because its quadratic form returns all ideal
   clock and ruler comparisons in the V11 regime;
2. a **measurement map**, because it converts effective event displacements
   into measured durations, spatial lengths, and causal classifications; and
3. a **constitutive response**, because the quadratic form depends on the
   physical medium state.

These are three aspects of one object, not three metrics.  It is emergent and
operationally physical (its intervals and cones are measurable), but it is not
ontologically fundamental and is not identical to \(q\), \(C\), or the medium.
Calling it merely a coordinate representation would be too weak: coordinate
components are conventional, whereas its interval and cone predictions are
invariant.

## 2. The admissible map

Let \(\rho_f\) be the induced action of a change of effective spacetime chart
\(f\) on the accepted state representation.  The admissible family is

\[
\begin{split}
\mathfrak G:=\{G\mid{}&G[q,C]\in\Gamma(\operatorname{Lor}(T^*M)),\\
 &G[\rho_fq,\rho_fC]=f^*G[q,C],\\
 &G[q,C]\text{ is gauge-basic and objective},\\
 &G[q_0,{\bf1}]_p\simeq\eta\text{ in every unloaded local V11 frame},\\
 &d\tau^2=-c^{-2}G[q,C]_{\mu\nu}dx^\mu dx^\nu
   \text{ for ideal timelike clocks},\\
 &G[q,C]_{\mu\nu}k^\mu k^\nu=0
   \text{ for universal V11 signal characteristics}\}.             \tag{M-002}
\end{split}
\]

Here \(\simeq\eta\) means equality after choosing a local inertial frame; it is
not equality of component arrays in every chart.  The constant \(c\) is the
already retained V11 conversion between clock and ruler units, not a new
coefficient.  Gauge-basic means that representatives of the same physical
\(q\) give the same effective tensor.  Objectivity means that a rigid
frame/material relabeling cannot create an observable deformation or metric
change.

### Dependence and regularity

The completely general map is a functional with kernel

\[
 \delta g^{\rm eff}_{\mu\nu}(x)=
 \int R_{\mu\nu A}(x,y)\,\delta q^A(y)\,d\mu(y),\qquad
 R:=D_q\bigl(G[q,C[q,q_0]]\bigr).                                  \tag{M-003}
\]

This includes the direct and \(C\)-mediated dependence by the chain rule,

\[
 R=D_1G+D_2G\circ D_qC .                                           \tag{M-004}
\]

No independent microscopic field has been added.  If locality is imposed in a
later constitutive closure, the admissible subclass is a finite-jet natural
operator

\[
 g^{\rm eff}(x)=\mathcal G\bigl(j_x^rq,j_x^rC[q,q_0];\mathcal D_x\bigr),       \tag{M-005}
\]

with finite but presently unselected order \(r\ge0\).  The accepted foundations
do not select ultralocal, finite-derivative, or causal nonlocal dependence.
Nonlocal maps must at least be covariant and causal in the resulting effective
cone; that condition is implicit until the cone is solved, so existence and
uniqueness are additional closure requirements.

Continuity is required on \({\cal Q}_{\rm adm}\) so nearby admissible states do
not produce discontinuous measurements.  A controlled weak limit requires
Fréchet (or appropriate Gateaux) differentiability at \(q_0\).  A quadratic
error estimate requires \(C^2\) there.  Global differentiability through a hard
elastic boundary is neither meaningful nor implied.  Every output is symmetric,
dimensionally compatible with the coordinate convention, nondegenerate, and
time-orientable on the operational domain.

The rank-three SPD tensor \(C\) supplies three principal ruler-deformation
channels.  It supplies no temporal eigenvalue, lapse, shift, clock calibration,
or synchronization rule.  Hence there is no function \(G(C)\) uniquely fixed by
the accepted data.  Completeness of \(q\) permits the missing operational content
to be a functional of \(q\), while DURATION-001 constrains that content through
clock measurements.  This proves why the admissible notation is \(G[q,C]\) and
why replacing it with a chosen conformal, disformal, coframe, or induced-metric
formula is not derivable here.

## 3. Operational meaning

For an effective displacement \(dx\):

\[
 \boxed{d\tau^2=-c^{-2}g^{\rm eff}(dx,dx)}                           \tag{M-006}
\]

is the duration accumulated by an ideal clock on a timelike path.  For an
observer represented by a future unit vector \(u\),

\[
 g^{\rm eff}(u,u)=-c^2,\qquad
 h^{(u)}_{\mu\nu}=g^{\rm eff}_{\mu\nu}+c^{-2}u_\mu u_\nu             \tag{M-007}
\]

is the positive spatial ruler metric on the observer's local rest space.
This \(u\) is an operational observer tangent, not a new microscopic medium
field.  Null directions satisfy

\[
 g^{\rm eff}_{\mu\nu}k^\mu k^\nu=0,                                \tag{M-008}
\]

and encode the universal signal cone.  Signals plus clock coincidences define
radar distance and a synchronization convention.  Simultaneity is therefore a
choice of spacelike hypersurfaces or local rest spaces, not an invariant extra
medium structure.  Clock proper duration and causal order are invariant;
coordinate time, one-way coordinate speed, and distant simultaneity are
chart/synchronization dependent.

Thus rulers measure propagation-defined spatial separation through the medium;
clocks count propagation-bearing cycles; and the metric is the unique effective
quadratic summary demanded by the accepted one-metric V11 limit.  It does not
claim that arbitrary microscopic processes are literally geodesic rods or
clocks outside that limit.

## 4. Weak-deformation recovery

Write \(q=q_0+\delta q\) in a local representative and
\(C={\bf1}+2\varepsilon+O(\varepsilon^2)\), as fixed by
DEFORMATION-001.  Differentiability gives

\[
 \boxed{g^{\rm eff}=\eta+h^{\rm eff}+O(\|\delta q\|^2),\qquad
 h^{\rm eff}=R_0[\delta q]
 =D_1G_0[\delta q]+2D_2G_0[\varepsilon].}                           \tag{M-009}
\]

For a nonlocal map, \(R_0[\delta q]\) has the integral form (M-003); for a
local map it is the corresponding local linear differential operator.  V11
compatibility is exactly the matching condition

\[
 [R_0\delta q]_{\rm gauge}=[h^{\rm V11}]_{\rm gauge},               \tag{M-010}
\]

together with one common null cone and (M-006) for all ideal matter clocks.
Consequently, in a local inertial chart, (M-006) reduces to
\(d\tau=dt\sqrt{1-|\mathbf v|^2/c^2}\), and (M-008) reduces to the V11
light cone.  Equations (M-009)--(M-010) demonstrate weak-deformation recovery
without selecting or fitting \(R_0\).  The foundations fix the target and
normalization, not which medium variations produce each allowed metric
perturbation.  No field equation follows from this kinematic match.

## 5. Near the finite elastic bound

Let \(C_n\) approach a boundary point of the accepted bounded SPD spectral
domain from its interior.  Admissibility requires, throughout the operational
domain,

\[
 \det G[q_n,C_n]\ne0,quad
 \operatorname{Inertia}(G[q_n,C_n])=(1,3),quad
 \text{a consistent time orientation}.                             \tag{M-011}
\]

There are only three structurally distinguishable possibilities consistent
with the accepted inputs:

1. **regular extension:** \(G[q_n,C_n]\) has a finite Lorentzian limit;
2. **operational boundary:** the map remains regular inside but the effective
   description ends at the elastic boundary; or
3. **loss of effective regularity before/on the boundary:** excluded from any
   domain on which V11 measurements are claimed valid, but not ruled out as a
   boundary of that domain.

A bounded \(C\) does not mathematically imply a bounded metric response or
bounded derivatives of \(G\).  Conversely, a hard or energetic elastic bound
does not force \(\det g\to0\), curvature blow-up, horizon formation, or any
named spacetime solution.  A globally regular V11 representation up to the
bound would additionally require uniform nondegeneracy, bounded \(G\) and
inverse \(G^{-1}\) in chosen local norms, and as many bounded derivatives as the
later equations require.  Those are closure conditions, not consequences of
the finite bound.

## 6. Covariance audit

Coordinate covariance belongs to the effective representation: (M-002) makes
\(G\) diffeomorphism-equivariant and all measured statements tensorial.  Local
Lorentz invariance is the frame symmetry of each output Lorentzian metric in
the V11 regime; it is a required operational symmetry, not evidence for a
fundamental Minkowski substrate.

The medium instead has the accepted objectivity/gauge invariance of \(q\), the
similarity covariance of \(C\), and reparametrization invariance of its ordered
histories.  These native symmetries need not be identified with spacetime
diffeomorphisms or Lorentz transformations until a realization of \(G\) supplies
that relation.  The map must descend to the quotient state and intertwine the
two representations; it may not turn a gauge/material relabeling into an
observable metric change.

DURATION-001 is respected because \(s\) never enters \(G\) as a clock and
(M-006) matches its invariant accumulated duration.  FOUNDATION-001 is respected
because the one medium remains fundamental, geometry and time remain emergent,
the complete state remains \(q\), no new constant appears, and the V11
Lorentzian operational limit is preserved.

## 7. Remaining closure gaps before FIELD-001

The following are not supplied by the metric map and remain mandatory:

1. a concrete local realization and function space for \(q\), including the
   event/material correspondence and boundary regularity;
2. selection of one \(G\in\mathfrak G\), including locality/kernel, derivative
   order, normalization away from \(q_0\), and proof of existence;
3. the clock/ruler universality theorem connecting all admissible physical
   processes to that one metric outside the infinitesimal regime;
4. selection of the accepted native action/kinetic functional and its measure;
5. the still-unselected stored-energy function within the accepted HYPER-001
   family and the implementation of the finite bound;
6. admissible variations, initial/boundary data, constraints, and gauge fixing;
7. the map from native constitutive response to effective stress-energy or
   geometric response, including matter coupling and source projection;
8. regularity, hyperbolicity/causality, stability, and well-posedness criteria;
9. the homogeneous reduction connecting local \(q,C,G\) to V11 background
   variables without double counting; and
10. a conservation/Noether compatibility proof and recovery tests for the full
    V11 weak-field and local-Lorentz regime.

These are prerequisites for governing equations; none is derived here.

## 8. Equation traceability and completion

| ID | Content | Accepted source / inference | Status |
|---|---|---|---|
| M-001 | typed medium-to-metric map | STATE-002 completeness; FP-2/FP-5; DURATION-001 | derived family |
| M-002 | definition of \(\mathfrak G\) | quotient objectivity; V11 one-metric limit; DU-009 | derived constraints |
| M-003 | general first variation/kernel | differentiability of M-001 | mathematical identity |
| M-004 | direct plus deformation chain rule | \(C=C[q,q_0]\) from S2-003/S2-008 | mathematical identity |
| M-005 | finite-jet local subclass | DYNAMICS-001 locality classification | conditional subclass |
| M-006 | ideal-clock duration | DURATION-001 DU-009 | accepted matching condition |
| M-007 | observer ruler metric | M-001/M-006; Lorentzian linear algebra | derived operational relation |
| M-008 | propagation cone | accepted V11 universal propagation | matching condition |
| M-009 | weak expansion | DEFORMATION-001 small strain; differentiable \(G\) | derived |
| M-010 | V11 perturbative recovery | FP-5; V11 fixed target | mandatory gate |
| M-011 | bound-domain signature constraints | HYPER-001/ENERGY-PRINCIPLE-001 finite bound | derived admissibility |

METRIC-001 is complete at the maximum strength authorized by the fixed inputs:
the map's type, admissible family, operational content, weak-limit matching,
finite-bound constraints, and covariance split are derived.  A unique formula
is underdetermined.  Presenting one would violate the prohibition on new
constitutive structure rather than complete the derivation.
