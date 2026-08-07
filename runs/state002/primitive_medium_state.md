# PBUF STATE-002 — Definition of the Primitive Medium State

## Canonical definition

FOUNDATION-001, DEFORMATION-001, HYPER-001, and ENERGY-PRINCIPLE-001 are held
fixed. The definition of \(q\) is an accepted premise.

Let \({\cal M}\) denote the one medium, \({\cal A}\) its admissible complete
configuration representatives, and let \({\cal G}\) contain only the descriptive
redundancies already accepted by DEFORMATION-001:

\[
{\cal G}=\operatorname{Diff}({\cal M})\ltimes{\cal G}_{\rm int}.             \tag{S2-001}
\]

The canonical mathematical type is

\[
\boxed{{\cal Q}_{\rm phys}:={\cal A}/{\cal G},\qquad
q=[\widehat q]_{\cal G}\in{\cal Q}_{\rm phys}.}                              \tag{S2-002}
\]

Thus \(q\) is a point of a quotient configuration space—globally, a gauge class
of admissible sections—not a coordinate tuple, metric, clock, embedding, or
representation of spacetime. The representative \(\widehat q\) is not an extra
field. Exactly one \(q\) is occupied at each order state; the state space is a
representational possibility set (FP-4).

Completeness means every later physical object is a functional of \(q\), or a
relational functional of \(q\) and a declared comparison state:

\[
{\cal O}={\cal O}[q],\qquad C=C[q,q_0].                                      \tag{S2-003}
\]

It does not claim that a component realization of \(q\) has been discovered.

## State space, regularity, and objectivity

The smallest admissible state space justified by the accepted chain is

\[
{\cal Q}_{\rm adm}:=\{q\in{\cal Q}_{\rm phys}:C[q,q_0]\in{\cal D}_C\},
\quad{\cal D}_C\subset\operatorname{Sym}^+(3),\quad{\bf1}\in{\cal D}_C.       \tag{S2-004}
\]

\({\cal D}_C\) is the path-connected, objective, permutation-invariant spectral
domain fixed by HYPER-001 and ENERGY-PRINCIPLE-001. If their accepted finite-bound
premise is active, its closure is compactly contained in the positive cone.
STATE-002 selects neither its boundary nor endpoint implementation.

The minimum topology is the initial topology induced locally by the comparison:

\[
q_n\to q\ \Longrightarrow\ C[q_n,q_0]\to C[q,q_0].                           \tag{S2-005}
\]

The map \(q\mapsto C\) must be continuous. It need only be \(C^1\) where a
classical energetic response is pulled back to state space, and \(C^2\) near
\(q_0\) where a pulled-back weak tangent is required. No stronger spatial
smoothness follows because the accepted energy contains no \(\nabla C\).
Specific Sobolev class, boundary regularity, and effective-metric
differentiability await a realization and dynamics.

Objectivity requires the comparison to descend through the quotient:

\[
C[g\widehat q,g\widehat q_0]
=R_g^{-1}C[\widehat q,\widehat q_0]R_g,\qquad g\in{\cal G}.                  \tag{S2-006}
\]

Consequently the unordered spectrum and \(I_1,I_2,I_3\) are
representative-independent. A superposed rigid spatial rotation cancels. This is
objectivity, not invariance of tensor components.

## Reference state

The canonical undeformed equilibrium state is

\[
\boxed{q_0=[\widehat q_0]_{\cal G}\in{\cal Q}_{\rm adm},\qquad
C[q_0,q_0]={\bf1},\quad W({\bf1})=0,\quad DW({\bf1})=0.}                     \tag{S2-007}
\]

It is nondegenerate and orientation/signature admissible, admits the local V11
Minkowski limit, and is homogeneous and isotropic when used as the V11
cosmological background. The energy equality fixes only the additive zero.

Representatives \(g\widehat q_0\) are one physical reference. Once an undeformed
reference prescription is fixed, \(q_0\) is unique modulo gauge. The accepted
inputs do not prove that \(W\) has no other minimizer and do not choose among the
fixed, instantaneous-natural, or parameter-dependent reference families left
open by DEFORMATION-001. Any permitted choice obeys (S2-007); it is not a second
primitive state variable.

## Relationship to deformation

For an admissible pair let \(F(q,q_0):V_0\to V_q\) be a representative relative
material/coframe isomorphism. With the unloaded structure defining the adjoint,

\[
\boxed{C(q,q_0):=F(q,q_0)^{\sharp_0}F(q,q_0):V_0\to V_0.}                   \tag{S2-008}
\]

This is exactly the abstract mapping fixed by DEFORMATION-001, and

\[
C=C^{\sharp_0}>0,\qquad\det C>0,\qquad C(q_0,q_0)={\bf1}.                   \tag{S2-009}
\]

Under \(F\mapsto QF\), \(Q^\sharp Q={\bf1}\), \(C\) is unchanged. Under a
reference relabelling \(R\),

\[
C\mapsto R^{-1}CR,\qquad
\operatorname{spec}C,\ I_1,\ I_2,\ I_3\ \hbox{unchanged}.                   \tag{S2-010}
\]

This includes \(C^I{}_J=B^{IK}\kappa_{KJ}\). An equivalent relative
symmetric-tensor realization does not make \(q\) a metric. The abstract mapping
is unique by authority; \(F\) is nonunique up to left orthogonal action.
\(E=(C-{\bf1})/2\) and \(H=(\log C)/2\) are reparametrizations, not new states.
The map need not be injective because \(q\) is complete while \(C\) records only
relative deformation. No extra distinguishing sector is introduced here.

On a relative tensor representative, the accepted weak limit is

\[
C={\bf1}+2\varepsilon+O(\varepsilon^2),\qquad
\delta C=q_0^{-1}\delta q+O(\delta q^2).                                    \tag{S2-011}
\]

This preserves DEFORMATION-001 and does not define a metric map.

## Evolution space

Let \((S,\prec)\) be an oriented totally ordered set, represented conveniently by
an interval. An admissible history and its equivalence space are

\[
\gamma:S\to{\cal Q}_{\rm adm},\quad s\mapsto q(s),                           \tag{S2-012}
\]
\[
\boxed{{\cal H}_{\rm adm}:=
C^0_{\rm ord}(S,{\cal Q}_{\rm adm})/\operatorname{Homeo}_+(S).}             \tag{S2-013}
\]

Strictly increasing relabellings preserve the history. Order and the curve in
state space are meaningful; the numerical scale, origin, and rate of \(s\) are
not. No fundamental time dimension, clock, lapse, or duration is introduced.

Minimum evolution is continuous, so \(C[q(s),q_0]\), and \(W\) wherever finite
and continuous, vary continuously. Piecewise \(C^1\) is optional when a tangent
is useful; no speed or second-order law follows. Jumps require a separately
authorized transition rule and are not in the minimal evolution space.

Each \(s\) labels one realized state. This does not imply deterministic
continuation:

\[
q(s_1)\ \not\Rightarrow\ \hbox{a unique }q(s_2),\qquad s_1\prec s_2.         \tag{S2-014}
\]

A segment is dynamically reversible only if an order-reversing bijection
\(\rho\), together with any required state involution \(\Theta\), gives

\[
q_{\rm rev}(s)=\Theta q(\rho(s))\in{\cal Q}_{\rm adm}                       \tag{S2-015}
\]

and the future evolution law admits that reversed history. The single-valued,
history-free \(W\) gives quasistatic path recovery on one elastic branch, not
dynamical reversibility, invertibility, or an arrow of order.

## Minimum degrees of freedom

At a material point the accepted rank-three tensor representative occupies

\[
\operatorname{Sym}^+(3)\simeq GL^+(3)/SO(3),\qquad\dim=9-3=6.                \tag{S2-016}
\]

These are six representative components. DEFORMATION-001 further fixes the
physical deformation data as the spectrum: material relabellings act by
similarity, so principal-axis components are descriptive and the quotient
\(\operatorname{Sym}^+(3)/SO(3)\) has three generic coordinates. They may be the
three positive eigenvalues or \(I_1,I_2,I_3\). This does not make \(C\) a scalar
or vector: the objective tensor is required to encode generic volume and shear,
even though its gauge-invariant local orbit is three-dimensional.

Because every admissible \(C\) must derive from \(q\), the comparison map must be
locally surjective onto the admissible spectral quotient. Therefore

\[
\boxed{\dim_{\rm local}{\cal Q}_{\rm phys}\ge3,\qquad
\dim_{\rm local}^{\rm minimal}=3\ \hbox{spectral degrees of freedom}.}       \tag{S2-017}
\]

Globally this is a field configuration—minimally three independent spectral
functions over the continuous medium after the local similarity quotient—so the
global space is infinite-dimensional. A tensor representative uses six component
functions before that quotient. Completeness may eventually require information
not visible in \(C\), but additional microscopic sectors are not authorized.
Three is therefore the proven physical lower bound and canonical minimum; it is
not a claim that \(q\) has been replaced by its three invariants.

## Dependency graph and compatibility

    FOUNDATION-001 + accepted STATE-002 definition
                         |
                         v
              q=[q-hat]_G in Q_adm
                         |
                 compare with q0
                         v
        C(q,q0)=F^sharp_0 F in D_C
                         |
                         v
          W(C)=Phi(I1,I2,I3) in F_EP
                         |
                         v
                    P_C=DW(C)

Thus

\[
\boxed{q\longrightarrow C[q,q_0]\longrightarrow W(C)
\longrightarrow P_C:=DW(C).}                                               \tag{S2-018}
\]

| Authority | Compatibility result |
|---|---|
| FP-1 | One quotient state belongs to the one medium; no substrate or second sector is added. |
| FP-2 | Deformation is available for a later gravity map; no field equation is derived. |
| FP-3 | (S2-012)--(S2-015) use order modulo monotone relabelling, not a time dimension. |
| FP-4 | Exactly one \(q(s)\) is occupied; state/history spaces remain representational. |
| FP-5 / V11 | \(q_0\) admits the local Minkowski limit and (S2-011) preserves symmetric weak-field kinematics. No preferred observable foliation or saturation-to-strain map is added. |
| FP-6 | No coefficient, coupling, or constant is introduced. |
| DEFORMATION-001 | (S2-008)--(S2-011) preserve the objective dimensionless SPD rank-three \(C\). |
| HYPER-001 | \(C\) feeds exactly \(W=\Phi(I_1,I_2,I_3)\) on the accepted domain. |
| ENERGY-PRINCIPLE-001 | No \(\Phi\), modulus, endpoint, or unique minimum is selected. |

V11 compatibility remains a downstream gate. The effective Lorentzian metric,
clocks/rulers, source map, homogeneous reduction, and dynamics are not derived.
The ordering label cannot define an observable preferred clock or foliation.

## Completion boundary

(S2-002), (S2-007), (S2-008), (S2-013), and (S2-017) are the canonical state,
reference, deformation map, evolution space, and minimum-DOF results. No field
equation, constitutive selection, photon rule, microscopic field, coupling,
fundamental constant, or V11 modification has been introduced.
