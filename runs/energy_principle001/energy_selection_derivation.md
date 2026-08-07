# PBUF ENERGY-PRINCIPLE-001 — Derivation of the Native Energy Selection Principle

## 0. Decision and scope

FOUNDATION-001, DEFORMATION-001, HYPER-001, and the V11 operational regime are
held fixed.  This milestone derives restrictions; it does not select a
constitutive formula.

The strongest justified result is

\[
 \boxed{W(C)=\Phi(I_1,I_2,I_3),\qquad C\in {\cal D}_C\subset
 \operatorname{Sym}^+(3)}                                                   \tag{EP-001}
\]

where \(I_1=\operatorname{tr}C\),
\(I_2=((\operatorname{tr}C)^2-\operatorname{tr}C^2)/2\), and
\(I_3=\det C\).  The domain is objective and permutation invariant, contains
\({\bf1}\), and—under the separately accepted finite-bound premise—has compact
closure inside the positive spectral cone.  Near \(i_0=(3,3,1)\), \(\Phi\) is
at least \(C^2\), is normalized and stationary there, and has a positive
definite quadratic form if strict local elastic stability is required.

These restrictions do **not** select a unique \(\Phi\).  A claim of uniqueness
from the accepted inputs would be an unauthorized constitutive assumption.

## 1. Logical-status rule

The source of each statement matters:

* **P** — postulate: an FP statement, not a theorem;
* **D** — mathematical consequence of cited accepted inputs;
* **A** — additional accepted premise already present in the mission/HYPER-001;
* **C** — compatibility gate, necessary only for a proposed completion;
* **S** — sufficient option, not mandatory.

In particular, a finite elastic bound does not occur in FP-1--FP-6.  It is a
fixed mission/HYPER-001 premise.  Likewise, continuity of the medium does not by
itself prove locality, homogeneity, isotropy, hyperelasticity, convexity, or a
unique ground state.

## 2. Restrictions attributable to the accepted foundations

### 2.1 FP-1: one continuous physical medium

FP-1 prohibits a sum of independently ontological substrate energies.  Given
the already accepted local isotropic family, there is one scalar energy density
for the medium, not one freely chosen law per representation:

\[
 W(R^\sharp C R)=W(C).                                                       \tag{EP-002}
\]

Equation (EP-002), however, follows mathematically from the accepted
objectivity/isotropy of HYPER-001, not from the word “one.”  FP-1 alone imposes
no derivative sign, convexity, endpoint behavior, or energy zero.

### 2.2 FP-2: emergent gravity

The energy must possess enough differentiability for an elastic response to be
defined:

\[
 P_C=DW(C)=\Phi_1{\bf1}+\Phi_2(I_1{\bf1}-C)+\Phi_3I_3C^{-1}.                 \tag{EP-003}
\]

Thus \(C^1\) interior regularity is required wherever a classical response is
used.  FP-2 does not determine the sign, magnitude, or functional form of that
response, nor the map from \(P_C\) to gravitational sources.  Calling
(EP-003) gravity before that map is supplied would exceed FP-2.

### 2.3 FP-3: emergent time

The accepted \(W\) is statewise: it depends on the current relative
configuration and contains no fundamental temporal eigenvalue, rate, or history
variable.  Along any monotone relabeling \(s\mapsto f(s)\) of the emergent order,

\[
 W[C(s)]=W[C(f^{-1}(f(s)))].                                                  \tag{EP-004}
\]

This is parameter independence of a state function, not an evolution law.
Rate independence and reversibility follow from the accepted hyperelastic
functional *provided the state remains on one branch*.  FP-3 alone does not
exclude future dissipation; the present functional class does.

### 2.4 FP-4: one complete physical configuration per state

At an order state \(s\), only \(C(s)\) is occupied.  This neither requires a
unique minimizer of \(W\) nor forbids multiple mathematical equilibria.  The
configuration space and reference configuration remain representational.
Consequently FP-4 produces no convexity or uniqueness theorem.

### 2.5 FP-5: V11 compatibility

FP-5 is a gate.  At the unloaded state the local effective description must
admit the V11 Minkowski/Lorentz regime, and its weak perturbations must be able
to reproduce the retained relativistic response.  Necessary constitutive
regularity is

\[
 \Phi\in C^2(U),\quad i_0\in U\subset{\cal D}_I,                              \tag{EP-005}
\]

so that a quadratic weak-field tangent exists.  The tangent has the inevitable
isotropic form

\[
 W={K\over2}(\operatorname{tr}\varepsilon)^2
   +\mu\,\varepsilon_{\rm TF}:\varepsilon_{\rm TF}
   +O(|\varepsilon|^3).                                                       \tag{EP-006}
\]

Compatibility with V11 requires that, after the missing metric/clock and source
maps are supplied, (EP-006) reproduce rather than alter V11's operational local
Lorentz and weak-field GR behavior.  No equation involving only \(\Phi\) can
state that matching condition yet.  V11's \(R_{\max}\), \(\epsilon_{0,T}\),
\(\alpha_T\), \(S(a)\), and \(\Omega_\sigma(a)\) do not define \({\cal D}_C\),
a principal stretch, or a derivative of \(\Phi\).

### 2.6 FP-6: no unauthorized constants

No coefficient, endpoint, exponent, modulus, or energy scale may be inserted as
an independently adjustable fundamental constant.  Symbols such as

\[
 \mu=-2(\Phi_2^0+\Phi_3^0),\qquad
 K=4\left[\Phi_2^0+\Phi_3^0+v^T H^\Phi v\right]+{2\over3}\mu,
 \quad v=(1,2,1)^T,                                                         \tag{EP-007}
\]

are derivative combinations, not newly derived numerical constants.  Their
values and the overall energy scale remain unclosed.  (Equivalently use the
HYPER-001 definitions \(\lambda=4[\Phi_2^0+\Phi_3^0+v^TH^\Phi v]\) and
\(K=\lambda+2\mu/3\).)

## 3. Mandatory energy-selection constraints

The unloaded reference is fixed by DEFORMATION-001 as \(C={\bf1}\).  Choosing
its energy zero is a harmless normalization:

\[
 \Phi(3,3,1)=0.                                                              \tag{EP-008}
\]

For an unloaded, stress-free reference,

\[
 DW({\bf1})=0
 \iff \Phi_1^0+2\Phi_2^0+\Phi_3^0=0.                                        \tag{EP-009}
\]

If “stored energy” includes stable unloaded equilibrium, nonnegativity makes it
a global minimum:

\[
 W(C)\ge0,qquad W(C)=0\text{ at }C={\bf1}.                                  \tag{EP-010}
\]

Neither uniqueness of that minimum nor global convexity follows.  Strict local
stability is the additional stability requirement

\[
 D^2W({\bf1})[H,H]>0\quad(H\ne0),
 \qquad \Longleftrightarrow\qquad \mu>0, K>0.                               \tag{EP-011}
\]

Global strict convexity on a convex domain would ensure at most one minimizer,
but is only sufficient.  Finite-deformation convexity in the native variable is
equivalent, where eigenvalues are distinct, to

\[
 [w_{,AB}]\succeq0,qquad
 {w_{,A}-w_{,B}\over\lambda_A-\lambda_B}\ge0,                               \tag{EP-012}
\]

with continuous divided-difference limits.  Objectivity does not imply
(EP-012).  Rank-one convexity, polyconvexity, and hyperbolicity require the
missing realization/kinetic completion and are not restrictions derived here.

Locality, frame indifference, isotropy, parity evenness, hyperelasticity, and
rate independence are mandatory because HYPER-001 is authoritative.  They are
not newly derived from FP-1--FP-6.  Reversibility is conditional on remaining
within this single-valued, history-free elastic branch.  Dissipation, fracture,
plasticity, and hysteresis would require unauthorized state variables.

## 4. Finite-bound audit

Let \(\Lambda=(\lambda_1,\lambda_2,\lambda_3)\).  The admitted set must be
path-connected (for elastic access from the reference), permutation invariant,
contain \((1,1,1)\), and satisfy

\[
 \overline{{\cal D}_\lambda}\Subset(0,\infty)^3,
 \quad\text{i.e.}\quad
 0<\underline\lambda\le\lambda_A\le\overline\lambda<\infty.               \tag{EP-013}
\]

The bounds in (EP-013) are existential, not material coefficients or a box
domain.  The invariant domain is the exact image

\[
 {\cal D}_I=\{(e_1(\Lambda),e_2(\Lambda),e_3(\Lambda)):
                   \Lambda\in{\cal D}_\lambda\},                            \tag{EP-014}
\]

and hence is bounded but is not an arbitrary subset of \(\mathbb R^3\); its
points must be roots of a cubic with three positive real roots.

There are only three endpoint implementations authorized by HYPER-001:

1. a hard constraint, represented by the lower-semicontinuous extension
   \(W=+\infty\) outside the closed admissible set;
2. an interior barrier, \(W(C_n)\to+\infty\) at each forbidden boundary;
3. a finite regular one-sided endpoint plus a separate state constraint.

A bound inaccessible under arbitrary finite energetic loading needs option 1
or 2.  A finite endpoint alone does not prevent continuation.  No saturation
shape, limiting spectrum, singular exponent, stress plateau, bounded energy,
or vanishing tangent follows.

## 5. Smallest remaining family

Let \({\cal F}_{EP}\) be all pairs \((\Phi,{\cal D}_I)\) satisfying:

\[
\begin{aligned}
{\cal F}_{EP}=\{(\Phi,{\cal D}_I):\;&{\cal D}_I=\iota({\cal D}_\lambda),\
 (1,1,1)\in{\cal D}_\lambda,\
 \overline{{\cal D}_\lambda}\Subset(0,\infty)^3,\\
&{\cal D}_\lambda\text{ path-connected and permutation invariant};\\
&\Phi\in C^1(\operatorname{int}{\cal D}_I)\cap C^2(U_{i_0}),\
 \Phi(i_0)=0,\\
&\Phi_1^0+2\Phi_2^0+\Phi_3^0=0,\
 \Phi\circ\iota\ge0,\ \mu>0,\ K>0;\\
&\text{one of the three endpoint classes is specified}\}.                  \tag{EP-015}
\end{aligned}
\]

Here \(\iota(\Lambda)=(e_1,e_2,e_3)\).  The nonnegativity condition is inherited
from HYPER-001's stable stored-energy requirement.  If only neutral rather than
strict stability is demanded, replace the strict tangent inequalities by
non-strict ones; FP-1--FP-6 alone do not choose between those stability grades.

The remaining freedoms are infinite-dimensional: the shape of \({\cal D}_I\),
the endpoint class, every higher-order derivative of \(\Phi\), the two tangent
derivative combinations and overall scale, possible additional minima, and
finite-deformation convexity.  FP-6 prevents treating any of them as fitted new
fundamental constants, but does not calculate them.

## 6. Constraint matrix

| Constraint | Domain | Symmetry | Regularity | Convexity/stability | Endpoint | Normalization/class | Status |
|---|---:|---:|---:|---:|---:|---:|---|
| one SPD rank-3 \(C\) | yes | — | — | — | — | one scalar density | D: DEFORMATION/HYPER |
| invariant dependence | — | yes | — | — | — | \(\Phi(I_1,I_2,I_3)\) | D: HYPER |
| state-local/rate-free | — | — | — | — | — | excludes history/gradients | A: HYPER |
| classical response | — | — | \(C^1\) interior | — | — | — | C: FP-2 |
| weak tangent | — | — | \(C^2\) near \(i_0\) | two isotropic channels | — | — | C: FP-5/V11 |
| finite bound | compact in SPD cone | permutation invariant | interior only | none implied | 3 classes | extended-valued allowed | A + D |
| unloaded reference | contains \(i_0\) | — | stationary | minimum if stable | — | zero additive constant | A/D |
| strict local stability | — | — | Hessian exists | \(\mu>0,K>0\) | — | — | C, not FP theorem |
| unique minimizer | — | — | — | strict convexity sufficient | — | — | not derived |
| no new constants | endpoints symbolic | — | derivatives symbolic | moduli symbolic | rates absent | forbids fitted choices | P: FP-6 |
| V11 saturation history | no local-domain inference | — | — | no derivative inference | no spectrum inference | downstream mapping needed | C: FP-5 |

## 7. FIELD-001 closure gates

No governing equation can yet be derived because all of the following are
missing:

1. the primitive realization \(C[q,q_0]\) chosen from the accepted kinematic
   family and the status/evolution of \(q_0\);
2. the one-metric and clock/ruler identification mapping medium variables to
   the effective Lorentzian metric used by V11;
3. an action, integration measure, kinetic/inertial term, and emergent-order
   variational prescription;
4. the map from \(P_C\) to spacetime stress-energy/geometric response;
5. matter variables, their coupling to the one metric, and source projection;
6. admissible variations, boundary/initial data, and treatment of the elastic
   boundary;
7. constitutive closure selecting \(\Phi\) and its scale, or a proof that the
   field equations depend only on already fixed derivative combinations;
8. a homogeneous-reduction map connecting \(C\) to V11's \(a\),
   \(\Omega_\sigma(a)\), and thermal lookup variables;
9. proof that the resulting local equations retain the V11 Lorentz/Einstein
   operational regime and possess a well-posed evolution.

These are dependency statements, not field equations.

## 8. Conclusion

The native **selection principle presently available is a restriction
principle, not a uniqueness principle**.  It fixes invariant dependence,
reference normalization, response regularity, the stable quadratic cone, and
finite-domain alternatives.  It leaves an infinite-dimensional admissible
family.  Selecting one member requires additional independently derived
physics; the accepted foundations do not contain it.
