# PBUF HYPER-001 — Derivation of the native stored-energy functional

## 0. Status and scope

**Result.** On the accepted rank-three material branch, the complete minimal local,
parity-even, isotropic hyperelastic family is

\[
 \boxed{\;W(C)=\Phi(I_1(C),I_2(C),I_3(C)),\qquad C\in\mathcal D_C\subset
 \operatorname{Sym}^+(3)\;}                                      \tag{H-001}
\]

with the extended-value completion (W=+\infty) outside the closed admissible
set when the elastic bound is imposed as a hard state constraint.  Here

\[
 I_1=\operatorname{tr}C,\qquad
 I_2={1\over2}\big[(\operatorname{tr}C)^2-\operatorname{tr}(C^2)\big],
 \qquad I_3=\det C.                                                \tag{H-002}
\]

The scalar function \(\Phi\), its admissible invariant domain, and its derivatives
are constitutive data.  FOUNDATION-001, DEFORMATION-001, and V11 do **not** select
a unique \(\Phi\), a limiting spectrum, or numerical moduli.  Thus (H-001), with
the restrictions derived below, is the strongest derivation authorized by the
inputs.  Named material laws, a saturation curve, field equations, and an
identification of cosmological variables with local strain are outside this
milestone.

The ontology premises used here are explicitly postulates, not results:
FOUNDATION-001 FP-1--FP-6.  Locality and isotropy are authorized by the mission.
The primitive rank-three realization and \(C\) are authorized by DEFORMATION-001.

## 1. Representation theorem and minimality

At a material point, objectivity has already removed a superposed spatial rotation:
if \(F\mapsto QF\), \(Q^\sharp Q=\mathbf1\), then

\[
 C=F^\sharp F\mapsto F^\sharp Q^\sharp QF=C.                       \tag{H-003}
\]

Material isotropy acts by orthogonal similarity, \(C\mapsto R^\sharp C R\).
Consequently an objective isotropic scalar obeys

\[
 W(R^\sharp C R)=W(C).                                             \tag{H-004}
\]

Every \(C\in\operatorname{Sym}^+(3)\) has three positive eigenvalues
\(\lambda_1,\lambda_2,\lambda_3\).  Equation (H-004) makes \(W\) a symmetric
function \(w(\lambda_1,\lambda_2,\lambda_3)\).  The elementary symmetric
polynomials separate unordered triples, so the isotropic representation theorem
gives (H-001)--(H-002).  Three generic scalar arguments are minimal: fewer cannot
separate all unordered spectra and would discard a volumetric or shear channel.
The spectrum itself, or \(J=\sqrt{I_3}\) plus two invariants of
\(\bar C=I_3^{-1/3}C\), is an equivalent coordinate system, not extra physics.

No \(\nabla C\), curvature, history variable, four-velocity, parity-odd scalar,
or rank-four clock/ruler component occurs.  Therefore the energy is state-local,
hyperelastic, and rate-independent.  A reference-volume energy is integrated as
\(\int_{\mathcal B_0}W(C)\,dV_0\); this is invariant under coordinate changes
because \(W\) is a material scalar and \(dV_0\) is the reference volume element.

### Regularity

For a well-defined stress, require \(\Phi\in C^1(\mathcal D_I)\).  For tangent
stiffness and the weak-deformation expansion, require \(\Phi\in C^2\) in a
neighborhood of \(i_0=(3,3,1)\).  A continuous tangent stiffness uses \(C^2\);
standard local existence analyses commonly require \(C^{2,1}\) or \(C^3\), but
that stronger condition is not derivable here.  At a hard boundary the proper
global requirement is lower semicontinuity of the extended-valued energy; smoothness
is required only in the interior.

Normalize the unloaded state without adding a parameter:

\[
 W(\mathbf1)=\Phi(3,3,1)=0,\qquad DW(\mathbf1)=0.                   \tag{H-005}
\]

Adding a constant does not change response, so the first equality fixes the energy
zero.  Nonnegativity and (H-005) make the unloaded state a global minimum.

## 2. Constitutive response

The invariant derivatives are

\[
 {\partial I_1\over\partial C}=\mathbf1,\quad
 {\partial I_2\over\partial C}=I_1\mathbf1-C,\quad
 {\partial I_3\over\partial C}=I_3C^{-1}.                          \tag{H-006}
\]

Writing \(\Phi_a=\partial\Phi/\partial I_a\), the energetic tensor conjugate to
\(C\) is therefore

\[
 \boxed{\;P_C:={\partial W\over\partial C}
 =\Phi_1\mathbf1+\Phi_2(I_1\mathbf1-C)+\Phi_3I_3C^{-1}.\;}          \tag{H-007}
\]

It is symmetric, coaxial with \(C\), and transforms covariantly under material
frame changes.  In principal axes,

\[
 p_A={\partial w\over\partial\lambda_A}
 =\Phi_1+\Phi_2(\lambda_B+\lambda_C)+\Phi_3\lambda_B\lambda_C,
 \quad\{A,B,C\}=\{1,2,3\}.                                       \tag{H-008}
\]

At the reference state, stress freedom is the single derivative constraint

\[
 \Phi_1^0+2\Phi_2^0+\Phi_3^0=0.                                   \tag{H-009}
\]

If a later realization uses \(C=F^\sharp F\), the conventional second
Piola-type response is \(2P_C\).  That factor is a choice of conjugate pair, not
a new coefficient.  No spacetime stress-energy or field equation is asserted here.

For completeness, the tangent map in invariant notation is

\[
 \mathbb A[H]=\sum_{a,b}\Phi_{ab}(DI_b[H])G_a
 +\Phi_2\big(\operatorname{tr}H\,\mathbf1-H\big)
 +\Phi_3I_3\big(\operatorname{tr}(C^{-1}H)C^{-1}-C^{-1}HC^{-1}\big), \tag{H-010}
\]

where \(G_1=\mathbf1\), \(G_2=I_1\mathbf1-C\),
\(G_3=I_3C^{-1}\), and \(DI_b[H]=G_b:H\).  This is the most general local
isotropic tangent response implied by (H-001).

## 3. Weak-deformation limit

Set \(C=\mathbf1+2\varepsilon\), \(t=\operatorname{tr}\varepsilon\), and
\(s_2=\operatorname{tr}(\varepsilon^2)\).  Through quadratic order,

\[
 I_1=3+2t,\qquad
 I_2=3+4t+2(t^2-s_2),\qquad
 I_3=1+2t+2(t^2-s_2)+O(|\varepsilon|^3).                            \tag{H-011}
\]

Let \(v=(1,2,1)^T\), \(g_a=\Phi_a(i_0)\), and
\(H^\Phi_{ab}=\Phi_{ab}(i_0)\).  Taylor expansion and (H-009) give

\[
 W=2\,[g_2+g_3+v^TH^\Phi v]t^2-2(g_2+g_3)s_2
 +O(|\varepsilon|^3).                                              \tag{H-012}
\]

Define, rather than empirically introduce, the two tangent moduli

\[
 \mu:=-2(g_2+g_3),\qquad
 \lambda:=4[g_2+g_3+v^TH^\Phi v].                                 \tag{H-013}
\]

Then the universal isotropic tangent form is

\[
 \boxed{\;W={\lambda\over2}(\operatorname{tr}\varepsilon)^2
 +\mu\operatorname{tr}(\varepsilon^2)+O(|\varepsilon|^3).\;}       \tag{H-014}
\]

Equivalently, with \(\varepsilon=\tfrac13t\mathbf1+arepsilon_{\rm TF}\),

\[
 W={K\over2}t^2+\mu\,\varepsilon_{\rm TF}:\varepsilon_{\rm TF}
 +O(|\varepsilon|^3),\qquad K=\lambda+{2\mu\over3}.               \tag{H-015}
\]

Thus the tangent law contains precisely the volumetric and shear channels fixed by
DEFORMATION-001.  For \(F=\mathbf1+\nabla u\),
\(\varepsilon=\operatorname{sym}\nabla u+O(|\nabla u|^2)\), so rigid rotations
drop out.  For the accepted relative representation,
\(\delta C=q_0^{-1}\delta q+O(\delta q^2)\); hence (H-014) is compatible with
V11 weak-field tensor kinematics and the retained local Minkowski/Lorentz limit.
This is a compatibility statement only: the absent one-metric map and dynamics
prevent HYPER-001 from deriving or modifying V11 equations.

## 4. Finite elastic bound

Let \(\Lambda=(\lambda_1,\lambda_2,\lambda_3)\).  The admissible spectral domain
must be a path-connected, permutation-invariant set \(\mathcal D_\lambda\) containing
\((1,1,1)\), and

\[
 \mathcal D_C=\{Q^\sharp\operatorname{diag}(\Lambda)Q:
 Q^\sharp Q=\mathbf1,\ \Lambda\in\mathcal D_\lambda\}.             \tag{H-016}
\]

A finite, nondegenerate elastic bound requires the closure to satisfy, for some
symbolic bounds not fixed by the inputs,

\[
 0<\underline\lambda\leq\lambda_A\leq\overline\lambda<\infty,
 \qquad A=1,2,3.                                                    \tag{H-017}
\]

The actual domain need not be a box; (H-017) only states compact containment in
the positive spectral cone.  It excludes loss of rank, orientation reversal,
unbounded stretch, and singular \(C^{-1}\) or \(\log C\).  Its image
\(\mathcal D_I\) under (H-002) is the bounded admissible invariant region and
must still obey the cubic-root positivity conditions implicit in
\(\mathcal D_\lambda\); arbitrary triples \((I_1,I_2,I_3)\) are not admissible.

There are exactly three mathematical endpoint classes consistent with a local
energy description:

1. **Hard constraint:** finite interior \(\Phi\), with the lower-semicontinuous
   extension \(W=+\infty\) outside the closed admissible set.
2. **Barrier:** \(W(C_n)\to+\infty\) as an interior sequence approaches every
   forbidden boundary component.  Its stress or tangent may diverge earlier; a
   particular rate is not derivable.
3. **Regular finite endpoint:** \(W\), and possibly \(P_C\), has a finite one-sided
   limit.  This prevents crossing only when an independent admissibility constraint
   is imposed; smooth continuation alone would not encode a finite bound.

Thus a state inaccessible under arbitrary finite energetic loading requires a hard
constraint or coercive barrier.  A stress plateau, bounded energy, vanishing
tangent, or specific singular exponent is **not** required by a finite bound and
would be an invented saturation law.  V11's \(R_{\max}\), activation function,
\(\epsilon_{0,T}\), and \(\Omega_\sigma(a)\) do not identify
\(\underline\lambda\), \(\overline\lambda\), or a path in \(\mathcal D_\lambda\).

## 5. Convexity and stability

### Necessary equilibrium conditions

Positive stored energy and a stable unloaded equilibrium require

\[
 W(C)\geq0,quad W(\mathbf1)=0,quad DW(\mathbf1)=0.                \tag{H-018}
\]

Strict local stability requires

\[
 D^2W(\mathbf1)[H,H]>0\quad\text{for every nonzero admissible }
 H\in\operatorname{Sym}(3).                                       \tag{H-019}
\]

Using (H-015), this is equivalent to

\[
 \boxed{\;\mu>0,qquad K=\lambda+2\mu/3>0\;}                       \tag{H-020}
\]

(positive shear and bulk tangent stiffness).  Semidefinite inequalities give
neutral modes, not strict stability.  The symbol \(\lambda\) here is the Lamé-type
tangent derivative defined by (H-013), not a cosmological constant and not an
independently introduced parameter.

### Finite-deformation tangent condition

For spectral energy \(w(\Lambda)\), at distinct eigenvalues the exact Hessian is

\[
 D^2W(C)[H,H]=\sum_{A,B}w_{,AB}H_{AA}H_{BB}
 +2\sum_{A<B}{w_{,A}-w_{,B}\over\lambda_A-\lambda_B}H_{AB}^2,       \tag{H-021}
\]

with continuous limiting divided differences when eigenvalues coincide.  Convexity
with respect to the accepted variable \(C\) is therefore equivalent to

\[
 [w_{,AB}]\succeq0,qquad
 {w_{,A}-w_{,B}\over\lambda_A-\lambda_B}\geq0                     \tag{H-022}
\]

throughout a convex spectral domain; replace \(\succeq\) by \(\succ\) for positive
tangent stiffness.  On a nonconvex admissible domain, (H-022) is interpreted
locally along admissible segments.  Global strict convexity of \(W(C)\) is a
sufficient uniqueness/stability condition, not forced by objectivity.  Convexity
in the raw \(F\), rank-one convexity, polyconvexity, and hyperbolicity belong to a
chosen \(F\)-based kinetic completion and are not silently substituted here.

### Unloading and perturbations

Because (H-001) is single-valued and contains no history/rate variables,
quasistatic loading and unloading follow the same path and recover the stored
energy; hysteresis and dissipation are absent.  They cannot be added without new
state structure.  For a body region \(B\), finite-energy perturbations satisfy

\[
 \Delta\mathcal E=\int_B[W(C+\delta C)-W(C)]\,dV_0<\infty.          \tag{H-023}
\]

Uniform positive tangent bounds on a compact interior subset give
\(c\|\delta C\|_{L^2}^2\leq2\Delta\mathcal E+o(\|\delta C\|^2)\)
for some symbolic \(c>0\).  Since the minimal energy has no gradients, it controls
an \(L^2\) deformation perturbation but supplies no \(H^1\) regularity or interfacial
penalty.  Boundary-reaching perturbations have infinite energy in the hard/barrier
classes; in the regular-endpoint class they must be excluded by the state constraint.

## 6. PBUF compatibility audit

| Requirement | Result |
|---|---|
| FOUNDATION-001 FP-1 | One local energy is assigned to the one continuous medium; no substrate is added. |
| FP-2 / emergent gravity | (H-007) supplies an elastic response that may feed a later emergence map; it is not mislabeled as a gravitational field equation. |
| FP-3 / emergent time | \(W\) is statewise and rate-free. An ordering parameter is needed only to label a succession of configurations; no temporal material eigenvalue is introduced. |
| FP-4 | \(C\) describes the one occupied configuration relative to the accepted reference representative; the configuration space is not reified. |
| FP-5 / V11 | The unloaded local state and symmetric tangent kinematics admit the retained Minkowski, Lorentz, and weak-field operational limits. No V11 equation is changed. |
| FP-6 | \(\Phi\), its derivatives, domain endpoints, \(\lambda,\mu,K\) remain symbolic constitutive data; no new independently adjusted constant is asserted. |
| Isotropy/objectivity | Equations (H-003)--(H-004) reduce dependence exactly to (H-002). |
| Covariance | \(W\) is a scalar of the material endomorphism; the integrated energy uses the covariant reference measure. |
| Rank-three ontology | Only three material eigenvalues/invariants occur; the rank-four clock/ruler representation is not primitive here. |

Compatibility does not equal closure.  V11 cannot fix \(\Phi\) from its homogeneous
saturation history without a reference-state map, metric identification, and
homogeneous reduction.  Accordingly this milestone ends at the functional and its
response.

## 7. Dependency graph toward FIELD-001

```text
FOUNDATION-001 FP-1--FP-6       V11 operational limits
              \                       /
               v                     v
       DEFORMATION-001: rank-3 objective C and R0
                            |
                            v
        HYPER-001: invariants (I1,I2,I3), domain D_C
                            |
                            v
      W=Phi(I1,I2,I3) + normalization/bound/stability gates
                            |
                            v
                P_C=dW/dC and tangent A
                            |
          +-----------------+------------------+
          |                                    |
          v                                    v
  missing metric/clock-ruler map       missing kinetic/action measure
          |                                    |
          +-----------------+------------------+
                            v
                         FIELD-001
       (variation, source projection, evolution/field equations)
```

The two “missing” nodes are downstream closure gates, not additions made by
HYPER-001.

## 8. Completion statement

The minimal invariant family is (H-001), its exact constitutive response is
(H-007), its linear elastic tangent is (H-014), its finite-bound alternatives are
classified in Section 4, and its stability conditions are (H-018)--(H-022).
No unique scalar function within this family can be selected from the accepted
inputs.  Claiming otherwise would introduce precisely the free constitutive
assumption prohibited by the mission.
