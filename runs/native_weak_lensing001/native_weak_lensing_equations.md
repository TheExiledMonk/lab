# PBUF NATIVE-WEAK-LENSING-001 — First Complete Native Forward Weak-Lensing Equation

## 0. Status, scope, and result

This milestone freezes no new ontology and changes no frozen constitutive
equation.  It selects one member of each previously open mathematical family,
explicitly as a closure hypothesis, so that the static forward chain from a
supplied baryonic mass distribution to weak-lensing observables is executable.

The completed chain is

\[
 \rho_b\longmapsto b\longmapsto y\longmapsto C
 \longmapsto \Phi\longmapsto g^{\rm eff}
 \longmapsto (x^\mu,k^\mu,{\cal D})
 \longmapsto(\kappa,\gamma_1,\gamma_2,\mu,g_1,g_2).
 \tag{NWL-001}
\]

Here \(y\), rather than \(C\), is the unknown in the native placement solve;
\(C\) is derived from it.  The potential \(\Phi\), metric, rays, and Jacobi map
are derived computational objects, not new ontological state variables or new
particles.  The model is static.  No kinetic law, cosmology, matter-medium
derivation, or interaction derivation is asserted.

All equations tagged **Frozen** are copied from, or are direct typed uses of,
the authoritative milestones named by the mission.  All added choices are
tagged **Closure Hypothesis**.  A closure is a predictive postulate here, not a
claim of uniqueness or derivation.

## 1. Inputs, unknowns, conventions, and domains

### 1.1 Required numerical inputs

The solver receives:

1. a nonnegative baryonic rest-mass density
   \(\rho_b\in L^1(\Omega)\), or a finite nonnegative mass measure, in physical
   reference coordinates \(X\in\Omega\subseteq\mathbb R^3\);
2. the frozen \(K_0>0\), \(\mu_0>0\), admissible domain
   \({\cal D}_C\Subset\operatorname{Sym}^+(3)\), and constants \(c,G_N\) in
   the retained operational unit system;
3. lens, source, and observer events, plus the observer tetrad and the source
   angular-size or ellipticity convention needed for the requested observable;
4. either the isolated boundary conditions below or explicitly supplied
   external displacement/traction and optical-background boundary data; and
5. a mesh, quadrature rule, solver tolerances, and a declared weak-lensing
   angular field.  These are numerical controls, not physical closures.

The reference and effective spatial coordinates are identified by the unloaded
Cartesian representative.  Indices \(I,J\) are reference spatial indices,
\(i,j\) are effective spatial indices, and \(\mu,\nu=0,1,2,3\).  We use
signature \((-+++ )\), \(x^0=ct\), and \(\Delta=\partial_i\partial_i\).

### 1.2 Forward unknowns and derived fields

The only native BVP unknown is the placement

\[
 y:\Omega\to\mathbb R^3,\qquad u(X):=y(X)-X.                 \tag{NWL-002}
\]

The optical stages solve successively for the scalar \(\chi\), scalar
\(\Phi\), ray phase-space variables \((x^\mu,k^\mu)\), and the \(2\times2\)
Jacobi map \({\cal D}_{AB}\).  These are outputs of declared maps, not
independent medium sectors.

For an isolated calculation the mathematical domain is \(\mathbb R^3\), with

\[
 u\to0,\quad C\to\mathbf1,\quad P_F\to0,\quad
 \chi\to0,\quad\Phi\to0\qquad(|X|\to\infty).                \tag{NWL-003}
\]

A finite-domain implementation must impose Dirichlet/traction data converged
against domain enlargement.  For \(\chi\), use the free-space Green operator
or its boundary trace induced by the same isolated exterior.  Arbitrary
zero-boundary truncation is an approximation and must be convergence-tested.

## 2. Frozen native mechanics

### 2.1 Kinematics and constitutive law — Frozen

On the authorized rank-three placement branch,

\[
 F=\operatorname{Grad}_0y=\mathbf1+\operatorname{Grad}_0u,
 \qquad C=F^TF,
 \qquad E={1\over2}(C-\mathbf1),                                \tag{NWL-004}
\]

\[
 t=\operatorname{tr}E,\qquad E_{\rm TF}=E-{t\over3}\mathbf1,    \tag{NWL-005}
\]

\[
 W_A(C)=
 \begin{cases}
 {K_0\over2}t^2+\mu_0 E_{\rm TF}:E_{\rm TF},
     &C\in\overline{{\cal D}_C},\\
 +\infty,&C\notin\overline{{\cal D}_C},
 \end{cases}                                                     \tag{NWL-006}
\]

\[
 P_C={K_0\over2}t\mathbf1+\mu_0E_{\rm TF},\qquad
 P_F=2FP_C=F\bigl(K_0t\mathbf1+2\mu_0E_{\rm TF}\bigr).         \tag{NWL-007}
\]

No coefficient or nonlinear remainder is added to this law.

### 2.2 Static balance — Frozen

For a regular native load, the strong equation is

\[
 -\operatorname{Div}_0P_F=b.                                    \tag{NWL-008}
\]

For distributional mass inputs the definitive equation is the weak form

\[
 \int_\Omega P_F(y):\operatorname{Grad}_0\eta\,dV_0
 =\langle b,\eta\rangle
 +\int_{\Gamma_N}\bar t\cdot\eta\,dA_0                         \tag{NWL-009}
\]

for every admissible zero-Dirichlet test field \(\eta\).  The solve is
restricted to the connected, uniformly strongly elliptic component containing
\(C=\mathbf1\).  A state reaching \(\partial{\cal D}_C\), losing
\(\det F>0\), or leaving that elliptic component is reported as outside this
solver's operational branch, not silently continued.

## 3. Closure hypotheses

### Closure Hypothesis CH-001 — Static isolated placement realization

**Mathematical definition.**  Use (NWL-002)--(NWL-003), quotient rigid modes by
the asymptotic normalization, and choose the locally unique equilibrium branch
continuously connected to \(y=X\) as \(\rho_b\to0\).

**Physical interpretation.**  The observed lens is a quasistatic loading of
the one medium relative to its unloaded local background.

**Reason for selection.**  It is the minimal realization already authorized
for the frozen local elastic problem and needs no kinetic, dissipative, or new
state variable.

**Expected observational consequence.**  The prediction is instantaneous and
history-free.  Rapidly evolving lenses, retardation, hysteresis, and wave
transients are outside the declared domain.

### Closure Hypothesis CH-002 — Baryonic rest-energy gradient load

**Mathematical definition.**  The supplied baryonic density is projected into
the placement dual by

\[
 \boxed{\ b=c^2\operatorname{Grad}_0\rho_b\ },                    \tag{NWL-010}
\]

in distributions.  Equivalently,

\[
 \boxed{\ \langle b,\eta\rangle
 =-c^2\int_\Omega\rho_b\,\operatorname{Div}_0\eta\,dV_0\ }.     \tag{NWL-011}
\]

Equation (NWL-011) is used for particles and discontinuous density maps and
does not numerically differentiate the data.  Supplied external traction is
added only through the separate boundary term in (NWL-009).

**Physical interpretation.**  A baryonic rest-energy density gradient applies
a purely longitudinal compressive native load.  It is one explicit
matter-to-load hypothesis, not a derived interaction law.

**Reason for selection.**  It is local, linear in the supplied mass measure,
translation/rotation covariant, gives a genuine placement covector, uses no
new coefficient, and has the mandatory units
\([c^2\nabla\rho_b]={\rm N\,m^{-3}}\).

**Expected observational consequence.**  At weak deformation it sources only
the longitudinal native mode.  Shape dependence still enters through the
three-dimensional mass distribution; there is no independent vector or
curl-load response.  A failure correlated with transverse material structure
would falsify CH-002.

### Closure Hypothesis CH-003 — Longitudinal deformation-to-potential map

Define the frozen longitudinal tangent modulus

\[
 M_0:=K_0+{4\mu_0\over3}.                                      \tag{NWL-012}
\]

**Mathematical definition.**  From the solved displacement compute

\[
 \Delta\chi=\operatorname{Div}_0u,\qquad \chi\to0
 \quad(|X|\to\infty),                                          \tag{NWL-013}
\]

or, identically in free space,

\[
 \chi(X)=-{1\over4\pi}\int_{\mathbb R^3}
 {\operatorname{Div}u(X')\over|X-X'|}\,d^3X'.                  \tag{NWL-014}
\]

The native optical potential is

\[
 \boxed{\ \Phi(X)=-{4\pi G_NM_0\over c^2}\,\chi(X)\ }.        \tag{NWL-015}
\]

Thus CH-003 is one support-controlled nonlocal member of the frozen admissible
metric-map family.  Its support and boundary prescription are explicit.

**Physical interpretation.**  Optical clock/ruler response is controlled by
the accumulated longitudinal deformation of the medium.  \(\chi\) is the
Helmholtz scalar computed from the native placement, not an additional
physical substance.

**Reason for selection.**  A scalar input cannot covariantly select a nonzero
vector load without a derivative, and the resulting longitudinal elastic
response must be accumulated once to recover a monopolar long-range optical
field.  The coefficient in (NWL-015) is fixed, not fitted, by the required
weak-field normalization in Section 7.1.

**Expected observational consequence.**  Linear weak lensing has the standard
inverse-distance monopole and is independent of \(K_0,\mu_0\) after their
cancellation.  Finite native deformation produces definite nonlinear,
shape-dependent departures because (NWL-008) remains geometrically nonlinear.

### Closure Hypothesis CH-004 — One exponential effective metric

**Mathematical definition.**  Select exactly one optical representation:

\[
 \boxed{
 ds_{\rm eff}^2=-e^{2\Phi/c^2}(dx^0)^2
                 +e^{-2\Phi/c^2}\delta_{ij}dx^idx^j .}
                                                                    \tag{NWL-016}
\]

This formula, with \(\Phi\) from (NWL-015), is the complete selected
medium-to-metric map.  No second metric, slip field, refractive tensor, or
independent lapse/shift is introduced.

**Physical interpretation.**  Longitudinal native deformation changes ideal
clock and isotropic ruler response oppositely while preserving a single
universal Lorentzian cone.

**Reason for selection.**  The metric is Lorentzian and nondegenerate for every
finite real \(\Phi\), reduces to Minkowski at the unloaded boundary, uses no
additional coefficient, and has the V11-compatible weak form

\[
 g_{00}=-1-{2\Phi\over c^2}+O(\Phi^2/c^4),\qquad
 g_{ij}=\left(1-{2\Phi\over c^2}\right)\delta_{ij}
          +O(\Phi^2/c^4).                                      \tag{NWL-017}
\]

**Expected observational consequence.**  The weak gravitational-slip ratio is
one and light responds to twice the single-potential contribution.  The exact
exponential form also makes higher-order lensing predictions; those terms are
closure predictions, not frozen PBUF consequences.

### Closure Hypothesis CH-005 — Universal null-geodesic bundle

**Mathematical definition.**  Light follows affinely parameterized null
geodesics of (NWL-016), and infinitesimal image distortion follows its Jacobi
equation, as written completely in Section 5.

**Physical interpretation.**  Continuous optical waves propagate along the
characteristics of the single effective metric; rays are the geometric-optics
characteristics of that continuous propagation.

**Reason for selection.**  This is the one-representation realization of the
frozen universal null-cone and duration conditions.  It adds no photon-medium
interaction law or particle species.

**Expected observational consequence.**  Lensing is achromatic in geometric
optics and obeys reciprocity for the static metric.  Frequency-dependent shear
would falsify this closure within its wavelength-resolution domain.

### Closure Hypothesis CH-006 — Observable convention

**Mathematical definition.**  The primary catalogue prediction is reduced
shear

\[
 (g_1,g_2)={1\over1-\kappa}(\gamma_1,\gamma_2),                  \tag{NWL-018}
\]

with magnification \(\mu=[(1-\kappa)^2-|\gamma|^2]^{-1}\).  For
an unbiased population of randomly oriented, intrinsically elliptical sources
in the usual complex-ellipticity convention, predict
\(\langle e_{\rm obs}\rangle=g_1+ig_2\) in the noncritical weak regime.

**Physical interpretation.**  Shape catalogues measure reduced shear rather
than shear directly.

**Reason for selection.**  It supplies the previously missing final
measurement map and is directly comparable to standard weak-lensing data.

**Expected observational consequence.**  The model predicts both image shape
and magnification consistently.  Near \(\det{\cal A}=0\), the weak catalogue
formula is invalid even though the ray/Jacobi solver remains defined.

## 4. Complete native field problem

Combining the frozen mechanics with CH-001 and CH-002 gives the closed native
equilibrium problem

\[
 \boxed{
 -\operatorname{Div}_0\!\left[
 F\left(K_0t\mathbf1+2\mu_0E_{\rm TF}\right)\right]
 =c^2\operatorname{Grad}_0\rho_b,}
                                                                    \tag{NWL-019}
\]

with (NWL-003)--(NWL-005), \(C(X)\in{\cal D}_C\), and the weak meaning
(NWL-009)--(NWL-011).  No symbol in (NWL-019) denotes an unspecified map.

The variational implementation minimizes

\[
 {\cal I}[y]=\int_\Omega W_A(F^TF)\,dV_0
 +c^2\int_\Omega\rho_b\,\operatorname{Div}_0u\,dV_0
 -\int_{\Gamma_N}\bar t\cdot y\,dA_0                         \tag{NWL-020}
\]

over the admissible placement class.  Its first variation is (NWL-009) with
(NWL-011).  Continuation in a load factor \(a\in[0,1]\), replacing
\(\rho_b\) by \(a\rho_b\), selects the branch specified by CH-001.

## 5. Native optical propagation

### 5.1 Ray equations

From (NWL-016), compute the inverse metric and Christoffel symbols

\[
 \Gamma^\mu{}_{\alpha\beta}
 ={1\over2}g^{\mu\nu}
 (\partial_\alpha g_{\nu\beta}+\partial_\beta g_{\nu\alpha}
  -\partial_\nu g_{\alpha\beta}).                               \tag{NWL-021}
\]

For each observed angular direction, initialize a past-directed null vector
in the supplied observer tetrad and integrate

\[
 {dx^\mu\over d\lambda}=k^\mu,\qquad
 {dk^\mu\over d\lambda}=-\Gamma^\mu{}_{\alpha\beta}k^\alpha k^\beta,
 \qquad g_{\mu\nu}k^\mu k^\nu=0.                               \tag{NWL-022}
\]

Affine rescaling of \(k\) is fixed numerically by
\(-g_{\mu\nu}u_o^\mu k^\nu=1\) at the observer.  This normalization does not
alter image distortion.

### 5.2 Screen transport and Jacobi map

Choose orthonormal screen vectors \(s_A^\mu\), \(A=1,2\), at the observer,
orthogonal to \(k^\mu\) and observer velocity \(u_o^\mu\).  Parallel transport
them:

\[
 {ds_A^\mu\over d\lambda}
 =-\Gamma^\mu{}_{\alpha\beta}k^\alpha s_A^\beta.                \tag{NWL-023}
\]

Compute curvature from the same metric,

\[
 R^\mu{}_{\nu\alpha\beta}
 =\partial_\alpha\Gamma^\mu{}_{\nu\beta}
 -\partial_\beta\Gamma^\mu{}_{\nu\alpha}
 +\Gamma^\mu{}_{\sigma\alpha}\Gamma^\sigma{}_{\nu\beta}
 -\Gamma^\mu{}_{\sigma\beta}\Gamma^\sigma{}_{\nu\alpha},     \tag{NWL-024}
\]

and the optical tidal matrix

\[
 {\cal R}_{AB}=R_{\mu\nu\alpha\beta}
 s_A^\mu k^\nu k^\alpha s_B^\beta.                             \tag{NWL-025}
\]

Integrate the Jacobi map

\[
 {d^2{\cal D}_{AB}\over d\lambda^2}
 =-{\cal R}_{AC}{\cal D}_{CB},\qquad
 {\cal D}_{AB}(0)=0,\qquad
 {d{\cal D}_{AB}\over d\lambda}(0)=\delta_{AB}.                \tag{NWL-026}
\]

The sign in (NWL-026) is tied to the curvature convention (NWL-024); changing
curvature convention requires changing both together.  Stop the integration
on the supplied source hypersurface.  These equations are the selected native
continuous-wave optical law in geometric-optics form.

### 5.3 Amplification, shear, convergence, and observable

Run the identical Jacobi integration in the declared unloaded/background
metric between the same operational endpoints, yielding scalar background
angular-diameter distance \(D_{A,0}\).  Define

\[
 {\cal A}:={{\cal D}(\lambda_s)\over D_{A,0}}.                   \tag{NWL-027}
\]

Remove any antisymmetric numerical/rotation part only for the standard
convergence-shear decomposition:

\[
 {\cal A}_{\rm sym}={1\over2}({\cal A}+{\cal A}^T)
 =\begin{pmatrix}
 1-\kappa-\gamma_1&-\gamma_2\\
 -\gamma_2&1-\kappa+\gamma_1
 \end{pmatrix}.                                                  \tag{NWL-028}
\]

Therefore

\[
 \boxed{
 \kappa=1-{A_{11}+A_{22}\over2},\quad
 \gamma_1={A_{22}-A_{11}\over2},\quad
 \gamma_2=-{A_{12}+A_{21}\over2}.}                              \tag{NWL-029}
\]

With \(|\gamma|^2=\gamma_1^2+\gamma_2^2\),

\[
 \boxed{
 \mu={1\over\det{\cal A}},\qquad
 (g_1,g_2)={1\over1-\kappa}(\gamma_1,\gamma_2).}                \tag{NWL-030}
\]

Equations (NWL-018) and (NWL-030) are the final predicted observables.

## 6. Weak-lensing fast path and validation identity

The exact production path is Sections 4--5.  In the Born, thin-lens,
small-potential limit, define physical line-of-sight coordinate \(z\), lens
plane coordinate \(\boldsymbol\xi=D_l\boldsymbol\theta\), and supplied
background distances \(D_l,D_s,D_{ls}\).  Then

\[
 \Sigma_b(\boldsymbol\xi)=\int\rho_b(\boldsymbol\xi,z)\,dz,
 \qquad
 \Sigma_{\rm crit}={c^2\over4\pi G_N}{D_s\over D_lD_{ls}},       \tag{NWL-031}
\]

\[
 \kappa(\boldsymbol\theta)
 ={\Sigma_b(D_l\boldsymbol\theta)\over\Sigma_{\rm crit}},       \tag{NWL-032}
\]

and, with \(\psi\) the dimensionless lensing potential,

\[
 \psi(\boldsymbol\theta)
 ={2D_{ls}\over c^2D_lD_s}
 \int\Phi(D_l\boldsymbol\theta,z)\,dz,                          \tag{NWL-033}
\]

\[
 \kappa={1\over2}(\psi_{,11}+\psi_{,22}),\quad
 \gamma_1={1\over2}(\psi_{,11}-\psi_{,22}),\quad
 \gamma_2=\psi_{,12}.                                           \tag{NWL-034}
\]

For a weak, interior-branch run, the exact solver must converge to
(NWL-031)--(NWL-034).  This is a regression test, not a second optical law.
Distances in (NWL-031)--(NWL-033) are supplied operational background data;
this milestone derives no cosmology.

## 7. Internal consistency demonstrations

### 7.1 Weak normalization

Linearize \(y=X+u\).  The frozen equation becomes

\[
 -\mu_0\Delta u-left(K_0+{\mu_0\over3}\right)
 \nabla(\nabla\cdot u)=c^2\nabla\rho_b.                          \tag{NWL-035}
\]

Taking the divergence and using (NWL-003) gives

\[
 \nabla\cdot u=-{c^2\over M_0}\rho_b.                           \tag{NWL-036}
\]

Equations (NWL-013) and (NWL-015) then imply exactly

\[
 \boxed{\ \Delta\Phi=4\pi G_N\rho_b\ }.                        \tag{NWL-037}
\]

Thus \(K_0\) and \(\mu_0\) govern the native deformation but cancel from the
linear optical normalization.  This cancellation is a prediction of the
selected closure chain, not a change to either modulus.

### 7.2 Dimensions

\[
 [C]=[E]=[t]=1,quad [K_0]=[\mu_0]=[P_F]={\rm J,m^{-3}},
\]

\[
 [b]=[c^2\nabla\rho_b]={\rm N\,m^{-3}},\quad [u]={\rm m},
\quad[\chi]={\rm m^2},
\]

\[
 [G_NM_0\chi/c^2]={\rm m^2s^{-2}}=[\Phi],\qquad
 [\Phi/c^2]=1.                                                    \tag{NWL-038}
\]

Hence every exponent and every metric component in (NWL-016) is
dimensionally valid.  \({\cal A},\kappa,\gamma,\mu,g\) are dimensionless.

### 7.3 Frozen locality

The source map (NWL-011) is finite-jet local and distributionally supported on
the baryonic mass.  The elastic constitutive operator is the frozen local
first-gradient operator.  Its elliptic solution depends on regional source and
boundary data, as allowed by LOCALITY-001.  CH-003 is explicitly nonlocal but
support-controlled: its single Green kernel and boundary condition are fully
specified in (NWL-013)--(NWL-014).  Optical propagation uses only the metric
and its derivatives on the ray tube.  Therefore no hidden global functional or
undefined support remains.

### 7.4 Frozen duration and continuous propagation

For timelike curves, (NWL-016) uses the frozen duration rule

\[
 d\tau^2=-c^{-2}g^{\rm eff}_{\mu\nu}dx^\mu dx^\nu.               \tag{NWL-039}
\]

For optical curves, \(d\tau=0\) and (NWL-022) propagates the common null cone.
The affine parameter is a ray parameter, not fundamental time.  The metric is
static only because CH-001 selects one quasistatic order state; no fundamental
time dimension is added to the medium ontology.

### 7.5 Frozen balance and constitutive law

The internal virtual work is exactly \(\int P_F:\nabla\eta\), the external
work is exactly the covector (NWL-011) plus declared boundary traction, and
(NWL-009) equates them.  The load never enters \(W\), \(P_C\), or \(P_F\) as
a constitutive argument.  Internal exchange is not double counted.  No energy,
momentum, or propagation conservation law beyond the frozen balance structure
is asserted.

### 7.6 Objectivity, covariance, and reference behavior

The native energy depends on \(C=F^TF\), so superposed rigid rotations cancel.
CH-002 is Euclidean-covariant in the selected unloaded representative.  The
free-space inverse Laplacian in CH-003 is invariant under translations and
rotations of that representative.  Equation (NWL-016) is a tensor once written
in that chart and transforms by pullback under effective coordinate changes.
At \(\rho_b=0\), the selected branch has
\(u=0,\chi=0,\Phi=0,g^{\rm eff}=\eta\), zero optical tidal matrix, and
\({\cal A}=\mathbf1\).

## 8. Executable computational sequence

1. **Ingest mass.**  Put \(\rho_b\) or the mass measure on the three-dimensional
   reference mesh.  Preserve total mass under deposition.
2. **Assemble load.**  Do not differentiate noisy density.  Assemble
   \(-c^2\int\rho_b\nabla\cdot\eta\) directly in the weak residual.
3. **Solve native equilibrium.**  Continue the load from zero to full strength;
   at each Newton step evaluate (NWL-004)--(NWL-007) and solve (NWL-009).
4. **Audit the branch.**  At all quadrature points require
   \(\det F>0\), \(C\in\operatorname{int}{\cal D}_C\), and positive lifted
   Legendre--Hadamard form.  Stop with an out-of-domain result if any fails.
5. **Extract longitudinal deformation.**  Project \(\nabla\cdot u\) into the
   scalar Poisson space and solve (NWL-013) with the isolated Green boundary.
6. **Build optical geometry.**  Evaluate (NWL-015)--(NWL-016), plus first and
   second spatial derivatives consistent with the FE/Poisson interpolation.
7. **Trace rays.**  For every observed direction integrate
   (NWL-021)--(NWL-023) from observer to source surface while monitoring the
   null constraint.
8. **Trace bundles.**  Integrate (NWL-024)--(NWL-026) on the same rays and
   compute the matched background \(D_{A,0}\).
9. **Form lensing fields.**  Evaluate (NWL-027)--(NWL-030) on the angular grid.
10. **Emit observables.**  Return \(\kappa,\gamma_1,\gamma_2,\mu,g_1,g_2\),
    quality flags, boundary-convergence results, native-domain margins, and
    null/Jacobi integration residuals.
11. **Validate.**  Reduce the mass amplitude successively and verify convergence
    to (NWL-031)--(NWL-037); refine mesh/ray steps and enlarge the domain until
    requested observables meet the declared numerical tolerance.

No computational block in this sequence calls an undefined physical map.

## 9. Full dependency graph

```text
rho_b(X), lens/source/observer data, K0, mu0, D_C, c, G_N
  |
  | CH-002: <b,eta> = -c^2 integral rho_b Div eta
  v
native load covector b + declared boundary traction
  |
  | FROZEN: static balance and W_A; CH-001 branch/boundary
  v
y(X) -> u -> F -> C -> E -> P_C -> P_F
  |
  | CH-003: Delta chi = Div u; Phi = -4 pi G_N M0 chi/c^2
  v
Phi(X)
  |
  | CH-004: one exponential Lorentzian metric
  v
g_eff
  |
  | CH-005: null geodesic + screen + Jacobi equations
  v
ray x(lambda), optical Jacobi map D_AB
  |
  | matched unloaded/background D_A,0
  v
A_AB
  |
  | algebraic decomposition
  v
kappa, gamma_1, gamma_2, magnification
  |
  | CH-006 measurement convention
  v
reduced shear g and mean observed ellipticity
```

## 10. Frozen assumptions versus closure assumptions

### 10.1 Frozen and authoritative

1. FOUNDATION-001 FP-1--FP-6: one continuous medium, emergent gravity and
   time, one complete configuration per state, V11 operational compatibility,
   and no unauthorized fundamental constants.
2. STATE-003: \(q\) is ontically complete but is neither identified with
   placement alone nor proven to be unique Cauchy data.
3. DEFORMATION-001: objective dimensionless rank-three relative deformation
   \(C=F^TF\), with unloaded value \(\mathbf1\).
4. HYPER-001 and CONSTITUTIVE-CONSTRUCTION-001: local isotropic hyperelastic
   response and the selected minimal hard quadratic \(W_A\), with frozen
   \(K_0,\mu_0,{\cal D}_C\).
5. BALANCE-001: local virtual-work/balance structure and no invented universal
   conservation charge.
6. DURATION-001: emergent proper duration and the one effective Lorentzian
   metric/null-cone matching condition.
7. LOCALITY-001, WEAK-LENSING-LOCALITY-001, and LOCAL-STATE-001: regional
   static elastic sufficiency, boundary replacement of the exterior, local
   constitutive support, and the requirement that any added source/metric map
   declare its support.

### 10.2 Closure hypotheses selected only for this forward solver

1. CH-001: quasistatic isolated placement branch and its boundary/gauge rule.
2. CH-002: \(b=c^2\nabla\rho_b\), with its distributional weak definition.
3. CH-003: inverse-Laplacian longitudinal deformation map and fixed
   \(G_N\)-normalized potential.
4. CH-004: the single exponential effective metric (NWL-016).
5. CH-005: null-geodesic and Jacobi optical propagation for continuous waves.
6. CH-006: reduced shear/magnification catalogue convention.

The closures are jointly testable and are not claimed unique.  In particular,
the standard linear lensing normalization follows by construction, whereas
finite-deformation corrections, exact exponential-metric corrections,
achromaticity, and absence of an independent transverse source channel are
genuine observational consequences of this selected closure set.

## 11. Completion statement

Equations (NWL-003)--(NWL-030), together with supplied numerical/problem data,
form a closed forward mathematical model.  They specify the mass-to-load map,
native equilibrium, deformation, deformation-to-metric map and its support,
duration-compatible optical metric, null propagation, bundle distortion, and
observable convention.  Equations (NWL-031)--(NWL-037) provide an independent
weak-limit implementation check.

An independent numerical physicist can implement this solver without another
theoretical milestone.  Empirical validation or rejection of CH-001--CH-006
is subsequent model-testing work; it does not reopen the frozen framework.
