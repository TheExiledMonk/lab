# PBUF CONSTITUTIVE-SELECTION-001 — Native selection audit

## Decision

The frozen architecture **does not uniquely select** any of the three remaining completions. It reduces Branch A to a constrained scalar invariant function and Branch C to a boundary/asymptotic condition, while Branch B remains a choice among fundamentally different communication operators. The governing-equation readiness outcome is therefore **C: multiple fundamentally different governing equations remain**.

This is an internal implication audit only. It introduces no ontology, constants, coefficients, observational fit, or weak-lensing/V11 change.

## Constitutive elimination table

| Branch | Status | Smallest surviving family | Exact unresolved content |
|---|---|---|---|
| A | Reduced | Phi on the frozen invariant domain, with fixed reference 2-jet and admissibility inequalities | all higher-order invariant dependence and boundary/asymptotic profile |
| B | Undetermined | one admissible communication operator L_comm: balance divergence, positive gradient variation, or symmetric causal nonlocal variation | operator type, differential order/kernel, and boundary data |
| C | Reduced | one boundary completion of Phi on the frozen admissible spectral domain: hard extended-value boundary, interior blow-up barrier, or stable continuation when the domain is unbounded | finite-boundary energy behavior or unbounded-domain asymptotic growth law |

No whole branch is **Eliminated**, and no branch is **Selected**. “Reduced” means the frozen constraints remove inadmissible members but leave inequivalent completions. “Undetermined” for B records that even the mechanism type is not selected.

## Branch A — stored-energy completion

On the frozen rank-three branch,

\[
 W(C)=\Phi(I_1,I_2,I_3),\qquad i_0=(3,3,1).
\]

The accepted inputs require interior differentiability sufficient for stress and tangent, \(\Phi(i_0)=0\), \(D W({\bf1})=0\), the frozen positive weak-field tangent, lower boundedness, and ellipticity/stability on the declared propagation domain. They do not impose global convexity: FOUNDATION-001's one occupied configuration is not a uniqueness-of-minimizer theorem.

Nonuniqueness is constructive. If \(\Phi_0\) is admissible near \(i_0\), then

\[
 \Phi=\Phi_0+R,\qquad R(i_0)=D R(i_0)=D^2R(i_0)=0,
\]

has the same frozen weak-field tangent. Infinitely many smooth invariant remainders exist; small compactly supported remainders preserve strict inequalities on a compact stable subdomain. Thus wave support and weak-field matching cannot determine nonlinear higher derivatives. Branch A is **Reduced**, not selected.

## Branch B — communication mechanism

For a local realization, variation of \(\int\Phi(C)\,dV_0\) and BALANCE-001 produce the balance-divergence operator \(\operatorname{{Div}}P\), which is already sufficient for neighbour communication. Two other frozen-admissible mechanisms are

\[
 {\delta\over\delta q}\int[\Phi(C)+\Psi(\nabla C)]dV_0,
 \qquad
 {\delta\over\delta q}{1\over4}\iint\Delta C(x):K(x,y):\Delta C(y)\,dxdy.
\]

A positive gradient sector and a symmetric positive, causally admissible kernel can support stable waves using only the existing state. Conversely, continuous wave propagation only constrains the linear symbol: required eigenvalues must be real/nonnegative and the resulting characteristics must match the effective V11 cone in its regime. It does not invert a symbol into a unique local, gradient, or integral operator. METRIC-001 explicitly leaves ultralocal, finite-jet, and causal nonlocal dependence open. Branch B is **Undetermined**.

No-new-constants does not prove balance-only communication: it forbids an independently adjustable length/coupling, but parameter-free operators or scales derived from already frozen data are not logically excluded. Selecting balance-only would require a separately frozen minimal-locality axiom.

## Branch C — large deformation and the Planck bound

Let \({\cal D}_C\) be the frozen admissible spectral domain. Finite elastic capacity asserts a boundary of admissible states, but does not imply a unique energy limit. At least two inequivalent completions implement the same capacity:

1. a hard constraint, with finite smooth \(\Phi\) in the interior and extended value \(+\infty\) outside the closed admissible set;
2. an interior barrier, \(\Phi(C_n)\to+\infty\) as \(C_n\to\partial{\cal D}_C\).

If the admissible domain is instead unbounded, coercivity
\(\|C\|+\|C^{-1}\|\to\infty\Rightarrow\Phi(C)\to\infty\) is a sufficient confinement condition, but it is neither a finite-capacity theorem nor equivalent to a finite barrier. Smooth hyperbolicity is required in the operational interior only. Branch C is therefore **Reduced** to a boundary/asymptotic choice, not selected.

## Frozen-constraint reports

### Ontology

One three-dimensional continuous medium fixes an objective isotropic response of the existing state and excludes independent fibres, particles, lattices, phases, or memory variables. Emergent gravity says that response must later feed the effective gravitational description, but supplies no source map or formula for \(\Phi\). Emergent time rules out inserting fundamental-time or rate dependence into the equilibrium energy. None selects a communication operator or boundary profile.

### Wave medium

Continuous waves require a differentiable tangent at the reference, positive acoustic response in required polarizations, spatial communication, and stability on the propagation domain. Convexity is sufficient in some realizations but is not necessary as a global condition and does not follow from wave existence. Gradient and integral mechanisms may introduce dispersion; admissibility requires stability and V11-compatible low-energy causal behavior, not zero dispersion at every scale.

### Emergent metric

The metric map may be a functional kernel or a finite-jet natural operator. It requires covariance, objectivity, causal consistency, Lorentzian signature, and nondegeneracy on its operational domain. Those are output constraints on a chosen closure, not a preference among balance, gradient, and integral communication.

### Duration

The constitutive evolution must permit positive, additive, reparametrization-invariant clock accumulation, stable propagation, and V11 proper-duration matching. Conservative \(\Phi\) permits path-independent recovery on one elastic branch. Dissipation is not selected, and irreversible memory would require unauthorized extra state unless derived from the complete existing \(q\). These restrictions again do not identify \(\Phi\), \(L_{{\rm comm}}\), or the boundary law.

## Minimal remaining freedom and readiness

The smallest exact closure still required is:

1. **one scalar function** \(\Phi\), or equivalently its remainder beyond the frozen reference 2-jet, together with one boundary/asymptotic condition on the already accepted domain; and
2. **one communication operator** \(L_{{\rm comm}}\), including the boundary/domain data needed to define it, selected from the surviving balance, gradient, or integral realizations.

Branch C can be encoded in the domain/extended-value definition of \(\Phi\), so it need not be counted as a third independent constitutive object. Nevertheless, different choices in both items change derivative order, boundary data, dispersion, and nonlinear response. A unique native governing equation does not yet follow. Future derivation therefore remains at outcome **C**, not A or B.

## Traceability and logical boundary

The machine-readable constraint catalogue cites the frozen source class for every restriction. Counterexamples are used only as mathematical constructions inside the already authorized invariant/operator families; no external constitutive theory is imported. Failure to select is the result, rather than a license to prefer a named model.
