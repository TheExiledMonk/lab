# PBUF MATERIAL-LAB-001 — Comparative Material Laboratory

## Scope and common procedure

The frozen ontology and milestones are inputs, not objects of review. The native deformation remains the objective tensor `C[q,q0]` in the frozen admissible SPD domain. For comparable algebra, each local model is restricted to a smooth reversible ray `C(e)` and tested using `w(e)=W(C(e))`, `sigma=w'`, and `k_t=w''`. Passing a ray test is necessary, not sufficient, for convexity, rank-one convexity, strong ellipticity, or hyperbolicity of a future tensor theory.

Every model receives the same five stages: (1) regularity, domain, positivity, convexity, stability, and equilibrium; (2) native response, energy, tangent, hardening, recovery, neighbor coupling, and waves; (3) frozen-milestone compatibility; (4) geometry, duration, propagation, recovery, finite energy, and stiffening consequences; and (5) downstream readiness. No observations or weak-lensing implementation enter the evaluation.

Common gates are `w(0)=sigma(0)=0`, `k_t(0)=k0>0`, positive energy away from equilibrium, positive tangent, and nondecreasing tangent for progressive hardening. A strict energy minimum gives a restoring tendency, but recovery in time requires an authorized evolution law. Coercive energies make finite-energy sublevels bounded; only a barrier model also has a finite mathematical endpoint.

## Model A — Progressive strain hardening

1. **Physical interpretation.** Resistance rises with accumulated reversible deformation.
2. **Constitutive assumptions.** hyperelastic ray restriction; k_t(e)>0 and k_t'(e)>=0; quartic representative has h>0.
3. **Derived equation.** `sigma=k0 e+h e^3`. For B, variation additionally gives `(1-ell^2 Laplacian)e=s/k0`; boundary conditions are required.
4. **Stored energy.** `w=k0 e^2/2+h e^4/4` on `signed e in R (representative); tensor lift restricted to the frozen SPD domain`.
5. **Stress and tangent.** `sigma=k0 e+h e^3`; `k_t=k0+3h e^2`.
6. **Progressive hardening.** pass. The hypothesis uniquely fixes only monotonicity inequalities unless the displayed representative is stated.
7. **Weak-field limit.** `w=k0 e^2/2+o(e^2)`, `sigma=k0 e+o(e)`.
8. **Recovery.** pass restoring tendency; evolution pending.
9. **Mathematical stability.** pass on ray; tensor proof pending; smoothness: pass.
10. **Strengths.** minimal smooth even hardening representative; strictly convex and coercive.
11. **Weaknesses.** hypothesis selects inequalities, not a unique law; local form supplies no neighbor coupling.
12. **Open questions.** Which invariant tensor lift?; What fixes the hardening scale?.

Neighbor interaction: none implied by the local hypothesis. Wave support: conditional on an authorized kinetic/gradient completion.

## Model B — Wave-equilibrium material

1. **Physical interpretation.** A single continuum rearranges collectively through a positive gradient-energy term.
2. **Constitutive assumptions.** one continuum only; ell>0 is a candidate correlation length; static scalar ray proxy.
3. **Derived equation.** `local sigma=k0 e; microstress xi=k0 ell^2 grad e`. For B, variation additionally gives `(1-ell^2 Laplacian)e=s/k0`; boundary conditions are required.
4. **Stored energy.** `W=int[k0 e^2/2+k0 ell^2 |grad e|^2/2-s e] dV` on `e in H1(Omega), with boundary data making the variational problem well posed`.
5. **Stress and tangent.** `local sigma=k0 e; microstress xi=k0 ell^2 grad e`; `local k_t=k0; Fourier Hessian=k0(1+ell^2|q|^2)`.
6. **Progressive hardening.** fail (amplitude tangent is constant). The hypothesis uniquely fixes only monotonicity inequalities unless the displayed representative is stated.
7. **Weak-field limit.** `w=k0 e^2/2+o(e^2)`, `sigma=k0 e+o(e)`.
8. **Recovery.** pass restoring tendency; evolution pending.
9. **Mathematical stability.** pass on ray; tensor proof pending; smoothness: pass.
10. **Strengths.** explicit neighbor coupling; positive Fourier stiffness; unique static equilibrium under standard boundary conditions.
11. **Weaknesses.** does not progressively harden in amplitude; ell and kinetic closure are underived.
12. **Open questions.** What frozen quantity fixes ell?; What tensor-gradient invariant and kinetic term are admissible?.

Neighbor interaction: intrinsic: Euler-Lagrange equilibrium (1-ell^2 Laplacian)e=s/k0. Wave support: yes conditionally: positive inertia gives omega^2=c0^2(q^2+ell^2 q^4); inertia is not selected here.

## Model C — Exponential hardening

1. **Physical interpretation.** Incremental work grows exponentially in deformation magnitude.
2. **Constitutive assumptions.** b>0; signed-symmetric hyperelastic representative.
3. **Derived equation.** `sigma=k0 sinh(b e)/b`. For B, variation additionally gives `(1-ell^2 Laplacian)e=s/k0`; boundary conditions are required.
4. **Stored energy.** `w=k0[cosh(b e)-1]/b^2` on `signed e in R`.
5. **Stress and tangent.** `sigma=k0 sinh(b e)/b`; `k_t=k0 cosh(b e)`.
6. **Progressive hardening.** pass. The hypothesis uniquely fixes only monotonicity inequalities unless the displayed representative is stated.
7. **Weak-field limit.** `w=k0 e^2/2+o(e^2)`, `sigma=k0 e+o(e)`.
8. **Recovery.** pass restoring tendency; evolution pending.
9. **Mathematical stability.** pass on ray; tensor proof pending; smoothness: pass.
10. **Strengths.** smooth, strictly convex, coercive; rapid progressive stiffening.
11. **Weaknesses.** no finite mathematical endpoint; exponential form and scale are not frozen consequences.
12. **Open questions.** What derives b?; Does a globally stable tensor lift exist on the full domain?.

Neighbor interaction: none implied locally. Wave support: conditional on an authorized kinetic/gradient completion.

## Model D — Finite-extensibility material

1. **Physical interpretation.** A divergent elastic barrier prevents reaching finite extensibility e_star.
2. **Constitutive assumptions.** e_star>0; Gent-type symmetric ray representative.
3. **Derived equation.** `sigma=k0 e/[1-(e/e_star)^2]`. For B, variation additionally gives `(1-ell^2 Laplacian)e=s/k0`; boundary conditions are required.
4. **Stored energy.** `w=-(k0 e_star^2/2) log[1-(e/e_star)^2]` on `|e|<e_star (normalized e_star=1)`.
5. **Stress and tangent.** `sigma=k0 e/[1-(e/e_star)^2]`; `k_t=k0[1+(e/e_star)^2]/[1-(e/e_star)^2]^2`.
6. **Progressive hardening.** pass. The hypothesis uniquely fixes only monotonicity inequalities unless the displayed representative is stated.
7. **Weak-field limit.** `w=k0 e^2/2+o(e^2)`, `sigma=k0 e+o(e)`.
8. **Recovery.** pass restoring tendency; evolution pending.
9. **Mathematical stability.** pass on ray; tensor proof pending; smoothness: pass on open admissible domain.
10. **Strengths.** finite domain with infinite-energy barrier; strict convexity and hardening on each loading branch.
11. **Weaknesses.** endpoint is not derived; ray result does not prove full tensor ellipticity.
12. **Open questions.** Which spectral boundary defines extensibility?; How are compression boundaries controlled?.

Neighbor interaction: none implied locally. Wave support: conditional on an authorized kinetic/gradient completion.

## Model E — Polynomial hardening

1. **Physical interpretation.** A general nonlinear polynomial supplies controllable hardening orders.
2. **Constitutive assumptions.** even energy for signed recovery; a_2m>=0 with at least one positive nonlinear coefficient.
3. **Derived equation.** `sigma=k0 e+sum a_2m e^(2m-1)`. For B, variation additionally gives `(1-ell^2 Laplacian)e=s/k0`; boundary conditions are required.
4. **Stored energy.** `w=k0 e^2/2+sum_(m=2)^N a_2m e^(2m)/(2m)` on `signed e in R for the nonnegative even-energy subclass`.
5. **Stress and tangent.** `sigma=k0 e+sum a_2m e^(2m-1)`; `k_t=k0+sum(2m-1)a_2m e^(2m-2)`.
6. **Progressive hardening.** pass. The hypothesis uniquely fixes only monotonicity inequalities unless the displayed representative is stated.
7. **Weak-field limit.** `w=k0 e^2/2+o(e^2)`, `sigma=k0 e+o(e)`.
8. **Recovery.** pass restoring tendency; evolution pending.
9. **Mathematical stability.** pass on ray; tensor proof pending; smoothness: pass.
10. **Strengths.** systematically extensible; stable for the displayed coefficient restrictions.
11. **Weaknesses.** mixed coefficients need separate domain proof; order and coefficients remain free.
12. **Open questions.** Which truncation is justified?; What constrains coefficients and tensor invariants?.

Neighbor interaction: none implied locally. Wave support: conditional on an authorized kinetic/gradient completion.

## Unweighted comparative matrix

| Property | A | B | C | D | E |
|---|---|---|---|---|---|
| Weak-field recovery | pass (common k0 tangent) | pass (common k0 tangent) | pass (common k0 tangent) | pass (common k0 tangent) | pass (common k0 tangent) |
| Progressive hardening | pass | fail (amplitude tangent is constant) | pass | pass | pass |
| Smooth constitutive response | pass | pass | pass | pass on open admissible domain | pass |
| Stable equilibrium | pass on ray; tensor proof pending | pass on ray; tensor proof pending | pass on ray; tensor proof pending | pass on ray; tensor proof pending | pass on ray; tensor proof pending |
| Wave-support capability | conditional (completion required) | pass (conditional kinetic closure) | conditional (completion required) | conditional (completion required) | conditional (completion required) |
| Recovery | pass restoring tendency; evolution pending | pass restoring tendency; evolution pending | pass restoring tendency; evolution pending | pass restoring tendency; evolution pending | pass restoring tendency; evolution pending |
| Finite-energy compatibility | pass | pass | pass | pass | pass |
| Emergent geometry compatibility | conditional on METRIC-001 map | conditional on METRIC-001 map | conditional on METRIC-001 map | conditional on METRIC-001 map | conditional on METRIC-001 map |
| Governing-equation readiness | conditional; closure slots explicit | conditional; closure slots explicit | conditional; closure slots explicit | conditional; closure slots explicit | conditional; closure slots explicit |
## Frozen PBUF compatibility

All candidates preserve FOUNDATION-001's one-medium ontology and STATE-002/DEFORMATION-001's native objective state. Their weak tangent can satisfy HYPER-001 and ENERGY-PRINCIPLE-001, subject to a full invariant tensor lift and spectral stability proof. They are statewise and introduce no fundamental time, so DURATION-001 is preserved. Emergent geometry remains conditional on a regular METRIC-001 map. Stress is an admissible BALANCE-001 closure ingredient, not a completed balance/evolution equation. CONSTITUTIVE-002 remains unmodified and supplies restrictions rather than a selected formula.

## Derived behavior and readiness

All five have a stable weak equilibrium and finite-energy compatibility. A, C, D, and E progressively stiffen in amplitude; B instead stiffens short wavelengths through neighbor coupling and is the only initial model with native spatial propagation structure. D alone imposes a finite mathematical deformation endpoint. None independently derives duration or the metric, and none closes governing equations without tensor, kinetic, source, boundary, and metric-map choices. Consequently lensing, cosmology, compact-object, and quantum work are only conditionally ready after those closures; no weak-lensing modification is warranted here.

## Comparative ranking without weights

There is no scalar score and no unique winner. The criterion-only result is a partial ordering by distinctive capability:

1. **Broad static hardening set:** A, C, D, and E pass all displayed ray-level hardening gates; none dominates the others absent a selected endpoint or growth principle.
2. **Spatial-interaction leader:** B alone passes native neighbor interaction and conditional wave-support gates, but does not pass amplitude hardening.
3. **Finite-endpoint leader:** D alone supplies finite extensibility through an energetic barrier; the frozen framework does not establish that this extra property must select the law.

This is a capability ranking, not a preference or declaration of physical correctness.

## Recommendations for subsequent development

Preserve all five as separate mechanism classes. Next derive invariant tensor lifts and run spectral Hessian/strong-ellipticity tests on the frozen SPD domain. For B, derive admissible tensor-gradient and kinetic terms and boundary conditions. Keep the weak tangent common, keep model parameters symbolic until independently derived, and only after governing closure evaluate later metric, lensing, cosmology, compact-object, or quantum consequences.
