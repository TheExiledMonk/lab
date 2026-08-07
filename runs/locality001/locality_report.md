# PBUF LOCALITY-001 — Native sufficiency of local constitutive communication

## Final decision

**Outcome B.** Local variational communication is sufficient for every frozen constitutive requirement. Gradient, integral, and hybrid mechanisms remain mathematically admissible enrichments, but none satisfies a requirement unavailable to local balance divergence. They are not necessary. The independent communication branch is closed for the minimal native equations, without declaring enriched models impossible.

No ontology, field, coefficient, length, kernel horizon, fit, V11 change, or weak-lensing change is introduced.

## 1. Local closure audit

For the placement realization authorized by DEFORMATION-001,

\[
F=\operatorname{Grad}y,\quad C=F^\sharp F,\quad
\mathcal E[y]=\int_{\mathcal B_0}\Phi(I_1,I_2,I_3)\,dV_0,
\]

HYPER-001 gives \(P_C=D_CW\). The chain rule gives \(P_F=2FP_C\), and

\[
\delta\mathcal E
=-\int_{\mathcal B_0}\operatorname{Div}P_F\cdot\delta y\,dV_0
+\int_{\partial\mathcal B_0}(P_FN)\cdot\delta y\,dA_0. \tag{L-001}
\]

A pointwise density in \(C\) is therefore not ultralocal in \(y\): \(C\) already contains a spatial derivative. Equation L-001 supplies neighbour coupling, stress transmission, and traction. BALANCE-001 supplies localization and the frozen kinetic/duration structure supplies the temporal part,

\[
\mathcal K_{\rm frozen}[y]+\operatorname{Div}P_F=b. \tag{L-002}
\]

Variation does not supply inertia or sources; those are separate frozen structures, not a missing communication law. Exactness \(P_C:dC=dW\) and the stable reference give conservative restoring response and quasistatic recovery; a trajectory still requires the frozen evolution closure.

### Necessary notation correction

The literal \(\operatorname{Div}(\partial\Phi/\partial C)\) is generally type-incomplete. \(P_C\) is conjugate to \(C\), whereas placement balance takes the divergence of \(P_F=2FP_C\). Generally the force is \((D_yC)^*P_C\), reducing to \(-\operatorname{Div}P_F\) here. The theorem is proved with this chain rule, not by identifying the two stresses.

## 2. Necessity test

| Mechanism | Classification | Unique frozen requirement |
|---|---|---|
| balance divergence | minimum sufficient | all communication requirements |
| positive gradient | optional enrichment; sometimes a kernel long-wave approximation | none |
| integral kernel | optional enrichment | none |
| hybrid | optional composition | none |

Gradient energy may regularize and disperse but raises differential order and needs extra boundary data; its usual length is unfrozen. A kernel requires unfrozen support, causal prescription, and horizon. Admissibility is not necessity. This sharpens CONSTITUTIVE-SELECTION-001: non-uniqueness in the broad admissible class is not independent necessary freedom under the reversed burden of proof.

## 3. Gradient audit

Intrinsic gradient dependence is **optional**, and redundant relative to the minimum for communication, balance, recovery, and wave existence. When selected it is a genuine higher-order enrichment, not identical to local elasticity. Some regular kernels admit it as a long-wave approximation, but not all and not exactly at finite wavelength. No frozen milestone requires \(\nabla C\), an intrinsic length, fourth-order equations, couple stress, or extra boundary conditions. CONSTITUTIVE-PRINCIPLES-001 explicitly classifies it as optional.

## 4. Nonlocal audit

- FOUNDATION-001: continuity of one medium does not imply two-point dependence.
- METRIC-001: admits local finite-jet and functional maps and selects neither.
- BALANCE-001: supplies a local divergence template, not an integral demand.
- DURATION-001: constrains clocks/order/causality, not a horizon or kernel.
- CONSTITUTIVE-PRINCIPLES-001: makes nonlocality optional and balance divergence minimal.

Thus no frozen principle requires integral nonlocality.

## 5. Wave audit

At a homogeneous reference, a plane wave in the linearization of L-002 obeys

\[
\rho\omega^2a=Q(n)a\,|k|^2,\qquad
Q_{ik}(n)=A_{iJkL}n_Jn_L. \tag{L-003}
\]

Positive frozen inertia and acoustic eigenvalues give continuous longitudinal/shear waves. Distinguish:

1. **existence:** inertia + divergence + an acoustic mode;
2. **stability:** acoustic positivity/strong ellipticity on the propagation domain;
3. **dispersion:** phase speed dependence on \(|k|\).

Local homogeneous elasticity gives \(\omega=c_\alpha(n)|k|\), so waves exist without dispersion. Gradient/nonlocal enrichments can add dispersion but it is not frozen; higher order also does not guarantee stability.

## 6. Metric audit

METRIC-001 requires covariance, Lorentzian signature, nondegeneracy, V11 clock/ruler matching, and cone agreement. Local waves supply characteristics for that gate. It derives no kernel, derivative order, or separate mechanical communication law. Its still-unselected metric map is a separate closure issue.

## 7. Planck-bound audit

Finite capacity restricts pointwise spectra and endpoint behavior. A hard extended-value domain or local barrier \(\Phi(C)\to+\infty\) enforces it. Large stress/tangent near the bound is amplitude dependence, not spatial-gradient dependence. No frozen theorem turns a state-space boundary into a length or kernel.

## 8. Minimal communication theorem

**Theorem.** Assume the accepted placement realization, local invariant energy, BALANCE-001 localization, and already-required kinetic/duration closure. If \(\Phi\) passes frozen reference stability and acoustic positivity, its variation supplies neighbour communication, traction/stress transmission, restoring response, balance, and waves. Metric and finite-capacity gates add no gradient/kernel premise. Hence no frozen constitutive requirement needs communication beyond local variational divergence.

**Proof.** L-001 proves communication/traction; frozen balance gives L-002; exactness plus the stable reference gives recovery; linearization and acoustic positivity give L-003 and stable waves with frozen inertia. The source audits show metric, duration, and capacity constrain outputs, evolution, or state domain, not operator type. This exhausts the frozen catalogue. ∎

The earliest implication would fail only if \(C\) were an independent point variable and \(D_qC\) contained no spatial derivative; then no divergence follows. The candidate explicitly assumes the placement/local-divergence realization, so this is a scope condition, not a frozen counterexample.

## 9. Counterexample catalogue

Explicit attempted failures and constructions are machine-readable in counterexample_catalogue.json. L-001, L-003, and the local hard-domain/barrier construction refute the proposed failures. The decision is based on implication proofs, not merely failure to find a counterexample.

## 10. Governing-equation implications

The **constitutive** part now depends only on selecting \(\Phi\), including its frozen domain/boundary completion. No independent \(L_{\rm comm}\) remains. The complete problem does not depend on \(\Phi\) alone: frozen inertia, sources, initial/boundary data, and METRIC-001's map remain separate.

## Closure

**Outcome B:** communication freedom is eliminated from minimal constitutive development; enriched architectures remain optional.
