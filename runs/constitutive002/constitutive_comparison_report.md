# PBUF CONSTITUTIVE-002 — Comparative Evaluation of Admissible Constitutive Laws

## 0. Decision and scope

**Comparative result:** the accepted principles eliminate globally linear elasticity as a complete finite-response law, but do not uniquely select one nonlinear formula. Polynomial, exponential, and finite-extensibility/barrier families can all satisfy CP-1--CP-7 after admissible parameter and domain restrictions. The most natural *class* is a smooth, strictly stable, progressively hardening hyperelastic law with the HYPER-001 weak tangent. If the separately frozen finite elastic bound is enforced energetically, the asymptotic barrier class is the strongest match; Gent-type finite extensibility is a standard representative, not a uniquely derived PBUF law. Polynomial and exponential laws impose no mathematical cutoff, yet finite total available energy makes every finite-energy sublevel physically bounded when their energy is coercive.

No family is declared native. Selection of a formula, deformation path, tensor extension, energy scale, or endpoint remains a later milestone. No governing equation, observation, V11 expression, ontology, microscopic constituent, or coupling is introduced.

## 1. Common mathematical test

PBUF's native variable remains the objective tensor

\[
C[q,q_0]\in {\cal D}_C\subset\operatorname{Sym}^+(3),\qquad C(q_0,q_0)={\bf1},
\]

and every admissible local isotropic energy has the form

\[
W(C)=\Phi(I_1,I_2,I_3),\quad W({\bf1})=0,\quad DW({\bf1})=0.
\tag{C2-001}
\]

To compare functional growth without inventing a scalar ontology, choose any admissible smooth deformation ray (C(e)), with (C(0)={\bf1}), and define

\[
w(e)=W(C(e)),\qquad \sigma(e)=w'(e),\qquad k_t(e)=\sigma'(e)=w''(e).
\tag{C2-002}
\]

Here (e\ge0) is accumulated deformation along that ray. The conclusions must hold on every physical ray for a tensor law to pass globally. The weak tangent is denoted (k_0>0); it is a ray projection of the already symbolic Lamé-type tangent, not a new fundamental constant. Near equilibrium,

\[
w(e)=\tfrac12k_0e^2+o(e^2),\quad \sigma(e)=k_0e+o(e),\quad k_t(0)=k_0.
\tag{C2-003}
\]

CP-5 means (k_t'(e)\ge0). Where (k_t>0), the reciprocal compliance is

\[
{de\over d\sigma}={1\over k_t(e)},\qquad {d\over de}\left({de\over d\sigma}\right)=-{k_t'(e)\over k_t(e)^2}\le0.
\tag{C2-004}
\]

Thus a monotonically *increasing* (de/d\sigma) would describe softening and contradict CP-5. This report evaluates the requested derivative but interprets progressive hardening by (C2-004).

For signed reversible strain, replace odd one-sided hardening laws by their even-energy extension. All formulas below are normalized to the same (k_0). Shape/end-point symbols describe candidate families only; choosing their values would require later authorization under FP-6.

## 2. Candidate definitions and derivations

### A — Linear elastic

\[
w_A(e)=\tfrac12k_0e^2,\qquad \sigma_A=k_0e,\qquad k_A=k_0,
\qquad {de\over d\sigma_A}=k_0^{-1}.
\tag{C2-A}
\]

It is positive, strictly convex, and coercive on an unbounded ray. It has the correct weak response and stable equilibrium. It does not progressively harden: resistance is constant. A finite energy budget (E) restricts reachability to (e\le\sqrt{2E/k_0}), although the mathematical domain remains unbounded. It is admissible only as the required tangent approximation.

### B — Polynomial hardening

The general one-sided family is

\[
\sigma_B(e)=k_0e+\sum_{m=2}^{N}a_m e^m,\qquad
w_B(e)=\tfrac12k_0e^2+\sum_{m=2}^{N}{a_m\over m+1}e^{m+1},
\tag{C2-B1}
\]

with

\[
k_B(e)=k_0+\sum_{m=2}^{N}m a_m e^{m-1}.
\tag{C2-B2}
\]

Nonnegative coefficients give positive tangent and monotone hardening for (e\ge0), provided at least one (a_m>0). The minimum one-sided nonlinear order is quadratic stress/cubic energy. For an unbiased signed strain whose stored energy is smooth and even, the first admissible correction is instead

\[
\sigma_B=k_0e+a_3e^3,\qquad
w_B=\tfrac12k_0e^2+\tfrac14a_3e^4,\quad a_3>0;
\tag{C2-B3}
\]

hence cubic stress/quartic energy is the minimum symmetry-compatible order. Negative leading coefficients cause loss of coercivity or tangent stability; mixed signs require a domain-specific proof. With a positive leading energy coefficient, (w_B\to\infty), so every finite-energy sublevel is bounded but no mathematical deformation endpoint exists.

### C — Exponential hardening

A signed-symmetric representative is

\[
w_C(e)={k_0\over b^2}\,[\cosh(be)-1],\qquad
\sigma_C(e)={k_0\over b}\sinh(be),\qquad
k_C(e)=k_0\cosh(be),
\tag{C2-C}
\]

for candidate shape parameter (b>0). It is positive, strictly convex and coercive; (w_C=\tfrac12k_0e^2+O(e^4)). On (e\ge0), stiffness increases and compliance (1/[k_0\cosh(be)]) decreases. Every finite-energy sublevel is bounded, while the mathematical domain is unbounded. Exponential growth satisfies CP-1--CP-7 naturally, but neither the exponential form nor (b) follows from the frozen architecture.

### D — Asymptotic hardening

Let (e_* >0) denote a candidate finite endpoint and (p>0). A normalized generic barrier is

\[
w_D(e)={k_0e_*^2\over p(p+1)}\left[(1-e/e_*)^{-p}-1-p(e/e_*)\right],
\quad 0\le e<e_*,
\tag{C2-D1}
\]

so that

\[
\sigma_D={k_0e_*\over p+1}\left[(1-e/e_*)^{-p-1}-1\right],
\quad
k_D=k_0(1-e/e_*)^{-p-2}.
\tag{C2-D2}
\]

It has the common weak tangent, positive strictly increasing stiffness, decreasing compliance, strict convexity, and (w_D\to\infty) as (e\uparrow e_*). Its open domain has a mathematical limit (e_*); any finite energy reaches only a smaller (e_E<e_*). This is compatible with HYPER-001's barrier implementation of a finite elastic bound. The endpoint and exponent are not derived, and a tensor extension must erect barriers at every forbidden spectral boundary, including compression/degeneracy boundaries.

### E — Standard finite-elastic comparison families

These are mathematical exemplars, not imported medium physics.

**Gent-type finite extensibility.** Along a normalized ray with (x=e/e_*\in[0,1)),

\[
w_{E,G}(e)=-\tfrac12k_0e_*^2\log(1-x^2),
\quad \sigma_{E,G}={k_0e\over1-x^2},
\quad k_{E,G}=k_0{1+x^2\over(1-x^2)^2}.
\tag{C2-E1}
\]

It has a quadratic weak limit, positive increasing tangent, and a logarithmically divergent barrier at (e_*). Finite energy cannot reach the endpoint. In its established tensor form the limiting invariant must be adapted to PBUF's full admissible spectral domain; a one-invariant Gent law alone does not guarantee volumetric or all-boundary control.

**Hencky/logarithmic strain elasticity.** With logarithmic strain (H=\tfrac12\log C), the standard quadratic form

\[
W_{E,H}(C)=\tfrac12K(\operatorname{tr}H)^2+\mu\,H_{\rm TF}:H_{\rm TF},
\quad K>0,\ \mu>0,
\tag{C2-E2}
\]

is objective, positive and has the required weak tangent. On a pure logarithmic-strain ray it is linear elastic in the chosen coordinate, so it does not progressively harden there. It is coercive as eigenvalues approach zero or infinity, but has no finite mathematical stretch endpoint and global convexity depends on the chosen native variable/domain. It is therefore a useful finite-strain parametrization, not a CP-5 winner.

**Logarithmic finite-extensibility barriers.** Energies of the form

\[
w_{E,L}(e)=-A\log[1-(e/e_*)^r]-\text{terms required to set }w(0)=w'(0)=0
\tag{C2-E3}
\]

can satisfy the same gates when their weak quadratic coefficient is (k_0/2), (r\ge2), and the tangent remains positive and nondecreasing. Gent is the clean (r=2) representative. The label “logarithmic” alone does not guarantee hardening; (C2-E2) and (C2-E3) have materially different endpoint behavior.

## 3. Stored energy, stability, and finite energy

For all families, recovery is an energetic statement, not a derived evolution equation. If external loading is removed and an authorized conservative or dissipative evolution follows decreasing/constant total energy on the same elastic branch, the unique strict minimum at (e=0) supplies a restoring force. The stored energy alone does **not** prove a relaxation rate, damping, asymptotic convergence, or time law. This respects DURATION-001 and EVOLUTION-001.

The identical static admissibility gates are:

1. (w(0)=w'(0)=0), (w(e)>0) for (e>0);
2. (k_t=w''>0) for strict tangent stability;
3. (k_t'\ge0) for CP-5 hardening;
4. coercivity or a hard/barrier domain to bound finite-energy sublevels;
5. for a full tensor law, the HYPER-001 spectral Hessian conditions, not merely ray convexity;
6. regularity of the metric map and effective Lorentzian signature throughout the selected domain.

Ray convexity is necessary but not sufficient for rank-one convexity, polyconvexity, strong ellipticity, hyperbolic dynamics, or well-posedness of an unselected field realization. Those properties cannot be claimed before the tensor extension, action, and kinetic closure are selected.

Let (E_{\rm av}<\infty) and (S_E=\{e:w(e)\le E_{\rm av}\}). Linear, positive-leading polynomial, and exponential energies have unbounded mathematical domains but bounded (S_E). Barrier families have both a finite mathematical endpoint and a stricter reachable endpoint (\sup S_E<e_*). A regular finite endpoint with finite energy would not be excluded by finite energy alone; it needs the independent state constraint already allowed by HYPER-001. “Finite-energy universe” therefore never by itself creates a mathematical cutoff.

## 4. Architecture compatibility

- **FOUNDATION-001:** all candidates are response functions of the one medium; none changes ontology. FP-6 blocks selecting their free shape/end-point data without derivation.
- **DEFORMATION-001 / STATE-002:** the scalar formulas are ray restrictions of (W(C)=\Phi(I_1,I_2,I_3)), not substitutes for objective tensor (C). A valid native lift must be permutation invariant on the SPD spectral domain.
- **HYPER-001 / ENERGY-PRINCIPLE-001:** every viable lift must obey positivity, stress-free reference, (K>0,\mu>0), spectral tangent gates, and the accepted endpoint class. Barrier laws implement, but do not derive, the finite bound.
- **DURATION-001:** all energies are statewise and history-free. None introduces fundamental time, viscosity, or a relaxation clock.
- **METRIC-001:** constitutive response may feed the still-unselected map (g^{\rm eff}=G[q,C;{\cal D}]). No formula alone ensures emergent geometry; regularity, one-metric universality, and the V11 linear map remain gates.
- **BALANCE-001:** (dW=P_C:dC) supplies the elastic covector only. Fluxes, sources, inertia, and boundary work remain closure slots.
- **V11 operational recovery:** every viable family has the same linear tangent and can pass the weak-field gate conditionally. None by itself derives Einstein dynamics or fixes normalization.

## 5. Unweighted comparison matrix

Legend: **pass** = intrinsically meets the stated scalar/ray criterion under displayed restrictions; **conditional** = needs a tensor/domain/dynamic completion; **fail** = contradicts the criterion as a complete law. The entries are categorical and deliberately unweighted.

| Family | weak field | progressive hardening | stored energy | finite-energy compatibility | recovery | mathematical stability | emergent geometry | V11 recovery | overall consequence |
|---|---|---|---|---|---|---|---|---|---|
| A Linear | pass | fail | pass | pass: reachable bound only | conditional | pass statically | conditional | pass conditionally | tangent limit only |
| B Polynomial | pass | pass with coefficient restrictions | pass with positive leading term | pass: reachable bound only | conditional | conditional globally | conditional | pass conditionally | viable, not selected |
| C Exponential | pass | pass | pass | pass: reachable bound only | conditional | pass on rays; tensor proof pending | conditional | pass conditionally | viable, not selected |
| D Asymptotic barrier | pass | pass | pass | pass: mathematical and reachable bounds distinguished | conditional | pass on rays; tensor proof pending | conditional | pass conditionally | strongest generic finite-bound match |
| E1 Gent-type | pass | pass | pass | pass: mathematical and reachable bounds distinguished | conditional | pass on displayed ray; full-domain proof pending | conditional | pass conditionally | strong standard barrier exemplar |
| E2 Hencky quadratic | pass | fail on log-strain rays | pass | pass: reachable bound only | conditional | conditional in native variable | conditional | pass conditionally | finite-strain coordinate model, not progressive winner |

There is no score total and no weighting. D and Gent-type E1 dominate only if energetic enforcement of the separately accepted finite bound is made a selection requirement. Without that extra choice, B, C, D, and E1 remain admissible competitors. CP-1--CP-7 alone do not distinguish polynomial from exponential versus barrier growth.

## 6. Dependency graph

```text
FOUNDATION-001 + STATE-002
             |
DEFORMATION-001: C in Sym+(3), objective spectrum
             |
HYPER-001 + ENERGY-PRINCIPLE-001
 W=Phi(I1,I2,I3), equilibrium/tangent/domain gates
             |
CONSTITUTIVE-002 (this comparison)
  |-- common ray test: w, sigma=w', kt=w''
  |-- A linear -> weak tangent only
  |-- B polynomial -> viable conditional family
  |-- C exponential -> viable conditional family
  |-- D barrier -> strongest generic finite-bound class
  `-- E standard forms -> Gent exemplar; Hencky partial match
             |
future native-law selection (not performed)
             |
tensor lift + metric/source/action closure
             |
DURATION-001 + METRIC-001 + BALANCE-001 gates
             |
V11 operational validation / governing equations (downstream)
```

## 7. Completion boundary

All requested families have been tested with identical definitions of stress, tangent, compliance, weak response, energy, reachability, recovery, and stability. The comparison supports a progressively hardening hyperelastic class and, conditionally on energetic enforcement of the frozen finite bound, favors an asymptotic barrier family with Gent-type laws as established exemplars. It does not justify a unique constitutive equation. The next milestone must choose or derive the invariant tensor lift and its permitted constitutive data before any governing-equation closure.
