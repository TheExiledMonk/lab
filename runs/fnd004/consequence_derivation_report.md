# PBUF FND-004 — Consequences of the three-dimensional microscopic ontology

## Result

Adopting the premises as axioms yields a small set of exact counting, equal-coupling, and identifiability consequences. It does **not** uniquely yield a continuum field equation, stiffness, propagation law, or weak-lensing signal. Those results become conditional only after the inherited symmetry, locality, stability, scalarization, dynamics, and photon-coupling assumptions are stated. No parameter was fitted and the frozen weak-lensing laboratory was not imported or modified.

## Axioms and explicit auxiliary conditions

- **A1:** Spacetime has three independent microscopic degrees of freedom corresponding to the three spatial dimensions.
- **A2:** Matter couples identically to each microscopic degree through g_dev=1/137.
- **A3:** Macroscopic continuum behaviour emerges by coarse graining the microscopic state.

- **C1:** q transforms as the defining spatial-vector representation and no background direction is present.
- **C2:** The long-wave theory is local, analytic, homogeneous, parity-even and rotationally invariant.
- **C3:** The reference state is stable and nondegenerate, and retained modes separate from microscopic scales.
- **C4:** A covariant scalarization/source map is specified and non-observed modes decouple.
- **C5:** Time dynamics (inertial or dissipative) and kinetic normalization are specified.
- **C6:** Light couples to a specified effective metric or ray-deflection functional.
- **C7:** The matter-coupling normalization is fixed independently of all other source coefficients.

The C-conditions are not new foundation axioms: they are clearly exposed hypotheses delimiting conditional effective-theory statements. A1's phrase “corresponding to spatial dimensions” does not alone specify a vector transformation law, so C1 remains necessary.

## Consequence and axiom traceability matrix

| ID | Prediction/boundary | Class | Minimum axioms | Extra conditions | Relation | Observable meaning |
|---|---|---|---|---|---|---|
| P01 | Exactly three microscopic component labels; generic N has N and a scalar has one. | exact | A1 | none | q=(q1,q2,q3) | Component count is physical only if components can be independently excited/read out. |
| P02 | Equal bare component-source vertices. | exact | A2 | fixed common normalization | g1=g2=g3=g_dev | No component-dependent coupling at the axiom scale. |
| P03 | The equal-coupling vector is the normalized singlet direction and two source-orthogonal combinations are unsourced by a common scalar load. | exact linear algebra | A1,A2 | linear common scalar source | q_parallel=(q1+q2+q3)/sqrt(3); J_perp=0 | Two dark combinations exist at linear source level, but their propagation is not fixed. |
| P04 | The coherent source strength carries a sqrt(3) amplitude (factor 3 in a quadratic response) relative to one equally normalized component. | conditional quantitative | A1,A2 | independent orthonormal components; linear response; same per-component susceptibility | \|g\|=sqrt(3) g_dev; \|g\|^2=3 g_dev^2 | N ontology replaces 3 by N; scalar has factor 1. |
| P05 | A rotationally invariant rank-two constitutive tensor in the vector sector has transverse/longitudinal form. | conditional symmetry | A1 | C1,C2 | Gamma_ab(k)=A(k^2) delta_ab+B(k^2) k_a k_b | One longitudinal and two degenerate transverse eigenmodes. |
| P06 | The leading static vector continuum equation has two gradient stiffnesses, not a uniquely fixed scalar Helmholtz law. | conditional continuum | A1,A3 | C1,C2,C3 | K q-a Laplacian(q)-b grad(div q)=J | G_T=a and G_L=a+b; positivity requires K>0, a>0, a+b>0. |
| P07 | If only the common/longitudinal scalar sector is retained, the CORE-type scalar equation follows in form. | conditional reduction | A1,A2,A3 | C2,C3,C4 | K u-div(G grad u)=s(rho) | The axioms fix neither K, G nor s(rho); the equation is not inevitable from A1-A3 alone. |
| P08 | Isotropy forbids component-dependent masses and directional stiffness at leading order. | conditional symmetry | A1,A2 | C1,C2 | K_ab=K delta_ab; equal transverse polarizations | Splitting signals broken isotropy, unequal coupling, or extra structure. |
| P09 | Propagation has one longitudinal and two transverse branches when inertial vector dynamics is supplied. | conditional propagation | A1,A3 | C1,C2,C3,C5 | omega_L^2=(K+(a+b)k^2)/M; omega_T^2=(K+a k^2)/M (double) | Three fixes transverse degeneracy at two; N internal components generally gives N-1 source-orthogonal modes, not spatial L/T modes unless it is a vector. |
| P10 | g_dev directly fixes the normalized common matter vertex, but A2 postulates rather than derives its numerical value. | corrected identifiability boundary | A2 | a completed calibrated response/readout chain for measurement | g_vec=g_dev(1,1,1); no independent coupling multiplier | Absolute calibrated response can be g_dev-sensitive; normalized component ratios cancel g_dev and cannot determine its magnitude. |
| P11 | No weak-lensing deflection law or amplitude follows from A1-A3. | non-derivation | A1,A2,A3 | C4,C6,C7 absent | none | A specified photon/effective-metric coupling is required before lensing is predicted. |
| P12 | Coarse graining preserves the three-component target space but does not by itself select locality, dynamics, coefficients, or a scalar observable. | exact boundary | A3 | none | Q_L=C_L[q] in R^3 | Kernel moments and scale separation determine finite-resolution corrections. |

## Continuum, constitutive, symmetry and propagation derivation

Under C1–C3, Fourier-space rotational covariance permits only `Gamma_ab(k)=A(k^2)delta_ab+B(k^2)k_a k_b`. Expanding analytically at small `k` and varying the quadratic energy gives `K q-a Laplacian(q)-b grad(div q)=J`. The longitudinal and transverse static stiffnesses are `G_L=a+b` and `G_T=a`; stability requires `K>0`, `G_L>0`, and `G_T>0`. Neither their magnitudes nor equality is fixed by three components or by `1/137`.

With an inertial normalization `M`, the branches are `omega_L^2=(K+G_L k^2)/M` and `omega_T^2=(K+G_T k^2)/M`, with the transverse branch doubled. Overdamped dynamics would instead produce relaxation rates, demonstrating that propagation is conditional on C5. A scalar-only closure has one stiffness and no transverse branch.

A normalized, isotropic coarse-graining operator preserves constants and rotational covariance, with higher kernel moments contributing resolution corrections. A3 alone does not guarantee such a kernel, scale separation, or locality. After C4, a single retained scalar may obey `K u-div(G grad u)=s(rho)`, making the *form* of the existing effective closure natural, but not its coefficients or nonlinear constitutive law.

## Alternative ontologies

| Feature | Three-component PBUF | Generic N | Scalar only | Three-specific finding |
|---|---|---|---|---|
| Microscopic count | 3 (axiom) | N | 1 | Exactly two combinations orthogonal to one common source. |
| Equal-coupling norm | sqrt(3) g_dev | sqrt(N) g_dev | g_dev | Factor 3 in normalized quadratic coherent response. |
| Spatial-vector polarizations | 1 longitudinal + 2 transverse | Not defined unless representation is specified | 1 scalar | Twofold transverse degeneracy, conditional on spatial-vector identity. |
| Isotropic static stiffnesses | K, G_L, G_T | Representation-dependent | K, G | Vector L/T split, not its coefficient values. |
| Scalar continuum closure | Possible after projection/decoupling | Possible after projection/decoupling | Direct | None at leading scalar order. |
| Weak lensing | Undetermined | Undetermined | Undetermined | None without a light-coupling law. |

For a generic equal common source, the internal source direction is `(1,...,1)/sqrt(N)` and there are `N-1` orthogonal combinations. The coherent amplitude is `sqrt(N) g_dev`. These are internal-space facts; spatial longitudinal/transverse language is valid only when the state actually carries the spatial-vector representation.

## Unique PBUF prediction catalogue

| ID | Claim | Status | Required experiment | Dependencies |
|---|---|---|---|---|
| U1 | Two source-orthogonal microscopic combinations for a single equal common source. | exact under linear common sourcing | Independently excite/read component combinations; search for two dark channels. | A1,A2 plus linear common source |
| U2 | Coherent quadratic response multiplicity exactly 3 rather than N or 1. | conditional quantitative | Compare calibrated single-component and coherent susceptibilities without fitting. | A1,A2 plus common normalization and equal susceptibility |
| U3 | One longitudinal and two degenerate transverse modes. | conditional representation prediction | Directional spectroscopy/propagation measurement; test transverse degeneracy. | A1,C1,C2,C3,C5 |
| U4 | No unique scalar-lensing signature follows merely from having three components. | exact negative result | Not directly testable; it is a model-identifiability constraint. | Logical comparison of A1-A3 |

The positive unique predictions are representation/counting signatures, not a weak-lensing curve. In particular, `g_dev=1/137` supplies a bare scale but remains observationally degenerate wherever another coefficient multiplies it. A normalized microscopic action or independent observable is required to expose the number itself.

## Weak-lensing implication

There is no derived weak-lensing amplitude, radial profile, or photon trajectory from A1–A3. At most, if C4 and C6 reduce the ontology to the same scalar continuum interface, the existing laboratory can consume that scalar field. Such compatibility is not a prediction and cannot discriminate three components from a scalar or generic N model when extra modes decouple.

## Next validation milestone

Proceed to **FND-005: covariant mode-and-coupling validation**. Specify a normalized quadratic microscopic action, the exact rotation representation, matter source tensor, and time kinetic term. Its preregistered no-fit tests should be: (1) measure or calculate the one-longitudinal/two-transverse spectrum and degeneracy; (2) test the direct `1/137` vertex with a calibrated response/readout chain; (3) derive the scalarization and photon coupling; and only then (4) run the unchanged weak-lensing laboratory as an out-of-sample consequence. Failure at steps 1–3 should stop lensing claims rather than trigger tuning.

## Completion checks

Every axiom has traced consequences and every prediction names its minimum axioms and auxiliary conditions. All automated checks pass: **True**.
