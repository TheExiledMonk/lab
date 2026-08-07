# PBUF FND-002 — Justification of the microscopic state

## Result

All nine CORE-001 assumptions are documented and assigned exactly one required label. The audit finds one universal leading-order derivation, five optional realizations, one interface-supported statement, and two new foundational postulates. No supplied PBUF equation derives either the number three or `1/137`; their appearance in CORE-001 establishes provenance, not mathematical support.

## Assumption audit matrix

| ID | Assumption | Label | Explicit? | Derivable? | Can weaken? | Can remove? | Observable consequence | What breaks? |
|---|---|---|---|---|---|---|---|---|
| FND-002-A01 | Exactly three microscopic degrees of freedom | NEW POSTULATE | Yes: CORE-001 calls it an existing PBUF working premise, while explicitly saying it is not derived. | No; neither the supplied equations nor symmetry data select dim(q)=3. | Yes: require only a microscopic state with at least one matter-coupled projection. | From the scalar closure, yes; from the stated three-degree PBUF ontology, no. | Extra components matter only if they have distinct sources, modes, or couplings; the frozen laboratory observes one projection. | CORE-001 internal O(3) interpretation changes, but its scalar long-wave equation need not. |
| FND-002-A02 | Characteristic coupling scale g_dev=1/137 | NEW POSTULATE | Yes: CORE-001 identifies 1/137 as an existing characteristic-scale premise, not a derivation. | No. Corrected CORE-001 places g_dev directly in the normalized matter vertex; its numerical value remains a PBUF premise rather than a consequence of the other assumptions. | No within the stated PBUF ontology: replacing g_dev by an effective coefficient would discard the fundamental-coupling premise. | No without removing the stated PBUF matter coupling; the zero-coupling limit removes matter loading. | The normalized microscopic source and conditional coarse source scale directly with g_dev; downstream observability still requires response and access maps. | The former source-normalization degeneracy is withdrawn; lack of a derivation of the numerical value remains. |
| FND-002-A03 | Isotropic lattice/regulator | OPTIONAL | No earlier supplied material requires a lattice; CORE-001 introduced it as a regulator. | No specific regulator is derivable. Isotropic leading-order propagation follows conditionally from spatial rotational symmetry. | Yes: any statistically homogeneous, isotropic local substrate with a valid long-wave limit suffices. | Yes: a continuum microfield or isotropic graph Laplacian gives the same effective operator. | Finite-spacing anisotropy and dispersion can appear beyond the long-wave regime. | Only regulator-level corrections change if the same isotropic continuum limit is retained. |
| FND-002-A04 | Linear matter-state interaction | OPTIONAL | No quantitative matter-state law exists in the supplied pre-CORE material. | Only as the first Taylor term of an analytic interaction about a reference state, assuming a nonzero linear coefficient. | Yes: require a differentiable local source whose leading perturbative term is linear. | Yes, but a purely nonlinear source changes the weak-loading closure. | It fixes superposition and leading response to density; nonlinear terms create amplitude-dependent response. | The exact Helmholtz form and density scaling generally cease to hold outside the linearized regime. |
| FND-002-A05 | Quadratic local recovery energy | DERIVED | CORE-001 defines it; earlier material asks for stiffness/recovery but does not fix the potential. | Yes as the universal leading nonconstant term of a smooth stable energy expanded about a nondegenerate equilibrium. | Yes: assume only smoothness, stability, and positive Hessian at equilibrium. | The exact quadratic model can be removed; some stabilizing recovery is required for finite static susceptibility. | Leading response is linear; higher powers produce nonlinear saturation away from equilibrium. | Without positive local curvature the unloaded state may be unstable or the Helmholtz mass term vanishes. |
| FND-002-A06 | Nearest-neighbour transmission | OPTIONAL | No; it is a CORE-001 discretization choice. | No. A gradient-squared term is the leading local isotropic spatial correction, but many microscopic couplings generate it. | Yes: finite-range or sufficiently decaying, symmetric couplings with a finite second moment. | Yes: graph, nonlocal, or direct continuum operators may share the same small-wavenumber expansion. | Higher-order dispersion and short-scale propagation depend on the coupling stencil. | Long-wave G is renormalized; nonlocal tails or anisotropy can change the continuum operator itself. |
| FND-002-A07 | Overdamped relaxation | OPTIONAL | No supplied theory selects first-order time evolution; CORE-001 chose it for static relaxation. | No. The static Euler-Lagrange equation is independent of whether relaxation is overdamped, inertial, or otherwise equilibrating. | Yes: require dynamics that converge to stationary extrema of the effective energy. | Yes for the static milestone and frozen laboratory. | Transient response, mode spectrum, causality, and damping differ. | Nothing static if the same equilibrium is reached; time-dependent predictions change completely. |
| FND-002-A08 | Choice of coarse-graining kernel | OPTIONAL | Coarse graining is required by CORE-001, but no supplied PBUF source fixes a kernel. | No unique kernel. Normalization and symmetry constraints follow from preserving constants and rotations. | Yes: any normalized, localized operator with vanishing first moment and controlled higher moments. | A distinct kernel can be removed in a direct continuum formulation; a scale-separation map remains conceptually necessary for discrete models. | Kernel moments control smoothing and finite-resolution corrections, not the leading constant/slow-field limit. | A nonnormalized or biased kernel fails constant preservation; anisotropic kernels imprint directional artifacts. |
| FND-002-A09 | Scalar continuum field u(x) | SUPPORTED | Yes at the existing continuum/frozen-laboratory interface; earlier discovery material still lists scalar versus tensor character as open microscopically. | No fundamental scalar ontology follows. A scalar effective field follows conditionally when the laboratory couples to one invariant projection and other modes decouple. | Yes: u may be the scalar observable sector of vector or tensor internal states. | Not without changing the frozen laboratory interface; it can be emergent rather than fundamental. | Additional unsuppressed tensor/vector sectors would generate responses absent from a single scalar field. | The current constitutive and gradient interfaces require revision if extra continuum modes couple appreciably. |

`DERIVED` means derivable only under the stated regularity and stability conditions, not derivable as an exact microscopic polynomial. `SUPPORTED` means required by the frozen interface and compatible with supplied PBUF material; it is not proof of a fundamental ontology. `OPTIONAL` marks a replaceable realization. `NEW POSTULATE` marks a claimed foundational fact that the supplied theory does not derive.

## Derivation matrix

| ID | Premises | Expansion/reduction | Consequence | Boundary |
|---|---|---|---|---|
| FND-002-D01 | Stable smooth local energy at q=0 | V(q)=V(0)+q^T H q/2+O(\|q\|^3), H>0 | Quadratic recovery at leading order | Fails at a degenerate/critical or nonsmooth equilibrium |
| FND-002-D02 | Analytic local matter interaction | I(eta,q)=I(0,0)-eta h.q+O(eta\|q\|^2,eta^2) | Linear loading at leading order | A symmetry may force h=0 |
| FND-002-D03 | Local isotropic spatial response | Gamma(k)=K+G\|k\|^2+O(\|k\|^4) | K u-G Laplacian(u)=s at long wavelength | Long-range kernels can be nonanalytic; anisotropy makes G a tensor |
| FND-002-D04 | Normalized centered localized kernel | C[u]=u+(mu_2/2d)Laplacian(u)+higher moments | Kernel-independent leading slow field | No scale separation or divergent moments |
| FND-002-D05 | One sourced light mode; others massive/decoupled | u=e.q; integrate out transverse modes | Scalar effective sector | Multiple sourced light modes require vector/tensor closure |
| FND-002-D06 | Any equilibrating dynamics for the same F | stationary state satisfies delta F/delta u=0 | Static equation independent of overdamped choice | Driven/non-equilibrium states need a dynamical postulate |

These reductions define a universality class. They derive the *leading continuum form*, not the exact microscopic substrate or coefficient values.

## Alternative models

| Model | Construction | Same continuum? | Assumption cost | Finding |
|---|---|---|---|---|
| Regular lattice | site state + symmetric finite differences | Yes | Adds regulator geometry and stencil choices | Representative, not preferred foundationally |
| Irregular graph/network | node state + weighted graph Laplacian | Yes, if graph homogenizes isotropically | Adds graph ensemble/weights | Equivalent universality class |
| Continuum microfield | local effective energy with gradient expansion | Yes, directly | Avoids lattice and kernel postulates; assumes locality/cutoff implicitly | Preferred minimal effective formulation |
| Tensor/vector internal state | one scalar sourced projection plus decoupled modes | Yes conditionally | Adds mode masses and couplings | No benefit until extra observables require it |
| Alternative coarse graining | normalized centered convolution, block average, spectral low-pass | Yes at leading order | Finite-scale corrections differ | Operator is not uniquely selectable here |

## Minimal postulate set and revised model

1. P1: A microscopic PBUF state exists and admits a stable, approximately local, rotationally invariant long-wavelength sector.
2. P2: Matter sources at least one dimensionless scalar observable sector u; other microscopic components are absent, massive, or decoupled at laboratory scales.
3. P3: The effective static response has positive finite local susceptibility and a positive leading spatial stiffness.
4. P4: PBUF's exactly-three-state and characteristic-1/137 statements are retained as ontology premises, but existing material supplies no derivation and the scalar closure cannot identify them separately.

The preferred revised model is therefore effective and regulator-independent:

`F_eff[u] = integral [K u^2/2 + G |grad u|^2/2 - s(rho)u + O(u^3, u grad^2 u, grad^4)] dx`, with `K>0`, `G>0`.

Its stationary leading-order equation is `K u-div(G grad u)=s(rho)`. A lattice, nearest-neighbour stencil, Gaussian kernel, and overdamped dynamics are examples that realize this equation; none is foundational. A multi-component microscopic state is permitted, but only its sourced scalar projection belongs in the current continuum interface. This revision changes no frozen weak-lensing code or parameter.

## Irreducible postulates

The irreducible content is P1–P3 for the effective theory. P4 is additionally irreducible only if the distinctive PBUF ontology—exactly three microscopic degrees and the `1/137` association—is demanded. In corrected CORE-001, `g_dev` directly normalizes matter loading and is no longer rescaling-degenerate with an auxiliary coupling. Its numerical value is still a premise rather than a result of this assumption audit, and the two transverse components have no independent observable in the frozen laboratory.

## Recommendation for FND-003

FND-003 should be a **symmetry and identifiability derivation**, not another regulator construction. It should specify the symmetry group and representation of the three microscopic degrees and define how stress-energy (not only static density) transforms and couples. Its decision tests should be: (1) derive three and `1/137` from those structures, or formally retain them as axioms; (2) determine whether extra modes decouple; and (3) derive coefficient signs and scaling without lensing fits. Dynamic, causal evolution should be deferred until the static representation and coupling are fixed.

## Completion checks

All checks pass: **{validation['all_checks_pass']}**. Label counts: `{validation['label_counts']}`. Scope remained theory-only.
