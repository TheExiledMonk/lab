# PBUF NATURE-001 — Natural Analog Constraints for the Spacetime Medium

## Decision

Natural media support a narrow set of **behavioral constraints**, not a material identity or a constitutive equation. The robust intersection is a stable, recoverable, energy-storing continuum with a linear tangent response near equilibrium, a nonlinear completion over a finite admissible domain, and causal wave-capable dynamics. Thermal/background state dependence and a tensorially adequate response must be retained. Nature does **not** decide whether the PBUF constitutive operator is local, gradient-enriched, or genuinely nonlocal, nor does it select a strain variable, stored-energy formula, saturation mechanism, or coefficient.

This audit therefore strengthens the hyperelastic/reversible-elastic classification from CONSTITUTIVE-CLASS-001 while preserving the unresolved choices identified by MEDIUM-001. It introduces no parameter, fit, equation, or change to a frozen validation experiment.

## 1. Scope and method

The analogies are used only as a cross-material behavior survey. A behavior is retained when it recurs across mechanically distinct systems and is independently compatible with established PBUF requirements. A mechanism is rejected when it depends on a particular substrate, geometry, dimensionality, microstructure, boundary, dissipation model, or named material law.

Classification legend:

- **A — Supported:** follows from an authoritative PBUF requirement and is reinforced by the natural behavior.
- **B — Compatible:** can occur without conflict, but is not selected or derived by current PBUF inputs.
- **C — Incompatible:** conflicts with an authoritative PBUF requirement as a native behavior.
- **D — Requires future derivation:** cannot be classified without the missing PBUF kinematics, action, coarse graining, or dynamics.

“Universal” below means universal enough to constrain the shortlisted PBUF elastic class; it does not mean that every material displays the behavior under every condition.

## 2. Behaviour comparison matrix

The detailed matrix is in `behaviour_comparison_matrix.csv`. Each extracted behavior is compared separately with the finite elastic bound, wave-supporting spacetime, stress-energy loading, thermal dependence, saturation, recovery, and cosmological evolution.

| Extracted behavior | Finite bound | Waves | Stress-energy load | Thermal | Saturation | Recovery | Cosmological evolution | Result |
|---|---:|---:|---:|---:|---:|---:|---:|---|
| Wave propagation | A | A | B | B | D | B | B | A |
| Energy storage | A | A | A | B | B | A | A | A |
| Recovery after unloading | A | B | A | B | B | A | B | A |
| Increasing resistance with deformation | B | B | B | B | B | A | B | B |
| Saturation or finite deformation limits | A | B | B | B | A | B | A | A |
| Stability near extreme loading | A | B | B | B | B | A | B | A |
| Energy release after unloading | A | A | A | B | B | A | B | A |
| Local versus nonlocal response | B | A | B | B | B | B | B | D |

Two boundaries are decisive. First, PBUF background saturation is not evidence for a local stress–strain plateau. It constrains the homogeneous history that the future action must reproduce. Second, spatially extended response from a Green function does not establish a genuinely nonlocal constitutive operator.

No surveyed universal behavior is classified C in the table because incompatible material-specific mechanisms have been removed before the PBUF comparison. The rejected list contains the C and unsupported cases.

## 3. Natural analogy audit

### Elastic solids

Elastic solids demonstrate storage, recovery, wave propagation, and—when their kinematics permit it—shear response. Their transferable content is a stable reference state and a positive tangent response. Crystal structure, preferred axes, defects, plastic yield, fracture, and terrestrial coefficients are not transferable.

### Hyperelastic materials

Hyperelasticity demonstrates that finite reversible response can be organized by a state function and that local linear elasticity can arise as the tangent theory. Some hyperelastic materials stiffen or approach limiting extension, but their chain statistics and named energy functions are mechanisms, not universal PBUF facts. This analogy supports the functional architecture only.

### Viscoelastic materials

Viscoelastic systems show causal memory, relaxation, attenuation, and delayed recovery. These behaviors warn that coarse graining can generate an open-system sector, but they require relaxation spectra, entropy production, or dissipative coefficients absent from the accepted ontology. Native viscoelasticity is therefore not selected. It may be reconsidered only if derived later.

### Water and wave-supporting fluids

Water illustrates finite-speed disturbance propagation, interference, reflection, and energy transport. It establishes that observable wave behavior depends on the dynamical completion rather than on static storage alone. Its molecular structure, free surface, viscosity, vorticity, and lack of static shear are analogy-specific. Water is usable for wave behavior only and is incompatible as the native material identity.

### Membranes and cloth

Membranes and cloth illustrate geometric nonlinearity, state-dependent tangent stiffness, wrinkling, anisotropy, and changing wave speeds under tension. The transferable lesson is that stability and incremental response can depend on the deformed state. Sheet dimensionality, weave, slack fibers, bending thickness, and preferred textile directions are rejected.

### Twisted fibres

Twisted fibres illustrate coupled deformation modes, increasing resistance, finite extension, and recoverable energy release. They show that nonlinear coupling can be conservative. Filaments, helicity, friction, knots, chirality, and polymer locking cannot be imported into PBUF.

### Additional natural continua

Crystalline continua reinforce the positive-Hessian and branch-stability tests but do not license a lattice, defect, or preferred-frame ontology. Evolving self-gravitating continua reinforce the distinction between background state evolution and local tangent response but do not supply a PBUF equation of state or phase-transition mechanism.

The row-level audit is in `natural_analogy_audit.csv`.

## 4. Universal mechanical properties shortlist

1. **Stable recoverable reference state.** The unloaded state must be a physical energy minimum with finite positive susceptibility on admissible perturbations.
2. **Conservative stored energy.** Native loading stores energy, and unloading returns it through allowed medium modes or consistent exchange with coupled sectors.
3. **Causal wave-capable dynamics.** The static energy must admit a well-posed kinetic completion; a static stress–strain relation alone is insufficient.
4. **Linear tangent plus nonlinear completion.** Linear elasticity is required near equilibrium, while the native finite-response theory must remain valid beyond that neighborhood.
5. **Finite admissible response.** The physical branch must respect the finite elastic bound through a regular endpoint, bounded domain, or saturation behavior; the mechanism is not selected.
6. **State and thermal dependence.** The future response must inherit the established thermal and cosmological state dependence.
7. **Tensorially adequate loading and response.** The eventual variables and metric map must accept the full variational stress-energy load and recover the effective-GR shear/tensor sector.
8. **Explicit locality class.** Local constitutive response, gradient enrichment, and genuine causal nonlocality must be distinguished and derived rather than inferred from spatially extended solutions.

The machine-readable shortlist is in `universal_mechanical_properties.json`.

Increasing resistance is retained as a compatible design tendency, not a mandatory global monotonicity rule. Natural systems may stiffen, soften, buckle, or change branch. PBUF requires stability on its claimed branch, not imitation of every natural failure mode.

## 5. Properties rejected as analogy-specific

The audit rejects molecular or atomic identity, water-like fluid identity, polymer-chain locking, named textbook hyperelastic functions, viscosity and relaxation spectra, plastic yield, fracture, damage, textile geometry, fibre chirality, crystallographic axes, defects, free surfaces, surface tension, and all analogy-derived coefficients or scales.

This separation matters because several rejected mechanisms can reproduce one desired behavior. For example, finite extension may result from chain locking, geometric constraints, or a restricted state domain. Observing the shared behavior cannot choose among those mechanisms. The complete rationale is in `rejected_analogy_specific_properties.csv`.

## 6. Recommended constraints on the future stored-energy functional

The future PBUF stored-energy functional/action should:

1. be written only after covariant deformation variables and the reference state are established;
2. possess a stable unloaded minimum and positive physical Hessian;
3. recover the linear-elastic tangent regime;
4. remain single-valued and reversible in the native sector;
5. respect the finite elastic bound through a regular admissible endpoint without selecting its mechanism prematurely;
6. avoid unstable softening, singular release, ghosts, and loss of hyperbolicity on the claimed branch;
7. inherit the authoritative thermal and cosmological dependence and reproduce homogeneous saturation without identifying it with a local plateau;
8. obtain the generalized load from the full stress-energy tensor through the eventual variational metric map;
9. admit causal wave dynamics and a controlled effective-GR limit with consistent energy-momentum exchange;
10. declare whether gradients or genuine nonlocal kernels occur and, if so, derive their covariance, positivity, and causal support; and
11. contain no analogy-selected material formula, coefficient, scale, relaxation time, limiting strain, or fitted parameter.

These are constraints, not a proposed functional. They leave the deformation variable, invariants, locality class, functional shape, coefficients, and any derived dissipative sector unresolved. The machine-readable record is in `stored_energy_constraints.json`.

## 7. Implication for the next constitutive milestone

Nature supports the reversible hyperelastic architecture but does not close MEDIUM-001’s missing constitutive principle. The next milestone must derive or postulate, with independent PBUF justification, the covariant deformation variables and admissible invariants before selecting a functional. Candidate functionals must then be tested against the stability, finite-bound, thermal inheritance, causal wave, tensor loading, conservation, homogeneous-history, and effective-GR gates listed here.

Weak-lensing artifacts remain downstream and were neither imported, executed, nor modified.

## Completion statement

All eight requested behaviors have been audited across the accepted analogies and compared with all seven stated PBUF constraints using the A–D classification. Universal behavioral constraints have been separated from substrate-, geometry-, and mechanism-specific properties. The result constrains a future native PBUF stored-energy functional without claiming a known-material identity, deriving an analogy-based equation, introducing parameters, fitting observations, or modifying frozen validation experiments.
