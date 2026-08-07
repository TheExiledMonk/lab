# PBUF CONSTITUTIVE-CLASS-001 — Native constitutive-class identification

## Decision

**Recommendation B: several realizations remain viable, inside a substantially narrowed family.** The established PBUF ontology is most compatible with a **nondissipative, action-based, recoverable elastic medium**, generically described at finite response by the hyperelastic umbrella and at small deformation by its linear-elastic tangent limit. The most natural covariant realization is a relativistic elastic solid/effective field theory of media.

PBUF does not yet discriminate among a local first-gradient realization, a stable gradient enrichment, and a causal nonlocal enrichment. A micromorphic/Cosserat branch remains conditional rather than excluded because PBUF has discussed multicomponent microscopic states, but no authoritative director, microrotation, couple stress, or extra continuum mode has been established. Thus no single classical textbook family is uniquely selected.

This is a classification result, not a constitutive law. No response functional, coefficient, kernel, fractional order, relaxation time, coupling, or weak-lensing input is introduced.

## 1. Authoritative ontology used

The selection inputs are limited to the supplied theory record:

- V11 section 1.1 preserves operational Lorentz invariance, standard GR, Einstein's equations, and standard quantum dynamics. It calls the spacetime picture a constitutive interpretation.
- V11 sections 2.2–2.3.3 introduce a homogeneous elastic stress/density sector, temperature-dependent rigidity inputs, and a scale-factor activation/saturation history. These background quantities are not a local stress–strain law.
- V11 section 2.4 says gravitational and electromagnetic signals propagate as modes of the same medium and imposes present-epoch speed agreement. It does not supply the medium kinetic operator or polarization content.
- MB-001 supplies only a conditional static scalar balance. Its local recovery and gradient terms are not an authoritative native law, and its scalar field cannot establish the tensorial shear sector.
- MEDIUM-001 requires a stable, causal covariant action/effective action, variational stress-energy loading, a retarded response, total Noether consistency, and a controlled effective-GR limit. It explicitly leaves the field type, locality, derivative order, memory, and genuine nonlocality unresolved.
- The later ontology work permits microscopic multicomponent structure and coarse graining but does not establish that the continuum field is a spatial vector, microrotation, or fundamental scalar. Historical CORE/FND regulator choices and the retired direct vertex are not promoted to V11 facts.

Two distinctions prevent false inferences. First, V11's background saturation factor is a cosmological history, not evidence for a particular local nonlinear stress–strain curve. Second, a local differential elastic operator can have a spatially extended Green-function response; that is not by itself genuine constitutive nonlocality.

## 2. Constitutive-family comparison matrix

Legend: **strong** means native compatibility; **conditional** means compatible only after missing PBUF structure is supplied; **limit** means useful only as an approximation; **poor** means it adds or conflicts with native behavior. The detailed cell-by-cell findings are in `constitutive_family_comparison_matrix.csv`.

| Family | Waves | Storage / recovery | Shear | Long range | Saturation / thermal | Stability / causality | GR-limit prospect | Classification |
|---|---|---|---|---|---|---|---|---|
| Linear elasticity | strong conditionally | strong near equilibrium | strong if tensorial | propagation only | limit / compatible | conditional | conditional | tangent limit only |
| Hyperelastic | strong conditionally | strong | strong if tensorial | enrichment needed | strong capability / compatible | conditional | conditional | core umbrella |
| Viscoelastic | attenuating | mixed/dissipative | timescale-dependent | temporal memory | capable / strong | demanding | possible only after decoupling | not native |
| Nonlocal elasticity | dispersive | strong if positive kernel | field-dependent | intrinsic | nonlinear extension / compatible | demanding retarded kernel | conditional | viable enrichment |
| Micromorphic/Cosserat | multiple branches | strong | strong | internal-length response | nonlinear extension / compatible | extra-mode burden | conditional | conditional branch |
| Gradient elasticity | dispersive | strong if positive | field-dependent | quasilocal | nonlinear extension / compatible | derivative/ghost burden | conditional | viable enrichment |
| Fractional models | anomalous | model-dependent | model-dependent | intrinsic scale-free memory | not implied / compatible | demanding | unsupported | eliminate now |
| Ordinary fluid | sound only | bulk only | no static shear | hydrodynamic | equation-of-state dependent | relativistic form possible | can source GR but wrong identity | eliminate as identity |
| Viscoelastic fluid | attenuating | mixed | transient only | memory | capable / strong | causal-relaxation burden | unsupported | eliminate as native |
| Relativistic elastic solid / media EFT | strong conditionally | strong | strong | extensible | strong capability / compatible | designed for covariant tests | best structural prospect | preferred realization |

Families in this table are not mutually exclusive at the same taxonomic level: hyperelasticity describes reversible finite-response behavior, while gradient and nonlocal elasticity describe spatial organization, and micromorphic/Cosserat theory describes additional kinematics. The shortlist is therefore hierarchical rather than a forced one-label choice.

## 3. Property-by-property compatibility

| Property | What PBUF requires or leaves open | Classification consequence |
|---|---|---|
| Wave propagation | Medium modes are required, but inertia, damping, and polarizations are absent. | Keep elastic families with a possible causal kinetic completion; do not infer viscoelasticity from propagation alone. |
| Energy storage | The native sector is described as elastic stress/rigidity and MEDIUM-001 calls for an action. | Prefer conservative stored-energy families. |
| Recovery | A stable unloaded state and finite susceptibility are required at least locally. | Prefer recoverable elasticity; exclude plastic/damage behavior as native. |
| Shear support | A full solid-like and tensor-wave completion needs distortional response, although V11 gives no modulus. | Exclude ordinary fluid identity and reject scalar-only ontology as established. |
| Long-range response | Local propagation, gradient response, and true nonlocal kernels are all open. | Retain local, gradient, and causal-nonlocal branches. |
| Saturation | Only homogeneous activation/saturation is authoritative; local nonlinear saturation is unknown. | Linear elasticity is a limit, while hyperelasticity supplies the appropriate nonlinear umbrella without selecting a function. |
| Thermal dependence | Thermal rigidity/background dependence is authoritative; a local thermodynamic law is not. | Any complete candidate must inherit, not redefine, the V11 thermal pipeline. This does not select one family. |
| Stability and causality | Positive static response, well-posed causal dynamics, retarded support, and operational Lorentz invariance are mandatory. | Exclude unstable, ghostly, instantaneous spacetime-nonlocal, or naive diffusive fundamental laws. |
| Ontology | Elastic physical medium is authoritative; field type, reference structure, locality, and memory are not. | Hyperelastic behavior fits; micromorphic and nonlocal structure remain conditional. |
| Effective GR limit | One physical Lorentzian metric, universal matter coupling, Einstein normalization, Bianchi/Noether consistency, and controlled extra modes are required. | Only covariant realizations qualify; no ordinary continuum label passes automatically. |

The full evidence and compatibility mapping is in `property_compatibility.csv`.

## 4. Elimination rationale

**Ordinary fluid models are incompatible as the native physical identity.** They support pressure and longitudinal sound but not static shear or recovery of shape. A homogeneous PBUF elastic density can be represented as an effective cosmological stress-energy component without turning the underlying medium into water or any other fluid. Water remains an analogy only.

**Viscoelastic and viscoelastic-fluid models are not supported as native.** They require constitutive memory, relaxation spectra, viscosity or entropy production. None is established in V11, MB-001, or MEDIUM-001. Such effects could later emerge as causal open-system corrections, but selecting them now would add ontology.

**Fractional models are presently inadmissible.** Fractional space/time response can summarize scale-free correlations or broad memory spectra, but PBUF derives neither a fractional operator nor an order. This family should be reconsidered only if coarse graining produces it without a fitted exponent.

**Globally linear elasticity is incomplete rather than false.** It is expected as the tangent theory about a stable unloaded state, but it cannot by itself classify finite deformation or local nonlinear response. It remains a required limit of a viable smooth elastic completion, not the native complete law.

**Plastic, damage, and active-media families are excluded as additional families.** Permanent deformation, damage history, yield surfaces, or internally injected power conflict with or exceed the established recoverable, stable elastic ontology.

**Micromorphic/Cosserat theory is not eliminated, but it is demoted to a conditional branch.** Multicomponent microscopic ontology alone does not establish continuum directors, independent rotations, or couple stress. These models also add modes that must be hidden or matched in the GR regime.

## 5. Shortlist

The admissible search space is:

1. **Core behavioral family:** covariant hyperelastic/reversible elastic medium with a stable unloaded state and a stored-action description.
2. **Preferred structural realization:** relativistic local elastic solid or effective field theory of media, provided it yields one universal Lorentzian metric and the Einstein limit.
3. **Viable spatial enrichments:** stable gradient elasticity or causal covariant nonlocal elasticity. Current inputs do not choose between them or prove either is necessary.
4. **Conditional kinematic enrichment:** micromorphic/Cosserat hyperelasticity, only if the missing definition of `chi` establishes independent microstructure and its additional modes pass GR matching.
5. **Required approximation:** linear elasticity as the small-deformation tangent of any smooth shortlisted model.

The shortlist is encoded in `shortlist.json`. It narrows the search without choosing a stored-energy function, metric map, kernel, or coefficients.

## 6. Why the recommendation is B rather than A or C

It is not A because “hyperelastic,” “gradient/nonlocal,” and “micromorphic” answer different classification questions, and PBUF does not fix the kinematics or spatial-response class needed to select one complete combination. Even the preferred relativistic-solid realization is not derivable from V11.

It is not C in the unrestricted sense because existing information does discriminate: it favors reversible stored-energy elasticity, requires a linear tangent regime, rules out an ordinary fluid identity, and supplies no basis for native viscous, fractional, plastic, damage, or active response. The evidence therefore reduces the search space materially while leaving several elastic realizations viable.

## 7. Gate for the next milestone

The next constitutive milestone should choose the **field/kinematic class before choosing a functional form**. It must establish from authoritative PBUF input whether the medium variable is tensor/strain, coframe/material-coordinate, scalar-plus-frame, or genuine microstructure; whether response is local, gradient, or intrinsically nonlocal; and whether dynamics is conservative or requires a derived dissipative sector. Only after that classification is fixed may a native action be derived and tested for thermal inheritance, stability, causality, shear/tensor modes, one-metric matter coupling, and the effective GR limit.

Weak lensing remains downstream validation and played no role in this selection.

## Completion statement

Every required candidate and all ten criteria have been evaluated. The native behavioral class is narrowed to covariant recoverable elasticity/hyperelasticity; linear elasticity is its limiting theory; gradient and causal-nonlocal forms remain viable enrichments; micromorphic/Cosserat structure is conditional; and fluid, native viscoelastic, fractional, plastic, damage, and active identities are excluded on current evidence. The milestone is complete with Recommendation B and no premature constitutive equation.
