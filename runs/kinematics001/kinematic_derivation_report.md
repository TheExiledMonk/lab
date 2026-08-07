# PBUF KINEMATICS-001 — Native deformation variables

## Decision

**Outcome D, with the restricted multiplicity of Outcome B/C.** The authoritative corpus does not uniquely identify a native deformation variable. It does, however, exclude a bare scalar as the complete kinematics and restrict viable descriptions to covariant, tensorially adequate relative-deformation frameworks with a stable unloaded configuration.

The exact missing kinematic principle is:

> **Reference-and-identification principle:** specify which covariant medium configuration carries physical clocks and rulers, which configuration is unloaded, and whether the comparison is a material map, a relative coframe, or a relative rank-two tensor; specify its internal gauge quotient and whether the comparison acts on three material directions or all four spacetime directions.

Until this principle is supplied, the strongest native notation is a dimensionless relative-deformation endomorphism `C`, defined only within each admissible realization. Its invariant eigenvalue data—not its basis components—are candidates for a future local isotropic stored-energy law. This report does not define such a law.

## 1. Authoritative boundary

V11 keeps operational Lorentz invariance, Einstein's equations, standard GR, and standard quantum dynamics intact. It supplies the homogeneous functions `alpha_T(a)`, `epsilon_0,T(a)`, and `Omega_sigma(a)`, including an activation/saturation history. It supplies no local medium field, material labels, reference metric/coframe, strain tensor, deformation gradient, gauge quotient, or map from deformation to the physical metric. Consequently `Omega_sigma(a)` and its saturation factor cannot be reinterpreted as a local strain or a limiting stretch.

MEDIUM-001 leaves the tensor type of `chi^A` and `g_eff[chi]` open. CONSTITUTIVE-CLASS-001 favors reversible covariant elasticity with a tensorially adequate small-deformation limit. NATURE-001 requires a recoverable reference state, finite admissible response, and a linear tangent regime, but explicitly does not select a strain variable. GEOMETRY-001 shows that a scalar cannot generate one Lorentzian metric without prior tensorial structure and leaves tensor, coframe, material-coordinate, and scalar-plus-frame realizations inequivalent.

Historical CORE/FND lattice variables and the conditional MB scalar `u` are not promoted to native kinematics: the former are non-authoritative conditional constructions, and the latter is at most a coarse-grained projection.

## 2. Kinematics, not dynamics

Kinematics consists only of:

1. a configuration field `q` and its gauge equivalence;
2. an unloaded configuration `q_0`;
3. a covariant relative comparison `C[q,q_0]` satisfying `C[q_0,q_0]=1`;
4. the admissible domain (orientation, rank, signature, and branch conditions); and
5. scalar invariants of `C` under diffeomorphisms and internal relabellings.

No balance equation, kinetic term, stress, modulus, response kernel, matter coupling coefficient, stored-energy density, or evolution rule is kinematic. Gradients of a strain field may be admitted as higher-kinematic data, but choosing to use them is a locality/derivative-order decision and is not fixed here.

## 3. Admissible frameworks

### 3.1 Material-coordinate (relativistic-solid) class

Let three scalar labels `phi^I(x)` (`I=1,2,3`) identify medium elements, and let `kappa_IJ(phi)` be the unloaded material metric. Given the physical metric retained by V11,

`B^IJ = g^mu nu (partial_mu phi^I)(partial_nu phi^J)`

is a covariant spacetime scalar and a material contravariant tensor. The dimensionless mixed relative deformation can be represented by

`C^I_J = B^IK kappa_KJ`,

with the convention for `kappa` chosen so that `C=1` in the unloaded state. Material relabellings act by similarity, so the eigenvalues or characteristic invariants of `C`, not its coordinate components, are objective. The rank-three condition and positive spatial spectrum select the physical solid branch.

This is the preferred *restricted family* because it gives shear, separates material relabelling from spacetime covariance, and has a standard infinitesimal-strain limit. It is not uniquely authorized: V11 does not supply `phi^I`, `kappa_IJ`, a rest congruence, or an explanation of how a material rest structure coexists with operational Lorentz invariance. It also treats `g` as an input to strain and therefore does not by itself solve GEOMETRY-001's emergence problem.

### 3.2 Relative-coframe class

Let an invertible coframe `E^a_mu` describe the current clock-and-ruler configuration and `Ebar^a_mu` the unloaded one. The relative map

`F^a_b = E^a_mu (Ebar^{-1})^mu_b`

is dimensionless. A deformation endomorphism is obtained from an internal-metric adjoint, schematically `C=F^sharp F`. Local Lorentz transformations common to both coframes and diffeomorphisms must be quotiented; otherwise six frame-orientation variables are falsely counted as strain.

This class can contain the full four-dimensional clock-and-ruler comparison and naturally maps to metric perturbations. Its unresolved points are decisive: V11 does not define a reference coframe, whether boosts are gauge or physical medium modes, or whether a Lorentzian `F^sharp F` has the positivity/spectral properties required for a single real principal branch.

### 3.3 Relative symmetric-tensor class

If a current symmetric Lorentzian tensor `q_mu nu` and unloaded tensor `qbar_mu nu` are independently physical medium configurations, define

`C^mu_nu = qbar^mu alpha q_alpha nu`.

It is diffeomorphism covariant and equals the identity when unloaded. When `q` is simply the already-effective metric, however, this is a parametrization of geometry rather than an independently derived medium deformation, contrary to the caution in MEDIUM-001. If `q` is independent, its relation to the one physical metric is precisely the missing GEOMETRY-001 map. Lorentzian relative tensors can also have non-positive or non-real spectra unless the admissible branch is separately postulated.

### 3.4 Micromorphic extension

Any of the above may be enlarged by an independent internal distortion/director field. Its relative map and invariants would be additional kinematics. This branch remains conditional because no authoritative PBUF director, microrotation, or couple-stress degree of freedom exists.

### 3.5 Scalar and vector fields

A scalar can encode volume/activation or another projection, and a vector can identify a congruence or preferred direction. Neither alone carries generic shear and tensor-wave kinematics or constructs a nondegenerate Lorentzian metric without additional structure. They are admissible only as parts or projections of a mixed system, not as the unique primitive deformation variable.

## 4. Reference state

The reference state is an equivalence class, not a preferred coordinate chart:

`R_0 = { q_0 modulo diffeomorphisms and internal gauge | C[q_0,q_0]=1 }`.

It must satisfy all of the following purely kinematic conditions:

- nondegenerate Lorentzian physical geometry and the standard local Minkowski form in a freely falling frame;
- homogeneous and isotropic representative when used for the V11 cosmological background;
- zero relative strain (`C=1`) rather than zero field components;
- admissible orientation, rank, and signature;
- compatibility with the V11 thermal/cosmological state, so a family `R_0(T,a)` is allowed if V11 microphysics changes the unloaded configuration; and
- no observable coordinate, frame, or material-label dependence.

V11 does not say whether the unloaded reference is Minkowski, the instantaneous FLRW background, a temperature-dependent natural configuration, or a fixed material configuration pulled along cosmological evolution. This is not a small convention: these choices assign background expansion either to reference evolution or to deformation and therefore change every invariant. The numerical values `alpha`, `epsilon_0`, `Omega_sigma`, and `Rmax` must retain their V11 meanings and do not define `C` or `R_0`.

## 5. Minimal invariant catalogue

For an isotropic, parity-even, local rank-three solid branch, let `lambda_A` be the three positive eigenvalues of `C`. A complete algebraically independent set is equivalently

`I1 = tr C`,

`I2 = 1/2[(tr C)^2-tr(C^2)]`,

`I3 = det C`.

At the unloaded state these take `(3,3,1)`. Any smooth symmetric local scalar of `C` can be expressed locally through these three invariants. A useful separation is volume `J=sqrt(I3)` (subject to the precise stretch convention) plus two independent isochoric/shear invariants formed from `Cbar=I3^(-1/3) C`. This is only a change of invariant coordinates.

For a genuine rank-four clock-and-ruler comparison, four characteristic invariants `e_n(C)`, `n=1,...,4`, are required in general. Reducing them to three requires a principle fixing or gauging the temporal eigenvalue. This is a principal unresolved choice.

If orientation or parity violation is physical, pseudoscalar invariants may be needed; current PBUF supplies no such premise, so they are excluded from the minimal catalogue. Curvature, `nabla C`, higher derivatives, and two-point combinations are admissible only in gradient or nonlocal enrichments and are not part of the minimal local set. Invariants involving a medium four-velocity are required only if that velocity is independently physical.

## 6. Finite-bound and weak-field variables

No variable can be proved regular at the PBUF finite elastic bound because V11 does not map that bound to a local principal stretch. Conditional on a finite, nonzero, orientation-preserving endpoint:

- `C`, its elementary invariants, Green-type strain `E=(C-1)/2`, and logarithmic strain `H=(1/2) log C` remain finite;
- inverse/Eulerian measures can diverge if an eigenvalue tends to zero;
- determinant-normalized shear variables become singular if `det C -> 0`; and
- a bounded saturation coordinate would require the missing bound and normalization and therefore cannot be selected now.

Thus the regularity-safe recommendation is to retain the principal spectrum or elementary invariants on an explicitly stated nondegenerate branch, without choosing a singular reciprocal parametrization.

For `C=1+delta C`, the infinitesimal strain

`varepsilon = (1/2) delta C + O(delta C^2)`

is the natural weak-field variable. Its trace is the scalar/volumetric channel and its trace-free part is the shear channel. In the tensor/coframe classes it is linearly related to the metric perturbation after gauge quotient; in the material class it is the symmetrized displacement-gradient form on the unloaded background. This recovers the *kinematic form* needed by a linear effective-GR limit, but Einstein dynamics and normalization remain downstream matching conditions.

## 7. Covariance and symmetry audit

Objectivity requires invariant dependence under spacetime diffeomorphisms and the internal gauge/relabeling group. Isotropy allows only symmetric functions of the principal values; it does not remove shear. Operational Lorentz invariance requires either that local Lorentz transformations are gauge or that any medium rest structure is unobservable in the effective regime. A fixed coordinate tensor, coordinate displacement `x'-x`, or scalar density contrast is not a native covariant strain.

Using a reference structure does not automatically violate covariance, but making it nondynamical or observationally preferred may violate V11's retained physics. Whether the reference is gauge bookkeeping, a physical second structure, or a natural configuration derived from the same medium is part of the missing principle.

## 8. Unresolved choices

1. Three material directions versus four clock-and-ruler directions.
2. Material map, coframe, or independent symmetric tensor as the primitive configuration.
3. Whether the physical metric is input, output, or one component of the medium configuration.
4. Fixed, evolving, or temperature-dependent unloaded state.
5. Internal symmetry/relabeling group and which apparent modes are gauge.
6. Local first-gradient kinematics versus gradient, micromorphic, or causal nonlocal enrichment.
7. The local meaning, if any, of the finite elastic bound.
8. Whether parity-odd or additional internal invariants exist.

Covariance, isotropy, finite response, and the weak-field limit do not resolve these choices. Multiple inequivalent frameworks pass them.

## 9. Recommendation for HYPER-001

HYPER-001 should be a **kinematic closure gate**, not a stored-energy selection exercise. Before any functional is written it should require an authoritative or Quantum-Engine-derived choice of:

1. the primitive configuration class and its dimensions;
2. the reference-state family and whether cosmological/thermal evolution changes it;
3. the three-versus-four-dimensional comparison and internal gauge group;
4. a dimensionless relative map `C` with `C=1` unloaded and a nondegenerate admissible branch;
5. the one-metric identification consistent with GEOMETRY-001;
6. the minimal characteristic invariants for that rank; and
7. the interpretation of the finite bound in those invariants.

Only after those gates close may HYPER-001 catalogue candidate functional dependence and test stability, causality, thermal inheritance, homogeneous V11 matching, and the effective-GR limit. It must not infer a local strain bound from `S(a)` or introduce a response coefficient to normalize the metric map.

## Completion statement

KINEMATICS-001 is complete under the stated criterion by isolating the precise missing kinematic principle. Current inputs select a restricted family of covariant relative-deformation descriptions and a unique rank-dependent invariant construction, but not one primitive field or reference state. The scientifically warranted outcome is D, accompanied by B/C for the surviving mathematical family.
