# PBUF INERTIA-001 — Physical Origin of Inertia in an Elastic Spacetime Medium

## 0. Verdict

Inertia in PBUF must represent the **momentum response of the complete medium
state to physically calibrated change**.  It is not an elastic restoring force,
not another part of the stored-energy law, and not a synonym for wave
propagation.

The frozen framework has two distinct results.

1. If “wave-supporting spacetime” is a required physical property, then some
   nondegenerate, energy-storing dynamical response is necessary.  A purely
   static or first-order relaxational medium does not support reversible waves.
2. Neither continuity, the frozen elastic energy, nor the bare existence of
   waves determines that response's normalization, tensor structure, temporal
   order, locality, or state dependence.  Those data remain a kinetic closure.

Consequently, the existence of an inertial sector is conditionally necessary,
but its constitutive geometry is not derived.  The slot \(\mathcal K_\tau\) is
more constrained than GOVERNING-EQUATION-001 stated only after the intended
meaning of *wave-supporting*, FP-3, and FP-5 are imposed together.  It cannot
yet be replaced by a unique formula or coefficient.

This is not a kinetic closure.  No final evolution equation is selected here.

## 1. Definitions and logical separation

For this audit:

- **elasticity** is the state-dependent storage law and its restoring
  covector, already fixed in PBUF by \(W(C)\) and \(D_qW\);
- **inertia** is the part of the dynamical law that maps physical history
  change into momentum change and allows energy to be stored in motion;
- **kinetic structure** is the geometry or functional from which momentum and
  inertial response may be obtained;
- **mass density** is one special scalar coefficient in a local Newtonian
  realization, not the definition of inertia;
- **wave support** means a well-posed propagating initial-value sector with
  nonzero real characteristic speeds, rather than diffusion, quasistatic
  adjustment, or a merely spatial sinusoidal equilibrium;
- **emergence** means derivability from previously frozen premises, not merely
  reinterpretation of an additional coefficient.

The physical question “why does the medium resist acceleration?” has two
levels.  At the structural level, resistance means that motion carries
momentum and kinetic energy, so changing it requires exchange with another part
of the closed system or a source.  At the normalization level, one must specify
how much momentum corresponds to a given calibrated rate.  Standard continuum
theory and frozen PBUF can motivate the former, but neither derives the latter
from an elastic potential alone.

## 2. Established-continuum survey: assumptions versus derivations

### 2.1 Classical elasticity

Classical elastodynamics combines three logically separate ingredients:

1. kinematics supplies deformation and acceleration;
2. balance of linear momentum supplies the material momentum rate (usually
   mass density times acceleration); and
3. an elastic constitutive law supplies stress from strain or stored energy.

Mass density and momentum balance are assumed/measured independently of the
elastic moduli.  Substitution of Hooke's law into momentum balance derives the
elastic wave equation and speeds proportional to square roots of
stiffness-to-density ratios.  It does **not** derive density from stiffness.
The same elastic energy consistently admits static equilibrium, arbitrary
positive density, or a quasistatic approximation in which inertia is omitted.

This separation is explicit in standard continuum presentations: balance laws
apply independently of material type, while constitutive equations distinguish
materials.  MIT's elasticity notes likewise display density/acceleration and
elastic moduli as separate terms, and seismic speeds depend on both modulus and
density ([NPTEL balance laws](https://archive.nptel.ac.in/content/storage2/courses/105106049/lecnotes/mainch5.html),
[MIT elasticity notes](https://ocw.mit.edu/courses/20-410j-molecular-cellular-and-tissue-biomechanics-be-410j-spring-2003/262e7bbdd3cacf750d2795b7519305c4_rev_ses_mid.pdf),
[MIT geophysics notes](https://ocw.mit.edu/courses/12-201-essentials-of-geophysics-fall-2004/ed326cbc7ba3cbbfb0c21fb6a1dba00a_ch1.pdf)).

**Origin classification:** independent balance/kinetic property at continuum
level; microscopically reducible to constituent momentum only in a richer
theory.

### 2.2 Fluids

For fluids, mass conservation and momentum balance supply inertia through the
material rate of momentum.  Pressure comes from an equation of state or
internal energy; viscous stress comes from a rate constitutive law.  Acoustic
waves arise when compressive storage and inertia are both retained.  Viscosity
dissipates them.  In the creeping-flow limit the same fluid constitutive law
can be used after neglecting inertia, producing an elliptic instantaneous
boundary-value problem rather than a wave problem.

Thus compressibility does not imply inertia, and inertia does not imply
compressibility.  MIT's continuum-electromechanics text treats conservation of
mass, momentum, compressibility, kinetic-energy storage, and viscous
dissipation as separate ingredients ([MIT Continuum Electromechanics,
chapter 7](https://ocw.mit.edu/courses/res-6-001-continuum-electromechanics-spring-2009/pages/open-textbook/)).

**Origin classification:** mass/momentum balance; independent of the pressure
and viscosity closures.

### 2.3 Viscoelastic and hereditary media

Viscoelasticity modifies the stress response by adding rates, internal
variables, or memory kernels.  A standard dynamic model still retains a
separate momentum balance and mass density.  Removing the inertial term yields
viscoelastic relaxation or creep, not reversible propagation.  Eliminating
unresolved internal oscillators can make the observed response frequency
dependent and operator-valued; it can place apparent inertia into a memory
kernel or apparent elasticity into a dynamic modulus.  That is redistribution
under coarse graining, not creation from static stored energy.

The separation is empirically nontrivial: spring-and-mass networks can possess
frequency-dependent effective elasticity while having zero effective mass
density, demonstrating that dynamic elasticity does not uniquely encode
inertia ([Milton and Seppecher
2011](https://arxiv.org/abs/1105.0941)).  Standard memory-kernel treatments also
retain density while the kernel interpolates viscous and elastic response
([MIT nonequilibrium statistical mechanics
notes](https://ocw.mit.edu/courses/5-72-non-equilibrium-statistical-mechanics-spring-2012/902154dd2d50becdc0878ee12052ef87_MIT5_72S12_master4.pdf)).

**Origin classification:** independent at a resolved Markovian level;
potentially operator-valued and history dependent after hidden kinetic degrees
of freedom are eliminated.

### 2.4 Electromagnetic continua

The electromagnetic field is not an elastic material with a separately
measured Newtonian mass density.  Its action, field equations, spacetime
symmetry, and stress-energy tensor jointly assign energy and momentum to field
configurations.  Radiation therefore has inertial response: changing field
momentum requires force, and confined field energy contributes to a composite
system's invariant mass.  In material media the division of total momentum
between field and matter is representation dependent; the conserved total is
the safe object.

This is an example of inertia emerging from **distributed field energy plus
relativistic spacetime symmetry and conservation**, not from energy density
alone.  The stress-energy construction embodies energy and momentum
conservation, while electromagnetic energy in matter can contribute to the
matter sector's inertia ([Crenshaw
2013](https://arxiv.org/abs/1302.6492), [Medina and Stephany
2017](https://arxiv.org/abs/1703.02109)).

**Origin classification:** emergent property of the complete relativistic
field action/stress-energy representation; subsystem attribution is not
unique.

### 2.5 Classical and relativistic field theories

In Lagrangian field theory, spatial-gradient and potential terms determine
restoring response while time-derivative terms determine canonical momentum
and kinetic response.  Both belong to the action, but one is not generally
derivable from the other.  Internal symmetry may constrain their tensor form.
Lorentz invariance can tie temporal and spatial derivative coefficients into a
single spacetime contraction and therefore fixes their *relative* form once a
field normalization and effective metric are chosen.  It still does not make
an arbitrary static potential determine the kinetic metric.

Noether's theorem derives conserved currents from symmetries of the **complete
action**, not from the potential alone.  Relativity further identifies total
rest energy with inertial mass for a closed system
([Einstein Online](https://www.einstein-online.info/en/emc/)).

**Origin classification:** an independent sector of a general action;
constrained, and sometimes partly unified with gradient energy, by spacetime
symmetry.

## 3. Independence theorem

### Theorem I — A rate-free elastic energy does not determine inertia

Let \(\mathcal Q\) be a nontrivial configuration space and let
\(W:\mathcal Q\rightarrow\mathbb R\cup\{+\infty\}\) be any fixed rate-free
stored energy.  Then \(W\) alone determines no unique inertial response.

**Proof.** The same \((\mathcal Q,W)\) is compatible with each of the
following mutually inequivalent histories:

- static stationarity, with no temporal operator;
- first-order gradient flow under any positive mobility;
- reversible second-order dynamics under any positive kinetic metric;
- constrained, nonlocal, or memory-bearing dynamics with the same static
  extrema.

Even within the local reversible second-order subclass, multiplying the
kinetic metric by any positive factor leaves \(W\), its equilibria, and its
elastic tangent unchanged while changing accelerations and every wave speed.
Therefore no functional of \(W\) alone can select the inertial map or its
normalization.  ∎

### Corollaries

1. Positive elastic moduli do not imply positive mass density.
2. A wave-speed observation determines at most a ratio/composite symbol unless
   one sector is independently known.
3. Calling stored energy “distributed” does not close the gap; one must still
   specify how boosts or history rates convert that energy into momentum.
4. Deformation history can produce memory, damping, or dynamic effective mass,
   but only after a history law or hidden variables with their own dynamics
   are added.

### Theorem II — Nontrivial reversible waves require kinetic storage

Suppose a linearized, closed PBUF branch supports nonzero-frequency,
nondissipative propagating modes and admits a positive conserved quadratic
energy.  Then its phase space must contain a nonzero rate/momentum sector, and
the quadratic conserved energy must contain a nonnegative kinetic part on each
physical propagating mode.

**Proof sketch.** A rate-free elastic quadratic form alone is a function of
configuration.  Along a nonstatic periodic orbit it alternately decreases and
increases.  Conservation therefore requires a complementary rate-dependent
storage term.  Positivity on physical modes prevents exponentially growing
negative-energy directions.  A purely first-order positive gradient flow
decreases \(W\) monotonically and cannot realize such an orbit.  ∎

This theorem establishes a necessary kinetic sector, not a particular
acceleration operator.  Gyroscopic first-order Hamiltonian formulations,
canonical pairs, and second-order configuration equations can represent the
same inertial physics.

## 4. Frozen-PBUF physical-origin analysis

### 4.1 What the ontology actually forces

| Frozen input | Consequence for inertia | What it does not supply |
|---|---|---|
| FP-1: one continuous medium | momentum/inertia, if present, belongs to the one medium or derived complete-state variables; no external ether/substrate may carry it | kinetic metric or density |
| FP-3: emergent time | native history dynamics must be invariant under monotone relabeling; acceleration with respect to the arbitrary label is not physical | clock calibration or ordinary quadratic kinetic energy before gauge fixing |
| FP-4: one complete configuration per state | inertia must act on the tangent/cotangent data of the complete physical state and respect gauge equivalence | proof that configuration alone is complete Cauchy data |
| FP-5: V11 limit | effective excitations must share the one Lorentzian causal structure; closed-system energy contributes to effective inertial mass | microscopic medium-to-metric map or native normalization |
| FP-6: no new free constants | a future normalization must be derived, inherited from already authorized scales, or observationally required | permission to guess a density/constant |
| frozen \(W(C)\) | restoring covector and local elastic tangent | momentum map, kinetic storage, temporal order |
| intended wave support | nondegenerate reversible kinetic storage is necessary if “wave” has its strict propagating meaning | scalar, tensor, local, or unique inertia |

Continuity is purely kinematic/topological.  Many continuous media are modeled
quasistatically or diffusively.  Therefore “continuous spacetime medium” alone
does not imply inertia.  “Continuous **wave-supporting** spacetime” implies
inertial response only if wave support is itself an accepted physical
requirement and waves are reversible propagating modes.  It cannot be used
circularly as a derivation from FP-1.

### 4.2 What physically resists acceleration

The strongest PBUF-safe interpretation is:

> The complete medium possesses physical momentum associated with calibrated
> change.  Acceleration changes that momentum.  In a closed system the change
> can occur only through stress-mediated redistribution among regions and
> sectors; an external acceleration requires momentum flux or source.  The
> inertial operator is the cotangent-valued accounting map for that momentum
> change.

This explains resistance structurally through conservation and kinetic
storage.  It does not yet explain the numerical magnitude of the momentum.
Deriving that magnitude requires a kinetic geometry, a complete
energy-momentum construction, or a microscopic coarse-graining rule.

### 4.3 Stored elastic energy and distributed energy

The frozen \(W\) is energy per reference volume associated with deformation.
Three distinct claims must not be conflated:

1. **Restoring claim:** \(D_qW\) drives return toward an unloaded state. This
   is derived.
2. **Composite relativistic claim:** a localized closed excitation's total
   energy contributes to its effective inertial mass. This is required in the
   V11 relativistic regime.
3. **Native kinetic claim:** the local value or Hessian of \(W\) determines
   the medium's inertial operator. This does not follow.

The relativistic composite claim concerns the total conserved stress-energy,
including stresses and kinetic/field sectors.  Applying \(E/c^2\) pointwise to
the elastic density would be an unauthorized local closure and can be wrong
for stressed open subsystems.  At most, FP-5 supplies a downstream matching
condition for isolated complete excitations after the effective metric and
stress-energy map are constructed.

### 4.4 Deformation history

The frozen baseline explicitly contains no memory constitutive variable or
hereditary law.  History therefore cannot currently generate inertia.  A
future elimination of unresolved medium modes could yield a causal convolution
or frequency-dependent effective inertia, but that derivation would presuppose
kinetic dynamics in the unresolved sector.  Memory can relocate inertia in a
reduced description; it cannot originate it from a static energy.

## 5. Emergent time and the form of inertia

Emergent time changes the admissible representation before and after clock
calibration.

### 5.1 Before duration calibration

The arbitrary order label \(s\) has no invariant rate.  Consequently neither
\(dq/ds\), \(d^2q/ds^2\), nor a coefficient multiplying them is independently
physical.  DYNAMICS-001 correctly requires a degree-one native history action.
Its Legendre map is radially degenerate and produces a constraint.  Native
inertia is therefore most generally a constrained cotangent response on the
unparameterized history, not a bare “mass times \(s\)-acceleration.”

### 5.2 After a physical duration gauge is selected

Once a monotonic relational duration \(\tau\) is calibrated, an ordinary
kinetic representation may be legal.  Reparametrization covariance still
requires that descriptions using another label transform to the same physical
history.  Thus \(\mathcal K_\tau\) must be built from physical duration rates,
or be the gauge-fixed image of the native constrained action.

Emergent time therefore forces:

- a native reparametrization identity/constraint;
- no physical dependence on the arbitrary order label;
- clock-gauge covariance;
- a static zero-rate limit; and
- compatibility of characteristics with the effective metric cone.

It does **not** force second temporal order.  A second-order configuration
form, a first-order Hamiltonian form on an enlarged complete state, or a causal
memory operator can encode inertia.  Temporal order is representation- and
state-choice dependent.

## 6. Structural theorem for any future PBUF inertial operator

### Theorem III — Minimal admissible inertial structure

Let \(\mathcal K_\tau\) close GE-001 on the physical state quotient and suppose
the selected branch is intended to be closed, reversible, stable, causal, and
wave supporting.  Then the following properties are necessary.

1. **Typed duality.** \(\mathcal K_\tau[q]\) is a covector (or covector-valued
   distribution) pairing with admissible variations.  A scalar “mass” is only
   a special representation.
2. **Gauge basicness and objectivity.** It descends to the physical quotient,
   annihilates pure gauge directions as required, and transforms covariantly
   under material relabeling and accepted frame actions.
3. **Positive kinetic storage.** The linearized kinetic form is nonnegative and
   nondegenerate on physical propagating modes, modulo genuine gauge or
   constraint null directions.
4. **Symmetric conservative principal part.** For a reversible variational
   branch, the velocity Hessian/Legendre derivative is symmetric (formally
   self-adjoint under the physical pairing).  Antisymmetric gyroscopic pieces
   may occur but cannot alone supply positive kinetic storage.
5. **Causality.** Its retarded response uses no future state; if nonlocal in
   history, its kernel has causal support.  Together with elasticity its
   characteristics lie on or within the single effective V11 cone.
6. **Conservation compatibility.** In the source-free closed system it admits
   a momentum balance and, with clock-translation symmetry, a total-energy
   balance.  Constraint reactions do no work along admissible directions.
7. **Reparametrization/clock covariance.** Before gauge fixing it respects the
   degree-one native action identity; after gauge fixing it is independent of
   the arbitrary label and transforms consistently between clock gauges.
8. **Static consistency.** It contributes no inertial force on a constant
   physical history, apart from separately identified constraint reactions.
9. **Regularity and well-posedness.** On the operating domain it has enough
   continuity/coercivity for the chosen weak formulation and preserves the
   admissible state/tangent cone.
10. **No unauthorized scale.** Its normalization is derived from accepted
    structure or explicitly justified observational input, not inserted by
    relabeling a new constant.

**Conditional properties, not theorem requirements:** spatial locality,
ultralocality, strict reversibility in all regimes, Markovianity, and a
second-order differential form.  Dissipation may be present in an effective
branch, but it must be separated into a passive causal part; it cannot be
called inertial merely because it opposes motion.

### Linearized diagnostic

Without selecting an equation, a stable homogeneous reference branch must
possess a kinetic symbol whose conservative part is Hermitian and positive on
physical propagating polarizations.  Combined with the frozen acoustic symbol,
the resulting characteristic frequencies must be real at leading lossless
order.  This is a candidate-rejection test, not a formula for
\(\mathcal K_\tau\).

## 7. Scalar, tensorial, operator-valued, or state dependent?

The answer is hierarchical.

| Classification | PBUF status | Reason |
|---|---|---|
| scalar coefficient | permitted special case, not derived | requires an isotropic local placement realization with a simple kinetic metric |
| rank-two/tensorial inertia | permitted and generically expected in coordinates | the kinetic geometry maps state tangents to cotangents and can mix components |
| operator-valued | most general safe classification | constraints, projections, spatial nonlocality, and memory can all enter the tangent-to-cotangent map |
| state dependent | permitted, not required | configuration-space kinetic geometry may vary with \(q\); the frozen ontology does not impose homogeneity |
| frequency dependent | permitted only as an effective causal representation | normally signals eliminated internal dynamics and must satisfy passivity/causality constraints |

Intrinsically, inertia is a (possibly nonlinear) bundle map or variational
operator from admissible histories/tangents into the cotangent bundle.  Its
linearization is a bilinear form/operator.  Calling it scalar before choosing
the realization would discard possibilities that the frozen framework has not
eliminated.

Isotropy at the unloaded reference can reduce a local linear kinetic tensor to
scalar blocks on irreducible physical representations.  It cannot fix their
values, prove that only one block exists, or eliminate gauge projections.

## 8. Compatibility audit and stronger constraints on \(\mathcal K_\tau\)

| Candidate property | Frozen status after INERTIA-001 | Audit result |
|---|---|---|
| nonzero on accelerating physical modes | conditionally required by strict wave support | stronger than an unconstrained slot |
| positive on physical propagating modes | required for stable lossless waves | stronger |
| covector valued | already implicit in GE-001 | confirmed |
| gauge/objective | required by STATE/DYNAMICS foundations | confirmed and made explicit |
| second order in \(\tau\) | not forced | remains free |
| local in space | not forced | remains free; causal nonlocal closure allowed |
| local in time/Markovian | not forced | remains free |
| reversible | required only for the conservative baseline/strict lossless waves | effective damping remains allowed |
| symmetric | required for the conservative kinetic Hessian, not every full response term | refined |
| momentum-conserving | conditional on full-action translations and closed boundaries | cannot be claimed from \(W\) alone |
| scalar density | not forced | remains free |
| normalization from \(K_0,\mu_0\) | not derivable | forbidden without an added length/time relation |
| normalization from elastic energy density via \(c^{-2}\) | only a total effective closed-system matching statement | not a native local derivation |
| cone compatibility | required by FP-5 | stronger rejection gate |
| native reparametrization constraint | required by FP-3 and DYNAMICS-001 | stronger rejection gate |

GOVERNING-EQUATION-001 was correct not to fill the slot.  INERTIA-001 narrows
its admissible family by turning positivity, quotient covariance,
clock-reparametrization covariance, conservation compatibility, and the V11
causal cone into explicit gates whenever the intended dynamics is stable and
wave supporting.  It does not derive a coefficient or final operator.

## 9. Minimal inertial requirements and remaining freedom

### Minimal requirements

Any future closure must:

1. identify the complete instantaneous state and its physical tangent/cotangent
   pairing;
2. define momentum from calibrated medium change;
3. provide positive, nondegenerate kinetic storage on physical wave modes;
4. respect gauge, objectivity, the hard admissible tangent cone, and order
   reparametrization;
5. be causal and pass the effective one-metric cone test;
6. support the appropriate closed-system momentum/energy balances;
7. vanish on constant physical histories; and
8. derive or independently justify its normalization without violating FP-6.

### Remaining freedom

Still undetermined are:

- the complete phase-space variables;
- the native kinetic functional or symplectic/Poisson structure;
- the physical-duration calibration entering rates;
- local versus spatially nonlocal response;
- Markovian versus hereditary response;
- first-order phase-space versus second-order configuration representation;
- scalar, block-tensor, or general operator structure;
- state dependence and nonlinear velocity dependence;
- reversible baseline versus permitted effective damping;
- normalization and its relation to existing PBUF scales;
- the native momentum-to-effective stress-energy map; and
- proof that the characteristic cone agrees with the V11 effective metric.

## 10. Required additional principle

Stored elastic energy is insufficient.  The minimum missing principle is a
**kinetic momentum principle**:

> Assign to every admissible calibrated physical history of the complete
> medium an objective, gauge-basic momentum covector (or equivalent
> symplectic/kinetic structure), positive and nondegenerate on physical
> propagation modes, whose source-free change is balanced by the already
> derived internal stress flux, whose native representation is invariant under
> order relabeling, and whose effective characteristics and total
> energy-momentum match the single V11 Lorentzian metric.

This principle may ultimately be derived from a deeper microscopic
energy-momentum construction.  In the present frozen ontology it is additional
to \(W\).  Its normalization is the specific unresolved part that prevents
\(\mathcal K_\tau\) from becoming numerical.

## 11. Final physical-origin statement

PBUF inertia is not “the elasticity of time” and not resistance created by the
deformation potential.  It is the medium's capacity to carry momentum and
store energy in physical change.  Wave propagation makes that capacity
necessary; conservation explains why changing it requires force; the
relativistic limit requires total energy to contribute to the inertia of
closed excitations.  But a kinetic momentum principle is still required to
map a calibrated rate of the complete state to momentum.  Until that map and
its normalization are derived, \(\mathcal K_\tau\) remains a constrained
operator family rather than a free arbitrary slot or a completed equation.

## 12. Completion statement

The continuum survey, independence proof, wave-support necessity result,
emergent-time analysis, structural theorem, PBUF compatibility audit, minimum
requirements, and remaining freedom are complete.  The milestone closes the
physical interpretation of inertia without selecting a kinetic equation.  Its
negative result is exact: elasticity alone cannot derive inertia.  Its positive
result is also exact: a stable reversible wave-supporting PBUF completion must
possess positive causal momentum-carrying kinetic structure on its physical
modes.
