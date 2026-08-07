# PBUF NONLINEARITY-001 — Physical Origins of the Nonlinear Constitutive Response

## Decision

**Outcome B: a reduced family survives; no unique physical origin is selected.**

Within the frozen architecture, every admissible nonlinear mechanism must reduce to

\[
W(C)=W_2(C)+R(C),\qquad C\in\mathcal D_C\Subset\operatorname{Sym}^+(3),
\tag{N-001}
\]

where \(W_2\) is the frozen quadratic weak-field representative and

\[
R({\bf1})=0,\qquad DR({\bf1})=0,\qquad D^2R({\bf1})=0. \tag{N-002}
\]

It must also be objective, isotropic, parity-even, single-valued, local, and
rate independent. N-002 is the exact weak-field preservation test. For a
sufficiently smooth remainder, its first possible nonzero Taylor term is cubic
in \(C-{\bf1}\).

All admissible physical stories reduce to three response classes:

1. **D1 — reversible state-local anharmonicity** in the interior;
2. **D2 — an energetic finite-capacity barrier** at
   \(\partial\mathcal D_C\); and
3. **D3 — a kinematic finite-capacity endpoint**, hard or regular constrained.

Progressive stiffening, geometric alignment, finite extensibility, asymptotic
stiffening, and barrier formation can instantiate these classes. Their
microscopic explanations are not observable in the frozen state \(C[q,q_0]\),
so PBUF cannot distinguish them without unauthorized structure. Plastic work
hardening, history-dependent jamming, and genuinely dynamical cooperative-wave
mechanisms are rejected. Bounded saturation on an extendible domain is also
rejected.

No governing equation, coefficient, microscopic constituent, auxiliary field,
fit, V11 object, metric map, or weak-lensing implementation is introduced or
modified.

## 1. Audit standard

### 1.1 Frozen gate

HYPER-001 permits only

\[
W(C)=\Phi(I_1,I_2,I_3),\quad
(I_1,I_2,I_3)=(\operatorname{tr}C,
\tfrac12[(\operatorname{tr}C)^2-\operatorname{tr}(C^2)],\det C).
\tag{N-003}
\]

The physical deformation is the unordered spectrum of the dimensionless
objective relative tensor \(C\). A proposed origin is admissible only if all of
its constitutive effect is a single-valued function of that spectrum. External
examples below establish occurrence in ordinary materials; they are analogies,
not claims that PBUF contains their microscopic constituents.

### 1.2 Weak-field lemma

Let \(H=C-{\bf1}\). The frozen reference value, stress, and Hessian are unchanged
exactly when N-002 holds, equivalently

\[
R({\bf1}+H)=o(\lVert H\rVert^2). \tag{N-004}
\]

If \(R\in C^3\) near the reference, \(R=O(\lVert H\rVert^3)\). This applies to
the tensor function, not separately to each invariant-coordinate derivative.
The frozen stress-free condition remains
\(\Phi_1^0+2\Phi_2^0+\Phi_3^0=0\).

Subtracting a candidate's 0-, 1-, and 2-jets is a mathematical projection onto
the allowed remainder space, not a derivation or permission to add fitting
coefficients. A mechanism naturally preserves the weak field only when its
effective contribution can begin beyond quadratic order without an extra
cancellation parameter.

### 1.3 Large-deformation vocabulary

| Property | Constitutive meaning |
|---|---|
| progressive hardening | tangent stiffness increases along a specified loading path |
| coercive growth | energy diverges on every escaping sequence of an unbounded domain |
| complete barrier | energy diverges on every approach to a forbidden finite-domain boundary |
| finite extensibility | continuation ends at finite deformation, energetically or kinematically |
| stress/energy saturation | response or energy approaches a finite limit |
| deformation saturation | state reaches a finite-domain endpoint; stress need not saturate |

Increasing stiffness does not prove a finite endpoint, and bounded energy does
not prevent finite loading from crossing an otherwise extendible state.

## 2. Mechanism catalogue

### A. Progressive strain hardening

**Interpretation, origin, and occurrence.** Reversible *strain stiffening* is an
increasing incremental modulus caused by alignment, a bending-to-stretching
transition, or nonlinear force extension. Metallurgical *work hardening* is an
irreversible evolution of material structure. Collagen, fibrin, actin, and other
cross-linked biopolymer networks exhibit strong stiffening through both
filament-level and network-geometric routes [Janmey et al.][R1],
[Storm et al.][R2], [Onck et al.][R3].

**Nonlinearity and decision.** Reversible stiffening can give positive growth of
directional second derivatives of \(W\), with polynomial, power-like, rapid, or
mixed-invariant \(R\). It **survives conditionally** if it obeys N-002 and stays
acoustically positive. Plastic work hardening is **rejected**: plastic strain,
damage, accumulated work, or history makes equal \(C\) states respond
differently, contradicting STATE-002, HYPER-001, and DURATION-001. No onset
strain, exponent, or hardening modulus is supplied by the mechanism name.

**Large deformation.** Progressive hardening; possibly coercive growth; only
conditionally a barrier. It does not imply saturation or finite extensibility.

### B. Finite extensibility

**Interpretation, origin, and occurrence.** In polymeric exemplars, initially
coiled or fluctuating strands lose configurational freedom as contour extension
is approached. Rubber-like and semiflexible networks show the associated
large-strain stiffening [MacKintosh et al.][R4], [Meng and Terentjev][R5].

**Nonlinearity and decision.** At continuum level this produces divergent
force/energy at a finite spectral boundary, or a finite endpoint enforced by an
admissibility constraint. Logarithmic, rational, and algebraic divergences are
possible but the physical label selects none. The endpoint property
**survives**, since the finite \(\mathcal D_C\) and its three authorized endpoint
completions are frozen already. A literal chain-locking explanation is not
native: chains, contour lengths, thermal scales, and network statistics would
be new ontology or constitutive data. Weak-field preservation additionally
requires regularity at \({\bf1}\) and N-002; named chain laws do not pass
automatically.

**Large deformation.** Finite extensibility with a complete barrier or
hard/regular constrained endpoint; usually asymptotic stiffening, not response
saturation.

### C. Geometric locking

**Interpretation, origin, and occurrence.** Networks can stiffen even when
elements are linearly elastic: rotation and alignment transfer load from soft
bending modes to stiff stretching modes. Dense packings can lose rearrangement
directions as contacts crowd. Fiber-network experiments and simulations
establish these geometric sources [Onck et al.][R3], [van Dillen et al.][R6].

**Nonlinearity and decision.** Equilibrium homogenization can produce
mixed-invariant stiffening, loss of compliance near a spectral boundary, or a
hard boundary. It **survives conditionally** only if the relevant geometry is
already encoded by \(q\) and \(C\), and the response is single-valued and
isotropic. Contact topology, fabric tensors, coordination, rearrangement
history, or grain/fiber scales violate STATE-002, isotropy, locality, or FP-6.
With those removed, locking is indistinguishable from D1, D2, or D3.

**Large deformation.** Progressive hardening, possible finite extensibility
and barrier; not necessarily coercive and not saturation.

### D. Anharmonic elasticity

**Interpretation, origin, and occurrence.** An anharmonic restoring landscape
is not exhausted by its quadratic expansion. Crystals exhibit higher-order
elastic response, thermal expansion, phonon interaction, and nonlinear acoustic
effects; nonlinear lattices are established mathematical-physics exemplars
[FPU][R7], [Dauxois et al.][R8]. They are analogies, not PBUF microstructure.

**Nonlinearity and decision.** This mechanism most directly generates \(R(C)\):
higher-order invariant terms change finite-deformation stress while leaving the
reference 2-jet intact. Polynomial, analytic, rapid-growth, logarithmic-
correction, and mixed-invariant behavior are possible. It **survives natively
as D1**, requiring no variable beyond \(C\), no history, and no length scale.
But “anharmonic” is a property of the unknown remainder, not a physical
selection rule; it fixes neither sign, shape, coefficient, ellipticity, nor
endpoint.

**Large deformation.** Hardening, stable-branch softening, coercive growth, or
barrier behavior depending on the member. Saturation remains subject to the
endpoint rules.

### E. Cooperative wave interaction

**Interpretation, origin, and occurrence.** Finite-amplitude waves in nonlinear
media exchange energy among modes, self-modulate, and form coherent structures
in nonlinear lattices, fluids, solids, and optics [Dauxois et al.][R8],
[Kartashov et al.][R9].

**Nonlinearity and decision.** Wave mixing is normally a consequence of a
nonlinear constitutive law and dynamics, not a state function explaining that
law. A wave population needs amplitudes, phases, occupations, gradients,
dispersion, or history. It is **rejected as an independent origin**:
mode variables violate STATE-002; gradients/kernels violate LOCALITY-001;
frequency/rate dependence violates DURATION-001 and the equilibrium state
function; dispersive scales/couplings violate FP-6. If all dynamic variables
are eliminated and the result is exactly a local invariant of \(C\), it survives
only by reduction to D1 and the wave story has no distinguishable constitutive
content.

**Large deformation.** No intrinsic classification; focusing, saturation, and
mode crowding depend on the rejected dynamical closure.

### F. Elastic saturation

**Interpretation and origin.** “Saturation” can mean a finite limit of stress,
a finite limit of energy, or exhaustion of admissible deformation. Soft
networks may show plateaus through rearrangement, damage, unfolding, or phase
conversion; limiting deformation can instead coexist with divergent stress.

**Nonlinearity and decision.**

- Deformation saturation **survives** as D3, the already-authorized hard or
  regular constrained endpoint.
- Asymptotic stiffening **survives** as D2 when it is a complete barrier; this
  is divergent resistance, not saturation of response.
- Bounded energy/stress saturation on an extendible domain is **rejected**. It
  supplies no energetic obstruction to continuation and may create a vanishing
  or negative tangent.
- A finite endpoint with a plateau survives only when the independent frozen
  state constraint forbids continuation and the operational interior remains
  stable. Damage- or rearrangement-driven plateaus requiring history are
  rejected.

### G. Barrier formation

**Interpretation and origin.** An energy barrier represents rapidly increasing
resistance as a state exhausts an admissible region. Continuum barriers enforce
orientation preservation, noninterpenetration, limiting stretch, or other
domain restrictions. Here the only authorized physical reading is exhaustion
of the frozen spectral state domain.

**Nonlinearity and decision.** A complete barrier obeys

\[
C_n\to\partial\mathcal D_C\quad\Longrightarrow\quad W(C_n)\to+\infty
\tag{N-005}
\]

for every forbidden boundary approach. Logarithmic, algebraic, rational, and
faster divergences can realize N-005. It **survives natively as D2** using only
existing \(C,W,\mathcal D_C\). A singularity inside the declared operational
domain is rejected because stress fails to be \(C^1\) and the domain is
effectively split; it can survive only by being declared boundary. A partial
barrier is insufficient unless every other forbidden boundary is hard
constrained. Position and profile remain unselected constitutive data.

### H. Phase switching and multiwell elasticity

Martensitic and other transformations show nonlinear multiwell response. A bare
invariant \(W(C)\) can be multiwell and preserve N-002, so it is mathematically
admissible branchwise. But phase fractions, interfaces, nucleation, hysteresis,
and gradient regularization require extra state or length scales. The literal
mechanism is therefore **not selected**; stripped of that physics it is merely a
D1 multiwell member and does not naturally ensure unique recovery or capacity.

### I. Contact, crowding, and jamming

Dense particulate and cellular systems stiffen as free volume and rearrangement
modes disappear. The literal mechanism is **rejected** because constituents,
contact/fabric state, dissipation, and rearrangement history are unauthorized.
An equilibrium isotropic coarse-graining depending only on \(C\) survives, but
reduces to D1/D2/D3 and supplies no distinct origin. A nonsmooth contact surface
inside the operational domain also fails the classical stress/acoustic gates.

## 3. Weak-field compatibility audit

| Mechanism | Can satisfy N-002? | Reference/Hessian verdict |
|---|---:|---|
| reversible strain stiffening | yes | survives only when its contribution begins beyond the frozen 2-jet |
| plastic/work hardening | no as a function of \(C\) alone | rejected for history/internal state before jet comparison |
| finite extensibility | yes | a remote boundary contribution can have zero 2-jet; named laws do not automatically |
| geometric alignment/locking | yes after state-local reduction | literal geometry variables are rejected |
| anharmonic elasticity | yes, directly | canonical smooth source of a higher-order remainder |
| cooperative waves | no independently | extra dynamics alter the state; complete local elimination reduces to D1 |
| deformation saturation | yes with frozen hard constraint | endpoint can leave the local jet unchanged |
| bounded constitutive saturation | sometimes locally | rejected globally on an extendible domain; jet preservation is insufficient |
| complete barrier | yes | boundary divergence is compatible with a regular zero-2-jet reference neighborhood |
| phase/multiwell | yes branchwise | literal transformation physics still needs extra state |
| contact/jamming | only after smooth coarse-graining | literal nonsmooth/history-dependent response is rejected |

No mechanism passes because its correction is merely “small.” Its value,
gradient, and Hessian at the reference must vanish exactly.

## 4. Frozen-framework compatibility

| Frozen authority | Mathematical requirement | Incompatibility |
|---|---|---|
| FOUNDATION-001 | one medium, no external substrate, new sector, or free constant | reifying chains, grains, lattices, defects, or wave modes as PBUF constituents; adding onset scales/couplings |
| STATE-002 | equal \(C[q,q_0]\) gives equal energy and stress | internal variables, phase fractions, damage, fabric, mode occupation, or history |
| DEFORMATION-001 | dimensionless objective SPD relative deformation; spectral physical content | raw rotation/component dependence, preferred direction, or non-equivalent strain variable |
| HYPER-001 | \(W=\Phi(I_1,I_2,I_3)\), local, isotropic, parity-even, \(C^1\) inside and \(C^2\) near reference | anisotropy, rate law, hysteresis, interior cusp/pole, or non-invariant formula |
| DURATION-001 | equilibrium energy is independent of history parametrization and fundamental rate/frequency | work-per-cycle, relaxation time, frequency stiffness, or wave evolution in \(R\) |
| METRIC-001 | Lorentzian/nondegenerate effective metric and frozen V11 weak limit; compatible characteristics | cone degeneracy, multiple incompatible signal cones, or changed weak map |
| BALANCE-001 | exact variational response \(P_C=DW\) and acoustic positivity on propagation domain | non-exact stress, negative acoustic mode, undefined stress, or dissipation disguised as storage |
| LOCALITY-001 | intrinsic constitutive freedom is pointwise in \(C\) | \(\nabla C\), kernels, horizons, dispersive lengths, or cooperative nonlocal response |

METRIC-001 compatibility is a gate, not a theorem supplied by a mechanism. Each
explicit \(R\) would still need member-by-member ellipticity and cone audits.

## 5. Constitutive simplicity

| Candidate | New variable/state? | New scale/parameter? | Disposition |
|---|---:|---:|---|
| state-local anharmonicity | no | not structurally; a formula still carries unselected higher-order data | survives |
| reversible stiffening phenotype | no after reduction | onset/profile unselected | survives as D1 |
| finite-capacity boundary | no; already authorized | domain shape remains frozen-but-unselected data | survives |
| complete barrier | no | divergence profile remains unselected | survives |
| regular constrained endpoint | no | endpoint response remains unselected | survives conditionally |
| work hardening/plasticity | yes | usually yes | rejected |
| literal chain extensibility | yes | contour and thermal scales | literal origin rejected; endpoint survives |
| literal crowding/jamming | yes | constituent and packing scales | literal origin rejected; phenotype reduces |
| cooperative wave population | yes | frequency, dispersion, coupling scales | rejected independently |
| damage/phase saturation | yes | thresholds/interfacial scales | rejected |
| unguarded bounded saturation | not necessarily | profile data | rejected for failure to enforce capacity |

“No new parameter” does not mean arbitrary \(R\) contains no constitutive
information. This milestone introduces none: higher derivatives, domain shape,
and endpoint profile remain unresolved data inherited from ENERGY-SEARCH-001.

## 6. Convergence analysis

| Asymptotic structure | Independent physical stories | Frozen status |
|---|---|---|
| smooth polynomial/power growth | anharmonic forces; reversible network stiffening; reduced alignment | admissible memberwise after zero-2-jet, positivity, and ellipticity checks |
| rapid/exponential-type growth | strong anharmonic stiffening; cooperative effects after complete local elimination | same D1 class; formula name has no priority |
| logarithmic correction without pole | anharmonic response; fluctuation-inspired analogy | admissible only if regular and zero-2-jet |
| logarithmic barrier | finite extensibility; state exhaustion; locking | D2 if every forbidden path is covered |
| algebraic/rational barrier | finite chain extension; crowding; limiting stretch | same D2 asymptotic class, though stresses differ |
| finite regular endpoint | deformation saturation; kinematic locking | D3 only with frozen state constraint |
| finite limit on extendible domain | response saturation; damage-like plateau | rejected as standalone capacity mechanism |
| multiwell/nonmonotone response | phase switching; competing configurations; anharmonicity | branchwise D1 survivor; literal phase physics rejected |

Two equivalence notions are required:

1. **Asymptotic equivalence:** different functions share polynomial, rapid,
   barrier, or finite-limit behavior. Their finite-deformation stresses can
   still differ.
2. **Frozen-observability equivalence:** physical origins yielding identical
   \(W(C)\) on \(\mathcal D_C\) are indistinguishable because no authorized
   hidden variable records their origin.

Thus chain locking, geometric locking, and state-space exhaustion can converge
to one complete barrier. Alignment, higher restoring forces, and fully reduced
collective response can converge to one smooth anharmonic hardening law. The
framework classifies responses but cannot infer their microscopic story.

## 7. Minimal reduction

| Reduced class | Includes | Excludes | PBUF content |
|---|---|---|---|
| **D1 — regular interior anharmonicity** | reversible stiffening, reduced alignment, higher-order elasticity, eliminated collective effects, branchwise multiwell behavior | history, rates, gradients, hidden modes | local invariant higher-order shape of \(R\) |
| **D2 — energetic finite-capacity barrier** | finite extensibility, limiting deformation, locking, state exhaustion with complete divergence | partial/interior poles and unguarded paths | divergence of the scalar energy at every forbidden boundary approach |
| **D3 — kinematic finite-capacity endpoint** | hard extension and finite regular endpoint, including deformation saturation | bounded saturation without constraint | already-authorized boundary of \(\mathcal D_C\) |

D2 and D3 remain physically distinct: D2 makes the boundary energetically
inaccessible, while D3 makes it inadmissible by state definition. D1 may be
combined with either endpoint class and is not an additional constitutive
object. Plastic, dissipative, dynamic-wave, and microstructure-dependent
mechanisms do not form survivor classes.

## 8. Native selection assessment

The result is **B, a reduced family survives**, with the stronger epistemic
conclusion that **no physical origin is uniquely selectable**.

The frozen ontology mildly favors descriptions intrinsic to the one medium:
state-local anharmonicity and exhaustion of its admissible domain require the
least translation. That economy is not a selection theorem. It supplies no
proposition distinguishing:

- smooth hardening from a remote barrier;
- energetic barrier from hard kinematic endpoint;
- polynomial from exponential interior growth; or
- chain-like, geometric, and collective analogies after reduction to the same
  scalar \(W(C)\).

Unique selection would require an additional authorized fact: boundary
accessibility/regularity, a finite-deformation response datum, a theorem mapping
complete \(q\) to configuration geometry, or a new symmetry restriction on
higher jets. None exists. Inferring an origin from “one medium,” “finite
capacity,” or “waves exist” would affirm the consequent because multiple
inequivalent mechanisms share those consequences.

## 9. Recommendations for the next milestone

1. Preserve \(R\) as an invariant function satisfying the exact zero-2-jet
   condition; do not derive a final law.
2. Select the next discriminator rather than a formula family: energetic D2
   versus kinematic D3 capacity, and the operational domain of strong
   ellipticity.
3. Use parameter-free property tests: monotonicity along volumetric and
   isochoric spectral paths, complete boundary coverage, lower semicontinuity,
   and lifted rank-one acoustic positivity.
4. Keep chains, lattices, fibers, grains, and wave modes as evidence that
   response classes are physically established, never as PBUF constituents.
5. If physical selection is required, first authorize evidence that separates
   D1–D3 while preserving V11, the weak 2-jet, locality, and existing state.
   Observational fitting and weak-lensing modification remain outside scope.

## 10. Traceability and references

Frozen sources: runs/foundation001/downstream_reference_contract.md;
runs/state002/primitive_medium_state.md;
runs/deformation001/deformation_measure_report.md;
runs/hyper001/stored_energy_derivation.md;
runs/energy_principle001/energy_selection_derivation.md;
runs/duration001/emergent_duration_derivation.md;
runs/metric001/effective_metric_derivation.md;
runs/balance001/native_balance_laws.md;
runs/locality001/locality_report.md; and
runs/energy_search001/energy_search_report.md.

Evidence for established analogues:

[R1]: https://doi.org/10.1016/j.ceb.2013.06.002 "Janmey et al., Effects of nonlinearity on cell–ECM interactions"
[R2]: https://doi.org/10.1038/nature04036 "Storm et al., Nonlinear elasticity in biological gels"
[R3]: https://doi.org/10.1103/PhysRevLett.95.178102 "Onck et al., Alternative explanation of stiffening in cross-linked semiflexible networks"
[R4]: https://doi.org/10.1103/PhysRevLett.75.4425 "MacKintosh et al., Elasticity of semiflexible biopolymer networks"
[R5]: https://doi.org/10.3390/polym8020052 "Meng and Terentjev, Theory of semiflexible filaments and networks"
[R6]: https://doi.org/10.1016/j.bpj.2008.01.009 "van Dillen et al., On the stiffening of semiflexible networks"
[R7]: https://doi.org/10.2172/4376203 "Fermi, Pasta, and Ulam, Studies of nonlinear problems"
[R8]: https://doi.org/10.1016/j.physrep.2004.09.001 "Dauxois et al., The Fermi–Pasta–Ulam nonlinear lattice"
[R9]: https://doi.org/10.1103/RevModPhys.83.247 "Kartashov et al., Solitons in nonlinear lattices"

## Completion

All requested groups and two additional established classes have been
catalogued. Their physical meaning, origin, occurrence, constitutive behavior,
weak-field jets, large-deformation asymptotics, frozen compatibility, and
structural costs were audited. Equivalent mechanisms reduce to D1–D3. The
result is Outcome B with no unique native physical selection.
