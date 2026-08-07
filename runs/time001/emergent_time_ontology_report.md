# PBUF TIME-001 — Emergent time from medium evolution

## Decision

**Conditional consistency; emergence is not established by the current corpus.** A PBUF ontology can consistently take the medium to be ontically three-dimensional and represent change by an ordered curve through a space of instantaneous medium states. The curve parameter need not be a coordinate or a fourth constituent of the medium. This establishes logical compatibility, not a derivation of physical time.

The exact missing principle is a **relational temporal-identification principle**:

> Specify a complete gauge-invariant instantaneous medium state, a lawful and causal relation among such states, and a monotonic internal clock functional whose intervals coincide in the effective regime with the proper time of the single Lorentzian metric retained by V11, independently of arbitrary parametrization and admissible slicing.

Current PBUF does not supply this principle. In particular, it does not define complete local medium dynamics, the clock functional, or the normalized medium-to-metric map. TIME-001 therefore does not promote an ordering parameter to physical time, change V11, or introduce a field equation.

## 1. Authority and scope

V11 is authoritative. It explicitly retains Special Relativity, General Relativity, operational Lorentz invariance, Einstein's equations, and standard quantum dynamics. Its calculations use scale factor `a`, expansion rate `H(a)`, thermal evolution, and relativistic propagation. Thus V11 already uses effective relativistic time operationally. V11's elastic-spacetime language is a constitutive interpretation; it does not give a local microscopic medium state, an evolution law for that state, or a derivation of clocks.

The later PBUF ontology takes the medium to be fundamentally three-dimensional. FND-003 shows that three-dimensionality alone does not select the representation or number of components of a microscopic state. GEOMETRY-001 finds no normalized map from a medium variable to one Lorentzian metric. KINEMATICS-001 finds no unique primitive configuration, unloaded reference, or choice between three material directions and four clock-and-ruler directions. These boundaries are inherited here.

This is an ontology and representation audit only. All equations below are definitions or structural conditions, not gravitational or medium field equations.

## 2. Ontology versus representation

The ontological claim is:

- at any occurrence, the medium has an instantaneous physical configuration on a three-dimensional spatial carrier `Sigma`;
- change is a relation among configurations; and
- physical duration, if emergent, is read from change rather than added as a fourth material direction.

A mathematical representation may nevertheless use `lambda`, four-index notation, a four-manifold, or a foliation. Such devices do not by themselves add a fourth ontic medium dimension. Conversely, merely writing `q(x,lambda)` does not show that time has emerged: `lambda` may be only an external Newtonian time renamed.

Let `C` be a candidate set of medium configurations on a three-manifold `Sigma`, and let `G` denote spatial-coordinate and internal gauge redundancies. The instantaneous physical configuration is an equivalence class

`[q] in C/G`.

If future dynamics require velocities, canonical momenta, memory variables, or other Cauchy data, the complete instantaneous state is instead

`X = [q, p, ...] in Q_phys`.

No current PBUF result determines the entries represented by the ellipsis. Therefore the minimum object justified now is **an equivalence class of sufficient Cauchy data on a three-dimensional carrier**, not specifically the historical three-component `q`, the scalar `u`, or a three-metric alone. A configuration is a state only if it suffices to determine all physically allowed continuations under the eventual law.

## 3. Ordering without a fundamental fourth dimension

A history can be represented as an oriented curve

`gamma: I -> Q_phys`, `lambda |-> X(lambda)`.

The ordering relation is inherited from the orientation of the interval `I`: `X_1` precedes `X_2` along a chosen history. To prevent `lambda` from becoming hidden absolute time, all strictly increasing reparametrizations `lambda -> f(lambda)` must represent the same physical history. Then only ordering and relational coincidences are invariant; the numerical value and rate of `lambda` are gauge.

This construction needs no fourth material coordinate. It does, however, assume more than a set of states: an orientation, an admissible-history relation, and enough structure to distinguish gauge change from physical change. An unordered set of configurations contains no time.

The construction also distinguishes three levels that must not be conflated:

1. **Sequence label:** arbitrary `lambda`; mathematically useful and physically unobservable.
2. **Relational clock:** a monotonic functional `T_C[X]` or a subsystem correlation used to compare changes; potentially physical but not yet relativistic time.
3. **Effective relativistic time:** local proper-time intervals supplied by the single Lorentzian metric to which all clocks and matter couple in the V11 regime.

Emergence requires a derivation linking levels 2 and 3. Replacing `t` by `lambda` supplies only level 1.

## 4. Emergent time versus coordinate time

Coordinate time is a chart label `x^0` on an effective spacetime description. It depends on coordinates or foliation and is not itself an observable duration. Proper time is a metric interval along a timelike worldline and is invariant under coordinate changes. In the proposed ontology, emergent time is the physical ordering-and-duration structure abstracted from medium change and realized operationally by internal clocks.

Consequently:

- emergent time must exist before a particular coordinate time is selected;
- many coordinate times may represent the same emergent history;
- a state-label `lambda` is neither coordinate time nor proper time until a clock/metric identification is supplied; and
- the effective four-dimensional spacetime of GR may remain the correct macroscopic representation even if it is not the fundamental ontology.

## 5. Comparison with standard 3+1 GR

The comparison is foundational, not adversarial. In standard 3+1 GR, a four-dimensional Lorentzian spacetime is foliated into spatial hypersurfaces for an initial-value representation. Lapse and shift encode how successive slices are related; different admissible slicings can describe the same geometry. GR does not normally privilege one global time coordinate, while proper time and causal structure are geometric.

PBUF's proposed ontology reverses the explanatory direction: it begins with three-dimensional medium states and seeks to derive the effective Lorentzian clock, ruler, and causal structure. The mathematical objects can look similar to a 3+1 decomposition, but their status differs. In GR, the slices are parts of a four-geometry; in emergent-time PBUF, the four-geometry would have to be reconstructed as an effective representation of state correlations.

This reversal is consistent only if it preserves the V11 relativistic limit. It must not introduce a detectable preferred foliation, an external universal clock, or inequivalent predictions for descriptions related by diffeomorphisms.

## 6. What is emergent and what remains fundamental

Under the candidate ontology, the following may be emergent:

- numerical duration and clock rate;
- the local proper-time functional;
- Lorentzian causal cones and time dilation;
- coordinate-time descriptions and foliation choices;
- the four-dimensional spacetime history as a macroscopic reconstruction.

The following must remain primitive unless a deeper derivation supplies them:

- the existence of distinguishable physical change;
- an orientation or precedence relation on admissible histories;
- the law or constraint selecting admissible histories;
- the separation of physical change from gauge transformation; and
- the conditions allowing a subsystem to function monotonically as a clock.

Thus the proposal can make **metric duration** emergent, but it cannot derive change from nothing. Calling all ordering emergent would be circular unless the non-temporal admissibility relation is independently defined.

## 7. Past-state reconstruction

Emergent time does **not** by itself prohibit reconstruction of arbitrary past global states. If the complete state and evolution relation are deterministic, invertible, and exactly known, a past state may be mathematically reconstructible even when time is emergent.

Non-reconstructibility follows only with additional conditions, any of which remains unestablished in PBUF:

- coarse graining is many-to-one and microscopic information is discarded;
- effective evolution is dissipative or stochastic;
- horizons or causal disconnection prevent access to global Cauchy data;
- gauge equivalence prevents a unique global representative;
- the clock is only locally monotonic; or
- no global foliation/Cauchy surface exists.

These conditions can make operational reconstruction impossible or non-unique. They do not prove destruction of fundamental information. TIME-001 therefore rejects a claimed automatic arrow of time or automatic loss of the past.

## 8. Recovery of V11 relativistic time

The emergent description recovers effective relativistic time used in V11 only if all of the following hold:

1. **Complete state:** a gauge-invariant Cauchy state on `Sigma` is defined.
2. **Relational evolution:** admissible histories and their orientation are defined without relying on an unobservable external time, and reparametrizations leave physics unchanged.
3. **One metric:** a normalized, diffeomorphism-equivariant medium-to-metric map yields one nondegenerate Lorentzian metric with the required signature.
4. **Universal clocks:** all admissible matter clocks and light couple to that same metric, so their relational readings agree with metric proper time in the effective regime.
5. **Causality and well-posedness:** the completed dynamics has a causal initial-value formulation compatible with the metric cones and does not introduce observable preferred slicing contrary to V11.
6. **GR limit:** the effective metric and total stress-energy obey the V11-retained Einstein/Noether/Bianchi architecture, including standard local Lorentz and weak-field behavior.
7. **Cosmological identification:** the scale factor and `H(a)` used by V11 are recovered as effective observables of the state history, not used circularly to define its microscopic evolution.

TIME-001 does not assert that these gates are satisfied. Items 3 and 4 overlap the missing principle identified by GEOMETRY-001; items 1, 2, and 5 are new temporal closure requirements exposed by this audit.

## 9. Answer to the research questions

1. **Minimum instantaneous object:** an equivalence class of sufficient Cauchy data on a three-dimensional carrier, `X in Q_phys`; its detailed field content is unresolved.
2. **Ordering parameter:** yes. An oriented, reparametrization-equivalent history orders states without making the parameter a fundamental dimension. This supplies order, not yet physical duration.
3. **Difference from coordinate time:** emergent time is invariant relational/clock structure; coordinate time is a chart or foliation label; proper time is the effective metric realization that must be derived.
4. **Past reconstruction:** no natural prohibition follows from emergence alone. Non-invertibility, coarse graining, causal limits, or incomplete data must be independently established.
5. **V11 recovery:** only through a universal one-metric clock identification, causal well-posed dynamics, foliation/gauge consistency, and the retained GR limit.

## 10. Conclusion

The proposed ontology is internally coherent as a **conditional architecture**: three-dimensional instantaneous states plus a reparametrization-invariant ordering relation can underlie an effective four-dimensional description. The current PBUF framework does not establish that architecture physically because it lacks the relational temporal-identification principle and its underlying complete dynamics. This exact missing principle, rather than a contradiction, is the TIME-001 result.

