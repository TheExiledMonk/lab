# PBUF CONSTITUTIVE-PRINCIPLES-001

## Scope and method

This is a convergence analysis of all 44 MATERIAL-DISCOVERY-001 candidates. No candidate was omitted: rejected and ontology-changing entries are used as failure controls. Frequency is descriptive, not a vote. A principle is **Necessary** only when it follows from the frozen design inputs and is not defeated by a viable mathematical counterexample; **Conditional** is required only after choosing a class or global requirement; **Optional** supplies a permissible specialization; **Incompatible** contradicts a necessary principle or the frozen state.

## Result

The common core is not a named material model. It is an objective, invariant, differentiable stored-energy sector on the frozen SPD deformation domain, normalized at a stable reference state, with the frozen positive weak-field tangent. Stress follows by variation. BALANCE-001 supplies spatial communication (or an optional positive gradient/nonlocal energetic extension does), and DURATION-001 plus stability and communication permits waves. Conservative recovery follows from the stable single-well energetic branch.

Progressive hardening is **conditional**, not universal. Finite extensibility, gradient dependence, integral nonlocality and dissipation are **optional**. Independent internal structures, irreversible completion and zero-shear fluid completion are **incompatible** with the frozen minimal constitutive core.

## Dependency graph

```text
P01 objective C ──> P02 invariant Phi ──> P04 stress/tangent ──┐
                 P03 reference minimum ────────────────────────┼──> P07 recovery
P03 + P04 ──> P05 positive weak tangent ──> P06 domain audit  │
                         P05 + P08 communication + duration ───┴──> P09 waves

P13 gradient ──> P08       P14 integral nonlocal ──> P08
P12 barrier --conditional/pathwise--> P10 hardening
P16 extra state -X- P01    P17 irreversibility -X- P07    P18 no shear -X- P05
```

The graph is conjunctive: wave support is not downstream of recovery, and neither energy nor neighbour communication alone implies waves.

## Equivalence and independence

No nontrivial exact equivalence survives. A finite barrier conditionally implies asymptotic hardening only along paths that approach the guarded boundary; it does not ensure hardening in all tensor modes. Gradient elasticity implies communication, but communication also arises from local stress divergence or an integral operator. A gradient model is merely the long-wave expansion of some kernels, not exactly equivalent to general nonlocality. Hardening is independent of communication; waves are independent of finite extensibility; recovery does not require hardening or dissipation.

## Family reduction

The 30 catalogue labels reduce to one required core family and four genuine specializations: (F1) local invariant energetic; (F2) finite-domain barrier; (F3) gradient energetic; (F4) integral nonlocal energetic; and (F5) duration-compatible dissipative completion. Polynomial, exponential, logarithmic and spectral forms are representations inside F1, not distinct principles. “Hybrid” is a composition of these slots, not a sixth primitive family.

## Minimal constitutive architecture

1. **State/domain:** frozen `C[q,q0]` in its admissible SPD domain; no additional independent state.
2. **Scalar potential:** `Phi(C)=phi(I1,I2,I3)`, with `Phi(I)=0` and `D Phi(I)=0`.
3. **Regularity and tangent:** `Phi` is at least `C^2` in the interior; its Hessian at `I` equals the frozen positive weak-field elasticity tensor.
4. **Admissibility gates:** bounded-below energy and positive incremental/acoustic response on the declared operating domain; a global completion must use either coercive growth or a complete boundary barrier, but this choice is not yet fixed.
5. **Response map:** stress is the authorized variational derivative of `Phi`; no final functional form is selected.
6. **Communication slot:** the already frozen balance divergence is the minimum. A positive objective gradient term or symmetric positive nonlocal kernel is optional and mutually alternative at minimum architecture level.
7. **Evolution slot:** the frozen duration/kinetic closure combines with stress and communication. Dissipation may only be an additive, nonnegative, duration-compatible later completion that preserves equilibrium energy.

This architecture introduces no new free parameter. In particular, no barrier limit, gradient length, kernel horizon, relaxation time or fitted coefficient is chosen.

## Governing-equation readiness

**Outcome B — Partially constrained.** The principles fix the variational skeleton, reference conditions, weak-field Hessian, stability gates, conservative recovery requirement and the admissible locations of communication and dissipation. They do **not** uniquely determine `phi`, decide coercive versus finite-barrier completion, or decide balance-only versus gradient versus integral interaction. Those alternatives yield inequivalent nonlinear stresses and dispersion relations while sharing all necessary principles. Consequently a unique governing equation cannot yet be derived without an additional authorized selection/derivation criterion; the problem is not unconstrained, but it is not unique.

## Deliverable map

The machine-readable catalogue, frequency and necessity tables, dependency/equivalence graphs, independence analysis, family reduction, principle-elimination report, architecture above, and readiness assessment together satisfy the milestone deliverables.
