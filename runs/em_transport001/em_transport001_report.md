# PBUF EM-TRANSPORT-001 -- Native Electromagnetic Transport of the Spacetime Medium

## 0. Decision

**Outcome B.** The neighbour-to-neighbour transport law required for
weak-lensing wavefront evolution is **not** mathematically identical
to, or derivable from, the local propagation equations implied by
the V11 electromagnetic microscopic structure.

The factor `alpha_resolved ~ 3 alpha_EM` is a numerical identity and
a motivating dimensional-counting argument in V11 section 2.3.1.
The CORE-001 formalization supplies a three-component microscopic
state `q in R^3` with `g_dev = 1/137` as the matter-vertex coupling,
but the microscopic free energy has a mass-like onsite term and a
scalar nearest-neighbour term of the form `kappa_1 |q_j - q_i|^2`,
not the gauge-invariant curl form `|curl A|^2` of a Maxwell field.
The CORE-001 local evolution `tau dq_i/dt = -d(F/epsilon_*)/dq_i`
is overdamped and first-order in time; the coarse-grained field
satisfies the time-independent Helmholtz equation
`K u - Div(G grad u) = s(rho)`.

The exact missing local physical principle is therefore **the
kinetic sector** that supplies positive momentum density (or an
equivalent symplectic structure).  It was already identified by
INERTIA-001 as the irreducible closure gap left by the static
elastic energy.  EM-TRANSPORT-001 confirms that no derivation of
that kinetic sector from `alpha_EM`, `alpha_resolved`, or `g_dev`
is available inside the V11 microscopic structure as frozen.

No ontology, field, coupling, length, kernel, fit, V11 change, or
weak-lensing change is introduced.

## 1. Inputs

The audit cites only frozen sources (see `validation.json`):
FOUNDATION-001 (FP-1, FP-5, FP-6), STATE-002, DEFORMATION-001,
HYPER-001, BALANCE-001, LOCALITY-001, INERTIA-001, DURATION-001,
DYNAMICS-001, EQUILIBRIUM-001, ENERGY-SEARCH-001, PHOTON-001,
CORE-001, the V11 preprint, and the V11-ALPHA-001 brief.

## 2. What V11 actually says about the microscopic structure

The full per-record inventory of what the frozen corpus asserts
about the V11 / CORE-001 microscopic structure is given in
`microscopic_structure_audit.csv`.  The decisive rows are:

* `M-V11-01`.  V11 equation (4) states
  `alpha_resolved ~ 3 alpha_EM = 3/137.036 ~ 0.0219`.  V11 itself
  classifies the factor of three as a "motivating consistency
  argument" with the QFT derivation deferred.
* `M-V11-05`.  V11 section 2.4 records the GW170817 multimessenger
  constraint that gravitational and electromagnetic waves propagate
  as wave modes of the same medium with `epsilon_0 ~ 1`.  This is a
  constraint on a parameter, not a derivation of a propagation law.
* `M-V11-06`.  V11 equation (16) writes `Omega_b0 = 2 alpha_resolved`
  and attributes the factor of two to the two transverse EM
  polarizations.  This is polarization counting; it is not a
  structural identification of the microscopic field with an EM
  vector potential.
* `M-CORE-03`.  CORE-001-E01 introduces
  `F = epsilon_* sum_i [kappa_0|q_i|^2/2 + kappa_1 sum_<ij>|q_j - q_i|^2/2
  - g_dev eta_i e.q_i]`.  The gradient term is a SCALAR gradient
  `|q_j - q_i|^2`, not the curl form that would be required for an
  EM vector potential.
* `M-CORE-04`.  CORE-001-E02 gives the local evolution
  `tau dq_i/dt = -d(F/epsilon_*)/dq_i + xi_i`.  This is first-order
  in time; it relaxes, it does not propagate.
* `M-CORE-06`.  CORE-001-E09 gives the coarse-grained field equation
  `K u - Div(G grad u) = s(rho)`.  This is a Helmholtz-type
  elliptic equation; it has no time derivative; no wavefront
  follows.

## 3. Per-mechanism audit of the four transport questions

The full per-mechanism classification is in `native_transport_audit.csv`.
The conclusions are:

1. **Local phase transfer.**  Not derivable as EM-like transport.
   The nearest-neighbour term `kappa_1|q_j - q_i|^2/2` couples
   amplitudes, but CORE-001-E02 evolves them overdamped.  A
   second-order kinetic sector, or a Maxwell-like first-order
   structure, is required.
2. **Local field rotation.**  Not derivable.  The triplet `q`
   admits rotations, but CORE-001-E02 does not propagate rotations
   coherently; a non-dissipative dynamics is required.
3. **Neighbour coupling.**  Present at the energy level
   (CORE-001-E01).  LOCALITY-001 already established that
   `Div(P_F)` supplies all required static communication without
   invoking this term.  Neighbour coupling is therefore a static
   modelling choice, not a transport law.
4. **Wavefront evolution.**  Not derived.  CORE-001-E02 is
   overdamped; CORE-001-E09 is elliptic.  The V11 numerical
   identity `alpha_resolved ~ 3 alpha_EM` does not supply a
   second-order time structure.  INERTIA-001 already identified
   the kinetic sector as the missing closure.

## 4. EM-local microscopic mechanism audit

The standard local mechanisms of electromagnetism are listed in
`em_local_microscopic_mechanism.csv`.  None of them is present in
the V11 / CORE-001 microscopic structure:

| EM mechanism           | Present? | Structural mismatch                              |
|------------------------|----------|--------------------------------------------------|
| Faraday induction      | no       | no antisymmetric pair, no curl, no time derivative on a field-strength |
| Ampere-Maxwell         | no       | no current, no curl operator                     |
| D'Alembertian / wave   | no       | no kinetic sector, no Lorentzian signature       |
| Gauge invariance       | no       | `kappa_0|q|^2` is a Proca-like mass term; `kappa_1|q_j - q_i|^2` is not the curl form |
| Two polarizations      | counting only (V11 eq. 16) | counting is not a derivation of the transport law |
| Dispersionless c       | constraint only (V11 sec. 2.4) | V11 uses GW170817 to fix `epsilon_0 ~ 1`; it does not derive `c = sqrt(G/K)` |

The mismatch is structural, not numerical.  CORE-001's `q in R^3`
plus `kappa_0|q|^2` plus `kappa_1|q_j - q_i|^2` is the form of a
massive scalar triplet, not an EM vector potential.

## 5. Wavefront evolution audit

The candidate laws for `u(x,t)` and their derivation status are
recorded in `wavefront_evolution_audit.csv`.  Of the six entries:

* `WE-002` and `WE-003` are present in the frozen corpus but are
  time-independent (overdamped relaxation and Helmholtz equilibrium).
* `WE-001` is the elastic wave equation accepted by LOCALITY-001
  L-003 but is derived only **after** a positive momentum density is
  supplied.  INERTIA-001 left that supply open.
* `WE-004` would be a Maxwell-like wave equation and is absent.
* `WE-005` would require an independent flow vector and is
  forbidden by FP-1 / FP-4.
* `WE-006` is the numerical identity `alpha_resolved ~ 3 alpha_EM`
  itself, which is not a differential law.

No row in the audit produces a wavefront from the V11 microscopic
structure alone.

## 6. The missing principle, identified precisely

`kinetic_closure_requirement.json` records the precise gap:

> A local, conservative, second-order-in-time kinetic sector for
> the medium that supplies positive momentum density (or an
> equivalent symplectic structure) and thereby turns the static
> constitutive chain into a wave-bearing evolution equation.
> Equivalently, a Maxwell-like first-order structure with a curl
> kinetic operator and a conserved field-strength pair, supplied
> without introducing an independent EM sector.

This is exactly the closure gap that INERTIA-001 left open.  The
present milestone re-derives it from the EM side of the
microscopic structure: `alpha_EM = 1/137` fixes only the static
matter-vertex coupling; `alpha_resolved = 3 alpha_EM` fixes only an
amplitude identity; `kappa_1 |q_j - q_i|^2` fixes only a static
elastic energy; none of them supplies a positive momentum density,
a symplectic structure, a curl operator, or a Lorentzian signature.
The kinetic closure is unavoidable.

## 7. Compliance with the milestone brief

| Constraint                                       | Status |
|--------------------------------------------------|--------|
| No new spacetime inertia                         | yes    |
| No phenomenological steering coefficients         | yes    |
| No metric ansaetze                                | yes    |
| No cosmology solution                             | yes    |
| No quantum mechanics solution                     | yes    |
| No independent EM sector separated from medium    | yes    |
| No free transport constants                       | yes    |
| No new ontology                                   | yes    |
| No fit to data                                    | yes    |
| No V11 modification                               | yes    |
| No weak-lensing modification                      | yes    |

The audit uses only frozen sources, traces every claim back to a
frozen artifact, and reports a single explicit decision (`decision.json`).
The validation record is in `validation.json`.

## 8. Closure

**Outcome B.**  The neighbour-to-neighbour transport law is not
contained in the V11 electromagnetic microscopic structure.
The missing local principle is precisely the kinetic closure
already flagged by INERTIA-001.  The decision is reported in
`decision.json` and the completion record in `validation.json`.
