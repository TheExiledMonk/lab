# DEV DOC — Native Medium Response Accumulation 001

## Purpose

Investigate the missing native PBUF stage between local medium loading and long-range accumulated medium response.

This document deliberately avoids constructing or naming a native gravity field. In PBUF, gravity is emergent. The native causal chain under investigation is:

`matter / stress-energy -> local medium loading -> spatially accumulated medium response -> emergent macroscopic gravity-like behavior`

The target of this development step is therefore the middle native bridge only:

`local medium loading -> long-range accumulated medium response`

## Current evidence

The preceding native-response fingerprint comparison established that `c_state` is the strongest current candidate for a local medium-loading/state variable.

Without using G, macroscopic benchmark amplitudes, native amplitude rescaling, fitting, lensing targets, Rmax, Quantum Engine, or Planck-scale inputs, `c_state` reproduced four of six frozen effective-response fingerprints:

- local loading scales linearly with source density: exponent approximately `+1`;
- response scales linearly with source mass at fixed radius: exponent approximately `+1`;
- far response scales linearly with source mass at fixed probe radius: exponent approximately `+1`;
- two-source `c_state` superposition is linear to numerical precision.

It did not reproduce the required spatial accumulation behavior:

- fixed-mass surface-radius scaling measured approximately `R^-2.821`, rather than the frozen effective-response fingerprint `R^-1`;
- fixed-source far-radius scaling measured approximately `r^-24.138`, rather than `r^-1`.

The signed M10 interface-vector pipeline showed weaker agreement:

- density exponent approximately `+1`;
- mass exponent approximately `+1`;
- far-mass exponent approximately `+1`;
- fixed-mass radius exponent approximately `R^-3.293`;
- far-radius exponent approximately `r^-18.769`;
- worst signed-component superposition residual approximately `0.975`.

These results support a working interpretation that `c_state` may encode local medium loading/state, while the currently tested M10 representation is not yet identified as the missing long-range accumulated-response variable.

This interpretation is structural only. It does not establish SI dimensions, absolute amplitude, microscopic substrate ontology, or the physical origin of the coupling.

## Scientific question

Does the existing native PBUF pipeline already contain a downstream operation, state, transform, relaxation stage, or accumulated quantity that converts the approximately local `c_state` loading into a long-range medium response with the frozen spatial fingerprints?

The required response structure is:

- source-density linearity retained;
- source-mass linearity retained;
- fixed-mass source-radius scaling approximately `R^-1`;
- fixed-source far-radius scaling approximately `r^-1`;
- far response remains linear in source mass;
- weak-field / low-loading superposition remains approximately additive.

No amplitude requirement is imposed in this step.

## Epistemic split

### Derived / testable native layer

`known native loading -> native downstream operation -> accumulated medium response fingerprint`

### Speculative layer

`why stress-energy couples to the underlying spacetime substrate in that particular way`

The speculative layer must remain explicitly separate. A successful response law does not identify what spacetime is made of.

## Gravity guardrail

Gravity must not be used as a native construction variable or native field label in this development step.

Allowed wording:

- `local medium loading`;
- `medium state`;
- `accumulated medium response`;
- `long-range medium response`;
- `effective response fingerprint`;
- `emergent macroscopic gravity-like behavior` only when discussing the later observable layer.

Avoid wording such as:

- `native gravity field`;
- `final gravitational deformation field`;
- `gravity source field`;
- any statement implying gravity is fundamental in PBUF.

## Frozen response fingerprint

The next implementation must use the already established response-shape targets only as comparison fingerprints, not as amplitude calibration inputs:

- `local_density_exponent = +1`;
- `surface_mass_exponent_fixed_R = +1`;
- `surface_radius_exponent_fixed_M = -1`;
- `far_mass_exponent_fixed_probe_r = +1`;
- `far_radius_exponent_fixed_source = -1`;
- low-loading additivity residual approximately zero.

The comparison window may be predeclared before execution, but must not be adjusted after seeing results.

## Primary mission

Create a fact-finding lab tentatively named:

`native_medium_response_accumulation001.py`

The lab must inventory and test existing downstream native operations reachable from `c_state` before proposing any new propagation law.

Priority order:

1. existing native scalar downstream states;
2. existing native vector/tensor/interface states;
3. existing relaxation, transport, iterative, projection-independent, divergence, potential-like, integral, or accumulated representations;
4. existing operators that mathematically integrate local loading over space;
5. only if none exist, report the missing operator explicitly.

Do not invent a new kernel merely because `1/r` is known to be required.

## Candidate inventory rule

For every tested candidate, record:

- exact repository symbol / function / file provenance;
- whether it is directly computed from `c_state` or from another downstream state;
- whether the transform is local, derivative, transport, iterative, integral, inverse-operator-like, projected, or nonlinear;
- whether any coefficient or normalization is embedded in the implementation;
- whether the candidate is scalar, vector, tensor, projected, or representation-only;
- whether it is appropriate for 3D native interpretation.

Projected observer-space quantities must not be promoted into native medium variables.

## Required synthetic tests

Use controlled synthetic sources so source loading and geometry are known by construction.

At minimum:

### 1. Density ladder

Fixed source radius, varying density/loading amplitude.

Measure candidate exponent against:

`response ~ rho^1`

### 2. Mass ladder at fixed radius

Fixed radius, vary native integrated source by varying loading amplitude.

Measure:

`response ~ M^1`

### 3. Radius ladder at fixed native mass

Change source radius while preserving integrated native source.

Measure candidate response near the source boundary against:

`response_surface ~ R^-1`

### 4. Far-radius ladder

Fixed source, vary probe radius outside the source.

Measure:

`response_far ~ r^-1`

### 5. Far-mass ladder

Fixed probe radius, vary source mass/loading.

Measure:

`response_far ~ M^1`

### 6. Two-source superposition

Construct two sources independently and together.

For scalar candidates compare:

`X_12` versus `X_1 + X_2`

For vector/tensor candidates compare signed components before taking magnitudes.

Do not use magnitude-only superposition tests when signed native components are available.

## Important distinction: local loading versus accumulation

A candidate that preserves density and mass linearity but decays too rapidly with source radius or probe radius should be classified as a local/short-range medium response, not rejected as physically meaningless.

A candidate that reproduces `R^-1` and `r^-1` while preserving mass linearity and additivity is a candidate accumulated-response representation.

Do not automatically identify such a candidate with metric strain or physical spacetime deformation. That requires a later bridge.

## Existing evidence to preserve

The previous result suggests:

`c_state -> local medium loading/state candidate`

This is not yet a final physical assignment. The next lab may strengthen, weaken, or qualify it.

The previous M10 result must also remain open. Its current vector form failed the complete fingerprint, but another downstream M10-derived quantity may still carry useful accumulation structure.

Do not discard a family merely because one representation failed.

## No-fit rule

This stage is structural fact-finding.

Forbidden:

- multiplying a candidate by a chosen constant to improve agreement;
- selecting a radial exponent or kernel after seeing the result;
- choosing a smoothing length to produce `1/r`;
- adjusting iteration count to match the target unless the iteration count already has independent repository provenance;
- solving for coefficients from the frozen response benchmark;
- rescaling native amplitude to G, h00, Newtonian potential, lensing, or any observer target.

A partial or null result is valid.

## Macroscopic benchmark rule

`G` is not needed in this lab and should not appear in the native calculation.

The frozen response exponents are shape/scaling constraints only.

Do not use measured gravity amplitude to normalize or rank candidates.

## Other hard guardrails

- gravity is emergent, not fundamental;
- no Rmax;
- no replacement future-turnaround/free activation parameter;
- no cosmology;
- no distance-redshift closure;
- no kappa, shear, HST morphology, or observer lensing targets;
- no legacy `0.18`;
- no Quantum Engine execution;
- no Planck-scale input;
- no random noise injection;
- no production-code modification merely to make a candidate pass;
- no run-directory creation unless independently required by an existing tested function; prefer stdout-only audit behavior;
- preserve historical untracked `runs/...` directories and stash state.

## Expected classifications

The implementation should be able to return one of the following without treating a non-match as execution failure:

### `EXISTING_NATIVE_ACCUMULATION_CANDIDATE_FOUND`

At least one existing 3D native candidate reproduces the frozen spatial accumulation fingerprint within predeclared tolerances while retaining source linearity and low-loading additivity.

This does not close amplitude or microscopic interpretation.

### `EXISTING_NATIVE_ACCUMULATION_PARTIAL_MATCH_ONLY`

One or more candidates reproduce some but not all required accumulation behavior.

Report exactly which fingerprints pass and fail.

### `NATIVE_ACCUMULATION_OPERATOR_NOT_FOUND`

No existing downstream native operation tested produces the required long-range spatial accumulation structure.

This is a legitimate derivation boundary. It licenses a later development step to derive or explicitly speculate about the minimal missing accumulation law, but does not itself justify inventing one.

## If no existing operator is found

Only after the repository inventory and native candidate tests are exhausted should the next development phase consider candidate mathematical accumulation mechanisms.

At that later stage, the question should be phrased:

> What minimal medium-response operator maps a localized, approximately linear native loading state into the independently required long-range accumulated-response fingerprint?

Possible mathematical classes may then be studied as hypotheses, but they must be labelled speculative unless independently derived from deeper PBUF constitutive physics.

A familiar inverse spatial operator or Green-function structure may be among the hypotheses because such operators can produce long-range responses from localized sources, but it must not be imported merely because standard gravity uses a similar form.

## Success criterion for this development step

The immediate goal is not to derive gravity.

The goal is to reduce the open native bridge from:

`matter -> unknown medium response -> emergent macroscopic behavior`

to either:

`matter -> c_state/local loading -> identified existing accumulated-response candidate -> later physical mapping`

or, if no candidate exists:

`matter -> c_state/local loading -> explicitly identified missing accumulation operator -> later speculative/constitutive derivation`.

That is the scientific boundary this dev doc is intended to establish.
