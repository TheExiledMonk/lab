# PBUF FOUNDATION — Native Response Fingerprint Match 001

## Purpose

Compare existing native PBUF response variables against the frozen effective weak-field fingerprint established by `known_source_inverse_response_fingerprint001`.

The purpose is not to invent a new spacetime law. It is to ask whether the current native implementation already contains variables with the correct structural behavior.

## Frozen effective fingerprint

The comparison target is fixed before the native run:

- local loading: `response ~ rho`;
- surface response at fixed radius: `response ~ M`;
- surface response at fixed mass: `response ~ R^-1`;
- far response at fixed probe radius: `response ~ M`;
- far response for a fixed source: `response ~ r^-1`;
- weak-field response: additive to first order.

These are scaling constraints only. No macroscopic amplitude is imported.

## Native candidates

The lab compares:

1. `c_state` — using center/local and shell/far samples;
2. M10 interface vector — using vector magnitude for scaling and signed vector components for superposition.

## Independent tests

### Density ladder

Fixed source radius, varying source density. Tests whether local/native response is linear in loading.

### Mass ladder at fixed radius

The integrated native source is varied while radius is held fixed. Tests the expected mass exponent of +1.

### Radius ladder at fixed native mass

Source radius varies while integrated native source is held constant by construction. Tests the required surface exponent of -1.

### Far-radius ladder

A fixed source is probed along exact grid locations outside the source. Tests whether the native field falls approximately as `r^-1`.

### Far-mass ladder

At a fixed far probe point, source mass is varied. Tests the expected +1 mass exponent.

### Two-source superposition

Two separated synthetic sources are run individually and together. `c_state` is compared directly; M10 is compared through signed vector components rather than magnitudes.

## Acceptance

Exponent comparisons use a predeclared ±0.20 structural window. Superposition uses a relative RMS tolerance of `1e-10`.

A full match means only that a native variable reproduces the effective weak-field scaling structure. It does **not** establish:

- SI amplitude;
- microscopic substrate identity;
- a derivation of `G`;
- a cosmological closure;
- a lensing normalization.

## Guardrails

- no `G` in the native calculation;
- no macroscopic response amplitude imported;
- no amplitude rescaling;
- no fitting or tuning;
- no microscopic substrate claim;
- no `Rmax` or replacement parameter;
- no cosmology;
- no lensing target, kappa, shear, morphology, or observer target;
- no historical `0.18`;
- no Quantum Engine or Planck input;
- stdout only; no new run directory.

## Possible statuses

`NATIVE_EFFECTIVE_RESPONSE_FINGERPRINT_MATCH_FOUND`

or

`NATIVE_EFFECTIVE_RESPONSE_FINGERPRINT_PARTIAL_MATCH_ONLY`

A partial result is scientifically useful: it identifies exactly which accumulation/propagation property is missing from the current native stage.

## Runner

```bash
PYTHONPATH=. python pbuf/labs/foundation/native_response_fingerprint_match001.py
```

Return branch/HEAD, exit code, complete raw stdout and stderr, `git status --short`, and preservation-stash confirmation. Do not modify, repair, rescale, tune, or reinterpret a partial/failing scientific result.
