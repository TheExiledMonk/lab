# PBUF FOUNDATION — Native Medium Response Accumulation 001

## Purpose

Execute the development target defined by `DEV_DOC_native_medium_response_accumulation001.md` directly in code.

Native causal target only:

`local medium loading -> long-range accumulated medium response`

Gravity remains emergent and is not used as a native construction variable or field label.

## Existing candidates tested

The lab tests existing repository operations only:

- `c_state`;
- `laplacian(c_state)` via `native_field_curvature_dimension_audit001._laplacian`;
- signed M10 interface vector via `m10_coverage_25pct_science001._interface_vector`;
- `divergence(M10)` via `native_field_curvature_dimension_audit001._divergence`.

No new `1/r` kernel, inverse operator, Green function, or propagation law is introduced.

## Frozen structural tests

Each candidate is tested against:

- local density exponent `+1`;
- surface mass exponent at fixed radius `+1`;
- surface radius exponent at fixed integrated native source `-1`;
- far mass exponent at fixed probe radius `+1`;
- far radius exponent for a fixed source `-1`;
- signed-component two-source additivity.

These are shape/scaling fingerprints only. No amplitude calibration is performed.

## Valid scientific classifications

- `EXISTING_NATIVE_ACCUMULATION_CANDIDATE_FOUND`
- `EXISTING_NATIVE_ACCUMULATION_PARTIAL_MATCH_ONLY`
- `NATIVE_ACCUMULATION_OPERATOR_NOT_FOUND`

A partial or null scientific result is not an execution failure.

## Guardrails

- gravity is not a native variable;
- no G;
- no macroscopic amplitude benchmark;
- no native amplitude rescaling;
- no fitting or tuning;
- no new accumulation kernel;
- no Rmax or replacement parameter;
- no cosmology;
- no lensing target;
- no legacy `0.18`;
- no Quantum Engine;
- no Planck input;
- stdout only;
- historical untracked `runs/...` and stash state must remain untouched.

## Runner

```bash
PYTHONPATH=. python pbuf/labs/foundation/native_medium_response_accumulation001.py
```

Return current branch/HEAD, exit code, complete raw stdout and stderr, `git status --short`, and `git stash list`. Do not repair, modify, rescale, tune, or reinterpret the result.
