# A8 Fast/Slow Eigenmode Decomposition Audit 001

## Question

Does the historical reduction from local `c_state` bond scale to A8 pair-transfer scale arise from the actual fast/slow channel structure rather than from an arbitrary scalar partition?

## Exact modal rewrite

The frozen terminal A8 channels are rewritten without approximation:

```text
c = (u_fast + u_slow)/2
d = (u_fast - u_slow)/2
u_fast = c + d
u_slow = c - d
```

The existing pair-transfer law

```text
A_ij = coef_fast * Delta u_fast + coef_slow * Delta u_slow
```

therefore becomes exactly

```text
A_ij = (coef_fast + coef_slow) * Delta c
     + (coef_fast - coef_slow) * Delta d
```

For the frozen constants this gives common-mode coefficient `0.033` and difference-mode coefficient `0.027`.

## Diagnostics

For each of the five canonical local benchmark clusters, the audit reports:

- RMS difference-mode bond relative to common-mode bond;
- correlation between difference- and common-mode bonds;
- descriptive projection of `Delta d` onto `Delta c`;
- exact common- and difference-mode contributions to pair transfer;
- resulting pair-transfer RMS relative to common-mode bond RMS;
- exact modal-formula reproduction error.

The projection coefficient is descriptive only. It is never used to construct a native channel split or tune a model.

## Guardrails

- historical unit-loading A8 dynamics unchanged;
- five canonical local benchmark FITS only;
- no network access;
- no observed lensing values;
- no native fast/slow reconstruction;
- no inferred ratio applied;
- no replacement strength scalar;
- no normalization/rescaling;
- no fitting/tuning;
- no G3D candidate;
- no GR/Weyl/LCDM/Rmax/QE/Planck input.

## Run

```bash
PYTHONPATH=. python pbuf/labs/foundation/a8_fast_slow_eigenmode_decomposition_audit001.py
```

## Status values

- `A8_FAST_SLOW_EIGENMODE_DECOMPOSITION_AUDIT_EXECUTED`
- `A8_FAST_SLOW_EIGENMODE_DECOMPOSITION_AUDIT_PARTIAL_EXECUTION`
- `A8_FAST_SLOW_EIGENMODE_DECOMPOSITION_AUDIT_NOT_ESTABLISHED`
