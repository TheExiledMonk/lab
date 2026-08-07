# PBUF FOUNDATION — REAL-CANDIDATE COVARIANCE LOCALIZATION LAB 001

**Lab ID**: PBUF-FOUNDATION-REAL-CANDIDATE-COVARIANCE-LOCALIZATION-001
**Head SHA**: `b5a5e12aa4bfb53db85bff1c885440681c0b9301`
**Branch**: `main`
**Cluster**: MACS0416
**Candidate**: PL1_PM1_PS2
**Nz / profile / stencil / boundary / strength / seed**: 9 / gaussian / N6 / reflective / 0.18 / 12345
**Native shape**: (9, 64, 64)
**Conventions version**: 1.1.0-correction001
**First-failure threshold**: 1e-08
**Duration**: 17.59 s

## Checkpoint covariance table

| RC | rho3d | u_slow | u_fast | c_state | eL | PT | Aij-oriented | Rij-oriented | endpoint | interface | first failure |
|---|---|---|---|---|---|---|---|---|---|---|---|
| RC0 | 0.000e+00 | 0.000e+00 | 0.000e+00 | 0.000e+00 | 0.000e+00 | 0.000e+00 | 0.000e+00 | 0.000e+00 | 0.000e+00 | 0.000e+00 | none |
| RC1 | 0.000e+00 | 5.557e-02 | 5.569e-02 | 5.419e-02 | 4.814e-01 | 3.050e-01 | 1.222e-01 | 1.272e+00 | 3.212e-01 | 3.464e-01 | u_slow |
| RC2 | 0.000e+00 | 5.580e-02 | 5.635e-02 | 5.460e-02 | 4.866e-01 | 3.055e-01 | 1.222e-01 | 1.138e+00 | 3.130e-01 | 3.349e-01 | u_slow |
| RC3 | 0.000e+00 | 5.553e-02 | 5.471e-02 | 5.374e-02 | 4.921e-01 | 3.052e-01 | 1.226e-01 | 1.130e+00 | 3.274e-01 | 3.485e-01 | u_slow |
| RC4 | 0.000e+00 | 5.656e-02 | 5.704e-02 | 5.536e-02 | 4.986e-01 | 3.051e-01 | 1.244e-01 | 1.372e+00 | 1.219e+00 | 3.465e-01 | u_slow |
| RC5 | 0.000e+00 | 5.580e-02 | 5.564e-02 | 5.429e-02 | 4.919e-01 | 3.045e-01 | 1.226e-01 | 1.253e+00 | 1.165e+00 | 3.425e-01 | u_slow |
| RC6 | 0.000e+00 | 5.638e-02 | 5.698e-02 | 5.519e-02 | 4.986e-01 | 3.082e-01 | 1.240e-01 | 1.303e+00 | 1.171e+00 | 3.492e-01 | u_slow |

## Wrong-control summary

| RC | WC1 scalar-only (endpoint) | WC2 ignore endpoint swap (pair-amp) | WC3 omit antisymmetric sign (pair-amp) |
|---|---|---|---|
| RC0 | 0.000e+00 | 0.000e+00 | 0.000e+00 |
| RC1 | 1.128e+00 | 9.270e-01 | 1.222e-01 |
| RC2 | 1.129e+00 | 1.192e+00 | 1.222e-01 |
| RC3 | 1.117e+00 | 1.195e+00 | 1.226e-01 |
| RC4 | 1.319e+00 | 1.196e+00 | 1.361e+00 |
| RC5 | 1.320e+00 | 1.191e+00 | 1.032e+00 |
| RC6 | 1.286e+00 | 9.970e-01 | 1.046e+00 |

## First-failure summary

| RC | first_failure_checkpoint | relative_error | previous_checkpoint | previous_relative_error |
|---|---|---|---|---|
| RC0 | none | n/a | interface | 0.000e+00 |
| RC1 | u_slow | 5.557e-02 | rho_3d | 0.000e+00 |
| RC2 | u_slow | 5.580e-02 | rho_3d | 0.000e+00 |
| RC3 | u_slow | 5.553e-02 | rho_3d | 0.000e+00 |
| RC4 | u_slow | 5.656e-02 | rho_3d | 0.000e+00 |
| RC5 | u_slow | 5.580e-02 | rho_3d | 0.000e+00 |
| RC6 | u_slow | 5.638e-02 | rho_3d | 0.000e+00 |

## Interpretation

**Outcome A — A8/T1 evolved-state covariance failure**

The first failing checkpoint for every non-trivial RC is u_slow / u_fast / c_state. The T1 evolution applied to the spatially-permuted rho produces a state that does not round-trip through the inverse spatial transform to machine precision. A8/T1 transformed evolution is not covariant under the seven-rotation set. Earlier checkpoints (rho_3d) and downstream checkpoints (eL, PT, pair amplitudes, pair responses, endpoint, interface) all inherit this dominant propagation.

## Per-pair-slot summary

| RC | n_pairs | n_endpoint_swap_required | n_orientation_sign_required | max_abs_diff | mean_abs_diff |
|---|---|---|---|---|---|
| RC0 | 105344 | 0 | 0 | 0.000e+00 | 0.000e+00 |
| RC1 | 105344 | 0 | 0 | 1.466e-05 | 2.426e-06 |
| RC2 | 105344 | 0 | 0 | 1.464e-05 | 2.435e-06 |
| RC3 | 105344 | 0 | 0 | 1.524e-05 | 2.436e-06 |
| RC4 | 105344 | 32768 | 32768 | 1.524e-05 | 2.469e-06 |
| RC5 | 105344 | 36288 | 36288 | 1.386e-05 | 2.444e-06 |
| RC6 | 105344 | 36288 | 36288 | 1.544e-05 | 2.470e-06 |

## Pair-direction transform table

| transform | source_direction | mapped_signed_direction | canonical_direction | endpoint_swap | orientation_sign |
|---|---|---|---|---|---|
| RC0 | xp | xp | xp | False | 1 |
| RC0 | yp | yp | yp | False | 1 |
| RC0 | zp | zp | zp | False | 1 |
| RC1 | xp | yp | yp | False | 1 |
| RC1 | yp | xp | xp | False | 1 |
| RC1 | zp | zp | zp | False | 1 |
| RC2 | xp | zp | zp | False | 1 |
| RC2 | yp | yp | yp | False | 1 |
| RC2 | zp | xp | xp | False | 1 |
| RC3 | xp | xp | xp | False | 1 |
| RC3 | yp | zp | zp | False | 1 |
| RC3 | zp | yp | yp | False | 1 |
| RC4 | xp | xp | xp | False | 1 |
| RC4 | yp | zp | zp | False | 1 |
| RC4 | zp | ym | yp | True | -1 |
| RC5 | xp | zm | zp | True | -1 |
| RC5 | yp | yp | yp | False | 1 |
| RC5 | zp | xp | xp | False | 1 |
| RC6 | xp | yp | yp | False | 1 |
| RC6 | yp | xm | xp | True | -1 |
| RC6 | zp | zp | zp | False | 1 |

## Hard rules

- NO SOURCE CHANGES
- NO TOLERANCE CHANGES
- NO SYNTHETIC SUBSTITUTE
- NO RAY TRACING
- NO JACOBIAN
- NO OBSERVATIONAL FITTING
- NO FIXING DURING EXECUTION
