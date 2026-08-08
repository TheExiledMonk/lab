# Homogeneous medium state derivation 001

## Purpose

Determine whether current PBUF already defines a physically justified homogeneous medium state variable `chi(a)` after retirement of `Rmax`.

This is a derivation-boundary audit. It does not choose a state variable because it is convenient or numerically well behaved.

## Retained Rmax-free inputs

- `alpha_T(a)`
- `epsilon0_T(a)`
- `k_max(a)=epsilon0_T(a)-alpha_T(a)`

## Candidate rule

A candidate may be selected as the homogeneous state only if existing PBUF sources establish all of the following together:

1. direct identity / definition;
2. physical role as the medium state rather than merely a coefficient, amplitude, response, or bound;
3. dimensional character / normalization;
4. an evolution or state equation.

Numeric similarity, monotonicity, convenient algebra, or cosmological usefulness are not evidence of identity.

## Candidates audited

- `chi(a)=alpha_T(a)`
- `chi(a)=epsilon0_T(a)`
- `chi(a)=k_max(a)`
- a separate `chi(a)` governed by retained thermal quantities

The audit does not promote any candidate automatically.

## Expected safe outcome

If no candidate satisfies all four requirements, return:

`HOMOGENEOUS_MEDIUM_STATE_NOT_YET_DERIVED`

Then constitutive-energy derivation remains blocked until the medium variable and its dynamics are independently derived or explicitly postulated from deeper constitutive microphysics.

## Guardrails

- `Rmax` remains retired;
- no replacement free parameter;
- no assumed future reversal;
- no cosmological target used to choose `chi`;
- no `Omega_sigma` construction;
- no flatness normalization;
- no cluster redshift promoted to expansion redshift;
- no kappa, shear, lensing morphology, or lensing amplitude;
- no LambdaCDM distance substitution;
- no measured-`G` backsolve;
- no historical `0.18`;
- no Quantum Engine execution;
- no Planck-scale input;
- stdout only; no run directory.

## Run

From repository root on this branch:

```bash
PYTHONPATH=. python pbuf/labs/foundation/homogeneous_medium_state_derivation001.py
```

## Runner contract

Return:

1. current HEAD SHA and branch name;
2. exit code;
3. complete raw stdout and stderr;
4. `git status --short` after the run;
5. confirmation that the preservation stash was untouched.

Do not alter the lab, source files, thermal cache, constants, or repository data. Do not select a candidate manually. Return an open/blocked result unchanged.

Do not delete, move, clean, modify, or commit historical untracked `runs/...` directories. Do not alter the preservation stash.
