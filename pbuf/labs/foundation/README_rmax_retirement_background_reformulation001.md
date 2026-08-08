# Rmax retirement background reformulation 001

## Purpose

Retire `Rmax` from active PBUF background physics and determine what remains of the historical V11 background when no future-reversal scale is allowed as an input.

`Rmax` is treated as a historical placeholder for a hypothetical future expansion reversal. It is not assigned a value, recovered from old code, fitted, inherited, or replaced by another free parameter.

## Core rule

```text
Rmax = RETIRED_FROM_ACTIVE_PBUF
```

A future turnaround may only re-enter as an output of completed dynamics, for example if an independently derived background eventually yields

```text
H(a_turn) = 0
```

No turnaround is assumed in advance.

## Historical relations under review

The audit classifies each relation separately:

```text
alpha_T(a), epsilon0_T(a)
  -> k_max(a)=epsilon0_T(a)-alpha_T(a)
  -> historical decay(a)=exp(-a/Rmax)
  -> historical S(a)
  -> historical Omega_sigma_raw(a)
  -> historical flat-today rescale
  -> Omega_sigma(a)
  -> E(a)
  -> H(a)
```

Expected classification logic:

- thermal `alpha_T(a)` and `epsilon0_T(a)`: retain if independently available;
- `k_max(a)`: retain because it is independent of `Rmax`;
- historical `decay(a)`: retire with `Rmax`;
- historical `S(a)`: retire in its current form because it inherits the retired activation factor;
- historical `Omega_sigma_raw(a)`: retire in its current form;
- historical `sigma_rescale`: retire because it normalises the retired raw construction;
- `Omega_sigma(a)`: concept may survive but requires independent Rmax-free derivation;
- `E(a)` and `H(a)`: structure-only until the elastic background is independently closed.

## Repository inventory

The lab also inventories:

1. remaining repository occurrences of `Rmax` for retirement planning;
2. textual `Omega_sigma` expressions that do not contain `Rmax` on the same line.

Textual candidates are never promoted automatically to physical laws. They are reported only for later source review.

## Required outcome

If no independently justified Rmax-free constitutive law for `Omega_sigma(a)` is already present, the correct result is:

```text
PRESENT_ERA_BACKGROUND_LAW_NOT_YET_DERIVED
```

This is not a failure. It identifies the exact next physics problem without inserting a replacement scale.

## Guardrails

- no numerical `Rmax`;
- no replacement free parameter;
- no assumption of future reversal;
- no cluster observed redshift promoted to expansion redshift;
- no kappa, shear, lensing target, morphology, or amplitude;
- no LambdaCDM distance substitution;
- no measured-`G` backsolve;
- no legacy `0.18`;
- no Quantum Engine execution;
- no Planck-scale input;
- no fabricated distance values;
- no production-code changes;
- stdout only; no run directory.

## Run

From repository root on the exact PR branch:

```bash
PYTHONPATH=. python pbuf/labs/foundation/rmax_retirement_background_reformulation001.py
```

## Runner contract

The runner is an executor only. Return:

1. current HEAD SHA and branch name;
2. exit code;
3. complete raw stdout and stderr;
4. `git status --short` after execution;
5. confirmation that the preservation stash is unchanged.

Do not modify equations, constants, thermal data, classifications, repository history, or historical untracked `runs/...` directories to change the result.
