# PBUF homogeneous elastic-energy closure audit 001

## Purpose

Test whether current PBUF already contains an independently justified, Rmax-free constitutive bridge from the retained thermal quantities to a homogeneous elastic energy density and then to `Omega_sigma(a)`.

Retained inputs:

- `alpha_T(a)`
- `epsilon0_T(a)`
- `k_max(a)=epsilon0_T(a)-alpha_T(a)`

`Rmax` remains retired from active PBUF.

## Required closure chain

A physical homogeneous elastic background requires all of the following:

```text
homogeneous medium state chi(a)
    -> normalized constitutive energy functional S_sigma / S_med
    -> physical rho_sigma(a) with absolute energy normalization
    -> Omega_sigma(a)
    -> E(a)
    -> H(a)
```

The audit checks whether these links already exist in the accepted foundation sources.

## Forbidden shortcuts

The lab must not:

- set `Omega_sigma(a)=alpha_T(a)`;
- set `Omega_sigma(a)=k_max(a)`;
- multiply thermal quantities merely to obtain a useful curve;
- normalize a candidate to flatness to hide missing amplitude;
- restore the retired `Rmax -> decay -> S -> Omega_sigma_raw` construction;
- introduce a replacement activation/amplitude parameter;
- fit cosmological or lensing observations.

Textual mentions of `Omega_sigma`, `S_sigma`, energy density, or medium actions are evidence for review only. They are not physical closure unless the source supplies the required normalized constitutive chain.

## Expected safe result

If the repository still lacks the normalized covariant medium action / constitutive energy functional, return:

```text
HOMOGENEOUS_ELASTIC_ENERGY_CLOSURE_NOT_YET_DERIVED
```

That is a valid result. Do not repair it by adding assumptions.

## Guardrails

- `Rmax` retired and numerically unused;
- no replacement free parameter;
- no assumed future reversal;
- no `alpha_T` or `k_max` promotion to `Omega_sigma`;
- no flatness normalization of an unclosed candidate;
- no observed cluster redshift promoted to expansion redshift;
- no kappa, shear, lensing morphology, or lensing amplitude;
- no LambdaCDM distance substitution;
- no measured-G backsolve;
- no legacy `0.18`;
- no Quantum Engine execution;
- no Planck-scale input;
- stdout only; no run directory.

## Run

From repository root on the exact PR branch:

```bash
PYTHONPATH=. python pbuf/labs/foundation/homogeneous_elastic_energy_closure_audit001.py
```

## Runner contract

The runner is an executor only. Do not modify the lab, source files, thermal cache, constants, classifications, or repository data.

Return:

1. current HEAD SHA and branch name;
2. exit code;
3. complete raw stdout and stderr;
4. `git status --short` after the run;
5. confirmation that the preservation stash was untouched.

Preserve all historical untracked `runs/...` directories and the existing preservation stash.
