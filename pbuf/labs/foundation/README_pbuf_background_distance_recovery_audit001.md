# PBUF background distance recovery audit 001

## Purpose

Recover the historical V11 PBUF background-expansion chain from current repository evidence before any cluster redshift is converted into a distance.

This is a **fact-finding recovery audit**. It does not fit cosmological parameters and it does not import LambdaCDM distances.

## Historical chain under audit

```text
alpha_T(a), epsilon0_T(a)
  -> k_max(a)
  -> decay(a)=exp(-a/Rmax)
  -> S(a)
  -> Omega_sigma_raw(a)
  -> flat-today rescale
  -> Omega_sigma(a)
  -> E(a)
  -> H(a)
  -> future radial-distance integral
```

The equations are taken from the existing V11 trace/provenance already present in the repository. The current thermal cache at `pbuf/data/quantum/thermal_table_cache.json` is inventoried and checked for numeric readability.

## What the audit must not do

The audit must not choose a historical or best-fit value merely because a symbol/value appears somewhere in the repository. In particular, it does not silently choose values for:

- `Rmax`
- `H0`
- `Omega_r0`
- `alpha_resolved`
- `BARYON_FRACTION`

String/numeric candidates are reported as provenance only until a later audit identifies an authoritative current source.

The thermal cache field `alpha_qm` is explicitly **not** promoted to `alpha_resolved` by assumption.

## Redshift boundary

This lab is background-only. It does not use the five observed cluster redshifts as expansion redshifts.

```text
z_observed != z_expansion by assumption
```

The recovered background may only be applied to clusters after the redshift decomposition is separately audited.

## Algebra control

A deterministic synthetic wiring test verifies the equation implementation:

```text
k_max = epsilon0_T - alpha_T
decay = exp(-a/Rmax)
S = 1 - (1-k_max) decay
Omega_sigma_raw = alpha_T (1-decay) S
sigma_rescale = Omega_sigma_target / Omega_sigma_raw(1)
Omega_sigma = sigma_rescale Omega_sigma_raw
E^2 = Omega_m0 a^-3 + Omega_r0 a^-4 + Omega_sigma
H = H0 E
```

The synthetic values have no physical meaning. They are used only to verify algebraic consistency.

## Guardrails

- no kappa or shear
- no lensing morphology or lensing amplitude
- no LambdaCDM distance substitution
- no fitted `Rmax`, `H0`, `Omega_r0`, `alpha`, or redshift correction
- no measured `G` backsolve
- no legacy `0.18`
- no Quantum Engine execution
- no Planck-scale input
- no production-code changes
- stdout only; no run directory

## Run

From the repository root on the exact PR branch:

```bash
PYTHONPATH=. python pbuf/labs/foundation/pbuf_background_distance_recovery_audit001.py
```

## Runner contract

The runner is an **executor only**.

Do not modify the lab, production modules, thermal cache, constants, thresholds, or repository data before or after execution.

Do not repair a failure. Return the raw failure unchanged.

Return exactly:

1. current HEAD SHA and branch name;
2. process exit code;
3. complete raw stdout and stderr;
4. `git status --short` after the run;
5. confirmation that the preservation stash was not altered.

Do not delete, move, clean, modify, or commit historical untracked `runs/...` directories.

Do not pop, apply, drop, rewrite, or otherwise alter the preservation stash.
