# PBUF FOUNDATION — Bounded-Strain Finite-Domain Convergence 001

## Purpose

Test the strongest remaining numerical alternative from the bounded-strain constitutive lab without changing the physics: determine whether the apparent far-field steepening is caused by the finite outer boundary.

The constitutive law is frozen:

`W(e)=-(K e_max^2/2) ln(1-(e/e_max)^2)`

`sigma(e)=K e/(1-(e/e_max)^2)`

with `K0=1`, `epsilon_max=1`, and load `8`.

Only the outer boundary is moved:

`R_B = 128, 256, 512, 1024`.

The same fixed physical probe radii are measured at every boundary.

## Question

Does the frozen bounded-strain response converge toward radial exponent `-1` at fixed probes as the artificial zero-response boundary is moved outward?

## Guardrails

- no constitutive-law change;
- no K0 change;
- no epsilon_max change;
- no load change;
- no inserted `1/r` response;
- no fitted radial law;
- no rescaling or tuning;
- no G;
- no Rmax;
- no cosmology;
- no lensing target;
- no Quantum Engine;
- no Planck input;
- stdout only;
- preserve historical `runs/...` and stash state.

## Valid statuses

- `BOUNDED_STRAIN_FINITE_DOMAIN_CONVERGENCE_SUPPORTED`
- `BOUNDED_STRAIN_FINITE_DOMAIN_CONVERGENCE_PARTIAL_SUPPORT`
- `BOUNDED_STRAIN_FINITE_DOMAIN_CONVERGENCE_NOT_SUPPORTED`

A partial or null outcome is a valid result.

## Runner

```bash
PYTHONPATH=. python pbuf/labs/foundation/bounded_strain_finite_domain_convergence001.py
```

Return branch/HEAD, exit code, complete raw stdout/stderr, `git status --short`, and `git stash list`. Do not modify, repair, tune, rescale, reinterpret, or merge anything.
