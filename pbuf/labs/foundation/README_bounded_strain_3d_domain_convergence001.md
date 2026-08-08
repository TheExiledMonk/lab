# Bounded-Strain 3D Domain Convergence 001

## Purpose

Test whether the remaining far-field slope miss in the full 3D bounded-strain nearest-neighbor network is caused by the finite cubic zero-response boundary.

The physics is frozen. Only the cubic domain size changes.

## Frozen model

- `K0 = 1.0` as the same dimensionless structural normalization used by the prior lab; it is not scanned or fitted.
- `epsilon_max = 1.0`.
- Constitutive stress: `sigma(e)=K e/(1-(e/e_max)^2)`.
- Six Cartesian nearest-neighbor bonds.
- Network equilibrium: discrete divergence of bounded-strain bond stresses equals the source.
- Zero-response boundary on the outer cubic faces.
- Source radius `3.5`.
- Integrated source load `2.0`.
- Fixed probe radii `6,7,8,9,10`.

No spherical equilibrium relation is used.

## Domain ladder

`N = 49, 65, 81, 97`

The physical source and probe positions remain unchanged while the boundary recedes.

## Primary question

Does the measured far-radius exponent, and the local exponent at fixed `r=8`, move monotonically toward `-1` as the domain grows?

## Guardrails

No G, no macroscopic amplitude calibration, no native rescaling, no fitted or tuned stiffness, no inserted `1/r` response, no spherical-equilibrium shortcut, no Rmax, no cosmology, no lensing target, no Quantum Engine, no Planck input.

## Run

```bash
PYTHONPATH=. python pbuf/labs/foundation/bounded_strain_3d_domain_convergence001.py
```

## Valid scientific statuses

- `BOUNDED_STRAIN_3D_DOMAIN_CONVERGENCE_SUPPORTED`
- `BOUNDED_STRAIN_3D_DOMAIN_CONVERGENCE_PARTIAL_SUPPORT`
- `BOUNDED_STRAIN_3D_DOMAIN_CONVERGENCE_NOT_SUPPORTED`

A partial or null result is scientifically valid. Do not modify or tune the model after seeing the result.

## Runner contract

Return branch and HEAD, exit code, complete raw stdout/stderr, `git status --short`, and `git stash list`. Do not modify, repair, tune, rescale, reinterpret, or merge anything. Preserve existing stash and historical untracked `runs/...` directories.
