# PBUF FOUNDATION — Bounded Strain Neighbor Constitutive 001

## Purpose

Test the working constitutive idea that neighbor separation becomes progressively harder as strain increases, with a possible finite limiting strain that requires unbounded energy to approach.

This is a direct implementation, not a paper argument. The model is intentionally structural and dimensionless.

## Constitutive hypothesis

The lab uses the fixed barrier energy

`W(e)=-(K e_max^2/2) ln(1-(e/e_max)^2)`

with restoring stress

`sigma(e)=K e/(1-(e/e_max)^2)`.

As `|e| -> e_max`, stored energy, restoring stress, and tangent stiffness diverge.

For weak strain the law reduces to the previous linear neighbor medium:

`sigma ~ K e`.

## Static accumulation test

For a localized spherical source, the lab solves the spherically symmetric equilibrium condition

`d/dr[r^2 sigma(e(r))]=0`

outside the source, numerically inverts the nonlinear constitutive stress, and integrates strain to obtain the accumulated response.

No `1/r` response is inserted. No response exponent is fitted into the model.

The lab reports:

- constitutive energy/stress/stiffness growth near the strain barrier;
- weak-strain recovery of the linear constitutive regime;
- weak-load source and far-response scaling;
- strong-load approach toward the strain limit;
- local radial response exponent as strain decreases outward.

## Guardrails

- gravity remains emergent and is not a native variable;
- no G;
- no macroscopic amplitude calibration;
- no native rescaling;
- no fitting or tuning;
- no inserted `1/r` response;
- no Rmax;
- no cosmology;
- no lensing target;
- no legacy `0.18`;
- no Quantum Engine;
- no Planck input;
- no EM-origin claim;
- stdout only.

## Valid scientific statuses

- `BOUNDED_STRAIN_NEIGHBOR_CONSTITUTIVE_STRUCTURE_SUPPORTED`
- `BOUNDED_STRAIN_NEIGHBOR_CONSTITUTIVE_PARTIAL_SUPPORT`
- `BOUNDED_STRAIN_NEIGHBOR_CONSTITUTIVE_NOT_SUPPORTED`

Partial/null support is a valid scientific outcome and is not permission to alter parameters after execution.

## Runner

Run exactly:

```bash
PYTHONPATH=. python pbuf/labs/foundation/bounded_strain_neighbor_constitutive001.py
```

Return branch/HEAD, exit code, complete raw stdout/stderr, `git status --short`, and `git stash list`. Do not modify, repair, tune, rescale, reinterpret, or merge anything.
