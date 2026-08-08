# PBUF FOUNDATION — Pretensioned Neighbor Medium Response 001

## Purpose

Test a minimal locally coupled, pre-tensioned medium as a candidate mechanism for the missing native accumulation behavior.

The model is intentionally simple and is derived from a nearest-neighbor coupling energy rather than from an inserted radial response law:

```text
E = (T/2) sum_<ij> (u_i - u_j)^2 - sum_i S_i u_i
```

Stationarity gives:

```text
T L[u]/DX^2 = S
```

The same local coupling has a source-free dynamical sector:

```text
mu u_tt = -T L[u]/DX^2
```

The lab asks whether this one minimal structure can produce both:

1. the required static long-range accumulated-response fingerprints; and
2. a wave-supporting dynamical sector whose speed follows the coupling/inertia scaling implied by the same equation.

No claim is made that the dynamical wave sector is electromagnetism. The point is only to test whether the same substrate-style local coupling can support both equilibrium redistribution and propagating disturbances.

## What is not inserted

- no `1/r` Green function;
- no inverse-square law;
- no G;
- no macroscopic amplitude benchmark;
- no native amplitude rescaling;
- no fitting or tuning to the target exponents;
- no Rmax;
- no cosmology;
- no lensing target;
- no legacy `0.18`;
- no Quantum Engine;
- no Planck input.

## Static tests

The discrete equilibrium is solved numerically with a matrix-free conjugate-gradient solve of the nearest-neighbor operator.

The runner measures:

- density exponent at fixed source radius;
- surface-response mass exponent at fixed radius;
- surface-response radius exponent at fixed integrated native load;
- far-response mass exponent at fixed probe radius;
- far-response radial exponent for a fixed source;
- two-source superposition residual.

The independently frozen shape targets remain:

```text
rho exponent = +1
surface mass exponent = +1
surface radius at fixed load = -1
far mass exponent = +1
far radius exponent = -1
additivity residual ~ 0
```

## Dynamic test

For a plane-wave mode of the same nearest-neighbor equation:

```text
omega^2 = (4 T / (mu DX^2)) sum_a sin^2(k_a DX/2)
```

The lab evaluates the long-wavelength branch and measures whether:

```text
wave speed ~ T^(+1/2)
wave speed ~ mu^(-1/2)
```

These exponents come from the candidate medium equation itself; they are not fitted to an electromagnetic target.

## Valid statuses

- `PRETENSIONED_NEIGHBOR_MEDIUM_STATIC_DYNAMIC_STRUCTURE_FOUND`
- `PRETENSIONED_NEIGHBOR_MEDIUM_PARTIAL_STRUCTURE_ONLY`
- `PRETENSIONED_NEIGHBOR_MEDIUM_STRUCTURE_NOT_SUPPORTED`

A partial or null scientific outcome is valid and must be returned unchanged.

## Runner

Run exactly:

```bash
PYTHONPATH=. python pbuf/labs/foundation/pretensioned_neighbor_medium_response001.py
```

Return HEAD/branch, exit code, complete raw stdout/stderr, `git status --short`, and `git stash list`.

Do not modify, repair, tune, rescale, reinterpret, or merge anything.
