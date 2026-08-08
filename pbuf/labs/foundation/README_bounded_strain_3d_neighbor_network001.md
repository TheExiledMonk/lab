# Bounded-Strain 3D Neighbor Network 001

## Purpose

Test the frozen bounded-strain nearest-neighbor constitutive law in a genuine 3D discrete network, without using the spherical equilibrium shortcut from the radial constitutive labs.

Each nearest-neighbor bond uses:

```text
W(e)=-(K e_max^2/2) ln(1-(e/e_max)^2)
sigma(e)=K e/(1-(e/e_max)^2)
```

with frozen `K0=1` and `epsilon_max=1`.

The 3D node equilibrium is the discrete divergence of the six nearest-neighbor bond stresses balancing the source. The nonlinear network is solved self-consistently by Picard iteration with matrix-free conjugate-gradient inner solves.

## Frozen fingerprint

The lab measures the same structural requirements:

1. density exponent `+1`
2. surface mass exponent at fixed radius `+1`
3. surface radius exponent at fixed load `-1`
4. far mass exponent at fixed probe `+1`
5. far radial exponent at fixed source `-1`
6. weak-regime additivity

It also applies one stronger load and reports the maximum bond-strain fraction to verify the finite strain barrier remains respected.

## Guardrails

- no spherical equilibrium relation `r^2 sigma(e)=constant`
- no inserted `1/r` response or inverse-square law
- no G
- no macroscopic amplitude calibration
- no native rescaling
- no parameter fitting/tuning
- no Rmax
- no cosmology
- no lensing target
- no legacy 0.18
- no Quantum Engine
- no Planck input

## Valid scientific outcomes

```text
BOUNDED_STRAIN_3D_NEIGHBOR_NETWORK_STRUCTURE_SUPPORTED
BOUNDED_STRAIN_3D_NEIGHBOR_NETWORK_PARTIAL_SUPPORT
BOUNDED_STRAIN_3D_NEIGHBOR_NETWORK_NOT_SUPPORTED
```

A partial or null result is valid and is not permission to modify or tune the model.

## Runner command

```bash
PYTHONPATH=. python pbuf/labs/foundation/bounded_strain_3d_neighbor_network001.py
```

Return the exact raw stdout/stderr, exit code, branch and HEAD, `git status --short`, and `git stash list`.

Do not modify, repair, tune, rescale, reinterpret, or merge anything.
