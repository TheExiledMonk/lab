# Constitutive Tangent Propagation Coupling 001

## Purpose

Test whether the already-frozen bounded-strain constitutive law closes the medium-to-propagation coupling by itself, without GR/LCDM potentials, observed lensing amplitudes, or a fitted propagation coefficient.

The supported native chain remains frozen:

`rho -> existing A8 transport -> raw c_state -> bounded-strain accumulated state u`

The bond law is

`sigma(e) = K0 e / (1 - (e/epsilon_max)^2)`

with exact tangent stiffness

`K_t(e) = K0 (1+q)/(1-q)^2`, where `q=(e/epsilon_max)^2`.

For the hypothesis that small propagating disturbances use the same local neighbor mode, `v^2 ~ K_t/mu`. With constant inertial density, the relative propagation index is fixed without an amplitude parameter:

`n/n0 = sqrt(K0/K_t)`.

The directional response tested is

`Delta k_x = integral ds partial_x ln(n/n0)`.

This is deliberately falsifiable. A wrong sign or wrong scaling is a rejection/localization of this propagation hypothesis, not permission to tune the accumulation bridge.

## Frozen checks

- exact constitutive tangent identity;
- zero source -> zero response;
- centered path -> zero transverse response;
- reflection antisymmetry;
- response points toward the source;
- weak response scales approximately as source mass^1;
- impact response magnitude scales approximately as b^-1.

## Guardrails

No `G`; no GR potential decomposition; no LCDM; no observed lensing target; no kappa/shear observation; no fitted propagation amplitude; no extra coupling coefficient; no native rescaling; no tuned `K0`; no inserted `1/r` or `1/b`; no spherical shortcut; no `Rmax`; no cosmology; no Quantum Engine; no Planck input.

## Run

```bash
PYTHONPATH=. python pbuf/labs/foundation/constitutive_tangent_propagation_coupling001.py
```

Return complete stdout/stderr and repository state. Do not repair or tune a partial/null scientific result.
