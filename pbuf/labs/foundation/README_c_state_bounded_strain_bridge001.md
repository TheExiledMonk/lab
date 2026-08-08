# C_STATE Bounded-Strain Bridge 001

## Purpose

Test whether the existing native `c_state` local-loading candidate can drive the bounded-strain 3D nearest-neighbor accumulation network and reproduce the frozen response fingerprint without amplitude fitting.

Native chain:

`rho -> existing A8 transport -> raw c_state`

Accumulation chain:

`raw c_state -> six-neighbor bounded-strain equilibrium -> accumulated response`

The raw `c_state` field is embedded into the larger accumulation grid without normalization or amplitude rescaling.

## Frozen constitutive law

`W(e)=-(K e_max^2/2) ln(1-(e/e_max)^2)`

`sigma(e)=K e/(1-(e/e_max)^2)`

with structural normalization `K0=1`, `epsilon_max=1`.

## Frozen six-check fingerprint

1. density exponent `+1`
2. surface mass exponent at fixed radius `+1`
3. surface radius exponent at fixed mass `-1`
4. far mass exponent at fixed probe `+1`
5. far radial exponent at fixed source `-1`
6. weak-regime additivity

## Guardrails

No G, no macroscopic amplitude calibration, no native amplitude rescaling, no fitted/tuned K, no inserted `1/r` response, no spherical-equilibrium shortcut, no Rmax, no cosmology, no lensing target, no Quantum Engine, no Planck input.

## Run

```bash
PYTHONPATH=. python pbuf/labs/foundation/c_state_bounded_strain_bridge001.py
```

Valid scientific statuses:

- `C_STATE_BOUNDED_STRAIN_BRIDGE_STRUCTURE_SUPPORTED`
- `C_STATE_BOUNDED_STRAIN_BRIDGE_PARTIAL_SUPPORT`
- `C_STATE_BOUNDED_STRAIN_BRIDGE_NOT_SUPPORTED`

Partial/null outcomes are scientifically valid. The runner must not modify, repair, tune, rescale, reinterpret, or merge anything.
