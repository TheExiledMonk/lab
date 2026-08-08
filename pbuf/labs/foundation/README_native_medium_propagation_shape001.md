# Native medium propagation shape 001

## Question

Does the already-supported native accumulation chain

`rho -> existing A8 transport -> raw c_state -> six-neighbor bounded-strain equilibrium -> accumulated medium state u`

produce a viable propagation-response **shape** directly from the medium state, without importing GR/LCDM potential machinery or fitting a physical propagation amplitude?

## Structural propagation hypothesis

The only new hypothesis in this lab is that a propagating disturbance crossing an inhomogeneous accumulated medium receives a transverse directional response proportional to the path-integrated transverse gradient of that accumulated state:

`Delta k_x ~ integral ds partial_x u`

The proportionality coefficient is intentionally left unspecified. Therefore this lab tests only coefficient-independent structure.

## Predeclared checks

1. zero source gives zero propagation response;
2. a centered propagation path gives zero transverse response;
3. reflection symmetry gives opposite transverse response on opposite sides of the source;
4. the response points toward the centered source;
5. weak-regime response scales as source mass^1;
6. response magnitude scales approximately as impact_parameter^-1.

## Frozen inputs

- existing A8 raw `c_state` generation;
- existing bounded-strain accumulation bridge from the supported prior lab;
- `K0=1` structural normalization;
- `epsilon_max=1`;
- source radius `3.5`;
- reference mass `2.0`;
- mass ladder `0.5,1,2,4,8`;
- impact parameters `6,7,8,9,10`;
- accumulation grid inherited unchanged from the supported bridge.

## Guardrails

No G. No GR potential decomposition. No LCDM. No physical propagation amplitude. No amplitude calibration. No native rescaling. No fitted/tuned K. No inserted `1/r` response. No inserted `1/b` response. No spherical-equilibrium shortcut. No Rmax. No cosmology. No observed lensing target. No kappa/shear observation. No Quantum Engine. No Planck input.

A failure is a scientific result and is not permission to alter the model.

## Run

```bash
PYTHONPATH=. python pbuf/labs/foundation/native_medium_propagation_shape001.py
```

## Valid statuses

- `NATIVE_MEDIUM_PROPAGATION_SHAPE_SUPPORTED`
- `NATIVE_MEDIUM_PROPAGATION_SHAPE_PARTIAL_SUPPORT`
- `NATIVE_MEDIUM_PROPAGATION_SHAPE_NOT_SUPPORTED`

The runner must return HEAD/branch, exit code, complete raw stdout/stderr, `git status --short`, and `git stash list`, and must not modify, repair, tune, rescale, reinterpret, or merge anything.
