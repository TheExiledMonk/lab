# G3D to 2D Observer Method Sweep 001

## Question

Can the poor final weak-lensing observable performance arise primarily from the chosen 3D-to-2D observer mapping rather than from the already-tested received G3D ray state?

## Frozen upstream lane

For each of the five canonical local clusters, generate exactly one current-native received G3D state:

`existing local benchmark source -> current native fast/slow transfer -> PM1/PS2 -> M10 -> LOS -> frozen G3D propagation -> received 3D ray state`

No observer candidate may alter any upstream state or rerun the trajectory differently.

## Observer sweep

Seven mappings consume the same received 3D state:

1. `xy_current` — current fixed x-y screen control.
2. `xz_control` — fixed x-z screen.
3. `yz_control` — fixed y-z screen.
4. `tangent_global` — detector plane perpendicular to the global mean received ray direction.
5. `tangent_local` — detector plane perpendicular to each source-bin mean received ray direction.
6. `gram_polar` — intrinsic 3D sheet deformation using `U = sqrt(J3^T J3)`.
7. `pca_global` — target-blind global PCA screen derived only from received 3D endpoints.

The full received differential state is

`J3 = d(xf,yf,zf) / d(x0,y0)`.

Each candidate produces a 2x2 observer mapping `A`. The same extraction convention is then used for all candidates:

- `kappa = 1 - det(A)`
- `gamma1 = 0.5 * (A00 - A11)`
- `gamma2 = 0.5 * (A01 + A10)`

## Diagnostics

For each method and each cluster report:

- finite coverage;
- kappa Pearson, Spearman and raw RMS amplitude ratio;
- total shear Pearson, Spearman and raw RMS amplitude ratio;
- gamma1 and gamma2 raw RMS amplitude ratios;
- predicted gamma2/gamma1 balance;
- shear orientation cosine against the local observed shear components;
- cross-cluster kappa-amplitude CV.

Aggregate results are descriptive only. No method is automatically promoted from benchmark agreement.

## Guardrails

- one received G3D state shared by every observer candidate within a cluster;
- current native transfer only;
- no historical strength 0.18;
- no unit-control replacement;
- no 0.02925 coefficient;
- no output normalization or rescaling;
- no fitting or tuning;
- no cluster-specific observer choice;
- no target-derived screen orientation;
- no modification of source, A8, PM1/PS2, M10, LOS, or G3D trajectory;
- observed kappa/gamma/gamma1/gamma2 are comparison targets only.

## Run

```bash
PYTHONPATH=. python pbuf/labs/foundation/g3d_to_2d_observer_method_sweep001.py
```

## Status values

- `G3D_TO_2D_OBSERVER_METHOD_SWEEP_EXECUTED`
- `G3D_TO_2D_OBSERVER_METHOD_SWEEP_PARTIAL_EXECUTION`
- `G3D_TO_2D_OBSERVER_METHOD_SWEEP_NOT_ESTABLISHED`
