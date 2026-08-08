# Native Observable Extraction Method Sweep 001

## Question

Does the remaining weak-lensing mismatch arise primarily from the algorithm used to convert an already-received detector-plane ray pattern into `kappa`, `gamma1`, and `gamma2`?

## Frozen upstream state

Each cluster generates exactly one current-native G3D received ray state:

`canonical local benchmark source -> current native fast/slow transfer -> PM1/PS2 -> M10 -> LOS -> frozen G3D -> received 3D rays`

The received rays are then mapped once onto one target-blind detector plane perpendicular to the global mean received ray direction. The detector first axis is anchored by projecting global +x into the plane so its orientation is deterministic and not benchmark-selected.

Every extraction method receives the exact same initial/final detector-plane coordinates.

## Extraction sweep

Eight methods are tested in parallel:

1. `histogram_density` — occupancy-ratio extraction.
2. `kernel_density` — Gaussian KDE extraction.
3. `jacobian_affine` — per-bin affine ray-bundle Jacobian.
4. `covariance_area` — finite covariance-area distortion.
5. `displacement_divergence` — gradient of the mean detector-plane displacement field.
6. `knn_density` — adaptive k-nearest-neighbour density.
7. `polar_jacobian` — rotation-free polar stretch of the per-bin affine Jacobian.
8. `covariance_transport` — symmetric positive-definite covariance transport between the initial and final ray bundles.

The first six reuse established extraction algorithms from `observable_lab001.py` where they remain applicable to the current two-dimensional detector screen. The last two explicitly test rotation removal and covariance-preserving bundle deformation.

## Diagnostics

For every method and every cluster report:

- finite coverage;
- kappa Pearson and Spearman;
- raw kappa RMS amplitude ratio;
- total-shear Pearson and Spearman;
- raw total-shear amplitude ratio;
- gamma1 and gamma2 amplitude ratios;
- predicted gamma2/gamma1 balance;
- shear-orientation cosine;
- cross-cluster kappa-amplitude CV.

Scientific agreement is reported only. It is not an execution gate.

## Guardrails

- one received current-native G3D state per cluster;
- one fixed target-blind tangent detector screen per cluster;
- extraction algorithm is the only swept layer;
- no historical `strength=0.18`;
- no unit-control replacement;
- no inferred `0.02925` coefficient;
- no fitting or tuning;
- no output normalization or rescaling;
- no target-derived screen orientation;
- no cluster-specific method selection;
- no modifications to source construction, A8, pair transfer, PM1/PS2, M10, LOS, or G3D propagation.

## Run

```bash
PYTHONPATH=. python pbuf/labs/foundation/native_observable_extraction_method_sweep001.py
```

## Status values

- `NATIVE_OBSERVABLE_EXTRACTION_METHOD_SWEEP_EXECUTED`
- `NATIVE_OBSERVABLE_EXTRACTION_METHOD_SWEEP_PARTIAL_EXECUTION`
- `NATIVE_OBSERVABLE_EXTRACTION_METHOD_SWEEP_NOT_ESTABLISHED`

No method is automatically promoted from benchmark agreement. The sweep is discriminatory evidence for which observable-extraction families deserve follow-up.