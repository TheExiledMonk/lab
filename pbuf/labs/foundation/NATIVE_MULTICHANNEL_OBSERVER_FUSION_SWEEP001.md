# Native Multichannel Observer Fusion Sweep 001

## Question

Does the received native ray field carry weak-lensing information in several distinct detector channels that should be layered together, rather than forcing one extraction operator to recover all observables?

## Frozen upstream state

For each of the five canonical local benchmark clusters, generate one and only one current-native received G3D ray state:

`established local benchmark source -> current native fast/slow transfer -> PM1/PS2 -> M10 -> LOS -> G3D -> received 3D rays`

Map those rays once onto one target-blind global tangent detector screen. Every raw channel and every stacked candidate consumes those exact same detector-plane ray coordinates.

Nothing upstream may differ between candidates.

## Raw detector channels

### Convergence-like channels

- `k_knn` — adaptive nearest-neighbour ray-density change.
- `k_kernel` — Gaussian KDE ray-density change.
- `k_area` — covariance-area change.
- `k_jacobian` — affine Jacobian determinant channel.
- `k_divergence` — displacement-divergence channel.

### Shear-like channels

- `g_jacobian` — affine Jacobian `(gamma1,gamma2)`.
- `g_covariance` — SPD covariance-transport `(gamma1,gamma2)`.
- `g_divergence` — displacement-gradient `(gamma1,gamma2)`.
- `g_knn` — KNN extractor displacement-derived `(gamma1,gamma2)`.

These channels are not fitted to the observations. They are existing geometrical/density measurements of the same frozen ray field.

## Predeclared stacks

Nine stacks are tested with no fitted coefficients:

1. `knn_plus_jacobian`
2. `knn_plus_covariance`
3. `knn_plus_divergence`
4. `knn_plus_knn`
5. `kernel_plus_jacobian`
6. `area_plus_jacobian`
7. `jacobian_all_control`
8. `covariance_all_control`
9. `divergence_all_control`

Each stack selects one complete raw convergence channel and one complete raw shear channel. No averaging, learned weight, amplitude matching, target normalization, or cluster-specific choice is allowed.

## Combined image diagnostic

For every stack, form a single geometric channel-overlay magnitude

`C = sqrt(kappa^2 + gamma1^2 + gamma2^2)`

and compare it against the same invariant formed from the observed products. This is a diagnostic of whether the layered channel state reconstructs the overall observed spatial structure more clearly than any one extractor alone.

## Metrics

For every cluster and stack report:

- kappa Pearson, Spearman, and raw RMS amplitude ratio;
- gamma Pearson, Spearman, and raw RMS amplitude ratio;
- gamma1 and gamma2 amplitudes;
- shear orientation cosine;
- combined-image Pearson, Spearman, and raw amplitude ratio;
- cross-cluster kappa-amplitude CV.

Raw channel metrics are also retained in the JSON output so the contribution of each detector channel remains visible.

## Guardrails

- one frozen current-native received G3D state per cluster;
- one frozen target-blind tangent detector screen per cluster;
- no historical `strength=0.18`;
- no unit-control replacement;
- no inferred `0.02925` coefficient;
- no output normalization or rescaling;
- no fitted channel weights;
- no fitting or tuning;
- no target-derived channel construction;
- no cluster-specific stack selection;
- observed lensing products used only for after-the-fact comparison.

## Run

```bash
PYTHONPATH=. python pbuf/labs/foundation/native_multichannel_observer_fusion_sweep001.py
```

## Status values

- `NATIVE_MULTICHANNEL_OBSERVER_FUSION_SWEEP_EXECUTED`
- `NATIVE_MULTICHANNEL_OBSERVER_FUSION_SWEEP_PARTIAL_EXECUTION`
- `NATIVE_MULTICHANNEL_OBSERVER_FUSION_SWEEP_NOT_ESTABLISHED`
