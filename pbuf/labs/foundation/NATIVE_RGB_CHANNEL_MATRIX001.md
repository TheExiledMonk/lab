# Native RGB-Like Channel Matrix 001

## Question

Does the frozen received native ray field behave like a multi-channel image in which different receiver operators carry different parts of the information, so that independently selecting and combining the kappa, gamma1, and gamma2 channels yields a clearer reconstructed observer state?

## Frozen upstream lane

For each of the five canonical local clusters, generate exactly one received state:

`local benchmark source -> current native fast/slow transfer -> PM1/PS2 -> M10 -> LOS -> G3D -> received 3D rays -> one target-blind global tangent detector screen`

Every candidate consumes the same detector-plane rays within a cluster.

## Full channel matrix

The audit evaluates every combination of:

- 5 convergence receivers;
- 4 gamma1 receivers;
- 4 gamma2 receivers.

Total:

`5 x 4 x 4 = 80` observer states.

Convergence receivers:

- `k_knn`
- `k_kernel`
- `k_area`
- `k_jacobian`
- `k_divergence`

Gamma1 and gamma2 may independently come from:

- `g_jacobian`
- `g_covariance`
- `g_divergence`
- `g_knn`

This explicitly tests whether the two shear components are most cleanly received by different operators.

## RGB-like overlays

For every 3-channel state, three fixed no-weight overlays are formed:

1. `vector_magnitude = sqrt(kappa^2 + gamma1^2 + gamma2^2)`
2. `absolute_overlay = |kappa| + |gamma1| + |gamma2|`
3. `signed_overlay = kappa + gamma1 + gamma2`

The same formula is applied to the observed channels for comparison.

These are observer diagnostics only. They are not fitted observables and do not change the underlying channel amplitudes.

## Diagnostics

For each of the 80 states and each cluster report:

- kappa Pearson, Spearman, and raw amplitude ratio;
- gamma1 Pearson, Spearman, and raw amplitude ratio;
- gamma2 Pearson, Spearman, and raw amplitude ratio;
- total gamma morphology/amplitude;
- shear orientation consistency;
- vector-magnitude morphology/amplitude;
- absolute-overlay morphology/amplitude;
- signed-overlay morphology/amplitude;
- cross-cluster amplitude stability.

Also print the top diagnostic states for each overlay metric. These rankings are descriptive only and do not automatically promote a state.

## Guardrails

- one received G3D state shared by all 80 states within a cluster;
- one target-blind tangent detector screen shared by all 80 states within a cluster;
- identity channel transfer only;
- no fitted channel weights;
- no output normalization or rescaling;
- no target-derived channel construction;
- no cluster-specific channel choice;
- no historical strength `0.18`;
- no unit-control replacement;
- no inferred `0.02925` coefficient;
- no modification of source, A8, PM1/PS2, M10, LOS, or G3D trajectory;
- observed products are comparison targets only.

## Run

```bash
PYTHONPATH=. python pbuf/labs/foundation/native_rgb_channel_matrix001.py
```

## Status values

- `NATIVE_RGB_CHANNEL_MATRIX_EXECUTED`
- `NATIVE_RGB_CHANNEL_MATRIX_PARTIAL_EXECUTION`
- `NATIVE_RGB_CHANNEL_MATRIX_NOT_ESTABLISHED`
