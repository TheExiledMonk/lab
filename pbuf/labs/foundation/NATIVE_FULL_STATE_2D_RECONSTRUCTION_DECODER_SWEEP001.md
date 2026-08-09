# PBUF FOUNDATION — Full-State 2D Reconstruction Decoder Sweep 001

## Purpose

PR #104 established that the frozen current-native received G3D ray state contains substantially more independent observer-side information than the conventional three-field collapse. This lab asks the next question without changing upstream physics:

> Can the full retained received-state bank be converted into coherent 2D observer maps by target-blind reconstruction rules before any conventional weak-lensing collapse is imposed?

This is an observer/readout audit only. It is not a new propagation mechanism and it does not fit a lensing law.

## Frozen upstream state

For each canonical cluster the lab builds exactly one current-native received G3D state using the existing source, zero-flux terminal fast/slow A8 channels, exact pair amplitude, c_state bond geometry, PM1/PS2, M10, LOS and G3D trajectory chain.

One target-blind global tangent detector screen is shared by all decoder candidates.

No upstream physics is changed.

## Full received-state bank

The lab reuses the full PR #104 decoded inventory before any reduction:

- 24 established 2D extraction fields: 8 extraction methods × convergence/shear_g1/shear_g2;
- 21 explicit 3D receipt fields containing displacement, direction, spread, covariance and local 3×2 differential information.

Fixed requested inventory:

\[
45\text{ channels}.
\]

Constant/all-missing fields remain part of the requested inventory but are excluded from standardized decoder matrices because they carry no variance.

## Target-blind reconstruction sweep

All reconstruction candidates are completed before observed kappa/gamma products are requested.

The predeclared sweep includes:

1. full-bank L2 energy;
2. full-bank L1 magnitude;
3. full-bank signed mean;
4. family-balanced L2/L1/signed fusion;
5. first eight full-bank principal-component maps;
6. PCA energy maps retaining 2, 3, 5, 8 and 10 components;
7. whitened-PCA energy maps retaining 2, 3, 5, 8 and 10 components;
8. explicit-3D-only controls;
9. established-2D-only controls;
10. no-depth controls;
11. L2/L1/signed/PC1 reconstructions for every decoded channel family.

Per-channel centering/scaling is target-blind. Missing cells are represented as zero deviation after centering. PCA signs are fixed deterministically by requiring the largest-magnitude loading in each component to be positive. No observation is used to orient a component.

## After-the-fact comparison only

Once every candidate image is frozen, the lab compares it descriptively against four local observed targets:

- kappa;
- absolute kappa;
- total gamma magnitude;
- combined observer norm \(\sqrt{\kappa^2+\gamma_1^2+\gamma_2^2}\).

For each candidate/target pair the existing comparison helper reports Pearson morphology, Spearman morphology, RMS amplitude ratio and coverage.

The ranking is descriptive only. It does not select or promote a decoder and is never fed back into candidate construction.

## Guardrails

- canonical five-cluster local benchmark only;
- one frozen received G3D state per cluster;
- one target-blind detector screen per cluster;
- full 45-channel bank exists before reconstruction;
- same candidate inventory for every cluster;
- observations accessed only after candidate construction;
- no target-derived channel selection;
- no observational regression;
- no fitted channel weights;
- no physical-output rescaling;
- no cluster-specific decoder choice;
- no historical strength 0.18;
- no inferred 0.02925 replacement coefficient;
- no source/A8/M10/LOS/G3D modification.

## Interpretation boundary

A successful candidate would establish only that a target-blind observer reconstruction can recover useful 2D morphology from the rich received state. It would not yet establish the final physical mapping to conventional \((\kappa,\gamma_1,\gamma_2)\).

The intended architecture remains:

\[
\boxed{
\mathcal R^{3D}_{\rm received}
\rightarrow
\mathcal I^{2D}_{\rm observer}
\rightarrow
(\kappa,\gamma_1,\gamma_2)\text{ for scientific comparison}
}
\]

with the conventional three-field representation delayed until the last stage rather than imposed at receipt.

## Run

```bash
PYTHONPATH=. python pbuf/labs/foundation/native_full_state_2d_reconstruction_decoder_sweep001.py
```

Expected success status:

`FULL_STATE_2D_RECONSTRUCTION_DECODER_SWEEP_EXECUTED`
