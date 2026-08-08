# Native Full Received-State Information Retention 001

## Question

Before collapsing the observer end to a small set of weak-lensing observables, how much independent information is already present in the frozen received 3D ray state, and how much is lost as the observer representation is reduced?

This is an observer/decoding audit only. It does not alter the PBUF propagation physics.

## Frozen upstream lane

For each of the five canonical local clusters, build exactly one frozen current-native received ray state:

`local benchmark source -> current native fast/slow transfer -> PM1/PS2 -> M10 -> LOS -> G3D -> received 3D rays`

That state is mapped once to one target-blind global tangent detector screen. Every decoded channel for that cluster consumes the same received rays and the same screen.

## Full decoded channel bank first

The audit deliberately retains a wide bank before any reduction.

### Established 2D decoder channels

All three output fields from all eight established extraction methods are retained independently:

- histogram density
- kernel density
- affine Jacobian
- covariance area
- displacement divergence
- KNN density
- polar Jacobian
- covariance transport

For each method retain:

- convergence-like field
- gamma1-like field
- gamma2-like field

This contributes 24 channels.

### Explicit received-3D channels

The audit additionally retains detector-frame quantities that exist because the ray trajectory was actually propagated in 3D:

- mean transverse displacements `du`, `dv`
- mean normal/depth displacement `dw`
- mean received direction components `t1`, `t2`, `tn`
- within-bin standard deviations of all six quantities
- displacement covariance cross terms
- the complete local 3x2 differential receipt map
  `d(e1,e2,n final) / d(u0,v0)`

This contributes 21 additional channels.

Total predeclared decoded inventory:

`24 + 21 = 45 channels`

This is a broad fixed inventory, not a claim that no other mathematical decoder could ever be constructed.

## Direct 3D-versus-2D information test

Before binning, compare the standardized per-ray received state:

`(du, dv, dw, t1, t2, tn)`

against its transverse-only reduction:

`(du, dv, t1, t2)`.

Report numerical rank, effective rank, participation ratio, and the incremental rank carried by depth/normal channels.

## Progressive observer reductions

Using target-blind column standardization for information-geometry diagnostics only, compare:

1. `full_decoded_bank` — all 45 channels;
2. `full_minus_explicit_depth3d` — remove explicit depth/normal family;
3. `all_established_2d_decoders` — retain only the 24 established 2D decoder outputs;
4. `canonical_jacobian_three_field` — affine-Jacobian convergence/gamma1/gamma2 only;
5. `density_three_channel_control` — three independent density convergence channels.

For every stage report:

- numerical rank;
- effective rank;
- participation ratio;
- top-five singular-value variance fraction;
- median absolute channel correlation;
- effective-rank fraction retained relative to the full decoded bank;
- numerical-rank fraction retained relative to the full decoded bank.

The standardization exists only so heterogeneous observer channels can be compared in an internal SVD/correlation audit. It is never fed into a physical observable and is not amplitude tuning.

## Channel-family diagnostics

Also report intrinsic information geometry for the predeclared families:

- density
- area
- differential shape
- 2D displacement
- 3D displacement
- direction
- depth/normal 3D
- 3D differential/J3

This helps distinguish genuinely complementary receiver channels from redundant alternate decoders.

## Observed lensing products

Observed kappa/gamma fields are not used to create, orient, normalize, select, weight, or reduce any channel.

Only after the full bank and all reduction stages have been constructed, each decoded channel is compared descriptively with observed `kappa`, `gamma1`, `gamma2`, and total `gamma`.

The lab prints the top observer-relevance channels only as an after-the-fact diagnostic. No winner is automatically promoted.

## Guardrails

- one frozen current-native G3D state per cluster;
- one target-blind global tangent screen per cluster;
- decode full channel bank before reductions;
- no historical `strength=0.18`;
- no inferred `0.02925` coefficient;
- no fitted channel weights;
- no observational regression;
- no target-derived channel construction;
- no physical-output rescaling;
- no cluster-specific decoder selection;
- no changes to source construction, A8, pair transfer, PM1/PS2, M10, LOS, G3D, or trajectory integration;
- observations are after-the-fact comparison targets only.

## Run

```bash
PYTHONPATH=. python pbuf/labs/foundation/native_full_received_state_information_retention001.py
```

## Interpretation target

The key discriminator is whether the full received-state representation has materially greater intrinsic dimensionality than its 2D/three-field reductions, and especially whether explicit depth/normal receipt channels add independent dimensions.

If information rank falls systematically as the observer representation is collapsed while the upstream state is unchanged, that supports an observer-side information-loss diagnosis rather than a missing propagation mechanism.
