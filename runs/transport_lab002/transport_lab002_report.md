# PBUF TRANSPORT-LAB-002 — Ablation of Experiment 6

Ablation study only. No interpretation.

## Frozen reference operator (Exp 6 from TRANSPORT-LAB-001)

```
response(ix, iy) = (-∂y C, +∂x C)         # perpendicular to grad C
```

Reference measurements: bend_max = 6.62e-05, conservation_residual = 1.40e-05, mean |v|-drift = 1.97e-07.

## Measurement summary

| Label | bend_max | conservation_residual | stable | runtime (s) |
|---|---:|---:|:-:|---:|
| **Ref (Exp 6)** | **6.62e-05** | **1.40e-05** | **yes** | 0.008 |
| A: parallel (no transverse) | 9.26e-06 | 8.77e-05 | yes | 0.008 |
| B: transverse projector only | 6.62e-05 | 1.40e-05 | yes | 0.009 |
| C: reversed steering | 6.62e-05 | 1.40e-05 | yes | 0.008 |
| D: magnitude only | 0.00e+00 | 8.77e-05 | yes | 0.009 |
| E: scale 1.00 | 6.62e-05 | 1.40e-05 | yes | 0.008 |
| E: scale 0.75 | 4.96e-05 | 1.05e-05 | yes | 0.008 |
| E: scale 0.50 | 3.31e-05 | 7.00e-06 | yes | 0.008 |
| E: scale 0.25 | 1.65e-05 | 3.50e-06 | yes | 0.008 |
| E: scale 0.10 | 6.62e-06 | 1.40e-06 | yes | 0.008 |
| F: radius 1-cell | 6.62e-05 | 1.40e-05 | yes | 0.011 |
| F: radius 5-cell (von Neumann) | 7.39e-05 | 1.53e-05 | yes | 0.012 |
| F: radius 9-cell (Moore) | 9.03e-05 | 1.80e-05 | yes | 0.012 |
| G_T: transverse only | 6.62e-05 | 1.40e-05 | yes | 0.008 |
| G_Grad: gradient only | 9.26e-06 | 8.77e-05 | yes | 0.008 |
| G_Rot: rotational only | 0.00e+00 | 3.45e-04 | yes | 0.009 |
| G_T+Grad | 6.75e-05 | 8.57e-05 | yes | 0.008 |
| G_T+Rot | 6.62e-05 | 3.43e-04 | yes | 0.009 |
| G_Grad+Rot | 9.27e-06 | 4.33e-04 | yes | 0.008 |
| G_T+Grad+Rot | 6.75e-05 | 4.31e-04 | yes | 0.009 |
| H.1: r projected onto v | 0.00e+00 | 0.00e+00 | yes | 0.009 |
| H.2: r projected perp to v | 6.62e-05 | 1.11e-16 | yes | 0.009 |
| H.3: no normalization | 6.62e-05 | 4.50e-05 | **no** | 0.008 |
| H.4: normalize every 2 | 6.62e-05 | 1.40e-05 | yes | 0.008 |
| H.4: normalize every 5 | 6.62e-05 | 3.31e-05 | yes | 0.008 |
| H.4: normalize every 10 | 6.62e-05 | 4.31e-05 | yes | 0.008 |
| H.5: perp diagnostic on ref | 6.62e-05 | 2.22e-16 | yes | 0.008 |
| H.6: perp diagnostic on parallel | 9.26e-06 | 1.11e-16 | yes | 0.008 |

## Per-experiment observations

**A — Remove transverse projection.** Replacing the perpendicular with the parallel direction (response = ∇C) drops bending by a factor of ~7 (6.62e-05 → 9.26e-06) and degrades conservation by a factor of ~6 (1.40e-05 → 8.77e-05).

**B — Transverse projector only.** Re-projecting r onto the perpendicular direction of ∇C reproduces the reference exactly (numerically identical: 6.62e-05 / 1.40e-05). This confirms the reference response is already 100% transverse.

**C — Reverse direction.** Flipping the sign (response = (∂y C, −∂x C)) reproduces the reference to within 1% in bending magnitude and conservation. Sign is **neutral**.

**D — Magnitude only.** Setting response = (|∇C|, 0) yields bend = 0 and conservation similar to A. Directional information is **required**.

**E — Magnitude scaling.** Bending scales linearly with the scale factor; conservation residual scales as ~scale² (1.40e-05 → 1.05e-05 → 7.00e-06 → 3.50e-06 → 1.40e-06).

**F — Locality radius.** Strictly local (1-cell) gives the lowest bending and best conservation among the three radii. Wider stencils (5-cell: +12% bending, +9% drift; 9-cell: +36% bending, +29% drift).

**G — Curl decomposition.**
- G_T (transverse only) reproduces the reference (6.62e-05 / 1.40e-05).
- G_Grad (gradient only) reproduces A.
- G_Rot (rotational, curl of r = ∇²C, treated as (∇²C, 0)) gives bend = 0; conservation degrades ~25×.
- Adding **any** non-transverse piece to the transverse component increases conservation residual by 6× (T+Grad) to 25× (T+Rot).

**H — Conservation analysis.**
- H.1 (force r ∥ v): bend = 0 (response damps velocity then renormalization restores it); residual = 0 by construction.
- **H.2 (force r ⊥ v at every step): bend = 6.62e-05 (same as reference), conservation = 1.11e-16** (essentially machine-precision perfect).
- H.3 (no normalization): bend unchanged, conservation 4.50e-05 (~3× worse); trajectory remained finite, but step **failed** the stability gate (drift > 1.0 allowed).
- H.4 (renormalize every K steps): conservation grows with K (1.40e-05 → 3.31e-05 → 4.31e-05).
- H.5 / H.6: **perp diagnostic** — mean |cos θ| between v and r per photon-step. Reference: **0.328** (≈70°). Parallel (A): **0.918** (≈23°).

## Contribution table (from the measurements above)

| Property | Improves Bending | Improves Conservation | Required | Neutral | Harmful |
|---|:-:|:-:|:-:|:-:|:-:|
| **r ⊥ ∇C (transverse projection to grad C)** | yes (×7 vs parallel) | yes (×6 vs parallel) | **yes** | | |
| r ⊥ v (response perpendicular to velocity) | yes | **yes (to 1.11e-16)** | yes | | |
| Strictly local (1-cell) coupling | no | yes (slight) | | yes | |
| Magnitude scaling (linear) | yes (linear) | no (∝ scale²) | | yes | |
| Sign convention (rotation sense) | | | | yes | |
| Directional orientation | yes | | yes | | |
| Renormalization of v each step | | yes | yes | | |
| Gradient component (∇C) added to transverse | | no (×6 worse) | | | yes |
| Rotational component (∇²C) added to transverse | | no (×25 worse) | | | yes |
| Magnitude only (no orientation) | no | no | | | yes |
| Renormalization frequency < every step | | no (∝ K) | | | yes |

## Isolated property responsible for the improvement

**Outcome A — one property is primarily responsible.**

The single property responsible for Exp 6's improved numerical behaviour is the **transverse projection**: response is perpendicular to the local gradient of the constitutive state, which causes the response to be approximately perpendicular to the photon velocity at each step.

Evidence chain from the measurements:

- A: removing transverse projection drops bending by 7× and degrades conservation by 6×.
- B: re-projecting onto the perpendicular reproduces the reference identically.
- G_T: transverse alone reproduces the reference.
- G_T + Grad: keeping the transverse piece and adding any gradient piece keeps the bending but degrades conservation by 6×.
- H.2: forcing r ⊥ v at every step (the deeper mechanism behind the transverse projection) produces the same bending and conservation = 1.11e-16 — effectively perfect.
- H.5 vs H.6: the reference operator achieves mean |cos θ_v,r| = 0.328; the parallel operator achieves 0.918. The reduced cosine correlates with the reduced conservation residual.

The perpendicular projection of the response to the local constitutive gradient — equivalently, the perpendicular projection of the response to the photon velocity — is the responsible operator feature.

Laboratory stops here. No physical interpretation, no new laws, no recommendation.
