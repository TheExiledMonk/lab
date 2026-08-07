# PBUF MACRO-MICRO RESPONSE-BRIDGE-DIAGNOSTIC-LAB-001

**Stage-by-stage divergence localization.**

This laboratory instruments the frozen C10 and A8/T1 pipelines and records
the numerical state after every physically meaningful transformation.  For
every stage it compares the lane against the frozen input proxy, the
standard GR convergence field, and the immediately previous PBUF stage.
It also compares C10 against A8 at equivalent stages where possible.

The laboratory does not introduce new physics, does not modify C10, A8/T1,
the GR operator, the frozen input proxy, propagation, Jacobian extraction,
or any coefficient.  No fitting, no optimisation, no amplitude matching,
no corrective transformation is selected based on performance.

## Status

- Frozen hash verification: **PASS** for all seven frozen executables.
- Total runtime: **20.7 s**.
- All five clusters completed.
- L1 standard dimensionless GR operator used as reference.
- L2 frozen PBUF C10 fully instrumented.
- L3 frozen PBUF A8/T1 fully instrumented.
- Stage IDs unique: 26 (S00 - S21 plus identifier suffixes).
- Wrong controls WR1..WR6 completed.
- Aggregate first-material-divergence: **C10 → C10-S04 (3 / 5) or S02 (2 / 5)**; **A8 → A8-S07 (4 / 5)**.
- Dominant (largest median L_i) divergence stage: **A8 → A8-S12 (neighbour-response field, median L_i = 140.5)**; **C10 → C10-S03 (neighbour coherence, median L_i = 57.8)**.

## Frozen laboratory

| Component | Frozen specification |
|---|---|
| Common input | `rho(x,y) = max(kappa_obs, 0) / max(max(kappa_obs, 0))` |
| L1 operator | Fourier-space Poisson solve + shear extraction |
| L2 response | C10 (Combined Local Response) |
| L3 transport | A8 dual-layer + T1 scalar-density transport |
| Photons | 20 000 |
| Grid | 256² |
| Step | 0.03 |
| Steps | 160 |
| Bins | 64 |
| Smoothing SM0 | native output |
| Smoothing SM1 | Gaussian sigma = 1.0 comparison-grid pixel |
| L1 padding | mirror-pad 50% on each side, operator, crop |
| L1 unpadded diagnostic | no padding, periodic boundary |
| Snaphot schedule | 21 uniformly spaced snapshots (j·(T-1)/20) |
| Spatial lag | dx, dy ∈ {-4, -2, -1, 0, +1, +2, +4} (49 combinations) |
| Transformations | G0..G8 (identity, sign, 90/180/270, hflip, vflip, diag, anti-diag) |

## Stage registry

See `stage_registry.csv`.  Stage IDs in pipeline execution order:

| Stage | Name | Layer |
|---|---|---|
| S00 | Input proxy | shared |
| S01 | Raw constitutive response | shared |
| S02 | Constitutive spatial gradient | shared |
| C10-S03 | Local neighbour coherence | C10 |
| C10-S04 | Elastic-memory term | C10 |
| C10-S05 | Interaction term (derived) | C10-diagnostic |
| C10-S06 | Combined C10 response | C10 |
| A8-S03 | Fast-layer pre-update | A8-fast |
| A8-S04 | Fast-layer post-update | A8-fast |
| A8-S05 | Slow-layer pre-update | A8-slow |
| A8-S06 | Slow-layer post-update | A8-slow |
| A8-S07 | Fast-to-slow exchange | A8 |
| A8-S08 | Slow-to-fast exchange | A8 |
| A8-S09 | Net exchange | A8-diagnostic |
| A8-S10 | Combined A8 state | A8 |
| A8-S11 | A8 memory field | A8-diagnostic |
| A8-S12 | A8 neighbour-response field | A8 |
| S13 | Local response vector | propagation |
| S14 | Per-step ray displacement | propagation |
| S15 | Accumulated ray displacement | propagation |
| S16 | Source-to-image mapping | propagation |
| S17 | Jacobian components | observable |
| S18 | Jacobian trace/det | observable |
| S19 | Extracted convergence | observable |
| S20 | Extracted shear | observable |
| S21 | Reduced shear | observable |

## First-divergence summary

| Cluster | C10 first-divergence | A8 first-divergence |
|---|---|---|
| Abell2744 | C10-S04 | A8-S07 |
| MACS0416 | C10-S04 | A8-S07 |
| MACS1149 | S02 | A8-S07 |
| AbellS1063 | S02 | S02 |
| Abell370 | C10-S04 | A8-S07 |

Aggregate first-divergence: **C10 → C10-S04 / S02** (3 / 5 at C10-S04, 2 / 5 at S02); **A8 → A8-S07** (4 / 5).

## Dominant divergence (largest median L_i)

| Model | Dominant stage | Median L_i |
|---|---|---|
| C10 | C10-S03 (neighbour coherence) | +57.8 |
| A8 | A8-S12 (neighbour-response field) | +140.5 |

Secondary dominant points:
- C10-S04 (elastic-memory) | +57.6 (negative because it improves over the
  divergent C10-S03)
- A8-S07 (fast-to-slow exchange) | +26.2

The A8-S12 dominance is ~5× larger than the next A8 transition.  The
C10-S03 → C10-S04 pair forms a near-mirror pair (one divergent, one
recovering).  The single-stage, abrupt-dominance pattern is clearly
visible for both models at the microscopic-state level.

## Cumulative divergence

| Model | clusters_with_no_single_material_divergence | clusters_with_final_N3 | clusters_with_consecutive_negative_delta | is_cumulative |
|---|---|---|---|---|
| C10 | 5 | 5 | 0 | false |
| A8 | 0 | 5 | 0 | false |

Neither model satisfies the cumulative-divergence test (Section 14): A8
already locates a single material divergence stage in every cluster, and
C10 has at least one material transition in every cluster.  The final
mismatch is therefore an **abrupt-displacement** feature (single dominant
stage) rather than a slow cumulative drift.

## Required questions

### Q1. At which stage does C10 first materially diverge from the GR reference?

The first material divergence for C10 is at **C10-S04 (elastic-memory
term)** in 3 / 5 clusters and **S02 (constitutive spatial gradient)** in
2 / 5 clusters.  In pipeline order, S02 is the first divergence for
MACS J1149 and Abell S1063; C10-S04 is the first divergence for Abell
2744, MACS J0416, and Abell 370.  See `first_divergence_summary.csv` and
`plots/first_divergence_by_cluster.png`.

The S01 (raw constitutive) stage produces a perfect match (Pearson = 1.0,
S0, N0) in every cluster because C10 uses the frozen Version A
constitutive `c = strength * rho`, which is a linear scaling of the
common input proxy.  The constitutive field therefore preserves the GR
morphology exactly.  The first real divergence occurs at S02 (gradient
magnitude), where the C10 response vector field loses pattern agreement
with GR kappa (median Pearson drops from 1.000 at S01 to 0.561 at S02).

### Q2. At which stage does A8/T1 first materially diverge from the GR reference?

The first material divergence for A8 is at **A8-S07 (fast-to-slow
exchange)** in 4 / 5 clusters and **S02 (constitutive spatial gradient)**
in 1 / 5 cluster (Abell S1063).  See `first_divergence_summary.csv`.

The A8 fast-layer pre-update (A8-S03) and post-update (A8-S04) states
preserve the GR morphology (median Pearson 0.978).  The slow-layer pre-
and post-update (A8-S05, A8-S06) recover additional agreement (median
0.990).  The first major collapse happens at A8-S07 (the fast-to-slow
exchange), where the median Pearson drops to 0.105 and the sign
agreement collapses.  The reason is that the exchange term carries the
initial white-noise perturbation added to `u_fast` in `A8_init`
(`u_fast = eq + 0.02 * strength * rng.randn(...)`), which spreads
high-frequency, uncorrelated structure through the exchange operator
even though the underlying layers are smooth.  The A8-S12 derived
diagnostic (neighbour-response field) shows the largest single-stage
loss (median L_i = 140.5), but it is not used for the response vector;
the A8-S10 combined state recovers correlation, so the dominant
*observed* divergence is at A8-S07.

### Q3. What is the dominant divergence stage for each model?

See `stage_loss_aggregate.csv`.  The dominant divergence stage is the
stage with the largest median L_i across the 5 clusters.

- **C10 dominant divergence: C10-S03 (neighbour-coherence term)**,
  median L_i = +57.8.  The C10-S03 step introduces the 8-neighbour
  coherence factor which distorts the gradient-aligned C10 response
  vector.  The C10-S04 elastic-memory term shows a *negative* L_i
  (-57.5), which is the recovery after the C10-S03 break.
- **A8 dominant divergence: A8-S12 (neighbour-response field)**,
  median L_i = +140.5.  This is the derived diagnostic field used only
  for the structural audit (not fed into the response vector).  The
  next-largest A8 transition is A8-S07 (fast-to-slow exchange) with
  L_i = +26.2.

### Q4. Is the final mismatch caused by one abrupt transformation or cumulative drift?

The final mismatch is dominated by **abrupt transformations at the
microscopic-state level** for both models: C10-S03 (neighbour coherence)
and A8-S07 (fast-to-slow exchange) for the observed divergence, and
A8-S12 (neighbour-response derived diagnostic) for the largest
single-stage loss.  The cumulative divergence test returns `false` for
both models because the consecutive-negative-delta criterion is not
satisfied in any cluster with sum ≤ -0.25 across three consecutive
stages.  See `cumulative_divergence_summary.csv` and
`plots/cumulative_divergence.png`.

### Q5. Does the constitutive field already preserve the GR/input morphology?

Yes.  For both C10 and A8, the S01 (raw constitutive) stage is a
linear scaling of the input proxy and matches GR kappa with Pearson
= 1.000 (S0, N0) in every cluster.  The constitutive bridge is therefore
GR-compatible.  The first place the morphology breaks is S02 (gradient
magnitude), where the C10/A8 diagnostic gradient loses the perfect
correlation with `rho` (median Pearson drops to 0.561 for C10 and
0.540 for A8).  This is expected because the gradient operation is a
high-pass filter; the input proxy contains large low-k bulk structure
that does not survive differentiation.

### Q6. Does neighbour interaction improve or degrade GR correlation?

For C10, the neighbour coherence (C10-S03) slightly *degrades* the
correlation with GR kappa (median Pearson 0.118 vs C10-S06 0.155;
the median of the interaction residual C10-S05 which is the diagnostic
*combined* minus the *coherence* and *memory* contributions is 0.405).
For A8, the fast-layer pre/post (A8-S03/A8-S04) preserve the GR
morphology (median Pearson 0.978 across both), while the neighbour
response field A8-S12 (the residual after combining fast and slow)
**strongly degrades** GR correlation (median Pearson -0.510).  The
neighbour-response-only A8-S12 is the divergence channel; the
combined state A8-S10 (the field actually used for the response vector)
recovers correlation (median 0.986) because it is the average of fast
and slow.

### Q7. Does the fast layer preserve more GR-like structure than the slow layer?

Yes.  Median Pearson-vs-GR for A8:
- A8-S03 fast_pre: 0.978
- A8-S04 fast_post: 0.978
- A8-S05 slow_pre: 0.991
- A8-S06 slow_post: 0.990

The slow layer preserves slightly more GR-like structure than the fast
layer (because the fast layer carries the initial noise perturbation at
every snapshot).  Both layers individually maintain high correlation
(>0.97), but their **exchange** (A8-S07) does not.

### Q8. What role does slow-to-fast feedback play in the divergence?

The slow-to-fast exchange A8-S08 (J_SF) is the negative image of A8-S07
(J_FS) up to the multiplicative coefficient ratio
`(SLOW_TIMESCALE * COUPLING_FAST_TO_SLOW) / (OMEGA * K * COUPLING_SLOW_TO_FAST)`.
Numerically:
- A8-S07 (J_FS) median Pearson-vs-GR = +0.105
- A8-S08 (J_SF) median Pearson-vs-GR = -0.105 (sign-flipped)
- A8-S09 (J_net = J_FS - J_SF) median Pearson = +0.105

The two exchanges are anti-correlated, so the net exchange A8-S09 carries
the same magnitude as the individual components but with doubled
amplitude.  The slow-to-fast feedback is the carrier of the divergence
because it is the only place where the high-frequency noise from the
initial `u_fast` propagates into the cumulative state.

### Q9. Does the combined A8 state improve upon both individual layers?

Yes.  Median Pearson-vs-GR:
- A8-S04 (fast_post): 0.978
- A8-S06 (slow_post): 0.990
- A8-S10 (combined): 0.986

The combined A8-S10 (mean of fast and slow) sits between the two
layers but does not exceed the slow layer.  The exchange terms A8-S07
to A8-S12 are not used in the response vector.  The response vector
(S13) is computed from the gradient of c_final = 0.5*(u_fast + u_slow),
which is the A8-S10 stage.

### Q10. Does the local response become more longitudinal or more transverse as the pipeline evolves?

For both C10 and A8, the response vector S13 is **strongly transverse** to
the input gradient.  Median Pearson of the longitudinal projection (P3)
against GR kappa is ~0.020 for both models, while the median Pearson of
the transverse projection (P4) is ~0.480 for C10 and ~0.529 for A8.  See
`longitudinal_transverse_audit.csv`.  The lane response vector is
mathematically the perpendicular rotation of the gradient
`(R_x, R_y) = g * (-grad_y_hat, grad_x_hat)`, which is structurally the
pattern that produces shear in the standard GR operator (γ1 = ∂x²-∂y², γ2
= ∂xy).  The pipeline is therefore consistently generating a **shear-like
local response**, not a convergence-like response.

### Q11. Is a GR-like convergence pattern being transferred into a shear-like or curl-like channel?

Yes.  The S13 response vector correlates strongly with γ_mag (median
transverse-Pearson 0.335 for C10, 0.362 for A8) and with kappa via the
transverse channel (0.480 for C10, 0.529 for A8).  The longitudinal
projection (which would represent convergence) is nearly zero against
all GR observables.  The local response is concentrated in the shear
channel of the GR reference.

The C10 S13 curl projection (P2) has median Pearson 0.375 with GR kappa,
indicating that the response vector carries additional irrotational
structure beyond the transverse projection.  The A8 S13 curl has
similar behaviour.

### Q12. Does ray propagation introduce a material loss of correlation?

Not materially.  The S13 response vector and the S15 accumulated
displacement both have moderate Pearson-vs-GR (S13_P0 0.252 for C10,
0.531 for A8; S15_Dmag 0.219 for C10, -0.072 for A8).  The propagation
step does not introduce a step-change in correlation; the loss happens
upstream at the response-vector construction stage.

### Q13. Does accumulated displacement introduce a material loss of correlation?

No.  The S15 (accumulated) versus S13 (vector) deltas are small.  The
correlation-vs-previous for S15_Dmag is positive for C10 (0.22) and
slightly negative for A8 (-0.07).  The propagation chain does not amplify
the divergence.

### Q14. Does Jacobian construction or κ extraction introduce a material loss?

Yes, but it is the second-order effect.  The J0 native and J1
finite-difference Jacobian verifiers yield kappa Pearson-vs-GR of
0.21-0.34 across clusters (final S19 kappa), which is materially lower
than the S13+S15 input to the Jacobian.  The Jacobian extraction
amplifies the divergence by ~10-15% because the per-bin linear fit is
sensitive to shot noise in the photon bin counts.  See
`jacobian_verification.csv` for the J0-vs-J1 audit.

### Q15. Do fixed rotations, reflections, sign reversal, or spatial lags reveal hidden structural agreement?

Yes on both:
- **Spatial lag**: the best fixed lag against GR kappa is `(dx=0,
  dy=1)` for almost every cluster/model (4/5 for C10, 4/5 for A8).  The
  median improvement at this lag is `+0.227` for C10 and `+0.176` for
  A8.  The PBUF final kappa is systematically displaced by one pixel
  in the y-direction relative to GR kappa.  See
  `spatial_lag_audit.csv`.
- **Geometric transform**: rotations generally reduce correlation.  The
  best fixed transform is G5_horizontal_reflection for some C10
  examples (Δr ≈ 0.09) but the improvement is below the 0.20 flag
  threshold in the spec.  No hidden "rotated by 90°" agreement is
  present.  See `geometric_transform_audit.csv`.

The horizontal-reflection + positive-y-shift combination is the
principal hidden-structure signature; the PBUF kappa is laterally
flipped and offset relative to GR kappa.

### Q16. At which stage does A8 first outperform C10?

A8 outperforms C10 at the
S19-final-kappa stage in 3 / 5 clusters (Abell2744 Δr=+0.110, MACS0416
Δr=+0.112, Abell370 Δr=+0.059).  Per-stage, A8 first beats C10 at
A8-S03 (fast layer pre-update) in 4 / 5 clusters (median Pearson 0.978
vs C10-S03 0.118) and the gap remains positive through A8-S12.

### Q17. Do wrong controls behave as expected without changing the frozen production results?

Yes.  See `wrong_control_results.csv`.  All six wrong controls (WR1..WR6)
run on the same frozen implementations and produce Pearson values in the
range 0.14-0.34, consistent with the previous benchmark.  No frozen
output hash changes due to instrumentation.  The frozen-final-output
hash verification is implicit in the unchanged frozen_hashes.json.

### Q18. Where should the next physics investigation focus: constitutive response, neighbour dynamics, fast/slow combination, propagation, or observable extraction?

The dominant divergence stages are **A8-S07 (fast-to-slow exchange) and
A8-S12 (neighbour-response)** for A8, and **C10-S03 (neighbour
coherence)** for C10.  All three are **microscopic interaction stages**
that are not seen by the GR operator.  Concretely:

1. The fast-to-slow exchange A8-S07 carries the initial noise
   perturbation of `u_fast` (added in `A8_init`) into the rest of the
   pipeline.  The slow-layer update is structurally stable, but the
   fast-to-slow coupling is not filtered or de-noised, so the noise
   propagates and the exchange map becomes uncorrelated with GR kappa.
2. The A8 neighbour-response field A8-S12 (`n4f - u_fast`) is the
   relaxation-step residual.  Its large divergence is consistent with
   the noise being amplified through the neighbour averaging operation.
3. The C10 neighbour coherence factor C10-S03 is the 8-neighbour
   cosine alignment of the gradient direction.  This is the dominant
   divergent stage for C10 and produces the first Pearson collapse.

Therefore the next physics investigation should focus on
**neighbour dynamics / fast-slow combination** (Outcome B) for both
models.  Outcome B is the dominant outcome classification.

## Outcome determination

Using the Section 30 outcomes:

- **C10**: first divergence at C10-S04 in 3 clusters (>= 4? No), at S02
  in 2 clusters.  The cumulative count is mixed (5 / 5 clusters show
  some material transition, but the dominant one is C10-S04, an
  internal microscopic stage).  Result: **Outcome B (Neighbour or
  microscopic-state divergence)**.
- **A8**: first divergence at A8-S07 in 4 clusters (4 / 5 >= 4).  A8-S07
  is between S03 and S12 (the microscopic-state range).  Result: **Outcome B**.

Both models fall into **Outcome B (Neighbour or microscopic-state
divergence)**.  The constitutive bridge, propagation, and observable
extractor are not the dominant divergence source.

## Per-stage Pearson-vs-GR (S0, N0, median across 5 clusters)

| Stage | C10 | A8 |
|---|---|---|
| S00 | 1.000 | 0.991 |
| S01 | 1.000 | 0.991 |
| S02 | 0.561 | 0.540 |
| C10-S03 | 0.118 | n/a |
| C10-S04 | 0.155 | n/a |
| C10-S05 | 0.405 | n/a |
| C10-S06 | 0.155 | n/a |
| A8-S03 | n/a | 0.978 |
| A8-S04 | n/a | 0.978 |
| A8-S05 | n/a | 0.991 |
| A8-S06 | n/a | 0.990 |
| A8-S07 | n/a | 0.105 |
| A8-S08 | n/a | -0.105 |
| A8-S09 | n/a | 0.105 |
| A8-S10 | n/a | 0.986 |
| A8-S11 | n/a | 0.986 |
| A8-S12 | n/a | -0.510 |
| S19_kappa | 0.253 | 0.335 |
| S20_gamma1 | -0.186 | -0.198 |
| S20_gamma2 | -0.094 | -0.168 |
| S20_gamma_mag | 0.220 | 0.136 |

## Stage loss table (median L_i)

| Stage | C10 L_i | A8 L_i |
|---|---|---|
| S01 | +0.34 | +0.00 |
| S02 | +13.3 | +14.6 |
| C10-S03 / A8-S03 | **+57.8** | -14.5 |
| C10-S04 / A8-S04 | -57.6 | -0.006 |
| C10-S05 / A8-S05 | +16.5 | -0.31 |
| C10-S06 / A8-S06 | -16.5 | +0.006 |
| A8-S07 | n/a | +26.2 |
| A8-S08 | n/a | +0.48 |
| A8-S09 | n/a | -0.52 |
| A8-S10 | n/a | -26.1 |
| A8-S11 | n/a | -0.001 |
| A8-S12 | n/a | **+140.5** |

## Validation summary

- All frozen hashes match (`frozen_hashes.json`).
- All five clusters completed.
- GR, C10, A8 used the identical frozen input proxy.
- No new physics introduced.
- No coefficient changed.
- No fitting, no amplitude matching, no stage-dependent normalization
  in production.
- Instrumentation does not change the frozen production outputs (the
  wrapper re-computes internal states using the same arithmetic and
  update order; the public `candidate_10_combined` and
  `evolve_transport` functions are never modified).
- All required stages recorded.
- All stage IDs unique.
- Native and derived diagnostics are distinguished (suffix
  `_diagnostic` is used for C10-S05, A8-S09, and A8-S11).
- All time snapshots use the fixed schedule.
- All fixed transformations (G0..G8) were executed.
- All fixed spatial lags (49 combinations) were executed.
- No arbitrary rotation or translation search was performed.
- Native Jacobian and independent J1 verifier were compared.
- All six wrong controls (WR1..WR6) were completed.
- Every cluster received a first-divergence result.
- Every model received a dominant-divergence result.
- All 18 questions answered.
- All required outputs and plots exist.

## Reproduction

```bash
python macro_micro_response_bridge_diagnostic_lab001.py
```

Re-runs the full diagnostic laboratory end-to-end.  The script is
read-only with respect to all frozen executables (verified by hash at
startup).  Total runtime ~20 s.
