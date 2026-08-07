# PBUF SAME-INPUT LCDM/GR-C10-A8 BENCHMARK-LAB-001

**Apples-to-Apples Standard-Operator Comparison.**

Same-input Bridge Class D comparison: three lanes (L1 standard dimensionless GR operator, L2 frozen PBUF C10, L3 frozen PBUF A8/T1) receive the exact same frozen dimensionless cluster input `rho(x,y) = max(kappa_obs, 0) / max(max(kappa_obs, 0))` and are compared on the same common grid, mask, and statistics.

No fitting.  No parameter optimisation.  No microscopic-equation changes.

## Status

- Frozen hash verification: **PASS** (all seven frozen executables verified byte-identical to the LAB-FREEZE-001 and MICROSCOPIC-TRANSPORT-EQUIVALENCE-LAB-001 registries).
- Total runtime: **16.8 s**.
- Bridge class: **D** (dimensionless same-input operator comparison).
- L1 standard dimensionless GR operator: **ran** (padded and unpadded diagnostics).
- L2 frozen PBUF C10: **ran**.
- L3 frozen PBUF A8/T1: **ran**.
- Wrong controls WR1..WR5: **completed**.
- Aggregate classification: C10 = **G3**, A8/T1 = **G3**.

## Frozen laboratory

| Component | Frozen specification |
|---|---|
| Common input | `rho(x,y) = max(kappa_obs, 0) / max(max(kappa_obs, 0))` |
| L1 operator | Fourier-space Poisson solve + shear extraction |
| L2 response | Candidate 10 / Combined Local Response |
| L3 transport | A8 dual-layer + T1 scalar-density transport |
| Photons | 20 000 |
| Grid | 256² |
| Step | 0.03 |
| Steps | 160 |
| Bins | 64 |
| Smoothing S0 | native output |
| Smoothing S1 | Gaussian sigma = 1.0 comparison-grid pixel |
| L1 padding | mirror-pad 50% on each side, operator, crop |
| L1 unpadded diagnostic | no padding, periodic boundary |

## Preflight: bridge classification (Section 3)

All five clusters are **Bridge Class D**: the existing frozen dimensionless matter proxy `rho(x,y) = max(kappa_obs, 0) / max(max(kappa_obs, 0))` is reclassified as a valid common controlled input for operator-response comparison.  L1 runs on this same dimensionless input; comparisons against L0 (the observational map from which the proxy was derived) are labelled `conditional_same_source` and never claimed as independent predictive tests.

## Per-cluster Pearson (kappa vs L1, S0, C1)

| Cluster | C10 | A8/T1 |
|---|---|---|
| Abell2744 | +0.1420 | +0.2520 |
| MACS0416 | +0.2528 | +0.3645 |
| MACS1149 | +0.2109 | +0.2366 |
| AbellS1063 | +0.3273 | +0.3228 |
| Abell370 | +0.3677 | +0.4268 |
| **median** | **+0.2528** | **+0.3228** |

## Neighbourhood classification (per cluster)

| Cluster | C10 | A8/T1 |
|---|---|---|
| Abell2744 | N3 | N3 |
| MACS0416 | N3 | N3 |
| MACS1149 | N3 | N3 |
| AbellS1063 | N3 | N3 |
| Abell370 | N3 | N3 |
| **aggregate** | **G3** | **G3** |

## A8 improvement test (Section 22)

A8 improves over C10 in a cluster when at least 3 of 4 conditions hold: Delta_r > 0, Delta_D_NRMS < 0, Delta_D_Q < 0, Delta_D_P < 0.

| Cluster | Delta_r | Delta_D_NRMS | Delta_D_Q | Delta_D_P | n | A8 improves? |
|---|---|---|---|---|---|---|
| Abell2744 | +0.1099 | -0.0186 | +3.6988 | -0.0295 | 3 | yes |
| MACS0416 | +0.1117 | +0.0654 | -0.1498 | -0.0487 | 3 | yes |
| MACS1149 | +0.0257 | +0.0778 | +0.0574 | -0.1253 | 2 | no |
| AbellS1063 | -0.0045 | +0.0703 | -0.9744 | -0.0959 | 2 | no |
| Abell370 | +0.0591 | +0.0538 | -18.5105 | -0.0840 | 3 | yes |

A8 improves over C10 in **3 / 5** clusters.

## Wrong controls (Section 24)

| Control | Mean RMSE kappa |
|---|---|
| WR1_rotated_matter_for_L1 | 0.0777 |
| WR2_phase_scrambled_matter_for_L1 | 0.1159 |
| WR3_radially_symmetrized_matter_for_L1 | 0.0528 |
| WR4_mismatched_cluster | 0.0826 |
| WR5_uniform_matter_for_L1 | 1.0000 |

Expected behaviour:
- WR1: amplitude retained, correlation reduced.
- WR2: power spectrum retained, morphology destroyed.
- WR3: radial profile retained, asymmetric substructure removed.
- WR4: mismatched-cluster correlation lower than matched.
- WR5: uniform input produces ~zero shear in the bulk.

## Required questions (Section 25)

### Q1.  Did L1 run on the exact same frozen proxy used by C10 and A8?

Yes.  All three lanes receive `rho(x,y) = max(kappa_obs, 0) / max(max(kappa_obs, 0))` constructed once per cluster from the frozen observation FITS; the per-cluster SHA-256 of this proxy is recorded in `proxy_statistics.csv`.

### Q2.  Were all three lanes processed with identical grids, masks, smoothing, and statistics?

Yes (Section 6 apples-to-apples rule).  All lanes use the 64x64 common grid on [-8, 8] in pipeline units, the same valid-pixel mask, the same S0/S1 Gaussian smoothing (sigma = 1 comparison-grid pixel), and the same metric implementations.

### Q3.  Does C10 lie in the same, adjacent, related, or different operator neighbourhood relative to standard GR?

C10 aggregate classification: **G3**.  Per-cluster: N3, N3, N3, N3, N3.

### Q4.  Does A8/T1 lie in the same, adjacent, related, or different operator neighbourhood relative to standard GR?

A8/T1 aggregate classification: **G3**.  Per-cluster: N3, N3, N3, N3, N3.

### Q5.  Is either PBUF lane formally classified as N3/Mars?

C10 receives N3 in **5 / 5** clusters; A8/T1 receives N3 in **5 / 5** clusters.

### Q6.  Does A8 improve on C10 relative to the standard operator?

A8 improves over C10 in **3 / 5** clusters under the 3-of-4 condition test.  See the table above for the per-cluster breakdown.

### Q7.  Are the differences primarily amplitude differences, morphology differences, or both?

Both.  The PBUF lanes differ from L1 in both amplitude (RMS amplitude ratio generally < 1 because the PBUF matter input is dimensionless and the source-plane launch only populates the central field) and morphology (different peak positions, different multipole spectrum shape, different power-spectrum slope).  See `operator_pair_metrics.csv` and the figures.

### Q8.  Do the models agree more strongly in the core, middle, or outer radial regions?

Inspection of `radial_profiles.csv` shows that the PBUF lanes produce finite predictions only in the central ~10-15 of the 20 radial bins (r/r_max ~ 0.5-0.75).  In the central region where all three lanes have finite values, the agreement is qualitative (positive, near-zero).  In the outer bins, the PBUF lanes are NaN so no comparison is possible.

### Q9.  Do C10 and A8 reproduce similar convergence peaks to the standard operator?

No.  L1 (the standard operator on a smooth dimensionless source) produces a single broad peak at the cluster centre; the PBUF lanes produce 2-5 sharp peaks offset from the centre (see `peak_statistics.csv` and `plots/peak_comparison_*.png`).

### Q10.  Do C10 and A8 reproduce similar multipole structure?

No.  L1 multipoles drop steeply with m (e.g., for Abell 2744: |Q1|=25, |Q2|=23, |Q3|=10, |Q4|=2).  Both PBUF lanes produce a roughly flat |Q_m| spectrum.  See `multipole_statistics.csv` and `plots/multipole_comparison_*.png`.

### Q11.  Do C10 and A8 reproduce similar spatial power spectra?

Only partially.  The PBUF lanes share the broad-scale power with L1 (low-k ratio) but differ in mid- and high-k power.  See `power_spectrum_statistics.csv` and `plots/power_spectrum_comparison_*.png`.

### Q12.  Does common smoothing materially change the neighbourhood classification?

Smoothing (S1, sigma=1 pixel) tightens the per-cluster Pearson values slightly but does not move any cluster between N0/N1/N2/N3 in this run.  See `operator_pair_metrics.csv` (compare S0 vs S1 rows).

### Q13.  Does Fourier padding materially change the L1 comparison?

The padded-vs-unpadded L1 difference is reported in `padding_diagnostics.csv`.  Periodic-boundary effects are small for these cluster maps because the proxy is already near-zero at the field edges; the Pearson between padded and unpadded L1 maps is high (>0.99) in every cluster.

### Q14.  Do wrong controls behave as expected?

Yes (see the table above).  WR1 reduces correlation while preserving amplitude; WR2 destroys morphology; WR3 produces the most radially-symmetric response; WR4 mismatched-cluster correlation is the lowest of the matched-cluster metrics; WR5 uniform input produces zero shear in the bulk.

### Q15.  Do any independently generated residual ratios recur near alpha, 3alpha, or 6alpha?

The alpha audit (`fundamental_constant_audit.csv`) shows a 6alpha-dominant distribution of nearest multiples.  The PBUF matter input itself is derived from `kappa_obs`, so every entry is flagged `alpha_input_dependency = indirect` and the audit is **passive**.
Nearest-multiple counts: 6alpha=71, alpha=3, 3alpha=1.

### Q16.  Are the current PBUF outputs broadly in the conventional weak-lensing operator neighbourhood when tested apples to apples?

C10 aggregate = **G3**; A8/T1 aggregate = **G3**.  Under identical dimensionless input, the PBUF responses differ from the standard GR operator in both amplitude and morphology, and the multipole and power-spectrum distances are non-trivial.  See `Outcome determination` below for the formal interpretation.

## Outcome determination (Section 29)

C10 -> G3 -> Outcome D
A8/T1 -> G3 -> Outcome D

## Reproduction

```bash
python same_input_lcdm_gr_c10_a8_benchmark_lab001.py
```

Re-runs the full benchmark end-to-end (L1 padded + unpadded, L2 C10, L3 A8/T1, all 5 clusters, all metrics, all plots, all CSVs, registry append, validation, and report).  The script is read-only with respect to all frozen executables (verified by hash at startup).

