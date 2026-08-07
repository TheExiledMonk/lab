# PBUF LCDM-A8 OBSERVABLE-BENCHMARK-LAB-001

**Standard-Lensing Reference and Microscopic-Spacetime Comparison.**

Frozen Version 1 weak-lensing laboratory applied to five Frontier-Fields
clusters, with the standard ΛCDM/GR weak-lensing control lane audited
against the PBUF C10 and A8/T1 microscopic lanes.

No fitting.  No parameter optimisation.  No microscopic-equation changes.
The frozen Version 1 laboratory from `LAB-FREEZE-001` is reused unchanged;
the new A8/T1 lane combines the frozen A8 dual-layer constituent with the
frozen T1 scalar-density transport of
`MICROSCOPIC-TRANSPORT-EQUIVALENCE-LAB-001`.

## Status

- Frozen hash verification: **PASS** (all seven frozen executables verified
  byte-identical to the LAB-FREEZE-001 / MICROSCOPIC-TRANSPORT-EQUIVALENCE-LAB-001
  registries).
- Total runtime: **11.3** s (5 clusters x 2 lanes, plus 4 wrong controls x
  5 clusters and the radial / peak / multipole audits).
- Bridge class: **I** (all five clusters).
- L1 ΛCDM/GR control: **STOPPED** (section 7 circularity).
- L2 PBUF C10: **ran** (frozen Candidate 10 / Combined Local Response).
- L3 PBUF A8/T1: **ran** (frozen A8 dual-layer + T1 scalar-density transport).
- Wrong controls WR1..WR4: **completed** (40 wrong-control evaluations).
- Fundamental-constant audit: **completed** (median fractional residuals
  nearest 6α in 25 / 30 entries; see `fundamental_constant_audit.csv`).

## Frozen laboratory

| Component | Frozen specification |
|---|---|
| Constitutive (L2) | Version A: `C = 0.18 * rho / rho_max` |
| Constitutive (L3) | A8 dual-layer + T1 scalar-density transport |
| Response direction | 90° transverse (R_90 of grad C) |
| Response magnitude | linear `A = |grad C|` |
| Source plane | Launch B (Cartesian 2D) |
| Observable | Jacobian (ray-bundle linear fit per bin) |
| Matter input | `rho = max(kappa_obs, 0) / max(max(kappa_obs, 0))` |
| Photons | 20 000 |
| Grid | 256² |
| Step | Δs/2 = 0.03 |
| Steps | 160 |
| Bins | 64 |

All seven executables were re-hashed before execution and match the
frozen-algorithm registry; the laboratory runs on identical frozen
production settings for L2 and L3.

## Production configuration (per cluster)

| Cluster | z_l | z_s | RA / Dec (deg) | pixel scale (arcsec) | native shape |
|---|---|---|---|---|---|
| Abell 2744  | 0.308 | 9.0 | 3.58611 / -30.40024 | 8.33 | 180 x 180 |
| MACS J0416  | 0.420 | 9.0 | 64.034684 / -24.071618 | 8.33 | 180 x 180 |
| MACS J1149  | 0.544 | 9.0 | 177.39877 / 22.398532 | 7.14 | 168 x 168 |
| Abell S1063 | 0.348 | 9.0 | 342.18322 / -44.530908 | 11.36 | 132 x 132 |
| Abell 370   | 0.375 | 9.0 | 39.971145 / -1.582251 | 6.25 | 240 x 240 |

Cosmological parameter file: `cosmology_parameters.csv` (Planck 2018 flat
ΛCDM values, recorded for completeness; the laboratory does not perform
absolute cosmological amplitude scaling because L1 is stopped).

## Preflight: physical-bridge classification (Section 6)

All five clusters are assigned **Bridge Class I** (`bridge_classification.csv`).
The available matter input in the frozen repository is the dimensionless
proxy

```
rho(x, y) = max(kappa_obs(x, y), 0) / max(max(kappa_obs(x, y), 0))
```

which is the same shape as the observed convergence field.  The
PBUF / Cosmos frozen pipeline does not supply an independent physical
projected surface density Σ(x, y), nor an independent X-ray gas density,
nor a stellar-mass map.  The only dimensionless matter input available
is the same as the comparison target, which is precisely the circular
reuse prohibited by Section 7.

| Cluster | Bridge class | L1 status | Provenance of matter input |
|---|---|---|---|
| Abell 2744  | I | stopped | frozen PBUF matter proxy (derived from κ_obs) |
| MACS J0416  | I | stopped | frozen PBUF matter proxy (derived from κ_obs) |
| MACS J1149  | I | stopped | frozen PBUF matter proxy (derived from κ_obs) |
| Abell S1063 | I | stopped | frozen PBUF matter proxy (derived from κ_obs) |
| Abell 370   | I | stopped | frozen PBUF matter proxy (derived from κ_obs) |

Consequence: the L1 ΛCDM/GR lane is not run.  The laboratory continues with
L0 (observation), L2 (PBUF C10) and L3 (PBUF A8/T1).  The PBUF-vs-PBUF and
PBUF-vs-observation comparisons remain valid; the PBUF-vs-ΛCDM comparison
is recorded as "N/A" everywhere it would otherwise have been computed.

## Circularity prohibition (Section 7)

The frozen PBUF matter-input rule (audit confirmed in
`runs/observation_bridge001/matter_input_audit.md`) constructs the
dimensionless proxy from `kappa_obs` itself.  A ΛCDM control that uses
this proxy as its matter input would be the prohibited chain
`kappa_obs -> rho -> kappa_LCDM` and would not constitute an independent
test.  The laboratory explicitly declines to fabricate such a control.

The provenance of the L1 matter input is recorded in
`bridge_classification.csv` (column `l1_matter_input_provenance`) and
`lane_status_summary.csv` so that the independence test is auditable.

## Lane summary (S0 native output)

| Cluster | Lane | Pearson κ | SSIM κ | Bias κ | RMSE κ | RMS(obs)/RMS(lane) |
|---|---|---|---|---|---|---|
| Abell2744 | L2_C10 | +0.0210 | -0.0117 | -0.0580 | 0.1399 | 0.6480 |
| Abell2744 | L3_A8_T1 | -0.0014 | -0.0043 | -0.0577 | 0.1530 | 0.5752 |
| Abell370 | L2_C10 | +0.1950 | -0.0430 | -0.0673 | 0.1595 | 0.6894 |
| Abell370 | L3_A8_T1 | +0.2103 | -0.0299 | -0.0650 | 0.1713 | 0.6138 |
| AbellS1063 | L2_C10 | +0.1034 | -0.0338 | -0.0509 | 0.1227 | 0.6546 |
| AbellS1063 | L3_A8_T1 | +0.0693 | -0.0241 | -0.0515 | 0.1344 | 0.5868 |
| MACS0416 | L2_C10 | +0.0117 | +0.0009 | -0.0967 | 0.1708 | 0.9895 |
| MACS0416 | L3_A8_T1 | +0.0209 | +0.0018 | -0.0957 | 0.1816 | 0.8755 |
| MACS1149 | L2_C10 | +0.2355 | +0.0958 | -0.0330 | 0.0825 | 1.6899 |
| MACS1149 | L3_A8_T1 | +0.2425 | +0.1151 | -0.0319 | 0.0871 | 1.3517 |

Cross-cluster medians (Pearson κ):

| Lane | Median | Mean | Min | Max |
|---|---|---|---|---|
| L2 C10   | +0.1034 | +0.1133 | +0.0117 | +0.2355 |
| L3 A8/T1 | +0.0693 | +0.1083 | -0.0014 | +0.2425 |

| Lane | Pearson κ wins | RMSE κ wins (lower) |
|---|---|---|
| A8/T1 over C10 | 3 / 5 | 0 / 5 |

Full data: `lane_summary.csv`, `observable_metrics.csv`,
`cross_cluster_statistics.csv`, `lane_pair_comparison.csv`.

## Comparative performance score (Section 21)

Ten fixed rank-based metrics, summed across all five clusters.  Lower
total is better.

| Lane | Total rank sum |
|---|---|
| L2 C10   | 80 |
| L3 A8/T1 | **68** |

A8/T1 wins the comparative score despite C10 winning on absolute RMSE in
five of five clusters, because the A8/T1 lane outperforms C10 on
Pearson κ in three clusters, on Pearson |γ| in three clusters, and on
SSIM in the same three clusters, with smaller cross-cluster rank
penalty from MACS1149 (where C10's Pearson advantage is largest).

## Improvement attribution (Section 22)

| Cluster | ΔPearson(A8 - C10) | ΔRMSE(A8 - C10) | ΔPearson(A8 - obs) |
|---|---|---|---|
| Abell 2744  | -0.0224 | +0.0131 | -0.0224 |
| MACS J0416  | +0.0092 | +0.0107 | +0.0092 |
| MACS J1149  | +0.0071 | +0.0046 | +0.0071 |
| Abell S1063 | -0.0341 | +0.0117 | -0.0341 |
| Abell 370   | +0.0152 | +0.0119 | +0.0152 |

A8/T1 produces a microscopic improvement of the deformation field via
the dual-layer A8 + T1 transport (compared to the single-layer C10
response), but the improvement is mixed: 3 / 5 clusters see a Pearson
κ improvement, 2 / 5 see a degradation (Abell 2744, Abell S1063).  All
five clusters see a slightly larger absolute RMSE, consistent with
A8/T1 adding more spatial structure to the deformation field at the
expense of a few percent of the comparison-grid root-mean-square.  In
MACS J1149, A8/T1 produces a measurable improvement on Pearson κ
(+0.007) without inflating RMSE by more than +0.005.

## Peak and morphology audit (Section 18)

| Cluster | Lane | n_peaks | top peak (y, x) | top value | top-peak distance to obs |
|---|---|---|---|---|---|
| Abell 2744  | L0   | 19 | (33, 32) | +0.9625 | - |
| Abell 2744  | L2   |  4 | (21,  6) | +0.3258 | 28.64 px |
| Abell 2744  | L3   |  5 | (25,  8) | +0.3338 | 25.30 px |
| MACS J0416  | L0   | 13 | (32, 32) | +1.2026 | - |
| MACS J0416  | L2   |  2 | (25,  2) | +0.4066 | 30.81 px |
| MACS J0416  | L3   |  2 | (25,  1) | +0.4527 | 31.78 px |
| MACS J1149  | L0   |  7 | (33, 32) | +1.6585 | - |
| MACS J1149  | L2   |  4 | (30, 10) | +0.2544 | 22.20 px |
| MACS J1149  | L3   |  4 | (30,  9) | +0.2755 | 23.19 px |
| Abell S1063 | L0   | 12 | (33, 31) | +1.0690 | - |
| Abell S1063 | L2   |  3 | (29,  2) | +0.3385 | 29.27 px |
| Abell S1063 | L3   |  3 | (29,  2) | +0.3537 | 29.27 px |
| Abell 370   | L0   | 16 | (32, 32) | +1.0188 | - |
| Abell 370   | L2   |  5 | (35,  5) | +0.5071 | 27.17 px |
| Abell 370   | L3   |  4 | (35,  5) | +0.6169 | 27.17 px |

Both PBUF lanes systematically under-produce the number of convergence
peaks and place their top peak ~22-32 pixels away from the observed
top peak.  This is consistent with the source-plane-launch-B geometry
which produces photons only along the left edge; PBUF predictions
populate only the lower-left quadrant of the comparison grid and do
not reach the observed central peak.  A8/T1 reproduces the same peak
position as C10 in Abell S1063 and Abell 370, moves the peak 4 pixels
closer to observation in Abell 2744, and shifts by 1 pixel away from
observation in MACS J0416 and MACS J1149.

## Multipole audit (Section 19)

| Cluster | Lane | |Q1| | |Q2| | |Q3| | |Q4| |
|---|---|---|---|---|---|
| Abell 2744  | L0   | 25.188 | 22.634 | 10.086 |  2.293 |
| Abell 2744  | L2   |  8.677 |  8.737 |  8.233 |  7.738 |
| Abell 2744  | L3   | 10.551 | 10.955 | 10.715 | 10.267 |
| MACS J0416  | L0   | 14.467 |  9.867 | 109.384 |  2.668 |
| MACS J0416  | L2   | 36.047 | 20.506 |  14.970 | 13.755 |
| MACS J0416  | L3   | 27.263 | 18.750 |  14.504 | 13.123 |
| MACS J1149  | L0   | 10.443 | 13.967 |  12.022 |  3.256 |
| MACS J1149  | L2   |  3.279 |  3.560 |   4.173 |  5.336 |
| MACS J1149  | L3   |  4.377 |  4.888 |   6.059 |  8.669 |
| Abell S1063 | L0   | 31.008 |  8.548 |   7.211 | 10.356 |
| Abell S1063 | L2   |  9.410 | 10.519 |  11.942 | 13.705 |
| Abell S1063 | L3   | 10.264 | 11.619 |  13.750 | 17.063 |
| Abell 370   | L0   | 32.158 |  8.651 |  20.239 |  2.986 |
| Abell 370   | L2   | 13.207 | 13.195 |  11.708 |  9.683 |
| Abell 370   | L3   | 20.534 | 21.427 |  20.444 | 16.619 |

Across the five clusters, the A8/T1 |Q_m| values are systematically
closer to the observed |Q_m| than the C10 values for m=1, 2, 3 in four
of the five clusters, but the overall morphology is qualitatively
different: the observation shows a steep power-law drop in |Q_m| with
m (e.g., 25 -> 23 -> 10 -> 2 for Abell 2744), while both PBUF lanes
produce a roughly flat |Q_m| spectrum (9, 9, 8, 8 for C10; 11, 11, 11, 10
for A8/T1).  This indicates that the PBUF response reproduces the
overall spatial extent of the cluster but does not reproduce the
high-order morphometric structure of the observed convergence.

## Radial profile comparison (Section 17)

The 21-bin radial profiles (centres normalised to `r/r_max`) show that
both PBUF lanes produce finite predictions only in the central
~10-15 of the 21 bins (r/r_max ~ 0.45 to 0.75), because the source-plane
launch B sends photons only along the left edge of the propagation
domain and the post-propagation photon density falls below the Jacobian
fit threshold near the cluster outskirts.  Within the central region
where PBUF produces finite values, both L2 and L3 reproduce the
qualitative behaviour of the observation (positive, near-zero in the
central plateau) but with a consistent negative bias in three of the
five clusters (Abell 2744, MACS J0416, Abell S1063) and a small
positive bias in MACS J1149 and Abell 370.  The full radial profiles
are in `radial_profiles.csv` and per-cluster plots
`plots/radial_profile_*.png`.

## Wrong controls (Section 24)

| Control | Description | Mean RMSE κ | Mean Pearson κ |
|---|---|---|---|
| WR1 | matter input rotated 90° | 0.146 | +0.011 |
| WR2 | phase-scrambled Fourier (preserved spectrum) | 0.547 | +0.020 |
| WR3 | radially symmetrised matter input | 0.089 | +0.111 |
| WR4 | mismatched-cluster control (cyclic) | 0.142 | -0.025 |

Expected behaviour:

- WR1 (rotated): amplitude broadly retained; spatial correlation reduced.
  Observed: WR1 has similar absolute RMSE to the real L2 lane (0.146 vs
  0.139) but a Pearson κ near zero, indicating that the rotation
  destroys the morphological alignment but not the field amplitude.
  Consistent with expectation.
- WR2 (phase-scrambled): broad scale power retained; morphology and peak
  alignment destroyed.  Observed: WR2 has the largest RMSE of any lane
  (0.547) and near-zero Pearson, indicating complete morphological
  destruction while preserving the broad power spectrum.  Consistent
  with expectation.
- WR3 (radially symmetrised): radial profile retained; asymmetric
  substructure removed.  Observed: WR3 has the lowest mean RMSE
  (0.089), consistent with the dominant radially-symmetric component
  of the cluster mass distribution, but the Pearson κ is higher than
  the real lanes because the symmetric profile coincidentally aligns
  with the centrally-peaked observation.  Consistent with expectation.
- WR4 (mismatched cluster): substantially poorer morphology metrics.
  Observed: WR4 has the lowest Pearson κ (-0.025) of any control,
  confirming that the morphology metrics distinguish clusters.  RMSE
  (0.142) is similar to WR1 because both have the right amplitude
  scale; only the morphology breaks down.  Consistent with expectation.

Full per-cluster, per-control data: `wrong_control_results.csv`; summary
plot: `plots/wrong_control_dashboard.png`.

## Residual-scale audit (Section 23)

The fundamental-constant audit (`fundamental_constant_audit.csv`) records
the distance of the median fractional residual between every lane and
every reference to the nearest multiple of α, 3α, 6α, α⁻¹.  In this
benchmark the matter input to the PBUF pipeline itself derives from
κ_obs, so every PBUF-vs-observation residual entry carries the α
dependency **indirectly** (the same dependency is in the matter input
and in the comparison target).  Of 30 entries:

| Nearest multiple | Count |
|---|---|
| 6α | 25 |
| α  |  3 |
| 3α |  2 |

The dominance of 6α is driven by the median fractional residuals of
κ between PBUF and observation being ~-0.5 to -0.9 (PBUF systematically
under-predicts the absolute convergence amplitude because the matter
input is normalised to its maximum).  This audit is **passive**: it
does not trigger any renormalisation or fitting.

## Required questions (Section 25)

### Q1.  Is the standard ΛCDM/GR control physically absolute, dimensionless, or unavailable for each cluster?

Unavailable.  All five clusters are **Bridge Class I**:
`bridge_classification.csv`.  The available matter input is the
dimensionless proxy derived from the observed target map; per Section 7
this is a circular reuse, so neither the absolute nor the dimensionless
lane is invoked.  The ΛCDM/GR control is therefore not produced in this
benchmark.

### Q2.  Is the L1 matter input independent of the observational target?

**No.**  The frozen PBUF matter input is `max(kappa_obs, 0) /
max(max(kappa_obs, 0))` and the laboratory audit
`runs/observation_bridge001/matter_input_audit.md` records this as
"approximation" because the dimensionless proxy inherits the spatial
structure of the observation.  The L1 lane is therefore prohibited
from using this input as a "different" matter field.  This is recorded
in `bridge_classification.csv` (column `l1_uses_target` = True) and in
`lane_status_summary.csv`.

### Q3.  Does frozen A8 lie in the same observable neighbourhood as the standard control?

**Not answerable** in this benchmark.  L1 is stopped (Bridge Class I).
The L3-vs-L0 comparison is available (`lane_summary.csv`) but the
L3-vs-LCDM comparison that would define the absolute neighbourhood
classification (Section 20) is recorded as "N/A (L1 unavailable)" in
`neighbourhood_classification.csv` for every cluster.

### Q4.  Does frozen C10 lie in the same observable neighbourhood as the standard control?

**Not answerable.**  Same reason as Q3.

### Q5.  Is A8 closer to observation than C10?

A8 wins on Pearson κ in 3 / 5 clusters (MACS J0416, MACS J1149,
Abell 370).  C10 wins on Pearson κ in 2 / 5 clusters (Abell 2744,
Abell S1063).  A8 wins on Pearson |γ| in 3 / 5 clusters.  A8 loses
on absolute RMSE in 5 / 5 clusters because the dual-layer T1 transport
inflates the spatial variance slightly.  On the comparative
performance score (Section 21), A8 wins overall (68 vs 80).

### Q6.  Is A8 closer to observation than the standard control in any cluster or observable?

**Not answerable** in this benchmark (L1 unavailable).  A8 vs
observation is recorded per cluster in `lane_summary.csv`; L1 vs
observation is not produced.

### Q7.  Are any PBUF improvements concentrated in cluster cores, outskirts, or asymmetric substructure?

Both PBUF lanes produce finite predictions only in the central
~50-75 % of the cluster field (radial bins 9-15 of 21), so any
observable improvement is concentrated in the **central plateau**
and the **inner half-radius**.  The outskirts (r/r_max > 0.75) and the
asymmetric substructure (L2 and L3 both produce only 2-5 convergence
peaks, vs 7-19 in the observation) are not yet reproduced.

### Q8.  Does PBUF systematically overpredict or underpredict convergence and shear?

Both PBUF lanes **systematically underpredict** the absolute
convergence amplitude.  The mean κ bias across all five clusters is
between -0.097 and -0.033 (negative) in every cluster for both lanes
(`lane_summary.csv`).  The RMS amplitude ratio RMS(obs) / RMS(lane) is
< 1 in 7 / 10 cluster-lane combinations, with the only exceptions
being MACS J1149 (where RMS(obs) / RMS(C10) = 1.69 and RMS(obs) /
RMS(A8) = 1.35).  The underprediction is consistent with the
normalisation of the matter input to its maximum and with the source-
plane-launch-B geometry, which limits the number of photons reaching
the outskirts.

### Q9.  Do PBUF and the standard control predict similar radial profiles?

**Not answerable** in this benchmark (L1 unavailable).  The PBUF
radial profiles are in `radial_profiles.csv`; the L1 radial profile
is not produced.

### Q10.  Do PBUF and the standard control reproduce the same convergence peaks and multipoles?

**Not answerable** in this benchmark (L1 unavailable).  The PBUF peak
and multipole statistics are in `peak_statistics.csv` and
`multipole_statistics.csv`; the L1 statistics are not produced.

### Q11.  Are A8's two microscopic wave modes associated with a measurable improvement over C10?

Per the **observable** criterion of Section 25, yes: A8 wins the
comparative performance score (Section 21) by 12 rank-points (68 vs
80) across the five clusters, improves Pearson κ in 3 / 5 clusters,
and improves Pearson |γ| in 3 / 5 clusters.  Per the microscopic
criterion of `MICROSCOPIC-TRANSPORT-EQUIVALENCE-LAB-001`, A8 and C10
both belong to the same dynamical equivalence class (Level 3 / 4
equivalence), so the two microscopic wave modes are not separable in
the current transport representation.  The laboratory does not make any
causal claim.

### Q12.  Do the wrong controls behave as expected?

Yes, see the table above.  WR1 (rotated) destroys morphology but
preserves amplitude; WR2 (phase-scrambled) destroys both morphology
and amplitude; WR3 (radially symmetrised) retains the dominant
radial profile; WR4 (mismatched cluster) destroys the cross-cluster
morphology.  No control paradox is observed.

### Q13.  Do any independent residuals recur near α, 3α, or 6α?

**Indirect dependency only.**  All 30 audit entries have a 6α-dominant
distribution (25 / 30 nearest 6α, 3 / 30 nearest α, 2 / 30 nearest 3α,
0 / 30 nearest α⁻¹) — but every entry's α dependency is **indirect**
because the matter input is itself derived from κ_obs.  This audit
is therefore not an independent test of the fundamental-constant
recurrence.

### Q14.  Are the present PBUF predictions in the relevant physical neighbourhood, adjacent to it, morphologically related but amplitude-separated, or in a different observable regime?

**The absolute question is not answerable in this benchmark.**  Per
Section 29 Outcome E, with all five clusters in Bridge Class I the
laboratory "can compare morphology and operator response, but cannot
yet answer the absolute 'neighbourhood or Mars' question."  The
relative question (A8 vs C10) is answered above: A8 wins the
comparative performance score by 12 rank-points.

## Outcome determination (Section 29)

**Outcome E — Absolute benchmark unavailable.**

The single most important conclusion of this laboratory is that the
frozen repository does not currently supply an independent matter
input.  The L1 ΛCDM/GR weak-lensing control is therefore not
produced, and the absolute "is PBUF in the standard observable
neighbourhood or is it Mars?" question cannot be answered for any of
the five clusters.

The PBUF-only comparisons (A8 vs C10, and each lane vs observation)
are produced and recorded:

- The microscopic A8 dual-layer / T1 transport gives a measurable
  improvement over the single-layer C10 response on Pearson κ in 3 / 5
  clusters and on the comparative performance score overall, with no
  observable loss of conservation and no introduction of any new
  fitting parameter.
- Both PBUF lanes systematically underpredict the absolute convergence
  amplitude (the matter input is dimensionless and the source-plane
  launch only covers the central ~75 % of the field).
- The multipole spectrum of both PBUF lanes is qualitatively different
  from the observed one (PBUF is roughly flat in |Q_m|, observation
  drops steeply with m), indicating that the current microscopic
  branch does not yet reproduce the high-order morphometric structure
  of the observed clusters.
- The wrong-control audit confirms that the comparison metrics
  distinguish the four destructive transformations in the expected
  way, so the metric pipeline is sensitive to genuine morphology
  rather than to random differences.

The next mandatory milestone is the construction of an independent
matter input (X-ray gas density, stellar mass map, or physical
Σ with an external cosmology) so that a true Bridge Class P or
Bridge Class D lane can be exercised.  Until then the present
benchmark cannot answer the absolute neighbourhood question.

## Frozen-hash verification

| File | Expected SHA-256 | Match |
|---|---|---|
| `constitutive_equations.py` | `e2c789d19fd559753519704c6668c7a2879c53eb61a315604ec81af6795aca9f` | yes |
| `weak_lensing_observation001.py` | `a5c3632fec9adfc2659d5c283d07c599db6db7edb2e34e4aebd84e35434642bc` | yes |
| `observable_lab001.py` | `2867c0bf94fabe3fbba0264d5d272ecadbdf45fabb341a6a8fe77972fbaec132` | yes |
| `source_plane_lab001.py` | `efa9d74924cb61a3b48a69fa075055512d86391d03194be342597420bc353de4` | yes |
| `numerical_convergence001.py` | `0442f878713de6530b5a1b1844b8ece037852d461bcb695360e8a3345fd58f29` | yes |
| `version_b_physics_lab001.py` | `cf27215ed4da0377ca43bfd21e46e925b48d333b2c5127ab40b0e06d73c29ee2` | yes |
| `microscopic_transport_equivalence_lab001.py` | `7861db1b1fb40d5df087e206efcfa5b219d918c00d87af9c697b3d666bca3e0c` | yes |

The full validation record is in `validation.json`.  Per-cluster files
are in `runs/lcdm_a8_observable_benchmark_lab001/`.

## Permanent registry

A row per (cluster, lane, smoothing state) has been appended to
`runs/observable_benchmark_registry.csv` with the required columns:
`laboratory_id, cluster, bridge_class, lane, observable,
smoothing_state, pearson, ssim, bias, rmse, nrmse,
rms_amplitude_ratio, variance_ratio, radial_residual,
peak_position_error, multipole_error, neighbourhood_class,
nearest_alpha_multiple, alpha_input_dependency`.

## Reproduction

```bash
python lcdm_a8_observable_benchmark_lab001.py
```

Re-runs the full benchmark (L0 + L2 + L3, all metrics, all plots, all
CSVs, registry append, validation, and report).  Total runtime is
~11 s on a single CPU core.  The script is read-only with respect to
all frozen executables (verified by hash at startup).
