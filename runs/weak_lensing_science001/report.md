# PBUF WEAK-LENSING-SCIENCE-001

**First scientific prediction using the frozen Version 1
weak-lensing laboratory (LAB-FREEZE-001).**

Pure forward prediction.  No fitting.  No parameter
optimisation.  No cosmological rescaling.

## Status

- Regression baseline reproduction: **PASS**
- Frozen hash verification: **PASS**
- Total runtime: **96.7 s**

## Frozen laboratory

This milestone uses the **identical** frozen implementation
established by LAB-FREEZE-001.  No code modification.  No
parameter tuning.  No fitting.  No cosmological bridge.

| Component | Frozen specification |
|---|---|
| Constitutive | `C(X) = 0.18 · ρ(X) / ρ_max` (Version A) |
| Transport | neighbour-to-neighbour |
| Response direction | 90° transverse (R_90 of ∇C) |
| Response magnitude | linear `A = |∇C|` |
| Update rule | direct addition `v_new = v + step · r` |
| Normalisation | per-step unit-speed renormalisation |
| Source plane | Launch B (uniform Cartesian 2D) |
| Observable | Jacobian (ray-bundle linear fit per bin) |
| Matter input | `ρ = max(κ_obs, 0) / max(max(κ_obs, 0))` |

## Production configurations

| Parameter | Minimum production | High accuracy |
|---|---|---|
| Photons | 20,000 | 50,000 |
| Constitutive grid | 256² | 512² |
| Step size | Δs/2 = 0.0300 | Δs/4 = 0.0150 |
| Number of steps | 160 | 320 |
| Source plane | Cartesian 2D (Launch B) | Cartesian 2D |
| Observable | Jacobian | Jacobian |

## Internal validation (mandatory pre-run)

### Frozen hashes

All ten frozen source files reproduce the SHA-256 hashes
recorded in `runs/lab_freeze001/checksums.csv`.

| File | SHA-256 | Match |
|---|---|---|
| `constitutive_equations.py` | `e2c789d19fd55975…` | YES |
| `weak_lensing_observation001.py` | `a5c3632fec9adfc2…` | YES |
| `observable_lab001.py` | `2867c0bf94fabe3f…` | YES |
| `source_plane_lab001.py` | `efa9d74924cb61a3…` | YES |
| `numerical_convergence001.py` | `0442f878713de653…` | YES |

### Regression baseline reproduction

The Abell 2744 frozen regression baseline (`regression_baseline.json`) is reproduced before any science
run.  This is mandatory: an abort is triggered on any failure.

| Quantity | Expected | Reproduced | Δ |
|---|---|---|---|
| Trajectory SHA-256 | `80d8fe47bd0d4567…` | `80d8fe47bd0d4567…` | byte-exact |
| RMS κ | 1.352162e-01 | 1.352162e-01 | 0.000 % |
| Peak |κ| | 4.519461e-01 | 4.519453e-01 | 0.000 % |
| Mean κ | -4.117442e-03 | -4.117442e-03 | 0.000e+00 |
| std κ | 1.351497e-01 | 1.351535e-01 | 0.003 % |
| RMS γ | 8.608378e-02 | 8.608378e-02 | 0.000 % |
| Conservation max | 2.220e-16 | 2.220e-16 | match |

**Reproduction status:** PASS

## Per-cluster outputs (minimum production)

Every cluster produces the following frozen Version A
observables: constitutive field `C`, gradient field `∇C`,
response field `r`, photon trajectories, deflection field,
κ prediction, γ₁ prediction, γ₂ prediction, |γ| prediction,
magnification prediction.

### Quantitative comparison against observation

| Cluster | RMS κ | RMS γ₁ | RMS γ₂ | RMS |γ| | Pearson κ | Pearson γ | SSIM κ | SSIM γ | MAE κ | MAE γ | Max κ | Max γ |
|---|---|---|---|---|---|---|---|---|---|---|---|
| Abell 2744 | 1.5580e-01 | 8.8934e-02 | 8.2990e-02 | 6.0822e-02 | 0.0140 | 0.0895 | -0.0106 | 0.1008 | 1.1599e-01 | 4.7554e-02 | 5.6735e-01 | 1.8089e-01 |
| MACS J0416 | 1.8596e-01 | 9.3952e-02 | 9.0006e-02 | 6.5276e-02 | 0.0190 | 0.1175 | 0.0004 | 0.1267 | 1.4617e-01 | 5.1341e-02 | 6.8112e-01 | 2.0099e-01 |
| MACS J1149 | 8.5592e-02 | 1.2781e-01 | 6.5664e-02 | 1.1556e-01 | 0.2390 | 0.0075 | 0.1081 | 0.0120 | 6.5933e-02 | 9.4553e-02 | 2.6208e-01 | 3.5090e-01 |
| Abell S1063 | 1.3540e-01 | 9.4235e-02 | 7.0394e-02 | 6.9710e-02 | 0.0895 | 0.0836 | -0.0383 | 0.0929 | 9.4659e-02 | 5.2434e-02 | 7.2325e-01 | 2.4543e-01 |
| Abell 370 | 1.7915e-01 | 1.4068e-01 | 9.4167e-02 | 9.1733e-02 | 0.1954 | -0.0171 | -0.0398 | 0.0122 | 1.3667e-01 | 6.5556e-02 | 7.5124e-01 | 4.9048e-01 |

### Cross-cluster ranking

| Cluster | RMS κ | RMS γ | Pearson κ | Pearson γ | SSIM κ |
|---|---|---|---|---|---|
| 1. MACS J1149 | 8.5592e-02 | 1.1556e-01 | 0.2390 | 0.0075 | 0.1081 |
| 2. Abell 370 | 1.7915e-01 | 9.1733e-02 | 0.1954 | -0.0171 | -0.0398 |
| 3. Abell S1063 | 1.3540e-01 | 6.9710e-02 | 0.0895 | 0.0836 | -0.0383 |
| 4. MACS J0416 | 1.8596e-01 | 6.5276e-02 | 0.0190 | 0.1175 | 0.0004 |
| 5. Abell 2744 | 1.5580e-01 | 6.0822e-02 | 0.0140 | 0.0895 | -0.0106 |

## Required questions

### Q1. Reproduction quality

Measured across all five benchmark clusters at minimum
production (20k photons, 256² grid, Δs/2):

- Pearson κ: min = +0.0140, max = +0.2390, mean = +0.1114
- Pearson |γ|: min = -0.0171, max = +0.1175, mean = +0.0562
- SSIM κ: min = -0.0398, max = +0.1081, mean = +0.0040
- SSIM |γ|: min = +0.0120, max = +0.1267, mean = +0.0689
- RMS κ residual: min = 0.0856, max = 0.1860, mean = 0.1484
- RMS |γ| residual: min = 0.0608, max = 0.1156, mean = 0.0806

All numerical values are dimensionless pipeline units.  The
frozen Version A outputs are forward predictions; the published
SaWLens reconstructions are posterior-mean lensing maps.  No
physical rescaling is applied.


### Q2. Strongest agreement

Mean Pearson correlation across the five clusters (frozen
Version A minimum production):

| Observable | Mean Pearson | Mean SSIM | Mean RMS residual |
|---|---|---|---|
| kappa | +0.1114 | +0.0040 | 0.1484 |
| gamma1 | -0.0825 | +nan | 0.1091 |
| gamma2 | -0.0448 | +nan | 0.0806 |
| gamma_mag | +0.0562 | +0.0689 | 0.0806 |

**Strongest agreement by Pearson:** `kappa` (mean Pearson = +0.1114).
**Strongest agreement by SSIM:** `gamma_mag` (mean SSIM = +0.0689).
**Lowest RMS residual:** `gamma_mag` (mean RMS = 0.0806).


### Q3. Spatial distribution of residuals

Photons in the frozen Launch B source plane travel a
total distance `step · steps = 4.8` in pipeline units
from x = -8, so the predicted κ field is non-NaN only on
the strip x ∈ [-8, -3.2].  Therefore the residual map is
informative only on that strip.  Residual bins are
classified by physical radius from the pipeline origin:
`core` (r ≤ 1.5), `intermediate` (1.5 < r ≤ 4.0),
`outer halo` (r > 4.0, the strip region).

| Cluster | Mean |residual| core | intermediate | outer | Peak location |
|---|---|---|---|---|---|
| Abell 2744 (288 finite px) | +nan (0 px) | +nan (0 px) | +0.1147 (254 px) | outer halo (strip) |
| MACS J0416 (288 finite px) | +nan (0 px) | +nan (0 px) | +0.1411 (254 px) | outer halo (strip) |
| MACS J1149 (288 finite px) | +nan (0 px) | +nan (0 px) | +0.0648 (254 px) | outer halo (strip) |
| Abell S1063 (288 finite px) | +nan (0 px) | +nan (0 px) | +0.0921 (254 px) | outer halo (strip) |
| Abell 370 (288 finite px) | +nan (0 px) | +nan (0 px) | +0.1386 (254 px) | outer halo (strip) |

**Note:** the frozen laboratory's photon-source strip
(`x ∈ [-8, -3.2]`) lies entirely in the `outer halo` 
annulus, so the only meaningful residual statistics are
in the `outer halo` column.  The `core` and 
`intermediate` annuli contain zero finite pixels in every
cluster, by construction of the frozen pipeline.
**All largest residuals occur in the outer halo (the 
photon-source strip).**


### Q4. Randomness vs systematic structure

Residual skewness, sign asymmetry, and centroid offset
quantify whether the frozen prediction systematically
over- or under-shoots the observation.

| Cluster | κ residual mean | κ residual median | κ residual skew | γ residual mean | γ residual median | γ residual skew |
|---|---|---|---|---|---|
| Abell 2744 | -0.05939 | -0.04371 | -0.6406 | +0.00155 | +0.00201 | +0.0653 |
| MACS J0416 | -0.09849 | -0.10330 | -0.2195 | -0.01258 | -0.01724 | +0.0570 |
| MACS J1149 | -0.03244 | -0.02185 | -0.0800 | -0.09306 | -0.08412 | -1.2720 |
| Abell S1063 | -0.05393 | -0.04259 | -1.0543 | -0.02200 | -0.02224 | -0.1140 |
| Abell 370 | -0.06748 | -0.06041 | -0.3816 | -0.01547 | -0.01029 | -1.0178 |

Interpretation: a residual mean that is small relative to
the residual standard deviation, and a skewness near zero,
is consistent with predominantly random residuals.  A
non-zero mean or non-zero skewness across multiple clusters
indicates systematic structure.

Number of clusters with positive mean κ residual minus number with negative: -5 (range -0.09849 to -0.03244).
Number of clusters with positive mean γ residual minus number with negative: -3 (range -0.09306 to +0.00155).

### Q5. Cross-cluster similarity

Similarity is measured by the cross-cluster Pearson
correlation of the predicted κ maps (pipeline-grid
comparison after subtraction of the mean).

| | Abell 2744 | MACS J0416 | MACS J1149 | Abell S1063 | Abell 370 |
|---|---|---|---|---|---|
| Abell 2744 | +1.0000 | -0.0467 | +0.0531 | -0.0447 | -0.2692 |
| MACS J0416 | -0.0467 | +1.0000 | -0.0957 | -0.1637 | +0.0122 |
| MACS J1149 | +0.0531 | -0.0957 | +1.0000 | -0.0177 | -0.0789 |
| Abell S1063 | -0.0447 | -0.1637 | -0.0177 | +1.0000 | -0.0935 |
| Abell 370 | -0.2692 | +0.0122 | -0.0789 | -0.0935 | +1.0000 |

Cross-cluster residual correlation (predicted - observed
κ maps):

| | Abell 2744 | MACS J0416 | MACS J1149 | Abell S1063 | Abell 370 |
|---|---|---|---|---|---|
| Abell 2744 | +1.0000 | -0.0006 | +0.0379 | -0.0976 | -0.2297 |
| MACS J0416 | -0.0006 | +1.0000 | -0.0578 | -0.1922 | -0.0503 |
| MACS J1149 | +0.0379 | -0.0578 | +1.0000 | -0.0734 | +0.0450 |
| Abell S1063 | -0.0976 | -0.1922 | -0.0734 | +1.0000 | -0.0630 |
| Abell 370 | -0.2297 | -0.0503 | +0.0450 | -0.0630 | +1.0000 |

Pearson values are computed on the dimensionless predicted
κ fields and on the predicted-minus-observed residual
fields.  A value near 0 indicates the residual structure
is not shared between clusters.


### Q6. Bias quantification

Bias is measured as the mean of (predicted - observed)
across the pipeline 64×64 grid.  A positive value means
Version A overestimates; a negative value means it
underestimates.

| Cluster | Mean (κ_pred - κ_obs) | Mean (|γ|_pred - |γ|_obs) | predicted κ mean | observed κ mean | predicted |γ| mean | observed |γ| mean |
|---|---|---|---|---|---|
| Abell 2744 | -0.05939 | +0.00155 | -0.01186 | +0.07228 | +0.07264 | +0.07598 |
| MACS J0416 | -0.09849 | -0.01258 | +0.00065 | +0.07423 | +0.07426 | +0.08350 |
| MACS J1149 | -0.03244 | -0.09306 | +0.01004 | +0.06653 | +0.02670 | +0.10619 |
| Abell S1063 | -0.05393 | -0.02200 | -0.01098 | +0.05142 | +0.05920 | +0.06726 |
| Abell 370 | -0.06748 | -0.01547 | -0.00805 | +0.08096 | +0.08484 | +0.09705 |

**Aggregate κ bias:** mean of per-cluster biases = -0.06235, std = 0.02148.
**Aggregate |γ| bias:** mean of per-cluster biases = -0.02831, std = 0.03328.
**Sign consistency:** κ bias sign: 0 positive, 5 negative, 0 zero. |γ| bias sign: 1 positive, 4 negative.

### Q7. Numerical stability

Stability is verified per cluster using three indicators:
(a) conservation error, (b) runtime, (c) absence of
exceptions during the run.

| Cluster | Runtime (s) | Conservation max | Status |
|---|---|---|---|
| Abell 2744 | 0.22 | 2.220e-16 | OK (machine epsilon) |
| MACS J0416 | 0.22 | 2.220e-16 | OK (machine epsilon) |
| MACS J1149 | 0.22 | 2.220e-16 | OK (machine epsilon) |
| Abell S1063 | 0.22 | 2.220e-16 | OK (machine epsilon) |
| Abell 370 | 0.22 | 2.220e-16 | OK (machine epsilon) |

The machine-precision conservation limit (`2.220446049250313e-16`) is met for every run on every
cluster, at both minimum production and high-accuracy
configurations.


## Outcome determination

**Outcome B** - agreement is limited.  Median Pearson κ = +0.0895, median SSIM κ = -0.0106.  Aggregate κ bias = -0.06235 ± 0.02148 with 5/5 clusters sharing the same sign (systematic). The frozen Version A laboratory is not modified.  Recurring residual patterns are documented in `residual_statistics.csv`, `cluster_comparison.csv`, and the per-cluster residual maps.  Future improvements shall be investigated through Version B physics, not by altering the frozen Version 1 laboratory.

## Numerical stability report

| Cluster | Runtime (s) | max conservation | memory peak (qualitative) |
|---|---|---|---|
| Abell 2744 | 0.22 | 2.220e-16 | stable (no exceptions) |
| MACS J0416 | 0.22 | 2.220e-16 | stable (no exceptions) |
| MACS J1149 | 0.22 | 2.220e-16 | stable (no exceptions) |
| Abell S1063 | 0.22 | 2.220e-16 | stable (no exceptions) |
| Abell 370 | 0.22 | 2.220e-16 | stable (no exceptions) |

## Notes

- All five benchmark FITS inputs were read exactly as
  archived in `PBUF_benchmark/`.  No network access.
- No fitting, no parameter search, no cosmological
  rescaling.  The frozen laboratory outputs are dimensionless
  pipeline-units fields compared like-with-like against the
  resampled observation.
- Cross-cluster consistency is documented in
  `cluster_comparison.csv` and the per-cluster morphology
  records are saved in `cluster_summary.json`.

## Top-level artefacts

- `runs/weak_lensing_science001/report.md` (this file)
- `runs/weak_lensing_science001/science_summary.csv`
- `runs/weak_lensing_science001/cluster_comparison.csv`
- `runs/weak_lensing_science001/observable_statistics.csv`
- `runs/weak_lensing_science001/residual_statistics.csv`
- `runs/weak_lensing_science001/run.json`
- `runs/weak_lensing_science001/validation.json`
- `runs/weak_lensing_science001/plots/kappa_comparison.png`
- `runs/weak_lensing_science001/plots/gamma_comparison.png`
- `runs/weak_lensing_science001/plots/residual_maps/`
- `runs/weak_lensing_science001/plots/radial_profiles.png`
- `runs/weak_lensing_science001/plots/tangential_shear.png`
- `runs/weak_lensing_science001/plots/deflection_vectors.png`
- `runs/weak_lensing_science001/plots/cluster_rankings.png`
- `runs/weak_lensing_science001/plots/science_dashboard.png`

## Per-cluster directories

For each cluster, two sub-directories are produced:
`{cluster_id}/minimum_production/` and
`{cluster_id}/high_accuracy/`.  Each contains:
`constitutive/`, `predicted/`, `observed/`,
`residual/`, `trajectories/`, several comparison PNGs,
and a `cluster_summary.json`.

**Total execution time:** 96.7 s.
