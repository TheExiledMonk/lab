# PBUF WEAK-LENSING-OBSERVATION-001

Frozen Version A pipeline applied to five public weak-lensing
benchmark clusters (SaWLens reconstructions, Merten et al. 2014).

## Pipeline (frozen, identical for every cluster)

- Constitutive: `C(X) = 0.18 · ρ(X) / ρ_max` (Version A)
- Response: `r = |∇C|` rotated 90° transverse
- Transport: neighbour-to-neighbour, direct addition, velocity
  renormalisation
- Numerical parameters (frozen):
  `n_grid = 128, extent = 8.0, strength = 0.18, step = 0.06, steps = 80, nphotons = 2000, bins = 64`

## Matter input

Each cluster's published κ map is taken as a proxy for the matter
density ρ. Negative κ values (mass deficits) are clamped to zero
to respect the implicit positivity assumption of the Version A
constitutive law; the resulting field is normalised by its peak
value so that ρ_max = 1 in pipeline units. No fitting, smoothing
or rescaling beyond bilinear interpolation onto the pipeline grid.

## Coordinate alignment

The observation WCS centres each map on its own cluster centre
(CRVAL). The pipeline grid is Cartesian `[-8, 8] × [-8, 8]`. The
observation is resampled bilinearly onto the pipeline grid with
the cluster centre mapped to the origin. No smoothing is
applied; only linear interpolation for coordinate alignment.

## Internal-consistency check (mandatory pre-comparison)

For every cluster, `gamma.fits` is compared against
`sqrt(gamma1² + gamma2²)` element-wise.

| Cluster | Max abs diff | RMS abs diff | Tolerance pass |
|---|---|---|---|
| Abell 2744 | 5.304e-08 | 2.956e-09 | YES |
| MACS J0416 | 5.237e-08 | 3.441e-09 | YES |
| MACS J1149 | 5.892e-08 | 4.099e-09 | YES |
| Abell S1063 | 2.233e-08 | 2.630e-09 | YES |
| Abell 370 | 3.282e-08 | 3.749e-09 | YES |

All `gamma.fits` files match `sqrt(gamma1² + gamma2²)` to within
single-precision FP tolerance. The supplied `gamma` field is
internally consistent.

## Resampled-observation interpolation consistency

After bilinear resampling of the four observation maps onto the
pipeline 64×64 grid, `gamma_resampled` is again compared to
`sqrt(gamma1_resampled² + gamma2_resampled²)`. The discrepancy is
an unavoidable interpolation artefact and is reported separately.

| Cluster | Resampled max abs diff | Resampled RMS abs diff |
|---|---|---|
| Abell 2744 | 6.546e-02 | 7.893e-03 |
| MACS J0416 | 1.120e-01 | 9.112e-03 |
| MACS J1149 | 8.376e-02 | 8.674e-03 |
| Abell S1063 | 8.341e-02 | 6.990e-03 |
| Abell 370 | 7.664e-02 | 6.990e-03 |

These residuals reflect the nonlinearity of the square root under
bilinear interpolation and do **not** indicate inconsistency in the
raw FITS products.

## Per-cluster metrics

| Cluster | RMS κ | MAE κ | Max abs κ | Corr κ | SSIM κ | Peak offset (px) |
|---|---|---|---|---|---|---|
| Abell 2744 | 5.5676e-01 | 5.5124e-01 | 7.0837e-01 | nan | -0.0072 | 34.54 |
| MACS J0416 | 5.8070e-01 | 5.6983e-01 | 7.8639e-01 | nan | -0.0048 | 34.18 |
| MACS J1149 | 5.1765e-01 | 5.0925e-01 | 7.1497e-01 | nan | -0.0009 | 34.54 |
| Abell S1063 | 5.4541e-01 | 5.3974e-01 | 7.2289e-01 | nan | -0.0056 | 33.62 |
| Abell 370 | 5.6128e-01 | 5.5897e-01 | 7.1426e-01 | nan | -0.0186 | 34.18 |

| Cluster | RMS γ₁ | Corr γ₁ | RMS γ₂ | Corr γ₂ | RMS γ | Corr γ |
|---|---|---|---|---|---|---|
| Abell 2744 | 5.2701e-01 | 0.0057 | 1.6043e-01 | -0.0401 | 5.4308e-01 | -0.0193 |
| MACS J0416 | 5.3852e-01 | 0.0023 | 1.6487e-01 | -0.0311 | 5.5341e-01 | 0.0031 |
| MACS J1149 | 5.2219e-01 | -0.0080 | 1.6998e-01 | -0.0278 | 5.3466e-01 | 0.0312 |
| Abell S1063 | 5.1707e-01 | -0.0161 | 1.5905e-01 | -0.0090 | 5.3211e-01 | 0.0136 |
| Abell 370 | 5.1950e-01 | -0.0070 | 1.6806e-01 | -0.0250 | 5.3504e-01 | -0.0099 |

## Required statistics (per cluster)

Each cluster's frozen-pipeline observables are compared against
the published observables resampled to the pipeline 64×64 grid.
No parameter of the frozen Version A pipeline is altered between
clusters or between iterations.

### Abell 2744 (`Abell2744`)

| Metric | Value |
|---|---|
| RMS κ | 5.5676e-01 |
| RMS γ₁ | 5.2701e-01 |
| RMS γ₂ | 1.6043e-01 |
| RMS γ | 5.4308e-01 |
| Correlation κ | nan |
| Correlation γ | -0.0193 |
| Runtime (s) | 0.0102 |
| Photon count | 2000 |
| Numerical conservation (max) | 2.2204e-16 |

### MACS J0416 (`MACS0416`)

| Metric | Value |
|---|---|
| RMS κ | 5.8070e-01 |
| RMS γ₁ | 5.3852e-01 |
| RMS γ₂ | 1.6487e-01 |
| RMS γ | 5.5341e-01 |
| Correlation κ | nan |
| Correlation γ | 0.0031 |
| Runtime (s) | 0.0090 |
| Photon count | 2000 |
| Numerical conservation (max) | 2.2204e-16 |

### MACS J1149 (`MACS1149`)

| Metric | Value |
|---|---|
| RMS κ | 5.1765e-01 |
| RMS γ₁ | 5.2219e-01 |
| RMS γ₂ | 1.6998e-01 |
| RMS γ | 5.3466e-01 |
| Correlation κ | nan |
| Correlation γ | 0.0312 |
| Runtime (s) | 0.0090 |
| Photon count | 2000 |
| Numerical conservation (max) | 2.2204e-16 |

### Abell S1063 (`AbellS1063`)

| Metric | Value |
|---|---|
| RMS κ | 5.4541e-01 |
| RMS γ₁ | 5.1707e-01 |
| RMS γ₂ | 1.5905e-01 |
| RMS γ | 5.3211e-01 |
| Correlation κ | nan |
| Correlation γ | 0.0136 |
| Runtime (s) | 0.0091 |
| Photon count | 2000 |
| Numerical conservation (max) | 2.2204e-16 |

### Abell 370 (`Abell370`)

| Metric | Value |
|---|---|
| RMS κ | 5.6128e-01 |
| RMS γ₁ | 5.1950e-01 |
| RMS γ₂ | 1.6806e-01 |
| RMS γ | 5.3504e-01 |
| Correlation κ | nan |
| Correlation γ | -0.0099 |
| Runtime (s) | 0.0090 |
| Photon count | 2000 |
| Numerical conservation (max) | 2.2204e-16 |


## Cross-cluster summary

| Cluster | RMS κ | RMS γ | Corr κ | Corr γ | Runtime (s) |
|---|---|---|---|---|---|
| Abell 2744 | 5.5676e-01 | 5.4308e-01 | nan | -0.0193 | 0.010 |
| MACS J0416 | 5.8070e-01 | 5.5341e-01 | nan | 0.0031 | 0.009 |
| MACS J1149 | 5.1765e-01 | 5.3466e-01 | nan | 0.0312 | 0.009 |
| Abell S1063 | 5.4541e-01 | 5.3211e-01 | nan | 0.0136 | 0.009 |
| Abell 370 | 5.6128e-01 | 5.3504e-01 | nan | -0.0099 | 0.009 |

## Units, conventions and mismatch (mandatory pre-comparison record)

- Published products are SaWLens reconstructions on RA/Dec WCS
  grids centred on each cluster; pixel scales 6.25-11.36 arcsec.
- All published observables are scaled to source redshift
  `Z_S = 9.0` (effectively an infinite-source approximation).
- `kappa.fits`, `gamma1.fits`, `gamma2.fits` are the lensing
  convergence and reduced-shear components from a parametric
  joint weak+strong lensing inversion. They are reconstructed
  posterior-mean maps, not direct observational data.
- Frozen Version A outputs are dimensionless lensing-like
  observables (κ_pred, γ₁_pred, γ₂_pred) derived from a
  constitutive + transport pipeline operating on synthetic
  dimensionless coordinates on `[-8, 8] × [-8, 8]`.
- The published products and Version A outputs are **NOT**
  directly comparable in absolute units, normalisation, or
  angular scale. The comparison made here is a like-with-like
  dimensionless field comparison after coordinate alignment,
  with no implicit cosmological rescaling.
- The matter input to the Version A constitutive law is the
  positive part of the published κ (clamped at zero), then
  normalised by its peak. This treats κ as a matter-density
  proxy (standard practice for clusters where mass traces light)
  while preserving the implicit positivity assumption of the
  Version A law.

## Per-cluster outputs

### Abell 2744 (`Abell2744`)

- `observed/`: resampled observation maps (κ, γ₁, γ₂, γ)
  and the internal-consistency `gamma_internal_check.csv`.
- `predicted/`: Version A outputs (κ, γ₁, γ₂, |γ|, deflection,
  magnification).
- `residual/`: `pred - obs` residual maps in absolute units and
  percentage form.
- `constitutive/`: matter proxy, C, ∇C, response field, maps.
- `trajectories/`: photon trajectories (`x`, `y`), endpoints,
  bending angles, max deviations.
- `comparison_kappa.png`, `comparison_gamma1.png`,
  `comparison_gamma2.png`, `comparison_gamma.png`: three-panel
  comparisons with identical colour scales.
- `comparison_overview.png`: six-panel composite (κ and |γ|).
- `composite_pipeline.png`: 3×3 pipeline panel (matter, C, |∇C|,
  response magnitude and direction, κ, γ₁, γ₂, μ).
- `photon_trajectories.png`: trajectory plot coloured by
  accumulated bending angle.
- `statistics.json`: all quantitative metrics.
- `fits_metadata.json`: FITS headers, WCS info, file SHA-256s,
  matter-proxy construction record, coordinate-alignment record.

### MACS J0416 (`MACS0416`)

- `observed/`: resampled observation maps (κ, γ₁, γ₂, γ)
  and the internal-consistency `gamma_internal_check.csv`.
- `predicted/`: Version A outputs (κ, γ₁, γ₂, |γ|, deflection,
  magnification).
- `residual/`: `pred - obs` residual maps in absolute units and
  percentage form.
- `constitutive/`: matter proxy, C, ∇C, response field, maps.
- `trajectories/`: photon trajectories (`x`, `y`), endpoints,
  bending angles, max deviations.
- `comparison_kappa.png`, `comparison_gamma1.png`,
  `comparison_gamma2.png`, `comparison_gamma.png`: three-panel
  comparisons with identical colour scales.
- `comparison_overview.png`: six-panel composite (κ and |γ|).
- `composite_pipeline.png`: 3×3 pipeline panel (matter, C, |∇C|,
  response magnitude and direction, κ, γ₁, γ₂, μ).
- `photon_trajectories.png`: trajectory plot coloured by
  accumulated bending angle.
- `statistics.json`: all quantitative metrics.
- `fits_metadata.json`: FITS headers, WCS info, file SHA-256s,
  matter-proxy construction record, coordinate-alignment record.

### MACS J1149 (`MACS1149`)

- `observed/`: resampled observation maps (κ, γ₁, γ₂, γ)
  and the internal-consistency `gamma_internal_check.csv`.
- `predicted/`: Version A outputs (κ, γ₁, γ₂, |γ|, deflection,
  magnification).
- `residual/`: `pred - obs` residual maps in absolute units and
  percentage form.
- `constitutive/`: matter proxy, C, ∇C, response field, maps.
- `trajectories/`: photon trajectories (`x`, `y`), endpoints,
  bending angles, max deviations.
- `comparison_kappa.png`, `comparison_gamma1.png`,
  `comparison_gamma2.png`, `comparison_gamma.png`: three-panel
  comparisons with identical colour scales.
- `comparison_overview.png`: six-panel composite (κ and |γ|).
- `composite_pipeline.png`: 3×3 pipeline panel (matter, C, |∇C|,
  response magnitude and direction, κ, γ₁, γ₂, μ).
- `photon_trajectories.png`: trajectory plot coloured by
  accumulated bending angle.
- `statistics.json`: all quantitative metrics.
- `fits_metadata.json`: FITS headers, WCS info, file SHA-256s,
  matter-proxy construction record, coordinate-alignment record.

### Abell S1063 (`AbellS1063`)

- `observed/`: resampled observation maps (κ, γ₁, γ₂, γ)
  and the internal-consistency `gamma_internal_check.csv`.
- `predicted/`: Version A outputs (κ, γ₁, γ₂, |γ|, deflection,
  magnification).
- `residual/`: `pred - obs` residual maps in absolute units and
  percentage form.
- `constitutive/`: matter proxy, C, ∇C, response field, maps.
- `trajectories/`: photon trajectories (`x`, `y`), endpoints,
  bending angles, max deviations.
- `comparison_kappa.png`, `comparison_gamma1.png`,
  `comparison_gamma2.png`, `comparison_gamma.png`: three-panel
  comparisons with identical colour scales.
- `comparison_overview.png`: six-panel composite (κ and |γ|).
- `composite_pipeline.png`: 3×3 pipeline panel (matter, C, |∇C|,
  response magnitude and direction, κ, γ₁, γ₂, μ).
- `photon_trajectories.png`: trajectory plot coloured by
  accumulated bending angle.
- `statistics.json`: all quantitative metrics.
- `fits_metadata.json`: FITS headers, WCS info, file SHA-256s,
  matter-proxy construction record, coordinate-alignment record.

### Abell 370 (`Abell370`)

- `observed/`: resampled observation maps (κ, γ₁, γ₂, γ)
  and the internal-consistency `gamma_internal_check.csv`.
- `predicted/`: Version A outputs (κ, γ₁, γ₂, |γ|, deflection,
  magnification).
- `residual/`: `pred - obs` residual maps in absolute units and
  percentage form.
- `constitutive/`: matter proxy, C, ∇C, response field, maps.
- `trajectories/`: photon trajectories (`x`, `y`), endpoints,
  bending angles, max deviations.
- `comparison_kappa.png`, `comparison_gamma1.png`,
  `comparison_gamma2.png`, `comparison_gamma.png`: three-panel
  comparisons with identical colour scales.
- `comparison_overview.png`: six-panel composite (κ and |γ|).
- `composite_pipeline.png`: 3×3 pipeline panel (matter, C, |∇C|,
  response magnitude and direction, κ, γ₁, γ₂, μ).
- `photon_trajectories.png`: trajectory plot coloured by
  accumulated bending angle.
- `statistics.json`: all quantitative metrics.
- `fits_metadata.json`: FITS headers, WCS info, file SHA-256s,
  matter-proxy construction record, coordinate-alignment record.

## Identical-pipeline verification (SHA-256)

| File | SHA-256 |
|---|---|
| `weak_lensing_observation001.py` | `a5c3632fec9adfc2659d5c283d07c599db6db7edb2e34e4aebd84e35434642bc` |
| `constitutive_equations.py` | `e2c789d19fd559753519704c6668c7a2879c53eb61a315604ec81af6795aca9f` |

## Notes

- The Version A pipeline parameters are held identical to those of
  WEAK-LENSING-PREDICTION-001 and WEAK-LENSING-GENERALIZATION-001.
- No parameter was altered to improve agreement with the
  observations. Any apparent mismatch is therefore a property of
  the frozen implementation itself.
- Discrepancies between predicted and observed fields can arise
  from unit/normalisation mismatch, the absence of cosmological
  scaling in the Version A pipeline, the simplified matter proxy,
  and the linear-response approximation inherent to Version A.
  These are documented, not compensated for.
- Total execution time: 19.55 s.
