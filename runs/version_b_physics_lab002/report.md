# PBUF VERSION-B PHYSICS-LAB-002

**Response family decomposition of C10 inside the frozen
Version 1 weak-lensing laboratory (LAB-FREEZE-001).**

C10 (Combined Local Response) from PHYSICS-LAB-001 is decomposed
into the two physical mechanisms it contains: Neighbour Coherence
(A) and Elastic Memory (B). Each mechanism and their combination
are evaluated inside the frozen laboratory without modifying any
other component.

## Status

- Frozen hash verification: **PASS**
- Variants tested: **4** (Control + 3 mechanism subsets of {A, B})
- Clusters: **5**

## Frozen laboratory

The Version 1 laboratory is used as the measurement instrument
without modification.  All frozen source files are verified by
SHA-256 against LAB-FREEZE-001.

| Component | Frozen specification |
|---|---|
| Constitutive | `C(X) = 0.18 * rho(X) / rho_max` (Version A) |
| Transport | neighbour-to-neighbour, direct addition, |
| | per-step unit-speed renormalisation |
| Response direction | 90 deg transverse (R_90 of grad C) |
| Source plane | Launch B (Cartesian 2D) |
| Observable | Jacobian (ray-bundle linear fit per bin) |
| Matter input | `rho = max(kappa_obs, 0) / max(max(kappa_obs, 0))` |

## Decomposed mechanisms in C10

C10 in PHYSICS-LAB-001 was implemented as

    r = coherence_factor * ((1-w) * R_90(g) + w * R_90(g_prev))

with `w = 0.5` and

    coherence_factor = 0.5 * (1 + mean_cos(theta_self, theta_8nn))

Therefore C10 contains exactly two physical mechanisms:

| Code | Mechanism | Implementation |
|---|---|---|
| A | Neighbour Coherence | Multiplicative factor over 8 nearest neighbours |
| B | Elastic Memory      | One-step persistence mix (w = 0.5) |

## Variants

Every subset of {A, B} is realised as a local response law.  No
new physics is introduced.

| Code | Mechanism A | Mechanism B | Description |
|---|---|---|---|
| CONTROL | no | no | Response = |grad C|; frozen Version A control. |
| C10-A | yes | no | Magnitude scaled by (1 + mean_cos)/2 over 8 neighbours. |
| C10-B | no | yes | r_new = (1-w)*R(g) + w*R(g_prev); w = 0.5. |
| C10-C | yes | yes | Coherence factor times elastic memory mix (original C10). |

All fixed parameters are documented in the candidate source
code.  No parameter is fitted.

## Production configuration

| Parameter | Value |
|---|---|
| Photons | 20,000 |
| Constitutive grid | 256^2 |
| Step size | Delta s / 2 = 0.0300 |
| Number of steps | 160 |
| Source plane | Cartesian 2D (Launch B) |
| Observable | Jacobian |

## Per-variant, per-cluster metrics

Computed metrics for every (variant, cluster) pair:

| Variant | Cluster | RMS k | RMS g | Pearson k | Pearson g | SSIM k | SSIM g | k bias | g bias | conservation | runtime (s) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| CONTROL Gradient (control) | Abell 2744 | 1.5580e-01 | 6.0822e-02 | +0.0140 | +0.0895 | -0.0106 | +0.1008 | -5.9392e-02 | +1.5525e-03 | 2.220e-16 | 0.233 |
| CONTROL Gradient (control) | MACS J0416 | 1.8596e-01 | 6.5276e-02 | +0.0190 | +0.1175 | +0.0004 | +0.1267 | -9.8489e-02 | -1.2579e-02 | 2.220e-16 | 0.230 |
| CONTROL Gradient (control) | MACS J1149 | 8.5592e-02 | 1.1556e-01 | +0.2390 | +0.0075 | +0.1081 | +0.0120 | -3.2439e-02 | -9.3065e-02 | 2.220e-16 | 0.231 |
| CONTROL Gradient (control) | Abell S1063 | 1.3540e-01 | 6.9710e-02 | +0.0895 | +0.0836 | -0.0383 | +0.0929 | -5.3933e-02 | -2.1999e-02 | 2.220e-16 | 0.231 |
| CONTROL Gradient (control) | Abell 370 | 1.7915e-01 | 9.1733e-02 | +0.1954 | -0.0171 | -0.0398 | +0.0122 | -6.7484e-02 | -1.5468e-02 | 2.220e-16 | 0.226 |
| C10-A Neighbour Coherence only | Abell 2744 | 1.4344e-01 | 5.7249e-02 | +0.0256 | +0.0860 | -0.0129 | +0.0973 | -5.8019e-02 | -5.3675e-03 | 2.220e-16 | 0.215 |
| C10-A Neighbour Coherence only | MACS J0416 | 1.7308e-01 | 6.3177e-02 | +0.0235 | +0.1136 | +0.0010 | +0.1239 | -9.7374e-02 | -1.9917e-02 | 2.220e-16 | 0.216 |
| C10-A Neighbour Coherence only | MACS J1149 | 8.3277e-02 | 1.1704e-01 | +0.2391 | +0.0082 | +0.0960 | +0.0112 | -3.3313e-02 | -9.5478e-02 | 2.220e-16 | 0.217 |
| C10-A Neighbour Coherence only | Abell S1063 | 1.2608e-01 | 6.8610e-02 | +0.0935 | +0.0810 | -0.0306 | +0.0869 | -5.0951e-02 | -2.6319e-02 | 2.220e-16 | 0.223 |
| C10-A Neighbour Coherence only | Abell 370 | 1.6048e-01 | 8.9565e-02 | +0.2063 | +0.0129 | -0.0292 | +0.0406 | -6.4601e-02 | -2.6507e-02 | 2.220e-16 | 0.215 |
| C10-B Elastic Memory only | Abell 2744 | 1.5551e-01 | 6.0280e-02 | +0.0088 | +0.0824 | -0.0093 | +0.0953 | -5.9732e-02 | -9.6621e-04 | 2.220e-16 | 0.234 |
| C10-B Elastic Memory only | MACS J0416 | 1.8648e-01 | 6.5138e-02 | +0.0105 | +0.1403 | +0.0001 | +0.1478 | -9.8932e-02 | -1.3677e-02 | 2.220e-16 | 0.220 |
| C10-B Elastic Memory only | MACS J1149 | 8.5737e-02 | 1.1673e-01 | +0.2328 | -0.0389 | +0.1053 | -0.0002 | -3.2457e-02 | -9.3967e-02 | 2.220e-16 | 0.228 |
| C10-B Elastic Memory only | Abell S1063 | 1.3382e-01 | 6.8874e-02 | +0.0946 | +0.0862 | -0.0403 | +0.0960 | -5.3841e-02 | -2.3095e-02 | 2.220e-16 | 0.216 |
| C10-B Elastic Memory only | Abell 370 | 1.7942e-01 | 9.0267e-02 | +0.1866 | +0.0012 | -0.0403 | +0.0298 | -6.7925e-02 | -1.9230e-02 | 2.220e-16 | 0.232 |
| C10-C Combined (original C10) | Abell 2744 | 1.3990e-01 | 5.6182e-02 | +0.0210 | +0.0824 | -0.0117 | +0.0940 | -5.8042e-02 | -9.3295e-03 | 2.220e-16 | 0.217 |
| C10-C Combined (original C10) | MACS J0416 | 1.7082e-01 | 6.3011e-02 | +0.0117 | +0.1395 | +0.0009 | +0.1469 | -9.6652e-02 | -2.2600e-02 | 2.220e-16 | 0.229 |
| C10-C Combined (original C10) | MACS J1149 | 8.2467e-02 | 1.1850e-01 | +0.2355 | -0.0208 | +0.0958 | +0.0047 | -3.3028e-02 | -9.7154e-02 | 2.220e-16 | 0.216 |
| C10-C Combined (original C10) | Abell S1063 | 1.2270e-01 | 6.7860e-02 | +0.1034 | +0.0865 | -0.0338 | +0.0910 | -5.0881e-02 | -2.8686e-02 | 2.220e-16 | 0.227 |
| C10-C Combined (original C10) | Abell 370 | 1.5946e-01 | 8.8799e-02 | +0.1950 | +0.0357 | -0.0430 | +0.0593 | -6.7297e-02 | -3.1552e-02 | 2.220e-16 | 0.223 |

## Cross-cluster evaluation

For every variant the following medians/means are taken across
the five benchmark clusters.

| Variant | Median Pearson k | Median Pearson g | Median SSIM k | Mean k Bias | Mean g Bias | Mean Pearson k | Conservation max | Runtime (s) |
|---|---|---|---|---|---|---|---|---|
| CONTROL Gradient (control) | +0.0895 | +0.0836 | -0.0106 | -6.2347e-02 | -2.8312e-02 | +0.1114 | 2.220e-16 | 0.231 |
| C10-A Neighbour Coherence only | +0.0935 | +0.0810 | -0.0129 | -6.0852e-02 | -3.4718e-02 | +0.1176 | 2.220e-16 | 0.216 |
| C10-B Elastic Memory only | +0.0946 | +0.0824 | -0.0093 | -6.2577e-02 | -3.0187e-02 | +0.1067 | 2.220e-16 | 0.228 |
| C10-C Combined (original C10) | +0.1034 | +0.0824 | -0.0117 | -6.1180e-02 | -3.7864e-02 | +0.1133 | 2.220e-16 | 0.223 |

## Family ranking by individual contribution

Variants ranked by `median Pearson kappa` minus the control
(individual contribution).  Combined contribution and
lost-when-removed are reported alongside.

| Rank | Code | Name | Family | Alone delta Pearson k | Combined delta | Lost when removed | Redundant |
|---|---|---|---|---|---|---|---|
| 1 | C10-C | Combined (original C10) | combined response | +0.01389 | +0.01389 | +0.00000 | YES |
| 2 | C10-B | Elastic Memory only | elastic memory | +0.00507 | +0.01389 | +0.00881 | NO |
| 3 | C10-A | Neighbour Coherence only | neighbour coherence | +0.00396 | +0.01389 | +0.00993 | NO |

## Contribution analysis

Decomposition of C10 improvement (vs frozen Version A control):

| Source | Delta Pearson k | Delta SSIM k | Delta kappa bias |
|---|---|---|---|
| A alone (Coherence) | +0.00396 | -0.00235 | +0.00150 |
| B alone (Memory)    | +0.00507 | +0.00129 | -0.00023 |
| Sum A + B           | +0.00903 | -0.00106 | +0.00127 |
| Combined (C10-C)    | +0.01389 | -0.00115 | +0.00117 |
| Interaction         | +0.00485 | -0.00009 | -0.00010 |

## Required questions

### Q1. Largest individual contribution

Control median Pearson kappa = +0.08951, median SSIM kappa = -0.01058, mean kappa bias = -0.06235.

Individual mechanism contributions (excluding the combined):

| Mechanism | Median Pearson k delta | Median SSIM k delta | Mean k bias delta |
|---|---|---|---|
| C10-B Elastic Memory only | +0.00507 | +0.00129 | -0.00023 |
| C10-A Neighbour Coherence only | +0.00396 | -0.00235 | +0.00150 |

For reference, the combined C10-C delta = +0.01389 (Pearson k), -0.00115 (SSIM k), +0.00117 (k bias).

**Largest individual contribution:** `C10-B` (Elastic Memory only) with median Pearson kappa delta = +0.00507.

### Q2. Single-mechanism explanation of C10

| Mechanism | Alone delta Pearson k | Combined delta Pearson k | Share |
|---|---|---|---|
| Coherence (A) | +0.00396 | +0.01389 | 28.5% |
| Memory (B)    | +0.00507 | +0.01389 | 36.5% |

Best single-mechanism explanation share = 36.5% (mechanism `Memory (B)`).
C10 improvement explained entirely by one mechanism: **NO**.

Additional metrics:

| Mechanism | SSIM k delta alone | SSIM k delta combined | Bias delta alone | Bias delta combined |
|---|---|---|---|---|
| Coherence (A) | -0.00235 | -0.00115 | +0.00150 | +0.00117 |
| Memory (B)    | +0.00129 | -0.00115 | -0.00023 | +0.00117 |

### Q3. Mutual reinforcement

Does the combined improvement exceed the best single-mechanism
improvement?

- Coherence alone delta: +0.00396
- Memory alone delta:    +0.00507
- Combined delta:        +0.01389
- Best single delta:     +0.00507

Reinforcement present (combined > best alone): **YES**.

### Q4. Redundancy check

Combined delta Pearson kappa = +0.01389.
Redundancy threshold = 0.01% of |combined delta|.

| Mechanism | Combined - alone | Redundant? |
|---|---|---|
| Coherence (A) | +0.00993 | NO |
| Memory (B)    | +0.00881 | NO |

### Q5. Universally sign-consistent improvement

| Variant | Median delta Pearson k | # clusters +ve | # clusters -ve | Sign consistent |
|---|---|---|---|---|
| C10-A Neighbour Coherence only | +0.00449 | 5 | 0 | YES |
| C10-C Combined (original C10) | -0.00031 | 2 | 3 | NO |
| C10-B Elastic Memory only | -0.00623 | 1 | 4 | NO |

Number of variants with sign-consistent improvement on all 5 clusters: 1/3.

Cross-check with VERSION-B PHYSICS-LAB-001 (Candidate 2 = C10-A):
  C2 / C10-A sign-consistent: **YES** (5 +ve, 0 -ve, median delta = +0.00449).

### Q6. Nature of the interaction

| Metric | A alone | B alone | A+B combined | Interaction |
|---|---|---|---|---|
| Pearson kappa delta | +0.00396 | +0.00507 | +0.01389 | +0.00485 |
| SSIM kappa delta    | -0.00235 | +0.00129 | -0.00115 | -0.00009 |
| Kappa bias delta    | +0.00150 | -0.00023 | +0.00117 | -0.00010 |

Interaction magnitude relative to combined improvement: 35.0%.

Verdict: **NONLINEAR / MIXED**.


## Outcome determination

**Outcome C** - the apparent improvement cannot be attributed to any individual mechanism. Best single-mechanism share = 36.5%, relative interaction = 35.0%. The combined behaviour is emergent and requires further investigation.

## Numerical stability report

| Variant | Median runtime (s) | Max conservation |
|---|---|---|
| CONTROL Gradient (control) | 0.231 | 2.220e-16 |
| C10-A Neighbour Coherence only | 0.216 | 2.220e-16 |
| C10-B Elastic Memory only | 0.228 | 2.220e-16 |
| C10-C Combined (original C10) | 0.223 | 2.220e-16 |

## Top-level artefacts

- `runs/version_b_physics_lab002/report.md` (this file)
- `runs/version_b_physics_lab002/component_contributions.csv`
- `runs/version_b_physics_lab002/interaction_matrix.csv`
- `runs/version_b_physics_lab002/leave_one_out.csv`
- `runs/version_b_physics_lab002/family_ranking.csv`
- `runs/version_b_physics_lab002/run.json`
- `runs/version_b_physics_lab002/validation.json`
- `runs/version_b_physics_lab002/plots/component_importance.png`
- `runs/version_b_physics_lab002/plots/interaction_heatmap.png`
- `runs/version_b_physics_lab002/plots/leave_one_out.png`
- `runs/version_b_physics_lab002/plots/synergy_matrix.png`
- `runs/version_b_physics_lab002/plots/family_contributions.png`

**Total execution time:** 5.4 s.
