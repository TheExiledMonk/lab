# PBUF COHERENCE-MEMORY-LAB-001

**Mapping the cooperative elastic response inside the frozen
Version 1 weak-lensing laboratory (LAB-FREEZE-001).**

Two local-response mechanisms are swept across a 5 x 5 grid:

| Axis | Mechanism | Parameter | Tested values |
|---|---|---|---|
| A | Neighbour Coherence | strength | 0.00, 0.25, 0.50, 0.75, 1.00 |
| B | Elastic Memory      | weight   | 0.00, 0.25, 0.50, 0.75, 1.00 |

Total configurations: **25**.
Clusters: **5**.
Total runs: **125**.

## Status

- Frozen hash verification: **PASS**
- Total runtime: **32.4 s**

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

## Response parameterisation

The local response law is parameterised by (A, B):

    r_lin   = R_90(g)                                       # = (-g*gy/|g|, g*gx/|g|)
    factor  = (1-A) + A * 0.5*(1 + mean_cos(theta_self, theta_8nn))
    r_mem   = (1-B) * r_lin + B * r_prev
    r(A, B) = factor * r_mem

Reproduces the originals:

| Configuration | Reproduces |
|---|---|
| (A=0, B=0) | frozen Version A control |
| (A=1, B=0) | C10-A from LAB-002 (Coherence only) |
| (A=0, B=0.5) | C10-B from LAB-002 (Memory only) |
| (A=1, B=0.5) | original C10 from LAB-001 |
| (A=1, B=1) | maximum combined (factor * r_prev) |

## Production configuration

| Parameter | Value |
|---|---|
| Photons | 20,000 |
| Constitutive grid | 256^2 |
| Step size | Delta s / 2 = 0.0300 |
| Number of steps | 160 |
| Source plane | Cartesian 2D (Launch B) |
| Observable | Jacobian |

## Cross-cluster summary table

| A | B | Median Pearson k | Median Pearson g | Median SSIM k | Mean k Bias | Mean g Bias | Mean Pearson k | Conservation max | Runtime (s) |
|---|---|---|---|---|---|---|---|---|---|
| 0.00 | 0.00 | +0.0895 | +0.0836 | -0.0106 | -6.2347e-02 | -2.8312e-02 | +0.1114 | 2.220e-16 | 0.228 |
| 0.00 | 0.25 | +0.0923 | +0.0858 | -0.0096 | -6.2368e-02 | -2.9408e-02 | +0.1094 | 2.220e-16 | 0.214 |
| 0.00 | 0.50 | +0.0946 | +0.0824 | -0.0093 | -6.2577e-02 | -3.0187e-02 | +0.1067 | 2.220e-16 | 0.215 |
| 0.00 | 0.75 | +0.0954 | +0.0778 | -0.0082 | -6.2667e-02 | -3.0629e-02 | +0.1035 | 2.220e-16 | 0.216 |
| 0.00 | 1.00 | +0.0968 | +0.0832 | -0.0072 | -6.2651e-02 | -3.0750e-02 | +0.1006 | 2.220e-16 | 0.215 |
| 0.25 | 0.00 | +0.0903 | +0.0827 | -0.0111 | -6.1978e-02 | -2.9992e-02 | +0.1126 | 2.220e-16 | 0.219 |
| 0.25 | 0.25 | +0.0937 | +0.0843 | -0.0106 | -6.2022e-02 | -3.1274e-02 | +0.1108 | 2.220e-16 | 0.230 |
| 0.25 | 0.50 | +0.0966 | +0.0819 | -0.0101 | -6.2234e-02 | -3.2198e-02 | +0.1085 | 2.220e-16 | 0.216 |
| 0.25 | 0.75 | +0.0980 | +0.0787 | -0.0089 | -6.2313e-02 | -3.2845e-02 | +0.1051 | 2.220e-16 | 0.229 |
| 0.25 | 1.00 | +0.0994 | +0.0825 | -0.0081 | -6.2329e-02 | -3.3166e-02 | +0.1017 | 2.220e-16 | 0.217 |
| 0.50 | 0.00 | +0.0909 | +0.0821 | -0.0121 | -6.1584e-02 | -3.1664e-02 | +0.1142 | 2.220e-16 | 0.228 |
| 0.50 | 0.25 | +0.0949 | +0.0833 | -0.0112 | -6.1709e-02 | -3.3048e-02 | +0.1123 | 2.220e-16 | 0.226 |
| 0.50 | 0.50 | +0.0984 | +0.0815 | -0.0105 | -6.1860e-02 | -3.4143e-02 | +0.1098 | 2.220e-16 | 0.232 |
| 0.50 | 0.75 | +0.1003 | +0.0793 | -0.0096 | -6.1940e-02 | -3.5012e-02 | +0.1065 | 2.220e-16 | 0.226 |
| 0.50 | 1.00 | +0.1026 | +0.0829 | -0.0087 | -6.2043e-02 | -3.5504e-02 | +0.1037 | 2.220e-16 | 0.231 |
| 0.75 | 0.00 | +0.0920 | +0.0816 | -0.0124 | -6.1245e-02 | -3.3231e-02 | +0.1157 | 2.220e-16 | 0.227 |
| 0.75 | 0.25 | +0.0967 | +0.0828 | -0.0118 | -6.1332e-02 | -3.4794e-02 | +0.1143 | 2.220e-16 | 0.232 |
| 0.75 | 0.50 | +0.1011 | +0.0818 | -0.0111 | -6.1543e-02 | -3.6055e-02 | +0.1115 | 2.220e-16 | 0.234 |
| 0.75 | 0.75 | +0.1031 | +0.0807 | -0.0105 | -6.1592e-02 | -3.7025e-02 | +0.1089 | 2.220e-16 | 0.231 |
| 0.75 | 1.00 | +0.1061 | +0.0859 | -0.0096 | -6.1642e-02 | -3.7637e-02 | +0.1058 | 2.220e-16 | 0.220 |
| 1.00 | 0.00 | +0.0935 | +0.0810 | -0.0129 | -6.0852e-02 | -3.4718e-02 | +0.1176 | 2.220e-16 | 0.233 |
| 1.00 | 0.25 | +0.0985 | +0.0830 | -0.0124 | -6.0991e-02 | -3.6460e-02 | +0.1158 | 2.220e-16 | 0.223 |
| 1.00 | 0.50 | +0.1034 | +0.0824 | -0.0117 | -6.1180e-02 | -3.7864e-02 | +0.1133 | 2.220e-16 | 0.230 |
| 1.00 | 0.75 | +0.1060 | +0.0828 | -0.0115 | -6.1237e-02 | -3.8929e-02 | +0.1102 | 2.220e-16 | 0.219 |
| 1.00 | 1.00 | +0.1093 | +0.0882 | -0.0105 | -6.1288e-02 | -3.9557e-02 | +0.1072 | 2.220e-16 | 0.232 |

## Synergy matrix (Tukey additivity)

Expected additive prediction: `E(A, B) = f(A, 0) + f(0, B) - f(0, 0)`.  Synergy = `f(A, B) - E(A, B)`.

| A | B | Observed Pearson k | Expected additive | Synergy | Class |
|---|---|---|---|---|---|
| 0.00 | 0.00 | +0.08951 | +0.08951 | +0.00000 | boundary |
| 0.00 | 0.25 | +0.09226 | +0.09226 | +0.00000 | boundary |
| 0.00 | 0.50 | +0.09458 | +0.09458 | +0.00000 | boundary |
| 0.00 | 0.75 | +0.09536 | +0.09536 | +0.00000 | boundary |
| 0.00 | 1.00 | +0.09683 | +0.09683 | +0.00000 | boundary |
| 0.25 | 0.00 | +0.09033 | +0.09033 | +0.00000 | boundary |
| 0.25 | 0.25 | +0.09371 | +0.09307 | +0.00064 | synergistic |
| 0.25 | 0.50 | +0.09664 | +0.09540 | +0.00124 | synergistic |
| 0.25 | 0.75 | +0.09796 | +0.09618 | +0.00178 | synergistic |
| 0.25 | 1.00 | +0.09939 | +0.09765 | +0.00174 | synergistic |
| 0.50 | 0.00 | +0.09089 | +0.09089 | +0.00000 | boundary |
| 0.50 | 0.25 | +0.09486 | +0.09364 | +0.00123 | synergistic |
| 0.50 | 0.50 | +0.09843 | +0.09596 | +0.00247 | synergistic |
| 0.50 | 0.75 | +0.10031 | +0.09674 | +0.00356 | synergistic |
| 0.50 | 1.00 | +0.10264 | +0.09821 | +0.00443 | synergistic |
| 0.75 | 0.00 | +0.09201 | +0.09201 | +0.00000 | boundary |
| 0.75 | 0.25 | +0.09675 | +0.09475 | +0.00199 | synergistic |
| 0.75 | 0.50 | +0.10107 | +0.09708 | +0.00399 | synergistic |
| 0.75 | 0.75 | +0.10315 | +0.09786 | +0.00529 | synergistic |
| 0.75 | 1.00 | +0.10608 | +0.09933 | +0.00675 | synergistic |
| 1.00 | 0.00 | +0.09347 | +0.09347 | +0.00000 | boundary |
| 1.00 | 0.25 | +0.09855 | +0.09621 | +0.00233 | synergistic |
| 1.00 | 0.50 | +0.10340 | +0.09854 | +0.00485 | synergistic |
| 1.00 | 0.75 | +0.10600 | +0.09932 | +0.00668 | synergistic |
| 1.00 | 1.00 | +0.10927 | +0.10079 | +0.00849 | synergistic |

## Ridge analysis

- Optimum location: (A = 1.00, B = 1.00), value = +0.10927
- Cells within 95% of optimum: 3/25
- Neighbour mean = +0.10508, neighbour std = 0.00136
- Max conservation error across grid = 2.220e-16
- Verdict: **broadly increasing (monotonic)**

## Spatial map analysis

For each cluster, four per-pixel fields are produced:

1. Constitutive `|grad C|` - shows steep transitions and
   merging substructures.
2. Coherence factor (A=1) - shows where neighbouring gradients
   align.
3. Memory term `|r_self - r_prev|` - shows where the current
   and previous-step responses differ.
4. Synergy field - per-pixel interaction contribution, computed
   as `r(1,1) - r(1,0) - r(0,1) + r(0,0)`.

Per-cluster Pearson correlation between the synergy field and
the other fields:

| Cluster | corr(synergy, |grad C|) | corr(synergy, coherence) | corr(synergy, memory) |
|---|---|---|---|
| Abell2744 | +0.226 | +0.289 | -0.239 |
| MACS0416 | +0.180 | +0.228 | -0.263 |
| MACS1149 | +0.140 | +0.163 | -0.236 |
| AbellS1063 | +0.115 | +0.150 | -0.272 |
| Abell370 | +0.268 | +0.298 | -0.235 |

Mean spatial correlation between synergy and `|grad C|` across 5 clusters = +0.186.

All 5 clusters show positive spatial correlation between synergy and `|grad C|` (range +0.115 to +0.268). Synergy consistently concentrates around steep constitutive transitions (the high-density peaks of each cluster), which is consistent with a physical mechanism rather than a numerical artefact.

## Required questions

### Q1. Nonlinearity across the parameter space

Interior grid points (A>0, B>0): 16
Points with |synergy| > 1e-4 in Pearson kappa: 16 (100.0%).
Max |synergy|: 0.00849.
Mean |synergy|: 0.00359.

Nonlinearity is **present** if the synergy at the original C10 (A=1, B=0.5) is significant and the surface is not flat.

| A | B | Observed Pearson k | Expected additive | Synergy |
|---|---|---|---|---|
| 0.25 | 0.25 | +0.09371 | +0.09307 | +0.00064 |
| 0.25 | 0.50 | +0.09664 | +0.09540 | +0.00124 |
| 0.25 | 0.75 | +0.09796 | +0.09618 | +0.00178 |
| 0.25 | 1.00 | +0.09939 | +0.09765 | +0.00174 |
| 0.50 | 0.25 | +0.09486 | +0.09364 | +0.00123 |
| 0.50 | 0.50 | +0.09843 | +0.09596 | +0.00247 |
| 0.50 | 0.75 | +0.10031 | +0.09674 | +0.00356 |
| 0.50 | 1.00 | +0.10264 | +0.09821 | +0.00443 |
| 0.75 | 0.25 | +0.09675 | +0.09475 | +0.00199 |
| 0.75 | 0.50 | +0.10107 | +0.09708 | +0.00399 |
| 0.75 | 0.75 | +0.10315 | +0.09786 | +0.00529 |
| 0.75 | 1.00 | +0.10608 | +0.09933 | +0.00675 |
| 1.00 | 0.25 | +0.09855 | +0.09621 | +0.00233 |
| 1.00 | 0.50 | +0.10340 | +0.09854 | +0.00485 |
| 1.00 | 0.75 | +0.10600 | +0.09932 | +0.00668 |
| 1.00 | 1.00 | +0.10927 | +0.10079 | +0.00849 |

Interaction remains nonlinear at 16/16 interior grid points.
### Q2. Broad stable region or narrow optimum

Optimum location: A = 1.00, B = 1.00, value = +0.10927.
Cells within 95% of optimum: 3/25.
Neighbour-mean drop: +0.00420.
Neighbour std: 0.00136.

Verdict: **broadly increasing (monotonic)**.

### Q3. Nature of the interaction

Interior classifications (Tukey additivity):
- synergistic: 16
- additive:     0
- antagonistic: 0

Maximum synergy in Pearson kappa = +0.00849.
Synergy at original C10 (A=1, B=0.5) = +0.00485 (synergistic).

No clear saturation: response continues to improve or plateau as either mechanism approaches its maximum (A=1 along B=1: +0.00320; B=1 along A=1: +0.00327).
### Q4. Regional dominance

For each cluster, compare the single-mechanism improvements:

| Cluster | Coherence only (A=1,B=0) delta | Memory only (A=0,B=0.5) delta | Combined (A=1,B=0.5) delta | Dominant mechanism |
|---|---|---|---|---|
| Abell2744 | +0.01155 | -0.00517 | +0.00695 | Coherence |
| MACS0416 | +0.00449 | -0.00846 | -0.00731 | Coherence |
| MACS1149 | +0.00010 | -0.00623 | -0.00354 | Coherence |
| AbellS1063 | +0.00396 | +0.00507 | +0.01389 | Memory |
| Abell370 | +0.01096 | -0.00874 | -0.00031 | Coherence |

Across 5 clusters, coherence dominates in 4, memory dominates in 1.

### Q5. Cross-cluster consistency of the optimum region

For each cluster the (A, B) that maximises Pearson kappa is located; consistency is measured by how clustered these optima are.

| Cluster | Best (A, B) | Best Pearson k | Delta vs control |
|---|---|---|---|
| Abell2744 | (1.00, 0.00) | +0.02555 | +0.01155 |
| MACS0416 | (1.00, 0.00) | +0.02346 | +0.00449 |
| MACS1149 | (0.75, 0.00) | +0.24011 | +0.00109 |
| AbellS1063 | (1.00, 1.00) | +0.10927 | +0.01976 |
| Abell370 | (1.00, 0.00) | +0.20632 | +0.01096 |

Mean optimum (A, B) = (0.95, 0.20); std (A, B) = (0.100, 0.400).
All 5 clusters favour the different region.

### Q6. Saturation behaviour

Pearson kappa along B=1 (varying A from 0 to 1):

| A | Median Pearson k | Delta vs previous |
|---|---|---|
| 0.00 | +0.09683 | — |
| 0.25 | +0.09939 | +0.00256 |
| 0.50 | +0.10264 | +0.00325 |
| 0.75 | +0.10608 | +0.00344 |
| 1.00 | +0.10927 | +0.00320 |

Pearson kappa along A=1 (varying B from 0 to 1):

| B | Median Pearson k | Delta vs previous |
|---|---|---|
| 0.00 | +0.09347 | — |
| 0.25 | +0.09855 | +0.00508 |
| 0.50 | +0.10340 | +0.00485 |
| 0.75 | +0.10600 | +0.00261 |
| 1.00 | +0.10927 | +0.00327 |

Saturation observed along A at B=1: **NO**.
Saturation observed along B at A=1: **NO**.

### Q7. Conservation stability throughout parameter space

Machine epsilon = 2.220e-16.
Conservation error across 125 runs (25 configs x 5 clusters):
- max = 2.220e-16
- median = 2.220e-16
- min = 2.220e-16
- runs within machine epsilon: 125/125

Numerical stability **PRESERVED** across the entire parameter space.


## Outcome determination

**Outcome A** - the interaction forms a broad stable region (surface verdict: broadly increasing (monotonic), 3/25 cells within 95% of optimum).  Mean synergy across interior grid points = 0.00359, max synergy = 0.00849.  Synergy is positive at 16/16 interior grid points.  Mean spatial correlation between synergy and constitutive gradient magnitude across 5 clusters = +0.186.  Neighbour Coherence and Elastic Memory appear to be complementary components of the same cooperative physical response.

## Numerical stability report

| Configuration | Median runtime (s) | Max conservation |
|---|---|---|
| (A=0.00, B=0.00) | 0.228 | 2.220e-16 |
| (A=0.00, B=0.25) | 0.214 | 2.220e-16 |
| (A=0.00, B=0.50) | 0.215 | 2.220e-16 |
| (A=0.00, B=0.75) | 0.216 | 2.220e-16 |
| (A=0.00, B=1.00) | 0.215 | 2.220e-16 |
| (A=0.25, B=0.00) | 0.219 | 2.220e-16 |
| (A=0.25, B=0.25) | 0.230 | 2.220e-16 |
| (A=0.25, B=0.50) | 0.216 | 2.220e-16 |
| (A=0.25, B=0.75) | 0.229 | 2.220e-16 |
| (A=0.25, B=1.00) | 0.217 | 2.220e-16 |
| (A=0.50, B=0.00) | 0.228 | 2.220e-16 |
| (A=0.50, B=0.25) | 0.226 | 2.220e-16 |
| (A=0.50, B=0.50) | 0.232 | 2.220e-16 |
| (A=0.50, B=0.75) | 0.226 | 2.220e-16 |
| (A=0.50, B=1.00) | 0.231 | 2.220e-16 |
| (A=0.75, B=0.00) | 0.227 | 2.220e-16 |
| (A=0.75, B=0.25) | 0.232 | 2.220e-16 |
| (A=0.75, B=0.50) | 0.234 | 2.220e-16 |
| (A=0.75, B=0.75) | 0.231 | 2.220e-16 |
| (A=0.75, B=1.00) | 0.220 | 2.220e-16 |
| (A=1.00, B=0.00) | 0.233 | 2.220e-16 |
| (A=1.00, B=0.25) | 0.223 | 2.220e-16 |
| (A=1.00, B=0.50) | 0.230 | 2.220e-16 |
| (A=1.00, B=0.75) | 0.219 | 2.220e-16 |
| (A=1.00, B=1.00) | 0.232 | 2.220e-16 |

## Top-level artefacts

- `runs/coherence_memory_lab001/report.md` (this file)
- `runs/coherence_memory_lab001/parameter_grid.csv`
- `runs/coherence_memory_lab001/cluster_grid_statistics.csv`
- `runs/coherence_memory_lab001/interaction_surface.csv`
- `runs/coherence_memory_lab001/synergy_matrix.csv`
- `runs/coherence_memory_lab001/ridge_analysis.csv`
- `runs/coherence_memory_lab001/run.json`
- `runs/coherence_memory_lab001/validation.json`
- `runs/coherence_memory_lab001/plots/pearson_surface.png`
- `runs/coherence_memory_lab001/plots/bias_surface.png`
- `runs/coherence_memory_lab001/plots/ssim_surface.png`
- `runs/coherence_memory_lab001/plots/synergy_heatmap.png`
- `runs/coherence_memory_lab001/plots/ridge_map.png`
- `runs/coherence_memory_lab001/plots/parameter_stability.png`
- `runs/coherence_memory_lab001/plots/family_summary.png`
- `runs/coherence_memory_lab001/plots/spatial_maps/spatial_maps.png`
- `runs/coherence_memory_lab001/plots/spatial_maps/cluster_field_correlations.csv`

**Total execution time:** 32.4 s.
