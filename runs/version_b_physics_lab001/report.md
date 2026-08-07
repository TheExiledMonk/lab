# PBUF VERSION-B PHYSICS-LAB-001

**Local response hypothesis survey using the frozen
Version 1 weak-lensing laboratory (LAB-FREEZE-001).**

No parameter fitting.  No optimisation.  Each candidate is
tested exactly once on each of the five benchmark clusters
with the frozen minimum production configuration.

## Status

- Frozen hash verification: **PASS**
- Total runtime: **12.5 s**
- Candidates tested: **10**
- Clusters: **5**

## Frozen laboratory

The Version 1 laboratory is used as the measurement
instrument without modification.

| Component | Frozen specification |
|---|---|
| Constitutive | `C(X) = 0.18 * rho(X) / rho_max` (Version A) |
| Transport | neighbour-to-neighbour, direct addition, |
| | per-step unit-speed renormalisation |
| Response direction | 90 deg transverse (R_90 of grad C) |
| Source plane | Launch B (Cartesian 2D) |
| Observable | Jacobian (ray-bundle linear fit per bin) |
| Matter input | `rho = max(kappa_obs, 0) / max(max(kappa_obs, 0))` |

## Production configuration

| Parameter | Value |
|---|---|
| Photons | 20,000 |
| Constitutive grid | 256^2 |
| Step size | Delta s / 2 = 0.0300 |
| Number of steps | 160 |
| Source plane | Cartesian 2D (Launch B) |
| Observable | Jacobian |

## Candidates

| # | Name | Family | Description |
|---|---|---|---|
| 1 | Gradient (control) | gradient | Response = |grad C|; frozen Version A control. |
| 2 | Local Neighbour Coherence | neighbour coherence | Magnitude scaled by (1 + mean_cos)/2 over 8 neighbours. |
| 3 | Cooperative Neighbour Response | cooperative response | Cell-averaged gradient (3x3 box). |
| 4 | Elastic Memory | elastic memory | r_new = (1-w)*R(g) + w*R(g_upstream); w = 0.5. |
| 5 | Gradient Curvature | gradient curvature | A = |grad C| + 0.5 * |Laplacian C|. |
| 6 | Phase-Coherent Response | phase coherence | A = |grad C| * mean_cos(phase differences). |
| 7 | Relaxation Response | relaxation | One Jacobi relaxation step toward neighbour mean of response. |
| 8 | Weak-Gradient Enhancement | weak-gradient enhancement | A = |grad C| + 0.05 * exp(-|grad C| / 0.05). |
| 9 | Constitutive Coupling | constitutive coupling | A = |grad C| * (1 + C). |
| 10 | Combined Local Response | combined response | Combine neighbour coherence (Cand 2) and elastic memory (Cand 4). |

All fixed parameters are documented in the candidate source
code.  No parameter is fitted.

## Per-candidate, per-cluster metrics

Computed metrics for every (candidate, cluster) pair:

| Candidate | Cluster | RMS k | RMS g | Pearson k | Pearson g | SSIM k | SSIM g | k bias | g bias | conservation | runtime (s) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| C1 Gradient (control) | Abell 2744 | 1.5580e-01 | 6.0822e-02 | +0.0140 | +0.0895 | -0.0106 | +0.1008 | -5.9392e-02 | +1.5525e-03 | 2.220e-16 | 0.232 |
| C1 Gradient (control) | MACS J0416 | 1.8596e-01 | 6.5276e-02 | +0.0190 | +0.1175 | +0.0004 | +0.1267 | -9.8489e-02 | -1.2579e-02 | 2.220e-16 | 0.230 |
| C1 Gradient (control) | MACS J1149 | 8.5592e-02 | 1.1556e-01 | +0.2390 | +0.0075 | +0.1081 | +0.0120 | -3.2439e-02 | -9.3065e-02 | 2.220e-16 | 0.233 |
| C1 Gradient (control) | Abell S1063 | 1.3540e-01 | 6.9710e-02 | +0.0895 | +0.0836 | -0.0383 | +0.0929 | -5.3933e-02 | -2.1999e-02 | 2.220e-16 | 0.229 |
| C1 Gradient (control) | Abell 370 | 1.7915e-01 | 9.1733e-02 | +0.1954 | -0.0171 | -0.0398 | +0.0122 | -6.7484e-02 | -1.5468e-02 | 2.220e-16 | 0.216 |
| C2 Local Neighbour Coherence | Abell 2744 | 1.4344e-01 | 5.7249e-02 | +0.0256 | +0.0860 | -0.0129 | +0.0973 | -5.8019e-02 | -5.3675e-03 | 2.220e-16 | 0.229 |
| C2 Local Neighbour Coherence | MACS J0416 | 1.7308e-01 | 6.3177e-02 | +0.0235 | +0.1136 | +0.0010 | +0.1239 | -9.7374e-02 | -1.9917e-02 | 2.220e-16 | 0.222 |
| C2 Local Neighbour Coherence | MACS J1149 | 8.3277e-02 | 1.1704e-01 | +0.2391 | +0.0082 | +0.0960 | +0.0112 | -3.3313e-02 | -9.5478e-02 | 2.220e-16 | 0.214 |
| C2 Local Neighbour Coherence | Abell S1063 | 1.2608e-01 | 6.8610e-02 | +0.0935 | +0.0810 | -0.0306 | +0.0869 | -5.0951e-02 | -2.6319e-02 | 2.220e-16 | 0.212 |
| C2 Local Neighbour Coherence | Abell 370 | 1.6048e-01 | 8.9565e-02 | +0.2063 | +0.0129 | -0.0292 | +0.0406 | -6.4601e-02 | -2.6507e-02 | 2.220e-16 | 0.215 |
| C3 Cooperative Neighbour Response | Abell 2744 | 1.4119e-01 | 5.3292e-02 | -0.0113 | +0.1334 | -0.0007 | +0.1449 | -5.3313e-02 | -5.8948e-03 | 2.220e-16 | 0.215 |
| C3 Cooperative Neighbour Response | MACS J0416 | 1.7109e-01 | 6.1086e-02 | +0.0083 | +0.1205 | -0.0002 | +0.1309 | -9.9727e-02 | -2.2080e-02 | 2.220e-16 | 0.216 |
| C3 Cooperative Neighbour Response | MACS J1149 | 8.5610e-02 | 1.1676e-01 | +0.1996 | -0.0569 | +0.0691 | -0.0026 | -3.4697e-02 | -9.4633e-02 | 2.220e-16 | 0.223 |
| C3 Cooperative Neighbour Response | Abell S1063 | 1.2120e-01 | 6.6337e-02 | +0.0939 | +0.0804 | -0.0399 | +0.0864 | -5.3185e-02 | -2.7592e-02 | 2.220e-16 | 0.216 |
| C3 Cooperative Neighbour Response | Abell 370 | 1.5836e-01 | 8.7458e-02 | +0.2168 | +0.0313 | -0.0321 | +0.0569 | -6.4800e-02 | -2.6289e-02 | 2.220e-16 | 0.216 |
| C4 Elastic Memory | Abell 2744 | 1.5551e-01 | 6.0280e-02 | +0.0088 | +0.0824 | -0.0093 | +0.0953 | -5.9732e-02 | -9.6621e-04 | 2.220e-16 | 0.220 |
| C4 Elastic Memory | MACS J0416 | 1.8648e-01 | 6.5138e-02 | +0.0105 | +0.1403 | +0.0001 | +0.1478 | -9.8932e-02 | -1.3677e-02 | 2.220e-16 | 0.214 |
| C4 Elastic Memory | MACS J1149 | 8.5737e-02 | 1.1673e-01 | +0.2328 | -0.0389 | +0.1053 | -0.0002 | -3.2457e-02 | -9.3967e-02 | 2.220e-16 | 0.217 |
| C4 Elastic Memory | Abell S1063 | 1.3382e-01 | 6.8874e-02 | +0.0946 | +0.0862 | -0.0403 | +0.0960 | -5.3841e-02 | -2.3095e-02 | 2.220e-16 | 0.215 |
| C4 Elastic Memory | Abell 370 | 1.7942e-01 | 9.0267e-02 | +0.1866 | +0.0012 | -0.0403 | +0.0298 | -6.7925e-02 | -1.9230e-02 | 2.220e-16 | 0.224 |
| C5 Gradient Curvature | Abell 2744 | 1.6100e-01 | 6.2657e-02 | +0.0107 | +0.0913 | -0.0094 | +0.1018 | -5.9274e-02 | +4.6482e-03 | 2.220e-16 | 0.218 |
| C5 Gradient Curvature | MACS J0416 | 1.9132e-01 | 6.6495e-02 | +0.0189 | +0.1176 | -0.0000 | +0.1259 | -9.9321e-02 | -9.5718e-03 | 2.220e-16 | 0.217 |
| C5 Gradient Curvature | MACS J1149 | 8.6722e-02 | 1.1488e-01 | +0.2365 | +0.0068 | +0.1117 | +0.0122 | -3.2072e-02 | -9.1954e-02 | 2.220e-16 | 0.219 |
| C5 Gradient Curvature | Abell S1063 | 1.3891e-01 | 7.0290e-02 | +0.0928 | +0.0834 | -0.0401 | +0.0937 | -5.4415e-02 | -1.9980e-02 | 2.220e-16 | 0.218 |
| C5 Gradient Curvature | Abell 370 | 1.8638e-01 | 9.2829e-02 | +0.1943 | -0.0158 | -0.0440 | +0.0123 | -6.8690e-02 | -1.1365e-02 | 2.220e-16 | 0.218 |
| C6 Phase-Coherent Response | Abell 2744 | 1.3457e-01 | 5.4379e-02 | +0.0303 | +0.0885 | -0.0098 | +0.0992 | -5.4669e-02 | -1.0335e-02 | 2.220e-16 | 0.215 |
| C6 Phase-Coherent Response | MACS J0416 | 1.6368e-01 | 6.1554e-02 | +0.0585 | +0.0966 | +0.0017 | +0.1077 | -9.7772e-02 | -2.4381e-02 | 2.220e-16 | 0.216 |
| C6 Phase-Coherent Response | MACS J1149 | 8.4099e-02 | 1.1686e-01 | +0.2181 | -0.0050 | +0.0750 | +0.0085 | -3.4610e-02 | -9.5548e-02 | 2.220e-16 | 0.215 |
| C6 Phase-Coherent Response | Abell S1063 | 1.1976e-01 | 6.7243e-02 | +0.0967 | +0.0757 | -0.0134 | +0.0815 | -4.6377e-02 | -2.8425e-02 | 2.220e-16 | 0.216 |
| C6 Phase-Coherent Response | Abell 370 | 1.5736e-01 | 8.8235e-02 | +0.1465 | +0.0311 | -0.0061 | +0.0556 | -6.0980e-02 | -3.1642e-02 | 2.220e-16 | 0.217 |
| C7 Relaxation Response | Abell 2744 | 1.4781e-01 | 5.6701e-02 | +0.0013 | +0.1107 | -0.0047 | +0.1225 | -5.6246e-02 | -2.4969e-03 | 2.220e-16 | 0.215 |
| C7 Relaxation Response | MACS J0416 | 1.7780e-01 | 6.2939e-02 | +0.0142 | +0.1177 | +0.0001 | +0.1283 | -9.8907e-02 | -1.7690e-02 | 2.220e-16 | 0.225 |
| C7 Relaxation Response | MACS J1149 | 8.5341e-02 | 1.1624e-01 | +0.2221 | -0.0202 | +0.0872 | +0.0050 | -3.3716e-02 | -9.4012e-02 | 2.220e-16 | 0.216 |
| C7 Relaxation Response | Abell S1063 | 1.2791e-01 | 6.8014e-02 | +0.0899 | +0.0796 | -0.0389 | +0.0868 | -5.3615e-02 | -2.5036e-02 | 2.220e-16 | 0.216 |
| C7 Relaxation Response | Abell 370 | 1.6807e-01 | 8.9543e-02 | +0.2048 | -0.0007 | -0.0359 | +0.0288 | -6.6043e-02 | -2.1303e-02 | 2.220e-16 | 0.215 |
| C8 Weak-Gradient Enhancement | Abell 2744 | 3.3985e-01 | 1.3327e-01 | -0.0600 | -0.1052 | +0.0095 | -0.0333 | -8.3702e-02 | +8.7493e-02 | 2.220e-16 | 0.218 |
| C8 Weak-Gradient Enhancement | MACS J0416 | 3.7056e-01 | 1.3802e-01 | -0.0531 | -0.0659 | -0.0083 | -0.0256 | -6.8712e-02 | +8.3984e-02 | 2.220e-16 | 0.220 |
| C8 Weak-Gradient Enhancement | MACS J1149 | 2.6924e-01 | 9.3041e-02 | +0.0677 | +0.2057 | +0.0068 | +0.2095 | -3.9525e-02 | -2.5691e-03 | 2.220e-16 | 0.217 |
| C8 Weak-Gradient Enhancement | Abell S1063 | 3.3532e-01 | 1.2790e-01 | +0.1060 | +0.0535 | -0.0361 | +0.0477 | -7.4845e-02 | +7.4901e-02 | 2.220e-16 | 0.216 |
| C8 Weak-Gradient Enhancement | Abell 370 | 2.9746e-01 | 1.1937e-01 | +0.0913 | -0.0273 | +0.0018 | -0.0070 | -5.9229e-02 | +4.1505e-02 | 2.220e-16 | 0.214 |
| C9 Constitutive Coupling | Abell 2744 | 1.5739e-01 | 6.1422e-02 | +0.0138 | +0.0913 | -0.0106 | +0.1021 | -5.9521e-02 | +2.4227e-03 | 2.220e-16 | 0.216 |
| C9 Constitutive Coupling | MACS J0416 | 1.8737e-01 | 6.5537e-02 | +0.0187 | +0.1185 | +0.0003 | +0.1275 | -9.8553e-02 | -1.1629e-02 | 2.220e-16 | 0.215 |
| C9 Constitutive Coupling | MACS J1149 | 8.5754e-02 | 1.1547e-01 | +0.2390 | +0.0068 | +0.1087 | +0.0119 | -3.2396e-02 | -9.2900e-02 | 2.220e-16 | 0.216 |
| C9 Constitutive Coupling | Abell S1063 | 1.3672e-01 | 7.0018e-02 | +0.0894 | +0.0842 | -0.0386 | +0.0939 | -5.4109e-02 | -2.1336e-02 | 2.220e-16 | 0.217 |
| C9 Constitutive Coupling | Abell 370 | 1.8127e-01 | 9.1833e-02 | +0.1945 | -0.0129 | -0.0407 | +0.0158 | -6.7769e-02 | -1.4233e-02 | 2.220e-16 | 0.216 |
| C10 Combined Local Response | Abell 2744 | 1.3990e-01 | 5.6182e-02 | +0.0210 | +0.0824 | -0.0117 | +0.0940 | -5.8042e-02 | -9.3295e-03 | 2.220e-16 | 0.226 |
| C10 Combined Local Response | MACS J0416 | 1.7082e-01 | 6.3011e-02 | +0.0117 | +0.1395 | +0.0009 | +0.1469 | -9.6652e-02 | -2.2600e-02 | 2.220e-16 | 0.216 |
| C10 Combined Local Response | MACS J1149 | 8.2467e-02 | 1.1850e-01 | +0.2355 | -0.0208 | +0.0958 | +0.0047 | -3.3028e-02 | -9.7154e-02 | 2.220e-16 | 0.231 |
| C10 Combined Local Response | Abell S1063 | 1.2270e-01 | 6.7860e-02 | +0.1034 | +0.0865 | -0.0338 | +0.0910 | -5.0881e-02 | -2.8686e-02 | 2.220e-16 | 0.216 |
| C10 Combined Local Response | Abell 370 | 1.5946e-01 | 8.8799e-02 | +0.1950 | +0.0357 | -0.0430 | +0.0593 | -6.7297e-02 | -3.1552e-02 | 2.220e-16 | 0.229 |

## Cross-cluster evaluation

For every candidate the following medians/means are taken
across the five benchmark clusters.

| Candidate | Median Pearson k | Median Pearson g | Median SSIM k | Mean k Bias | Mean g Bias | Mean Pearson k | Conservation max | Runtime (s) |
|---|---|---|---|---|---|---|---|---|
| C1 Gradient (control) | +0.0895 | +0.0836 | -0.0106 | -6.2347e-02 | -2.8312e-02 | +0.1114 | 2.220e-16 | 0.230 |
| C2 Local Neighbour Coherence | +0.0935 | +0.0810 | -0.0129 | -6.0852e-02 | -3.4718e-02 | +0.1176 | 2.220e-16 | 0.215 |
| C3 Cooperative Neighbour Response | +0.0939 | +0.0804 | -0.0007 | -6.1144e-02 | -3.5298e-02 | +0.1015 | 2.220e-16 | 0.216 |
| C4 Elastic Memory | +0.0946 | +0.0824 | -0.0093 | -6.2577e-02 | -3.0187e-02 | +0.1067 | 2.220e-16 | 0.217 |
| C5 Gradient Curvature | +0.0928 | +0.0834 | -0.0094 | -6.2755e-02 | -2.5645e-02 | +0.1106 | 2.220e-16 | 0.218 |
| C6 Phase-Coherent Response | +0.0967 | +0.0757 | -0.0061 | -5.8882e-02 | -3.8066e-02 | +0.1100 | 2.220e-16 | 0.216 |
| C7 Relaxation Response | +0.0899 | +0.0796 | -0.0047 | -6.1705e-02 | -3.2107e-02 | +0.1065 | 2.220e-16 | 0.216 |
| C8 Weak-Gradient Enhancement | +0.0677 | -0.0273 | +0.0018 | -6.5203e-02 | +5.7063e-02 | +0.0304 | 2.220e-16 | 0.217 |
| C9 Constitutive Coupling | +0.0894 | +0.0842 | -0.0106 | -6.2470e-02 | -2.7535e-02 | +0.1111 | 2.220e-16 | 0.216 |
| C10 Combined Local Response | +0.1034 | +0.0824 | -0.0117 | -6.1180e-02 | -3.7864e-02 | +0.1133 | 2.220e-16 | 0.226 |

## Candidate ranking (by median Pearson kappa)

| Rank | Candidate | Family | Median Pearson k | Median SSIM k | Mean k Bias | Delta Pearson vs control | Delta SSIM vs control |
|---|---|---|---|---|---|---|---|
| 1 | C10 Combined Local Response | combined response | +0.1034 | -0.0117 | -6.1180e-02 | +0.0139 | -0.0012 |
| 2 | C6 Phase-Coherent Response | phase coherence | +0.0967 | -0.0061 | -5.8882e-02 | +0.0072 | +0.0045 |
| 3 | C4 Elastic Memory | elastic memory | +0.0946 | -0.0093 | -6.2577e-02 | +0.0051 | +0.0013 |
| 4 | C3 Cooperative Neighbour Response | cooperative response | +0.0939 | -0.0007 | -6.1144e-02 | +0.0044 | +0.0099 |
| 5 | C2 Local Neighbour Coherence | neighbour coherence | +0.0935 | -0.0129 | -6.0852e-02 | +0.0040 | -0.0024 |
| 6 | C5 Gradient Curvature | gradient curvature | +0.0928 | -0.0094 | -6.2755e-02 | +0.0033 | +0.0012 |
| 7 | C7 Relaxation Response | relaxation | +0.0899 | -0.0047 | -6.1705e-02 | +0.0004 | +0.0059 |
| 8 | C1 Gradient (control) | gradient | +0.0895 | -0.0106 | -6.2347e-02 | +0.0000 | +0.0000 |
| 9 | C9 Constitutive Coupling | constitutive coupling | +0.0894 | -0.0106 | -6.2470e-02 | -0.0001 | -0.0000 |
| 10 | C8 Weak-Gradient Enhancement | weak-gradient enhancement | +0.0677 | +0.0018 | -6.5203e-02 | -0.0218 | +0.0124 |

## Required questions

### Q1. Reduction of the systematic kappa underprediction

Control (C1) mean kappa bias = -0.06235.  A negative bias means the
frozen Version A laboratory underpredicts kappa.  A reduction
of this bias (toward zero or positive) is the success criterion.

| Candidate | Mean kappa bias | Delta vs control |
|---|---|---|
| C6 Phase-Coherent Response | -0.05888 | +0.00347 |
| C2 Local Neighbour Coherence | -0.06085 | +0.00150 |
| C3 Cooperative Neighbour Response | -0.06114 | +0.00120 |
| C10 Combined Local Response | -0.06118 | +0.00117 |
| C7 Relaxation Response | -0.06171 | +0.00064 |
| C9 Constitutive Coupling | -0.06247 | -0.00012 |
| C4 Elastic Memory | -0.06258 | -0.00023 |
| C5 Gradient Curvature | -0.06275 | -0.00041 |
| C8 Weak-Gradient Enhancement | -0.06520 | -0.00286 |

Best improvement: C6 Phase-Coherent Response (mean kappa bias = -0.05888, delta = +0.00347).
Number of candidates that *reduce* the magnitude of the bias: 5/9.

### Q2. Largest improvement

Control (C1) median Pearson kappa = +0.0895, median SSIM kappa = -0.0106, mean kappa bias = -0.06235.

**Largest median Pearson kappa:** C10 Combined Local Response (+0.1034, delta = +0.0139).
**Largest median SSIM kappa:** C8 Weak-Gradient Enhancement (+0.0018, delta = +0.0124).
**Lowest |mean kappa bias|:** C6 Phase-Coherent Response (-0.05888, delta = +0.00347).

Top-3 by median Pearson kappa:

| Rank | Candidate | Median Pearson kappa | Median SSIM kappa | Mean kappa bias |
|---|---|---|---|---|
| 1 | C10 Combined Local Response | +0.1034 | -0.0117 | -0.06118 |
| 2 | C6 Phase-Coherent Response | +0.0967 | -0.0061 | -0.05888 |
| 3 | C4 Elastic Memory | +0.0946 | -0.0093 | -0.06258 |

### Q3. Numerical stability impact

Control (C1) maximum conservation error = 2.220e-16, median runtime = 0.230 s.
Machine epsilon = 2.220446049250313e-16.

| Candidate | Max conservation | Runtime (s) | Status |
|---|---|---|---|
| C1 Gradient (control) | 2.220e-16 | 0.230 | machine-epsilon |
| C2 Local Neighbour Coherence | 2.220e-16 | 0.215 | machine-epsilon |
| C3 Cooperative Neighbour Response | 2.220e-16 | 0.216 | machine-epsilon |
| C4 Elastic Memory | 2.220e-16 | 0.217 | machine-epsilon |
| C5 Gradient Curvature | 2.220e-16 | 0.218 | machine-epsilon |
| C6 Phase-Coherent Response | 2.220e-16 | 0.216 | machine-epsilon |
| C7 Relaxation Response | 2.220e-16 | 0.216 | machine-epsilon |
| C8 Weak-Gradient Enhancement | 2.220e-16 | 0.217 | machine-epsilon |
| C9 Constitutive Coupling | 2.220e-16 | 0.216 | machine-epsilon |
| C10 Combined Local Response | 2.220e-16 | 0.226 | machine-epsilon |

Number of candidates that preserve the machine-precision conservation bound: 10/10.

No candidate exceeds the machine-precision conservation bound observed for the frozen control.

### Q4. Machine-precision conservation preservation

Machine epsilon = 2.220e-16.
Number of candidates that satisfy `max conservation <= machine epsilon`: 10/10.

| Candidate | Max conservation | Preserves |
|---|---|---|
| C1 Gradient (control) | 2.220e-16 | YES |
| C2 Local Neighbour Coherence | 2.220e-16 | YES |
| C3 Cooperative Neighbour Response | 2.220e-16 | YES |
| C4 Elastic Memory | 2.220e-16 | YES |
| C5 Gradient Curvature | 2.220e-16 | YES |
| C6 Phase-Coherent Response | 2.220e-16 | YES |
| C7 Relaxation Response | 2.220e-16 | YES |
| C8 Weak-Gradient Enhancement | 2.220e-16 | YES |
| C9 Constitutive Coupling | 2.220e-16 | YES |
| C10 Combined Local Response | 2.220e-16 | YES |

### Q5. Cross-cluster consistency of improvement

Improvement is consistent across clusters if every cluster
exhibits a positive delta in Pearson kappa versus the control.

| Candidate | Median delta Pearson k | # clusters +ve | # clusters -ve | Sign consistent |
|---|---|---|---|---|
| C6 Phase-Coherent Response | +0.0072 | 3 | 2 | NO |
| C2 Local Neighbour Coherence | +0.0045 | 5 | 0 | YES |
| C9 Constitutive Coupling | -0.0002 | 0 | 5 | NO |
| C10 Combined Local Response | -0.0003 | 2 | 3 | NO |
| C5 Gradient Curvature | -0.0011 | 1 | 4 | NO |
| C7 Relaxation Response | -0.0047 | 2 | 3 | NO |
| C4 Elastic Memory | -0.0062 | 1 | 4 | NO |
| C3 Cooperative Neighbour Response | -0.0106 | 2 | 3 | NO |
| C8 Weak-Gradient Enhancement | -0.0740 | 1 | 4 | NO |

Number of candidates with sign-consistent improvement on all 5 clusters: 1/9.

### Q6. Best-performing physical family

| Rank | Family | Best candidate | Median Pearson k | Median SSIM k | Mean k bias |
|---|---|---|---|---|---|
| 1 | combined response | C10 Combined Local Response | +0.1034 | -0.0117 | -0.06118 |
| 2 | phase coherence | C6 Phase-Coherent Response | +0.0967 | -0.0061 | -0.05888 |
| 3 | elastic memory | C4 Elastic Memory | +0.0946 | -0.0093 | -0.06258 |
| 4 | cooperative response | C3 Cooperative Neighbour Response | +0.0939 | -0.0007 | -0.06114 |
| 5 | neighbour coherence | C2 Local Neighbour Coherence | +0.0935 | -0.0129 | -0.06085 |
| 6 | gradient curvature | C5 Gradient Curvature | +0.0928 | -0.0094 | -0.06275 |
| 7 | relaxation | C7 Relaxation Response | +0.0899 | -0.0047 | -0.06171 |
| 8 | gradient | C1 Gradient (control) | +0.0895 | -0.0106 | -0.06235 |
| 9 | constitutive coupling | C9 Constitutive Coupling | +0.0894 | -0.0106 | -0.06247 |
| 10 | weak-gradient enhancement | C8 Weak-Gradient Enhancement | +0.0677 | +0.0018 | -0.06520 |

**Best family:** `combined response` (best candidate: C10 Combined Local Response, median Pearson kappa = +0.1034).


## Outcome determination

**Outcome A (weak)** - 1 candidate(s) show sign-consistent improvement on all 5 clusters while preserving machine-precision conservation, but the median delta is below 0.01 in absolute value.  No candidate produces a strong, consistent improvement.  The frozen laboratory is not modified; the result suggests the missing physics lies outside the tested local-response families.

## Numerical stability report

| Candidate | Mean runtime (s) | Max conservation |
|---|---|---|
| C1 Gradient (control) | 0.230 | 2.220e-16 |
| C2 Local Neighbour Coherence | 0.215 | 2.220e-16 |
| C3 Cooperative Neighbour Response | 0.216 | 2.220e-16 |
| C4 Elastic Memory | 0.217 | 2.220e-16 |
| C5 Gradient Curvature | 0.218 | 2.220e-16 |
| C6 Phase-Coherent Response | 0.216 | 2.220e-16 |
| C7 Relaxation Response | 0.216 | 2.220e-16 |
| C8 Weak-Gradient Enhancement | 0.217 | 2.220e-16 |
| C9 Constitutive Coupling | 0.216 | 2.220e-16 |
| C10 Combined Local Response | 0.226 | 2.220e-16 |

## Top-level artefacts

- `runs/version_b_physics_lab001/report.md` (this file)
- `runs/version_b_physics_lab001/candidate_summary.csv`
- `runs/version_b_physics_lab001/cross_cluster_statistics.csv`
- `runs/version_b_physics_lab001/candidate_ranking.csv`
- `runs/version_b_physics_lab001/run.json`
- `runs/version_b_physics_lab001/validation.json`
- `runs/version_b_physics_lab001/plots/candidate_rankings.png`
- `runs/version_b_physics_lab001/plots/bias_comparison.png`
- `runs/version_b_physics_lab001/plots/pearson_comparison.png`
- `runs/version_b_physics_lab001/plots/ssim_comparison.png`
- `runs/version_b_physics_lab001/plots/cluster_performance.png`
- `runs/version_b_physics_lab001/plots/response_family_summary.png`

**Total execution time:** 12.5 s.
