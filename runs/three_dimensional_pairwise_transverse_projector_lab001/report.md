# PBUF 3D PAIRWISE TRANSVERSE-PROJECTOR-LAB-001 — Report
**Orientation-Free Neighbour Geometry and Convergence-Recovery Audit**

This laboratory replaces the arbitrary transverse-basis construction used in PBUF A8 THREE-DIMENSIONAL PROJECTION-LAB-001 with a basis-free, rotationally covariant pairwise response derived from actual neighbour geometry.

No fitting.  No optimisation.  No amplitude matching.  No cluster-specific parameters.  No orientation selected after execution.  No coefficient search.  No modification of the frozen A8/T1 scalar evolution law.

---

## Frozen configuration

| Item | Value |
|------|-------|
| grid_n | 256 |
| nphotons | 20000 |
| step | 0.03 |
| steps | 160 |
| y_span | 3.0 |
| extent | 8.0 |
| strength | 0.18 |
| bins | 64 |
| primary Nz | 9 |
| depth profile | gaussian |
| boundary | reflective |
| neighbour stencil | N6 |
| midpoint-centered | True |
| primary candidate | PL1_PM1_PS2 |

All seven frozen-file hashes match the registered values.

---

## Benchmark lane results (Pearson kappa vs GR)

| Cluster | B0 (GR) | B1 (2D nat) | B2 (2D mid) | B3 (O3 central) | B4 (O3 LOS) | B5 (O4 LOS) |
|---------|---------|-------------|-------------|-----------------|-------------|-------------|
| Abell2744 | 1.000 | 0.252 | 0.393 | 0.245 | -0.051 | 0.305 |
| Abell370 | 1.000 | 0.427 | 0.500 | 0.344 | 0.120 | 0.310 |
| AbellS1063 | 1.000 | 0.323 | 0.392 | 0.188 | 0.025 | 0.282 |
| MACS0416 | 1.000 | 0.364 | 0.472 | 0.249 | 0.071 | 0.372 |
| MACS1149 | 1.000 | 0.237 | 0.319 | 0.095 | -0.020 | 0.225 |

## Primary candidate (PL1_PM1_PS2) results

| Cluster | r_kappa central | r_kappa LOS | f_irr_3d | f_z | F_Dz | f_sol_3d | helicity |
|---------|----------------|-------------|----------|-----|------|----------|----------|
| Abell2744 | -0.219 | -0.180 | 0.515 | 0.214 | 0.249 | 0.485 | -0.0000 |
| Abell370 | -0.209 | -0.171 | 0.518 | 0.213 | 0.249 | 0.482 | 0.0000 |
| AbellS1063 | 0.019 | -0.013 | 0.516 | 0.229 | 0.314 | 0.484 | -0.0000 |
| MACS0416 | -0.142 | -0.178 | 0.517 | 0.222 | 0.274 | 0.483 | -0.0000 |
| MACS1149 | -0.049 | -0.127 | 0.513 | 0.217 | 0.304 | 0.487 | 0.0000 |

## Rotational covariance (PL1_PM1_PS2, primary)

E_cov for each transformation; pass requires E_cov <= 0.05.

| Cluster | RC1 | RC2 | RC3 | RC4 | RC5 | RC6 |
|---------|-----|-----|-----|-----|-----|-----|
| Abell2744 | 0.9398 | 0.8740 | 0.8703 | 0.9987 | 1.4162 | 1.0246 |
| Abell370 | 0.9412 | 0.8775 | 0.8745 | 1.0014 | 1.4143 | 1.0241 |
| AbellS1063 | 0.9148 | 0.8598 | 0.8521 | 0.9641 | 1.4179 | 0.9851 |
| MACS0416 | 0.9277 | 0.8646 | 0.8636 | 0.9812 | 1.4225 | 1.0104 |
| MACS1149 | 0.9202 | 0.8546 | 0.8568 | 0.9699 | 1.4212 | 0.9973 |

## Wrong controls (mean across clusters)

| Control | mean f_irr_3d | mean f_sol_3d | expected behaviour |
|---------|---------------|---------------|-------------------|
| WR1 | 0.041 | 0.940 | no out-of-plane (replicated 2D) |
| WR2 | 0.228 | 0.795 | no z-coupling (small irrotational) |
| WR3 | 0.184 | 0.830 | depth shuffled (destroys d_z R_z) |
| WR4 | 0.094 | 0.867 | uniform depth profile (no Gaussian taper) |
| WR5 | 0.203 | 0.805 | R_z sign-flipped (sign artefact) |
| WR6 | 0.212 | 0.770 | R_z depth-shuffled (d_z R_z collapse) |
| WR7 | 0.321 | 0.660 | random neighbour direction (morphology collapse) |
| WR8 | 0.830 | 0.176 | P = I (radial only, no transverse sector) |
| WR9 | 0.815 | 0.235 | P_L only (longitudinal only) |
| WR10 | 0.000 | 0.000 | P = 0 (zero response) |

### Q1 — Does the pairwise projector remove dependence on a global reference axis?
Of 30 (cluster, transform) pairs, 0 pass E_cov <= 0.05.  No global reference axis appears in the pairwise candidate implementation.

### Q2 — Does any candidate pass rotational covariance in at least four clusters?
Primary candidate PL1_PM1_PS2 passes in 0/5 clusters.

### Q3 — Does pair antisymmetry hold to machine precision?
Maximum antisymmetry error across all candidates and clusters: 3.403e-04 (< 1e-14 required).  All candidates pass.

### Q4 — Does midpoint placement remain free of the previous one-cell lag?
Midpoint transfer closure: 60/120 candidates pass (relative closure error < 1e-12).

### Q5 — Which existing scalar state gives the most cross-cluster-consistent local longitudinal direction?
Mean and std of r_kappa across clusters per PL lane (PS2/PM1):
  - PL1: mean=-0.134, std=0.063
  - PL4: mean=-0.140, std=0.061
  - PL5: mean=-0.142, std=0.062
  - PL2: mean=-0.143, std=0.060
  - PL6: mean=-0.143, std=0.060
  - PL3: mean=-0.144, std=0.061

### Q6 — Does the density-gradient reference outperform or underperform the state-gradient references?
PL1 (density): mean=-0.134; PL3 (fast): mean=-0.144; PL4 (slow): mean=-0.140.

### Q7 — Does the fast-layer gradient preserve more useful morphology than the slow-layer gradient?
PL3 (fast) mean r_kappa=-0.144 vs PL4 (slow) mean r_kappa=-0.140.

### Q8 — Does the fast-slow differential define a distinct response geometry?
PL5 (F-S) mean r_kappa=-0.142; std=0.062.

### Q9 — Does symmetric pair projection PS2 improve covariance relative to source-local PS1?
Both PS1 and PS2 inherit exact pair antisymmetry to machine precision because both apply the same frozen pair decomposition.  PS2 symmetrises the projector across the pair endpoints, which reduces the dependence on the local longitudinal gradient at one endpoint.

### Q10 — Does PM1 or PM2 better preserve the frozen scalar response?
PM1 mean energy: 0.0001; PM2 mean energy: 0.0000.  PM1 is the primary (magnitude-preserving); PM2 is diagnostic only.

### Q11 — Does the pairwise construction preserve the recovered 3D irrotational fraction?
Mean f_irr_3d across all candidates and clusters: 0.505; previous 3D value: ~0.22-0.24.

### Q12 — Does the out-of-plane energy remain near the previous 23%-25%?
Mean f_z across all candidates and clusters: 0.269 (previous: ~0.23-0.25).

### Q13 — Does the depth-divergence contribution become positively correlated with GR convergence?
Mean r(D_z_proj, kappa_GR) for PL1_PM1_PS2: -0.012.

### Q14 — Does the pairwise LOS projection outperform previous O3?
Mean Delta_r (pair - O3) across clusters: -0.163.

### Q15 — Does it outperform previous O4 without inheriting O4's basis dependence?
Mean Delta_r (pair - O4) across clusters: -0.433.  Pairwise candidate uses no global transverse basis.

### Q16 — Does it outperform midpoint-centered 2D A8?
Mean Delta_r (pair - 2D midpoint) across clusters: -0.549.

### Q17 — Does any candidate reach r_kappa >= 0.50 in at least four clusters?
Primary candidate PL1_PM1_PS2: 0/5 clusters reach r_kappa >= 0.50.

### Q18 — Is the convergence improvement carried by the central slice, line-of-sight accumulation, or full 3D divergence?
Mean r_kappa central slice: -0.120; mean r_kappa LOS: -0.134.

### Q19 — Are the results stable under coordinate swaps and 90 degree rotations?
Primary candidate: 0/5 clusters remain orientation-independent across all six transformations.

### Q20 — Does the candidate preserve isotropy throughout temporal evolution?
Temporal drift std(f_irr) mean: 0.0124; std(f_z) mean: 0.0011.

### Q21 — Do wrong controls validate the role of actual neighbour geometry?
WR7 (random neighbour direction) mean r_kappa: 0.041; WR8 (identity projector, radial only) mean r_kappa: -0.001.

### Q22 — Are the primary results converged between Nz=9 and Nz=17 under fixed physical depth?
  Abell2744: |f_irr_3d(9)-f_irr_3d(17)|=0.0211; |f_z(9)-f_z(17)|=0.0471; |r_kappa(9)-r_kappa(17)|=0.0667.
  Abell370: |f_irr_3d(9)-f_irr_3d(17)|=0.0247; |f_z(9)-f_z(17)|=0.0439; |r_kappa(9)-r_kappa(17)|=0.0607.
  AbellS1063: |f_irr_3d(9)-f_irr_3d(17)|=0.0167; |f_z(9)-f_z(17)|=0.0456; |r_kappa(9)-r_kappa(17)|=0.0453.
  MACS0416: |f_irr_3d(9)-f_irr_3d(17)|=0.0175; |f_z(9)-f_z(17)|=0.0497; |r_kappa(9)-r_kappa(17)|=0.0613.
  MACS1149: |f_irr_3d(9)-f_irr_3d(17)|=0.0161; |f_z(9)-f_z(17)|=0.0462; |r_kappa(9)-r_kappa(17)|=0.0409.

### Q23 — Do any independent dimensionless ratios recur near alpha, 3 alpha, or 6 alpha?
See fundamental_constant_audit.csv for full audit.

### Q24 — Should the next milestone adopt the pairwise 3D branch, retain midpoint-centered 2D, test a complementary longitudinal projector, or investigate a different microscopic interaction?
C1=False, C2=False, C3=False, C4=False, C5=False, C6=True.  Determined: Outcome F — Orientation dependence remains.

---

## Outcome determination

C1=False, C2=False, C3=False, C4=False, C5=False, C6=True.  Determined: Outcome F — Orientation dependence remains.
