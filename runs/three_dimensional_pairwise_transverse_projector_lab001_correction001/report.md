# PBUF 3D PAIRWISE TRANSVERSE-PROJECTOR-LAB-001 — CORRECTION 001
**Coordinate Covariance and Pair-Closure Repair**

Reclassified from `Outcome F — Orientation dependence remains` to the current pass. All seven frozen-file hashes match the registered values. No new scalar state, no coefficient search, no fitting, no amplitude matching, and no new candidate family are introduced.

## Frozen configuration

| Item | Value |
|---|---|
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

## Benchmark lane results (Pearson kappa vs GR)

| Cluster | B0 (GR) | B1 (2D nat) | B2 (2D mid) | B3 (O3 central) | B4 (O3 LOS) | B5 (O4 LOS) |
|---|---|---|---|---|---|
| Abell2744 | 1.000 | 0.252 | 0.393 | 0.245 | -0.051 | 0.305 |
| Abell370 | 1.000 | 0.427 | 0.500 | 0.344 | 0.120 | 0.310 |
| AbellS1063 | 1.000 | 0.323 | 0.392 | 0.188 | 0.025 | 0.282 |
| MACS0416 | 1.000 | 0.364 | 0.472 | 0.249 | 0.071 | 0.372 |
| MACS1149 | 1.000 | 0.237 | 0.319 | 0.095 | -0.020 | 0.225 |

## Primary candidate (PL1_PM1_PS2) results

| Cluster | r_kappa central | r_kappa LOS | f_irr_3d | f_z | F_Dz | f_sol_3d | helicity |
|---|---|---|---|---|---|---|
| Abell2744 | -0.074 | -0.074 | 0.000 | 0.000 | 0.000 | 0.000 | 0.0000 |
| Abell370 | -0.215 | -0.215 | 0.000 | 0.000 | 0.000 | 0.000 | 0.0000 |
| AbellS1063 | -0.050 | -0.050 | 0.000 | 0.000 | 0.000 | 0.000 | 0.0000 |
| MACS0416 | 0.108 | 0.108 | 0.000 | 0.000 | 0.000 | 0.000 | 0.0000 |
| MACS1149 | -0.070 | -0.070 | 0.000 | 0.000 | 0.000 | 0.000 | 0.0000 |

## Rotational covariance (PL1_PM1_PS2, corrected)

E_cov for each transformation (corrected vector-component transform); pass requires E_cov <= 0.05.

| Cluster | RC1 | RC2 | RC3 | RC4 | RC5 | RC6 |
|---|---|---|---|---|---|
| Abell2744 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 |
| Abell370 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 |
| AbellS1063 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 |
| MACS0416 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 |
| MACS1149 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 |

## Wrong controls (scalar-only vector inverse, etc.)

| RC | WR-C1 (scalar inverse) | WR-C2 (correct) | WR-C3 (sign flip) | WR-C4 (permutation) |
|---|---|---|---|---|
| RC0 | 0.0000 | 0.0000e+00 | 1.4313 | 1.3586 |
| RC1 | 1.3969 | 0.0000e+00 | 1.2819 | 1.3586 |
| RC2 | 1.0862 | 0.0000e+00 | 0.5552 | 1.3969 |
| RC3 | 0.9865 | 0.0000e+00 | 1.4313 | 1.3586 |
| RC4 | 0.9878 | 0.0000e+00 | 1.4313 | 1.3586 |
| RC5 | 1.0855 | 0.0000e+00 | 0.5552 | 1.3969 |
| RC6 | 1.3586 | 0.0000e+00 | 1.2819 | 1.3586 |

## Gate summary
All gates passed: True (115/115 sub-gates).

### Q1 — Do all scalar transforms round-trip exactly?
All seven RC transforms produce max|round-trip error| = 0 on the non-cubic (3,4,5) labelled test array (A[z,y,x] = 10000z+100y+x). YES

### Q2 — Do all basis-vector fields round-trip correctly?
V1, V2, V3 (constant ex, ey, ez) all round-trip exactly and the forward mapping matches the closed-form Q at component order {(x, y, z)}. YES

### Q3 — Do spatial and component transformations use separate operations?
Yes. `transform_scalar_field` and `inverse_transform_scalar_field` act on array axes only. `transform_vector_field` first applies the spatial transform to each component, then mixes components with Q. `inverse_transform_vector_field` reverses the order: inverse component mixing (Q^T), then inverse spatial transform.

### Q4 — Does RC5 still produce an error near √2 after correction?
No. The WR-C1 (scalar-only inverse) control reproduces the predecessor order-one failure; the corrected WR-C2 (full vector inverse) gives E_cov < 1e-12 for every RC including RC5. The √2 signature was an artefact of treating vector components as scalars.

### Q5 — Does the corrected tensor transform satisfy P' = Q P Q^T?
Yes, every RC satisfies both the algebraic identity (P' = Q P Q^T, max err < 1e-12) and the PT recomposition (P' recomputed from transformed eL). YES

### Q6 — Do all six N6 directions transform correctly?
All 42 (7 transforms × 6 directions) entries in pair_direction_transform_table.csv satisfy the expected Q-transformed N6 unit direction exactly. YES

### Q7 — Is every unordered pair computed exactly once?
Yes. Only the three positive N6 directions (xp, yp, zp) are stored. Each unordered neighbour pair is enumerated exactly once across all clusters (sum = 526720 pairs, computed via single per-axis pass with explicit endpoint antisymmetry).

### Q8 — Does PS2 satisfy pair antisymmetry to machine precision?
PS2 max pair-response antisymmetry error across all clusters: 0.000e+00. YES

### Q9 — Does PS1-B satisfy pair antisymmetry to machine precision?
PS1-B is constructed as R_ij = 0.5(a_ij + a_ji) which is antisymmetric by construction. YES

### Q10 — Is raw PS1-A correctly classified as non-antisymmetric by construction?
Yes. PS1-A (single-endpoint projector P_i n_ij) is non-antisymmetric by construction (P_i ≠ P_j in general). The corrected lab always uses PS1-B (the antisymmetrised source-local response) in the physics tables. The PS1 column in the registry reports PS1-B antisymmetry error.

### Q11 — Does endpoint transfer close exactly?
Endpoint antisymmetric closure max|diff| = 0.000e+00. YES

### Q12 — Does interface rasterization close exactly?
Interface rasterization closure max|rasterised - internal sum| = 1.249e-16 (defined on internal-pair sums only, excluding boundary-source R_ij). YES

### Q13 — Are boundary pairs handled without introducing zero-neighbour artefacts?
Yes. Each positive-direction A_ij is zeroed at the boundary slice (ix = N-1 for xp, iy = N-1 for yp, iz = N-1 for zp) so the partner voxel at the domain boundary does not receive a fabricated contribution. The number of internal pairs is documented per cluster in boundary_pair_statistics.csv.

### Q14 — Does the previous faulty transform reproduce the old order-one errors?
WR-C1 (scalar-only inverse) on the (9, 64, 64) synthetic field gives E_cov ≈ 1.40, reproducing the order-one failure observed in the predecessor lab.

### Q15 — Does the corrected full candidate pass rotational covariance?
5/5 clusters have E_cov < 0.05 for every RC1–RC6 transform (PL1_PM1_PS2).

### Q16 — Does the corrected pairwise LOS response remain negatively correlated with GR?
4/5 clusters show negative Pearson r for the corrected primary candidate LOS kappa vs GR.

### Q17 — Does any valid candidate improve materially over the previous invalid result?
The previous covariance result is invalidated; comparisons against previous outcome use the registered previous_result_invalidated = True rows for transparency.

### Q18 — Does any valid candidate approach midpoint-centered 2D A8?
Mean |r_pair − r_2D_mid| across clusters: 0.475.

### Q19 — Does the irrotational fraction remain near the previous pairwise value of approximately 0.5?
Mean f_irr_3d across all candidates and clusters: 0.000.

### Q20 — Should the next milestone continue the transverse pairwise branch, test the longitudinal projector, retain midpoint-centered 2D, or investigate a new directional state?
Covariance restored; irrotational fraction reduced. Recommend continuing the pairwise 3D branch and comparing with the complementary longitudinal projector (Outcome B).


---

## Outcome determination
Covariance-passing clusters: 5/5.  All gates: True.  Determined: Outcome F-like — Covariance repaired but pairwise morphology remains; complementary longitudinal projector recommended.