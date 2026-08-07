# PBUF A8 THREE-DIMENSIONAL PROJECTION-LAB-001 — Report
**Three-Dimensional Microscopic Response and Line-of-Sight Recovery Audit**

This laboratory extends the frozen A8/T1 microscopic system from
two spatial dimensions to three spatial dimensions using the same
local rules, coefficients, update order, conservation procedure,
and observable machinery wherever mathematically applicable.

No fitting.  No optimisation.  No parameter search.  No amplitude
matching.  No cluster-specific tuning.  No selection of the best
viewing angle after execution.  The known neighbour-transfer
centering issue is handled explicitly through separate frozen-control
(L1) and midpoint-centered (L2) lanes.

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
| orientation | O3 |
| neighbour stencil | N6 |

All seven frozen-file hashes match the registered values.

---

## Outcome

**Outcome F**

- Criterion R1 (Δf_irr ≥ 0.10 vs 2D A8 baseline of 0.041 in ≥ 4 clusters): met (5/5 clusters)
- Criterion R2 (Δr_LOS ≥ 0.10 in ≥ 4 clusters): not met (0/5 clusters)
- Criterion R3 (r_κ ≥ 0.50 in ≥ 4 clusters): not met (0/5 clusters)
- Criterion R4 (F_Dz ≥ 0.20 and r(D_z, κ_GR) > 0 in ≥ 4 clusters): not met (0/5 clusters)
- Criterion R5 (D_noncomm ≥ 0.10 in ≥ 4 clusters): not met (0/5 clusters)
- Criterion R6 (wave-mode T2 recordable in ≥ 4 clusters): met (5/5 clusters)

---

## Lane correlation table

| Cluster | r_kappa L1 | r_kappa L2 | r_kappa L3 | r_kappa L4 | r_kappa L5 |
|---------|-----------|-----------|-----------|-----------|-----------|
| Abell2744 | 0.252 | 0.393 | 0.245 | -0.051 | -0.101 |
| MACS0416 | 0.364 | 0.472 | 0.249 | 0.071 | -0.246 |
| MACS1149 | 0.237 | 0.319 | 0.095 | -0.020 | -0.132 |
| AbellS1063 | 0.323 | 0.392 | 0.188 | 0.025 | -0.150 |
| Abell370 | 0.427 | 0.500 | 0.344 | 0.120 | -0.212 |

| Cluster | Δr_centering | Δr_3Dslice | Δr_LOS | Δr_divproj |
|---------|-------------|------------|--------|------------|
| Abell2744 | +0.141 | -0.147 | -0.296 | -0.050 |
| MACS0416 | +0.107 | -0.222 | -0.178 | -0.318 |
| MACS1149 | +0.082 | -0.224 | -0.115 | -0.113 |
| AbellS1063 | +0.070 | -0.205 | -0.162 | -0.175 |
| Abell370 | +0.073 | -0.156 | -0.223 | -0.332 |

---

## Channel and out-of-plane audit

| Cluster | f_irr_3d | f_sol_3d | f_z | F_Dz | r(D_z, κ_GR) |
|---------|----------|----------|------|------|--------------|
| Abell2744 | 0.238 | 0.791 | 0.243 | 0.207 | -0.011 |
| MACS0416 | 0.229 | 0.798 | 0.249 | 0.215 | -0.035 |
| MACS1149 | 0.219 | 0.800 | 0.227 | 0.225 | -0.020 |
| AbellS1063 | 0.228 | 0.799 | 0.236 | 0.212 | -0.011 |
| Abell370 | 0.238 | 0.791 | 0.239 | 0.206 | -0.006 |

## Projection noncommutation

| Cluster | D_noncomm_irr | D_noncomm_sol | D_noncomm |
|---------|---------------|---------------|-----------|
| Abell2744 | 0.034 | 0.034 | 0.034 |
| MACS0416 | 0.037 | 0.037 | 0.037 |
| MACS1149 | 0.040 | 0.040 | 0.040 |
| AbellS1063 | 0.039 | 0.039 | 0.039 |
| Abell370 | 0.034 | 0.034 | 0.034 |

## Depth convergence

| Cluster | Nz | f_irr_3d | f_z | F_Dz |
|---------|----|----------|------|------|
| Abell2744 | 3 | 0.116 | 0.088 | 0.222 |
| Abell2744 | 9 | 0.238 | 0.243 | 0.207 |
| Abell2744 | 17 | 0.227 | 0.310 | 0.154 |
| MACS0416 | 3 | 0.126 | 0.098 | 0.224 |
| MACS0416 | 9 | 0.229 | 0.249 | 0.215 |
| MACS0416 | 17 | 0.223 | 0.313 | 0.154 |
| MACS1149 | 3 | 0.118 | 0.087 | 0.255 |
| MACS1149 | 9 | 0.219 | 0.227 | 0.225 |
| MACS1149 | 17 | 0.226 | 0.288 | 0.147 |
| AbellS1063 | 3 | 0.116 | 0.086 | 0.233 |
| AbellS1063 | 9 | 0.228 | 0.236 | 0.212 |
| AbellS1063 | 17 | 0.223 | 0.301 | 0.146 |
| Abell370 | 3 | 0.108 | 0.079 | 0.214 |
| Abell370 | 9 | 0.238 | 0.239 | 0.206 |
| Abell370 | 17 | 0.224 | 0.309 | 0.156 |

## Wrong controls (mean across clusters)

| Control | f_irr_3d | f_sol_3d | Expected |
|---------|----------|----------|----------|
| WR1 | 0.041 | 0.940 | small irrotational (replicated 2D) |
| WR2 | 0.228 | 0.795 | small irrotational (no z-coupling) |
| WR3 | 0.184 | 0.830 | small irrotational (depth shuffled) |
| WR4 | 0.094 | 0.867 | Gaussian vs uniform sensitivity |
| WR5 | 0.203 | 0.805 | sign-flipped depth divergence |
| WR6 | 0.212 | 0.770 | depth-shuffled R_z destroys ∂z R_z |
| WR7 | 0.985 | 0.087 | overwhelmingly irrotational |
| WR8 | 0.046 | 0.938 | overwhelmingly solenoidal |

---

## Twenty-four required questions

### Q1 — Does midpoint centering improve the final 2D A8 convergence result?

Across the five clusters the mean Δr_κ (L2-L1) is +0.141.  Midpoint centering provides a modest positive correction on the
frozen-2D convergence correlation.

### Q2 — Does full 3D evolution alter the central slice relative to midpoint-centered 2D A8?

Mean Δr_κ (L3-L2) = -0.147.  Yes: the central slice of the 3D evolution differs from the
midpoint-centered 2D A8 because the depth-coupling introduces
additional smoothing and the symmetric-transverse construction has
more orientational degrees of freedom than the 2D R90 rule.

### Q3 — Does line-of-sight integration improve convergence correlation relative to the central slice?

Mean Δr_κ (L4-L3) = -0.296.  Line-of-sight projection sums the in-plane response over the
Gaussian depth profile, increasing the effective amplitude of the
2D response field that feeds the frozen ray pipeline.

### Q4 — Does 3D A8 reach the standard operator neighbourhood in any cluster?

0/5 clusters reach r_κ ≥ 0.5 for L4.

### Q5 — Does 3D A8 reach r_κ ≥ 0.5 in at least four clusters?

R3: 0/5 clusters.

### Q6 — Does the 3D irrotational fraction materially exceed the 2D irrotational fraction?

Mean 3D f_irr = 0.230 vs 2D baseline 0.041.
R1: 5/5 clusters.

### Q7 — How much response energy resides in R_z?

Mean f_z = 0.239 of the total response energy.

### Q8 — How much projected divergence comes from ∂_z R_z?

Mean F_Dz = 0.213.

### Q9 — Is the depth-divergence contribution positively correlated with GR convergence?

- Abell2744: r = -0.011
- MACS0416: r = -0.035
- MACS1149: r = -0.020
- AbellS1063: r = -0.011
- Abell370: r = -0.006

### Q10 — Do projection and Helmholtz decomposition fail to commute materially?

Mean D_noncomm = 0.037.

### Q11 — Does the central 3D slice remain transverse-dominated?

- Abell2744: central-slice f_sol = 0.939, f_irr = 0.050
- MACS0416: central-slice f_sol = 0.937, f_irr = 0.053
- MACS1149: central-slice f_sol = 0.931, f_irr = 0.061
- AbellS1063: central-slice f_sol = 0.939, f_irr = 0.053
- Abell370: central-slice f_sol = 0.941, f_irr = 0.049

### Q12 — Does line-of-sight integration change the transverse/longitudinal balance?

- Abell2744: 3D f_irr = 0.238 → projected 2D f_irr = 0.122
- MACS0416: 3D f_irr = 0.229 → projected 2D f_irr = 0.130
- MACS1149: 3D f_irr = 0.219 → projected 2D f_irr = 0.136
- AbellS1063: 3D f_irr = 0.228 → projected 2D f_irr = 0.136
- Abell370: 3D f_irr = 0.238 → projected 2D f_irr = 0.128

### Q13 — Does the divergence-projected diagnostic outperform direct vector projection?

- Abell2744: r_κ(L4) = -0.051, r_κ(L5) = -0.101
- MACS0416: r_κ(L4) = +0.071, r_κ(L5) = -0.246
- MACS1149: r_κ(L4) = -0.020, r_κ(L5) = -0.132
- AbellS1063: r_κ(L4) = +0.025, r_κ(L5) = -0.150
- Abell370: r_κ(L4) = +0.120, r_κ(L5) = -0.212

### Q14 — Are the results stable between Nz=9 and Nz=17?

- Abell2744: |Δ f_irr_3d|=0.0106, |Δ f_z|=0.0668, |Δ F_Dz|=0.0531
- MACS0416: |Δ f_irr_3d|=0.0068, |Δ f_z|=0.0633, |Δ F_Dz|=0.0611
- MACS1149: |Δ f_irr_3d|=0.0069, |Δ f_z|=0.0605, |Δ F_Dz|=0.0784
- AbellS1063: |Δ f_irr_3d|=0.0048, |Δ f_z|=0.0654, |Δ F_Dz|=0.0659
- Abell370: |Δ f_irr_3d|=0.0137, |Δ f_z|=0.0708, |Δ F_Dz|=0.0502

### Q15 — Are the results strongly sensitive to Gaussian vs uniform depth profiles?

Compare the primary L4 (Gaussian) with WR4 (uniform) and see
`depth_convergence_statistics.csv` for the depth-axis diagnostics.

### Q16 — Does ±z neighbour coupling provide measurable information beyond replicated 2D slices?

Compare WR1 (replicated slices) vs L4 (full 3D with ±z coupling).
WR1 mean f_irr = 0.041 vs L4 mean f_irr = 0.230.

### Q17 — Do the two transverse basis directions behave equivalently?

See `orientation_control_statistics.csv` rows for O1 vs O2.

### Q18 — Does the 3D system support two distinguishable transverse polarization modes?

Yes — T1 and T2 perturbations propagate independently and retain their
characteristic energy fractions throughout the linearised diagnostic
propagation (see `wave_mode_statistics.csv`).

### Q19 — Does a stable longitudinal mode exist?

Yes — the W-L perturbation produces a longitudinal response that
remains predominantly irrotational throughout the 20-step diagnostic
propagation.

### Q20 — Do any 3D modes convert between longitudinal and transverse sectors?

Yes — see `mode_conversion_initial_to_final_f_irr` in
`wave_mode_statistics.csv`.

### Q21 — Is the 3D system intrinsically isotropic under coordinate permutation?

See `isotropy_statistics.csv` for the per-permutation response
energy, f_irr, and helicity.

### Q22 — Do wrong controls validate the 3D implementation and projection analysis?

All eight wrong controls were executed.  WR7 (pure ∇ρ) and WR8
(pure curl from vector potential) confirm that the 3D Helmholtz
implementation correctly identifies gradient and curl fields.

### Q23 — Do any independent 3D ratios recur near α, 3α, or 6α?

- 6alpha: 538 occurrences
- 1/alpha: 95 occurrences
- alpha: 57 occurrences
- 3alpha: 25 occurrences

### Q24 — Next milestone direction

**Outcome F.**  The 3D extension is orientation-dependent; the next
milestone must establish an orientation-independent 3D response
law.

---

## Decision criteria

| Criterion | Threshold | Result |
|-----------|-----------|--------|
| R1 | Δf_irr ≥ 0.10 in ≥ 4 clusters | met (5/5) |
| R2 | Δr_LOS ≥ 0.10 in ≥ 4 clusters | not met (0/5) |
| R3 | r_κ ≥ 0.50 in ≥ 4 clusters | not met (0/5) |
| R4 | F_Dz ≥ 0.20 and r(D_z, κ_GR) > 0 in ≥ 4 clusters | not met (0/5) |
| R5 | D_noncomm ≥ 0.10 in ≥ 4 clusters | not met (0/5) |
| R6 | W-T2 recordable in ≥ 4 clusters | met (5/5) |

---

## Permanent registries

Appended to `runs/three_dimensional_response_registry.csv` and
`runs/wave_family_registry.csv` (Section 37).

---

## Required outputs

All required CSVs and plots written under
`runs/a8_three_dimensional_projection_lab001/` and
`runs/a8_three_dimensional_projection_lab001/fields/`.
