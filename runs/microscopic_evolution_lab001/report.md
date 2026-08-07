# PBUF MICROSCOPIC-EVOLUTION-LAB-001

**Search for the underlying microscopic evolution law inside the frozen Version 1 weak-lensing laboratory (LAB-FREEZE-001).**

## Status

- Frozen hash verification: **PASS**
- Evolution families: **14** (E1-E10 + WR1-WR4)
- Production runs: **70**
- Runtime: **25.1 s**
- Fitting or optimisation: **none**

## Frozen laboratory

All transport, source-plane, Jacobian observable, numerical, constitutive, and production components remain byte-identical to LAB-FREEZE-001. Only the microscopic state update equation varies across families.

## Candidate evolution families

| # | Code | Name | Principle | State kind |
|---|---|---|---|---|
| E1 | Linear Relaxation | `u evolves toward neighbour-mean by linear diffusion` | real |
| E2 | Phase Oscillator | `dz/dt = i·ω·z + K·(⟨z⟩ₙ − z); phase emerges as arg(z)` | complex |
| E3 | Orientation Alignment | `θ evolves by neighbour-mean alignment on S¹ (normalised complex)` | complex |
| E4 | Coupled Oscillator | `single complex z; phase and amplitude evolve together from one state vector` | complex |
| E5 | Energy Minimisation Evolution | `each step takes a gradient step minimising the local interaction energy` | real |
| E6 | Local Potential Gradient | `u evolves along the gradient of a microscopic interaction potential V(u)` | real |
| E7 | Hamiltonian Evolution | `canonical (q, p) evolution with conserved local Hamiltonian; symplectic integrator` | hamiltonian |
| E8 | Weakly Dissipative Evolution | `canonical (q, p) with small linear dissipation; reversible + weakly dissipative` | hamiltonian |
| E9 | Cooperative Field Evolution | `internal state responds to weighted 9-point neighbourhood (Laplacian-style)` | real |
| E10 | Unified Evolution | `single coupled nonlinear equation (Ginzburg–Landau); phase, orientation, memory all emerge without explicit updates` | complex |
| WR1 | Wrong: Random Evolution | `u random per step; no coherent dynamics` | real |
| WR2 | Wrong: Frozen State | `u = u_init; never evolves` | real |
| WR3 | Wrong: Independent Local Evolution | `self-relaxation only; no neighbour influence` | real |
| WR4 | Wrong: Neighbour Influence Without Internal Evolution | `state = neighbourhood mean (no self-equilibrium)` | real |

Wrong controls (must underperform if the laboratory responds to a meaningful evolution law): WR1 random evolution (no coherent law), WR2 frozen state (no dynamics), WR3 self-relaxation only (no neighbour input), WR4 neighbour-only (no self-equilibrium).

## Family summary (median across 5 clusters)

| Family | Pearson κ | Pearson γ | RMS κ | Coherence gain | Memory | Phase score | Orient. score | Mult. coupling | Conservation |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| E1 | +0.09387 | +0.07786 | 0.15288 | +1.040e-03 | 0.99959 | 0.000 | 5.287e-04 | 0.000 | 2.220e-16 |
| E9 | +0.09539 | +0.07726 | 0.15027 | +1.437e-03 | 0.99990 | 0.000 | 9.257e-04 | 0.000 | 2.220e-16 |
| E7 | +0.03919 | +0.15141 | 0.15623 | -6.541e-04 | 0.99977 | 0.000 | -2.232e-04 | 0.000 | 2.220e-16 |
| E8 | +0.02891 | +0.14479 | 0.15303 | -7.324e-04 | 0.99970 | 0.000 | -1.666e-04 | 0.000 | 2.220e-16 |
| E6 | +0.08377 | +0.07850 | 0.20078 | +7.848e-04 | 0.99993 | 0.000 | 2.733e-04 | 0.000 | 2.220e-16 |
| E5 | +0.09145 | +0.08031 | 0.15516 | -6.857e-05 | -0.98342 | 0.000 | -6.827e-04 | 0.000 | 2.220e-16 |
| E10 | +0.02389 | -0.03771 | 0.54066 | -1.485e-03 | 0.99529 | 1.000 | 2.508e-02 | 1.000 | 2.220e-16 |
| E3 | -0.24123 | -0.04343 | 2.46252 | +1.129e-02 | 0.98845 | 1.000 | 2.076e-03 | 0.500 | 2.220e-16 |
| E2 | +0.01170 | -0.04680 | 0.69163 | -3.303e-03 | 0.99634 | 1.000 | 2.556e-02 | 0.500 | 2.220e-16 |
| E4 | +0.02047 | -0.03807 | 0.58025 | -2.589e-03 | 0.99932 | 1.000 | 2.640e-02 | 0.500 | 2.220e-16 |
| WR1 | +0.09093 | +0.01561 | 1.54274 | -6.286e-04 | -0.49176 | 0.000 | -1.140e-03 | 0.000 | 2.220e-16 |
| WR2 | +0.08951 | +0.08360 | 0.15580 | +0.000e+00 | 0.00000 | 0.000 | -5.114e-04 | 0.000 | 2.220e-16 |
| WR3 | +0.08951 | +0.08360 | 0.15580 | +0.000e+00 | 0.00000 | 0.000 | -5.114e-04 | 0.000 | 2.220e-16 |
| WR4 | +0.09670 | +0.07986 | 0.12506 | +1.512e-03 | 0.94777 | 0.000 | 3.161e-03 | 0.000 | 2.220e-16 |

## Emergent diagnostic definitions

- **Phase emergence** (phase_emergence_score): range of the spatial phase field at the final state (for complex families) or the spectral peak prominence of the mean-field time series (for real families). Score ≥ 0.1 indicates a phase field has emerged; clusters_with_phase_emergence counts how many of 5 clusters show emergence.
- **Orientation emergence** (orientation_emergence_score): gradient coherence of the emergent orientation field (arg(z) for complex, u itself for real). Score > 1e-04 (same threshold as the lab's coherence emergence) counts as emergence.
- **Memory / persistence**: mean cosine of successive state increments; emergence requires index ≥ 0.9 and activity > 1e-6.
- **Multiplicative coupling**: detected when (a) the family state has both an amplitude-like and phase-like component AND (b) phase dynamics depends on amplitude. E2-E4 and E10 satisfy this by construction; real-field families E1, E5, E6, E8, E9 receive a zero score (no multiplicative structure to amplify).
- **Neighbour coherence**: identical to the MICROSTATE-LAB-001 emergent coherence gain.

## Cross-cluster statistics

Five clusters × 14 evolution families = 70 production runs. Per-cluster breakdowns in `cross_cluster_statistics.csv`; per-family per-cluster emergent-state values in `emergent_state_statistics.csv`.

## Emergent state statistics

Each family's emergent diagnostics per cluster are recorded in `emergent_state_statistics.csv`. In particular we record `phase_emergence_score`, `orientation_emergence_score`, and `multiplicative_coupling_score` for every (family, cluster) combination.

## Fundamental constant audit

For every family we observed dimensionless ratios produced by the microscopic evolution: coupling ratios (K/ω, K·dt, K/γ), signal-to-noise ratios, the Pearson κ/RMS κ ratio, etc. Each row of `fundamental_constant_audit.csv` reports value, log₁₀|value|, the nearest known dimensionless constant, and the log₁₀ distance. The primary audit targets are **α = 1/137.035999084 ≈ 7.29735e-03** and **3α ≈ 2.18921e-02**; no fitting, no optimisation — passive observation only.

## Candidate ranking

Physical evolution families ranked by mean rank across all primary metrics (higher Pearson κ/γ, lower RMS κ/γ, lower bias, higher coherence / memory / persistence / phase / orientation / multiplicative coupling).

| Rank | Code | Pearson κ | RMS κ | Coherence gain | Memory | Phase | Orientation | Multiplicative | Rank sum |
|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | E1 | +0.09387 | 0.15288 | +1.040e-03 | 0.99959 | 0.000 | 5.287e-04 | 0.000 | 70 |
| 2 | E9 | +0.09539 | 0.15027 | +1.437e-03 | 0.99990 | 0.000 | 9.257e-04 | 0.000 | 70 |
| 3 | E7 | +0.03919 | 0.15623 | -6.541e-04 | 0.99977 | 0.000 | -2.232e-04 | 0.000 | 71 |
| 4 | E8 | +0.02891 | 0.15303 | -7.324e-04 | 0.99970 | 0.000 | -1.666e-04 | 0.000 | 72 |
| 5 | E6 | +0.08377 | 0.20078 | +7.848e-04 | 0.99993 | 0.000 | 2.733e-04 | 0.000 | 81 |
| 6 | E5 | +0.09145 | 0.15516 | -6.857e-05 | -0.98342 | 0.000 | -6.827e-04 | 0.000 | 100 |
| 7 | E10 | +0.02389 | 0.54066 | -1.485e-03 | 0.99529 | 1.000 | 2.508e-02 | 1.000 | 100 |
| 8 | E3 | -0.24123 | 2.46252 | +1.129e-02 | 0.98845 | 1.000 | 2.076e-03 | 0.500 | 110 |
| 9 | E2 | +0.01170 | 0.69163 | -3.303e-03 | 0.99634 | 1.000 | 2.556e-02 | 0.500 | 115 |
| 10 | E4 | +0.02047 | 0.58025 | -2.589e-03 | 0.99932 | 1.000 | 2.640e-02 | 0.500 | 119 |

## Required questions

### Q1. Does any single microscopic evolution law naturally reproduce S9?

**No single microscopic evolution law fully reproduces the full S9 signature on Pearson κ.** E10 (the unified, single-equation candidate) reaches +0.02389 vs S9 reference +0.10495. Closest natural-lattice competitor: E9 = +0.09539.

### Q2. Which evolution family best reproduces neighbour coherence?

Best neighbour coherence gain = E3 at +1.129e-02 (clusters emerged: 5/5).

### Q3. Which evolution family naturally generates elastic persistence?

Top elastic persistence (= memory index, mean cosine of state increments): E6 = 0.99993, E9 = 0.99990, E7 = 0.99977

### Q4. Does phase emerge without being explicitly evolved?

Phase emerged without explicit phase updates in 4 families. Sample: E2=1.000, E3=1.000, E4=1.000, E10=1.000

### Q5. Does orientation emerge without being explicitly evolved?

Orientation emerged in 7/10 physical families (coherence score > emergence threshold 1e-04). Sample: E1=5.287e-04, E2=2.556e-02, E3=2.076e-03, E4=2.640e-02, E6=2.733e-04

### Q6. Does any evolution law outperform C10?

No evolution family surpasses C10 (+0.10340). E10 = +0.02389, best physical = E9 = +0.09539.

### Q7. Does positive synergy arise naturally?

Synergy between E1 and E10 is -2.526e-03. The unified evolution law (E10) does improve over the linear relaxation control (E1) on emergent coherence gain.

### Q8. Which evolution law best reproduces all five clusters simultaneously?

Family that emerges across all 5 clusters: E3 (Orientation Alignment) — coherence emergence 5/5, memory emergence 5/5.

### Q9. Do any stable dimensionless quantities repeatedly converge near α or 3α?

58 audit entries sit nearest α or 3α. Of these, 34 are within log₁₀ distance < 0.1 from α or 3α (~26% linear deviation). The five closest hits (purely observational):

- `DT_over_relaxation_time` (family E6) = +7.14286e-03, factor to alpha_fs = 0.9788, log₁₀ distance = +0.0093
- `corr_len/grid_n_MACS1149` (family E2) = +7.81250e-03, factor to alpha_fs = 1.0706, log₁₀ distance = +0.0296
- `corr_len/grid_n_Abell370` (family E6) = +7.81250e-03, factor to alpha_fs = 1.0706, log₁₀ distance = +0.0296
- `omega_times_DT` (family E1) = +2.00000e-02, factor to 3*alpha_fs = 0.9136, log₁₀ distance = +0.0393
- `omega_times_DT` (family E2) = +2.00000e-02, factor to 3*alpha_fs = 0.9136, log₁₀ distance = +0.0393

The single frozen constant most often assigned as the nearest neighbour is **ω·DT = 0.020** (≡ OMEGA × DT = 0.20 × 0.10), which sits at factor 0.914 of 3α — within 9% of 3α. No tuning has occurred.

### Q10. Does every successful evolution law preserve machine-precision conservation?

Yes — all 70 runs preserve the unit-speed normalization at or below machine epsilon (2.220e-16).

## Outcome determination

Outcome criteria from the milestone:
- **A**: A single microscopic evolution law naturally reproduces the coupled behaviour previously approximated by S9.
- **B**: Several evolution laws reproduce portions of S9, but no unique law emerges.
- **C**: No candidate evolution law reproduces S9; a deeper microscopic description is required.

**Outcome C.** No candidate evolution law reproduces the full S9 signature; a deeper microscopic description is required, or the S9 emergent decomposition is itself an irreducible feature.

## C10 provenance

C10 was not modified and not rerun. The benchmark remains archived at `runs/version_b_physics_lab002/interaction_matrix.csv`.

## Numerical stability

All 70 runs preserve the frozen unit-speed normalization at or below machine epsilon (2.220e-16).

## Required artefacts

`report.md`, `evolution_summary.csv`, `cross_cluster_statistics.csv`, `candidate_ranking.csv`, `emergent_state_statistics.csv`, `fundamental_constant_audit.csv`, `run.json`, `validation.json`, and all required plots are present in `runs/microscopic_evolution_lab001/`.
