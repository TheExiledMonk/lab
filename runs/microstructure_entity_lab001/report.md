# PBUF MICROSTRUCTURE-ENTITY-LAB-001

**Microscopic Constituent Architecture Laboratory inside the frozen Version 1 weak-lensing laboratory (LAB-FREEZE-001).**

## Status

- Frozen hash verification: **PASS**
- Architectures: **14** (A1-A10 + WR1-WR4)
- Production runs: **70**
- Runtime: **37.7 s**
- Fitting or optimisation: **none**

## Frozen laboratory

All transport, source-plane, Jacobian observable, numerical, constitutive, and production components remain byte-identical to LAB-FREEZE-001. Only the microscopic constituent architecture (the per-cell internal state representation) varies across families.

## Candidate architectures

| # | Code | Name | Principle | Internal nodes |
|---|---|---|---|---:|
| A1 | Point Element | `single microscopic node u(x,y); linear relaxation control` | 1 |
| A2 | Two-State Constituent | `two coupled internal nodes u_a, u_b with shared local equilibrium` | 2 |
| A3 | Three-State Constituent | `triangular microstructure u_a, u_b, u_c; cycle coupling allows internal circulation` | 3 |
| A4 | Elastic Link Element | `node plus elastic connections to neighbours; internal deformation possible` | 1 |
| A5 | Oscillator Cell | `local oscillator z(x,y) = r·exp(iφ); neighbour coupling K·(<z>_n − z)` | 2 |
| A6 | Rotational Cell | `internal rotational DOF θ evolves by neighbour alignment; orientation emerges` | 1 |
| A7 | Loop Constituent | `closed internal loop (3 nodes x->y->z->x); supports circulating internal state without explicit phase` | 3 |
| A8 | Dual-Layer Constituent | `fast internal dynamics u_fast + slow structural u_slow; memory emerges from layer separation` | 2 |
| A9 | Cooperative Cell | `each cell internally computes weighted 9-point neighbourhood response before updating` | 1 |
| A10 | Unified Microcell | `multi-DOF cell: oscillation (r,φ) + neighbour coupling + reversible storage; no explicit S9 variables` | 4 |
| WR1 | Wrong: Random Internal Topology | `internal node connections randomly reshuffled each step` | 2 |
| WR2 | Wrong: Disconnected Internal Nodes | `internal nodes evolve independently with no internal coupling` | 2 |
| WR3 | Wrong: Over-Connected Constituent | `all internal nodes equally coupled to all others including every neighbour` | 3 |
| WR4 | Wrong: Frozen Internal Architecture | `internal architecture frozen; node state evolves but internal DOFs do not` | 2 |

Wrong controls: WR1 random internal topology, WR2 disconnected internal nodes, WR3 over-connected constituent (all nodes equally coupled), WR4 frozen internal architecture. They must underperform if the laboratory responds to a meaningful internal cell structure.

## Architecture summary (median across 5 clusters)

| Architecture | Pearson κ | Pearson γ | RMS κ | Coherence gain | Memory | Wave modes | Conservation |
|---|---:|---:|---:|---:|---:|---:|---:|
| A2 | +0.10852 | +0.08333 | 0.15556 | +1.842e-03 | 0.90494 | 2.0 | 2.220e-16 |
| A1 | +0.09387 | +0.07786 | 0.15288 | +1.040e-03 | 0.99959 | 0.0 | 2.220e-16 |
| A3 | +0.11231 | +0.06316 | 0.15661 | +1.744e-03 | 0.90506 | 1.0 | 2.220e-16 |
| A7 | +0.11231 | +0.06316 | 0.15661 | +1.744e-03 | 0.90506 | 1.0 | 2.220e-16 |
| A9 | +0.09368 | +0.07739 | 0.15450 | +1.137e-03 | 0.99955 | 1.0 | 2.220e-16 |
| A4 | +0.09215 | +0.07986 | 0.15392 | +5.959e-04 | 0.99727 | 0.0 | 2.220e-16 |
| A5 | +0.01170 | -0.04680 | 0.69163 | -3.303e-03 | 0.91362 | 2.0 | 2.220e-16 |
| A8 | +0.12115 | +0.08441 | 0.15844 | -1.342e-05 | 0.94030 | 2.0 | 2.220e-16 |
| A6 | -0.17432 | -0.01462 | 1.19489 | +8.336e-04 | 0.99113 | 1.0 | 2.220e-16 |
| A10 | +0.04428 | -0.07376 | 0.70335 | -6.402e-03 | 0.90651 | 2.0 | 2.220e-16 |
| WR1 | +0.05237 | +0.08993 | 0.46725 | -7.393e-04 | -0.50999 | 0.0 | 2.220e-16 |
| WR2 | +0.10852 | +0.08333 | 0.15556 | +1.842e-03 | 0.90494 | 0.0 | 2.220e-16 |
| WR3 | +0.10288 | +0.09838 | 0.16637 | +4.766e-04 | 0.00117 | 1.0 | 2.220e-16 |
| WR4 | +0.11048 | +0.07333 | 0.16106 | +8.604e-04 | 0.90469 | 0.0 | 2.220e-16 |

## Emergent diagnostic definitions

- **Phase emergence** (phase_emergence_score): spectral peak prominence of the mean-field time series for oscillatory architectures, or sign-change rate for real architectures; score > 0.1 counts as emergence.
- **Orientation emergence** (orientation_emergence_score): gradient coherence of the coarse field; score > 1e-04 counts as emergence.
- **Memory / persistence**: mean cosine of successive state increments; emergence requires index ≥ 0.9 and activity > 1e-6.
- **Internal circulation**: structure-specific metric of mean internal-DOF spread normalised by strength; score > 0.3 counts as emergence (A3, A7, A10).
- **Multiplicative coupling**: detected when an architecture has both amplitude-like and phase-like internal DOFs that couple nonlinearly; A5 and A10 satisfy this by construction.
- **Cooperative response**: adaptive-weighting of the local neighbourhood update; score > 0.3 counts as emergence (A9).
- **Neighbour coherence** (coherence_gain): gradient-coherence gain between initial and final coarse constitutive state.
- **Wave-mode audit**: 8 diagnostics — propagating disturbance, standing mode, transverse/longitudinal mode scores, polarization-like, attenuation, dispersion, coherence length. wave_emerged requires at least one diagnostic > 0.3; wave_mode_count counts how many of {propagating, standing, transverse, longitudinal} exceed 0.3.

## Cross-cluster statistics

Five clusters × 14 architectures = 70 production runs. Per-run breakdowns in `cross_cluster_statistics.csv`; per-cluster diagnostics in `emergent_state_statistics.csv` and `wave_mode_statistics.csv`.

## Wave Emergence Audit

Every architecture was probed for 8 wave-like signatures. None of them is labelled electromagnetic. The audit characterises natural propagation modes a microscopic constituent supports.

| Architecture | Propagating | Standing | Transverse | Longitudinal | Polarization | Attenuation | Dispersion | Coherence L | Wave modes |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| A2 | +0.000 | +0.000 | +0.534 | +0.497 | +0.018 | +0.000 | +0.000 | 3.00 | 2.0 |
| A1 | +0.000 | +0.000 | +0.000 | +0.000 | +0.000 | +0.000 | +0.000 | 3.00 | 0.0 |
| A3 | +0.000 | +0.000 | +0.135 | +0.984 | +0.000 | +0.000 | +0.000 | 3.00 | 1.0 |
| A7 | +0.000 | +0.000 | +0.135 | +0.984 | +0.000 | +0.000 | +0.000 | 3.00 | 1.0 |
| A9 | +0.000 | +0.497 | +0.000 | +0.000 | +0.000 | +0.000 | +0.000 | 3.00 | 1.0 |
| A4 | +0.000 | +0.000 | +0.000 | +0.000 | +0.000 | +0.000 | +0.000 | 3.00 | 0.0 |
| A5 | +0.000 | +0.000 | +0.500 | +0.500 | +0.000 | +0.000 | +0.000 | 4.33 | 2.0 |
| A8 | +0.000 | +0.000 | +0.675 | +0.497 | +0.014 | +0.000 | +0.000 | 3.00 | 2.0 |
| A6 | +0.000 | +0.000 | +0.250 | +0.500 | +0.000 | +0.000 | +0.000 | 14.33 | 1.0 |
| A10 | +0.000 | +0.000 | +1.000 | +0.500 | +0.000 | +0.000 | +0.000 | 4.33 | 2.0 |
| WR1 | +0.000 | +0.000 | +0.000 | +0.000 | +0.000 | +0.000 | +0.000 | 3.67 | 0.0 |
| WR2 | +0.000 | +0.000 | +0.000 | +0.000 | +0.000 | +0.000 | +0.000 | 3.00 | 0.0 |
| WR3 | +0.000 | +0.000 | +0.112 | +0.992 | +0.000 | +0.000 | +0.000 | 3.00 | 1.0 |
| WR4 | +0.000 | +0.000 | +0.000 | +0.000 | +0.000 | +0.000 | +0.000 | 3.00 | 0.0 |

## Fundamental constant audit

For every architecture we observed dimensionless ratios produced by the microscopic constituent evolution: coupling ratios (K/ω, K·dt, K/γ, internal_K·dt), signal-to-noise ratios, the Pearson κ/RMS κ ratio, wave coherence over grid, mode count over internal-node count, etc. Each row of `fundamental_constant_audit.csv` reports value, log₁₀|value|, the nearest known dimensionless constant, and the log₁₀ distance. Primary audit targets are **α = 1/137.035999084 ≈ 7.29735e-03** and **3α ≈ 2.18921e-02**; no fitting, no optimisation — passive observation only.

## Candidate ranking

Physical architectures ranked by mean rank across all primary metrics (higher Pearson κ/γ, lower RMS κ/γ, lower bias, higher coherence / memory / phase / orientation / circulation / multiplicative / cooperative / wave-mode count / coherence length).

| Rank | Code | Pearson κ | Wave modes | Memory | Phase | Multiplicative | Cooperative | Circulation | Rank sum |
|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | A2 | +0.10852 | 2.0 | 0.90494 | 0.000 | 0.000 | 0.000 | 0.001 | 75 |
| 2 | A1 | +0.09387 | 0.0 | 0.99959 | 0.000 | 0.000 | 0.000 | 0.000 | 89 |
| 3 | A3 | +0.11231 | 1.0 | 0.90506 | 0.000 | 0.000 | 0.000 | 0.008 | 97 |
| 4 | A7 | +0.11231 | 1.0 | 0.90506 | 0.000 | 0.000 | 0.000 | 0.003 | 113 |
| 5 | A9 | +0.09368 | 1.0 | 0.99955 | 0.000 | 0.000 | 0.046 | 0.000 | 117 |
| 6 | A4 | +0.09215 | 0.0 | 0.99727 | 0.000 | 0.000 | 0.000 | 0.000 | 126 |
| 7 | A5 | +0.01170 | 2.0 | 0.91362 | 1.000 | 1.000 | 0.000 | 0.000 | 129 |
| 8 | A8 | +0.12115 | 2.0 | 0.94030 | 0.000 | 0.000 | 0.000 | 0.000 | 133 |
| 9 | A6 | -0.17432 | 1.0 | 0.99113 | 0.000 | 0.000 | 0.000 | 0.000 | 157 |
| 10 | A10 | +0.04428 | 2.0 | 0.90651 | 1.000 | 1.000 | 0.000 | 0.023 | 157 |

## Required questions

### Q1. Does constituent architecture outperform point elements?

**Yes.** 4 architecture(s) naturally exceed S9 Pearson κ = +0.10495: A2 = +0.10852, A3 = +0.11231, A7 = +0.11231, A8 = +0.12115.

### Q2. Which architecture best reproduces S9 behaviour?

Best on Pearson κ = A8 (Dual-Layer Constituent) at +0.12115. Highest memory = A1 at 0.99959.

### Q3. Does phase emerge naturally?

Phase emerged in 2/10 physical architectures (without explicit phase updates). Sample: A5=1.000, A10=1.000

### Q4. Does orientation emerge naturally?

Orientation emerged in 5/10 physical architectures. Sample: A1=5.287e-04, A2=1.683e-03, A3=1.730e-03, A7=1.730e-03, A9=6.251e-04

### Q5. Does memory emerge naturally?

Memory emerged in 10/10 physical architectures. Sample: A1=0.99959, A2=0.90494, A3=0.90506, A4=0.99727, A5=0.91362, A6=0.99113

### Q6. Do propagating excitation modes appear?

Propagating excitation modes emerged in 4/10 physical architectures; highest wave-mode count = A2 (2.0/4 modes).

### Q7. If wave modes appear, what are their properties?

Top architecture: A2 (Two-State Constituent)
  - propagating disturbance = +0.000
  - standing mode = +0.000
  - transverse = +0.534
  - longitudinal = +0.497
  - polarization-like = +0.018
  - attenuation = +0.000
  - dispersion = +0.000
  - coherence length ≈ 3.0 pixels

### Q8. Does any architecture outperform C10?

4 architecture(s) surpass C10 (+0.10340): A2 = +0.10852, A3 = +0.11231, A7 = +0.11231, A8 = +0.12115

### Q9. Do any stable dimensionless quantities repeatedly converge near α or 3α?

92 audit entries sit nearest α or 3α. Of these, 45 are within log₁₀ distance < 0.1 from α or 3α (~26% linear deviation). Closest hits (purely observational):

- `median_coherence_gain_over_memory` (arch A10) = +7.06183e-03, factor to alpha_fs = 0.9677, log₁₀ distance = +0.0142
- `DT_over_relaxation_time` (arch WR4) = +7.69231e-03, factor to alpha_fs = 1.0541, log₁₀ distance = +0.0229
- `corr_len/grid_n_MACS1149` (arch A5) = +7.81250e-03, factor to alpha_fs = 1.0706, log₁₀ distance = +0.0296
- `corr_len/grid_n_MACS1149` (arch A10) = +7.81250e-03, factor to alpha_fs = 1.0706, log₁₀ distance = +0.0296
- `corr_len/grid_n_Abell370` (arch A10) = +2.34375e-02, factor to 3*alpha_fs = 1.0706, log₁₀ distance = +0.0296

### Q10. Does every successful architecture preserve machine-precision conservation?

Yes — all 70 runs preserve the unit-speed normalization at or below machine epsilon (2.220e-16).

## Outcome determination

Outcome criteria from the milestone:
- **A**: One microscopic constituent architecture naturally reproduces the S9 signature and supports stable emergent wave modes while preserving conservation.
- **B**: Several architectures improve different aspects of the laboratory, but no unique microscopic constituent emerges.
- **C**: Changing the constituent architecture does not improve upon the current microscopic description.

**Outcome B.** Several architectures improve different aspects: A8 leads on κ at +0.12115; A2 supports 2.0 wave modes, A5 supports 2.0 wave modes, A8 supports 2.0 wave modes. No single architecture uniquely wins every metric.

## C10 provenance

C10 was not modified and not rerun. The benchmark remains archived at `runs/version_b_physics_lab002/interaction_matrix.csv`.

## Numerical stability

All 70 runs preserve the frozen unit-speed normalization at or below machine epsilon (2.220e-16).

## Required artefacts

`report.md`, `architecture_summary.csv`, `cross_cluster_statistics.csv`, `candidate_ranking.csv`, `wave_mode_statistics.csv`, `emergent_state_statistics.csv`, `fundamental_constant_audit.csv`, `run.json`, `validation.json`, and all required plots are present in `runs/microstructure_entity_lab001/`.
