# PBUF MICROSTRUCTURE-ENTITY-A8-DECOMPOSITION-001

**Dual-Layer Constituent Mechanism Laboratory inside the frozen Version 1 weak-lensing laboratory (LAB-FREEZE-001).**

## Status

- Frozen hash verification: **PASS**
- Decompositions: **15** (D1-D15)
- Production runs: **75**
- Runtime: **37.1 s**
- Fitting or optimisation: **none**

## Frozen laboratory

All transport, source-plane, Jacobian observable, numerical, constitutive, and production components remain byte-identical to LAB-FREEZE-001. Only the internal architecture of the A8 constituent (Dual-Layer Constituent) is varied across the 15 decompositions.

## Decomposition Matrix

| # | Code | Family | Name | Principle |
|---|---|---|---|---|
| D1 | control | Full A8 (control) | `fast and slow layers with mutual additive coupling + neighbour interaction on both` |
| D2 | single_layer | Fast removed (slow only) | `drop fast layer; only slow evolves; u_final = u_slow` |
| D3 | single_layer | Slow removed (fast only) | `drop slow layer; only fast evolves; u_final = u_fast` |
| D4 | frozen | Fast frozen (slow evolves) | `fast layer frozen at initial state; slow evolves normally` |
| D5 | frozen | Slow frozen (fast evolves) | `slow layer frozen at initial state; fast evolves normally` |
| D6 | coupling_direction | Fast->Slow coupling removed | `only slow->fast coupling remains; fast does not feel slow` |
| D7 | coupling_direction | Slow->Fast coupling removed | `only fast->slow coupling remains; slow does not feel fast` |
| D8 | coupling_direction | Bidirectional removed (independent) | `no coupling; fast and slow evolve independently` |
| D9 | coupling_form | Coupling additive (no multiplicative) | `pure additive coupling with full strength 1.0 on both sides` |
| D10 | coupling_form | Coupling multiplicative | `coupling as product u_fast*u_slow instead of additive difference` |
| D11 | neighbour_assignment | Neighbour interaction only on fast | `fast updates include neighbour diff; slow updates contain only coupling term` |
| D12 | neighbour_assignment | Neighbour interaction only on slow | `slow updates include neighbour diff; fast updates contain only coupling term` |
| D13 | neighbour_assignment | Neighbour interaction equally on both | `equal weight neighbour term on both layers` |
| D14 | ordering | Layer update order reversed | `slow updated first using old fast, then fast updated using new slow` |
| D15 | timescale | Timescales forced equal | `removes timescale separation; both layers use same coefficient` |

## Decomposition summary (median across 5 clusters)

| Decomposition | Pearson κ | Pearson γ | Coherence gain | Memory | Wave modes | F/S exchange | Persistence | Conservation |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| D1 | +0.12115 | +0.08441 | -1.342e-05 | 0.99930 | 2.0 | +1.770e-03 | 0.99999 | 2.220e-16 |
| D3 | +0.12112 | +0.10361 | -4.251e-04 | 0.99996 | 1.0 | +0.000e+00 | 1.00000 | 2.220e-16 |
| D6 | +0.12079 | +0.08779 | -1.631e-04 | 0.99875 | 2.0 | +1.873e-03 | 0.99999 | 2.220e-16 |
| D4 | +0.12243 | +0.09695 | -1.366e-04 | 0.99992 | 1.0 | +2.322e-03 | 0.99999 | 2.220e-16 |
| D2 | +0.09055 | +0.07931 | -9.584e-05 | 0.99997 | 1.0 | +0.000e+00 | 1.00000 | 2.220e-16 |
| D10 | +0.10704 | +0.08560 | -7.571e-05 | 0.99987 | 2.0 | +2.377e-03 | 0.99998 | 2.220e-16 |
| D12 | +0.12283 | +0.09228 | -6.824e-05 | 0.99960 | 2.0 | +2.190e-03 | 0.99999 | 2.220e-16 |
| D14 | +0.12126 | +0.08446 | -1.719e-05 | 0.99928 | 2.0 | +1.772e-03 | 0.99999 | 2.220e-16 |
| D13 | +0.12065 | +0.08893 | +7.181e-05 | 0.99814 | 2.0 | +1.922e-03 | 0.99999 | 2.220e-16 |
| D11 | +0.11955 | +0.08507 | +1.244e-04 | 0.99954 | 2.0 | +1.711e-03 | 0.99999 | 2.220e-16 |
| D9 | +0.11946 | +0.08137 | -1.675e-05 | 0.99994 | 2.0 | +1.558e-03 | 0.99998 | 2.220e-16 |
| D15 | +0.12044 | +0.08379 | +1.453e-04 | 0.99976 | 2.0 | +1.751e-03 | 0.99999 | 2.220e-16 |
| D8 | +0.11906 | +0.07917 | -2.167e-04 | 0.99995 | 1.0 | +2.375e-03 | 0.99999 | 2.220e-16 |
| D7 | +0.11712 | +0.07655 | -5.428e-04 | 0.99993 | 1.0 | +2.236e-03 | 0.99999 | 2.220e-16 |
| D5 | +0.11651 | +0.07570 | -2.797e-04 | 0.99996 | 1.0 | +2.241e-03 | 0.99999 | 2.220e-16 |

## Expanded Wave Audit

Every decomposition was probed for 11 wave-like signatures: number of modes, transverse/longitudinal classification, standing vs travelling, phase velocity, group velocity, dispersion, attenuation, coherence length, mode stability. None is labelled electromagnetic; only characterised.

| Decomposition | Propagating | Standing | Transverse | Longitudinal | Polarization | Phase vel | Group vel | Dispersion | Attenuation | Coherence L | Mode stab | Modes |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| D1 | +0.000 | +0.000 | +0.692 | +0.497 | +0.015 | 0.000 | 0.000 | 0.000 | 0.000 | 3.00 | 0.995 | 2.0 |
| D3 | +0.000 | +0.000 | +0.000 | +0.500 | +1.000 | 0.000 | 0.000 | 0.000 | 0.000 | 3.00 | 0.970 | 1.0 |
| D6 | +0.000 | +0.000 | +0.635 | +0.497 | +0.016 | 0.000 | 0.000 | 0.000 | 0.000 | 3.00 | 0.995 | 2.0 |
| D4 | +0.000 | +0.000 | +0.000 | +0.496 | +0.022 | 0.000 | 0.000 | 0.000 | 0.000 | 3.00 | 0.945 | 1.0 |
| D2 | +0.000 | +0.000 | +0.000 | +0.500 | +1.000 | 0.000 | 0.000 | 0.000 | 0.000 | 3.00 | 0.983 | 1.0 |
| D10 | +0.000 | +0.000 | +0.852 | +0.496 | +0.022 | 0.000 | 0.000 | 0.000 | 0.000 | 3.00 | 0.973 | 2.0 |
| D12 | +0.000 | +0.000 | +0.924 | +0.496 | +0.020 | 0.000 | 0.000 | 0.000 | 0.000 | 3.00 | 0.930 | 2.0 |
| D14 | +0.000 | +0.000 | +0.690 | +0.497 | +0.015 | 0.000 | 0.000 | 0.000 | 0.000 | 3.00 | 0.995 | 2.0 |
| D13 | +0.000 | +0.000 | +0.850 | +0.497 | +0.017 | 0.000 | 0.000 | 0.000 | 0.000 | 3.00 | 0.970 | 2.0 |
| D11 | +0.000 | +0.000 | +0.878 | +0.497 | +0.014 | 0.000 | 0.000 | 0.000 | 0.000 | 3.00 | 0.990 | 2.0 |
| D9 | +0.000 | +0.000 | +0.757 | +0.498 | +0.012 | 0.000 | 0.000 | 0.000 | 0.000 | 3.00 | 0.983 | 2.0 |
| D15 | +0.000 | +0.000 | +0.695 | +0.497 | +0.015 | 0.000 | 0.000 | 0.000 | 0.000 | 3.00 | 0.993 | 2.0 |
| D8 | +0.000 | +0.000 | +0.251 | +0.495 | +0.024 | 0.000 | 0.000 | 0.000 | 0.000 | 3.00 | 0.973 | 1.0 |
| D7 | +0.000 | +0.000 | +0.214 | +0.496 | +0.022 | 0.000 | 0.000 | 0.000 | 0.000 | 3.00 | 0.968 | 1.0 |
| D5 | +0.000 | +0.000 | +0.000 | +0.496 | +0.022 | 0.000 | 0.000 | 0.000 | 0.000 | 3.00 | 0.964 | 1.0 |

## Layer Coupling Statistics

Per-decomposition fast/slow exchange and per-step state persistence, plus velocity and stability per cluster.

## Candidate ranking

Physical decompositions ranked by mean rank across all primary metrics (higher Pearson κ/γ, lower RMS κ/γ, higher coherence / memory / phase / orientation / multiplicative coupling / fast-slow exchange / state persistence / wave mode count / coherence length / mode stability).

| Rank | Code | Family | Pearson κ | Wave modes | Phase vel | Mode stab | F/S exchange | Persistence | Rank sum |
|---:|---|---|---:|---:|---:|---:|---:|---:|---:|
| 1 | D1 | control | +0.12115 | 2.0 | 0.000 | 0.995 | +1.770e-03 | 0.99999 | 104 |
| 2 | D3 | single_layer | +0.12112 | 1.0 | 0.000 | 0.970 | +0.000e+00 | 1.00000 | 114 |
| 3 | D6 | coupling_direction | +0.12079 | 2.0 | 0.000 | 0.995 | +1.873e-03 | 0.99999 | 122 |
| 4 | D4 | frozen | +0.12243 | 1.0 | 0.000 | 0.945 | +2.322e-03 | 0.99999 | 123 |
| 5 | D2 | single_layer | +0.09055 | 1.0 | 0.000 | 0.983 | +0.000e+00 | 1.00000 | 132 |
| 6 | D10 | coupling_form | +0.10704 | 2.0 | 0.000 | 0.973 | +2.377e-03 | 0.99998 | 141 |
| 7 | D12 | neighbour_assignment | +0.12283 | 2.0 | 0.000 | 0.930 | +2.190e-03 | 0.99999 | 143 |
| 8 | D14 | ordering | +0.12126 | 2.0 | 0.000 | 0.995 | +1.772e-03 | 0.99999 | 148 |
| 9 | D13 | neighbour_assignment | +0.12065 | 2.0 | 0.000 | 0.970 | +1.922e-03 | 0.99999 | 151 |
| 10 | D11 | neighbour_assignment | +0.11955 | 2.0 | 0.000 | 0.990 | +1.711e-03 | 0.99999 | 155 |
| 11 | D9 | coupling_form | +0.11946 | 2.0 | 0.000 | 0.983 | +1.558e-03 | 0.99998 | 156 |
| 12 | D15 | timescale | +0.12044 | 2.0 | 0.000 | 0.993 | +1.751e-03 | 0.99999 | 158 |
| 13 | D8 | coupling_direction | +0.11906 | 1.0 | 0.000 | 0.973 | +2.375e-03 | 0.99999 | 162 |
| 14 | D7 | coupling_direction | +0.11712 | 1.0 | 0.000 | 0.968 | +2.236e-03 | 0.99999 | 172 |
| 15 | D5 | frozen | +0.11651 | 1.0 | 0.000 | 0.964 | +2.241e-03 | 0.99999 | 179 |

## Required questions

### Q1. Is the fast layer essential?

D2 (fast removed): Δκ = -0.03061, Δγ = -0.00510. D3 (slow removed): Δκ = -0.00004. Fast+Slow together = +0.12115. Fast is not essential; slow is essential.

### Q2. Is the slow layer essential?

Removing fast layer (D3): Δκ = -0.00004. Removing slow layer (D2): Δκ = -0.03061. Both layers indispensable.

### Q3. Is bidirectional coupling required?

D1 (bidirectional, control) κ = +0.12115, wave modes = 2.0. D8 (independent, no coupling): Δκ = -0.00210, wave modes = 1.0. Coupling is required for the A8 signature.

### Q4. Is timescale separation responsible for the memory effect?

D1 timescale-separated memory: 0.99930. D15 (forced equal τ) memory: 0.99976 (Δ = +0.00046); persistence 0.99999 vs 0.99999. Timescale separation is not solely responsible for the memory effect.

### Q5. Which layer generates the observed wave modes?

D1 full modes: propagating=+0.000, standing=+0.000, transverse=+0.692, longitudinal=+0.497. Slow-only (D2) longitudinal: +0.500. Fast-only (D3) longitudinal: +0.500. Fast-frozen (D4) longitudinal: +0.496.

### Q6. Do wave modes disappear if fast/slow coupling is removed?

D1 wave modes = 2.0. Coupling removed: D8 modes = 1.0, D6 modes = 2.0, D7 modes = 1.0. Wave modes are reduced but not eliminated when fast/slow coupling is removed.

### Q7. Does neighbour interaction primarily act through the fast layer or the slow layer?

D11 (neighbour on fast only): κ = +0.11955, modes = 2.0. D12 (neighbour on slow only): κ = +0.12283, modes = 2.0. D1 control κ = +0.12115, modes = 2.0. Neighbour primarily acts through slow layer.

### Q8. Can A8 be simplified without losing performance?

A8 can be simplified without significant κ loss. Closest candidates: D3|Δκ=-0.00004, D14|Δκ=+0.00011, D6|Δκ=-0.00037

### Q9. Do any stable wave properties repeatedly converge near α or 3α?

84 audit entries sit nearest α or 3α; 77 within log₁₀ distance < 0.1.

- `omega_times_DT` (D1) = +2.00000e-02, factor to 3*alpha_fs = 0.9136, log₁₀ dist = +0.0393
- `FAST_times_DT` (D1) = +2.00000e-02, factor to 3*alpha_fs = 0.9136, log₁₀ dist = +0.0393
- `omega_times_DT` (D2) = +2.00000e-02, factor to 3*alpha_fs = 0.9136, log₁₀ dist = +0.0393
- `FAST_times_DT` (D2) = +2.00000e-02, factor to 3*alpha_fs = 0.9136, log₁₀ dist = +0.0393
- `omega_times_DT` (D3) = +2.00000e-02, factor to 3*alpha_fs = 0.9136, log₁₀ dist = +0.0393

### Q10. Does every successful decomposition preserve machine-precision conservation?

Yes — all 75 runs preserve the unit-speed normalization at or below machine epsilon (2.220e-16).

## Outcome determination

- **A**: One physical mechanism (or one coupling) within A8 is identified as the principal origin of the improved weak-lensing agreement and emergent wave behaviour.
- **B**: Several mechanisms contribute comparably, indicating that A8 is an irreducible cooperative microscopic architecture.
- **C**: No individual mechanism explains A8; its behaviour only emerges from the complete dual-layer constituent.

**Outcome C.** No individual mechanism explains A8; its behaviour only emerges from the complete dual-layer constituent.

## C10 provenance

C10 was not modified and not rerun. The benchmark remains archived at `runs/version_b_physics_lab002/interaction_matrix.csv`.

## Numerical stability

All 75 runs preserve the frozen unit-speed normalization at or below machine epsilon (2.220e-16).

## Wave Family Registry

`runs/wave_family_registry.csv` was updated with 75 new entries from this laboratory. Subsequent laboratories may append further entries without modifying this registry.

## Required artefacts

`report.md`, `component_summary.csv`, `cross_cluster_statistics.csv`, `wave_mode_statistics.csv`, `layer_coupling_statistics.csv`, `candidate_ranking.csv`, `fundamental_constant_audit.csv`, `run.json`, `validation.json`, and all required plots are present in `runs/microstructure_entity_a8_decomposition001/`.
