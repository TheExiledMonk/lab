# PBUF MICROSTRUCTURE-LAB-002

**Collective Neighbourhood Interaction Laboratory in the frozen Version 1 weak-lensing laboratory (LAB-FREEZE-001).**

## Status

- Frozen hash verification: **PASS**
- Production runs: **80**
- Runtime: **58.9 s**
- Fitting or optimisation: **none**

## Frozen laboratory

Only the microscopic neighbourhood interaction rule varies. The constitutive pipeline, transport, source plane, Jacobian observable, and numerical configuration remain byte-identical.

## Interaction families

All step counts, time step, weights, decay scales, radii, and cross-coupling constants were fixed a priori. They are dimensionless or set by the matter field; no sweep or fit was performed.

| Family | Law | Principle |
|---|---|---|
| N1 | First Neighbour Only (Control) | `Ring 1 weight 1.0` |
| N2 | Two-Ring Neighbourhood | `Ring 1: 1.0, Ring 2: 0.5` |
| N3 | Three-Ring Neighbourhood | `Ring 1: 1.0, Ring 2: 0.5, Ring 3: 0.25` |
| N4 | Exponential Decay | `w(d) = exp(-d/lambda), lambda=1.0` |
| N5 | Gaussian Local Field | `w(d) = exp(-d^2 / (2 sigma^2)), sigma=1.2` |
| N6 | Inverse Distance | `w(d) = 1/d` |
| N7 | Inverse Square | `w(d) = 1/d^2` |
| N8 | Finite Radius | `Equal weight within radius R=3` |
| N9 | Orientation-Weighted | `weight *= (1 + 0.5 cos(theta_diff))` |
| N10 | Adaptive Local Field | `weight *= (1 + 0.5 * |grad u|)` |
| N11 | Cooperative Relaxation | `F = (mean_w * neighbours) - u, with damping` |
| N12 | Mean Local Potential | `F = mean_potential of neighbourhood` |
| WR1 | Wrong: Random Weights | `weight = random per step` |
| WR2 | Wrong: Negative Weights | `weight = -|w|` |
| WR3 | Wrong: Extremely Long-Range | `Equal weight across full grid` |
| WR4 | Wrong: Shuffled Topology | `NN-only with shuffled neighbour indices` |

Wrong controls: WR1 random weights, WR2 negative (repulsive) weights, WR3 uniform across full grid (destroys locality), WR4 nearest-neighbour with shuffled topology (tests geometry). They must underperform if the laboratory responds to a meaningful local interaction.

## Emergent index definitions

Emergent Coherence Index = constitutive-gradient-magnitude-weighted mean cosine alignment with 4 neighbours; emergence requires gain > `0.0001`. Emergent Memory Index = mean cosine of successive microscopic increments; persistence requires index >= `0.9` and activity > `1e-06`. Both are computed before photon launch.

## Family summary

| Family | Pearson k | Pearson g | SSIM k | RMS k | Coherence gain | Memory index | Coherence | Memory | Effective radius | Conservation |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| N1 | +0.08475 | +0.10633 | -0.00262 | 0.16441 | +1.251e-03 | 1.00000 | 3/5 | 5/5 | 40.0 | 2.220e-16 |
| N11 | +0.07456 | +0.09347 | -0.00155 | 0.15070 | +3.346e-04 | 1.00000 | 3/5 | 5/5 | 40.0 | 2.220e-16 |
| N9 | +0.08193 | +0.08599 | -0.00073 | 0.15797 | +8.380e-04 | 0.99110 | 3/5 | 5/5 | 40.0 | 2.220e-16 |
| N8 | +0.07705 | +0.09120 | -0.00210 | 0.15170 | +2.405e-03 | 1.00000 | 5/5 | 5/5 | 40.0 | 2.220e-16 |
| N6 | +0.07979 | +0.08329 | -0.00194 | 0.15586 | +1.246e-03 | 1.00000 | 4/5 | 5/5 | 40.0 | 2.220e-16 |
| N4 | +0.08197 | +0.08074 | -0.00187 | 0.15781 | +7.284e-04 | 1.00000 | 3/5 | 5/5 | 40.0 | 2.220e-16 |
| N7 | +0.08602 | +0.07641 | -0.00177 | 0.15811 | +5.365e-04 | 1.00000 | 3/5 | 5/5 | 40.0 | 2.220e-16 |
| N3 | +0.08115 | +0.08026 | -0.00199 | 0.15658 | +1.123e-03 | 1.00000 | 4/5 | 5/5 | 40.0 | 2.220e-16 |
| N2 | +0.08939 | +0.07909 | -0.00176 | 0.16143 | -6.968e-04 | 1.00000 | 2/5 | 5/5 | 40.0 | 2.220e-16 |
| N5 | +0.08832 | +0.07879 | -0.00178 | 0.15968 | -1.946e-04 | 1.00000 | 2/5 | 5/5 | 40.0 | 2.220e-16 |
| N10 | +0.08115 | +0.08026 | -0.00199 | 0.15658 | +1.123e-03 | 1.00000 | 4/5 | 5/5 | 40.0 | 2.220e-16 |
| N12 | +0.06356 | +0.13986 | -0.00948 | 0.19288 | -1.288e-03 | 1.00000 | 1/5 | 5/5 | 40.0 | 2.220e-16 |
| WR3 | +0.08452 | +0.10302 | -0.00160 | 0.13894 | +3.571e-03 | 1.00000 | 4/5 | 5/5 | 140.0 | 2.220e-16 |
| WR1 | +0.06573 | +0.08522 | -0.00142 | 0.15195 | +2.055e-03 | 0.47544 | 5/5 | 0/5 | 40.0 | 2.220e-16 |
| WR4 | +0.08475 | +0.10633 | -0.00262 | 0.16441 | +1.251e-03 | 1.00000 | 3/5 | 5/5 | 4.0 | 2.220e-16 |
| WR2 | +0.06680 | +0.14604 | -0.00992 | 0.19293 | -1.290e-03 | 1.00000 | 1/5 | 5/5 | 40.0 | 2.220e-16 |

## Emergent synergy (Tukey-style, N3 - N1 - N2)

`synergy = N3 - N1 - N2` (three-ring minus first-ring minus two-ring).

- Pearson-kappa synergy: **-0.081196**
- Nonlinear synergy emerged: **YES**

## Cross-cluster validation

| Family | Clusters improving Pearson k | Coherence emergence | Memory emergence | Spatial L | Temporal T | Effective radius |
|---|---:|---:|---:|---:|---:|---:|
| N1 | 0/5 | 3/5 | 5/5 | 3.00 | 4.00 | 40.0 |
| N2 | 4/5 | 2/5 | 5/5 | 3.00 | 4.00 | 40.0 |
| N3 | 3/5 | 4/5 | 5/5 | 3.00 | 4.00 | 40.0 |
| N4 | 3/5 | 3/5 | 5/5 | 3.00 | 4.00 | 40.0 |
| N5 | 3/5 | 2/5 | 5/5 | 3.00 | 4.00 | 40.0 |
| N6 | 3/5 | 4/5 | 5/5 | 3.00 | 4.00 | 40.0 |
| N7 | 3/5 | 3/5 | 5/5 | 3.00 | 4.00 | 40.0 |
| N8 | 2/5 | 5/5 | 5/5 | 3.00 | 4.00 | 40.0 |
| N9 | 2/5 | 3/5 | 5/5 | 3.00 | 4.00 | 40.0 |
| N10 | 3/5 | 4/5 | 5/5 | 3.00 | 4.00 | 40.0 |
| N11 | 3/5 | 3/5 | 5/5 | 3.00 | 4.00 | 40.0 |
| N12 | 1/5 | 1/5 | 5/5 | 3.00 | 5.00 | 40.0 |
| WR1 | 1/5 | 5/5 | 0/5 | 3.00 | 4.00 | 40.0 |
| WR2 | 1/5 | 1/5 | 5/5 | 3.00 | 0.00 | 40.0 |
| WR3 | 2/5 | 4/5 | 5/5 | 3.00 | 4.00 | 140.0 |
| WR4 | 0/5 | 3/5 | 5/5 | 3.00 | 4.00 | 4.0 |

## Required questions

### Q1. Does collective neighbourhood interaction outperform pairwise interaction?

The pairwise control N1 reaches median Pearson kappa +0.08475; the best physical neighbourhood family reaches +0.08939 (N2). Collective interaction outperforms the pairwise control.

### Q2. Which interaction decay law performs best?

Among pure-decay families (N2-N8), the best is **N2 — Two-Ring Neighbourhood** at Pearson kappa +0.08939, RMS kappa 0.16143.

### Q3. Does anisotropy remain important once neighbourhood effects are included?

N9 (orientation-weighted) reaches Pearson kappa +0.08193, spatial correlation 3.00; anisotropy still contributes positively relative to the N3 baseline.

### Q4. Does memory emerge naturally from neighbourhood relaxation?

N11 (cooperative relaxation with damping) reaches emergent memory index 1.00000 and activity 0.05835. Memory emerges naturally on all five clusters.

### Q5. Does positive synergy return?

Nonlinear synergy N3 - N1 - N2 = -0.081196. NO — synergy is not positive; the previously observed positive synergy does not return under collective neighbourhood interaction.

### Q6. Does any neighbourhood interaction outperform C10?

No. C10 remains at Pearson kappa +0.10340; no physical neighbourhood family simultaneously exceeds both primary metrics.

### Q7. Is the interaction fundamentally local or semi-local?

Median effective interaction radius across physical families: 40.0 neighbours; median spatial correlation length: 3.00 pixels. The interaction is semi-local, extending beyond the nearest ring.

### Q8. Does performance improve consistently across all five clusters?

Per-family improvement counts — see cross-cluster table. Top improvers: N2=4/5, N3=3/5, N4=3/5.

### Q9. Does every successful interaction preserve machine-precision conservation?

Yes. All 12 physical families have maximum speed-normalisation error <= 2.220e-16.

## Wrong-control diagnostics

| Wrong family | Pearson k | Coherence | Memory | Conservation |
|---|---:|---:|---:|---:|
| WR1 — Wrong: Random Weights | +0.06573 | 5/5 | 0/5 | 2.220e-16 |
| WR2 — Wrong: Negative Weights | +0.06680 | 1/5 | 5/5 | 2.220e-16 |
| WR3 — Wrong: Extremely Long-Range | +0.08452 | 4/5 | 5/5 | 2.220e-16 |
| WR4 — Wrong: Shuffled Topology | +0.08475 | 3/5 | 5/5 | 2.220e-16 |

These deliberately wrong interactions are included to verify that the laboratory is not merely responding to added complexity or destroyed locality.

## Outcome determination

**Outcome B.** Neighbourhood interactions outperform pairwise interactions but remain inferior to C10.

## C10 provenance

Archived reference: `runs/version_b_physics_lab002/interaction_matrix.csv`, SHA-256 `6f0d83d691296e07727d721c8a080e850c743afd6e7455a430b7eb17073116a9`. Not rerun or modified.

## Numerical stability

All 80 runs preserve the frozen unit-speed normalization at or below machine epsilon (2.220e-16).

## Required artefacts

`interaction_summary.csv`, `cross_cluster_statistics.csv`, `interaction_radius.csv`, `correlation_statistics.csv`, `synergy_statistics.csv`, `candidate_ranking.csv`, `run.json`, `validation.json`, and all eight requested plots are present in `runs/microstructure_lab002/`.
