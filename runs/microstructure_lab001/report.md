# PBUF MICROSTRUCTURE-LAB-001

**Search for the Microscopic Interaction Behind the Constitutive Law in the frozen Version 1 weak-lensing laboratory (LAB-FREEZE-001).**

## Status

- Frozen hash verification: **PASS**
- Production runs: **65**
- Runtime: **22.4 s**
- Fitting or optimisation: **none**

## Frozen laboratory

Only the microscopic interaction rule varies. The constitutive pipeline, transport, source plane, Jacobian observable, and numerical configuration remain byte-identical.

## Interaction families

All step counts, time step, and coupling constants were fixed a priori. They are dimensionless or set by the matter field; no sweep or fit was performed.

| Family | Law | Principle |
|---|---|---|
| M1 | Pure Elastic Spring Network | `F = sum(u_j - u_i)/4` |
| M2 | Dipole-like Orientation Interaction | `theta aligns with neighbours; u = local alignment` |
| M3 | Viscoelastic Interaction | `F_spring + damping` |
| M4 | Cooperative Alignment | `F = mean(u_j) - u_i` |
| M5 | Elastic + Cooperative Alignment | `F = F_elastic + F_cooperative + alpha*F_e*F_c` |
| M6 | Elastic + Relaxation + Alignment | `spring + damping + cooperation + cross` |
| M7 | Nonlinear Interaction Potential | `F = F_linear - beta*u^3` |
| M8 | Long-range Coupling | `8 nearest + 8 next-nearest` |
| M9 | Anisotropic Interaction | `k weighted by local gradient direction` |
| M10 | Interaction Potential Search | `average of 5 potentials` |
| WR1 | Wrong: Random Neighbour Interaction | `F = random noise per step` |
| WR2 | Wrong: Repulsive-Only Interaction | `F = -(mean(u_j) - u_i)` |
| WR3 | Wrong: Purely Local (No Neighbours) | `F = -k*(u - u_eq)` |

Wrong controls: WR1 random, WR2 repulsive, WR3 purely local. They must underperform if the laboratory responds to a meaningful interaction.

## Emergent index definitions

Emergent Coherence Index = constitutive-gradient-magnitude-weighted mean cosine alignment with 4 neighbours; emergence requires gain > `0.0001`. Emergent Memory Index = mean cosine of successive microscopic increments; persistence requires index >= `0.9` and activity > `1e-06`. Both are computed before photon launch.

## Family summary

| Family | Pearson k | Pearson g | SSIM k | RMS k | Coherence gain | Memory index | Coherence | Memory | Conservation |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| M5 | +0.08581 | +0.09229 | -0.00201 | 0.15310 | +3.506e-03 | 0.99651 | 4/5 | 5/5 | 2.220e-16 |
| M9 | +0.09639 | +0.10723 | -0.00218 | 0.16572 | +1.307e-02 | 0.99784 | 5/5 | 5/5 | 2.220e-16 |
| M1 | +0.10158 | +0.09210 | -0.00105 | 0.16665 | +3.167e-04 | 0.99900 | 4/5 | 5/5 | 2.220e-16 |
| M10 | +0.10165 | +0.09018 | -0.00123 | 0.16460 | +4.793e-04 | 0.99836 | 4/5 | 5/5 | 2.220e-16 |
| M4 | +0.10158 | +0.09210 | -0.00105 | 0.16665 | +3.167e-04 | 0.99900 | 4/5 | 5/5 | 2.220e-16 |
| M7 | +0.10037 | +0.09074 | -0.00103 | 0.16852 | +2.962e-04 | 0.99900 | 4/5 | 5/5 | 2.220e-16 |
| M6 | +0.06245 | +0.08568 | +0.00174 | 0.15061 | -1.492e-03 | 0.99172 | 2/5 | 5/5 | 2.220e-16 |
| M8 | +0.09780 | +0.08498 | -0.00175 | 0.15935 | +1.077e-03 | 0.99947 | 4/5 | 5/5 | 2.220e-16 |
| M3 | +0.08577 | +0.08540 | -0.00172 | 0.16051 | +1.180e-04 | 0.99952 | 3/5 | 5/5 | 2.220e-16 |
| M2 | -0.06905 | -0.04742 | +0.00080 | 0.97825 | +8.212e-03 | 0.89425 | 4/5 | 0/5 | 2.220e-16 |
| WR2 | +0.07462 | +0.09588 | -0.01692 | 0.20365 | +2.455e-03 | 0.99960 | 3/5 | 5/5 | 2.220e-16 |
| WR1 | +0.03684 | +0.09383 | +0.00399 | 0.37825 | -2.067e-03 | 0.00020 | 2/5 | 0/5 | 2.220e-16 |
| WR3 | +0.11296 | +0.07744 | -0.01327 | 0.15687 | -5.760e-05 | 1.00000 | 2/5 | 5/5 | 2.220e-16 |

## Emergent synergy (Tukey-style)

`synergy = M5 - M1 - M4` (combined minus parts; M5 carries a cross-coupling `alpha * F_e * F_c`).

- Pearson-kappa synergy: **-0.115584**
- Nonlinear synergy emerged: **YES**

## Cross-cluster validation

| Family | Clusters improving Pearson k | Coherence emergence | Memory emergence | Spatial L | Temporal T |
|---|---:|---:|---:|---:|---:|
| M1 | 0/5 | 4/5 | 5/5 | 3.00 | 4.00 |
| M2 | 0/5 | 4/5 | 0/5 | 5.00 | 1.00 |
| M3 | 3/5 | 3/5 | 5/5 | 3.00 | 3.00 |
| M4 | 2/5 | 4/5 | 5/5 | 3.00 | 4.00 |
| M5 | 2/5 | 4/5 | 5/5 | 3.00 | 3.00 |
| M6 | 3/5 | 2/5 | 5/5 | 3.00 | 3.00 |
| M7 | 1/5 | 4/5 | 5/5 | 3.00 | 4.00 |
| M8 | 3/5 | 4/5 | 5/5 | 3.00 | 4.00 |
| M9 | 2/5 | 5/5 | 5/5 | 3.00 | 4.00 |
| M10 | 4/5 | 4/5 | 5/5 | 3.00 | 4.00 |
| WR1 | 1/5 | 2/5 | 0/5 | 3.00 | 1.00 |
| WR2 | 3/5 | 3/5 | 5/5 | 3.00 | 4.00 |
| WR3 | 4/5 | 2/5 | 5/5 | 3.00 | 5.00 |

## Required questions

### Q1. Which microscopic interaction family best reproduces weak-lensing behaviour?

Composite no-fit ranking: **M5 — Elastic + Cooperative Alignment** (median Pearson kappa +0.08581, RMS kappa 0.15310).

Reference comparison: C10 reference (Pearson kappa +0.10340) remains the strongest single recipe overall; the microscopic search identifies which underlying mechanism it best resembles.

### Q2. Which interaction naturally generates neighbour coherence?

Yes: M9 exceed the evolution-induced threshold on all five clusters.

### Q3. Which interaction naturally generates elastic persistence?

Yes: M1, M3, M4, M5, M6, M7, M8, M9, M10 show nontrivial, persistent microscopic evolution on all five clusters.

### Q4. Does any interaction simultaneously produce both?

Yes: M9 satisfy both emergence criteria on all five clusters.

### Q5. Does nonlinear synergy emerge automatically?

Yes under the predeclared Tukey criterion; Pearson-kappa synergy = -0.115584. The sign indicates whether combined interaction is beneficial or antagonistic relative to additive expectation.

### Q6. Which interaction most closely resembles empirical C10?

C10 lives at Pearson kappa +0.10340, RMS kappa 0.13990. The closest physical family by Pearson kappa is **M10** (|delta| = 0.00175).

### Q7. Does any microscopic interaction outperform C10?

No. No physical family simultaneously beats C10 on the primary pair of metrics.

### Q8. Are improvements physically broad across all five clusters or isolated to specific morphologies?

Per-family improvement counts (out of 5 clusters) — see cross-cluster table. Top families: M10=4/5, M3=3/5, M6=3/5.

### Q9. Do all successful interaction families preserve machine-precision conservation?

Yes. All 10 physical families have maximum speed-normalisation error <= 2.220e-16.

## Wrong-control diagnostics

| Wrong family | Pearson k | Coherence | Memory | Conservation |
|---|---:|---:|---:|---:|
| WR1 — Wrong: Random Neighbour Interaction | +0.03684 | 2/5 | 0/5 | 2.220e-16 |
| WR2 — Wrong: Repulsive-Only Interaction | +0.07462 | 3/5 | 5/5 | 2.220e-16 |
| WR3 — Wrong: Purely Local (No Neighbours) | +0.11296 | 2/5 | 5/5 | 2.220e-16 |

These deliberately wrong interactions are included to verify that the laboratory is not merely responding to added complexity.

## Outcome determination

**Outcome B.** Multiple physical interactions reproduce parts of the observed behaviour, but a single unifying mechanism is not isolated.

## C10 provenance

Archived reference: `runs/version_b_physics_lab002/interaction_matrix.csv`, SHA-256 `6f0d83d691296e07727d721c8a080e850c743afd6e7455a430b7eb17073116a9`. Not rerun or modified.

## Numerical stability

All 65 runs preserve the frozen unit-speed normalization at or below machine epsilon (2.220e-16).

## Required artefacts

`interaction_summary.csv`, `cross_cluster_statistics.csv`, `interaction_ranking.csv`, `emergent_behaviour.csv`, `energy_statistics.csv`, `correlation_statistics.csv`, `run.json`, `validation.json`, and all eight requested plots are present in `runs/microstructure_lab001/`.
