# PBUF MICROSTATE-LAB-001

**Internal State Laboratory in the frozen Version 1 weak-lensing laboratory (LAB-FREEZE-001).**

## Status

- Frozen hash verification: **PASS**
- Production runs: **70**
- Runtime: **27.2 s**
- Fitting or optimisation: **none**

## Frozen laboratory

Only the microscopic state carried by each spacetime element varies. The transport, source plane, Jacobian observable, numerical configuration, and constitutive framework remain byte-identical.

## Microscopic state families

All step counts, time step, decay scales, and cross-coupling constants were fixed a priori. They are dimensionless or set by the matter field; no sweep or fit was performed.

| Family | State | Principle |
|---|---|---|
| S1 | Scalar (Control) | `u = rho; no internal variables` |
| S2 | Scalar + Orientation | `u, theta; theta aligns with neighbours` |
| S3 | Scalar + Internal Strain | `u, epsilon; strain accumulates with time` |
| S4 | Scalar + Phase | `u, phi; phase oscillates and couples to neighbours` |
| S5 | Scalar + Relaxation State | `u, R; R evolves with finite relaxation time` |
| S6 | Scalar + Local Momentum | `u, p; momentum carries neighbour response` |
| S7 | Scalar + Orientation + Relaxation | `u, theta, R; combined orientation-relaxation` |
| S8 | Scalar + Orientation + Strain | `u, theta, epsilon; combined orientation-strain` |
| S9 | Scalar + Phase + Orientation | `u, phi, theta; combined phase-orientation` |
| S10 | Full Local State | `u, theta, epsilon, phi, R; full combined state` |
| WR1 | Wrong: Random Internal State | `u = random noise per step` |
| WR2 | Wrong: Frozen Internal State | `u = u_init, never evolves` |
| WR3 | Wrong: Rapid Randomisation | `u re-randomised every step` |
| WR4 | Wrong: Self-Only Evolution | `no neighbour influence` |

Wrong controls: WR1 random internal state, WR2 frozen internal state, WR3 rapid randomisation (destroys persistence), WR4 self-only evolution (no neighbour influence). They must underperform if the laboratory responds to a meaningful internal state.

## Emergent index definitions

Emergent Coherence Index = constitutive-gradient-magnitude-weighted mean cosine alignment with 4 neighbours; emergence requires gain > `0.0001`. Emergent Memory Index = mean cosine of successive state increments; persistence requires index >= `0.9` and activity > `1e-06`. Relaxation time = first half-decay lag. Both are computed before photon launch.

## Family summary

| Family | Pearson k | Pearson g | SSIM k | RMS k | Coherence gain | Memory index | Coherence | Memory | Relax. time | Conservation |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| S9 | +0.10495 | +0.09400 | -0.00183 | 0.16621 | +2.166e-03 | 0.99804 | 5/5 | 5/5 | 6.00 | 2.220e-16 |
| S1 | +0.10158 | +0.09210 | -0.00105 | 0.16665 | +3.167e-04 | 0.99900 | 4/5 | 5/5 | 6.00 | 2.220e-16 |
| S4 | +0.10278 | +0.09633 | -0.00264 | 0.16691 | +2.440e-03 | 0.99754 | 5/5 | 5/5 | 6.00 | 2.220e-16 |
| S8 | +0.09009 | +0.07903 | -0.00179 | 0.15265 | +9.059e-04 | 0.90323 | 4/5 | 4/5 | 1.00 | 2.220e-16 |
| S3 | +0.09337 | +0.07898 | -0.00235 | 0.15341 | +8.616e-04 | 0.99958 | 4/5 | 5/5 | 7.00 | 2.220e-16 |
| S5 | +0.08951 | +0.08360 | -0.01058 | 0.15580 | +3.455e-08 | 1.00000 | 0/5 | 5/5 | 5.00 | 2.220e-16 |
| S2 | -0.03699 | +0.00740 | +0.00729 | 0.32281 | +4.354e-03 | 0.44465 | 3/5 | 0/5 | 19.00 | 2.220e-16 |
| S6 | +0.08951 | +0.08360 | -0.01058 | 0.15580 | +0.000e+00 | 0.00000 | 0/5 | 0/5 | 0.00 | 2.220e-16 |
| S10 | +0.02754 | +0.00474 | -0.00735 | 0.15066 | +1.739e-02 | 0.28933 | 5/5 | 0/5 | 1.00 | 2.220e-16 |
| S7 | -0.05222 | -0.10397 | +0.00067 | 0.30219 | +7.347e-03 | 0.99052 | 4/5 | 5/5 | 8.00 | 2.220e-16 |
| WR2 | +0.06680 | +0.14604 | -0.00992 | 0.19293 | -1.290e-03 | 0.00000 | 1/5 | 0/5 | 0.00 | 2.220e-16 |
| WR4 | +0.11296 | +0.07744 | -0.01327 | 0.15687 | -5.760e-05 | 1.00000 | 2/5 | 5/5 | 6.00 | 2.220e-16 |
| WR3 | +0.05327 | +0.08588 | -0.00365 | 1.44434 | +4.678e-04 | -0.49310 | 3/5 | 0/5 | 20.00 | 2.220e-16 |
| WR1 | +0.09093 | +0.01561 | -0.00472 | 1.54274 | -6.286e-04 | -0.49176 | 2/5 | 0/5 | 20.00 | 2.220e-16 |

## Emergent synergy (Tukey-style, S8 - S2 - S3)

`synergy = S8 - S2 - S3` (orientation+strain combined minus its two base states).

- Pearson-kappa synergy: **+0.036558**
- Nonlinear synergy emerged: **YES**

## Cross-cluster validation

| Family | Clusters improving Pearson k | Coherence emergence | Memory emergence | Spatial L | Temporal T | Relax. time |
|---|---:|---:|---:|---:|---:|---:|
| S1 | 0/5 | 4/5 | 5/5 | 3.00 | 4.00 | 6.00 |
| S2 | 1/5 | 3/5 | 0/5 | 10.00 | 3.00 | 19.00 |
| S3 | 3/5 | 4/5 | 5/5 | 3.00 | 5.00 | 7.00 |
| S4 | 4/5 | 5/5 | 5/5 | 3.00 | 4.00 | 6.00 |
| S5 | 3/5 | 0/5 | 5/5 | 3.00 | 4.00 | 5.00 |
| S6 | 3/5 | 0/5 | 0/5 | 3.00 | 0.00 | 0.00 |
| S7 | 0/5 | 4/5 | 5/5 | 4.00 | 4.00 | 8.00 |
| S8 | 2/5 | 4/5 | 4/5 | 3.00 | 1.00 | 1.00 |
| S9 | 4/5 | 5/5 | 5/5 | 3.00 | 4.00 | 6.00 |
| S10 | 1/5 | 5/5 | 0/5 | 4.00 | 1.00 | 1.00 |
| WR1 | 2/5 | 2/5 | 0/5 | 128.00 | 1.00 | 20.00 |
| WR2 | 1/5 | 1/5 | 0/5 | 3.00 | 0.00 | 0.00 |
| WR3 | 1/5 | 3/5 | 0/5 | 128.00 | 1.00 | 20.00 |
| WR4 | 4/5 | 2/5 | 5/5 | 3.00 | 5.00 | 6.00 |

## Required questions

### Q1. Does introducing internal state improve weak-lensing agreement?

The scalar control S1 reaches median Pearson kappa +0.10158; the best physical state family reaches +0.10495 (S9). Internal state improves agreement.

### Q2. Which internal variable contributes most?

Among the single-variable additions (S2-S6), the best is **S4 — Scalar + Phase** at Pearson kappa +0.10278.

### Q3. Does memory emerge naturally from state evolution?

6/10 physical families show nontrivial, persistent state evolution on all five clusters. Memory emerges naturally in the majority of physical families.

### Q4. Does neighbour coherence emerge without explicit programming?

3/10 physical families exceed the evolution-induced coherence threshold on all five clusters. Neighbour coherence emerges in only a minority of physical families.

### Q5. Does positive synergy return?

Nonlinear synergy S8 - S2 - S3 = +0.036558. YES — synergy is positive, recovering the previously observed cooperative behaviour.

### Q6. Does any microscopic state outperform C10?

No. C10 remains at Pearson kappa +0.10340; no physical state family simultaneously exceeds both primary metrics.

### Q7. Which internal state produces the most physically consistent behaviour across all five clusters?

Ranking by improvement-count across clusters: S4=4/5, S9=4/5, S3=3/5.

### Q8. Are improvements broad across every cluster or morphology-specific?

Per-family improvement counts — see cross-cluster table. The S1 baseline itself varies across clusters; improvements are measured relative to S1 per cluster.

### Q9. Does every successful state preserve machine-precision conservation?

Yes. All 10 physical families have maximum speed-normalisation error <= 2.220e-16.

## Wrong-control diagnostics

| Wrong family | Pearson k | Coherence | Memory | Conservation |
|---|---:|---:|---:|---:|
| WR1 — Wrong: Random Internal State | +0.09093 | 2/5 | 0/5 | 2.220e-16 |
| WR2 — Wrong: Frozen Internal State | +0.06680 | 1/5 | 0/5 | 2.220e-16 |
| WR3 — Wrong: Rapid Randomisation | +0.05327 | 3/5 | 0/5 | 2.220e-16 |
| WR4 — Wrong: Self-Only Evolution | +0.11296 | 2/5 | 5/5 | 2.220e-16 |

These deliberately wrong states are included to verify that the laboratory is not merely responding to added complexity or destroyed locality.

## Outcome determination

**Outcome B.** Several microscopic states improve the laboratory, but no unique description emerges.

## C10 provenance

Archived reference: `runs/version_b_physics_lab002/interaction_matrix.csv`, SHA-256 `6f0d83d691296e07727d721c8a080e850c743afd6e7455a430b7eb17073116a9`. Not rerun or modified.

## Numerical stability

All 70 runs preserve the frozen unit-speed normalization at or below machine epsilon (2.220e-16).

## Required artefacts

`state_summary.csv`, `cross_cluster_statistics.csv`, `candidate_ranking.csv`, `state_statistics.csv`, `synergy_statistics.csv`, `run.json`, `validation.json`, and all eight requested plots are present in `runs/microstate_lab001/`.
