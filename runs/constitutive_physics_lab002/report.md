# PBUF CONSTITUTIVE-PHYSICS-LAB-002

**Search for the Constitutive Law in the frozen Version 1 weak-lensing laboratory.**

## Status

- Frozen hash verification: **PASS**
- Production runs: **35**
- Runtime: **10.2 s**
- Fitting or optimisation: **none**

## Frozen laboratory

Transport, source plane, Jacobian observable, numerical configuration, and validation components remain byte-identical. Only the scalar constitutive state supplied to the frozen transverse response is evolved.

## Constitutive laws

All iteration counts and dimensionless coefficients were fixed before inspection of results. They were not swept or fitted.

| Family | Law | Principle |
|---|---|---|
| F1 | Instantaneous Constitutive Response | `C=Ceq` |
| F2 | Relaxation Constitutive Law | `dC/ds=(Ceq-C)/tau` |
| F3 | Local Constitutive Evolution | `dC/ds=Laplacian(C)` |
| F4 | Constitutive Energy Functional | `min integral[(C-Ceq)^2+|grad C|^2]` |
| F5 | Gradient-Driven Constitutive Evolution | `dC/ds=div(g(|grad C|) grad C)` |
| F6 | Relaxation + Neighbour Evolution | `dC/ds=(Ceq-C)/tau+Laplacian(C)` |
| F7 | Variational Constitutive Law | `gradient flow of convex local elastic energy with continuity` |

The fixed discretizations are: F2 eight steps at 0.25; F3 eight four-neighbour diffusion steps at 0.20; F4 twelve Jacobi minimization steps with equal fidelity and continuity weights; F5 eight edge-stopping gradient-flow steps at 0.20 with scale fixed by the initial lattice-gradient RMS; F6 eight unified steps with relaxation 0.25 and neighbour evolution 0.15; F7 twelve convex variational-flow steps at 0.15 with continuity weight 0.50 and normalized quartic weight 0.25.

## Emergent-index definitions

The Emergent Coherence Index is the constitutive-gradient-magnitude-weighted mean cosine alignment with four neighbours. Emergence requires final-minus-initial gain > `0.0001`. The Emergent Memory Index is the mean cosine similarity of successive constitutive-state increments. Persistence requires index >= `0.9` and normalized evolution activity > `1e-06`. Both are computed before photon launch; absolute input smoothness is not counted as emergence.

## Family summary

| Rank | Family | Pearson k | Pearson g | SSIM k | RMS k | Coherence gain | Memory index | Improved clusters | Conservation |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | F3 | +0.09371 | +0.07931 | -0.00113 | 0.14066 | +7.770e-02 | 0.99804 | 2/5 | 2.220e-16 |
| 2 | F7 | +0.08857 | +0.08411 | -0.00553 | 0.13584 | +4.707e-02 | 0.99995 | 1/5 | 2.220e-16 |
| 3 | F4 | +0.09053 | +0.07945 | -0.00445 | 0.14616 | +6.452e-02 | 0.93850 | 2/5 | 2.220e-16 |
| 4 | F5 | +0.09037 | +0.08499 | -0.01012 | 0.15567 | +7.717e-03 | 0.99828 | 3/5 | 2.220e-16 |
| 5 | F6 | +0.08868 | +0.08271 | -0.00519 | 0.13977 | +5.085e-02 | 0.99995 | 1/5 | 2.220e-16 |
| 6 | F1 | +0.08951 | +0.08360 | -0.01058 | 0.15580 | +0.000e+00 | 0.00000 | 0/5 | 2.220e-16 |
| 7 | F2 | +0.08822 | +0.08268 | -0.00866 | 0.14401 | +7.197e-10 | 1.00000 | 2/5 | 2.220e-16 |

## Emergent synergy

For the unified F6 law, Tukey nonadditivity is `F6 - F2 - F3 + F1`. This is evaluated in final lensing metrics and directly in the constitutive state before propagation.

- Pearson kappa interaction: **-0.003743**
- Pearson gamma interaction: **+0.004330**
- Median constitutive nonadditivity index: **3.999880e-03**
- Nonlinear synergy emerged: **YES**

## Cross-cluster validation

| Family | Clusters improving Pearson kappa | Coherence emergence | Memory emergence |
|---|---:|---:|---:|
| F1 | 0/5 | 0/5 | 0/5 |
| F2 | 2/5 | 0/5 | 5/5 |
| F3 | 2/5 | 5/5 | 5/5 |
| F4 | 2/5 | 5/5 | 4/5 |
| F5 | 3/5 | 5/5 | 5/5 |
| F6 | 1/5 | 5/5 | 5/5 |
| F7 | 1/5 | 5/5 | 5/5 |

## Required questions

### Q1. Does any constitutive family naturally reproduce neighbour coherence?

Yes: F3, F4, F5, F6, F7 exceed the evolution-induced threshold on all five clusters.

### Q2. Does any constitutive family naturally reproduce elastic persistence?

Yes: F2, F3, F5, F6, F7 show nontrivial, persistent constitutive evolution on all five clusters.

### Q3. Does nonlinear synergy emerge without explicitly programming it?

Yes for nonlinear nonadditivity under the predeclared dual criterion. The Pearson-kappa interaction is -0.003743, so it is antagonistic rather than automatically equivalent to the previously observed positive synergy.

### Q4. Which constitutive family gives the greatest improvement over Version A?

The composite no-fit ranking selects **F3 — Local Constitutive Evolution**, with median Pearson-kappa change +0.00420 and improvement on 2/5 clusters.

### Q5. Which family best explains the previous coherence-memory interaction?

**F6** is the closest structural explanation because relaxation and neighbour evolution coexist in one state equation and are nonadditive. However, its Pearson-kappa interaction has the opposite sign from the previous positive synergy, so it is not a complete explanation.

### Q6. Does any family outperform manually combined C10?

No. C10 remains at Pearson kappa +0.10340 and RMS kappa 0.13990; no constitutive family beats both.

### Q7. Do all successful families preserve machine-precision conservation?

Yes. All 1 successful families have maximum speed-normalisation error <= 2.220e-16.

## Outcome determination

**Outcome A.** At least one successful constitutive family naturally reproduces both predeclared behaviours.

## C10 provenance

Archived reference: `runs/version_b_physics_lab002/interaction_matrix.csv`, SHA-256 `6f0d83d691296e07727d721c8a080e850c743afd6e7455a430b7eb17073116a9`. It was not rerun or modified.

## Numerical stability

All 35 runs preserve the frozen unit-speed normalization at or below machine epsilon (2.220e-16).

## Required artefacts

`family_summary.csv`, `cross_cluster_statistics.csv`, `emergent_behaviour.csv`, `constitutive_ranking.csv`, `run.json`, `validation.json`, and all six requested plots are present in `runs/constitutive_physics_lab002/`.
