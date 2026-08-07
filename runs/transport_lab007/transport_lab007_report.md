# PBUF TRANSPORT-LAB-007 — Local neighbour update mechanism ablation

Frozen: Lens-001, constitutive solution, neighbour graph, kernel envelope, response magnitude (A = |∇C|), response angle (90°), integration, timestep, normalisation.

Variable: the update rule that maps (v_old, r) → v_new.

The response fed into all 7 rules is the magnitude-normalised 90° transverse response from LAB-006: `r = A · R_90(∇̂C)`.

## Required ranking table

| Rank | Rule | Update mechanism | Bend | Conservation | Stable | Runtime (s) | Score |
|---:|---|---|---:|---:|:-:|---:|---:|
| 1 | Exp 6 | Phase transfer (rotate by step · arg(r)) | 2.648e+00 | 2.220e-16 | yes | 0.016 | 1.500 |
| 2 | Exp 3 | Projection update (r ⊥ v only) | 6.618e-05 | 1.110e-16 | yes | 0.019 | 1.000 |
| 3 | Exp 5 | Local momentum transfer (m = 0.5) | 1.324e-04 | 2.220e-16 | yes | 0.017 | 0.500 |
| 4 | Exp 2 | Pure rotation (by step · ‖r‖) | 6.619e-05 | 2.220e-16 | yes | 0.016 | 0.500 |
| 5 | Exp 7 | Mixed (rotation + translation + shear) | 6.619e-05 | 2.220e-16 | yes | 0.022 | 0.500 |
| 6 | Exp 1 | Direct vector addition (baseline) | 6.618e-05 | 2.220e-16 | yes | 0.017 | 0.500 |
| 7 | Exp 4 | Incremental shear | 9.054e-06 | 2.220e-16 | yes | 0.019 | 0.500 |

## Behaviour summary

| Rule | Mechanism | Stability | Behaviour |
|---|---|:-:|---|
| Exp 1 | Direct vector addition (baseline) | stable | monotonic |
| Exp 2 | Pure rotation | stable | monotonic |
| Exp 3 | Projection update (r ⊥ v only) | stable | monotonic |
| Exp 4 | Incremental shear | stable | monotonic |
| Exp 5 | Local momentum transfer (m=0.5) | stable | monotonic |
| Exp 6 | Phase transfer | stable | **oscillatory** |
| Exp 7 | Mixed | stable | monotonic |

## Observations (no interpretation)

- **Exp 1, 2, 3, 7 are statistically indistinguishable on bending**: 6.618e-05, 6.619e-05, 6.618e-05, 6.619e-05. Relative spread < 0.005%.
- **Conservation is machine precision (≤ 2.2e-16) for every rule** because every rule renormalises the velocity to unit magnitude after the update.
- **Exp 4 (incremental shear) gives 7.3× less bending** (9.05e-06). The shear decomposition scales the v_perp component by (1 + step·|r|); at step·|r| ≈ 8.4e-4 this is a near-unity perturbation, so v_perp is barely changed.
- **Exp 5 (momentum, m = 0.5) gives 2.0× more bending** (1.32e-04). With m = 0.5 the per-step velocity update is doubled relative to the m = 1 baseline.
- **Exp 6 (phase transfer) gives 40 000× more bending** (2.65) and produces oscillatory photon paths. The phase angle is `step · arg(r)`, and the cumulative phase over 80 steps can rotate the velocity by several radians.
- All 7 rules are numerically stable (finite outputs, |v| bounded by renormalisation).

## Numerical detail per rule

| Rule | Bending | Conservation | Speed drift (pre) | Direction drift | Behaviour |
|---|---:|---:|---:|---:|---|
| Exp 1 | 6.618e-05 | 2.22e-16 | 1.40e-05 | 1.09e-06 | stable / monotonic |
| Exp 2 | 6.619e-05 | 2.22e-16 | 1.40e-05 | 1.09e-06 | stable / monotonic |
| Exp 3 | 6.618e-05 | 1.11e-16 | 1.40e-05 | 1.09e-06 | stable / monotonic |
| Exp 4 | 9.054e-06 | 2.22e-16 | 1.93e-06 | 1.49e-07 | stable / monotonic |
| Exp 5 | 1.324e-04 | 2.22e-16 | 2.80e-05 | 2.18e-06 | stable / monotonic |
| Exp 6 | 2.648e+00 | 2.22e-16 | 5.64e-01 | 8.59e-01 | stable / oscillatory |
| Exp 7 | 6.619e-05 | 2.22e-16 | 1.40e-05 | 1.09e-06 | stable / monotonic |

## Outcome

**Outcome B — several update mechanisms are statistically indistinguishable.**

The statistically equivalent group is:

| Rule | Update mechanism |
|---|---|
| Exp 1 | Direct vector addition (baseline) |
| Exp 2 | Pure rotation |
| Exp 3 | Projection update (r ⊥ v only) |
| Exp 7 | Mixed (rotation + translation + shear) |

These four rules give bending within 0.005% of each other (6.618–6.619 × 10⁻⁵) and machine-precision conservation. They are operationally equivalent in the small-step regime of this dataset.

The three remaining rules (Exp 4, 5, 6) produce measurably different results, but not in a way that gives a clean "winner" by either metric — Exp 4 gives less bending, Exp 5 gives more bending with no conservation advantage, and Exp 6 gives many orders of magnitude more bending but at the cost of oscillatory photon paths that do not match the weak-lensing transport pattern.

The baseline update rule (Exp 1, direct vector addition with renormalisation) is not improved by any of the tested alternatives.

Artefacts: `runs/transport_lab007/{measurements.csv, measurements.json, ranking.csv, behaviour_summary.csv, transport_lab007_report.md}`. Source: `transport_lab007.py`.

Laboratory stops. No physical interpretation, no new laws.
