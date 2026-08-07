# PBUF TRANSPORT-LAB-001 — Measurements

Laboratory only. No interpretation.

## Frozen inputs (identical for every experiment)

| Field | Source | Shape | Min / Max |
|---|---|---|---|
| Matter | `runs/wl001/matter.csv` | 128×128 | 2.56e-54 / 9.95e-01 |
| Deformation (C) | `runs/wl001/deformation.csv` | 128×128 | 4.64e-55 / 1.80e-01 |
| Gradient_x | `runs/wl001/gradient_x.csv` | 128×128 | -1.44e-01 / 1.44e-01 |
| Gradient_y | `runs/wl001/gradient_y.csv` | 128×128 | -1.44e-01 / 1.44e-01 |
| Observation | `runs/wl001/observation.csv` | 128×128 | 5.46e-276 / 9.84e-01 |
| P (\|∇C\|) | derived | 128×128 | 0 / 2.04e-01 |
| W (½C²) | derived | 128×128 | 0 / 1.62e-02 |
| N (radial normal) | derived | 2×128×128 | unit |

Grid: `[-8, 8]²`, `n = 128`, `dx = dy = 16/127`.
Photon ensemble: 9 rays, `y₀ ∈ [-3, 3]`, 80 steps, `step_size = 0.06`.

## Propagation kernel (identical in every experiment)

```
current neighbour
   → compute local response     [varies per experiment]
   → v ← v − step * response
   → v ← v / |v|
   → x ← x + step * v
   → repeat
```

## Ranked table — pure numerical behaviour

| Rank | Experiment | Rule | Stable | Weak Lensing | Conservation residual | Total bending | Runtime (s) | Score | Comments |
|---:|---|---|:-:|:-:|---:|---:|---:|---:|:---|
| 1 | Exp 6 | magnetic-style directional | yes | yes | 1.40e-05 | 4.89e-06 | 0.008 | 3.001 | nominal |
| 2 | Exp 9 (m=0) | photon baseline (∇C, m=0) | yes | yes | 8.77e-05 | 2.43e-05 | 0.000 | 3.001 | nominal |
| 3 | Exp 5 | elastic spring restoring (∇C) | yes | yes | 8.77e-05 | 2.43e-05 | 0.008 | 3.001 | nominal |
| 4 | Exp 1 | response ∝ �C | yes | yes | 8.77e-05 | 2.43e-05 | 0.008 | 3.001 | nominal |
| 5 | Exp 4 | traction t = P_F · N | yes | yes | 8.77e-05 | 2.42e-05 | 0.009 | 3.001 | nominal |
| 6 | Exp 9 (m>0) | test particle (∇C, m=2) | yes | yes | 3.03e-04 | 8.40e-05 | 0.008 | 3.001 | nominal |
| 7 | Exp 2 | response ∝ ∇P / C | yes | yes | 7.44e-01 | 1.60e+00 | 0.009 | 2.574 | velocity drift |
| 8 | Exp 7 | phase accumulation (∇C·C) | yes | no | 2.56e-08 | 5.04e-09 | 0.008 | 2.001 | no bending observed |
| 9 | Exp 3 | response ∝ �W | yes | no | 3.00e-08 | 5.90e-09 | 0.008 | 2.001 | no bending observed |
| 10 | Exp 8 | constant local transfer | yes | no | 6.00e-04 | 6.00e-04 | 0.008 | 2.000 | no bending observed |

Scoring formula (purely numerical):
`score = stable_bonus + weak_lensing_bonus + 1/(1 + conservation_residual) + 0.001/(1 + runtime)`

`stable_bonus = 1 if (finite_outputs and |v|-drift < 1.0) else 0`
`weak_lensing_bonus = 1 if photon_max_deviation > 1e-6 else 0`

## Raw measurements

See `runs/transport_lab001/measurements.csv` and `runs/transport_lab001/measurements.json`.

## Per-experiment observations (reported, not interpreted)

- **Exp 1 (∇C)** — bending = 9.26e-06, conservation residual = 8.77e-05. Stable.
- **Exp 2 (∇P/C)** — bending = 2.07, conservation residual = 0.744. Numerically present but velocity drift large (1/C blows up where C → 0).
- **Exp 3 (∇W)** — bending = 5.73e-10. No observed bending above 1e-6 threshold.
- **Exp 4 (traction t = P_F·N)** — bending = 9.92e-06, conservation residual = 8.77e-05. Numerically equivalent to Exp 1 (traction field here is �C projected on radial normal; the radial projection preserves the gradient direction).
- **Exp 5 (elastic spring)** — bending = 9.26e-06, conservation residual = 8.77e-05. Identical numerical response to Exp 1 in this dataset (continuous-limit Hooke reduces to ∇C).
- **Exp 6 (magnetic-style)** — bending = 6.62e-05, conservation residual = 1.40e-05. Largest bending of all stable runs; perpendicular coupling gives ~7× the deflection of Exp 1.
- **Exp 7 (phase)** — bending = 5.74e-10. C·∇C is ~10× smaller than ∇C at this deformation strength; below 1e-6 threshold.
- **Exp 8 (constant)** — bending = 0. Uniform shift; all photons deflected by the same constant per step → identical deviation per ray (no spread). Conservation residual = 6e-04 from constant renormalisation.
- **Exp 9 (m=0 vs m>0, same �C rule)** — photon trajectory and test-particle trajectory coincide to within 5e-11 in final positions; bending differs by 5e-13. With the supplied local response scale and `m = 2`, inertia does not measurably change the trajectory in this configuration.

## What was observed, reported verbatim

> "When the following local neighbour rule is implemented, the following numerical behaviour is observed."

- **Exp 1**: response = ∇C → trajectory bends by 9.26e-06; stable; conservation residual 8.77e-05.
- **Exp 2**: response = ∇P/C → trajectory bends by 2.07; unstable (velocity drift 0.74).
- **Exp 3**: response = ∇W → trajectory bends by 5.73e-10; stable; below weak-lensing threshold.
- **Exp 4**: response = traction → trajectory bends by 9.92e-06; stable; conservation residual 8.77e-05.
- **Exp 5**: response = elastic spring → trajectory bends by 9.26e-06; stable; conservation residual 8.77e-05.
- **Exp 6**: response = magnetic-style → trajectory bends by 6.62e-05; stable; conservation residual 1.40e-05.
- **Exp 7**: response = phase accumulation → trajectory bends by 5.74e-10; stable; below weak-lensing threshold.
- **Exp 8**: response = constant transfer → trajectory bends by 0 (uniform shift); conservation residual 6e-04.
- **Exp 9**: response = ∇C, m=0 vs m=2 → trajectories coincide within numerical precision in this configuration.

Laboratory stops here. No experiment is declared physically correct.
