# PBUF TRANSPORT-LAB-004 — Continuous angular sweep

The local response is parameterised as

```
r(θ) = cos(θ) · ĝ + sin(θ) · ĝ_⊥
```

with response magnitude ‖∇C‖ and the same kernel as TRANSPORT-LAB-003. Only the angle θ varies.

## Required table

### Coarse sweep (0° – 180°, step 15°)

| Angle (deg) | Bend | Conservation | Speed Drift | Stable |
|---:|---:|---:|---:|:-:|
|   0 | 9.26e-06 | 8.77e-05 | 8.77e-05 | yes |
|  15 | 1.84e-05 | 8.42e-05 | 8.42e-05 | yes |
|  30 | 3.42e-05 | 7.50e-05 | 7.50e-05 | yes |
|  45 | 4.77e-05 | 6.06e-05 | 6.06e-05 | yes |
|  60 | 5.80e-05 | 4.21e-05 | 4.21e-05 | yes |
|  75 | 6.43e-05 | 2.84e-05 | 2.84e-05 | yes |
| **90** | **6.62e-05** | **1.40e-05** | **1.40e-05** | yes |
| 105 | 6.36e-05 | 2.61e-05 | 2.61e-05 | yes |
| 120 | 5.67e-05 | 4.55e-05 | 4.55e-05 | yes |
| 135 | 4.59e-05 | 6.34e-05 | 6.34e-05 | yes |
| 150 | 3.20e-05 | 7.69e-05 | 7.69e-05 | yes |
| 165 | 2.02e-05 | 8.52e-05 | 8.52e-05 | yes |
| 180 | 9.26e-06 | 8.77e-05 | 8.77e-05 | yes |

### Fine sweep (80° – 100°, step 1°)

| Angle (deg) | Bend | Conservation | Speed Drift | Stable |
|---:|---:|---:|---:|:-:|
|  80 | 6.540e-05 | 2.375e-05 | 2.375e-05 | yes |
|  81 | 6.557e-05 | 2.280e-05 | 2.280e-05 | yes |
|  82 | 6.572e-05 | 2.184e-05 | 2.184e-05 | yes |
|  83 | 6.584e-05 | 2.088e-05 | 2.088e-05 | yes |
|  84 | 6.595e-05 | 1.991e-05 | 1.991e-05 | yes |
|  85 | 6.604e-05 | 1.894e-05 | 1.894e-05 | yes |
|  86 | 6.611e-05 | 1.796e-05 | 1.796e-05 | yes |
|  87 | 6.616e-05 | 1.697e-05 | 1.697e-05 | yes |
|  88 | 6.618e-05 | 1.598e-05 | 1.598e-05 | yes |
|  89 | 6.619e-05 | 1.499e-05 | 1.499e-05 | yes |
| **90** | **6.618e-05** | **1.399e-05** | **1.399e-05** | yes |
|  91 | 6.615e-05 | 1.484e-05 | 1.484e-05 | yes |
|  92 | 6.610e-05 | 1.568e-05 | 1.568e-05 | yes |
|  93 | 6.602e-05 | 1.651e-05 | 1.651e-05 | yes |
|  94 | 6.593e-05 | 1.734e-05 | 1.734e-05 | yes |
|  95 | 6.582e-05 | 1.817e-05 | 1.817e-05 | yes |
|  96 | 6.568e-05 | 1.898e-05 | 1.898e-05 | yes |
|  97 | 6.553e-05 | 1.980e-05 | 1.980e-05 | yes |
|  98 | 6.536e-05 | 2.060e-05 | 2.060e-05 | yes |
|  99 | 6.516e-05 | 2.140e-05 | 2.140e-05 | yes |
| 100 | 6.495e-05 | 2.220e-05 | 2.220e-05 | yes |

## Plots

`runs/transport_lab004/angle_sweep_coarse.png` and `runs/transport_lab004/angle_sweep_fine.png` show:

- Bending: peaks at θ = 89° – 90° (fine), 90° (coarse). Profile tracks sin(θ) within the half-period [0°, 180°], with a small non-zero residual at θ = 0° and θ = 180° from the offset between v and ∇C along the photon path.
- Conservation: monotonically improves from θ = 0° to 90°, then monotonically degrades symmetrically to 180°. Sharp minimum at θ = 90°.
- Speed drift: numerically identical to conservation residual (they measure the same pre-normalization quantity here).
- Direction drift: tracks the bending profile (largest direction change at the angles that produce the largest bending).

## Optima (coarse, 5%-tolerance plateaus)

| Quantity | Value | Plateau within 5% |
|---|---:|---|
| **Maximum bending angle** | **θ = 90°** | 75° – 105° |
| **Conservation optimum** | **θ = 90°** | 90° – 90° (sharp) |
| **Combined optimum** (bending + conservation, normalised sum) | **θ = 90°** | — |

Fine sweep (80° – 100°, step 1°) confirms the location of the bending peak at θ = 89° – 90° (6.6193e-05 vs 6.6181e-05) and the conservation minimum exactly at θ = 90° (1.3994e-05).

## Numerical observations

- The bending curve is approximately proportional to sin(θ) in [0°, 180°]: small at θ = 0° (9.26e-06), rises monotonically to 6.62e-05 at 90°, falls symmetrically back to 9.26e-06 at 180°. The non-zero residual at θ = 0° (and θ = 180°) arises from the offset between the photon velocity and ∇C along the photon path.
- The conservation curve is approximately proportional to cos²(θ): maximum at θ = 0°/180° (8.77e-05), minimum at θ = 90° (1.40e-05), factor 6.3× ratio.
- All 34 measurements (13 coarse + 21 fine) were numerically stable.

## Outcome

**Outcome A — a unique optimum angle exists.**

| | |
|---|---:|
| Maximum bending angle | **θ = 90°** |
| Conservation optimum | **θ = 90°** |
| Combined optimum | **θ = 90°** |

The 90° response is the unique simultaneous optimum for both metrics. The bending has a broad plateau (75° – 105° within 5% of the maximum), but the conservation optimum is sharp at exactly 90°.

Laboratory stops. No physical interpretation, no new laws.
