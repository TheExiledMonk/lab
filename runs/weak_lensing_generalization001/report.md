# PBUF WEAK-LENSING-GENERALIZATION-001

Frozen pipeline: matter → C = 0.18 · ρ/ρ_max → ∇C → |∇C| response → 90° transverse response → direct addition + renormalisation → photon propagation → observables.

## Identical pipeline check (SHA-256)

| File | SHA-256 |
|---|---|
| weak_lensing_generalization001.py | `7669454c09045569ae676a1bcb047bd7d749ae34efcb1362893520bd62c9d7bd` |
| constitutive_equations.py | `e2c789d19fd559753519704c6668c7a2879c53eb61a315604ec81af6795aca9f` |

## Cross-dataset summary

| Dataset | Runtime (s) | Max bend | RMS bend | Max κ | RMS κ | Max γ | RMS γ | Cons. |
|---|---|---|---|---|---|---|---|---|
| Isolated cluster | 0.007 | 6.618e-05 | 3.126e-05 | 5.000e-01 | 5.000e-01 | 4.740e+00 | 5.442e-01 | 2.22e-16 |
| Binary cluster | 0.007 | 9.473e-03 | 4.323e-03 | 5.000e-01 | 5.000e-01 | 4.740e+00 | 5.442e-01 | 2.22e-16 |
| Elongated cluster | 0.007 | 4.376e-02 | 1.513e-02 | 5.000e-01 | 5.000e-01 | 4.740e+00 | 5.442e-01 | 2.22e-16 |
| Asymmetric cluster | 0.007 | 3.819e-07 | 1.744e-07 | 5.000e-01 | 5.000e-01 | 4.740e+00 | 5.442e-01 | 1.11e-16 |
| Sparse field | 0.007 | 9.865e-02 | 8.092e-02 | 5.000e-01 | 5.000e-01 | 4.740e+00 | 5.442e-01 | 2.22e-16 |
| Dense field | 0.007 | 0.000e+00 | 0.000e+00 | 5.000e-01 | 5.000e-01 | 4.740e+00 | 5.442e-01 | 0.00e+00 |

## Per-dataset products

### Isolated cluster (`isolated`)

- CSVs: `matter.csv`, `constitutive.csv`, `gradient_*.csv`, `response_*.csv`, `convergence.csv`, `shear_g1.csv`, `shear_g2.csv`, `deflection_*.csv`, `magnification.csv`, `photon_*.csv`
- Maps: `matter_map.png`, `constitutive_map.png`, `gradient_map.png`, `response_magnitude_map.png`, `response_direction_map.png`, `convergence_map.png`, `shear_g1_map.png`, `shear_g2_map.png`, `shear_magnitude_map.png`, `deflection_magnitude_map.png`, `magnification_map.png`
- Trajectories: `photon_trajectories.png`; composite: `composite_observables.png`
- Statistics: `dataset_statistics.json`

### Binary cluster (`binary`)

- CSVs: `matter.csv`, `constitutive.csv`, `gradient_*.csv`, `response_*.csv`, `convergence.csv`, `shear_g1.csv`, `shear_g2.csv`, `deflection_*.csv`, `magnification.csv`, `photon_*.csv`
- Maps: `matter_map.png`, `constitutive_map.png`, `gradient_map.png`, `response_magnitude_map.png`, `response_direction_map.png`, `convergence_map.png`, `shear_g1_map.png`, `shear_g2_map.png`, `shear_magnitude_map.png`, `deflection_magnitude_map.png`, `magnification_map.png`
- Trajectories: `photon_trajectories.png`; composite: `composite_observables.png`
- Statistics: `dataset_statistics.json`

### Elongated cluster (`elongated`)

- CSVs: `matter.csv`, `constitutive.csv`, `gradient_*.csv`, `response_*.csv`, `convergence.csv`, `shear_g1.csv`, `shear_g2.csv`, `deflection_*.csv`, `magnification.csv`, `photon_*.csv`
- Maps: `matter_map.png`, `constitutive_map.png`, `gradient_map.png`, `response_magnitude_map.png`, `response_direction_map.png`, `convergence_map.png`, `shear_g1_map.png`, `shear_g2_map.png`, `shear_magnitude_map.png`, `deflection_magnitude_map.png`, `magnification_map.png`
- Trajectories: `photon_trajectories.png`; composite: `composite_observables.png`
- Statistics: `dataset_statistics.json`

### Asymmetric cluster (`asymmetric`)

- CSVs: `matter.csv`, `constitutive.csv`, `gradient_*.csv`, `response_*.csv`, `convergence.csv`, `shear_g1.csv`, `shear_g2.csv`, `deflection_*.csv`, `magnification.csv`, `photon_*.csv`
- Maps: `matter_map.png`, `constitutive_map.png`, `gradient_map.png`, `response_magnitude_map.png`, `response_direction_map.png`, `convergence_map.png`, `shear_g1_map.png`, `shear_g2_map.png`, `shear_magnitude_map.png`, `deflection_magnitude_map.png`, `magnification_map.png`
- Trajectories: `photon_trajectories.png`; composite: `composite_observables.png`
- Statistics: `dataset_statistics.json`

### Sparse field (`sparse`)

- CSVs: `matter.csv`, `constitutive.csv`, `gradient_*.csv`, `response_*.csv`, `convergence.csv`, `shear_g1.csv`, `shear_g2.csv`, `deflection_*.csv`, `magnification.csv`, `photon_*.csv`
- Maps: `matter_map.png`, `constitutive_map.png`, `gradient_map.png`, `response_magnitude_map.png`, `response_direction_map.png`, `convergence_map.png`, `shear_g1_map.png`, `shear_g2_map.png`, `shear_magnitude_map.png`, `deflection_magnitude_map.png`, `magnification_map.png`
- Trajectories: `photon_trajectories.png`; composite: `composite_observables.png`
- Statistics: `dataset_statistics.json`

### Dense field (`dense`)

- CSVs: `matter.csv`, `constitutive.csv`, `gradient_*.csv`, `response_*.csv`, `convergence.csv`, `shear_g1.csv`, `shear_g2.csv`, `deflection_*.csv`, `magnification.csv`, `photon_*.csv`
- Maps: `matter_map.png`, `constitutive_map.png`, `gradient_map.png`, `response_magnitude_map.png`, `response_direction_map.png`, `convergence_map.png`, `shear_g1_map.png`, `shear_g2_map.png`, `shear_magnitude_map.png`, `deflection_magnitude_map.png`, `magnification_map.png`
- Trajectories: `photon_trajectories.png`; composite: `composite_observables.png`
- Statistics: `dataset_statistics.json`

## Notes

All datasets were processed by the identical frozen Version A pipeline. No parameter was altered between datasets. No comparison with ΛCDM, GR, or observations was performed at this milestone.