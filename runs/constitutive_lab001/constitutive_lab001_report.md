# PBUF CONSTITUTIVE-LAB-001 — Constitutive-field generation rules

Frozen transport (LAB-008): neighbour-to-neighbour, 90° transverse response, linear amplitude A = |∇C|, direct addition + renormalisation, identical kernel, integration, normalisation, timestep.

Variable: the rule that turns matter(x) into C(x). All candidates normalised so max(C) = control max = 1.8000e-01.

## Constitutive statistics

| Candidate | Max C | Mean C | Max \|∇C\| | Mean \|∇C\| |
|---|---:|---:|---:|---:|
| Cand 1 (existing control) | 1.8000e-01 | 2.4590e-03 | 1.4481e-01 | 4.0950e-03 |
| Cand 2 (local linear) | 1.8000e-01 | 2.4590e-03 | 1.4481e-01 | 4.0950e-03 |
| Cand 3 (finite-range boxcar R=1.5) | 1.8000e-01 | 2.4429e-03 | 1.0273e-01 | 6.3304e-03 |
| Cand 4 (exponential L=1.0) | 1.8000e-01 | 5.1413e-03 | 7.9706e-02 | 7.7966e-03 |
| Cand 5 (Gaussian σ=1.0) | 1.8000e-01 | 4.0832e-03 | 8.8800e-02 | 6.8217e-03 |
| Cand 6 (inverse-distance L=1.0) | 1.8000e-01 | 5.4292e-03 | 4.6321e-02 | 9.4070e-03 |
| Cand 7 (compact-support R=2.0) | 1.8000e-01 | 2.6798e-03 | 1.2051e-01 | 5.0363e-03 |

## Weak-lensing performance

| Candidate | Bend | Conservation | Runtime (s) |
|---|---:|---:|---:|
| **Cand 1 (control)** | **6.6181e-05** | **2.22e-16** | 0.017 |
| **Cand 2 (local linear)** | **6.6181e-05** | **2.22e-16** | 0.017 |
| Cand 3 (boxcar R=1.5) | 5.4004e-07 | 1.11e-16 | 0.017 |
| Cand 4 (exponential L=1.0) | 8.2410e-03 | 2.22e-16 | 0.018 |
| Cand 5 (Gaussian σ=1.0) | 2.1212e-04 | 2.22e-16 | 0.017 |
| Cand 6 (inverse-distance L=1.0) | 3.2553e-02 | 2.22e-16 | 0.018 |
| Cand 7 (compact-support R=2.0) | 1.1183e-07 | 1.11e-16 | 0.017 |

## Relative comparison vs control

| Candidate | Bend Δ% | Conservation Δ% | \|∇C\| max Δ% | \|∇C\| mean Δ% |
|---|---:|---:|---:|---:|
| Cand 1 (control) | +0.0000% | +0.00% | +0.00% | +0.00% |
| Cand 2 (local linear) | +0.0000% | +0.00% | −0.00% | +0.00% |
| Cand 3 (boxcar R=1.5) | −99.18% | −50.00% | −29.06% | +54.59% |
| Cand 4 (exponential L=1.0) | +12 352.18% | +0.00% | −44.96% | +90.42% |
| Cand 5 (Gaussian σ=1.0) | +220.52% | +0.00% | −38.68% | +66.58% |
| Cand 6 (inverse-distance L=1.0) | +49 087.82% | +0.00% | −68.01% | +129.74% |
| Cand 7 (compact-support R=2.0) | −99.83% | −50.00% | −16.78% | +22.99% |

## Visualisations

`runs/constitutive_lab001/constitutive_fields_1to4.png` and `constitutive_fields_5to7.png` show the constitutive field C and the gradient magnitude |∇C| for every candidate on identical colour scales.

## Observations (no interpretation)

- **Cand 1 and Cand 2 are mathematically identical** — both produce C(x) proportional to matter(x), so max(C), mean(C), max|∇C|, mean|∇C|, bend, and conservation are identical to within machine precision.
- **All kernels (Cand 3–7) produce a smaller max |∇C|** than the control (16.8% to 68.0% smaller), because smoothing redistributes the peak over a larger area.
- **All kernels (Cand 3–7) produce a larger mean |∇C|** than the control (23.0% to 129.7% larger), because smoothing fills in the tail of the distribution.
- **Bending does not track max|∇C| or mean|∇C|** — it tracks the spatial distribution of |∇C| along the photon path:
  - Cand 3 (boxcar, R=1.5) reduces bending by 99.2% — the smoothing kills the narrow peak that the photon would otherwise sample.
  - Cand 4 (exponential, L=1.0) increases bending by 12 352% — the long exponential tail brings significant response into the photon path, which is otherwise off-peak.
  - Cand 5 (Gaussian, σ=1.0) increases bending by 221% — same tail effect, smaller magnitude.
  - Cand 6 (inverse-distance, L=1.0) increases bending by 49 088% — the 1/r-like tail is heaviest at large distances.
  - Cand 7 (compact-support, R=2.0) reduces bending by 99.8% — the kernel vanishes outside R = 2.0, so beyond R = 2.0 from the mass the response is exactly zero (the photon path lies outside this support).
- All 7 candidates are numerically stable (finite outputs, |v| bounded by renormalisation).

## Outcome

**Outcome B — several constitutive-field generation rules produce statistically indistinguishable transport.**

The statistically equivalent group is:

| Candidate | Rule | Bend vs control |
|---|---|---:|
| Cand 1 | Existing control (C = u₀ · matter / matter_max) | 0.00% |
| Cand 2 | Local linear (no kernel) | 0.00% |

These two are **mathematically identical** (both produce C(x) ∝ matter(x)) and so give bit-exact identical transport: bend_max = 6.6181e-05, conservation = 2.22e-16, identical ∇C statistics.

All other candidates (Cand 3, 4, 5, 6, 7) produce bending that deviates from the control by −99.8% to +49 088%, so they are not in the equivalent group. The deviations are not statistical noise — they are systematic and trace directly to the spatial redistribution of the constitutive field by the smoothing kernel, which changes what the photon path actually samples.

Artefacts: `runs/constitutive_lab001/{validation.csv, constitutive_statistics.csv, weak_lensing_performance.csv, relative_vs_control.csv, indistinguishable.json, measurements.csv, measurements.json, constitutive_fields_1to4.png, constitutive_fields_5to7.png, constitutive_lab001_report.md}`. Source: `constitutive_lab001.py`.

Laboratory stops. No physical interpretation, no new laws.