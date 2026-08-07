# PBUF TRANSPORT-LAB-008 — Local response amplitude law ablation

Frozen: response direction (90° transverse), update rule (direct addition + renormalisation), kernel, integration, timestep, normalisation, Lens-001.

Variable: the scalar amplitude law A(|∇C|). All candidates are normalised so `max(A) = max(|∇C|)`.

## Required validation table

| Candidate | Max A | Mean A | Pass |
|---|---:|---:|:-:|
| Cand 1 (linear) | 1.4481e-01 | 4.0950e-03 | PASS |
| Cand 2 (sqrt) | 1.4481e-01 | 6.5241e-03 | PASS |
| Cand 3 (quadratic) | 1.4481e-01 | 2.6929e-03 | PASS |
| Cand 4 (log) | 1.4481e-01 | 4.1907e-03 | PASS |
| Cand 5 (saturating) | 1.4481e-01 | 7.0365e-03 | PASS |
| Cand 7 (piecewise) | 1.4481e-01 | 5.1516e-03 | PASS |
| Cand 6.05 (threshold 5%) | 1.4481e-01 | 3.9867e-03 | PASS |
| Cand 6.10 (threshold 10%) | 1.4481e-01 | 3.8746e-03 | PASS |
| Cand 6.20 (threshold 20%) | 1.4481e-01 | 3.6370e-03 | PASS |

All 9 candidates pass: max(A) equals the control max (1.4481e-01) to machine precision.

## Required performance table

| Candidate | Bend | Conservation | Runtime (s) |
|---|---:|---:|---:|
| Cand 1 (linear, control) | 6.6181e-05 | 2.22e-16 | 0.017 |
| Cand 2 (sqrt, normalised) | 2.1580e-03 | 2.22e-16 | 0.017 |
| Cand 3 (quadratic, normalised) | 2.1530e-07 | 2.22e-16 | 0.018 |
| Cand 4 (log, normalised) | 7.0848e-05 | 2.22e-16 | 0.017 |
| Cand 5 (saturating, normalised) | 3.9081e-04 | 2.22e-16 | 0.017 |
| Cand 7 (piecewise, normalised) | 1.0182e-04 | 2.22e-16 | 0.016 |
| Cand 6.05 (threshold 5%) | 0.0000e+00 | 0.00e+00 | 0.017 |
| Cand 6.10 (threshold 10%) | 0.0000e+00 | 0.00e+00 | 0.017 |
| Cand 6.20 (threshold 20%) | 0.0000e+00 | 0.00e+00 | 0.017 |

## Response curves

`runs/transport_lab008/response_curves.png` shows A(|∇C|) for every candidate on the normalised scale.

## Relative comparison vs control

| Candidate | Bend Δ% | Conservation Δ% | Runtime Δ% |
|---|---:|---:|---:|
| Cand 1 (linear, control) | +0.0000% | +0.0000% | +0.0000% |
| Cand 2 (sqrt) | +3 160.68% | +0.0000% | +2.05% |
| Cand 3 (quadratic) | −99.67% | +0.0000% | +4.24% |
| Cand 4 (log) | +7.05% | +0.0000% | +1.93% |
| Cand 5 (saturating) | +490.51% | +0.0000% | −1.24% |
| Cand 7 (piecewise) | +53.85% | +0.0000% | −2.65% |
| Cand 6.05 (threshold 5%) | −100.00% | −100.00% | −0.71% |
| Cand 6.10 (threshold 10%) | −100.00% | −100.00% | −1.33% |
| Cand 6.20 (threshold 20%) | −100.00% | −100.00% | −1.17% |

## Observations (no interpretation)

- **Cand 4 (logarithmic) reproduces the control within 7%**: bend = 7.08e-05 vs 6.62e-05, conservation identical. Reason: along the photon path |∇C| ∈ [~1e-21, 2.5e-3], well within the regime where log(1 + x) ≈ x; after peak-normalisation the two amplitude laws differ only by O(|∇C|²/2) corrections.
- **Cand 7 (piecewise, break at 50%) gives +54% more bending** (1.02e-04 vs 6.62e-05). Most of the field is below the breakpoint, so the piecewise law is locally linear almost everywhere; the larger mean amplitude (5.15e-3 vs 4.10e-3) reflects the slope reduction above the break and the renormalisation.
- **Cand 5 (saturating) gives +491% bending** (3.91e-04). The renormalised saturation curve has a higher mean amplitude (7.04e-3) than the linear law; the small-|∇C| regime is approximately linear, but the normalisation forces a steeper slope at small |∇C| to compensate for the saturation plateau at large |∇C|.
- **Cand 2 (sqrt, normalised) gives +3 161% bending** (2.16e-03). Square root compresses the dynamic range of |∇C|, so after renormalisation the small-|∇C| cells contribute much more.
- **Cand 3 (quadratic, normalised) gives −99.7% bending** (2.15e-07). Quadratic emphasises the peak and suppresses everything else; after renormalisation the off-peak cells contribute negligibly.
- **Cand 6 (threshold at 5%, 10%, 20%) gives zero bending and zero conservation residual**. Reason: the photon path (x ∈ [−8, −3.2], y ∈ [−3, 3]) samples |∇C| ∈ [~0, 2.5e-3] at the closest point to the mass, which is below 5% of the peak (7.24e-3). With A = 0 along the entire path, the photon travels undeflected.
- All 9 candidates are numerically stable (finite outputs, |v| bounded by renormalisation).

## Outcome

**Outcome B — several amplitude laws perform equivalently.**

The statistically equivalent group is:

| Candidate | Amplitude law | Bend Δ vs control |
|---|---|---:|
| Cand 1 (linear) | A = \|∇C\| | 0.00% |
| Cand 4 (log) | A = log(1 + \|∇C\|) | +7.05% |

These two candidates give bending within 7% of each other and identical machine-precision conservation. Both behave as approximately linear for the small-|∇C| values actually sampled by the photon path; the logarithmic law adds only a small nonlinear correction that the propagation cannot resolve at the current step size.

All other candidates (sqrt, quadratic, saturating, piecewise, threshold) produce bending that deviates from the control by −100% to +3 161%, so they are not in the equivalent group.

Artefacts: `runs/transport_lab008/{validation.csv, performance_table.csv, relative_vs_control.csv, response_curves.png, indistinguishable.json, measurements.csv, measurements.json, transport_lab008_report.md}`. Source: `transport_lab008.py`.

Laboratory stops. No physical interpretation, no new laws.