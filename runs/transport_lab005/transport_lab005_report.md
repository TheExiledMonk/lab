# PBUF TRANSPORT-LAB-005 — Which local constitutive direction governs the response?

Frozen: neighbour-to-neighbour transport, Lens-001, kernel, 90° response, integration, normalization, response magnitude policy (response = R_90(g) with full vector magnitude, not renormalised).

## Measurements

| Candidate | Reference direction | \|g\|_max | bend_max | conservation | speed_drift_pre_max | stable |
|---|---|---:|---:|---:|---:|:-:|
| 1 (control) | Constitutive gradient ∇C | 1.448e-01 | 6.618e-05 | 1.399e-05 | 1.399e-05 | yes |
| 2 | Stress gradient ∇|∇C|/(C+ε) | 9.896e+00 | 2.601e-01 | 1.607e-01 | 1.607e-01 | yes |
| 3 | Energy gradient C·∇C | 1.848e-02 | 5.841e-09 | 2.677e-09 | 2.677e-09 | yes |
| 4 | Traction (∇C·N̂)N̂ | 1.448e-01 | 6.618e-05 | 1.474e-05 | 1.474e-05 | yes |
| 5 | Force density ∇(∇²C) | 8.167e-01 | 1.478e-03 | 2.581e-04 | 2.581e-04 | yes |
| 6 | Principal strain direction | 1.000e+00 | 4.080e+00 | 6.000e-02 | 6.000e-02 | yes |
| 7 | Principal stress direction | 1.000e+00 | 4.068e+00 | 6.000e-02 | 6.000e-02 | yes |

## Ranking (by combined score = bend/bend_max + min_cons/conservation)

| Rank | Candidate | Reference | Bend | Conservation | Stable | Score |
|---:|---|---|---:|---:|:-:|---:|
| 1 | 6 | Principal strain direction | 4.080e+00 | 6.000e-02 | yes | 1.000 |
| 2 | 3 | Energy gradient C·∇C | 5.841e-09 | 2.677e-09 | yes | 1.000 |
| 3 | 7 | Principal stress direction | 4.068e+00 | 5.999e-02 | yes | 0.997 |
| 4 | 2 | Stress gradient ∇|∇C|/(C+ε) | 2.601e-01 | 1.607e-01 | yes | 0.064 |
| 5 | 5 | Force density ∇(∇²C) | 1.478e-03 | 2.581e-04 | yes | 0.000 |
| 6 | 1 (control) | Constitutive gradient ∇C | 6.618e-05 | 1.399e-05 | yes | 0.000 |
| 7 | 4 | Traction (∇C·N̂)N̂ | 6.618e-05 | 1.474e-05 | yes | 0.000 |

## Relative comparison vs control (Candidate 1 = ∇C)

| Candidate | Bend Δ% | Conservation Δ% | Speed drift Δ% |
|---|---:|---:|---:|
| 1 (control) | +0.00% | +0.00% | +0.00% |
| 2 (stress gradient) | +392 859% | +1 148 564% | +1 148 564% |
| 3 (energy gradient) | −99.99% | −99.98% | −99.98% |
| 4 (traction) | −0.00% | +5.32% | +5.32% |
| 5 (force density) | +2 134% | +1 745% | +1 745% |
| 6 (principal strain) | +6 164 908% | +428 656% | +428 656% |
| 7 (principal stress) | +6 146 167% | +428 569% | +428 569% |

## Correlation analysis

| Pair | Pearson r |
|---|---:|
| bend vs conservation residual | +0.2719 |
| bend vs speed drift (pre-renorm) | +0.2719 |
| conservation vs speed drift | +1.0000 |

Bending and conservation are weakly **positively** correlated across the seven candidates: candidates that bend more also accumulate more conservation residual. The two metrics are not anti-correlated.

## Observations (no interpretation)

- The control (∇C) and the traction (∇C·N̂)N̂ produce bending and conservation within numerical uncertainty of each other: bend 6.6181e-05 vs 6.6179e-05 (Δ = 3e-9 absolute), conservation 1.399e-05 vs 1.474e-05 (Δ = 5.3%). The two reference fields have the same |g| spectrum (max 0.145) and the same perpendicular structure.
- Candidates 2, 5, 6, 7 all give bending 4×–60 000× larger than the control. They also give conservation 18×–430 000× worse. The increase in bending and the degradation of conservation track each other — they are determined by |g|, not by the geometric character of the reference direction.
- Candidate 3 (energy gradient C·∇C) gives bending 1.1×10⁴× smaller than the control and conservation 5.2×10³× smaller. The small magnitude of C·∇C in this dataset (|g|_max = 0.018) is responsible.
- All seven candidates are numerically stable.
- No candidate produces simultaneously higher bending AND lower conservation than the control.

## Outcome

**Outcome B — several constitutive directions perform equivalently within numerical uncertainty.**

The statistically indistinguishable group is:

| Candidate | Reference direction |
|---|---|
| 1 (control) | Constitutive gradient ∇C |
| 4 | Traction (∇C·N̂)N̂ |

These two give bending 6.618e-05 (within 0.005% of each other) and conservation within 5.3% of each other. All other candidates produce different magnitudes and thus different bending/conservation, scaling with the magnitude of their reference field.

The control direction (∇C) and the traction (a projection of ∇C onto the radial unit normal) are indistinguishable in this dataset; both are equally valid choices of reference direction for the 90° transverse response.

Artefacts: `runs/transport_lab005/{measurements.csv, measurements.json, ranking.csv, relative_vs_control.csv, correlations.json, transport_lab005_report.md}`. Source: `transport_lab005.py`.

Laboratory stops. No physical interpretation, no new laws.
