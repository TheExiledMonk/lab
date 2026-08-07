# PBUF TRANSPORT-LAB-006 — Pure direction ablation (magnitude-normalised)

Critical fix vs LAB-005: every candidate produces a response of identical magnitude A = |∇C|. Only the response direction varies.

## Validation (pre-propagation)

| Candidate | Mean |r − A| | Max |r − A| | Pass/Fail |
|---|---:|---:|:-:|
| 1 (∇C) | 2.78e-17 | 2.78e-17 | **PASS** |
| 2 (∇P/C) | 2.78e-17 | 2.78e-17 | **PASS** |
| 3 (∇W) | 2.78e-17 | 2.78e-17 | **PASS** |
| 4 (traction) | 2.78e-17 | 2.78e-17 | **PASS** |
| 5 (Div P_F) | 2.78e-17 | 2.78e-17 | **PASS** |
| 6 (principal strain) | 2.78e-17 | 2.78e-17 | **PASS** |
| 7 (principal stress) | 2.78e-17 | 2.78e-17 | **PASS** |

All seven candidates pass: |r| = A = |∇C| at every cell to machine precision.

## Performance table

| Candidate | Bend | Conservation | Speed Drift | Stable |
|---|---:|---:|---:|:-:|
| 1 (∇C) | 6.6181e-05 | 1.3994e-05 | 1.3994e-05 | yes |
| 2 (∇P/C) | 6.6180e-05 | 1.4357e-05 | 1.4357e-05 | yes |
| 3 (∇W) | 6.6181e-05 | 1.3994e-05 | 1.3994e-05 | yes |
| 4 (traction) | 6.6179e-05 | 1.4759e-05 | 1.4759e-05 | yes |
| 5 (Div P_F) | 6.6180e-05 | 1.4452e-05 | 1.4452e-05 | yes |
| 6 (principal strain) | 6.6178e-05 | 1.4972e-05 | 1.4972e-05 | yes |
| 7 (principal stress) | 6.6175e-05 | 1.5895e-05 | 1.5895e-05 | yes |

## Relative comparison vs control (Candidate 1 = ∇C)

| Candidate | Bend Δ% | Conservation Δ% | Speed drift Δ% |
|---|---:|---:|---:|
| 1 (control) | +0.0000% | +0.0000% | +0.0000% |
| 2 (∇P/C) | −0.0016% | +2.5953% | +2.5953% |
| 3 (∇W) | +0.0000% | +0.0000% | +0.0000% |
| 4 (traction) | −0.0030% | +5.4646% | +5.4646% |
| 5 (Div P_F) | −0.0015% | +3.2725% | +3.2725% |
| 6 (principal strain) | −0.0045% | +6.9888% | +6.9888% |
| 7 (principal stress) | −0.0091% | +13.5871% | +13.5871% |

## Indistinguishable group (within 5% on BOTH bend and conservation)

| Candidate | Bend Δ% | Conservation Δ% |
|---|---:|---:|
| 1 (∇C) | 0.00% | 0.00% |
| 2 (∇P/C) | −0.0016% | +2.60% |
| 3 (∇W) | 0.00% | 0.00% |
| 5 (Div P_F) | −0.0015% | +3.27% |

(Not in the 5% group: 4, 6, 7 — they exceed 5% on conservation.)

## Observations (no interpretation)

- All seven candidates produce bending within 0.01% of each other (range 6.6175e-05 to 6.6181e-05). With magnitude fixed, the response direction has negligible effect on bending.
- Conservation residuals vary across candidates from 1.40e-05 to 1.59e-05 — a 13.6% spread, but all within an order of magnitude.
- Candidate 1 (∇C, control) and Candidate 3 (∇W = C·∇C) give **identical** numerical results (6.6181e-05 / 1.3994e-05) because C·∇C is a positive scalar multiple of ∇C in this dataset (C ≥ 0 everywhere), so their unit directions are equal.
- Candidate 7 (principal stress direction) is the worst on conservation (1.59e-05) — 13.6% above the control. The principal eigenvector differs slightly from the gradient direction where the Hessian is non-negligible relative to the gradient-magnitude term.
- All seven are numerically stable.

## Outcome

**Outcome C — the apparent superiority observed in LAB-005 disappears completely after magnitude normalization.**

Comparison of LAB-005 (raw magnitude) vs LAB-006 (magnitude-normalised) for each candidate:

| Candidate | LAB-005 bend | LAB-006 bend | LAB-005 cons | LAB-006 cons |
|---|---:|---:|---:|---:|
| 1 (∇C) | 6.62e-05 | 6.62e-05 | 1.40e-05 | 1.40e-05 |
| 2 (∇P/C) | 2.60e-01 | 6.62e-05 | 1.61e-01 | 1.44e-05 |
| 3 (∇W) | 5.84e-09 | 6.62e-05 | 2.68e-09 | 1.40e-05 |
| 4 (traction) | 6.62e-05 | 6.62e-05 | 1.47e-05 | 1.48e-05 |
| 5 (Div P_F) | 1.48e-03 | 6.62e-05 | 2.58e-04 | 1.45e-05 |
| 6 (principal strain) | 4.08e+00 | 6.62e-05 | 6.00e-02 | 1.50e-05 |
| 7 (principal stress) | 4.07e+00 | 6.62e-05 | 6.00e-02 | 1.59e-05 |

In LAB-005, bending ranged over 9 orders of magnitude (5.84e-09 to 4.08e+00) and conservation over 5 orders of magnitude (2.68e-09 to 1.61e-01). In LAB-006, bending is constant to within 0.01% and conservation varies by only 13.6%. The reference direction has effectively no influence on the transport once the response magnitude is fixed.

Artefacts: `runs/transport_lab006/{validation.csv, performance_table.csv, relative_vs_control.csv, indistinguishable.json, status.json, measurements.csv, measurements.json, transport_lab006_report.md}`. Source: `transport_lab006.py`.

Laboratory stops. No physical interpretation, no new laws.
