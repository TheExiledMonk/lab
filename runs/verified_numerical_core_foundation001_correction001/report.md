# PBUF FOUNDATION-001-CORRECTION-001 — Verified Numerical Core

**Correction ID**: FOUNDATION-001-CORRECTION-001
**Conventions version**: 1.1.0-correction001
**Duration**: 0.3s

## Status reset (§2)

- Experimental modules: M04_tensor_transforms, M08_pair_transfer, M10_midpoint_rasterization, M12_differential_operators, M13_helmholtz_3d, M16_observable_extraction
- Pending contract review: M01_conventions, M02_coordinate_transforms, M03_vector_transforms, M05_pair_enumeration, M11_field_diagnostics, M14_los_projection, M15_ray_interface

## Module results

- `M01_conventions`: unit_verified_pending_contract_review, max_error=0.000e+00, 19/19 tests pass
- `M02_coordinate_transforms`: unit_verified_pending_contract_review, max_error=0.000e+00, 36/36 tests pass
- `M03_vector_transforms`: unit_verified_pending_contract_review, max_error=0.000e+00, 42/42 tests pass
- `M04_tensor_transforms`: experimental, max_error=0.000e+00, 42/42 tests pass
- `M05_pair_enumeration`: unit_verified_pending_contract_review, max_error=0.000e+00, 39/39 tests pass
- `M08_pair_transfer`: experimental, max_error=0.000e+00, 3/3 tests pass
- `M10_midpoint_rasterization`: experimental, max_error=4.988e+00, 5/5 tests pass
- `M11_field_diagnostics`: unit_verified_pending_contract_review, max_error=0.000e+00, 4/4 tests pass
- `M12_differential_operators`: experimental, max_error=4.441e-16, 9/9 tests pass
- `M13_helmholtz_3d`: experimental, max_error=0.000e+00, 8/8 tests pass
- `M14_los_projection`: unit_verified_pending_contract_review, max_error=0.000e+00, 8/8 tests pass
- `M15_ray_interface`: unit_verified_pending_contract_review, max_error=0.000e+00, 6/6 tests pass
- `M16_observable_extraction`: experimental, max_error=0.000e+00, 15/15 tests pass

## Wrong controls (§19)

- `WC1 old-rasterization-omits-final-pair`: passes=True
- `WC2 wrong-midpoint-geometry`: passes=True
- `WC3 wrong-curl-reference`: passes=True
- `WC4 duplicate-helmholtz-marked`: passes=True
- `WC5 old-tied-rank-vs-new`: passes=True
- `WC6 structured-small-acceptance`: passes=True

## Stage R1 — synthetic integration

Shape: (4, 5, 6), n_pairs: 286, passes: True

- Endpoint closure (random): 4.070e-15
- Endpoint energy (random): 5.720e+02
- Interface energy (random): 1.316e+02
- Interface RMS (random): 1.047e+00

## Stage R2 — MACS0416 restricted recovery

Cluster: MACS0416, candidate: PL1_PM1_PS2, shape: (9, 32, 32), passes: True

- Endpoint energy: 4.364e-04
- Endpoint closure norm: 4.064e-17
- Interface energy: 2.952e-04
- LOS Rx RMS: 3.812e-05
- LOS Ry RMS: 3.868e-05
- Central Rx RMS: 1.505e-05
- Central Ry RMS: 1.357e-05
- Ray input classification: structured_normal
- Helmholtz f_irr_partition: 0.365207804627653
- Helmholtz f_sol_partition: 0.6347921953723471
- Helmholtz f_irr_native: 0.36511805338202197
- Helmholtz f_sol_native: 0.6346361927088509

## Stage R3 — covariance confirmation

All RC0..RC6 within 0.05: True

- RC0: E_cov_correct=0.000e+00, E_cov_wrong=0.000e+00, passes=True
- RC1: E_cov_correct=0.000e+00, E_cov_wrong=1.363e+00, passes=True
- RC2: E_cov_correct=0.000e+00, E_cov_wrong=1.093e+00, passes=True
- RC3: E_cov_correct=0.000e+00, E_cov_wrong=1.002e+00, passes=True
- RC4: E_cov_correct=0.000e+00, E_cov_wrong=9.845e-01, passes=True
- RC5: E_cov_correct=0.000e+00, E_cov_wrong=1.083e+00, passes=True
- RC6: E_cov_correct=0.000e+00, E_cov_wrong=1.363e+00, passes=True

## Outcome determination

Outcome A (all modules requalified): CORRECTION-001 stabilises the verified core. Stage R1, R2, and R3 all pass; covariance remains restored; the restricted recovery remains nontrivial. Full five-cluster restricted physics rerun MAY proceed once a second review accepts this correction output (review_status = pending_second_review at the time of writing).

Physics results from the previous FOUNDATION-001 run are PROVISIONAL and require requalification against the corrected core before downstream interpretation.

## Final report questions (§22)

1. The previous rasterization used `[:-2]` for the valid source slice, omitting the LAST valid internal pair adjacent to the upper boundary. For shape (4, 5, 6) this omitted 20 + 24 + 30 = 74 pairs in total (4 per axis-slice).
2. Yes. The omission materially altered interface energy because the boundary-adjacent pairs were skipped; the corrected closure identity (sum_i R_interface = sum R_ij) no longer holds for the predecessor.
3. Yes. Pair midpoints now use the generic `0.5*(i+j)` rule and pass M05-C1 (exact identity), C2 (fixed-axis coordinates), C3 (direction displacement) and C4 (noncubic grid).
4. Yes. The corrected reference curl agrees with production to 1e-12 on three independent curl fixtures.
5. Yes. Curl fixtures 1, 2, 3 exercise Cx, Cy, Cz independently.
6. Algebraically equivalent BEFORE magnitude normalisation; distinct AFTER PM1 on a spatially varying projector. The candidate registry marks them as distinct after PM1.
7. Yes. PS1-A is restricted to diagnostic-only (raw single-endpoint, no antisymmetry). pair_antisymmetry_expected = false; physics_candidate = false.
8. Yes. All RC0..RC6 inverse round-trips enforced; the previous waiver has been removed and a tensor-transform einsum bug fixed.
9. Yes. Helmholtz now uses physical spacing (K = 2π·fftfreq(n,d)).
10. Yes. Padded and cropped closures are reported separately.
11. Yes. Pure-longitudinal, pure-transverse, and mixed-mode analytic Fourier fixtures recover exact fractions (f_irr=1, f_sol=1, f_irr=0.2 respectively).
12. Yes. LOS projection exposes `project_vector_to_image_plane` with metadata (los_axis, depth_array_axis, image_component_1, image_component_2, output_plane_axis_order).
13. Yes. Spearman uses scipy.stats.rankdata(method='average') semantics; tied plateau produces different results from double-argsort.
14. Yes. `extract_jacobian_observables` accepts an optional `reference_kappa` argument; when None, no GR correlation is computed.
15. Yes. `array_fingerprint` exposes both `raw_sha256` (dtype+shape+raw bytes) and `canonical_float64_sha256` (content only).
16. Yes. The corrected ray interface accepts structured_small fields (variance below 1e-12 but spatial variation exists).
17. Yes. Endpoint energy and interface energy remain > 0 on MACS0416 with the corrected core.
18. Yes. RC0..RC6 round-trips restore E_cov < 0.05.
19. All values affected by the interface rasterization off-by-one correction (interface energy, RMS, central RMS, LOS projection) changed. The Helmholtz fractions are also affected by the spacing/padding correction.
20. Zero modules are now legitimately frozen. All thirteen modules either remain `experimental` (6) or `unit_verified_pending_contract_review` (7) until the second review is performed.