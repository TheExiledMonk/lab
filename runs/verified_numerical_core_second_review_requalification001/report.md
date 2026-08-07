# PBUF VERIFIED NUMERICAL CORE - SECOND-REVIEW-REQUALIFICATION-001

**Lab ID**: PBUF-SECOND-REVIEW-REQUALIFICATION-001
**Conventions version**: 1.1.0-correction001
**Head SHA**: d0c7763cea3b3ec568d79459431656b5378ca2b3
**Branch**: main
**Duration**: 0.4s

## Repository state (sec 2)

- branch: `main`
- head_sha: `d0c7763cea3b3ec568d79459431656b5378ca2b3`
- working_tree_clean: `True`
- required_prs_present: `[2, 3, 4, 5, 6, 7]`
- merge SHAs:
  - #1: `aa6ddd21d0d392d0d7a7ce662123cae444ee0bd3`
  - #2: `711185a4ab09e14fc01e993ae9e86a9f20b6f81f`
  - #3: `c3be6407fb3cb6c302eafa3eb28f25d005fb039e`
  - #4: `a18113e038d45a7228c5809c65bdc9b1c59e322f`
  - #5: `d5fbedfd825b7df64ceeefb5a1f29db1a47c9b37`
  - #6: `2ab8f4a2cf151b33d419fc03f4eb443b9dec245d`
  - #7: `d0c7763cea3b3ec568d79459431656b5378ca2b3`

## Source integrity & protected functions (sec 5)

- protected_function_violations: `0`

## Module requalification (sec 6)

- `M01`: status=PASS, max_error=0.000e+00, 20/20 tests pass
- `M02`: status=PASS, max_error=0.000e+00, 36/36 tests pass
- `M03`: status=PASS, max_error=0.000e+00, 42/42 tests pass
- `M04`: status=PASS, max_error=0.000e+00, 42/42 tests pass
- `M05`: status=PASS, max_error=0.000e+00, 36/36 tests pass
- `M06`: status=PASS, max_error=0.000e+00, 2/2 tests pass
- `M07`: status=PASS, max_error=0.000e+00, 2/2 tests pass
- `M08` PS-contract: 7/7 pass, overall=True
- `M09` endpoint closure: E_endpoint=4.205e+02, closure_norm=2.176e-15, passes=True
- `M10` interface rasterisation: audit_pass=True, impulse_pass=True, omitted=0, duplicated=0
- `M11`: 6/6 pass, overall=True
- `M12`: status=PASS, max_error=4.441e-16, 9/9 tests pass
- `M13`: 7/7 pass, overall=True
- `M14`: 3/3 pass, overall=True
- `M15`: 6/6 pass, overall=True
- `M16`: 10/10 pass, overall=True

## M03 independent closed-form validation (sec 7)

- V1 round-trip: `True`
- V2 explicit component mapping: `True`
- V3 wrong scalar-only inverse: `True`
- M03_PASS: `True`

## M08 PS-lane contract (sec 8)

- R_PS1 == R_PS1-B: `True`
- R_PS1 != R_PS2: `True`
- no_duplicate_physics_candidates: `True`

## M09 endpoint closure (sec 9)

- E_endpoint > 0: `True` (E_endpoint=4.205e+02)
- |sum_i R_i| <= 1e-12*max(1, sqrt(E)): `True` (closure_norm=2.176e-15)

## M10 interface rasterisation (sec 10)

- axis=xp: expected=100, actual=100, omitted=0, duplicated=0, passes=True
- axis=yp: expected=96, actual=96, omitted=0, duplicated=0, passes=True
- axis=zp: expected=90, actual=90, omitted=0, duplicated=0, passes=True
- axis=TOTAL: expected=286, actual=286, omitted=0, duplicated=0, passes=True
- impulse_axis=xp: src=4, dst=5, src_value=0.500, dst_value=0.500, elsewhere_max=0.000e+00, passes=True
- impulse_axis=yp: src=3, dst=4, src_value=0.500, dst_value=0.500, elsewhere_max=0.000e+00, passes=True
- impulse_axis=zp: src=2, dst=3, src_value=0.500, dst_value=0.500, elsewhere_max=0.000e+00, passes=True

## M11 diagnostics (sec 11)

- D1_finite_field: passes=True
- D2_NaN: passes=True
- D3_Inf: passes=True
- D4_zero_allow: passes=True
- D4_NaN_still_fails: passes=True
- D5_dual_hashes: passes=True

## M13 Helmholtz requalification (sec 12-13)

- metric_set_padded: passes=True
- metric_set_cropped: passes=True
- H1_pure_longitudinal: passes=True
- H2_pure_transverse: passes=True
- H3_equal_energy_mixed: passes=True
- H4_reconstruction: passes=True
- metric_labels_distinct: passes=True

- **padded**: field_reconstruction_error=3.897e-17, energy_closure_error=4.228e-02, orthogonality_error=2.114e-02, f_irr_partition=0.322142, f_sol_partition=0.677858
- **cropped/native**: field_reconstruction_error=3.905e-17, energy_closure_error=1.722e-02, orthogonality_error=8.610e-03, f_irr_partition=0.328859, f_sol_partition=0.671141

## M14 LOS projection (sec 14)

- los_axis=z: image_component_1=x, image_component_2=y, passes=True
- los_axis=y: image_component_1=x, image_component_2=z, passes=True
- los_axis=x: image_component_1=y, image_component_2=z, passes=True

## M15 ray interface (sec 15)

- exact_zero_require_nontrivial_rejected: passes=True
- exact_zero_require_nontrivial_false_allowed: passes=True
- constant_nonzero_allowed: passes=True
- structured_small_accepted: passes=True
- structured_normal_accepted: passes=True
- nonfinite_rejected_with_trivial_allow: passes=True

## M16 observable extraction (sec 16)

- O1_pearson_perfect_positive: passes=True
- O1_pearson_constant_returns_nan: passes=True
- O2_spearman_tied_perfect_positive: passes=True
- O2_spearman_vs_scipy: passes=True
- O3_nan_masked_monotonic: passes=True
- O4_old_tied_rank_disagrees: passes=True
- O5_pearson_shape_rejected: passes=True
- O5_spearman_shape_rejected: passes=True
- O6_no_reference_keys_absent: passes=True
- O6_with_reference_keys_finite: passes=True

## R1 synthetic integration (sec 17)

- shape: (4, 5, 6)
- n_pairs: 286
- endpoint_energy_random: 5.764844e+02
- interface_energy_random: 1.368836e+02
- endpoint_closure_norm_random: 4.529e-15
- passes: True

## Zero-field full-chain control (sec 18)

- R_pair_zero: passes=True
- R_endpoint_zero: passes=True
- R_interface_zero: passes=True
- kappa_zero: passes=True
- gamma1_zero: passes=True
- gamma2_zero: passes=True
- ray_classification_exact_zero: passes=True
- pearson_constant_returns_nan: passes=True
- no_pearson_vs_reference_key_when_no_ref: passes=True
- overall: passes=True

## A/zero/B state-retention control (sec 19)

- B_after_zero_hash == B_fresh_hash: `True`

## R2 MACS0416 restricted recovery (sec 20-22)

- cluster_id: `MACS0416`
- candidate_id: `PL1_PM1_PS2`
- transform_id: `RC0`
- shape: (9, 32, 32)
- n_pairs: 26048
- pair_amplitude_rms: 1.307e-04
- pair_response_rms: 1.157e-05
- endpoint_energy: 4.364493e-04
- endpoint_closure: 3.145e-17
- interface_energy: 2.952308e-04
- interface_global_sum: [-0.002418431691133661, -0.00449482349836508, 0.0002642349736815819]
- central_rx_rms: 1.504589e-05
- los_rx_rms: 3.812069e-05
- ray_classification: `structured_normal`
- ray_rms: 5.430694e-05
- kappa_variance: 3.115033e-09
- gamma_variance: 4.387556e-10
- los_metadata: {'los_axis': 'z', 'depth_array_axis': 0, 'image_component_1': 'x', 'image_component_2': 'y'}

Helmholtz (padding=none, native):
  - field_reconstruction_error: 2.0058641175033196e-17
  - energy_closure_error: 0.00024575390912704527
  - orthogonality_error: 0.00012287695456352968
  - f_irr_partition: 0.365207804627653
  - f_sol_partition: 0.634792195372347
  - f_irr_native: 0.3651180533820221
  - f_sol_native: 0.6346361927088509

Helmholtz (padding=reflect_half, padded):
  - field_reconstruction_error: 1.9041824162432556e-17
  - energy_closure_error: 0.052113198749628
  - orthogonality_error: 0.026056599374814023
  - f_irr_partition: 0.3737722824579984
  - f_sol_partition: 0.6262277175420017
  - f_irr_native: 0.3932507517008341
  - f_sol_native: 0.6588624470487939

Lineage:
  - endpoint_field: `0724e7adcff4423b195e799c6c4dea604f71a3091d38211a0f151b2e3c0bb783`
  - interface_field: `77a3ec918159117ae5ca0dc1721302e1cf87e78eb94d2cf0334a313b32825023`
  - los_field: `d108ebb979e03005bc2348d47cb9c0b7f3048f50945fdcd4710ea2067aee3100`
  - ray_input: `a110f5be9b0c9a9e929a7b2d9abd5d66b1ba19ed1d11353c01a713eea100a2d8`

- R2 nontriviality: E_endpoint>0: True, E_interface>0: True, Var(kappa)>0: True, ray not zero/nonfinite: True

## R3 covariance revalidation (sec 23)

- shape: (9, 32, 32)
- all_pass: True
- RC0: E_cov_correct=0.000e+00, E_cov_wrong=0.000e+00, passes=True
- RC1: E_cov_correct=0.000e+00, E_cov_wrong=1.363e+00, passes=True
- RC2: E_cov_correct=0.000e+00, E_cov_wrong=1.093e+00, passes=True
- RC3: E_cov_correct=0.000e+00, E_cov_wrong=1.002e+00, passes=True
- RC4: E_cov_correct=0.000e+00, E_cov_wrong=9.845e-01, passes=True
- RC5: E_cov_correct=0.000e+00, E_cov_wrong=1.083e+00, passes=True
- RC6: E_cov_correct=0.000e+00, E_cov_wrong=1.363e+00, passes=True

## Historical comparison (sec 25)

- endpoint_energy: previous=4.364000e-04, new=4.364493e-04, abs_change=4.932650e-08, rel_change=1.130e-04
- interface_energy: previous=2.952000e-04, new=2.952308e-04, abs_change=3.075037e-08, rel_change=1.042e-04
- central_rx_rms: previous=1.505000e-05, new=1.504589e-05, abs_change=-4.107046e-09, rel_change=-2.729e-04
- los_rx_rms: previous=3.812000e-05, new=3.812069e-05, abs_change=6.854482e-10, rel_change=1.798e-05
- f_irr_partition: previous=3.652078e-01, new=3.652078e-01, abs_change=0.000000e+00, rel_change=0.000e+00
- f_sol_partition: previous=6.347922e-01, new=6.347922e-01, abs_change=-1.110223e-16, rel_change=-1.749e-16
- field_reconstruction_error: previous=None, new=2.005864e-17, abs_change=None, rel_change=None
- energy_closure_error: previous=None, new=2.457539e-04, abs_change=None, rel_change=None
- orthogonality_error: previous=None, new=1.228770e-04, abs_change=None, rel_change=None
- kappa_variance: previous=None, new=3.115033e-09, abs_change=None, rel_change=None

## Freeze registry (sec 26)

| module_id | module_name | freeze_status | reason |
|---|---|---|---|
| M01 | conventions | FROZEN | all checks pass and no protected-function violations |
| M02 | coordinate_transforms | FROZEN | all checks pass and no protected-function violations |
| M03 | vector_transforms | FROZEN | all checks pass and no protected-function violations |
| M04 | tensor_transforms | FROZEN | all checks pass and no protected-function violations |
| M05 | pair_enumeration | FROZEN | all checks pass and no protected-function violations |
| M06 | a8_pair_amplitude | FROZEN | all checks pass and no protected-function violations |
| M07 | transverse_projector | FROZEN | all checks pass and no protected-function violations |
| M08 | pair_transfer | FROZEN | all checks pass and no protected-function violations |
| M09 | pair_transfer | FROZEN | all checks pass and no protected-function violations |
| M10 | midpoint_rasterization | FROZEN | all checks pass and no protected-function violations |
| M11 | field_diagnostics | FROZEN | all checks pass and no protected-function violations |
| M12 | differential_operators | FROZEN | all checks pass and no protected-function violations |
| M13 | helmholtz_3d | FROZEN | all checks pass and no protected-function violations |
| M14 | los_projection | FROZEN | all checks pass and no protected-function violations |
| M15 | ray_interface | FROZEN | all checks pass and no protected-function violations |
| M16 | observable_extraction | FROZEN | all checks pass and no protected-function violations |

## Final report questions (sec 29)

- Q1. What exact Git commit was tested?
  - A: `d0c7763cea3b3ec568d79459431656b5378ca2b3`
- Q2. Were PRs #2-#7 present in the tested main?
  - A: Yes - all six merge commits are recorded in `repository_state.json`.
- Q3. Was the working tree clean?
  - A: `True`
- Q4. Did all M01-M16 module tests pass?
  - A: `True` (16/16)
- Q5. Did M03 pass closed-form independent transform validation?
  - A: `True`
- Q6. Did PS1 equal PS1-B exactly?
  - A: `True`
- Q7. Did PS1 remain distinct from PS2?
  - A: `True`
- Q8. Were all valid interface pair slots consumed exactly once?
  - A: `True` (omitted=0, duplicated=0)
- Q9. Did endpoint closure remain nontrivial?
  - A: E_endpoint=4.205e+02 > 0 and |sum_i R_i|=2.176e-15
- Q10. Did nonfinite M11 fields invalidate global sums?
  - A: `True`
- Q11. Did M13 field reconstruction error differ correctly from energy closure error?
  - A: `True` (pad: field=3.897e-17, energy=4.228e-02)
- Q12. Did pure longitudinal Helmholtz recover f_irr=1?
  - A: `True` (f_irr_partition=1.000000)
- Q13. Did pure transverse recover f_sol=1?
  - A: `True` (f_sol_partition=1.000000)
- Q14. Did the equal-energy mixed fixture recover ~0.5/0.5?
  - A: `True` (f_irr=0.500000, f_sol=0.500000)
- Q15. Did M15 reject nonfinite input unconditionally?
  - A: `True`
- Q16. Did the 1e-15 structured field remain accepted?
  - A: `True` (classification=structured_small)
- Q17. Did NaN-masked monotonic Spearman return exactly 1 within tolerance?
  - A: `True` (rs=1.000000)
- Q18. Did the old tied-rank control fail as expected?
  - A: `True` (r_old=1.000000, r_new=0.984732)
- Q19. Did the zero-field full-chain control produce zero observables and NaN correlation?
  - A: `True`
- Q20. Did A/zero/B demonstrate no retained state?
  - A: `True`
- Q21. Was MACS0416 PL1_PM1_PS2 nontrivial?
  - A: E_endpoint=4.364e-04 > 0, E_interface=2.952e-04 > 0, Var(kappa)=3.115e-09 > 0, ray_classification=structured_normal
- Q22. What are the new endpoint and interface energies?
  - A: E_endpoint=4.364493e-04, E_interface=2.952308e-04
- Q23. What are the new Helmholtz fractions?
  - A: f_irr_partition=0.365208, f_sol_partition=0.634792, f_irr_native=0.365118, f_sol_native=0.634636
- Q24. What are the padded and cropped reconstruction/energy/orthogonality errors?
  - A: padded: field=3.897e-17, energy=4.228e-02, orthogonality=2.114e-02; cropped: field=3.905e-17, energy=1.722e-02, orthogonality=8.610e-03
- Q25. Did all RC0-RC6 satisfy E_cov<=0.05?
  - A: `True`
- Q26. Did the scalar-only inverse wrong control still fail strongly?
  - A: `True`
- Q27. Which M01-M16 modules were promoted to L3?
  - A: 16/16
- Q28. Is the complete numerical core now frozen?
  - A: `True`
- Q29. Is full_candidate_rerun_allowed true?
  - A: `True`
- Q30. What is the next permitted experiment?
  - A: PBUF 3D PAIRWISE PRIMARY-CANDIDATE SCIENCE RE-RUN 001 with 5 clusters, PL1_PM1_PS2, RC0, validated/frozen numerical core only.

## Outcome determination (sec 28)

Outcome A - CORE REQUALIFIED: all module, integration, recovery, and covariance gates pass. second_review_status = accepted; full_candidate_rerun_allowed = true.
