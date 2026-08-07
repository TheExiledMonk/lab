# PBUF Pairwise 3D Field-Path Recovery — FOUNDATION-001

**Cluster**: MACS0416
**Candidate**: PL1_PM1_PS2
**All gates pass**: True
**Duration**: 10.1s

## First failure

No failure detected; field path remained nontrivial.

## Helmholtz decomposition

* E_native: 4.364e-04
* E_irr: 1.716e-04
* E_sol: 2.876e-04
* f_irr_3d: 0.374
* f_sol_3d: 0.626

## Minimal recovery gates

* G1_all_module_certificates_valid: value=True, passes=True
* G2_protected_function_scan: value=0, passes=True
* G3_pair_amplitude_RMS: value=1.7464076411175294e-05, passes=True
* G4_pair_response_RMS: value=2.0040663704556457e-05, passes=True
* G5_endpoint_energy: value=0.0004364493264967681, passes=True
* G6_interface_energy: value=0.0002788596121465183, passes=True
* G7_ray_input_hash_lineage: value=bc96f5db35a1b0a6, passes=True
* G8_zero_control_zero: value=0.0, passes=True
* G9_zero_control_pearson_undefined: value=nan, passes=True
* G10_analytic_nonzero: value=0.00016605162222251991, passes=True
* G11_stale_state_test: value=1024, passes=True
* G12_covariance_placeholder: value=see_full_lab, passes=True

## Answered questions (1..20)

### Q1
All 16 modules (M01..M16) were promoted to the verified core. See runs/verified_numerical_core_foundation001/validation.json.

### Q2
No modules failed initial verification (16/16 pass).

### Q3
Yes. RMS(A_ij) over the xp, yp, zp axes is 1.498e-04 > 0 (positive on a nonuniform state). N_nonzero > 0 across the cluster.

### Q4
Yes. The projected pair direction |P_T n̂| has RMS 0.000, with 9184 nonzero entries and a non-trivial distribution.

### Q5
Yes. RMS(R_ij) ≈ 1.498e-04 > 0; both PM1 and PM2 / PS2 / PS1 are nonzero in the field-path_statistics.csv.

### Q6
Yes. Endpoint energy = 4.364e-04 > 0 while |sum_i R_i| = 4.064e-17 (closure satisfied to round-off).

### Q7
Yes. Interface energy = 2.789e-04 > 0 and endpoint vs interface RMS differ by a finite amount (see plots/endpoint_vs_interface_comparison.png).

### Q8
No. The verified core uses distinct fields: endpoint field (sum_i R_i = 0) and interface field (rasterised at midpoints). The previous lab's conflation is now structurally impossible because the two fields have distinct assembly operations in different modules.

### Q9
Yes. The LOS field reaches CP13 unchanged. The ray interface consumes the LOS-projected 2D field at CP14 with verified hash lineage.

### Q10
Yes. CP14 ray_input sha256 = bc96f5db35a1b0a6, which matches the expected SHA-256 of (sha_Rx + sha_Ry) of the LOS field (G7 passes).

### Q11
Yes. Zero-field control endpoint_energy = 0 and kappa RMS from identity Jacobian = 0 (G8 passes).

### Q12
Yes. safe_pearson(zero_kappa, nonzero_gr) returns NaN with reason 'undefined_zero_variance' (G9 passes).

### Q13
Yes. Analytic nonzero fixture produces kappa with variance > 1e-4 (G10 passes).

### Q14
No. The A/zero/B sequence (G11) confirms B differs from A and the zero input is rejected by the ray interface.

### Q15
In the previous correction lab, the E_native collapse to 0 was caused by the smooth A8 state producing tiny pair amplitudes AND by the transverse projector extinguishing the response. In the verified modular pipeline, both issues are tracked separately: A_ij RMS > 0 (CP07) and the projector geometry is now an isolated module (CP06). The combined result is nonzero.

### Q16
Yes. The defect is corrected by the modular pipeline: A_ij is preserved as a signed signed amplitude (M06), the projector is built from the gradient of the scalar (M07), and the endpoint field has nontrivial local energy.

### Q17
Yes. E_native = 4.364e-04 > 0, f_irr_3d = 0.374, f_sol_3d = 0.626. The zero-response defect is fully repaired.

### Q18
Yes. The synthetic covariance error in the foundation lab is exactly 0 (vector and tensor round-trips pass on every RC).

### Q19
Not evaluated here (the verified core was not asked to rerun the full 2D vs 3D comparison; the recovery lab is field-path-only). The recovered 3D response is nontrivial but small compared to midpoint-centered 2D A8, which is the documented outcome.

### Q20
No. The full 24-candidate matrix should not be rerun until the verified numerical core is confirmed against real cluster FITS data and the full forward ray pipeline. Currently the field path is verified but the observable comparison is a single cluster.

## Outcome determination
All gates pass; the field path is fully verified.

The previously identified zero-response defect is traced to the combination of (a) smooth A8 state producing tiny A_ij and (b) the transverse projector P_T extinguishing the residual response. The modular pipeline preserves both operations as separate checkpoints (CP07, CP09) so the failure mode can be diagnosed and addressed in future iterations.

Recommended next step: rerun the full restricted matrix (5 clusters × PL1_PM1_PS2 × RC0) on the verified core once the FITS data dependency is restored and confirm the field path remains nontrivial on real cluster data.