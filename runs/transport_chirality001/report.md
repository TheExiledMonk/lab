# PBUF TRANSPORT-CHIRALITY-001

## Frozen conditions
- Constitutive: `u = 0.18 · ρ/ρ_max` (Version A)
- Transport magnitude: `|∇C|`
- Response angle: 90° (the angle is frozen; only the sign of the transverse rotation varies)
- Update rule: direct addition + renormalisation (`transport_lab007.upd_direct_addition`)
- Timestep, normalisation, neighbour transport, Lens-001 dataset: unchanged

## Candidates
| ID | Rule | Selection |
|---|---|---|
| Cand 1 | Global +90 (control) | always `R_{+90}(∇̂C)` |
| Cand 2 | Global −90 | always `R_{−90}(∇̂C)` |
| Cand 3 | Local Laplacian | `laplacian ≥ 0` → +90, else −90 |
| Cand 4 | Symmetric dual (diagnostic) | propagate both, average final positions |
| Cand 5 | Local centre-seeking | response · (mass_center − position) ≥ 0 → +90, else −90 |

## Comparison table

| Candidate | Mirror Δ | Mirror rel. | Bend | Bend. angle | Cons. | Runtime (s) |
|---|---|---|---|---|---|---|
| Cand 1 — Global +90 (control) | 1.324e-04 | 2.00 | 6.618e-05 | 3.028e-04 | 1.1e-16 | 0.0020 |
| Cand 2 — Global −90 | 1.324e-04 | 2.00 | 6.618e-05 | 3.028e-04 | 2.2e-16 | 0.0019 |
| Cand 3 — Local Laplacian | 1.324e-04 | 2.00 | 6.618e-05 | 3.028e-04 | 1.1e-16 | 0.0021 |
| Cand 4 — Symmetric dual (diag.) | 7.174e-10 | 1.08e-05 | 6.618e-05 | 3.028e-04 | 2.2e-16 | 0.0037 |
| Cand 5 — Local centre-seeking | 1.324e-04 | 2.00 | 6.618e-05 | 3.028e-04 | 2.2e-16 | 0.0022 |

## Symmetry summary

| Candidate | Repeatability | Translation | Rotation | Mirror |
|---|---|---|---|---|
| Cand 1 — Global +90 (control) | PASS | PASS | FAIL | **FAIL** |
| Cand 2 — Global −90 | PASS | PASS | FAIL | **FAIL** |
| Cand 3 — Local Laplacian | PASS | PASS | FAIL | **FAIL** |
| Cand 4 — Symmetric dual (diag.) | PASS | PASS | PASS | PASS |
| Cand 5 — Local centre-seeking | PASS | FAIL | FAIL | **FAIL** |

## Required validation: mirror test

The complete mirror test from WEAK-LENSING-VALIDATION-001 was repeated for every candidate.

- **Mirror error (max trajectory delta under y → −y transform):** 1.324e-04 for Cand 1/2/3/5; 7.174e-10 for Cand 4.
- **Trajectory difference (per-photon path delta):** identical to mirror error.
- **Bending difference (`|bend − bend_mirror|`):** 0.0 for every candidate (the bending magnitude is preserved because the mirror only flips the sign of the transverse kick).

## Visualisation

`trajectory_overlay.png` shows, for every candidate, the baseline trajectories, the mirrored trajectories, and the per-step overlay difference, all on identical plotting scales.

## Outcome

**Outcome B.** No local transverse-selection rule among Cand 1–3 and Cand 5 restores mirror symmetry. The relative mirror error is exactly 2.0 in every case, matching the control's intrinsic chirality. The frozen 90° transverse response is intrinsically one-handed: the handedness is set by the choice of `R_{+90}` vs `R_{−90}` at the moment the frozen transport law was selected, and no local geometric criterion can remove it without modifying the frozen transport. Only the diagnostic symmetric-dual propagation (Cand 4), which is not a physical model, produces a mirror-symmetric observable by averaging the two chiral trajectories.

The rotation test fails for the same reason: a 90° rotation of the lens exchanges the role of `x` and `y` in the frozen response, and the rotation cannot be compensated by a local sign choice without future information.

**No further implementation changes.** The chirality of the frozen transport is an intrinsic property of the formulation, not a numerical artefact.
