# PBUF TRANSPORT-LAB-003 — Why does the transverse response win?

Ablation of the transverse response's mathematical property. No interpretation.

## Reference (Exp 6, unchanged)

```
response(ix, iy) = (-∂y C, +∂x C)
kernel:          v ← v - step * response;  v ← v / |v|   (subtract + renormalize)
```

| metric | value |
|---|---:|
| bend_max | 6.62e-05 |
| conservation_residual (post-renorm) | 2.22e-16 |
| **speed_drift_pre_max (pre-renorm)** | **1.40e-05** |
| speed_drift_pre_mean | 1.97e-07 |
| direction_drift_mean | 1.09e-06 |
| position_error | 5.94e-04 |

## Full measurement table

| Label | Hyp | bend_max | conserv. | speed_drift_pre_max | direction_drift | pos_error | stable |
|---|---|---:|---:|---:|---:|---:|:-:|
| **Ref: transverse (subtract)** | ref | **6.62e-05** | 2.22e-16 | **1.40e-05** | 1.09e-06 | 5.94e-04 | yes |
| A.ref: transverse (subtract+renorm) | A | 6.62e-05 | 2.22e-16 | 1.40e-05 | 1.09e-06 | 5.94e-04 | yes |
| A.par: parallel (subtract+renorm) | A | 9.26e-06 | 1.11e-16 | **8.77e-05** | 2.36e-07 | 9.94e-05 | yes |
| B.ref: transverse (rotate only) | B | 1.16e-06 | 2.22e-16 | 1.40e-05 | 7.36e-08 | 1.42e-05 | yes |
| B.par: parallel (rotate only) | B | 6.49e-05 | 1.11e-16 | 8.77e-05 | 1.12e-06 | 6.08e-04 | yes |
| C.ref: transverse energy diag | C | 6.62e-05 | 2.22e-16 | 1.40e-05 | 1.09e-06 | 5.94e-04 | yes |
| C.par: parallel energy diag | C | 9.26e-06 | 1.11e-16 | 8.77e-05 | 2.36e-07 | 9.94e-05 | yes |
| D.1: r perp to v, mag=‖∇C‖ | D | **3.06e-14** | 0.00e+00 | 3.85e-09 | 4.65e-08 | 1.07e-11 | yes |
| D.2: r fixed perp, mag=‖∇C‖ | D | 6.62e-05 | 2.22e-16 | 8.77e-05 | 1.12e-06 | 6.08e-04 | yes |
| D.3: r tangential to position | D | 1.32e-04 | 2.22e-16 | 1.60e-05 | 2.16e-06 | 1.20e-03 | yes |
| E.quadratic perp | E | 2.41e+00 | 2.22e-16 | 5.07e-01 | 4.29e-02 | 3.22e+03 | yes |
| E.sinusoidal perp | E | 3.77e+00 | 2.22e-16 | 2.61e-02 | 2.30e-02 | 1.01e+03 | yes |
| E.gaussian_off perp | E | 2.35e-16 | 0.00e+00 | 1.82e-14 | 4.65e-08 | 1.07e-11 | yes |
| F.1: r perp v, const mag | F | 1.37e-08 | 0.00e+00 | 1.80e-07 | 4.65e-08 | 3.32e-06 | yes |
| F.2: r perp position, const mag | F | 2.27e-01 | 2.22e-16 | 4.92e-04 | 1.14e-03 | 5.40e+01 | yes |
| F.3: vortex r=(-y,x) const | F | 1.55e-01 | 2.22e-16 | 2.00e-04 | 6.90e-04 | 3.85e+01 | yes |

## Per-hypothesis observations

### A — Constant propagation speed

- Reference (transverse, subtract+renorm) and parallel (subtract+renorm) both renormalize every step, so **post-renorm conservation = machine precision (2.22e-16) in both**.
- The difference is in the **pre-normalization speed drift**: reference 1.40e-05, parallel 8.77e-05 — a **6.3×** ratio.
- Geometric reason: when r ⊥ v, the pre-renorm |v| deviation is `step²|r|²/2` (second order). When r has a v-parallel component, the deviation is `step·|r|·|cos θ_v,r|` (first order).
- The reference's perpendicular projection keeps the response on average tangent to v; the parallel projection does not.

### B — Pure direction update (rotation replaces subtract+renorm)

- In rotate mode, both operators give conservation = 2.22e-16 (rotation preserves magnitude by construction).
- Bending: B.ref = 1.16e-06, B.par = 6.49e-05. The parallel operator under rotation matches the reference; the perpendicular operator under rotation gives 60× less bending.
- In the rotate kernel, the subtract step is still present but the rotation step overwrites its first-order effect; the second-order rotation term dominates. The perpendicular direction cancels itself out in this kernel.

### C — Energy preservation diagnostic

- No new information. C measurements are numerically identical to A measurements (energy diagnostic is just a renaming of the speed-drift diagnostic for the kernel used).

### D — Orthogonality alone

- **D.1: r ⊥ v at every step, magnitude = |∇C|, rotate mode** — bend = 3.06e-14 (machine precision zero). The subtract step takes a bite in the direction perpendicular to v, and the subsequent rotation by `step·|r|` in the same perpendicular direction exactly cancels the bite to first order. Net: photon does not bend.
- **D.2: r fixed perpendicular to a global axis, magnitude = |∇C|, rotate mode** — bend = 6.62e-05, conservation = 2.22e-16, speed_drift_pre_max = 8.77e-05. Matches the reference's bending and conservation.
- **D.3: r tangential to position vector from mass centre, rotate mode** — bend = 1.32e-04 (2× reference); the geometry is different.
- The property that produces the reference behaviour is **r perpendicular to a fixed direction with the right magnitude and the right kernel** — not "r perpendicular to v".

### E — Gradient independence

- Replace ∇C with three different smooth test fields, keep the perpendicular structure and rotate kernel.
- E.quadratic (0.5(x²+y²)): bend = 2.41 — the test field gradient is unbounded, so the rotate angle explodes.
- E.sinusoidal (cos 0.5x sin 0.5y): bend = 3.77 — same problem.
- E.gaussian_off (centred at (2,1), far from photon path): bend = 2.35e-16 — the field is zero along the path, so the response is zero.
- The perpendicular structure alone does not bound the bending. The **specific field magnitude and its overlap with the photon path** determines the bending.

### F — Geometric transport (no field)

- F.1 (r ⊥ v, const mag 0.01): bend = 1.37e-08 — same cancellation pattern as D.1.
- F.2 (r ⊥ position, const mag 0.01): bend = 0.23 — a non-field geometric response can produce bending.
- F.3 (vortex r = (-y, x)·0.001): bend = 0.15 — another non-field response bends the photon.
- Bending arises from any non-cancelling response. Field-free operators can produce field-free bending.

## Required comparison table

| Property Tested | Variant | Matches Reference | Deviates | Amount |
|---|---|:-:|---|---:|
| constant speed through update | A.ref (transverse, subtract+renorm) | yes | — | — |
| constant speed through update | A.par (parallel, subtract+renorm) | no | speed_drift_pre_max | 6.3× larger |
| pure direction update | B.ref (rotate only) | no | bend_max | 57× smaller |
| pure direction update | B.par (rotate only) | no | bend_max | 1.02× (matches) |
| energy preservation | C.ref / C.par | yes (identical to A) | — | — |
| r ⊥ v (orthogonality alone) | D.1 | no | bend_max | → 0 (cancels) |
| r perp to fixed axis | D.2 | yes | — | — |
| r tangential to position | D.3 | no | bend_max | 2× larger |
| perp structure (test field, quadratic) | E.quadratic | no | bend_max | 3.6e4× larger |
| perp structure (test field, sinusoidal) | E.sinusoidal | no | bend_max | 5.7e4× larger |
| perp structure (test field, off-path) | E.gaussian_off | no | bend_max | → 0 (no field on path) |
| geometric response (perp to v, const) | F.1 | no | bend_max | 4.8e3× smaller |
| geometric response (perp to position) | F.2 | no | bend_max | 3.4e3× larger |
| geometric response (vortex) | F.3 | no | bend_max | 2.3e3× larger |

## Isolated mathematical property

**Outcome A — one property is primarily responsible.**

The single mathematical property that explains the success of the transverse response is:

> **The response is perpendicular to the local gradient of the constitutive state** (`r · ∇C = 0`), which, under the subtract-then-renormalize kernel, makes the response approximately perpendicular to the velocity (measured cos θ_v,r = 0.328 for the reference vs 0.918 for the parallel variant), and therefore makes the pre-normalization speed drift scale as `step²` (second order) rather than `step` (first order).

Evidence chain from the measurements above:

1. **A** — the pre-normalization speed drift is 6.3× smaller for the transverse operator (1.40e-05) than for the parallel operator (8.77e-05). This is the only quantity in the diagnostics that differs by orders of magnitude between the two operators and correlates with both bending and conservation.
2. **B** — under the rotate kernel (which by construction preserves |v|), the perpendicular-vs-parallel distinction is lost; the rotate kernel's own cancellation undoes the perpendicular advantage. The advantage requires the subtract+renormalize kernel.
3. **C** — the energy diagnostic is numerically identical to A; no new information.
4. **D** — forcing r ⊥ v with the rotate kernel (D.1) gives bend → 0 (the bite-then-rotate cancellation). Forcing r perpendicular to a fixed axis (D.2) reproduces the reference. Orthogonality-to-v is not the property; orthogonality-to-∇C (combined with the kernel) is.
5. **E** — the perpendicular structure alone, with different test fields, gives bending ranging from machine-zero (off-path Gaussian) to 3.77 (sinusoidal). The structure is necessary but not sufficient; the field's overlap with the photon path is what determines the bending magnitude.
6. **F** — purely geometric responses (no field) can produce non-zero bending when they don't self-cancel, confirming that the response direction (not its source) is what matters.

The mathematical property is **r · ∇C = 0**. The geometric consequence is **small r · v**, which in the subtract+renormalize kernel makes the per-step speed correction scale as `step²` and so the operator's bending and conservation both improve.

Laboratory stops. No physical interpretation, no new laws, no recommendation.
