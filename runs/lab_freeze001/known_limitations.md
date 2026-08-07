# PBUF LAB-FREEZE-001 — Known Laboratory Limitations

This document records every limitation of the Version 1 weak-lensing
laboratory that has been **experimentally established** by the completed
validation programme. No speculative items are listed.

Every limitation below is supported by an artifact in `runs/` and is
reproduced from the original report. The laboratory does not attempt to
hide these limitations; they are part of the freeze record.

---

## L1. No cosmological bridge

**Status:** established (OBSERVATION-BRIDGE-001).

The laboratory's κ is a local photon-density distortion (dimensionless
Cartesian units). It is not `Σ/Σ_crit` and depends on no cosmological
distance ratio, no lens redshift, and no source redshift. The published
Frontier-Fields κ maps, in contrast, are surface-mass-density ratios
`Σ/Σ_crit(z_l, z_s)` and depend explicitly on the source redshift
`z_s = 9` in the SaWLens reconstruction.

A direct numerical comparison between laboratory κ and published κ is
therefore only an order-of-magnitude comparison; the symbols coincide,
the physics does not.

**Source:** `runs/observation_bridge001/unit_table.md`.

## L2. Dimensionless constitutive field

**Status:** established (WEAK-LENSING-OBSERVATION-001).

The constitutive field `C(X) = 0.18 · ρ(X) / ρ_max` is dimensionless. No
physical length, mass, or energy is attached to `C`, `∇C`, or `r`. The
photon step `Δs = 0.06` is also dimensionless; the total photon travel
`step · steps = 4.80` is a dimensionless length on the same grid.

**Source:** `runs/observation_bridge001/version_a_chain.md`.

## L3. No Σ_crit

**Status:** established (OBSERVATION-BRIDGE-001).

The critical surface-mass density `Σ_crit(z_l, z_s) = c² / (4πG) ·
D_s / (D_l · D_ls)` is never computed by the laboratory. There is no
input that supplies `z_l` or `z_s`. Therefore no published κ can be
scaled to the laboratory's internal κ; the comparison is symbolic only.

**Source:** `runs/observation_bridge001/unit_table.md`.

## L4. No source-redshift dependence

**Status:** established (OBSERVATION-BRIDGE-001).

Photon initial velocity is `(1, 0)` for every photon; the launch plane
is Launch B (uniform Cartesian 2D on `x ∈ [-extent, -extent + y_span]`,
`y ∈ [-y_span, y_span]`). There is no parameter that supplies a source
redshift, and no term in any laboratory equation contains `z_s`.

**Source:** `runs/source_plane_lab001/report.md`, `runs/observation_bridge001/unit_table.md`.

## L5. No physical angular scale

**Status:** established (OBSERVATION-BRIDGE-001).

The pipeline grid is dimensionless Cartesian on `[-8, 8]²`. There is no
WCS, no RA/Dec, no arcsec-per-pixel. The Frontier-Fields benchmark maps
have CDELT in deg/pixel ranging from `6.25e-5` to `1.13e-5` deg/pixel
depending on cluster; the laboratory has no equivalent quantity.

**Source:** `runs/observation_bridge001/coordinate_audit.md`.

## L6. Benchmark comparison pending (no quantitative agreement)

**Status:** established (INPUT-LAB-001, INPUT-LAB-002, WEAK-LENSING-OBSERVATION-001).

Across all 14 input candidates, all 5 clusters (Abell 2744, MACS J0416,
MACS J1149, Abell S1063, Abell 370), and all 30 propagation configurations
tested in INPUT-LAB-002, the laboratory κ and γ have not demonstrated
quantitative agreement with the published Frontier-Fields κ and γ:

- Pearson(κ) is `nan` in every frozen run (constant κ = -0.5 from the
  histogram method masks any signal).
- Pearson(γ) is at most +0.059 (INPUT-LAB-002, steps=120) and never
  exceeds `|0.1|` under any tested configuration.
- RMS κ is constant across all 12 (steps, step) combinations in the
  INPUT-LAB-002 saturation sweep.

This is an intrinsic property of the κ-extraction rule combined with the
finite photon travel distance, not a tuning failure. See INPUT-LAB-002 Q5.

**Source:** `runs/input_lab002/report.md`, `runs/input_lab001/cluster_statistics.csv`.

## L7. Version A physics only

**Status:** established (constitutive_equations.py).

The laboratory uses Version A exclusively. Versions B (quadratic loading),
C (matter-dependent rigidity), and D (Helmholtz-propagated quadratic
loading) are present in `constitutive_equations.py` but are NOT frozen
and are NOT used by the laboratory. Any change to a different version
constitutes a new laboratory.

**Source:** `constitutive_equations.py`.

## L8. Mirror asymmetry is a property of the transport

**Status:** established (WEAK-LENSING-VALIDATION-001).

The 90° transverse response carries a right-handed rotation `R_90(∇̂C)`.
A mirror reflection flips the handedness of the plane, so the mirrored
response is `R_{-90}(∇̂C)` rather than `R_90`. The mirror-symmetry test
fails by construction with `max transformed trajectory delta = 1.324e-04`
and `relative = 2.00e+00`.

This is documented as a known property of the frozen transport, not an
implementation bug. Removing it would require modifying the transport.

**Source:** `runs/weak_lensing_validation001/report.md`.

## L9. Histogram method is degenerate by construction

**Status:** established (OBSERVABLE-LAB-001, INPUT-LAB-002).

The `histogram` (occupancy) extraction method evaluates κ only at bins
with `N_initial > 0`. Photons launched from `x = -8` leave the launch
column entirely within `step · steps = 4.80` units, so `N_final` at the
launch column is zero and κ reduces to the constant `-0.5` everywhere.

This is a property of the convergence extraction rule, not of the
constitutive law. The histogram method is retained only as a historical
reference and MUST NOT be used as a primary observable.

**Source:** `runs/observable_lab001/report.md`, `runs/input_lab002/report.md`.

## L10. No domain rescaling

**Status:** established (NUMERICAL-CONVERGENCE-001 Group D).

The Group D audit varies the domain half-extent `L ∈ {8, 12, 16, 24}`.
The RMS κ changes monotonically with `L` because the FITS matter field
is rescaled to fill the entire domain. This is a consistency check, not
a refinement study. The frozen laboratory uses `L = 8` exclusively.

**Source:** `runs/numerical_convergence001/report.md`.

## L11. kNN Jacobian neighbourhood sensitivity

**Status:** established (NUMERICAL-CONVERGENCE-001 Group E).

The kNN Jacobian (Group E audit, separate implementation that does NOT
modify the frozen `method_jacobian`) converges with neighbourhood size
(`p_obs ≈ 3.9` for the mean κ, `p_obs ≈ 1.1` for the field RMS). The
primary frozen observable (linear-fit `method_jacobian`) is independent
of the kNN neighbourhood parameter.

**Source:** `runs/numerical_convergence001/report.md`.

## L12. Bending saturates and then collapses

**Status:** established (INPUT-LAB-002).

The travelling-distance sweep (Groups A and B of INPUT-LAB-002) covers
theoretical distances from `2.40` to `38.40`. RMS γ first increases
with distance (photons traverse more of the high-response region),
reaches a maximum around `theoretical_distance ≈ 14`, and then collapses
to a near-zero value once photons exit the domain. This is loss of
signal, not convergence.

**Source:** `runs/input_lab002/report.md`.

## L13. Conservation is at machine precision (a property, not a limitation)

**Status:** established (all validation milestones).

The maximum deviation of photon speed from 1 is `2.220446049250313e-16`
(= 2⁻⁵², IEEE-754 double-precision machine epsilon) for every frozen
run. This is a property of the renormalisation step (per-step
unit-speed rescaling); the laboratory does not "conserve" anything in
the physical sense.

**Source:** all `runs/*/validation.json` and `runs/*/run.json`.

## L14. Single-cluster reproduction

**Status:** established (WEAK-LENSING-OBSERVATION-001).

The regression baseline is for the Abell 2744 cluster only. The
frozen trajectory SHA-256 (`80d8fe47bd0d45672e16c659fc1b095114106eacd6a52634fccde1fc81b41b3f`)
corresponds to this specific input. Other clusters (MACS J0416, MACS
J1149, Abell S1063, Abell 370) are validated against the pipeline but
do not have a frozen trajectory checksum.

**Source:** `runs/source_plane_lab001/run.json`.

---

— PBUF LAB-FREEZE-001