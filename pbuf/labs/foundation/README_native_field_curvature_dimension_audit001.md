# Native field curvature-dimension audit 001

## Purpose

Determine how the existing native PBUF medium fields scale with synthetic source density and source radius before assigning any of them the physical role of the mass-induced spacetime response.

This is a **fact-finding scaling audit**, not a lensing run and not a calibration run.

The preceding mass-spacetime-response audit established that nature supplies a measurable macroscopic mass-to-spacetime response while leaving the microscopic PBUF origin and native-field mapping open. This lab deliberately does **not** use `G` at all. It tests the current PBUF implementation itself.

## Synthetic experiment

The lab constructs centered uniform 3D spheres on a fixed `33 x 33 x 33` grid and runs two frozen ladders:

1. **Radius ladder**
   - fixed source density `0.05`
   - nominal radii `{2.5, 3.5, 4.5, 5.5, 6.5}` grid units
   - effective radius is calculated from the actual voxelized sphere volume

2. **Density ladder**
   - fixed nominal radius `4.5` grid units
   - densities `{0.0125, 0.025, 0.05, 0.10, 0.20}`

Initialization is noise-free:

```text
u_slow0 = rho3
u_fast0 = rho3
```

The historical `strength=0.18` coefficient is not used.

## Native quantities audited

- `c_state`
- `|grad c_state|`
- `|laplacian c_state|`
- M10 interface-vector magnitude
- `|div M10|`
- Frobenius norm of the symmetrized spatial gradient of the M10 vector
- M14 image-plane magnitude only as a finite representation check; there is no observer or lensing comparison

## Predeclared scaling diagnostics

At fixed radius, a linear source response should scale approximately as density to the first power.

At fixed density, the measured radius exponent is compared against these predeclared structural classes:

```text
R^0  local source / curvature-like scaling
R^1  one-length integrated / gradient-connection-like scaling
R^2  two-length integrated / strain-metric-deformation-like scaling
```

These names are **diagnostic labels only**. A numerical field landing near one exponent does not prove that it already has the corresponding SI physical units. The current frozen transport contains a grid scale, so the audit identifies scaling behavior/rank, not a completed metric interpretation.

The synthetic source itself supplies controls:

- integrated source must scale exactly as effective-radius cubed at fixed density;
- integrated source must scale linearly with density at fixed radius.

## Hard guardrails

- no Newton `G`
- gravity is not declared fundamental in PBUF
- no kappa
- no shear
- no HST data
- no lens morphology
- no observer comparison
- no `strength=0.18`
- no fit, tuning, optimization, or solved coefficient
- no Quantum Engine
- no Planck-scale input
- no random injection noise
- no accepted run may hit the frozen state clip
- do not change production code
- stdout only; the lab must not create a run directory

## Run

From the repository root, on the exact PR branch:

```bash
PYTHONPATH=. python pbuf/labs/foundation/native_field_curvature_dimension_audit001.py
```

## Coder / runner contract

The coder is an **executor only** for this run.

Do not modify the lab, constants, source ladders, thresholds, production modules, or existing run artifacts before or after execution.

Do not attempt to repair a failing run. If the program exits nonzero, return the complete raw failure exactly as produced so the scientific/code issue can be reviewed separately.

Return exactly:

1. current HEAD SHA and branch name;
2. process exit code;
3. complete raw stdout and stderr, with no summarization or reinterpretation;
4. `git status --short` after the run;
5. confirmation that the preservation stash was not altered.

Do not delete, move, clean, commit, or modify historical untracked `runs/...` directories.

Do not pop, apply, drop, rewrite, or otherwise alter the preservation stash.
