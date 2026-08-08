# Mass-spacetime response strength audit 001

## Purpose

Establish the absolute **macroscopic** mass/stress-energy -> spacetime-response scale that any microscopic PBUF constitutive law must eventually reproduce, while preserving the PBUF premise that gravity is emergent rather than fundamental.

This lab does **not** derive Newton's constant and does **not** insert Newton/GR into the native PBUF propagation chain. Measured `G` is used only as an experimental macroscopic response anchor.

The key empirical quantities are:

- stress-energy -> curvature compliance
  - `C = 8*pi*G_measured/c^4`
- inverse effective stiffness representation
  - `K = c^4/(8*pi*G_measured)`
- nonrelativistic rest-mass-density curvature source
  - `q_rho = 8*pi*G_measured*rho/c^2`
- external dimensionless weak-field response scale
  - `eta = G_measured*M/(r*c^2)`
  - `|h00| = 2*eta`

`K` is an empirical effective stiffness representation of the measured macroscopic response. The lab must **not** describe it as the fundamental microscopic PBUF stiffness.

## What this lab is testing

1. Compute the empirical absolute stress-energy -> curvature compliance implied by measured gravity.
2. Express its inverse as an effective spacetime stiffness scale.
3. Evaluate clean analytic response scales for Earth and Sun without any lensing data.
4. Verify the uniform-sphere identity
   - `q_rho * R^2 = 6*GM/(R*c^2) = 3*|h00|`.
5. Inventory the historical `STRENGTH=0.18` value only to prove that it is **not** used as physical input here.
6. Leave the microscopic PBUF constitutive origin and native-field mapping explicitly open.

## Interpretation rules

The coder/runner must preserve these distinctions exactly:

- **PBUF claim:** gravity is not assumed fundamental.
- **Measured G:** an experimentally determined macroscopic response anchor.
- **`c^4/(8*pi*G)`:** an empirical effective stiffness representation, not a first-principles PBUF derivation.
- **Legacy `0.18`:** historical dimensionless trajectory/amplitude diagnostic only. It must not be used to infer, scale, tune, or guess the physical response strength.
- **Earth/Sun:** analytic scale checks only; they are not fit targets.

## Hard guardrails

- no fit or tuning
- no optimization or parameter search
- no kappa or shear
- no HST pixel values
- no lens morphology or lens amplitude target
- no weak-lensing comparison in this lab
- no use of `0.18` in a physical equation
- no solving a replacement strength from `0.18`
- no Quantum Engine input
- no Planck length or Planck-unit input
- do not rename measured `G` as a fundamental PBUF coupling
- do not inject Newtonian potential or `h00` into the PBUF observer pipeline
- stdout only; do not create a run directory

## Run

From repository root:

```bash
PYTHONPATH=. python pbuf/labs/foundation/mass_spacetime_response_strength_audit001.py
```

## Runner contract

The coder/runner is an executor only.

Return, in this exact order:

1. current `HEAD` SHA
2. process exit code
3. complete raw stdout from the command above
4. `git status --short` after the run

Do not modify the lab, constants, equations, inputs, or guardrails before or after execution. Do not repair or reinterpret a failed run. Do not create helper patches because an output looks surprising. Do not delete historical untracked `runs/...` directories. Do not touch, pop, apply, drop, or otherwise modify the preservation stash.

If the lab fails, return the failure exactly as produced so the scientific code can be reviewed separately.

## Expected scientific outcome class

This audit is allowed to conclude only that the **empirical macroscopic response scale is known** while the **microscopic PBUF origin and native-field mapping remain open**.

The intended next step, if all checks pass, is a separate native-field dimensional/mapping audit asking which existing PBUF medium variable can carry the physical curvature/loading quantity without promoting measured `G` to a fundamental microscopic law.
