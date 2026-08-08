# Native length / physical-scale mapping audit 001

## Purpose

Test the interpretation from the preceding native-field scaling audit: the missing absolute normalization appears to enter at the local matter -> `c_state` boundary, while M10 already carries additional spatial accumulation.

The lab asks whether the current audited PBUF source/state implementation already contains enough **independently normalized physical scale information** to map SI baryonic density into the dimensionless native state.

## Dimensional closure under test

The established macroscopic response anchor is

```text
q_rho = (8*pi*G/c^2) * rho_SI            [1/m^2]
```

If the local native state is dimensionless, a coarse-grained mapping can only have the schematic form

```text
c_state = q_rho * L_cg^2
```

Write the unknown SI density normalization as

```text
rho_SI = rho_native * RHO0
```

and measure the native transfer directly from the frozen PBUF dynamics:

```text
c_state = T_native * rho_native
```

Then consistency requires

```text
T_native = (8*pi*G/c^2) * RHO0 * L_cg^2
```

This is the central test: native response linearity constrains only the product `RHO0 * L_cg^2`. It cannot separately determine the physical density normalization and physical coarse-graining length.

## What the lab does

1. Re-runs the clean fixed-radius density ladder from the preceding audit with noise-free unit loading.
2. Verifies that `c_state` remains linear in native density and measures `T_native = c/rho_native`.
3. Inventories the scales actually active in the present source/state path (`A8_INIT_DT`, `A8_INIT_K`, step count, grid `DX`, and simulation extent) and refuses to promote numerical/grid quantities to SI physical scales.
4. Inventories native source/state files for explicit SI density/length markers.
5. Demonstrates the `RHO0` / `L_cg` degeneracy using several deliberately arbitrary example lengths. These are identity controls, **not physical candidates**.
6. Reports `SCALE_CLOSURE_NOT_YET_AVAILABLE` unless independently normalized SI density and length scales already exist upstream.

## Guardrails

- gravity is not fundamental in PBUF;
- measured `G` is a macroscopic response anchor only;
- no kappa, shear, HST, lens morphology, or observer target;
- no historical `strength=0.18`;
- no fitting, tuning, optimization, or backward solving from lensing;
- no Quantum Engine;
- no Planck-scale input;
- numerical grid spacing and timestep may not be interpreted as meters/seconds without independent provenance;
- no production-code modification;
- stdout only; no run directory.

## Run

From repository root on the exact PR branch:

```bash
PYTHONPATH=. python pbuf/labs/foundation/native_length_scale_mapping_audit001.py
```

## Coder / runner contract

The coder is an executor only. Do not edit the lab or any production input before or after execution.

Return exactly:

1. current HEAD SHA and branch name;
2. process exit code;
3. complete raw stdout and stderr;
4. `git status --short` after the run;
5. confirmation that the preservation stash was untouched.

If the run fails, return the raw failure. Do not change constants, thresholds, scales, or equations to make it pass.

Do not delete, clean, move, alter, or commit historical untracked `runs/...` directories. Do not pop, apply, drop, or modify the preservation stash.
