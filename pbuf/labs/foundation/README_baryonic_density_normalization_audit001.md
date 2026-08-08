# Baryonic density normalization audit 001

## Purpose

Test the opposite side of the native-scale relation by asking whether the current independent HST/F160W source can be converted into an **absolute baryonic mass-density field in SI units** without using any weak-lensing target.

The prior native-side audit established a stable native transfer but showed that

```text
T_native = (8*pi*G/c^2) * RHO0 * L_cg^2
```

constrains only the product `RHO0 * L_cg^2` unless one side is independently known.

This lab therefore tries to close `RHO0` from source physics alone. Only if `RHO0` closes independently is the opposite-side prediction

```text
L_cg = sqrt(T_native*c^2/(8*pi*G*RHO0))
```

allowed to be evaluated.

## What is audited

For each existing independent HST/F160W cluster source, the lab inventories whether the source path contains:

1. absolute detector-to-flux calibration metadata;
2. source redshift / luminosity-distance information;
3. independently justified stellar mass-to-light or stellar-population conversion;
4. diffuse/hot-gas baryonic mass information;
5. physical pixel area / angular-diameter-distance conversion;
6. line-of-sight depth or deprojection information needed for kg/m^3;
7. preservation of absolute amplitude into `rho2` / `rho3`.

The current common-footprint source is explicitly checked for the existing rule:

```text
rho2 = luminous / max(luminous)
```

which preserves morphology but erases absolute photometric amplitude.

## Guardrails

- no observed kappa pixel values;
- no shear or lensing morphology;
- no lensing amplitude fit;
- no historical `strength=0.18`;
- no fitted mass-to-light ratio;
- no gas normalization inferred from lensing;
- no Quantum Engine;
- no Planck-scale input;
- measured `G` is a macroscopic response anchor only and is not used unless the source-side SI density normalization independently closes;
- gravity is not declared fundamental in PBUF;
- stdout only; no run directory;
- do not alter production code or historical run artifacts.

A clean negative result is expected and valid if the current source does not contain enough independent information. The lab must report `BARYONIC_DENSITY_NORMALIZATION_NOT_YET_CLOSED` rather than invent missing astrophysical conversions.

## Run

From the repository root on the exact PR branch:

```bash
PYTHONPATH=. python pbuf/labs/foundation/baryonic_density_normalization_audit001.py
```

## Coder / runner contract

The runner is an executor only.

Do not modify the lab, constants, source pipeline, thresholds, cluster inputs, or production modules. Do not add catalog values, mass-to-light ratios, redshifts, gas fractions, physical depth assumptions, or cosmological conversions by hand.

If the run fails, return the raw failure exactly as produced. Do not repair or reinterpret it.

Return exactly:

1. current HEAD SHA and branch name;
2. process exit code;
3. complete raw stdout and stderr;
4. `git status --short` after the run;
5. confirmation that the preservation stash was not altered.

Do not delete, move, clean, modify, or commit historical untracked `runs/...` directories.

Do not pop, apply, drop, rewrite, or otherwise alter the preservation stash.
