# Cluster distance / redshift closure audit 001

## Purpose

Audit the first still-open physical link in the five-cluster baryonic source chain: measured spectral redshift to physical distance geometry.

The audit deliberately separates the observable

```text
z_observed
```

from the interpretation

```text
z_expansion
```

because the measured wavelength/frequency shift can in principle contain expansion, peculiar-motion, local/gravitational, and possible PBUF-medium propagation contributions. No numerical correction is invented in this lab.

## Historical PBUF background

PBUF V11 historically specifies

```text
E(a)^2 = Omega_m0 a^-3 + Omega_r0 a^-4 + Omega_sigma(a)
H(a)   = H0 E(a)
```

with `Omega_sigma(a)` constructed from thermal-table `alpha_T(a)`, `epsilon0_T(a)`, `kmax(a)`, activation, and flat-today normalization.

Those equations are recorded as provenance only. The current weak-lensing foundation does not expose an audited complete distance API plus the thermal/LUT inputs needed to reproduce that background exactly, so this lab does not recreate it from guessed constants.

## Closure requirements

The distance link is closed only after all of the following are physically justified:

1. `z_observed -> z_expansion` decomposition or a demonstrated negligible correction;
2. exact audited PBUF `H(z)` / `E(a)` with all required inputs;
3. comoving/radial distance integration;
4. angular-diameter distance;
5. luminosity distance;
6. validation of distance duality / reciprocity under PBUF photon-medium propagation.

The last item matters because `D_L=(1+z)^2 D_A` must not be silently assumed if the same medium can alter photon frequency or amplitude during propagation.

## Guardrails

- no kappa pixels;
- no shear;
- no lensing morphology or amplitude;
- no LambdaCDM distance silently substituted for PBUF distance;
- no fitted redshift correction;
- no G backsolve;
- no stellar M/L or gas normalization;
- no Quantum Engine or Planck-scale input;
- gravity is not fundamental in PBUF;
- stdout only; no run directory.

## Run

From repository root on the exact PR branch:

```bash
PYTHONPATH=. python pbuf/labs/foundation/cluster_distance_redshift_closure_audit001.py
```

## Runner contract

The runner is an executor only.

Do not modify the lab, constants, classifications, source redshifts, historical equations, or production modules before or after execution.

Do not add cosmological parameters, distance values, redshift corrections, peculiar velocities, gravitational-redshift corrections, medium-redshift corrections, or a distance-duality rule by hand.

If the lab reports the distance geometry is still open, return that result unchanged.

Return exactly:

1. current HEAD SHA and branch name;
2. process exit code;
3. complete raw stdout and stderr;
4. `git status --short` after the run;
5. confirmation that the preservation stash was untouched.

Do not delete, clean, move, modify, or commit historical untracked `runs/...` directories.

Do not pop, apply, drop, rewrite, or otherwise alter the preservation stash.
