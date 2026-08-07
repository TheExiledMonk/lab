# Radial Wave Origin Audit 001

Fact-finding continuation of `spatial_coherence_kernel_audit001.py`.

## Question

Does the previously interesting inverse-r spatial amplitude arise generically from radial wave propagation in three spatial dimensions, rather than remaining an arbitrary test kernel?

## What the lab does

- derives the radial amplitude exponent from shell-flux geometry in 1D, 2D, 3D, and 4D;
- verifies that 3D gives amplitude proportional to `1/r`;
- verifies the source-free radial 3D Helmholtz solution `psi = cos(k r)/r` analytically through its residual;
- defines a bounded wave-origin kernel with unit core and exact `L/r` exterior;
- compares its uniform-sphere average to the prior `1/(1+r/L)` kernel at large `R/L`;
- propagates the same seven microscopic candidate amplitudes through the same target-blind mass-radius cases;
- explicitly audits the circularity that the conventional numerical Planck length already embeds Newton's G via `G = c^3 l_P^2 / hbar`.

## What the lab does not do

It does not claim a deeper substrate exists, does not derive electromagnetism, does not derive gravity, does not load an external numerical value of G, does not use GR/Newtonian force or potential laws, does not use lensing/HST/kappa/shear targets, does not fit or tune anything, and does not select a microscopic candidate.

## Runner contract

After merge, sync to latest `main` and run exactly:

```bash
PYTHONPATH=. python pbuf/labs/foundation/radial_wave_origin_audit001.py
```

Return current HEAD SHA, exit code, complete raw stdout, and `git status --short` after the run.

Do not repair, modify, reinterpret, optimize, change inputs, or rerun with altered inputs if it fails. Do not alter or pop the preservation stash.
