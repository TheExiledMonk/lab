# C_STATE Bridge Domain Convergence 001

## Purpose

Verify that the leading native accumulation bridge candidate survives a stricter finite-domain test without changing the mechanism:

`rho -> existing A8 transport -> raw c_state -> six-neighbor bounded-strain equilibrium -> accumulated response`

Only the accumulation-domain size is varied.

## Frozen ingredients

- native grid: `33^3`
- source radius: `3.5`
- source mass: `2.0`
- probe radii: `6, 7, 8, 9, 10`
- constitutive stiffness normalization: `K0=1`
- bounded strain: `epsilon_max=1`
- constitutive law:
  - `W(e)=-(K e_max^2/2) ln(1-(e/e_max)^2)`
  - `sigma(e)=K e/(1-(e/e_max)^2)`
- raw `c_state` amplitude is preserved exactly; no normalization or amplitude rescaling
- same Picard/CG solver tolerances and damping

## Domain ladder

`N = 49, 65, 81, 97`

The same native `c_state` field is embedded at the center of each larger accumulation grid.

## Main question

Does the frozen native bridge move toward the asymptotic `-1` long-range radial exponent as only the artificial outer accumulation boundary recedes?

The lab records:

- far radial fit over fixed physical probes;
- local logarithmic exponent at `r=8`;
- surface response;
- source-region maximum strain;
- native `rho` and `c_state` integrals;
- nonlinear convergence diagnostics.

## Predeclared structural checks

- far fit moves monotonically toward `-1` with increasing domain;
- fixed-`r=8` local exponent moves monotonically toward `-1`;
- largest-domain far fit is within `0.20` of `-1`;
- largest-domain local exponent is within `0.20` of `-1`;
- source strain is stable across domains to relative span `<= 5e-4`;
- raw native `c_state` integral remains equal to the source integral to relative error `<= 1e-12`.

These thresholds are frozen before execution and must not be changed after seeing the result.

## Guardrails

No G, no macroscopic amplitude calibration, no native amplitude rescaling, no fitted/tuned K, no inserted `1/r` response, no spherical-equilibrium shortcut, no Rmax, no cosmology, no lensing target, no Quantum Engine, no Planck input.

## Valid outcomes

- `C_STATE_BRIDGE_DOMAIN_CONVERGENCE_SUPPORTED`
- `C_STATE_BRIDGE_DOMAIN_CONVERGENCE_PARTIAL_SUPPORT`
- `C_STATE_BRIDGE_DOMAIN_CONVERGENCE_NOT_SUPPORTED`

Partial or null support is a valid scientific result and is not permission to repair or tune the lab.

## Runner

```bash
PYTHONPATH=. python pbuf/labs/foundation/c_state_bridge_domain_convergence001.py
```

Return branch/HEAD, exit code, complete raw stdout/stderr, `git status --short`, and `git stash list`. Do not modify, repair, tune, rescale, reinterpret, or merge anything.
