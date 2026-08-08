# PBUF FOUNDATION — Known-Source Inverse Response Fingerprint 001

## Purpose

Attack the mass→spacetime problem from the opposite end.

Instead of assuming what the spacetime substrate is or how it couples microscopically, begin with sources whose mass and geometry are specified independently by construction and derive the effective weak-field response fingerprint that any candidate PBUF medium law must reproduce.

## Epistemic split

Derived/effective side:

`known source -> measured/macroscopic weak-field response -> required response fingerprint`

Speculative side:

`stress-energy -> unknown microscopic substrate coupling -> effective spacetime response`

The lab does not claim to identify the substrate or microscopic coupling.

## Controlled sources

Uniform spheres are synthetic controlled inputs. Their `M`, `R`, and `rho` are defined by construction rather than inferred from orbital/gravitational behavior.

The macroscopic benchmark uses measured `G` only in its allowed role as an empirical response coefficient. `G` is not treated as fundamental PBUF microphysics and is not backsolved into a microscopic stiffness.

## Fingerprints tested

- local source curvature versus density;
- surface response versus mass at fixed radius;
- surface response versus radius at fixed mass;
- far-field response versus mass at fixed probe radius;
- fixed-density size scaling;
- the uniform-sphere identity `q_rho R^2 = 3 h00_surface`;
- first-order weak-field additivity as an effective benchmark.

## Forbidden moves

- no fit or tuning;
- no microscopic substrate claim;
- no promotion of `G` to a fundamental PBUF constant;
- no `Rmax`;
- no replacement free parameter;
- no cosmology or distance closure;
- no lensing target, kappa, shear, or morphology;
- no legacy `0.18`;
- no Quantum Engine or Planck-scale input.

## Expected status

`EFFECTIVE_RESPONSE_FINGERPRINT_DERIVED_MICROSCOPIC_CAUSE_SPECULATIVE`

This status means the effective response constraints are established while the causal substrate-level interpretation remains explicitly speculative.

## Runner

```bash
PYTHONPATH=. python pbuf/labs/foundation/known_source_inverse_response_fingerprint001.py
```

Return branch/HEAD, exit code, complete raw stdout and stderr, `git status --short`, and preservation-stash confirmation. Do not repair, rescale, reinterpret, or tune a failing/open result.
