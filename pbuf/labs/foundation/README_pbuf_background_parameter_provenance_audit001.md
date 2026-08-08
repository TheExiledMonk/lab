# PBUF background parameter provenance audit 001

## Purpose

Recover the numerical provenance of the five background parameters left open by `PBUF-FOUNDATION-BACKGROUND-DISTANCE-RECOVERY-AUDIT-001`:

- `Rmax`
- `alpha_resolved`
- `Omega_r0`
- `BARYON_FRACTION`
- `H0`

This is a **fact-finding provenance audit**, not a parameter fit and not a background reconstruction.

## Method

The lab scans repository text/data sources for each target symbol and aliases, records nearby numeric candidates and context, and distinguishes explicit symbol/value bindings from mere mentions.

A candidate is not automatically authoritative just because it exists in the repository. Historical files, optimisation output, tests, examples, synthetic controls, benchmarks, and unrelated symbols remain evidence only.

The lab may report a single explicit repository binding, multiple conflicting bindings, or no audited binding. Even a single binding is marked for provenance review before it can enter the physical PBUF background.

## Hard guardrails

- `alpha_qm` from the thermal cache is **not** identified with `alpha_resolved` by numerical closeness.
- no kappa, shear, or lensing target;
- no LambdaCDM distance substitution;
- no fitting or tuning;
- no measured-G backsolve;
- no historical `0.18`;
- no Quantum Engine execution;
- no Planck-scale input;
- no `E(a)` or `H(a)` reconstruction from unreviewed candidates;
- no production-code changes;
- stdout only; no run directory.

## Scientific interpretation

The prior recovery audit already established that the thermal table and V11 equation wiring are reproducible. The present question is narrower:

> Can the repository prove the authoritative numerical value and role of every remaining background parameter?

A value should only become physical input after its exact source, role, and current validity are established.

## Run

From repository root on the exact PR branch:

```bash
PYTHONPATH=. python pbuf/labs/foundation/pbuf_background_parameter_provenance_audit001.py
```

## Runner contract

The runner is an **executor only**.

Do not modify the lab, constants, repository data, production modules, thermal cache, or candidate parameter files.

Do not select or insert values manually for `Rmax`, `alpha_resolved`, `Omega_r0`, `BARYON_FRACTION`, or `H0`.

Do not repair a failure or reinterpret an open result.

Return exactly:

1. current HEAD SHA and branch name;
2. process exit code;
3. complete raw stdout and stderr;
4. `git status --short` after the run;
5. confirmation that the preservation stash was untouched.

Do not delete, move, clean, modify, or commit historical untracked `runs/...` directories.

Do not pop, apply, drop, rewrite, or otherwise alter the preservation stash.
