# Rmax-forbidden PBUF background closure audit 001

## Purpose

Determine where the current PBUF cosmological background becomes underdetermined when `Rmax` is physically unknown and unavailable as an input.

`Rmax` is not treated as a missing constant to recover. It represents the unknown future expansion-reversal scale, so this audit assigns no numerical value to it.

## Policy

```text
Rmax = UNKNOWN
Rmax numeric use = NOT ALLOWED
Rmax fitting or inheritance = NOT ALLOWED
```

## Dependency chain

```text
alpha_T(a), epsilon0_T(a)
  -> k_max(a)
  -> decay(a)=exp(-a/Rmax)
  -> S(a)
  -> Omega_sigma_raw(a)
  -> flat-today rescale
  -> Omega_sigma(a)
  -> E(a)
  -> H(a)
  -> distances
```

The audit propagates this dependency symbolically and reports the first link blocked by unknown `Rmax`.

The thermal relation

```text
k_max(a)=epsilon0_T(a)-alpha_T(a)
```

can remain available because it does not depend on `Rmax`.

The historical activation factor

```text
decay(a)=exp(-a/Rmax)
```

cannot be evaluated physically while `Rmax` is unknown. Downstream background relations remain unresolved unless an independent present-era law that does not use the future reversal scale is established.

## Redshift boundary

```text
z_observed != z_expansion by assumption
```

No cluster redshift is used in this audit.

## Scientific guardrails

- no numerical `Rmax`;
- no lensing target or morphology input;
- no LambdaCDM distance substitution;
- no measured-G backsolve;
- no historical `0.18`;
- no Quantum Engine execution;
- no Planck-scale input;
- no production-code changes;
- no fabricated distances;
- stdout only, no run directory.

## Run

```bash
PYTHONPATH=. python pbuf/labs/foundation/rmax_forbidden_background_closure_audit001.py
```

Return the branch and HEAD, exit code, complete stdout/stderr, `git status --short`, and confirmation that the preservation stash was unchanged. Preserve all historical untracked `runs/...` directories. If the audit reports that the background remains blocked by unknown `Rmax`, return that result unchanged.
