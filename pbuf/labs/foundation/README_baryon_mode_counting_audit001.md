# PBUF FOUNDATION — Baryon Mode Counting Audit 001

## Purpose

Fact-finding only. This lab asks whether the exponent 18 seen in the prior G-free microscopic-scale audit follows independently from the existing PBUF structural hypothesis of three spatial dimensions with two baryonic source modes per dimension.

It does **not** assume that 18 is correct and does **not** use measured `G` or the conventional Planck length to construct a candidate exponent.

## Structural lanes

With `d=3`, `modes_per_dimension=2`, total source modes are `6`.

The lab keeps several possible response-count semantics separate:

- six source modes;
- six longitudinal response components;
- twelve transverse response components;
- eighteen full-vector response components (`6 x 3`), which requires the additional assumption that every source mode drives all three spatial response components independently;
- fifteen unordered distinct mode pairs;
- thirty ordered distinct mode pairs.

The key guardrail is that `6` source modes alone do not derive `18`. If `18` appears, its extra channel assumption is printed explicitly.

For each predeclared structural count `N`, the G-free candidate algebra is

```text
alpha_G,candidate = alpha_EM^N
G_eff = (hbar*c/m_p^2) * alpha_EM^N
L0 = (hbar/(m_p*c)) * alpha_EM^(N/2)
```

Only after all rows are constructed does the lab load conventional `G` and Planck length for post-hoc ratios. Those reference values cannot influence the candidate count or exponent.

## Hard guardrails

- no numerical `G` upstream;
- no conventional Planck length upstream;
- no exponent solved from `G`;
- no fitted or fractional exponent selected to improve agreement;
- no candidate selected or ranked;
- no kappa, shear, HST, or lens benchmark input;
- no GR/Newtonian force, potential, deflection, or calibration law;
- no Quantum Engine;
- electromagnetic `alpha_EM` is not PBUF alpha;
- existing PBUF baryon/alpha structure is not modified;
- stdout only; no run directory.

## Runner contract

After this PR is merged, sync to latest main and run exactly:

```bash
PYTHONPATH=. python pbuf/labs/foundation/baryon_mode_counting_audit001.py
```

Return exactly:

1. Current HEAD SHA
2. Exit code
3. Complete raw stdout
4. `git status --short` after the run

Do not modify, repair, reinterpret, optimize, or change inputs. Do not rerun with altered inputs if it fails. Do not create or edit files. Do not pop, remove, or alter the preservation stash.
