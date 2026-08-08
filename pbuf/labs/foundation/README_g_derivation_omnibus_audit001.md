# G derivation omnibus audit 001

## Purpose

Collect the currently plausible PBUF routes to an effective Newton coupling into one fact-finding audit without allowing the measured value of `G` to choose or tune any upstream route.

The measured Newton constant is treated as a valid experimental constant and enters only after candidate construction for comparison. If no independent route closes, the explicit fallback is to retain measured `G` as an external physical constant and defer first-principles derivation.

## Routes

1. **G-free microscopic-length route**
   - `G_eff(L0)=c^3 L0^2/hbar`
   - proton reduced-Compton, electron reduced-Compton, and proton charge-radius seed lengths
   - systematic `alpha_EM^n`, `n=0..12`

2. **G-free structural dimensionless-coupling route**
   - `G_eff=(hbar*c/m_p^2)*alpha_EM^N`
   - frozen structural counts `N={6,12,15,18,30}` from the prior counting audit

3. **Algebraic equivalence audit**
   - verifies that the proton reduced-Compton length route with `n=N/2` is the same mathematics as the dimensionless-coupling route for even `N`
   - equivalent formulas are not counted as independent evidence

4. **Conventional Planck identity control**
   - `G=c^3*l_P^2/hbar`
   - evaluated only post-hoc
   - explicitly circular because conventional `l_P` contains `G`

5. **Absolute constitutive-stiffness route**
   - retained as a legitimate future route
   - current audited foundation does not yet provide an independently normalized absolute stiffness/metric-response map
   - no stiffness is solved from measured `G`

6. **Medium-state dressing route**
   - symbolic form `G_eff(a,state)=G_bare*R_medium(a,state)`
   - retained as a legitimate possibility
   - current foundation does not derive `R_medium`, so no missing multiplier is inserted

7. **Growth / sigma8 route**
   - classified as a diagnostic of a supplied `G_eff(a)` history, not an absolute derivation of `G`
   - growth can be used after a coupling law is derived, not to manufacture the coupling upstream

## Guardrails

- no fit or tuning
- no candidate selection/ranking
- no exponent solved from `G`
- no fractional exponent chosen for agreement
- no `L0`, stiffness, or medium response solved from measured `G`
- no kappa, shear, HST, or lens target input
- no Quantum Engine input
- `alpha_EM` is the electromagnetic fine-structure constant, not PBUF alpha
- measured `G` and conventional Planck length enter post-hoc only
- stdout only; no run directory

## Run

```bash
PYTHONPATH=. python pbuf/labs/foundation/g_derivation_omnibus_audit001.py
```

Return current HEAD SHA, exit code, complete raw stdout, and `git status --short` after the run. Do not modify inputs or repair/reinterpret a failed run. Do not alter the preservation stash.
