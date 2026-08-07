# PS Contract Resolution — FOUNDATION-001-CORRECTION-001

This document records the explicit resolution of the four PS (pair
symmetrisation) lanes for the verified numerical core.

## Definitions (CORRECTION-001 §8)

For a pair `(i, j)` with `j = i + hat_axis`, the per-pair response
is `R_ij = A_ij * v_ij` where `A_ij` is the pair amplitude and `v_ij`
is a 3-vector derived from the projector tensor field.

The four PS lanes are DECLARED DISTINCT (per §8.2):

### PS1-A — raw single-endpoint directional diagnostic

```
v_ij = P_i @ n̂_ij
```

at the source voxel `i`. NOT antisymmetrised. NOT used for physics.
This is the "raw" diagnostic lane and is restricted to diagnostic-only
use.

```python
pair_antisymmetry_expected = False
physics_candidate = False
diagnostic_only = True
```

### PS1 — antisymmetrised source-local

```
v_ij = 0.5 * (v_i - v_j)
```

where `v_i = P_i @ n̂_ij` and `v_j = P_j @ n̂_ij` (partner-side projector
at `j`).

### PS1-B — midpoint antisymmetrised

```
v_ij = 0.5 * (v_i - v_j)
```

Algebraically identical to PS1 (the unscaled `v_ij` is the same).

### PS2 — midpoint-symmetrised projector

```
P̄ = 0.5 * (P_i + P_j)
v_ij = P̄ @ n̂_ij
```

## Algebraic equivalence (§8.4)

BEFORE magnitude normalisation (PM2), PS1, PS1-B, and PS2 produce
different vectors:

- PS1 / PS1-B produce `0.5 (v_i - v_j)`.
- PS2 produces `0.5 (P_i + P_j) n̂_ij = 0.5 (v_i + v_j)` (when the
  projector is evaluated component-by-component).

These differ because `(v_i - v_j) ≠ (v_i + v_j)` in general. The
underlying UN-NORMALISED vectors for PS1 and PS1-B are identical,
which is why they share the equivalence class `PS1-B_EQ_PS1`.

AFTER magnitude normalisation (PM1), the three lanes generally
produce DIFFERENT normalised directions (because the inputs to the
normalisation differ). The candidate registry marks them as distinct
when PM1 is selected.

## Production implementation

The implementation in `pbuf/core/pair_transfer.py` declares all four
lanes as named constants:

```python
PS1_A = "PS1-A"
PS1 = "PS1"
PS1_B = "PS1-B"
PS2 = "PS2"
```

Each has a dedicated code path inside
`_build_pair_response_per_axis` and inside
`build_pair_responses_reference` (the explicit-loop reference).

## Test coverage

The test `_PS_lanes_distinct_test` verifies that the production code
produces distinct R_ij fields for at least PS1-A vs PS2 (which is the
only unambiguously distinct lane — raw vs symmetrised). The test
`_PS1B_PS2_equivalence_class_test` reports the equivalence class
between PS1-B and PS2 (algebraically equivalent before PM, distinct
after PM1 on a spatially-varying projector).

## Resolution

| lane    | physics?  | antisymmetrised? | distinct from PS2? |
|---------|-----------|------------------|---------------------|
| PS1-A   | no (diag) | no               | yes                 |
| PS1     | yes       | yes              | depends on PM       |
| PS1-B   | yes       | yes              | depends on PM       |
| PS2     | yes       | no               | baseline             |

* **PS1-A**: diagnostic-only, raw single-endpoint, never placed under
  strict pair-response antisymmetry.
* **PS1, PS1-B**: produce `0.5 (v_i - v_j)` (algebraically identical to
  each other); classified as `PS1_B_EQ_PS1`.
* **PS2**: midpoint-symmetrised, produces `0.5 (v_i + v_j)`.
* After PM1, PS1 and PS1-B still produce the same result (they share
  the same code path); PS2 differs unless `v_i + v_j` happens to
  parallel `v_i - v_j`.
* The candidate matrix therefore has **three** distinct PS lanes
  (PS1-A, PS1-B=PS1, PS2) at the level of physics candidates.

## Predecessor

The predecessor routed PS1-A, PS1, PS1-B, and PS2 through the SAME
averaging branch `0.5 (v_i + v_j)`. This made the four labels
algebraically identical in production.

CORRECTION-001 splits the four lanes into distinct code paths:
* PS1-A: raw `P_i n̂`.
* PS1 / PS1-B: `0.5 (v_i - v_j)` (same code path).
* PS2: `0.5 (P_i + P_j) n̂` via the midpoint-symmetrised projector.
