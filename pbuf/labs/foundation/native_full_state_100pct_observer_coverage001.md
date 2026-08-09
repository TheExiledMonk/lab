# PBUF FOUNDATION — Full-state 100% observer coverage 001

## Purpose

Test whether the richer observer-side information seen at 25% source-plane coverage becomes clearer when the same frozen native propagation and the exact PR #105 decoder inventory are evaluated over the full 100% source plane.

## Experimental change

Only source-plane sampling coverage changes.

- 25% control: existing 267×267 Cartesian launch over the established 8×8 rectangle.
- 100% lane: deterministic 534×534 Cartesian launch over the full `[-extent,+extent]^2` plane.
- The full lane therefore uses exactly 4× the ray count for 4× the area, preserving approximately the same rays per supported source bin.

Both lanes share the same cluster source and one native M10 interface construction. No source, A8, fast/slow transfer, c_state geometry, PM1/PS2, M10, LOS field, propagation coefficient, step size, propagation length, or observer-decoder weight changes.

## Observer decode

For each lane the lab builds:

1. the received G3D ray state;
2. one target-blind global tangent screen;
3. the full 45-channel receipt bank from PR #104;
4. the exact predeclared PR #105 reconstruction inventory.

All candidate maps for both lanes are completed before observed weak-lensing products are requested.

## Carry-forward diagnostics

Three PR #105 candidates are frozen before this 100% run and reported explicitly:

- `full_whitened_pca_energy_8`
- `full_l2`
- `established2d_l2`

The complete PR #105 candidate inventory is still evaluated so the test remains a wide-net coverage audit.

## Launch-aware receipt repair

The first run exposed a hidden 25%-coverage assumption in the reused PR #104 receipt helper: its initial 3D ray positions were reconstructed internally from the historical 25% launcher. The new 100% lane therefore attempted to combine 285156 terminal rays with 71289 initial rays and failed before any science comparison.

The repair runner `native_full_state_100pct_observer_coverage_launch_fix001.py` changes only that bookkeeping boundary. It supplies the actual lane-specific `x0/y0` launch coordinates to an otherwise equivalent receipt-binning implementation. Physics, propagation, source coverage, channel definitions, decoder inventory, target use, and all guardrails remain unchanged.

## What this tests

The audit asks whether going from 25% to 100% receipt:

- increases spatial support as intended;
- retains or increases full 3D received-state effective rank;
- changes the effective rank of the full 45-channel bank;
- improves morphology correlation of the frozen full-state reconstruction candidates;
- improves the broader target-blind decoder family without changing physics.

## Guardrails

- same native M10 shared by the 25% and 100% lanes;
- exactly 4× rays for exactly 4× source-plane area;
- same decoder inventory in both lanes;
- lane-specific launch coordinates preserved through receipt extraction;
- no observational regression;
- no fitted decoder weights;
- no target-derived channel orientation or selection;
- no physical-output rescaling;
- no cluster-specific decoder choice;
- no historical `strength=0.18`;
- no inferred replacement coefficient;
- no upstream physics changes.

## Run

```bash
PYTHONPATH=. python pbuf/labs/foundation/native_full_state_100pct_observer_coverage_launch_fix001.py
```

Expected success status:

`FULL_STATE_100PCT_OBSERVER_COVERAGE_EXECUTED`
