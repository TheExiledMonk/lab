# PBUF FOUNDATION — Full-state 100% observer coverage fix 001

## Failure repaired

The first PR #106 execution reached the 100% lane and failed before science comparison with:

`ValueError: operands could not be broadcast together with shapes (285156,) (71289,)`

The 25% lane carries 71289 rays and the 100% lane carries 285156 rays. The reused PR #104 helper `RET._binned_received_3d(...)` contained a hidden historical assumption: it reconstructed initial 3D ray positions by calling the 25% launcher internally rather than using the launch coordinates belonging to the current lane.

## Repair

`native_full_state_100pct_observer_coverage_fix001.py` leaves the merged PR #106 experiment intact and overrides only the observer-side receipt-binning boundary with a launch-aware equivalent.

For each lane the exact `x0/y0` arrays used by G3D propagation are now also used to construct the initial 3D positions needed for `du,dv,dw,t1,t2,tn` and the full 45-channel receipt bank.

An explicit ray-count identity gate checks `x0`, `y0`, detector coordinates, terminal positions and terminal velocity arrays before any receipt arithmetic.

## Scientific guardrails

No source, native response, fast/slow transfer, pair amplitude, PM1/PS2, M10, LOS, G3D propagation, ray count, source coverage, decoder inventory, fitting, target use, weighting, or physical-output scaling changes.

The original #106 main routine, paired 25%/100% aggregation, frozen #105 carry-forward candidates, checks, and output format are preserved.

## Run

```bash
PYTHONPATH=. python pbuf/labs/foundation/native_full_state_100pct_observer_coverage_fix001.py
```

Expected success status:

`FULL_STATE_100PCT_OBSERVER_COVERAGE_EXECUTED`
