# A8 zero-flux native benchmark fix 001

## Finding

The five-cluster local benchmark run exposed a repeatable 1.1–1.3% change in the raw native `c_state` integral. The cause is numerical: the historical A8 helper named `reflective` uses NumPy `mode="reflect"`, which mirrors the next interior voxel at the finite-box edge. For this nearest-neighbour transfer operator that is not a conservative no-through-flow boundary when source support reaches the box edge.

Earlier compact centered source tests barely sampled the edge and therefore hid the defect.

## Correction

The historical `boundary="reflective"` implementation is retained unchanged for exact legacy reproducibility.

A new explicit `boundary="zero_flux"` mode uses edge replication for the missing outside neighbour. This keeps the six-neighbour operator symmetric and preserves the global discrete integral in the absence of clipping/source terms.

Only the native raw-`c_state` lane of the local benchmark audit uses `zero_flux`.

No A8 coefficient, source amplitude, K0, epsilon_max, Picard/CG setting, LOS/G3D setting, observed benchmark, or propagation amplitude is changed.

## Run

```bash
PYTHONPATH=. python pbuf/labs/foundation/native_accumulated_full_lensing_local_benchmark001.py
```

## Interpretation

This is an implementation correction to the finite-box boundary condition, not model tuning. Historical controls remain frozen on the historical boundary so their prior results remain directly comparable.