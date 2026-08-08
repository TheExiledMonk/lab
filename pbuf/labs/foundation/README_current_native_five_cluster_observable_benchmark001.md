# Current Native Five-Cluster Observable Benchmark 001

## Goal

Measure the current native propagation path on the established five-cluster local benchmark without reintroducing historical amplitude assistance.

## Data path

All local weak-lensing products are read through `pbuf.core.benchmark_data`:

- `kappa`
- `gamma`
- `gamma1`
- `gamma2`

The benchmark source construction is the existing frozen local route already used by recent foundation audits:

`kappa -> construct_common_proxy -> construct_rho_3d`

No new source format or external source dataset is invented.

## Current native model path

`rho3 -> zero-flux terminal fast/slow channels -> exact A8 pair transfer -> native terminal c_state geometry -> PM1/PS2 -> M10 -> LOS -> photon propagation -> Jacobian -> kappa/gamma1/gamma2`

The fast and slow transfer coefficients are computed from the current frozen A8 constants at runtime. They are not entered as benchmark correction numbers.

## Excluded

- historical `strength=0.18`
- unit-control lane
- replacement strength/amplitude scalar
- inferred `~0.02925` coefficient
- output normalization or rescaling to observed weak lensing
- fitting or tuning
- cluster-specific correction
- HST/F160W source logic
- invented independent 3-D source FITS requirement

## Interpretation

This is the current-model score on the established benchmark setup. Because the existing benchmark source construction begins from kappa, the kappa comparison is not an independent matter-source prediction. The separately loaded gamma/gamma1/gamma2 products are reported as observable comparisons and no target product changes the model response.

## Run

```bash
PYTHONPATH=. python pbuf/labs/foundation/current_native_five_cluster_observable_benchmark001.py
```
