# Clean Current Weak-Lensing Benchmark 001

## Goal

Measure the current frozen native propagation lane without any benchmark-assisted source construction, HST/F160W proxy, historical control lane, replacement strength, inferred transfer coefficient, fitted normalization, or target-dependent tuning.

## Source discipline

The five canonical Merten v1 kappa FITS files under `PBUF_benchmark` are targets only.

A cluster is eligible only if its local benchmark directory contains exactly one independent 3-D PBUF native-loading FITS cube. The source cube must explicitly declare native loading in FITS metadata (`PBUFROLE=INDEPENDENT_NATIVE_LOADING` / `NATIVE_MATTER_LOADING`, or `BUNIT=PBUF_NATIVE_LOADING` / `PBUF_NATIVE_RHO`). Filename guessing is not used.

If no such source exists, the lab reports `CLEAN_CURRENT_WEAK_LENSING_SOURCE_NOT_AVAILABLE`. It does not construct a source from kappa, normalize light, fetch HST data, or invent an SI-to-native conversion.

## Frozen native lane

```text
independent local native rho3
 -> zero-flux A8 terminal u_fast/u_slow
 -> A_ij = 0.03 Delta u_fast + 0.003 Delta u_slow
 -> terminal c_state geometry
 -> PM1/PS2/M10
 -> LOS
 -> existing G3D observer
```

The 0.03 and 0.003 coefficients are read from the existing frozen A8 implementation. They are not fitted or selected by this benchmark.

## Target reveal

Only after the full prediction is complete is the local kappa target loaded. Target morphology is deterministically sampled onto the observer grid for Pearson/Spearman measurement only. It never enters the source or propagation lane.

## Explicit exclusions

- no kappa-derived source proxy;
- no HST/F160W;
- no network;
- no legacy `strength=0.18` lane;
- no unit-control lane;
- no replacement strength scalar;
- no inferred `~0.02925` coefficient;
- no output normalization to kappa;
- no fit/tuning;
- no GR/Weyl/LCDM/Rmax/QE/Planck input.

## Run

```bash
PYTHONPATH=. python pbuf/labs/foundation/clean_current_weak_lensing_benchmark001.py
```

## Status values

- `CLEAN_CURRENT_WEAK_LENSING_BENCHMARK_EXECUTED`
- `CLEAN_CURRENT_WEAK_LENSING_SOURCE_NOT_AVAILABLE`
- `CLEAN_CURRENT_WEAK_LENSING_BENCHMARK_PARTIAL_EXECUTION`
- `CLEAN_CURRENT_WEAK_LENSING_BENCHMARK_NOT_ESTABLISHED`

`SOURCE_NOT_AVAILABLE` is an informative clean result: it means the current five local weak-lensing targets alone are insufficient for an independent present-day benchmark without fabricating or calibrating a source from the target.