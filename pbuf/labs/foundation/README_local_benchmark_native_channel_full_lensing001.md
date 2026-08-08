# Local Benchmark Native Channel Full Lensing 001

## Purpose

Run the repaired native propagation interface from PR #93 through the complete existing LOS/G3D observer stack on the five canonical local weak-lensing benchmark FITS files.

This lab deliberately returns to the same **local benchmark-assisted source construction** used throughout the recent bridge audits. It does not use the retired HST/F160W independent-source experiment.

## Local-only source lane

For each canonical cluster:

```text
PBUF_benchmark/<cluster>/hlsp_frontier_model_*_merten_v1_kappa.fits
 -> frozen construct_common_proxy
 -> frozen construct_rho_3d
 -> zero-flux A8 terminal u_fast/u_slow
 -> A_ij = 0.03*Delta u_fast + 0.003*Delta u_slow
 -> terminal c_state geometry
 -> PM1 / PS2 / M10
 -> LOS
 -> existing G3D / angular observer
```

The historical unit-loading lane is run from the same benchmark-assisted proxy as an internal control.

## Scientific role

This is an end-to-end **benchmark-assisted structural response test**, not an independent weak-lensing prediction, because the benchmark kappa morphology is used to construct the source proxy and is also the morphology used for end-stage response comparison.

The purpose is to answer the immediate engineering/scientific question: does the newly repaired native channel-transfer bridge run consistently through the complete weak-lensing machinery for all five canonical local files, with no fitted amplitude correction?

## Explicit retirements for this lane

This lab contains no:

- HST F160W source discovery;
- STScI/archive URL discovery;
- download cache;
- urllib/network code;
- network fallback;
- replacement strength scalar;
- inferred 0.02925 coefficient;
- normalization to the benchmark response;
- fitting or tuning.

If a canonical local FITS file is missing, execution fails with `FileNotFoundError` through `pbuf.core.benchmark_data`.

## Outputs

Per cluster the lab reports:

- native/unit LOS RMS and ratio;
- native/unit final G3D angular RMS and ratio;
- Pearson and Spearman response correlation against the frozen benchmark proxy;
- terminal common/history identity;
- G3D unit-speed validity;
- exact local benchmark path used.

## Run

```bash
PYTHONPATH=. python pbuf/labs/foundation/local_benchmark_native_channel_full_lensing001.py
```

## Status values

- `LOCAL_BENCHMARK_NATIVE_CHANNEL_FULL_LENSING_EXECUTED`
- `LOCAL_BENCHMARK_NATIVE_CHANNEL_FULL_LENSING_PARTIAL_EXECUTION`
- `LOCAL_BENCHMARK_NATIVE_CHANNEL_FULL_LENSING_NOT_ESTABLISHED`
