# Native accumulated full lensing — local benchmark 001

## Purpose

Correct the data-loading mistake in the first full-lensing audit and standardize access to the five benchmark FITS files already present in the repository.

This PR adds `pbuf.core.benchmark_data` as the canonical local loader for:

- `WL-001_Abell2744`
- `WL-002_MACS0416`
- `WL-003_MACS1149`
- `WL-004_AbellS1063`
- `WL-005_Abell370`

The loader performs **local filesystem I/O only**. It has no URL discovery or download fallback.

## Strength-removal audit

The corrected lab uses the same local benchmark-assisted source construction used by earlier foundation lensing diagnostics and runs three lanes:

1. Legacy `strength=0.18` control.
2. Unit-loading diagnostic control.
3. Native no-strength lane: `rho3 -> raw c_state -> bounded-strain accumulation -> -grad(u) -> existing LOS/G3D -> observer`.

This is intentionally a benchmark-assisted diagnostic of the historical strength factor and end-to-end propagation. Because local kappa morphology supplies the source proxy, it is **not** labeled an independent prediction.

## Guardrails

- no network access;
- no HST discovery or download;
- all five canonical local FITS files must exist;
- no replacement strength scalar;
- no native normalization/rescaling;
- no amplitude matching between lanes;
- no fitting or optimization;
- existing LOS/G3D propagation remains unchanged;
- no GR/Weyl/LCDM machinery;
- no Rmax;
- no Quantum Engine;
- no Planck input.

## Run

```bash
PYTHONPATH=. python pbuf/labs/foundation/native_accumulated_full_lensing_local_benchmark001.py
```

## Valid statuses

- `LOCAL_BENCHMARK_NATIVE_FULL_LENSING_EXECUTED`
- `LOCAL_BENCHMARK_NATIVE_FULL_LENSING_PARTIAL_EXECUTION`
- `LOCAL_BENCHMARK_NATIVE_FULL_LENSING_NOT_ESTABLISHED`

A partial or null result is valid. Do not tune or repair based on scientific output.
