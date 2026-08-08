# Native Channel Full Weak-Lensing Five-Cluster 001

## Question

Does the now-closed native terminal fast/slow propagation interface survive the complete existing weak-lensing pipeline on all five canonical clusters without reintroducing a fitted amplitude factor?

## Native prediction lane

```text
independent rho3 source proxy
 -> frozen zero-flux A8 terminal u_fast/u_slow
 -> A_ij = 0.03*Delta u_fast + 0.003*Delta u_slow
 -> native c_state geometry
 -> frozen PM1 / PS2 / M10
 -> frozen LOS projection
 -> existing G3D propagation / angular observer
```

The historical `strength=0.18` lane and the historical unit-loading lane are retained only as controls.

## Benchmark discipline

For every cluster, all prediction lanes are completed before observed weak-lensing kappa values are loaded. The observed map is therefore an end-of-chain external comparison, not an input to construction, normalization, tuning, or candidate selection.

The current independent source remains the normalized HST/F160W luminous-structure proxy used by the existing benchmark stack. It is not an absolute baryonic mass map, so this test is an end-to-end structural/observational weak-lensing test rather than an absolute SI-amplitude closure claim.

## Reported diagnostics

For each of the five clusters, report:

- native/legacy/unit LOS amplitudes;
- final G3D angular RMS amplitude;
- Pearson and Spearman correlation of final angular RMS-angle magnitude with observed kappa;
- the corresponding historical 0.18 and unit-loading control correlations;
- G3D unit-speed and terminal-channel identity gates.

No target correlation is used as an execution gate. Benchmark agreement is reported, not optimized.

## Guardrails

- five canonical local weak-lensing benchmarks only;
- no network access;
- no replacement `strength` scalar;
- no inferred ~0.02925 coefficient;
- no normalization or rescaling;
- no fitting or tuning;
- frozen zero-flux A8 coefficients unchanged;
- frozen PM1/PS2/M10 unchanged;
- frozen LOS/G3D/observer unchanged;
- observed kappa revealed only after predictions;
- no GR/Weyl/LCDM/Rmax/QE/Planck input.

## Run

```bash
PYTHONPATH=. python pbuf/labs/foundation/native_channel_full_weak_lensing_five_cluster001.py
```

## Status values

- `NATIVE_CHANNEL_FULL_WEAK_LENSING_FIVE_CLUSTER_EXECUTED`
- `NATIVE_CHANNEL_FULL_WEAK_LENSING_FIVE_CLUSTER_PARTIAL_EXECUTION`
- `NATIVE_CHANNEL_FULL_WEAK_LENSING_FIVE_CLUSTER_NOT_ESTABLISHED`
