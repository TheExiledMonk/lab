# Propagation Interface Wide-Net Sweep 001

## Goal

Replace serial one-candidate-at-a-time interface tests with one broad discriminatory audit over many already-available native and historical local observables.

The internal reference is the historical **unit-loading M10 interface** only. Observed lensing values are never used.

## Candidate families

The sweep evaluates 18 candidates across all five canonical local benchmark clusters, including:

- native accumulated-state gradient;
- native centered face-bond vector;
- bounded-strain face traction;
- bounded-strain net traction imbalance;
- native `c_state` gradient;
- gradients of native accumulated/state Laplacians;
- frozen fast/slow/combined traction probes;
- one actual frozen fast/slow A8 update from the native equilibrium state;
- historical common, difference, fast, and slow bond controls;
- native-vs-historical state-representation residuals.

## Diagnostics

Every candidate is evaluated on the same source set for:

1. amplitude ratio to historical M10;
2. amplitude-ratio coefficient of variation across sources;
3. component-wise morphology correlation with M10;
4. mean local vector cosine with M10;
5. positive-direction/sign fraction;
6. construction/locality class.

A composite score is printed only for triage. It does **not** fit, rescale, promote, or apply a candidate.

## Guardrails

- five canonical local `PBUF_benchmark` FITS only;
- canonical benchmark loader;
- no network access;
- no observed lensing values;
- no replacement `strength` scalar;
- no normalization/rescaling to M10;
- no fitting/tuning;
- no candidate fed into G3D;
- native zero-flux accumulation unchanged;
- bounded-strain `K0` and `epsilon_max` unchanged;
- historical A8/M10 machinery unchanged;
- no GR/Weyl/LCDM/Rmax/QE/Planck input.

## Run

```bash
PYTHONPATH=. python pbuf/labs/foundation/propagation_interface_wide_net_sweep001.py
```

## Status values

- `PROPAGATION_INTERFACE_WIDE_NET_SWEEP_EXECUTED`
- `PROPAGATION_INTERFACE_WIDE_NET_SWEEP_PARTIAL_EXECUTION`
- `PROPAGATION_INTERFACE_WIDE_NET_SWEEP_NOT_ESTABLISHED`

The intended use is to eliminate weak hypotheses in bulk and reserve follow-up PRs for only the strongest one or two structural survivors.
