# Native Channel Transfer Closure Sweep 001

## Question

Does the actual frozen native zero-flux fast/slow A8 channel evolution generate the transfer needed to turn the now-strong native `c_state` geometry into the correct local propagation-interface scale?

PR #92 showed that native `c_state` bond geometry is a strong structural match to the historical interface after PM1/PS2/M10 and LOS representation, but its unweighted amplitude is roughly 35x too large. This sweep tests several forward channel-transfer constructions at once rather than one hypothesis per PR.

## Frozen native lane

```text
rho3
 -> equal unit native A8 initialization
 -> frozen A8 evolution with boundary=zero_flux
 -> terminal u_fast, u_slow
 -> common mode c=(u_fast+u_slow)/2
 -> difference mode d=(u_fast-u_slow)/2
 -> candidate frozen pair transfer
 -> native c_state geometry
 -> frozen PM1/PS2
 -> frozen M10
 -> frozen LOS
```

The independently produced native `c_state` from the bounded-strain bridge is required to match the terminal common mode from the explicit zero-flux A8 evolution. This is an execution/identity gate, not a normalization.

## Candidates

The run tests in parallel:

- unweighted native `c_state` bond control;
- native `c_state` bond times the frozen common coefficient `0.033`;
- terminal common-mode transfer only;
- terminal difference-mode transfer only;
- terminal full common+difference modal transfer;
- terminal exact direct fast+slow pair law;
- terminal fast-channel-only transfer;
- terminal slow-channel-only transfer;
- bounded-strain slow-traction control;
- historical exact full-transfer positive control.

No inferred effective coefficient such as the descriptive ~0.02925 from earlier audits is inserted.

## Diagnostics

Each candidate is transformed through the same frozen interface machinery and compared M10-to-M10 and LOS-to-LOS across all five canonical local benchmark sources:

- amplitude ratio and cross-source CV;
- component correlation;
- vector cosine;
- positive-direction fraction;
- composite triage score.

The ranking is diagnostic only.

## Guardrails

- five canonical local benchmark FITS only;
- no network access;
- no observed lensing values;
- no replacement `strength` scalar;
- no inferred transfer coefficient applied;
- no normalization or rescaling;
- no fitting or tuning;
- native zero-flux A8 transport unchanged;
- native bounded-strain accumulation unchanged;
- frozen A8 coefficients unchanged;
- frozen PM1/PS2/M10 unchanged;
- frozen LOS projection unchanged;
- no candidate fed into G3D;
- no GR/Weyl/LCDM/Rmax/QE/Planck input.

## Run

```bash
PYTHONPATH=. python pbuf/labs/foundation/native_channel_transfer_closure_sweep001.py
```

## Status values

- `NATIVE_CHANNEL_TRANSFER_CLOSURE_SWEEP_EXECUTED`
- `NATIVE_CHANNEL_TRANSFER_CLOSURE_SWEEP_PARTIAL_EXECUTION`
- `NATIVE_CHANNEL_TRANSFER_CLOSURE_SWEEP_NOT_ESTABLISHED`

A promising result must survive amplitude, morphology, direction, sign, and cross-source stability together. A numerical amplitude coincidence alone is not sufficient.
