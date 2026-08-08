# Interface-to-Interface Survivor Sweep 001

## Question

Which native/local survivor quantities remain plausible once they are transformed through the same frozen PM1/PS2/M10 representation as the historical internal reference?

PR #91 showed that direct raw-field morphology comparisons are not representation-fair: even historical bond ingredients correlate only weakly with final M10 before pair response and midpoint rasterisation. This audit therefore compares interface to interface, then LOS to LOS.

## Frozen reference lane

```text
rho3
 -> historical unit-loading A8 state
 -> historical pair amplitudes / PM1 / PS2
 -> historical M10
 -> frozen LOS projection
```

This is an internal structural reference only. No observed lensing values are used.

## Candidate lane

```text
rho3
 -> zero-flux raw c_state
 -> bounded-strain accumulated equilibrium u
 -> candidate local pair amplitudes
 -> frozen PM1 / PS2
 -> frozen M10
 -> frozen LOS projection
```

Twelve survivors are tested in one run:

- native accumulated bond with native `c_state` geometry;
- native accumulated bond with accumulated-`u` geometry;
- native bounded-strain traction with both geometry choices;
- frozen slow traction transfer with both geometry choices;
- frozen fast traction transfer;
- frozen fast+slow traction transfer;
- native `c_state` bond;
- actual one-step fast-update bond;
- actual one-step slow-update bond;
- historical common-bond control.

No candidate is normalized or rescaled to historical M10.

## Diagnostics

For each candidate and each of the five canonical local sources:

- M10 amplitude ratio;
- M10 component correlation;
- M10 local vector cosine and positive-direction fraction;
- LOS amplitude ratio;
- LOS component correlation;
- LOS vector cosine and positive-direction fraction;
- source-to-source amplitude CV.

A composite score is used only for triage. It does not modify a field and does not establish a physical law.

## Guardrails

- five canonical local `PBUF_benchmark` FITS only;
- canonical benchmark loader;
- no network access;
- no observed lensing values;
- no replacement `strength` scalar;
- no normalization or amplitude rescaling;
- no fitting or tuning;
- native zero-flux accumulation unchanged;
- bounded-strain `K0` and `epsilon_max` unchanged;
- historical A8/PM1/PS2/M10 unchanged;
- frozen LOS projection unchanged;
- no candidate fed into G3D;
- no GR/Weyl/LCDM/Rmax/QE/Planck input.

## Run

```bash
PYTHONPATH=. python pbuf/labs/foundation/interface_to_interface_survivor_sweep001.py
```

## Status values

- `INTERFACE_TO_INTERFACE_SURVIVOR_SWEEP_EXECUTED`
- `INTERFACE_TO_INTERFACE_SURVIVOR_SWEEP_PARTIAL_EXECUTION`
- `INTERFACE_TO_INTERFACE_SURVIVOR_SWEEP_NOT_ESTABLISHED`

The ranking is diagnostic triage only. Follow-up work should focus only on candidates that survive both interface-level and LOS-level morphology/direction checks without fitted amplitude correction.
