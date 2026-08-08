# Native Local Dynamic Response Audit 001

## Question

Does the large static-native-vs-M10 interface mismatch arise because the bounded-strain solver provides a full equilibrium deformation/traction while the historical propagation interface represents only a local incremental transfer?

PR #86 established that the historical M10 scale is generated primarily before M10 rasterisation by the frozen local T1 pair-transfer equation, while native nearest-neighbour bond strain/traction remains roughly 144–159x larger than M10.

This audit does **not** create a replacement interface law. It measures whether the mismatch can be understood as a local transfer-scale issue.

## Compared quantities

Historical unit route:

```text
rho3
 -> unit A8 state
 -> local T1 pair amplitude
 -> PM1/PS2
 -> M10 interface
```

Native route:

```text
rho3
 -> zero-flux c_state
 -> bounded-strain equilibrium u
 -> positive-N6 bond strain
 -> bounded-strain bond traction
```

The native bond traction is projected through the already-existing frozen A8 update factors only as diagnostics:

```text
fast diagnostic     = (dt * omega * K) * traction
slow diagnostic     = (dt * tau_slow) * traction
combined diagnostic = (fast coefficient + slow coefficient) * traction
```

No diagnostic is fed into G3D.

The audit also reports inverse-source ratios:

```text
historical pair-amplitude RMS / native traction RMS
historical M10 RMS            / native traction RMS
```

These are used only to test source-to-source stability. They are not fitted coefficients and are not inserted into the model.

## Guardrails

- five canonical local `PBUF_benchmark` FITS only;
- canonical benchmark loader;
- no network access;
- no observed lensing values;
- no replacement strength scalar;
- no amplitude normalization or rescaling;
- no fitting or tuning;
- no G3D propagation of a new candidate;
- frozen A8 coefficients unchanged;
- bounded-strain `K0` and `epsilon_max` unchanged;
- no division by 360;
- no GR/Weyl/LCDM/Rmax/QE/Planck input.

## Run

```bash
PYTHONPATH=. python pbuf/labs/foundation/native_local_dynamic_response_audit001.py
```

## Status values

- `NATIVE_LOCAL_DYNAMIC_RESPONSE_AUDIT_EXECUTED`
- `NATIVE_LOCAL_DYNAMIC_RESPONSE_AUDIT_PARTIAL_EXECUTION`
- `NATIVE_LOCAL_DYNAMIC_RESPONSE_AUDIT_NOT_ESTABLISHED`

A stable inferred transfer ratio would motivate a later derivation of the medium's local dynamic response law. It would not by itself establish that law.
