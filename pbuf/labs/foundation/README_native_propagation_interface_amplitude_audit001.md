# Native propagation interface amplitude audit 001

## Question

Where does the approximately 10^2 native/unit amplitude mismatch first appear in the already-established five-cluster full-lensing pipeline?

## Frozen comparison

Unit historical diagnostic:

`rho3 -> unit A8 c_state -> M10 interface vector -> LOS -> existing G3D`

Native route:

`rho3 -> zero-flux raw c_state -> bounded-strain accumulated u -> -grad(u) -> LOS -> existing G3D`

The audit records RMS and selected max amplitudes at each stage, plus the native/unit ratio at:

1. A8 `c_state`;
2. propagation-input 3-vector (`M10` versus `-grad(u)`);
3. LOS-projected 2-vector;
4. final G3D angular response.

It also reports the native accumulated deformation amplitude and the internal conversion ratios `native gradient / native c_state` and `unit M10 vector / unit c_state`.

A `>=10x` marker is used only to identify the first stage where a large mismatch becomes visible. It is diagnostic only and is not a fit, gate, normalization target, or physical threshold.

## Guardrails

- five existing local `PBUF_benchmark` FITS datasets only;
- canonical `pbuf.core.benchmark_data` loader;
- no network access;
- no observed lensing values used;
- no historical `strength=0.18` lane required for this audit;
- unit loading is diagnostic only;
- no replacement strength scalar;
- no normalization or rescaling;
- no fitting or tuning;
- native zero-flux `c_state` boundary preserved;
- bounded-strain solver preserved;
- M10 historical unit route preserved;
- LOS/G3D machinery preserved;
- no GR/Weyl/LCDM/Rmax/QE/Planck inputs.

## Run

```bash
PYTHONPATH=. python pbuf/labs/foundation/native_propagation_interface_amplitude_audit001.py
```

Valid statuses:

- `NATIVE_PROPAGATION_INTERFACE_AMPLITUDE_AUDIT_EXECUTED`
- `NATIVE_PROPAGATION_INTERFACE_AMPLITUDE_AUDIT_PARTIAL_EXECUTION`
- `NATIVE_PROPAGATION_INTERFACE_AMPLITUDE_AUDIT_NOT_ESTABLISHED`
