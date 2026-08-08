# M10 Local Interface Decomposition Audit 001

## Question

Where does the historical unit-loading M10 propagation interface acquire its approximately `0.008 * c_state` scale?

PR #85 established that native and unit `c_state` RMS are essentially identical, while the propagation-input vectors differ by roughly 264–292x. This audit decomposes the old M10 route before attempting any new propagation mapping.

## Historical unit route

```text
rho3
 -> unit A8 state
 -> positive-N6 neighbour differences
 -> frozen T1 pair amplitude
      A_ij = (dt * omega * K) * Delta u_fast
           + (dt * tau_slow) * Delta u_slow
 -> PM1/PS2 pair response
 -> M10 midpoint rasterisation
```

Frozen coefficients are read directly from `pbuf.models.a8_state`:

- `coef_fast = A8_INIT_DT * A8_INIT_OMEGA * A8_INIT_K`
- `coef_slow = A8_INIT_DT * A8_INIT_SLOW_TIMESCALE`

The audit records the scale reduction at every stage relative to `c_state` and to the preceding stage.

## Native comparison

```text
rho3
 -> zero-flux raw c_state
 -> bounded-strain accumulated equilibrium u
 -> positive-N6 bond strain Delta u
 -> bounded-strain bond traction
      sigma = K0 * Delta u / (1 - (Delta u / epsilon_max)^2)
```

These native bond quantities are diagnostics only. They are not fed into G3D and are not declared to be the correct propagation interface.

## Guardrails

- five canonical local `PBUF_benchmark` FITS only;
- canonical benchmark loader;
- no network access;
- no observed lensing values used;
- no division by 360 or angular redistribution assumption;
- no replacement `strength` scalar;
- no normalization or rescaling;
- no fitting or tuning;
- no change to A8 coefficients;
- no change to PM1/PS2/M10;
- no change to bounded-strain `K0` or `epsilon_max`;
- no G3D propagation in the native bond comparison;
- no GR/Weyl/LCDM/Rmax/QE/Planck input.

## Run

```bash
PYTHONPATH=. python pbuf/labs/foundation/m10_local_interface_decomposition_audit001.py
```

## Status values

- `M10_LOCAL_INTERFACE_DECOMPOSITION_AUDIT_EXECUTED`
- `M10_LOCAL_INTERFACE_DECOMPOSITION_AUDIT_PARTIAL_EXECUTION`
- `M10_LOCAL_INTERFACE_DECOMPOSITION_AUDIT_NOT_ESTABLISHED`

The result is an audit of the existing equations only. It does not select or fit a replacement interface law.
