# Native Transfer Factor Decomposition Audit 001

## Question

Can the source-stable inverse transfer scale found in the native local dynamic-response audit be explained by already-frozen local factors, without fitting or introducing a replacement strength coefficient?

The previous audit measured approximately:

- pair-amplitude / native-traction ≈ 0.0087;
- M10 / native-traction ≈ 0.0066.

These remain diagnostic inverse-source ratios only.

## Frozen factors tested

This audit constructs a small, predeclared set of physically identifiable combinations from quantities already present in the code:

- A8 fast coefficient `dt * omega * K`;
- A8 slow coefficient `dt * tau_slow`;
- existing fast/slow coupling factors;
- the exact M10 midpoint half-share;
- the N6 one-sixth neighbour share.

No arbitrary search over real-valued coefficients is performed.

## Guardrails

- five canonical local benchmark FITS only;
- no network access;
- no observed lensing values;
- no replacement strength scalar;
- no normalization/rescaling;
- no fitting/tuning;
- no candidate applied to the native model;
- no candidate fed into G3D;
- no division by 360;
- no GR/Weyl/LCDM/Rmax/QE/Planck input.

## Run

```bash
PYTHONPATH=. python pbuf/labs/foundation/native_transfer_factor_decomposition_audit001.py
```

## Status values

- `NATIVE_TRANSFER_FACTOR_DECOMPOSITION_AUDIT_EXECUTED`
- `NATIVE_TRANSFER_FACTOR_DECOMPOSITION_AUDIT_PARTIAL_EXECUTION`
- `NATIVE_TRANSFER_FACTOR_DECOMPOSITION_AUDIT_NOT_ESTABLISHED`

A numerically close frozen-factor combination is not promoted to a physical law by this audit. It only identifies which existing structural factors deserve a direct forward derivation next.
