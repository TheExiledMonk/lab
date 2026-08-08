# Native accumulated full lensing 001

## Purpose

Run the existing real-cluster G3D lensing/observer stack with the historical `STRENGTH=0.18` source-amplitude training wheel isolated to a control lane, while the prediction lane uses the newly supported native chain:

`rho3 -> raw c_state -> bounded-strain accumulated response -> -grad(u) -> existing LOS/G3D -> observer`

No replacement strength scalar is introduced.

## Three lanes

1. **Legacy control** — exact historical `strength=0.18` A8/M10/G3D route.
2. **Unit diagnostic control** — same historical route with unit loading only to expose the old scalar's amplitude effect; this is not a physical prediction.
3. **Native prediction** — raw `c_state` with no strength multiplier, bounded-strain six-neighbor accumulation, deformation gradient, existing LOS/G3D ray tracker and existing observer.

The independent source is the existing common-footprint HST/F160W luminous-structure proxy. Benchmark kappa pixel values are withheld until all three prediction/control lanes are complete, then used only for end-of-chain comparison.

## Interpretation of strength

The historical code applies `STRENGTH=0.18` at initialization:

- `u_slow0 = strength * rho3`
- `u_fast0 = strength * rho3 + strength * injection_noise`

It is therefore an initial source/loading amplitude multiplier, not a ray-bending coefficient.

The native lane tests whether its role can disappear entirely because source loading and long-range response are now supplied by the derived native bridge.

## Important limitation

The current independent HST/F160W source is normalized luminous structure, not an absolute baryonic mass map. Therefore this is a full lensing propagation/comparison audit, but any absolute observational amplitude agreement remains diagnostic rather than an SI-normalized prediction.

## Guardrails

- no replacement `strength` scalar
- no native amplitude normalization or rescaling
- no fitting or optimization against kappa
- no benchmark pixels before prediction lanes are complete
- no GR potential decomposition
- no LCDM
- no Rmax
- no Quantum Engine
- no Planck input
- legacy `0.18` confined to historical control
- unit lane explicitly diagnostic only
- existing G3D ray tracker preserved

## Run

```bash
PYTHONPATH=. python pbuf/labs/foundation/native_accumulated_full_lensing001.py
```

## Valid statuses

- `NATIVE_ACCUMULATED_FULL_LENSING_EXECUTED`
- `NATIVE_ACCUMULATED_FULL_LENSING_PARTIAL_EXECUTION`
- `NATIVE_ACCUMULATED_FULL_LENSING_NOT_ESTABLISHED`

These statuses describe execution/bridge viability. Observed-lensing correlation is reported but deliberately not used as a pass/fail tuning gate.

## Runner contract

Return:

1. current branch and HEAD SHA
2. exit code
3. complete raw stdout
4. complete raw stderr
5. `git status --short`
6. `git stash list`

Do not modify, repair, tune, rescale, normalize, reinterpret, or merge anything. A partial or unsuccessful scientific result is valid and is not permission to change the code.
