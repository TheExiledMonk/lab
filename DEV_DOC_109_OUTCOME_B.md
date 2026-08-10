# Dev Doc 109 — Outcome B

Vulkan float64 G3D propagation is operational and repeatable, but full
downstream observer parity is not established.

- worst CPU/Vulkan ray-state difference: approximately `1.78e-15`
- screen parity: passed
- same-device repeatability: passed
- 100% ray count: `285156`
- 100% source support: `4096`
- warm Vulkan speedup: `3.68x`

The machine-scale propagation difference and downstream observer sensitivity
are preserved as historical evidence. No CPU/GPU output correction is used.
