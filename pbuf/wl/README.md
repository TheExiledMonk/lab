# Weak-lensing pipeline

## Propagation backend

`CpuReferenceBackend` is the frozen scientific reference. `VulkanBackend` is
the compute-acceleration backend and must pass CPU parity before scientific
use. Callers may pass either backend instance to `run_wl_pipeline`; CPU remains
the default and Vulkan failures never fall back silently.

Vulkan changes execution, not physics. It owns only the canonical G3D
ray-step kernel, device/runtime setup, buffer transfer, dispatch, checkpoint
capture, and timing. Source construction, native response, M10/LOS, launch,
screen/received-state construction, decoding, reconstruction, and comparisons
remain backend-neutral.

The Vulkan path requires a compute-capable device with shader `float64`, a
Vulkan loader and development headers, a C compiler, and either
`glslangValidator` or `glslc`. It compiles the checked-in GLSL and small
headless Vulkan bridge into a versioned temporary cache; it downloads nothing.
Set `PBUF_VULKAN_DEVICE_INDEX` to select a valid device from the deterministic
compute+float64 candidate list. Use `vulkan_diagnostics()` for non-raising
availability details.

100% coverage is the primary science workload. 25% coverage remains a parity
control and regression lane.
