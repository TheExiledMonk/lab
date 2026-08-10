"""Exact CPU reference and Vulkan backend selection for frozen pairwise KDE."""

import numpy as np
import observable_lab001 as OLD
from .vulkan_kde_runtime import VulkanKDERuntime


class CpuExactKDE:
    name = "cpu"
    def evaluate(self, u, v, *, values=None, config=None):
        data = np.vstack((np.asarray(u, dtype=np.float64), np.asarray(v, dtype=np.float64)))
        bandwidth = OLD._gaussian_kde_bandwidth(data) if config is None else np.asarray(config, dtype=np.float64)
        result = OLD._diag_kde(data, bandwidth)(data)
        if values is not None:
            # Frozen KDE has no value-weighted mode; rejecting prevents a formula change.
            raise ValueError("frozen exact KDE does not support values")
        return result


class VulkanExactKDE:
    name = "vulkan"
    def __init__(self, workgroup_size=256, runtime=None):
        self._owns_runtime = runtime is None
        self.runtime = runtime or VulkanKDERuntime(workgroup_size)
        self.last_timing = {}
    def evaluate(self, u, v, *, values=None, config=None):
        if values is not None: raise ValueError("frozen exact KDE does not support values")
        u = np.ascontiguousarray(u, dtype=np.float64); v = np.ascontiguousarray(v, dtype=np.float64)
        if u.ndim != 1 or v.shape != u.shape or not u.size: raise ValueError("u and v must be nonempty matching vectors")
        h = OLD._gaussian_kde_bandwidth(np.vstack((u, v))) if config is None else np.asarray(config, dtype=np.float64)
        result, self.last_timing = self.runtime.evaluate(u, v, h)
        return result
    def close(self):
        if self._owns_runtime and self.runtime is not None: self.runtime.close(); self.runtime = None
    def __enter__(self): return self
    def __exit__(self, *_): self.close()


def make_kde_backend(name):
    if name == "cpu": return CpuExactKDE()
    if name == "vulkan": return VulkanExactKDE()
    raise ValueError(f"unsupported KDE backend: {name}")
