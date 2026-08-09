import os
import unittest
from unittest import mock

import numpy as np

from pbuf.labs.foundation import los_consistent_ray_geometry001 as GEO
from pbuf.wl.backends import VulkanBackend, VulkanRuntime, vulkan_available
from pbuf.wl.backends.vulkan_runtime import VulkanUnavailableError, _device_index
from pbuf.wl.config import CHECKPOINT, EXTENT, PROPAGATION_STEP, PROPAGATION_STEPS
from pbuf.wl.launch import RayLaunch
from pbuf.wl.propagation import PropagationConfig


class VulkanBackendTests(unittest.TestCase):
    def setUp(self):
        grid = np.linspace(-EXTENT, EXTENT, 64)
        gx, gy = np.meshgrid(grid, grid)
        self.field = {"xgrid": grid, "ygrid": grid,
                      "rx": 1e-3*gx + 2e-4*gy, "ry": -3e-4*gx + 7e-4*gy}
        x = np.array([.010, .011, .012, .013, .014, .015, .016])
        y = np.array([.010, .012, .014, .016, .018, .019, .020])
        self.launch = RayLaunch(x, y, np.zeros(7), np.zeros(7), "synthetic", 1)
        self.config = PropagationConfig(PROPAGATION_STEP, PROPAGATION_STEPS, CHECKPOINT)

    def test_device_index_default_explicit_and_invalid_text(self):
        with mock.patch.dict(os.environ, {}, clear=True): self.assertEqual(_device_index(), -1)
        with mock.patch.dict(os.environ, {"PBUF_VULKAN_DEVICE_INDEX": "0"}): self.assertEqual(_device_index(), 0)
        with mock.patch.dict(os.environ, {"PBUF_VULKAN_DEVICE_INDEX": "bad"}):
            with self.assertRaisesRegex(VulkanUnavailableError, "must be an integer"): _device_index()

    def test_ray_count_mismatch(self):
        bad = RayLaunch(self.launch.x0, self.launch.y0[:-1], self.launch.vx0,
                        self.launch.vy0, "bad", 1)
        with self.assertRaisesRegex(ValueError, "ray-count mismatch"):
            VulkanBackend._validated_inputs(self.field, bad)

    @unittest.skipUnless(vulkan_available(), "float64 Vulkan compute unavailable")
    def test_synthetic_parity_tail_and_repeatability(self):
        cpu_cp, _ = GEO._propagate_g3d(self.field, self.config.step, self.config.steps,
                                       self.launch.x0, self.launch.y0)
        with VulkanBackend(workgroup_size=64) as backend:
            first = backend.propagate(self.field, self.launch, self.config)
            second = backend.propagate(self.field, self.launch, self.config)
        for step in GEO.CHECKPOINTS:
            for name in ("x", "y", "z", "vx", "vy", "vz"):
                np.testing.assert_allclose(first["checkpoints"][step][name], cpu_cp[step][name],
                                           rtol=1e-10, atol=1e-12)
                self.assertTrue(np.array_equal(first["checkpoints"][step][name],
                                               second["checkpoints"][step][name]))

    @unittest.skipUnless(vulkan_available(), "float64 Vulkan compute unavailable")
    def test_workgroup_size_parity(self):
        with VulkanBackend(64) as a, VulkanBackend(256) as b:
            aa = a.propagate(self.field, self.launch, self.config)["final_snapshot"]
            bb = b.propagate(self.field, self.launch, self.config)["final_snapshot"]
        for name in ("x", "y", "z", "vx", "vy", "vz"):
            self.assertTrue(np.array_equal(aa[name], bb[name]))


if __name__ == "__main__": unittest.main()
