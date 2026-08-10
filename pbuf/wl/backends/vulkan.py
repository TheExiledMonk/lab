"""Vulkan compute implementation of the frozen canonical G3D propagation."""

from __future__ import annotations

import time
import numpy as np

from pbuf.labs.foundation import g3d_angular_received_distribution001 as ANG
from pbuf.labs.foundation import g3d_native_angular_detector_image001 as DET
from pbuf.labs.foundation import los_consistent_ray_geometry001 as GEO
from ..config import UNIT_SPEED_TOL
from .vulkan_runtime import VulkanRuntime


class VulkanBackend:
    """Float64 Vulkan backend; Vulkan changes execution, not physics."""

    def __init__(self, workgroup_size: int = 256, runtime: VulkanRuntime | None = None):
        self._owns_runtime = runtime is None
        self.runtime = runtime or VulkanRuntime(workgroup_size)
        self.last_timing: dict[str, float] = {}

    def close(self):
        if getattr(self, "_owns_runtime", False) and getattr(self, "runtime", None) is not None:
            self.runtime.close()
            self.runtime = None

    def __enter__(self): return self
    def __exit__(self, *_): self.close()
    def __del__(self): self.close()

    @staticmethod
    def _validated_inputs(field, launch):
        lengths = {name: len(np.asarray(getattr(launch, name))) for name in ("x0", "y0", "vx0", "vy0")}
        if len(set(lengths.values())) != 1:
            detail = ", ".join(f"{key}={value}" for key, value in lengths.items())
            raise ValueError(f"Vulkan propagation ray-count mismatch: {detail}")
        n = lengths["x0"]
        if n == 0:
            raise ValueError("Vulkan propagation ray-count mismatch: ray count must be positive")
        xgrid = np.asarray(field["xgrid"], dtype=np.float64)
        ygrid = np.asarray(field["ygrid"], dtype=np.float64)
        rx = np.asarray(field["rx"], dtype=np.float64)
        ry = np.asarray(field["ry"], dtype=np.float64)
        if rx.shape != (ygrid.size, xgrid.size) or ry.shape != rx.shape:
            raise ValueError(f"Vulkan field shape mismatch: xgrid={xgrid.size}, ygrid={ygrid.size}, rx={rx.shape}, ry={ry.shape}")
        zeros = np.zeros(n, dtype=np.float64)
        # The frozen CPU kernel defines initial transverse velocity as zero and vz as one.
        arrays = [np.asarray(launch.x0, dtype=np.float64), np.asarray(launch.y0, dtype=np.float64),
                  zeros, zeros, zeros, np.ones(n, dtype=np.float64),
                  xgrid, ygrid, rx.ravel(), ry.ravel()]
        return n, xgrid, ygrid, rx, ry, arrays

    def propagate(self, field, launch, config) -> dict:
        n, xgrid, ygrid, rx, ry, arrays = self._validated_inputs(field, launch)
        checkpoints_wanted = tuple(int(k) for k in GEO.CHECKPOINTS)
        if max(checkpoints_wanted) >= config.steps:
            raise RuntimeError(f"checkpoint >= steps: {max(checkpoints_wanted)} >= {config.steps}")
        started = time.perf_counter()
        output = self.runtime.propagate(arrays, n, xgrid.size, ygrid.size,
                                        float(config.step), int(config.steps), checkpoints_wanted)
        self.last_timing = {"warm_total_seconds": time.perf_counter() - started}
        names = ("x", "y", "z", "vx", "vy", "vz")
        checkpoints = {}
        for slot, step_index in enumerate(checkpoints_wanted):
            snap = {name: output[j][slot].copy() for j, name in enumerate(names)}
            sampled_rx, sampled_ry = GEO._sample(field, snap["x"], snap["y"])
            snap["rx_sample"] = sampled_rx.copy(); snap["ry_sample"] = sampled_ry.copy()
            checkpoints[step_index] = snap
        returned_lengths = {name: output[j].shape[1] for j, name in enumerate(names)}
        if set(returned_lengths.values()) != {n}:
            detail = ", ".join(f"{key}={value}" for key, value in returned_lengths.items())
            raise ValueError(f"Vulkan propagation ray-count mismatch: expected={n}, {detail}")
        final_snapshot = checkpoints[int(config.checkpoint)]
        max_unit_error = float(np.max(output[6]))
        if max_unit_error > UNIT_SPEED_TOL:
            raise RuntimeError(f"G3D unit-speed gate failed: {max_unit_error}")
        g3d = {name: final_snapshot[name] for name in names}
        g3d["max_unit_speed_error"] = max_unit_error
        groups = GEO._source_groups(launch.x0, launch.y0)
        if len(groups) != launch.expected_support_bins:
            raise RuntimeError(f"expected {launch.expected_support_bins} source bins, got {len(groups)}")
        los_mag = np.hypot(rx, ry)
        first = GEO._first_step_geometry(field, launch.x0, launch.y0, checkpoints[1], np.zeros_like(los_mag), los_mag)
        if not first["first_step_exact_pass"]:
            raise RuntimeError("first-step exact geometry gate failed")
        final_ang = ANG._angular_distribution_fields(final_snapshot, groups)
        gates = ANG._moment_gates(final_ang)
        if gates["second_moment_equals_cov_plus_centroid_outer_relative_rms_error"] > ANG.MOMENT_IDENTITY_TOL:
            raise RuntimeError("angular second-moment identity failed")
        if not gates["covariance_psd_pass"]: raise RuntimeError("angular covariance PSD gate failed")
        if not gates["direction_mean_vector_bound_pass"]: raise RuntimeError("angular direction-mean bound failed")
        if float(np.min(np.abs(final_snapshot["vz"]))) <= DET.VZ_MIN:
            raise RuntimeError("final tangent projection vz too small")
        return {"checkpoints": checkpoints, "g3d": g3d, "final_snapshot": final_snapshot,
                "groups": groups, "first_step": first, "angular_gates": gates, "final_ang": final_ang,
                "backend_metadata": {"device": self.runtime.device,
                    "workgroup_size": self.runtime.workgroup_size, "timing": self.last_timing}}

    def propagate_final_snapshot(self, field, launch, config) -> dict:
        """Propagate one streaming tile and return only the frozen checkpoint.

        This uses the identical shader and float64 inputs as :meth:`propagate`,
        but omits global observer diagnostics which are not mathematically
        composable per tile.
        """
        n, xgrid, ygrid, _rx, _ry, arrays = self._validated_inputs(field, launch)
        output = self.runtime.propagate(arrays, n, xgrid.size, ygrid.size,
                                        float(config.step), int(config.steps),
                                        (int(config.checkpoint),))
        names = ("x", "y", "z", "vx", "vy", "vz")
        snapshot = {name: output[j][0].copy() for j, name in enumerate(names)}
        max_unit_error = float(np.max(output[6]))
        if max_unit_error > UNIT_SPEED_TOL:
            raise RuntimeError(f"G3D unit-speed gate failed: {max_unit_error}")
        snapshot["max_unit_speed_error"] = max_unit_error
        return snapshot
