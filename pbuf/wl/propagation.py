"""Propagation backend contract and bit-preserving CPU reference backend."""

from dataclasses import dataclass
from typing import Protocol
import numpy as np

from pbuf.labs.foundation import g3d_angular_received_distribution001 as ANG
from pbuf.labs.foundation import g3d_native_angular_detector_image001 as DET
from pbuf.labs.foundation import los_consistent_ray_geometry001 as GEO
from .config import UNIT_SPEED_TOL
from .launch import RayLaunch


@dataclass(frozen=True)
class PropagationConfig:
    step: float
    steps: int
    checkpoint: object


class PropagationBackend(Protocol):
    def propagate(self, field: dict, launch: RayLaunch, config: PropagationConfig) -> dict: ...


class CpuReferenceBackend:
    def propagate(self, field: dict, launch: RayLaunch, config: PropagationConfig) -> dict:
        groups = GEO._source_groups(launch.x0, launch.y0)
        if len(groups) != launch.expected_support_bins:
            raise RuntimeError(f"expected {launch.expected_support_bins} source bins, got {len(groups)}")
        checkpoints, g3d = GEO._propagate_g3d(field, config.step, config.steps, launch.x0, launch.y0)
        if g3d["max_unit_speed_error"] > UNIT_SPEED_TOL:
            raise RuntimeError(f"G3D unit-speed gate failed: {g3d['max_unit_speed_error']}")
        los_mag = np.hypot(field["rx"], field["ry"])
        first = GEO._first_step_geometry(field, launch.x0, launch.y0, checkpoints[1], np.zeros_like(los_mag), los_mag)
        if not first["first_step_exact_pass"]:
            raise RuntimeError("first-step exact geometry gate failed")
        final_snapshot = checkpoints[config.checkpoint]
        final_ang = ANG._angular_distribution_fields(final_snapshot, groups)
        gates = ANG._moment_gates(final_ang)
        if gates["second_moment_equals_cov_plus_centroid_outer_relative_rms_error"] > ANG.MOMENT_IDENTITY_TOL:
            raise RuntimeError("angular second-moment identity failed")
        if not gates["covariance_psd_pass"]:
            raise RuntimeError("angular covariance PSD gate failed")
        if not gates["direction_mean_vector_bound_pass"]:
            raise RuntimeError("angular direction-mean bound failed")
        if float(np.min(np.abs(final_snapshot["vz"]))) <= DET.VZ_MIN:
            raise RuntimeError("final tangent projection vz too small")
        return {"checkpoints": checkpoints, "g3d": g3d, "final_snapshot": final_snapshot,
                "groups": groups, "first_step": first, "angular_gates": gates, "final_ang": final_ang}
