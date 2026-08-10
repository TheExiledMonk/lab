"""Passive Dev140 propagation-parameter and physical-unit audit helpers."""
from __future__ import annotations

from dataclasses import asdict, dataclass
from statistics import median
from typing import Iterable

C_SI = 299_792_458.0
TIME_CLASSES = {
    "NATIVE_TIME_EXPLICIT", "NATIVE_TIME_IMPLICIT_FROM_PROPAGATION_STEP",
    "AFFINE_PATH_PARAMETER_ONLY", "SPATIAL_STEP_ONLY", "NO_NATIVE_TIME_STATE",
}


@dataclass(frozen=True)
class PropagationParameterAudit:
    name: str
    classification: str
    units_provenance: str
    update_equation: str
    relation_to_spatial_displacement: str
    relation_to_direction_update: str
    monotonic: bool
    physical_time_interpretation_established: bool

    def __post_init__(self):
        if self.classification not in TIME_CLASSES:
            raise ValueError("invalid propagation-parameter classification")
        if self.physical_time_interpretation_established and self.classification not in {
            "NATIVE_TIME_EXPLICIT", "NATIVE_TIME_IMPLICIT_FROM_PROPAGATION_STEP"
        }:
            raise ValueError("path/solver parameter cannot be promoted to physical time")

    def to_dict(self):
        return asdict(self)


def current_propagation_parameter_audit():
    return PropagationParameterAudit(
        name="PropagationConfig.step / loop index",
        classification="SPATIAL_STEP_ONLY",
        units_provenance="native coordinate displacement; no clock state",
        update_equation="position_next = position + step * unit_direction + response correction",
        relation_to_spatial_displacement="direct native-coordinate integration increment",
        relation_to_direction_update="indexes successive spatial tangent updates",
        monotonic=True,
        physical_time_interpretation_established=False,
    )


def reject_solver_iteration_as_time(parameter_name: str):
    if "iter" in parameter_name.lower() or "solver" in parameter_name.lower():
        raise ValueError("NUMERICAL_SOLVER_ITERATION_USED_AS_PHYSICAL_TIME")


def l0_over_t0(native_speed: float, audit: PropagationParameterAudit, *, speed_constructed_from_c=False):
    if not audit.physical_time_interpretation_established:
        return {"outcome": "NATIVE_TIME_NOT_ESTABLISHED", "value": None}
    if speed_constructed_from_c:
        return {"outcome": "DIMENSIONAL_MAPPING_AMBIGUOUS", "value": None, "circular": True}
    if native_speed <= 0:
        raise ValueError("native speed must be positive")
    return {"outcome": "L0_OVER_T0_ESTABLISHED", "value": C_SI / float(native_speed), "units": "m/s per (native_length/native_time)"}


def speed_statistics(samples: Iterable[float]):
    x = sorted(float(v) for v in samples)
    if not x or any(v <= 0 for v in x):
        raise ValueError("positive native speed samples required")
    m = median(x)
    dev = sorted(abs(v-m) for v in x)
    def pct(q):
        j=(len(x)-1)*q; lo=int(j); hi=min(lo+1,len(x)-1); return x[lo]+(x[hi]-x[lo])*(j-lo)
    return {"median": m, "MAD": median(dev), "P05": pct(.05), "P95": pct(.95)}


def compare_messenger_speeds(photon_speed, gw_speed, *, shared_operator=False, tolerance=1e-12):
    if photon_speed is None:
        return {"classification": "PHOTON_NATIVE_DYNAMIC_STATE_UNAVAILABLE"}
    if gw_speed is None:
        return {"classification": "GW_NATIVE_DYNAMIC_STATE_UNAVAILABLE", "photon_native_speed": photon_speed}
    ratio=float(photon_speed)/float(gw_speed)
    classification = "ALGEBRAICALLY_IDENTICAL_SPEED_CONSTRAINT" if shared_operator else "INDEPENDENT_COMMON_SPEED_CONSTRAINT"
    return {"classification": classification, "ratio": ratio, "consistent": abs(ratio-1) <= tolerance}


def native_travel_state(path_length_native, audit, native_speed=None):
    timed=audit.physical_time_interpretation_established and native_speed is not None
    return {"contract": "PBUF_NATIVE_TRAVEL_STATE_V1", "path_length_native": float(path_length_native),
            "travel_parameter_native": float(path_length_native/native_speed) if timed else None,
            "native_speed": native_speed, "L0_over_T0_relation": None,
            "physical_seconds_established": False, "physical_metres_established": False}
