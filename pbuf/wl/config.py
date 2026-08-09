"""Frozen execution configuration owned by the established WL labs."""

from pbuf.labs.foundation import m10_coverage_25pct_science001 as BASE
from pbuf.labs.foundation import native_accumulated_full_lensing001 as G3D
from pbuf.labs.foundation import native_full_state_100pct_observer_coverage001 as LEGACY

OBS_BINS = int(BASE.OBS_BINS)
EXTENT = float(BASE.CFG["extent"])
PROPAGATION_STEP = float(BASE.CFG["step"])
PROPAGATION_STEPS = int(BASE.CFG["steps"])
CHECKPOINT = G3D.CHECKPOINT
UNIT_SPEED_TOL = float(G3D.UNIT_SPEED_TOL)

COVERAGE_25PCT = {
    "label": "coverage_25pct",
    "ray_count": int(LEGACY.N25),
    "expected_support_bins": int(LEGACY.EXPECTED_SUPPORT25),
}
COVERAGE_100PCT = {
    "label": "coverage_100pct",
    "ray_count": int(LEGACY.N100),
    "expected_support_bins": int(LEGACY.EXPECTED_SUPPORT100),
}
