"""Explicit ray-launch definitions."""

from dataclasses import dataclass
import numpy as np

from pbuf.labs.foundation import m10_coverage_25pct_science001 as BASE
from pbuf.labs.foundation import native_full_state_100pct_observer_coverage001 as LEGACY
from .config import COVERAGE_25PCT, COVERAGE_100PCT


@dataclass(frozen=True)
class RayLaunch:
    x0: np.ndarray
    y0: np.ndarray
    vx0: np.ndarray
    vy0: np.ndarray
    coverage_label: str
    expected_support_bins: int


def _make(values, cfg: dict) -> RayLaunch:
    return RayLaunch(*values, cfg["label"], cfg["expected_support_bins"])


def launch_25pct() -> RayLaunch:
    return _make(BASE._launch_expanded_25pct(), COVERAGE_25PCT)


def launch_100pct() -> RayLaunch:
    return _make(LEGACY._launch_full_100pct(), COVERAGE_100PCT)
