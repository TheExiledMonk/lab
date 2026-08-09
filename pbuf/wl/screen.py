"""Target-blind global tangent detector screen."""

from pbuf.labs.foundation import native_observable_extraction_method_sweep001 as EX
from .launch import RayLaunch


def build_detector_screen(launch: RayLaunch, propagation: dict) -> dict:
    return EX._screen_coordinates(launch.x0, launch.y0, propagation["final_snapshot"])
