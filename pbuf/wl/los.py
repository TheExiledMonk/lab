"""Line-of-sight projection for the canonical interface vector."""

import numpy as np

from pbuf.core import los_projection as M14
from .config import EXTENT


def project_interface_to_los(vector: tuple[np.ndarray, np.ndarray, np.ndarray]) -> dict:
    projected = M14.project_vector_to_image_plane(*vector, los_axis="z")
    rx = np.asarray(projected["comp_1"], dtype=np.float64)
    ry = np.asarray(projected["comp_2"], dtype=np.float64)
    grid = np.linspace(-EXTENT, EXTENT, rx.shape[0])
    return {
        "Rx": rx,
        "Ry": ry,
        "los_mag": np.hypot(rx, ry),
        "grid": grid,
        "field": {"xgrid": grid, "ygrid": grid, "rx": rx, "ry": ry},
    }
