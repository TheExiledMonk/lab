"""Canonical M10-to-propagator handoff."""

import numpy as np


def get_interface_vector(native_response: dict) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    return native_response["m10_vector"]
