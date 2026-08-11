"""Mechanical DEV168 receipt adapter for the frozen historical observer.

This module only renames/repackages native receipt fields.  It deliberately
contains no propagation, steering, normalization, or observational logic.
"""
from __future__ import annotations

import numpy as np

from .native_finite_receipt import NativeReceivedState
from pbuf.labs.foundation import native_observable_extraction_method_sweep001 as EX
from pbuf.labs.foundation import native_full_received_state_information_retention001 as RET
from pbuf.labs.foundation.native_full_state_100pct_observer_coverage_fix001 import (
    _binned_received_3d_with_launch,
)


def adapt_native_receipt(receipt: NativeReceivedState) -> dict:
    """Return the exact historical observer primitives plus native metadata."""
    launch = np.asarray(receipt.source_positions, dtype=np.float64)
    received = np.asarray(receipt.received_positions, dtype=np.float64)
    direction = np.asarray(receipt.directions, dtype=np.float64)
    snap = {
        "x": received[:, 1], "y": received[:, 2], "z": received[:, 0],
        "vx": direction[:, 1], "vy": direction[:, 2], "vz": direction[:, 0],
    }
    screen = EX._screen_coordinates(launch[:, 1], launch[:, 2], snap)
    return {
        "launch_coordinates_3d": launch,
        "received_position_3d": received,
        "direction_3d": direction,
        "deposition_weight": np.asarray(receipt.weights, dtype=np.float64),
        "source_lineage": np.asarray(receipt.native_cell_ids),
        "progression_step": np.asarray(receipt.progression_steps),
        "finite_support_state": {
            "local_displacement": np.asarray(receipt.local_displacement),
            "local_momentum": np.asarray(receipt.local_momentum),
            "local_flux": np.asarray(receipt.local_flux),
            "local_content_candidates": np.asarray(receipt.local_content_candidates),
        },
        "screen": screen,
        "snapshot": snap,
    }


def execute_frozen_observer(adapted: dict, *, bins: int = 6) -> tuple[dict, dict]:
    """Execute the unchanged 24 two-dimensional and 21 3D observer channels."""
    launch = adapted["launch_coordinates_3d"]
    screen = adapted["screen"]
    values = np.concatenate((screen["u0"], screen["v0"], screen["uf"], screen["vf"]))
    extent = max(float(np.max(np.abs(values))) * 1.001, 1.0)
    extracted = EX._extract_all(screen, extent, bins)
    receipt3d = _binned_received_3d_with_launch(
        screen, adapted["snapshot"], extent, bins, launch[:, 1], launch[:, 2]
    )
    bank, family = RET._decoded_bank(extracted, receipt3d)
    if len(bank) != 45:
        raise RuntimeError(f"frozen observer produced {len(bank)} channels, expected 45")
    return bank, {"family": family, "extent": extent, "bins": bins,
                  "primary_channel": "histogram_density__convergence"}
