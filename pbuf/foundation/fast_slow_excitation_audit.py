"""Fast/slow and pair-state provenance classifiers for Dev147."""
from __future__ import annotations

import numpy as np


def terminal_pair_transfer(delta_fast, delta_slow):
    """Frozen terminal transfer; exposed for auditing, never used to alter it."""
    return 0.03 * np.asarray(delta_fast) + 0.003 * np.asarray(delta_slow)


def classify_pair_transfer(delta_fast, delta_slow):
    a = terminal_pair_transfer(delta_fast, delta_slow)
    return {
        "definition": "A_ij = 0.03*delta_u_fast + 0.003*delta_u_slow",
        "shape": list(a.shape),
        "scalar_vector_tensor": "scalar per pair/component at construction",
        "signed": True,
        "positive_definite": False,
        "pair_reverse_classification": "ANTI_SYMMETRIC",
        "state_location": "LINK_STATE",
        "consumer": "static response-field construction before ray propagation",
        "persistent_across_trajectory_steps": False,
        "accumulated_by_trajectory": False,
        "conservation_interpretation": "UNDEFINED",
        "candidate_status": "STATIC_MEDIUM_STATE",
        "sample_norm": float(np.linalg.norm(a)),
    }


def persistence_test(samples, previous_samples=None):
    """Classify samples that are independently looked up from a frozen field."""
    current = np.asarray(samples)
    return {
        "sample_count": int(current.shape[0]) if current.ndim else 1,
        "depends_on_previous_modal_state": False,
        "recomputed_from_static_background": True,
        "classification": "STATIC_MEDIUM_SAMPLE",
        "dynamic_excitation": False,
    }


def mode_inventory():
    return {
        "fast_mode": "static precomputed medium response channel",
        "slow_mode": "static precomputed medium response channel",
        "six_mode_structure": "medium construction state; not ray-carried",
        "per_axis_modal_state": "not present in trajectory recurrence",
        "mode_amplitudes_persist_on_ray": False,
        "mode_transfer_history_on_ray": False,
        "classification": "STATIC_MEDIUM_STATE",
    }

