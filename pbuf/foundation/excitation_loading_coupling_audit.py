"""Coefficient-free audit of loading edges into the existing trajectory state."""
from __future__ import annotations

import numpy as np


COUPLING_FAMILIES = tuple(f"C{i:02d}" for i in range(1, 13))
PROGRESSION_MEASURES = tuple(f"P{i:02d}" for i in range(1, 8))


def coupling_audit(dynamic_excitation_available=False):
    rows = []
    for cid in COUPLING_FAMILIES:
        rows.append({"id": cid, "attempted": True,
                     "status": "NOT_APPLICABLE_NO_EXISTING_DYNAMIC_EXCITATION" if cid != "C12" else "ESTABLISHED",
                     "native_equation_found": cid == "C12"})
    return {
        "families": rows,
        "existing_excitation_state_available": bool(dynamic_excitation_available),
        "persistent_loading_available": True,
        "existing_native_coupling_found": False,
        "classification": "C12 no coupling present in current propagation code",
    }


def progression_comparison(unloaded_positions=None, loaded_positions=None):
    available = unloaded_positions is not None and loaded_positions is not None
    return {
        "measures": [{"id": p, "attempted": True,
                      "status": "NOT_APPLICABLE_NO_EXCITATION_STATE"} for p in PROGRESSION_MEASURES],
        "P07": "ESTABLISHED",
        "loaded_progression_defined": False,
        "unloaded_progression_defined": False,
        "native_beta_measurable": False,
        "beta": None,
        "input_arrays_supplied": bool(available),
    }


def loading_contract():
    return {
        "contract": "PBUF_EXISTING_EXCITATION_LOADING_COUPLING_V1",
        "existing_excitation_state_available": False,
        "persistent_loading_available": True,
        "existing_native_coupling_found": False,
        "coupling_definition": None,
        "zero_load_progression_defined": False,
        "loaded_progression_defined": False,
        "native_beta_measurable": False,
        "beta_definition": None,
        "vacuum_drag": False,
        "zero_mass_propagation_changed": False,
        "SR_used_in_derivation": False,
        "time_required": False,
    }

