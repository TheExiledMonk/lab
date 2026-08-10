"""Mixed loaded/excited state construction and coefficient-free progression."""
from __future__ import annotations
import numpy as np
from .native_neighbor_loaded_excitation import LOADS, EXCITATIONS, loading_profile, excitation, frames_from_loading
from .native_neighbor_frame_transport import transport


def construct_case(load_index: int, excitation_index: int, n: int = 64, amplitude: float = 1.0):
    if not 0 <= load_index < len(LOADS) or not 0 <= excitation_index < len(EXCITATIONS):
        raise IndexError("LOAD/EX index out of range")
    L = loading_profile(load_index, n)
    X = amplitude * excitation(excitation_index, n)
    return {"load_id": f"LOAD{load_index:02d}", "excitation_id": f"EX{excitation_index+1:02d}",
            "L": L, "X": X, "frames": frames_from_loading(L)}


def progress_case(case: dict, steps: int = 16, frame_candidate: str = "F04") -> dict:
    """Nearest-neighbour packet progression through loaded local frames."""
    L, frames = np.asarray(case["L"]), np.asarray(case["frames"])
    x = np.asarray(case["X"], float).copy()
    history = [np.column_stack((L, x))]
    for _ in range(int(steps)):
        nxt = np.zeros_like(x)
        for i in range(len(x)):
            j = (i + 1) % len(x)
            nxt[j] = transport(x[i], frames[i], frames[j], frame_candidate)
        x = nxt
        history.append(np.column_stack((L, x)))
    h = np.asarray(history)
    norm = np.sum(h[..., 1:] ** 2, axis=(1, 2))
    return {"history": h, "excitation_norm": norm,
            "relative_norm_drift": float(np.ptp(norm) / max(norm[0], 1e-30)),
            "loading_backreaction": 0.0, "new_interaction_coefficients": 0}


def full_matrix(n=64, steps=16, frame_candidate="F04"):
    return [(construct_case(li, ei, n), progress_case(construct_case(li, ei, n), steps, frame_candidate))
            for li in range(9) for ei in range(8)]

