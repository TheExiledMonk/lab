"""Moving-source schedules and forced use of the unchanged Dev156 F03 map."""
from __future__ import annotations

import numpy as np

from pbuf.excitation.native_relational_dynamics import f03_step
from .native_source_state import NativeSourceState
from .native_source_medium_interaction import source_medium_response, stationary_response


def integer_schedule(start: tuple[int, int, int], axis: int, moves: int,
                     dwell: int = 1) -> list[tuple[int, int, int]]:
    if axis not in (0, 1, 2) or moves < 0 or dwell < 1:
        raise ValueError("invalid integer movement schedule")
    positions = []
    p = list(start)
    for _ in range(moves + 1):
        positions.extend([tuple(p)] * dwell)
        p[axis] += 1
    return positions


def forced_f03_step(q: np.ndarray, retained: np.ndarray,
                    source: NativeSourceState | None) -> tuple[np.ndarray, np.ndarray]:
    """F03 plus an external equilibrium kick; f03_step itself is frozen."""
    q1, r1 = f03_step(q, retained)
    if source is not None:
        r1 = r1 + source_medium_response(q.shape, source)
        q1 = q1 + source_medium_response(q.shape, source)
    return q1, r1


def evolve_schedule(shape: tuple[int, int, int], amplitude: float,
                    positions: list[tuple[int, int, int]], q0: np.ndarray | None = None,
                    retained0: np.ndarray | None = None) -> dict[str, np.ndarray]:
    q = np.zeros(shape) if q0 is None else np.array(q0, dtype=float, copy=True)
    retained = np.zeros(shape) if retained0 is None else np.array(retained0, dtype=float, copy=True)
    states, residuals, constraints = [], [], []
    for position in positions:
        source = NativeSourceState(position, amplitude)
        q, retained = forced_f03_step(q, retained, source)
        equilibrium = stationary_response(shape, source)
        states.append(q.copy())
        residuals.append((q - equilibrium).copy())
        constraints.append(source_medium_response(shape, source))
    return {"states": np.asarray(states), "retained": retained,
            "dynamic_residual": np.asarray(residuals), "constraints": np.asarray(constraints)}


def release(q: np.ndarray, retained: np.ndarray, steps: int) -> np.ndarray:
    states = []
    for _ in range(steps):
        q, retained = f03_step(q, retained)
        states.append(q.copy())
    return np.asarray(states)
