"""DEV199 native-only bond-state cross-event identities.

This module intentionally imports no EM/QED material.  It is a direct,
coefficient-free reading of the frozen DEV167 pair law.
"""
from __future__ import annotations
import numpy as np
from pbuf.excitation.native_vector_pair_dynamics import bounded_stress, positive_relations


def bond_state(displacement: np.ndarray) -> dict[str, np.ndarray]:
    """Positive-N6 bond strain, unit relation and DEV167 force."""
    relation = positive_relations(displacement)
    length = np.linalg.norm(relation, axis=-1)
    strain = length - 1.0
    unit = relation / length[..., None]
    stress = bounded_stress(strain)
    # Preserve DEV167's multiplication/division ordering for bitwise force
    # reconstruction, while retaining the unit relation as native state.
    return {"strain": strain, "unit": unit, "stress": stress,
            "force": stress[..., None] * relation / length[..., None]}


def four_state_cross_term(background, a_only, b_only, ab) -> dict[str, np.ndarray]:
    """Exact inclusion--exclusion term for four existing native bond states."""
    states = [bond_state(x) for x in (background, a_only, b_only, ab)]
    s0, sa, sb, sab = states
    force = sab["force"] - sa["force"] - sb["force"] + s0["force"]
    strain = sab["strain"] - sa["strain"] - sb["strain"] + s0["strain"]
    # Exact telescoping split.  The first term changes only scalar stress at
    # the background unit relation; the remainder is the orientation change.
    constitutive = (sab["stress"] - sa["stress"] - sb["stress"] + s0["stress"])[..., None] * s0["unit"]
    geometric = force - constitutive
    return {"background": s0, "A": sa, "B": sb, "AB": sab,
            "force_cross": force, "strain_cross": strain,
            "constitutive_cross": constitutive, "geometric_cross": geometric}


def sigma_prime(epsilon):
    e = np.asarray(epsilon, dtype=np.float64)
    return (1.0 + e * e) / (1.0 - e * e) ** 2


def sigma_second(epsilon):
    e = np.asarray(epsilon, dtype=np.float64)
    return 2.0 * e * (3.0 + e * e) / (1.0 - e * e) ** 3
