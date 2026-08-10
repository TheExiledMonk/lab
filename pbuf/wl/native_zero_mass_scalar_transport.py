"""Coefficient-free candidate registry and passive scalar transport machinery."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Mapping, Sequence

import numpy as np

from pbuf.wl.native_zero_mass_scalar import ScalarStep, ZeroMassScalarState


STATUSES = {"ESTABLISHED", "DERIVABLE", "STRUCTURALLY_VALID", "RELATION_ONLY",
            "MISSING_CONSTITUTIVE_COEFFICIENT", "MISSING_NATIVE_STATE", "NON_DIMENSIONLESS",
            "CIRCULAR", "ENDPOINT_ONLY", "NONREVERSIBLE_UNEXPLAINED", "Q0_DEPENDENT",
            "TRAJECTORY_CONTAMINATING", "NUMERICALLY_UNSTABLE", "REDUNDANT", "NOT_APPLICABLE"}


@dataclass(frozen=True)
class TransportCandidate:
    id: str
    name: str
    lane: str
    driver: str | None
    status: str
    provenance: str
    classification: str
    reversible: str
    dimensionless: bool
    free_coefficients: int = 0


_SPECS = [
 ("U01","identity transport","multiplicative",None,"STRUCTURALLY_VALID","scalar carriage infrastructure","IDENTITY","REVERSIBLE",True),
 ("U02","additive local-state difference","additive","u","MISSING_CONSTITUTIVE_COEFFICIENT","native u difference","LOCAL_REVERSIBLE","UNDETERMINED",False),
 ("U03","additive directional gradient","additive","directional_grad_u","MISSING_CONSTITUTIVE_COEFFICIENT","native directional gradient","PATH_ACCUMULATED","UNDETERMINED",False),
 ("U04","multiplicative local-state ratio","multiplicative","u_positive","ENDPOINT_ONLY","positive native state ratio","ENDPOINT_ONLY","REVERSIBLE",True),
 ("U05","multiplicative response-gradient transfer","multiplicative","grad_u","MISSING_CONSTITUTIVE_COEFFICIENT","native response gradient","PATH_ACCUMULATED","UNDETERMINED",False),
 ("U06","multiplicative strain-gradient transfer","multiplicative","grad_epsilon","MISSING_CONSTITUTIVE_COEFFICIENT","native strain gradient","PATH_ACCUMULATED","UNDETERMINED",False),
 ("U07","fast-channel local transfer","multiplicative","delta_u_fast","MISSING_CONSTITUTIVE_COEFFICIENT","native fast channel","LOCAL_REVERSIBLE","UNDETERMINED",False),
 ("U08","slow-channel local transfer","multiplicative","delta_u_slow","MISSING_CONSTITUTIVE_COEFFICIENT","native slow channel","LOCAL_REVERSIBLE","UNDETERMINED",False),
 ("U09","combined fast/slow local transfer","multiplicative","fast_slow","NON_DIMENSIONLESS","native terminal trajectory combination","LOCAL_REVERSIBLE","UNDETERMINED",False),
 ("U10","neighbor-state transition transfer","multiplicative","neighbor_response","MISSING_CONSTITUTIVE_COEFFICIENT","native neighbor pair","LOCAL_REVERSIBLE","UNDETERMINED",False),
 ("U11","entry/exit boundary transfer","multiplicative","boundary","MISSING_CONSTITUTIVE_COEFFICIENT","native boundary states","BOUNDARY_ONLY","UNDETERMINED",False),
 ("U12","bounded-strain energy-difference transfer","multiplicative","W_positive","ENDPOINT_ONLY","native bounded-strain W ratio","ENDPOINT_ONLY","REVERSIBLE",True),
 ("U13","convex-excitation-difference transfer","multiplicative","excitation_positive","ENDPOINT_ONLY","native positive excitation ratio","ENDPOINT_ONLY","REVERSIBLE",True),
 ("U14","trajectory-curvature-coupled transfer","additive","curvature","TRAJECTORY_CONTAMINATING","trajectory curvature only","PATH_ACCUMULATED","UNDETERMINED",False),
 ("U15","path-projected response transfer","additive","projected_response","MISSING_CONSTITUTIVE_COEFFICIENT","native projected response","PATH_ACCUMULATED","UNDETERMINED",False),
 ("U16","bundle-expansion transfer","multiplicative","bundle_area","MISSING_NATIVE_STATE","optional bundle diagnostic","PATH_ACCUMULATED","UNDETERMINED",True),
 ("U17","local Hessian transfer","additive","hessian","MISSING_CONSTITUTIVE_COEFFICIENT","native local Hessian","PATH_ACCUMULATED","UNDETERMINED",False),
 ("U18","exact-state-ratio transfer","multiplicative","positive_state","ENDPOINT_ONLY","generic positive native state ratio","ENDPOINT_ONLY","REVERSIBLE",True),
 ("U19","reversible pair-state transfer","multiplicative","pair_ratio","ENDPOINT_ONLY","positive pair-state ratio","LOCAL_REVERSIBLE","REVERSIBLE",True),
 ("U20","conservative endpoint-potential transfer","multiplicative","potential_positive","ENDPOINT_ONLY","positive endpoint-potential ratio","ENDPOINT_ONLY","REVERSIBLE",True),
]


CANDIDATES = {s[0]: TransportCandidate(*s, free_coefficients=(1 if s[4] == "MISSING_CONSTITUTIVE_COEFFICIENT" else 0)) for s in _SPECS}


def candidate_registry() -> dict[str, TransportCandidate]:
    return dict(CANDIDATES)


def ratio_transfer(x_i: float, x_j: float) -> float:
    if not np.isfinite(x_i) or not np.isfinite(x_j) or x_i <= 0.0 or x_j <= 0.0:
        raise ValueError("ratio transfer requires positive finite states")
    return float(x_j / x_i)


def inverse_pair_factor(factor: float) -> float:
    if not np.isfinite(factor) or factor == 0.0:
        raise ValueError("pair factor must be finite and nonzero")
    return 1.0 / float(factor)


def apply_factors(q0: float, factors: Sequence[float], *, candidate_id="U01",
                  path_positions: Sequence[float] | None = None) -> ZeroMassScalarState:
    state = ZeroMassScalarState(q0, candidate_id=candidate_id)
    positions = np.arange(len(factors) + 1, dtype=float) if path_positions is None else np.asarray(path_positions, float)
    if len(positions) != len(factors) + 1:
        raise ValueError("path_positions must have len(factors)+1 entries")
    for i, raw in enumerate(factors):
        factor = float(raw)
        if not np.isfinite(factor) or factor <= 0.0:
            raise ValueError("multiplicative factor must be positive and finite")
        before = state.q_receive; after = before * factor
        state.history.append(ScalarStep(i, float(positions[i+1]), before, after, factor,
                                        float(np.log(factor)), factor, {}, {}))
        state.q_scalar = after
    return state


def apply_state_ratio(q0: float, positive_state: Sequence[float], *, candidate_id="U18") -> ZeroMassScalarState:
    x = np.asarray(positive_state, float)
    return apply_factors(q0, [ratio_transfer(a, b) for a, b in zip(x[:-1], x[1:])], candidate_id=candidate_id)


def apply_additive(q0: float, increments: Sequence[float], *, candidate_id="U02") -> ZeroMassScalarState:
    state = ZeroMassScalarState(q0, candidate_id=candidate_id)
    for i, inc in enumerate(np.asarray(increments, float)):
        before = state.q_receive; after = before + float(inc)
        ratio = after / before if before else np.nan
        state.history.append(ScalarStep(i, float(i+1), before, after, ratio,
                                        float(np.log(ratio)) if ratio > 0 else np.nan, float(inc), {}, {}))
        state.q_scalar = after
    return state


def transport_on_frozen_trajectory(q0: float, positions: np.ndarray, directions: np.ndarray,
                                   factors: Sequence[float], *, candidate_id="U01"):
    """Return unchanged trajectory arrays and an independently evolved scalar."""
    p = np.asarray(positions, float).copy(); d = np.asarray(directions, float).copy()
    return p, d, apply_factors(q0, factors, candidate_id=candidate_id)


def telescoping_test(states: Sequence[float], *, rtol=1e-12) -> dict[str, object]:
    x = np.asarray(states, float)
    factors = x[1:] / x[:-1]
    product = float(np.prod(factors)); endpoint = float(x[-1] / x[0])
    return {"TELESCOPING": bool(np.isclose(product, endpoint, rtol=rtol, atol=rtol)),
            "PATH_MEMORY": False, "product": product, "endpoint_ratio": endpoint}
