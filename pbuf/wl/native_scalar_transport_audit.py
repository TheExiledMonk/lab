"""Reusable numerical controls for neutral scalar candidate audits."""
from __future__ import annotations
import numpy as np
from pbuf.wl.native_zero_mass_scalar_transport import apply_factors


def coefficient_of_variation(values) -> float:
    a = np.asarray(values, float); mean = float(np.mean(a))
    return 0.0 if np.all(a == mean) else float(np.std(a) / abs(mean)) if mean else float("inf")


def uniform_identity(q0=1.0, steps=64) -> dict[str, object]:
    s = apply_factors(q0, np.ones(steps))
    ratios = np.array([1.0] + [h.q_after/q0 for h in s.history])
    return {"max_abs_ratio_error": float(np.max(np.abs(ratios-1))), "pass": bool(np.max(np.abs(ratios-1)) <= 1e-10)}


def q0_fractional_invariance(factors, q0_values=(.25,.5,1,2,4,8)) -> dict[str, object]:
    curves = np.array([[1.0]+[h.q_after/q for h in apply_factors(q, factors).history] for q in q0_values])
    cv = np.array([coefficient_of_variation(curves[:, i]) for i in range(curves.shape[1])])
    return {"q0": list(q0_values), "max_cv": float(cv.max()), "strong": bool(cv.max() <= .01), "curves": curves}


def reverse_path(factors, q0=1.0) -> dict[str, object]:
    forward = apply_factors(q0, factors); reverse = apply_factors(forward.q_receive, 1/np.asarray(factors)[::-1])
    error = abs(reverse.q_receive/q0-1)
    return {"forward_ratio": forward.q_ratio, "recovered_ratio": reverse.q_receive/q0,
            "relative_error": float(error), "pass": bool(error <= 1e-8)}


def path_memory(factors_a, factors_b, q0=1.0, *, tolerance=1e-10) -> dict[str, object]:
    ra=apply_factors(q0,factors_a).q_ratio; rb=apply_factors(q0,factors_b).q_ratio
    return {"ratio_a":ra,"ratio_b":rb,"PATH_MEMORY":bool(abs(ra-rb)>tolerance)}


def refinement_cv(state_function, densities=(.5,1,2,4)) -> dict[str, object]:
    ratios=[float(state_function(d)) for d in densities]; cv=coefficient_of_variation(ratios)
    return {"sampling_density":list(densities),"ratios":ratios,"cv":cv,
            "classification":"STRONG" if cv<=.02 else "MODERATE" if cv<=.05 else "PATH_DISCRETIZATION_UNSTABLE"}
