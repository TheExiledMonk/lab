"""Exact vector tangent/remainder diagnostics for the frozen DEV167 law."""
from __future__ import annotations
import numpy as np
from pbuf.excitation.native_vector_pair_dynamics import positive_relations, bounded_stress
from pbuf.observer.local_state_cross_event import sigma_prime


def positive_force(displacement: np.ndarray) -> np.ndarray:
    r = positive_relations(displacement); length = np.linalg.norm(r, axis=-1)
    return bounded_stress(length - 1.0)[..., None] * r / length[..., None]


def tangent_force(background: np.ndarray, perturbed: np.ndarray) -> np.ndarray:
    """Full geometric Frechet tangent of sigma(|r|-1) rhat on positive bonds."""
    r0 = positive_relations(background); r1 = positive_relations(perturbed)
    length = np.linalg.norm(r0, axis=-1); unit = r0 / length[..., None]
    eps = length - 1.0; stress = bounded_stress(eps); dr = r1-r0
    radial = np.sum(unit*dr, axis=-1)
    transverse = dr-unit*radial[..., None]
    return (stress[..., None]*unit + sigma_prime(eps)[..., None]*radial[..., None]*unit
            + (stress/length)[..., None]*transverse)


def summary(values: np.ndarray) -> dict:
    a=np.asarray(values); a=a[np.isfinite(a)]
    if not a.size: return {"count":0,"min":None,"median":None,"mean":None,"max":None,"quantiles":{},"l2":0.0}
    return {"count":int(a.size),"min":float(a.min()),"median":float(np.median(a)),"mean":float(a.mean()),"max":float(a.max()),"quantiles":{str(q):float(np.quantile(a,q)) for q in (.01,.1,.25,.75,.9,.99)},"l2":float(np.linalg.norm(a))}
