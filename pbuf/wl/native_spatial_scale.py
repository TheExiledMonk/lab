"""Native spatial-scale candidates and scale-invariance diagnostics.

This module is deliberately independent of observational arrays and does not
alter propagation.  It supplies provenance checks used by Dev136 and tests.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Iterable
import numpy as np

PHYSICAL_LENGTH_UNITS = {"m", "meter", "metre", "km", "pc", "kpc", "Mpc"}
ANGULAR_UNITS = {"rad", "radian", "arcsec", "degree"}
FORBIDDEN = {"rmax": "FORBIDDEN_RMAX", "strength=0.18": "FORBIDDEN_0P18",
             "planck": "FORBIDDEN_PLANCK_ASSUMPTION", "lcdm": "FORBIDDEN_LCDM"}


@dataclass(frozen=True)
class PhysicalScaleCandidate:
    candidate_id: str
    source_stage: str
    source_quantity: str
    value: float | None
    units: str
    derivation: str
    assumptions: tuple[str, ...] = ()
    target_dependency: bool = False
    dimensional_validity: bool = False
    coordinate_lineage_validity: bool = False
    status: str = "UNRESOLVED"
    rejection_reason: str | None = None

    def to_dict(self) -> dict:
        value = asdict(self); value["assumptions"] = list(self.assumptions); return value


def classify_candidate(candidate: PhysicalScaleCandidate) -> PhysicalScaleCandidate:
    text = " ".join((candidate.source_stage, candidate.source_quantity,
                     candidate.derivation, *candidate.assumptions)).lower()
    reason = next((reason for token, reason in FORBIDDEN.items() if token in text), None)
    if reason:
        status = "REJECTED"
    elif candidate.target_dependency:
        status, reason = "TARGET_CONTAMINATED", "TARGET_CONTAMINATED"
    elif candidate.units not in PHYSICAL_LENGTH_UNITS | ANGULAR_UNITS:
        status, reason = "NUMERICAL_ONLY", "NUMERICAL_ONLY"
    elif not candidate.dimensional_validity:
        status, reason = "DIMENSIONALLY_INVALID", "DIMENSION_MISMATCH"
    elif not candidate.coordinate_lineage_validity:
        status, reason = "PROVENANCE_INCOMPLETE", "BROKEN_LINEAGE"
    else:
        status, reason = "AUTHORITATIVE", None
    data = candidate.to_dict(); data.update(status=status, rejection_reason=reason)
    data["assumptions"] = tuple(data["assumptions"])
    return PhysicalScaleCandidate(**data)


def recover_scale(spacing: float, units: str, transforms: Iterable[float],
                  lineage_complete: bool = True) -> PhysicalScaleCandidate:
    factor = float(np.prod(tuple(transforms)))
    raw = PhysicalScaleCandidate("synthetic_spacing", "source_grid", "cell_spacing",
        float(spacing) / factor, units, f"source spacing / coordinate factor {factor:g}",
        (), False, units in PHYSICAL_LENGTH_UNITS, lineage_complete)
    return classify_candidate(raw)


def scale_invariance_control(alphas=(.5, 1., 2., 4.)) -> dict:
    """Controlled affine-ray rescaling; compares dimensionless observables."""
    t = np.linspace(0., 1., 33)
    base = np.column_stack((t, .08*np.sin(2*np.pi*t), t*t*.1))
    rows = []
    reference = None
    for alpha in alphas:
        path = base * float(alpha); delta = np.diff(path, axis=0)
        lengths = np.linalg.norm(delta, axis=1); directions = delta / lengths[:, None]
        endpoint = (path[-1]-path[0]) / np.linalg.norm(path[-1]-path[0])
        turns = np.linalg.norm(np.diff(directions, axis=0), axis=1)
        covariance = np.cov((path/np.ptp(path, axis=0)).T)
        metrics = np.r_[endpoint, turns.mean(), turns.max(), lengths/lengths.sum(),
                        np.linalg.eigvalsh(covariance)]
        if reference is None: reference = metrics
        error = float(np.max(np.abs(metrics-reference)))
        rows.append({"alpha": float(alpha), "normalized_endpoint_morphology": endpoint.tolist(),
                     "relative_direction_morphology": float(turns.mean()),
                     "bundle_anisotropy": float(np.linalg.cond(covariance)),
                     "relative_path_ratios": (lengths/lengths.sum()).tolist(),
                     "condition_number_structure": float(np.linalg.cond(covariance)),
                     "max_dimensionless_error": error})
    passed = all(r["max_dimensionless_error"] < 1e-12 for r in rows)
    return {"alphas": rows, "tolerance": 1e-12, "passed": passed,
            "outcome": "GLOBAL_SCALE_DEGENERACY_SUPPORTED" if passed else "GLOBAL_SCALE_DEGENERACY_REJECTED",
            "propagation_physics_modified": False}
