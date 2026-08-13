"""DEV191 observable-boundary schema and representation classification.

This is deliberately descriptive.  It contains no image values, conventional
shear calibration, source reconstruction, or mapping from an astronomical
quantity to ``TransferCentroidJacobian/v1``.
"""
from __future__ import annotations

REQUIRED_CANDIDATE_FIELDS = (
    "candidate_id", "family", "observable_name", "measurement_level",
    "direct_or_derived", "PSF_dependency", "noise_dependency",
    "intrinsic_source_dependency", "population_prior_dependency",
    "absolute_scale_dependency", "orientation_dependency", "GR_dependency",
    "LCDM_dependency", "lens_model_dependency", "calibration_dependency",
    "native_counterpart", "transformation_class",
    "coefficient_free_bridge_status", "current_dataset_availability",
    "historical_PBUF_relation", "status", "reason",
)

TRANSFORMATION_CLASSES = {
    "scalar": {"spin": "spin-0", "rotation": "invariant", "reflection": "even", "axis_swap": "invariant", "isotropic_scaling": "depends on definition"},
    "vector": {"spin": "spin-1", "rotation": "R v", "reflection": "R v", "axis_swap": "S v", "isotropic_scaling": "linear for a position"},
    "symmetric_tensor": {"spin": "mixed spin-0 + spin-2", "rotation": "R Q R^T", "reflection": "R Q R^T", "axis_swap": "S Q S^T", "isotropic_scaling": "quadratic for a second moment"},
    "tracefree_tensor": {"spin": "spin-2", "rotation": "R Q_TF R^T", "reflection": "R Q_TF R^T", "axis_swap": "S Q_TF S^T", "isotropic_scaling": "quadratic before normalization"},
    "third_moment": {"spin": "spin-1 + spin-3", "rotation": "R⊗R⊗R", "reflection": "covariant", "axis_swap": "covariant", "isotropic_scaling": "cubic"},
    "rotation": {"spin": "spin-0 pseudoscalar", "rotation": "invariant", "reflection": "sign reversal", "axis_swap": "sign reversal", "isotropic_scaling": "invariant"},
    "mixed": {"spin": "mixed", "rotation": "basis-dependent coefficient action", "reflection": "basis-dependent coefficient action", "axis_swap": "basis-dependent coefficient action", "isotropic_scaling": "basis-scale dependent"},
}

def validate_candidate(candidate: dict) -> None:
    missing = [key for key in REQUIRED_CANDIDATE_FIELDS if key not in candidate]
    if missing:
        raise ValueError(f"candidate {candidate.get('candidate_id', '<unknown>')} missing {missing}")
    if candidate["transformation_class"] not in TRANSFORMATION_CLASSES:
        raise ValueError(f"unknown transformation class {candidate['transformation_class']}")
