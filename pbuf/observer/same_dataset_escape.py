"""DEV192 source-elimination vocabulary; no image measurement or fitting."""
from __future__ import annotations

REQUIRED_ESCAPE_FIELDS=("candidate_id","parent_DEV191_candidate","observable","source_term","source_elimination_method","required_assumptions","required_dataset_fields","native_counterpart","exact_formula","scale_dependency","orientation_dependency","selection_dependency","PSF_dependency","population_dependency","symmetry_prediction","null_prediction","status","rejection_reason")
ELIMINATION_METHODS={"EXACT_ALGEBRAIC_CANCELLATION","SYMMETRY_EXPECTATION_CANCELLATION","ENSEMBLE_EXPECTATION_CANCELLATION","CONDITIONAL_CANCELLATION","NULL_CHANNEL_CANCELLATION","NOT_ELIMINABLE"}
STATUSES={"EXACT_SOURCE_FREE_BRIDGE","ENSEMBLE_SOURCE_CANCELLATION_CANDIDATE","NULL_TEST_CANDIDATE","RELATIVE_GEOMETRY_CANDIDATE","FIELD_STATISTIC_CANDIDATE","HIGHER_ORDER_CANDIDATE","STRUCTURAL_ONLY","SOURCE_PRIOR_REQUIRED","SOURCE_NOT_ELIMINABLE","MEASUREMENT_UNAVAILABLE","MODEL_CONTAMINATED","REJECTED","UNRESOLVED"}
def validate_escape(row):
    missing=[x for x in REQUIRED_ESCAPE_FIELDS if x not in row]
    if missing: raise ValueError(f"missing fields: {missing}")
    if row['source_elimination_method'] not in ELIMINATION_METHODS: raise ValueError('invalid elimination method')
    if row['status'] not in STATUSES: raise ValueError('invalid status')
