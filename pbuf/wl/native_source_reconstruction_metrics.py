"""Truth-scoring and aggregation metrics for the Dev139 blind sweep.

This module is deliberately separate from reconstruction.  Functions accepting
truth are only called after the prediction manifest has been durably frozen.
"""
from __future__ import annotations

from collections import defaultdict
import math
import numpy as np


def ambiguity_area(surface, threshold_fraction=.05):
    q = np.asarray(surface, dtype=float)
    if q.size == 0 or not np.isfinite(q).any():
        return math.nan
    cut = np.nanmin(q) + threshold_fraction * (np.nanmax(q) - np.nanmin(q))
    return float(np.mean(q <= cut))


def information_gain(area, baseline_area):
    return float(1.0 - area / baseline_area) if baseline_area > 0 else math.nan


def fractional_error(prediction, truth):
    return float(abs(prediction - truth) / abs(truth)) if truth != 0 else math.nan


def depth_error(prediction, truth, lens_depth=1.0):
    return float(abs(prediction - truth) / abs(truth - lens_depth))


def coefficient_of_variation(values):
    values = np.asarray(values, dtype=float)
    mean = float(np.mean(values))
    return float(np.std(values) / abs(mean)) if mean else math.nan


def score_depth(prediction, truth, lens_depth=1.0):
    candidates = tuple(float(x) for x in prediction.get("depth_candidates", ()))
    support = prediction.get("support_interval")
    err = depth_error(float(prediction["primary_depth"]), truth, lens_depth)
    contains = any(abs(x - truth) / abs(truth - lens_depth) <= .10 for x in candidates)
    unique = len(candidates) == 1
    false_unique = bool(unique and support and not support[0] <= truth <= support[1])
    if unique and err <= .10 and not false_unique:
        classification = "UNIQUE_DEPTH_SUCCESS"
    elif not unique and contains:
        classification = "CORRECT_MULTIVALUED"
    elif false_unique:
        classification = "FALSE_UNIQUE_DEPTH"
    else:
        classification = "FAIL"
    return {"depth_error": err, "absolute_depth_error": abs(float(prediction["primary_depth"])-truth),
            "classification": classification, "unique_success": classification == "UNIQUE_DEPTH_SUCCESS",
            "correct_multivalued": classification == "CORRECT_MULTIVALUED", "false_unique": false_unique}


def hierarchical_summary(rows, value, families=("morphology", "lens_family", "source_size", "source_depth", "response_regime")):
    """Equal-cell summary: duplicating a row within one cell cannot reweight peers."""
    cells = defaultdict(list)
    for row in rows:
        cells[tuple(row.get(k) for k in families)].append(float(row[value]))
    cell_values = np.array([np.median(v) for v in cells.values()], dtype=float)
    if not len(cell_values):
        return {"median": math.nan, "mad": math.nan, "p10": math.nan, "p90": math.nan, "cell_count": 0}
    med = float(np.median(cell_values))
    return {"median": med, "mad": float(np.median(np.abs(cell_values-med))),
            "p10": float(np.percentile(cell_values, 10)), "p90": float(np.percentile(cell_values, 90)),
            "cell_count": int(len(cell_values))}


def outcome_gates(summary):
    depth = summary["unique_success_rate"] + summary["correct_multivalued_rate"]
    gates = {
        "relative_depth_established": depth >= .70 and summary["median_depth_error"] <= .10 and summary["false_unique_rate"] <= .10,
        "depth_size_degeneracy_reduced": summary["C4_ambiguity_area"] <= .5*summary["C1_ambiguity_area"] and summary["triplet_reduced_fraction"] >= .70 and summary["C4_triplet_accuracy"]-summary["C1_triplet_accuracy"] >= .20,
        "scale_free_geometry_established": summary["distance_ratio_error"] <= .10 and summary["size_ratio_error"] <= .10 and summary["coordinate_rescaling_CV"] <= .05 and summary["resolution_CV"] <= .10,
        "roundtrip_reconstruction_established": summary["rich_roundtrip"] < summary["position_roundtrip"] and summary["rich_roundtrip"] < summary["straight_roundtrip"],
    }
    gates["established"] = all(gates.values())
    return gates
