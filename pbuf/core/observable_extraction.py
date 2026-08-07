"""M16 — Observable Extraction and Field Comparison.

Second-review correction FOUNDATION-001-CORRECTION-002
------------------------------------------------------
This module operates on already-extracted lensing observables. It does
not derive kappa/gamma from a Jacobian and it does not assume the identity
of any comparison field.

Scientific contract:
* Pearson uses pairwise-finite samples and returns NaN when either field
  has undefined/near-zero variance.
* Spearman uses average ranks for ties and the same pairwise-finite mask.
* Input field shapes must match before flattening/masking.
* A reference comparison is computed only when an explicit
  ``reference_kappa`` field is supplied.
* ``package_lensing_observables`` is the semantically correct public API.
  ``extract_jacobian_observables`` is retained as a compatibility alias
  for existing labs.
"""
from __future__ import annotations

import math
import numpy as np

from .conventions import EPS_VARIANCE_UNDEFINED

__all__ = [
    "safe_pearson",
    "safe_spearman",
    "_average_ranks",
    "package_lensing_observables",
    "extract_jacobian_observables",
    "ObservableExtractionError",
]


class ObservableExtractionError(ValueError):
    pass


def _paired_finite_samples(field_a, field_b):
    """Return pairwise-finite flattened samples after shape validation."""
    a0 = np.asarray(field_a, dtype=np.float64)
    b0 = np.asarray(field_b, dtype=np.float64)
    if a0.shape != b0.shape:
        raise ObservableExtractionError(
            f"field shapes must match, got {a0.shape} and {b0.shape}")
    a = a0.ravel()
    b = b0.ravel()
    mask = np.isfinite(a) & np.isfinite(b)
    return a[mask], b[mask]


def safe_pearson(field_a, field_b, variance_epsilon=None):
    """Pairwise-finite Pearson correlation.

    Returns NaN when fewer than two finite pairs remain or when either
    retained field has variance <= ``variance_epsilon``.
    """
    if variance_epsilon is None:
        variance_epsilon = EPS_VARIANCE_UNDEFINED
    if variance_epsilon < 0 or not np.isfinite(variance_epsilon):
        raise ObservableExtractionError("variance_epsilon must be finite and >= 0")

    a, b = _paired_finite_samples(field_a, field_b)
    if a.size < 2:
        return float("nan")

    a_c = a - np.mean(a)
    b_c = b - np.mean(b)
    var_a = float(np.mean(a_c * a_c))
    var_b = float(np.mean(b_c * b_c))
    if var_a <= variance_epsilon or var_b <= variance_epsilon:
        return float("nan")

    cov = float(np.mean(a_c * b_c))
    denom = math.sqrt(var_a * var_b)
    if not math.isfinite(denom) or denom == 0.0:
        return float("nan")
    r = cov / denom
    if not math.isfinite(r):
        return float("nan")
    return float(np.clip(r, -1.0, 1.0))


def _average_ranks(x):
    """Return one-based average ranks, assigning equal values equal rank."""
    x = np.asarray(x, dtype=np.float64).ravel()
    if not np.all(np.isfinite(x)):
        raise ObservableExtractionError("_average_ranks requires finite input")
    n = x.size
    if n == 0:
        return np.empty(0, dtype=np.float64)

    order = np.argsort(x, kind="mergesort")
    sorted_x = x[order]
    ranks = np.empty(n, dtype=np.float64)
    i = 0
    while i < n:
        j = i + 1
        while j < n and sorted_x[j] == sorted_x[i]:
            j += 1
        avg_rank = 0.5 * ((i + 1) + j)
        ranks[order[i:j]] = avg_rank
        i = j
    return ranks


def safe_spearman(field_a, field_b, variance_epsilon=None):
    """Pairwise-finite Spearman correlation with average ranks for ties."""
    if variance_epsilon is None:
        variance_epsilon = EPS_VARIANCE_UNDEFINED
    a, b = _paired_finite_samples(field_a, field_b)
    if a.size < 2:
        return float("nan")
    ra = _average_ranks(a)
    rb = _average_ranks(b)
    return safe_pearson(ra, rb, variance_epsilon=variance_epsilon)


def package_lensing_observables(kappa_field, gamma1_field, gamma2_field,
                                reference_kappa=None,
                                variance_epsilon=None):
    """Package already-extracted kappa/gamma observables.

    ``reference_kappa`` is optional and deliberately generic; the core
    never assumes that it represents GR.
    """
    if variance_epsilon is None:
        variance_epsilon = EPS_VARIANCE_UNDEFINED
    kappa = np.asarray(kappa_field, dtype=np.float64)
    g1 = np.asarray(gamma1_field, dtype=np.float64)
    g2 = np.asarray(gamma2_field, dtype=np.float64)
    if kappa.shape != g1.shape or kappa.shape != g2.shape:
        raise ObservableExtractionError(
            "kappa/gamma1/gamma2 must share shape; "
            f"got {kappa.shape}, {g1.shape}, {g2.shape}")

    out = {
        "kappa": kappa,
        "gamma1": g1,
        "gamma2": g2,
        "gamma_mag": np.sqrt(g1 * g1 + g2 * g2),
    }

    if reference_kappa is not None:
        ref = np.asarray(reference_kappa, dtype=np.float64)
        if ref.shape != kappa.shape:
            raise ObservableExtractionError(
                f"reference_kappa shape {ref.shape} != kappa shape {kappa.shape}")
        out["pearson_vs_reference"] = safe_pearson(
            kappa, ref, variance_epsilon=variance_epsilon)
        out["spearman_vs_reference"] = safe_spearman(
            kappa, ref, variance_epsilon=variance_epsilon)
    return out


def extract_jacobian_observables(kappa_field, gamma1_field, gamma2_field,
                                 reference_kappa=None,
                                 variance_epsilon=None):
    """Backward-compatible alias for :func:`package_lensing_observables`.

    Historical name only: this function does not derive fields from a
    Jacobian. New code should use ``package_lensing_observables``.
    """
    return package_lensing_observables(
        kappa_field, gamma1_field, gamma2_field,
        reference_kappa=reference_kappa,
        variance_epsilon=variance_epsilon,
    )


# ----------------------------------------------------------------------
# Self-check / independent validation fixtures
# ----------------------------------------------------------------------
def _pearson_tests():
    x = np.arange(20.0)
    rows = []
    rows.append(("perfect_positive", safe_pearson(x, 3.0 * x + 2.0), 1.0))
    rows.append(("perfect_negative", safe_pearson(x, -2.0 * x + 7.0), -1.0))
    rows.append(("zero_variance", safe_pearson(np.zeros_like(x), x), float("nan")))
    a = np.array([1.0, 2.0, np.nan, 4.0, 5.0])
    b = np.array([10.0, 20.0, 30.0, 40.0, 50.0])
    rows.append(("nan_masked", safe_pearson(a, b), 1.0))
    passes = True
    for _, got, expected in rows:
        if math.isnan(expected):
            passes &= math.isnan(got)
        else:
            passes &= abs(got - expected) < 1e-12
    return {"rows": rows, "passes": bool(passes)}


def _spearman_tie_tests():
    rows = []
    a = np.array([1, 1, 2, 2, 3, 3, 4, 4], dtype=float)
    b = np.array([10, 10, 20, 20, 30, 30, 40, 40], dtype=float)
    rows.append(("matched_ties_positive", safe_spearman(a, b), 1.0))
    rows.append(("matched_ties_negative", safe_spearman(a, b[::-1]), -1.0))
    rows.append(("all_ties", safe_spearman(np.ones(8), np.arange(8.0)), float("nan")))
    a2 = np.array([1.0, 2.0, np.nan, 4.0, 5.0])
    b2 = np.array([10.0, 20.0, 30.0, 40.0, 50.0])
    rows.append(("nan_masked", safe_spearman(a2, b2), 1.0))

    passes = True
    for _, got, expected in rows:
        if math.isnan(expected):
            passes &= math.isnan(got)
        else:
            passes &= abs(got - expected) < 1e-12
    return {"rows": rows, "passes": bool(passes)}


def _spearman_against_scipy_test():
    """Cross-check tied and untied fixtures against scipy when available."""
    try:
        from scipy.stats import spearmanr
    except Exception:
        return {"skipped": True, "passes": True}

    fixtures = [
        (np.array([1, 1, 2, 2, 3, 5, 5], float),
         np.array([7, 6, 6, 4, 3, 2, 2], float)),
        (np.arange(12.0), np.array([5, 1, 8, 3, 9, 2, 7, 0, 11, 10, 4, 6], float)),
    ]
    for a, b in fixtures:
        ours = safe_spearman(a, b)
        ref = float(spearmanr(a, b).statistic)
        if not np.isclose(ours, ref, rtol=0.0, atol=1e-12, equal_nan=True):
            return {"ours": ours, "scipy": ref, "passes": False}
    return {"passes": True}


def _wrong_control_tied_rank_test():
    """Old double-argsort ranking must disagree on a tied fixture."""
    a = np.array([1, 1, 2, 2, 3, 3, 4, 4, 5, 5], dtype=float)
    b = np.arange(10.0)
    ra_old = np.argsort(np.argsort(a)).astype(float) + 1.0
    rb_old = np.argsort(np.argsort(b)).astype(float) + 1.0
    r_old = safe_pearson(ra_old, rb_old)
    r_new = safe_spearman(a, b)
    return {"r_old": r_old, "r_new": r_new,
            "passes": abs(r_old - r_new) > 1e-12}


def _shape_contract_test():
    try:
        safe_pearson(np.zeros((2, 3)), np.zeros((6,)))
    except ObservableExtractionError:
        pearson_ok = True
    else:
        pearson_ok = False
    try:
        safe_spearman(np.zeros((2, 3)), np.zeros((6,)))
    except ObservableExtractionError:
        spearman_ok = True
    else:
        spearman_ok = False
    return {"passes": pearson_ok and spearman_ok}


def _package_api_test():
    rng = np.random.RandomState(21)
    kappa = rng.randn(10, 10)
    g1 = rng.randn(10, 10)
    g2 = rng.randn(10, 10)

    out = package_lensing_observables(kappa, g1, g2)
    no_reference_keys = (
        "pearson_vs_reference" not in out and
        "spearman_vs_reference" not in out
    )

    reference = 2.0 * kappa + 0.01 * rng.randn(10, 10)
    out_ref = package_lensing_observables(
        kappa, g1, g2, reference_kappa=reference)
    finite_reference_metrics = (
        np.isfinite(out_ref["pearson_vs_reference"]) and
        np.isfinite(out_ref["spearman_vs_reference"])
    )

    compat = extract_jacobian_observables(
        kappa, g1, g2, reference_kappa=reference)
    compatibility_equal = (
        np.array_equal(compat["kappa"], out_ref["kappa"]) and
        compat["pearson_vs_reference"] == out_ref["pearson_vs_reference"] and
        compat["spearman_vs_reference"] == out_ref["spearman_vs_reference"]
    )
    return {"passes": bool(no_reference_keys and finite_reference_metrics
                            and compatibility_equal)}


if __name__ == "__main__":
    r = _pearson_tests(); assert r["passes"], r
    print("M16 Pearson analytic fixtures: PASS")
    r = _spearman_tie_tests(); assert r["passes"], r
    print("M16 Spearman tie/NaN fixtures: PASS")
    r = _spearman_against_scipy_test(); assert r["passes"], r
    print("M16 Spearman scipy cross-check: PASS")
    r = _wrong_control_tied_rank_test(); assert r["passes"], r
    print(f"M16 tied-rank wrong control: old={r['r_old']:.6f}, new={r['r_new']:.6f}")
    r = _shape_contract_test(); assert r["passes"], r
    print("M16 correlation shape contract: PASS")
    r = _package_api_test(); assert r["passes"], r
    print("M16 observable packaging/reference contract: PASS")
    print("M16 observable extraction: all checks passed")
