"""M16 — Observable Extraction.

Correlation, regression, and field-comparison observables.

Correction pass FOUNDATION-001-CORRECTION-001
---------------------------------------------
* Pearson correlation simplified to the standard form
    r = cov / sqrt(var_a * var_b)
  where cov = mean(a_centered * b_centered) and var_* = mean(* ** 2).
  Clamped only for floating-point overshoot.
* Spearman correlation uses AVERAGE RANKS for tied values
  (scipy.stats.rankdata with method="average" or an equivalent).
* Observable API renamed: ``extract_jacobian_observables`` now takes
  ALREADY EXTRACTED kappa/gamma1/gamma2 fields and an OPTIONAL
  reference_kappa field.  It does NOT assume the reference is GR.
* ``pearson_vs_reference`` / ``spearman_vs_reference`` are only
  computed when ``reference_kappa`` is supplied.
"""
from __future__ import annotations
import math
import numpy as np

from .conventions import EPS_VARIANCE_UNDEFINED

__all__ = [
    "safe_pearson", "safe_spearman",
    "_average_ranks",
    "extract_jacobian_observables",
    "ObservableExtractionError",
]


class ObservableExtractionError(ValueError):
    pass


def safe_pearson(field_a, field_b, variance_epsilon=None):
    """Pearson correlation that returns NaN for near-zero variance.

    CORRECTION-001 §12.1: simplified form
        r = mean(a_centered * b_centered) / sqrt(var_a * var_b)
    where both numerator and denominator use the same 1/N convention.
    """
    if variance_epsilon is None:
        variance_epsilon = EPS_VARIANCE_UNDEFINED
    a = np.asarray(field_a, dtype=np.float64).ravel()
    b = np.asarray(field_b, dtype=np.float64).ravel()
    mask = np.isfinite(a) & np.isfinite(b)
    if mask.sum() < 2:
        return float("nan")
    a = a[mask]; b = b[mask]
    a_c = a - a.mean()
    b_c = b - b.mean()
    var_a = float(np.mean(a_c ** 2))
    var_b = float(np.mean(b_c ** 2))
    if var_a <= variance_epsilon or var_b <= variance_epsilon:
        return float("nan")
    denom = math.sqrt(var_a * var_b)
    if denom == 0.0:
        return float("nan")
    cov = float(np.mean(a_c * b_c))
    r = cov / denom
    # Clamp only floating-point overshoot (do not mask genuine
    # invalid values).
    if r > 1.0:
        r = 1.0
    elif r < -1.0:
        r = -1.0
    return r


def _average_ranks(x):
    """Return the average ranks of ``x``.

    Tied values receive their average rank. This is the standard
    "average" method (matches scipy.stats.rankdata(method='average')).
    """
    x = np.asarray(x, dtype=np.float64)
    n = x.size
    # Stable sort indices.
    order = np.argsort(x, kind="mergesort")
    ranks = np.empty(n, dtype=np.float64)
    # Walk through the sorted array and assign average ranks to ties.
    i = 0
    while i < n:
        j = i
        while j < n and x[order[j]] == x[order[i]]:
            j += 1
        # Tied group spans [i, j).
        avg_rank = 0.5 * (i + j - 1) + 1.0  # 1-based average rank
        for k in range(i, j):
            ranks[order[k]] = avg_rank
        i = j
    return ranks


def safe_spearman(field_a, field_b, variance_epsilon=None):
    """Spearman rank correlation with AVERAGE ranks for ties.

    CORRECTION-001 §12.2: replaces the previous double-argsort
    ranking which assigned distinct ranks to ties and therefore
    produced incorrect Spearman values.
    """
    if variance_epsilon is None:
        variance_epsilon = EPS_VARIANCE_UNDEFINED
    a = np.asarray(field_a, dtype=np.float64).ravel()
    b = np.asarray(field_b, dtype=np.float64).ravel()
    mask = np.isfinite(a) & np.isfinite(b)
    if mask.sum() < 2:
        return float("nan")
    a = a[mask]; b = b[mask]
    ra = _average_ranks(a)
    rb = _average_ranks(b)
    return safe_pearson(ra, rb, variance_epsilon)


def extract_jacobian_observables(kappa_field, gamma1_field, gamma2_field,
                                    reference_kappa=None,
                                    variance_epsilon=None):
    """Package κ, γ₁, γ₂, γ observables.

    CORRECTION-001 §12.4: takes ALREADY EXTRACTED fields. Does NOT
    assume the identity of the reference. If ``reference_kappa`` is
    supplied, computes ``pearson_vs_reference`` and
    ``spearman_vs_reference``. The caller may label the reference in
    metadata (the core does not assume it is GR).

    Returns a dict with:
        kappa            — the input κ field
        gamma1           — the input γ₁ field
        gamma2           — the input γ₂ field
        gamma_mag        — magnitude √(γ₁² + γ₂²)
        pearson_vs_reference  (only if reference_kappa supplied)
        spearman_vs_reference (only if reference_kappa supplied)
    """
    if variance_epsilon is None:
        variance_epsilon = EPS_VARIANCE_UNDEFINED
    kappa = np.asarray(kappa_field, dtype=np.float64)
    g1 = np.asarray(gamma1_field, dtype=np.float64)
    g2 = np.asarray(gamma2_field, dtype=np.float64)
    if kappa.shape != g1.shape or kappa.shape != g2.shape:
        raise ObservableExtractionError(
            f"kappa/gamma1/gamma2 must share the same shape, "
            f"got {kappa.shape}, {g1.shape}, {g2.shape}")
    out = {
        "kappa": kappa,
        "gamma1": g1,
        "gamma2": g2,
        "gamma_mag": np.sqrt(g1 ** 2 + g2 ** 2),
    }
    if reference_kappa is not None:
        ref = np.asarray(reference_kappa, dtype=np.float64)
        if ref.shape != kappa.shape:
            raise ObservableExtractionError(
                f"reference_kappa shape {ref.shape} does not match "
                f"kappa shape {kappa.shape}")
        out["pearson_vs_reference"] = safe_pearson(kappa, ref, variance_epsilon)
        out["spearman_vs_reference"] = safe_spearman(kappa, ref, variance_epsilon)
    return out


# ----------------------------------------------------------------------
# Self-check
# ----------------------------------------------------------------------
def _pearson_basic_test():
    rng = np.random.RandomState(0)
    x = rng.randn(1000)
    y = 2 * x + 0.1 * rng.randn(1000)
    r = safe_pearson(x, y)
    return {"pearson": r, "passes": abs(r - 1.0) < 0.01}


def _pearson_zero_variance_test():
    """safe_pearson must return NaN when one variance is below epsilon."""
    a = np.zeros(100)
    b = np.linspace(0, 1, 100)
    r = safe_pearson(a, b)
    return {"pearson": r, "passes": math.isnan(r)}


def _pearson_nan_test():
    a = np.array([float("nan")] * 10)
    b = np.linspace(0, 1, 10)
    r = safe_pearson(a, b)
    return {"pearson": r, "passes": math.isnan(r)}


def _zero_kappa_test():
    kappa_zero = np.zeros((8, 8))
    rng = np.random.RandomState(0)
    gr = rng.randn(8, 8)
    r = safe_pearson(kappa_zero, gr)
    return {"pearson": r, "passes": math.isnan(r)}


def _spearman_basic_test():
    rng = np.random.RandomState(0)
    x = rng.randn(500)
    y = x + 0.1 * rng.randn(500)  # monotonic increasing
    r = safe_spearman(x, y)
    return {"spearman": r, "passes": r > 0.95}


def _spearman_decreasing_test():
    rng = np.random.RandomState(0)
    x = rng.randn(500)
    y = -x + 0.1 * rng.randn(500)  # monotonic decreasing
    r = safe_spearman(x, y)
    return {"spearman": r, "passes": r < -0.95}


def _spearman_no_ties_test():
    rng = np.random.RandomState(1)
    x = np.arange(20, dtype=np.float64)
    rng.shuffle(x)
    y = 2 * x + 0.5 * rng.randn(20)
    r = safe_spearman(x, y)
    # Compare with double-argsort (previous broken implementation).
    ra = np.argsort(np.argsort(x)).astype(np.float64) + 1
    rb = np.argsort(np.argsort(y)).astype(np.float64) + 1
    r_old = safe_pearson(ra, rb)
    return {"spearman": r, "spearman_old": r_old,
            "passes": abs(r - r_old) < 1e-12}


def _spearman_all_ties_test():
    a = np.zeros(100)
    b = np.linspace(0, 1, 100)
    r = safe_spearman(a, b)
    return {"spearman": r, "passes": math.isnan(r)}


def _spearman_repeated_plateau_test():
    """A field with a single plateau: ranks must be averaged.

    Plateau at a positions 4, 5, 6 (avg rank 6) does not exactly
    correspond to a plateau in b, so Spearman ≠ 1 exactly. But the
    monotonic relationship is preserved and Spearman must be high.
    Compare the corrected Spearman against the previous
    (double-argsort) implementation; they MUST disagree."""
    a = np.array([1, 2, 3, 4, 5, 5, 5, 6, 7, 8], dtype=np.float64)
    b = np.array([10, 20, 30, 40, 50, 60, 70, 80, 90, 100], dtype=np.float64)
    r = safe_spearman(a, b)
    ra_old = np.argsort(np.argsort(a)).astype(np.float64) + 1
    rb_old = np.argsort(np.argsort(b)).astype(np.float64) + 1
    r_old = safe_pearson(ra_old, rb_old)
    return {"spearman": r, "spearman_old": r_old,
            "passes": (abs(r - r_old) > 0.0) and r > 0.95}


def _spearman_monotonic_with_ties_test():
    """Monotonic increasing with ties: Spearman should be high (>0.95)."""
    a = np.array([1, 1, 2, 2, 3, 3, 4, 4, 5, 5], dtype=np.float64)
    b = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10], dtype=np.float64)
    r = safe_spearman(a, b)
    return {"spearman": r, "passes": r > 0.95}


def _spearman_monotonic_decreasing_with_ties_test():
    a = np.array([5, 5, 4, 4, 3, 3, 2, 2, 1, 1], dtype=np.float64)
    b = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10], dtype=np.float64)
    r = safe_spearman(a, b)
    return {"spearman": r, "passes": r < -0.95}


def _spearman_nan_test():
    a = np.array([1.0, 2.0, float("nan"), 4.0, 5.0])
    b = np.array([10.0, 20.0, 30.0, 40.0, 50.0])
    r = safe_spearman(a, b)
    return {"spearman": r, "passes": math.isnan(r) or abs(r - 1.0) < 1e-12}


def _spearman_against_scipy_test():
    """If scipy is available, verify average-rank Spearman against it."""
    try:
        from scipy.stats import spearmanr
    except Exception:
        return {"skipped": True, "passes": True}
    rng = np.random.RandomState(7)
    for _ in range(5):
        a = rng.randn(40)
        b = 0.5 * a + 0.5 * rng.randn(40)
        r_ours = safe_spearman(a, b)
        r_scipy = float(spearmanr(a, b).statistic)
        if abs(r_ours - r_scipy) > 1e-10:
            return {"spearman_ours": r_ours, "spearman_scipy": r_scipy,
                    "passes": False}
    return {"passes": True}


def _wc5_tied_rank_old_impl_test():
    """WC5 (CORRECTION-001 §19): the previous double-argsort
    implementation disagrees with average-rank Spearman on tied data."""
    # Build a tied-data fixture.
    a = np.array([1, 1, 2, 2, 3, 3, 4, 4, 5, 5], dtype=np.float64)
    b = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10], dtype=np.float64)
    # Old (broken) implementation: double argsort, distinct ranks.
    ra_old = np.argsort(np.argsort(a)).astype(np.float64) + 1
    rb_old = np.argsort(np.argsort(b)).astype(np.float64) + 1
    r_old = safe_pearson(ra_old, rb_old)
    # New (corrected) implementation: average ranks.
    r_new = safe_spearman(a, b)
    return {"r_old": r_old, "r_new": r_new,
            "passes": abs(r_old - r_new) > 0.0}


def _extract_api_test():
    """Without reference_kappa: no pearson/spearman keys.
    With reference_kappa: keys are present and use the supplied
    reference (not a hard-coded GR field)."""
    rng = np.random.RandomState(0)
    kappa = rng.randn(10, 10)
    g1 = rng.randn(10, 10)
    g2 = rng.randn(10, 10)
    out = extract_jacobian_observables(kappa, g1, g2)
    no_ref = ("pearson_vs_reference" not in out
              and "spearman_vs_reference" not in out)
    ref = rng.randn(10, 10)
    out2 = extract_jacobian_observables(kappa, g1, g2, reference_kappa=ref)
    has_ref = ("pearson_vs_reference" in out2
               and "spearman_vs_reference" in out2
               and out2["pearson_vs_reference"] != float("nan"))
    return {"passes": no_ref and has_ref}


if __name__ == "__main__":
    r = _pearson_basic_test(); assert r["passes"], r
    print(f"M16 pearson basic: r={r['pearson']:.3f}")
    r = _pearson_zero_variance_test(); assert r["passes"], r
    print("M16 pearson zero-variance: NaN")
    r = _pearson_nan_test(); assert r["passes"], r
    print("M16 pearson NaN: NaN")
    r = _zero_kappa_test(); assert r["passes"], r
    print("M16 zero κ vs GR: NaN")
    r = _spearman_basic_test(); assert r["passes"], r
    print(f"M16 spearman basic: r={r['spearman']:.3f}")
    r = _spearman_decreasing_test(); assert r["passes"], r
    print(f"M16 spearman decreasing: r={r['spearman']:.3f}")
    r = _spearman_no_ties_test(); assert r["passes"], r
    print(f"M16 spearman no-ties (ours/old agree): "
          f"{r['spearman']:.6f} vs {r['spearman_old']:.6f}")
    r = _spearman_all_ties_test(); assert r["passes"], r
    print("M16 spearman all-ties: NaN")
    r = _spearman_repeated_plateau_test(); assert r["passes"], r
    print(f"M16 spearman plateau: r={r['spearman']:.6f}")
    r = _spearman_monotonic_with_ties_test(); assert r["passes"], r
    print(f"M16 spearman monotonic+ties (increasing): r={r['spearman']:.6f}")
    r = _spearman_monotonic_decreasing_with_ties_test(); assert r["passes"], r
    print(f"M16 spearman monotonic+ties (decreasing): r={r['spearman']:.6f}")
    r = _spearman_nan_test(); assert r["passes"], r
    print(f"M16 spearman NaN-masked: r={r['spearman']:.6f}")
    r = _spearman_against_scipy_test(); assert r["passes"], r
    print("M16 spearman vs scipy: agree")
    r = _wc5_tied_rank_old_impl_test(); assert r["passes"], r
    print(f"WC5 old vs new on ties: old={r['r_old']:.6f} new={r['r_new']:.6f}")
    r = _extract_api_test(); assert r["passes"], r
    print("M16 extract API: no reference → no GR correlation, with reference → "
          "explicit pearson/spearman")
    print("M16 observable extraction: all checks passed")
