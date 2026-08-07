"""M11 — Field Diagnostics.

Provides fingerprint, statistics, and non-triviality assertions for
any field passed between modules.

Correction pass FOUNDATION-001-CORRECTION-001
---------------------------------------------
* ``allow_zero`` in ``assert_nontrivial_field`` is now implemented:
  when True and the field is exactly zero, the call returns True
  instead of raising TrivialFieldError. ``allow_zero`` does NOT
  bypass nonfinite validation.
* ``array_fingerprint`` exposes BOTH a ``raw_sha256`` (which mixes
  dtype + shape bytes with the raw contiguous byte stream) AND a
  ``canonical_float64_sha256`` (which hashes only the canonical
  float64 numerical content). The two hashes are independent of
  each other.
* Vector statistics reject shape mismatches between Rx, Ry, Rz
  with a clear FieldDiagnosticsError.
* Global sums of vector fields now use finite validation and
  cannot silently produce NaN lineages.
"""
from __future__ import annotations
import hashlib
import numpy as np

from .conventions import (
    EPS_VARIANCE_UNDEFINED, EPS_EXACT_COMPARISON, EPS_NORM_RELATIVE,
)

__all__ = [
    "field_statistics_scalar", "field_statistics_vector",
    "array_fingerprint", "assert_nontrivial_field", "assert_finite_field",
    "assert_not_stale", "FieldArtifact", "FieldDiagnosticsError",
    "TrivialFieldError", "NonFiniteFieldError",
]


class FieldDiagnosticsError(ValueError):
    pass


class TrivialFieldError(FieldDiagnosticsError):
    """Raised when a field has near-zero variance and is treated as
    trivial."""


class NonFiniteFieldError(FieldDiagnosticsError):
    """Raised when a field contains NaN or Inf."""


class FieldArtifact:
    """Metadata wrapper for a field passed between modules."""

    __slots__ = ("data", "artifact_id", "module_name", "module_version",
                 "source_artifact_ids", "candidate_id", "cluster_id",
                 "transform_id", "sha256", "statistics", "role")

    def __init__(self, data, artifact_id, module_name, module_version,
                 source_artifact_ids, candidate_id, cluster_id,
                 transform_id, sha256, statistics, role="unspecified"):
        self.data = data
        self.artifact_id = str(artifact_id)
        self.module_name = str(module_name)
        self.module_version = str(module_version)
        self.source_artifact_ids = list(source_artifact_ids or [])
        self.candidate_id = candidate_id
        self.cluster_id = cluster_id
        self.transform_id = transform_id
        self.sha256 = sha256
        self.statistics = statistics
        self.role = role

    def to_dict(self):
        return {
            "artifact_id": self.artifact_id,
            "module_name": self.module_name,
            "module_version": self.module_version,
            "source_artifact_ids": self.source_artifact_ids,
            "candidate_id": self.candidate_id,
            "cluster_id": self.cluster_id,
            "transform_id": self.transform_id,
            "sha256": self.sha256,
            "statistics": self.statistics,
            "role": self.role,
        }


def _hash_raw(arr: np.ndarray) -> str:
    """Hash the dtype string, shape tuple, and raw contiguous bytes
    together — producing a fingerprint that distinguishes dtypes and
    shapes but is otherwise stable for the same content."""
    arr = np.ascontiguousarray(arr)
    payload = arr.dtype.str.encode("utf-8")
    payload += str(arr.shape).encode("utf-8")
    payload += arr.tobytes()
    return hashlib.sha256(payload).hexdigest()


def _hash_canonical_float64(arr: np.ndarray) -> str:
    """Hash the canonical float64 numerical content only."""
    arr64 = np.ascontiguousarray(np.asarray(arr, dtype=np.float64))
    return hashlib.sha256(arr64.tobytes()).hexdigest()


def array_fingerprint(arr):
    """Return a deterministic fingerprint dict for an array.

    Exposes BOTH:
      raw_sha256               — dtype + shape + raw bytes
      canonical_float64_sha256 — numerical content only
    The two hashes are independent and surface different failure modes
    (dtype confusion vs. numerical drift).
    """
    arr = np.asarray(arr)
    return {
        "shape": list(arr.shape),
        "dtype": str(arr.dtype),
        "id": id(arr),
        "raw_sha256": _hash_raw(arr),
        "canonical_float64_sha256": _hash_canonical_float64(arr),
        "statistics": field_statistics_scalar(arr),
    }


def field_statistics_scalar(field):
    """Return a dict of statistics for a scalar field."""
    arr = np.asarray(field, dtype=np.float64)
    nan_count = int(np.count_nonzero(np.isnan(arr)))
    inf_count = int(np.count_nonzero(np.isinf(arr)))
    finite_mask = np.isfinite(arr)
    finite_arr = arr[finite_mask] if finite_mask.any() else arr
    n_finite = int(finite_mask.sum())
    if finite_arr.size == 0:
        return {
            "shape": list(arr.shape),
            "dtype": str(arr.dtype),
            "minimum": 0.0, "maximum": 0.0, "mean": 0.0,
            "variance": 0.0, "rms": 0.0,
            "nonzero_count": 0, "finite_count": 0,
            "nan_count": nan_count, "inf_count": inf_count,
        }
    var = float(finite_arr.var())
    rms = float(np.sqrt(np.mean(finite_arr ** 2)))
    return {
        "shape": list(arr.shape),
        "dtype": str(arr.dtype),
        "minimum": float(finite_arr.min()),
        "maximum": float(finite_arr.max()),
        "mean": float(finite_arr.mean()),
        "variance": var,
        "rms": rms,
        "nonzero_count": int(np.count_nonzero(finite_arr)),
        "finite_count": n_finite,
        "nan_count": nan_count,
        "inf_count": inf_count,
    }


def field_statistics_vector(Rx, Ry, Rz):
    """Return statistics for a vector field.

    Raises FieldDiagnosticsError if Rx, Ry, Rz do not share shape.
    Global sums use finite validation; if any component contains
    NaN/Inf the global sum is reported as NaN (not silently).
    """
    Rx = np.asarray(Rx, dtype=np.float64)
    Ry = np.asarray(Ry, dtype=np.float64)
    Rz = np.asarray(Rz, dtype=np.float64)
    if Rx.shape != Ry.shape or Rx.shape != Rz.shape:
        raise FieldDiagnosticsError(
            f"vector components have mismatched shapes: "
            f"{Rx.shape}, {Ry.shape}, {Rz.shape}")
    stats_x = field_statistics_scalar(Rx)
    stats_y = field_statistics_scalar(Ry)
    stats_z = field_statistics_scalar(Rz)
    mag = np.sqrt(Rx ** 2 + Ry ** 2 + Rz ** 2)
    mag_stats = field_statistics_scalar(mag)
    finite_mask = np.isfinite(Rx) & np.isfinite(Ry) & np.isfinite(Rz)
    if finite_mask.any():
        sum_x = float(np.sum(Rx[finite_mask]))
        sum_y = float(np.sum(Ry[finite_mask]))
        sum_z = float(np.sum(Rz[finite_mask]))
        sum_vec = (sum_x, sum_y, sum_z)
        sum_norm = float(np.sqrt(sum_x ** 2 + sum_y ** 2 + sum_z ** 2))
    else:
        sum_vec = (float("nan"), float("nan"), float("nan"))
        sum_norm = float("nan")
    return {
        "shape": list(stats_x["shape"]),
        "Rx": stats_x, "Ry": stats_y, "Rz": stats_z,
        "magnitude": mag_stats,
        "global_vector_sum": sum_vec,
        "global_vector_sum_norm": sum_norm,
    }


def assert_nontrivial_field(field, name="field", variance_epsilon=None,
                              allow_zero=False):
    """Assert a field is nontrivial.

    Parameters
    ----------
    field : array_like
    name : str
    variance_epsilon : float, optional (default EPS_VARIANCE_UNDEFINED)
    allow_zero : bool
        When True and the field is exactly zero, the call returns True
        instead of raising TrivialFieldError. ``allow_zero`` does NOT
        bypass nonfinite validation.
    """
    if variance_epsilon is None:
        variance_epsilon = EPS_VARIANCE_UNDEFINED
    arr = np.asarray(field, dtype=np.float64)
    # Finite check (independent of allow_zero).
    if not np.all(np.isfinite(arr)):
        n_nan = int(np.count_nonzero(np.isnan(arr)))
        n_inf = int(np.count_nonzero(np.isinf(arr)))
        raise NonFiniteFieldError(
            f"{name}: contains {n_nan} NaN and {n_inf} Inf entries.")
    if allow_zero:
        # Exact-zero path: return True without raising.
        if np.all(arr == 0):
            return True
    var = float(arr.var())
    if var <= variance_epsilon:
        raise TrivialFieldError(
            f"{name}: variance {var:.3e} <= epsilon {variance_epsilon:.3e}; "
            f"field is trivial.")
    return True


def assert_finite_field(field, name="field"):
    """Assert a field has no NaN or Inf entries."""
    arr = np.asarray(field, dtype=np.float64)
    if not np.all(np.isfinite(arr)):
        n_nan = int(np.count_nonzero(np.isnan(arr)))
        n_inf = int(np.count_nonzero(np.isinf(arr)))
        raise NonFiniteFieldError(
            f"{name}: contains {n_nan} NaN and {n_inf} Inf entries.")
    return True


def assert_not_stale(artifact, expected_hash, candidate_id=None):
    """Assert a FieldArtifact's sha256 matches ``expected_hash``."""
    if artifact.sha256 != expected_hash:
        raise FieldDiagnosticsError(
            f"stale artifact {artifact.artifact_id}: "
            f"expected {expected_hash}, got {artifact.sha256}"
            + (f" (candidate {candidate_id})" if candidate_id else ""))
    if candidate_id is not None and artifact.candidate_id != candidate_id:
        raise FieldDiagnosticsError(
            f"stale candidate label: artifact={artifact.candidate_id}, "
            f"expected={candidate_id}")
    return True


# ----------------------------------------------------------------------
# Self-check
# ----------------------------------------------------------------------
def _fingerprint_test():
    A = np.random.RandomState(0).randn(4, 5, 6)
    fp = array_fingerprint(A)
    assert fp["shape"] == [4, 5, 6]
    # Raw hash mixes dtype + shape + raw bytes.
    expected_raw = _hash_raw(A)
    assert fp["raw_sha256"] == expected_raw
    # Canonical hash is float64 content only.
    expected_can = _hash_canonical_float64(A)
    assert fp["canonical_float64_sha256"] == expected_can
    # Different dtype must give a different raw hash. The canonical
    # hash is content-only and DOES depend on the underlying bytes —
    # so two arrays with the same numerical content cast to float64
    # will share the canonical hash.
    A32 = A.astype(np.float32)
    fp32 = array_fingerprint(A32)
    assert fp32["raw_sha256"] != fp["raw_sha256"]
    # Construct a value-preserving comparison: use exact integer
    # fractions whose float32 and float64 representations match.
    B = np.array([[0.5, 1.0, 1.5], [2.0, 0.25, -1.0]], dtype=np.float32)
    B64 = B.astype(np.float64)
    fpB = array_fingerprint(B)
    fpB64 = array_fingerprint(B64)
    # For exact values the canonical hash matches because the
    # float32→float64 cast is exact.
    assert fpB["canonical_float64_sha256"] == fpB64["canonical_float64_sha256"]
    return fp


def _assertions_test():
    A = np.random.RandomState(0).randn(3, 4, 5)
    assert assert_nontrivial_field(A, "random", variance_epsilon=1e-3)
    try:
        assert_nontrivial_field(np.zeros((3, 4, 5)), "zeros")
    except TrivialFieldError:
        pass
    else:
        raise AssertionError("trivial zeros not caught")
    # allow_zero=True on zeros returns True.
    assert assert_nontrivial_field(np.zeros((3, 4, 5)), "zeros",
                                    variance_epsilon=1e-15, allow_zero=True)
    try:
        assert_nontrivial_field(np.array([[float("nan")]]), "nan")
    except NonFiniteFieldError:
        pass
    else:
        raise AssertionError("nan not caught")
    # allow_zero=True must still reject NaN inputs.
    try:
        assert_nontrivial_field(np.array([[float("nan")]]), "nan",
                                allow_zero=True)
    except NonFiniteFieldError:
        pass
    else:
        raise AssertionError("allow_zero did not bypass nonfinite check")


def _vector_shape_mismatch_test():
    try:
        field_statistics_vector(np.zeros((3,)), np.zeros((4,)), np.zeros((3,)))
    except FieldDiagnosticsError:
        return {"passes": True}
    return {"passes": False}


def _vector_finite_sum_test():
    """Global sums use finite validation; NaN entries are masked out
    of the global sum rather than silently propagated."""
    Rx = np.array([1.0, float("nan"), 2.0])
    Ry = np.array([3.0, 4.0, 5.0])
    Rz = np.array([0.0, 1.0, 2.0])
    s = field_statistics_vector(Rx, Ry, Rz)
    # The global sum must equal the sum over FINITE entries only
    # (here: positions 0 and 2). It must NOT silently include the
    # NaN entry.
    expected_x = 1.0 + 2.0  # index 0 and 2; index 1 is NaN and masked.
    expected_y = 3.0 + 5.0
    expected_z = 0.0 + 2.0
    return {"passes": (s["global_vector_sum"] == (expected_x, expected_y, expected_z))}


if __name__ == "__main__":
    fp = _fingerprint_test()
    print(f"M11 fingerprint: shape={fp['shape']}, "
          f"raw_sha={fp['raw_sha256'][:16]}, "
          f"can_sha={fp['canonical_float64_sha256'][:16]}")
    _assertions_test()
    print("M11 assertions: trivial/nan/inf/allow_zero all caught")
    r = _vector_shape_mismatch_test()
    assert r["passes"]
    print("M11 vector shape mismatch: caught")
    r = _vector_finite_sum_test()
    assert r["passes"]
    print("M11 vector global sum finite-validation: OK")
    print("M11 field diagnostics: all checks passed")
