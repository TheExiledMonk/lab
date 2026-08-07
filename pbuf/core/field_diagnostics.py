"""M11 — Field Diagnostics.

Provides fingerprint, statistics, and non-triviality assertions for
any field passed between modules.

Second-review repair FOUNDATION-001-CORRECTION-002
--------------------------------------------------
* ``allow_zero`` is implemented without bypassing nonfinite checks.
* ``array_fingerprint`` exposes both a raw artifact SHA-256 and a
  canonical float64 numerical SHA-256.
* Vector statistics reject component shape mismatches.
* Vector global sums are now STRICT: if any component contains NaN or
  Inf anywhere, ``field_is_finite`` is False and the reported global
  vector sum/norm are NaN.  No partial finite-mask closure can be
  mistaken for a valid conservation result.
"""
from __future__ import annotations
import hashlib
import numpy as np

from .conventions import EPS_VARIANCE_UNDEFINED

__all__ = [
    "field_statistics_scalar", "field_statistics_vector",
    "array_fingerprint", "assert_nontrivial_field", "assert_finite_field",
    "assert_not_stale", "FieldArtifact", "FieldDiagnosticsError",
    "TrivialFieldError", "NonFiniteFieldError",
]


class FieldDiagnosticsError(ValueError):
    pass


class TrivialFieldError(FieldDiagnosticsError):
    """Raised when a field has near-zero variance and is treated as trivial."""


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
    """Hash dtype, shape, and raw contiguous bytes together."""
    arr = np.ascontiguousarray(arr)
    payload = arr.dtype.str.encode("utf-8")
    payload += str(arr.shape).encode("utf-8")
    payload += arr.tobytes()
    return hashlib.sha256(payload).hexdigest()


def _hash_canonical_float64(arr: np.ndarray) -> str:
    """Hash canonical float64 numerical content only."""
    arr64 = np.ascontiguousarray(np.asarray(arr, dtype=np.float64))
    return hashlib.sha256(arr64.tobytes()).hexdigest()


def array_fingerprint(arr):
    """Return deterministic raw and canonical fingerprints for an array."""
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
    """Return descriptive statistics for a scalar field.

    Statistics are computed over finite entries only, while nonfinite
    counts are always reported explicitly.  ``field_is_finite`` tells
    callers whether the statistics represent the complete field.
    """
    arr = np.asarray(field, dtype=np.float64)
    nan_count = int(np.count_nonzero(np.isnan(arr)))
    inf_count = int(np.count_nonzero(np.isinf(arr)))
    finite_mask = np.isfinite(arr)
    finite_arr = arr[finite_mask]
    n_finite = int(finite_mask.sum())
    field_is_finite = bool(np.all(finite_mask))

    if finite_arr.size == 0:
        return {
            "shape": list(arr.shape),
            "dtype": str(arr.dtype),
            "minimum": float("nan"),
            "maximum": float("nan"),
            "mean": float("nan"),
            "variance": float("nan"),
            "rms": float("nan"),
            "nonzero_count": 0,
            "finite_count": 0,
            "nan_count": nan_count,
            "inf_count": inf_count,
            "field_is_finite": False,
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
        "field_is_finite": field_is_finite,
    }


def field_statistics_vector(Rx, Ry, Rz):
    """Return statistics for a vector field.

    Raises ``FieldDiagnosticsError`` if component shapes differ.

    Conservation/closure semantics are strict: if *any* component has
    any NaN or Inf entry, the field cannot be treated as globally valid.
    In that case ``field_is_finite`` is False and both
    ``global_vector_sum`` and ``global_vector_sum_norm`` are NaN.  The
    routine never masks bad entries and reports a plausible-looking
    partial sum.
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

    field_is_finite = bool(
        stats_x["field_is_finite"]
        and stats_y["field_is_finite"]
        and stats_z["field_is_finite"]
    )

    if field_is_finite:
        sum_x = float(np.sum(Rx))
        sum_y = float(np.sum(Ry))
        sum_z = float(np.sum(Rz))
        sum_vec = (sum_x, sum_y, sum_z)
        sum_norm = float(np.sqrt(sum_x ** 2 + sum_y ** 2 + sum_z ** 2))
    else:
        nan = float("nan")
        sum_vec = (nan, nan, nan)
        sum_norm = nan

    return {
        "shape": list(stats_x["shape"]),
        "Rx": stats_x,
        "Ry": stats_y,
        "Rz": stats_z,
        "magnitude": mag_stats,
        "field_is_finite": field_is_finite,
        "global_vector_sum": sum_vec,
        "global_vector_sum_norm": sum_norm,
    }


def assert_nontrivial_field(field, name="field", variance_epsilon=None,
                              allow_zero=False):
    """Assert that a scalar field is finite and nontrivial.

    ``allow_zero=True`` permits an exactly zero field but never permits
    NaN or Inf.  A constant nonzero field still has zero variance and is
    therefore trivial under this scalar-field contract.
    """
    if variance_epsilon is None:
        variance_epsilon = EPS_VARIANCE_UNDEFINED
    arr = np.asarray(field, dtype=np.float64)
    if not np.all(np.isfinite(arr)):
        n_nan = int(np.count_nonzero(np.isnan(arr)))
        n_inf = int(np.count_nonzero(np.isinf(arr)))
        raise NonFiniteFieldError(
            f"{name}: contains {n_nan} NaN and {n_inf} Inf entries.")
    if allow_zero and np.all(arr == 0):
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
    """Assert a FieldArtifact's sha256 and optional candidate label."""
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
# Self-check / unit tests
# ----------------------------------------------------------------------
def _fingerprint_test():
    A = np.random.RandomState(0).randn(4, 5, 6)
    fp = array_fingerprint(A)
    assert fp["shape"] == [4, 5, 6]
    assert fp["raw_sha256"] == _hash_raw(A)
    assert fp["canonical_float64_sha256"] == _hash_canonical_float64(A)

    A32 = A.astype(np.float32)
    fp32 = array_fingerprint(A32)
    assert fp32["raw_sha256"] != fp["raw_sha256"]

    B = np.array([[0.5, 1.0, 1.5], [2.0, 0.25, -1.0]], dtype=np.float32)
    B64 = B.astype(np.float64)
    assert (array_fingerprint(B)["canonical_float64_sha256"] ==
            array_fingerprint(B64)["canonical_float64_sha256"])
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

    assert assert_nontrivial_field(
        np.zeros((3, 4, 5)), "zeros", variance_epsilon=1e-15,
        allow_zero=True)

    for allow_zero in (False, True):
        try:
            assert_nontrivial_field(
                np.array([[float("nan")]]), "nan", allow_zero=allow_zero)
        except NonFiniteFieldError:
            pass
        else:
            raise AssertionError("nonfinite field escaped validation")


def _vector_shape_mismatch_test():
    try:
        field_statistics_vector(
            np.zeros((3,)), np.zeros((4,)), np.zeros((3,)))
    except FieldDiagnosticsError:
        return {"passes": True}
    return {"passes": False}


def _vector_finite_sum_test():
    """Finite vectors must report the exact full-field global sum."""
    Rx = np.array([1.0, 2.0, 3.0])
    Ry = np.array([4.0, 5.0, 6.0])
    Rz = np.array([-1.0, 0.0, 1.0])
    s = field_statistics_vector(Rx, Ry, Rz)
    expected = (6.0, 15.0, 0.0)
    return {
        "passes": (
            s["field_is_finite"] is True
            and s["global_vector_sum"] == expected
            and np.isfinite(s["global_vector_sum_norm"])
        )
    }


def _vector_nonfinite_sum_test():
    """Any NaN/Inf must invalidate the complete global vector sum.

    This specifically guards against the predecessor behavior that
    masked bad entries and returned a partial sum that could falsely
    appear to satisfy closure.
    """
    fixtures = [
        (
            np.array([1.0, float("nan"), 2.0]),
            np.array([3.0, 4.0, 5.0]),
            np.array([0.0, 1.0, 2.0]),
            "nan_in_Rx",
        ),
        (
            np.array([1.0, 2.0, 3.0]),
            np.array([3.0, float("inf"), 5.0]),
            np.array([0.0, 1.0, 2.0]),
            "inf_in_Ry",
        ),
    ]
    rows = []
    for Rx, Ry, Rz, label in fixtures:
        s = field_statistics_vector(Rx, Ry, Rz)
        vec = s["global_vector_sum"]
        passes = (
            s["field_is_finite"] is False
            and all(np.isnan(v) for v in vec)
            and np.isnan(s["global_vector_sum_norm"])
        )
        rows.append({"fixture": label, "passes": passes})
    return rows


def _scalar_all_nonfinite_test():
    """An all-nonfinite scalar field must not report fake zero stats."""
    s = field_statistics_scalar(np.array([float("nan"), float("inf")]))
    return {
        "passes": (
            s["field_is_finite"] is False
            and s["finite_count"] == 0
            and np.isnan(s["mean"])
            and np.isnan(s["variance"])
            and np.isnan(s["rms"])
        )
    }


if __name__ == "__main__":
    fp = _fingerprint_test()
    print(f"M11 fingerprint: shape={fp['shape']}, "
          f"raw_sha={fp['raw_sha256'][:16]}, "
          f"can_sha={fp['canonical_float64_sha256'][:16]}")

    _assertions_test()
    print("M11 assertions: trivial/nonfinite/allow_zero contract passed")

    r = _vector_shape_mismatch_test()
    assert r["passes"]
    print("M11 vector shape mismatch: caught")

    r = _vector_finite_sum_test()
    assert r["passes"]
    print("M11 finite vector global sum: exact full-field sum")

    rows = _vector_nonfinite_sum_test()
    assert all(r["passes"] for r in rows)
    print("M11 nonfinite vector global sum: invalidated as NaN")

    r = _scalar_all_nonfinite_test()
    assert r["passes"]
    print("M11 all-nonfinite scalar statistics: NaN, not fake zeros")

    print("M11 field diagnostics: all checks passed")
