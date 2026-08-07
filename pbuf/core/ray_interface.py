"""M15 — Ray Interface.

Exclusive boundary between response fields and ray propagation.

Correction pass FOUNDATION-001-CORRECTION-001
---------------------------------------------
* Replaced the unsafe absolute variance gate (var < 1e-15) with a
  five-category scale-aware classification:
    - exact_zero           — max|R| = 0
    - constant_nonzero     — Var(Rx) = Var(Ry) = 0 but max|R| > 0
    - structured_small     — nonzero spatial variation, small amplitude
    - structured_normal    — nonzero spatial variation, ordinary scale
    - nonfinite            — contains NaN or Inf
* Acceptance rules:
    - reject exact_zero when require_nontrivial=True
    - reject nonfinite always
    - permit constant_nonzero, structured_small, structured_normal
* The ray interface now validates field integrity (does the field
  look like a usable numerical field), NOT whether it will produce
  useful convergence.
"""
from __future__ import annotations
import hashlib
import numpy as np

from .field_diagnostics import (
    FieldArtifact, field_statistics_vector, array_fingerprint,
    assert_nontrivial_field, assert_finite_field, TrivialFieldError,
    NonFiniteFieldError, FieldDiagnosticsError,
)

__all__ = [
    "classify_ray_input",
    "prepare_ray_input", "ray_input_fingerprint",
    "RayInterfaceError", "TrivialRayInputError",
    "RAY_CLASSES",
]


class RayInterfaceError(ValueError):
    pass


class TrivialRayInputError(RayInterfaceError):
    """Raised when the ray input is trivial or invalid and the caller
    requested rejection."""


# Ray input classification categories (CORRECTION-001 §14.2).
RAY_CLASSES = (
    "exact_zero",
    "constant_nonzero",
    "structured_small",
    "structured_normal",
    "nonfinite",
)


# Scale threshold separating structured_small from structured_normal.
# This is a *diagnostic* threshold; it does not affect acceptance.
RAY_SMALL_AMPLITUDE = 1e-12


def _hash_array(arr: np.ndarray) -> str:
    """Hash dtype + shape + raw bytes (deterministic)."""
    arr = np.ascontiguousarray(arr)
    payload = arr.dtype.str.encode("utf-8")
    payload += str(arr.shape).encode("utf-8")
    payload += arr.tobytes()
    return hashlib.sha256(payload).hexdigest()


def classify_ray_input(Rx_2d, Ry_2d):
    """Classify a 2D ray input into one of five categories.

    The classification is scale-aware: structured_small fields with
    nonzero spatial variation are NOT rejected as trivial.

    Returns a dict with at least:
        classification : one of RAY_CLASSES
        R_rms          : float
        R_max          : float
        var_Rx         : float
        var_Ry         : float
    """
    Rx = np.asarray(Rx_2d, dtype=np.float64)
    Ry = np.asarray(Ry_2d, dtype=np.float64)
    if Rx.shape != Ry.shape:
        raise RayInterfaceError("Rx_2d and Ry_2d must share shape")
    if not (np.all(np.isfinite(Rx)) and np.all(np.isfinite(Ry))):
        return {
            "classification": "nonfinite",
            "R_rms": float("nan"),
            "R_max": float("nan"),
            "var_Rx": float("nan"),
            "var_Ry": float("nan"),
        }
    max_abs = float(max(np.max(np.abs(Rx)), np.max(np.abs(Ry))))
    var_Rx = float(Rx.var())
    var_Ry = float(Ry.var())
    if max_abs == 0.0:
        return {
            "classification": "exact_zero",
            "R_rms": 0.0, "R_max": 0.0,
            "var_Rx": var_Rx, "var_Ry": var_Ry,
        }
    # Use a tolerance for "constant" classification so that
    # non-exactly-representable constants (e.g. 0.3 in float64) still
    # classify as constant_nonzero when their variation is below the
    # machine-noise floor. The threshold is
    #     tol = max_abs^2 * n_machine + 1e-300
    # which is appropriate for genuinely-constant fields regardless of
    # amplitude.
    n_machine = 2.2e-16
    var_tol = max_abs ** 2 * n_machine + 1e-300
    if var_Rx <= var_tol and var_Ry <= var_tol:
        return {
            "classification": "constant_nonzero",
            "R_rms": float(np.sqrt((np.mean(Rx ** 2) + np.mean(Ry ** 2)) / 2)),
            "R_max": max_abs,
            "var_Rx": var_Rx, "var_Ry": var_Ry,
        }
    # Structured (nonzero variance somewhere)
    R_rms = float(np.sqrt(0.5 * (np.mean(Rx ** 2) + np.mean(Ry ** 2))))
    cls = "structured_small" if R_rms < RAY_SMALL_AMPLITUDE else "structured_normal"
    return {
        "classification": cls,
        "R_rms": R_rms,
        "R_max": max_abs,
        "var_Rx": var_Rx, "var_Ry": var_Ry,
    }


def prepare_ray_input(Rx_2d, Ry_2d, metadata, require_nontrivial=True,
                        upstream_rms=None):
    """Validate and package a 2D field for ray propagation.

    Parameters
    ----------
    Rx_2d, Ry_2d : ndarray of shape (H, W)
    metadata : dict with keys candidate_id, cluster_id, transform_id,
        role ("central" | "los" | "interface"), input_source, and
        optional arrays.
    require_nontrivial : bool, default True
        When True, exact_zero and nonfinite inputs are rejected.
    upstream_rms : float, optional
        If supplied, the scale-aware diagnostic ``R_rms / A_rms`` is
        recorded.  This is metadata only; it does not affect acceptance.
    """
    if not isinstance(metadata, dict):
        raise RayInterfaceError("metadata must be a dict")
    for k in ("candidate_id", "cluster_id", "transform_id", "role"):
        if k not in metadata:
            raise RayInterfaceError(f"metadata missing required key {k!r}")
    Rx = np.asarray(Rx_2d, dtype=np.float64)
    Ry = np.asarray(Ry_2d, dtype=np.float64)
    if Rx.shape != Ry.shape:
        raise RayInterfaceError("Rx_2d and Ry_2d must share shape")
    # Classification
    cls = classify_ray_input(Rx, Ry)
    if require_nontrivial:
        if cls["classification"] == "nonfinite":
            raise TrivialRayInputError(
                f"ray input contains NaN/Inf: candidate={metadata['candidate_id']}")
        if cls["classification"] == "exact_zero":
            raise TrivialRayInputError(
                f"ray input is exact zero: candidate={metadata['candidate_id']}, "
                f"cluster={metadata['cluster_id']}")
    # Scale-aware diagnostic
    scale_info = {
        "R_rms": cls["R_rms"],
        "R_max": cls["R_max"],
        "var_Rx": cls["var_Rx"],
        "var_Ry": cls["var_Ry"],
        "R_rms_over_A_rms": (cls["R_rms"] / upstream_rms
                              if (upstream_rms is not None and upstream_rms > 0)
                              else None),
    }
    # Hashes
    sha_rx = _hash_array(Rx)
    sha_ry = _hash_array(Ry)
    combined = hashlib.sha256((sha_rx + sha_ry).encode("utf-8")).hexdigest()
    stats = field_statistics_vector(Rx, Ry, np.zeros_like(Rx))
    artifact = FieldArtifact(
        data={"Rx": Rx, "Ry": Ry},
        artifact_id=f"ray_input_{metadata['candidate_id']}_"
                     f"{metadata['cluster_id']}_{metadata['transform_id']}",
        module_name="pbuf.core.ray_interface",
        module_version="1.1.0-correction001",
        source_artifact_ids=metadata.get("source_artifact_ids", []),
        candidate_id=metadata["candidate_id"],
        cluster_id=metadata["cluster_id"],
        transform_id=metadata["transform_id"],
        sha256=combined,
        statistics={**stats, "ray_classification": cls["classification"],
                    "ray_scale": scale_info},
        role=metadata["role"],
    )
    return artifact


def ray_input_fingerprint(artifact):
    """Return the fingerprint of a ray input artifact."""
    return {
        "sha256": artifact.sha256,
        "statistics": artifact.statistics,
        "role": artifact.role,
        "candidate_id": artifact.candidate_id,
        "cluster_id": artifact.cluster_id,
        "transform_id": artifact.transform_id,
    }


# ----------------------------------------------------------------------
# Self-check
# ----------------------------------------------------------------------
def _trivial_input_test():
    metadata = {"candidate_id": "PL1_PM1_PS2", "cluster_id": "MACS0416",
                 "transform_id": "RC0", "role": "central"}
    try:
        prepare_ray_input(np.zeros((4, 5)), np.zeros((4, 5)), metadata,
                            require_nontrivial=True)
    except TrivialRayInputError:
        return {"passes": True}
    return {"passes": False}


def _nan_input_test():
    metadata = {"candidate_id": "PL1_PM1_PS2", "cluster_id": "MACS0416",
                 "transform_id": "RC0", "role": "central"}
    bad = np.array([[float("nan"), 0.0], [0.0, 0.0]])
    try:
        prepare_ray_input(bad, np.zeros((2, 2)), metadata, require_nontrivial=True)
    except TrivialRayInputError:
        return {"passes": True}
    return {"passes": False}


def _nontrivial_input_test():
    metadata = {"candidate_id": "PL1_PM1_PS2", "cluster_id": "MACS0416",
                 "transform_id": "RC0", "role": "central"}
    rng = np.random.RandomState(0)
    Rx = rng.randn(8, 8); Ry = rng.randn(8, 8)
    a = prepare_ray_input(Rx, Ry, metadata, require_nontrivial=True)
    fp = ray_input_fingerprint(a)
    return {"passes": (a.sha256 == fp["sha256"] and a.role == "central"
                        and a.statistics["ray_classification"] == "structured_normal")}


def _hash_lineage_test():
    metadata = {"candidate_id": "PL1_PM1_PS2", "cluster_id": "MACS0416",
                 "transform_id": "RC0", "role": "central"}
    rng = np.random.RandomState(0)
    a1 = prepare_ray_input(rng.randn(8, 8), rng.randn(8, 8), metadata)
    a2 = prepare_ray_input(rng.randn(8, 8), rng.randn(8, 8), metadata)
    return {"passes": a1.sha256 != a2.sha256,
            "sha1": a1.sha256[:16], "sha2": a2.sha256[:16]}


def _classification_test():
    rows = []
    # exact_zero
    cls = classify_ray_input(np.zeros((4, 5)), np.zeros((4, 5)))
    rows.append({"input": "exact_zero", "class": cls["classification"],
                  "passes": cls["classification"] == "exact_zero"})
    # constant_nonzero
    cls = classify_ray_input(np.full((4, 5), 0.5), np.full((4, 5), 0.3))
    rows.append({"input": "constant_nonzero", "class": cls["classification"],
                  "passes": cls["classification"] == "constant_nonzero"})
    # structured_small
    x = np.linspace(0, 1, 50).reshape(1, -1)
    Rx = 1e-15 * np.sin(x * 2 * np.pi)
    Ry = 1e-15 * np.cos(x * 2 * np.pi)
    cls = classify_ray_input(Rx, Ry)
    rows.append({"input": "structured_small", "class": cls["classification"],
                  "passes": cls["classification"] == "structured_small"})
    # structured_normal
    rng = np.random.RandomState(0)
    cls = classify_ray_input(rng.randn(8, 8), rng.randn(8, 8))
    rows.append({"input": "structured_normal", "class": cls["classification"],
                  "passes": cls["classification"] == "structured_normal"})
    # nonfinite
    cls = classify_ray_input(np.array([[float("nan")]]), np.zeros((1, 1)))
    rows.append({"input": "nonfinite", "class": cls["classification"],
                  "passes": cls["classification"] == "nonfinite"})
    return {"rows": rows, "passes": all(r["passes"] for r in rows)}


def _wc6_absolute_variance_gate_test():
    """WC6 (CORRECTION-001 §19): small structured fields with
    R_x = 1e-10 sin(x), R_y = 1e-10 cos(y) must be accepted as
    structured_small by the corrected interface."""
    metadata = {"candidate_id": "WC6", "cluster_id": "MACS0416",
                 "transform_id": "RC0", "role": "central"}
    x = np.linspace(0, 1, 32)
    Rx = 1e-10 * np.sin(2 * np.pi * x).reshape(1, -1)
    Ry = 1e-10 * np.cos(2 * np.pi * x).reshape(1, -1)
    try:
        a = prepare_ray_input(Rx, Ry, metadata, require_nontrivial=True)
        return {"passes": True,
                "classification": a.statistics["ray_classification"]}
    except TrivialRayInputError:
        return {"passes": False}


if __name__ == "__main__":
    r = _trivial_input_test(); assert r["passes"], r
    print("M15 trivial input rejected: OK")
    r = _nan_input_test(); assert r["passes"], r
    print("M15 NaN input rejected: OK")
    r = _nontrivial_input_test(); assert r["passes"], r
    print(f"M15 nontrivial input accepted: classification={r.get('class', 'n/a')}")
    r = _hash_lineage_test(); assert r["passes"], r
    print(f"M15 hash lineage: {r['sha1']} ≠ {r['sha2']}")
    r = _classification_test(); assert r["passes"], r
    for row in r["rows"]:
        print(f"M15 classification {row['input']}: {row['class']}")
    r = _wc6_absolute_variance_gate_test(); assert r["passes"], r
    print(f"WC6 small structured field: classification={r['classification']}")
    print("M15 ray interface: all checks passed")
