"""M15 — Ray Interface.

Exclusive boundary between response fields and ray propagation.

FOUNDATION-001 second-review correction
---------------------------------------
This module classifies 2D ray-input vector fields without conflating
small amplitude with numerical triviality.

Classes:
  * exact_zero
  * constant_nonzero
  * structured_small
  * structured_normal
  * nonfinite

Scientific contract:
  * nonfinite inputs are ALWAYS rejected by ``prepare_ray_input``;
  * exact-zero inputs are rejected only when ``require_nontrivial=True``;
  * constant nonzero and structured inputs are accepted;
  * the small/normal threshold is diagnostic only and never controls
    acceptance;
  * only two-dimensional image-plane fields are valid ray inputs.
"""
from __future__ import annotations

import hashlib
import numpy as np

from .field_diagnostics import FieldArtifact, field_statistics_vector

__all__ = [
    "classify_ray_input",
    "prepare_ray_input",
    "ray_input_fingerprint",
    "RayInterfaceError",
    "TrivialRayInputError",
    "RAY_CLASSES",
]


class RayInterfaceError(ValueError):
    pass


class TrivialRayInputError(RayInterfaceError):
    """Raised when an exact-zero ray input is forbidden by contract."""


RAY_CLASSES = (
    "exact_zero",
    "constant_nonzero",
    "structured_small",
    "structured_normal",
    "nonfinite",
)

# Diagnostic threshold only. It MUST NOT affect field acceptance.
RAY_SMALL_AMPLITUDE = 1e-12


def _validate_pair_shape(Rx, Ry):
    if Rx.shape != Ry.shape:
        raise RayInterfaceError("Rx_2d and Ry_2d must share shape")
    if Rx.ndim != 2:
        raise RayInterfaceError(
            f"ray inputs must be 2D image-plane fields, got ndim={Rx.ndim}"
        )


def _hash_array(arr: np.ndarray) -> str:
    """Hash dtype + shape + raw contiguous bytes."""
    arr = np.ascontiguousarray(arr)
    payload = arr.dtype.str.encode("utf-8")
    payload += str(arr.shape).encode("utf-8")
    payload += arr.tobytes()
    return hashlib.sha256(payload).hexdigest()


def classify_ray_input(Rx_2d, Ry_2d):
    """Classify a 2D image-plane ray input.

    Classification is scale-aware but acceptance-neutral: a tiny field
    with real spatial structure remains a valid structured field.
    """
    Rx = np.asarray(Rx_2d, dtype=np.float64)
    Ry = np.asarray(Ry_2d, dtype=np.float64)
    _validate_pair_shape(Rx, Ry)

    finite = np.all(np.isfinite(Rx)) and np.all(np.isfinite(Ry))
    if not finite:
        return {
            "classification": "nonfinite",
            "R_rms": float("nan"),
            "R_max": float("nan"),
            "var_Rx": float("nan"),
            "var_Ry": float("nan"),
        }

    max_abs = float(max(np.max(np.abs(Rx)), np.max(np.abs(Ry))))
    var_Rx = float(np.var(Rx))
    var_Ry = float(np.var(Ry))
    R_rms = float(np.sqrt(0.5 * (np.mean(Rx ** 2) + np.mean(Ry ** 2))))

    if max_abs == 0.0:
        cls = "exact_zero"
    else:
        # Floating-point floor for a mathematically constant field.
        # This is relative to field scale and does not impose a physical
        # amplitude cutoff.
        eps = np.finfo(np.float64).eps
        var_tol = max_abs ** 2 * eps + np.finfo(np.float64).tiny
        if var_Rx <= var_tol and var_Ry <= var_tol:
            cls = "constant_nonzero"
        else:
            cls = (
                "structured_small"
                if R_rms < RAY_SMALL_AMPLITUDE
                else "structured_normal"
            )

    return {
        "classification": cls,
        "R_rms": R_rms,
        "R_max": max_abs,
        "var_Rx": var_Rx,
        "var_Ry": var_Ry,
    }


def prepare_ray_input(
    Rx_2d,
    Ry_2d,
    metadata,
    require_nontrivial=True,
    upstream_rms=None,
):
    """Validate and package a 2D field for ray propagation.

    Nonfinite inputs are invalid regardless of ``require_nontrivial``.
    ``require_nontrivial=False`` only permits the deliberate zero-field
    control lane.
    """
    if not isinstance(metadata, dict):
        raise RayInterfaceError("metadata must be a dict")
    for key in ("candidate_id", "cluster_id", "transform_id", "role"):
        if key not in metadata:
            raise RayInterfaceError(f"metadata missing required key {key!r}")

    Rx = np.asarray(Rx_2d, dtype=np.float64)
    Ry = np.asarray(Ry_2d, dtype=np.float64)
    _validate_pair_shape(Rx, Ry)

    cls = classify_ray_input(Rx, Ry)

    # Nonfinite is always invalid. This is an integrity failure, not a
    # triviality choice.
    if cls["classification"] == "nonfinite":
        raise RayInterfaceError(
            f"ray input contains NaN/Inf: candidate={metadata['candidate_id']}, "
            f"cluster={metadata['cluster_id']}"
        )

    if require_nontrivial and cls["classification"] == "exact_zero":
        raise TrivialRayInputError(
            f"ray input is exact zero: candidate={metadata['candidate_id']}, "
            f"cluster={metadata['cluster_id']}"
        )

    if upstream_rms is not None:
        upstream_rms = float(upstream_rms)
        if not np.isfinite(upstream_rms) or upstream_rms < 0.0:
            raise RayInterfaceError("upstream_rms must be finite and >= 0")

    scale_info = {
        "R_rms": cls["R_rms"],
        "R_max": cls["R_max"],
        "var_Rx": cls["var_Rx"],
        "var_Ry": cls["var_Ry"],
        "R_rms_over_A_rms": (
            cls["R_rms"] / upstream_rms
            if upstream_rms is not None and upstream_rms > 0.0
            else None
        ),
    }

    sha_rx = _hash_array(Rx)
    sha_ry = _hash_array(Ry)
    combined = hashlib.sha256((sha_rx + sha_ry).encode("utf-8")).hexdigest()

    stats = field_statistics_vector(Rx, Ry, np.zeros_like(Rx))
    artifact = FieldArtifact(
        data={"Rx": Rx, "Ry": Ry},
        artifact_id=(
            f"ray_input_{metadata['candidate_id']}_"
            f"{metadata['cluster_id']}_{metadata['transform_id']}"
        ),
        module_name="pbuf.core.ray_interface",
        module_version="1.2.0-second-review",
        source_artifact_ids=metadata.get("source_artifact_ids", []),
        candidate_id=metadata["candidate_id"],
        cluster_id=metadata["cluster_id"],
        transform_id=metadata["transform_id"],
        sha256=combined,
        statistics={
            **stats,
            "ray_classification": cls["classification"],
            "ray_scale": scale_info,
        },
        role=metadata["role"],
    )
    return artifact


def ray_input_fingerprint(artifact):
    return {
        "sha256": artifact.sha256,
        "statistics": artifact.statistics,
        "role": artifact.role,
        "candidate_id": artifact.candidate_id,
        "cluster_id": artifact.cluster_id,
        "transform_id": artifact.transform_id,
    }


# ----------------------------------------------------------------------
# Self-check / validation fixtures
# ----------------------------------------------------------------------
def _metadata(candidate="PL1_PM1_PS2"):
    return {
        "candidate_id": candidate,
        "cluster_id": "MACS0416",
        "transform_id": "RC0",
        "role": "central",
    }


def _classification_test():
    rows = []

    cases = [
        (
            "exact_zero",
            np.zeros((4, 5)),
            np.zeros((4, 5)),
            "exact_zero",
        ),
        (
            "constant_nonzero",
            np.full((4, 5), 0.5),
            np.full((4, 5), 0.3),
            "constant_nonzero",
        ),
    ]

    x = np.linspace(0.0, 1.0, 50).reshape(1, -1)
    cases.append((
        "structured_small",
        1e-15 * np.sin(2 * np.pi * x),
        1e-15 * np.cos(2 * np.pi * x),
        "structured_small",
    ))

    rng = np.random.RandomState(0)
    cases.append((
        "structured_normal",
        rng.randn(8, 8),
        rng.randn(8, 8),
        "structured_normal",
    ))

    for name, Rx, Ry, expected in cases:
        got = classify_ray_input(Rx, Ry)["classification"]
        rows.append({"input": name, "class": got, "passes": got == expected})

    got = classify_ray_input(
        np.array([[float("nan")]]), np.zeros((1, 1))
    )["classification"]
    rows.append({"input": "nonfinite", "class": got,
                 "passes": got == "nonfinite"})

    return {"rows": rows, "passes": all(r["passes"] for r in rows)}


def _acceptance_policy_test():
    # exact zero: rejected normally, permitted only for deliberate control.
    zero = np.zeros((4, 5))
    try:
        prepare_ray_input(zero, zero, _metadata(), require_nontrivial=True)
    except TrivialRayInputError:
        zero_rejected = True
    else:
        zero_rejected = False

    zero_artifact = prepare_ray_input(
        zero, zero, _metadata("ZERO_CONTROL"), require_nontrivial=False
    )
    zero_control_ok = (
        zero_artifact.statistics["ray_classification"] == "exact_zero"
    )

    # Nonfinite must be rejected even if zero/trivial controls are allowed.
    bad = np.zeros((2, 2)); bad[0, 0] = np.nan
    try:
        prepare_ray_input(
            bad, np.zeros((2, 2)), _metadata("NONFINITE"),
            require_nontrivial=False,
        )
    except RayInterfaceError:
        nonfinite_rejected = True
    else:
        nonfinite_rejected = False

    return {
        "passes": zero_rejected and zero_control_ok and nonfinite_rejected,
        "zero_rejected": zero_rejected,
        "zero_control_ok": zero_control_ok,
        "nonfinite_rejected_even_when_trivial_allowed": nonfinite_rejected,
    }


def _structured_small_acceptance_test():
    # Deliberately below RAY_SMALL_AMPLITUDE but spatially structured.
    x = np.linspace(0.0, 2 * np.pi, 64, endpoint=False)
    Rx = 1e-15 * np.sin(x).reshape(1, -1)
    Ry = 1e-15 * np.cos(x).reshape(1, -1)
    art = prepare_ray_input(Rx, Ry, _metadata("SMALL"), True)
    return {
        "passes": art.statistics["ray_classification"] == "structured_small",
        "classification": art.statistics["ray_classification"],
    }


def _wc6_old_threshold_fixture_test():
    """The historical 1e-10 fixture is not 'structured_small'.

    Since RAY_SMALL_AMPLITUDE=1e-12, 1e-10 is correctly classified as
    structured_normal. The important scientific requirement is that it
    is accepted, not that it be mislabeled small.
    """
    x = np.linspace(0.0, 2 * np.pi, 64, endpoint=False)
    Rx = 1e-10 * np.sin(x).reshape(1, -1)
    Ry = 1e-10 * np.cos(x).reshape(1, -1)
    art = prepare_ray_input(Rx, Ry, _metadata("WC6"), True)
    return {
        "passes": art.statistics["ray_classification"] == "structured_normal",
        "classification": art.statistics["ray_classification"],
    }


def _hash_lineage_test():
    metadata = _metadata()
    rng = np.random.RandomState(1)
    Rx = rng.randn(8, 8); Ry = rng.randn(8, 8)
    a1 = prepare_ray_input(Rx, Ry, metadata)
    a2 = prepare_ray_input(Rx.copy(), Ry.copy(), metadata)
    Ry2 = Ry.copy(); Ry2[0, 0] += 1e-9
    a3 = prepare_ray_input(Rx, Ry2, metadata)
    return {
        "passes": a1.sha256 == a2.sha256 and a1.sha256 != a3.sha256,
    }


def _dimension_contract_test():
    try:
        classify_ray_input(np.zeros((2, 2, 2)), np.zeros((2, 2, 2)))
    except RayInterfaceError:
        return {"passes": True}
    return {"passes": False}


if __name__ == "__main__":
    r = _classification_test(); assert r["passes"], r
    for row in r["rows"]:
        print(f"M15 classification {row['input']}: {row['class']}")

    r = _acceptance_policy_test(); assert r["passes"], r
    print("M15 acceptance policy: zero-control and nonfinite rules PASS")

    r = _structured_small_acceptance_test(); assert r["passes"], r
    print(f"M15 structured-small accepted: {r['classification']}")

    r = _wc6_old_threshold_fixture_test(); assert r["passes"], r
    print(f"M15 historical 1e-10 fixture: {r['classification']}")

    r = _hash_lineage_test(); assert r["passes"], r
    print("M15 hash lineage: deterministic and change-sensitive")

    r = _dimension_contract_test(); assert r["passes"], r
    print("M15 2D image-plane contract: enforced")

    print("M15 ray interface: all checks passed")
