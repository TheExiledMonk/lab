"""Stage-by-stage parity utilities for the legacy and canonical WL paths."""

import hashlib
import numpy as np

from pbuf.labs.foundation import native_full_state_100pct_observer_coverage001 as LEGACY
from pbuf.labs.foundation import native_full_state_100pct_observer_coverage_fix001 as FIX
from pbuf.labs.foundation import current_native_five_cluster_observable_benchmark001 as CUR
from pbuf.labs.foundation import m10_coverage_25pct_science001 as BASE
from pbuf.core import los_projection as M14
from .pipeline import compare_with_observations, run_wl_pipeline

_LEGACY_CACHE: dict[tuple[str, str], dict] = {}
_CANONICAL_CACHE: dict[tuple[str, str], object] = {}


def sha256_float64_array(a: np.ndarray) -> str:
    arr = np.ascontiguousarray(np.asarray(a, dtype=np.float64))
    return hashlib.sha256(arr.tobytes()).hexdigest()


def relative_rms_error(old, new) -> float:
    a, b = np.asarray(old, dtype=np.float64), np.asarray(new, dtype=np.float64)
    finite = np.isfinite(a) & np.isfinite(b)
    if not np.any(finite):
        return 0.0
    d = a[finite] - b[finite]
    denom = max(float(np.sqrt(np.mean(a[finite] ** 2))), 1e-30)
    return float(np.sqrt(np.mean(d ** 2)) / denom)


def array_report(old, new, *, rtol: float, atol: float, exact: bool = False) -> dict:
    a, b = np.asarray(old), np.asarray(new)
    same = np.array_equal(a, b, equal_nan=True) if exact else np.allclose(a, b, rtol=rtol, atol=atol, equal_nan=True)
    finite = np.isfinite(a) & np.isfinite(b)
    return {
        "pass": bool(same),
        "max_abs_error": float(np.max(np.abs(a[finite] - b[finite]))) if np.any(finite) else 0.0,
        "relative_rms_error": relative_rms_error(a, b),
    }


def run_legacy_lane(cluster: dict, coverage: str) -> dict:
    key = (cluster["id"], coverage)
    if key in _LEGACY_CACHE:
        return _LEGACY_CACHE[key]
    data = CUR.local_cluster(cluster)
    vector, channel = CUR.current_native_m10(data["rho3"])
    launch = BASE._launch_expanded_25pct() if coverage == "25pct" else LEGACY._launch_full_100pct()
    expected = LEGACY.EXPECTED_SUPPORT25 if coverage == "25pct" else LEGACY.EXPECTED_SUPPORT100
    chain = LEGACY._run_g3d_with_launch(vector, launch[0], launch[1], expected)
    lane = FIX._decode_lane_launch_aware(data, channel, chain, launch[0], launch[1], f"coverage_{coverage}")
    result = {"data": data, "vector": vector, "launch": launch, "chain": chain, **lane}
    # Deep parity and metric parity reuse Abell2744. Avoid retaining the very
    # large propagation snapshots for every cluster in the five-cluster run.
    if cluster["id"] == "Abell2744":
        _LEGACY_CACHE[key] = result
    return result


def _run_canonical(cluster: dict, coverage: str):
    key = (cluster["id"], coverage)
    if key in _CANONICAL_CACHE:
        return _CANONICAL_CACHE[key]
    result = run_wl_pipeline(cluster, coverage)
    if cluster["id"] == "Abell2744":
        _CANONICAL_CACHE[key] = result
    return result


def deep_parity(cluster: dict, coverage: str) -> dict:
    old = run_legacy_lane(cluster, coverage)
    new = _run_canonical(cluster, coverage)
    exact = {}
    for name, a, b in zip(("m10_x", "m10_y", "m10_z"), old["vector"], new.native_response["m10_vector"]):
        exact[name] = array_report(a, b, rtol=0.0, atol=0.0, exact=True)
    old_projected = M14.project_vector_to_image_plane(*old["vector"], los_axis="z")
    old_los = old["chain"]
    for name, a, b in (("los_Rx", old_projected["comp_1"], new.los["Rx"]),
                       ("los_Ry", old_projected["comp_2"], new.los["Ry"]),
                       ("los_mag", old_los["los_mag"], new.los["los_mag"])):
        exact[name] = array_report(a, b, rtol=0.0, atol=0.0, exact=True)
    for name, a, b in zip(("launch_x0", "launch_y0", "launch_vx0", "launch_vy0"), old["launch"],
                          (new.launch.x0, new.launch.y0, new.launch.vx0, new.launch.vy0)):
        exact[name] = array_report(a, b, rtol=0.0, atol=0.0, exact=True)
    final = {k: array_report(old["chain"]["checkpoints"][LEGACY.G3D.CHECKPOINT][k], new.propagation["final_snapshot"][k], rtol=1e-13, atol=1e-13)
             for k in ("x", "y", "z", "vx", "vy", "vz")}
    screen = {k: array_report(old["screen"][k], new.screen[k], rtol=1e-13, atol=1e-13)
              for k in old["screen"] if np.issubdtype(np.asarray(old["screen"][k]).dtype, np.number)}
    names = list(old["bank"])
    channels = {k: array_report(old["bank"][k], new.channel_bank[k], rtol=1e-12, atol=1e-12) for k in names if k in new.channel_bank}
    candidate_names = list(old["candidates"])
    candidates = {k: array_report(old["candidates"][k], new.reconstruction_candidates[k], rtol=1e-12, atol=1e-12)
                  for k in candidate_names if k in new.reconstruction_candidates}
    hashes = {}
    hash_arrays = {"m10_x": new.native_response["m10_vector"][0], "m10_y": new.native_response["m10_vector"][1],
                   "m10_z": new.native_response["m10_vector"][2], "los_Rx": new.los["Rx"], "los_Ry": new.los["Ry"],
                   "launch_x0": new.launch.x0, "launch_y0": new.launch.y0,
                   **{f"final_{k}": new.propagation["final_snapshot"][k] for k in ("x", "y", "z", "vx", "vy", "vz")}}
    hashes.update({k: sha256_float64_array(v) for k, v in hash_arrays.items()})
    return {"coverage": coverage, "exact": exact, "final": final, "screen": screen, "channels": channels,
            "channel_names_exact": names == list(new.channel_bank), "candidates": candidates,
            "candidate_names_exact": candidate_names == list(new.reconstruction_candidates), "hashes": hashes}


def metric_parity(cluster: dict, coverage: str) -> dict:
    old = run_legacy_lane(cluster, coverage)
    new = _run_canonical(cluster, coverage)
    targets = LEGACY.DEC._targets_after_decoding(old["data"])
    old_metrics = LEGACY.DEC._compare_candidates(old["candidates"], targets)
    new_metrics = compare_with_observations(new, new.source["data"])
    errors = []
    for candidate in old_metrics:
        for target in old_metrics[candidate]:
            for metric in ("pearson", "spearman"):
                a, b = old_metrics[candidate][target][metric], new_metrics[candidate][target][metric]
                if not (np.isnan(a) and np.isnan(b)):
                    errors.append(abs(float(a) - float(b)))
    return {"pass": bool(list(old_metrics) == list(new_metrics) and max(errors, default=0.0) <= 1e-12),
            "max_metric_abs_error": max(errors, default=0.0)}
