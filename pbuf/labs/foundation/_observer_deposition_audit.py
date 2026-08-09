"""Shared target-blind machinery for Dev Doc 110 observer audit labs."""

from __future__ import annotations

import time
import numpy as np

from pbuf.core import observable_extraction as M16
from pbuf.labs.foundation import native_full_received_state_information_retention001 as RET
from pbuf.labs.foundation import native_full_state_2d_reconstruction_decoder_sweep001 as DEC
from pbuf.wl.channels import decode_full_channel_bank
from pbuf.wl.config import EXTENT, OBS_BINS
from pbuf.wl.deposition import METHODS
from pbuf.wl.received_state import build_received_state
from pbuf.wl.reconstruction import build_reconstruction_candidates
from pbuf.wl.screen import build_detector_screen
from pbuf.wl.observer_cache import ObserverPrimitiveCache, ObserverStateId

EPSILONS = (1e-16, 1e-15, 1e-14, 1e-13, 1e-12)
DIRECTIONS = {"+u": (1, 0), "-u": (-1, 0), "+v": (0, 1), "-v": (0, -1),
              "+(u,v)": (1, 1), "+(u,-v)": (1, -1)}


def error(a, b):
    x, y = np.asarray(a, dtype=np.float64), np.asarray(b, dtype=np.float64)
    finite = np.isfinite(x) & np.isfinite(y)
    if not np.any(finite):
        return {"max_abs_error": 0.0, "relative_rms_error": 0.0, "finite": True}
    delta = x[finite] - y[finite]
    scale = max(float(np.sqrt(np.mean(x[finite] ** 2))), 1e-30)
    return {"max_abs_error": float(np.max(np.abs(delta))),
            "relative_rms_error": float(np.sqrt(np.mean(delta ** 2)) / scale),
            "finite": bool(np.all(np.isfinite(x)) and np.all(np.isfinite(y)))}


def bank_error(reference, candidate):
    rows = {name: error(reference[name], candidate[name]) for name in reference}
    rel = np.array([row["relative_rms_error"] for row in rows.values()])
    return {"max_abs_error": max(row["max_abs_error"] for row in rows.values()),
            "mean_relative_rms_error": float(np.mean(rel)),
            "median_relative_rms_error": float(np.median(rel)),
            "max_relative_rms_error": float(np.max(rel)),
            "number_channels_changed": sum(row["max_abs_error"] > 0 for row in rows.values()),
            "number_channels_failed_1e-9": sum(row["relative_rms_error"] > 1e-9 for row in rows.values()),
            "all_finite": all(row["finite"] for row in rows.values())}


def decode(prepared, propagation, screen, method, *, reconstruct=True, cache=None,
           state_id=None, kde_backend=None):
    received = build_received_state(prepared["launch"], propagation, screen)
    t0 = time.perf_counter()
    decoded = decode_full_channel_bank(screen, received, method, cache=cache,
                                       state_id=state_id, kde_backend=kde_backend)
    decode_seconds = time.perf_counter() - t0
    candidates, meta, reconstruction_seconds = {}, {}, 0.0
    if reconstruct:
        t0 = time.perf_counter()
        candidates, meta = build_reconstruction_candidates(decoded["bank"], decoded["family"])
        reconstruction_seconds = time.perf_counter() - t0
    return {**decoded, "received": received, "candidates": candidates, "meta": meta,
            "decode_seconds": decode_seconds, "reconstruction_seconds": reconstruction_seconds}


def perturb(screen, epsilon, direction):
    du, dv = DIRECTIONS[direction]
    out = dict(screen)
    out["uf"] = np.asarray(screen["uf"], dtype=np.float64) + du * epsilon
    out["vf"] = np.asarray(screen["vf"], dtype=np.float64) + dv * epsilon
    return out


def boundary_proximity(screen):
    edges = np.linspace(-EXTENT, EXTENT, OBS_BINS + 1, dtype=np.float64)
    u, v = np.asarray(screen["uf"]), np.asarray(screen["vf"])
    width = edges[1] - edges[0]
    du = np.minimum(np.mod(u + EXTENT, width), width - np.mod(u + EXTENT, width))
    dv = np.minimum(np.mod(v + EXTENT, width), width - np.mod(v + EXTENT, width))
    mask = (du <= 1e-12) | (dv <= 1e-12)
    return mask, {"number_of_boundary_near_rays": int(mask.sum()),
                  "fraction_of_total_rays": float(mask.mean())}


def conservation(method, screen):
    u, v = screen["uf"], screen["vf"]
    valid = (np.isfinite(u) & np.isfinite(v) & (u >= -EXTENT) & (u <= EXTENT) &
             (v >= -EXTENT) & (v <= EXTENT))
    deposited = method.deposit(u, v, None, bins=OBS_BINS, extent=EXTENT)
    expected = int(valid.sum())
    absolute_error = abs(float(deposited.sum()) - expected)
    return {"deposited_weight": float(deposited.sum()), "valid_rays": expected,
            "absolute_error": absolute_error,
            "pass": bool(np.isclose(deposited.sum(), expected, rtol=1e-12, atol=1e-12))}


def information(bank):
    X, _ = RET._standardized_matrix(bank, list(bank))
    singular = np.linalg.svd(X, compute_uv=False) if X.size else np.array([])
    variance = singular ** 2
    total = float(variance.sum())
    p = variance / total if total else variance
    return {"numerical_rank": int(np.linalg.matrix_rank(X)),
            "effective_rank": float(np.exp(-np.sum(p[p > 0] * np.log(p[p > 0])))) if p.size else 0.0,
            "participation_ratio": float(total ** 2 / np.sum(variance ** 2)) if np.any(variance) else 0.0,
            "top5_variance_fraction": float(p[:5].sum()), "top10_variance_fraction": float(p[:10].sum())}


def morphology(control, bank):
    rows = []
    for name in control:
        x, y = np.asarray(control[name]).ravel(), np.asarray(bank[name]).ravel()
        finite = np.isfinite(x) & np.isfinite(y)
        if finite.sum() < 2 or np.std(x[finite]) == 0 or np.std(y[finite]) == 0:
            pearson = spearman = 1.0 if np.array_equal(x, y, equal_nan=True) else 0.0
        else:
            pearson = float(M16.safe_pearson(x[finite], y[finite]))
            spearman = float(M16.safe_spearman(x[finite], y[finite]))
        rows.append((pearson, spearman, error(x, y)["relative_rms_error"]))
    return {"median_channel_morphology_pearson": float(np.median([x[0] for x in rows])),
            "median_channel_morphology_spearman": float(np.median([x[1] for x in rows])),
            "median_channel_relative_rms": float(np.median([x[2] for x in rows]))}


def ray_difference(cpu_screen, gpu_screen, cpu_prop, gpu_prop):
    def rr(names, a, b):
        x = np.concatenate([np.asarray(a[n], dtype=np.float64) for n in names])
        y = np.concatenate([np.asarray(b[n], dtype=np.float64) for n in names])
        return error(x, y)
    pos = rr(("x", "y", "z"), cpu_prop["final_snapshot"], gpu_prop["final_snapshot"])
    vel = rr(("vx", "vy", "vz"), cpu_prop["final_snapshot"], gpu_prop["final_snapshot"])
    return {"max_abs_position_difference": pos["max_abs_error"],
            "max_abs_velocity_difference": vel["max_abs_error"],
            "relative_rms_position_difference": pos["relative_rms_error"],
            "relative_rms_velocity_difference": vel["relative_rms_error"]}


def audit_lane(prepared, cpu_prop, gpu_prop):
    cpu_screen = build_detector_screen(prepared["launch"], cpu_prop)
    gpu_screen = build_detector_screen(prepared["launch"], gpu_prop)
    target_blind = {}; cache = ObserverPrimitiveCache()
    cpu_id = ObserverStateId(f"cpu_{prepared['launch'].coverage_label}_base", backend="cpu")
    gpu_id = ObserverStateId(f"vulkan_{prepared['launch'].coverage_label}_base", backend="vulkan")
    for method in METHODS:
        base = decode(prepared, cpu_prop, cpu_screen, method, cache=cache, state_id=cpu_id)
        gpu = decode(prepared, gpu_prop, gpu_screen, method, cache=cache, state_id=gpu_id)
        curves = {}
        machine_ok = True
        for eps in EPSILONS:
            directions = {}
            for direction in DIRECTIONS:
                changed = decode(prepared, cpu_prop, perturb(cpu_screen, eps, direction), method,
                                 reconstruct=False, cache=cache,
                                 state_id=ObserverStateId(cpu_id.base_state, f"translate_{direction}_{eps}", "cpu"))
                directions[direction] = bank_error(base["bank"], changed["bank"])
                if eps <= 1e-14:
                    row = directions[direction]
                    machine_ok &= row["all_finite"] and row["median_relative_rms_error"] <= 1e-9 and row["max_relative_rms_error"] <= 1e-6
            curves[str(eps)] = directions
        target_blind[method.name] = {"base": base, "gpu": gpu,
            "conservation": conservation(method, cpu_screen), "perturbation": curves,
            "machine_scale_stable": bool(machine_ok), "cpu_vulkan": bank_error(base["bank"], gpu["bank"]),
            "information": information(base["bank"])}
    control = target_blind["hard_bin_current"]["base"]["bank"]
    for row in target_blind.values():
        row["morphology"] = morphology(control, row["base"]["bank"])
    # Binding boundary: targets are first accessed only after all target-blind builds.
    targets = DEC._targets_after_decoding(prepared["source"]["data"])
    for row in target_blind.values():
        row["observational"] = DEC._compare_candidates(row["base"]["candidates"], targets)
    mask, proximity = boundary_proximity(cpu_screen)
    proximity["deposition_probe"] = {m.name: {
        str(eps): {direction: error(
            m.deposit(cpu_screen["uf"][mask], cpu_screen["vf"][mask], None, bins=OBS_BINS, extent=EXTENT),
            m.deposit(cpu_screen["uf"][mask] + DIRECTIONS[direction][0] * eps,
                      cpu_screen["vf"][mask] + DIRECTIONS[direction][1] * eps,
                      None, bins=OBS_BINS, extent=EXTENT))
                   for direction in DIRECTIONS} for eps in EPSILONS}
        for m in METHODS}
    return target_blind, ray_difference(cpu_screen, gpu_screen, cpu_prop, gpu_prop), proximity


def survivors(rows):
    hard_erank = rows["hard_bin_current"]["information"]["effective_rank"]
    order = ("hard_bin_half_open", "nearest_center", "bilinear_cic", "tsc_3x3", "gaussian_sigma_half_cell")
    return [name for name in order if rows[name]["machine_scale_stable"] and
            rows[name]["conservation"]["pass"] and rows[name]["cpu_vulkan"]["number_channels_failed_1e-9"] == 0 and
            len(rows[name]["base"]["bank"]) == 45 and len(rows[name]["base"]["candidates"]) == 68 and
            rows[name]["information"]["effective_rank"] >= .8 * hard_erank and
            rows[name]["morphology"]["median_channel_morphology_pearson"] >= .95]
