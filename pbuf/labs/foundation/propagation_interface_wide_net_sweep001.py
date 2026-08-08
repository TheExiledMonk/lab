#!/usr/bin/env python3
"""PBUF FOUNDATION — PROPAGATION INTERFACE WIDE-NET SWEEP 001.

Purpose
-------
Replace serial one-hypothesis-at-a-time interface audits with one broad,
discriminatory sweep over many already-available local observables.

The sweep uses the same five canonical local benchmark sources and compares
candidate native/local interface observables against the historical unit-loading
M10 interface ONLY as an internal structural reference. No observed lensing
values are used. No candidate is fitted, rescaled, normalised to M10, or fed to
G3D.

Candidates are ranked on several independent diagnostics:
  * amplitude-ratio stability across sources;
  * morphology/component correlation with historical M10;
  * vector-direction cosine/sign agreement;
  * locality / construction class;
  * source-scaling stability.

The ranking is diagnostic. It does not promote a candidate to a physical law.
"""
from __future__ import annotations

import json
import math
import statistics
import subprocess
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))

from pbuf.core import benchmark_data as BENCH
from pbuf.models import a8_state as A8
import pbuf.labs.foundation.m10_coverage_25pct_science001 as BASE
import pbuf.labs.foundation.m10_local_interface_decomposition_audit001 as DEC
import pbuf.labs.foundation.native_accumulated_full_lensing_local_benchmark001 as LOCAL
import pbuf.labs.foundation.native_accumulated_full_lensing001 as FULL

LAB_ID = "PBUF-FOUNDATION-PROPAGATION-INTERFACE-WIDE-NET-SWEEP-001"
EXPECTED_CLUSTER_IDS = ("Abell2744", "MACS0416", "MACS1149", "AbellS1063", "Abell370")
EPS = 1.0e-30


def git(*args: str) -> str:
    return subprocess.check_output(["git", *args], cwd=ROOT, text=True).strip()


def repo_state() -> dict:
    return {
        "repository": "TheExiledMonk/lab",
        "branch": git("rev-parse", "--abbrev-ref", "HEAD"),
        "head_sha": git("rev-parse", "HEAD"),
        "tracked_changes": git("diff", "--name-only"),
        "staged_changes": git("diff", "--name-only", "--cached"),
    }


def rms(a) -> float:
    x = np.asarray(a, dtype=np.float64)
    return float(np.sqrt(np.mean(x * x))) if x.size else float("nan")


def vrms(v) -> float:
    x, y, z = (np.asarray(q, dtype=np.float64) for q in v)
    return float(np.sqrt(np.mean(x*x + y*y + z*z)))


def cv(values) -> float:
    vals = [float(v) for v in values if math.isfinite(float(v))]
    if len(vals) < 2:
        return float("nan")
    m = statistics.fmean(vals)
    return float(statistics.pstdev(vals) / max(abs(m), EPS))


def safe_corr(a, b) -> float:
    x = np.asarray(a, dtype=np.float64).ravel()
    y = np.asarray(b, dtype=np.float64).ravel()
    mask = np.isfinite(x) & np.isfinite(y)
    if int(mask.sum()) < 2:
        return float("nan")
    x = x[mask]; y = y[mask]
    sx = float(np.std(x)); sy = float(np.std(y))
    if sx <= EPS or sy <= EPS:
        return 0.0
    return float(np.corrcoef(x, y)[0, 1])


def vector_component_corr(a, b) -> float:
    vals = [safe_corr(x, y) for x, y in zip(a, b)]
    vals = [v for v in vals if math.isfinite(v)]
    return float(statistics.fmean(vals)) if vals else float("nan")


def vector_cosine(a, b) -> tuple[float, float]:
    ax, ay, az = (np.asarray(q, dtype=np.float64) for q in a)
    bx, by, bz = (np.asarray(q, dtype=np.float64) for q in b)
    dot = ax*bx + ay*by + az*bz
    am = np.sqrt(ax*ax + ay*ay + az*az)
    bm = np.sqrt(bx*bx + by*by + bz*bz)
    mask = (am > EPS) & (bm > EPS) & np.isfinite(dot)
    if not np.any(mask):
        return float("nan"), float("nan")
    cos = dot[mask] / (am[mask] * bm[mask])
    return float(np.mean(cos)), float(np.mean(cos > 0.0))


def grad3(s):
    gz, gy, gx = np.gradient(np.asarray(s, dtype=np.float64), edge_order=1)
    return (gx, gy, gz)


def laplacian(s):
    f = np.asarray(s, dtype=np.float64)
    p = np.pad(f, ((1,1),(1,1),(1,1)), mode="edge")
    return (
        p[1:-1,1:-1,2:] + p[1:-1,1:-1,:-2]
        + p[1:-1,2:,1:-1] + p[1:-1,:-2,1:-1]
        + p[2:,1:-1,1:-1] + p[:-2,1:-1,1:-1] - 6.0*f
    )


def face_bond_vector(s):
    """Centered local face-difference vector from a scalar state."""
    f = np.asarray(s, dtype=np.float64)
    p = np.pad(f, ((1,1),(1,1),(1,1)), mode="edge")
    dx = 0.5 * (p[1:-1,1:-1,2:] - p[1:-1,1:-1,:-2])
    dy = 0.5 * (p[1:-1,2:,1:-1] - p[1:-1,:-2,1:-1])
    dz = 0.5 * (p[2:,1:-1,1:-1] - p[:-2,1:-1,1:-1])
    return dx, dy, dz


def traction_face_vector(u):
    """Centered bounded-strain face-traction vector, no fitted coefficient."""
    f = np.asarray(u, dtype=np.float64)
    p = np.pad(f, ((1,1),(1,1),(1,1)), mode="edge")
    xp = p[1:-1,1:-1,2:] - f; xm = f - p[1:-1,1:-1,:-2]
    yp = p[1:-1,2:,1:-1] - f; ym = f - p[1:-1,:-2,1:-1]
    zp = p[2:,1:-1,1:-1] - f; zm = f - p[:-2,1:-1,1:-1]

    def sigma(q):
        frac2 = (q / FULL.EPSILON_MAX) ** 2
        if np.any(frac2 >= 1.0):
            raise RuntimeError("bond exceeded bounded-strain domain")
        return FULL.K0 * q / (1.0 - frac2)

    return (
        0.5 * (sigma(xp) + sigma(xm)),
        0.5 * (sigma(yp) + sigma(ym)),
        0.5 * (sigma(zp) + sigma(zm)),
    )


def net_traction_imbalance_vector(u):
    """Difference of positive/negative face tractions on each axis."""
    f = np.asarray(u, dtype=np.float64)
    p = np.pad(f, ((1,1),(1,1),(1,1)), mode="edge")
    xp = p[1:-1,1:-1,2:] - f; xm = f - p[1:-1,1:-1,:-2]
    yp = p[1:-1,2:,1:-1] - f; ym = f - p[1:-1,:-2,1:-1]
    zp = p[2:,1:-1,1:-1] - f; zm = f - p[:-2,1:-1,1:-1]

    def sigma(q):
        frac2 = (q / FULL.EPSILON_MAX) ** 2
        if np.any(frac2 >= 1.0):
            raise RuntimeError("bond exceeded bounded-strain domain")
        return FULL.K0 * q / (1.0 - frac2)

    return sigma(xp) - sigma(xm), sigma(yp) - sigma(ym), sigma(zp) - sigma(zm)


def historical_fields(rho3):
    rng = np.random.RandomState(BASE.SEED)
    eq = np.asarray(rho3, dtype=np.float64)
    initial = {
        "rho_3d": rho3.copy(),
        "u_slow0": eq.copy(),
        "u_fast0": eq + A8.A8_INIT_INJECTION_NOISE * rng.randn(*rho3.shape),
    }
    state = BASE._evolve(initial)
    candidate = BASE._candidate(state)
    m10 = tuple(np.asarray(v, dtype=np.float64) for v in BASE._interface_vector(candidate))
    uf = np.asarray(state["u_fast"], dtype=np.float64)
    us = np.asarray(state["u_slow"], dtype=np.float64)
    c = np.asarray(state["c_state"], dtype=np.float64)
    d = 0.5 * (uf - us)
    return {"uf": uf, "us": us, "c": c, "d": d, "m10": m10}


def candidate_catalog(rho3, hist, native):
    c_native = np.asarray(native["c_state"], dtype=np.float64)
    u = np.asarray(native["accumulated"], dtype=np.float64)
    coef_fast = A8.A8_INIT_DT * A8.A8_INIT_OMEGA * A8.A8_INIT_K
    coef_slow = A8.A8_INIT_DT * A8.A8_INIT_SLOW_TIMESCALE

    gu = grad3(u)
    gc = grad3(c_native)
    lap_u = laplacian(u)
    lap_c = laplacian(c_native)
    traction = traction_face_vector(u)
    imbalance = net_traction_imbalance_vector(u)
    bond_u = face_bond_vector(u)
    bond_c = face_bond_vector(c_native)

    h_c = face_bond_vector(hist["c"])
    h_d = face_bond_vector(hist["d"])
    h_uf = face_bond_vector(hist["uf"])
    h_us = face_bond_vector(hist["us"])

    # One actual frozen A8 update from equal native channels. This is a state
    # diagnostic, not a candidate law: it asks what the update operator itself
    # produces before any G3D use.
    n = A8.neighbours6_face_zero_flux_3d(u)
    dfast = A8.A8_INIT_DT * A8.A8_INIT_OMEGA * A8.A8_INIT_K * (n - u)
    dslow = A8.A8_INIT_DT * A8.A8_INIT_SLOW_TIMESCALE * (n - u)

    return {
        "native_grad_u": (gu, "local_gradient"),
        "native_face_bond_u": (bond_u, "local_bond"),
        "native_face_traction": (traction, "local_constitutive_bond"),
        "native_net_traction_imbalance": (imbalance, "local_force_balance"),
        "native_grad_c_state": (gc, "local_source_state_gradient"),
        "native_grad_laplacian_u": (grad3(lap_u), "local_curvature_gradient"),
        "native_grad_laplacian_c": (grad3(lap_c), "local_source_curvature_gradient"),
        "native_fast_scaled_traction": (tuple(coef_fast*x for x in traction), "frozen_fast_transfer_probe"),
        "native_slow_scaled_traction": (tuple(coef_slow*x for x in traction), "frozen_slow_transfer_probe"),
        "native_sum_scaled_traction": (tuple((coef_fast+coef_slow)*x for x in traction), "frozen_combined_transfer_probe"),
        "native_fast_update_gradient": (grad3(dfast), "actual_one_step_fast_update"),
        "native_slow_update_gradient": (grad3(dslow), "actual_one_step_slow_update"),
        "historical_common_bond": (h_c, "historical_internal_control"),
        "historical_difference_bond": (h_d, "historical_internal_control"),
        "historical_fast_bond": (h_uf, "historical_internal_control"),
        "historical_slow_bond": (h_us, "historical_internal_control"),
        "native_minus_historical_common_bond": (tuple(a-b for a,b in zip(bond_u,h_c)), "state_representation_residual"),
        "native_c_minus_historical_c_bond": (tuple(a-b for a,b in zip(bond_c,h_c)), "source_state_representation_residual"),
    }


def run_cluster(cluster):
    rho3 = DEC.local_rho3(cluster)
    hist = historical_fields(rho3)
    native = LOCAL.native_accumulated_vector_zero_flux(rho3)
    ref = hist["m10"]
    ref_rms = vrms(ref)
    rows = []
    for name, (vec, locality) in candidate_catalog(rho3, hist, native).items():
        cand_rms = vrms(vec)
        corr = vector_component_corr(vec, ref)
        cos, sign = vector_cosine(vec, ref)
        rows.append({
            "cluster_id": cluster["id"],
            "candidate": name,
            "locality": locality,
            "candidate_rms": cand_rms,
            "historical_m10_rms": ref_rms,
            "amplitude_ratio_to_m10": cand_rms / max(ref_rms, EPS),
            "component_correlation_to_m10": corr,
            "mean_vector_cosine_to_m10": cos,
            "positive_direction_fraction": sign,
            "native_c_integral_relative_error": float(native["c_state_integral_relative_error"]),
            "native_accumulation_converged": bool(native["converged"]),
        })
    return rows


def aggregate(rows):
    names = sorted(set(r["candidate"] for r in rows))
    out = []
    for name in names:
        rs = [r for r in rows if r["candidate"] == name]
        ratios = [r["amplitude_ratio_to_m10"] for r in rs]
        corrs = [r["component_correlation_to_m10"] for r in rs]
        cosines = [r["mean_vector_cosine_to_m10"] for r in rs]
        signs = [r["positive_direction_fraction"] for r in rs]
        ratio_mean = statistics.fmean(ratios)
        ratio_cv = cv(ratios)
        corr_mean = statistics.fmean(corrs)
        cosine_mean = statistics.fmean(cosines)
        sign_mean = statistics.fmean(signs)
        # Diagnostic score only. No candidate field is altered. Amplitude enters
        # through log distance to unity so orders-of-magnitude misses are penalised
        # symmetrically; morphology/direction and source stability are independent.
        amp_score = math.exp(-abs(math.log10(max(abs(ratio_mean), EPS))))
        stability_score = 1.0 / (1.0 + max(ratio_cv, 0.0))
        morphology_score = 0.5 * (max(corr_mean, -1.0) + 1.0)
        direction_score = 0.5 * (max(cosine_mean, -1.0) + 1.0)
        score = 0.30*amp_score + 0.20*stability_score + 0.25*morphology_score + 0.15*direction_score + 0.10*sign_mean
        out.append({
            "candidate": name,
            "locality": rs[0]["locality"],
            "amplitude_ratio_mean": ratio_mean,
            "amplitude_ratio_cv": ratio_cv,
            "component_correlation_mean": corr_mean,
            "vector_cosine_mean": cosine_mean,
            "positive_direction_fraction_mean": sign_mean,
            "diagnostic_score": score,
        })
    out.sort(key=lambda r: r["diagnostic_score"], reverse=True)
    for i, r in enumerate(out, 1):
        r["rank"] = i
    return out


def main() -> int:
    state = repo_state()
    clusters = list(BENCH.clusters())
    inventory = BENCH.inventory()
    ids = tuple(c["id"] for c in clusters)
    rows = []
    failures = []

    if ids == EXPECTED_CLUSTER_IDS and len(inventory) == 5 and all(x["exists"] for x in inventory):
        for cluster in clusters:
            try:
                rows.extend(run_cluster(cluster))
            except Exception as exc:
                failures.append({"cluster_id": cluster["id"], "error": f"{type(exc).__name__}: {exc}"})

    ranking = aggregate(rows) if rows else []
    expected_candidates = 18
    checks = {
        "canonical_five_cluster_inventory": ids == EXPECTED_CLUSTER_IDS,
        "all_five_local_fits_present": bool(len(inventory) == 5 and all(x["exists"] for x in inventory)),
        "all_five_clusters_completed": bool(len(rows) == 5 * expected_candidates and not failures),
        "broad_candidate_count_at_least_15": bool(ranking and len(ranking) >= 15),
        "all_measured_values_finite": bool(rows and all(
            all(math.isfinite(float(r[k])) for k in (
                "candidate_rms", "historical_m10_rms", "amplitude_ratio_to_m10",
                "component_correlation_to_m10", "mean_vector_cosine_to_m10",
                "positive_direction_fraction",
            )) for r in rows
        )),
        "native_c_state_integral_preserved": bool(rows and all(r["native_c_integral_relative_error"] <= 1e-12 for r in rows)),
        "native_accumulation_converged": bool(rows and all(r["native_accumulation_converged"] for r in rows)),
        "no_tracked_or_staged_changes": bool(not state["tracked_changes"] and not state["staged_changes"]),
    }
    gate = bool(all(checks.values()))
    if gate:
        status = "PROPAGATION_INTERFACE_WIDE_NET_SWEEP_EXECUTED"
    elif rows:
        status = "PROPAGATION_INTERFACE_WIDE_NET_SWEEP_PARTIAL_EXECUTION"
    else:
        status = "PROPAGATION_INTERFACE_WIDE_NET_SWEEP_NOT_ESTABLISHED"

    result = {
        "lab_id": LAB_ID,
        "status": status,
        "repo_state": state,
        "model": {
            "reference": "historical unit-loading M10 internal interface only",
            "observed_lensing_values_used": False,
            "network_access_used": False,
            "replacement_strength_scalar": None,
            "normalization_or_rescaling": False,
            "fit_or_tuning": False,
            "candidate_fed_to_G3D": False,
            "ranking_role": "diagnostic triage only; not a physical-law selector",
        },
        "candidate_count": len(ranking),
        "ranking": ranking,
        "rows": rows,
        "failures": failures,
        "checks": checks,
        "execution_gate_pass": gate,
    }

    print(LAB_ID)
    print(f"status={status}")
    print(f"head_sha={state['head_sha']}")
    print("network_access_used=false")
    print("observed_lensing_values_used=false")
    print("replacement_strength_scalar=none")
    print("native_response_rescaled=false")
    print("candidate_fed_to_G3D=false")
    print("fit_or_tuning=false")
    print(f"candidate_count={len(ranking)}")
    print()
    print("RANKING")
    for r in ranking:
        print(
            f"rank={r['rank']} candidate={r['candidate']} locality={r['locality']} "
            f"amp_ratio={r['amplitude_ratio_mean']:.12g} amp_cv={r['amplitude_ratio_cv']:.12g} "
            f"corr={r['component_correlation_mean']:.12g} cosine={r['vector_cosine_mean']:.12g} "
            f"sign={r['positive_direction_fraction_mean']:.12g} score={r['diagnostic_score']:.12g}"
        )
    print()
    print("CHECKS")
    for k,v in checks.items():
        print(f"{k}={str(v).lower()}")
    print(f"execution_gate_pass={str(gate).lower()}")
    print(json.dumps(result, sort_keys=True, separators=(",", ":"), default=str))
    return 0 if gate else 1


if __name__ == "__main__":
    raise SystemExit(main())
