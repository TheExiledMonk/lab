#!/usr/bin/env python3
"""PBUF FOUNDATION — INTERFACE-TO-INTERFACE SURVIVOR SWEEP 001.

Purpose
-------
Take the most plausible native/local survivor quantities from the wide-net sweep,
pass them through the SAME frozen PM1/PS2/M10 interface construction, and compare
interface-to-interface and LOS-to-LOS against the historical unit-loading M10
lane. This corrects the representation mismatch exposed by PR #91, where even
known historical bond ingredients correlated only weakly with final M10 before
pair/interface construction.

No observed lensing values are loaded. No candidate is fitted, normalized,
rescaled to M10, or fed into G3D. Historical M10 and LOS are internal structural
references only.
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
from pbuf.core import pair_enumeration as M05
from pbuf.core import pair_transfer as M08
from pbuf.core import midpoint_rasterization as M10
from pbuf.core import los_projection as M14
from pbuf.models import a8_state as A8
from pbuf.models import a8_pair_amplitude as PAIR
from pbuf.models import transverse_projector as PROJ
import pbuf.labs.foundation.m10_coverage_25pct_science001 as BASE
import pbuf.labs.foundation.m10_local_interface_decomposition_audit001 as DEC
import pbuf.labs.foundation.native_accumulated_full_lensing_local_benchmark001 as LOCAL
import pbuf.labs.foundation.native_accumulated_full_lensing001 as FULL

LAB_ID = "PBUF-FOUNDATION-INTERFACE-TO-INTERFACE-SURVIVOR-SWEEP-001"
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
    comps = [np.asarray(q, dtype=np.float64) for q in v]
    mag2 = np.zeros_like(comps[0])
    for q in comps:
        mag2 += q * q
    return float(np.sqrt(np.mean(mag2)))


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
    if float(np.std(x)) <= EPS or float(np.std(y)) <= EPS:
        return 0.0
    return float(np.corrcoef(x, y)[0, 1])


def component_corr(a, b) -> float:
    vals = [safe_corr(x, y) for x, y in zip(a, b)]
    vals = [v for v in vals if math.isfinite(v)]
    return float(statistics.fmean(vals)) if vals else float("nan")


def vector_cosine(a, b) -> tuple[float, float]:
    aa = [np.asarray(x, dtype=np.float64) for x in a]
    bb = [np.asarray(x, dtype=np.float64) for x in b]
    dot = sum(x*y for x, y in zip(aa, bb))
    am = np.sqrt(sum(x*x for x in aa))
    bm = np.sqrt(sum(x*x for x in bb))
    mask = (am > EPS) & (bm > EPS) & np.isfinite(dot)
    if not np.any(mask):
        return float("nan"), float("nan")
    c = dot[mask] / (am[mask] * bm[mask])
    return float(np.mean(c)), float(np.mean(c > 0.0))


def historical_state_and_interface(rho3: np.ndarray) -> dict:
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
    los = M14.project_vector_to_image_plane(*m10, los_axis="z")
    los_vec = (np.asarray(los["comp_1"], dtype=np.float64), np.asarray(los["comp_2"], dtype=np.float64))
    return {
        "u_fast": np.asarray(state["u_fast"], dtype=np.float64),
        "u_slow": np.asarray(state["u_slow"], dtype=np.float64),
        "c_state": np.asarray(state["c_state"], dtype=np.float64),
        "m10": m10,
        "los": los_vec,
    }


def positive_bond_amplitudes(field: np.ndarray, scale: float = 1.0) -> dict:
    f = np.asarray(field, dtype=np.float64)
    A_xp = np.zeros_like(f); A_yp = np.zeros_like(f); A_zp = np.zeros_like(f)
    A_xp[:, :, :-1] = scale * (f[:, :, 1:] - f[:, :, :-1])
    A_yp[:, :-1, :] = scale * (f[:, 1:, :] - f[:, :-1, :])
    A_zp[:-1, :, :] = scale * (f[1:, :, :] - f[:-1, :, :])
    return {"A_xp": A_xp, "A_yp": A_yp, "A_zp": A_zp}


def positive_traction_amplitudes(u: np.ndarray, scale: float = 1.0) -> dict:
    f = np.asarray(u, dtype=np.float64)
    dx = f[:, :, 1:] - f[:, :, :-1]
    dy = f[:, 1:, :] - f[:, :-1, :]
    dz = f[1:, :, :] - f[:-1, :, :]

    def sigma(q):
        frac2 = (q / FULL.EPSILON_MAX) ** 2
        if np.any(frac2 >= 1.0):
            raise RuntimeError("candidate bond exceeded bounded-strain domain")
        return FULL.K0 * q / (1.0 - frac2)

    A_xp = np.zeros_like(f); A_yp = np.zeros_like(f); A_zp = np.zeros_like(f)
    A_xp[:, :, :-1] = scale * sigma(dx)
    A_yp[:, :-1, :] = scale * sigma(dy)
    A_zp[:-1, :, :] = scale * sigma(dz)
    return {"A_xp": A_xp, "A_yp": A_yp, "A_zp": A_zp}


def m10_from_amplitudes(amplitudes: dict, geometry_scalar: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    shape = tuple(np.asarray(geometry_scalar).shape)
    pairs = M05.enumerate_internal_pairs(shape)
    ex, ey, ez, valid, gmag = PROJ.build_longitudinal_direction(np.asarray(geometry_scalar, dtype=np.float64))
    projector = PROJ.build_transverse_projector(ex, ey, ez)
    response = M08.build_pair_responses(
        pairs, amplitudes, projector,
        magnitude_formulation="PM1", pair_symmetrization="PS2",
    )
    interface = M10.rasterize_interface_field(response, shape)
    return (
        np.asarray(interface["Rx_3d_interface"], dtype=np.float64),
        np.asarray(interface["Ry_3d_interface"], dtype=np.float64),
        np.asarray(interface["Rz_3d_interface"], dtype=np.float64),
    )


def native_one_step_fields(u: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    # Frozen A8 zero-flux neighbor operator, used only to generate actual one-step
    # state increments for candidate pair amplitudes. No inferred coefficient.
    n = A8.neighbours6_face_zero_flux_3d(np.asarray(u, dtype=np.float64))
    fast = (A8.A8_INIT_DT * A8.A8_INIT_OMEGA * A8.A8_INIT_K) * (n - u)
    slow = (A8.A8_INIT_DT * A8.A8_INIT_SLOW_TIMESCALE) * (n - u)
    return np.asarray(fast), np.asarray(slow)


def candidate_catalog(native: dict, hist: dict) -> dict:
    c = np.asarray(native["c_state"], dtype=np.float64)
    u = np.asarray(native["accumulated"], dtype=np.float64)
    coef_fast = float(A8.A8_INIT_DT * A8.A8_INIT_OMEGA * A8.A8_INIT_K)
    coef_slow = float(A8.A8_INIT_DT * A8.A8_INIT_SLOW_TIMESCALE)
    fast_step, slow_step = native_one_step_fields(u)

    # Each entry is (pair amplitudes, geometry scalar, construction class).
    # The amplitude operation and geometry operation are kept explicit so that
    # candidate ranking does not silently normalize or fit either component.
    return {
        "native_u_bond_cgeom": (positive_bond_amplitudes(u), c, "native_bond_native_c_geometry"),
        "native_u_bond_ugeom": (positive_bond_amplitudes(u), u, "native_bond_native_u_geometry"),
        "native_traction_cgeom": (positive_traction_amplitudes(u), c, "native_traction_native_c_geometry"),
        "native_traction_ugeom": (positive_traction_amplitudes(u), u, "native_traction_native_u_geometry"),
        "native_slow_traction_cgeom": (positive_traction_amplitudes(u, coef_slow), c, "frozen_slow_traction_native_c_geometry"),
        "native_slow_traction_ugeom": (positive_traction_amplitudes(u, coef_slow), u, "frozen_slow_traction_native_u_geometry"),
        "native_fast_traction_cgeom": (positive_traction_amplitudes(u, coef_fast), c, "frozen_fast_traction_native_c_geometry"),
        "native_sum_traction_cgeom": (positive_traction_amplitudes(u, coef_fast + coef_slow), c, "frozen_sum_traction_native_c_geometry"),
        "native_c_bond_cgeom": (positive_bond_amplitudes(c), c, "native_c_bond_native_c_geometry"),
        "native_fast_step_bond_cgeom": (positive_bond_amplitudes(fast_step), c, "actual_fast_step_bond_native_c_geometry"),
        "native_slow_step_bond_cgeom": (positive_bond_amplitudes(slow_step), c, "actual_slow_step_bond_native_c_geometry"),
        "historical_common_control": (positive_bond_amplitudes(hist["c_state"]), hist["c_state"], "historical_common_bond_control"),
    }


def run_cluster(cluster: dict) -> list[dict]:
    rho3 = DEC.local_rho3(cluster)
    hist = historical_state_and_interface(rho3)
    native = LOCAL.native_accumulated_vector_zero_flux(rho3)
    ref_m10 = hist["m10"]
    ref_los = hist["los"]
    ref_m10_rms = vrms(ref_m10)
    ref_los_rms = vrms(ref_los)
    rows = []
    for name, (amps, geom, locality) in candidate_catalog(native, hist).items():
        m10 = m10_from_amplitudes(amps, geom)
        los_raw = M14.project_vector_to_image_plane(*m10, los_axis="z")
        los = (np.asarray(los_raw["comp_1"], dtype=np.float64), np.asarray(los_raw["comp_2"], dtype=np.float64))
        m10_cos, m10_sign = vector_cosine(m10, ref_m10)
        los_cos, los_sign = vector_cosine(los, ref_los)
        rows.append({
            "cluster_id": cluster["id"],
            "candidate": name,
            "locality": locality,
            "m10_rms": vrms(m10),
            "historical_m10_rms": ref_m10_rms,
            "m10_amplitude_ratio": vrms(m10) / max(ref_m10_rms, EPS),
            "m10_component_correlation": component_corr(m10, ref_m10),
            "m10_vector_cosine": m10_cos,
            "m10_positive_direction_fraction": m10_sign,
            "los_rms": vrms(los),
            "historical_los_rms": ref_los_rms,
            "los_amplitude_ratio": vrms(los) / max(ref_los_rms, EPS),
            "los_component_correlation": component_corr(los, ref_los),
            "los_vector_cosine": los_cos,
            "los_positive_direction_fraction": los_sign,
            "native_c_integral_relative_error": float(native["c_state_integral_relative_error"]),
            "native_accumulation_converged": bool(native["converged"]),
        })
    return rows


def aggregate(rows: list[dict]) -> list[dict]:
    names = sorted(set(r["candidate"] for r in rows))
    out = []
    for name in names:
        rs = [r for r in rows if r["candidate"] == name]
        mrat = [r["m10_amplitude_ratio"] for r in rs]
        lrat = [r["los_amplitude_ratio"] for r in rs]
        mcorr = statistics.fmean(r["m10_component_correlation"] for r in rs)
        mcos = statistics.fmean(r["m10_vector_cosine"] for r in rs)
        msign = statistics.fmean(r["m10_positive_direction_fraction"] for r in rs)
        lcorr = statistics.fmean(r["los_component_correlation"] for r in rs)
        lcos = statistics.fmean(r["los_vector_cosine"] for r in rs)
        lsign = statistics.fmean(r["los_positive_direction_fraction"] for r in rs)
        mmean = statistics.fmean(mrat); lmean = statistics.fmean(lrat)
        mcv = cv(mrat); lcv = cv(lrat)

        amp_score = 0.5 * math.exp(-abs(math.log10(max(abs(mmean), EPS)))) + 0.5 * math.exp(-abs(math.log10(max(abs(lmean), EPS))))
        stability = 0.5/(1.0 + max(mcv, 0.0)) + 0.5/(1.0 + max(lcv, 0.0))
        morphology = 0.25*(mcorr+1.0) + 0.25*(lcorr+1.0)
        direction = 0.25*(mcos+1.0) + 0.25*(lcos+1.0)
        sign = 0.5*(msign + lsign)
        score = 0.25*amp_score + 0.20*stability + 0.25*morphology + 0.20*direction + 0.10*sign
        out.append({
            "candidate": name,
            "locality": rs[0]["locality"],
            "m10_amplitude_ratio_mean": mmean,
            "m10_amplitude_ratio_cv": mcv,
            "m10_correlation_mean": mcorr,
            "m10_cosine_mean": mcos,
            "m10_sign_mean": msign,
            "los_amplitude_ratio_mean": lmean,
            "los_amplitude_ratio_cv": lcv,
            "los_correlation_mean": lcorr,
            "los_cosine_mean": lcos,
            "los_sign_mean": lsign,
            "diagnostic_score": score,
        })
    out.sort(key=lambda x: x["diagnostic_score"], reverse=True)
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
    checks = {
        "canonical_five_cluster_inventory": ids == EXPECTED_CLUSTER_IDS,
        "all_five_local_fits_present": bool(len(inventory) == 5 and all(x["exists"] for x in inventory)),
        "all_five_clusters_completed": bool(len(rows) == 5 * 12 and not failures),
        "candidate_count_at_least_10": bool(len(ranking) >= 10),
        "all_measured_values_finite": bool(rows and all(
            all(math.isfinite(float(v)) for k, v in r.items() if k not in ("cluster_id", "candidate", "locality") and not isinstance(v, bool))
            for r in rows
        )),
        "native_c_state_integral_preserved": bool(rows and all(r["native_c_integral_relative_error"] <= 1.0e-12 for r in rows)),
        "native_accumulation_converged": bool(rows and all(r["native_accumulation_converged"] for r in rows)),
        "no_tracked_or_staged_changes": bool(not state["tracked_changes"] and not state["staged_changes"]),
    }
    execution_gate_pass = bool(all(checks.values()))
    if execution_gate_pass:
        status = "INTERFACE_TO_INTERFACE_SURVIVOR_SWEEP_EXECUTED"
    elif rows:
        status = "INTERFACE_TO_INTERFACE_SURVIVOR_SWEEP_PARTIAL_EXECUTION"
    else:
        status = "INTERFACE_TO_INTERFACE_SURVIVOR_SWEEP_NOT_ESTABLISHED"

    result = {
        "lab_id": LAB_ID,
        "status": status,
        "repo_state": state,
        "candidate_count": len(ranking),
        "model": {
            "reference": "historical unit-loading PM1/PS2/M10 interface and its frozen LOS projection only",
            "candidate_path": "native/local pair amplitude -> frozen PM1/PS2 -> frozen M10 -> frozen LOS projection",
            "observed_lensing_values_used": False,
            "replacement_strength_scalar": None,
            "normalization_or_rescaling": False,
            "fit_or_tuning": False,
            "candidate_fed_to_G3D": False,
        },
        "ranking": ranking,
        "rows": rows,
        "failures": failures,
        "checks": checks,
        "execution_gate_pass": execution_gate_pass,
        "interpretation_rule": "Rank interface-to-interface and LOS-to-LOS structural survivors only. Do not promote, fit, rescale, or feed any candidate into G3D in this audit.",
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
    print("RANKING")
    for r in ranking:
        print(
            f"rank={r['rank']} candidate={r['candidate']} locality={r['locality']} "
            f"m10_amp={r['m10_amplitude_ratio_mean']:.12g} m10_cv={r['m10_amplitude_ratio_cv']:.12g} "
            f"m10_corr={r['m10_correlation_mean']:.12g} m10_cos={r['m10_cosine_mean']:.12g} m10_sign={r['m10_sign_mean']:.12g} "
            f"los_amp={r['los_amplitude_ratio_mean']:.12g} los_cv={r['los_amplitude_ratio_cv']:.12g} "
            f"los_corr={r['los_correlation_mean']:.12g} los_cos={r['los_cosine_mean']:.12g} los_sign={r['los_sign_mean']:.12g} "
            f"score={r['diagnostic_score']:.12g}"
        )
    print("CHECKS")
    for k, v in checks.items():
        print(f"{k}={str(v).lower()}")
    print(f"execution_gate_pass={str(execution_gate_pass).lower()}")
    for f in failures:
        print(f"failure_cluster={f['cluster_id']} error={f['error']}")
    print(json.dumps(result, sort_keys=True, separators=(",", ":"), default=str))
    return 0 if execution_gate_pass else 1


if __name__ == "__main__":
    raise SystemExit(main())
