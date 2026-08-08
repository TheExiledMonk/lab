#!/usr/bin/env python3
"""PBUF FOUNDATION — NATIVE CHANNEL TRANSFER CLOSURE SWEEP 001.

Purpose
-------
Test several forward native fast/slow transfer constructions at once around the
strong PR #92 result that native c_state geometry closely reproduces historical
interface morphology while the unweighted c-bond amplitude is ~35x too large.

No inferred effective coefficient is inserted.  The frozen A8 coefficients and
zero-flux native transport are used directly.  Every candidate is passed through
the same frozen PM1/PS2/M10 and LOS representation before comparison with the
historical unit-loading internal reference.
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
from pbuf.core import los_projection as M14
from pbuf.models import a8_state as A8
import pbuf.labs.foundation.m10_local_interface_decomposition_audit001 as DEC
import pbuf.labs.foundation.interface_to_interface_survivor_sweep001 as S92
import pbuf.labs.foundation.native_accumulated_full_lensing_local_benchmark001 as LOCAL

LAB_ID = "PBUF-FOUNDATION-NATIVE-CHANNEL-TRANSFER-CLOSURE-SWEEP-001"
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


def cv(values) -> float:
    vals = [float(v) for v in values if math.isfinite(float(v))]
    if len(vals) < 2:
        return float("nan")
    m = statistics.fmean(vals)
    return float(statistics.pstdev(vals) / max(abs(m), EPS))


def combine_amplitudes(a: dict, b: dict) -> dict:
    return {k: np.asarray(a[k], dtype=np.float64) + np.asarray(b[k], dtype=np.float64)
            for k in ("A_xp", "A_yp", "A_zp")}


def scale_amplitudes(a: dict, scale: float) -> dict:
    return {k: float(scale) * np.asarray(a[k], dtype=np.float64)
            for k in ("A_xp", "A_yp", "A_zp")}


def native_terminal_channels(rho3: np.ndarray) -> dict:
    """Run the frozen A8 transport with equal unit native initialization and zero flux."""
    u0 = np.asarray(rho3, dtype=np.float64).copy()
    us, uf, history = A8.evolve_a8_transport_3d(
        u0.copy(), u0.copy(), stencil="N6", boundary="zero_flux"
    )
    us = np.asarray(us, dtype=np.float64)
    uf = np.asarray(uf, dtype=np.float64)
    c = 0.5 * (uf + us)
    d = 0.5 * (uf - us)
    return {"u_slow": us, "u_fast": uf, "c": c, "d": d,
            "history_c": np.asarray(history[-1], dtype=np.float64)}


def candidate_catalog(native: dict, channels: dict, hist: dict) -> dict:
    c = np.asarray(native["c_state"], dtype=np.float64)
    u = np.asarray(native["accumulated"], dtype=np.float64)
    uf = channels["u_fast"]
    us = channels["u_slow"]
    cm = channels["c"]
    dm = channels["d"]

    cf = float(A8.A8_INIT_DT * A8.A8_INIT_OMEGA * A8.A8_INIT_K)
    cs = float(A8.A8_INIT_DT * A8.A8_INIT_SLOW_TIMESCALE)
    cc = cf + cs
    cd = cf - cs

    b_c = S92.positive_bond_amplitudes(c)
    b_cm = S92.positive_bond_amplitudes(cm)
    b_dm = S92.positive_bond_amplitudes(dm)
    b_uf = S92.positive_bond_amplitudes(uf)
    b_us = S92.positive_bond_amplitudes(us)

    common = scale_amplitudes(b_cm, cc)
    difference = scale_amplitudes(b_dm, cd)
    full_modal = combine_amplitudes(common, difference)
    full_direct = combine_amplitudes(scale_amplitudes(b_uf, cf), scale_amplitudes(b_us, cs))

    hist_full_direct = combine_amplitudes(
        scale_amplitudes(S92.positive_bond_amplitudes(hist["u_fast"]), cf),
        scale_amplitudes(S92.positive_bond_amplitudes(hist["u_slow"]), cs),
    )

    # Wide set in one run.  All factors below are already frozen A8 coefficients;
    # no diagnostic/inferred 0.02925-type coefficient is applied.
    return {
        "native_c_unweighted_control": (b_c, c, "native_c_bond_control"),
        "native_c_common_coef": (scale_amplitudes(b_c, cc), c, "native_c_times_frozen_common_coefficient"),
        "native_terminal_common_only": (common, c, "zero_flux_terminal_common_mode_only"),
        "native_terminal_difference_only": (difference, c, "zero_flux_terminal_difference_mode_only"),
        "native_terminal_full_modal": (full_modal, c, "zero_flux_terminal_common_plus_difference"),
        "native_terminal_full_direct": (full_direct, c, "zero_flux_terminal_exact_fast_plus_slow_pair_law"),
        "native_terminal_fast_only": (scale_amplitudes(b_uf, cf), c, "zero_flux_terminal_fast_channel_only"),
        "native_terminal_slow_only": (scale_amplitudes(b_us, cs), c, "zero_flux_terminal_slow_channel_only"),
        "native_traction_slow_control": (S92.positive_traction_amplitudes(u, cs), c, "bounded_traction_slow_control"),
        "historical_full_transfer_positive_control": (hist_full_direct, hist["c_state"], "historical_exact_pair_transfer_control"),
    }


def run_cluster(cluster: dict) -> tuple[list[dict], dict]:
    rho3 = DEC.local_rho3(cluster)
    hist = S92.historical_state_and_interface(rho3)
    native = LOCAL.native_accumulated_vector_zero_flux(rho3)
    channels = native_terminal_channels(rho3)

    # The independently returned native c_state should be the terminal common mode
    # of the same frozen zero-flux A8 evolution.  This is a gate, not a rescaling.
    c_rel = S92.rms(channels["c"] - native["c_state"]) / max(S92.rms(native["c_state"]), EPS)
    history_rel = S92.rms(channels["history_c"] - channels["c"]) / max(S92.rms(channels["c"]), EPS)

    ref_m10 = hist["m10"]
    ref_los = hist["los"]
    ref_m10_rms = S92.vrms(ref_m10)
    ref_los_rms = S92.vrms(ref_los)
    rows = []

    for name, (amps, geom, locality) in candidate_catalog(native, channels, hist).items():
        m10 = S92.m10_from_amplitudes(amps, geom)
        lp = M14.project_vector_to_image_plane(*m10, los_axis="z")
        los = (np.asarray(lp["comp_1"], dtype=np.float64), np.asarray(lp["comp_2"], dtype=np.float64))
        mc, ms = S92.vector_cosine(m10, ref_m10)
        lc, ls = S92.vector_cosine(los, ref_los)
        rows.append({
            "cluster_id": cluster["id"], "candidate": name, "locality": locality,
            "m10_amplitude_ratio": S92.vrms(m10) / max(ref_m10_rms, EPS),
            "m10_component_correlation": S92.component_corr(m10, ref_m10),
            "m10_vector_cosine": mc, "m10_positive_direction_fraction": ms,
            "los_amplitude_ratio": S92.vrms(los) / max(ref_los_rms, EPS),
            "los_component_correlation": S92.component_corr(los, ref_los),
            "los_vector_cosine": lc, "los_positive_direction_fraction": ls,
            "native_c_terminal_common_relative_rms_error": c_rel,
            "terminal_history_common_relative_rms_error": history_rel,
            "native_c_integral_relative_error": float(native["c_state_integral_relative_error"]),
            "native_accumulation_converged": bool(native["converged"]),
        })
    return rows, {"cluster_id": cluster["id"], "c_terminal_relerr": c_rel,
                  "history_terminal_relerr": history_rel}


def aggregate(rows: list[dict]) -> list[dict]:
    out = []
    for name in sorted(set(r["candidate"] for r in rows)):
        rs = [r for r in rows if r["candidate"] == name]
        mrat = [r["m10_amplitude_ratio"] for r in rs]
        lrat = [r["los_amplitude_ratio"] for r in rs]
        mcorr = statistics.fmean(r["m10_component_correlation"] for r in rs)
        lcorr = statistics.fmean(r["los_component_correlation"] for r in rs)
        mcos = statistics.fmean(r["m10_vector_cosine"] for r in rs)
        lcos = statistics.fmean(r["los_vector_cosine"] for r in rs)
        msign = statistics.fmean(r["m10_positive_direction_fraction"] for r in rs)
        lsign = statistics.fmean(r["los_positive_direction_fraction"] for r in rs)
        mm = statistics.fmean(mrat); lm = statistics.fmean(lrat)
        mcv = cv(mrat); lcv = cv(lrat)
        amp_score = 0.5 * (math.exp(-abs(math.log10(max(abs(mm), EPS)))) + math.exp(-abs(math.log10(max(abs(lm), EPS)))))
        stability = 0.5 * (1.0/(1.0+mcv) + 1.0/(1.0+lcv))
        morphology = 0.25 * (mcorr + lcorr + 2.0)
        direction = 0.25 * (mcos + lcos + 2.0)
        sign = 0.5 * (msign + lsign)
        score = 0.30*amp_score + 0.15*stability + 0.25*morphology + 0.20*direction + 0.10*sign
        out.append({
            "candidate": name, "locality": rs[0]["locality"],
            "m10_amplitude_ratio_mean": mm, "m10_amplitude_ratio_cv": mcv,
            "m10_component_correlation_mean": mcorr, "m10_vector_cosine_mean": mcos,
            "m10_positive_direction_fraction_mean": msign,
            "los_amplitude_ratio_mean": lm, "los_amplitude_ratio_cv": lcv,
            "los_component_correlation_mean": lcorr, "los_vector_cosine_mean": lcos,
            "los_positive_direction_fraction_mean": lsign,
            "diagnostic_score": score,
        })
    out.sort(key=lambda r: r["diagnostic_score"], reverse=True)
    for i, r in enumerate(out, 1): r["rank"] = i
    return out


def main() -> int:
    state = repo_state()
    clusters = list(BENCH.clusters())
    inventory = BENCH.inventory()
    ids = tuple(c["id"] for c in clusters)
    rows, channel_checks, failures = [], [], []

    if ids == EXPECTED_CLUSTER_IDS and len(inventory) == 5 and all(x["exists"] for x in inventory):
        for cluster in clusters:
            try:
                rr, cc = run_cluster(cluster); rows.extend(rr); channel_checks.append(cc)
            except Exception as exc:
                failures.append({"cluster_id": cluster["id"], "error": f"{type(exc).__name__}: {exc}"})

    ranking = aggregate(rows) if rows else []
    checks = {
        "canonical_five_cluster_inventory": ids == EXPECTED_CLUSTER_IDS,
        "all_five_local_fits_present": bool(len(inventory) == 5 and all(x["exists"] for x in inventory)),
        "all_five_clusters_completed": bool(len(channel_checks) == 5 and not failures),
        "broad_candidate_count_at_least_8": bool(len(ranking) >= 8),
        "all_measured_values_finite": bool(rows and all(all(math.isfinite(float(v)) for k,v in r.items() if k not in ("cluster_id","candidate","locality") and not isinstance(v,bool)) for r in rows)),
        "native_terminal_common_matches_native_c": bool(channel_checks and all(x["c_terminal_relerr"] <= 1.0e-12 for x in channel_checks)),
        "terminal_history_matches_terminal_common": bool(channel_checks and all(x["history_terminal_relerr"] <= 1.0e-12 for x in channel_checks)),
        "native_c_state_integral_preserved": bool(rows and all(r["native_c_integral_relative_error"] <= 1.0e-12 for r in rows)),
        "native_accumulation_converged": bool(rows and all(r["native_accumulation_converged"] for r in rows)),
        "no_tracked_or_staged_changes": bool(not state["tracked_changes"] and not state["staged_changes"]),
    }
    passed = bool(all(checks.values()))
    status = "NATIVE_CHANNEL_TRANSFER_CLOSURE_SWEEP_EXECUTED" if passed else ("NATIVE_CHANNEL_TRANSFER_CLOSURE_SWEEP_PARTIAL_EXECUTION" if rows else "NATIVE_CHANNEL_TRANSFER_CLOSURE_SWEEP_NOT_ESTABLISHED")

    result = {
        "lab_id": LAB_ID, "status": status, "repo_state": state,
        "model": {
            "question": "Does the actual frozen native zero-flux fast/slow channel evolution generate the transfer needed by the native c_state interface geometry?",
            "reference": "historical unit-loading M10 and LOS internal structural controls only",
            "inferred_effective_coefficient_applied": False,
            "replacement_strength_scalar": None,
            "normalization_or_rescaling": False,
            "fit_or_tuning": False,
            "observed_lensing_values_used": False,
            "candidate_fed_to_G3D": False,
        },
        "frozen_coefficients": {"fast": float(A8.A8_INIT_DT*A8.A8_INIT_OMEGA*A8.A8_INIT_K), "slow": float(A8.A8_INIT_DT*A8.A8_INIT_SLOW_TIMESCALE)},
        "candidate_count": len(ranking), "ranking": ranking, "rows": rows,
        "channel_checks": channel_checks, "failures": failures, "checks": checks,
        "execution_gate_pass": passed,
        "interpretation_rule": "Rank forward native channel-transfer paths only. Do not promote, fit, rescale, or feed any candidate into G3D in this audit.",
    }

    print(LAB_ID); print(f"status={status}"); print(f"head_sha={state['head_sha']}")
    print("network_access_used=false"); print("observed_lensing_values_used=false")
    print("replacement_strength_scalar=none"); print("inferred_effective_coefficient_applied=false")
    print("native_response_rescaled=false"); print("candidate_fed_to_G3D=false"); print("fit_or_tuning=false")
    print(f"candidate_count={len(ranking)}")
    print("RANKING")
    for r in ranking:
        print(f"rank={r['rank']} candidate={r['candidate']} locality={r['locality']} m10_amp={r['m10_amplitude_ratio_mean']:.12g} m10_cv={r['m10_amplitude_ratio_cv']:.12g} m10_corr={r['m10_component_correlation_mean']:.12g} m10_cos={r['m10_vector_cosine_mean']:.12g} m10_sign={r['m10_positive_direction_fraction_mean']:.12g} los_amp={r['los_amplitude_ratio_mean']:.12g} los_cv={r['los_amplitude_ratio_cv']:.12g} los_corr={r['los_component_correlation_mean']:.12g} los_cos={r['los_vector_cosine_mean']:.12g} los_sign={r['los_positive_direction_fraction_mean']:.12g} score={r['diagnostic_score']:.12g}")
    print("CHECKS")
    for k,v in checks.items(): print(f"{k}={str(v).lower()}")
    print(f"execution_gate_pass={str(passed).lower()}")
    print(json.dumps(result, sort_keys=True, separators=(",",":"), default=str))
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
