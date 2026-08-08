#!/usr/bin/env python3
"""PBUF FOUNDATION — NATIVE FORWARD TRANSFER OPERATOR AUDIT 001.

Purpose
-------
Test the next forward-derivation step after the frozen-factor decomposition.
Rather than multiplying already-frozen factors by inspection, take the native
bounded-strain accumulated equilibrium and pass it through the actual frozen
A8 pair-transfer + PM1/PS2 + M10 machinery.

This audit makes exactly one explicit structural embedding assumption:

    u_fast = u_slow = u_accum

for the purpose of asking what the frozen A8 local transfer law would produce
from the native accumulated scalar field if both historical response channels
sample the same settled medium configuration.

The assumption is diagnostic only. It is not promoted to PBUF physics, is not
fed into G3D, and is not fitted to the historical M10 amplitude.

The exact A8 pair law is then used unchanged:

    A_ij = (dt*omega*K) * Delta u_fast
         + (dt*tau_slow) * Delta u_slow

followed by the existing PM1/PS2 response and M10 midpoint rasterisation.
No hand-applied 0.5 coupling factor is inserted. If the frozen fast/slow
couplings do not enter this forward construction algebraically, the audit must
show that rather than forcing the previous 0.00825 factor combination.
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
from pbuf.models import a8_pair_amplitude as PAIR
from pbuf.models import a8_state as A8
from pbuf.models import transverse_projector as PROJ
import pbuf.labs.foundation.m10_local_interface_decomposition_audit001 as DEC
import pbuf.labs.foundation.native_accumulated_full_lensing_local_benchmark001 as LOCAL

LAB_ID = "PBUF-FOUNDATION-NATIVE-FORWARD-TRANSFER-OPERATOR-AUDIT-001"
EXPECTED_CLUSTER_IDS = (
    "Abell2744", "MACS0416", "MACS1149", "AbellS1063", "Abell370"
)


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


def ratio(a: float, b: float) -> float:
    return float(a / max(abs(b), 1.0e-30))


def vector_rms(vector) -> float:
    comps = [np.asarray(v, dtype=np.float64) for v in vector]
    mag2 = np.zeros_like(comps[0])
    for c in comps:
        mag2 += c * c
    return float(np.sqrt(np.mean(mag2)))


def concat_pair_amplitudes(amps: dict) -> np.ndarray:
    return np.concatenate((
        np.asarray(amps["A_xp"], dtype=np.float64)[:, :, :-1].ravel(),
        np.asarray(amps["A_yp"], dtype=np.float64)[:, :-1, :].ravel(),
        np.asarray(amps["A_zp"], dtype=np.float64)[:-1, :, :].ravel(),
    ))


def native_forward_operator(rho3: np.ndarray) -> dict:
    build = LOCAL.native_accumulated_vector_zero_flux(rho3)
    c = np.asarray(build["c_state"], dtype=np.float64)
    u = np.asarray(build["accumulated"], dtype=np.float64)
    shape = tuple(u.shape)

    # Explicit diagnostic embedding: both frozen historical transfer channels
    # sample the same settled native scalar field. No noise, strength or extra
    # coupling multiplier is introduced.
    us = u.copy()
    uf = u.copy()

    pairs = M05.enumerate_internal_pairs(shape)
    ex, ey, ez, valid, _gmag = PROJ.build_longitudinal_direction(c)
    projector = PROJ.build_transverse_projector(ex, ey, ez)

    amps = PAIR.compute_a8_pair_amplitudes(us, uf, c, pairs)
    response = M08.build_pair_responses(
        pairs, amps, projector,
        magnitude_formulation="PM1", pair_symmetrization="PS2",
    )
    interface = M10.rasterize_interface_field(response, shape)
    m10_vec = (
        np.asarray(interface["Rx_3d_interface"], dtype=np.float64),
        np.asarray(interface["Ry_3d_interface"], dtype=np.float64),
        np.asarray(interface["Rz_3d_interface"], dtype=np.float64),
    )

    pair_amp = concat_pair_amplitudes(amps)
    bonds = DEC.concat_positive_bonds(u)
    coef_fast = float(A8.A8_INIT_DT * A8.A8_INIT_OMEGA * A8.A8_INIT_K)
    coef_slow = float(A8.A8_INIT_DT * A8.A8_INIT_SLOW_TIMESCALE)
    exact = (coef_fast + coef_slow) * bonds
    exact_relerr = ratio(rms(pair_amp - exact), rms(exact))

    # Native bounded-strain traction is retained only as an amplitude reference.
    frac = np.abs(bonds) / DEC.FULL.EPSILON_MAX
    if np.any(frac >= 1.0):
        raise RuntimeError("native bond exceeded bounded-strain domain")
    traction = DEC.FULL.K0 * bonds / (1.0 - frac * frac)

    return {
        "c_rms": rms(c),
        "u_rms": rms(u),
        "bond_rms": rms(bonds),
        "traction_rms": rms(traction),
        "pair_amp_rms": rms(pair_amp),
        "m10_rms": vector_rms(m10_vec),
        "pair_formula_relative_rms_error": exact_relerr,
        "coef_fast": coef_fast,
        "coef_slow": coef_slow,
        "coef_sum": coef_fast + coef_slow,
        "valid_longitudinal_count": int(np.count_nonzero(valid)),
        "c_state_integral_relative_error": float(build["c_state_integral_relative_error"]),
        "accumulation_converged": bool(build["converged"]),
    }


def run_cluster(cluster: dict) -> dict:
    rho3 = DEC.local_rho3(cluster)
    historical = DEC.unit_decomposition(rho3)
    native = native_forward_operator(rho3)

    return {
        "cluster_id": cluster["id"],
        "native_traction_rms": native["traction_rms"],
        "forward_pair_amp_rms": native["pair_amp_rms"],
        "forward_m10_rms": native["m10_rms"],
        "historical_pair_amp_rms": historical["pair_amp_rms"],
        "historical_m10_rms": historical["m10_vector_rms"],
        "forward_pair_over_native_traction": ratio(native["pair_amp_rms"], native["traction_rms"]),
        "forward_m10_over_native_traction": ratio(native["m10_rms"], native["traction_rms"]),
        "forward_pair_over_historical_pair": ratio(native["pair_amp_rms"], historical["pair_amp_rms"]),
        "forward_m10_over_historical_m10": ratio(native["m10_rms"], historical["m10_vector_rms"]),
        "forward_m10_over_forward_pair": ratio(native["m10_rms"], native["pair_amp_rms"]),
        "native_bond_over_traction": ratio(native["bond_rms"], native["traction_rms"]),
        "coef_fast": native["coef_fast"],
        "coef_slow": native["coef_slow"],
        "coef_sum": native["coef_sum"],
        "pair_formula_relative_rms_error": native["pair_formula_relative_rms_error"],
        "native_c_state_integral_relative_error": native["c_state_integral_relative_error"],
        "native_accumulation_converged": native["accumulation_converged"],
    }


def cv(values: list[float]) -> float:
    if len(values) < 2:
        return float("nan")
    m = statistics.fmean(values)
    return float(statistics.pstdev(values) / max(abs(m), 1.0e-30))


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
                rows.append(run_cluster(cluster))
            except Exception as exc:
                failures.append({"cluster_id": cluster["id"], "error": f"{type(exc).__name__}: {exc}"})

    pair_forward = [r["forward_pair_over_native_traction"] for r in rows]
    m10_forward = [r["forward_m10_over_native_traction"] for r in rows]
    pair_hist_ratio = [r["forward_pair_over_historical_pair"] for r in rows]
    m10_hist_ratio = [r["forward_m10_over_historical_m10"] for r in rows]

    checks = {
        "canonical_five_cluster_inventory": ids == EXPECTED_CLUSTER_IDS,
        "all_five_local_fits_present": bool(len(inventory) == 5 and all(x["exists"] for x in inventory)),
        "all_five_clusters_completed": bool(len(rows) == 5 and not failures),
        "all_measured_values_finite": bool(rows and all(
            all(math.isfinite(float(v)) for k, v in r.items() if k != "cluster_id" and not isinstance(v, bool))
            for r in rows
        )),
        "exact_forward_pair_formula_reproduced": bool(rows and all(r["pair_formula_relative_rms_error"] <= 1.0e-14 for r in rows)),
        "native_c_state_integral_preserved": bool(rows and all(r["native_c_state_integral_relative_error"] <= 1.0e-12 for r in rows)),
        "native_accumulation_converged": bool(rows and all(r["native_accumulation_converged"] for r in rows)),
        "no_tracked_or_staged_changes": bool(not state["tracked_changes"] and not state["staged_changes"]),
    }
    execution_gate_pass = bool(all(checks.values()))

    if execution_gate_pass:
        status = "NATIVE_FORWARD_TRANSFER_OPERATOR_AUDIT_EXECUTED"
    elif rows:
        status = "NATIVE_FORWARD_TRANSFER_OPERATOR_AUDIT_PARTIAL_EXECUTION"
    else:
        status = "NATIVE_FORWARD_TRANSFER_OPERATOR_AUDIT_NOT_ESTABLISHED"

    summary = {
        "forward_pair_over_native_traction_mean": statistics.fmean(pair_forward) if pair_forward else float("nan"),
        "forward_pair_over_native_traction_cv": cv(pair_forward),
        "forward_m10_over_native_traction_mean": statistics.fmean(m10_forward) if m10_forward else float("nan"),
        "forward_m10_over_native_traction_cv": cv(m10_forward),
        "forward_pair_over_historical_pair_mean": statistics.fmean(pair_hist_ratio) if pair_hist_ratio else float("nan"),
        "forward_m10_over_historical_m10_mean": statistics.fmean(m10_hist_ratio) if m10_hist_ratio else float("nan"),
        "frozen_coef_fast": float(A8.A8_INIT_DT * A8.A8_INIT_OMEGA * A8.A8_INIT_K),
        "frozen_coef_slow": float(A8.A8_INIT_DT * A8.A8_INIT_SLOW_TIMESCALE),
        "frozen_coef_sum": float(A8.A8_INIT_DT * A8.A8_INIT_OMEGA * A8.A8_INIT_K + A8.A8_INIT_DT * A8.A8_INIT_SLOW_TIMESCALE),
    }

    result = {
        "lab_id": LAB_ID,
        "status": status,
        "repo_state": state,
        "model": {
            "explicit_embedding_assumption": "u_fast=u_slow=u_accum for diagnostic forward A8 transfer only",
            "pair_law": "A_ij=(dt*omega*K)*Delta_u_fast+(dt*tau_slow)*Delta_u_slow",
            "PM1_PS2_M10": "existing frozen machinery unchanged",
            "hand_applied_coupling_factor": False,
            "replacement_strength_scalar": None,
            "normalization_or_rescaling": False,
            "observed_lensing_values_used": False,
            "candidate_fed_to_G3D": False,
            "fit_or_tuning": False,
        },
        "summary": summary,
        "rows": rows,
        "failures": failures,
        "checks": checks,
        "execution_gate_pass": execution_gate_pass,
        "interpretation_rule": "Determine what the exact frozen A8/PM1/PS2/M10 operator produces from the native accumulated field under the stated equal-channel embedding. Do not add factors to improve agreement.",
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
    print("embedding_assumption=u_fast_equals_u_slow_equals_u_accum")
    print("hand_applied_coupling_factor=false")
    print(f"coef_fast={summary['frozen_coef_fast']:.12g}")
    print(f"coef_slow={summary['frozen_coef_slow']:.12g}")
    print(f"coef_sum={summary['frozen_coef_sum']:.12g}")
    print()
    print("CLUSTERS")
    for r in rows:
        print(
            f"cluster={r['cluster_id']} "
            f"forward_pair_over_traction={r['forward_pair_over_native_traction']:.12g} "
            f"forward_m10_over_traction={r['forward_m10_over_native_traction']:.12g} "
            f"forward_pair_over_historical={r['forward_pair_over_historical_pair']:.12g} "
            f"forward_m10_over_historical={r['forward_m10_over_historical_m10']:.12g} "
            f"m10_over_pair={r['forward_m10_over_forward_pair']:.12g}"
        )
    for f in failures:
        print(f"failure_cluster={f['cluster_id']} error={f['error']}")
    print()
    print("SUMMARY")
    for k, v in summary.items():
        print(f"{k}={v:.12g}")
    print("CHECKS")
    for k, v in checks.items():
        print(f"{k}={str(v).lower()}")
    print(f"execution_gate_pass={str(execution_gate_pass).lower()}")
    print(json.dumps(result, sort_keys=True, separators=(",", ":"), default=str))
    return 0 if execution_gate_pass else 1


if __name__ == "__main__":
    raise SystemExit(main())
