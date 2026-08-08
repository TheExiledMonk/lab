#!/usr/bin/env python3
"""PBUF FOUNDATION — NATIVE LOCAL DYNAMIC RESPONSE AUDIT 001.

Purpose
-------
Test the specific hypothesis raised after the M10 local-interface decomposition:
the bounded-strain solver gives a static accumulated equilibrium, while the
historical propagation interface is built from a local per-update transfer.

This lab does NOT create a new propagation law.  It compares the already-derived
native local bond traction with the already-frozen historical transfer scales and
reports the source-by-source effective transfer coefficient that would map the
native local traction onto the historical pair-amplitude and M10 scales.

Nothing inferred here is fed back into the model.  In particular, the inferred
coefficients are diagnostics, not fitted parameters.
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
import pbuf.labs.foundation.m10_local_interface_decomposition_audit001 as DEC
import pbuf.labs.foundation.native_accumulated_full_lensing001 as FULL

LAB_ID = "PBUF-FOUNDATION-NATIVE-LOCAL-DYNAMIC-RESPONSE-AUDIT-001"
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


def ratio(a: float, b: float) -> float:
    return float(a / max(abs(b), 1.0e-30))


def cv(values: list[float]) -> float:
    vals = [float(v) for v in values if math.isfinite(float(v))]
    if len(vals) < 2:
        return float("nan")
    mean = statistics.fmean(vals)
    return float(statistics.pstdev(vals) / max(abs(mean), 1.0e-30))


def run_cluster(cluster: dict) -> dict:
    rho3 = DEC.local_rho3(cluster)
    unit = DEC.unit_decomposition(rho3)
    native = DEC.native_bond_decomposition(rho3)

    traction = float(native["bond_traction_rms"])
    pair_amp = float(unit["pair_amp_rms"])
    m10 = float(unit["m10_vector_rms"])

    coef_fast = float(A8.A8_INIT_DT * A8.A8_INIT_OMEGA * A8.A8_INIT_K)
    coef_slow = float(A8.A8_INIT_DT * A8.A8_INIT_SLOW_TIMESCALE)
    coef_sum = coef_fast + coef_slow

    # These are diagnostics only: what one local native traction increment would
    # look like if multiplied by the already-existing frozen A8 update factors.
    fast_increment = coef_fast * traction
    slow_increment = coef_slow * traction
    combined_increment = coef_sum * traction

    # Inverse-source diagnostics.  They are reported but never used to modify
    # any field or produce a propagation prediction.
    inferred_pair_transfer_coefficient = ratio(pair_amp, traction)
    inferred_m10_transfer_coefficient = ratio(m10, traction)

    return {
        "cluster_id": cluster["id"],
        "native_c_rms": float(native["c_rms"]),
        "native_bond_rms": float(native["bond_diff_rms"]),
        "native_traction_rms": traction,
        "native_traction_over_bond": float(native["traction_over_bond_diff"]),
        "unit_pair_amp_rms": pair_amp,
        "unit_m10_rms": m10,
        "coef_fast": coef_fast,
        "coef_slow": coef_slow,
        "coef_fast_plus_slow": coef_sum,
        "fast_increment_rms": fast_increment,
        "slow_increment_rms": slow_increment,
        "combined_increment_rms": combined_increment,
        "fast_increment_over_pair_amp": ratio(fast_increment, pair_amp),
        "slow_increment_over_pair_amp": ratio(slow_increment, pair_amp),
        "combined_increment_over_pair_amp": ratio(combined_increment, pair_amp),
        "fast_increment_over_m10": ratio(fast_increment, m10),
        "slow_increment_over_m10": ratio(slow_increment, m10),
        "combined_increment_over_m10": ratio(combined_increment, m10),
        "inferred_pair_transfer_coefficient": inferred_pair_transfer_coefficient,
        "inferred_m10_transfer_coefficient": inferred_m10_transfer_coefficient,
        "inferred_pair_over_fast_coef": ratio(inferred_pair_transfer_coefficient, coef_fast),
        "inferred_pair_over_slow_coef": ratio(inferred_pair_transfer_coefficient, coef_slow),
        "inferred_m10_over_fast_coef": ratio(inferred_m10_transfer_coefficient, coef_fast),
        "inferred_m10_over_slow_coef": ratio(inferred_m10_transfer_coefficient, coef_slow),
        "native_c_state_integral_relative_error": float(native["c_state_integral_relative_error"]),
        "native_accumulation_converged": bool(native["accumulation_converged"]),
        "pair_amplitude_exact_formula_relative_rms_error": float(unit["pair_amp_exact_formula_relative_rms_error"]),
    }


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
                failures.append({
                    "cluster_id": cluster["id"],
                    "error": f"{type(exc).__name__}: {exc}",
                })

    inferred_pair = [r["inferred_pair_transfer_coefficient"] for r in rows]
    inferred_m10 = [r["inferred_m10_transfer_coefficient"] for r in rows]
    pair_mean = statistics.fmean(inferred_pair) if inferred_pair else float("nan")
    m10_mean = statistics.fmean(inferred_m10) if inferred_m10 else float("nan")

    checks = {
        "canonical_five_cluster_inventory": ids == EXPECTED_CLUSTER_IDS,
        "all_five_local_fits_present": bool(len(inventory) == 5 and all(x["exists"] for x in inventory)),
        "all_five_clusters_completed": bool(len(rows) == 5 and not failures),
        "all_measured_values_finite": bool(rows and all(
            all(math.isfinite(float(v)) for k, v in r.items()
                if k != "cluster_id" and not isinstance(v, bool))
            for r in rows
        )),
        "historical_pair_formula_reproduced": bool(rows and all(
            r["pair_amplitude_exact_formula_relative_rms_error"] <= 1.0e-14 for r in rows
        )),
        "native_c_state_integral_preserved": bool(rows and all(
            r["native_c_state_integral_relative_error"] <= 1.0e-12 for r in rows
        )),
        "native_accumulation_converged": bool(rows and all(r["native_accumulation_converged"] for r in rows)),
        "no_tracked_or_staged_changes": bool(not state["tracked_changes"] and not state["staged_changes"]),
    }
    execution_gate_pass = bool(all(checks.values()))

    if execution_gate_pass:
        status = "NATIVE_LOCAL_DYNAMIC_RESPONSE_AUDIT_EXECUTED"
    elif rows:
        status = "NATIVE_LOCAL_DYNAMIC_RESPONSE_AUDIT_PARTIAL_EXECUTION"
    else:
        status = "NATIVE_LOCAL_DYNAMIC_RESPONSE_AUDIT_NOT_ESTABLISHED"

    summary = {
        "inferred_pair_transfer_coefficient_mean": pair_mean,
        "inferred_pair_transfer_coefficient_cv": cv(inferred_pair),
        "inferred_m10_transfer_coefficient_mean": m10_mean,
        "inferred_m10_transfer_coefficient_cv": cv(inferred_m10),
        "frozen_fast_coefficient": float(A8.A8_INIT_DT * A8.A8_INIT_OMEGA * A8.A8_INIT_K),
        "frozen_slow_coefficient": float(A8.A8_INIT_DT * A8.A8_INIT_SLOW_TIMESCALE),
    }

    result = {
        "lab_id": LAB_ID,
        "status": status,
        "repo_state": state,
        "model": {
            "question": "Does the static-native-vs-local-M10 mismatch reduce to a local incremental-transfer scale rather than an over-strong equilibrium deformation?",
            "native_static_quantity": "bounded-strain positive-N6 bond traction",
            "historical_local_quantity": "A8 pair transfer then PM1/PS2/M10 interface",
            "inferred_coefficients_role": "diagnostic inverse-source ratios only; not fitted or applied",
            "replacement_strength_scalar": None,
            "normalization_or_rescaling": False,
            "observed_lensing_values_used": False,
            "native_candidate_fed_to_G3D": False,
            "fit_or_tuning": False,
        },
        "summary": summary,
        "rows": rows,
        "failures": failures,
        "checks": checks,
        "execution_gate_pass": execution_gate_pass,
        "interpretation_rule": "Compare source-stability of the inferred local transfer coefficient and the no-fit projections using already-frozen A8 update factors. Do not promote any coefficient to a native law in this audit.",
    }

    print(LAB_ID)
    print(f"status={status}")
    print(f"head_sha={state['head_sha']}")
    print("benchmark_loader=pbuf.core.benchmark_data")
    print("network_access_used=false")
    print("observed_lensing_values_used=false")
    print("replacement_strength_scalar=none")
    print("native_response_rescaled=false")
    print("native_candidate_fed_to_G3D=false")
    print("fit_or_tuning=false")
    print(f"coef_fast={summary['frozen_fast_coefficient']:.12g}")
    print(f"coef_slow={summary['frozen_slow_coefficient']:.12g}")
    print()
    print("CLUSTERS")
    for r in rows:
        print(
            f"cluster={r['cluster_id']} "
            f"traction_over_c={r['native_traction_rms']/max(r['native_c_rms'],1e-30):.12g} "
            f"fast_over_m10={r['fast_increment_over_m10']:.12g} "
            f"slow_over_m10={r['slow_increment_over_m10']:.12g} "
            f"combined_over_m10={r['combined_increment_over_m10']:.12g} "
            f"inferred_pair_coef={r['inferred_pair_transfer_coefficient']:.12g} "
            f"inferred_m10_coef={r['inferred_m10_transfer_coefficient']:.12g}"
        )
    for failure in failures:
        print(f"failure_cluster={failure['cluster_id']} error={failure['error']}")
    print()
    print("SOURCE_STABILITY")
    print(f"inferred_pair_coef_mean={pair_mean:.12g}")
    print(f"inferred_pair_coef_cv={summary['inferred_pair_transfer_coefficient_cv']:.12g}")
    print(f"inferred_m10_coef_mean={m10_mean:.12g}")
    print(f"inferred_m10_coef_cv={summary['inferred_m10_transfer_coefficient_cv']:.12g}")
    print("CHECKS")
    for key, value in checks.items():
        print(f"{key}={str(value).lower()}")
    print(f"execution_gate_pass={str(execution_gate_pass).lower()}")
    print(json.dumps(result, sort_keys=True, separators=(",", ":"), default=str))
    return 0 if execution_gate_pass else 1


if __name__ == "__main__":
    raise SystemExit(main())
