#!/usr/bin/env python3
"""PBUF FOUNDATION — NATIVE TRANSFER FACTOR DECOMPOSITION AUDIT 001.

Purpose
-------
Test whether the source-stable inverse transfer scale found in the previous audit
can be explained by already-frozen local factors, without fitting or applying a
new coefficient.

The target quantities are diagnostics only:
  pair target = historical pair-amplitude RMS / native bond-traction RMS
  M10 target  = historical M10 RMS / native bond-traction RMS

Candidate factors are fixed before comparison and are built only from quantities
already present in the implementation: A8 fast/slow update coefficients, the two
0.5 fast/slow coupling coefficients, the 0.5 M10 midpoint share, and the N6
one-sixth neighbour share. No candidate is fed into G3D or promoted to a law.
"""
from __future__ import annotations

import json
import math
import statistics
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))

from pbuf.core import benchmark_data as BENCH
from pbuf.models import a8_state as A8
import pbuf.labs.foundation.native_local_dynamic_response_audit001 as PREV

LAB_ID = "PBUF-FOUNDATION-NATIVE-TRANSFER-FACTOR-DECOMPOSITION-AUDIT-001"
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


def mean(xs):
    return statistics.fmean(float(x) for x in xs) if xs else float("nan")


def cv(xs):
    vals = [float(x) for x in xs]
    m = mean(vals)
    return statistics.pstdev(vals) / max(abs(m), 1e-30) if len(vals) >= 2 else float("nan")


def relerr(value, target):
    return abs(float(value) - float(target)) / max(abs(float(target)), 1e-30)


def frozen_factors() -> dict[str, float]:
    fast = float(A8.A8_INIT_DT * A8.A8_INIT_OMEGA * A8.A8_INIT_K)
    slow = float(A8.A8_INIT_DT * A8.A8_INIT_SLOW_TIMESCALE)
    f2s = float(A8.A8_INIT_COUP_F2S)
    s2f = float(A8.A8_INIT_COUP_S2F)
    midpoint = 0.5  # exact M10 contract: +R_ij/2 to each adjacent voxel
    n6_share = 1.0 / 6.0
    return {
        "fast": fast,
        "slow": slow,
        "fast_plus_slow": fast + slow,
        "fast_times_f2s": fast * f2s,
        "fast_times_s2f": fast * s2f,
        "fast_times_midpoint": fast * midpoint,
        "fast_times_f2s_times_midpoint": fast * f2s * midpoint,
        "fast_plus_slow_times_midpoint": (fast + slow) * midpoint,
        "fast_plus_slow_times_f2s_times_midpoint": (fast + slow) * f2s * midpoint,
        "fast_over_N6": fast * n6_share,
        "fast_plus_slow_over_N6": (fast + slow) * n6_share,
        "slow_plus_fast_times_midpoint": slow + fast * midpoint,
        "slow_plus_fast_times_f2s_times_midpoint": slow + fast * f2s * midpoint,
    }


def main() -> int:
    state = repo_state()
    clusters = list(BENCH.clusters())
    inventory = BENCH.inventory()
    ids = tuple(c["id"] for c in clusters)

    rows, failures = [], []
    if ids == EXPECTED_CLUSTER_IDS and len(inventory) == 5 and all(x["exists"] for x in inventory):
        for cluster in clusters:
            try:
                rows.append(PREV.run_cluster(cluster))
            except Exception as exc:
                failures.append({"cluster_id": cluster["id"], "error": f"{type(exc).__name__}: {exc}"})

    pair_vals = [r["inferred_pair_transfer_coefficient"] for r in rows]
    m10_vals = [r["inferred_m10_transfer_coefficient"] for r in rows]
    pair_target = mean(pair_vals)
    m10_target = mean(m10_vals)

    factors = frozen_factors()
    comparisons = []
    for name, value in factors.items():
        comparisons.append({
            "candidate": name,
            "value": value,
            "relative_error_to_pair_target": relerr(value, pair_target),
            "relative_error_to_m10_target": relerr(value, m10_target),
        })
    pair_ranked = sorted(comparisons, key=lambda x: x["relative_error_to_pair_target"])
    m10_ranked = sorted(comparisons, key=lambda x: x["relative_error_to_m10_target"])

    checks = {
        "canonical_five_cluster_inventory": ids == EXPECTED_CLUSTER_IDS,
        "all_five_local_fits_present": bool(len(inventory) == 5 and all(x["exists"] for x in inventory)),
        "all_five_clusters_completed": bool(len(rows) == 5 and not failures),
        "all_targets_finite": bool(math.isfinite(pair_target) and math.isfinite(m10_target)),
        "all_frozen_candidates_finite": bool(all(math.isfinite(v) for v in factors.values())),
        "no_tracked_or_staged_changes": bool(not state["tracked_changes"] and not state["staged_changes"]),
    }
    execution_gate_pass = bool(all(checks.values()))
    if execution_gate_pass:
        status = "NATIVE_TRANSFER_FACTOR_DECOMPOSITION_AUDIT_EXECUTED"
    elif rows:
        status = "NATIVE_TRANSFER_FACTOR_DECOMPOSITION_AUDIT_PARTIAL_EXECUTION"
    else:
        status = "NATIVE_TRANSFER_FACTOR_DECOMPOSITION_AUDIT_NOT_ESTABLISHED"

    result = {
        "lab_id": LAB_ID,
        "status": status,
        "repo_state": state,
        "targets": {
            "pair_transfer_mean": pair_target,
            "pair_transfer_cv": cv(pair_vals),
            "m10_transfer_mean": m10_target,
            "m10_transfer_cv": cv(m10_vals),
            "role": "inverse-source diagnostics only; not fitted or applied",
        },
        "frozen_inputs": {
            "A8_INIT_DT": A8.A8_INIT_DT,
            "A8_INIT_OMEGA": A8.A8_INIT_OMEGA,
            "A8_INIT_K": A8.A8_INIT_K,
            "A8_INIT_SLOW_TIMESCALE": A8.A8_INIT_SLOW_TIMESCALE,
            "A8_INIT_COUP_F2S": A8.A8_INIT_COUP_F2S,
            "A8_INIT_COUP_S2F": A8.A8_INIT_COUP_S2F,
            "M10_midpoint_share": 0.5,
            "N6_neighbour_share": 1.0/6.0,
        },
        "candidate_comparisons": comparisons,
        "best_pair_candidate": pair_ranked[0] if pair_ranked else None,
        "best_m10_candidate": m10_ranked[0] if m10_ranked else None,
        "rows": rows,
        "failures": failures,
        "checks": checks,
        "execution_gate_pass": execution_gate_pass,
        "guardrails": {
            "replacement_strength_scalar": None,
            "normalization_or_rescaling": False,
            "fit_or_tuning": False,
            "candidate_applied_to_native_model": False,
            "candidate_fed_to_G3D": False,
            "observed_lensing_values_used": False,
            "degrees_or_360_factor_used": False,
        },
    }

    print(LAB_ID)
    print(f"status={status}")
    print(f"head_sha={state['head_sha']}")
    print("replacement_strength_scalar=none")
    print("native_response_rescaled=false")
    print("candidate_applied_to_native_model=false")
    print("candidate_fed_to_G3D=false")
    print("observed_lensing_values_used=false")
    print("degrees_or_360_factor_used=false")
    print("fit_or_tuning=false")
    print()
    print("TARGETS")
    print(f"pair_transfer_mean={pair_target:.12g}")
    print(f"pair_transfer_cv={cv(pair_vals):.12g}")
    print(f"m10_transfer_mean={m10_target:.12g}")
    print(f"m10_transfer_cv={cv(m10_vals):.12g}")
    print()
    print("FROZEN_CANDIDATES")
    for item in comparisons:
        print(
            f"candidate={item['candidate']} value={item['value']:.12g} "
            f"pair_relerr={item['relative_error_to_pair_target']:.12g} "
            f"m10_relerr={item['relative_error_to_m10_target']:.12g}"
        )
    if pair_ranked:
        b = pair_ranked[0]
        print(f"best_pair_candidate={b['candidate']} value={b['value']:.12g} relerr={b['relative_error_to_pair_target']:.12g}")
    if m10_ranked:
        b = m10_ranked[0]
        print(f"best_m10_candidate={b['candidate']} value={b['value']:.12g} relerr={b['relative_error_to_m10_target']:.12g}")
    print("CHECKS")
    for k, v in checks.items():
        print(f"{k}={str(v).lower()}")
    print(f"execution_gate_pass={str(execution_gate_pass).lower()}")
    print(json.dumps(result, sort_keys=True, separators=(",", ":"), default=str))
    return 0 if execution_gate_pass else 1


if __name__ == "__main__":
    raise SystemExit(main())
