#!/usr/bin/env python3
"""PBUF FOUNDATION — A8 FAST/SLOW EIGENMODE DECOMPOSITION AUDIT 001.

Purpose
-------
Determine whether the historical suppression from local c_state bond scale to
A8 pair-transfer scale is explained by the actual fast/slow channel structure,
rather than by an arbitrary scalar partition.

The frozen historical channels are rewritten exactly as

    c = (u_fast + u_slow)/2
    d = (u_fast - u_slow)/2

so that

    u_fast = c + d
    u_slow = c - d

and the frozen pair law becomes

    A_ij = coef_fast * Delta u_fast + coef_slow * Delta u_slow
         = (coef_fast + coef_slow) * Delta c
         + (coef_fast - coef_slow) * Delta d.

This is an audit of the existing dynamics only.  No native channel split is
constructed, no inferred ratio is applied, and nothing is fed into G3D.
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
import pbuf.labs.foundation.m10_coverage_25pct_science001 as BASE

LAB_ID = "PBUF-FOUNDATION-A8-FAST-SLOW-EIGENMODE-DECOMPOSITION-AUDIT-001"
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


def cv(values) -> float:
    vals = [float(v) for v in values if math.isfinite(float(v))]
    if len(vals) < 2:
        return float("nan")
    mean = statistics.fmean(vals)
    return float(statistics.pstdev(vals) / max(abs(mean), 1.0e-30))


def corr(a, b) -> float:
    x = np.asarray(a, dtype=np.float64).ravel()
    y = np.asarray(b, dtype=np.float64).ravel()
    if x.size < 2 or y.size != x.size:
        return float("nan")
    sx = float(np.std(x)); sy = float(np.std(y))
    if sx <= 1.0e-30 or sy <= 1.0e-30:
        return float("nan")
    return float(np.corrcoef(x, y)[0, 1])


def historical_state(rho3: np.ndarray) -> dict:
    """Reconstruct the exact historical unit-loading A8 terminal state."""
    rng = np.random.RandomState(BASE.SEED)
    eq = np.asarray(rho3, dtype=np.float64)
    noise = A8.A8_INIT_INJECTION_NOISE * rng.randn(*rho3.shape)
    initial = {
        "rho_3d": rho3.copy(),
        "u_slow0": eq.copy(),
        "u_fast0": eq + noise,
    }
    return BASE._evolve(initial)


def run_cluster(cluster: dict) -> dict:
    rho3 = DEC.local_rho3(cluster)
    state = historical_state(rho3)
    us = np.asarray(state["u_slow"], dtype=np.float64)
    uf = np.asarray(state["u_fast"], dtype=np.float64)
    c = np.asarray(state["c_state"], dtype=np.float64)

    common = 0.5 * (uf + us)
    difference = 0.5 * (uf - us)
    common_identity_err = rms(common - c) / max(rms(c), 1.0e-30)

    dc = DEC.concat_positive_bonds(common)
    dd = DEC.concat_positive_bonds(difference)
    df = DEC.concat_positive_bonds(uf)
    ds = DEC.concat_positive_bonds(us)

    coef_fast = float(A8.A8_INIT_DT * A8.A8_INIT_OMEGA * A8.A8_INIT_K)
    coef_slow = float(A8.A8_INIT_DT * A8.A8_INIT_SLOW_TIMESCALE)
    coef_common = coef_fast + coef_slow
    coef_difference = coef_fast - coef_slow

    common_contrib = coef_common * dc
    difference_contrib = coef_difference * dd
    modal_pair = common_contrib + difference_contrib
    direct_pair = coef_fast * df + coef_slow * ds

    exact_modal_relerr = rms(modal_pair - direct_pair) / max(rms(direct_pair), 1.0e-30)

    # Descriptive projection only: how much of the difference-mode bond follows
    # the common-mode bond in the actual frozen terminal state.  This is never
    # used to construct or tune a model.
    denom = float(np.dot(dc, dc))
    beta_d_on_c = float(np.dot(dd, dc) / denom) if denom > 0.0 else float("nan")
    projected_effective_coef = coef_common + coef_difference * beta_d_on_c

    return {
        "cluster_id": cluster["id"],
        "coef_fast": coef_fast,
        "coef_slow": coef_slow,
        "coef_common": coef_common,
        "coef_difference": coef_difference,
        "common_identity_relative_rms_error": common_identity_err,
        "common_bond_rms": rms(dc),
        "difference_bond_rms": rms(dd),
        "difference_over_common_bond_rms": ratio(rms(dd), rms(dc)),
        "difference_common_bond_correlation": corr(dd, dc),
        "difference_on_common_projection_beta": beta_d_on_c,
        "common_contribution_rms": rms(common_contrib),
        "difference_contribution_rms": rms(difference_contrib),
        "difference_common_contribution_correlation": corr(difference_contrib, common_contrib),
        "pair_rms": rms(direct_pair),
        "pair_over_common_bond_rms": ratio(rms(direct_pair), rms(dc)),
        "pair_over_common_only_prediction": ratio(rms(direct_pair), rms(common_contrib)),
        "projected_effective_coefficient": projected_effective_coef,
        "modal_pair_exact_formula_relative_rms_error": exact_modal_relerr,
        "c_state_rms": rms(c),
        "pair_over_c_state_rms": ratio(rms(direct_pair), rms(c)),
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
                failures.append({"cluster_id": cluster["id"], "error": f"{type(exc).__name__}: {exc}"})

    betas = [r["difference_on_common_projection_beta"] for r in rows]
    effs = [r["projected_effective_coefficient"] for r in rows]
    pair_bond = [r["pair_over_common_bond_rms"] for r in rows]

    checks = {
        "canonical_five_cluster_inventory": ids == EXPECTED_CLUSTER_IDS,
        "all_five_local_fits_present": bool(len(inventory) == 5 and all(x["exists"] for x in inventory)),
        "all_five_clusters_completed": bool(len(rows) == 5 and not failures),
        "all_measured_values_finite": bool(rows and all(
            all(math.isfinite(float(v)) for k, v in r.items() if k != "cluster_id") for r in rows
        )),
        "common_mode_equals_c_state": bool(rows and all(r["common_identity_relative_rms_error"] <= 1.0e-14 for r in rows)),
        "exact_modal_pair_formula_reproduced": bool(rows and all(r["modal_pair_exact_formula_relative_rms_error"] <= 1.0e-14 for r in rows)),
        "no_tracked_or_staged_changes": bool(not state["tracked_changes"] and not state["staged_changes"]),
    }
    execution_gate_pass = bool(all(checks.values()))

    if execution_gate_pass:
        status = "A8_FAST_SLOW_EIGENMODE_DECOMPOSITION_AUDIT_EXECUTED"
    elif rows:
        status = "A8_FAST_SLOW_EIGENMODE_DECOMPOSITION_AUDIT_PARTIAL_EXECUTION"
    else:
        status = "A8_FAST_SLOW_EIGENMODE_DECOMPOSITION_AUDIT_NOT_ESTABLISHED"

    summary = {
        "difference_on_common_projection_beta_mean": statistics.fmean(betas) if betas else float("nan"),
        "difference_on_common_projection_beta_cv": cv(betas),
        "projected_effective_coefficient_mean": statistics.fmean(effs) if effs else float("nan"),
        "projected_effective_coefficient_cv": cv(effs),
        "pair_over_common_bond_rms_mean": statistics.fmean(pair_bond) if pair_bond else float("nan"),
        "pair_over_common_bond_rms_cv": cv(pair_bond),
        "coef_common": float(A8.A8_INIT_DT * A8.A8_INIT_OMEGA * A8.A8_INIT_K + A8.A8_INIT_DT * A8.A8_INIT_SLOW_TIMESCALE),
        "coef_difference": float(A8.A8_INIT_DT * A8.A8_INIT_OMEGA * A8.A8_INIT_K - A8.A8_INIT_DT * A8.A8_INIT_SLOW_TIMESCALE),
    }

    result = {
        "lab_id": LAB_ID,
        "status": status,
        "repo_state": state,
        "model": {
            "common_mode": "c=(u_fast+u_slow)/2",
            "difference_mode": "d=(u_fast-u_slow)/2",
            "exact_pair_modal_law": "A=(coef_fast+coef_slow)*Delta_c+(coef_fast-coef_slow)*Delta_d",
            "native_channel_split_constructed": False,
            "inferred_ratio_applied": False,
            "candidate_fed_to_G3D": False,
            "normalization_or_rescaling": False,
            "fit_or_tuning": False,
            "observed_lensing_values_used": False,
        },
        "summary": summary,
        "rows": rows,
        "failures": failures,
        "checks": checks,
        "execution_gate_pass": execution_gate_pass,
        "interpretation_rule": "Determine whether actual historical fast/slow difference-mode structure cancels part of the common-mode pair transfer. Projection coefficients are descriptive diagnostics only and are never applied.",
    }

    print(LAB_ID)
    print(f"status={status}")
    print(f"head_sha={state['head_sha']}")
    print("network_access_used=false")
    print("observed_lensing_values_used=false")
    print("native_channel_split_constructed=false")
    print("inferred_ratio_applied=false")
    print("candidate_fed_to_G3D=false")
    print("fit_or_tuning=false")
    print(f"coef_common={summary['coef_common']:.12g}")
    print(f"coef_difference={summary['coef_difference']:.12g}")
    print()
    print("CLUSTERS")
    for r in rows:
        print(
            f"cluster={r['cluster_id']} "
            f"d_over_c_bond={r['difference_over_common_bond_rms']:.12g} "
            f"corr_d_c={r['difference_common_bond_correlation']:.12g} "
            f"beta_d_on_c={r['difference_on_common_projection_beta']:.12g} "
            f"pair_over_dc={r['pair_over_common_bond_rms']:.12g} "
            f"pair_over_common_only={r['pair_over_common_only_prediction']:.12g} "
            f"projected_eff_coef={r['projected_effective_coefficient']:.12g}"
        )
    for failure in failures:
        print(f"failure_cluster={failure['cluster_id']} error={failure['error']}")
    print()
    print("SUMMARY")
    for key, value in summary.items():
        print(f"{key}={value:.12g}")
    print("CHECKS")
    for key, value in checks.items():
        print(f"{key}={str(value).lower()}")
    print(f"execution_gate_pass={str(execution_gate_pass).lower()}")
    print(json.dumps(result, sort_keys=True, separators=(",", ":"), default=str))
    return 0 if execution_gate_pass else 1


if __name__ == "__main__":
    raise SystemExit(main())
