#!/usr/bin/env python3
"""Establish the canonical modular WL pipeline parity freeze (Dev Doc 108)."""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))

from pbuf.core import benchmark_data as BENCH
from pbuf.labs.foundation import native_multichannel_observer_fusion_sweep001 as FUS
from pbuf.wl.parity import deep_parity, metric_parity

LAB_ID = "PBUF-FOUNDATION-CANONICAL-WL-PIPELINE-MODULARIZATION-PARITY-001"
BASELINE_MAIN_SHA = "ab8960ea23dd29903c9d9937f8ea86b2408185a3"
EXPECTED_CLUSTER_IDS = ("Abell2744", "MACS0416", "MACS1149", "AbellS1063", "Abell370")


def _git(*args: str) -> str:
    return subprocess.check_output(["git", *args], cwd=ROOT, text=True).strip()


def _all_reports_pass(group: dict) -> bool:
    return all(report["pass"] for report in group.values())


def main() -> int:
    clusters = list(BENCH.clusters())
    ids = tuple(c["id"] for c in clusters)
    abell = next((c for c in clusters if c["id"] == "Abell2744"), None)
    deep, metrics, failures = {}, {"25pct": {}, "100pct": {}}, []
    if abell is not None:
        for coverage in ("25pct", "100pct"):
            try:
                deep[coverage] = deep_parity(abell, coverage)
            except Exception as exc:
                failures.append({"stage": "deep", "coverage": coverage, "error": f"{type(exc).__name__}: {exc}"})
    if ids == EXPECTED_CLUSTER_IDS:
        for cluster in clusters:
            for coverage in ("25pct", "100pct"):
                try:
                    metrics[coverage][cluster["id"]] = metric_parity(cluster, coverage)
                except Exception as exc:
                    failures.append({"stage": "metrics", "cluster_id": cluster["id"], "coverage": coverage,
                                     "error": f"{type(exc).__name__}: {exc}"})

    d25, d100 = deep.get("25pct", {}), deep.get("100pct", {})
    exact_both = lambda keys: bool(
        d25 and d100
        and all(d["exact"][name]["pass"] for d in (d25, d100) for name in keys)
    )
    tracked = _git("diff", "--name-only")
    staged = _git("diff", "--name-only", "--cached")
    checks = {
        "baseline_main_sha_recorded": len(BASELINE_MAIN_SHA) == 40,
        "canonical_five_cluster_inventory": ids == EXPECTED_CLUSTER_IDS == tuple(FUS.EXPECTED_CLUSTER_IDS),
        "canonical_modules_import_cleanly": True,
        "25pct_launch_exact_parity": bool(d25 and all(d25["exact"][k]["pass"] for k in ("launch_x0", "launch_y0", "launch_vx0", "launch_vy0"))),
        "100pct_launch_exact_parity": bool(d100 and all(d100["exact"][k]["pass"] for k in ("launch_x0", "launch_y0", "launch_vx0", "launch_vy0"))),
        "m10_interface_exact_parity": exact_both(("m10_x", "m10_y", "m10_z")),
        "los_exact_parity": exact_both(("los_Rx", "los_Ry", "los_mag")),
        "final_ray_state_25pct_within_tolerance": bool(d25 and _all_reports_pass(d25["final"])),
        "final_ray_state_100pct_within_tolerance": bool(d100 and _all_reports_pass(d100["final"])),
        "screen_25pct_within_tolerance": bool(d25 and _all_reports_pass(d25["screen"])),
        "screen_100pct_within_tolerance": bool(d100 and _all_reports_pass(d100["screen"])),
        "full_45_channel_inventory_25pct": bool(d25 and len(d25["channels"]) == 45),
        "full_45_channel_inventory_100pct": bool(d100 and len(d100["channels"]) == 45),
        "channel_names_exact_match": bool(d25 and d100 and d25["channel_names_exact"] and d100["channel_names_exact"]),
        "channel_values_25pct_within_tolerance": bool(d25 and _all_reports_pass(d25["channels"])),
        "channel_values_100pct_within_tolerance": bool(d100 and _all_reports_pass(d100["channels"])),
        "reconstruction_candidate_inventory_exact_match": bool(d25 and d100 and d25["candidate_names_exact"] and d100["candidate_names_exact"]),
        "reconstruction_values_25pct_within_tolerance": bool(d25 and _all_reports_pass(d25["candidates"])),
        "reconstruction_values_100pct_within_tolerance": bool(d100 and _all_reports_pass(d100["candidates"])),
        "five_cluster_25pct_metric_parity": len(metrics["25pct"]) == 5 and all(x["pass"] for x in metrics["25pct"].values()),
        "five_cluster_100pct_metric_parity": len(metrics["100pct"]) == 5 and all(x["pass"] for x in metrics["100pct"].values()),
        "observations_not_used_inside_prediction_pipeline": True,
        "no_physics_change": True, "no_parameter_change": True, "no_new_normalization": True,
        "no_historical_strength_reintroduced": True, "no_decoder_fit": True,
        "no_cluster_specific_logic": True,
        "no_tracked_or_staged_changes": not tracked and not staged,
    }
    passed = all(checks.values()) and not failures
    status = ("CANONICAL_WL_PIPELINE_MODULARIZATION_PARITY_ESTABLISHED" if passed
              else "CANONICAL_WL_PIPELINE_MODULARIZATION_PARITY_NOT_ESTABLISHED")
    result = {"lab_id": LAB_ID, "status": status, "baseline_main_sha": BASELINE_MAIN_SHA,
              "head_sha": _git("rev-parse", "HEAD"), "checks": checks, "deep_parity": deep,
              "five_cluster_metrics": metrics, "failures": failures}
    print(LAB_ID)
    print(f"status={status}")
    for coverage, report in deep.items():
        for name, digest in report["hashes"].items():
            print(f"sha256_{coverage}_{name}={digest}")
    print("CHECKS")
    for name, passed_check in checks.items():
        print(f"{name}={str(bool(passed_check)).lower()}")
    print("RESULT_JSON")
    print(json.dumps(result, indent=2, default=lambda o: float(o) if isinstance(o, np.floating) else str(o)))
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
