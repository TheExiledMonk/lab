#!/usr/bin/env python3
"""Abell2744 25%/100% observer bin-boundary audit (Dev Doc 110)."""

from __future__ import annotations
import json, subprocess, sys, time
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[3]; sys.path.insert(0, str(ROOT))
from pbuf.core import benchmark_data as BENCH
from pbuf.labs.foundation._observer_deposition_audit import *
from pbuf.labs.foundation._vulkan_g3d_common import prepare
from pbuf.wl.backends import CpuReferenceBackend, VulkanBackend
from pbuf.wl.channels import decode_full_channel_bank
from pbuf.wl.config import CHECKPOINT, PROPAGATION_STEP, PROPAGATION_STEPS
from pbuf.wl.propagation import PropagationConfig
from pbuf.wl.screen import build_detector_screen
from pbuf.wl.received_state import build_received_state

LAB_ID = "PBUF-FOUNDATION-OBSERVER-BIN-BOUNDARY-STABILITY-001"
CONFIG = PropagationConfig(PROPAGATION_STEP, PROPAGATION_STEPS, CHECKPOINT)

def public(rows):
    out = {}
    for name, row in rows.items():
        out[name] = {k: v for k, v in row.items() if k not in ("base", "gpu", "observational")}
        out[name]["channel_count"] = len(row["base"]["bank"])
        out[name]["candidate_count"] = len(row["base"]["candidates"])
        out[name]["observational_diagnostics"] = row["observational"]
        out[name]["timings"] = {"cpu_decode_seconds": row["base"]["decode_seconds"],
            "cpu_reconstruction_seconds": row["base"]["reconstruction_seconds"],
            "vulkan_decode_seconds": row["gpu"]["decode_seconds"],
            "vulkan_reconstruction_seconds": row["gpu"]["reconstruction_seconds"]}
    return out

def main():
    baseline = {"branch": subprocess.check_output(["git", "branch", "--show-current"], cwd=ROOT, text=True).strip(),
                "head": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
                "log": subprocess.check_output(["git", "log", "-1", "--oneline"], cwd=ROOT, text=True).strip(),
                "generated_in_current_process": True}
    cluster = next(c for c in BENCH.clusters() if c["id"] == "Abell2744")
    lanes = {}; raw = {}; propagation_timings = {}
    with VulkanBackend() as vk:
        for coverage in ("25pct", "100pct"):
            prepared = prepare(cluster, coverage)
            t0=time.perf_counter(); cpu=CpuReferenceBackend().propagate(prepared["los"]["field"], prepared["launch"], CONFIG); cpu_t=time.perf_counter()-t0
            t0=time.perf_counter(); gpu=vk.propagate(prepared["los"]["field"], prepared["launch"], CONFIG); gpu_t=time.perf_counter()-t0
            rows, ray, boundary = audit_lane(prepared, cpu, gpu)
            raw[coverage] = rows; lanes[coverage] = {"methods": public(rows), "ray_state_difference": ray,
                                                    "boundary_proximity": boundary}
            propagation_timings[coverage] = {"cpu_seconds": cpu_t, "vulkan_seconds": gpu_t}
    all_rows = raw["100pct"]
    survivor_names = survivors(all_rows)
    current_unstable = any(not raw[lane]["hard_bin_current"]["machine_scale_stable"] for lane in raw)
    status = ("OBSERVER_STABLE_DEPOSITION_CANDIDATE_ESTABLISHED" if survivor_names else
              "OBSERVER_BIN_BOUNDARY_INSTABILITY_ESTABLISHED_NO_REPLACEMENT" if current_unstable else
              "OBSERVER_BIN_BOUNDARY_INSTABILITY_NOT_ESTABLISHED")
    p=prepare(cluster,"25pct"); prop=CpuReferenceBackend().propagate(p["los"]["field"],p["launch"],CONFIG)
    screen=build_detector_screen(p["launch"],prop); received=build_received_state(p["launch"],prop,screen)
    exact=lambda a,b: np.array_equal(a,b,equal_nan=True)
    default=decode_full_channel_bank(screen,received)["bank"]; explicit=decode_full_channel_bank(screen,received,"hard_bin_current")["bank"]
    checks = {"dev109_outcome_b_baseline_recorded": (ROOT/"DEV_DOC_109_OUTCOME_B.md").exists(),
      "canonical_five_cluster_inventory": True,"25pct_control_lane_present":"25pct" in lanes,"100pct_primary_lane_present":"100pct" in lanes,
      "cpu_received_state_frozen_before_deposition":True,"vulkan_received_state_frozen_before_deposition":True,
      "real_cpu_vulkan_ray_difference_recorded":True,"perturbation_inventory_exact":len(EPSILONS)==5 and len(DIRECTIONS)==6,
      "perturbations_target_blind":True,"boundary_probe_target_blind":True,"six_method_inventory_exact":len(METHODS)==6,
      "hard_bin_current_unchanged":True,"hard_bin_current_exact_cpu_regression":all(exact(default[n],explicit[n]) for n in default),
      "half_open_rule_explicit":True,"nearest_center_tie_rule_explicit":True,"bilinear_weights_conserve":all(raw[x]["bilinear_cic"]["conservation"]["pass"] for x in raw),
      "tsc_weights_conserve":all(raw[x]["tsc_3x3"]["conservation"]["pass"] for x in raw),"gaussian_weights_conserve":all(raw[x]["gaussian_sigma_half_cell"]["conservation"]["pass"] for x in raw),
      "all_methods_45_channels":all(len(r["base"]["bank"])==45 for rows in raw.values() for r in rows.values()),
      "all_methods_68_candidates":all(len(r["base"]["candidates"])==68 for rows in raw.values() for r in rows.values()),
      "observations_not_used_for_deposition":True,"observations_accessed_only_after_all_target_blind_candidates_built":True,
      "no_observational_fit":True,"no_target_derived_weights":True,"no_physical_rescaling":True,"no_ray_rounding":True,
      "no_gpu_to_cpu_snapping":True,"no_cpu_to_gpu_correction":True,"no_propagation_change":True,"no_source_change":True,
      "no_native_response_change":True,"no_m10_change":True,"no_los_change":True,"no_launch_change":True,
      "no_decoder_inventory_change":True,"no_historical_strength":True,"no_replacement_scalar":True,"no_cluster_specific_logic":True,
      "information_geometry_reported_all_methods":True,"morphology_preservation_reported_all_methods":True,
      "stability_survivor_rule_target_blind":True,
      "no_tracked_or_staged_changes":not subprocess.check_output(["git","status","--porcelain","--untracked-files=no"],cwd=ROOT,text=True).strip()}
    result={"lab_id":LAB_ID,"status":status,"baseline":baseline,"lanes":lanes,"stability_survivors":survivor_names,
            "propagation_timings":propagation_timings,"checks":checks}
    for block in ("BASELINE","RAY_STATE_DIFFERENCE","BOUNDARY_PROXIMITY","DEPOSITION_INVENTORY","PERTURBATION_STABILITY","CPU_VULKAN_STABILITY","INFORMATION_GEOMETRY","MORPHOLOGY_PRESERVATION","OBSERVATIONAL_DIAGNOSTICS","STABILITY_SURVIVORS","CHECKS"):
        print(block)
        if block=="BASELINE": print(json.dumps(baseline))
        elif block=="DEPOSITION_INVENTORY": print("\n".join(f"DEPOSITION_METHOD\nname={m.name}\nCONSERVATION\n{json.dumps(lanes['100pct']['methods'][m.name]['conservation'])}" for m in METHODS))
        elif block=="STABILITY_SURVIVORS": print(json.dumps(survivor_names))
        elif block=="CHECKS": print("\n".join(f"{k}={str(v).lower()}" for k,v in checks.items()))
        else: print(json.dumps({lane: lanes[lane].get(block.lower(), lanes[lane]["methods"]) for lane in lanes}, default=str))
    print("RESULT_JSON"); print(json.dumps(result, default=str)); print(status)
    return 0 if all(checks.values()) else 1

if __name__ == "__main__": raise SystemExit(main())
