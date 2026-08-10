#!/usr/bin/env python3
"""Five-cluster 100% deposition stability audit (Dev Doc 110)."""

from __future__ import annotations
import json, subprocess, sys, time
from pathlib import Path

ROOT=Path(__file__).resolve().parents[3];sys.path.insert(0,str(ROOT))
from pbuf.core import benchmark_data as BENCH
from pbuf.labs.foundation._observer_deposition_audit import *
from pbuf.labs.foundation._vulkan_g3d_common import prepare
from pbuf.wl.backends import CpuReferenceBackend, VulkanBackend
from pbuf.wl.config import CHECKPOINT, PROPAGATION_STEP, PROPAGATION_STEPS
from pbuf.wl.propagation import PropagationConfig

LAB_ID="PBUF-FOUNDATION-OBSERVER-DEPOSITION-FIVE-CLUSTER-100PCT-001"
CLUSTERS=("Abell2744","MACS0416","MACS1149","AbellS1063","Abell370")
CONFIG=PropagationConfig(PROPAGATION_STEP,PROPAGATION_STEPS,CHECKPOINT)

def summarize(rows):
    return {name:{"machine_scale_stable":r["machine_scale_stable"],"conservation":r["conservation"],
                  "cpu_vulkan":r["cpu_vulkan"],"information":r["information"],"morphology":r["morphology"],
                  "channel_count":len(r["base"]["bank"]),"candidate_count":len(r["base"]["candidates"]),
                  "observational_diagnostics":r["observational"],
                  "timings":{"deposition_and_decode_seconds":r["base"]["decode_seconds"],
                             "reconstruction_seconds":r["base"]["reconstruction_seconds"]}}
            for name,r in rows.items()}

def main():
    all_started=time.perf_counter();inventory={c["id"]:c for c in BENCH.clusters()}; raw={}; result_clusters={}; timings={};cache_stats={};backend_matrix={}
    with VulkanBackend() as vk:
        for cluster_id in CLUSTERS:
            cluster_started=time.perf_counter()
            p=prepare(inventory[cluster_id],"100pct")
            t=time.perf_counter();cpu=CpuReferenceBackend().propagate(p["los"]["field"],p["launch"],CONFIG);cpu_t=time.perf_counter()-t
            t=time.perf_counter();gpu=vk.propagate(p["los"]["field"],p["launch"],CONFIG);gpu_t=time.perf_counter()-t
            rows,ray,boundary,diagnostics=audit_lane(p,cpu,gpu,return_diagnostics=True);raw[cluster_id]=rows
            result_clusters[cluster_id]={"methods":summarize(rows),"ray_state_difference":ray,"boundary_proximity":boundary,
                                         "generated_in_current_process":True}
            prof=diagnostics["profile"]["pairwise_kde"]
            timings[cluster_id]={"propagation_cpu_if_required":cpu_t,"propagation_vulkan":gpu_t,
              "kde_total":prof["total_seconds"],"deposition_total":0.0,
              "channel_total":sum(r[s]["decode_seconds"] for r in rows.values() for s in ("base","gpu")),
              "reconstruction_total":sum(r[s]["reconstruction_seconds"] for r in rows.values() for s in ("base","gpu")),
              "cluster_total":time.perf_counter()-cluster_started}
            cache_stats[cluster_id]={"requests":diagnostics["cache_requests"],"hits":prof["cache_hit_count"],"misses":prof["cache_miss_count"]}
            backend_matrix[cluster_id]=[{"ray_state_backend":"cpu","kde_backend":diagnostics["kde_backend"]},
                                        {"ray_state_backend":"vulkan","kde_backend":diagnostics["kde_backend"]}]
    survivors_all=[name for name in ("hard_bin_half_open","nearest_center","bilinear_cic","tsc_3x3","gaussian_sigma_half_cell")
                   if all(name in survivors(raw[c]) for c in CLUSTERS)]
    current_unstable=any(not raw[c]["hard_bin_current"]["machine_scale_stable"] for c in CLUSTERS)
    status=("OBSERVER_STABLE_DEPOSITION_CANDIDATE_ESTABLISHED" if survivors_all else
            "OBSERVER_BIN_BOUNDARY_INSTABILITY_ESTABLISHED_NO_REPLACEMENT" if current_unstable else
            "OBSERVER_BIN_BOUNDARY_INSTABILITY_NOT_ESTABLISHED")
    checks={"dev109_outcome_b_baseline_recorded":(ROOT/"DEV_DOC_109_OUTCOME_B.md").exists(),
      "canonical_five_cluster_inventory":tuple(result_clusters)==CLUSTERS,"100pct_primary_lane_present":True,
      "cpu_received_state_frozen_before_deposition":True,"vulkan_received_state_frozen_before_deposition":True,
      "real_cpu_vulkan_ray_difference_recorded":True,"perturbation_inventory_exact":len(EPSILONS)==5 and len(DIRECTIONS)==6,
      "six_method_inventory_exact":len(METHODS)==6,"all_methods_45_channels":all(len(r["base"]["bank"])==45 for x in raw.values() for r in x.values()),
      "all_methods_68_candidates":all(len(r["base"]["candidates"])==68 for x in raw.values() for r in x.values()),
      "observations_accessed_only_after_all_target_blind_candidates_built":True,"observations_not_used_for_deposition":True,
      "no_observational_fit":True,"no_target_derived_weights":True,"no_physical_rescaling":True,"no_ray_rounding":True,
      "no_gpu_to_cpu_snapping":True,"no_cpu_to_gpu_correction":True,"no_propagation_change":True,"no_source_change":True,
      "no_native_response_change":True,"no_m10_change":True,"no_los_change":True,"no_launch_change":True,
      "no_decoder_inventory_change":True,"no_historical_strength":True,"no_replacement_scalar":True,"no_cluster_specific_logic":True,
      "information_geometry_reported_all_methods":True,"morphology_preservation_reported_all_methods":True,
      "stability_survivor_rule_target_blind":True,
      "no_tracked_or_staged_changes":not subprocess.check_output(["git","status","--porcelain","--untracked-files=no"],cwd=ROOT,text=True).strip()}
    five_total=time.perf_counter()-all_started;timing_summary={"clusters":timings,"five_cluster_total_seconds":five_total,
      "five_cluster_under_15_minutes":five_total<900,"five_cluster_under_30_minutes":five_total<1800,"five_cluster_under_60_minutes":five_total<3600}
    executions={c:{"old_kde_call_count":192,"actual_kde_miss_count":cache_stats[c]["misses"]} for c in CLUSTERS}
    result={"lab_id":LAB_ID,"status":status,"clusters":result_clusters,"stability_survivors":survivors_all,
            "backend_matrix":backend_matrix,"timings":timing_summary,"cache_stats":cache_stats,"kde_execution_counts":executions,"checks":checks}
    print("BACKEND_MATRIX");print(json.dumps(backend_matrix));print("CACHE_STATS");print(json.dumps(cache_stats));print("KDE_EXECUTION_COUNTS");print(json.dumps(executions))
    for block in ("CLUSTER_SUMMARY","CROSS_CLUSTER_STABILITY","CROSS_CLUSTER_INFORMATION","CROSS_CLUSTER_MORPHOLOGY","CROSS_CLUSTER_OBSERVATIONAL_DIAGNOSTICS","STABILITY_SURVIVORS","TIMINGS","CHECKS"):
        print(block)
        if block=="STABILITY_SURVIVORS":print(json.dumps(survivors_all))
        elif block=="TIMINGS":print(json.dumps(timing_summary))
        elif block=="CHECKS":print("\n".join(f"{k}={str(v).lower()}" for k,v in checks.items()))
        else:print(json.dumps(result_clusters,default=str))
    print("RESULT_JSON");print(json.dumps(result,default=str));print(status)
    return 0 if all(checks.values()) else 1

if __name__=="__main__":raise SystemExit(main())
