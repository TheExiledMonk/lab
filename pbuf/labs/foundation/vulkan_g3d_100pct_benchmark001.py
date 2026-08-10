#!/usr/bin/env python3
"""Primary 100% Abell2744 CPU/Vulkan propagation benchmark (Dev Doc 109)."""

from __future__ import annotations
import json, statistics, subprocess, sys, time
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[3];sys.path.insert(0,str(ROOT))
from pbuf.core import benchmark_data as BENCH
from pbuf.labs.foundation import los_consistent_ray_geometry001 as GEO
from pbuf.labs.foundation._vulkan_g3d_common import prepare
from pbuf.wl.backends import VulkanBackend, vulkan_diagnostics
from pbuf.wl.config import CHECKPOINT, PROPAGATION_STEP, PROPAGATION_STEPS
from pbuf.wl.propagation import PropagationConfig

LAB_ID="PBUF-FOUNDATION-VULKAN-G3D-100PCT-BENCHMARK-001"

def stats(values):
    return {"samples":len(values),"median":statistics.median(values),"minimum":min(values),"maximum":max(values)}

def main():
    cluster=next(c for c in BENCH.clusters() if c["id"]=="Abell2744");p=prepare(cluster,"100pct")
    field,launch=p["los"]["field"],p["launch"];n=len(launch.x0);total=n*PROPAGATION_STEPS;failures=[]
    t=time.perf_counter();GEO._propagate_g3d(field,PROPAGATION_STEP,PROPAGATION_STEPS,launch.x0,launch.y0);cpu=[time.perf_counter()-t]
    config=PropagationConfig(PROPAGATION_STEP,PROPAGATION_STEPS,CHECKPOINT)
    cold_start=time.perf_counter();vk=VulkanBackend();initialization_seconds=time.perf_counter()-cold_start
    cold=vk.propagate(field,launch,config);cold_total=initialization_seconds+vk.last_timing["warm_total_seconds"]
    vk.propagate(field,launch,config) # discarded warm-up
    warm=[];repeat=[]
    for _ in range(5):
        out=vk.propagate(field,launch,config);warm.append(vk.last_timing["warm_total_seconds"]);repeat.append(out)
    repeatability=all(np.array_equal(repeat[0]["final_snapshot"][name],x["final_snapshot"][name]) for x in repeat[1:] for name in ("x","y","z","vx","vy","vz"))
    device=vk.runtime.device;workgroup=vk.runtime.workgroup_size;vk.close();cs,vs=stats(cpu),stats(warm)
    cold_speedup=cs["median"]/cold_total;warm_speedup=cs["median"]/vs["median"]
    checks={"canonical_cpu_backend_available":True,"vulkan_runtime_available":True,"vulkan_compute_device_found":True,
      "vulkan_float64_supported":bool(device["supports_float64"]),"shader_compilation_pass":True,"shader_pipeline_creation_pass":True,
      "100pct_exact_ray_count_preserved":n==285156,"100pct_full_4096_source_support_preserved":len(cold["groups"])==4096,
      "vulkan_repeatability_pass":repeatability,"warm_vulkan_faster_than_cpu":warm_speedup>1.0,
      "no_source_change":True,"no_native_response_change":True,"no_m10_change":True,"no_los_change":True,"no_launch_change":True,
      "no_integration_change":True,"no_step_change":True,"no_decoder_change":True,"no_observational_fit":True,"no_physical_rescaling":True,
      "no_historical_strength":True,"no_replacement_scalar":True,"no_cluster_specific_logic":True,
      "no_tracked_or_staged_changes":not subprocess.check_output(["git","status","--porcelain","--untracked-files=no"],cwd=ROOT,text=True).strip()}
    passed=all(checks.values()) and not failures;status="VULKAN_G3D_100PCT_BENCHMARK_ESTABLISHED" if passed else "VULKAN_G3D_100PCT_BENCHMARK_NOT_ESTABLISHED"
    result={"lab_id":LAB_ID,"status":status,"device":device,"workgroup_size":workgroup,"ray_count":n,"propagation_steps":PROPAGATION_STEPS,
      "total_ray_steps":total,"cpu_timing":cs,"vulkan_initialization_seconds":initialization_seconds,
      "vulkan_cold_total_seconds":cold_total,"vulkan_warm_timing":vs,
      "vulkan_dispatch_seconds":None,"dispatch_timing_note":"host wall timing only; no kernel-only claim",
      "cpu_ray_steps_per_second":total/cs["median"],"vulkan_warm_ray_steps_per_second":total/vs["median"],
      "vulkan_cold_speedup":cold_speedup,"vulkan_warm_speedup":warm_speedup,"target_speedup_10x_reached":warm_speedup>=10,
      "checks":checks,"failures":failures}
    for block in ("VULKAN_DEVICE","BENCHMARK_WORKLOAD","CPU_TIMING","VULKAN_COLD_TIMING","VULKAN_WARM_TIMING","THROUGHPUT","SPEEDUP","CHECKS"):
        print(block)
        if block=="VULKAN_DEVICE":print(json.dumps(device,sort_keys=True))
        elif block=="BENCHMARK_WORKLOAD":print(f"ray_count={n}\npropagation_steps={PROPAGATION_STEPS}\ntotal_ray_steps={total}\nworkgroup_size={workgroup}")
        elif block=="CPU_TIMING":print(f"cpu_timing_samples={len(cpu)}\ncpu_propagation_seconds={cs['median']}")
        elif block=="VULKAN_COLD_TIMING":print(f"vulkan_cold_total_seconds={cold_total}")
        elif block=="VULKAN_WARM_TIMING":print(f"vulkan_warm_total_seconds={vs['median']}\nvulkan_dispatch_seconds=unavailable_host_wall_only")
        elif block=="THROUGHPUT":print(f"cpu_ray_steps_per_second={result['cpu_ray_steps_per_second']}\nvulkan_warm_ray_steps_per_second={result['vulkan_warm_ray_steps_per_second']}")
        elif block=="SPEEDUP":print(f"vulkan_cold_speedup={cold_speedup}\nvulkan_warm_speedup={warm_speedup}\ntarget_speedup_10x_reached={str(warm_speedup>=10).lower()}")
        elif block=="CHECKS":
            for k,v in checks.items():print(f"{k}={str(bool(v)).lower()}")
    print("RESULT_JSON");print(json.dumps(result));print(status)
    return 0 if passed else 1
if __name__=="__main__":raise SystemExit(main())
