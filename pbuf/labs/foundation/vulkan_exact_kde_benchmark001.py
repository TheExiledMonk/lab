#!/usr/bin/env python3
"""Production-size exact KDE benchmark (Dev Doc 112)."""
from __future__ import annotations
import json, os, time
import numpy as np
from pbuf.core import benchmark_data as BENCH
from pbuf.labs.foundation._vulkan_g3d_common import prepare
from pbuf.wl.backends import CpuReferenceBackend
from pbuf.wl.backends.vulkan_kde import CpuExactKDE, VulkanExactKDE
from pbuf.wl.config import CHECKPOINT, PROPAGATION_STEP, PROPAGATION_STEPS
from pbuf.wl.propagation import PropagationConfig
from pbuf.wl.screen import build_detector_screen

CONFIG=PropagationConfig(PROPAGATION_STEP,PROPAGATION_STEPS,CHECKPOINT)

def main():
    full=os.environ.get("PBUF_KDE_BENCHMARK_SCALE")=="full"
    if full:
        saved=os.environ.get("PBUF_KDE_CPU_SECONDS")
        if saved is None: raise RuntimeError("full benchmark reuses the one parity CPU timing via PBUF_KDE_CPU_SECONDS")
        cpu=float(saved); cluster=next(c for c in BENCH.clusters() if c["id"]=="Abell2744")
        p=prepare(cluster,"100pct"); prop=CpuReferenceBackend().propagate(p["los"]["field"],p["launch"],CONFIG)
        screen=build_detector_screen(p["launch"],prop);u=np.asarray(screen["u0"]);v=np.asarray(screen["v0"])
    else:
        n=int(os.environ.get("PBUF_KDE_BENCHMARK_RAYS","1000"));r=np.random.default_rng(111);u=r.normal(size=n);v=r.normal(size=n)
        t=time.perf_counter();CpuExactKDE().evaluate(u,v);cpu=time.perf_counter()-t
    n=u.size; init=time.perf_counter();gpu=VulkanExactKDE();context_pipeline=time.perf_counter()-init
    try:
        gpu.evaluate(u,v) # discarded warmup
        warm=[]
        for _ in range(5):
            t=time.perf_counter();gpu.evaluate(u,v);warm.append(time.perf_counter()-t)
        memory=dict(gpu.last_timing);device=gpu.runtime.device
    finally:gpu.close()
    median=float(np.median(warm));speed=cpu/median;pairs=n*n
    memory.update({"ray_count":n,"input_bytes":int(u.nbytes+v.nbytes),"output_bytes":int(n*8),
      "temporary_global_buffer_bytes":int(memory.get("temporary_buffer_bytes",0)),"workgroup_shared_bytes":0,
      "estimated_peak_gpu_bytes":int(memory["estimated_total_gpu_bytes"])})
    result={"ray_count":n,"logical_pair_evaluations":pairs,"cpu_timing_samples":1,"vulkan_warm_samples":5,
      "cpu_seconds":cpu,"vulkan_warm_median_seconds":median,"vulkan_warm_min_seconds":min(warm),
      "vulkan_warm_max_seconds":max(warm),"warm_speedup":speed,"device":device,"memory":memory,
      "timing_classification":"TIMING_CLASSIFICATION_UNAVAILABLE"}
    blocks=(("WORKLOAD",{"cluster":"Abell2744" if full else "synthetic","coverage":"100pct" if full else "fast","ray_count":n}),
      ("CPU_TIMING",{"cpu_timing_samples":1,"cpu_seconds":cpu,"source":"full parity verifier" if full else "current process"}),
      ("VULKAN_INITIALIZATION",{"vulkan_context_and_pipeline_init_seconds":context_pipeline}),
      ("VULKAN_TRANSFER",{"available":False}),("VULKAN_DISPATCH",{"dispatch_host_wait_seconds":warm}),
      ("VULKAN_WARM_TIMING",{"vulkan_warm_samples":5,"minimum":min(warm),"median":median,"maximum":max(warm)}),
      ("PAIR_THROUGHPUT",{"logical_pair_evaluations":pairs,"cpu_pair_evaluations_per_second":pairs/cpu,
        "vulkan_pair_evaluations_per_second":pairs/median}),
      ("SPEEDUP",{"warm_speedup":speed,"target_10x_reached":speed>=10,"target_25x_reached":speed>=25,
        "target_50x_reached":speed>=50,"target_100x_reached":speed>=100,
        "single_100pct_kde_under_60s":median<60,"single_100pct_kde_under_30s":median<30,"single_100pct_kde_under_10s":median<10}),
      ("TIMING_CLASSIFICATION",{"value":"TIMING_CLASSIFICATION_UNAVAILABLE","reason":"device dispatch timing unavailable"}),
      ("MEMORY",memory),("CHECKS",{"fullscale_workload":not full or n==285156,"cpu_reference_not_recomputed":full,
       "vulkan_kde_faster_than_cpu":speed>1,"fullscale_no_n_squared_allocation":memory["no_n_squared_buffer_allocation"]}),
      ("RESULT_JSON",result))
    for name,value in blocks:print(name);print(json.dumps(value,sort_keys=True))
    print("FULLSCALE_VULKAN_KDE_ACCELERATION_ESTABLISHED" if speed>1 else "FULLSCALE_VULKAN_KDE_PARITY_ESTABLISHED_ACCELERATION_NOT_ESTABLISHED")
    return 0 if speed>1 else 1
if __name__=="__main__":raise SystemExit(main())
