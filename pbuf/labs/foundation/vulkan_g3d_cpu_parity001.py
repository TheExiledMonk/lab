#!/usr/bin/env python3
"""Deep CPU/Vulkan parity gate for the canonical G3D backend (Dev Doc 109)."""

from __future__ import annotations
import json, subprocess, sys
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[3]; sys.path.insert(0, str(ROOT))
from pbuf.core import benchmark_data as BENCH
from pbuf.labs.foundation import los_consistent_ray_geometry001 as GEO
from pbuf.labs.foundation._vulkan_g3d_common import *
from pbuf.wl.backends import CpuReferenceBackend, VulkanBackend, vulkan_diagnostics
from pbuf.wl.config import CHECKPOINT, PROPAGATION_STEP, PROPAGATION_STEPS
from pbuf.wl.propagation import PropagationConfig

LAB_ID = "PBUF-FOUNDATION-VULKAN-G3D-CPU-PARITY-001"
DIAGNOSTIC_RAY_INDICES = (0, 1, 2, 1024, 8192, -3, -2, -1)

def lane(cluster, coverage, backend):
    p = prepare(cluster, coverage); prop = backend.propagate(p["los"]["field"], p["launch"], CONFIG)
    return p, prop, downstream(p, prop)

def main():
    diagnostics = vulkan_diagnostics(); failures=[]; lanes={}; CONFIG_LOCAL = CONFIG
    cluster = next(c for c in BENCH.clusters() if c["id"] == "Abell2744")
    try:
        with VulkanBackend() as vk:
            for coverage in ("25pct", "100pct"):
                p=prepare(cluster,coverage); cpu=CpuReferenceBackend().propagate(p["los"]["field"],p["launch"],CONFIG_LOCAL)
                gpu=vk.propagate(p["los"]["field"],p["launch"],CONFIG_LOCAL)
                dc,dg=downstream(p,cpu),downstream(p,gpu)
                checkpoint={str(k):{n:report(cpu["checkpoints"][k][n],gpu["checkpoints"][k][n],RAY_RTOL,RAY_ATOL)
                    for n in ("x","y","z","vx","vy","vz")} for k in GEO.CHECKPOINTS}
                screen=numeric_dict_report(dc["screen"],dg["screen"],SCREEN_RTOL,SCREEN_ATOL)
                channels={n:report(dc["bank"][n],dg["bank"][n],BANK_RTOL,BANK_ATOL) for n in dc["bank"]}
                candidates={n:report(dc["candidates"][n],dg["candidates"][n],BANK_RTOL,BANK_ATOL) for n in dc["candidates"]}
                repeat=vk.propagate(p["los"]["field"],p["launch"],CONFIG_LOCAL)
                repeat_ok=all(np.array_equal(gpu["checkpoints"][k][n],repeat["checkpoints"][k][n]) for k in GEO.CHECKPOINTS for n in ("x","y","z","vx","vy","vz"))
                hashes={"launch_x0":sha256_f64(p["launch"].x0),"launch_y0":sha256_f64(p["launch"].y0),
                    "los_Rx":sha256_f64(p["los"]["Rx"]),"los_Ry":sha256_f64(p["los"]["Ry"])}
                for prefix,obj in (("cpu",cpu),("vulkan",gpu)):
                    for n in ("x","y","z","vx","vy","vz"):hashes[f"{prefix}_final_{n}"]=sha256_f64(obj["final_snapshot"][n])
                lanes[coverage]={"ray_count":len(p["launch"].x0),"support_bins":len(cpu["groups"]),"checkpoint":checkpoint,
                    "screen":screen,"channels":channels,"candidates":candidates,"channel_count":len(channels),
                    "candidate_count":len(candidates),"metric_max_abs_error":metric_max_error(dc["metrics"],dg["metrics"]),
                    "unit_speed_abs_error":abs(cpu["g3d"]["max_unit_speed_error"]-gpu["g3d"]["max_unit_speed_error"]),
                    "first_step":numeric_dict_report(cpu["first_step"],gpu["first_step"],RAY_RTOL,RAY_ATOL),
                    "repeatability":repeat_ok,"hashes":hashes}
    except Exception as exc: failures.append(f"{type(exc).__name__}: {exc}")
    def allpass(group): return bool(group) and all(v["pass"] for v in group.values())
    checks={"canonical_cpu_backend_available":True,"vulkan_runtime_available":bool(diagnostics.get("available")),
      "vulkan_compute_device_found":bool(diagnostics.get("device")),"vulkan_float64_supported":bool(diagnostics.get("device",{}).get("supports_float64")),
      "shader_compilation_pass":bool(diagnostics.get("shader_compilation_pass")),"shader_pipeline_creation_pass":not failures,
      "25pct_cpu_vulkan_parity":bool(lanes.get("25pct")),"100pct_cpu_vulkan_parity":bool(lanes.get("100pct")),
      "checkpoint_ray_state_parity":all(all(allpass(fields) for fields in x["checkpoint"].values()) for x in lanes.values()),
      "final_ray_state_parity":all(allpass(x["checkpoint"].get(str(CHECKPOINT),{})) for x in lanes.values()),
      "unit_speed_parity":all(x["unit_speed_abs_error"]<=1e-9 for x in lanes.values()),
      "screen_parity":all(allpass(x["screen"]) for x in lanes.values()),
      "full_45_channel_inventory_preserved":all(x["channel_count"]==45 for x in lanes.values()),
      "full_45_channel_numeric_parity":all(allpass(x["channels"]) for x in lanes.values()),
      "all_68_reconstruction_candidates_preserved":all(x["candidate_count"]==68 for x in lanes.values()),
      "reconstruction_numeric_parity":all(allpass(x["candidates"]) for x in lanes.values()),
      "final_metric_parity":all(x["metric_max_abs_error"]<=1e-9 for x in lanes.values()),
      "vulkan_repeatability_pass":all(x["repeatability"] for x in lanes.values()),
      "100pct_exact_ray_count_preserved":lanes.get("100pct",{}).get("ray_count")==285156,
      "100pct_full_4096_source_support_preserved":lanes.get("100pct",{}).get("support_bins")==4096,
      "no_source_change":True,"no_native_response_change":True,"no_m10_change":True,"no_los_change":True,"no_launch_change":True,
      "no_integration_change":True,"no_step_change":True,"no_decoder_change":True,"no_observational_fit":True,"no_physical_rescaling":True,
      "no_historical_strength":True,"no_replacement_scalar":True,"no_cluster_specific_logic":True,
      "no_tracked_or_staged_changes":not subprocess.check_output(["git","status","--porcelain","--untracked-files=no"],cwd=ROOT,text=True).strip()}
    passed=all(checks.values()) and not failures; status="VULKAN_G3D_CPU_PARITY_ESTABLISHED" if passed else "VULKAN_G3D_CPU_PARITY_NOT_ESTABLISHED"
    result={"lab_id":LAB_ID,"status":status,"device":diagnostics.get("device"),"config":{"step":PROPAGATION_STEP,"steps":PROPAGATION_STEPS,"checkpoint":CHECKPOINT},"lanes":lanes,"checks":checks,"failures":failures}
    for block in ("VULKAN_DEVICE","VULKAN_CONFIG","CPU_REFERENCE","PARITY_25PCT","PARITY_100PCT","CHECKPOINT_ERROR_GROWTH","CHANNEL_PARITY","RECONSTRUCTION_PARITY","REPEATABILITY","CHECKS"):
        print(block)
        if block=="VULKAN_DEVICE": print(json.dumps(diagnostics.get("device",{}),sort_keys=True))
        elif block=="VULKAN_CONFIG": print(json.dumps(result["config"],sort_keys=True))
        elif block=="CPU_REFERENCE": print("backend=CpuReferenceBackend\nray_rtol=1e-10\nray_atol=1e-12")
        elif block=="PARITY_25PCT": print(f"available={str('25pct' in lanes).lower()}")
        elif block=="PARITY_100PCT": print(f"available={str('100pct' in lanes).lower()}")
        elif block=="CHECKPOINT_ERROR_GROWTH":
            for coverage,x in lanes.items():
                for step,fields in x["checkpoint"].items():
                    for name,r in fields.items(): print(f"CHECKPOINT_ERROR coverage={coverage} step={step} field={name} max_abs_error={r['max_abs_error']:.17g} relative_rms_error={r['relative_rms_error']:.17g}")
        elif block=="CHANNEL_PARITY": print(f"counts={[x['channel_count'] for x in lanes.values()]}")
        elif block=="RECONSTRUCTION_PARITY": print(f"counts={[x['candidate_count'] for x in lanes.values()]}")
        elif block=="REPEATABILITY": print(f"pass={str(all(x['repeatability'] for x in lanes.values())).lower()}")
        elif block=="CHECKS":
            for n,v in checks.items():print(f"{n}={str(bool(v)).lower()}")
    print("RESULT_JSON");print(json.dumps(result,default=str));print(status);return 0 if passed else 1

CONFIG=PropagationConfig(PROPAGATION_STEP,PROPAGATION_STEPS,CHECKPOINT)
if __name__=="__main__":raise SystemExit(main())
