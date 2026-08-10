#!/usr/bin/env python3
"""Five-cluster 100% Vulkan smoke gate run after deep Abell2744 parity."""

from __future__ import annotations
import json, sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[3];sys.path.insert(0,str(ROOT))
from pbuf.core import benchmark_data as BENCH
from pbuf.labs.foundation._vulkan_g3d_common import downstream, prepare
from pbuf.wl.backends import VulkanBackend
from pbuf.wl.config import CHECKPOINT, PROPAGATION_STEP, PROPAGATION_STEPS
from pbuf.wl.propagation import PropagationConfig

EXPECTED=("Abell2744","MACS0416","MACS1149","AbellS1063","Abell370")

def main():
    rows=[];failures=[];config=PropagationConfig(PROPAGATION_STEP,PROPAGATION_STEPS,CHECKPOINT)
    with VulkanBackend() as backend:
        for cluster in BENCH.clusters():
            try:
                p=prepare(cluster,"100pct");prop=backend.propagate(p["los"]["field"],p["launch"],config);decoded=downstream(p,prop)
                rows.append({"cluster_id":cluster["id"],"ray_count":len(p["launch"].x0),"support_bins":len(prop["groups"]),
                             "channel_count":len(decoded["bank"]),"candidate_count":len(decoded["candidates"])})
            except Exception as exc:failures.append({"cluster_id":cluster["id"],"error":f"{type(exc).__name__}: {exc}"})
    checks={"five_cluster_inventory":tuple(x["cluster_id"] for x in rows)==EXPECTED,
            "all_100pct_ray_counts":all(x["ray_count"]==285156 for x in rows),
            "all_4096_source_support":all(x["support_bins"]==4096 for x in rows),
            "all_45_channels":all(x["channel_count"]==45 for x in rows),
            "all_68_candidates":all(x["candidate_count"]==68 for x in rows)}
    passed=all(checks.values()) and not failures
    print("VULKAN_G3D_FIVE_CLUSTER_100PCT_SMOKE")
    for k,v in checks.items():print(f"{k}={str(bool(v)).lower()}")
    print("RESULT_JSON");print(json.dumps({"rows":rows,"checks":checks,"failures":failures}))
    print("VULKAN_G3D_FIVE_CLUSTER_SMOKE_ESTABLISHED" if passed else "VULKAN_G3D_FIVE_CLUSTER_SMOKE_NOT_ESTABLISHED")
    return 0 if passed else 1
if __name__=="__main__":raise SystemExit(main())
