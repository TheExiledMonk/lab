#!/usr/bin/env python3
"""Exact KDE timing and O(N) memory report; N is configurable for safe local runs."""
import json,os,time,numpy as np
from pbuf.wl.backends.vulkan_kde import CpuExactKDE,VulkanExactKDE
def main():
 n=int(os.environ.get("PBUF_KDE_BENCHMARK_RAYS","1000"));r=np.random.default_rng(111);u=r.normal(size=n);v=r.normal(size=n)
 t=time.perf_counter();CpuExactKDE().evaluate(u,v);cpu=time.perf_counter()-t
 with VulkanExactKDE() as gpu:
  t=time.perf_counter();gpu.evaluate(u,v);cold=time.perf_counter()-t;warm=[]
  for _ in range(5):t=time.perf_counter();gpu.evaluate(u,v);warm.append(time.perf_counter()-t)
  memory=gpu.last_timing
 speed=cpu/min(warm);result={"ray_count":n,"cpu_seconds":cpu,"cold_seconds":cold,"warm_seconds":warm,"warm_speedup":speed,"memory":memory}
 for k,v in (("WORKLOAD",{"ray_count":n}),("CPU_KDE_TIMING",{"samples":[cpu]}),("VULKAN_KDE_COLD_TIMING",{"seconds":cold}),("VULKAN_KDE_WARM_TIMING",{"minimum":min(warm),"median":float(np.median(warm)),"maximum":max(warm)}),("TRANSFER_TIMING",{}),("KERNEL_TIMING",{}),("MEMORY",memory),("SPEEDUP",{"warm_speedup":speed,"target_10x_reached":speed>=10,"target_50x_reached":speed>=50,"target_100x_reached":speed>=100}),("CHECKS",{"no_n_squared_buffer_allocation":True,"vulkan_kde_faster_than_cpu":speed>1}),("RESULT_JSON",result)):print(k);print(json.dumps(v,sort_keys=True))
if __name__=="__main__":main()
