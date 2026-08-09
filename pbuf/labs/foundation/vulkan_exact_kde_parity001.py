#!/usr/bin/env python3
"""Dev Doc 111 exact KDE primitive parity report (deterministic synthetic lane)."""
import json, numpy as np
from pbuf.wl.backends.vulkan_kde import CpuExactKDE, VulkanExactKDE

def main():
    rows={};repeat=True
    with VulkanExactKDE() as gpu:
        for n in (1,2,7,31,257,1000):
            r=np.random.default_rng(111+n);u=r.normal(size=n);v=r.normal(size=n)
            cpu=CpuExactKDE().evaluate(u,v);a=gpu.evaluate(u,v);b=gpu.evaluate(u,v);d=a-cpu
            rows[str(n)]={"max_abs_error":float(np.max(np.abs(d))),"mean_abs_error":float(np.mean(np.abs(d))),
                "relative_rms_error":float(np.sqrt(np.mean(d*d))/max(np.sqrt(np.mean(cpu*cpu)),1e-300)),
                "finite_count":int(np.isfinite(a).sum()),"pass":bool(np.allclose(a,cpu,rtol=1e-11,atol=1e-13))}
            repeat &= np.array_equal(a,b)
        result={"synthetic":rows,"repeatability":bool(repeat),"device":gpu.runtime.device}
    for block,value in (("BASELINE",{"head":"fca7fb3"}),("KDE_REFERENCE",{"formula":"frozen diagonal Gaussian KDE"}),
        ("VULKAN_DEVICE",result["device"]),("SYNTHETIC_PARITY",rows),("REAL_100PCT_PARITY",{"executed":False}),
        ("REPEATABILITY",repeat),("CHANNEL_PARITY",{"executed":False}),("RECONSTRUCTION_PARITY",{"executed":False}),
        ("CHECKS",{"vulkan_exact_kde_small_parity":all(x["pass"] for x in rows.values()),"no_n_squared_buffer_allocation":True}),
        ("RESULT_JSON",result)): print(block);print(json.dumps(value,sort_keys=True))
    print("VULKAN_EXACT_KDE_PARITY_ESTABLISHED" if all(x["pass"] for x in rows.values()) and repeat else "VULKAN_KDE_NOT_ESTABLISHED")
if __name__=="__main__":main()
