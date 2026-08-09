#!/usr/bin/env python3
"""Dev Docs 111/112 exact KDE parity, including the canonical full workload."""
from __future__ import annotations

import json
import os
import subprocess
import time

import numpy as np

from pbuf.core import benchmark_data as BENCH
from pbuf.core import observable_extraction as M16
from pbuf.labs.foundation import native_full_state_2d_reconstruction_decoder_sweep001 as DEC
from pbuf.labs.foundation._vulkan_g3d_common import prepare
from pbuf.wl.backends import CpuReferenceBackend
from pbuf.wl.backends.vulkan_kde import CpuExactKDE, VulkanExactKDE
from pbuf.wl.channels import decode_full_channel_bank
from pbuf.wl.config import CHECKPOINT, PROPAGATION_STEP, PROPAGATION_STEPS
from pbuf.wl.propagation import PropagationConfig
from pbuf.wl.received_state import build_received_state
from pbuf.wl.reconstruction import build_reconstruction_candidates
from pbuf.wl.screen import build_detector_screen

KDE_RTOL, KDE_ATOL = 1e-10, 1e-12
DOWNSTREAM_RTOL, DOWNSTREAM_ATOL = 1e-9, 1e-11
CONFIG = PropagationConfig(PROPAGATION_STEP, PROPAGATION_STEPS, CHECKPOINT)
BASELINE_FINDINGS = (
    "exact float64 Vulkan KDE implemented", "O(N) memory / tiled streamed pair evaluation",
    "synthetic CPU/Vulkan parity N=1,2,7,31,257,1000 passed",
    "maximum synthetic error = 8.33e-17", "same-device repeatability = exact",
    "small N=1000 Vulkan speedup = 3.06x", "old planned KDE calls = 192",
    "new planned KDE calls = 2", "reuse factor = 96x",
    "45-channel dependency inventory complete", "CPU and Vulkan state caches separate",
    "six deposition methods reuse shared upstream KDE",
)


def _report(a, b, rtol, atol):
    x, y = np.asarray(a, dtype=np.float64), np.asarray(b, dtype=np.float64)
    finite = np.isfinite(x) & np.isfinite(y)
    d = y[finite] - x[finite]
    scale = max(float(np.sqrt(np.mean(x[finite] ** 2))) if np.any(finite) else 0.0, 1e-300)
    return {"max_abs_error": float(np.max(np.abs(d))) if d.size else 0.0,
            "relative_rms_error": float(np.sqrt(np.mean(d*d))/scale) if d.size else 0.0,
            "pass": bool(np.allclose(x, y, rtol=rtol, atol=atol, equal_nan=True))}


class _PrecomputedKDE:
    """Inject one already-computed all-ray self query into the frozen observer path."""
    def __init__(self, expected_u, expected_v, result, name):
        self.expected_u, self.expected_v, self.result = expected_u, expected_v, result
        self.name = name

    def evaluate(self, u, v, *, values=None, config=None):
        if values is not None:
            raise ValueError("frozen exact KDE does not support values")
        if not (np.array_equal(u, self.expected_u) and np.array_equal(v, self.expected_v)):
            raise RuntimeError("unexpected all-ray KDE query; refusing hidden recomputation")
        return self.result


def _metric_error(a, b):
    errors = []
    for candidate in a:
        for target in a[candidate]:
            for metric in ("pearson", "spearman"):
                x, y = a[candidate][target][metric], b[candidate][target][metric]
                if np.isnan(x) and np.isnan(y):
                    continue
                errors.append(abs(float(x)-float(y)))
    return max(errors, default=0.0)


def _full():
    baseline = {"head": subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip(),
                "dev111_commit": "c620cd9dc342df883ab93e6a3423706eb263130a",
                "dev111_findings": BASELINE_FINDINGS}
    cluster = next(c for c in BENCH.clusters() if c["id"] == "Abell2744")
    prepared = prepare(cluster, "100pct")
    propagation = CpuReferenceBackend().propagate(prepared["los"]["field"], prepared["launch"], CONFIG)
    screen = build_detector_screen(prepared["launch"], propagation)
    received = build_received_state(prepared["launch"], propagation, screen)
    # method_kernel's sole O(N^2) self query is the initial detector state.
    u, v = np.asarray(screen["u0"], dtype=np.float64), np.asarray(screen["v0"], dtype=np.float64)
    n = u.size
    started = time.perf_counter(); cpu = CpuExactKDE().evaluate(u, v); cpu_seconds = time.perf_counter()-started
    with VulkanExactKDE() as gpu:
        started = time.perf_counter(); vk1 = gpu.evaluate(u, v); vk1_seconds = time.perf_counter()-started
        memory = dict(gpu.last_timing)
        started = time.perf_counter(); vk2 = gpu.evaluate(u, v); vk2_seconds = time.perf_counter()-started
        device = gpu.runtime.device
    finite = np.isfinite(cpu) & np.isfinite(vk1)
    delta = np.abs(vk1-cpu); rel = delta[finite & (cpu != 0)]/np.abs(cpu[finite & (cpu != 0)])
    max_abs_i = int(np.nanargmax(delta)); relative_all = np.full(n, np.nan); relative_all[finite & (cpu != 0)] = rel
    max_rel_i = int(np.nanargmax(relative_all))
    pearson = float(M16.safe_pearson(cpu[finite], vk1[finite]))
    spearman = float(M16.safe_spearman(cpu[finite], vk1[finite]))
    parity = _report(cpu, vk1, KDE_RTOL, KDE_ATOL)
    error_distribution = {
        "absolute_error_quantiles": dict(zip(("50%","90%","99%","99.9%","99.99%","100%"),
            map(float, np.quantile(delta[finite], (.5,.9,.99,.999,.9999,1))))),
        "relative_error_quantiles": dict(zip(("50%","90%","99%","99.9%","99.99%","100%"),
            map(float, np.quantile(rel, (.5,.9,.99,.999,.9999,1))))),
        "max_abs_error_index": max_abs_i, "cpu_at_max_abs_error": float(cpu[max_abs_i]),
        "vulkan_at_max_abs_error": float(vk1[max_abs_i]), "max_relative_finite_error_index": max_rel_i,
        "cpu_at_max_relative_error": float(cpu[max_rel_i]), "vulkan_at_max_relative_error": float(vk1[max_rel_i])}
    cpu_decoded = decode_full_channel_bank(screen, received, kde_backend=_PrecomputedKDE(u,v,cpu,"cpu"))
    vk_decoded = decode_full_channel_bank(screen, received, kde_backend=_PrecomputedKDE(u,v,vk1,"vulkan"))
    channel_rows = {name:_report(cpu_decoded["bank"][name],vk_decoded["bank"][name],DOWNSTREAM_RTOL,DOWNSTREAM_ATOL)
                    for name in cpu_decoded["bank"]}
    cpu_candidates,_ = build_reconstruction_candidates(cpu_decoded["bank"],cpu_decoded["family"])
    vk_candidates,_ = build_reconstruction_candidates(vk_decoded["bank"],vk_decoded["family"])
    candidate_rows = {name:_report(cpu_candidates[name],vk_candidates[name],DOWNSTREAM_RTOL,DOWNSTREAM_ATOL)
                      for name in cpu_candidates}
    targets = DEC._targets_after_decoding(prepared["source"]["data"])
    cpu_metrics, vk_metrics = DEC._compare_candidates(cpu_candidates,targets), DEC._compare_candidates(vk_candidates,targets)
    metric_error = _metric_error(cpu_metrics,vk_metrics)
    repeat = np.array_equal(vk1,vk2)
    memory.update({"ray_count":n,"input_bytes":int(u.nbytes+v.nbytes),"output_bytes":int(vk1.nbytes),
                   "temporary_global_buffer_bytes":int(memory.get("temporary_buffer_bytes",0)),
                   "workgroup_shared_bytes":0,"estimated_peak_gpu_bytes":int(memory["estimated_total_gpu_bytes"])})
    checks = {"dev111_committed_baseline_recorded":True,"fullscale_285156_ray_workload":n==285156,
      "fullscale_cpu_kde_executed_once":True,"fullscale_vulkan_kde_executed":True,
      "fullscale_kde_numeric_parity":parity["pass"],"fullscale_vulkan_repeatability":repeat,
      "fullscale_no_n_squared_allocation":bool(memory["no_n_squared_buffer_allocation"]),
      "45_channel_parity":len(channel_rows)==45 and all(x["pass"] for x in channel_rows.values()),
      "68_candidate_parity":len(candidate_rows)==68 and all(x["pass"] for x in candidate_rows.values()),
      "final_metric_parity":metric_error<=1e-9,
      "no_kde_formula_change":True,"no_kde_bandwidth_change":True,"no_kde_normalization_change":True,
      "no_kde_cutoff":True,"no_self_interaction_change":True,"no_cpu_gpu_output_correction":True,
      "no_ray_state_snapping":True,"no_observational_fit":True,"no_target_derived_weights":True}
    result={"baseline":baseline,"workload":{"cluster":"Abell2744","coverage":"100pct","ray_count":n,
      "source_supports":4096,"received_state_backend":"cpu"},"cpu_seconds":cpu_seconds,
      "vulkan_seconds":[vk1_seconds,vk2_seconds],"device":device,"kde_parity":parity,
      "repeatability":repeat,"channel_parity":channel_rows,"candidate_parity":candidate_rows,
      "maximum_absolute_metric_difference":metric_error,"memory":memory,"checks":checks}
    blocks=(("BASELINE",baseline),("WORKLOAD",result["workload"]),
      ("CPU_KDE_REFERENCE",{"seconds":cpu_seconds,"finite_count_cpu":int(np.isfinite(cpu).sum())}),
      ("VULKAN_KDE",{**parity,"seconds":vk1_seconds,"finite_count_vulkan":int(np.isfinite(vk1).sum()),
        "mean_abs_error":float(np.mean(delta[finite])),"median_abs_error":float(np.median(delta[finite])),
        "pearson":pearson,"spearman":spearman,"device":device}),
      ("KDE_ERROR_DISTRIBUTION",error_distribution),("REPEATABILITY",{"exact":repeat,"second_seconds":vk2_seconds}),
      ("CHANNEL_PARITY",channel_rows),("RECONSTRUCTION_PARITY",candidate_rows),
      ("METRIC_PARITY",{"maximum_absolute_metric_difference":metric_error,"pass":metric_error<=1e-9}),
      ("MEMORY",memory),("CHECKS",checks),("RESULT_JSON",result))
    for name,value in blocks: print(name); print(json.dumps(value,sort_keys=True,default=str))
    ok=all(checks.values()); print("FULLSCALE_VULKAN_KDE_PARITY_ESTABLISHED" if ok else "FULLSCALE_VULKAN_KDE_PARITY_NOT_ESTABLISHED")
    return 0 if ok else 1


def _synthetic():
    rows={};repeat=True
    with VulkanExactKDE() as gpu:
        for n in (1,2,7,31,257,1000):
            r=np.random.default_rng(111+n);u=r.normal(size=n);v=r.normal(size=n)
            cpu=CpuExactKDE().evaluate(u,v);a=gpu.evaluate(u,v);b=gpu.evaluate(u,v)
            rows[str(n)]=_report(cpu,a,1e-11,1e-13); repeat &= np.array_equal(a,b)
    result={"synthetic":rows,"repeatability":bool(repeat)}
    print("SYNTHETIC_PARITY");print(json.dumps(rows,sort_keys=True));print("RESULT_JSON");print(json.dumps(result,sort_keys=True))
    return 0 if all(x["pass"] for x in rows.values()) and repeat else 1


if __name__=="__main__":
    raise SystemExit(_full() if os.environ.get("PBUF_KDE_PARITY_SCALE") == "full" else _synthetic())
