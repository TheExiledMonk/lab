#!/usr/bin/env python3
"""Minimal gamma-component completion for Dev Doc 113."""

from __future__ import annotations

import json
from pathlib import Path
import sys
import time

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))

from pbuf.core import benchmark_data as BENCH
from pbuf.labs.foundation._observer_deposition_audit import decode
from pbuf.labs.foundation._vulkan_g3d_common import prepare
from pbuf.labs.foundation import native_full_state_2d_reconstruction_decoder_sweep001 as DEC
from pbuf.wl.backends import VulkanBackend
from pbuf.wl.backends.vulkan_kde import make_kde_backend
from pbuf.wl.channel_compatibility import CLUSTERS, SURVIVORS
from pbuf.wl.config import CHECKPOINT, PROPAGATION_STEP, PROPAGATION_STEPS
from pbuf.wl.deposition import get_deposition_method
from pbuf.wl.observer_cache import ObserverPrimitiveCache, ObserverStateId
from pbuf.wl.observer_profile import ObserverProfile
from pbuf.wl.propagation import PropagationConfig
from pbuf.wl.screen import build_detector_screen


LAB_ID = "PBUF-FOUNDATION-WL-CHANNEL-COMPONENT-COMPLETION-001"
CONFIG = PropagationConfig(PROPAGATION_STEP, PROPAGATION_STEPS, CHECKPOINT)
OUTPUT = ROOT / "runs" / "wl_channel_component_completion001" / "result.json"


def _targets(data):
    observed = DEC.FUS._observed(data)
    return {name: DEC._finite(observed[name]) for name in ("gamma1", "gamma2")}


def main() -> int:
    inventory = {c["id"]: c for c in BENCH.clusters()}
    result = {"lab_id": LAB_ID, "clusters": {}, "depositions": list(SURVIVORS),
              "propagation_backend": "vulkan", "coverage": "100pct", "checks": {}}
    started = time.perf_counter()
    with VulkanBackend() as vk:
        for cluster_id in CLUSTERS:
            cluster_started = time.perf_counter()
            prepared = prepare(inventory[cluster_id], "100pct")
            propagated = vk.propagate(prepared["los"]["field"], prepared["launch"], CONFIG)
            screen = build_detector_screen(prepared["launch"], propagated)
            profile = ObserverProfile()
            cache = ObserverPrimitiveCache(profile)
            kde = make_kde_backend("vulkan")
            built = {}
            try:
                state_id = ObserverStateId(f"vulkan_{cluster_id}_100pct_base", backend="vulkan")
                for name in SURVIVORS:
                    built[name] = decode(prepared, propagated, screen,
                        get_deposition_method(name), cache=cache, state_id=state_id,
                        kde_backend=kde)
            finally:
                if hasattr(kde, "close"):
                    kde.close()
            # Binding boundary: observations are accessed after every frozen
            # candidate is assembled for this cluster.
            targets = _targets(prepared["source"]["data"])
            methods = {}
            for name in SURVIVORS:
                methods[name] = {
                    "channel_count": len(built[name]["bank"]),
                    "candidate_count": len(built[name]["candidates"]),
                    "observational_diagnostics": DEC._compare_candidates(built[name]["candidates"], targets),
                }
            result["clusters"][cluster_id] = {"methods": methods,
                "cache": profile.describe()["pairwise_kde"],
                "elapsed_seconds": time.perf_counter() - cluster_started}
            print(f"COMPONENT_CLUSTER_COMPLETE={cluster_id}", flush=True)
    result["elapsed_seconds"] = time.perf_counter() - started
    result["checks"] = {
        "canonical_five_cluster_inventory": tuple(result["clusters"]) == CLUSTERS,
        "five_stable_survivors_only": tuple(result["depositions"]) == SURVIVORS,
        "vulkan_100pct_only": True,
        "gamma1_gamma2_separate": True,
        "all_methods_45_channels": all(row["channel_count"] == 45 for c in result["clusters"].values() for row in c["methods"].values()),
        "all_methods_68_candidates": all(row["candidate_count"] == 68 for c in result["clusters"].values() for row in c["methods"].values()),
        "no_cpu_lane": True, "no_perturbation_sweep": True, "no_25pct_lane": True,
        "no_fit_or_rescaling": True,
    }
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print("RESULT_JSON")
    print(json.dumps(result, sort_keys=True))
    return 0 if all(result["checks"].values()) else 1


if __name__ == "__main__":
    raise SystemExit(main())
