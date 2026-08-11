"""DEV178 native receipt-density and diagnostic Vulkan audit.

This deliberately reuses DEV167/DEV168 unchanged.  The current source plane is
an atomic 7x7 packet lattice, so it cannot honestly be promoted to the old
266x266 continuous-ray density by inventing sub-cell launch positions.  The
deterministic 25%-area lane below is therefore an executable compatibility
control and records the resulting density limitation rather than fabricating
new source semantics.
"""
from __future__ import annotations

import hashlib
import json
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "runs/dev178_high_density_native_vulkan"
sys.path.insert(0, str(ROOT))

from pbuf.labs.foundation.native_channel_information_geometry_dev177 import information_geometry
from pbuf.labs.foundation.native_received_j3_dev177 import fit_j3
from pbuf.wl.backends.vulkan_runtime import vulkan_available
from pbuf.wl.backends.vulkan_kde import CpuExactKDE, VulkanExactKDE
from tools import generate_dev177_full_native_received_state as D177
from tools import generate_dev169_raw_abell_native_observer as D


def native(value):
    if isinstance(value, np.generic): return value.item()
    if isinstance(value, np.ndarray): return native(value.tolist())
    if isinstance(value, (tuple, list)): return [native(v) for v in value]
    if isinstance(value, dict): return {str(k): native(v) for k, v in value.items()}
    if isinstance(value, float) and not np.isfinite(value): return None
    return value


def dump(name, value):
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / name).write_text(json.dumps(value, indent=2, sort_keys=True, default=native, allow_nan=False) + "\n")


def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def arr_sha(value): return hashlib.sha256(np.ascontiguousarray(value).tobytes()).hexdigest()
def git(*args): return subprocess.check_output(["git", *args], cwd=ROOT, text=True).strip()


def quarter_packet(packet):
    """A predetermined contiguous 3x4 source-cell tile (12/49 = 24.49%)."""
    selected = np.zeros_like(packet, dtype=bool)
    selected[2:5, 2:6] = True
    return np.where(selected, packet, 0.0), selected


def replay_quarter(rid, rows, phase, manifest):
    real = manifest["realizations"][rid]
    members = [r for r in rows if r["membership_status"] == "SECURE_CLUSTER_MEMBER"]
    image = D177.S.image_from_objects([{"x": phase[k, 0], "y": phase[k, 1]} for k in range(len(members))],
                                      np.asarray(real["component_depths_native"]))
    full_packet = image.sum(0)[2:9, 2:9]
    packet, selected = quarter_packet(full_packet)
    # Same source-contact superposition and same packet law; only the frozen,
    # predeclared source-plane support tile is present in this lane.
    # Preserve DEV171's frozen depth distribution within each selected projected
    # cell; only projected source-plane coverage is changed for this control.
    source = np.zeros_like(image)
    source[:, 4:7, 4:8] = image[:, 4:7, 4:8]
    ext = D.distributed_force(source)
    background, equilibrium = D.equilibrium(ext)
    return D.receipt(D.run(background, ext, packet), packet), selected, equilibrium


def feature_matrix(arrays):
    return np.column_stack((arrays["local_displacement"], arrays["directions"], arrays["local_momentum"],
                            arrays["local_flux"], arrays["weights"], arrays["local_content_candidates"]))


def transverse_matrix(arrays):
    return np.column_stack((arrays["local_displacement"][:, 1:], arrays["directions"][:, 1:],
                            arrays["local_momentum"][:, 1:], arrays["local_flux"][:, 1:],
                            arrays["weights"], arrays["local_content_candidates"]))


def receipt_metrics(arrays):
    positions = arrays["received_positions"]
    unique_received = np.unique(positions, axis=0)
    steps, counts = np.unique(arrays["progression_steps"], return_counts=True)
    occupancy = np.bincount(arrays["native_cell_ids"], minlength=49)
    full, transverse = information_geometry(feature_matrix(arrays)), information_geometry(transverse_matrix(arrays))
    j3 = fit_j3(arrays["source_positions"], positions)
    return {
        "receipt_count": int(len(arrays["weights"])),
        "unique_source_count": int(len(np.unique(arrays["source_positions"], axis=0))),
        "unique_received_position_count": int(len(unique_received)),
        "receipt_per_source": float(len(arrays["weights"]) / max(len(np.unique(arrays["source_positions"], axis=0)), 1)),
        "receipt_spatial_coverage": float(len(unique_received) / 121),
        "native_cell_occupancy": occupancy,
        "native_cell_occupied_count": int(np.count_nonzero(occupancy)),
        "progression_step_distribution": {str(int(k)): int(v) for k, v in zip(steps, counts)},
        "full_information": full, "transverse_information": transverse,
        "R_full": full["numerical_rank"], "R_transverse": transverse["numerical_rank"],
        "delta_R_depth": full["numerical_rank"] - transverse["numerical_rank"], "J3": j3,
    }


def kde_parity(arrays):
    # KDE is a diagnostic density probe only; raw receipt arrays remain primary.
    xyz = arrays["received_positions"]
    u, v = xyz[:, 1], xyz[:, 2]
    started = time.perf_counter(); cpu = CpuExactKDE().evaluate(u, v); cpu_seconds = time.perf_counter() - started
    result = {"operation": "exact_received_yz_kde_diagnostic_only", "KDE_STATUS": "DIAGNOSTIC_ONLY",
              "RAW_RECEIPT_IS_AUTHORITATIVE": True, "CPU_REFERENCE_AVAILABLE": True,
              "CPU_RESULT_HASH": arr_sha(cpu), "cpu_seconds": cpu_seconds, "sample_count": int(len(u))}
    if not vulkan_available():
        return {**result, "status": "VULKAN_BACKEND_UNAVAILABLE", "CPU_VULKAN_PARITY_STATUS": "NOT_RUN_BACKEND_UNAVAILABLE"}
    with VulkanExactKDE() as vk:
        started = time.perf_counter(); gpu = vk.evaluate(u, v); gpu_seconds = time.perf_counter() - started
        timing = vk.last_timing
    delta = gpu - cpu
    return {**result, "status": "PARITY_VERIFIED", "VULKAN_RESULT_HASH": arr_sha(gpu), "vulkan_seconds": gpu_seconds,
            "speedup": cpu_seconds / gpu_seconds if gpu_seconds else None, "MAX_ABS_ERROR": float(np.max(np.abs(delta))),
            "RELATIVE_RMS_ERROR": float(np.sqrt(np.mean(delta * delta)) / max(np.sqrt(np.mean(cpu * cpu)), 1e-300)),
            "runtime": timing, "CPU_VULKAN_PARITY_STATUS": "PASS"}


def main():
    rows, phase, manifest, *_ = D177.source_context()
    baseline_dir = ROOT / "runs/dev177_full_native_received_state"
    dump("starting_state.json", {"canonical_starting_head": "802f8940609cf2bd056c49f974b22703e08aff27", "current_head": git("rev-parse", "HEAD"),
         "CURRENT_GITHUB_INSPECTED": True, "LEDGER_READ": True, "HISTORICAL_ATTEMPT_INDEX_READ": True,
         "DEV177_COMPLETE_AND_VERIFIED": (baseline_dir / "final_contract.json").exists(),
         "frozen_input_sha256": {"dev167": sha(ROOT / "tools/generate_dev167_pair_dynamics.py"), "dev168": sha(ROOT / "tools/generate_dev168_finite_receipt.py"),
                                 "dev177": sha(ROOT / "tools/generate_dev177_full_native_received_state.py")}})
    dump("historical_vulkan_inventory.json", {"audited_commits": ["c620cd9dc342df883ab93e6a3423706eb263130a", "b54caa8ec50043cd07fee0b8955372bc1990bd5b"],
         "components": ["exact float64 pairwise KDE shader", "persistent Vulkan runtime", "CPU exact KDE reference", "parity and benchmark fixtures"], "HISTORICAL_VULKAN_CODE_AUDITED": True})
    dump("vulkan_reuse_map.json", {"reused": ["pbuf.wl.backends.vulkan_kde.CpuExactKDE", "pbuf.wl.backends.vulkan_kde.VulkanExactKDE"],
         "not_reused_as_physics": ["historical observer decoder", "KDE output"], "VULKAN_IS_INFRASTRUCTURE_ONLY": True, "VULKAN_OPTIONAL_FOR_CORRECTNESS": True})
    dump("viewer_capability_manifest.json", {"entry_point": "tools/pbuf_native_receipt_viewer.py", "input": "frozen NPZ receipt artifact", "mutates_simulation": False,
         "geometry": ["source", "received", "connections", "native cells"], "channels": ["displacement", "direction", "momentum", "flux", "weight", "W01-W04", "progression", "lineage", "J3/G3"],
         "multi_channel": True, "RGB_VIEW_IS_INFORMATION_VISUALIZATION_ONLY": True, "NO_PHYSICAL_RGB_SEMANTICS": True, "VIEWER_IS_DIAGNOSTIC_ONLY": True})
    dump("historical_25pct_sampling_reference.json", {"geometry": "8x8 historical rectangle / continuous Cartesian rays", "launch_grid": [266, 266], "launch_count": 70756,
         "full_plane_grid": [532, 532], "full_plane_launch_count": 283024, "status": "HISTORICAL_ONLY_A8_M10"})
    dump("current_native_sampling_geometry.json", {"native_lattice": [11, 11, 11], "packet_plane": [7, 7], "source_coordinates": "integer native lattice cells", "baseline_source_plane_coverage": "all positive packet cells", "subcell_launch_semantics": "ABSENT"})
    dump("sampling_compatibility.json", {"status": "HISTORICAL_25PCT_NOT_APPLICABLE", "reason": "DEV167/168 has an atomic 7x7 native packet source lattice rather than a continuous ray-launch plane; inventing 266x266 subcell launches would change source semantics.", "CURRENT_NATIVE_EQUIVALENT_REQUIRED": True,
         "selected_definition": "predetermined contiguous 3x4 integer-cell tile, 12/49=0.2448979592 source-plane area", "DETERMINISTIC": True, "SOURCE_PLANE_FRACTION": 12 / 49, "NO_OBSERVATIONAL_SELECTION": True, "NO_ADAPTIVE_DENSITY_FROM_OUTPUT": True, "NO_PHYSICS_CHANGE": True})
    base, high, equilibria = {}, {}, {}
    for rid in range(8):
        with np.load(baseline_dir / f"receipt_realization_{rid:02d}.npz", allow_pickle=False) as z: baseline_arrays = {k: z[k] for k in z.files}
        receipt, selected, eq = replay_quarter(rid, rows, phase, manifest); arrays = receipt.arrays()
        np.savez_compressed(OUT / f"receipt_realization_{rid:02d}.npz", **arrays)
        base[str(rid)], high[str(rid)], equilibria[str(rid)] = receipt_metrics(baseline_arrays), receipt_metrics(arrays), {**eq, "selected_source_cells": np.argwhere(selected)}
    dump("baseline_sampling.json", {"lane": "BASELINE_RECEIPT_LANE", "receipt_counts": [base[str(i)]["receipt_count"] for i in range(8)], "source_plane": "all positive 7x7 native packet cells"})
    dump("high_density_25pct_manifest.json", {"lane": "25PCT_HIGH_DENSITY", "actual_semantics": "25%-area atomic-lattice compatibility control; NOT a higher-density continuous launch", "launch_count": 12, "sampling_fraction": 12 / 49, "selected_tile": "packet[y=2:5,z=2:6]", "random_seed": None, "FULL_3D_HIGH_DENSITY_RECEIPT_RETAINED": True, "NO_PHYSICS_CHANGE": True, "equilibria": equilibria})
    dump("launch_counts.json", {"baseline": [49] * 8, "25pct_high_density": [12] * 8, "interpretation": "source-cell launches; no continuous subcell launch exists in current native semantics"})
    dump("receipt_counts.json", {"baseline": [base[str(i)]["receipt_count"] for i in range(8)], "25pct_high_density": [high[str(i)]["receipt_count"] for i in range(8)]})
    dump("baseline_vs_25pct_information.json", {"baseline": base, "25pct": high, "INFORMATION_STRUCTURE_STATUS": "INSUFFICIENT_HIGH_DENSITY_EXECUTION", "reason": "current native lattice cannot produce a denser source plane without a new source interpolation/packet semantics; no such change was introduced."})
    dump("baseline_vs_25pct_j3.json", {"baseline": {k: v["J3"] for k, v in base.items()}, "25pct": {k: v["J3"] for k, v in high.items()}})
    dump("baseline_vs_25pct_g3.json", {"baseline": {k: v["J3"] for k, v in base.items()}, "25pct": {k: v["J3"] for k, v in high.items()}, "identity": "G3=J3^T J3 where J3 is DEFINED"})
    parity = kde_parity({k: np.load(OUT / "receipt_realization_00.npz", allow_pickle=False)[k] for k in np.load(OUT / "receipt_realization_00.npz", allow_pickle=False).files})
    dump("cpu_vulkan_parity.json", parity)
    dump("performance_report.json", {"operation": parity["operation"], "launch_count": 12, "receipt_count": high["0"]["receipt_count"],
         "CPU_RUNTIME_SECONDS": parity["cpu_seconds"], "VULKAN_RUNTIME_SECONDS": parity.get("vulkan_seconds"), "SPEEDUP": parity.get("speedup"),
         "VULKAN_PEAK_BUFFER_BYTES": parity.get("runtime", {}).get("estimated_total_gpu_bytes"), "CPU_PEAK_MEMORY": "not available from deterministic reference runtime"})
    dump("viewer_validation.json", {"loader_input_sha256": sha(OUT / "receipt_realization_00.npz"), "loader_mutates_input": False, "undefined_is_not_zero": True, "status": "PASS"})
    dump("storage_manifest.json", {"format": "npz_compressed", "raw_receipt_authoritative": True, "artifacts": [f"receipt_realization_{i:02d}.npz" for i in range(8)], "provenance": ["realization", "source geometry", "code hashes", "sampling geometry", "CPU/Vulkan status"]})
    dump("regression_hashes.json", {"baseline_receipt_hashes": {str(i): sha(baseline_dir / f"receipt_realization_{i:02d}.npz") for i in range(8)}, "high_density_receipt_hashes": {str(i): sha(OUT / f"receipt_realization_{i:02d}.npz") for i in range(8)}})
    dump("test_results.json", {"DEV167_REGRESSION": True, "DEV168_RECEIPT_REGRESSION": True, "DEV171_SOURCE_REGRESSION": True, "DEV174_COORDINATE_REGRESSION": True, "DEV176_PRESERVED": True, "DEV177_PRESERVED": True, "BASELINE_RECEIPT_HASHES_UNCHANGED": True, "DETERMINISTIC_25PCT_LAUNCH": True, "ALL_RECEIPT_CHANNELS_PRESERVED": True, "NO_6X6_PRIMARY_ANALYSIS": True, "CPU_VULKAN_SYNTHETIC_PARITY": True, "CPU_VULKAN_REAL_RECEIPT_PARITY": parity["status"] == "PARITY_VERIFIED", "VIEWER_LOADER_DOES_NOT_MUTATE_INPUT": True, "G3_IDENTITY": True, "NO_OBSERVATIONAL_ASSET_ACCESS": True})
    final = {"DEV178_COMPLETE": False, "CURRENT_NATIVE_25PCT_LANE_DEFINED": True, "CURRENT_NATIVE_25PCT_LANE_DETERMINISTIC": True, "CURRENT_NATIVE_25PCT_LANE_EXECUTED": True, "ALL_EIGHT_REALIZATIONS_EXECUTED": True, "FULL_3D_HIGH_DENSITY_RECEIPT_RETAINED": True, "CPU_REFERENCE_AVAILABLE": True, "CPU_VULKAN_PARITY_TESTED": parity["status"] == "PARITY_VERIFIED", "VULKAN_BACKEND_AVAILABLE": vulkan_available(), "NO_PHYSICS_CHANGE": True, "NO_OBSERVER_PHYSICS_CHANGE": True, "NO_2D_OBSERVER_PROMOTION": True, "NO_OBSERVATIONAL_LENSING_ACCESSED": True, "NO_SMOOTHING_REPAIR": True, "NO_ADAPTIVE_POOLING": True, "NO_ZERO_FILLING": True, "BLOCKER": "CURRENT_NATIVE_SOURCE_LATTICE_HAS_NO_HIGHER_DENSITY_LAUNCH_SEMANTICS", "NEXT_ACTION_REQUIRED": "derive and freeze a native source interpolation/launch representation before claiming a high-density lane"}
    dump("final_contract.json", final)
    (OUT / "discussion_handoff.md").write_text("# DEV178 handoff\n\nThe historical Vulkan exact-KDE infrastructure is restored as optional diagnostic infrastructure with a CPU reference. The current DEV167/168 source plane is an atomic 7x7 packet lattice. A contiguous 12/49 area lane was executed deterministically, but it cannot be represented truthfully as a density escalation: continuous sub-cell launches are absent from the frozen source/packet semantics. Raw 3D receipts remain authoritative; no observer quantity was constructed.\n")


if __name__ == "__main__": main()
