#!/usr/bin/env python3
"""PBUF FOUNDATION — INDEPENDENT SOURCE TRAINING-WHEELS-OFF 001 COMMON-FOOTPRINT FIX.

Corrects the original training-wheels-off audit after the first cluster exposed a
real field-of-view mismatch: the HST F160W main-cluster mosaic overlapped only a
small fraction of the much larger kappa reconstruction grid (1672/32400 pixels in
the failed run). The original 25% whole-grid overlap gate was therefore not a
valid scientific requirement.

This correction does NOT lower the gate and then pretend the fields cover the
same sky area. Instead it derives, using WCS metadata only, the actual common HST
/ benchmark sky footprint and constructs a fixed 64x64 comparison/source grid
inside that common footprint.

Ordering remains strict:

  HST image + HST WCS
      + benchmark FITS HEADER/SHAPE only
      -> common-footprint geometry
      -> independent HST luminous proxy on 64x64 common grid
      -> frozen PBUF chain COMPLETE

  ONLY THEN:
      load benchmark kappa pixel values
      -> sample them onto the exact same common sky grid
      -> build assisted control on that same footprint
      -> end-of-chain comparisons

Thus observed kappa pixel morphology is still absent from the independent PBUF
lane, while the assisted and independent lanes now occupy the same physical sky
footprint instead of comparing a small HST pointing against the full kappa map.

No conventional gravitational law, fitting, smoothing, morphology matching,
translation, rotation, cluster-specific tuning, or benchmark-pixel-dependent grid
selection is introduced.
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from scipy.ndimage import map_coordinates

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))

import pbuf.labs.foundation.independent_source_training_wheels_off001 as LAB
from a8_three_dimensional_projection_lab001 import construct_rho_3d

LAB_ID = "PBUF-FOUNDATION-INDEPENDENT-SOURCE-TRAINING-WHEELS-OFF-001-COMMON-FOOTPRINT-FIX"
OUT = ROOT / "runs" / "independent_source_training_wheels_off001_common_footprint_fix"
DOWNLOADS = OUT / "downloads"
COMMON_N = 64
MIN_RAW_OVERLAP_PIXELS = 100
MIN_COMMON_GRID_VALID_FRACTION = 0.50


def _positive_normalized(field: np.ndarray, valid: np.ndarray) -> np.ndarray:
    x = np.asarray(field, dtype=np.float64)
    valid = np.asarray(valid, dtype=bool)
    out = np.zeros_like(x)
    use = valid & np.isfinite(x)
    if not np.any(use):
        raise RuntimeError("common-footprint field has no finite valid pixels")
    pos = np.maximum(x[use], 0.0)
    mx = float(np.max(pos)) if pos.size else 0.0
    if not np.isfinite(mx) or mx <= 0.0:
        raise RuntimeError("common-footprint positive normalized proxy is empty")
    out[use] = pos / mx
    return out


def _common_grid_from_wcs(hst_data: np.ndarray, hst_hdr, target_hdr, target_shape):
    """Derive a 64x64 common sky footprint using WCS metadata only.

    No benchmark pixel values are read here. The full benchmark pixel grid is
    projected into the HST WCS to discover which target pixels are geometrically
    covered by the HST image. A bounding box of that metadata-only overlap is then
    sampled on a fixed 64x64 grid. Pixels in that box that still fall outside the
    exact HST footprint remain masked rather than extrapolated.
    """
    hst = np.asarray(np.squeeze(hst_data), dtype=np.float64)
    if hst.ndim != 2:
        raise RuntimeError(f"HST image must be 2-D after squeeze, got {hst.shape}")
    nyh, nxh = hst.shape
    hw = WCS(hst_hdr)
    tw = WCS(target_hdr)

    nyt, nxt = target_shape
    yy, xx = np.indices((nyt, nxt), dtype=np.float64)
    world = tw.pixel_to_world_values(xx, yy)
    hx, hy = hw.world_to_pixel_values(*world)
    raw_overlap = (
        np.isfinite(hx) & np.isfinite(hy)
        & (hx >= 0.0) & (hx <= nxh - 1.0)
        & (hy >= 0.0) & (hy <= nyh - 1.0)
    )
    raw_count = int(np.count_nonzero(raw_overlap))
    if raw_count < MIN_RAW_OVERLAP_PIXELS:
        raise RuntimeError(
            f"insufficient geometric HST/benchmark WCS overlap even for common-footprint construction: "
            f"{raw_count}/{raw_overlap.size}"
        )

    rr, cc = np.where(raw_overlap)
    r0, r1 = int(rr.min()), int(rr.max())
    c0, c1 = int(cc.min()), int(cc.max())
    if r1 <= r0 or c1 <= c0:
        raise RuntimeError("degenerate HST/benchmark common-footprint bounding box")

    gy = np.linspace(float(r0), float(r1), COMMON_N)
    gx = np.linspace(float(c0), float(c1), COMMON_N)
    target_x, target_y = np.meshgrid(gx, gy, indexing="xy")

    world_grid = tw.pixel_to_world_values(target_x, target_y)
    hgrid_x, hgrid_y = hw.world_to_pixel_values(*world_grid)
    hst_valid = (
        np.isfinite(hgrid_x) & np.isfinite(hgrid_y)
        & (hgrid_x >= 0.0) & (hgrid_x <= nxh - 1.0)
        & (hgrid_y >= 0.0) & (hgrid_y <= nyh - 1.0)
    )
    valid_fraction = float(np.mean(hst_valid))
    if valid_fraction < MIN_COMMON_GRID_VALID_FRACTION:
        raise RuntimeError(
            f"common-footprint 64x64 grid remains too sparsely covered by HST: "
            f"{int(np.count_nonzero(hst_valid))}/{hst_valid.size} ({valid_fraction:.6f})"
        )

    finite_native = np.isfinite(hst)
    if not np.any(finite_native):
        raise RuntimeError("HST image contains no finite pixels")
    native_fill = float(np.nanmedian(hst[finite_native]))
    src = np.where(finite_native, hst, native_fill)
    sampled = map_coordinates(
        src,
        np.array([hgrid_y, hgrid_x], dtype=np.float64),
        order=1,
        mode="constant",
        cval=np.nan,
    )
    valid = hst_valid & np.isfinite(sampled)
    if float(np.mean(valid)) < MIN_COMMON_GRID_VALID_FRACTION:
        raise RuntimeError("finite HST samples do not adequately cover common 64x64 grid")

    vals = sampled[valid]
    background = float(np.median(vals))
    luminous = np.zeros_like(sampled)
    luminous[valid] = np.maximum(sampled[valid] - background, 0.0)
    maxv = float(np.max(luminous))
    if not np.isfinite(maxv) or maxv <= 0.0:
        raise RuntimeError("positive F160W luminous proxy is empty on common footprint")
    rho2 = luminous / maxv

    geometry = {
        "target_pixel_x": target_x,
        "target_pixel_y": target_y,
        "valid_mask": valid,
    }
    diagnostics = {
        "target_full_pixel_count": int(raw_overlap.size),
        "raw_wcs_overlap_pixel_count": raw_count,
        "raw_wcs_overlap_fraction": float(np.mean(raw_overlap)),
        "common_bbox_row_min": r0,
        "common_bbox_row_max": r1,
        "common_bbox_col_min": c0,
        "common_bbox_col_max": c1,
        "common_grid_n": COMMON_N,
        "common_grid_pixel_count": int(valid.size),
        "common_grid_valid_pixel_count": int(np.count_nonzero(valid)),
        "common_grid_valid_fraction": float(np.mean(valid)),
        "fixed_background_estimator": "median_of_finite_HST_samples_on_metadata_only_common_WCS_footprint",
        "background_value": background,
        "positive_luminous_common_max": maxv,
        "positive_luminous_common_rms": LAB._rms(luminous),
    }
    return rho2, luminous, geometry, diagnostics


def _independent_source(cluster: dict) -> dict:
    # Header/shape only from the benchmark at this stage; NO kappa pixels.
    kpath, target_hdr, target_shape = LAB._target_header_and_shape(cluster)
    url, candidates = LAB._discover_f160w_url(cluster["id"])
    local = DOWNLOADS / cluster["id"] / url.split("/")[-1]
    LAB._download(url, local)
    with fits.open(local, memmap=True) as hdul:
        hst = np.asarray(hdul[0].data, dtype=np.float64)
        hst_hdr = hdul[0].header.copy()

    rho2, luminous, geometry, align = _common_grid_from_wcs(
        hst, hst_hdr, target_hdr, target_shape
    )
    rho3 = construct_rho_3d(rho2, LAB.BASE.NZ, profile=LAB.BASE.PROFILE)
    return {
        "rho2": rho2,
        "rho3": rho3,
        "luminous_common": luminous,
        "geometry": geometry,
        "kappa_path_for_header_only": str(kpath),
        "hst_url": url,
        "hst_discovered_candidates": candidates,
        "hst_local_path": str(local),
        "hst_sha256": LAB._sha_file(local),
        "alignment": align,
        "source_role": "independent_HST_F160W_positive_luminous_structure_proxy_on_metadata_only_common_sky_footprint_not_mass_map",
        "observed_kappa_pixel_values_used": False,
    }


def _benchmark_on_common_grid(kappa_native: np.ndarray, geometry: dict):
    gx = np.asarray(geometry["target_pixel_x"], dtype=np.float64)
    gy = np.asarray(geometry["target_pixel_y"], dtype=np.float64)
    valid = np.asarray(geometry["valid_mask"], dtype=bool)
    sampled = map_coordinates(
        np.asarray(kappa_native, dtype=np.float64),
        np.array([gy, gx], dtype=np.float64),
        order=1,
        mode="constant",
        cval=np.nan,
    )
    comparison = np.where(valid & np.isfinite(sampled), sampled, np.nan)
    source_values = np.where(valid & np.isfinite(sampled), sampled, 0.0)
    assisted_rho2 = _positive_normalized(source_values, valid & np.isfinite(sampled))
    return comparison, assisted_rho2


def _run_cluster(cluster: dict):
    cid = cluster["id"]

    # INDEPENDENT LANE FIRST. No observed kappa pixel values are loaded.
    independent_source = _independent_source(cluster)
    independent_chain = LAB._run_chain_from_rho3(
        independent_source["rho3"], observed_for_first_step=None
    )

    # ONLY NOW load kappa pixel values and sample the exact metadata-defined common grid.
    kpath = LAB._kappa_path(cluster)
    with fits.open(kpath) as hdul:
        kappa_native = np.asarray(hdul[0].data, dtype=np.float64)
    observed, assisted_rho2 = _benchmark_on_common_grid(
        kappa_native, independent_source["geometry"]
    )
    assisted_rho3 = construct_rho_3d(
        assisted_rho2, LAB.BASE.NZ, profile=LAB.BASE.PROFILE
    )
    assisted_chain = LAB._run_chain_from_rho3(
        assisted_rho3, observed_for_first_step=observed
    )

    indep_metrics = LAB._lane_vs_benchmark(
        "independent", independent_source["rho2"], independent_chain, observed
    )
    assisted_metrics = LAB._lane_vs_benchmark(
        "assisted", assisted_rho2, assisted_chain, observed
    )

    stage_pairs = {
        "source_rho2": (assisted_rho2, independent_source["rho2"]),
        "m10_los_mag": (assisted_chain["los_mag"], independent_chain["los_mag"]),
        "final_angular_centroid_mag": (
            assisted_chain["final_ang"]["angular_centroid_mag"],
            independent_chain["final_ang"]["angular_centroid_mag"],
        ),
        "final_angular_spread_rms": (
            assisted_chain["final_ang"]["angular_spread_rms"],
            independent_chain["final_ang"]["angular_spread_rms"],
        ),
        "final_angular_rms_angle_mag": (
            assisted_chain["final_ang"]["angular_rms_angle_mag"],
            independent_chain["final_ang"]["angular_rms_angle_mag"],
        ),
    }
    pair_metrics = {}
    for name, (a, i) in stage_pairs.items():
        p, s, n = LAB._corr(a, i)
        pair_metrics[f"assisted_vs_independent_{name}_pearson"] = p
        pair_metrics[f"assisted_vs_independent_{name}_spearman"] = s
        pair_metrics[f"assisted_vs_independent_{name}_count"] = n

    det = LAB._paired_detector(independent_chain, assisted_chain)
    summary = {
        "cluster_id": cid,
        "candidate_id": LAB.BASE.CANDIDATE_ID,
        "physical_source_representation": LAB.BASE.PHYSICAL_SOURCE,
        "assisted_source_role": "observed_kappa_positive_normalized_proxy_on_same_metadata_defined_common_sky_footprint_training_wheels_control",
        "independent_source_role": independent_source["source_role"],
        "independent_source_filter": "HST_WFC3_F160W",
        "independent_source_background_rule": independent_source["alignment"]["fixed_background_estimator"],
        "independent_source_observed_kappa_pixel_values_used": False,
        "benchmark_values_loaded_after_independent_lane_complete": True,
        "benchmark_role": "external_morphology_comparison_only_for_independent_lane",
        "independent_source_limit": "luminous_structure_proxy_not_baryonic_or_total_mass_map",
        "common_footprint_role": "WCS_metadata_only_geometry_shared_by_both_lanes_no_kappa_pixel_values_used_to_select_crop",
        "geometry_lane": "G3D_LOS_consistent_diagnostic",
        "observer_readout": "validated_per_source_angular_moments_plus_native_global_angular_histogram",
        "n_photons": int(len(independent_chain["tx"])),
        "source_supported_bins": int(len(independent_chain["groups"])),
        "hst_url": independent_source["hst_url"],
        "hst_sha256": independent_source["hst_sha256"],
        "hst_candidate_count": len(independent_source["hst_discovered_candidates"]),
        **{f"hst_{k}": v for k, v in independent_source["alignment"].items()},
        "independent_g3d_unit_speed_max_error": float(independent_chain["g3d"]["max_unit_speed_error"]),
        "independent_g3d_unit_speed_pass": bool(independent_chain["g3d"]["max_unit_speed_error"] <= LAB.UNIT_SPEED_TOL),
        "independent_first_step_exact_max_vector_error": independent_chain["first"]["first_step_exact_max_vector_error"],
        "independent_first_step_exact_pass": independent_chain["first"]["first_step_exact_pass"],
        "assisted_g3d_unit_speed_max_error": float(assisted_chain["g3d"]["max_unit_speed_error"]),
        "assisted_g3d_unit_speed_pass": bool(assisted_chain["g3d"]["max_unit_speed_error"] <= LAB.UNIT_SPEED_TOL),
        "assisted_first_step_exact_max_vector_error": assisted_chain["first"]["first_step_exact_max_vector_error"],
        "assisted_first_step_exact_pass": assisted_chain["first"]["first_step_exact_pass"],
        **assisted_metrics,
        **indep_metrics,
        **pair_metrics,
        "paired_detector_resolution": det["detector_resolution"],
        "paired_detector_half_range": det["detector_half_range"],
        "paired_detector_total_variation_distance": det["total_variation_distance"],
        "paired_detector_jensen_shannon_divergence_nats": det["jensen_shannon_divergence_nats"],
        "paired_detector_normalized_l2_distance": det["normalized_l2_distance"],
        "paired_detector_image_pearson": det["image_pearson"],
        "paired_detector_image_spearman": det["image_spearman"],
        "paired_detector_assisted_occupied_bins": det["assisted_occupied_bins"],
        "paired_detector_independent_occupied_bins": det["independent_occupied_bins"],
        "paired_detector_occupied_bin_intersection": det["occupied_bin_intersection"],
    }
    arrays = {
        "observed_kappa_reference_only_common_footprint": observed,
        "common_footprint_valid_mask": independent_source["geometry"]["valid_mask"],
        "common_footprint_target_pixel_x": independent_source["geometry"]["target_pixel_x"],
        "common_footprint_target_pixel_y": independent_source["geometry"]["target_pixel_y"],
        "independent_luminous_common": independent_source["luminous_common"],
        "assisted_source_rho2": assisted_rho2,
        "independent_source_rho2": independent_source["rho2"],
        "assisted_m10_los_mag": assisted_chain["los_mag"],
        "independent_m10_los_mag": independent_chain["los_mag"],
        "assisted_final_angular_centroid_mag": assisted_chain["final_ang"]["angular_centroid_mag"],
        "independent_final_angular_centroid_mag": independent_chain["final_ang"]["angular_centroid_mag"],
        "assisted_final_angular_spread_rms": assisted_chain["final_ang"]["angular_spread_rms"],
        "independent_final_angular_spread_rms": independent_chain["final_ang"]["angular_spread_rms"],
        "assisted_final_angular_rms_angle_mag": assisted_chain["final_ang"]["angular_rms_angle_mag"],
        "independent_final_angular_rms_angle_mag": independent_chain["final_ang"]["angular_rms_angle_mag"],
        "paired_detector_assisted_hist_r64": det["assisted_hist"],
        "paired_detector_independent_hist_r64": det["independent_hist"],
    }
    return summary, arrays


def _mean(summaries, key):
    vals = np.asarray([s[key] for s in summaries], dtype=np.float64)
    return float(np.nanmean(vals)) if np.any(np.isfinite(vals)) else float("nan")


def main() -> int:
    started = time.perf_counter()
    OUT.mkdir(parents=True, exist_ok=True)
    DOWNLOADS.mkdir(parents=True, exist_ok=True)
    # Reuse downloaded files from failed original run when available without moving them.
    LAB.DOWNLOADS = LAB.OUT / "downloads"

    repo = LAB._repo_state()
    LAB._write_json(OUT / "repository_state.json", repo)
    if repo["branch"] != "main" or repo["tracked_changes"] or repo["staged_changes"]:
        v = {"lab_id": LAB_ID, "outcome": "REPOSITORY_GATE_FAILURE", "head_sha": repo["head_sha"]}
        LAB._write_json(OUT / "validation.json", v)
        print(json.dumps(v, indent=2))
        return 2

    summaries = []
    failures = []
    for cluster in LAB.BASE.CLUSTERS:
        cid = cluster["id"]
        print(f"[{cid}] TRAINING WHEELS OFF COMMON FOOTPRINT: HST F160W -> frozen PBUF chain")
        try:
            summary, arrays = _run_cluster(cluster)
            summaries.append(summary)
            cdir = OUT / "clusters" / cid
            cdir.mkdir(parents=True, exist_ok=True)
            LAB._write_json(cdir / "independent_source_summary.json", summary)
            np.savez_compressed(cdir / "independent_source_fields.npz", **arrays)
        except Exception as exc:
            failures.append({"cluster_id": cid, "error": repr(exc)})
            LAB._write_json(OUT / "cluster_failures.json", failures)
            raise

    LAB._write_csv(OUT / "independent_source_summary.csv", summaries)
    LAB._write_json(OUT / "cluster_failures.json", failures)

    validation = {
        "lab_id": LAB_ID,
        "outcome": "Outcome A — INDEPENDENT SOURCE TRAINING-WHEELS-OFF COMMON-FOOTPRINT AUDIT COMPLETE",
        "head_sha": repo["head_sha"],
        "candidate_id": LAB.BASE.CANDIDATE_ID,
        "physical_source_representation": LAB.BASE.PHYSICAL_SOURCE,
        "assisted_source_role": "observed_kappa_positive_normalized_proxy_on_same_metadata_defined_common_sky_footprint_training_wheels_control",
        "independent_source_role": "HST_F160W_positive_luminous_structure_proxy_common_footprint_no_kappa_pixel_values",
        "independent_source_limit": "luminous_structure_proxy_not_baryonic_or_total_mass_map",
        "common_footprint_rule": "metadata_only_HST_and_benchmark_WCS_overlap_bbox_resampled_to_fixed_64x64_before_any_kappa_pixels_are_loaded",
        "benchmark_role": "external_morphology_comparison_only_for_independent_lane",
        "benchmark_values_loaded_after_independent_lane_complete": True,
        "cluster_count_expected": len(LAB.BASE.CLUSTERS),
        "cluster_count_completed": len(summaries),
        "all_independent_source_observed_kappa_pixel_values_unused": bool(all(not s["independent_source_observed_kappa_pixel_values_used"] for s in summaries)),
        "all_common_grid_valid_fraction_pass": bool(all(s["hst_common_grid_valid_fraction"] >= MIN_COMMON_GRID_VALID_FRACTION for s in summaries)),
        "all_independent_g3d_unit_speed_pass": bool(all(s["independent_g3d_unit_speed_pass"] for s in summaries)),
        "all_independent_first_step_exact_pass": bool(all(s["independent_first_step_exact_pass"] for s in summaries)),
        "all_assisted_g3d_unit_speed_pass": bool(all(s["assisted_g3d_unit_speed_pass"] for s in summaries)),
        "all_assisted_first_step_exact_pass": bool(all(s["assisted_first_step_exact_pass"] for s in summaries)),
        "mean_assisted_source_rho2_vs_observed_pearson": _mean(summaries, "assisted_source_rho2_vs_observed_pearson"),
        "mean_independent_source_rho2_vs_observed_pearson": _mean(summaries, "independent_source_rho2_vs_observed_pearson"),
        "mean_assisted_vs_independent_source_rho2_pearson": _mean(summaries, "assisted_vs_independent_source_rho2_pearson"),
        "mean_assisted_m10_los_mag_vs_observed_pearson": _mean(summaries, "assisted_m10_los_mag_vs_observed_pearson"),
        "mean_independent_m10_los_mag_vs_observed_pearson": _mean(summaries, "independent_m10_los_mag_vs_observed_pearson"),
        "mean_assisted_vs_independent_m10_los_mag_pearson": _mean(summaries, "assisted_vs_independent_m10_los_mag_pearson"),
        "mean_assisted_final_angular_rms_angle_mag_vs_observed_pearson": _mean(summaries, "assisted_final_angular_rms_angle_mag_vs_observed_pearson"),
        "mean_independent_final_angular_rms_angle_mag_vs_observed_pearson": _mean(summaries, "independent_final_angular_rms_angle_mag_vs_observed_pearson"),
        "mean_assisted_vs_independent_final_angular_rms_angle_mag_pearson": _mean(summaries, "assisted_vs_independent_final_angular_rms_angle_mag_pearson"),
        "mean_paired_detector_image_pearson": _mean(summaries, "paired_detector_image_pearson"),
        "mean_paired_detector_total_variation_distance": _mean(summaries, "paired_detector_total_variation_distance"),
        "physics_change_authorized": False,
        "candidate_change_authorized": False,
        "production_geometry_change_authorized": False,
        "independent_source_as_final_mass_model_authorized": False,
        "observable_selection_authorized": False,
        "training_wheels_interpretation_required": True,
        "continue_pbuf_authorized": False,
        "stop_pbuf_authorized": False,
        "next_experiment_authorized": False,
        "duration_seconds": float(time.perf_counter() - started),
    }
    LAB._write_json(OUT / "validation.json", validation)
    LAB._write_json(OUT / "run.json", {"validation": validation, "summaries": summaries})
    print(json.dumps(validation, indent=2, default=LAB._json_default))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
