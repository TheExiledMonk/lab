#!/usr/bin/env python3
"""PBUF FOUNDATION — INDEPENDENT SOURCE TRAINING-WHEELS-OFF 001.

First deliberate removal of the observed-kappa morphology from the PBUF source lane.

The current validated diagnostic chain used a 3-D matter proxy constructed directly
from the observed weak-lensing kappa map. That was useful for localizing information
loss, but it cannot establish independent predictive power and can mask defects by
feeding target morphology into the source.

This audit freezes the downstream PBUF chain and changes only the upstream morphology
source in a paired experiment:

  ASSISTED control:
      observed kappa -> historical positive normalized proxy -> frozen PBUF chain

  INDEPENDENT lane:
      public HST Frontier Fields F160W cluster image -> WCS alignment only ->
      fixed positive luminous-structure proxy -> same frozen PBUF chain

The independent lane never receives observed kappa pixel values. The kappa FITS header
and shape are used only as an astrometric target grid so both products occupy the same
sky coordinates. Observed kappa pixel values are loaded only after the complete
independent PBUF lane has been constructed, and are used solely for end-of-chain
comparison.

Important limitation:
  F160W surface brightness is an independent luminous-structure tracer, NOT a baryonic
  mass map and NOT a total-matter map. It may contain foreground/background sources and
  lensed background light. Therefore failure here would expose dependence on the old
  kappa-assisted source but would not by itself falsify PBUF. Success would be much
  stronger because no kappa morphology entered the independent source lane.

No conventional gravitational law is introduced. No fitting, smoothing, morphology
matching, rotation, translation, cluster-specific tuning, or benchmark-dependent
normalization is allowed.
"""
from __future__ import annotations

import csv
import hashlib
import html.parser
import json
import subprocess
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from scipy.ndimage import map_coordinates

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))

import pbuf.labs.foundation.m10_coverage_25pct_science001 as BASE
import pbuf.labs.foundation.los_consistent_ray_geometry001 as GEO
import pbuf.labs.foundation.g3d_angular_received_distribution001 as ANG
import pbuf.labs.foundation.g3d_native_angular_detector_image001 as DET
from pbuf.core import los_projection as M14
from pbuf.core import observable_extraction as M16
from a8_three_dimensional_projection_lab001 import construct_rho_3d

LAB_ID = "PBUF-FOUNDATION-INDEPENDENT-SOURCE-TRAINING-WHEELS-OFF-001"
OUT = ROOT / "runs" / "independent_source_training_wheels_off001"
DOWNLOADS = OUT / "downloads"
BENCHMARK = ROOT / "PBUF_benchmark"

CHECKPOINT = GEO.CHECKPOINTS[-1]
EXPECTED_SUPPORT = GEO.EXPECTED_SUPPORT
UNIT_SPEED_TOL = GEO.UNIT_SPEED_TOL
FIRST_STEP_TOL = GEO.FIRST_STEP_TOL
DETECTOR_RESOLUTION = 64
RANGE_PAD = 1.0 + 1e-12

# This mapping is archive naming only. It contains no lensing information.
ARCHIVE_SLUG = {
    "Abell2744": "abell2744",
    "MACS0416": "macs0416",
    "MACS1149": "macs1149",
    "AbellS1063": "abells1063",
    "Abell370": "abell370",
}


class _HrefParser(html.parser.HTMLParser):
    def __init__(self):
        super().__init__()
        self.hrefs: list[str] = []

    def handle_starttag(self, tag, attrs):
        if tag.lower() != "a":
            return
        for k, v in attrs:
            if k.lower() == "href" and v:
                self.hrefs.append(v)


def _git(*args: str) -> str:
    return subprocess.check_output(["git", *args], cwd=str(ROOT), text=True).strip()


def _repo_state() -> dict:
    return {
        "repository": "TheExiledMonk/lab",
        "branch": _git("rev-parse", "--abbrev-ref", "HEAD"),
        "head_sha": _git("rev-parse", "HEAD"),
        "tracked_changes": _git("diff", "--name-only"),
        "staged_changes": _git("diff", "--name-only", "--cached"),
    }


def _json_default(o):
    if isinstance(o, np.floating): return float(o)
    if isinstance(o, np.integer): return int(o)
    if isinstance(o, np.bool_): return bool(o)
    if isinstance(o, Path): return str(o)
    return str(o)


def _write_json(path: Path, obj) -> None:
    path.write_text(json.dumps(obj, indent=2, default=_json_default))


def _write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        path.write_text("")
        return
    keys: list[str] = []
    for r in rows:
        for k in r:
            if k not in keys:
                keys.append(k)
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader(); w.writerows(rows)


def _sha_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(8 * 1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def _corr(a, b) -> tuple[float, float, int]:
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    if a.shape != b.shape:
        raise RuntimeError(f"correlation shape mismatch: {a.shape} vs {b.shape}")
    m = np.isfinite(a) & np.isfinite(b)
    n = int(m.sum())
    if n < 2:
        return float("nan"), float("nan"), n
    return float(M16.safe_pearson(a[m], b[m])), float(M16.safe_spearman(a[m], b[m])), n


def _rms(a) -> float:
    a = np.asarray(a, dtype=np.float64)
    m = np.isfinite(a)
    return float(np.sqrt(np.mean(a[m] * a[m]))) if np.any(m) else float("nan")


def _url_text(url: str) -> str:
    req = urllib.request.Request(url, headers={"User-Agent": "PBUF-foundation-audit/1.0"})
    with urllib.request.urlopen(req, timeout=60) as resp:
        return resp.read().decode("utf-8", errors="replace")


def _hrefs(url: str) -> list[str]:
    p = _HrefParser(); p.feed(_url_text(url))
    return p.hrefs


def _discover_f160w_url(cluster_id: str) -> tuple[str, list[str]]:
    slug = ARCHIVE_SLUG[cluster_id]
    root = f"https://archive.stsci.edu/pub/hlsp/frontier/{slug}/images/hst/"
    dirs = []
    for href in _hrefs(root):
        h = urllib.parse.unquote(href)
        if "v1.0" in h.lower() and h.endswith("/"):
            dirs.append(urllib.parse.urljoin(root, href))
    dirs = sorted(set(dirs))
    if not dirs:
        raise RuntimeError(f"{cluster_id}: no v1.0 HST directories discovered at {root}")

    candidates: list[str] = []
    needle = f"_{slug}_f160w_"
    for d in dirs:
        for href in _hrefs(d):
            name = urllib.parse.unquote(href).split("/")[-1]
            low = name.lower()
            if (
                "hlsp_frontier_hst_wfc3-60mas_" in low
                and needle in low
                and low.endswith("_drz.fits")
                and "hffpar" not in low
                and "-par" not in low
            ):
                candidates.append(urllib.parse.urljoin(d, href))
    candidates = sorted(set(candidates))
    if not candidates:
        raise RuntimeError(f"{cluster_id}: no main-field WFC3 60mas F160W drizzled FITS discovered")

    # Deterministic preference: exact /v1.0/ directory first, then epoch1, epoch2,
    # then lexical fallback. This uses archive version naming only, never benchmark data.
    def rank(u: str):
        lu = u.lower()
        if "/v1.0/" in lu: p = 0
        elif "v1.0-epoch1" in lu: p = 1
        elif "v1.0-epoch2" in lu: p = 2
        else: p = 3
        return (p, u)

    return sorted(candidates, key=rank)[0], candidates


def _download(url: str, path: Path) -> None:
    if path.exists() and path.stat().st_size > 0:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".part")
    req = urllib.request.Request(url, headers={"User-Agent": "PBUF-foundation-audit/1.0"})
    with urllib.request.urlopen(req, timeout=120) as resp, tmp.open("wb") as fh:
        while True:
            chunk = resp.read(8 * 1024 * 1024)
            if not chunk: break
            fh.write(chunk)
    tmp.replace(path)


def _kappa_path(cluster: dict) -> Path:
    return BENCHMARK / cluster["directory"] / f"hlsp_frontier_model_{cluster['slug']}_merten_v1_kappa.fits"


def _target_header_and_shape(cluster: dict):
    path = _kappa_path(cluster)
    if not path.exists(): raise FileNotFoundError(path)
    hdr = fits.getheader(path, 0)
    shape = (int(hdr["NAXIS2"]), int(hdr["NAXIS1"]))
    return path, hdr, shape


def _sample_hst_onto_target(hst_data: np.ndarray, hst_hdr, target_hdr, target_shape) -> tuple[np.ndarray, dict]:
    hst_data = np.asarray(np.squeeze(hst_data), dtype=np.float64)
    if hst_data.ndim != 2:
        raise RuntimeError(f"HST image must be 2-D after squeeze, got {hst_data.shape}")
    finite = np.isfinite(hst_data)
    if not np.any(finite):
        raise RuntimeError("HST image contains no finite pixels")
    fill = float(np.nanmedian(hst_data[finite]))
    src = np.where(finite, hst_data, fill)

    tw = WCS(target_hdr)
    hw = WCS(hst_hdr)
    ny, nx = target_shape
    yy, xx = np.indices((ny, nx), dtype=np.float64)
    world = tw.pixel_to_world_values(xx, yy)
    hx, hy = hw.world_to_pixel_values(*world)
    coords = np.array([hy, hx], dtype=np.float64)
    sampled = map_coordinates(src, coords, order=1, mode="constant", cval=np.nan)
    overlap = np.isfinite(sampled)
    if int(overlap.sum()) < max(100, int(0.25 * sampled.size)):
        raise RuntimeError(f"insufficient HST/kappa WCS overlap: {int(overlap.sum())}/{sampled.size}")

    vals = sampled[overlap]
    background = float(np.median(vals))
    luminous = np.zeros_like(sampled)
    luminous[overlap] = np.maximum(sampled[overlap] - background, 0.0)
    maxv = float(np.max(luminous))
    if not np.isfinite(maxv) or maxv <= 0.0:
        raise RuntimeError("positive F160W luminous-structure proxy is empty after fixed median subtraction")

    return luminous, {
        "target_pixel_count": int(sampled.size),
        "wcs_overlap_pixel_count": int(overlap.sum()),
        "wcs_overlap_fraction": float(overlap.mean()),
        "fixed_background_estimator": "median_of_finite_wcs_overlap",
        "background_value": background,
        "positive_luminous_native_max": maxv,
        "positive_luminous_native_rms": _rms(luminous),
    }


def _independent_source(cluster: dict) -> dict:
    # Only target HEADER/SHAPE are read here. No observed kappa pixel values.
    kpath, target_hdr, target_shape = _target_header_and_shape(cluster)
    url, candidates = _discover_f160w_url(cluster["id"])
    local = DOWNLOADS / cluster["id"] / url.split("/")[-1]
    _download(url, local)
    with fits.open(local, memmap=True) as hdul:
        hst = np.asarray(hdul[0].data, dtype=np.float64)
        hst_hdr = hdul[0].header.copy()
    luminous_native, align = _sample_hst_onto_target(hst, hst_hdr, target_hdr, target_shape)
    rho2 = BASE.construct_common_proxy(luminous_native, bins=BASE.OBS_BINS, extent=BASE.CFG["extent"])
    rho3 = construct_rho_3d(rho2, BASE.NZ, profile=BASE.PROFILE)
    return {
        "rho2": rho2,
        "rho3": rho3,
        "luminous_native": luminous_native,
        "kappa_path_for_header_only": str(kpath),
        "hst_url": url,
        "hst_discovered_candidates": candidates,
        "hst_local_path": str(local),
        "hst_sha256": _sha_file(local),
        "alignment": align,
        "source_role": "independent_HST_F160W_positive_luminous_structure_proxy_not_mass_map",
        "observed_kappa_pixel_values_used": False,
    }


def _run_chain_from_rho3(rho3: np.ndarray, observed_for_first_step=None) -> dict:
    state = BASE._evolve(BASE._initial_state(rho3))
    candidate = BASE._candidate(state)
    vector = BASE._interface_vector(candidate)
    los = M14.project_vector_to_image_plane(*vector, los_axis="z")
    Rx = np.asarray(los["comp_1"], dtype=np.float64)
    Ry = np.asarray(los["comp_2"], dtype=np.float64)
    los_mag = np.hypot(Rx, Ry)
    grid = np.linspace(-BASE.CFG["extent"], BASE.CFG["extent"], Rx.shape[0])
    field = {"xgrid": grid, "ygrid": grid, "rx": Rx, "ry": Ry}

    x0, y0, _, _ = BASE._launch_expanded_25pct()
    groups = GEO._source_groups(x0, y0)
    if len(groups) != EXPECTED_SUPPORT:
        raise RuntimeError(f"expected {EXPECTED_SUPPORT} source bins, got {len(groups)}")
    checkpoints, g3d = GEO._propagate_g3d(field, BASE.CFG["step"], BASE.CFG["steps"], x0, y0)
    if g3d["max_unit_speed_error"] > UNIT_SPEED_TOL:
        raise RuntimeError(f"G3D unit speed gate failed: {g3d['max_unit_speed_error']}")

    # The exact first-step vector identity itself does not require benchmark morphology.
    # GEO's helper also computes benchmark diagnostics, so use a zero array when no
    # observed benchmark is allowed in this lane.
    control_obs = np.zeros_like(los_mag) if observed_for_first_step is None else observed_for_first_step
    first = GEO._first_step_geometry(field, x0, y0, checkpoints[1], control_obs, los_mag)
    if not first["first_step_exact_pass"]:
        raise RuntimeError("first-step exact geometry gate failed")

    final_ang = ANG._angular_distribution_fields(checkpoints[CHECKPOINT], groups)
    gates = ANG._moment_gates(final_ang)
    if gates["second_moment_equals_cov_plus_centroid_outer_relative_rms_error"] > ANG.MOMENT_IDENTITY_TOL:
        raise RuntimeError("angular second-moment identity failed")
    if not gates["covariance_psd_pass"]:
        raise RuntimeError("angular covariance PSD gate failed")
    if not gates["direction_mean_vector_bound_pass"]:
        raise RuntimeError("angular mean-direction bound failed")

    snap = checkpoints[CHECKPOINT]
    vz = np.asarray(snap["vz"], dtype=np.float64)
    if float(np.min(np.abs(vz))) <= DET.VZ_MIN:
        raise RuntimeError("final tangent projection vz too small")
    tx = np.asarray(snap["vx"], dtype=np.float64) / vz
    ty = np.asarray(snap["vy"], dtype=np.float64) / vz

    return {
        "candidate": candidate,
        "los_mag": los_mag,
        "checkpoints": checkpoints,
        "final_ang": final_ang,
        "tx": tx,
        "ty": ty,
        "g3d": g3d,
        "first": first,
        "angular_gates": gates,
        "groups": groups,
    }


def _lane_vs_benchmark(prefix: str, rho2: np.ndarray, chain: dict, observed: np.ndarray) -> dict:
    fields = {
        "source_rho2": rho2,
        "m10_los_mag": chain["los_mag"],
        "final_angular_centroid_mag": chain["final_ang"]["angular_centroid_mag"],
        "final_angular_spread_rms": chain["final_ang"]["angular_spread_rms"],
        "final_angular_rms_angle_mag": chain["final_ang"]["angular_rms_angle_mag"],
    }
    out = {}
    for name, field in fields.items():
        p, s, n = _corr(field, observed)
        out[f"{prefix}_{name}_vs_observed_pearson"] = p
        out[f"{prefix}_{name}_vs_observed_spearman"] = s
        out[f"{prefix}_{name}_vs_observed_count"] = n
        out[f"{prefix}_{name}_rms"] = _rms(field)
    return out


def _paired_detector(ind: dict, assisted: dict) -> dict:
    half = max(
        float(np.max(np.abs(ind["tx"]))), float(np.max(np.abs(ind["ty"]))),
        float(np.max(np.abs(assisted["tx"]))), float(np.max(np.abs(assisted["ty"]))),
    ) * RANGE_PAD
    if half <= 0.0: raise RuntimeError("paired detector half-range is zero")
    hi = DET._detector_hist(ind["tx"], ind["ty"], DETECTOR_RESOLUTION, half)
    ha = DET._detector_hist(assisted["tx"], assisted["ty"], DETECTOR_RESOLUTION, half)
    # Reuse detector metrics; key names say moment_control internally, so rename.
    m = DET._image_metrics(ha, hi)
    return {
        "detector_resolution": DETECTOR_RESOLUTION,
        "detector_half_range": half,
        "assisted_flux": m["empirical_flux"],
        "independent_flux": m["moment_control_flux"],
        "flux_relative_error": m["flux_relative_error"],
        "total_variation_distance": m["total_variation_distance"],
        "jensen_shannon_divergence_nats": m["jensen_shannon_divergence_nats"],
        "normalized_l2_distance": m["normalized_l2_distance"],
        "image_pearson": m["image_pearson"],
        "image_spearman": m["image_spearman"],
        "assisted_occupied_bins": m["empirical_occupied_bins"],
        "independent_occupied_bins": m["moment_control_occupied_bins"],
        "occupied_bin_intersection": m["occupied_bin_intersection"],
        "assisted_hist": ha,
        "independent_hist": hi,
    }


def _run_cluster(cluster: dict) -> tuple[dict, dict]:
    cid = cluster["id"]

    # ---------------------------------------------------------------
    # INDEPENDENT LANE FIRST. No observed kappa pixel values are loaded.
    # ---------------------------------------------------------------
    independent_source = _independent_source(cluster)
    independent_chain = _run_chain_from_rho3(independent_source["rho3"], observed_for_first_step=None)

    # Only now load benchmark pixel values. From this line onward they are comparison
    # data and the assisted-control source, never input to the already-complete
    # independent lane.
    kpath = _kappa_path(cluster)
    with fits.open(kpath) as hdul:
        kappa_native = np.asarray(hdul[0].data, dtype=np.float64)
    observed = BASE.resample_to_grid(kappa_native, BASE.OBS_BINS, BASE.CFG["extent"])

    # Historical assisted control reconstructed exactly as before.
    assisted_rho2 = BASE.construct_common_proxy(kappa_native, bins=BASE.OBS_BINS, extent=BASE.CFG["extent"])
    assisted_rho3 = construct_rho_3d(assisted_rho2, BASE.NZ, profile=BASE.PROFILE)
    assisted_chain = _run_chain_from_rho3(assisted_rho3, observed_for_first_step=observed)

    indep_metrics = _lane_vs_benchmark("independent", independent_source["rho2"], independent_chain, observed)
    assisted_metrics = _lane_vs_benchmark("assisted", assisted_rho2, assisted_chain, observed)

    # Stage-by-stage direct lane similarity exposes where the independent lane first
    # departs from the assisted diagnostic baseline.
    stage_pairs = {
        "source_rho2": (assisted_rho2, independent_source["rho2"]),
        "m10_los_mag": (assisted_chain["los_mag"], independent_chain["los_mag"]),
        "final_angular_centroid_mag": (assisted_chain["final_ang"]["angular_centroid_mag"], independent_chain["final_ang"]["angular_centroid_mag"]),
        "final_angular_spread_rms": (assisted_chain["final_ang"]["angular_spread_rms"], independent_chain["final_ang"]["angular_spread_rms"]),
        "final_angular_rms_angle_mag": (assisted_chain["final_ang"]["angular_rms_angle_mag"], independent_chain["final_ang"]["angular_rms_angle_mag"]),
    }
    pair_metrics = {}
    for name, (a, i) in stage_pairs.items():
        p, s, n = _corr(a, i)
        pair_metrics[f"assisted_vs_independent_{name}_pearson"] = p
        pair_metrics[f"assisted_vs_independent_{name}_spearman"] = s
        pair_metrics[f"assisted_vs_independent_{name}_count"] = n

    det = _paired_detector(independent_chain, assisted_chain)

    summary = {
        "cluster_id": cid,
        "candidate_id": BASE.CANDIDATE_ID,
        "physical_source_representation": BASE.PHYSICAL_SOURCE,
        "assisted_source_role": "observed_kappa_positive_normalized_proxy_training_wheels_control",
        "independent_source_role": independent_source["source_role"],
        "independent_source_filter": "HST_WFC3_F160W",
        "independent_source_background_rule": independent_source["alignment"]["fixed_background_estimator"],
        "independent_source_observed_kappa_pixel_values_used": independent_source["observed_kappa_pixel_values_used"],
        "benchmark_values_loaded_after_independent_lane_complete": True,
        "benchmark_role": "external_morphology_comparison_only_for_independent_lane",
        "independent_source_limit": "luminous_structure_proxy_not_baryonic_or_total_mass_map",
        "geometry_lane": "G3D_LOS_consistent_diagnostic",
        "observer_readout": "validated_per_source_angular_moments_plus_native_global_angular_histogram",
        "n_photons": int(len(independent_chain["tx"])),
        "source_supported_bins": int(len(independent_chain["groups"])),
        "hst_url": independent_source["hst_url"],
        "hst_sha256": independent_source["hst_sha256"],
        "hst_candidate_count": len(independent_source["hst_discovered_candidates"]),
        **{f"hst_{k}": v for k, v in independent_source["alignment"].items()},
        "independent_g3d_unit_speed_max_error": float(independent_chain["g3d"]["max_unit_speed_error"]),
        "independent_g3d_unit_speed_pass": bool(independent_chain["g3d"]["max_unit_speed_error"] <= UNIT_SPEED_TOL),
        "independent_first_step_exact_max_vector_error": independent_chain["first"]["first_step_exact_max_vector_error"],
        "independent_first_step_exact_pass": independent_chain["first"]["first_step_exact_pass"],
        "assisted_g3d_unit_speed_max_error": float(assisted_chain["g3d"]["max_unit_speed_error"]),
        "assisted_g3d_unit_speed_pass": bool(assisted_chain["g3d"]["max_unit_speed_error"] <= UNIT_SPEED_TOL),
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
        "observed_kappa_reference_only": observed,
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


def _mean(summaries: list[dict], key: str) -> float:
    vals = np.asarray([s[key] for s in summaries], dtype=np.float64)
    return float(np.nanmean(vals)) if np.any(np.isfinite(vals)) else float("nan")


def main() -> int:
    started = time.perf_counter()
    OUT.mkdir(parents=True, exist_ok=True); DOWNLOADS.mkdir(parents=True, exist_ok=True)
    repo = _repo_state(); _write_json(OUT / "repository_state.json", repo)
    if repo["branch"] != "main" or repo["tracked_changes"] or repo["staged_changes"]:
        v = {"lab_id": LAB_ID, "outcome": "REPOSITORY_GATE_FAILURE", "head_sha": repo["head_sha"]}
        _write_json(OUT / "validation.json", v); print(json.dumps(v, indent=2)); return 2

    summaries = []; failures = []
    for cluster in BASE.CLUSTERS:
        cid = cluster["id"]
        print(f"[{cid}] TRAINING WHEELS OFF: HST F160W independent morphology -> frozen PBUF chain")
        try:
            summary, arrays = _run_cluster(cluster)
            summaries.append(summary)
            cdir = OUT / "clusters" / cid; cdir.mkdir(parents=True, exist_ok=True)
            _write_json(cdir / "independent_source_summary.json", summary)
            np.savez_compressed(cdir / "independent_source_fields.npz", **arrays)
        except Exception as exc:
            failures.append({"cluster_id": cid, "error": repr(exc)})
            _write_json(OUT / "cluster_failures.json", failures)
            raise

    _write_csv(OUT / "independent_source_summary.csv", summaries)
    _write_json(OUT / "cluster_failures.json", failures)

    validation = {
        "lab_id": LAB_ID,
        "outcome": "Outcome A — INDEPENDENT SOURCE TRAINING-WHEELS-OFF AUDIT COMPLETE",
        "head_sha": repo["head_sha"],
        "candidate_id": BASE.CANDIDATE_ID,
        "physical_source_representation": BASE.PHYSICAL_SOURCE,
        "assisted_source_role": "observed_kappa_positive_normalized_proxy_training_wheels_control",
        "independent_source_role": "HST_F160W_positive_luminous_structure_proxy_WCS_aligned_no_kappa_pixel_values",
        "independent_source_limit": "luminous_structure_proxy_not_baryonic_or_total_mass_map",
        "benchmark_role": "external_morphology_comparison_only_for_independent_lane",
        "benchmark_values_loaded_after_independent_lane_complete": True,
        "cluster_count_expected": len(BASE.CLUSTERS),
        "cluster_count_completed": len(summaries),
        "all_independent_source_observed_kappa_pixel_values_unused": bool(all(not s["independent_source_observed_kappa_pixel_values_used"] for s in summaries)),
        "all_independent_g3d_unit_speed_pass": bool(all(s["independent_g3d_unit_speed_pass"] for s in summaries)),
        "all_independent_first_step_exact_pass": bool(all(s["independent_first_step_exact_pass"] for s in summaries)),
        "all_assisted_g3d_unit_speed_pass": bool(all(s["assisted_g3d_unit_speed_pass"] for s in summaries)),
        "all_assisted_first_step_exact_pass": bool(all(s["assisted_first_step_exact_pass"] for s in summaries)),
        "mean_assisted_source_rho2_vs_observed_pearson": _mean(summaries, "assisted_source_rho2_vs_observed_pearson"),
        "mean_independent_source_rho2_vs_observed_pearson": _mean(summaries, "independent_source_rho2_vs_observed_pearson"),
        "mean_assisted_m10_los_mag_vs_observed_pearson": _mean(summaries, "assisted_m10_los_mag_vs_observed_pearson"),
        "mean_independent_m10_los_mag_vs_observed_pearson": _mean(summaries, "independent_m10_los_mag_vs_observed_pearson"),
        "mean_assisted_final_angular_rms_angle_mag_vs_observed_pearson": _mean(summaries, "assisted_final_angular_rms_angle_mag_vs_observed_pearson"),
        "mean_independent_final_angular_rms_angle_mag_vs_observed_pearson": _mean(summaries, "independent_final_angular_rms_angle_mag_vs_observed_pearson"),
        "mean_assisted_vs_independent_source_rho2_pearson": _mean(summaries, "assisted_vs_independent_source_rho2_pearson"),
        "mean_assisted_vs_independent_m10_los_mag_pearson": _mean(summaries, "assisted_vs_independent_m10_los_mag_pearson"),
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
    _write_json(OUT / "validation.json", validation)
    _write_json(OUT / "run.json", {"validation": validation, "summaries": summaries})
    print(json.dumps(validation, indent=2, default=_json_default))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
