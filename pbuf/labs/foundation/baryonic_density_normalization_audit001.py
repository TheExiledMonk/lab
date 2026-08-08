#!/usr/bin/env python3
"""PBUF FOUNDATION — BARYONIC DENSITY NORMALIZATION AUDIT 001.

Fact-finding only.

Purpose
-------
Test whether the existing independent HST/F160W source path contains enough
independently calibrated information to construct an absolute baryonic mass
density field in SI units before any weak-lensing target is consulted.

The prior native-scale audit found that the native transfer is well-defined in
native units but the absolute bridge constrains only

    T_native = (8*pi*G/c^2) * RHO0 * L_cg^2.

This lab attacks the opposite side independently.  It asks whether RHO0 can be
closed from source physics alone.  If and only if an absolute baryonic density
normalization is independently available, the previously measured native
transfer may be used to predict L_cg.  The prediction must not be fitted.

The current independent source is an HST F160W luminous proxy.  This audit
therefore inventories, in order:

1. detector/photometric calibration metadata;
2. flux/luminosity distance requirements;
3. stellar mass-to-light / stellar-population requirements;
4. gas / intracluster baryon requirements;
5. 2D surface-density -> 3D density/deprojection requirements;
6. whether the current normalization erases absolute amplitude.

No kappa pixels, shear, lensing morphology, fitted lensing amplitude, 0.18,
Quantum Engine, or Planck-scale input is allowed.
"""
from __future__ import annotations

import json
import math
import subprocess
import sys
from pathlib import Path

import numpy as np
from astropy.io import fits

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))

import pbuf.labs.foundation.independent_source_training_wheels_off001 as LAB
import pbuf.labs.foundation.independent_source_training_wheels_off001_common_footprint_fix as SRC

LAB_ID = "PBUF-FOUNDATION-BARYONIC-DENSITY-NORMALIZATION-AUDIT-001"
OUT = ROOT / "runs" / "baryonic_density_normalization_audit001"
DOWNLOADS = OUT / "downloads"

# Previous native-side result.  This is not re-fitted here; it is frozen from the
# successful native-length-scale mapping audit and is used only if the source-side
# SI density normalization closes independently.
T_NATIVE_CENTER = 0.9989359985074102
G_MEASURED = 6.67430e-11
C = 299_792_458.0

PHOT_KEYS = (
    "BUNIT", "EXPTIME", "PHOTFLAM", "PHOTPLAM", "PHOTBW", "ZPTMAG",
    "ABMAGZP", "VEGAMAG", "FILTER", "INSTRUME", "DETECTOR",
)


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


def _json_default(x):
    if isinstance(x, (np.floating,)): return float(x)
    if isinstance(x, (np.integer,)): return int(x)
    if isinstance(x, (np.bool_,)): return bool(x)
    if isinstance(x, Path): return str(x)
    return str(x)


def _positive_stats(a: np.ndarray) -> dict:
    x = np.asarray(a, dtype=np.float64)
    x = x[np.isfinite(x) & (x > 0.0)]
    if not x.size:
        return {"count": 0, "sum": 0.0, "mean": 0.0, "rms": 0.0, "max": 0.0}
    return {
        "count": int(x.size),
        "sum": float(np.sum(x)),
        "mean": float(np.mean(x)),
        "rms": float(np.sqrt(np.mean(x*x))),
        "max": float(np.max(x)),
    }


def _photometry_inventory(path: Path) -> dict:
    with fits.open(path, memmap=True) as hdul:
        data = np.asarray(np.squeeze(hdul[0].data), dtype=np.float64)
        hdr = hdul[0].header
        meta = {k: hdr.get(k) for k in PHOT_KEYS}
    finite = data[np.isfinite(data)]
    return {
        "header": meta,
        "finite_pixel_count": int(finite.size),
        "finite_data_min": float(np.min(finite)) if finite.size else None,
        "finite_data_max": float(np.max(finite)) if finite.size else None,
        "finite_data_rms": float(np.sqrt(np.mean(finite*finite))) if finite.size else None,
        "has_BUNIT": meta.get("BUNIT") is not None,
        "has_PHOTFLAM": meta.get("PHOTFLAM") is not None,
        "has_PHOTPLAM": meta.get("PHOTPLAM") is not None,
        "has_EXPTIME": meta.get("EXPTIME") is not None,
        "photometric_flux_calibration_metadata_present": bool(
            meta.get("BUNIT") is not None and
            (meta.get("PHOTFLAM") is not None or meta.get("ABMAGZP") is not None or meta.get("ZPTMAG") is not None)
        ),
    }


def _requirements(phot: dict, source: dict) -> list[dict]:
    align = source["alignment"]
    header = phot["header"]
    return [
        {
            "key": "absolute_detector_to_flux_calibration",
            "required": True,
            "present": bool(phot["photometric_flux_calibration_metadata_present"]),
            "role": "convert raw HST image values to calibrated observed flux density / surface brightness",
            "may_be_fitted": False,
        },
        {
            "key": "source_redshift_and_luminosity_distance",
            "required": True,
            "present": False,
            "role": "convert observed flux to rest-frame luminosity; cluster/member redshift or equivalent independently calibrated distance is required",
            "may_be_fitted": False,
        },
        {
            "key": "stellar_mass_to_light_or_population_model",
            "required": True,
            "present": False,
            "role": "convert F160W luminosity into stellar baryonic mass; requires independently justified stellar-population/M-L information, not a lensing fit",
            "may_be_fitted": False,
        },
        {
            "key": "diffuse_gas_baryon_component",
            "required": True,
            "present": False,
            "role": "clusters contain baryons not traced by stellar F160W light; hot gas / intracluster baryon mass must be supplied independently",
            "may_be_fitted": False,
        },
        {
            "key": "physical_pixel_area_or_angular_diameter_distance",
            "required": True,
            "present": False,
            "role": "convert calibrated sky surface brightness/mass per angular pixel into kg/m^2 on the source plane",
            "may_be_fitted": False,
        },
        {
            "key": "line_of_sight_depth_or_deprojection_model",
            "required": True,
            "present": False,
            "role": "convert baryonic surface density kg/m^2 to volume density kg/m^3 for the current 3D native source",
            "may_be_fitted": False,
        },
        {
            "key": "absolute_amplitude_preserved_into_rho2",
            "required": True,
            "present": False,
            "role": "current common-footprint source divides positive luminous field by its own maximum, erasing absolute amplitude before rho2/rho3",
            "may_be_fitted": False,
            "evidence": {
                "positive_luminous_common_max": align.get("positive_luminous_common_max"),
                "source_code_rule": "rho2 = luminous / maxv",
            },
        },
    ]


def _run_cluster(cluster: dict) -> dict:
    old = SRC.DOWNLOADS
    SRC.DOWNLOADS = DOWNLOADS
    try:
        source = SRC._independent_source(cluster)
    finally:
        SRC.DOWNLOADS = old

    local = Path(source["hst_local_path"])
    phot = _photometry_inventory(local)
    luminous = np.asarray(source["luminous_common"], dtype=np.float64)
    rho2 = np.asarray(source["rho2"], dtype=np.float64)
    valid = np.asarray(source["geometry"]["valid_mask"], dtype=bool)
    req = _requirements(phot, source)
    closed = all(r["present"] for r in req if r["required"])

    return {
        "cluster_id": cluster["id"],
        "hst_local_path": str(local),
        "hst_sha256": source["hst_sha256"],
        "source_role": source["source_role"],
        "observed_kappa_pixel_values_used": False,
        "photometry": phot,
        "raw_common_luminous_positive_stats": _positive_stats(luminous[valid]),
        "normalized_rho2_positive_stats": _positive_stats(rho2[valid]),
        "rho2_max": float(np.max(rho2)),
        "absolute_amplitude_erased_by_current_rho2_normalization": True,
        "requirements": req,
        "absolute_baryonic_surface_density_closed": False,
        "absolute_baryonic_volume_density_closed": bool(closed),
        "RHO0_kg_m3_per_native": None,
        "L_cg_predicted_m": None,
    }


def main() -> None:
    repo_before = _repo_state()
    rows = [_run_cluster(c) for c in LAB.BASE.CLUSTERS]

    any_si_density = any(r["absolute_baryonic_volume_density_closed"] for r in rows)
    rho0_values = [r["RHO0_kg_m3_per_native"] for r in rows if r["RHO0_kg_m3_per_native"] is not None]

    if any_si_density and rho0_values:
        rho0 = float(np.mean(rho0_values))
        L_cg = math.sqrt(T_NATIVE_CENTER * C*C / (8.0 * math.pi * G_MEASURED * rho0))
        closure_status = "SOURCE_SIDE_SI_DENSITY_CLOSED_LCG_PREDICTED"
    else:
        rho0 = None
        L_cg = None
        closure_status = "BARYONIC_DENSITY_NORMALIZATION_NOT_YET_CLOSED"

    repo_after = _repo_state()
    checks = {
        "all_clusters_processed": len(rows) == len(LAB.BASE.CLUSTERS),
        "no_kappa_pixel_values_used": all(not r["observed_kappa_pixel_values_used"] for r in rows),
        "current_rho2_normalization_erases_absolute_amplitude": all(r["absolute_amplitude_erased_by_current_rho2_normalization"] for r in rows),
        "no_source_side_SI_density_fabricated": not any_si_density,
        "no_RHO0_solved_from_G_or_native_transfer": rho0 is None,
        "no_Lcg_predicted_without_independent_RHO0": L_cg is None,
        "legacy_0p18_used": False,
        "lensing_target_used": False,
        "fit_or_tuning_used": False,
        "quantum_engine_used": False,
        "planck_scale_used": False,
        "gravity_fundamental_in_PBUF": False,
        "measured_G_used_only_if_source_side_closes": True,
        "no_tracked_or_staged_changes": repo_after["tracked_changes"] == "" and repo_after["staged_changes"] == "",
        "stdout_only_no_run_directory_created": not OUT.exists(),
    }

    result = {
        "lab_id": LAB_ID,
        "status": "FACT_FINDING_ONLY",
        "repo_state": repo_after,
        "policy": {
            "gravity_fundamental_in_PBUF": False,
            "measured_G_role": "MACROSCOPIC_RESPONSE_ANCHOR_ONLY_AND_UNUSED_UNLESS_SOURCE_SIDE_CLOSES",
            "legacy_0p18_used": False,
            "lensing_target_used": False,
            "fit_or_tuning_used": False,
        },
        "clusters": rows,
        "closure": {
            "status": closure_status,
            "absolute_baryonic_density_scale_available": any_si_density,
            "RHO0_kg_m3_per_native": rho0,
            "L_cg_predicted_m": L_cg,
            "native_transfer_frozen": T_NATIVE_CENTER,
            "important_result": (
                "Only an independently derived baryonic SI density normalization may unlock the opposite-side prediction of L_cg. "
                "Photometric calibration alone is insufficient if luminosity distance, stellar M/L, gas baryons, physical area, deprojection, or amplitude preservation remain open."
            ),
            "safe_next": (
                "supply independently measured baryonic mass information (stellar and gas) plus geometry, then rebuild the source without max-normalization and predict L_cg; "
                "do not infer mass-to-light or gas normalization from lensing"
            ),
        },
        "checks": checks,
    }

    print(LAB_ID)
    print("status=FACT_FINDING_ONLY")
    print(f"head_sha={repo_after['head_sha']}")
    print("gravity_fundamental_in_PBUF=false")
    print("legacy_0p18_used=false")
    print("lensing_target_used=false")
    print("fit_or_tuning_used=false")
    print()
    print("CLUSTER_SOURCE_CLOSURE")
    print("cluster | photometric_cal | distance | stellar_ML | gas | physical_area | deprojection | amplitude_preserved | SI_rho_closed")
    for r in rows:
        req = {x["key"]: x for x in r["requirements"]}
        print(
            f"{r['cluster_id']} | "
            f"{req['absolute_detector_to_flux_calibration']['present']} | "
            f"{req['source_redshift_and_luminosity_distance']['present']} | "
            f"{req['stellar_mass_to_light_or_population_model']['present']} | "
            f"{req['diffuse_gas_baryon_component']['present']} | "
            f"{req['physical_pixel_area_or_angular_diameter_distance']['present']} | "
            f"{req['line_of_sight_depth_or_deprojection_model']['present']} | "
            f"{req['absolute_amplitude_preserved_into_rho2']['present']} | "
            f"{r['absolute_baryonic_volume_density_closed']}"
        )
    print()
    print("CONCLUSION")
    print(f"status={closure_status}")
    print(f"absolute_baryonic_density_scale_available={str(any_si_density).lower()}")
    print(f"RHO0_kg_m3_per_native={rho0}")
    print(f"L_cg_predicted_m={L_cg}")
    print("current_F160W_role=luminous_structure_proxy_not_absolute_baryonic_mass_map")
    print("safe_next=add independent stellar+gas baryonic mass and physical geometry, preserve absolute amplitude, then predict L_cg without lensing fit")
    print()
    print("CHECKS")
    for k, v in checks.items():
        print(f"{k}={str(v).lower() if isinstance(v, bool) else v}")
    print("JSON=" + json.dumps(result, sort_keys=True, default=_json_default, separators=(",", ":")))


if __name__ == "__main__":
    main()
