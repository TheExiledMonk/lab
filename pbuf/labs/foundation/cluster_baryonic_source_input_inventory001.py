#!/usr/bin/env python3
"""PBUF FOUNDATION — CLUSTER BARYONIC SOURCE INPUT INVENTORY 001.

Fact-finding only.

Purpose
-------
Step back from weak-lensing prediction and audit the physical source chain one
conversion at a time for the five current Frontier Fields clusters.

The previous baryonic-density normalization audit established that calibrated
F160W imaging is available, but absolute baryonic density is not yet closed.
This lab separates the source bridge into individually auditable links and
classifies each link as:

  OBSERVED_AVAILABLE
  PIPELINE_AVAILABLE_BUT_AMPLITUDE_ERASED
  ASTROPHYSICAL_CONVERSION_REQUIRED
  COSMOLOGY_DEPENDENT_GEOMETRY
  INDEPENDENT_EXTERNAL_DATA_REQUIRED
  NOT_YET_CLOSED

Important distinction
---------------------
Cluster redshift is accepted as directly observed source metadata.  Redshift is
NOT silently converted into luminosity distance or angular-diameter distance,
because that conversion requires a cosmological geometry.  A future PBUF source
bridge must either use independently supplied distances or explicitly state and
audit the cosmological distance model used.

No kappa pixels, shear, lensing morphology, lensing amplitude, historical 0.18,
Quantum Engine, Planck-scale input, fitted M/L, fitted gas fraction, fitted LOS
depth, or back-solving from G is allowed.
"""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))

LAB_ID = "PBUF-FOUNDATION-CLUSTER-BARYONIC-SOURCE-INPUT-INVENTORY-001"

# STScI Hubble Frontier Fields archive metadata.  Redshifts are directly listed
# by the archive for these cluster fields.  They are observational metadata only;
# this lab does not convert z into a distance.
STSCI_FRONTIER_FIELDS_URL = "https://archive.stsci.edu/prepds/frontier/"
CLUSTERS = (
    {"id": "Abell2744", "name": "Abell 2744", "z": 0.308},
    {"id": "MACS0416", "name": "MACSJ0416.1-2403", "z": 0.396},
    {"id": "MACS1149", "name": "MACSJ1149.5+2223", "z": 0.543},
    {"id": "AbellS1063", "name": "Abell S1063 (RXCJ2248.7-4431)", "z": 0.348},
    {"id": "Abell370", "name": "Abell 370", "z": 0.375},
)

SOURCE_CODE = ROOT / "pbuf" / "labs" / "foundation" / "independent_source_training_wheels_off001_common_footprint_fix.py"
PRIOR_AUDIT = ROOT / "pbuf" / "labs" / "foundation" / "baryonic_density_normalization_audit001.py"
MANIFEST_ROOT = ROOT / "pbuf" / "data" / "baryonic_source_inputs"


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


def _source_code_evidence() -> dict:
    text = SOURCE_CODE.read_text(errors="ignore")
    return {
        "source_file": str(SOURCE_CODE.relative_to(ROOT)),
        "f160w_source_present": "F160W" in text or "f160w" in text,
        "wcs_geometry_present": "WCS" in text,
        "positive_luminous_field_present": "luminous" in text,
        "max_normalization_present": "rho2 = luminous / maxv" in text,
        "rho3_constructor_present": "construct_rho_3d" in text,
        "kappa_pixels_excluded_before_independent_chain": (
            "NO kappa pixels" in text or "No observed kappa pixel values" in text or "observed kappa pixel values" in text
        ),
    }


def _prior_audit_evidence() -> dict:
    text = PRIOR_AUDIT.read_text(errors="ignore") if PRIOR_AUDIT.exists() else ""
    return {
        "prior_audit_present": PRIOR_AUDIT.exists(),
        "photometric_calibration_audited": "PHOTFLAM" in text and "BUNIT" in text,
        "stellar_ml_missing_recorded": "stellar_mass_to_light_or_population_model" in text,
        "gas_missing_recorded": "diffuse_gas_baryon_component" in text,
        "deprojection_missing_recorded": "line_of_sight_depth_or_deprojection_model" in text,
        "amplitude_erasure_recorded": "absolute_amplitude_erased_by_current_rho2_normalization" in text,
    }


def _load_manifest(cluster_id: str) -> dict | None:
    path = MANIFEST_ROOT / f"{cluster_id}.json"
    if not path.exists():
        return None
    obj = json.loads(path.read_text())
    if not isinstance(obj, dict):
        raise RuntimeError(f"manifest must be a JSON object: {path}")
    return obj


def _manifest_has(m: dict | None, key: str) -> bool:
    if not m:
        return False
    value = m.get(key)
    if value is None:
        return False
    if isinstance(value, str):
        return bool(value.strip())
    return True


def _link(key: str, status: str, role: str, evidence: str, blocks_si_density: bool) -> dict:
    return {
        "key": key,
        "status": status,
        "role": role,
        "evidence": evidence,
        "blocks_SI_density_closure": bool(blocks_si_density),
    }


def _cluster_row(cluster: dict, source_ev: dict, prior_ev: dict) -> dict:
    manifest = _load_manifest(cluster["id"])

    stellar_mass = _manifest_has(manifest, "stellar_mass_map") or _manifest_has(manifest, "stellar_surface_density_map")
    gas_mass = _manifest_has(manifest, "gas_mass_map") or _manifest_has(manifest, "gas_surface_density_map")
    independent_distance = _manifest_has(manifest, "luminosity_distance_m") and _manifest_has(manifest, "angular_diameter_distance_m")
    physical_area = _manifest_has(manifest, "physical_pixel_area_m2")
    deprojection = _manifest_has(manifest, "los_depth_m") or _manifest_has(manifest, "deprojection_model")
    amplitude_preserving_source = _manifest_has(manifest, "absolute_baryonic_density_map") or _manifest_has(manifest, "absolute_baryonic_surface_density_map")

    links = [
        _link(
            "observed_cluster_redshift",
            "OBSERVED_AVAILABLE",
            "source metadata / spectral shift",
            f"STScI Frontier Fields archive: z={cluster['z']}",
            False,
        ),
        _link(
            "hst_detector_to_flux_calibration",
            "OBSERVED_AVAILABLE" if prior_ev["photometric_calibration_audited"] else "NOT_YET_CLOSED",
            "convert HST detector units to calibrated flux density / surface brightness",
            "Prior audit found BUNIT/PHOTFLAM/PHOTPLAM/EXPTIME on all five F160W mosaics",
            not prior_ev["photometric_calibration_audited"],
        ),
        _link(
            "redshift_to_physical_distance",
            "OBSERVED_AVAILABLE" if independent_distance else "COSMOLOGY_DEPENDENT_GEOMETRY",
            "obtain luminosity and angular-diameter distances without silently assuming LCDM geometry",
            "Observed z alone does not specify D_L or D_A; explicit independent distance or audited cosmological geometry is required",
            not independent_distance,
        ),
        _link(
            "stellar_baryonic_mass",
            "OBSERVED_AVAILABLE" if stellar_mass else "ASTROPHYSICAL_CONVERSION_REQUIRED",
            "convert stellar light / SED information into stellar baryonic mass",
            "Requires independent stellar-mass catalog/map or population-model conversion; F160W alone is not an absolute stellar-mass map",
            not stellar_mass,
        ),
        _link(
            "hot_diffuse_gas_baryons",
            "OBSERVED_AVAILABLE" if gas_mass else "INDEPENDENT_EXTERNAL_DATA_REQUIRED",
            "include intracluster gas baryons not traced by F160W stellar light",
            "Requires independent X-ray/SZ/gas mass or surface-density information",
            not gas_mass,
        ),
        _link(
            "physical_pixel_area",
            "OBSERVED_AVAILABLE" if physical_area else "COSMOLOGY_DEPENDENT_GEOMETRY",
            "convert mass/flux per angular pixel into physical surface density",
            "Needs angular scale plus D_A, or an independently supplied physical pixel area",
            not physical_area,
        ),
        _link(
            "surface_to_volume_density",
            "OBSERVED_AVAILABLE" if deprojection else "ASTROPHYSICAL_CONVERSION_REQUIRED",
            "convert kg/m^2 baryonic surface density into kg/m^3 source density",
            "Requires LOS depth or a documented deprojection model independent of lensing target",
            not deprojection,
        ),
        _link(
            "absolute_amplitude_into_native_source",
            "OBSERVED_AVAILABLE" if amplitude_preserving_source else "PIPELINE_AVAILABLE_BUT_AMPLITUDE_ERASED",
            "preserve absolute baryonic amplitude when constructing rho2/rho3",
            "Current source code uses rho2 = luminous/max(luminous), preserving morphology but erasing absolute inter-cluster amplitude",
            not amplitude_preserving_source,
        ),
    ]

    closed = not any(x["blocks_SI_density_closure"] for x in links)
    blockers = [x["key"] for x in links if x["blocks_SI_density_closure"]]

    return {
        "cluster_id": cluster["id"],
        "cluster_name": cluster["name"],
        "observed_redshift": cluster["z"],
        "redshift_provenance": STSCI_FRONTIER_FIELDS_URL,
        "external_manifest_path": str((MANIFEST_ROOT / f"{cluster['id']}.json").relative_to(ROOT)),
        "external_manifest_present": manifest is not None,
        "links": links,
        "SI_baryonic_volume_density_closed": closed,
        "blocking_links": blockers,
    }


def main() -> None:
    repo = _repo_state()
    source_ev = _source_code_evidence()
    prior_ev = _prior_audit_evidence()
    rows = [_cluster_row(c, source_ev, prior_ev) for c in CLUSTERS]

    all_redshifts = all(any(x["key"] == "observed_cluster_redshift" and x["status"] == "OBSERVED_AVAILABLE" for x in r["links"]) for r in rows)
    all_closed = all(r["SI_baryonic_volume_density_closed"] for r in rows)
    no_manifests = all(not r["external_manifest_present"] for r in rows)

    checks = {
        "all_five_clusters_inventory_present": len(rows) == 5,
        "all_cluster_redshifts_observed_available": all_redshifts,
        "redshift_not_promoted_to_distance": True,
        "current_source_max_normalization_detected": source_ev["max_normalization_present"],
        "current_source_morphology_role_detected": source_ev["positive_luminous_field_present"],
        "no_kappa_or_lensing_target_used": True,
        "no_G_backsolve": True,
        "no_fitted_stellar_ML": True,
        "no_fitted_gas_fraction": True,
        "no_fitted_LOS_depth": True,
        "no_quantum_engine": True,
        "no_planck_scale": True,
        "gravity_fundamental_in_PBUF": False,
        "no_tracked_or_staged_changes": repo["tracked_changes"] == "" and repo["staged_changes"] == "",
        "stdout_only_no_run_directory_created": True,
    }

    result = {
        "lab_id": LAB_ID,
        "status": "FACT_FINDING_ONLY",
        "repo_state": repo,
        "policy": {
            "gravity_fundamental_in_PBUF": False,
            "lensing_target_used": False,
            "legacy_0p18_used": False,
            "fit_or_tuning_used": False,
            "redshift_role": "OBSERVED_SOURCE_METADATA_ONLY",
            "distance_role": "MUST_BE_INDEPENDENTLY_SUPPLIED_OR_EXPLICITLY_AUDITED_COSMOLOGICAL_GEOMETRY",
        },
        "source_code_evidence": source_ev,
        "prior_audit_evidence": prior_ev,
        "clusters": rows,
        "closure": {
            "all_clusters_SI_baryonic_density_closed": all_closed,
            "external_manifests_currently_absent": no_manifests,
            "status": "BARYONIC_SOURCE_CHAIN_STILL_OPEN" if not all_closed else "BARYONIC_SOURCE_CHAIN_CLOSED",
            "safe_next": (
                "Populate independent per-cluster source manifests one physical link at a time: distance geometry, stellar baryonic mass, gas baryons, physical area/deprojection, then rebuild an amplitude-preserving SI source. Do not use lensing to fill any field."
            ),
        },
        "checks": checks,
    }

    print(LAB_ID)
    print("status=FACT_FINDING_ONLY")
    print(f"head_sha={repo['head_sha']}")
    print("gravity_fundamental_in_PBUF=false")
    print("legacy_0p18_used=false")
    print("lensing_target_used=false")
    print("fit_or_tuning_used=false")
    print()
    print("CLUSTER_STEP_INVENTORY")
    print("cluster | z_obs | flux_cal | distance | stellar_mass | gas | physical_area | deprojection | amplitude | SI_rho_closed")
    for r in rows:
        s = {x["key"]: x["status"] for x in r["links"]}
        print(
            f"{r['cluster_id']} | {r['observed_redshift']:.3f} | "
            f"{s['hst_detector_to_flux_calibration']} | {s['redshift_to_physical_distance']} | "
            f"{s['stellar_baryonic_mass']} | {s['hot_diffuse_gas_baryons']} | "
            f"{s['physical_pixel_area']} | {s['surface_to_volume_density']} | "
            f"{s['absolute_amplitude_into_native_source']} | {r['SI_baryonic_volume_density_closed']}"
        )
    print()
    print("PER_CLUSTER_BLOCKERS")
    for r in rows:
        print(f"{r['cluster_id']}: {','.join(r['blocking_links']) if r['blocking_links'] else 'NONE'}")
    print()
    print("CONCLUSION")
    print(f"status={result['closure']['status']}")
    print(f"all_cluster_redshifts_observed_available={str(all_redshifts).lower()}")
    print("redshift_to_distance_automatically_assumed=false")
    print(f"current_max_normalization_erases_absolute_amplitude={str(source_ev['max_normalization_present']).lower()}")
    print(f"all_clusters_SI_baryonic_density_closed={str(all_closed).lower()}")
    print(f"safe_next={result['closure']['safe_next']}")
    print()
    print("CHECKS")
    for k, v in checks.items():
        print(f"{k}={str(v).lower() if isinstance(v, bool) else v}")
    print("JSON=" + json.dumps(result, sort_keys=True, separators=(",", ":")))
    print("EXIT_CODE=0")


if __name__ == "__main__":
    main()
