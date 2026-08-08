#!/usr/bin/env python3
"""PBUF FOUNDATION — CLUSTER DISTANCE / REDSHIFT CLOSURE AUDIT 001.

Fact-finding only.

Purpose
-------
Close, or precisely localize, the first open physical link in the five-cluster
baryonic source chain: observed spectral redshift -> physical distance geometry.

Critical distinction
--------------------
Observed redshift is an observable wave-frequency/wavelength ratio.  This audit
must NOT silently identify it with pure expansion redshift.  We therefore keep

    z_obs != z_expansion   by assumption

until peculiar, local/gravitational and possible medium/propagation contributions
are shown negligible or are independently modelled.

Historical PBUF V11 provenance
------------------------------
The V11 preprint specifies an elastic background

    E(a)^2 = Omega_m0 a^-3 + Omega_r0 a^-4 + Omega_sigma(a)
    H(a)   = H0 E(a)

where Omega_sigma(a) is built from thermal-table alpha_T(a), epsilon0_T(a),
kmax(a), an activation decay, and a flat-today rescaling.  The current weak-
lensing lab repository does not expose a complete audited distance API or the
thermal LUT required to reproduce that background here.  Historical equations
are therefore recorded as provenance, not promoted into a new implementation.

Distance closure requirements
-----------------------------
A reproducible PBUF distance mapping needs, at minimum:
1. a justified mapping from z_obs to z_expansion;
2. an audited PBUF H(z) / E(a) implementation and all of its inputs;
3. a comoving/line-of-sight distance integration rule;
4. an angular-diameter-distance rule;
5. a luminosity-distance rule;
6. validation that the usual distance-duality relation remains valid under the
   PBUF photon/medium propagation law, rather than assuming it automatically.

No kappa, shear, lensing amplitude, G backsolve, stellar M/L, gas fraction,
Quantum Engine, Planck-scale input, or fitted redshift correction is allowed.
"""
from __future__ import annotations

import json
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
LAB_ID = "PBUF-FOUNDATION-CLUSTER-DISTANCE-REDSHIFT-CLOSURE-AUDIT-001"

CLUSTERS = (
    ("Abell2744", "Abell 2744", 0.308),
    ("MACS0416", "MACSJ0416.1-2403", 0.396),
    ("MACS1149", "MACSJ1149.5+2223", 0.543),
    ("AbellS1063", "Abell S1063 (RXCJ2248.7-4431)", 0.348),
    ("Abell370", "Abell 370", 0.375),
)

# Structural repository search terms only. Presence is not sufficient for closure.
SEARCH_TERMS = {
    "pbuf_background_expansion": ("Omega_sigma", "alpha_T", "epsilon0_T", "kmax"),
    "hubble_function": ("H_of_z", "H_of_a", "E_of_z", "E_of_a", "hubble"),
    "comoving_distance": ("comoving_distance", "line_of_sight_distance", "chi_of_z"),
    "angular_diameter_distance": ("angular_diameter_distance", "D_A", "DA_of_z"),
    "luminosity_distance": ("luminosity_distance", "D_L", "DL_of_z"),
    "distance_duality": ("distance_duality", "etherington", "reciprocity_relation"),
    "redshift_decomposition": ("z_expansion", "z_peculiar", "z_gravitational", "z_medium"),
}

SKIP_PARTS = {".git", "runs", "__pycache__", ".venv", "venv"}
TEXT_SUFFIXES = {".py", ".md", ".txt", ".json", ".toml", ".yaml", ".yml"}


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


def _iter_text_files():
    for p in ROOT.rglob("*"):
        if not p.is_file() or p.suffix.lower() not in TEXT_SUFFIXES:
            continue
        if any(part in SKIP_PARTS for part in p.parts):
            continue
        if p.name == Path(__file__).name:
            continue
        yield p


def _repo_inventory() -> dict:
    hits = {k: [] for k in SEARCH_TERMS}
    files = list(_iter_text_files())
    for p in files:
        try:
            text = p.read_text(errors="ignore")
        except Exception:
            continue
        rel = str(p.relative_to(ROOT))
        for key, terms in SEARCH_TERMS.items():
            found = [t for t in terms if t in text]
            if found:
                hits[key].append({"path": rel, "terms": found})
    return hits


def _closed_by_repo_implementation(hits: dict, key: str) -> bool:
    """Conservative closure rule: marker presence alone never closes physics."""
    # This audit is intentionally conservative.  A later lab may whitelist a
    # verified implementation once its inputs and equations have been audited.
    return False


def _cluster_row(cid: str, name: str, z_obs: float, hits: dict) -> dict:
    links = [
        {
            "key": "observed_spectral_redshift",
            "status": "OBSERVED_AVAILABLE",
            "value": z_obs,
            "role": "measured total wavelength/frequency shift",
            "blocks_distance_closure": False,
        },
        {
            "key": "redshift_decomposition",
            "status": "PHYSICAL_DECOMPOSITION_OPEN",
            "role": "separate expansion, peculiar, local/gravitational and possible PBUF-medium propagation contributions",
            "blocks_distance_closure": True,
            "guardrail": "z_obs is not identified with z_expansion by assumption",
        },
        {
            "key": "pbuf_background_Hz",
            "status": "HISTORICAL_V11_EQUATION_PRESENT_CURRENT_WL_IMPLEMENTATION_NOT_AUDITED",
            "role": "provide E(a) and H(a) from the elastic background with all thermal/LUT inputs",
            "blocks_distance_closure": True,
        },
        {
            "key": "comoving_distance_rule",
            "status": "NOT_CLOSED_IN_CURRENT_WL_FOUNDATION",
            "role": "integrate the audited PBUF expansion history into a physical radial distance",
            "blocks_distance_closure": True,
        },
        {
            "key": "angular_diameter_distance_rule",
            "status": "NOT_CLOSED_IN_CURRENT_WL_FOUNDATION",
            "role": "map observed angular separations/pixel scales to transverse physical size",
            "blocks_distance_closure": True,
        },
        {
            "key": "luminosity_distance_rule",
            "status": "NOT_CLOSED_IN_CURRENT_WL_FOUNDATION",
            "role": "map calibrated flux to luminosity",
            "blocks_distance_closure": True,
        },
        {
            "key": "distance_duality_under_pbuf_propagation",
            "status": "REQUIRES_PROPAGATION_VALIDATION",
            "role": "test whether D_L=(1+z)^2 D_A remains valid if the medium changes wave frequency/amplitude",
            "blocks_distance_closure": True,
        },
    ]
    return {
        "cluster_id": cid,
        "cluster_name": name,
        "z_observed": z_obs,
        "z_expansion": None,
        "D_comoving_Mpc": None,
        "D_A_Mpc": None,
        "D_L_Mpc": None,
        "kpc_per_arcsec": None,
        "distance_geometry_closed": False,
        "links": links,
    }


def main() -> None:
    before = _repo_state()
    hits = _repo_inventory()
    rows = [_cluster_row(*c, hits) for c in CLUSTERS]
    after = _repo_state()

    historical = {
        "role": "provenance_only_not_reimplemented_here",
        "v11_background_equation": "E(a)^2=Omega_m0*a^-3+Omega_r0*a^-4+Omega_sigma(a); H(a)=H0*E(a)",
        "v11_elastic_inputs": ["alpha_T(a)", "epsilon0_T(a)", "kmax(a)", "Rmax", "flat_today_rescale"],
        "reason_not_promoted": "current WL foundation lacks an audited complete distance implementation and thermal/LUT provenance required to reproduce V11 background exactly",
    }

    blocking = [
        "redshift_decomposition",
        "audited_current_pbuf_Hz",
        "comoving_distance_rule",
        "angular_diameter_distance_rule",
        "luminosity_distance_rule",
        "distance_duality_under_pbuf_propagation",
    ]

    checks = {
        "all_five_observed_redshifts_present": len(rows) == 5 and all(r["z_observed"] > 0 for r in rows),
        "z_observed_not_promoted_to_z_expansion": all(r["z_expansion"] is None for r in rows),
        "no_distance_values_fabricated": all(r["D_A_Mpc"] is None and r["D_L_Mpc"] is None for r in rows),
        "historical_v11_background_recorded_as_provenance_only": True,
        "distance_duality_not_silently_assumed": True,
        "no_lcdm_distance_imported": True,
        "no_kappa_or_lensing_target_used": True,
        "no_G_backsolve": True,
        "no_fitted_redshift_correction": True,
        "no_quantum_engine": True,
        "no_planck_scale": True,
        "gravity_fundamental_in_PBUF": False,
        "no_tracked_or_staged_changes": after["tracked_changes"] == "" and after["staged_changes"] == "",
        "stdout_only_no_run_directory_created": True,
    }

    result = {
        "lab_id": LAB_ID,
        "status": "FACT_FINDING_ONLY",
        "repo_state": after,
        "policy": {
            "redshift_role": "OBSERVED_TOTAL_SPECTRAL_SHIFT_NOT_ASSUMED_PURE_EXPANSION",
            "lcdm_distance_role": "NOT_IMPORTED",
            "historical_v11_role": "PROVENANCE_ONLY_UNTIL_CURRENT_IMPLEMENTATION_AUDITED",
            "lensing_target_used": False,
            "fit_or_tuning_used": False,
            "gravity_fundamental_in_PBUF": False,
        },
        "repository_inventory": hits,
        "historical_pbuf_distance_provenance": historical,
        "clusters": rows,
        "closure": {
            "status": "DISTANCE_REDSHIFT_GEOMETRY_NOT_YET_CLOSED",
            "blocking_links": blocking,
            "safe_next": (
                "recover or rebuild the exact current PBUF background H(a) with its thermal/LUT inputs, then audit the distance integral and distance-duality relation; "
                "keep z_obs separate from z_expansion until medium/local contributions are quantified"
            ),
        },
        "checks": checks,
    }

    print(LAB_ID)
    print("status=FACT_FINDING_ONLY")
    print(f"head_sha={after['head_sha']}")
    print("gravity_fundamental_in_PBUF=false")
    print("z_observed_assumed_pure_expansion=false")
    print("lcdm_distance_imported=false")
    print("lensing_target_used=false")
    print("fit_or_tuning_used=false")
    print()
    print("CLUSTER_DISTANCE_REDSHIFT_STATUS")
    print("cluster | z_obs | z_expansion | PBUF_Hz | D_A | D_L | distance_duality | closed")
    for r in rows:
        print(
            f"{r['cluster_id']} | {r['z_observed']:.3f} | UNRESOLVED | "
            "HISTORICAL_ONLY_CURRENT_WL_NOT_AUDITED | UNRESOLVED | UNRESOLVED | REQUIRES_VALIDATION | False"
        )
    print()
    print("REPOSITORY_IMPLEMENTATION_INVENTORY")
    for key in SEARCH_TERMS:
        print(f"{key}_marker_hits={len(hits[key])}")
    print()
    print("HISTORICAL_V11_PROVENANCE")
    print(historical["v11_background_equation"])
    print("thermal_inputs=" + ",".join(historical["v11_elastic_inputs"]))
    print("promoted_to_current_distance_law=false")
    print()
    print("OPEN_PHYSICS")
    print("z_obs=(expansion)*(peculiar)*(local/gravitational)*(possible_medium_propagation) schematically; components not solved here")
    print("distance_duality_D_L_equals_(1+z)^2_D_A_assumed=false")
    print()
    print("CONCLUSION")
    print("status=DISTANCE_REDSHIFT_GEOMETRY_NOT_YET_CLOSED")
    print("blocking_links=" + ",".join(blocking))
    print("safe_next=" + result["closure"]["safe_next"])
    print()
    print("CHECKS")
    for k, v in checks.items():
        print(f"{k}={str(v).lower()}")
    print("JSON=" + json.dumps(result, sort_keys=True, separators=(",", ":")))


if __name__ == "__main__":
    main()
