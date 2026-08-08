#!/usr/bin/env python3
"""PBUF foundation audit: homogeneous elastic-energy closure after Rmax retirement.

Mission
-------
Determine whether the current repository already contains an independently
justified, Rmax-free constitutive law that maps retained thermal/medium inputs

    alpha_T(a), epsilon0_T(a), k_max(a)=epsilon0_T(a)-alpha_T(a)

into a homogeneous elastic energy density rho_sigma(a), and from there into an
elastic cosmological contribution Omega_sigma(a).

This is fact-finding only.  It does not invent a medium action, choose a modulus,
identify alpha_T or k_max with Omega_sigma, normalize a candidate to flatness,
fit observations, or import LCDM distances.  Rmax remains retired.
"""
from __future__ import annotations

import hashlib
import json
import re
import subprocess
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[3]
LAB_ID = "PBUF-FOUNDATION-HOMOGENEOUS-ELASTIC-ENERGY-CLOSURE-AUDIT-001"

THERMAL_PATH = ROOT / "pbuf/data/quantum/thermal_table_cache.json"
RMAX_RETIREMENT_PATH = ROOT / "pbuf/labs/foundation/rmax_retirement_background_reformulation001.py"
SOURCE_PATHS = [
    ROOT / "matter001_derivation.py",
    ROOT / "medium001_derivation.py",
    ROOT / "pbuf/labs/foundation/medium_dynamics_no_fit_closure001.py",
    ROOT / "constitutive_principles001.py",
    ROOT / "em_transport001.py",
    ROOT / "transport_research001.py",
    ROOT / "v11_alpha_audit.py",
]

# A source closes the requested bridge only if it supplies all four links with
# actual equations/definitions rather than prose references to a missing action.
REQUIRED_LINKS = [
    "homogeneous_medium_state",
    "constitutive_energy_density",
    "absolute_energy_normalization",
    "omega_sigma_conversion",
]

FORBIDDEN_SHORTCUTS = {
    "alpha_T_as_Omega_sigma": "Omega_sigma(a)=alpha_T(a)",
    "kmax_as_Omega_sigma": "Omega_sigma(a)=k_max(a)",
    "flatness_normalization": "choose/rescale amplitude so Omega totals unity",
    "retired_Rmax_chain": "reuse decay/S/Omega_sigma_raw/sigma_rescale",
    "replacement_parameter": "introduce a new free activation or amplitude parameter",
}


def git_text(*args: str) -> str:
    try:
        return subprocess.check_output(["git", *args], cwd=ROOT, text=True).strip()
    except Exception:
        return ""


def sha256(path: Path) -> str | None:
    if not path.is_file():
        return None
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def load_thermal() -> dict[str, Any]:
    if not THERMAL_PATH.is_file():
        return {"present": False, "readable": False}
    try:
        data = json.loads(THERMAL_PATH.read_text())
    except Exception as exc:
        return {"present": True, "readable": False, "error": repr(exc)}
    rows = data.get("rows") if isinstance(data, dict) else None
    valid = []
    if isinstance(rows, list):
        for row in rows:
            if not isinstance(row, dict):
                continue
            try:
                a = float(row["a"])
                alpha_t = float(row["alpha_T"])
                eps = float(row["epsilon0_T"])
            except Exception:
                continue
            valid.append((a, alpha_t, eps, eps - alpha_t))
    return {
        "present": True,
        "readable": True,
        "row_count": len(rows) if isinstance(rows, list) else 0,
        "valid_row_count": len(valid),
        "a_min": min((r[0] for r in valid), default=None),
        "a_max": max((r[0] for r in valid), default=None),
        "covers_a_equals_1": any(abs(r[0] - 1.0) <= 1e-12 for r in valid),
        "alpha_T_range": [min((r[1] for r in valid), default=None), max((r[1] for r in valid), default=None)],
        "epsilon0_T_range": [min((r[2] for r in valid), default=None), max((r[2] for r in valid), default=None)],
        "kmax_range": [min((r[3] for r in valid), default=None), max((r[3] for r in valid), default=None)],
        "sha256": sha256(THERMAL_PATH),
    }


def source_record(path: Path) -> dict[str, Any]:
    if not path.is_file():
        return {"path": str(path.relative_to(ROOT)), "present": False}
    text = path.read_text(errors="replace")
    low = text.lower()

    # Evidence is deliberately split into positive structural markers and
    # explicit missing/conditional language.  Textual occurrence alone never
    # counts as closure.
    markers = {
        "mentions_medium_action": bool(re.search(r"S_(?:med|sigma)|medium action|elastic action", text, re.I)),
        "mentions_energy_density": bool(re.search(r"energy density|rho[_ ]?sigma|rho_?\{?sigma\}?|Omega_sigma", text, re.I)),
        "mentions_homogeneous": "homogeneous" in low,
        "mentions_action_hessian": bool(re.search(r"hessian|second variation|quadratic kernel", text, re.I)),
        "mentions_absolute_normalization": bool(re.search(r"absolute normalization|dimensionally normalized|normalization", text, re.I)),
        "mentions_alpha_T": "alpha_T" in text,
        "mentions_epsilon0_T": "epsilon0_T" in text or "epsilon_0" in text,
        "mentions_kmax": "k_max" in text or "kmax" in low,
        "explicit_missing_action_language": bool(re.search(
            r"missing.*(?:action|constitutive)|unknown.*(?:action|constitutive)|as-yet unknown local completion|not fixed|not derive|cannot reconstruct|requires.*action",
            text,
            re.I | re.S,
        )),
        "explicit_conditional_only_language": bool(re.search(r"conditional|schematic|family of admissible|not a closed prediction|additional assumption", text, re.I)),
    }

    # Strong closure candidates must be equations tying a medium state/strain to
    # an energy density and then to Omega_sigma with an absolute normalization.
    equation_lines = []
    for no, line in enumerate(text.splitlines(), 1):
        if re.search(r"(?:rho[_ ]?sigma|Omega_sigma|energy density|S_sigma|S_med|W\s*=|psi\s*=)", line, re.I):
            equation_lines.append({"line": no, "text": line.strip()[:500]})

    return {
        "path": str(path.relative_to(ROOT)),
        "present": True,
        "sha256": sha256(path),
        "markers": markers,
        "candidate_lines": equation_lines[:80],
        "candidate_line_count": len(equation_lines),
    }


def evaluate_closure(records: list[dict[str, Any]]) -> dict[str, Any]:
    # Current accepted foundation documents explicitly say the normalized local
    # medium action/constitutive Hessian is missing.  Therefore textual formulas
    # are not promoted unless a source unambiguously supplies every required
    # link and does not simultaneously mark it conditional/missing.
    viable_sources = []
    for record in records:
        if not record.get("present"):
            continue
        m = record["markers"]
        structural = (
            m["mentions_medium_action"]
            and m["mentions_energy_density"]
            and m["mentions_homogeneous"]
            and m["mentions_absolute_normalization"]
        )
        blocked = m["explicit_missing_action_language"] or m["explicit_conditional_only_language"]
        if structural and not blocked:
            viable_sources.append(record["path"])

    links = {
        "homogeneous_medium_state": {
            "closed": False,
            "reason": "No authoritative current source defines the homogeneous dynamical medium variable chi(a) with dimensions and evolution law.",
        },
        "constitutive_energy_density": {
            "closed": False,
            "reason": "No normalized local S_sigma/S_med or constitutive energy functional is supplied from which rho_sigma(a) can be varied/derived.",
        },
        "absolute_energy_normalization": {
            "closed": False,
            "reason": "The thermal quantities are dimensionless and current foundation work does not supply the absolute constitutive energy scale required for rho_sigma.",
        },
        "omega_sigma_conversion": {
            "closed": False,
            "reason": "Without physical rho_sigma(a) and an independently closed critical-density normalization, Omega_sigma(a) cannot be evaluated physically.",
        },
    }

    if viable_sources:
        # Presence is surfaced for manual review only; this audit is not allowed
        # to silently override the explicit missing-action findings already in
        # the accepted foundation chain.
        for link in links.values():
            link["candidate_sources_for_manual_review"] = viable_sources

    return {
        "required_links": REQUIRED_LINKS,
        "links": links,
        "all_links_closed": all(v["closed"] for v in links.values()),
        "viable_sources_for_manual_review": viable_sources,
    }


def main() -> int:
    thermal = load_thermal()
    sources = [source_record(path) for path in SOURCE_PATHS]
    closure_map = evaluate_closure(sources)

    rmax_retirement_present = RMAX_RETIREMENT_PATH.is_file()
    status = (
        "HOMOGENEOUS_ELASTIC_ENERGY_CLOSURE_DERIVED"
        if closure_map["all_links_closed"]
        else "HOMOGENEOUS_ELASTIC_ENERGY_CLOSURE_NOT_YET_DERIVED"
    )

    closure = {
        "status": status,
        "Rmax_status": "RETIRED_FROM_ACTIVE_PBUF",
        "replacement_free_parameter_introduced": False,
        "retained_inputs": ["alpha_T(a)", "epsilon0_T(a)", "k_max(a)=epsilon0_T(a)-alpha_T(a)"],
        "rho_sigma_physically_closed": closure_map["links"]["constitutive_energy_density"]["closed"],
        "Omega_sigma_physically_closed": closure_map["all_links_closed"],
        "E_of_a_physically_closed": False,
        "H_of_a_physically_closed": False,
        "distance_integration_allowed": False,
        "future_turnaround_role": "OUTPUT_ONLY_IF_FUTURE_DYNAMICS_YIELD_H_EQUALS_ZERO",
        "safe_next": (
            "The current repository does not yet supply the normalized covariant medium action/constitutive energy functional needed to derive homogeneous rho_sigma(a). "
            "Do not identify alpha_T or k_max with Omega_sigma and do not restore the retired Rmax activation chain. "
            "The next physical closure must derive the medium variable, energy functional, and absolute constitutive normalization together; until then keep Omega_sigma(a), E(a), H(a), and distances open."
        ),
    }

    checks = {
        "Rmax_retired_from_active_PBUF": True,
        "Rmax_numeric_value_used": False,
        "replacement_free_parameter_introduced": False,
        "future_reversal_assumed": False,
        "thermal_cache_readable": bool(thermal.get("readable")),
        "alpha_T_not_promoted_to_Omega_sigma": True,
        "kmax_not_promoted_to_Omega_sigma": True,
        "flatness_not_used_to_normalize_candidate": True,
        "retired_Rmax_chain_not_reused": True,
        "no_lcdm_distance_imported": True,
        "no_cluster_redshift_used_as_expansion_redshift": True,
        "no_kappa_or_lensing_target_used": True,
        "no_G_backsolve": True,
        "legacy_0p18_used": False,
        "fit_or_tuning_used": False,
        "quantum_engine_executed": False,
        "planck_scale_used": False,
        "gravity_fundamental_in_PBUF": False,
        "no_tracked_or_staged_changes": git_text("diff", "--name-only") == "" and git_text("diff", "--cached", "--name-only") == "",
        "stdout_only_no_run_directory_created": True,
    }

    result = {
        "status": "FACT_FINDING_ONLY",
        "lab_id": LAB_ID,
        "policy": {
            "Rmax_role": "RETIRED_HISTORICAL_PLACEHOLDER",
            "replacement_parameter_allowed": False,
            "future_reversal_role": "NOT_ASSUMED_MAY_EMERGE_AS_OUTPUT",
            "gravity_fundamental_in_PBUF": False,
            "lensing_target_used": False,
            "fit_or_tuning_used": False,
            "forbidden_shortcuts": FORBIDDEN_SHORTCUTS,
        },
        "repo_state": {
            "repository": "TheExiledMonk/lab",
            "branch": git_text("branch", "--show-current"),
            "head_sha": git_text("rev-parse", "HEAD"),
            "tracked_changes": git_text("diff", "--name-only"),
            "staged_changes": git_text("diff", "--cached", "--name-only"),
        },
        "sources": {
            "thermal_cache": thermal,
            "rmax_retirement_audit_present": rmax_retirement_present,
            "rmax_retirement_audit_sha256": sha256(RMAX_RETIREMENT_PATH),
            "constitutive_source_inventory": sources,
        },
        "closure_map": closure_map,
        "closure": closure,
        "checks": checks,
    }

    print(LAB_ID)
    print("status=FACT_FINDING_ONLY")
    print(f"head_sha={result['repo_state']['head_sha']}")
    print("gravity_fundamental_in_PBUF=false")
    print("Rmax_status=RETIRED_FROM_ACTIVE_PBUF")
    print("replacement_free_parameter_introduced=false")
    print("future_reversal_assumed=false")
    print("alpha_T_promoted_to_Omega_sigma=false")
    print("kmax_promoted_to_Omega_sigma=false")
    print("flatness_normalization_used=false")
    print("lensing_target_used=false")
    print("fit_or_tuning_used=false")
    print()
    print("HOMOGENEOUS_ELASTIC_ENERGY_REQUIRED_LINKS")
    print("key | closed | reason")
    for key in REQUIRED_LINKS:
        row = closure_map["links"][key]
        print(f"{key} | {str(row['closed']).lower()} | {row['reason']}")
    print()
    print("SOURCE_SUMMARY")
    print("path | present | candidate_lines | missing_or_conditional")
    for record in sources:
        if not record.get("present"):
            print(f"{record['path']} | false | 0 | true")
            continue
        m = record["markers"]
        blocked = m["explicit_missing_action_language"] or m["explicit_conditional_only_language"]
        print(f"{record['path']} | true | {record['candidate_line_count']} | {str(blocked).lower()}")
    print()
    print("CONCLUSION")
    print(f"status={closure['status']}")
    print(f"rho_sigma_physically_closed={str(closure['rho_sigma_physically_closed']).lower()}")
    print(f"Omega_sigma_physically_closed={str(closure['Omega_sigma_physically_closed']).lower()}")
    print("E_of_a_physically_closed=false")
    print("H_of_a_physically_closed=false")
    print("distance_integration_allowed=false")
    print(f"safe_next={closure['safe_next']}")
    print()
    print("CHECKS")
    for key, value in checks.items():
        print(f"{key}={str(value).lower() if isinstance(value, bool) else value}")
    print("JSON=" + json.dumps(result, sort_keys=True, separators=(",", ":")))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
