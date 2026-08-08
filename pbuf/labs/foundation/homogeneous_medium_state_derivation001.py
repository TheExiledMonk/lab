#!/usr/bin/env python3
"""PBUF foundation audit: homogeneous medium state derivation 001.

Purpose
-------
Test whether current PBUF definitions already justify a homogeneous medium
state variable chi(a) without importing cosmological targets, fitting a new
parameter, or simply renaming alpha_T, epsilon0_T, or k_max for convenience.

The accepted Rmax-free thermal inputs are:
    alpha_T(a)
    epsilon0_T(a)
    k_max(a) = epsilon0_T(a) - alpha_T(a)

A candidate homogeneous state survives only if repository evidence establishes
its physical role, dimensional character, and evolution/state meaning. Numeric
correlation or convenient algebra is insufficient.

Fact-finding / derivation-boundary audit only. No lensing target, no LCDM
substitution, no G backsolve, no Quantum Engine execution, no Planck input,
and no run directory.
"""
from __future__ import annotations

import hashlib
import json
import re
import subprocess
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[3]
THERMAL_PATH = ROOT / "pbuf/data/quantum/thermal_table_cache.json"

LAB_ID = "PBUF-FOUNDATION-HOMOGENEOUS-MEDIUM-STATE-DERIVATION-001"

SOURCE_PATHS = [
    "matter001_derivation.py",
    "medium001_derivation.py",
    "constitutive_principles001.py",
    "em_transport001.py",
    "transport_research001.py",
    "v11_alpha_audit.py",
    "pbuf/labs/foundation/medium_dynamics_no_fit_closure001.py",
    "pbuf/labs/foundation/homogeneous_elastic_energy_closure_audit001.py",
]

CANDIDATES = {
    "alpha_T": {
        "expression": "chi(a)=alpha_T(a)",
        "required_role": "authoritative identification of alpha_T as the homogeneous medium state variable",
    },
    "epsilon0_T": {
        "expression": "chi(a)=epsilon0_T(a)",
        "required_role": "authoritative identification of epsilon0_T as the homogeneous medium state variable",
    },
    "k_max": {
        "expression": "chi(a)=k_max(a)=epsilon0_T(a)-alpha_T(a)",
        "required_role": "authoritative identification of k_max as state rather than bound/response/coefficient",
    },
    "separate_chi": {
        "expression": "chi(a) is a separate variable governed by retained thermal quantities",
        "required_role": "authoritative definition of chi, dimensions, and an evolution/state equation",
    },
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
                alpha = float(row["alpha_T"])
                eps = float(row["epsilon0_T"])
            except Exception:
                continue
            valid.append((a, alpha, eps, eps - alpha))
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


def source_inventory() -> list[dict[str, Any]]:
    rows = []
    patterns = {
        "chi_definition": re.compile(r"\bchi(?:\^?[A-Za-z0-9_]+)?\b.{0,80}(?:define|defined|state|field|variable|dimension)", re.I),
        "alpha_role": re.compile(r"alpha_T.{0,120}(?:state|field|variable|coefficient|amplitude|thermal|response)", re.I),
        "epsilon_role": re.compile(r"epsilon(?:0|_0)[,_]?T.{0,120}(?:state|field|variable|coefficient|thermal|response|bound)", re.I),
        "kmax_role": re.compile(r"k_max.{0,120}(?:state|field|variable|coefficient|response|bound|saturation|max)", re.I),
        "missing_chi": re.compile(r"(?:missing|does not|not supplied|unknown).{0,100}\bchi\b|\bchi\b.{0,100}(?:missing|not supplied|unknown)", re.I),
        "dimensions_chi": re.compile(r"\bchi\b.{0,100}(?:dimension|units)|(?:dimension|units).{0,100}\bchi\b", re.I),
        "evolution_chi": re.compile(r"\bchi\b.{0,120}(?:evolution|dynamics|equation of motion|equation|dot\{?chi|dchi)", re.I),
    }
    for rel in SOURCE_PATHS:
        path = ROOT / rel
        if not path.is_file():
            rows.append({"path": rel, "present": False})
            continue
        text = path.read_text(errors="replace")
        lines = text.splitlines()
        matches = []
        for i, line in enumerate(lines, 1):
            reasons = [name for name, pat in patterns.items() if pat.search(line)]
            if reasons:
                matches.append({"line": i, "text": line[:500], "reasons": reasons})
        rows.append({
            "path": rel,
            "present": True,
            "sha256": sha256(path),
            "matches": matches,
            "markers": {name: bool(pat.search(text)) for name, pat in patterns.items()},
        })
    return rows


def candidate_assessment(inventory: list[dict[str, Any]]) -> list[dict[str, Any]]:
    corpus = "\n".join(
        (ROOT / rel).read_text(errors="replace")
        for rel in SOURCE_PATHS
        if (ROOT / rel).is_file()
    )

    explicit_missing_chi = bool(re.search(
        r"(?:missing|not supplied|does not define|unknown).{0,150}\bchi\b|\bchi\b.{0,150}(?:missing|not supplied|unknown)",
        corpus,
        re.I | re.S,
    ))

    assessments = []
    for key, meta in CANDIDATES.items():
        if key == "alpha_T":
            direct_identity = bool(re.search(r"\bchi\b\s*(?:=|:=|is)\s*alpha_T\b|alpha_T\b\s*(?:=|:=|is)\s*\bchi\b", corpus, re.I))
            role_support = bool(re.search(r"alpha_T.{0,100}(?:thermal|amplitude|coefficient)", corpus, re.I | re.S))
            state_support = bool(re.search(r"alpha_T.{0,100}(?:homogeneous )?state variable", corpus, re.I | re.S))
        elif key == "epsilon0_T":
            direct_identity = bool(re.search(r"\bchi\b\s*(?:=|:=|is)\s*epsilon(?:0|_0)[,_]?T\b|epsilon(?:0|_0)[,_]?T\b\s*(?:=|:=|is)\s*\bchi\b", corpus, re.I))
            role_support = bool(re.search(r"epsilon(?:0|_0)[,_]?T.{0,100}(?:thermal|response|bound|coefficient)", corpus, re.I | re.S))
            state_support = bool(re.search(r"epsilon(?:0|_0)[,_]?T.{0,100}(?:homogeneous )?state variable", corpus, re.I | re.S))
        elif key == "k_max":
            direct_identity = bool(re.search(r"\bchi\b\s*(?:=|:=|is)\s*k_max\b|k_max\b\s*(?:=|:=|is)\s*\bchi\b", corpus, re.I))
            role_support = bool(re.search(r"k_max.{0,100}(?:response|bound|saturation|max|difference)", corpus, re.I | re.S))
            state_support = bool(re.search(r"k_max.{0,100}(?:homogeneous )?state variable", corpus, re.I | re.S))
        else:
            direct_identity = bool(re.search(r"\bchi\b.{0,100}(?:physical medium variable|medium variable|deformation variable|field)", corpus, re.I | re.S))
            role_support = direct_identity
            state_support = bool(re.search(r"\bchi\b.{0,160}(?:homogeneous|cosmolog).{0,120}(?:evolution|dynamics|equation|state)", corpus, re.I | re.S))

        has_dimensions = bool(re.search(r"\bchi\b.{0,120}(?:dimension|units)|(?:dimension|units).{0,120}\bchi\b", corpus, re.I | re.S))
        has_evolution = bool(re.search(r"\bchi\b.{0,160}(?:evolution law|equation of motion|dynamics|dchi|dot\{?chi)", corpus, re.I | re.S))

        if direct_identity and state_support and has_dimensions and has_evolution and not explicit_missing_chi:
            status = "DERIVED_CANDIDATE"
            reason = "repository supplies direct identity, state role, dimensions, and evolution law"
        elif key == "separate_chi" and direct_identity and not (has_dimensions and has_evolution):
            status = "CONCEPT_PRESENT_BUT_NOT_CLOSED"
            reason = "chi appears as a generic medium/deformation field but its normalized homogeneous state and evolution are not supplied"
        else:
            status = "NOT_DERIVED"
            reason = "no authoritative chain establishes identity, state role, dimensions, and evolution together"

        assessments.append({
            "candidate": key,
            "expression": meta["expression"],
            "status": status,
            "reason": reason,
            "direct_identity_found": direct_identity,
            "role_context_found": role_support,
            "explicit_state_role_found": state_support,
            "chi_dimensions_closed": has_dimensions,
            "chi_evolution_closed": has_evolution,
            "repository_explicitly_reports_missing_chi": explicit_missing_chi,
        })
    return assessments


def main() -> int:
    thermal = load_thermal()
    inventory = source_inventory()
    candidates = candidate_assessment(inventory)

    derived = [c for c in candidates if c["status"] == "DERIVED_CANDIDATE"]
    chi_closed = len(derived) == 1

    if chi_closed:
        closure_status = "HOMOGENEOUS_MEDIUM_STATE_DERIVED"
        selected = derived[0]["candidate"]
        safe_next = "With one uniquely derived homogeneous medium state, derive its normalized constitutive energy functional without cosmological fitting."
    else:
        closure_status = "HOMOGENEOUS_MEDIUM_STATE_NOT_YET_DERIVED"
        selected = None
        safe_next = (
            "Current PBUF does not uniquely define a homogeneous chi(a) with physical role, dimensions, and evolution law. "
            "Do not rename alpha_T, epsilon0_T, or k_max as chi by convenience. The next closure must derive or postulate "
            "the medium variable and its dynamics from constitutive microphysics before an elastic energy density can be built."
        )

    checks = {
        "Rmax_retired_from_active_PBUF": True,
        "Rmax_numeric_value_used": False,
        "replacement_free_parameter_introduced": False,
        "future_reversal_assumed": False,
        "thermal_cache_readable": bool(thermal.get("readable")),
        "alpha_T_not_auto_promoted_to_chi": not any(c["candidate"] == "alpha_T" and c["status"] == "DERIVED_CANDIDATE" for c in candidates),
        "epsilon0_T_not_auto_promoted_to_chi": not any(c["candidate"] == "epsilon0_T" and c["status"] == "DERIVED_CANDIDATE" for c in candidates),
        "kmax_not_auto_promoted_to_chi": not any(c["candidate"] == "k_max" and c["status"] == "DERIVED_CANDIDATE" for c in candidates),
        "no_cosmological_target_used": True,
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
        "status": "FACT_FINDING_DERIVATION_BOUNDARY",
        "lab_id": LAB_ID,
        "policy": {
            "Rmax_role": "RETIRED_FROM_ACTIVE_PBUF",
            "candidate_selection_rule": "identity+physical_role+dimensions+evolution_required",
            "numeric_similarity_or_convenience_is_insufficient": True,
            "replacement_parameter_allowed": False,
            "future_reversal_role": "OUTPUT_ONLY_IF_LATER_DYNAMICS_YIELD_H_EQUALS_ZERO",
            "gravity_fundamental_in_PBUF": False,
        },
        "repo_state": {
            "repository": "TheExiledMonk/lab",
            "branch": git_text("branch", "--show-current"),
            "head_sha": git_text("rev-parse", "HEAD"),
            "tracked_changes": git_text("diff", "--name-only"),
            "staged_changes": git_text("diff", "--cached", "--name-only"),
        },
        "thermal_inputs": thermal,
        "source_inventory": inventory,
        "candidate_assessments": candidates,
        "closure": {
            "status": closure_status,
            "homogeneous_chi_physically_closed": chi_closed,
            "selected_candidate": selected,
            "constitutive_energy_derivation_allowed": chi_closed,
            "Omega_sigma_derivation_allowed": False,
            "safe_next": safe_next,
        },
        "checks": checks,
    }

    print(LAB_ID)
    print("status=FACT_FINDING_DERIVATION_BOUNDARY")
    print(f"head_sha={result['repo_state']['head_sha']}")
    print("gravity_fundamental_in_PBUF=false")
    print("Rmax_status=RETIRED_FROM_ACTIVE_PBUF")
    print("replacement_free_parameter_introduced=false")
    print("future_reversal_assumed=false")
    print("cosmological_target_used=false")
    print("fit_or_tuning_used=false")
    print()
    print("HOMOGENEOUS_MEDIUM_STATE_CANDIDATES")
    print("candidate | status | direct_identity | state_role | dimensions | evolution")
    for c in candidates:
        print(
            f"{c['candidate']} | {c['status']} | {str(c['direct_identity_found']).lower()} | "
            f"{str(c['explicit_state_role_found']).lower()} | {str(c['chi_dimensions_closed']).lower()} | "
            f"{str(c['chi_evolution_closed']).lower()}"
        )
    print()
    print("CONCLUSION")
    print(f"status={closure_status}")
    print(f"homogeneous_chi_physically_closed={str(chi_closed).lower()}")
    print(f"selected_candidate={selected}")
    print(f"constitutive_energy_derivation_allowed={str(chi_closed).lower()}")
    print("Omega_sigma_derivation_allowed=false")
    print(f"safe_next={safe_next}")
    print()
    print("CHECKS")
    for key, value in checks.items():
        print(f"{key}={str(value).lower() if isinstance(value, bool) else value}")
    print("JSON=" + json.dumps(result, sort_keys=True, separators=(",", ":")))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
