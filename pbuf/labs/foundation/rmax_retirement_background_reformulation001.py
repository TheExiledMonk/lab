#!/usr/bin/env python3
"""PBUF foundation audit: retire Rmax and test an Rmax-free background reformulation.

Rmax is treated as a historical placeholder for a hypothetical future expansion
reversal, not as a physical input to current PBUF.  This audit does not assign,
recover, fit, inherit, or replace Rmax.  It classifies the historical V11
background equations into:

- independently retained relations;
- relations retired because their physical content depends on Rmax;
- structural relations that may survive only after the elastic background is
  independently rederived;
- open external normalisations.

The audit also inventories repository evidence for any independently derived
present-era Omega_sigma(a) law that does not use Rmax.  Presence of a textual
candidate is not promoted to physical closure.

Fact-finding only.  No lensing target, no LCDM distance substitution, no
measured-G backsolve, no Quantum Engine execution, no replacement free
parameter, and no run directory.
"""
from __future__ import annotations

import hashlib
import json
import re
import subprocess
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[3]
TRACE_PATH = ROOT / "v11_alpha_audit.py"
THERMAL_PATH = ROOT / "pbuf/data/quantum/thermal_table_cache.json"
PRIOR_RMAX_AUDIT = ROOT / "pbuf/labs/foundation/rmax_forbidden_background_closure_audit001.py"

LAB_ID = "PBUF-FOUNDATION-RMAX-RETIREMENT-BACKGROUND-REFORMULATION-001"

# Historical V11 relations.  The classification is a physics/provenance
# classification, not a numerical rewrite.
RELATIONS = [
    {
        "key": "thermal_alpha_T_epsilon0_T",
        "expression": "alpha_T(a), epsilon0_T(a)",
        "classification": "RETAIN_INDEPENDENT",
        "reason": "direct thermal/LUT inputs; no Rmax dependency",
    },
    {
        "key": "kmax",
        "expression": "k_max(a)=epsilon0_T(a)-alpha_T(a)",
        "classification": "RETAIN_INDEPENDENT",
        "reason": "thermal response difference; algebraically independent of Rmax",
    },
    {
        "key": "decay_activation",
        "expression": "decay(a)=exp(-a/Rmax)",
        "classification": "RETIRE_WITH_RMAX",
        "reason": "contains the retired future-reversal scale explicitly",
    },
    {
        "key": "saturation_S",
        "expression": "S(a)=1-(1-k_max(a))*decay(a)",
        "classification": "RETIRE_WITH_RMAX",
        "reason": "its historical form inherits the retired decay/Rmax factor",
    },
    {
        "key": "omega_sigma_raw",
        "expression": "Omega_sigma_raw(a)=alpha_T(a)*(1-decay(a))*S(a)",
        "classification": "RETIRE_WITH_RMAX",
        "reason": "historical elastic-background amplitude is built from retired activation terms",
    },
    {
        "key": "omega_sigma_target",
        "expression": "Omega_sigma_target=1-Omega_m0-Omega_r0-alpha_resolved",
        "classification": "HISTORICAL_NORMALISATION_REVIEW",
        "reason": "does not contain Rmax, but belongs to the historical flat-today normalisation and requires independent review",
    },
    {
        "key": "sigma_rescale",
        "expression": "sigma_rescale=Omega_sigma_target/Omega_sigma_raw(a=1)",
        "classification": "RETIRE_WITH_RMAX",
        "reason": "normalises the retired Omega_sigma_raw construction",
    },
    {
        "key": "omega_sigma",
        "expression": "Omega_sigma(a)=sigma_rescale*Omega_sigma_raw(a)",
        "classification": "REQUIRES_INDEPENDENT_REDERIVATION",
        "reason": "the elastic background concept may survive, but the historical formula is retired with its Rmax-dependent construction",
    },
    {
        "key": "E_of_a",
        "expression": "E(a)^2=Omega_m0*a^-3+Omega_r0*a^-4+Omega_sigma(a)",
        "classification": "STRUCTURE_ONLY_PENDING_OMEGA_SIGMA_REDERIVATION",
        "reason": "background structure cannot be physically evaluated until an Rmax-free Omega_sigma(a) is derived",
    },
    {
        "key": "H_of_a",
        "expression": "H(a)=H0*E(a)",
        "classification": "STRUCTURE_ONLY_PENDING_E_AND_H0",
        "reason": "kinematic scaling remains formal, but physical H(a) requires closed E(a) and H0 provenance",
    },
]

# Audit-generated files and historical provenance files are not evidence for a
# new constitutive law.  They are excluded from candidate promotion.
EXCLUDED_CANDIDATE_NAMES = {
    "rmax_retirement_background_reformulation001.py",
    "rmax_forbidden_background_closure_audit001.py",
    "pbuf_background_distance_recovery_audit001.py",
    "pbuf_background_parameter_provenance_audit001.py",
    "cluster_distance_redshift_closure_audit001.py",
    "v11_alpha_audit.py",
}

TEXT_SUFFIXES = {".py", ".md", ".txt", ".json", ".yaml", ".yml", ".toml"}


def sha256(path: Path) -> str | None:
    if not path.is_file():
        return None
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def git_text(*args: str) -> str:
    try:
        return subprocess.check_output(["git", *args], cwd=ROOT, text=True).strip()
    except Exception:
        return ""


def load_thermal() -> dict[str, Any]:
    if not THERMAL_PATH.is_file():
        return {"present": False, "readable": False}
    try:
        data = json.loads(THERMAL_PATH.read_text())
    except Exception as exc:
        return {"present": True, "readable": False, "error": repr(exc)}
    rows = data.get("rows") if isinstance(data, dict) else None
    valid: list[tuple[float, float, float]] = []
    if isinstance(rows, list):
        for row in rows:
            if not isinstance(row, dict):
                continue
            try:
                valid.append((float(row["a"]), float(row["alpha_T"]), float(row["epsilon0_T"])))
            except Exception:
                continue
    return {
        "present": True,
        "readable": True,
        "row_count": len(rows) if isinstance(rows, list) else 0,
        "valid_row_count": len(valid),
        "a_min": min((v[0] for v in valid), default=None),
        "a_max": max((v[0] for v in valid), default=None),
        "covers_a_equals_1": any(abs(v[0] - 1.0) <= 1e-12 for v in valid),
        "sha256": sha256(THERMAL_PATH),
    }


def trace_markers() -> dict[str, Any]:
    if not TRACE_PATH.is_file():
        return {"present": False, "markers": {}}
    text = TRACE_PATH.read_text(errors="replace")
    markers = {
        "kmax": "k_max(a)=epsilon_0,T(a)-alpha_T(a)" in text,
        "omega_sigma_raw": "Omega_sigma_raw(a)=alpha_T(a)(1-decay(a))S(a)" in text,
        "flat_today": "Omega_sigma_target=1-Omega_m0-Omega_r0-alpha_resolved" in text,
        "omega_sigma_rmax_model": bool(re.search(r"Omega_sigma\(a\)=alpha\(1-exp\(-a/R_max\)\)S\(a\)", text)),
        "background_to_E_H": '"target":"E(a), H(a)"' in text,
    }
    return {
        "present": True,
        "sha256": sha256(TRACE_PATH),
        "markers": markers,
        "all_required_trace_markers_present": all(markers.values()),
    }


def _is_text_candidate(path: Path) -> bool:
    if path.suffix.lower() not in TEXT_SUFFIXES:
        return False
    parts = set(path.parts)
    if ".git" in parts or "runs" in parts:
        return False
    if path.name in EXCLUDED_CANDIDATE_NAMES:
        return False
    return True


def inventory_rmax_occurrences() -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    pattern = re.compile(r"\bR_?max\b|\brmax\b", re.I)
    for path in sorted(ROOT.rglob("*")):
        if not path.is_file() or not _is_text_candidate(path):
            continue
        try:
            text = path.read_text(errors="replace")
        except Exception:
            continue
        if not pattern.search(text):
            continue
        hits = []
        for idx, line in enumerate(text.splitlines(), start=1):
            if pattern.search(line):
                hits.append({"line": idx, "context": line.strip()[:500]})
        if hits:
            rows.append({"path": str(path.relative_to(ROOT)), "hit_count": len(hits), "hits": hits[:20]})
    return rows


def inventory_rmax_free_omega_sigma_candidates() -> list[dict[str, Any]]:
    """Find textual Omega_sigma expressions not containing Rmax on the same line.

    These are candidates for human/source review only.  A hit is not evidence of
    an independently derived constitutive law.
    """
    rows: list[dict[str, Any]] = []
    omega = re.compile(r"Omega_sigma|omega_sigma", re.I)
    rmax = re.compile(r"R_?max|rmax", re.I)
    assignment = re.compile(r"Omega_sigma\s*\(|Omega_sigma\s*=|omega_sigma\s*=", re.I)
    for path in sorted(ROOT.rglob("*")):
        if not path.is_file() or not _is_text_candidate(path):
            continue
        try:
            lines = path.read_text(errors="replace").splitlines()
        except Exception:
            continue
        for idx, line in enumerate(lines, start=1):
            if not omega.search(line) or rmax.search(line):
                continue
            context = " ".join(x.strip() for x in lines[max(0, idx-2): min(len(lines), idx+1)])[:900]
            if assignment.search(line):
                rows.append({
                    "path": str(path.relative_to(ROOT)),
                    "line": idx,
                    "expression_line": line.strip()[:500],
                    "context": context,
                    "status": "TEXTUAL_CANDIDATE_NOT_PHYSICAL_CLOSURE",
                })
    return rows


def relation_counts() -> dict[str, int]:
    counts: dict[str, int] = {}
    for row in RELATIONS:
        counts[row["classification"]] = counts.get(row["classification"], 0) + 1
    return counts


def main() -> int:
    thermal = load_thermal()
    trace = trace_markers()
    rmax_occurrences = inventory_rmax_occurrences()
    omega_candidates = inventory_rmax_free_omega_sigma_candidates()

    retained = [r for r in RELATIONS if r["classification"] == "RETAIN_INDEPENDENT"]
    retired = [r for r in RELATIONS if r["classification"] == "RETIRE_WITH_RMAX"]
    rederive = [r for r in RELATIONS if "REDERIVATION" in r["classification"] or "PENDING" in r["classification"]]

    # A repository text hit cannot establish a physical law.  This audit only
    # closes the reformulation if an independently derived law is already
    # explicitly designated as such by an authoritative current source.  No
    # such promotion is performed automatically here.
    independently_derived_rmax_free_omega_sigma_found = False

    closure = {
        "status": "PRESENT_ERA_BACKGROUND_LAW_NOT_YET_DERIVED",
        "Rmax_status": "RETIRED_FROM_ACTIVE_PBUF",
        "Rmax_numeric_value": None,
        "replacement_free_parameter_introduced": False,
        "future_reversal_assumed": False,
        "thermal_inputs_retained": bool(thermal.get("readable")),
        "kmax_retained": True,
        "historical_decay_activation_retired": True,
        "historical_saturation_S_retired": True,
        "historical_omega_sigma_raw_retired": True,
        "historical_sigma_rescale_retired": True,
        "omega_sigma_concept_requires_rederivation": True,
        "independently_derived_Rmax_free_omega_sigma_found": independently_derived_rmax_free_omega_sigma_found,
        "E_of_a_physically_closed": False,
        "H_of_a_physically_closed": False,
        "distance_integration_allowed": False,
        "future_turnaround_role": "OUTPUT_ONLY_IF_DYNAMICS_EVENTUALLY_YIELD_H_EQUALS_ZERO",
        "safe_next": (
            "Derive an Rmax-free present-era elastic background from independently justified constitutive "
            "medium physics. Retain alpha_T(a), epsilon0_T(a), and k_max(a) as available inputs; do not "
            "reuse the retired decay/S/Omega_sigma_raw construction and do not introduce a replacement "
            "free parameter. If no constitutive Omega_sigma(a) law can yet be derived, keep E(a), H(a), "
            "and distance integration open."
        ),
    }

    checks = {
        "Rmax_retired_from_active_PBUF": closure["Rmax_status"] == "RETIRED_FROM_ACTIVE_PBUF",
        "Rmax_numeric_value_is_null": closure["Rmax_numeric_value"] is None,
        "no_replacement_free_parameter": not closure["replacement_free_parameter_introduced"],
        "future_reversal_not_assumed": not closure["future_reversal_assumed"],
        "thermal_cache_readable": bool(thermal.get("readable")),
        "v11_trace_present": bool(trace.get("all_required_trace_markers_present")),
        "kmax_retained_without_Rmax": closure["kmax_retained"],
        "retired_relations_not_used_to_build_E_or_H": not closure["E_of_a_physically_closed"] and not closure["H_of_a_physically_closed"],
        "no_cluster_redshift_used_as_expansion_redshift": True,
        "no_distance_values_fabricated": True,
        "no_lcdm_distance_imported": True,
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
            "Rmax_role": "RETIRED_HISTORICAL_FUTURE_REVERSAL_PLACEHOLDER",
            "Rmax_numeric_use_allowed": False,
            "replacement_parameter_allowed": False,
            "future_reversal_role": "NOT_ASSUMED_MAY_EMERGE_AS_OUTPUT",
            "cluster_redshift_role": "OBSERVED_TOTAL_SHIFT_NOT_EXPANSION_BY_ASSUMPTION",
            "fit_or_tuning_used": False,
            "lensing_target_used": False,
            "gravity_fundamental_in_PBUF": False,
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
            "v11_trace": trace,
            "prior_rmax_audit_present": PRIOR_RMAX_AUDIT.is_file(),
            "prior_rmax_audit_sha256": sha256(PRIOR_RMAX_AUDIT),
        },
        "relation_classification": RELATIONS,
        "relation_classification_counts": relation_counts(),
        "retained_relations": retained,
        "retired_relations": retired,
        "relations_requiring_rederivation_or_open_structure": rederive,
        "repository_Rmax_occurrences": rmax_occurrences,
        "Rmax_free_Omega_sigma_textual_candidates": omega_candidates,
        "closure": closure,
        "checks": checks,
    }

    print(LAB_ID)
    print("status=FACT_FINDING_ONLY")
    print(f"head_sha={result['repo_state']['head_sha']}")
    print("gravity_fundamental_in_PBUF=false")
    print("Rmax_status=RETIRED_FROM_ACTIVE_PBUF")
    print("Rmax_numeric_value=None")
    print("replacement_free_parameter_introduced=false")
    print("future_reversal_assumed=false")
    print("cluster_redshift_used_as_expansion_redshift=false")
    print("lcdm_distance_imported=false")
    print("lensing_target_used=false")
    print("fit_or_tuning_used=false")
    print()
    print("RMAX_RETIREMENT_RELATION_CLASSIFICATION")
    print("key | classification | reason")
    for row in RELATIONS:
        print(f"{row['key']} | {row['classification']} | {row['reason']}")
    print()
    print("REPOSITORY_INVENTORY")
    print(f"files_with_Rmax_occurrences={len(rmax_occurrences)}")
    print(f"Rmax_free_Omega_sigma_textual_candidates={len(omega_candidates)}")
    print("textual_candidates_promoted_to_physical_law=false")
    print()
    print("CONCLUSION")
    print(f"status={closure['status']}")
    print("Rmax_retired_from_active_PBUF=true")
    print("kmax_retained=true")
    print("historical_decay_activation_retired=true")
    print("historical_saturation_S_retired=true")
    print("historical_omega_sigma_raw_retired=true")
    print("historical_sigma_rescale_retired=true")
    print("omega_sigma_concept_requires_rederivation=true")
    print("E_of_a_physically_closed=false")
    print("H_of_a_physically_closed=false")
    print("distance_integration_allowed=false")
    print(f"future_turnaround_role={closure['future_turnaround_role']}")
    print(f"safe_next={closure['safe_next']}")
    print()
    print("CHECKS")
    for key, value in checks.items():
        print(f"{key}={str(value).lower() if isinstance(value, bool) else value}")
    print("JSON=" + json.dumps(result, sort_keys=True, separators=(",", ":")))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
