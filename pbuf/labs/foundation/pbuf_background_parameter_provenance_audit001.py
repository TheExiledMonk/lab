#!/usr/bin/env python3
"""PBUF FOUNDATION — BACKGROUND PARAMETER PROVENANCE AUDIT 001.

Fact-finding only.

Purpose
-------
Recover authoritative/current numerical provenance for the five background
parameters left open by PBUF-BACKGROUND-DISTANCE-RECOVERY-AUDIT-001:

    Rmax
    alpha_resolved
    Omega_r0
    BARYON_FRACTION
    H0

The audit searches repository text/data sources, extracts nearby numeric
candidates, records file provenance and context, and applies deliberately
conservative promotion rules. A value is promoted only when current repository
material explicitly binds the target symbol/role to a single numeric value and
no incompatible current candidate is found.

Historical equations, optimisation outputs, cached metadata, comments, tests,
and unrelated symbols are evidence but are not automatically authoritative.
In particular thermal-cache alpha_qm must not be silently identified with
alpha_resolved merely because the values are numerically close.

No lensing target, fitting, LambdaCDM distance substitution, measured-G
backsolve, Quantum Engine execution, or Planck-scale input is permitted.
"""
from __future__ import annotations

import json
import math
import re
import subprocess
from dataclasses import dataclass
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
LAB_ID = "PBUF-FOUNDATION-BACKGROUND-PARAMETER-PROVENANCE-AUDIT-001"

TARGETS = {
    "Rmax": {
        "aliases": [r"\bRmax\b", r"\bR_max\b", r"\bRMAX\b"],
        "role": "elastic activation scale in decay(a)=exp(-a/Rmax)",
    },
    "alpha_resolved": {
        "aliases": [r"\balpha_resolved\b", r"\balpha resolved\b", r"\bα_resolved\b"],
        "role": "resolved elastic amplitude used by V11 background normalization",
    },
    "Omega_r0": {
        "aliases": [r"\bOmega_r0\b", r"\bΩ_r0\b", r"\bomega_r0\b"],
        "role": "present radiation density parameter",
    },
    "BARYON_FRACTION": {
        "aliases": [r"\bBARYON_FRACTION\b", r"\bbaryon_fraction\b", r"\bf_b\b"],
        "role": "rule/value connecting Omega_b0 to Omega_m0",
    },
    "H0": {
        "aliases": [r"\bH0\b", r"\bH_0\b", r"\bHubble[_ ]?constant\b"],
        "role": "absolute background expansion / distance scale",
    },
}

TEXT_SUFFIXES = {".py", ".md", ".txt", ".json", ".yaml", ".yml", ".toml", ".ini", ".cfg", ".csv"}
EXCLUDE_DIRS = {".git", "runs", "__pycache__", ".pytest_cache", ".mypy_cache", "node_modules"}
MAX_BYTES = 4_000_000
NUMBER_RE = re.compile(r"(?<![A-Za-z0-9_])[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[eE][-+]?\d+)?")

CURRENT_HINTS = (
    "authoritative", "current", "production", "default", "config", "fixed",
    "v11", "equation", "parameter", "constant", "flat_today", "background",
)
WEAK_HINTS = (
    "example", "synthetic", "test", "toy", "placeholder", "candidate",
    "historical", "archive", "benchmark", "lcdm", "comparison", "deprecated",
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


def _iter_files():
    for p in ROOT.rglob("*"):
        if not p.is_file() or p.suffix.lower() not in TEXT_SUFFIXES:
            continue
        rel = p.relative_to(ROOT)
        if any(part in EXCLUDE_DIRS for part in rel.parts):
            continue
        try:
            if p.stat().st_size > MAX_BYTES:
                continue
        except OSError:
            continue
        yield p, rel


def _read_text(path: Path) -> str | None:
    try:
        return path.read_text(errors="replace")
    except Exception:
        return None


def _context_score(rel: Path, line: str, before: str, after: str) -> tuple[int, list[str]]:
    blob = " ".join([str(rel).lower(), before.lower(), line.lower(), after.lower()])
    score = 0
    reasons = []
    for token in CURRENT_HINTS:
        if token in blob:
            score += 1
            reasons.append(f"positive:{token}")
    for token in WEAK_HINTS:
        if token in blob:
            score -= 2
            reasons.append(f"negative:{token}")
    if str(rel).startswith("pbuf/labs/foundation/"):
        score -= 1
        reasons.append("negative:foundation_audit_context")
    if rel.name == "v11_alpha_audit.py":
        score += 2
        reasons.append("positive:v11_trace_source")
    if rel.as_posix() == "pbuf/data/quantum/thermal_table_cache.json":
        reasons.append("special:thermal_cache_metadata")
    return score, reasons


def _extract_candidates(target: str, aliases: list[str]) -> list[dict]:
    alias_re = re.compile("|".join(f"(?:{a})" for a in aliases), re.I)
    hits = []
    for path, rel in _iter_files():
        text = _read_text(path)
        if text is None or not alias_re.search(text):
            continue
        lines = text.splitlines()
        for i, line in enumerate(lines):
            if not alias_re.search(line):
                continue
            before = lines[i-1] if i > 0 else ""
            after = lines[i+1] if i + 1 < len(lines) else ""
            nums = []
            for s in (line, after):
                nums.extend(NUMBER_RE.findall(s))
            values = []
            for n in nums:
                try:
                    v = float(n)
                except ValueError:
                    continue
                if math.isfinite(v):
                    values.append(v)
            score, reasons = _context_score(rel, line, before, after)
            hits.append({
                "path": rel.as_posix(),
                "line": i + 1,
                "context": line.strip()[:500],
                "next_line": after.strip()[:500],
                "numeric_candidates": values[:12],
                "context_score": score,
                "context_reasons": reasons,
            })
    return hits


def _explicit_binding_values(target: str, hits: list[dict]) -> list[dict]:
    out = []
    # Conservative: numeric value must appear on the same line as the target,
    # and the line must look like assignment/configuration/equality.
    binding_tokens = ("=", ":", "default", "fixed", "value", "set to")
    for h in hits:
        line = h["context"]
        if not h["numeric_candidates"]:
            continue
        if not any(tok in line.lower() for tok in binding_tokens):
            continue
        # Exclude obvious equation identifiers/page labels and factor-counting
        # noise by requiring at least neutral context score.
        if h["context_score"] < 0:
            continue
        for v in h["numeric_candidates"]:
            out.append({"value": v, "path": h["path"], "line": h["line"], "context": line, "score": h["context_score"]})
    return out


def _unique_values(bindings: list[dict]) -> list[float]:
    vals = []
    for b in bindings:
        v = b["value"]
        if not any(math.isclose(v, x, rel_tol=1e-12, abs_tol=1e-15) for x in vals):
            vals.append(v)
    return vals


def _classify(target: str, hits: list[dict]) -> dict:
    bindings = _explicit_binding_values(target, hits)
    vals = _unique_values(bindings)
    selected = None
    status = "NO_AUDITED_NUMERIC_BINDING_FOUND"
    if len(vals) == 1:
        selected = vals[0]
        status = "SINGLE_EXPLICIT_REPOSITORY_BINDING_FOUND_REVIEW_PROVENANCE"
    elif len(vals) > 1:
        status = "MULTIPLE_INCOMPATIBLE_EXPLICIT_BINDINGS_FOUND"

    # alpha_resolved has an extra hard guardrail: alpha_qm is not accepted as a
    # binding unless the exact alpha_resolved symbol is present in that context.
    if target == "alpha_resolved" and selected is not None:
        source_contexts = " ".join(b["context"].lower() for b in bindings if math.isclose(b["value"], selected, rel_tol=1e-12, abs_tol=1e-15))
        if "alpha_qm" in source_contexts and "alpha_resolved" not in source_contexts:
            selected = None
            status = "THERMAL_ALPHA_QM_NOT_PROMOTED_TO_ALPHA_RESOLVED"

    return {
        "status": status,
        "selected_value": selected,
        "all_explicit_values": vals,
        "explicit_bindings": bindings,
        "all_hits": hits,
    }


def main() -> None:
    repo_before = _repo_state()
    result_rows = {}
    for key, spec in TARGETS.items():
        hits = _extract_candidates(key, spec["aliases"])
        c = _classify(key, hits)
        c["role"] = spec["role"]
        result_rows[key] = c

    selected = {k: v["selected_value"] for k, v in result_rows.items()}
    all_closed = all(v is not None for v in selected.values())

    # Do not build E(a)/H(a) here. This lab is provenance only. Even a single
    # explicit binding is flagged for review rather than silently promoted into
    # the physical background.
    closure = {
        "status": "ALL_FIVE_HAVE_SINGLE_EXPLICIT_REPOSITORY_BINDINGS_REVIEW_REQUIRED" if all_closed else "BACKGROUND_PARAMETER_PROVENANCE_STILL_OPEN",
        "all_five_single_bindings_found": all_closed,
        "selected_values_for_provenance_review_only": selected,
        "physical_background_reconstruction_allowed": False,
        "safe_next": (
            "Review each single-binding source against authoritative V11/current provenance. "
            "Promote only proven values into a dedicated background reconstruction; resolve any multiple/missing bindings first."
        ),
    }

    repo_after = _repo_state()
    checks = {
        "no_lensing_target_used": True,
        "no_lcdm_distance_imported": True,
        "no_G_backsolve": True,
        "legacy_0p18_used": False,
        "fit_or_tuning_used": False,
        "quantum_engine_executed": False,
        "planck_scale_used": False,
        "gravity_fundamental_in_PBUF": False,
        "thermal_alpha_qm_not_auto_promoted": True,
        "no_E_or_H_reconstructed_from_unreviewed_values": True,
        "no_tracked_or_staged_changes": repo_after["tracked_changes"] == "" and repo_after["staged_changes"] == "",
        "stdout_only_no_run_directory_created": True,
    }

    payload = {
        "lab_id": LAB_ID,
        "status": "FACT_FINDING_ONLY",
        "repo_state": repo_after,
        "parameters": result_rows,
        "closure": closure,
        "checks": checks,
        "policy": {
            "gravity_fundamental_in_PBUF": False,
            "measured_G_role": "NOT_USED",
            "lensing_target_used": False,
            "fit_or_tuning_used": False,
            "alpha_qm_role": "THERMAL_METADATA_ONLY_NOT_ALPHA_RESOLVED_BY_ASSUMPTION",
        },
    }

    print(LAB_ID)
    print("status=FACT_FINDING_ONLY")
    print(f"head_sha={repo_after['head_sha']}")
    print("gravity_fundamental_in_PBUF=false")
    print("lensing_target_used=false")
    print("fit_or_tuning_used=false")
    print("thermal_alpha_qm_auto_promoted=false")
    print()
    print("PARAMETER_PROVENANCE")
    print("parameter | status | selected_value | hit_count | explicit_binding_count")
    for key, row in result_rows.items():
        print(f"{key} | {row['status']} | {row['selected_value']} | {len(row['all_hits'])} | {len(row['explicit_bindings'])}")
    print()
    print("CONCLUSION")
    print(f"status={closure['status']}")
    print(f"all_five_single_bindings_found={str(all_closed).lower()}")
    print("physical_background_reconstruction_allowed=false")
    print(f"safe_next={closure['safe_next']}")
    print()
    print("CHECKS")
    for k, v in checks.items():
        print(f"{k}={str(v).lower() if isinstance(v, bool) else v}")
    print("JSON=" + json.dumps(payload, sort_keys=True, separators=(",", ":")))


if __name__ == "__main__":
    main()
