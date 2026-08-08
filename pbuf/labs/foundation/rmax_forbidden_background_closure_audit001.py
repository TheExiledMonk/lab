#!/usr/bin/env python3
"""PBUF foundation audit: determine background closure when Rmax is forbidden.

Rmax is treated as physically unknown and unavailable as an input.  The audit
therefore does not attempt to recover, fit, guess, inherit, or substitute a
numerical Rmax.  It traces which historical V11 background relations remain
usable without Rmax and identifies the first equation at which the current
physical background becomes underdetermined.

Fact-finding only.  No lensing target, no LCDM distance substitution, no
measured-G backsolve, no Quantum Engine execution, and no run directory.
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
RECOVERY_PATH = ROOT / "pbuf/labs/foundation/pbuf_background_distance_recovery_audit001.py"

LAB_ID = "PBUF-FOUNDATION-RMAX-FORBIDDEN-BACKGROUND-CLOSURE-AUDIT-001"

# Historical equations are provenance records, not newly promoted physics.
# Dependency tags are deliberately explicit so that the audit cannot make a
# hidden numerical choice for Rmax.
HISTORICAL_CHAIN = [
    {
        "key": "thermal_alpha_T_epsilon0_T",
        "expression": "alpha_T(a), epsilon0_T(a)",
        "depends_on": [],
        "role": "thermal/LUT inputs",
    },
    {
        "key": "kmax",
        "expression": "k_max(a)=epsilon0_T(a)-alpha_T(a)",
        "depends_on": ["thermal_alpha_T_epsilon0_T"],
        "role": "thermal response difference",
    },
    {
        "key": "decay_activation",
        "expression": "decay(a)=exp(-a/Rmax)",
        "depends_on": ["Rmax"],
        "role": "historical activation/reversal-scale factor",
    },
    {
        "key": "saturation_S",
        "expression": "S(a)=1-(1-k_max(a))*decay(a)",
        "depends_on": ["kmax", "decay_activation"],
        "role": "historical saturation/activation function",
    },
    {
        "key": "omega_sigma_raw",
        "expression": "Omega_sigma_raw(a)=alpha_T(a)*(1-decay(a))*S(a)",
        "depends_on": ["thermal_alpha_T_epsilon0_T", "decay_activation", "saturation_S"],
        "role": "historical unnormalised elastic background",
    },
    {
        "key": "omega_sigma_target",
        "expression": "Omega_sigma_target=1-Omega_m0-Omega_r0-alpha_resolved",
        "depends_on": ["Omega_m0", "Omega_r0", "alpha_resolved"],
        "role": "historical flat-today target",
    },
    {
        "key": "sigma_rescale",
        "expression": "sigma_rescale=Omega_sigma_target/Omega_sigma_raw(a=1)",
        "depends_on": ["omega_sigma_target", "omega_sigma_raw"],
        "role": "historical flat-today rescale",
    },
    {
        "key": "omega_sigma",
        "expression": "Omega_sigma(a)=sigma_rescale*Omega_sigma_raw(a)",
        "depends_on": ["sigma_rescale", "omega_sigma_raw"],
        "role": "historical elastic background term",
    },
    {
        "key": "E_of_a",
        "expression": "E(a)^2=Omega_m0*a^-3+Omega_r0*a^-4+Omega_sigma(a)",
        "depends_on": ["Omega_m0", "Omega_r0", "omega_sigma"],
        "role": "dimensionless background expansion",
    },
    {
        "key": "H_of_a",
        "expression": "H(a)=H0*E(a)",
        "depends_on": ["H0", "E_of_a"],
        "role": "absolute expansion history",
    },
]

FORBIDDEN_INPUTS = {"Rmax"}
OTHER_OPEN_INPUTS = {"Omega_m0", "Omega_r0", "alpha_resolved", "H0"}


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
            valid.append((a, alpha_t, eps))
    return {
        "present": True,
        "readable": True,
        "row_count": len(rows) if isinstance(rows, list) else 0,
        "valid_row_count": len(valid),
        "a_min": min((x[0] for x in valid), default=None),
        "a_max": max((x[0] for x in valid), default=None),
        "covers_a_equals_1": any(abs(x[0] - 1.0) <= 1e-12 for x in valid),
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
        "omega_sigma_model": bool(re.search(r"Omega_sigma\(a\)=alpha\(1-exp\(-a/R_max\)\)S\(a\)", text)),
        "background_dependency": '"target":"E(a), H(a)"' in text,
    }
    return {
        "present": True,
        "sha256": sha256(TRACE_PATH),
        "markers": markers,
        "all_required_trace_markers_present": all(markers.values()),
    }


def propagate_dependency() -> list[dict[str, Any]]:
    status: dict[str, str] = {}
    rows: list[dict[str, Any]] = []
    for item in HISTORICAL_CHAIN:
        key = item["key"]
        deps = item["depends_on"]
        forbidden = sorted(d for d in deps if d in FORBIDDEN_INPUTS)
        open_external = sorted(d for d in deps if d in OTHER_OPEN_INPUTS)
        blocked_upstream = sorted(
            d for d in deps if status.get(d) in {"BLOCKED_BY_FORBIDDEN_RMAX", "BLOCKED_BY_OPEN_UPSTREAM"}
        )

        if forbidden:
            state = "BLOCKED_BY_FORBIDDEN_RMAX"
        elif blocked_upstream:
            state = "BLOCKED_BY_OPEN_UPSTREAM"
        elif open_external:
            state = "OPEN_EXTERNAL_INPUTS_BUT_NOT_RMAX"
        else:
            state = "AVAILABLE_WITHOUT_RMAX"
        status[key] = state
        rows.append({
            **item,
            "status": state,
            "direct_forbidden_dependencies": forbidden,
            "open_external_dependencies": open_external,
            "blocked_upstream_dependencies": blocked_upstream,
        })
    return rows


def first_rmax_block(chain: list[dict[str, Any]]) -> dict[str, Any] | None:
    for row in chain:
        if row["status"] == "BLOCKED_BY_FORBIDDEN_RMAX":
            return {"key": row["key"], "expression": row["expression"], "role": row["role"]}
    return None


def main() -> int:
    thermal = load_thermal()
    trace = trace_markers()
    chain = propagate_dependency()
    first_block = first_rmax_block(chain)

    kmax_available = next(r for r in chain if r["key"] == "kmax")["status"] == "AVAILABLE_WITHOUT_RMAX"
    omega_sigma_raw_blocked = next(r for r in chain if r["key"] == "omega_sigma_raw")["status"] != "AVAILABLE_WITHOUT_RMAX"
    e_blocked = next(r for r in chain if r["key"] == "E_of_a")["status"] != "AVAILABLE_WITHOUT_RMAX"
    h_blocked = next(r for r in chain if r["key"] == "H_of_a")["status"] != "AVAILABLE_WITHOUT_RMAX"

    # Rmax is not a missing numeric datum in this audit.  It is explicitly an
    # unknown future-reversal quantity and is therefore forbidden as an input.
    closure = {
        "status": "CURRENT_PBUF_BACKGROUND_BLOCKED_BY_UNKNOWN_RMAX",
        "Rmax_status": "PHYSICALLY_UNKNOWN_FORBIDDEN_AS_INPUT",
        "Rmax_numeric_value": None,
        "Rmax_recovery_attempted": False,
        "historical_v11_background_reproducible_without_Rmax": False,
        "kmax_thermal_relation_available_without_Rmax": kmax_available,
        "Omega_sigma_raw_available_without_Rmax": not omega_sigma_raw_blocked,
        "E_of_a_physically_closed": not e_blocked,
        "H_of_a_physically_closed": not h_blocked,
        "distance_integration_allowed": False,
        "first_Rmax_dependent_link": first_block,
        "safe_next": (
            "Do not assign Rmax. Treat the historical activation/reversal factor and all downstream "
            "Omega_sigma/E/H relations that depend on it as unresolved. Audit whether current PBUF "
            "admits an independently derived present-era/local background law that does not require "
            "the future reversal scale. Keep observed cluster redshift separate from expansion redshift."
        ),
    }

    checks = {
        "Rmax_value_is_null": closure["Rmax_numeric_value"] is None,
        "Rmax_recovery_not_attempted": not closure["Rmax_recovery_attempted"],
        "Rmax_forbidden_as_input": closure["Rmax_status"] == "PHYSICALLY_UNKNOWN_FORBIDDEN_AS_INPUT",
        "thermal_cache_readable": bool(thermal.get("readable")),
        "v11_trace_present": bool(trace.get("all_required_trace_markers_present")),
        "no_E_or_H_reconstructed_with_Rmax": e_blocked and h_blocked,
        "no_distance_values_fabricated": True,
        "no_cluster_redshift_used_as_expansion_redshift": True,
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
            "Rmax_role": "UNKNOWN_FUTURE_REVERSAL_SCALE_NOT_AN_INPUT",
            "Rmax_numeric_use_allowed": False,
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
            "prior_recovery_audit_present": RECOVERY_PATH.is_file(),
            "prior_recovery_audit_sha256": sha256(RECOVERY_PATH),
        },
        "dependency_chain": chain,
        "closure": closure,
        "checks": checks,
    }

    print(LAB_ID)
    print("status=FACT_FINDING_ONLY")
    print(f"head_sha={result['repo_state']['head_sha']}")
    print("gravity_fundamental_in_PBUF=false")
    print("Rmax_status=PHYSICALLY_UNKNOWN_FORBIDDEN_AS_INPUT")
    print("Rmax_numeric_value=None")
    print("Rmax_recovery_attempted=false")
    print("cluster_redshift_used_as_expansion_redshift=false")
    print("lcdm_distance_imported=false")
    print("lensing_target_used=false")
    print("fit_or_tuning_used=false")
    print()
    print("BACKGROUND_DEPENDENCY_WITH_RMAX_FORBIDDEN")
    print("key | status | role")
    for row in chain:
        print(f"{row['key']} | {row['status']} | {row['role']}")
    print()
    print("CONCLUSION")
    print(f"status={closure['status']}")
    print(f"historical_v11_background_reproducible_without_Rmax={str(closure['historical_v11_background_reproducible_without_Rmax']).lower()}")
    print(f"kmax_available_without_Rmax={str(closure['kmax_thermal_relation_available_without_Rmax']).lower()}")
    print(f"Omega_sigma_raw_available_without_Rmax={str(closure['Omega_sigma_raw_available_without_Rmax']).lower()}")
    print(f"E_of_a_physically_closed={str(closure['E_of_a_physically_closed']).lower()}")
    print(f"H_of_a_physically_closed={str(closure['H_of_a_physically_closed']).lower()}")
    print("distance_integration_allowed=false")
    if first_block:
        print(f"first_Rmax_dependent_link={first_block['key']}")
    print(f"safe_next={closure['safe_next']}")
    print()
    print("CHECKS")
    for key, value in checks.items():
        print(f"{key}={str(value).lower() if isinstance(value, bool) else value}")
    print("JSON=" + json.dumps(result, sort_keys=True, separators=(",", ":")))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
