#!/usr/bin/env python3
"""PBUF FOUNDATION — BACKGROUND DISTANCE RECOVERY AUDIT 001.

Fact-finding only.

Recover as much of the historical V11 PBUF background-expansion chain as the
current repository can reproduce without importing LambdaCDM or filling missing
parameters by hand.

Historical V11 chain under audit:

    k_max(a) = epsilon0_T(a) - alpha_T(a)
    decay(a) = exp(-a/Rmax)
    S(a) = 1 - (1-k_max(a))*decay(a)
    Omega_sigma_raw(a) = alpha_T(a)*(1-decay(a))*S(a)
    Omega_sigma_target = 1 - Omega_m0 - Omega_r0 - alpha_resolved
    sigma_rescale = Omega_sigma_target/Omega_sigma_raw(1)
    Omega_sigma(a) = sigma_rescale*Omega_sigma_raw(a)
    E(a)^2 = Omega_m0*a^-3 + Omega_r0*a^-4 + Omega_sigma(a)
    H(a) = H0*E(a)

The audit does not attach this background to the cluster redshifts.  The prior
guardrail remains: observed spectral redshift is not identified with expansion
redshift by assumption.
"""
from __future__ import annotations

import hashlib
import json
import math
import re
import subprocess
import sys
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))

LAB_ID = "PBUF-FOUNDATION-BACKGROUND-DISTANCE-RECOVERY-AUDIT-001"
THERMAL_PATH = ROOT / "pbuf" / "data" / "quantum" / "thermal_table_cache.json"
V11_AUDIT_PATH = ROOT / "v11_alpha_audit.py"
ALG_TOL = 1.0e-12


def _git(*args: str) -> str:
    return subprocess.check_output(["git", *args], cwd=str(ROOT), text=True).strip()


def _repo_state() -> dict[str, Any]:
    return {
        "repository": "TheExiledMonk/lab",
        "branch": _git("rev-parse", "--abbrev-ref", "HEAD"),
        "head_sha": _git("rev-parse", "HEAD"),
        "tracked_changes": _git("diff", "--name-only"),
        "staged_changes": _git("diff", "--name-only", "--cached"),
    }


def _sha256(path: Path) -> str | None:
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else None


def _load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text())


def _thermal_inventory() -> dict[str, Any]:
    if not THERMAL_PATH.exists():
        return {"present": False, "status": "MISSING", "path": str(THERMAL_PATH.relative_to(ROOT))}

    payload = _load_json(THERMAL_PATH)
    meta = payload.get("metadata", {})
    rows = payload.get("rows", [])
    required = ("a", "T_K", "epsilon0_T", "alpha_T")
    valid = []
    for row in rows:
        try:
            vals = [float(row[k]) for k in required]
        except (KeyError, TypeError, ValueError):
            continue
        if all(np.isfinite(vals)) and vals[0] > 0.0:
            valid.append(row)

    a = np.asarray([float(r["a"]) for r in valid], dtype=np.float64)
    eps = np.asarray([float(r["epsilon0_T"]) for r in valid], dtype=np.float64)
    alpha_t = np.asarray([float(r["alpha_T"]) for r in valid], dtype=np.float64)
    order = np.argsort(a) if a.size else np.array([], dtype=int)
    a = a[order] if a.size else a
    eps = eps[order] if eps.size else eps
    alpha_t = alpha_t[order] if alpha_t.size else alpha_t

    return {
        "present": True,
        "status": "FOUND_AND_NUMERICALLY_READABLE" if valid else "FOUND_BUT_NO_VALID_ROWS",
        "path": str(THERMAL_PATH.relative_to(ROOT)),
        "sha256": _sha256(THERMAL_PATH),
        "row_count": len(rows),
        "valid_row_count": len(valid),
        "a_min": float(a.min()) if a.size else None,
        "a_max": float(a.max()) if a.size else None,
        "covers_a_equals_1": bool(a.size and a[0] <= 1.0 <= a[-1]),
        "epsilon0_T_min": float(eps.min()) if eps.size else None,
        "epsilon0_T_max": float(eps.max()) if eps.size else None,
        "alpha_T_min": float(alpha_t.min()) if alpha_t.size else None,
        "alpha_T_max": float(alpha_t.max()) if alpha_t.size else None,
        "metadata": {
            k: meta.get(k) for k in (
                "mode", "table_version", "method_version", "regulator", "field_content",
                "alpha_qm", "eps0_base", "micro_hash", "micro_source", "generated_at"
            )
        },
    }


def _v11_trace_inventory() -> dict[str, Any]:
    if not V11_AUDIT_PATH.exists():
        return {"present": False, "status": "MISSING", "markers": {}}
    text = V11_AUDIT_PATH.read_text(errors="ignore")
    markers = {
        "kmax": "k_max(a)=epsilon_0,T(a)-alpha_T(a)" in text,
        "omega_sigma_raw": "Omega_sigma_raw(a)=alpha_T(a)(1-decay(a))S(a)" in text,
        "flat_today": "Omega_sigma_target=1-Omega_m0-Omega_r0-alpha_resolved" in text,
        "background": "Omega_sigma(a), Omega_m0, Omega_r0" in text and "E(a), H(a)" in text,
        "baryon_identity": "Omega_b0=2 alpha_resolved" in text,
    }
    return {
        "present": True,
        "status": "FOUND_WITH_V11_TRACE_MARKERS",
        "path": str(V11_AUDIT_PATH.relative_to(ROOT)),
        "sha256": _sha256(V11_AUDIT_PATH),
        "markers": markers,
        "all_required_trace_markers_present": all(markers.values()),
    }


def _iter_text_files() -> list[Path]:
    skip = {".git", "runs", "__pycache__", ".venv", "venv"}
    return [
        p for p in ROOT.rglob("*")
        if p.is_file()
        and not any(part in skip for part in p.parts)
        and p.suffix.lower() in {".py", ".json", ".md", ".txt", ".toml", ".yaml", ".yml"}
    ]


def _numeric_candidates(text: str, symbol: str) -> list[float]:
    pattern = re.compile(rf"\b{re.escape(symbol)}\b\s*(?:=|:)\s*([-+]?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?)")
    vals = []
    for match in pattern.finditer(text):
        try:
            value = float(match.group(1))
        except ValueError:
            continue
        if np.isfinite(value):
            vals.append(value)
    return vals


def _parameter_inventory() -> dict[str, Any]:
    aliases = {
        "Rmax": ("Rmax", "R_max"),
        "H0": ("H0", "H_0"),
        "Omega_r0": ("Omega_r0", "OMEGA_R0"),
        "alpha_resolved": ("alpha_resolved", "ALPHA_RESOLVED"),
        "BARYON_FRACTION": ("BARYON_FRACTION",),
    }
    hits: dict[str, list[dict[str, Any]]] = {k: [] for k in aliases}
    for path in _iter_text_files():
        try:
            text = path.read_text(errors="ignore")
        except Exception:
            continue
        rel = str(path.relative_to(ROOT))
        for key, names in aliases.items():
            if not any(name in text for name in names):
                continue
            values = []
            for name in names:
                values.extend(_numeric_candidates(text, name))
            hits[key].append({"path": rel, "numeric_candidates": sorted(set(values))})

    thermal_alpha_qm = None
    if THERMAL_PATH.exists():
        thermal_alpha_qm = _load_json(THERMAL_PATH).get("metadata", {}).get("alpha_qm")

    # Inventory does not select among candidates.  A value becomes physical only
    # after a later audit establishes its authoritative/current provenance.
    return {
        "hits": hits,
        "thermal_alpha_qm": thermal_alpha_qm,
        "thermal_alpha_qm_promoted_to_alpha_resolved": False,
        "selected_values": {
            "Rmax": None,
            "H0": None,
            "Omega_r0": None,
            "alpha_resolved": None,
            "BARYON_FRACTION": None,
        },
    }


def _algebra_control() -> dict[str, Any]:
    """Synthetic wiring test only; values below have no physical meaning."""
    a = np.array([0.25, 0.5, 1.0, 2.0], dtype=np.float64)
    eps = np.array([0.91, 0.92, 0.93, 0.94], dtype=np.float64)
    alpha_t = np.array([0.020, 0.021, 0.022, 0.023], dtype=np.float64)
    rmax, om_m, om_r, alpha_resolved, h0 = 3.5, 0.31, 9.0e-5, 0.0219, 70.0

    kmax = eps - alpha_t
    decay = np.exp(-a / rmax)
    sat = 1.0 - (1.0 - kmax) * decay
    raw = alpha_t * (1.0 - decay) * sat
    target = 1.0 - om_m - om_r - alpha_resolved
    today = int(np.where(a == 1.0)[0][0])
    rescale = target / raw[today]
    omega_sigma = rescale * raw
    e2 = om_m / a**3 + om_r / a**4 + omega_sigma
    hz = h0 * np.sqrt(e2)

    direct_raw = np.array([
        alpha_t[i]
        * (1.0 - math.exp(-a[i] / rmax))
        * (1.0 - (1.0 - (eps[i] - alpha_t[i])) * math.exp(-a[i] / rmax))
        for i in range(a.size)
    ])
    direct_e2 = np.array([
        om_m * a[i] ** -3 + om_r * a[i] ** -4 + rescale * direct_raw[i]
        for i in range(a.size)
    ])

    raw_err = float(np.max(np.abs(raw - direct_raw)))
    e2_err = float(np.max(np.abs(e2 - direct_e2)))
    flat_err = float(abs(omega_sigma[today] - target))
    passed = bool(
        raw_err <= ALG_TOL
        and e2_err <= ALG_TOL
        and flat_err <= ALG_TOL
        and np.all(np.isfinite(hz))
        and np.all(hz > 0.0)
    )
    return {
        "synthetic_only_not_physical_model": True,
        "raw_vector_vs_scalar_max_abs_error": raw_err,
        "e2_vector_vs_scalar_max_abs_error": e2_err,
        "flat_today_elastic_value_error": flat_err,
        "H_finite_positive": bool(np.all(np.isfinite(hz)) and np.all(hz > 0.0)),
        "pass": passed,
    }


def _dependency_ledger(thermal: dict[str, Any], trace: dict[str, Any], params: dict[str, Any]) -> list[dict[str, Any]]:
    selected = params["selected_values"]
    return [
        {
            "key": "thermal_table_alpha_T_epsilon0_T",
            "status": "FOUND_AND_REPRODUCIBLE" if thermal.get("valid_row_count", 0) > 0 else "MISSING",
            "blocks_background_closure": not (thermal.get("valid_row_count", 0) > 0),
        },
        {
            "key": "v11_equation_trace",
            "status": "FOUND_AND_REPRODUCIBLE" if trace.get("all_required_trace_markers_present") else "INCOMPLETE",
            "blocks_background_closure": not bool(trace.get("all_required_trace_markers_present")),
        },
        {
            "key": "Rmax_activation_scale",
            "status": "MISSING_AUDITED_NUMERIC_VALUE" if selected["Rmax"] is None else "FOUND_AND_REPRODUCIBLE",
            "blocks_background_closure": selected["Rmax"] is None,
        },
        {
            "key": "alpha_resolved",
            "status": "MISSING_AUDITED_NUMERIC_VALUE" if selected["alpha_resolved"] is None else "FOUND_AND_REPRODUCIBLE",
            "blocks_background_closure": selected["alpha_resolved"] is None,
            "guardrail": "thermal alpha_qm not silently identified with alpha_resolved",
        },
        {
            "key": "Omega_r0",
            "status": "MISSING_AUDITED_NUMERIC_VALUE" if selected["Omega_r0"] is None else "FOUND_AND_REPRODUCIBLE",
            "blocks_background_closure": selected["Omega_r0"] is None,
        },
        {
            "key": "BARYON_FRACTION_to_Omega_m0",
            "status": "MISSING_AUDITED_NUMERIC_VALUE" if selected["BARYON_FRACTION"] is None else "FOUND_AND_REPRODUCIBLE",
            "blocks_background_closure": selected["BARYON_FRACTION"] is None,
        },
        {
            "key": "H0_absolute_distance_scale",
            "status": "MISSING_AUDITED_NUMERIC_VALUE" if selected["H0"] is None else "FOUND_AND_REPRODUCIBLE",
            "blocks_background_closure": selected["H0"] is None,
        },
        {
            "key": "flat_today_normalization_rule",
            "status": "HISTORICAL_EQUATION_TRACE_PRESENT" if trace.get("markers", {}).get("flat_today") else "MISSING",
            "blocks_background_closure": not bool(trace.get("markers", {}).get("flat_today")),
        },
    ]


def main() -> None:
    thermal = _thermal_inventory()
    trace = _v11_trace_inventory()
    params = _parameter_inventory()
    algebra = _algebra_control()
    ledger = _dependency_ledger(thermal, trace, params)
    blocking = [row["key"] for row in ledger if row["blocks_background_closure"]]
    background_closed = (not blocking) and algebra["pass"]
    status = (
        "PBUF_BACKGROUND_DEPENDENCIES_RECOVERED_READY_FOR_DISTANCE_IMPLEMENTATION"
        if background_closed
        else "PBUF_BACKGROUND_RECOVERY_PARTIAL"
    )
    repo = _repo_state()

    checks = {
        "thermal_cache_readable": bool(thermal.get("valid_row_count", 0) > 0),
        "thermal_alpha_qm_not_promoted_to_alpha_resolved": not params["thermal_alpha_qm_promoted_to_alpha_resolved"],
        "v11_equation_algebra_pass": bool(algebra["pass"]),
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
        "no_tracked_or_staged_changes": repo["tracked_changes"] == "" and repo["staged_changes"] == "",
        "stdout_only_no_run_directory_created": True,
    }

    result = {
        "lab_id": LAB_ID,
        "status": "FACT_FINDING_ONLY",
        "repo_state": repo,
        "thermal_cache": thermal,
        "v11_trace": trace,
        "parameter_inventory": params,
        "algebra_control": algebra,
        "dependency_ledger": ledger,
        "closure": {
            "status": status,
            "background_dependencies_closed": background_closed,
            "blocking_dependencies": blocking,
            "E_of_a_recovered_physical": bool(background_closed),
            "H_of_a_recovered_physical": bool(background_closed),
            "distance_integration_allowed": bool(background_closed),
            "safe_next": (
                "Recover authoritative/current numerical provenance for each remaining background parameter one at a time. "
                "Only after every dependency closes should a background-only E(a), H(a), and radial-distance integral be implemented; "
                "do not apply observed cluster redshift until redshift decomposition is separately audited."
            ),
        },
        "checks": checks,
    }

    print(LAB_ID)
    print("status=FACT_FINDING_ONLY")
    print(f"head_sha={repo['head_sha']}")
    print("gravity_fundamental_in_PBUF=false")
    print("cluster_redshift_used_as_expansion_redshift=false")
    print("lcdm_distance_imported=false")
    print("lensing_target_used=false")
    print("fit_or_tuning_used=false")
    print()

    print("DEPENDENCY_LEDGER")
    print("key | status | blocks_background_closure")
    for row in ledger:
        print(f"{row['key']} | {row['status']} | {row['blocks_background_closure']}")
    print()

    print("THERMAL_CACHE")
    print(f"status={thermal.get('status')}")
    print(f"valid_row_count={thermal.get('valid_row_count')}")
    print(f"a_range={thermal.get('a_min')}..{thermal.get('a_max')}")
    print(f"covers_a_equals_1={thermal.get('covers_a_equals_1')}")
    print(f"thermal_alpha_qm={params.get('thermal_alpha_qm')}")
    print("thermal_alpha_qm_promoted_to_alpha_resolved=false")
    print()

    print("ALGEBRA_CONTROL")
    print("synthetic_only_not_physical_model=true")
    print(f"raw_vector_vs_scalar_max_abs_error={algebra['raw_vector_vs_scalar_max_abs_error']:.17e}")
    print(f"e2_vector_vs_scalar_max_abs_error={algebra['e2_vector_vs_scalar_max_abs_error']:.17e}")
    print(f"flat_today_elastic_value_error={algebra['flat_today_elastic_value_error']:.17e}")
    print(f"pass={algebra['pass']}")
    print()

    print("CONCLUSION")
    print(f"status={status}")
    print(f"background_dependencies_closed={background_closed}")
    print("blocking_dependencies=" + (",".join(blocking) if blocking else "NONE"))
    print(f"distance_integration_allowed={background_closed}")
    print("safe_next=" + result["closure"]["safe_next"])
    print()

    print("CHECKS")
    for key, value in checks.items():
        print(f"{key}={str(value).lower() if isinstance(value, bool) else value}")
    print("JSON=" + json.dumps(result, sort_keys=True, separators=(",", ":")))


if __name__ == "__main__":
    main()
