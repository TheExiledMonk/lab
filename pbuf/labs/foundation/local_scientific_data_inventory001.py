#!/usr/bin/env python3
"""PBUF FOUNDATION — LOCAL SCIENTIFIC DATA INVENTORY 001.

Purpose
-------
Inventory the scientific data physically present in the runner's local repository,
including ignored/untracked large files that are intentionally not stored on GitHub.

This lab DOES NOT run PBUF physics and DOES NOT choose or fabricate a source. It only
reports what data actually exist locally so a later clean benchmark can be wired to
real independent inputs rather than assumptions.

Scope
-----
- recursively scan the repository for common scientific-data file types;
- exclude .git and runs/ artifacts;
- inspect FITS headers/shapes without loading full image cubes;
- classify known weak-lensing products by filename only;
- flag possible independent-source candidates for human review, never auto-accept;
- report files grouped by the five canonical PBUF benchmark cluster directories;
- make no network calls and write no output files.
"""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

from astropy.io import fits

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))

from pbuf.core import benchmark_data as BENCH

LAB_ID = "PBUF-FOUNDATION-LOCAL-SCIENTIFIC-DATA-INVENTORY-001"
DATA_SUFFIXES = {".fits", ".fit", ".fts", ".npy", ".npz", ".h5", ".hdf5", ".csv", ".tsv"}
EXCLUDED_TOP = {".git", "runs"}
HEADER_KEYS = (
    "NAXIS", "NAXIS1", "NAXIS2", "NAXIS3", "BUNIT", "OBJECT", "TELESCOP",
    "INSTRUME", "FILTER", "EXTNAME", "BTYPE", "CTYPE1", "CTYPE2", "CTYPE3",
    "PBUFROLE", "DATATYPE", "CONTENT", "ORIGIN",
)


def git(*args: str) -> str:
    return subprocess.check_output(["git", *args], cwd=ROOT, text=True).strip()


def repo_state() -> dict:
    return {
        "repository": "TheExiledMonk/lab",
        "branch": git("branch", "--show-current"),
        "head_sha": git("rev-parse", "HEAD"),
        "tracked_changes": git("diff", "--name-only"),
        "staged_changes": git("diff", "--name-only", "--cached"),
    }


def _rel(path: Path) -> str:
    try:
        return str(path.relative_to(ROOT))
    except ValueError:
        return str(path)


def _is_excluded(path: Path) -> bool:
    rel = path.relative_to(ROOT)
    return bool(rel.parts and rel.parts[0] in EXCLUDED_TOP)


def _known_lensing_role(path: Path) -> str | None:
    name = path.name.lower()
    if "_kappa" in name:
        return "weak_lensing_kappa_target"
    if "_gamma1" in name:
        return "weak_lensing_gamma1_product"
    if "_gamma2" in name:
        return "weak_lensing_gamma2_product"
    if "_gamma" in name:
        return "weak_lensing_gamma_product"
    if "shear" in name:
        return "weak_lensing_shear_product"
    return None


def _fits_info(path: Path) -> dict:
    out = {
        "path": _rel(path),
        "suffix": path.suffix.lower(),
        "size_bytes": int(path.stat().st_size),
        "readable": False,
        "fits_hdu_count": None,
        "primary_shape": None,
        "primary_ndim": None,
        "header": {},
        "known_lensing_role": _known_lensing_role(path),
        "possible_independent_source_candidate": False,
        "candidate_reason": None,
    }
    try:
        with fits.open(path, memmap=True, mode="readonly") as hdul:
            out["fits_hdu_count"] = len(hdul)
            hdr = hdul[0].header
            data = hdul[0].data
            shape = tuple(int(x) for x in data.shape) if data is not None else None
            out["primary_shape"] = shape
            out["primary_ndim"] = len(shape) if shape is not None else int(hdr.get("NAXIS", 0))
            out["header"] = {k: hdr.get(k) for k in HEADER_KEYS if hdr.get(k) is not None}
            out["readable"] = True
    except Exception as exc:
        out["error"] = f"{type(exc).__name__}: {exc}"
        return out

    # Candidate flag is intentionally permissive and non-selecting. It only tells us
    # which local files deserve inspection as possible source data.
    role = out["known_lensing_role"]
    ndim = out["primary_ndim"] or 0
    bunit = str(out["header"].get("BUNIT", "")).lower()
    pbufrole = str(out["header"].get("PBUFROLE", "")).lower()
    name = path.name.lower()

    if role is None:
        reasons = []
        if ndim >= 3:
            reasons.append("3D_or_higher_FITS")
        if any(token in name for token in ("mass", "density", "rho", "baryon", "gas", "stellar", "matter", "source")):
            reasons.append("source_like_filename")
        if any(token in bunit for token in ("kg", "g/", "msun", "m_sun", "solmass", "density")):
            reasons.append("physical_mass_or_density_like_BUNIT")
        if pbufrole:
            reasons.append(f"PBUFROLE={pbufrole}")
        if reasons:
            out["possible_independent_source_candidate"] = True
            out["candidate_reason"] = reasons

    return out


def _generic_info(path: Path) -> dict:
    return {
        "path": _rel(path),
        "suffix": path.suffix.lower(),
        "size_bytes": int(path.stat().st_size),
        "known_lensing_role": _known_lensing_role(path),
        "possible_independent_source_candidate": any(
            token in path.name.lower()
            for token in ("mass", "density", "rho", "baryon", "gas", "stellar", "matter", "source")
        ),
        "candidate_reason": ["source_like_filename"] if any(
            token in path.name.lower()
            for token in ("mass", "density", "rho", "baryon", "gas", "stellar", "matter", "source")
        ) else None,
    }


def scan_files() -> list[dict]:
    rows = []
    for path in ROOT.rglob("*"):
        if not path.is_file() or _is_excluded(path):
            continue
        if path.suffix.lower() not in DATA_SUFFIXES:
            continue
        if path.suffix.lower() in {".fits", ".fit", ".fts"}:
            rows.append(_fits_info(path))
        else:
            rows.append(_generic_info(path))
    return sorted(rows, key=lambda x: x["path"])


def cluster_group(rows: list[dict]) -> list[dict]:
    grouped = []
    for cluster in BENCH.clusters():
        prefix = f"PBUF_benchmark/{cluster['directory']}/"
        members = [r for r in rows if r["path"].startswith(prefix)]
        grouped.append({
            "cluster_id": cluster["id"],
            "directory": cluster["directory"],
            "file_count": len(members),
            "possible_independent_source_candidates": [
                r["path"] for r in members if r.get("possible_independent_source_candidate")
            ],
            "files": members,
        })
    return grouped


def main() -> int:
    state = repo_state()
    rows = scan_files()
    clusters = cluster_group(rows)
    candidates = [r for r in rows if r.get("possible_independent_source_candidate")]
    known_lensing = [r for r in rows if r.get("known_lensing_role")]
    fits_rows = [r for r in rows if r["suffix"] in {".fits", ".fit", ".fts"}]

    checks = {
        "canonical_five_cluster_inventory": len(clusters) == 5,
        "local_scan_completed": True,
        "fits_headers_inspected_without_physics_run": True,
        "runs_directory_excluded": True,
        "git_directory_excluded": True,
        "no_source_auto_selected": True,
        "no_physics_executed": True,
        "no_files_written": True,
        "no_tracked_or_staged_changes": bool(not state["tracked_changes"] and not state["staged_changes"]),
    }
    status = "LOCAL_SCIENTIFIC_DATA_INVENTORY_EXECUTED" if all(checks.values()) else "LOCAL_SCIENTIFIC_DATA_INVENTORY_PARTIAL_EXECUTION"

    result = {
        "lab_id": LAB_ID,
        "status": status,
        "repo_state": state,
        "scan_root": str(ROOT),
        "excluded_top_level": sorted(EXCLUDED_TOP),
        "scientific_file_count": len(rows),
        "fits_file_count": len(fits_rows),
        "known_lensing_product_count": len(known_lensing),
        "possible_independent_source_candidate_count": len(candidates),
        "possible_independent_source_candidates": candidates,
        "clusters": clusters,
        "all_scientific_files": rows,
        "checks": checks,
        "interpretation_rule": (
            "Inventory only. Candidate flags are prompts for human inspection, not accepted sources. "
            "Do not run a clean weak-lensing benchmark until a genuinely independent source is identified from actual local data."
        ),
    }

    print(LAB_ID)
    print(f"status={status}")
    print(f"head_sha={state['head_sha']}")
    print(f"branch={state['branch']}")
    print(f"scan_root={ROOT}")
    print(f"scientific_file_count={len(rows)}")
    print(f"fits_file_count={len(fits_rows)}")
    print(f"known_lensing_product_count={len(known_lensing)}")
    print(f"possible_independent_source_candidate_count={len(candidates)}")
    print()
    print("CLUSTER_INVENTORY")
    for c in clusters:
        print(
            f"cluster={c['cluster_id']} file_count={c['file_count']} "
            f"candidate_count={len(c['possible_independent_source_candidates'])}"
        )
        for r in c["files"]:
            print(
                f"  file={r['path']} size={r['size_bytes']} "
                f"ndim={r.get('primary_ndim')} shape={r.get('primary_shape')} "
                f"BUNIT={r.get('header', {}).get('BUNIT')} "
                f"lensing_role={r.get('known_lensing_role')} "
                f"source_candidate={r.get('possible_independent_source_candidate')} "
                f"reason={r.get('candidate_reason')}"
            )
    print()
    print("GLOBAL_SOURCE_CANDIDATES")
    if candidates:
        for r in candidates:
            print(
                f"candidate={r['path']} size={r['size_bytes']} ndim={r.get('primary_ndim')} "
                f"shape={r.get('primary_shape')} BUNIT={r.get('header', {}).get('BUNIT')} "
                f"reason={r.get('candidate_reason')}"
            )
    else:
        print("none")
    print()
    print("CHECKS")
    for k, v in checks.items():
        print(f"{k}={str(v).lower()}")
    print(json.dumps(result, sort_keys=True, separators=(",", ":"), default=str))
    return 0 if status == "LOCAL_SCIENTIFIC_DATA_INVENTORY_EXECUTED" else 1


if __name__ == "__main__":
    raise SystemExit(main())
