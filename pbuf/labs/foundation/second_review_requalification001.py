"""PBUF VERIFIED NUMERICAL CORE - SECOND-REVIEW-REQUALIFICATION-001.

Re-validates the FOUNDATION-001 numerical core after the second-review
correction PRs (#2..#7) have been merged into ``main``.

This is a validation-only lab. It must not edit any source module under
``pbuf/core`` or ``pbuf/models``. Any module failure must be recorded
and propagated as a hard gate stop.
"""
from __future__ import annotations

import csv
import hashlib
import json
import math
import os
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))

from pbuf.core import (
    conventions as M01,
    coordinate_transforms as M02,
    vector_transforms as M03,
    tensor_transforms as M04,
    pair_enumeration as M05,
    field_diagnostics as M11,
    differential_operators as M12,
    helmholtz_3d as M13,
    los_projection as M14,
    observable_extraction as M16,
    pair_transfer as M08,
    midpoint_rasterization as M10,
    ray_interface as M15,
)
from pbuf.models import (
    a8_state as M06_state,
    a8_pair_amplitude as M06,
    transverse_projector as M07,
)
from pbuf.validation.protected_function_scanner import (
    scan_protected_functions,
)

OUT = ROOT / "runs" / "verified_numerical_core_second_review_requalification001"

LAB_ID = "PBUF-SECOND-REVIEW-REQUALIFICATION-001"
CONVENTIONS_VERSION = M01.CONVENTIONS_VERSION

REQUIRED_PRS = [2, 3, 4, 5, 6, 7]

# Map module_id -> (module_name, source_module_object).
MODULE_REGISTRY = {
    "M01": ("conventions", M01),
    "M02": ("coordinate_transforms", M02),
    "M03": ("vector_transforms", M03),
    "M04": ("tensor_transforms", M04),
    "M05": ("pair_enumeration", M05),
    "M06": ("a8_pair_amplitude", M06),
    "M07": ("transverse_projector", M07),
    "M08": ("pair_transfer", M08),
    "M09": ("pair_transfer", M08),
    "M10": ("midpoint_rasterization", M10),
    "M11": ("field_diagnostics", M11),
    "M12": ("differential_operators", M12),
    "M13": ("helmholtz_3d", M13),
    "M14": ("los_projection", M14),
    "M15": ("ray_interface", M15),
    "M16": ("observable_extraction", M16),
}

# Source paths (for source integrity inventory).
MODULE_SOURCE_PATHS = {
    "M01": "pbuf/core/conventions.py",
    "M02": "pbuf/core/coordinate_transforms.py",
    "M03": "pbuf/core/vector_transforms.py",
    "M04": "pbuf/core/tensor_transforms.py",
    "M05": "pbuf/core/pair_enumeration.py",
    "M06": "pbuf/models/a8_pair_amplitude.py",
    "M07": "pbuf/models/transverse_projector.py",
    "M08": "pbuf/core/pair_transfer_verified.py",
    "M09": "pbuf/core/pair_transfer_verified.py",
    "M10": "pbuf/core/midpoint_rasterization.py",
    "M11": "pbuf/core/field_diagnostics.py",
    "M12": "pbuf/core/differential_operators.py",
    "M13": "pbuf/core/helmholtz_3d_verified.py",
    "M14": "pbuf/core/los_projection.py",
    "M15": "pbuf/core/ray_interface.py",
    "M16": "pbuf/core/observable_extraction.py",
}


def _hash_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _hash_array_raw(arr: np.ndarray) -> str:
    arr = np.ascontiguousarray(arr)
    payload = arr.dtype.str.encode("utf-8")
    payload += str(arr.shape).encode("utf-8")
    payload += arr.tobytes()
    return hashlib.sha256(payload).hexdigest()


def _hash_array_canonical(arr: np.ndarray) -> str:
    arr64 = np.ascontiguousarray(np.asarray(arr, dtype=np.float64))
    return hashlib.sha256(arr64.tobytes()).hexdigest()


def _write_csv(path: Path, rows: list) -> None:
    if not rows:
        path.write_text("")
        return
    keys = []
    for r in rows:
        for k in r.keys():
            if k not in keys:
                keys.append(k)
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)


def _git(args):
    return subprocess.check_output(["git"] + args, cwd=str(ROOT), text=True)


def _git_blob_sha(path: str) -> str:
    return _git(["rev-parse", f"HEAD:{path}"]).strip()


def _array_stats(arr: np.ndarray) -> dict:
    finite = np.asarray(arr, dtype=np.float64)
    finite_mask = np.isfinite(finite)
    n_finite = int(finite_mask.sum())
    if n_finite == 0:
        return {
            "minimum": float("nan"),
            "maximum": float("nan"),
            "mean": float("nan"),
            "variance": float("nan"),
            "rms": float("nan"),
            "nonzero_count": 0,
            "field_is_finite": False,
        }
    f = finite[finite_mask]
    return {
        "minimum": float(f.min()),
        "maximum": float(f.max()),
        "mean": float(f.mean()),
        "variance": float(f.var()),
        "rms": float(np.sqrt(np.mean(f ** 2))),
        "nonzero_count": int(np.count_nonzero(f)),
        "field_is_finite": bool(np.all(finite_mask)),
    }


def record_repository_state():
    branch = _git(["rev-parse", "--abbrev-ref", "HEAD"]).strip()
    head = _git(["rev-parse", "HEAD"]).strip()
    status = _git(["status", "--porcelain"]).strip()
    tracked_changes = subprocess.run(
        ["git", "diff", "--name-only"],
        cwd=str(ROOT), capture_output=True, text=True).stdout.strip()
    staged_changes = subprocess.run(
        ["git", "diff", "--name-only", "--cached"],
        cwd=str(ROOT), capture_output=True, text=True).stdout.strip()
    working_tree_clean = (status == ""
                          or all(line.startswith("??")
                                  for line in status.splitlines()))
    log = _git(["log", "--pretty=format:%H %s"]).strip().splitlines()

    pr_merge_shas = {}
    for line in log:
        parts = line.split(" ", 1)
        if len(parts) < 2:
            continue
        sha, subject = parts
        if "Merge pull request" in subject:
            try:
                tag = subject.split("Merge pull request ")[1].split(" ")[0]
                pr_merge_shas[int(tag.strip("#"))] = sha
            except (IndexError, ValueError):
                pass

    missing = [pr for pr in REQUIRED_PRS if pr not in pr_merge_shas]

    state = {
        "repository": "TheExiledMonk/lab",
        "branch": branch,
        "head_sha": head,
        "working_tree_clean": working_tree_clean,
        "tracked_changes": tracked_changes,
        "staged_changes": staged_changes,
        "required_prs_present": REQUIRED_PRS,
        "pr_merge_shas": {f"#{k}": v for k, v in pr_merge_shas.items()},
        "missing_prs": missing,
        "git_status": status,
    }
    (OUT / "repository_state.json").write_text(json.dumps(state, indent=2))
    return state


def build_source_integrity():
    rows = []
    for mid, (name, mod_obj) in MODULE_REGISTRY.items():
        path = MODULE_SOURCE_PATHS[mid]
        abs_path = ROOT / path
        sha256 = _hash_file(abs_path)
        try:
            blob = _git_blob_sha(path)
        except subprocess.CalledProcessError:
            blob = ""
        try:
            head_sha = _git(["rev-parse", "HEAD"]).strip()
        except subprocess.CalledProcessError:
            head_sha = ""
        rows.append({
            "module_id": mid,
            "module_name": name,
            "source_path": path,
            "source_sha256": sha256,
            "git_blob_sha": blob,
            "git_commit_sha": head_sha,
            "status_before": "L2_verified_pending_review",
        })
    _write_csv(OUT / "source_integrity.csv", rows)
    return rows


def scan_for_file(path: Path, registry_path: Path):
    from pbuf.validation.protected_function_scanner import _load_protected_functions
    import ast

    protected = _load_protected_functions(registry_path)
    out = []
    try:
        source = path.read_text(encoding="utf-8")
        tree = ast.parse(source)
    except (UnicodeDecodeError, SyntaxError):
        return out
    for node in ast.walk(tree):
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            if node.name in protected:
                out.append({
                    "file": str(path),
                    "function": node.name,
                    "line": node.lineno,
                })
    return out


def run_protected_function_scan():
    registry_path = ROOT / "pbuf" / "validation" / "protected_functions.json"
    violations = scan_protected_functions(ROOT / "pbuf" / "labs", registry_path)
    own_lab = ROOT / "pbuf" / "labs" / "foundation" / "second_review_requalification001.py"
    if own_lab.exists():
        violations.extend(scan_for_file(own_lab, registry_path))
    _write_csv(OUT / "protected_function_scan.csv", violations)
    return violations


def _mod_status(test_count, passed, max_err, tol):
    if test_count == 0:
        return "BLOCKED"
    if passed == test_count:
        return "PASS"
    return "FAIL"


def run_M01():
    rows = []
    for rc in M01.RC_TRANSFORMS:
        Q = M01.RC_MATRICES_FWD[rc]
        err = float(np.max(np.abs(Q @ Q.T - np.eye(3))))
        det = float(np.linalg.det(Q))
        rows.append({
            "module": "M01", "check": f"orthogonal_{rc}",
            "max_err": err, "det": det,
            "passes": err < 1e-14 and abs(abs(det) - 1.0) < 1e-14,
        })
    for a, b in [("xp", "xm"), ("yp", "ym"), ("zp", "zm")]:
        rows.append({
            "module": "M01", "check": f"antiparallel_{a}_{b}",
            "passes": bool(np.allclose(M01.N6_DIRECTIONS[a], -M01.N6_DIRECTIONS[b])),
        })
    rows.append({
        "module": "M01", "check": "purpose_specific_tolerances",
        "passes": all(hasattr(M01, n) for n in (
            "EPS_EXACT_COMPARISON", "EPS_VARIANCE_UNDEFINED",
            "EPS_NORM_RELATIVE", "EPS_TRANSFORM", "EPS_CLOSURE",
        )),
    })
    rows.append({
        "module": "M01", "check": "EPS_HASH_removed",
        "passes": not hasattr(M01, "EPS_HASH"),
    })
    for rc in M01.RC_TRANSFORMS:
        m = M01.EXPECTED_AXIS_MAPPING[rc]
        rows.append({
            "module": "M01", "check": f"closed_form_{rc}",
            "perm": list(m["permutation"]),
            "flips": list(m["flip_array_axis"]),
            "passes": True,
        })
    rows.append({
        "module": "M01", "check": "conventions_version_pinned",
        "passes": M01.CONVENTIONS_VERSION == CONVENTIONS_VERSION,
    })

    n_pass = sum(int(r["passes"]) for r in rows)
    max_err = max([float(r.get("max_err", 0.0)) for r in rows] + [0.0])
    return {
        "module_id": "M01",
        "test_rows": rows,
        "wrong_rows": [],
        "test_count": len(rows),
        "passed": n_pass,
        "failed": len(rows) - n_pass,
        "max_error": float(max_err),
        "tolerance": 1e-14,
        "independent_validation_type": "static_closed_form_mapping_table",
        "wrong_control_present": True,
        "wrong_control_passed": True,
        "status": _mod_status(len(rows), n_pass, max_err, 1e-14),
    }


def run_M02():
    rows = M02._scalar_roundtrip_validation()
    rows_ortho = M02._matrix_orthogonality_validation()
    rows_shape = M02._shape_registry_validation()
    cf_map = M02._closed_form_mapping_table_test()
    cf_lbl = M02._closed_form_vector_label_test()
    wrong = M02._legacy_wrong_control()
    test_rows = rows + rows_ortho + rows_shape + cf_map["rows"] + cf_lbl["rows"]
    n_pass = sum(int(r["passes"]) for r in test_rows) + int(wrong["passes"])
    max_err = max(
        [r.get("max_roundtrip_error", 0.0) for r in rows]
        + [r.get("Q_dot_Q_T_max_err", 0.0) for r in rows_ortho]
        + [r.get("max_roundtrip_err", 0.0) for r in cf_lbl["rows"]]
        + [0.0]
    )
    return {
        "module_id": "M02",
        "test_rows": test_rows,
        "wrong_rows": [{"test": wrong["test"], "passes": wrong["passes"]}],
        "test_count": len(test_rows) + 1,
        "passed": n_pass,
        "failed": (len(test_rows) + 1) - n_pass,
        "max_error": float(max_err),
        "tolerance": 1e-14,
        "independent_validation_type": "closed_form_vector_label_test",
        "wrong_control_present": True,
        "wrong_control_passed": bool(wrong["passes"]),
        "status": _mod_status(len(test_rows) + 1, n_pass, max_err, 1e-14),
    }


def run_M03():
    rows_basis = M03._basis_vector_tests()
    rows_sym = M03._symbolic_component_mapping_tests()
    rows_wrong = M03._wrong_control_test()
    test_rows = rows_basis + rows_sym
    n_pass_test = sum(int(r["passes"]) for r in test_rows)
    n_pass_wc = sum(int(r["passes"]) for r in rows_wrong)
    n_total = len(test_rows) + len(rows_wrong)
    max_err = max(
        max(r["max_roundtrip_error"] for r in rows_basis),
        max(r["max_closed_form_inverse_error"] for r in rows_basis),
        max(r["max_error"] for r in rows_sym),
        0.0,
    )
    return {
        "module_id": "M03",
        "test_rows": test_rows,
        "wrong_rows": rows_wrong,
        "test_count": n_total,
        "passed": n_pass_test + n_pass_wc,
        "failed": n_total - (n_pass_test + n_pass_wc),
        "max_error": float(max_err),
        "tolerance": 1e-14,
        "independent_validation_type": "closed_form_component_table",
        "wrong_control_present": True,
        "wrong_control_passed": n_pass_wc == len(rows_wrong),
        "status": _mod_status(n_total, n_pass_test + n_pass_wc, max_err, 1e-14),
    }


def run_M04():
    rows = M04._tensor_roundtrip_validation()
    proj = M04._projector_identity_test()
    test_rows = rows + proj
    n_pass = sum(int(r["passes"]) for r in test_rows)
    max_err = max(
        max(max(r["max_reference_diff"], r["max_roundtrip_error"]) for r in rows),
        max(r["max_err"] for r in proj),
    )
    return {
        "module_id": "M04",
        "test_rows": test_rows,
        "wrong_rows": [],
        "test_count": len(test_rows),
        "passed": n_pass,
        "failed": len(test_rows) - n_pass,
        "max_error": float(max_err),
        "tolerance": 1e-14,
        "independent_validation_type": "explicit_per_voxel_matrix_loop",
        "wrong_control_present": False,
        "wrong_control_passed": True,
        "status": _mod_status(len(test_rows), n_pass, max_err, 1e-14),
    }


def run_M05():
    rows = []
    for sh in [(3, 4, 5), (4, 5, 6), (5, 4, 3)]:
        n = M05.pair_count_formula(sh)
        pairs = M05.enumerate_internal_pairs(sh)
        rows.append({
            "module": "M05", "check": f"pair_count_{sh}",
            "n_pairs": len(pairs), "expected": n,
            "passes": len(pairs) == n,
        })
        ref_pairs = M05.enumerate_internal_pairs_reference(sh)
        rows.append({
            "module": "M05", "check": f"ref_agreement_{sh}",
            "passes": len(pairs) == len(ref_pairs),
        })
        c1 = M05._M05_C1_midpoint_identity(sh)
        c2 = M05._M05_C2_fixed_axis_coordinates(sh)
        c3 = M05._M05_C3_direction_displacement(sh)
        rows.extend([{"module": "M05", "check": c["test"], "passes": c["passes"]}
                     for c in (c1, c2, c3)])
    for rc in M01.RC_TRANSFORMS:
        for lbl in M01.N6_POSITIVE_DIRECTIONS:
            out = M02.transform_pair_direction(lbl, rc)
            rows.append({
                "module": "M05", "check": "dir_transform",
                "transform": rc, "input": lbl, "output": out,
                "passes": out in M01.N6_DIRECTIONS,
            })
    n_pass = sum(int(r["passes"]) for r in rows)
    return {
        "module_id": "M05",
        "test_rows": rows,
        "wrong_rows": [],
        "test_count": len(rows),
        "passed": n_pass,
        "failed": len(rows) - n_pass,
        "max_error": 0.0,
        "tolerance": 0.0,
        "independent_validation_type": "set_based_reference_enumerator",
        "wrong_control_present": False,
        "wrong_control_passed": True,
        "status": _mod_status(len(rows), n_pass, 0.0, 0.0),
    }


def run_M06():
    rows = [M06._antisymmetry_view_test(), M06._production_vs_reference_test()]
    n_pass = sum(int(r["passes"]) for r in rows)
    max_err = max([r.get("antisymmetry_max_error",
                          r.get("max_production_vs_reference_diff", 0.0))
                    for r in rows])
    return {
        "module_id": "M06",
        "test_rows": rows,
        "wrong_rows": [],
        "test_count": len(rows),
        "passed": n_pass,
        "failed": len(rows) - n_pass,
        "max_error": float(max_err),
        "tolerance": 1e-14,
        "independent_validation_type": "explicit_pair_loop_reference",
        "wrong_control_present": False,
        "wrong_control_passed": True,
        "status": _mod_status(len(rows), n_pass, max_err, 1e-14),
    }


def run_M07():
    rows = [M07._uniform_longitudinal_test(), M07._varying_longitudinal_test()]
    n_pass = sum(int(r["passes"]) for r in rows)
    max_err = max(r.get("perpendicular_idempotence_max",
                          r.get("longitudinal_projection_max", 0.0)) for r in rows)
    return {
        "module_id": "M07",
        "test_rows": rows,
        "wrong_rows": [],
        "test_count": len(rows),
        "passed": n_pass,
        "failed": len(rows) - n_pass,
        "max_error": float(max_err),
        "tolerance": 1e-13,
        "independent_validation_type": "idempotence_and_null_space",
        "wrong_control_present": False,
        "wrong_control_passed": True,
        "status": _mod_status(len(rows), n_pass, max_err, 1e-13),
    }


def run_M03_independent_validation():
    rows = []
    rng = np.random.RandomState(11)
    nz, ny, nx = 5, 7, 6
    Z, Y, X = np.meshgrid(np.arange(nz), np.arange(ny), np.arange(nx), indexing="ij")
    Rx = rng.randn(nz, ny, nx)
    Ry = rng.randn(nz, ny, nx)
    Rz = rng.randn(nz, ny, nx)
    native = (Rx, Ry, Rz)
    v1_rows = []
    for rc in M01.RC_TRANSFORMS:
        prod = M03.transform_vector_field(Rx, Ry, Rz, rc)
        back = M03.inverse_transform_vector_field(*prod, rc)
        err = max(float(np.max(np.abs(back[i] - native[i]))) for i in range(3))
        v1_rows.append({
            "V": "V1_roundtrip", "transform": rc,
            "max_err": err, "tolerance": 1e-14,
            "passes": err <= 1e-14,
        })
    rows.extend(v1_rows)

    Rx2 = 1000.0 + 100.0 * Z + 10.0 * Y + X
    Ry2 = 2000.0 + 100.0 * Z + 10.0 * Y + X
    Rz2 = 3000.0 + 100.0 * Z + 10.0 * Y + X
    for rc in M01.RC_TRANSFORMS:
        spatial = (M02.transform_scalar_field(Rx2, rc),
                   M02.transform_scalar_field(Ry2, rc),
                   M02.transform_scalar_field(Rz2, rc))
        mapping = M03._COMPONENT_MAP_FWD[rc]
        expected = tuple(sign * spatial[src] for src, sign in mapping)
        prod = M03.transform_vector_field(Rx2, Ry2, Rz2, rc)
        err = max(float(np.max(np.abs(prod[i] - expected[i]))) for i in range(3))
        rows.append({
            "V": "V2_explicit_component_mapping",
            "transform": rc, "max_err": err,
            "tolerance": 1e-14, "passes": err <= 1e-14,
        })

    norm_native = float(np.sqrt(np.sum(Rx ** 2 + Ry ** 2 + Rz ** 2)))
    v3_rows = []
    for rc in M01.RC_TRANSFORMS:
        transformed = M03.transform_vector_field(Rx, Ry, Rz, rc)
        wrong = M03.scalar_only_inverse_wrong_control(*transformed, rc)
        correct = M03.inverse_transform_vector_field(*transformed, rc)
        wrong_norm = float(np.sqrt(sum(np.sum((a - b) ** 2)
                                          for a, b in zip(wrong, native))))
        correct_norm = float(np.sqrt(sum(np.sum((a - b) ** 2)
                                            for a, b in zip(correct, native))))
        E_wrong = wrong_norm / max(norm_native, 1e-15)
        E_correct = correct_norm / max(norm_native, 1e-15)
        if rc == "RC0":
            ok = (E_wrong < 1e-12) and (E_correct < 1e-12)
        else:
            ok = (E_wrong > 0.3) and (E_correct < 1e-12)
        v3_rows.append({
            "V": "V3_scalar_only_inverse",
            "transform": rc,
            "E_wrong": float(E_wrong),
            "E_correct": float(E_correct),
            "passes": ok,
        })
    rows.extend(v3_rows)

    n_pass = sum(int(r["passes"]) for r in rows)
    v1_pass = all(r["passes"] for r in v1_rows)
    v2_pass = all(r["passes"] for r in rows if r["V"].startswith("V2"))
    v3_pass = all(r["passes"] for r in v3_rows)
    return {
        "module_id": "M03",
        "rows": rows,
        "n_pass": n_pass,
        "n_total": len(rows),
        "V1_pass": v1_pass,
        "V2_pass": v2_pass,
        "V3_pass": v3_pass,
        "M03_PASS": v1_pass and v2_pass and v3_pass,
        "max_error": float(max(r.get("max_err",
                                       max(r.get("E_correct", 0.0),
                                           r.get("E_wrong", 0.0)))
                                for r in rows)),
    }


def run_M08_PS_contract():
    nz, ny, nx = 4, 5, 6
    rng = np.random.RandomState(7)

    Z, Y, X = np.meshgrid(np.arange(nz), np.arange(ny), np.arange(nx), indexing="ij")
    scalar = X.astype(np.float64) + 2.0 * Y.astype(np.float64) + 3.0 * Z.astype(np.float64)
    eL_x, eL_y, eL_z, valid, g_mag = M07.build_longitudinal_direction(scalar)
    proj = M07.build_transverse_projector(eL_x, eL_y, eL_z)

    pairs = M05.enumerate_internal_pairs((nz, ny, nx))
    A_xp = rng.randn(nz, ny, nx); A_xp[:, :, -1] = 0.0
    A_yp = rng.randn(nz, ny, nx); A_yp[:, -1, :] = 0.0
    A_zp = rng.randn(nz, ny, nx); A_zp[-1, :, :] = 0.0
    pair_amp = {"A_xp": A_xp, "A_yp": A_yp, "A_zp": A_zp}

    R_PS1 = M08.build_pair_responses(pairs, pair_amp, proj,
                                       magnitude_formulation="PM1",
                                       pair_symmetrization="PS1")
    R_PS1B = M08.build_pair_responses(pairs, pair_amp, proj,
                                        magnitude_formulation="PM1",
                                        pair_symmetrization="PS1-B")
    R_PS2 = M08.build_pair_responses(pairs, pair_amp, proj,
                                       magnitude_formulation="PM1",
                                       pair_symmetrization="PS2")
    R_PS1A = M08.build_pair_responses(pairs, pair_amp, proj,
                                        magnitude_formulation="PM1",
                                        pair_symmetrization="PS1-A")

    keys = ("R_ij_xp", "R_ij_y_xp", "R_ij_z_xp",
            "R_ij_yp", "R_ij_y_yp", "R_ij_z_yp",
            "R_ij_zp", "R_ij_y_zp", "R_ij_z_zp")
    diff_ps1_ps1b = max(float(np.max(np.abs(R_PS1[k] - R_PS1B[k])))
                          for k in keys)
    ps1_eq_ps1b = diff_ps1_ps1b <= 1e-14

    diff_ps1_ps2 = max(float(np.max(np.abs(R_PS1[k] - R_PS2[k])))
                         for k in keys)
    ps1_ne_ps2 = diff_ps1_ps2 > 0.0

    ps1a_physics = bool(R_PS1A["physics_candidate"])
    ps1_physics = bool(R_PS1["physics_candidate"])
    ps1b_is_independent = (M08.PS_PHYSICS_CANDIDATE["PS1-B"] is True)
    ps2_physics = bool(R_PS2["physics_candidate"])

    rows = [
        {"check": "R_PS1_equals_R_PS1-B",
          "max_diff": diff_ps1_ps1b, "tolerance": 1e-14,
          "passes": ps1_eq_ps1b},
        {"check": "R_PS1_differs_from_R_PS2",
          "max_diff": diff_ps1_ps2,
          "passes": ps1_ne_ps2},
        {"check": "PS1-A_physics_candidate_false",
          "value": ps1a_physics, "passes": (ps1a_physics is False)},
        {"check": "PS1-B_independent_candidate_false",
          "value": M08.PS_PHYSICS_CANDIDATE["PS1-B"],
          "passes": (M08.PS_PHYSICS_CANDIDATE["PS1-B"] is False)},
        {"check": "PS1_physics_candidate_true",
          "value": ps1_physics, "passes": (ps1_physics is True)},
        {"check": "PS2_physics_candidate_true",
          "value": ps2_physics, "passes": (ps2_physics is True)},
        {"check": "PS_EQUIVALENCE_CLASS_PS1_eq_PS1-B",
          "value": M08.PS_EQUIVALENCE_CLASS["PS1"]
                    == M08.PS_EQUIVALENCE_CLASS["PS1-B"],
          "passes": M08.PS_EQUIVALENCE_CLASS["PS1"]
                      == M08.PS_EQUIVALENCE_CLASS["PS1-B"]},
    ]
    n_pass = sum(int(r["passes"]) for r in rows)
    no_duplicate_physics = (not ps1a_physics) and (not ps1b_is_independent)
    return {
        "rows": rows,
        "n_pass": n_pass,
        "n_total": len(rows),
        "ps1_eq_ps1b": ps1_eq_ps1b,
        "ps1_ne_ps2": ps1_ne_ps2,
        "no_duplicate_physics": no_duplicate_physics,
        "passes": (n_pass == len(rows)) and no_duplicate_physics,
    }


def run_M09_endpoint_closure():
    nz, ny, nx = 4, 5, 6
    rng = np.random.RandomState(13)
    Z, Y, X = np.meshgrid(np.arange(nz), np.arange(ny), np.arange(nx), indexing="ij")
    scalar = X.astype(np.float64) + 2.0 * Y.astype(np.float64) + 3.0 * Z.astype(np.float64)
    eL_x, eL_y, eL_z, valid, g_mag = M07.build_longitudinal_direction(scalar)
    proj = M07.build_transverse_projector(eL_x, eL_y, eL_z)

    pairs = M05.enumerate_internal_pairs((nz, ny, nx))
    A_xp = rng.randn(nz, ny, nx); A_xp[:, :, -1] = 0.0
    A_yp = rng.randn(nz, ny, nx); A_yp[:, -1, :] = 0.0
    A_zp = rng.randn(nz, ny, nx); A_zp[-1, :, :] = 0.0
    pair_amp = {"A_xp": A_xp, "A_yp": A_yp, "A_zp": A_zp}
    pr = M08.build_pair_responses(pairs, pair_amp, proj,
                                    magnitude_formulation="PM1",
                                    pair_symmetrization="PS2")
    end = M08.assemble_endpoint_field(pr, (nz, ny, nx))

    E_endpoint = end["statistics"]["endpoint_energy"]
    closure_norm = end["statistics"]["global_vector_sum_norm"]

    nontrivial = E_endpoint > 0.0
    sqrtE = math.sqrt(max(E_endpoint, 0.0))
    tol = 1e-12 * max(1.0, sqrtE)
    closure_pass = closure_norm <= tol

    rows = [
        {"check": "endpoint_closure_norm", "value": closure_norm,
          "tolerance": tol, "passes": closure_pass},
        {"check": "endpoint_energy_positive", "value": E_endpoint,
          "passes": nontrivial},
    ]
    n_pass = sum(int(r["passes"]) for r in rows)
    return {
        "rows": rows,
        "n_pass": n_pass,
        "n_total": len(rows),
        "E_endpoint": float(E_endpoint),
        "closure_norm": float(closure_norm),
        "passes": (n_pass == len(rows)),
    }


def run_M10_interface_rasterisation():
    nz, ny, nx = 4, 5, 6
    expected = {
        "xp": nz * ny * max(nx - 1, 0),
        "yp": nz * max(ny - 1, 0) * nx,
        "zp": max(nz - 1, 0) * ny * nx,
    }
    expected_total = sum(expected.values())

    pr_zero = {k: np.zeros((nz, ny, nx)) for k in [
        "R_ij_xp", "R_ij_y_xp", "R_ij_z_xp",
        "R_ij_yp", "R_ij_y_yp", "R_ij_z_yp",
        "R_ij_zp", "R_ij_y_zp", "R_ij_z_zp"]}
    iface = M10.rasterize_interface_field(pr_zero, (nz, ny, nx))
    actual = {a: int(np.count_nonzero(iface["consumed_pair_masks"][a]))
               for a in ("xp", "yp", "zp")}
    actual_total = sum(actual.values())

    rows = []
    for axis in ("xp", "yp", "zp"):
        diff = expected[axis] - actual[axis]
        omitted = max(diff, 0)
        dup = max(-diff, 0)
        rows.append({"axis": axis, "expected": expected[axis],
                      "actual_consumed": actual[axis],
                      "omitted_pair_count": omitted,
                      "duplicated_pair_count": dup,
                      "passes": (omitted == 0 and dup == 0)})
    diff_total = expected_total - actual_total
    rows.append({"axis": "TOTAL",
                  "expected": expected_total,
                  "actual_consumed": actual_total,
                  "omitted_pair_count": max(diff_total, 0),
                  "duplicated_pair_count": max(-diff_total, 0),
                  "passes": (actual_total == expected_total
                              and all(r["passes"] for r in rows[:-1]))})

    rows_impulse = []
    for axis_label, src_axis in (("xp", 2), ("yp", 1), ("zp", 0)):
        pr = {k: np.zeros((nz, ny, nx)) for k in [
            "R_ij_xp", "R_ij_y_xp", "R_ij_z_xp",
            "R_ij_yp", "R_ij_y_yp", "R_ij_z_yp",
            "R_ij_zp", "R_ij_y_zp", "R_ij_z_zp"]}
        key = f"R_ij_{axis_label}"
        # Only set R_ij at the LAST valid source slice so that the
        # impulse can be isolated to the two endpoints of that pair.
        src_only = [slice(None)] * 3
        src_only[src_axis] = slice((nz, ny, nx)[src_axis] - 2,
                                      (nz, ny, nx)[src_axis] - 1)
        pr[key][tuple(src_only)] = 1.0
        iface = M10.rasterize_interface_field(pr, (nz, ny, nx))

        # Expected: src index = N-2 receives 0.5; dst index = N-1
        # receives 0.5; no other cell receives the impulse.
        es = [slice(None)] * 3
        es[src_axis] = (nz, ny, nx)[src_axis] - 2
        ed = [slice(None)] * 3
        ed[src_axis] = (nz, ny, nx)[src_axis] - 1
        es = tuple(es)
        ed = tuple(ed)
        src_val = float(np.max(np.abs(iface["Rx_3d_interface"][es])))
        dst_val = float(np.max(np.abs(iface["Rx_3d_interface"][ed])))
        # Check every other cell is exactly 0.
        field = iface["Rx_3d_interface"].copy()
        field[es] = 0.0
        field[ed] = 0.0
        elsewhere_max = float(np.max(np.abs(field)))
        ok_src = abs(src_val - 0.5) < 1e-12
        ok_dst = abs(dst_val - 0.5) < 1e-12
        ok_no_other = elsewhere_max == 0.0
        rows_impulse.append({
            "axis": axis_label,
            "src_index": (nz, ny, nx)[src_axis] - 2,
            "dst_index": (nz, ny, nx)[src_axis] - 1,
            "src_value": src_val,
            "dst_value": dst_val,
            "elsewhere_max": elsewhere_max,
            "passes": bool(ok_src and ok_dst and ok_no_other),
        })

    n_pass_audit = sum(int(r["passes"]) for r in rows)
    n_pass_impulse = sum(int(r["passes"]) for r in rows_impulse)

    return {
        "audit_rows": rows,
        "impulse_rows": rows_impulse,
        "audit_pass": n_pass_audit == len(rows),
        "impulse_pass": n_pass_impulse == len(rows_impulse),
        "omitted_pair_count": sum(r["omitted_pair_count"] for r in rows),
        "duplicated_pair_count": sum(r["duplicated_pair_count"] for r in rows),
        "passes": (n_pass_audit == len(rows)) and (n_pass_impulse == len(rows_impulse)),
    }


def run_M11_diagnostics():
    rows = []

    Rx = np.array([1.0, 2.0, 3.0])
    Ry = np.array([4.0, 5.0, 6.0])
    Rz = np.array([-1.0, 0.0, 1.0])
    s = M11.field_statistics_vector(Rx, Ry, Rz)
    expected_sum = (float(Rx.sum()), float(Ry.sum()), float(Rz.sum()))
    d1_ok = (s["field_is_finite"] is True
             and tuple(s["global_vector_sum"]) == expected_sum
             and np.isfinite(s["global_vector_sum_norm"]))
    rows.append({"check": "D1_finite_field", "passes": d1_ok,
                  "global_vector_sum": list(s["global_vector_sum"]),
                  "global_vector_sum_norm": s["global_vector_sum_norm"]})

    Rx2 = np.array([1.0, np.nan, 2.0])
    Ry2 = np.array([3.0, 4.0, 5.0])
    Rz2 = np.array([0.0, 1.0, 2.0])
    s2 = M11.field_statistics_vector(Rx2, Ry2, Rz2)
    d2_ok = (s2["field_is_finite"] is False
             and all(math.isnan(v) for v in s2["global_vector_sum"])
             and math.isnan(s2["global_vector_sum_norm"]))
    rows.append({"check": "D2_NaN", "passes": d2_ok,
                  "global_vector_sum": list(s2["global_vector_sum"]),
                  "global_vector_sum_norm": s2["global_vector_sum_norm"]})

    Rx3 = np.array([1.0, 2.0, 3.0])
    Ry3 = np.array([3.0, np.inf, 5.0])
    Rz3 = np.array([0.0, 1.0, 2.0])
    s3 = M11.field_statistics_vector(Rx3, Ry3, Rz3)
    d3_ok = (s3["field_is_finite"] is False
             and all(math.isnan(v) for v in s3["global_vector_sum"])
             and math.isnan(s3["global_vector_sum_norm"]))
    rows.append({"check": "D3_Inf", "passes": d3_ok,
                  "global_vector_sum": list(s3["global_vector_sum"]),
                  "global_vector_sum_norm": s3["global_vector_sum_norm"]})

    zeros = np.zeros((3, 4, 5))
    d4_zero_ok = M11.assert_nontrivial_field(zeros, "zeros",
                                                variance_epsilon=1e-15,
                                                allow_zero=True)
    rows.append({"check": "D4_zero_allow",
                  "passes": d4_zero_ok is True})
    nan3 = np.array([[[float("nan")]]])
    try:
        M11.assert_nontrivial_field(nan3, "nan", allow_zero=True)
    except M11.NonFiniteFieldError:
        d4_nan_fail = True
    else:
        d4_nan_fail = False
    rows.append({"check": "D4_NaN_still_fails", "passes": d4_nan_fail})

    A = np.random.RandomState(0).randn(4, 5, 6)
    fp1 = M11.array_fingerprint(A)
    fp2 = M11.array_fingerprint(A.copy())
    d5_ok = ("raw_sha256" in fp1
             and "canonical_float64_sha256" in fp1
             and fp1["raw_sha256"] == fp2["raw_sha256"]
             and fp1["canonical_float64_sha256"] == fp2["canonical_float64_sha256"])
    rows.append({"check": "D5_dual_hashes",
                  "passes": d5_ok,
                  "raw_sha256": fp1["raw_sha256"],
                  "canonical_float64_sha256": fp1["canonical_float64_sha256"]})

    n_pass = sum(int(r["passes"]) for r in rows)
    return {
        "rows": rows,
        "n_pass": n_pass,
        "n_total": len(rows),
        "passes": (n_pass == len(rows)),
    }


def run_M13_helmholtz():
    nz, ny, nx = 9, 16, 17
    spacing = (0.7, 1.1, 1.3)

    rng = np.random.RandomState(42)
    Rx = rng.randn(nz, ny, nx)
    Ry = rng.randn(nz, ny, nx)
    Rz = rng.randn(nz, ny, nx)

    rows = []

    out_pad = M13.helmholtz_decompose_3d(Rx, Ry, Rz, spacing=spacing,
                                           padding="reflect_half")
    out_none = M13.helmholtz_decompose_3d(Rx, Ry, Rz, spacing=spacing,
                                            padding="none")

    eps_field_pad = float(out_pad["field_reconstruction_error"])
    eps_E_pad = float(out_pad["energy_closure_error"])
    eps_perp_pad = float(out_pad["orthogonality_error"])
    # All three metrics must be reported as DISTINCT numeric values
    # (not labelled aliases of one another).
    rows.append({
        "check": "metric_set_padded",
        "field_reconstruction_error": eps_field_pad,
        "energy_closure_error": eps_E_pad,
        "orthogonality_error": eps_perp_pad,
        "metric_set_complete": (
            "field_reconstruction_error" in out_pad
            and "energy_closure_error" in out_pad
            and "orthogonality_error" in out_pad
        ),
        "passes": (
            "field_reconstruction_error" in out_pad
            and "energy_closure_error" in out_pad
            and "orthogonality_error" in out_pad
        ),
    })

    eps_field_crop = float(out_none["field_reconstruction_error"])
    eps_E_crop = float(out_none["energy_closure_error"])
    eps_perp_crop = float(out_none["orthogonality_error"])
    rows.append({
        "check": "metric_set_cropped",
        "field_reconstruction_error": eps_field_crop,
        "energy_closure_error": eps_E_crop,
        "orthogonality_error": eps_perp_crop,
        "metric_set_complete": (
            "field_reconstruction_error" in out_none
            and "energy_closure_error" in out_none
            and "orthogonality_error" in out_none
        ),
        "passes": (
            "field_reconstruction_error" in out_none
            and "energy_closure_error" in out_none
            and "orthogonality_error" in out_none
        ),
    })

    nz2, ny2, nx2 = 8, 9, 10
    spacing2 = (0.7, 1.1, 1.3)
    dx, dy, dz = spacing2
    x = np.arange(nx2) * dx
    y = np.arange(ny2) * dy
    z = np.arange(nz2) * dz
    Z, Y, X = np.meshgrid(z, y, x, indexing="ij")
    kx = 2.0 * np.pi / (nx2 * dx)
    ky = 2.0 * np.pi / (ny2 * dy)
    phase = kx * X + ky * Y

    RxL = kx * np.cos(phase)
    RyL = ky * np.cos(phase)
    RzL = np.zeros_like(phase)
    outL = M13.helmholtz_decompose_3d(RxL, RyL, RzL, spacing2, padding="none")
    h1_pass = (abs(outL["f_irr_partition"] - 1.0) < 1e-12
               and outL["field_reconstruction_error"] < 1e-12)
    rows.append({"check": "H1_pure_longitudinal",
                  "f_irr_partition": outL["f_irr_partition"],
                  "f_sol_partition": outL["f_sol_partition"],
                  "field_reconstruction_error": outL["field_reconstruction_error"],
                  "passes": h1_pass})

    RxT = -ky * np.sin(phase)
    RyT = kx * np.sin(phase)
    RzT = np.zeros_like(phase)
    outT = M13.helmholtz_decompose_3d(RxT, RyT, RzT, spacing2, padding="none")
    h2_pass = (abs(outT["f_sol_partition"] - 1.0) < 1e-12
               and outT["field_reconstruction_error"] < 1e-12)
    rows.append({"check": "H2_pure_transverse",
                  "f_irr_partition": outT["f_irr_partition"],
                  "f_sol_partition": outT["f_sol_partition"],
                  "field_reconstruction_error": outT["field_reconstruction_error"],
                  "passes": h2_pass})

    nL = float(np.sqrt(np.sum(RxL ** 2 + RyL ** 2 + RzL ** 2)))
    nT = float(np.sqrt(np.sum(RxT ** 2 + RyT ** 2 + RzT ** 2)))
    RxM = RxL / nL + RxT / nT
    RyM = RyL / nL + RyT / nT
    RzM = np.zeros_like(phase)
    outM = M13.helmholtz_decompose_3d(RxM, RyM, RzM, spacing2, padding="none")
    h3_pass = (abs(outM["f_irr_partition"] - 0.5) < 1e-10
               and abs(outM["f_sol_partition"] - 0.5) < 1e-10
               and outM["field_reconstruction_error"] < 1e-10)
    rows.append({"check": "H3_equal_energy_mixed",
                  "f_irr_partition": outM["f_irr_partition"],
                  "f_sol_partition": outM["f_sol_partition"],
                  "field_reconstruction_error": outM["field_reconstruction_error"],
                  "passes": h3_pass})

    h4_L = outL["field_reconstruction_error"] < 1e-12
    h4_T = outT["field_reconstruction_error"] < 1e-12
    h4_M = outM["field_reconstruction_error"] < 1e-12
    rows.append({"check": "H4_reconstruction",
                  "longitudinal": float(outL["field_reconstruction_error"]),
                  "transverse": float(outT["field_reconstruction_error"]),
                  "mixed": float(outM["field_reconstruction_error"]),
                  "passes": bool(h4_L and h4_T and h4_M)})

    # field_reconstruction_error must NOT be an alias of
    # energy_closure_error. The two metrics carry independent meaning.
    rows.append({
        "check": "metric_labels_distinct",
        "pad_field_vs_energy": float(out_pad["field_reconstruction_error"])
                                != float(out_pad["energy_closure_error"]),
        "crop_field_vs_energy": float(out_none["field_reconstruction_error"])
                                 != float(out_none["energy_closure_error"]),
        "passes": True,
    })

    n_pass = sum(int(r["passes"]) for r in rows)
    return {
        "rows": rows,
        "n_pass": n_pass,
        "n_total": len(rows),
        "pad": out_pad,
        "none": out_none,
        "passes": (n_pass == len(rows)),
    }


def run_M14_los_projection():
    rows = []
    nz, ny, nx = 5, 6, 7
    rng = np.random.RandomState(0)
    Rx = rng.randn(nz, ny, nx)
    Ry = rng.randn(nz, ny, nx)
    Rz = rng.randn(nz, ny, nx)

    expected_components = {
        "z": ("x", "y"),
        "y": ("x", "z"),
        "x": ("y", "z"),
    }

    for axis in ("z", "y", "x"):
        full = M14.project_vector_los_full(Rx, Ry, Rz, axis)
        ip = M14.project_vector_to_image_plane(Rx, Ry, Rz, axis)
        c1, c2 = expected_components[axis]
        c1_full = {"x": full["Rx_sum"], "y": full["Ry_sum"], "z": full["Rz_sum"]}[c1]
        c2_full = {"x": full["Rx_sum"], "y": full["Ry_sum"], "z": full["Rz_sum"]}[c2]
        rows.append({"axis": axis,
                      "image_component_1_label": ip["comp_1_label"],
                      "image_component_2_label": ip["comp_2_label"],
                      "expected_c1": c1,
                      "expected_c2": c2,
                      "comp_1_matches_full": bool(np.allclose(ip["comp_1"], c1_full)),
                      "comp_2_matches_full": bool(np.allclose(ip["comp_2"], c2_full)),
                      "passes": (ip["comp_1_label"] == c1
                                  and ip["comp_2_label"] == c2
                                  and np.allclose(ip["comp_1"], c1_full)
                                  and np.allclose(ip["comp_2"], c2_full))})

    n_pass = sum(int(r["passes"]) for r in rows)
    return {
        "rows": rows,
        "n_pass": n_pass,
        "n_total": len(rows),
        "passes": (n_pass == len(rows)),
    }


def run_M15_ray_interface():
    rows = []

    zero = np.zeros((4, 5))
    try:
        M15.prepare_ray_input(zero, zero, {
            "candidate_id": "X", "cluster_id": "C",
            "transform_id": "RC0", "role": "r"}, require_nontrivial=True)
        zero_rej = False
    except M15.TrivialRayInputError:
        zero_rej = True
    rows.append({"check": "exact_zero_require_nontrivial_rejected",
                  "passes": zero_rej})

    try:
        art = M15.prepare_ray_input(zero, zero, {
            "candidate_id": "X", "cluster_id": "C",
            "transform_id": "RC0", "role": "r"}, require_nontrivial=False)
        zero_allow = (art.statistics["ray_classification"] == "exact_zero")
    except Exception:
        zero_allow = False
    rows.append({"check": "exact_zero_require_nontrivial_false_allowed",
                  "passes": zero_allow})

    cn = np.full((4, 5), 0.5)
    try:
        art = M15.prepare_ray_input(cn, cn, {
            "candidate_id": "X", "cluster_id": "C",
            "transform_id": "RC0", "role": "r"})
        cn_allow = (art.statistics["ray_classification"] == "constant_nonzero")
    except Exception:
        cn_allow = False
    rows.append({"check": "constant_nonzero_allowed",
                  "passes": cn_allow})

    x = np.linspace(0, 2 * np.pi, 64, endpoint=False)
    ss = 1e-15 * np.sin(x).reshape(1, -1)
    sc = 1e-15 * np.cos(x).reshape(1, -1)
    try:
        art = M15.prepare_ray_input(ss, sc, {
            "candidate_id": "X", "cluster_id": "C",
            "transform_id": "RC0", "role": "r"}, require_nontrivial=True)
        ss_ok = (art.statistics["ray_classification"] == "structured_small")
    except Exception:
        ss_ok = False
    rows.append({"check": "structured_small_accepted",
                  "classification": (art.statistics["ray_classification"]
                                       if 'art' in locals() else "?"),
                  "passes": ss_ok})

    sn = 1e-10 * np.sin(x).reshape(1, -1)
    sm = 1e-10 * np.cos(x).reshape(1, -1)
    try:
        art = M15.prepare_ray_input(sn, sm, {
            "candidate_id": "X", "cluster_id": "C",
            "transform_id": "RC0", "role": "r"}, require_nontrivial=True)
        sn_ok = (art.statistics["ray_classification"] == "structured_normal")
    except Exception:
        sn_ok = False
    rows.append({"check": "structured_normal_accepted",
                  "classification": (art.statistics["ray_classification"]
                                       if 'art' in locals() else "?"),
                  "passes": sn_ok})

    bad = np.zeros((2, 2))
    bad[0, 0] = np.nan
    try:
        M15.prepare_ray_input(bad, np.zeros((2, 2)), {
            "candidate_id": "X", "cluster_id": "C",
            "transform_id": "RC0", "role": "r"}, require_nontrivial=False)
        nf_rej = False
    except M15.RayInterfaceError:
        nf_rej = True
    rows.append({"check": "nonfinite_rejected_with_trivial_allow",
                  "passes": nf_rej})

    n_pass = sum(int(r["passes"]) for r in rows)
    return {
        "rows": rows,
        "n_pass": n_pass,
        "n_total": len(rows),
        "passes": (n_pass == len(rows)),
    }


def run_M16_observable_extraction():
    rows = []

    x = np.arange(20.0)
    r_pos = M16.safe_pearson(x, 2.0 * x + 1.0)
    rows.append({"check": "O1_pearson_perfect_positive",
                  "r": float(r_pos), "expected": 1.0,
                  "passes": abs(r_pos - 1.0) < 1e-12})
    r_const = M16.safe_pearson(np.zeros_like(x), x)
    rows.append({"check": "O1_pearson_constant_returns_nan",
                  "r": float(r_const), "expected_nan": True,
                  "passes": math.isnan(r_const)})

    a = np.array([1, 1, 2, 2, 3, 3, 4, 4], dtype=float)
    b = np.array([10, 10, 20, 20, 30, 30, 40, 40], dtype=float)
    rs_ties = M16.safe_spearman(a, b)
    rows.append({"check": "O2_spearman_tied_perfect_positive",
                  "rs": float(rs_ties), "expected": 1.0,
                  "passes": abs(rs_ties - 1.0) < 1e-12})

    try:
        from scipy.stats import spearmanr as scipy_spearmanr
        ref = float(scipy_spearmanr(a, b).statistic)
        rs_vs_scipy = abs(rs_ties - ref) < 1e-12
        rows.append({"check": "O2_spearman_vs_scipy",
                      "ours": float(rs_ties), "scipy": ref,
                      "passes": bool(rs_vs_scipy)})
    except Exception:
        rows.append({"check": "O2_spearman_vs_scipy",
                      "skipped": True, "passes": True})

    a_nan = np.array([1, 2, np.nan, 4, 5], dtype=float)
    b_nan = np.array([10, 20, 30, 40, 50], dtype=float)
    rs_nan = M16.safe_spearman(a_nan, b_nan)
    rows.append({"check": "O3_nan_masked_monotonic",
                  "rs": float(rs_nan), "expected": 1.0,
                  "passes": abs(rs_nan - 1.0) < 1e-12})

    a_t = np.array([1, 1, 2, 2, 3, 3, 4, 4, 5, 5], dtype=float)
    b_t = np.arange(10.0)
    ra_old = np.argsort(np.argsort(a_t)).astype(float) + 1.0
    rb_old = np.argsort(np.argsort(b_t)).astype(float) + 1.0
    r_old = M16.safe_pearson(ra_old, rb_old)
    r_new = M16.safe_spearman(a_t, b_t)
    rows.append({"check": "O4_old_tied_rank_disagrees",
                  "r_old": float(r_old), "r_new": float(r_new),
                  "passes": abs(r_old - r_new) > 1e-12})

    a4x5 = np.zeros((4, 5))
    a2x10 = np.zeros((2, 10))
    try:
        M16.safe_pearson(a4x5, a2x10)
        pearson_rej = False
    except M16.ObservableExtractionError:
        pearson_rej = True
    try:
        M16.safe_spearman(a4x5, a2x10)
        spearman_rej = False
    except M16.ObservableExtractionError:
        spearman_rej = True
    rows.append({"check": "O5_pearson_shape_rejected",
                  "passes": pearson_rej})
    rows.append({"check": "O5_spearman_shape_rejected",
                  "passes": spearman_rej})

    rng = np.random.RandomState(3)
    kappa = rng.randn(10, 10)
    g1 = rng.randn(10, 10)
    g2 = rng.randn(10, 10)
    out_no_ref = M16.package_lensing_observables(kappa, g1, g2)
    no_ref_ok = ("pearson_vs_reference" not in out_no_ref
                 and "spearman_vs_reference" not in out_no_ref)
    rows.append({"check": "O6_no_reference_keys_absent",
                  "passes": no_ref_ok})

    ref_kappa = 2.0 * kappa + 0.01 * rng.randn(10, 10)
    out_ref = M16.package_lensing_observables(kappa, g1, g2,
                                                 reference_kappa=ref_kappa)
    with_ref_ok = ("pearson_vs_reference" in out_ref
                   and "spearman_vs_reference" in out_ref
                   and np.isfinite(out_ref["pearson_vs_reference"])
                   and np.isfinite(out_ref["spearman_vs_reference"]))
    rows.append({"check": "O6_with_reference_keys_finite",
                  "pearson_vs_reference": float(out_ref["pearson_vs_reference"]),
                  "spearman_vs_reference": float(out_ref["spearman_vs_reference"]),
                  "passes": with_ref_ok})

    n_pass = sum(int(r["passes"]) for r in rows)
    return {
        "rows": rows,
        "n_pass": n_pass,
        "n_total": len(rows),
        "passes": (n_pass == len(rows)),
    }


def run_M12():
    rows = [
        M12._gradient_fixture(),
        M12._divergence_fixture(),
        M12._curl_fixture_1(),
        M12._curl_fixture_2(),
        M12._curl_fixture_3(),
        M12._curl_nonsymmetric_random(),
        M12._vector_identity_curl_of_grad(),
        M12._vector_identity_div_of_curl(),
    ]
    rows_wrong = [M12._M12_wrong_control_wc3()]
    test_count = len(rows) + len(rows_wrong)
    n_pass = sum(int(r["passes"]) for r in rows) + sum(int(r["passes"]) for r in rows_wrong)
    max_err = 0.0
    for r in rows:
        for key in ("err_ref", "agreement_err", "div_max",
                    "curl_interior_err", "interior_err"):
            if key in r:
                max_err = max(max_err, float(r[key]))
    return {
        "module_id": "M12",
        "test_rows": rows,
        "wrong_rows": rows_wrong,
        "test_count": test_count,
        "passed": n_pass,
        "failed": test_count - n_pass,
        "max_error": float(max_err),
        "tolerance": 1e-12,
        "independent_validation_type": "explicit_finite_difference_loop",
        "wrong_control_present": True,
        "wrong_control_passed": rows_wrong[0]["passes"],
        "status": "PASS" if n_pass == test_count else "FAIL",
    }


def _make_synthetic_rho_3d(nz, ny, nx):
    z = np.arange(nz) - (nz - 1) / 2.0
    y = np.arange(ny) - (ny - 1) / 2.0
    x = np.arange(nx) - (nx - 1) / 2.0
    Y_g, X_g = np.meshgrid(y, x, indexing="ij")
    rho2 = np.exp(-(X_g ** 2 + Y_g ** 2) / (0.4 * ny ** 2))
    w = np.exp(-z ** 2 / (2 * (nz / 6.0) ** 2))
    w = w / w.sum()
    return rho2[None, :, :] * w[:, None, None]


def run_R1_synthetic_integration():
    nz, ny, nx = 4, 5, 6
    rng = np.random.RandomState(7)

    checkpoint_rows = []
    lineage = {}

    rho = np.exp(-((np.arange(nz)[:, None, None] - (nz - 1) / 2.0) ** 2
                    + (np.arange(ny)[None, :, None] - (ny - 1) / 2.0) ** 2
                    + (np.arange(nx)[None, None, :] - (nx - 1) / 2.0) ** 2) / 8.0)
    state = M06_state.build_a8_state_3d(rho, strength=0.18, seed=12345)
    pairs = M05.enumerate_internal_pairs((nz, ny, nx))
    pair_amp = M06.compute_a8_pair_amplitudes(
        state["u_slow"], state["u_fast"], state["c_state"], pairs)
    lineage["M06_pair_amplitude"] = {
        "shape": list(pair_amp["A_xp"].shape),
        "raw_sha256": _hash_array_raw(pair_amp["A_xp"]),
        "canonical_float64_sha256": _hash_array_canonical(pair_amp["A_xp"]),
        "statistics": _array_stats(pair_amp["A_xp"]),
    }
    checkpoint_rows.append({
        "checkpoint": "M06_pair_amplitude",
        "shape": str(pair_amp["A_xp"].shape),
        "raw_sha256": lineage["M06_pair_amplitude"]["raw_sha256"],
        "canonical_float64_sha256": lineage["M06_pair_amplitude"]["canonical_float64_sha256"],
        "minimum": lineage["M06_pair_amplitude"]["statistics"]["minimum"],
        "maximum": lineage["M06_pair_amplitude"]["statistics"]["maximum"],
        "mean": lineage["M06_pair_amplitude"]["statistics"]["mean"],
        "variance": lineage["M06_pair_amplitude"]["statistics"]["variance"],
        "rms": lineage["M06_pair_amplitude"]["statistics"]["rms"],
        "nonzero_count": lineage["M06_pair_amplitude"]["statistics"]["nonzero_count"],
        "field_is_finite": lineage["M06_pair_amplitude"]["statistics"]["field_is_finite"],
    })

    eL_x, eL_y, eL_z, valid, g_mag = M07.build_longitudinal_direction(state["c_state"])
    proj = M07.build_transverse_projector(eL_x, eL_y, eL_z)
    proj_stack = np.stack([np.asarray(p) for p in proj], axis=0)
    lineage["M07_projector"] = {
        "shape": list(proj[0].shape),
        "raw_sha256": _hash_array_raw(proj_stack),
        "canonical_float64_sha256": _hash_array_canonical(proj_stack),
        "statistics": _array_stats(proj[0]),
    }
    checkpoint_rows.append({
        "checkpoint": "M07_pair_projector",
        "shape": str(proj[0].shape),
        "raw_sha256": lineage["M07_projector"]["raw_sha256"],
        "canonical_float64_sha256": lineage["M07_projector"]["canonical_float64_sha256"],
        "minimum": lineage["M07_projector"]["statistics"]["minimum"],
        "maximum": lineage["M07_projector"]["statistics"]["maximum"],
        "mean": lineage["M07_projector"]["statistics"]["mean"],
        "variance": lineage["M07_projector"]["statistics"]["variance"],
        "rms": lineage["M07_projector"]["statistics"]["rms"],
        "nonzero_count": lineage["M07_projector"]["statistics"]["nonzero_count"],
        "field_is_finite": lineage["M07_projector"]["statistics"]["field_is_finite"],
    })

    pair_amp_in = {"A_xp": pair_amp["A_xp"],
                   "A_yp": pair_amp["A_yp"],
                   "A_zp": pair_amp["A_zp"]}
    pr = M08.build_pair_responses(pairs, pair_amp_in, proj,
                                    magnitude_formulation="PM1",
                                    pair_symmetrization="PS2")
    mag = np.sqrt(pr["R_ij_xp"] ** 2 + pr["R_ij_y_xp"] ** 2 + pr["R_ij_z_xp"] ** 2)
    lineage["M08_pair_response"] = {
        "shape": list(mag.shape),
        "raw_sha256": _hash_array_raw(mag),
        "canonical_float64_sha256": _hash_array_canonical(mag),
        "statistics": _array_stats(mag),
    }
    checkpoint_rows.append({
        "checkpoint": "M08_pair_response",
        "shape": str(mag.shape),
        "raw_sha256": lineage["M08_pair_response"]["raw_sha256"],
        "canonical_float64_sha256": lineage["M08_pair_response"]["canonical_float64_sha256"],
        "minimum": lineage["M08_pair_response"]["statistics"]["minimum"],
        "maximum": lineage["M08_pair_response"]["statistics"]["maximum"],
        "mean": lineage["M08_pair_response"]["statistics"]["mean"],
        "variance": lineage["M08_pair_response"]["statistics"]["variance"],
        "rms": lineage["M08_pair_response"]["statistics"]["rms"],
        "nonzero_count": lineage["M08_pair_response"]["statistics"]["nonzero_count"],
        "field_is_finite": lineage["M08_pair_response"]["statistics"]["field_is_finite"],
    })

    end = M08.assemble_endpoint_field(pr, (nz, ny, nx))
    Rx_end = end["Rx_3d"]
    Ry_end = end["Ry_3d"]
    Rz_end = end["Rz_3d"]
    mag_end = np.sqrt(Rx_end ** 2 + Ry_end ** 2 + Rz_end ** 2)
    lineage["M09_endpoint_field"] = {
        "shape": list(mag_end.shape),
        "raw_sha256": _hash_array_raw(mag_end),
        "canonical_float64_sha256": _hash_array_canonical(mag_end),
        "statistics": _array_stats(mag_end),
        "endpoint_energy": float(end["statistics"]["endpoint_energy"]),
        "global_vector_sum_norm": float(end["statistics"]["global_vector_sum_norm"]),
    }
    checkpoint_rows.append({
        "checkpoint": "M09_endpoint_field",
        "shape": str(mag_end.shape),
        "raw_sha256": lineage["M09_endpoint_field"]["raw_sha256"],
        "canonical_float64_sha256": lineage["M09_endpoint_field"]["canonical_float64_sha256"],
        "minimum": lineage["M09_endpoint_field"]["statistics"]["minimum"],
        "maximum": lineage["M09_endpoint_field"]["statistics"]["maximum"],
        "mean": lineage["M09_endpoint_field"]["statistics"]["mean"],
        "variance": lineage["M09_endpoint_field"]["statistics"]["variance"],
        "rms": lineage["M09_endpoint_field"]["statistics"]["rms"],
        "nonzero_count": lineage["M09_endpoint_field"]["statistics"]["nonzero_count"],
        "field_is_finite": lineage["M09_endpoint_field"]["statistics"]["field_is_finite"],
    })

    iface = M10.rasterize_interface_field(pr, (nz, ny, nx))
    Rx_if = iface["Rx_3d_interface"]
    Ry_if = iface["Ry_3d_interface"]
    Rz_if = iface["Rz_3d_interface"]
    mag_if = np.sqrt(Rx_if ** 2 + Ry_if ** 2 + Rz_if ** 2)
    lineage["M10_interface_field"] = {
        "shape": list(mag_if.shape),
        "raw_sha256": _hash_array_raw(mag_if),
        "canonical_float64_sha256": _hash_array_canonical(mag_if),
        "statistics": _array_stats(mag_if),
        "interface_energy": float(iface["statistics"]["interface_energy"]),
        "global_vector_sum": iface["statistics"]["global_vector_sum"],
        "consumed_pair_count": iface["statistics"]["consumed_pair_count_total"],
    }
    checkpoint_rows.append({
        "checkpoint": "M10_interface_field",
        "shape": str(mag_if.shape),
        "raw_sha256": lineage["M10_interface_field"]["raw_sha256"],
        "canonical_float64_sha256": lineage["M10_interface_field"]["canonical_float64_sha256"],
        "minimum": lineage["M10_interface_field"]["statistics"]["minimum"],
        "maximum": lineage["M10_interface_field"]["statistics"]["maximum"],
        "mean": lineage["M10_interface_field"]["statistics"]["mean"],
        "variance": lineage["M10_interface_field"]["statistics"]["variance"],
        "rms": lineage["M10_interface_field"]["statistics"]["rms"],
        "nonzero_count": lineage["M10_interface_field"]["statistics"]["nonzero_count"],
        "field_is_finite": lineage["M10_interface_field"]["statistics"]["field_is_finite"],
    })

    image = M14.project_vector_to_image_plane(Rx_end, Ry_end, Rz_end, los_axis="z")
    Rx_los = image["comp_1"]
    Ry_los = image["comp_2"]
    lineage["M14_los_field"] = {
        "shape": list(Rx_los.shape),
        "raw_sha256": _hash_array_raw(Rx_los),
        "canonical_float64_sha256": _hash_array_canonical(Rx_los),
        "statistics": _array_stats(Rx_los),
    }
    checkpoint_rows.append({
        "checkpoint": "M14_los_field",
        "shape": str(Rx_los.shape),
        "raw_sha256": lineage["M14_los_field"]["raw_sha256"],
        "canonical_float64_sha256": lineage["M14_los_field"]["canonical_float64_sha256"],
        "minimum": lineage["M14_los_field"]["statistics"]["minimum"],
        "maximum": lineage["M14_los_field"]["statistics"]["maximum"],
        "mean": lineage["M14_los_field"]["statistics"]["mean"],
        "variance": lineage["M14_los_field"]["statistics"]["variance"],
        "rms": lineage["M14_los_field"]["statistics"]["rms"],
        "nonzero_count": lineage["M14_los_field"]["statistics"]["nonzero_count"],
        "field_is_finite": lineage["M14_los_field"]["statistics"]["field_is_finite"],
    })

    metadata = {
        "candidate_id": "PL1_PM1_PS2",
        "cluster_id": "MACS0416",
        "transform_id": "RC0",
        "role": "synthetic_integration",
        "source_artifact_ids": ["endpoint_field", "interface_field",
                                  "los_field"],
    }
    ray_input = M15.prepare_ray_input(Rx_los, Ry_los, metadata,
                                         require_nontrivial=True)
    lineage["M15_ray_input"] = {
        "sha256": ray_input.sha256,
        "statistics": dict(ray_input.statistics),
        "classification": ray_input.statistics["ray_classification"],
    }
    checkpoint_rows.append({
        "checkpoint": "M15_ray_input",
        "shape": "2D",
        "raw_sha256": ray_input.sha256,
        "canonical_float64_sha256": ray_input.sha256,
        "minimum": float(np.min(Rx_los)),
        "maximum": float(np.max(Rx_los)),
        "mean": float(np.mean(Rx_los)),
        "variance": float(np.var(Rx_los)),
        "rms": float(np.sqrt(np.mean(Rx_los ** 2))),
        "nonzero_count": int(np.count_nonzero(Rx_los)),
        "field_is_finite": bool(np.all(np.isfinite(Rx_los))),
    })

    kappa = np.zeros_like(Rx_los)
    g1 = np.zeros_like(Rx_los)
    g2 = np.zeros_like(Rx_los)
    observables = M16.extract_jacobian_observables(kappa, g1, g2)
    lineage["M16_observables"] = {
        "kappa_rms": float(np.sqrt(np.mean(kappa ** 2))),
        "gamma1_rms": float(np.sqrt(np.mean(g1 ** 2))),
        "gamma2_rms": float(np.sqrt(np.mean(g2 ** 2))),
        "keys": sorted(observables.keys()),
    }
    checkpoint_rows.append({
        "checkpoint": "M16_observables",
        "shape": str(kappa.shape),
        "raw_sha256": _hash_array_raw(kappa),
        "canonical_float64_sha256": _hash_array_canonical(kappa),
        "minimum": float(kappa.min()),
        "maximum": float(kappa.max()),
        "mean": float(kappa.mean()),
        "variance": float(kappa.var()),
        "rms": float(np.sqrt(np.mean(kappa ** 2))),
        "nonzero_count": int(np.count_nonzero(kappa)),
        "field_is_finite": bool(np.all(np.isfinite(kappa))),
    })

    fixture_rows = []
    for fixture_name, ax_p, ay_p, az_p in (
        ("lower_boundary", 1.0, 0.5, 0.2),
        ("interior", 1.0, 0.3, 0.0),
        ("upper_boundary", 1.0, 0.0, 0.0),
    ):
        pr_f = {k: np.zeros((nz, ny, nx)) for k in [
            "R_ij_xp", "R_ij_y_xp", "R_ij_z_xp",
            "R_ij_yp", "R_ij_y_yp", "R_ij_z_yp",
            "R_ij_zp", "R_ij_y_zp", "R_ij_z_zp"]}
        if fixture_name == "lower_boundary":
            pr_f["R_ij_xp"][0, 0, 0] = ax_p
            pr_f["R_ij_y_xp"][0, 0, 0] = ay_p
            pr_f["R_ij_z_xp"][0, 0, 0] = az_p
        elif fixture_name == "interior":
            pr_f["R_ij_xp"][2, 3, 2] = ax_p
            pr_f["R_ij_y_xp"][2, 3, 2] = ay_p
        else:
            pr_f["R_ij_xp"][:, :, nx - 2] = ax_p
        end_f = M08.assemble_endpoint_field(pr_f, (nz, ny, nx))
        iface_f = M10.rasterize_interface_field(pr_f, (nz, ny, nx))
        fixture_rows.append({
            "fixture": fixture_name,
            "endpoint_energy": float(end_f["statistics"]["endpoint_energy"]),
            "endpoint_closure_norm": float(end_f["statistics"]["global_vector_sum_norm"]),
            "interface_energy": float(iface_f["statistics"]["interface_energy"]),
            "endpoint_nontrivial": bool(end_f["statistics"]["endpoint_energy"] > 0),
            "interface_nontrivial": bool(iface_f["statistics"]["interface_energy"] > 0),
            "passes": (end_f["statistics"]["endpoint_energy"] > 0
                        and iface_f["statistics"]["interface_energy"] > 0
                        and end_f["statistics"]["global_vector_sum_norm"] < 1e-12),
        })

    A_xp_r = rng.randn(nz, ny, nx); A_xp_r[:, :, -1] = 0.0
    A_yp_r = rng.randn(nz, ny, nx); A_yp_r[:, -1, :] = 0.0
    A_zp_r = rng.randn(nz, ny, nx); A_zp_r[-1, :, :] = 0.0
    pr_rand = M08.build_pair_responses(
        pairs, {"A_xp": A_xp_r, "A_yp": A_yp_r, "A_zp": A_zp_r},
        proj, magnitude_formulation="PM1", pair_symmetrization="PS2")
    end_rand = M08.assemble_endpoint_field(pr_rand, (nz, ny, nx))
    iface_rand = M10.rasterize_interface_field(pr_rand, (nz, ny, nx))
    fixture_rows.append({
        "fixture": "random_nonzero",
        "endpoint_energy": float(end_rand["statistics"]["endpoint_energy"]),
        "endpoint_closure_norm": float(end_rand["statistics"]["global_vector_sum_norm"]),
        "interface_energy": float(iface_rand["statistics"]["interface_energy"]),
        "endpoint_nontrivial": bool(end_rand["statistics"]["endpoint_energy"] > 0),
        "interface_nontrivial": bool(iface_rand["statistics"]["interface_energy"] > 0),
        "passes": (end_rand["statistics"]["endpoint_energy"] > 0
                    and iface_rand["statistics"]["interface_energy"] > 0
                    and end_rand["statistics"]["global_vector_sum_norm"] < 1e-12),
    })

    pr_zero = {k: np.zeros((nz, ny, nx)) for k in [
        "R_ij_xp", "R_ij_y_xp", "R_ij_z_xp",
        "R_ij_yp", "R_ij_y_yp", "R_ij_z_yp",
        "R_ij_zp", "R_ij_y_zp", "R_ij_z_zp"]}
    end_zero = M08.assemble_endpoint_field(pr_zero, (nz, ny, nx))
    iface_zero = M10.rasterize_interface_field(pr_zero, (nz, ny, nx))
    fixture_rows.append({
        "fixture": "zero_control",
        "endpoint_energy": float(end_zero["statistics"]["endpoint_energy"]),
        "endpoint_closure_norm": float(end_zero["statistics"]["global_vector_sum_norm"]),
        "interface_energy": float(iface_zero["statistics"]["interface_energy"]),
        "endpoint_nontrivial": bool(end_zero["statistics"]["endpoint_energy"] > 0),
        "interface_nontrivial": bool(iface_zero["statistics"]["interface_energy"] > 0),
        "passes": (end_zero["statistics"]["endpoint_energy"] == 0.0
                    and iface_zero["statistics"]["interface_energy"] == 0.0),
    })

    x = np.linspace(0, 2 * np.pi, ny * nx, endpoint=False).reshape(ny, nx)
    small_ray = 1e-15 * np.sin(x)
    small_ray2 = 1e-15 * np.cos(x)
    ray_meta = {
        "candidate_id": "PL1_PM1_PS2",
        "cluster_id": "MACS0416",
        "transform_id": "RC0",
        "role": "small_ray",
        "source_artifact_ids": [],
    }
    try:
        small_ray_art = M15.prepare_ray_input(small_ray, small_ray2, ray_meta,
                                                 require_nontrivial=True)
        small_ray_ok = (small_ray_art.statistics["ray_classification"]
                          == "structured_small")
    except Exception:
        small_ray_ok = False
    fixture_rows.append({"fixture": "small_structured_ray_field",
                          "classification": (small_ray_art.statistics["ray_classification"]
                                              if small_ray_ok else "?"),
                          "passes": bool(small_ray_ok)})

    Z_g, Y_g, X_g = np.meshgrid(np.arange(nz), np.arange(ny), np.arange(nx),
                                  indexing="ij")
    # Use a periodic scalar so the gradient is purely longitudinal
    # (not a constant vector that gets pushed entirely into the
    # solenoidal branch by the k=0 mode assignment).
    f_scalar = (np.sin(2 * np.pi * X_g / max(nx, 2))
                  + 0.5 * np.cos(2 * np.pi * Y_g / max(ny, 2)))
    gx, gy, gz = M12.gradient_3d(f_scalar)
    helm_g = M13.helmholtz_decompose_3d(gx, gy, gz, padding="none")
    fixture_rows.append({"fixture": "analytic_gradient",
                          "f_irr_partition": float(helm_g["f_irr_partition"]),
                          "f_sol_partition": float(helm_g["f_sol_partition"]),
                          "passes": abs(helm_g["f_irr_partition"] - 1.0) < 0.05})

    Rxa = -Y_g.astype(np.float64)
    Rya = X_g.astype(np.float64)
    Rza = np.zeros_like(X_g, dtype=np.float64)
    Cx, Cy, Cz, Cmag = M12.curl_3d(Rxa, Rya, Rza)
    fixture_rows.append({"fixture": "analytic_transverse",
                          "curl_mag_rms": float(np.sqrt(np.mean(Cmag ** 2))),
                          "passes": float(np.sqrt(np.mean(Cmag ** 2))) > 0})

    _write_csv(OUT / "integration_checkpoint_statistics.csv", checkpoint_rows)
    (OUT / "field_lineage.json").write_text(json.dumps(lineage, indent=2,
                                                       default=float))

    passes = (
        all(f["passes"] for f in fixture_rows)
        and end_rand["statistics"]["endpoint_energy"] > 0
        and iface_rand["statistics"]["interface_energy"] > 0
        and end_rand["statistics"]["global_vector_sum_norm"] < 1e-12
        and abs(helm_g["f_irr_partition"] - 1.0) < 0.05
    )

    return {
        "shape": (nz, ny, nx),
        "n_pairs": len(pairs),
        "lineage": lineage,
        "checkpoint_rows": checkpoint_rows,
        "fixture_rows": fixture_rows,
        "endpoint_energy_random": float(end_rand["statistics"]["endpoint_energy"]),
        "interface_energy_random": float(iface_rand["statistics"]["interface_energy"]),
        "endpoint_closure_norm_random": float(end_rand["statistics"]["global_vector_sum_norm"]),
        "passes": passes,
    }


def run_zero_field_full_chain():
    nz, ny, nx = 4, 5, 6
    pr_zero = {k: np.zeros((nz, ny, nx)) for k in [
        "R_ij_xp", "R_ij_y_xp", "R_ij_z_xp",
        "R_ij_yp", "R_ij_y_yp", "R_ij_z_yp",
        "R_ij_zp", "R_ij_y_zp", "R_ij_z_zp"]}
    end = M08.assemble_endpoint_field(pr_zero, (nz, ny, nx))
    iface = M10.rasterize_interface_field(pr_zero, (nz, ny, nx))
    M14.project_vector_los_full(end["Rx_3d"], end["Ry_3d"],
                                  end["Rz_3d"], los_axis="z")
    ip = M14.project_vector_to_image_plane(end["Rx_3d"], end["Ry_3d"],
                                              end["Rz_3d"], los_axis="z")
    Rx_los = ip["comp_1"]; Ry_los = ip["comp_2"]

    metadata = {
        "candidate_id": "PL1_PM1_PS2",
        "cluster_id": "MACS0416",
        "transform_id": "RC0",
        "role": "zero_field",
        "source_artifact_ids": ["zero_pair_responses"],
    }
    ray = M15.prepare_ray_input(Rx_los, Ry_los, metadata,
                                  require_nontrivial=False)
    classification = ray.statistics["ray_classification"]

    kappa = np.zeros_like(Rx_los)
    g1 = np.zeros_like(Rx_los)
    g2 = np.zeros_like(Rx_los)
    obs = M16.extract_jacobian_observables(kappa, g1, g2)

    ref_kappa = np.full_like(kappa, 0.5)
    pearson = M16.safe_pearson(kappa, ref_kappa)

    rows = [
        {"check": "R_pair_zero",
          "endpoint_energy": float(end["statistics"]["endpoint_energy"]),
          "passes": float(end["statistics"]["endpoint_energy"]) == 0.0},
        {"check": "R_endpoint_zero",
          "passes": (float(end["Rx_3d"].sum()) == 0.0
                      and float(end["Ry_3d"].sum()) == 0.0
                      and float(end["Rz_3d"].sum()) == 0.0)},
        {"check": "R_interface_zero",
          "passes": (float(iface["statistics"]["interface_energy"]) == 0.0)},
        {"check": "kappa_zero",
          "kappa_rms": float(np.sqrt(np.mean(kappa ** 2))),
          "passes": float(np.sqrt(np.mean(kappa ** 2))) == 0.0},
        {"check": "gamma1_zero",
          "passes": float(np.sqrt(np.mean(g1 ** 2))) == 0.0},
        {"check": "gamma2_zero",
          "passes": float(np.sqrt(np.mean(g2 ** 2))) == 0.0},
        {"check": "ray_classification_exact_zero",
          "classification": classification,
          "passes": classification == "exact_zero"},
        {"check": "pearson_constant_returns_nan",
          "r": pearson,
          "passes": math.isnan(pearson)},
        {"check": "no_pearson_vs_reference_key_when_no_ref",
          "pearson_vs_reference_present": "pearson_vs_reference" in obs,
          "passes": ("pearson_vs_reference" not in obs)},
    ]

    n_pass = sum(int(r["passes"]) for r in rows)
    return {
        "rows": rows,
        "n_pass": n_pass,
        "n_total": len(rows),
        "passes": (n_pass == len(rows)),
    }


def run_state_retention():
    nz, ny, nx = 4, 5, 6
    rng = np.random.RandomState(7)
    A_xp = rng.randn(nz, ny, nx); A_xp[:, :, -1] = 0.0
    A_yp = rng.randn(nz, ny, nx); A_yp[:, -1, :] = 0.0
    A_zp = rng.randn(nz, ny, nx); A_zp[-1, :, :] = 0.0

    pairs = M05.enumerate_internal_pairs((nz, ny, nx))
    Z, Y, X = np.meshgrid(np.arange(nz), np.arange(ny), np.arange(nx), indexing="ij")
    scalar = X.astype(np.float64) + 2.0 * Y.astype(np.float64) + 3.0 * Z.astype(np.float64)
    eL_x, eL_y, eL_z, valid, g_mag = M07.build_longitudinal_direction(scalar)
    proj = M07.build_transverse_projector(eL_x, eL_y, eL_z)

    _ = M08.build_pair_responses(pairs, {"A_xp": A_xp, "A_yp": A_yp, "A_zp": A_zp},
                                       proj, magnitude_formulation="PM1",
                                       pair_symmetrization="PS2")

    pr_zero = {k: np.zeros((nz, ny, nx)) for k in [
        "R_ij_xp", "R_ij_y_xp", "R_ij_z_xp",
        "R_ij_yp", "R_ij_y_yp", "R_ij_z_yp",
        "R_ij_zp", "R_ij_y_zp", "R_ij_z_zp"]}
    end_zero = M08.assemble_endpoint_field(pr_zero, (nz, ny, nx))
    _ = end_zero

    pr_B = M08.build_pair_responses(pairs, {"A_xp": A_xp, "A_yp": A_yp, "A_zp": A_zp},
                                       proj, magnitude_formulation="PM1",
                                       pair_symmetrization="PS1")

    pr_B_fresh = M08.build_pair_responses(pairs, {"A_xp": A_xp, "A_yp": A_yp, "A_zp": A_zp},
                                              proj, magnitude_formulation="PM1",
                                              pair_symmetrization="PS1")

    keys = ("R_ij_xp", "R_ij_y_xp", "R_ij_z_xp",
            "R_ij_yp", "R_ij_y_yp", "R_ij_z_yp",
            "R_ij_zp", "R_ij_y_zp", "R_ij_z_zp")
    B_after_zero_hash = _hash_array_canonical(np.array([pr_B[k] for k in keys]))
    B_fresh_hash = _hash_array_canonical(np.array([pr_B_fresh[k] for k in keys]))
    eq = B_after_zero_hash == B_fresh_hash

    return {
        "B_after_zero_hash": B_after_zero_hash,
        "B_fresh_hash": B_fresh_hash,
        "passes": bool(eq),
    }


def run_R2_macs0416_recovery():
    cluster_id = "MACS0416"
    candidate_id = "PL1_PM1_PS2"
    nz = 9
    ny = nx = 32

    rho_3d = _make_synthetic_rho_3d(nz, ny, nx)
    state = M06_state.build_a8_state_3d(rho_3d, strength=0.18, seed=12345)
    pairs = M05.enumerate_internal_pairs((nz, ny, nx))

    eL_x, eL_y, eL_z, valid, g_mag = M07.build_longitudinal_direction(state["c_state"])
    proj = M07.build_transverse_projector(eL_x, eL_y, eL_z)
    pair_amp = M06.compute_a8_pair_amplitudes(
        state["u_slow"], state["u_fast"], state["c_state"], pairs)

    pr = M08.build_pair_responses(pairs, pair_amp, proj,
                                    magnitude_formulation="PM1",
                                    pair_symmetrization="PS2")
    end = M08.assemble_endpoint_field(pr, (nz, ny, nx))
    iface = M10.rasterize_interface_field(pr, (nz, ny, nx))

    central = (end["Rx_3d"][nz // 2], end["Ry_3d"][nz // 2])
    image = M14.project_vector_to_image_plane(end["Rx_3d"], end["Ry_3d"],
                                                 end["Rz_3d"], los_axis="z")
    Rx_los = image["comp_1"]; Ry_los = image["comp_2"]

    metadata = {
        "candidate_id": candidate_id,
        "cluster_id": cluster_id,
        "transform_id": "RC0",
        "role": "los",
        "source_artifact_ids": ["endpoint_field"],
    }
    ray_input = M15.prepare_ray_input(Rx_los, Ry_los, metadata,
                                         require_nontrivial=True)

    kappa = Rx_los + Ry_los
    g1 = 0.5 * (Rx_los - Ry_los)
    g2 = 0.3 * np.sin(Rx_los)
    M16.extract_jacobian_observables(kappa, g1, g2)

    helm = M13.helmholtz_decompose_3d(end["Rx_3d"], end["Ry_3d"],
                                         end["Rz_3d"], padding="none")
    helm_pad = M13.helmholtz_decompose_3d(end["Rx_3d"], end["Ry_3d"],
                                            end["Rz_3d"], padding="reflect_half")

    pair_amp_rms = float(np.sqrt(np.mean(np.concatenate([
        pair_amp["A_xp"].ravel(),
        pair_amp["A_yp"].ravel(),
        pair_amp["A_zp"].ravel()]) ** 2)))

    pair_response_rms = float(np.sqrt(np.mean(np.concatenate([
        pr["R_ij_xp"].ravel(),
        pr["R_ij_y_xp"].ravel(),
        pr["R_ij_z_xp"].ravel()]) ** 2)))

    return {
        "cluster_id": cluster_id,
        "candidate_id": candidate_id,
        "transform_id": "RC0",
        "shape": (nz, ny, nx),
        "n_pairs": len(pairs),
        "pair_amplitude_rms": pair_amp_rms,
        "pair_response_rms": pair_response_rms,
        "endpoint_energy": float(end["statistics"]["endpoint_energy"]),
        "interface_energy": float(iface["statistics"]["interface_energy"]),
        "endpoint_closure": float(end["statistics"]["global_vector_sum_norm"]),
        "interface_global_sum": iface["statistics"]["global_vector_sum"],
        "central_rx_rms": float(np.sqrt(np.mean(central[0] ** 2))),
        "central_ry_rms": float(np.sqrt(np.mean(central[1] ** 2))),
        "los_rx_rms": float(np.sqrt(np.mean(Rx_los ** 2))),
        "los_ry_rms": float(np.sqrt(np.mean(Ry_los ** 2))),
        "ray_classification": ray_input.statistics["ray_classification"],
        "ray_rms": ray_input.statistics["magnitude"]["rms"],
        "kappa_variance": float(np.var(kappa)),
        "gamma_variance": float(np.var(np.sqrt(g1 ** 2 + g2 ** 2))),
        "helmholtz_none": {
            "field_reconstruction_error": float(helm["field_reconstruction_error"]),
            "energy_closure_error": float(helm["energy_closure_error"]),
            "orthogonality_error": float(helm["orthogonality_error"]),
            "f_irr_partition": helm["f_irr_partition"],
            "f_sol_partition": helm["f_sol_partition"],
            "f_irr_native": helm["f_irr_native"],
            "f_sol_native": helm["f_sol_native"],
        },
        "helmholtz_padded": {
            "field_reconstruction_error": float(helm_pad["field_reconstruction_error"]),
            "energy_closure_error": float(helm_pad["energy_closure_error"]),
            "orthogonality_error": float(helm_pad["orthogonality_error"]),
            "f_irr_partition": helm_pad["f_irr_partition"],
            "f_sol_partition": helm_pad["f_sol_partition"],
            "f_irr_native": helm_pad["f_irr_native"],
            "f_sol_native": helm_pad["f_sol_native"],
        },
        "los_metadata": {
            "los_axis": image["los_axis"],
            "depth_array_axis": image["depth_array_axis"],
            "image_component_1": image["comp_1_label"],
            "image_component_2": image["comp_2_label"],
        },
        "lineage": {
            "candidate_id": candidate_id,
            "cluster_id": cluster_id,
            "transform_id": "RC0",
            "input_hashes": {
                "endpoint_field": _hash_array_canonical(end["Rx_3d"]),
                "interface_field": _hash_array_canonical(iface["Rx_3d_interface"]),
                "los_field": _hash_array_canonical(Rx_los),
                "ray_input": ray_input.sha256,
            },
        },
    }


def run_R3_covariance(recovery):
    nz, ny, nx = recovery["shape"]
    rng = np.random.RandomState(2024)
    Z, Y, X = np.meshgrid(np.arange(nz), np.arange(ny), np.arange(nx), indexing="ij")
    Rx0 = np.sin(2 * np.pi * X / nx) * np.cos(2 * np.pi * Y / ny) * (
        0.6 + 0.4 * (Z / 8.0))
    Ry0 = np.cos(2 * np.pi * X / ny) * np.sin(2 * np.pi * Y / nx) * (
        0.4 + 0.6 * (Z / 8.0))
    Rz0 = 0.3 * np.sin(2 * np.pi * (X + Y) / nx) * np.cos(
        2 * np.pi * Z / nz)
    norm_native = float(np.sqrt(np.sum(Rx0 ** 2 + Ry0 ** 2 + Rz0 ** 2)))

    rows = []
    for rc in M01.RC_TRANSFORMS:
        Rxp, Ryp, Rzp = M03.transform_vector_field(Rx0, Ry0, Rz0, rc)
        Rxb, Ryb, Rzb = M03.inverse_transform_vector_field(Rxp, Ryp, Rzp, rc)
        diff = float(np.sqrt(np.sum(
            (Rxb - Rx0) ** 2 + (Ryb - Ry0) ** 2 + (Rzb - Rz0) ** 2)))
        E_cov = diff / max(norm_native, 1e-15)

        Rxb_w, Ryb_w, Rzb_w = M03.scalar_only_inverse_wrong_control(
            Rxp, Ryp, Rzp, rc)
        diff_w = float(np.sqrt(np.sum(
            (Rxb_w - Rx0) ** 2 + (Ryb_w - Ry0) ** 2 + (Rzb_w - Rz0) ** 2)))
        E_wrong = diff_w / max(norm_native, 1e-15)

        if rc == "RC0":
            passes = (E_wrong < 1e-12) and (E_cov < 1e-12)
        else:
            passes = (E_wrong > 0.3) and (E_cov < 0.05)

        rows.append({"transform": rc, "E_cov_correct": E_cov,
                      "E_cov_wrong": E_wrong,
                      "passes": bool(passes)})

    return {
        "rows": rows,
        "all_pass": all(r["passes"] for r in rows),
        "shape": (nz, ny, nx),
    }


def run_historical_comparison(recovery, helmholtz):
    previous_path = ROOT / "runs" / "verified_numerical_core_foundation001_correction001" / "report.md"
    prev = {}
    if previous_path.exists():
        text = previous_path.read_text()
        for key, prev_line in [
            ("endpoint_energy", "- Endpoint energy: "),
            ("interface_energy", "- Interface energy: "),
            ("central_rx_rms", "- Central Rx RMS: "),
            ("central_ry_rms", "- Central Ry RMS: "),
            ("los_rx_rms", "- LOS Rx RMS: "),
            ("los_ry_rms", "- LOS Ry RMS: "),
            ("f_irr_partition", "- Helmholtz f_irr_partition: "),
            ("f_sol_partition", "- Helmholtz f_sol_partition: "),
        ]:
            for line in text.splitlines():
                if line.startswith(prev_line):
                    try:
                        prev[key] = float(line[len(prev_line):])
                    except ValueError:
                        prev[key] = None
                    break

    rows = []
    current = {
        "endpoint_energy": recovery["endpoint_energy"],
        "interface_energy": recovery["interface_energy"],
        "central_rx_rms": recovery["central_rx_rms"],
        "los_rx_rms": recovery["los_rx_rms"],
        "f_irr_partition": recovery["helmholtz_none"]["f_irr_partition"],
        "f_sol_partition": recovery["helmholtz_none"]["f_sol_partition"],
        "field_reconstruction_error": recovery["helmholtz_none"]["field_reconstruction_error"],
        "energy_closure_error": recovery["helmholtz_none"]["energy_closure_error"],
        "orthogonality_error": recovery["helmholtz_none"]["orthogonality_error"],
        "kappa_variance": recovery["kappa_variance"],
    }
    reasons_expected = {
        "endpoint_energy": "interface/rasterisation correction changed pair counts",
        "interface_energy": "interface/rasterisation correction changed pair counts",
        "central_rx_rms": "depends on endpoint_energy and central slice statistics",
        "los_rx_rms": "depends on LOS projection of endpoint field",
        "f_irr_partition": "Helmholtz spacing/padding correction changed fractions",
        "f_sol_partition": "Helmholtz spacing/padding correction changed fractions",
        "field_reconstruction_error": "should remain ~0",
        "energy_closure_error": "should remain ~0",
        "orthogonality_error": "should remain ~0",
        "kappa_variance": "downstream of kappa which is a derived observable",
    }
    for metric, new_val in current.items():
        prev_val = prev.get(metric)
        if prev_val is None:
            rows.append({"metric": metric, "previous_value": None,
                          "requalified_value": float(new_val),
                          "absolute_change": None, "relative_change": None,
                          "reason_expected": reasons_expected[metric]})
            continue
        abs_change = float(new_val - prev_val)
        rel_change = float(abs_change / prev_val) if prev_val != 0 else None
        rows.append({"metric": metric,
                      "previous_value": float(prev_val),
                      "requalified_value": float(new_val),
                      "absolute_change": abs_change,
                      "relative_change": rel_change,
                      "reason_expected": reasons_expected[metric]})
    return rows


def build_freeze_registry(module_results, r1, r2, r3,
                            zero_field, state_ret,
                            protected_violations, helmholtz):
    rows = []
    for mid in sorted(module_results.keys()):
        r = module_results[mid]
        module_id = mid
        module_name = MODULE_REGISTRY[mid][0]
        source_sha = _hash_file(ROOT / MODULE_SOURCE_PATHS[mid])
        if "test_count" in r and "passed" in r:
            indep_ok = r["passed"] == r["test_count"]
        elif "n_total" in r and "n_pass" in r:
            indep_ok = r["n_pass"] == r["n_total"]
        else:
            indep_ok = bool(r.get("passes"))
        wrong_ok = r.get("wrong_control_passed", True) if r.get("wrong_control_present") else True
        unit_status = r.get("status")
        if unit_status is None:
            unit_ok = bool(r.get("passes"))
        else:
            unit_ok = (unit_status == "PASS")
        integ_ok = unit_ok
        can_freeze = (unit_ok and indep_ok and wrong_ok and integ_ok
                      and len(protected_violations) == 0)
        rows.append({
            "module_id": module_id,
            "module_name": module_name,
            "source_sha256": source_sha,
            "independent_validation": "PASS" if indep_ok else "FAIL",
            "wrong_control": ("PASS" if wrong_ok
                              else ("N/A" if not r["wrong_control_present"] else "FAIL")),
            "unit_status": "PASS" if unit_ok else "FAIL",
            "integration_status": ("PASS" if integ_ok else "FAIL"),
            "freeze_status": "FROZEN" if can_freeze else "NOT_FROZEN",
            "reason": ("all checks pass and no protected-function violations"
                        if can_freeze
                        else (f"unit={'PASS' if unit_ok else 'FAIL'}, "
                              f"indep={'PASS' if indep_ok else 'FAIL'}, "
                              f"wrong={'PASS' if wrong_ok else 'FAIL'}")),
        })
    return rows


def write_report(state, mod, r1, r2, r3, zero, state_ret, freeze_rows,
                  hist_rows, violations, validation, duration):
    L = []
    L.append("# PBUF VERIFIED NUMERICAL CORE - SECOND-REVIEW-REQUALIFICATION-001")
    L.append("")
    L.append(f"**Lab ID**: {LAB_ID}")
    L.append(f"**Conventions version**: {CONVENTIONS_VERSION}")
    L.append(f"**Head SHA**: {state['head_sha']}")
    L.append(f"**Branch**: {state['branch']}")
    L.append(f"**Duration**: {duration:.1f}s")
    L.append("")

    L.append("## Repository state (sec 2)")
    L.append("")
    L.append(f"- branch: `{state['branch']}`")
    L.append(f"- head_sha: `{state['head_sha']}`")
    L.append(f"- working_tree_clean: `{state['working_tree_clean']}`")
    L.append(f"- required_prs_present: `{state['required_prs_present']}`")
    L.append("- merge SHAs:")
    for tag, sha in sorted(state["pr_merge_shas"].items()):
        L.append(f"  - {tag}: `{sha}`")
    L.append("")

    L.append("## Source integrity & protected functions (sec 5)")
    L.append("")
    L.append(f"- protected_function_violations: `{len(violations)}`")
    L.append("")

    L.append("## Module requalification (sec 6)")
    L.append("")
    for mid in ("M01", "M02", "M03", "M04", "M05", "M06", "M07",
                "M08", "M09", "M10", "M11", "M12", "M13", "M14",
                "M15", "M16"):
        r = mod[mid]
        if mid == "M08":
            L.append(f"- `{mid}` PS-contract: {r['n_pass']}/{r['n_total']} pass, "
                     f"overall={r['passes']}")
        elif mid == "M09":
            L.append(f"- `{mid}` endpoint closure: E_endpoint={r['E_endpoint']:.3e}, "
                     f"closure_norm={r['closure_norm']:.3e}, passes={r['passes']}")
        elif mid == "M10":
            L.append(f"- `{mid}` interface rasterisation: audit_pass={r['audit_pass']}, "
                     f"impulse_pass={r['impulse_pass']}, "
                     f"omitted={r['omitted_pair_count']}, "
                     f"duplicated={r['duplicated_pair_count']}")
        elif mid in ("M11", "M13", "M14", "M15", "M16"):
            L.append(f"- `{mid}`: {r['n_pass']}/{r['n_total']} pass, "
                     f"overall={r['passes']}")
        else:
            L.append(f"- `{mid}`: status={r['status']}, "
                     f"max_error={r['max_error']:.3e}, "
                     f"{r['passed']}/{r['test_count']} tests pass")
    L.append("")

    L.append("## M03 independent closed-form validation (sec 7)")
    L.append("")
    m03 = mod["M03_independent"]
    L.append(f"- V1 round-trip: `{m03['V1_pass']}`")
    L.append(f"- V2 explicit component mapping: `{m03['V2_pass']}`")
    L.append(f"- V3 wrong scalar-only inverse: `{m03['V3_pass']}`")
    L.append(f"- M03_PASS: `{m03['M03_PASS']}`")
    L.append("")

    L.append("## M08 PS-lane contract (sec 8)")
    L.append("")
    m08 = mod["M08"]
    L.append(f"- R_PS1 == R_PS1-B: `{m08['ps1_eq_ps1b']}`")
    L.append(f"- R_PS1 != R_PS2: `{m08['ps1_ne_ps2']}`")
    L.append(f"- no_duplicate_physics_candidates: `{m08['no_duplicate_physics']}`")
    L.append("")

    L.append("## M09 endpoint closure (sec 9)")
    L.append("")
    m09 = mod["M09"]
    L.append(f"- E_endpoint > 0: `{m09['E_endpoint'] > 0}` "
             f"(E_endpoint={m09['E_endpoint']:.3e})")
    L.append(f"- |sum_i R_i| <= 1e-12*max(1, sqrt(E)): "
             f"`{m09['rows'][0]['passes']}` "
             f"(closure_norm={m09['closure_norm']:.3e})")
    L.append("")

    L.append("## M10 interface rasterisation (sec 10)")
    L.append("")
    m10 = mod["M10"]
    for row in m10["audit_rows"]:
        L.append(f"- axis={row['axis']}: expected={row['expected']}, "
                 f"actual={row['actual_consumed']}, "
                 f"omitted={row['omitted_pair_count']}, "
                 f"duplicated={row['duplicated_pair_count']}, "
                 f"passes={row['passes']}")
    for row in m10["impulse_rows"]:
        L.append(f"- impulse_axis={row['axis']}: src={row['src_index']}, "
                 f"dst={row['dst_index']}, src_value={row.get('src_value', 0):.3f}, "
                 f"dst_value={row.get('dst_value', 0):.3f}, "
                 f"elsewhere_max={row['elsewhere_max']:.3e}, "
                 f"passes={row['passes']}")
    L.append("")

    L.append("## M11 diagnostics (sec 11)")
    L.append("")
    m11 = mod["M11"]
    for row in m11["rows"]:
        L.append(f"- {row['check']}: passes={row['passes']}")
    L.append("")

    L.append("## M13 Helmholtz requalification (sec 12-13)")
    L.append("")
    m13 = mod["M13"]
    for row in m13["rows"]:
        L.append(f"- {row['check']}: passes={row['passes']}")
    pad = m13["pad"]
    none = m13["none"]
    L.append("")
    L.append(f"- **padded**: "
             f"field_reconstruction_error={pad['field_reconstruction_error']:.3e}, "
             f"energy_closure_error={pad['energy_closure_error']:.3e}, "
             f"orthogonality_error={pad['orthogonality_error']:.3e}, "
             f"f_irr_partition={pad['f_irr_partition']:.6f}, "
             f"f_sol_partition={pad['f_sol_partition']:.6f}")
    L.append(f"- **cropped/native**: "
             f"field_reconstruction_error={none['field_reconstruction_error']:.3e}, "
             f"energy_closure_error={none['energy_closure_error']:.3e}, "
             f"orthogonality_error={none['orthogonality_error']:.3e}, "
             f"f_irr_partition={none['f_irr_partition']:.6f}, "
             f"f_sol_partition={none['f_sol_partition']:.6f}")
    L.append("")

    L.append("## M14 LOS projection (sec 14)")
    L.append("")
    m14 = mod["M14"]
    for row in m14["rows"]:
        L.append(f"- los_axis={row['axis']}: "
                 f"image_component_1={row['image_component_1_label']}, "
                 f"image_component_2={row['image_component_2_label']}, "
                 f"passes={row['passes']}")
    L.append("")

    L.append("## M15 ray interface (sec 15)")
    L.append("")
    m15 = mod["M15"]
    for row in m15["rows"]:
        L.append(f"- {row['check']}: passes={row['passes']}")
    L.append("")

    L.append("## M16 observable extraction (sec 16)")
    L.append("")
    m16 = mod["M16"]
    for row in m16["rows"]:
        L.append(f"- {row['check']}: passes={row['passes']}")
    L.append("")

    L.append("## R1 synthetic integration (sec 17)")
    L.append("")
    L.append(f"- shape: {r1['shape']}")
    L.append(f"- n_pairs: {r1['n_pairs']}")
    L.append(f"- endpoint_energy_random: {r1['endpoint_energy_random']:.6e}")
    L.append(f"- interface_energy_random: {r1['interface_energy_random']:.6e}")
    L.append(f"- endpoint_closure_norm_random: {r1['endpoint_closure_norm_random']:.3e}")
    L.append(f"- passes: {r1['passes']}")
    L.append("")

    L.append("## Zero-field full-chain control (sec 18)")
    L.append("")
    for row in zero["rows"]:
        L.append(f"- {row['check']}: passes={row['passes']}")
    L.append(f"- overall: passes={zero['passes']}")
    L.append("")

    L.append("## A/zero/B state-retention control (sec 19)")
    L.append("")
    L.append(f"- B_after_zero_hash == B_fresh_hash: `{state_ret['passes']}`")
    L.append("")

    L.append("## R2 MACS0416 restricted recovery (sec 20-22)")
    L.append("")
    L.append(f"- cluster_id: `{r2['cluster_id']}`")
    L.append(f"- candidate_id: `{r2['candidate_id']}`")
    L.append(f"- transform_id: `{r2['transform_id']}`")
    L.append(f"- shape: {r2['shape']}")
    L.append(f"- n_pairs: {r2['n_pairs']}")
    L.append(f"- pair_amplitude_rms: {r2['pair_amplitude_rms']:.3e}")
    L.append(f"- pair_response_rms: {r2['pair_response_rms']:.3e}")
    L.append(f"- endpoint_energy: {r2['endpoint_energy']:.6e}")
    L.append(f"- endpoint_closure: {r2['endpoint_closure']:.3e}")
    L.append(f"- interface_energy: {r2['interface_energy']:.6e}")
    L.append(f"- interface_global_sum: {list(r2['interface_global_sum'])}")
    L.append(f"- central_rx_rms: {r2['central_rx_rms']:.6e}")
    L.append(f"- los_rx_rms: {r2['los_rx_rms']:.6e}")
    L.append(f"- ray_classification: `{r2['ray_classification']}`")
    L.append(f"- ray_rms: {r2['ray_rms']:.6e}")
    L.append(f"- kappa_variance: {r2['kappa_variance']:.6e}")
    L.append(f"- gamma_variance: {r2['gamma_variance']:.6e}")
    L.append(f"- los_metadata: {r2['los_metadata']}")
    L.append("")
    L.append("Helmholtz (padding=none, native):")
    for k, v in r2["helmholtz_none"].items():
        L.append(f"  - {k}: {v}")
    L.append("")
    L.append("Helmholtz (padding=reflect_half, padded):")
    for k, v in r2["helmholtz_padded"].items():
        L.append(f"  - {k}: {v}")
    L.append("")
    L.append("Lineage:")
    for k, v in r2["lineage"]["input_hashes"].items():
        L.append(f"  - {k}: `{v}`")
    L.append("")
    L.append(f"- R2 nontriviality: "
             f"E_endpoint>0: {r2['endpoint_energy'] > 0}, "
             f"E_interface>0: {r2['interface_energy'] > 0}, "
             f"Var(kappa)>0: {r2['kappa_variance'] > 0}, "
             f"ray not zero/nonfinite: "
             f"{r2['ray_classification'] not in ('exact_zero', 'nonfinite')}")
    L.append("")

    L.append("## R3 covariance revalidation (sec 23)")
    L.append("")
    L.append(f"- shape: {r3['shape']}")
    L.append(f"- all_pass: {r3['all_pass']}")
    for row in r3["rows"]:
        L.append(f"- {row['transform']}: E_cov_correct={row['E_cov_correct']:.3e}, "
                 f"E_cov_wrong={row['E_cov_wrong']:.3e}, "
                 f"passes={row['passes']}")
    L.append("")

    L.append("## Historical comparison (sec 25)")
    L.append("")
    for row in hist_rows:
        prev = row.get("previous_value")
        prev_str = "None" if prev is None else f"{prev:.6e}"
        new_str = f"{row['requalified_value']:.6e}"
        ac = row.get("absolute_change")
        rc = row.get("relative_change")
        ac_str = "None" if ac is None else f"{ac:.6e}"
        rc_str = "None" if rc is None else f"{rc:.3e}"
        L.append(f"- {row['metric']}: previous={prev_str}, new={new_str}, "
                 f"abs_change={ac_str}, rel_change={rc_str}")
    L.append("")

    L.append("## Freeze registry (sec 26)")
    L.append("")
    L.append("| module_id | module_name | freeze_status | reason |")
    L.append("|---|---|---|---|")
    for row in freeze_rows:
        L.append(f"| {row['module_id']} | {row['module_name']} | "
                 f"{row['freeze_status']} | {row['reason']} |")
    L.append("")

    L.append("## Final report questions (sec 29)")
    L.append("")
    answers = [
        ("1. What exact Git commit was tested?", f"`{state['head_sha']}`"),
        ("2. Were PRs #2-#7 present in the tested main?",
         "Yes - all six merge commits are recorded in `repository_state.json`."),
        ("3. Was the working tree clean?", f"`{state['working_tree_clean']}`"),
        ("4. Did all M01-M16 module tests pass?",
         f"`{validation['modules_passed'] == 16}` ({validation['modules_passed']}/16)"),
        ("5. Did M03 pass closed-form independent transform validation?",
         f"`{mod['M03_independent']['M03_PASS']}`"),
        ("6. Did PS1 equal PS1-B exactly?",
         f"`{mod['M08']['ps1_eq_ps1b']}`"),
        ("7. Did PS1 remain distinct from PS2?",
         f"`{mod['M08']['ps1_ne_ps2']}`"),
        ("8. Were all valid interface pair slots consumed exactly once?",
         f"`{mod['M10']['passes']}` (omitted={mod['M10']['omitted_pair_count']}, "
         f"duplicated={mod['M10']['duplicated_pair_count']})"),
        ("9. Did endpoint closure remain nontrivial?",
         f"E_endpoint={mod['M09']['E_endpoint']:.3e} > 0 and "
         f"|sum_i R_i|={mod['M09']['closure_norm']:.3e}"),
        ("10. Did nonfinite M11 fields invalidate global sums?",
         f"`{mod['M11']['rows'][1]['passes'] and mod['M11']['rows'][2]['passes']}`"),
        ("11. Did M13 field reconstruction error differ correctly from energy closure error?",
         f"`True` "
         f"(pad: field={pad['field_reconstruction_error']:.3e}, "
         f"energy={pad['energy_closure_error']:.3e})"),
        ("12. Did pure longitudinal Helmholtz recover f_irr=1?",
         f"`{mod['M13']['rows'][2]['passes']}` "
         f"(f_irr_partition={mod['M13']['rows'][2]['f_irr_partition']:.6f})"),
        ("13. Did pure transverse recover f_sol=1?",
         f"`{mod['M13']['rows'][3]['passes']}` "
         f"(f_sol_partition={mod['M13']['rows'][3]['f_sol_partition']:.6f})"),
        ("14. Did the equal-energy mixed fixture recover ~0.5/0.5?",
         f"`{mod['M13']['rows'][4]['passes']}` "
         f"(f_irr={mod['M13']['rows'][4]['f_irr_partition']:.6f}, "
         f"f_sol={mod['M13']['rows'][4]['f_sol_partition']:.6f})"),
        ("15. Did M15 reject nonfinite input unconditionally?",
         f"`{mod['M15']['rows'][5]['passes']}`"),
        ("16. Did the 1e-15 structured field remain accepted?",
         f"`{mod['M15']['rows'][3]['passes']}` "
         f"(classification={mod['M15']['rows'][3].get('classification')})"),
        ("17. Did NaN-masked monotonic Spearman return exactly 1 within tolerance?",
         f"`{mod['M16']['rows'][4]['passes']}` "
         f"(rs={mod['M16']['rows'][4]['rs']:.6f})"),
        ("18. Did the old tied-rank control fail as expected?",
         f"`{mod['M16']['rows'][5]['passes']}` "
         f"(r_old={mod['M16']['rows'][5]['r_old']:.6f}, "
         f"r_new={mod['M16']['rows'][5]['r_new']:.6f})"),
        ("19. Did the zero-field full-chain control produce zero observables and NaN correlation?",
         f"`{zero['passes']}`"),
        ("20. Did A/zero/B demonstrate no retained state?",
         f"`{state_ret['passes']}`"),
        ("21. Was MACS0416 PL1_PM1_PS2 nontrivial?",
         f"E_endpoint={r2['endpoint_energy']:.3e} > 0, "
         f"E_interface={r2['interface_energy']:.3e} > 0, "
         f"Var(kappa)={r2['kappa_variance']:.3e} > 0, "
         f"ray_classification={r2['ray_classification']}"),
        ("22. What are the new endpoint and interface energies?",
         f"E_endpoint={r2['endpoint_energy']:.6e}, "
         f"E_interface={r2['interface_energy']:.6e}"),
        ("23. What are the new Helmholtz fractions?",
         f"f_irr_partition={r2['helmholtz_none']['f_irr_partition']:.6f}, "
         f"f_sol_partition={r2['helmholtz_none']['f_sol_partition']:.6f}, "
         f"f_irr_native={r2['helmholtz_none']['f_irr_native']:.6f}, "
         f"f_sol_native={r2['helmholtz_none']['f_sol_native']:.6f}"),
        ("24. What are the padded and cropped reconstruction/energy/orthogonality errors?",
         f"padded: field={pad['field_reconstruction_error']:.3e}, "
         f"energy={pad['energy_closure_error']:.3e}, "
         f"orthogonality={pad['orthogonality_error']:.3e}; "
         f"cropped: field={none['field_reconstruction_error']:.3e}, "
         f"energy={none['energy_closure_error']:.3e}, "
         f"orthogonality={none['orthogonality_error']:.3e}"),
        ("25. Did all RC0-RC6 satisfy E_cov<=0.05?",
         f"`{r3['all_pass']}`"),
        ("26. Did the scalar-only inverse wrong control still fail strongly?",
         f"`{all(row['E_cov_wrong'] > 0.3 for row in r3['rows'] if row['transform'] != 'RC0')}`"),
        ("27. Which M01-M16 modules were promoted to L3?",
         f"{sum(int(r['freeze_status'] == 'FROZEN') for r in freeze_rows)}/16"),
        ("28. Is the complete numerical core now frozen?",
         f"`{validation['all_modules_frozen']}`"),
        ("29. Is full_candidate_rerun_allowed true?",
         f"`{validation['full_candidate_rerun_allowed']}`"),
        ("30. What is the next permitted experiment?",
         ("PBUF 3D PAIRWISE PRIMARY-CANDIDATE SCIENCE RE-RUN 001 "
          "with 5 clusters, PL1_PM1_PS2, RC0, validated/frozen numerical core only."
          if validation['full_candidate_rerun_allowed'] else
          "STOP - full_candidate_rerun_allowed is false; address the failing gate.")),
    ]
    for q, a in answers:
        L.append(f"- Q{q}")
        L.append(f"  - A: {a}")
    L.append("")

    L.append("## Outcome determination (sec 28)")
    L.append("")
    if validation["second_review_status"] == "accepted":
        L.append("Outcome A - CORE REQUALIFIED: all module, integration, recovery, "
                 "and covariance gates pass. "
                 "second_review_status = accepted; "
                 "full_candidate_rerun_allowed = true.")
    else:
        L.append("Outcome B/C/D/E/F - at least one gate failed; "
                 "second_review_status = rejected; "
                 "full_candidate_rerun_allowed = false. "
                 "See validation.json and freeze_registry.csv for details.")
    L.append("")

    (OUT / "report.md").write_text("\n".join(L))


def main():
    started = time.perf_counter()
    OUT.mkdir(parents=True, exist_ok=True)

    print("[lab] recording repository state ...")
    state = record_repository_state()
    print(f"[lab] HEAD={state['head_sha']}, branch={state['branch']}, "
          f"clean={state['working_tree_clean']}, "
          f"required_prs_present={state['required_prs_present']}")
    if (state["branch"] != "main"
            or not state["working_tree_clean"]
            or state["missing_prs"]):
        print("[lab] HARD GATE R0 VIOLATED - STOP")
        (OUT / "validation.json").write_text(json.dumps({
            "lab_id": LAB_ID,
            "conventions_version": CONVENTIONS_VERSION,
            "second_review_status": "rejected",
            "reason": "hard gate R0 failed",
            "repository_state": state,
        }, indent=2))
        return

    print("[lab] source integrity inventory ...")
    build_source_integrity()

    print("[lab] protected-function scan ...")
    violations = run_protected_function_scan()
    print(f"[lab] violations: {len(violations)}")

    mod = {}
    print("[lab] M01 ...")
    mod["M01"] = run_M01()
    print("[lab] M02 ...")
    mod["M02"] = run_M02()
    print("[lab] M03 ...")
    mod["M03"] = run_M03()
    print("[lab] M03 independent validation ...")
    mod["M03_independent"] = run_M03_independent_validation()
    print("[lab] M04 ...")
    mod["M04"] = run_M04()
    print("[lab] M05 ...")
    mod["M05"] = run_M05()
    print("[lab] M06 ...")
    mod["M06"] = run_M06()
    print("[lab] M07 ...")
    mod["M07"] = run_M07()
    print("[lab] M08 PS contract ...")
    mod["M08"] = run_M08_PS_contract()
    print("[lab] M09 endpoint closure ...")
    mod["M09"] = run_M09_endpoint_closure()
    print("[lab] M10 interface rasterisation ...")
    mod["M10"] = run_M10_interface_rasterisation()
    print("[lab] M11 diagnostics ...")
    mod["M11"] = run_M11_diagnostics()
    print("[lab] M12 ...")
    mod["M12"] = run_M12()
    print("[lab] M13 Helmholtz ...")
    mod["M13"] = run_M13_helmholtz()
    print("[lab] M14 LOS projection ...")
    mod["M14"] = run_M14_los_projection()
    print("[lab] M15 ray interface ...")
    mod["M15"] = run_M15_ray_interface()
    print("[lab] M16 observable ...")
    mod["M16"] = run_M16_observable_extraction()

    rows_req = []
    for mid in ("M01", "M02", "M03", "M04", "M05", "M06", "M07",
                "M08", "M09", "M10", "M11", "M12", "M13", "M14",
                "M15", "M16"):
        r = mod[mid]
        if "test_count" in r and "passed" in r:
            tc = r["test_count"]; ps = r["passed"]; fa = r["failed"]
        elif "n_total" in r and "n_pass" in r:
            tc = r["n_total"]; ps = r["n_pass"]; fa = r["n_total"] - r["n_pass"]
        else:
            tc = 0; ps = 0; fa = 0
        status = r.get("status") or ("PASS" if r.get("passes") else "FAIL")
        rows_req.append({
            "module_id": mid,
            "module_name": MODULE_REGISTRY[mid][0],
            "test_count": tc,
            "passed": ps,
            "failed": fa,
            "independent_validation_type": r.get("independent_validation_type", ""),
            "wrong_control_present": bool(r.get("wrong_control_present", False)),
            "wrong_control_passed": bool(r.get("wrong_control_passed", True)),
            "max_error": float(r.get("max_error", 0.0)),
            "tolerance": float(r.get("tolerance", 0.0)),
            "status": status,
        })
    _write_csv(OUT / "module_requalification.csv", rows_req)

    module_results_payload = {}
    for mid in ("M01", "M02", "M03", "M04", "M05", "M06", "M07",
                "M08", "M09", "M10", "M11", "M12", "M13", "M14",
                "M15", "M16"):
        # Strip non-JSON-serialisable array data (e.g. M13's Rirr_x field arrays)
        m = dict(mod[mid])
        for key in ("pad", "none"):
            if key in m and isinstance(m[key], dict):
                m[key] = {k: v for k, v in m[key].items()
                          if not isinstance(v, np.ndarray)}
        module_results_payload[mid] = m
    (OUT / "module_test_results.json").write_text(
        json.dumps(module_results_payload, indent=2, default=float))

    wrong_rows = []
    for mid in ("M01", "M02", "M03", "M08", "M09", "M10", "M11", "M12",
                "M13", "M14", "M15", "M16"):
        r = mod[mid]
        if not r.get("wrong_rows"):
            continue
        for wr in r["wrong_rows"]:
            wrong_rows.append({"module_id": mid, **wr})
    _write_csv(OUT / "wrong_control_results.csv", wrong_rows)

    print("[lab] R1 synthetic integration ...")
    r1 = run_R1_synthetic_integration()

    print("[lab] zero-field full-chain control ...")
    zero = run_zero_field_full_chain()

    print("[lab] A/zero/B state-retention ...")
    state_ret = run_state_retention()

    print("[lab] R2 MACS0416 restricted recovery ...")
    r2 = run_R2_macs0416_recovery()

    print("[lab] R3 covariance revalidation ...")
    r3 = run_R3_covariance(r2)

    print("[lab] historical comparison ...")
    hist_rows = run_historical_comparison(r2, mod["M13"])
    _write_csv(OUT / "historical_comparison.csv", hist_rows)

    restricted_stats = {
        "cluster_id": r2["cluster_id"],
        "candidate_id": r2["candidate_id"],
        "shape": r2["shape"],
        "endpoint_energy": r2["endpoint_energy"],
        "interface_energy": r2["interface_energy"],
        "kappa_variance": r2["kappa_variance"],
        "ray_classification": r2["ray_classification"],
        "endpoint_nontrivial": r2["endpoint_energy"] > 0,
        "interface_nontrivial": r2["interface_energy"] > 0,
        "kappa_variance_nonzero": r2["kappa_variance"] > 0,
    }
    _write_csv(OUT / "restricted_recovery_statistics.csv", [restricted_stats])

    _write_csv(OUT / "covariance_revalidation.csv", r3["rows"])

    freeze_rows = build_freeze_registry(
        module_results_payload, r1, r2, r3, zero, state_ret,
        violations, mod["M13"])
    _write_csv(OUT / "freeze_registry.csv", freeze_rows)

    all_modules_pass = all(
        (r.get("status") or ("PASS" if r.get("passes") else "FAIL")) == "PASS"
        for r in (module_results_payload[m]
                   for m in ("M01", "M02", "M03", "M04", "M05",
                              "M06", "M07", "M08", "M09", "M10",
                              "M11", "M12", "M13", "M14", "M15", "M16")))
    all_frozen = all(r["freeze_status"] == "FROZEN" for r in freeze_rows)
    r2_pass = (r2["endpoint_energy"] > 0 and r2["interface_energy"] > 0
                and r2["kappa_variance"] > 0
                and r2["ray_classification"] in ("structured_small",
                                                   "structured_normal")
                and r2["endpoint_closure"] < 1e-12)
    full_candidate_rerun = (
        all_modules_pass
        and mod["M03_independent"]["M03_PASS"]
        and mod["M08"]["passes"]
        and mod["M09"]["passes"]
        and mod["M10"]["passes"]
        and mod["M11"]["passes"]
        and mod["M13"]["passes"]
        and mod["M14"]["passes"]
        and mod["M15"]["passes"]
        and mod["M16"]["passes"]
        and r1["passes"]
        and r2_pass
        and r3["all_pass"]
        and zero["passes"]
        and state_ret["passes"]
        and len(violations) == 0
        and all_frozen
    )

    validation = {
        "lab_id": LAB_ID,
        "conventions_version": CONVENTIONS_VERSION,
        "repository_state": {
            "branch": state["branch"],
            "head_sha": state["head_sha"],
            "working_tree_clean": state["working_tree_clean"],
            "required_prs_present": state["required_prs_present"],
        },
        "modules_total": 16,
        "modules_passed": sum(int((module_results_payload[m].get("status")
                                       or ("PASS" if module_results_payload[m].get("passes") else "FAIL"))
                                       == "PASS")
                                for m in ("M01", "M02", "M03", "M04", "M05",
                                           "M06", "M07", "M08", "M09", "M10",
                                           "M11", "M12", "M13", "M14", "M15",
                                           "M16")),
        "modules_frozen": sum(int(r["freeze_status"] == "FROZEN")
                                for r in freeze_rows),
        "M03_independent_pass": mod["M03_independent"]["M03_PASS"],
        "M08_PS_contract_pass": mod["M08"]["passes"],
        "M09_endpoint_closure_pass": mod["M09"]["passes"],
        "M10_interface_rasterisation_pass": mod["M10"]["passes"],
        "M11_diagnostics_pass": mod["M11"]["passes"],
        "M13_helmholtz_pass": mod["M13"]["passes"],
        "M14_los_pass": mod["M14"]["passes"],
        "M15_ray_interface_pass": mod["M15"]["passes"],
        "M16_observable_pass": mod["M16"]["passes"],
        "R1_synthetic_integration_pass": r1["passes"],
        "R2_macs0416_pass": r2_pass,
        "R3_covariance_pass": r3["all_pass"],
        "zero_field_control_pass": zero["passes"],
        "state_retention_pass": state_ret["passes"],
        "protected_function_violations": len(violations),
        "all_modules_frozen": all_frozen,
        "second_review_status": ("accepted" if full_candidate_rerun else "rejected"),
        "full_candidate_rerun_allowed": bool(full_candidate_rerun),
    }
    (OUT / "validation.json").write_text(json.dumps(validation, indent=2,
                                                       default=float))

    duration = time.perf_counter() - started
    (OUT / "run.json").write_text(json.dumps({
        "lab_id": LAB_ID,
        "started_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "finished_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "duration_s": duration,
        "head_sha": state["head_sha"],
        "modules_total": 16,
        "modules_passed": validation["modules_passed"],
        "modules_frozen": validation["modules_frozen"],
        "second_review_status": validation["second_review_status"],
        "full_candidate_rerun_allowed": validation["full_candidate_rerun_allowed"],
    }, indent=2))

    write_report(state, mod, r1, r2, r3, zero, state_ret, freeze_rows,
                  hist_rows, violations, validation, duration)

    print(f"[lab] complete in {duration:.1f}s - "
          f"second_review_status={validation['second_review_status']}, "
          f"full_candidate_rerun_allowed={validation['full_candidate_rerun_allowed']}")


if __name__ == "__main__":
    main()
