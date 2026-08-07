#!/usr/bin/env python3
"""PBUF FOUNDATION — REAL-CANDIDATE COVARIANCE LOCALIZATION LAB 001.

Diagnostic-only lab. Uses real MACS0416 input and the reviewed M01-M10
numerical core. It localizes the first covariance failure before LOS/rays.

No source mutation, no fitting, no tolerance changes, no synthetic substitute.
"""
from __future__ import annotations

import csv
import hashlib
import json
import math
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
from astropy.io import fits

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))

from a8_three_dimensional_projection_lab001 import (
    CLUSTERS,
    PRODUCTION,
    construct_common_proxy,
    construct_rho_3d,
)

from pbuf.core import conventions as M01
from pbuf.core import coordinate_transforms as M02
from pbuf.core import vector_transforms as M03
from pbuf.core import tensor_transforms as M04
from pbuf.core import pair_enumeration as M05
from pbuf.core import pair_transfer as M08
from pbuf.core import midpoint_rasterization as M10
from pbuf.models import a8_state as M06_state
from pbuf.models import a8_pair_amplitude as M06
from pbuf.models import transverse_projector as M07

OUT = ROOT / "runs" / "real_candidate_covariance_localization001"
BENCHMARK = ROOT / "PBUF_benchmark"
LAB_ID = "PBUF-FOUNDATION-REAL-CANDIDATE-COVARIANCE-LOCALIZATION-001"

CLUSTER_ID = "MACS0416"
CANDIDATE_ID = "PL1_PM1_PS2"
NZ = 9
PROFILE = "gaussian"
STRENGTH = 0.18
SEED = 12345
FIRST_FAILURE_TOL = 1.0e-8

CFG = dict(PRODUCTION)


def _git(*args: str) -> str:
    return subprocess.check_output(["git", *args], cwd=str(ROOT), text=True).strip()


def _write_json(name: str, obj) -> None:
    (OUT / name).write_text(
        json.dumps(
            obj,
            indent=2,
            default=lambda o: (
                float(o) if isinstance(o, np.floating)
                else int(o) if isinstance(o, np.integer)
                else bool(o) if isinstance(o, np.bool_)
                else list(o) if isinstance(o, tuple)
                else str(o)
            ),
        )
    )


def _write_csv(name: str, rows: list[dict]) -> None:
    path = OUT / name
    if not rows:
        path.write_text("")
        return
    keys: list[str] = []
    for row in rows:
        for key in row:
            if key not in keys:
                keys.append(key)
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)


def _sha_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _norm_scalar(a: np.ndarray) -> float:
    return float(np.linalg.norm(np.asarray(a, dtype=np.float64).ravel()))


def _err_scalar(ref: np.ndarray, test: np.ndarray) -> float:
    ref = np.asarray(ref, dtype=np.float64)
    test = np.asarray(test, dtype=np.float64)
    if ref.shape != test.shape:
        return float("inf")
    return float(np.linalg.norm((test - ref).ravel()) / max(_norm_scalar(ref), 1e-15))


def _norm_vector(v: tuple[np.ndarray, np.ndarray, np.ndarray]) -> float:
    x, y, z = (np.asarray(a, dtype=np.float64) for a in v)
    return float(np.sqrt(np.sum(x * x + y * y + z * z)))


def _err_vector(ref, test) -> float:
    rx, ry, rz = (np.asarray(a, dtype=np.float64) for a in ref)
    tx, ty, tz = (np.asarray(a, dtype=np.float64) for a in test)
    if rx.shape != tx.shape or ry.shape != ty.shape or rz.shape != tz.shape:
        return float("inf")
    d = float(np.sqrt(np.sum((tx-rx)**2 + (ty-ry)**2 + (tz-rz)**2)))
    return d / max(_norm_vector((rx, ry, rz)), 1e-15)


def _norm_tensor(P) -> float:
    Pxx, Pxy, Pxz, Pyy, Pyz, Pzz = (np.asarray(a, dtype=np.float64) for a in P)
    return float(np.sqrt(np.sum(
        Pxx**2 + Pyy**2 + Pzz**2 + 2.0*(Pxy**2 + Pxz**2 + Pyz**2)
    )))


def _err_tensor(ref, test) -> float:
    if any(np.asarray(a).shape != np.asarray(b).shape for a, b in zip(ref, test)):
        return float("inf")
    num = 0.0
    weights = (1.0, 2.0, 2.0, 1.0, 2.0, 1.0)
    for w, a, b in zip(weights, ref, test):
        num += w * float(np.sum((np.asarray(b)-np.asarray(a))**2))
    return math.sqrt(num) / max(_norm_tensor(ref), 1e-15)


def _classify(e: float) -> str:
    if e <= 1e-12:
        return "machine_precision"
    if e <= 1e-8:
        return "small"
    if e <= 1e-4:
        return "warning"
    return "failure"


def _repo_state() -> dict:
    return {
        "repository": "TheExiledMonk/lab",
        "branch": _git("rev-parse", "--abbrev-ref", "HEAD"),
        "head_sha": _git("rev-parse", "HEAD"),
        "tracked_changes": _git("diff", "--name-only"),
        "staged_changes": _git("diff", "--name-only", "--cached"),
    }


def _load_real_input() -> dict:
    cluster = next(c for c in CLUSTERS if c["id"] == CLUSTER_ID)
    fits_path = (
        BENCHMARK
        / cluster["directory"]
        / f"hlsp_frontier_model_{cluster['slug']}_merten_v1_kappa.fits"
    )
    if not fits_path.exists():
        raise FileNotFoundError(fits_path)
    with fits.open(fits_path) as hdul:
        kappa_native = np.asarray(hdul[0].data, dtype=np.float64)
    rho2 = construct_common_proxy(
        kappa_native,
        bins=CFG["bins"],
        extent=CFG["extent"],
    )
    rho3 = construct_rho_3d(rho2, NZ, profile=PROFILE)
    return {
        "cluster": cluster,
        "fits_path": fits_path,
        "fits_sha256": _sha_file(fits_path),
        "kappa_native": kappa_native,
        "rho2": rho2,
        "rho3": rho3,
    }


def _native_initial_state(rho3: np.ndarray) -> dict:
    rng = np.random.RandomState(SEED)
    eq = STRENGTH * rho3
    u_slow0 = eq.copy()
    noise = M06_state.A8_INIT_INJECTION_NOISE * STRENGTH * rng.randn(*rho3.shape)
    u_fast0 = eq + noise
    return {"rho_3d": rho3.copy(), "u_slow0": u_slow0, "u_fast0": u_fast0}


def _evolve_from_initial(initial: dict) -> dict:
    us, uf, history = M06_state.evolve_a8_transport_3d(
        initial["u_slow0"].copy(),
        initial["u_fast0"].copy(),
        stencil="N6",
        boundary="reflective",
    )
    return {
        "rho_3d": initial["rho_3d"].copy(),
        "u_slow": us,
        "u_fast": uf,
        "c_state": history[-1],
    }


def _candidate(state: dict) -> dict:
    shape = tuple(state["c_state"].shape)
    pairs = M05.enumerate_internal_pairs(shape)
    eLx, eLy, eLz, valid, gmag = M07.build_longitudinal_direction(state["c_state"])
    P = M07.build_transverse_projector(eLx, eLy, eLz)
    amp = M06.compute_a8_pair_amplitudes(
        state["u_slow"], state["u_fast"], state["c_state"], pairs
    )
    resp = M08.build_pair_responses(
        pairs,
        amp,
        P,
        magnitude_formulation="PM1",
        pair_symmetrization="PS2",
    )
    end = M08.assemble_endpoint_field(resp, shape)
    iface = M10.rasterize_interface_field(resp, shape)
    return {
        "shape": shape,
        "pairs": pairs,
        "eL": (eLx, eLy, eLz),
        "valid": valid,
        "gmag": gmag,
        "P": P,
        "amp": amp,
        "resp": resp,
        "end": end,
        "iface": iface,
    }


def _transform_initial(initial: dict, rc: str) -> dict:
    return {
        "rho_3d": M02.transform_scalar_field(initial["rho_3d"], rc),
        "u_slow0": M02.transform_scalar_field(initial["u_slow0"], rc),
        "u_fast0": M02.transform_scalar_field(initial["u_fast0"], rc),
    }


def _index_map(native_shape: tuple[int, int, int], rc: str):
    labels = np.arange(np.prod(native_shape), dtype=np.int64).reshape(native_shape)
    transformed = M02.transform_scalar_field(labels, rc)
    coords = np.indices(transformed.shape, dtype=np.int64).reshape(3, -1).T
    native_to_transformed = np.empty((labels.size, 3), dtype=np.int64)
    native_to_transformed[transformed.ravel().astype(np.int64)] = coords

    def map_index(idx):
        flat = np.ravel_multi_index(tuple(idx), native_shape)
        return tuple(int(x) for x in native_to_transformed[flat])

    return map_index


def _pair_lookup(pairs):
    return {
        frozenset((tuple(p.i_index), tuple(p.j_index))): p
        for p in pairs
    }


def _amp_at(amp: dict, pair) -> float:
    return float(amp[f"A_{pair.axis}"][tuple(pair.i_index)])


def _response_at(resp: dict, pair) -> np.ndarray:
    i = tuple(pair.i_index)
    if pair.axis == "xp":
        keys = ("R_ij_xp", "R_ij_y_xp", "R_ij_z_xp")
    elif pair.axis == "yp":
        keys = ("R_ij_yp", "R_ij_y_yp", "R_ij_z_yp")
    elif pair.axis == "zp":
        keys = ("R_ij_zp", "R_ij_y_zp", "R_ij_z_zp")
    else:
        raise ValueError(pair.axis)
    return np.array([float(resp[k][i]) for k in keys], dtype=np.float64)


def _pair_covariance(native: dict, transformed: dict, rc: str) -> dict:
    map_index = _index_map(native["shape"], rc)
    lookup = _pair_lookup(transformed["pairs"])
    Q = np.asarray(M01.get_coordinate_matrix(rc, inverse=False), dtype=np.float64)

    a_ref = []
    a_good = []
    a_no_swap = []
    a_wrong_sign = []
    r_ref = []
    r_good = []
    r_no_swap = []
    direction_rows = []
    missing = 0
    reversed_count = 0

    for p0 in native["pairs"]:
        i0 = tuple(p0.i_index)
        j0 = tuple(p0.j_index)
        it = map_index(i0)
        jt = map_index(j0)
        pt = lookup.get(frozenset((it, jt)))
        if pt is None:
            missing += 1
            continue

        same_orientation = tuple(pt.i_index) == it and tuple(pt.j_index) == jt
        sign = 1.0 if same_orientation else -1.0
        if sign < 0:
            reversed_count += 1

        a0 = _amp_at(native["amp"], p0)
        at = _amp_at(transformed["amp"], pt)

        a_ref.append(a0)
        a_good.append(sign * at)
        a_no_swap.append(at)
        a_wrong_sign.append(-sign * at)

        r0 = _response_at(native["resp"], p0)
        rt = _response_at(transformed["resp"], pt)
        rback = Q.T @ (sign * rt)
        rback_no_swap = Q.T @ rt

        r_ref.extend(r0.tolist())
        r_good.extend(rback.tolist())
        r_no_swap.extend(rback_no_swap.tolist())

        if len(direction_rows) < 18:
            n0 = np.asarray(p0.direction_xyz, dtype=np.float64)
            nt = Q @ n0
            direction_rows.append({
                "transform": rc,
                "source_direction": p0.axis,
                "mapped_vector_x": float(nt[0]),
                "mapped_vector_y": float(nt[1]),
                "mapped_vector_z": float(nt[2]),
                "canonical_transformed_direction": pt.axis,
                "endpoint_swap": bool(sign < 0),
                "orientation_sign": int(sign),
            })

    def rel_arrays(ref, test):
        ref = np.asarray(ref, dtype=np.float64)
        test = np.asarray(test, dtype=np.float64)
        return float(np.linalg.norm(test-ref) / max(np.linalg.norm(ref), 1e-15))

    return {
        "E_pair_amplitude_oriented": rel_arrays(a_ref, a_good),
        "E_pair_amplitude_no_swap_wrong_control": rel_arrays(a_ref, a_no_swap),
        "E_pair_amplitude_wrong_sign_control": rel_arrays(a_ref, a_wrong_sign),
        "E_pair_response_oriented": rel_arrays(r_ref, r_good),
        "E_pair_response_no_swap_wrong_control": rel_arrays(r_ref, r_no_swap),
        "mapped_pair_count": len(a_ref),
        "missing_pair_count": missing,
        "endpoint_swap_count": reversed_count,
        "direction_rows": direction_rows,
    }


def _scalar_checkpoint(rows, rc, name, ref, transformed):
    back = M02.inverse_transform_scalar_field(transformed, rc)
    e = _err_scalar(ref, back)
    rows.append({
        "transform": rc,
        "checkpoint": name,
        "field_type": "scalar",
        "relative_error": e,
        "classification": _classify(e),
        "passes_1e12": e <= 1e-12,
        "passes_1e8": e <= 1e-8,
    })
    return e


def _vector_checkpoint(rows, rc, name, ref, transformed):
    back = M03.inverse_transform_vector_field(*transformed, rc)
    e = _err_vector(ref, back)
    rows.append({
        "transform": rc,
        "checkpoint": name,
        "field_type": "vector",
        "relative_error": e,
        "classification": _classify(e),
        "passes_1e12": e <= 1e-12,
        "passes_1e8": e <= 1e-8,
    })
    return e


def _tensor_checkpoint(rows, rc, name, ref, transformed):
    back = M04.inverse_transform_symmetric_tensor_field(*transformed, rc)
    e = _err_tensor(ref, back)
    rows.append({
        "transform": rc,
        "checkpoint": name,
        "field_type": "tensor",
        "relative_error": e,
        "classification": _classify(e),
        "passes_1e12": e <= 1e-12,
        "passes_1e8": e <= 1e-8,
    })
    return e


def main() -> int:
    started = time.perf_counter()
    OUT.mkdir(parents=True, exist_ok=True)

    repo = _repo_state()
    _write_json("repository_state.json", repo)
    if repo["branch"] != "main" or repo["tracked_changes"] or repo["staged_changes"]:
        validation = {
            "lab_id": LAB_ID,
            "outcome": "REPOSITORY_GATE_FAILURE",
            "head_sha": repo["head_sha"],
        }
        _write_json("validation.json", validation)
        print(json.dumps(validation, indent=2))
        return 2

    real = _load_real_input()
    _write_json("input_provenance.json", {
        "cluster_id": CLUSTER_ID,
        "candidate_id": CANDIDATE_ID,
        "fits_path": str(real["fits_path"].relative_to(ROOT)),
        "fits_sha256": real["fits_sha256"],
        "native_fits_shape": list(real["kappa_native"].shape),
        "rho2_shape": list(real["rho2"].shape),
        "rho3_shape": list(real["rho3"].shape),
        "nz": NZ,
        "profile": PROFILE,
        "strength": STRENGTH,
        "seed": SEED,
    })

    initial0 = _native_initial_state(real["rho3"])
    state0 = _evolve_from_initial(initial0)
    cand0 = _candidate(state0)

    checkpoint_rows = []
    pair_summary_rows = []
    direction_rows = []
    wrong_rows = []
    first_failure = {}

    for rc in M01.RC_TRANSFORMS:
        initial_t = _transform_initial(initial0, rc)
        state_t = _evolve_from_initial(initial_t)
        cand_t = _candidate(state_t)

        errors = {}
        errors["rho_3d"] = _scalar_checkpoint(
            checkpoint_rows, rc, "rho_3d", state0["rho_3d"], state_t["rho_3d"]
        )
        errors["u_slow"] = _scalar_checkpoint(
            checkpoint_rows, rc, "u_slow", state0["u_slow"], state_t["u_slow"]
        )
        errors["u_fast"] = _scalar_checkpoint(
            checkpoint_rows, rc, "u_fast", state0["u_fast"], state_t["u_fast"]
        )
        errors["c_state"] = _scalar_checkpoint(
            checkpoint_rows, rc, "c_state", state0["c_state"], state_t["c_state"]
        )
        errors["eL"] = _vector_checkpoint(
            checkpoint_rows, rc, "eL", cand0["eL"], cand_t["eL"]
        )
        errors["PT"] = _tensor_checkpoint(
            checkpoint_rows, rc, "PT", cand0["P"], cand_t["P"]
        )

        pc = _pair_covariance(cand0, cand_t, rc)
        errors["pair_amplitude_oriented"] = pc["E_pair_amplitude_oriented"]
        errors["pair_response_oriented"] = pc["E_pair_response_oriented"]

        for name in ("pair_amplitude_oriented", "pair_response_oriented"):
            e = errors[name]
            checkpoint_rows.append({
                "transform": rc,
                "checkpoint": name,
                "field_type": "pair_oriented",
                "relative_error": e,
                "classification": _classify(e),
                "passes_1e12": e <= 1e-12,
                "passes_1e8": e <= 1e-8,
            })

        end_t = (
            cand_t["end"]["Rx_3d"],
            cand_t["end"]["Ry_3d"],
            cand_t["end"]["Rz_3d"],
        )
        end0 = (
            cand0["end"]["Rx_3d"],
            cand0["end"]["Ry_3d"],
            cand0["end"]["Rz_3d"],
        )
        errors["endpoint"] = _vector_checkpoint(
            checkpoint_rows, rc, "endpoint", end0, end_t
        )

        iface_t = (
            cand_t["iface"]["Rx_3d_interface"],
            cand_t["iface"]["Ry_3d_interface"],
            cand_t["iface"]["Rz_3d_interface"],
        )
        iface0 = (
            cand0["iface"]["Rx_3d_interface"],
            cand0["iface"]["Ry_3d_interface"],
            cand0["iface"]["Rz_3d_interface"],
        )
        errors["interface"] = _vector_checkpoint(
            checkpoint_rows, rc, "interface", iface0, iface_t
        )

        wrong_scalar = M03.scalar_only_inverse_wrong_control(*end_t, rc)
        wrong_scalar_e = _err_vector(end0, wrong_scalar)
        wrong_rows.extend([
            {
                "transform": rc,
                "control": "scalar_only_inverse_endpoint",
                "relative_error": wrong_scalar_e,
            },
            {
                "transform": rc,
                "control": "pair_amplitude_ignore_endpoint_swap",
                "relative_error": pc["E_pair_amplitude_no_swap_wrong_control"],
            },
            {
                "transform": rc,
                "control": "pair_amplitude_wrong_sign_after_swap",
                "relative_error": pc["E_pair_amplitude_wrong_sign_control"],
            },
            {
                "transform": rc,
                "control": "pair_response_ignore_endpoint_swap",
                "relative_error": pc["E_pair_response_no_swap_wrong_control"],
            },
        ])

        pair_summary_rows.append({
            "transform": rc,
            "mapped_pair_count": pc["mapped_pair_count"],
            "missing_pair_count": pc["missing_pair_count"],
            "endpoint_swap_count": pc["endpoint_swap_count"],
            "E_pair_amplitude_oriented": pc["E_pair_amplitude_oriented"],
            "E_pair_response_oriented": pc["E_pair_response_oriented"],
        })
        direction_rows.extend(pc["direction_rows"])

        order = [
            "rho_3d",
            "u_slow",
            "u_fast",
            "c_state",
            "eL",
            "PT",
            "pair_amplitude_oriented",
            "pair_response_oriented",
            "endpoint",
            "interface",
        ]
        previous = None
        first = None
        for name in order:
            if errors[name] > FIRST_FAILURE_TOL:
                first = name
                break
            previous = name
        first_failure[rc] = {
            "first_failure_checkpoint": first,
            "relative_error": errors[first] if first is not None else None,
            "previous_checkpoint": previous,
            "previous_relative_error": errors[previous] if previous is not None else None,
        }

    _write_csv("checkpoint_covariance.csv", checkpoint_rows)
    _write_csv("pair_slot_covariance.csv", pair_summary_rows)
    _write_csv("pair_direction_covariance.csv", direction_rows)
    _write_csv("wrong_control_covariance.csv", wrong_rows)
    _write_json("first_failure.json", first_failure)

    non_rc0 = [rc for rc in M01.RC_TRANSFORMS if rc != "RC0"]
    firsts = [first_failure[rc]["first_failure_checkpoint"] for rc in non_rc0]
    if any(x in ("u_slow", "u_fast", "c_state") for x in firsts):
        outcome = "Outcome A — A8 EVOLUTION COVARIANCE FAILURE"
    elif any(x in ("eL", "PT") for x in firsts):
        outcome = "Outcome B — LONGITUDINAL/PROJECTOR COVARIANCE FAILURE"
    elif any(x == "pair_amplitude_oriented" for x in firsts):
        outcome = "Outcome D — PAIR-AMPLITUDE COVARIANCE FAILURE"
    elif any(x == "pair_response_oriented" for x in firsts):
        outcome = "Outcome E — PAIR-RESPONSE COVARIANCE FAILURE"
    elif any(x in ("endpoint", "interface") for x in firsts):
        outcome = "Outcome F — ASSEMBLY COVARIANCE FAILURE"
    elif all(x is None for x in firsts):
        outcome = "Outcome G — NO LOCALIZED FAILURE"
    else:
        outcome = "Outcome C — PAIR-ORIENTATION/OTHER LOCALIZED FAILURE"

    validation = {
        "lab_id": LAB_ID,
        "head_sha": repo["head_sha"],
        "outcome": outcome,
        "first_failure_threshold": FIRST_FAILURE_TOL,
        "first_failure": first_failure,
        "full_candidate_rerun_allowed": False,
        "duration_seconds": time.perf_counter() - started,
    }
    _write_json("validation.json", validation)
    _write_json("run.json", {
        "lab_id": LAB_ID,
        "head_sha": repo["head_sha"],
        "cluster_id": CLUSTER_ID,
        "candidate_id": CANDIDATE_ID,
        "nz": NZ,
        "profile": PROFILE,
        "duration_seconds": validation["duration_seconds"],
    })

    report = [
        f"# {LAB_ID}",
        "",
        f"- tested commit: `{repo['head_sha']}`",
        f"- outcome: **{outcome}**",
        f"- threshold: `{FIRST_FAILURE_TOL:.1e}`",
        "",
        "## First failing checkpoint",
        "",
        "| RC | first failure | error | previous | previous error |",
        "|---|---|---:|---|---:|",
    ]
    for rc in M01.RC_TRANSFORMS:
        ff = first_failure[rc]
        report.append(
            f"| {rc} | {ff['first_failure_checkpoint']} | "
            f"{ff['relative_error'] if ff['relative_error'] is not None else '—'} | "
            f"{ff['previous_checkpoint']} | "
            f"{ff['previous_relative_error'] if ff['previous_relative_error'] is not None else '—'} |"
        )
    report += [
        "",
        "This lab is diagnostic only. It does not modify core physics modules and does not authorize a full candidate rerun.",
    ]
    (OUT / "report.md").write_text("\n".join(report) + "\n")

    print(json.dumps(validation, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
