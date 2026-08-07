#!/usr/bin/env python3
"""PBUF FOUNDATION — INTERFACE-SOURCE REQUALIFICATION 001.

Requalifies the real MACS0416 PL1_PM1_PS2 chain after the endpoint/interface
role audit established that M09 endpoint storage is orientation-dependent
bookkeeping while M10 interface storage is the coordinate-safe physical
pair-field representation.

This lab does not modify M08/M09/M10. It changes only the downstream source
used by the requalification experiment:

    M10 interface -> M14 native-z LOS -> M15 -> photons -> Jacobian -> kappa/gamma

R3 uses one fixed stochastic A8 initial realization, transforms that physical
initial state under RC0..RC6, reruns A8/T1 and candidate construction, then
checks M10 interface covariance and native-z LOS covariance.
"""
from __future__ import annotations

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
from weak_lensing_observation001 import propagate as wl_propagate, resample_to_grid
import observable_lab001 as obs_lab
import source_plane_lab001 as src_lab

from pbuf.core import coordinate_transforms as M02
from pbuf.core import vector_transforms as M03
from pbuf.core import pair_enumeration as M05
from pbuf.core import pair_transfer as M08
from pbuf.core import midpoint_rasterization as M10
from pbuf.core import helmholtz_3d as M13
from pbuf.core import los_projection as M14
from pbuf.core import ray_interface as M15
from pbuf.core import observable_extraction as M16
from pbuf.core import conventions as M01
from pbuf.models import a8_state as M06_state
from pbuf.models import a8_pair_amplitude as M06
from pbuf.models import transverse_projector as M07

LAB_ID = "PBUF-FOUNDATION-INTERFACE-SOURCE-REQUALIFICATION-001"
OUT = ROOT / "runs" / "interface_source_requalification001"
BENCHMARK = ROOT / "PBUF_benchmark"

CLUSTER_ID = "MACS0416"
CANDIDATE_ID = "PL1_PM1_PS2"
NZ = 9
PROFILE = "gaussian"
STRENGTH = 0.18
SEED = 12345
COV_TOL = 1e-8
CFG = dict(PRODUCTION)


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


def _write_json(name: str, obj) -> None:
    (OUT / name).write_text(json.dumps(obj, indent=2, default=_json_default))


def _write_csv(name: str, rows: list[dict]) -> None:
    import csv
    p = OUT / name
    if not rows:
        p.write_text("")
        return
    fields = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    with p.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)


def _json_default(o):
    if isinstance(o, np.floating):
        return float(o)
    if isinstance(o, np.integer):
        return int(o)
    if isinstance(o, np.bool_):
        return bool(o)
    if isinstance(o, Path):
        return str(o)
    if isinstance(o, tuple):
        return list(o)
    return str(o)


def _sha_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _sha_arr(a) -> str:
    return hashlib.sha256(
        np.ascontiguousarray(np.asarray(a, dtype=np.float64)).tobytes()
    ).hexdigest()


def _vec_sha(v) -> str:
    h = hashlib.sha256()
    for a in v:
        h.update(np.ascontiguousarray(np.asarray(a, dtype=np.float64)).tobytes())
    return h.hexdigest()


def _rms(a) -> float:
    x = np.asarray(a, dtype=np.float64)
    return float(np.sqrt(np.mean(x * x)))


def _err_vector(ref, test) -> float:
    r0, r1, r2 = (np.asarray(x, dtype=np.float64) for x in ref)
    t0, t1, t2 = (np.asarray(x, dtype=np.float64) for x in test)
    num = float(np.sqrt(np.sum((t0-r0)**2 + (t1-r1)**2 + (t2-r2)**2)))
    den = float(np.sqrt(np.sum(r0*r0 + r1*r1 + r2*r2)))
    return num / max(den, 1e-15)


def _err_plane(ref, test) -> float:
    a0, a1 = (np.asarray(x, dtype=np.float64) for x in ref)
    b0, b1 = (np.asarray(x, dtype=np.float64) for x in test)
    num = float(np.sqrt(np.sum((b0-a0)**2 + (b1-a1)**2)))
    den = float(np.sqrt(np.sum(a0*a0 + a1*a1)))
    return num / max(den, 1e-15)


def _load_real_input() -> dict:
    cluster = next(c for c in CLUSTERS if c["id"] == CLUSTER_ID)
    path = BENCHMARK / cluster["directory"] / (
        f"hlsp_frontier_model_{cluster['slug']}_merten_v1_kappa.fits"
    )
    if not path.exists():
        raise FileNotFoundError(path)
    with fits.open(path) as hdul:
        kappa = np.asarray(hdul[0].data, dtype=np.float64)
        hdr = hdul[0].header
    rho2 = construct_common_proxy(kappa, bins=CFG["bins"], extent=CFG["extent"])
    rho3 = construct_rho_3d(rho2, NZ, profile=PROFILE)
    return {
        "cluster": cluster,
        "path": path,
        "kappa": kappa,
        "rho2": rho2,
        "rho3": rho3,
        "provenance": {
            "input_kind": "observed_frontier_fields_fits",
            "cluster_id": CLUSTER_ID,
            "fits_path": str(path.relative_to(ROOT)),
            "fits_sha256": _sha_file(path),
            "fits_shape": list(kappa.shape),
            "proxy_shape": list(rho2.shape),
            "rho3d_shape": list(rho3.shape),
            "proxy_sha256": _sha_arr(rho2),
            "rho3d_sha256": _sha_arr(rho3),
            "Z_L": float(hdr["Z_L"]) if "Z_L" in hdr else None,
            "Z_S": float(hdr["Z_S"]) if "Z_S" in hdr else None,
            "nz": NZ,
            "depth_profile": PROFILE,
            "strength": STRENGTH,
            "seed": SEED,
        },
    }


def _native_initial_state(rho3: np.ndarray) -> dict:
    rng = np.random.RandomState(SEED)
    eq = STRENGTH * rho3
    u_slow0 = eq.copy()
    noise = M06_state.A8_INIT_INJECTION_NOISE * STRENGTH * rng.randn(*rho3.shape)
    u_fast0 = eq + noise
    return {
        "rho_3d": rho3.copy(),
        "u_slow0": u_slow0,
        "u_fast0": u_fast0,
    }


def _transform_initial(initial: dict, rc: str) -> dict:
    return {
        "rho_3d": M02.transform_scalar_field(initial["rho_3d"], rc),
        "u_slow0": M02.transform_scalar_field(initial["u_slow0"], rc),
        "u_fast0": M02.transform_scalar_field(initial["u_fast0"], rc),
    }


def _evolve(initial: dict) -> dict:
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
    ex, ey, ez, valid, gmag = M07.build_longitudinal_direction(state["c_state"])
    projector = M07.build_transverse_projector(ex, ey, ez)
    amp = M06.compute_a8_pair_amplitudes(
        state["u_slow"], state["u_fast"], state["c_state"], pairs
    )
    response = M08.build_pair_responses(pairs, amp, projector, "PM1", "PS2")
    endpoint = M08.assemble_endpoint_field(response, shape)
    interface = M10.rasterize_interface_field(response, shape)
    return {
        "shape": shape,
        "pairs": pairs,
        "valid_count": int(np.count_nonzero(valid)),
        "gradient_rms": _rms(gmag),
        "amp": amp,
        "response": response,
        "endpoint": endpoint,
        "interface": interface,
    }


def _interface_vector(candidate: dict):
    x = candidate["interface"]
    return (
        x["Rx_3d_interface"],
        x["Ry_3d_interface"],
        x["Rz_3d_interface"],
    )


def _real_ray(Rx: np.ndarray, Ry: np.ndarray, real: dict) -> dict:
    meta = {
        "candidate_id": CANDIDATE_ID,
        "cluster_id": CLUSTER_ID,
        "transform_id": "RC0",
        "role": "los",
        "physical_source_representation": "M10_interface_field",
        "source_artifact_ids": ["real_interface_field"],
    }
    artifact = M15.prepare_ray_input(Rx, Ry, meta, require_nontrivial=True)
    if Rx.shape != Ry.shape or Rx.ndim != 2 or Rx.shape[0] != Rx.shape[1]:
        raise RuntimeError("invalid ray image-plane shape")

    n = Rx.shape[0]
    grid = np.linspace(-CFG["extent"], CFG["extent"], n)
    field = {"xgrid": grid, "ygrid": grid, "rx": Rx, "ry": Ry}
    x0, y0, vx0, vy0 = src_lab.launch_B_cartesian(CFG["nphotons"])
    photons = wl_propagate(field, CFG["step"], CFG["steps"], x0, y0, vx0, vy0)
    photons["x0"] = x0
    photons["y0"] = y0

    jac = obs_lab.method_jacobian(
        x0, y0, photons["x"], photons["y"], CFG["extent"], CFG["bins"]
    )
    observed_kappa = resample_to_grid(real["kappa"], CFG["bins"], CFG["extent"])
    observable = M16.package_lensing_observables(
        jac["convergence"], jac["shear_g1"], jac["shear_g2"],
        reference_kappa=observed_kappa,
    )

    kappa = np.asarray(observable["kappa"])
    finite = np.isfinite(kappa)
    displacement = np.hypot(photons["x"] - x0, photons["y"] - y0)
    trajectory_hash = hashlib.sha256()
    for key in ("xs", "ys", "x", "y", "conservation"):
        trajectory_hash.update(
            np.ascontiguousarray(np.asarray(photons[key], dtype=np.float64)).tobytes()
        )

    metrics = {
        "n_photons": int(len(x0)),
        "trajectory_sha256": trajectory_hash.hexdigest(),
        "mean_endpoint_displacement": float(np.mean(displacement)),
        "max_endpoint_displacement": float(np.max(displacement)),
        "conservation_max": float(np.max(photons["conservation"])),
        "kappa_finite_count": int(finite.sum()),
        "kappa_total_count": int(kappa.size),
        "kappa_finite_fraction": float(finite.mean()),
        "kappa_variance_finite": float(np.var(kappa[finite])) if finite.sum() >= 2 else float("nan"),
        "kappa_rms_finite": _rms(kappa[finite]) if finite.any() else float("nan"),
        "pearson_vs_observed": observable.get("pearson_vs_reference", float("nan")),
        "spearman_vs_observed": observable.get("spearman_vs_reference", float("nan")),
    }
    return {"artifact": artifact, "observable": observable, "metrics": metrics}


def _run_R2(real: dict, initial: dict) -> dict:
    state = _evolve(initial)
    candidate = _candidate(state)
    vector = _interface_vector(candidate)
    los = M14.project_vector_to_image_plane(*vector, los_axis="z")
    Rx, Ry = los["comp_1"], los["comp_2"]
    ray = _real_ray(Rx, Ry, real)

    helmholtz_none = M13.helmholtz_decompose_3d(*vector, padding="none")
    helmholtz_padded = M13.helmholtz_decompose_3d(*vector, padding="reflect_half")

    interface_energy = float(candidate["interface"]["statistics"]["interface_energy"])
    endpoint_closure = float(candidate["endpoint"]["statistics"]["global_vector_sum_norm"])
    endpoint_energy = float(candidate["endpoint"]["statistics"]["endpoint_energy"])
    consumed = int(candidate["interface"]["statistics"]["consumed_pair_count_total"])
    expected = int(M08.expected_interface_pair_count(candidate["shape"]))
    km = ray["metrics"]

    passes = bool(
        real["provenance"]["input_kind"] == "observed_frontier_fields_fits"
        and interface_energy > 0.0
        and consumed == expected == len(candidate["pairs"])
        and ray["artifact"].statistics["ray_classification"]
            in ("structured_small", "structured_normal", "constant_nonzero")
        and km["kappa_finite_count"] >= 2
        and np.isfinite(km["kappa_variance_finite"])
        and km["kappa_variance_finite"] > 0.0
        and km["mean_endpoint_displacement"] > 0.0
        and np.isfinite(km["conservation_max"])
    )

    hkeys = (
        "field_reconstruction_error", "energy_closure_error", "orthogonality_error",
        "f_irr_partition", "f_sol_partition", "f_irr_native", "f_sol_native",
    )
    metrics = {
        "cluster_id": CLUSTER_ID,
        "candidate_id": CANDIDATE_ID,
        "transform_id": "RC0",
        "physical_source_representation": "M10_interface_field",
        "endpoint_role": "conservation_bookkeeping",
        "shape": list(candidate["shape"]),
        "n_pairs": int(len(candidate["pairs"])),
        "interface_energy": interface_energy,
        "interface_consumed_pair_count": consumed,
        "interface_expected_pair_count": expected,
        "endpoint_energy_bookkeeping": endpoint_energy,
        "endpoint_closure_bookkeeping": endpoint_closure,
        "los_rx_rms": _rms(Rx),
        "los_ry_rms": _rms(Ry),
        "ray_classification": ray["artifact"].statistics["ray_classification"],
        **km,
        "helmholtz_none": {k: helmholtz_none[k] for k in hkeys},
        "helmholtz_padded": {k: helmholtz_padded[k] for k in hkeys},
    }
    lineage = {
        "fits": real["provenance"]["fits_sha256"],
        "rho2": real["provenance"]["proxy_sha256"],
        "rho3": real["provenance"]["rho3d_sha256"],
        "u_slow": _sha_arr(state["u_slow"]),
        "u_fast": _sha_arr(state["u_fast"]),
        "c_state": _sha_arr(state["c_state"]),
        "interface_3d": _vec_sha(vector),
        "los": _vec_sha((Rx, Ry, np.zeros_like(Rx))),
        "ray_input": ray["artifact"].sha256,
        "trajectory": km["trajectory_sha256"],
        "kappa": _sha_arr(np.nan_to_num(ray["observable"]["kappa"], nan=0.0)),
    }
    return {
        "passes": passes,
        "state": state,
        "candidate": candidate,
        "metrics": metrics,
        "lineage": lineage,
    }


def _run_R3(initial: dict, r2: dict) -> dict:
    ref_vector = _interface_vector(r2["candidate"])
    ref_los_obj = M14.project_vector_to_image_plane(*ref_vector, los_axis="z")
    ref_los = (ref_los_obj["comp_1"], ref_los_obj["comp_2"])
    rows = []

    for rc in M01.RC_TRANSFORMS:
        state_t = _evolve(_transform_initial(initial, rc))
        candidate_t = _candidate(state_t)
        vector_t = _interface_vector(candidate_t)
        back = M03.inverse_transform_vector_field(*vector_t, rc)
        E_interface = _err_vector(ref_vector, back)

        los_back_obj = M14.project_vector_to_image_plane(*back, los_axis="z")
        los_back = (los_back_obj["comp_1"], los_back_obj["comp_2"])
        E_los = _err_plane(ref_los, los_back)

        wrong = M03.scalar_only_inverse_wrong_control(*vector_t, rc)
        E_wrong = _err_vector(ref_vector, wrong)

        endpoint = candidate_t["endpoint"]
        rows.append({
            "transform": rc,
            "shape": list(candidate_t["shape"]),
            "pair_count": int(len(candidate_t["pairs"])),
            "E_cov_interface": E_interface,
            "E_cov_native_z_los": E_los,
            "E_wrong_scalar_only_interface": E_wrong,
            "interface_energy": float(candidate_t["interface"]["statistics"]["interface_energy"]),
            "endpoint_closure_bookkeeping": float(endpoint["statistics"]["global_vector_sum_norm"]),
            "passes": bool(E_interface <= COV_TOL and E_los <= COV_TOL),
        })

    return {
        "rows": rows,
        "passes": all(row["passes"] for row in rows),
        "audit_kind": "fixed-realization-full-A8-interface-field-and-native-z-LOS-covariance",
        "covariance_tolerance": COV_TOL,
    }


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
            "full_candidate_rerun_allowed": False,
        }
        _write_json("validation.json", validation)
        print(json.dumps(validation, indent=2))
        return 2

    real = _load_real_input()
    initial = _native_initial_state(real["rho3"])
    _write_json("input_provenance.json", real["provenance"])

    print("[R2] M10 interface -> native-z LOS -> real photons -> Jacobian")
    r2 = _run_R2(real, initial)
    _write_json("field_lineage.json", r2["lineage"])
    _write_json("real_interface_ray_statistics.json", r2["metrics"])

    if r2["passes"]:
        print("[R3] fixed physical initial realization -> RC0..RC6 -> M10 interface covariance")
        r3 = _run_R3(initial, r2)
    else:
        r3 = {
            "rows": [],
            "passes": False,
            "audit_kind": "not_run_due_to_R2_failure",
            "covariance_tolerance": COV_TOL,
        }
    _write_csv("interface_covariance_revalidation.csv", r3["rows"])

    if r2["passes"] and r3["passes"]:
        outcome = "Outcome A — INTERFACE-SOURCE CORE REQUALIFIED"
        allowed = True
    elif r2["passes"]:
        outcome = "Outcome E — INTERFACE COVARIANCE FAILURE"
        allowed = False
    else:
        outcome = "Outcome D — REAL INTERFACE-RAY RECOVERY FAILURE"
        allowed = False

    validation = {
        "lab_id": LAB_ID,
        "outcome": outcome,
        "head_sha": repo["head_sha"],
        "physical_source_representation": "M10_interface_field",
        "endpoint_role": "conservation_bookkeeping_only",
        "R2_real_cluster_interface_ray_pass": bool(r2["passes"]),
        "R3_interface_and_native_z_los_covariance_pass": bool(r3["passes"]),
        "covariance_tolerance": COV_TOL,
        "second_review_status": "accepted" if allowed else "blocked",
        "full_candidate_rerun_allowed": allowed,
        "next_permitted_experiment": (
            "PBUF 3D PAIRWISE PRIMARY-CANDIDATE SCIENCE RE-RUN 001 — M10 INTERFACE SOURCE"
            if allowed else None
        ),
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
        "strength": STRENGTH,
        "seed": SEED,
        "config": CFG,
        "duration_seconds": validation["duration_seconds"],
    })
    print(json.dumps(validation, indent=2))
    return 0 if allowed else 1


if __name__ == "__main__":
    raise SystemExit(main())
