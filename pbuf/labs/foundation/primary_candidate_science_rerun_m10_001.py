#!/usr/bin/env python3
"""PBUF 3D PAIRWISE PRIMARY-CANDIDATE SCIENCE RE-RUN 001 — M10 INTERFACE SOURCE.

Five-cluster science rerun for the frozen PL1_PM1_PS2 primary candidate after
PBUF-FOUNDATION-INTERFACE-SOURCE-REQUALIFICATION-001 authorized the M10
interface field as the coordinate-safe downstream physical source.

Scientific routing, fixed for this run:

    real cluster kappa
      -> common positive-kappa proxy
      -> Nz=9 Gaussian depth reconstruction
      -> frozen A8/T1 evolution
      -> eL / PT
      -> A_ij
      -> PM1 + PS2 pair response
      -> M10 interface field
      -> native z LOS (M14)
      -> M15 ray input validation
      -> frozen photon propagation
      -> ray-bundle Jacobian
      -> kappa / gamma1 / gamma2
      -> comparison to the observed kappa map

M09 endpoint storage is retained only for conservation/source-sink bookkeeping.
It is NOT used as the physical source for LOS or rays.

No fitting, optimisation, amplitude matching, cluster-specific tuning, viewing
angle selection, or parameter search is performed.
"""
from __future__ import annotations

import csv
import hashlib
import json
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

from pbuf.core import pair_enumeration as M05
from pbuf.core import pair_transfer as M08
from pbuf.core import midpoint_rasterization as M10
from pbuf.core import helmholtz_3d as M13
from pbuf.core import los_projection as M14
from pbuf.core import ray_interface as M15
from pbuf.core import observable_extraction as M16
from pbuf.models import a8_state as M06_state
from pbuf.models import a8_pair_amplitude as M06
from pbuf.models import transverse_projector as M07

LAB_ID = "PBUF-3D-PAIRWISE-PRIMARY-CANDIDATE-SCIENCE-RE-RUN-001-M10"
OUT = ROOT / "runs" / "primary_candidate_science_rerun_m10_001"
BENCHMARK = ROOT / "PBUF_benchmark"

CANDIDATE_ID = "PL1_PM1_PS2"
PHYSICAL_SOURCE = "M10_interface_field"
ENDPOINT_ROLE = "conservation_bookkeeping_only"
NZ = 9
PROFILE = "gaussian"
STRENGTH = 0.18
SEED = 12345
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


def _json_default(obj):
    if isinstance(obj, np.floating):
        return float(obj)
    if isinstance(obj, np.integer):
        return int(obj)
    if isinstance(obj, np.bool_):
        return bool(obj)
    if isinstance(obj, Path):
        return str(obj)
    if isinstance(obj, tuple):
        return list(obj)
    return str(obj)


def _write_json(path: Path, obj) -> None:
    path.write_text(json.dumps(obj, indent=2, default=_json_default))


def _write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        path.write_text("")
        return
    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


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
    if x.size == 0:
        return float("nan")
    return float(np.sqrt(np.mean(x * x)))


def _load_cluster(cluster: dict) -> dict:
    path = BENCHMARK / cluster["directory"] / (
        f"hlsp_frontier_model_{cluster['slug']}_merten_v1_kappa.fits"
    )
    if not path.exists():
        raise FileNotFoundError(path)
    with fits.open(path) as hdul:
        kappa = np.asarray(hdul[0].data, dtype=np.float64)
        header = hdul[0].header
    rho2 = construct_common_proxy(kappa, bins=CFG["bins"], extent=CFG["extent"])
    rho3 = construct_rho_3d(rho2, NZ, profile=PROFILE)
    return {
        "cluster": cluster,
        "path": path,
        "kappa": kappa,
        "rho2": rho2,
        "rho3": rho3,
        "provenance": {
            "cluster_id": cluster["id"],
            "cluster_label": cluster["label"],
            "input_kind": "observed_frontier_fields_fits",
            "fits_path": str(path.relative_to(ROOT)),
            "fits_sha256": _sha_file(path),
            "fits_shape": list(kappa.shape),
            "proxy_shape": list(rho2.shape),
            "rho3d_shape": list(rho3.shape),
            "proxy_sha256": _sha_arr(rho2),
            "rho3d_sha256": _sha_arr(rho3),
            "Z_L": float(header["Z_L"]) if "Z_L" in header else None,
            "Z_S": float(header["Z_S"]) if "Z_S" in header else None,
        },
    }


def _initial_state(rho3: np.ndarray) -> dict:
    rng = np.random.RandomState(SEED)
    eq = STRENGTH * rho3
    u_slow0 = eq.copy()
    noise = M06_state.A8_INIT_INJECTION_NOISE * STRENGTH * rng.randn(*rho3.shape)
    u_fast0 = eq + noise
    return {"rho_3d": rho3.copy(), "u_slow0": u_slow0, "u_fast0": u_fast0}


def _evolve(initial: dict) -> dict:
    u_slow, u_fast, history = M06_state.evolve_a8_transport_3d(
        initial["u_slow0"].copy(),
        initial["u_fast0"].copy(),
        stencil="N6",
        boundary="reflective",
    )
    return {
        "rho_3d": initial["rho_3d"].copy(),
        "u_slow": u_slow,
        "u_fast": u_fast,
        "c_state": history[-1],
    }


def _candidate(state: dict) -> dict:
    shape = tuple(state["c_state"].shape)
    pairs = M05.enumerate_internal_pairs(shape)
    ex, ey, ez, valid, gmag = M07.build_longitudinal_direction(state["c_state"])
    projector = M07.build_transverse_projector(ex, ey, ez)
    amplitudes = M06.compute_a8_pair_amplitudes(
        state["u_slow"], state["u_fast"], state["c_state"], pairs
    )
    responses = M08.build_pair_responses(
        pairs, amplitudes, projector,
        magnitude_formulation="PM1",
        pair_symmetrization="PS2",
    )
    endpoint = M08.assemble_endpoint_field(responses, shape)
    interface = M10.rasterize_interface_field(responses, shape)
    return {
        "shape": shape,
        "pairs": pairs,
        "valid_count": int(np.count_nonzero(valid)),
        "gradient_rms": _rms(gmag),
        "endpoint": endpoint,
        "interface": interface,
    }


def _interface_vector(candidate: dict):
    interface = candidate["interface"]
    return (
        interface["Rx_3d_interface"],
        interface["Ry_3d_interface"],
        interface["Rz_3d_interface"],
    )


def _ray_and_observable(cluster_id: str, vector, real: dict) -> dict:
    los = M14.project_vector_to_image_plane(*vector, los_axis="z")
    Rx, Ry = los["comp_1"], los["comp_2"]
    metadata = {
        "candidate_id": CANDIDATE_ID,
        "cluster_id": cluster_id,
        "transform_id": "RC0",
        "role": "los",
        "physical_source_representation": PHYSICAL_SOURCE,
        "source_artifact_ids": [f"{cluster_id}_interface_field"],
    }
    artifact = M15.prepare_ray_input(Rx, Ry, metadata, require_nontrivial=True)
    if Rx.shape != Ry.shape or Rx.ndim != 2 or Rx.shape[0] != Rx.shape[1]:
        raise RuntimeError(f"{cluster_id}: invalid ray image-plane shape {Rx.shape}/{Ry.shape}")

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

    kappa = np.asarray(observable["kappa"], dtype=np.float64)
    g1 = np.asarray(observable["gamma1"], dtype=np.float64)
    g2 = np.asarray(observable["gamma2"], dtype=np.float64)
    finite = np.isfinite(kappa)
    displacement = np.hypot(photons["x"] - x0, photons["y"] - y0)

    trajectory_hash = hashlib.sha256()
    for key in ("xs", "ys", "x", "y", "conservation"):
        trajectory_hash.update(
            np.ascontiguousarray(np.asarray(photons[key], dtype=np.float64)).tobytes()
        )

    return {
        "los": los,
        "artifact": artifact,
        "photons": photons,
        "observable": observable,
        "observed_kappa": observed_kappa,
        "metrics": {
            "los_rx_rms": _rms(Rx),
            "los_ry_rms": _rms(Ry),
            "ray_classification": artifact.statistics["ray_classification"],
            "n_photons": int(len(x0)),
            "mean_endpoint_displacement": float(np.mean(displacement)),
            "max_endpoint_displacement": float(np.max(displacement)),
            "conservation_max": float(np.max(photons["conservation"])),
            "kappa_finite_count": int(finite.sum()),
            "kappa_total_count": int(kappa.size),
            "kappa_finite_fraction": float(finite.mean()),
            "kappa_variance_finite": float(np.var(kappa[finite])) if finite.sum() >= 2 else float("nan"),
            "kappa_rms_finite": _rms(kappa[finite]) if finite.any() else float("nan"),
            "gamma1_rms_finite": _rms(g1[np.isfinite(g1)]) if np.isfinite(g1).any() else float("nan"),
            "gamma2_rms_finite": _rms(g2[np.isfinite(g2)]) if np.isfinite(g2).any() else float("nan"),
            "pearson_vs_observed": observable.get("pearson_vs_reference", float("nan")),
            "spearman_vs_observed": observable.get("spearman_vs_reference", float("nan")),
            "trajectory_sha256": trajectory_hash.hexdigest(),
        },
    }


def _run_cluster(cluster: dict) -> dict:
    cluster_id = cluster["id"]
    real = _load_cluster(cluster)
    initial = _initial_state(real["rho3"])
    state = _evolve(initial)
    candidate = _candidate(state)
    vector = _interface_vector(candidate)

    expected_pairs = int(M08.expected_interface_pair_count(candidate["shape"]))
    consumed_pairs = int(candidate["interface"]["statistics"]["consumed_pair_count_total"])
    interface_energy = float(candidate["interface"]["statistics"]["interface_energy"])
    endpoint_closure = float(candidate["endpoint"]["statistics"]["global_vector_sum_norm"])
    endpoint_energy = float(candidate["endpoint"]["statistics"]["endpoint_energy"])

    ray = _ray_and_observable(cluster_id, vector, real)
    hn = M13.helmholtz_decompose_3d(*vector, padding="none")
    hp = M13.helmholtz_decompose_3d(*vector, padding="reflect_half")
    hkeys = (
        "field_reconstruction_error", "energy_closure_error", "orthogonality_error",
        "f_irr_partition", "f_sol_partition", "f_irr_native", "f_sol_native",
    )

    metrics = {
        "cluster_id": cluster_id,
        "cluster_label": cluster["label"],
        "candidate_id": CANDIDATE_ID,
        "physical_source_representation": PHYSICAL_SOURCE,
        "endpoint_role": ENDPOINT_ROLE,
        "shape": list(candidate["shape"]),
        "n_pairs": int(len(candidate["pairs"])),
        "interface_expected_pair_count": expected_pairs,
        "interface_consumed_pair_count": consumed_pairs,
        "interface_pair_count_ok": bool(expected_pairs == consumed_pairs == len(candidate["pairs"])),
        "interface_energy": interface_energy,
        "endpoint_energy_bookkeeping": endpoint_energy,
        "endpoint_closure_bookkeeping": endpoint_closure,
        "valid_longitudinal_count": int(candidate["valid_count"]),
        "gradient_rms": float(candidate["gradient_rms"]),
        **ray["metrics"],
        "helmholtz_none_field_reconstruction_error": hn["field_reconstruction_error"],
        "helmholtz_none_energy_closure_error": hn["energy_closure_error"],
        "helmholtz_none_orthogonality_error": hn["orthogonality_error"],
        "helmholtz_none_f_irr_partition": hn["f_irr_partition"],
        "helmholtz_none_f_sol_partition": hn["f_sol_partition"],
        "helmholtz_none_f_irr_native": hn["f_irr_native"],
        "helmholtz_none_f_sol_native": hn["f_sol_native"],
        "helmholtz_padded_field_reconstruction_error": hp["field_reconstruction_error"],
        "helmholtz_padded_energy_closure_error": hp["energy_closure_error"],
        "helmholtz_padded_orthogonality_error": hp["orthogonality_error"],
        "helmholtz_padded_f_irr_partition": hp["f_irr_partition"],
        "helmholtz_padded_f_sol_partition": hp["f_sol_partition"],
        "helmholtz_padded_f_irr_native": hp["f_irr_native"],
        "helmholtz_padded_f_sol_native": hp["f_sol_native"],
    }

    kappa = np.asarray(ray["observable"]["kappa"], dtype=np.float64)
    integrity_pass = bool(
        real["provenance"]["input_kind"] == "observed_frontier_fields_fits"
        and interface_energy > 0.0
        and metrics["interface_pair_count_ok"]
        and ray["artifact"].statistics["ray_classification"]
            in ("structured_small", "structured_normal", "constant_nonzero")
        and metrics["kappa_finite_count"] >= 2
        and np.isfinite(metrics["kappa_variance_finite"])
        and metrics["kappa_variance_finite"] > 0.0
        and metrics["mean_endpoint_displacement"] > 0.0
        and np.isfinite(metrics["conservation_max"])
    )
    metrics["integrity_pass"] = integrity_pass

    lineage = {
        "fits": real["provenance"]["fits_sha256"],
        "rho2": real["provenance"]["proxy_sha256"],
        "rho3": real["provenance"]["rho3d_sha256"],
        "u_slow": _sha_arr(state["u_slow"]),
        "u_fast": _sha_arr(state["u_fast"]),
        "c_state": _sha_arr(state["c_state"]),
        "interface_3d": _vec_sha(vector),
        "los": _vec_sha((ray["los"]["comp_1"], ray["los"]["comp_2"], np.zeros_like(ray["los"]["comp_1"]))),
        "ray_input": ray["artifact"].sha256,
        "trajectory": ray["metrics"]["trajectory_sha256"],
        "kappa": _sha_arr(np.nan_to_num(kappa, nan=0.0)),
    }

    cluster_dir = OUT / "clusters" / cluster_id
    cluster_dir.mkdir(parents=True, exist_ok=True)
    _write_json(cluster_dir / "input_provenance.json", real["provenance"])
    _write_json(cluster_dir / "statistics.json", metrics)
    _write_json(cluster_dir / "field_lineage.json", lineage)
    np.savez_compressed(
        cluster_dir / "observables.npz",
        kappa=kappa,
        gamma1=np.asarray(ray["observable"]["gamma1"], dtype=np.float64),
        gamma2=np.asarray(ray["observable"]["gamma2"], dtype=np.float64),
        observed_kappa=np.asarray(ray["observed_kappa"], dtype=np.float64),
        los_comp_1=np.asarray(ray["los"]["comp_1"], dtype=np.float64),
        los_comp_2=np.asarray(ray["los"]["comp_2"], dtype=np.float64),
    )
    return {"metrics": metrics, "lineage": lineage}


def main() -> int:
    started = time.perf_counter()
    OUT.mkdir(parents=True, exist_ok=True)
    repo = _repo_state()
    _write_json(OUT / "repository_state.json", repo)

    if repo["branch"] != "main" or repo["tracked_changes"] or repo["staged_changes"]:
        validation = {
            "lab_id": LAB_ID,
            "outcome": "REPOSITORY_GATE_FAILURE",
            "head_sha": repo["head_sha"],
            "science_run_complete": False,
        }
        _write_json(OUT / "validation.json", validation)
        print(json.dumps(validation, indent=2))
        return 2

    _write_json(OUT / "run_config.json", {
        "lab_id": LAB_ID,
        "head_sha": repo["head_sha"],
        "candidate_id": CANDIDATE_ID,
        "physical_source_representation": PHYSICAL_SOURCE,
        "endpoint_role": ENDPOINT_ROLE,
        "clusters": [c["id"] for c in CLUSTERS],
        "nz": NZ,
        "profile": PROFILE,
        "strength": STRENGTH,
        "seed": SEED,
        "config": CFG,
        "fitting": False,
        "optimisation": False,
        "cluster_specific_tuning": False,
        "amplitude_matching": False,
        "viewing_angle_search": False,
    })

    rows: list[dict] = []
    failures: list[dict] = []
    for cluster in CLUSTERS:
        cid = cluster["id"]
        print(f"[{cid}] M10 interface -> z LOS -> photons -> Jacobian -> observables")
        try:
            result = _run_cluster(cluster)
            rows.append(result["metrics"])
        except Exception as exc:
            failures.append({
                "cluster_id": cid,
                "exception_type": type(exc).__name__,
                "exception": str(exc),
            })
            raise

    _write_csv(OUT / "cluster_science_summary.csv", rows)
    _write_json(OUT / "cluster_failures.json", failures)

    integrity_all = bool(len(rows) == len(CLUSTERS) and all(bool(r["integrity_pass"]) for r in rows))
    finite_corr = [float(r["pearson_vs_observed"]) for r in rows if np.isfinite(r["pearson_vs_observed"])]
    finite_spear = [float(r["spearman_vs_observed"]) for r in rows if np.isfinite(r["spearman_vs_observed"])]

    validation = {
        "lab_id": LAB_ID,
        "outcome": (
            "Outcome A — FIVE-CLUSTER M10 PRIMARY-CANDIDATE SCIENCE RUN COMPLETE"
            if integrity_all else
            "Outcome D — FIVE-CLUSTER M10 PRIMARY-CANDIDATE INTEGRITY FAILURE"
        ),
        "head_sha": repo["head_sha"],
        "candidate_id": CANDIDATE_ID,
        "physical_source_representation": PHYSICAL_SOURCE,
        "endpoint_role": ENDPOINT_ROLE,
        "cluster_count_expected": len(CLUSTERS),
        "cluster_count_completed": len(rows),
        "all_cluster_integrity_pass": integrity_all,
        "pearson_mean_finite": float(np.mean(finite_corr)) if finite_corr else float("nan"),
        "pearson_median_finite": float(np.median(finite_corr)) if finite_corr else float("nan"),
        "spearman_mean_finite": float(np.mean(finite_spear)) if finite_spear else float("nan"),
        "spearman_median_finite": float(np.median(finite_spear)) if finite_spear else float("nan"),
        "science_interpretation_required": True,
        "duration_seconds": time.perf_counter() - started,
    }
    _write_json(OUT / "validation.json", validation)

    report = [
        f"# {LAB_ID}",
        "",
        f"- tested commit: `{repo['head_sha']}`",
        f"- outcome: **{validation['outcome']}**",
        f"- physical source: `{PHYSICAL_SOURCE}`",
        f"- endpoint role: `{ENDPOINT_ROLE}`",
        "",
        "| cluster | integrity | r(kappa,obs) | rho_s | kappa RMS | kappa var | LOS Rx RMS | LOS Ry RMS | f_irr | f_sol | interface energy |",
        "|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in rows:
        report.append(
            f"| {row['cluster_id']} | {row['integrity_pass']} | "
            f"{row['pearson_vs_observed']:.6g} | {row['spearman_vs_observed']:.6g} | "
            f"{row['kappa_rms_finite']:.6g} | {row['kappa_variance_finite']:.6g} | "
            f"{row['los_rx_rms']:.6g} | {row['los_ry_rms']:.6g} | "
            f"{row['helmholtz_none_f_irr_partition']:.6g} | {row['helmholtz_none_f_sol_partition']:.6g} | "
            f"{row['interface_energy']:.6g} |"
        )
    report += [
        "",
        "This run is a science measurement, not a fit. Correlations are reported exactly as obtained and are not acceptance gates.",
        "Earlier endpoint-source primary-candidate results are superseded as physical-source results but remain historical controls.",
    ]
    (OUT / "report.md").write_text("\n".join(report) + "\n")

    print(json.dumps(validation, indent=2))
    return 0 if integrity_all else 1


if __name__ == "__main__":
    raise SystemExit(main())
