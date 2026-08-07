#!/usr/bin/env python3
"""PBUF FOUNDATION — G3D ANGULAR HIGHER-MOMENT SUFFICIENCY 001.

Diagnostic-only follow-up to the validated angular received-distribution and
angular moment-decomposition labs.

Frozen science:
- candidate PL1_PM1_PS2
- physical source M10_interface_field
- C25 source geometry and photon count
- normalized G3D LOS-consistent propagation
- fixed observer normal n=(0,0,1)
- tangent-plane coordinates tx=vx/vz, ty=vy/vz only after full 3D arrival
- checkpoints 0,1,5,10,20,40,80,120,159

No conventional gravitational law is introduced. Observed kappa remains an
external morphology benchmark only.

Question tested here:

  Is the observer-received angular distribution adequately described, at the
  information-content level, by its validated first two moments (mu,C), or do
  third/fourth central moments retain substantial additional structured content?

This lab does NOT choose an observable and does NOT assume a Gaussian detector
response. It measures higher central angular moments directly from the per-ray
arrival directions. Gaussian reference values are reported only as descriptive
controls and are never execution gates.
"""
from __future__ import annotations

import csv
import json
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))

import pbuf.labs.foundation.m10_coverage_25pct_science001 as BASE
import pbuf.labs.foundation.los_consistent_ray_geometry001 as GEO
import pbuf.labs.foundation.g3d_angular_received_distribution001 as ANG
from pbuf.core import los_projection as M14
from pbuf.core import observable_extraction as M16

LAB_ID = "PBUF-FOUNDATION-G3D-ANGULAR-HIGHER-MOMENT-SUFFICIENCY-001"
OUT = ROOT / "runs" / "g3d_angular_higher_moment_sufficiency001"
CHECKPOINTS = GEO.CHECKPOINTS
EXPECTED_SUPPORT = GEO.EXPECTED_SUPPORT
UNIT_SPEED_TOL = GEO.UNIT_SPEED_TOL
FIRST_STEP_TOL = GEO.FIRST_STEP_TOL
IDENTITY_TOL = 1e-12
PARITY_TOL = 1e-14
ENERGY_EPS = 1e-30
VZ_MIN = 1e-12


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


def _json_default(o):
    if isinstance(o, np.floating): return float(o)
    if isinstance(o, np.integer): return int(o)
    if isinstance(o, np.bool_): return bool(o)
    if isinstance(o, Path): return str(o)
    return str(o)


def _write_json(path: Path, obj) -> None:
    path.write_text(json.dumps(obj, indent=2, default=_json_default))


def _write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        path.write_text("")
        return
    keys: list[str] = []
    for row in rows:
        for key in row:
            if key not in keys:
                keys.append(key)
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)


def _corr(a, b):
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    if a.shape != b.shape:
        raise RuntimeError(f"correlation shape mismatch {a.shape} vs {b.shape}")
    mask = np.isfinite(a) & np.isfinite(b)
    n = int(mask.sum())
    if n < 2:
        return float("nan"), float("nan"), n
    return (
        float(M16.safe_pearson(a[mask], b[mask])),
        float(M16.safe_spearman(a[mask], b[mask])),
        n,
    )


def _rms(a) -> float:
    x = np.asarray(a, dtype=np.float64)
    return float(np.sqrt(np.mean(x*x))) if x.size else float("nan")


def _safe_rel_rms(diff, reference) -> float:
    return ANG._safe_rel_rms(diff, reference)


def _empty_map() -> np.ndarray:
    return np.full((BASE.OBS_BINS, BASE.OBS_BINS), np.nan, dtype=np.float64)


def _higher_moment_fields(snap: dict, groups: dict[int, np.ndarray]) -> dict:
    vx = np.asarray(snap["vx"], dtype=np.float64)
    vy = np.asarray(snap["vy"], dtype=np.float64)
    vz = np.asarray(snap["vz"], dtype=np.float64)
    if float(np.min(np.abs(vz))) <= VZ_MIN:
        raise RuntimeError("observer tangent projection has |vz| <= VZ_MIN")

    tx = vx / vz
    ty = vy / vz

    names = (
        "central_second_trace",
        "central_third_xxx",
        "central_third_xxy",
        "central_third_xyy",
        "central_third_yyy",
        "central_third_tensor_frobenius_mag",
        "standardized_third_tensor_mag",
        "central_fourth_xxxx",
        "central_fourth_xxxy",
        "central_fourth_xxyy",
        "central_fourth_xyyy",
        "central_fourth_yyyy",
        "central_fourth_tensor_frobenius_mag",
        "central_fourth_radial_moment",
        "standardized_fourth_radial_ratio",
        "gaussian_reference_fourth_radial_excess",
        "standardized_fourth_tensor_mag",
        "higher_moment_standardized_combined_mag",
        "photon_count_per_bin",
    )
    out = {name: _empty_map() for name in names}

    for q, idx in groups.items():
        r, c = divmod(q, BASE.OBS_BINS)
        x = tx[idx]
        y = ty[idx]
        dx = x - np.mean(x)
        dy = y - np.mean(y)

        cxx = float(np.mean(dx*dx))
        cyy = float(np.mean(dy*dy))
        s2 = cxx + cyy

        m300 = float(np.mean(dx**3))
        m210 = float(np.mean((dx**2)*dy))
        m120 = float(np.mean(dx*(dy**2)))
        m030 = float(np.mean(dy**3))
        t3f = float(np.sqrt(m300*m300 + 3.0*m210*m210 + 3.0*m120*m120 + m030*m030))

        m400 = float(np.mean(dx**4))
        m310 = float(np.mean((dx**3)*dy))
        m220 = float(np.mean((dx**2)*(dy**2)))
        m130 = float(np.mean(dx*(dy**3)))
        m040 = float(np.mean(dy**4))
        t4f = float(np.sqrt(m400*m400 + 4.0*m310*m310 + 6.0*m220*m220 + 4.0*m130*m130 + m040*m040))
        radial4 = m400 + 2.0*m220 + m040

        if s2 > ENERGY_EPS:
            std3 = t3f / (s2 ** 1.5)
            std4rad = radial4 / (s2*s2)
            std4tensor = t4f / (s2*s2)
            gexcess = std4rad - 2.0
            combined = float(np.sqrt(std3*std3 + gexcess*gexcess))
        else:
            std3 = float("nan")
            std4rad = float("nan")
            std4tensor = float("nan")
            gexcess = float("nan")
            combined = float("nan")

        out["central_second_trace"][r,c] = s2
        out["central_third_xxx"][r,c] = m300
        out["central_third_xxy"][r,c] = m210
        out["central_third_xyy"][r,c] = m120
        out["central_third_yyy"][r,c] = m030
        out["central_third_tensor_frobenius_mag"][r,c] = t3f
        out["standardized_third_tensor_mag"][r,c] = std3
        out["central_fourth_xxxx"][r,c] = m400
        out["central_fourth_xxxy"][r,c] = m310
        out["central_fourth_xxyy"][r,c] = m220
        out["central_fourth_xyyy"][r,c] = m130
        out["central_fourth_yyyy"][r,c] = m040
        out["central_fourth_tensor_frobenius_mag"][r,c] = t4f
        out["central_fourth_radial_moment"][r,c] = radial4
        out["standardized_fourth_radial_ratio"][r,c] = std4rad
        out["gaussian_reference_fourth_radial_excess"][r,c] = gexcess
        out["standardized_fourth_tensor_mag"][r,c] = std4tensor
        out["higher_moment_standardized_combined_mag"][r,c] = combined
        out["photon_count_per_bin"][r,c] = float(len(idx))

    return out


def _gates(base: dict, higher: dict) -> dict:
    second_trace_parity = _safe_rel_rms(
        higher["central_second_trace"] - base["angular_cov_trace"],
        base["angular_cov_trace"],
    )
    radial_from_components = (
        higher["central_fourth_xxxx"]
        + 2.0*higher["central_fourth_xxyy"]
        + higher["central_fourth_yyyy"]
    )
    radial_identity = _safe_rel_rms(
        higher["central_fourth_radial_moment"] - radial_from_components,
        higher["central_fourth_radial_moment"],
    )

    s2 = higher["central_second_trace"]
    nonzero = np.isfinite(s2) & (s2 > ENERGY_EPS)
    zero = np.isfinite(s2) & ~nonzero
    std3 = higher["standardized_third_tensor_mag"]
    std4 = higher["standardized_fourth_radial_ratio"]
    std4t = higher["standardized_fourth_tensor_mag"]

    nonzero_defined_pass = bool(
        np.all(np.isfinite(std3[nonzero]))
        and np.all(np.isfinite(std4[nonzero]))
        and np.all(np.isfinite(std4t[nonzero]))
    )
    zero_undefined_pass = bool(
        np.all(~np.isfinite(std3[zero]))
        and np.all(~np.isfinite(std4[zero]))
        and np.all(~np.isfinite(std4t[zero]))
    ) if np.any(zero) else True

    third_nonnegative = bool(np.all(higher["central_third_tensor_frobenius_mag"][np.isfinite(higher["central_third_tensor_frobenius_mag"])] >= -IDENTITY_TOL))
    fourth_nonnegative = bool(np.all(higher["central_fourth_radial_moment"][np.isfinite(higher["central_fourth_radial_moment"])] >= -IDENTITY_TOL))

    return {
        "second_trace_vs_prior_cov_trace_relative_rms_error": second_trace_parity,
        "fourth_radial_component_identity_relative_rms_error": radial_identity,
        "standardized_higher_moment_defined_count": int(np.sum(nonzero)),
        "zero_spread_standardized_higher_moments_undefined_pass": zero_undefined_pass,
        "nonzero_spread_standardized_higher_moments_defined_pass": nonzero_defined_pass,
        "third_tensor_frobenius_nonnegative_pass": third_nonnegative,
        "fourth_radial_nonnegative_pass": fourth_nonnegative,
    }


def _require_gates(cid: str, k: int, prior: dict, gates: dict) -> None:
    if prior["second_moment_equals_cov_plus_centroid_outer_relative_rms_error"] > ANG.MOMENT_IDENTITY_TOL:
        raise RuntimeError(f"{cid}: inherited angular moment identity failed at step {k}")
    if not prior["covariance_psd_pass"]:
        raise RuntimeError(f"{cid}: inherited covariance PSD failed at step {k}")
    if not prior["direction_mean_vector_bound_pass"]:
        raise RuntimeError(f"{cid}: inherited direction mean-vector bound failed at step {k}")
    if gates["second_trace_vs_prior_cov_trace_relative_rms_error"] > PARITY_TOL:
        raise RuntimeError(f"{cid}: second-trace parity failed at step {k}")
    if gates["fourth_radial_component_identity_relative_rms_error"] > IDENTITY_TOL:
        raise RuntimeError(f"{cid}: fourth radial component identity failed at step {k}")
    if not gates["zero_spread_standardized_higher_moments_undefined_pass"]:
        raise RuntimeError(f"{cid}: zero-spread undefined semantics failed at step {k}")
    if not gates["nonzero_spread_standardized_higher_moments_defined_pass"]:
        raise RuntimeError(f"{cid}: nonzero-spread standardized moments undefined at step {k}")
    if not gates["third_tensor_frobenius_nonnegative_pass"]:
        raise RuntimeError(f"{cid}: third tensor norm negative at step {k}")
    if not gates["fourth_radial_nonnegative_pass"]:
        raise RuntimeError(f"{cid}: fourth radial moment negative at step {k}")


def _add_correlations(row: dict, fields: dict, observed, los_mag) -> None:
    names = (
        "central_second_trace",
        "central_third_tensor_frobenius_mag",
        "standardized_third_tensor_mag",
        "central_fourth_tensor_frobenius_mag",
        "central_fourth_radial_moment",
        "standardized_fourth_radial_ratio",
        "gaussian_reference_fourth_radial_excess",
        "standardized_fourth_tensor_mag",
        "higher_moment_standardized_combined_mag",
        "photon_count_per_bin",
    )
    for name in names:
        arr = fields[name]
        p,s,n = _corr(arr, observed)
        p2,s2,n2 = _corr(arr, los_mag)
        finite = arr[np.isfinite(arr)]
        row[f"{name}_vs_observed_pearson"] = p
        row[f"{name}_vs_observed_spearman"] = s
        row[f"{name}_vs_observed_count"] = n
        row[f"{name}_vs_los_mag_pearson"] = p2
        row[f"{name}_vs_los_mag_spearman"] = s2
        row[f"{name}_vs_los_mag_count"] = n2
        row[f"{name}_rms"] = _rms(finite) if finite.size else float("nan")


def _run_cluster(cluster):
    cid = cluster["id"]
    real = BASE._load_cluster(cluster)
    state = BASE._evolve(BASE._initial_state(real["rho3"]))
    candidate = BASE._candidate(state)
    vector = BASE._interface_vector(candidate)
    los = M14.project_vector_to_image_plane(*vector, los_axis="z")
    Rx = np.asarray(los["comp_1"], dtype=np.float64)
    Ry = np.asarray(los["comp_2"], dtype=np.float64)
    los_mag = np.hypot(Rx,Ry)
    observed = np.asarray(real["observed_kappa"], dtype=np.float64)
    grid = np.linspace(-BASE.CFG["extent"], BASE.CFG["extent"], Rx.shape[0])
    field = {"xgrid":grid,"ygrid":grid,"rx":Rx,"ry":Ry}

    x0,y0,_,_ = BASE._launch_expanded_25pct()
    groups = GEO._source_groups(x0,y0)
    if len(groups) != EXPECTED_SUPPORT:
        raise RuntimeError(f"{cid}: expected {EXPECTED_SUPPORT} source bins, got {len(groups)}")

    checkpoints,g3d = GEO._propagate_g3d(field, BASE.CFG["step"], BASE.CFG["steps"], x0, y0)
    if g3d["max_unit_speed_error"] > UNIT_SPEED_TOL:
        raise RuntimeError(f"{cid}: G3D unit-speed gate failed: {g3d['max_unit_speed_error']}")
    first = GEO._first_step_geometry(field,x0,y0,checkpoints[1],observed,los_mag)
    if not first["first_step_exact_pass"]:
        raise RuntimeError(f"{cid}: first-step exact geometry gate failed")

    rows=[]; all_fields={}
    for k in CHECKPOINTS:
        base = ANG._angular_distribution_fields(checkpoints[k], groups)
        prior = ANG._moment_gates(base)
        higher = _higher_moment_fields(checkpoints[k], groups)
        gates = _gates(base, higher)
        _require_gates(cid, k, prior, gates)
        row = {
            "cluster_id": cid,
            "step_index": int(k),
            "propagation_distance": float(k*BASE.CFG["step"]),
            **gates,
        }
        _add_correlations(row, higher, observed, los_mag)
        rows.append(row)
        all_fields[f"step_{k}"] = higher

    final = rows[-1]
    summary = {
        "cluster_id":cid,
        "candidate_id":BASE.CANDIDATE_ID,
        "physical_source_representation":BASE.PHYSICAL_SOURCE,
        "geometry_lane":"G3D_LOS_consistent_diagnostic",
        "observer_readout_lane":"per_ray_tangent_distribution_higher_central_moment_sufficiency",
        "observer_coordinates":"fixed_tangent_plane_tx_vx_over_vz_ty_vy_over_vz_after_full_3d_arrival_state",
        "observer_information_used":"third_and_fourth_central_moments_beyond_validated_mu_and_covariance",
        "benchmark_role":"external_morphology_comparison_only",
        "gaussian_reference_role":"descriptive_control_only_not_model_or_gate",
        "n_photons":int(len(x0)),
        "source_supported_bins":int(len(groups)),
        "checkpoint_count":len(rows),
        "g3d_unit_speed_max_error":float(g3d["max_unit_speed_error"]),
        "g3d_unit_speed_pass":bool(g3d["max_unit_speed_error"] <= UNIT_SPEED_TOL),
        "first_step_exact_max_vector_error":first["first_step_exact_max_vector_error"],
        "first_step_exact_pass":first["first_step_exact_pass"],
        "los_mag_vs_observed_pearson":_corr(los_mag,observed)[0],
        "los_mag_vs_observed_spearman":_corr(los_mag,observed)[1],
    }
    for key,val in final.items():
        if key not in ("cluster_id","step_index","propagation_distance"):
            summary[f"final_{key}"] = val
    return summary, rows, all_fields, {"los_mag":los_mag,"observed_benchmark":observed}


def _nanmean_key(summaries, key):
    vals=np.asarray([r[key] for r in summaries],dtype=np.float64)
    return float(np.nanmean(vals)) if np.any(np.isfinite(vals)) else float("nan")


def main() -> int:
    started=time.perf_counter()
    OUT.mkdir(parents=True,exist_ok=True)
    repo=_repo_state()
    _write_json(OUT/"repository_state.json",repo)
    if repo["branch"] != "main" or repo["tracked_changes"] or repo["staged_changes"]:
        validation={"lab_id":LAB_ID,"outcome":"REPOSITORY_GATE_FAILURE","head_sha":repo["head_sha"]}
        _write_json(OUT/"validation.json",validation)
        print(json.dumps(validation,indent=2)); return 2

    summaries=[]; checkpoint_rows=[]; failures=[]
    for cluster in BASE.CLUSTERS:
        cid=cluster["id"]
        print(f"[{cid}] higher angular moment sufficiency audit")
        try:
            summary,rows,fields,final_arrays=_run_cluster(cluster)
            summaries.append(summary); checkpoint_rows.extend(rows)
            cdir=OUT/"clusters"/cid; cdir.mkdir(parents=True,exist_ok=True)
            _write_json(cdir/"angular_higher_moment_summary.json",summary)
            _write_csv(cdir/"angular_higher_moment_checkpoints.csv",rows)
            npz={}
            for step_name,fd in fields.items():
                for name,arr in fd.items(): npz[f"{step_name}__{name}"]=arr
            npz.update(final_arrays)
            np.savez_compressed(cdir/"angular_higher_moment_checkpoint_fields.npz",**npz)
        except Exception as exc:
            failures.append({"cluster_id":cid,"error":repr(exc)})
            _write_json(OUT/"cluster_failures.json",failures)
            raise

    _write_csv(OUT/"angular_higher_moment_summary.csv",summaries)
    _write_csv(OUT/"angular_higher_moment_checkpoint_summary.csv",checkpoint_rows)
    _write_json(OUT/"cluster_failures.json",failures)

    validation={
        "lab_id":LAB_ID,
        "outcome":"Outcome A — G3D ANGULAR HIGHER-MOMENT SUFFICIENCY AUDIT COMPLETE",
        "head_sha":repo["head_sha"],
        "candidate_id":BASE.CANDIDATE_ID,
        "physical_source_representation":BASE.PHYSICAL_SOURCE,
        "geometry_lane":"G3D_LOS_consistent_diagnostic",
        "observer_readout_lane":"per_ray_tangent_distribution_higher_central_moment_sufficiency",
        "observer_coordinates":"fixed_tangent_plane_tx_vx_over_vz_ty_vy_over_vz_after_full_3d_arrival_state",
        "observer_information_used":"third_and_fourth_central_moments_beyond_validated_mu_and_covariance",
        "benchmark_role":"external_morphology_comparison_only",
        "gaussian_reference_role":"descriptive_control_only_not_model_or_gate",
        "cluster_count_expected":len(BASE.CLUSTERS),
        "cluster_count_completed":len(summaries),
        "checkpoint_steps":list(CHECKPOINTS),
        "identity_tolerance":IDENTITY_TOL,
        "parity_tolerance":PARITY_TOL,
        "all_cluster_g3d_unit_speed_pass":bool(all(r["g3d_unit_speed_pass"] for r in summaries)),
        "all_cluster_first_step_exact_pass":bool(all(r["first_step_exact_pass"] for r in summaries)),
        "all_cluster_final_nonzero_spread_higher_moments_defined_pass":bool(all(r["final_nonzero_spread_standardized_higher_moments_defined_pass"] for r in summaries)),
        "mean_los_mag_vs_observed_pearson":float(np.mean([r["los_mag_vs_observed_pearson"] for r in summaries])),
        "mean_final_central_third_tensor_frobenius_mag_vs_observed_pearson":_nanmean_key(summaries,"final_central_third_tensor_frobenius_mag_vs_observed_pearson"),
        "mean_final_standardized_third_tensor_mag_vs_observed_pearson":_nanmean_key(summaries,"final_standardized_third_tensor_mag_vs_observed_pearson"),
        "mean_final_central_fourth_radial_moment_vs_observed_pearson":_nanmean_key(summaries,"final_central_fourth_radial_moment_vs_observed_pearson"),
        "mean_final_standardized_fourth_radial_ratio_vs_observed_pearson":_nanmean_key(summaries,"final_standardized_fourth_radial_ratio_vs_observed_pearson"),
        "mean_final_gaussian_reference_fourth_radial_excess_vs_observed_pearson":_nanmean_key(summaries,"final_gaussian_reference_fourth_radial_excess_vs_observed_pearson"),
        "mean_final_higher_moment_standardized_combined_mag_vs_observed_pearson":_nanmean_key(summaries,"final_higher_moment_standardized_combined_mag_vs_observed_pearson"),
        "mean_final_standardized_third_tensor_mag_rms":float(np.nanmean([r["final_standardized_third_tensor_mag_rms"] for r in summaries])),
        "mean_final_gaussian_reference_fourth_radial_excess_rms":float(np.nanmean([r["final_gaussian_reference_fourth_radial_excess_rms"] for r in summaries])),
        "physics_change_authorized":False,
        "candidate_change_authorized":False,
        "production_geometry_change_authorized":False,
        "observer_projection_change_authorized":False,
        "observer_distribution_truncation_authorized":False,
        "observable_selection_authorized":False,
        "higher_moment_interpretation_required":True,
        "native_image_formation_authorized":False,
        "next_experiment_authorized":False,
        "duration_seconds":float(time.perf_counter()-started),
    }
    _write_json(OUT/"validation.json",validation)
    _write_json(OUT/"run.json",{
        "lab_id":LAB_ID,
        "head_sha":repo["head_sha"],
        "output_directory":str(OUT.relative_to(ROOT)),
        "duration_seconds":validation["duration_seconds"],
    })
    print(json.dumps(validation,indent=2,default=_json_default))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
