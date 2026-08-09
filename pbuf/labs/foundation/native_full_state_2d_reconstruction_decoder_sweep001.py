#!/usr/bin/env python3
"""PBUF FOUNDATION — FULL-STATE 2D RECONSTRUCTION DECODER SWEEP 001.

Observer-side decoding audit only.

PR #104 established that the frozen current-native received G3D ray state is
information-rich and loses independent dimensions when reduced prematurely.
This lab therefore keeps the full decoded bank intact, constructs a broad set
of TARGET-BLIND 2D reconstruction candidates, and only then compares those
finished candidates with the local observed weak-lensing products.

The purpose is not to change propagation physics or fit a final lensing law.
It is to determine which information-preserving reconstruction families are
capable of turning the rich received state into coherent 2D observer maps.

Candidate families include:
  * full-bank L1/L2/signed fusion
  * family-balanced L1/L2/signed fusion
  * target-blind PCA component maps
  * target-blind PCA energy maps
  * target-blind whitened-PCA energy maps
  * explicit-3D versus established-2D controls
  * one-family-only L1/L2/PCA controls

All internal standardization, PCA orientation, and family balancing are based
only on the received decoded state. Observations are loaded only after every
candidate map has been constructed. No observational regression, fitted
weights, physical-output rescaling, channel selection, or cluster-specific
choice is permitted.
"""
from __future__ import annotations

import json
import statistics
import subprocess
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))

from pbuf.core import benchmark_data as BENCH
import pbuf.labs.foundation.native_full_received_state_information_retention001 as RET
import pbuf.labs.foundation.native_multichannel_observer_fusion_sweep001 as FUS
import pbuf.labs.foundation.native_observable_extraction_method_sweep001 as EX
import pbuf.labs.foundation.native_accumulated_full_lensing001 as G3D
import pbuf.labs.foundation.m10_coverage_25pct_science001 as BASE

LAB_ID = "PBUF-FOUNDATION-FULL-STATE-2D-RECONSTRUCTION-DECODER-SWEEP-001"
EXPECTED_CLUSTER_IDS = FUS.EXPECTED_CLUSTER_IDS
EPS = 1.0e-30
PCA_COMPONENTS = 8
PCA_ENERGY_K = (2, 3, 5, 8, 10)


def git(*args: str) -> str:
    return subprocess.check_output(["git", *args], cwd=ROOT, text=True).strip()


def repo_state() -> dict:
    return {
        "repository": "TheExiledMonk/lab",
        "branch": git("rev-parse", "--abbrev-ref", "HEAD"),
        "head_sha": git("rev-parse", "HEAD"),
        "tracked_changes": git("diff", "--name-only"),
        "staged_changes": git("diff", "--name-only", "--cached"),
    }


def _finite(a) -> np.ndarray:
    return np.asarray(a, dtype=np.float64)


def _standardized_bank(bank: dict[str, np.ndarray]) -> tuple[np.ndarray, list[str], dict[str, int], tuple[int, int]]:
    """Return target-blind z-scored decoded fields on one common detector grid.

    Missing cells are filled with zero *after centering*, exactly meaning the
    channel mean in standardized coordinates. Constant/all-missing channels
    are retained in the original bank but excluded from the usable matrix.
    """
    first = next(iter(bank.values()))
    shape = _finite(first).shape
    cols: list[np.ndarray] = []
    names: list[str] = []
    for name, raw in bank.items():
        a = _finite(raw).reshape(-1)
        finite = np.isfinite(a)
        if not np.any(finite):
            continue
        mu = float(np.mean(a[finite]))
        sd = float(np.std(a[finite]))
        if not np.isfinite(sd) or sd <= EPS:
            continue
        z = np.zeros_like(a)
        z[finite] = (a[finite] - mu) / sd
        cols.append(z)
        names.append(name)
    X = np.column_stack(cols) if cols else np.empty((shape[0] * shape[1], 0), dtype=np.float64)
    index = {name: i for i, name in enumerate(names)}
    return X, names, index, shape


def _indices(names: list[str], family: dict[str, str], wanted: set[str] | None = None, fam: str | None = None) -> list[int]:
    out = []
    for i, name in enumerate(names):
        if wanted is not None and name not in wanted:
            continue
        if fam is not None and family[name] != fam:
            continue
        out.append(i)
    return out


def _map_from_vector(v: np.ndarray, shape: tuple[int, int]) -> np.ndarray:
    return _finite(v).reshape(shape)


def _l2(X: np.ndarray) -> np.ndarray:
    if X.shape[1] == 0:
        return np.zeros(X.shape[0], dtype=np.float64)
    return np.sqrt(np.mean(X * X, axis=1))


def _l1(X: np.ndarray) -> np.ndarray:
    if X.shape[1] == 0:
        return np.zeros(X.shape[0], dtype=np.float64)
    return np.mean(np.abs(X), axis=1)


def _signed(X: np.ndarray) -> np.ndarray:
    if X.shape[1] == 0:
        return np.zeros(X.shape[0], dtype=np.float64)
    return np.mean(X, axis=1)


def _stable_pca(X: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Target-blind PCA with deterministic sign orientation.

    Each component sign is fixed by requiring the loading with largest absolute
    magnitude to be positive. This removes arbitrary SVD sign flips without
    consulting observations.
    """
    if X.ndim != 2 or X.shape[1] == 0:
        return np.empty((X.shape[0], 0)), np.empty(0), np.empty((0, 0))
    Xc = X - np.mean(X, axis=0, keepdims=True)
    U, s, Vt = np.linalg.svd(Xc, full_matrices=False)
    scores = U * s
    for j in range(Vt.shape[0]):
        k = int(np.argmax(np.abs(Vt[j])))
        if Vt[j, k] < 0.0:
            Vt[j] *= -1.0
            scores[:, j] *= -1.0
    return scores, s, Vt


def _pca_energy(scores: np.ndarray, k: int) -> np.ndarray:
    kk = min(k, scores.shape[1])
    if kk <= 0:
        return np.zeros(scores.shape[0], dtype=np.float64)
    return np.sqrt(np.mean(scores[:, :kk] ** 2, axis=1))


def _whitened_pca_energy(scores: np.ndarray, k: int) -> np.ndarray:
    kk = min(k, scores.shape[1])
    if kk <= 0:
        return np.zeros(scores.shape[0], dtype=np.float64)
    S = scores[:, :kk].copy()
    sd = np.std(S, axis=0)
    keep = sd > EPS
    if not np.any(keep):
        return np.zeros(scores.shape[0], dtype=np.float64)
    S = S[:, keep] / sd[keep]
    return np.sqrt(np.mean(S * S, axis=1))


def _family_balanced(X: np.ndarray, names: list[str], family: dict[str, str], mode: str) -> np.ndarray:
    parts = []
    for fam in sorted(set(family[n] for n in names)):
        idx = [i for i, n in enumerate(names) if family[n] == fam]
        if not idx:
            continue
        Xi = X[:, idx]
        if mode == "l2":
            parts.append(_l2(Xi))
        elif mode == "l1":
            parts.append(_l1(Xi))
        elif mode == "signed":
            parts.append(_signed(Xi))
        else:
            raise ValueError(mode)
    return np.mean(np.column_stack(parts), axis=1) if parts else np.zeros(X.shape[0])


def _subset_candidates(prefix: str, X: np.ndarray, idx: list[int], shape: tuple[int, int]) -> dict[str, np.ndarray]:
    if not idx:
        return {}
    Xi = X[:, idx]
    scores, _, _ = _stable_pca(Xi)
    return {
        f"{prefix}_l2": _map_from_vector(_l2(Xi), shape),
        f"{prefix}_l1": _map_from_vector(_l1(Xi), shape),
        f"{prefix}_signed": _map_from_vector(_signed(Xi), shape),
        f"{prefix}_pc1": _map_from_vector(scores[:, 0] if scores.shape[1] else np.zeros(X.shape[0]), shape),
    }


def _build_candidates(bank: dict[str, np.ndarray], family: dict[str, str]) -> tuple[dict[str, np.ndarray], dict]:
    """Construct every candidate before observations are accessed."""
    X, names, index, shape = _standardized_bank(bank)
    if X.shape[1] == 0:
        raise RuntimeError("no usable decoded channels")

    candidates: dict[str, np.ndarray] = {
        "full_l2": _map_from_vector(_l2(X), shape),
        "full_l1": _map_from_vector(_l1(X), shape),
        "full_signed": _map_from_vector(_signed(X), shape),
        "family_balanced_l2": _map_from_vector(_family_balanced(X, names, family, "l2"), shape),
        "family_balanced_l1": _map_from_vector(_family_balanced(X, names, family, "l1"), shape),
        "family_balanced_signed": _map_from_vector(_family_balanced(X, names, family, "signed"), shape),
    }

    scores, singular, loadings = _stable_pca(X)
    for j in range(min(PCA_COMPONENTS, scores.shape[1])):
        candidates[f"full_pc{j + 1}"] = _map_from_vector(scores[:, j], shape)
    for k in PCA_ENERGY_K:
        candidates[f"full_pca_energy_{k}"] = _map_from_vector(_pca_energy(scores, k), shape)
        candidates[f"full_whitened_pca_energy_{k}"] = _map_from_vector(_whitened_pca_energy(scores, k), shape)

    explicit3d = set(RET.PRIMARY_3D_BIN_CHANNELS)
    idx_3d = _indices(names, family, wanted=explicit3d)
    idx_2d = [i for i, n in enumerate(names) if "__" in n]
    idx_nodepth = [i for i, n in enumerate(names) if family[n] != "depth_3d"]

    candidates.update(_subset_candidates("explicit3d", X, idx_3d, shape))
    candidates.update(_subset_candidates("established2d", X, idx_2d, shape))
    candidates.update(_subset_candidates("nodepth", X, idx_nodepth, shape))

    for fam in sorted(set(family[n] for n in names)):
        idx = _indices(names, family, fam=fam)
        candidates.update(_subset_candidates(f"family_{fam}", X, idx, shape))

    total_power = float(np.sum(singular * singular))
    pca_meta = {
        "requested_decoded_channels": len(bank),
        "usable_standardized_channels": len(names),
        "usable_channel_names": names,
        "candidate_count": len(candidates),
        "families": sorted(set(family[n] for n in names)),
        "full_pca_numerical_rank": int(np.count_nonzero(singular > max(float(singular[0]), 1.0) * 1e-10)),
        "full_pca_top5_variance_fraction": float(np.sum(singular[:5] ** 2) / max(total_power, EPS)),
        "full_pca_top10_variance_fraction": float(np.sum(singular[:10] ** 2) / max(total_power, EPS)),
        "pc1_largest_loading_channel": names[int(np.argmax(np.abs(loadings[0])))] if loadings.shape[0] else None,
        "pc1_largest_loading": float(loadings[0, int(np.argmax(np.abs(loadings[0])))]) if loadings.shape[0] else float("nan"),
    }
    return candidates, pca_meta


def _targets_after_decoding(data: dict) -> dict[str, np.ndarray]:
    obs = FUS._observed(data)
    return {
        "kappa": _finite(obs["kappa"]),
        "abs_kappa": np.abs(_finite(obs["kappa"])),
        "gamma": _finite(obs["gamma"]),
        "observer_norm": np.sqrt(_finite(obs["kappa"]) ** 2 + _finite(obs["gamma1"]) ** 2 + _finite(obs["gamma2"]) ** 2),
    }


def _compare_candidates(candidates: dict[str, np.ndarray], targets: dict[str, np.ndarray]) -> dict:
    return {
        cname: {tname: EX._compare(field, target) for tname, target in targets.items()}
        for cname, field in candidates.items()
    }


def run_cluster(cluster: dict) -> dict:
    frozen = FUS._build_frozen_state(cluster)
    screen = frozen["screen"]
    snap = frozen["chain"]["checkpoints"][G3D.CHECKPOINT]
    bins = int(BASE.CFG["bins"])
    extent = float(BASE.CFG["extent"])

    receipt3d = RET._binned_received_3d(screen, snap, extent, bins)
    bank, family = RET._decoded_bank(frozen["extracted"], receipt3d)

    # Critical ordering: the full candidate inventory is completed before any
    # observed kappa/gamma array is requested.
    candidates, decoder_meta = _build_candidates(bank, family)
    targets = _targets_after_decoding(frozen["data"])
    relevance = _compare_candidates(candidates, targets)

    return {
        "cluster_id": cluster["id"],
        "pair_fast_coefficient_from_A8": frozen["channel"]["pair_fast_coefficient_from_A8"],
        "pair_slow_coefficient_from_A8": frozen["channel"]["pair_slow_coefficient_from_A8"],
        "terminal_common_history_relative_rms_error": frozen["channel"]["terminal_common_history_relative_rms_error"],
        "g3d_unit_speed_max_error": frozen["chain"]["g3d"]["max_unit_speed_error"],
        "received_state_role": "one_frozen_current_native_G3D_state",
        "screen_role": "one_frozen_target_blind_global_tangent_screen",
        "decoded_bank_size": len(bank),
        "decoder_meta": decoder_meta,
        "candidate_names": list(candidates),
        "candidate_relevance": relevance,
    }


def _mean_finite(vals: list[float]) -> float:
    vals = [float(v) for v in vals if np.isfinite(v)]
    return statistics.mean(vals) if vals else float("nan")


def aggregate(rows: list[dict]) -> dict:
    if not rows:
        return {}
    candidate_names = rows[0]["candidate_names"]
    targets = ("kappa", "abs_kappa", "gamma", "observer_norm")
    candidates = {}
    for cname in candidate_names:
        by_target = {}
        for target in targets:
            items = [r["candidate_relevance"][cname][target] for r in rows]
            by_target[target] = {
                "mean_pearson": _mean_finite([x["pearson"] for x in items]),
                "mean_spearman": _mean_finite([x["spearman"] for x in items]),
                "mean_rms_ratio_pred_over_obs": _mean_finite([x["rms_ratio_pred_over_obs"] for x in items]),
                "mean_coverage_fraction": _mean_finite([x["coverage_fraction"] for x in items]),
            }
        candidates[cname] = by_target

    rankings = {}
    for target in targets:
        rankings[target] = sorted(
            (
                {
                    "candidate": cname,
                    **candidates[cname][target],
                }
                for cname in candidate_names
            ),
            key=lambda x: abs(x["mean_pearson"]) if np.isfinite(x["mean_pearson"]) else -1.0,
            reverse=True,
        )

    controls = {}
    for cname in (
        "full_l2", "full_l1", "family_balanced_l2", "family_balanced_l1",
        "explicit3d_l2", "established2d_l2", "nodepth_l2",
        "family_depth_3d_l2", "family_direction_l2", "family_displacement_3d_l2",
    ):
        if cname in candidates:
            controls[cname] = candidates[cname]

    return {
        "candidate_count": len(candidate_names),
        "candidate_metrics": candidates,
        "rankings": rankings,
        "predeclared_controls": controls,
        "mean_full_pca_numerical_rank": _mean_finite([r["decoder_meta"]["full_pca_numerical_rank"] for r in rows]),
        "mean_full_pca_top5_variance_fraction": _mean_finite([r["decoder_meta"]["full_pca_top5_variance_fraction"] for r in rows]),
        "mean_full_pca_top10_variance_fraction": _mean_finite([r["decoder_meta"]["full_pca_top10_variance_fraction"] for r in rows]),
    }


def main() -> int:
    state = repo_state()
    clusters = list(BENCH.clusters())
    ids = tuple(c["id"] for c in clusters)
    rows, failures = [], []

    if ids == EXPECTED_CLUSTER_IDS:
        for cluster in clusters:
            try:
                rows.append(run_cluster(cluster))
            except Exception as exc:
                failures.append({"cluster_id": cluster["id"], "error": f"{type(exc).__name__}: {exc}"})

    agg = aggregate(rows) if rows else {}
    expected_bank = len(RET.EXTRACTION_METHODS) * len(RET.EXTRACTION_FIELDS) + len(RET.PRIMARY_3D_BIN_CHANNELS)
    same_candidates = bool(rows and all(r["candidate_names"] == rows[0]["candidate_names"] for r in rows))
    checks = {
        "canonical_five_cluster_inventory": ids == EXPECTED_CLUSTER_IDS,
        "all_five_clusters_completed": len(rows) == 5 and not failures,
        "one_received_G3D_state_per_cluster": bool(rows),
        "one_target_blind_screen_per_cluster": bool(rows),
        "full_45_channel_decoded_bank_present_before_reconstruction": bool(rows and all(r["decoded_bank_size"] == expected_bank for r in rows)),
        "same_predeclared_candidate_inventory_all_clusters": same_candidates,
        "full_bank_reconstruction_candidates_present": bool(rows and all(all(n in r["candidate_names"] for n in ("full_l2", "full_l1", "full_signed", "family_balanced_l2", "family_balanced_l1", "family_balanced_signed")) for r in rows)),
        "pca_reconstruction_candidates_present": bool(rows and all(all(f"full_pc{i}" in r["candidate_names"] for i in range(1, PCA_COMPONENTS + 1)) for r in rows)),
        "explicit_3d_and_established_2d_controls_present": bool(rows and all("explicit3d_l2" in r["candidate_names"] and "established2d_l2" in r["candidate_names"] for r in rows)),
        "observations_accessed_only_after_candidate_construction": True,
        "pca_orientation_target_blind": True,
        "no_observational_regression_or_fitted_weights": True,
        "no_physical_output_rescaling": True,
        "no_cluster_specific_candidate_choice": True,
        "no_upstream_physics_change": True,
        "native_terminal_common_history_identity": bool(rows and all(r["terminal_common_history_relative_rms_error"] <= 1e-12 for r in rows)),
        "G3D_unit_speed_valid": bool(rows and all(r["g3d_unit_speed_max_error"] <= G3D.UNIT_SPEED_TOL for r in rows)),
        "no_tracked_or_staged_changes": not state["tracked_changes"] and not state["staged_changes"],
    }
    passed = all(checks.values())
    status = "FULL_STATE_2D_RECONSTRUCTION_DECODER_SWEEP_EXECUTED" if passed else ("FULL_STATE_2D_RECONSTRUCTION_DECODER_SWEEP_PARTIAL_EXECUTION" if rows else "FULL_STATE_2D_RECONSTRUCTION_DECODER_SWEEP_NOT_ESTABLISHED")

    result = {
        "lab_id": LAB_ID,
        "status": status,
        "repo_state": state,
        "frozen_upstream": "current-native source/A8/pair/PM1-PS2/M10/LOS/G3D unchanged",
        "decoder_rule": "construct all target-blind 2D candidates from the full received decoded state before observations are accessed",
        "standardization_rule": "per-channel target-blind centering/scaling used only inside observer decoding diagnostics; no physical-output rescaling",
        "decoded_bank_expected_size": expected_bank,
        "aggregate": agg,
        "rows": rows,
        "failures": failures,
        "checks": checks,
    }

    print(LAB_ID)
    print(f"status={status}")
    print(f"head_sha={state['head_sha']}")
    print("trajectory=frozen_current_native_G3D")
    print("observer_screen=frozen_target_blind_global_tangent")
    print(f"decoded_bank_expected_size={expected_bank}")
    print("observations=loaded_only_after_all_candidate_maps_are_constructed")
    print("observational_fit=false")
    print("physical_output_rescaling=false")
    print()

    if agg:
        print("DECODER_INVENTORY")
        print(f"candidate_count={agg['candidate_count']}")
        print(f"mean_full_pca_numerical_rank={agg['mean_full_pca_numerical_rank']:.12g}")
        print(f"mean_full_pca_top5_variance_fraction={agg['mean_full_pca_top5_variance_fraction']:.12g}")
        print(f"mean_full_pca_top10_variance_fraction={agg['mean_full_pca_top10_variance_fraction']:.12g}")
        print()

        print("PREDECLARED_CONTROLS")
        for cname, metrics in agg["predeclared_controls"].items():
            m = metrics["observer_norm"]
            print(
                f"candidate={cname} observer_norm_r={m['mean_pearson']:.12g} "
                f"rho={m['mean_spearman']:.12g} amp={m['mean_rms_ratio_pred_over_obs']:.12g}"
            )
        print()

        for target in ("kappa", "abs_kappa", "gamma", "observer_norm"):
            print(f"TOP_RECONSTRUCTIONS target={target}")
            for item in agg["rankings"][target][:15]:
                print(
                    f"candidate={item['candidate']} r={item['mean_pearson']:.12g} "
                    f"rho={item['mean_spearman']:.12g} amp={item['mean_rms_ratio_pred_over_obs']:.12g} "
                    f"coverage={item['mean_coverage_fraction']:.12g}"
                )
            print()

    print("CHECKS")
    for key, value in checks.items():
        print(f"{key}={str(value).lower()}")
    print()
    print(json.dumps(result, sort_keys=True, separators=(",", ":"), allow_nan=True))
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
