#!/usr/bin/env python3
"""PBUF FOUNDATION — MATTER LOADING PHYSICAL CLOSURE 001.

Purpose
-------
Establish the no-fit physical bridge requirements between measured matter
stress-energy and the PBUF medium.  This lab deliberately does NOT run lensing,
does NOT use observed kappa, and does NOT invent a replacement for the legacy
`strength = 0.18` coefficient.

The action-level architecture under audit is

    S = S_EH[g] + S_sigma[g, chi; PBUF microphysics] + S_m[g, Psi]

and standard minimal matter coupling gives the medium source functional

    J_A(x) = -(1/(2 sqrt(-g(x)))) \int d^4y sqrt(-g(y))
             T^{mu nu}(y) delta G_{mu nu}(y)/delta chi^A(x)

or, for a local algebraic metric-medium map,

    J_A = -(1/2) T^{mu nu} partial G_{mu nu}/partial chi^A.

This is important because it leaves no independent matter-coupling knob to fit.
The stress-energy T^{mu nu} is measurable/physical, while the response kernel
R_{mu nu A} = delta G_{mu nu}/delta chi^A must come from the PBUF metric-medium
map.  The medium dynamics must come from S_sigma.  Only after those are closed
may a coarse-grained scalar/vector loading field be derived.

The lab therefore:

1. verifies the local source contraction algebra with deterministic synthetic
   tensors (algebra test only; not a physical model);
2. verifies linearity in T and in the metric-medium response kernel;
3. inventories the current repository for the known legacy normalized-source
   and strength-amplitude construction;
4. checks whether the required physical closure objects are present;
5. emits a machine-readable bridge specification for the next implementation;
6. explicitly refuses to derive mass from HST light or use any fitted scalar.

No GR/Newtonian lensing law is injected into the PBUF response.  Standard
stress-energy coupling is used only to define what a physically closed PBUF
matter source must derive from its own metric-medium map.
"""
from __future__ import annotations

import hashlib
import json
import math
import subprocess
import sys
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))

LAB_ID = "PBUF-FOUNDATION-MATTER-LOADING-PHYSICAL-CLOSURE-001"
OUT = ROOT / "runs" / "matter_loading_physical_closure001"

LEGACY_STRENGTH = 0.18
ALG_TOL = 1.0e-14

# Files that establish the currently active placeholder architecture.
AUDIT_FILES = (
    ROOT / "constitutive_equations.py",
    ROOT / "pbuf" / "labs" / "foundation" / "m10_coverage_25pct_science001.py",
    ROOT / "pbuf" / "labs" / "foundation" / "strength_factorization_physical_bridge001.py",
)

# Search terms are deliberately structural.  Their presence alone is not proof
# of closure; they are used only for a transparent repository inventory.
CLOSURE_MARKERS = {
    "metric_medium_map": (
        "delta_g_delta_chi", "metric_medium_map", "metric_response_kernel",
        "delta G", "dG_dchi",
    ),
    "medium_action": (
        "S_sigma", "medium_action", "elastic_action", "vacuum_action",
    ),
    "stress_energy_source_functional": (
        "J_A", "stress_energy_source", "matter_source_functional",
    ),
    "coarse_graining_projection": (
        "coarse_graining", "coarse_grain", "source_projection", "projection_operator",
    ),
    "physical_mass_density_si": (
        "kg_m3", "kg/m^3", "mass_density_si", "surface_density_kg_m2",
    ),
}


@dataclass(frozen=True)
class BridgeRequirement:
    key: str
    role: str
    required: bool
    present_in_current_foundation: bool
    may_be_fitted: bool
    derivation_rule: str


@dataclass(frozen=True)
class SourceFunctionalSpec:
    action_architecture: str
    exact_nonlocal_source: str
    exact_local_source: str
    stress_energy_role: str
    metric_response_role: str
    independent_matter_coupling_constant_allowed: bool


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


def _json_default(x):
    if isinstance(x, (np.floating,)): return float(x)
    if isinstance(x, (np.integer,)): return int(x)
    if isinstance(x, (np.bool_,)): return bool(x)
    if isinstance(x, Path): return str(x)
    return str(x)


def _write_json(path: Path, obj) -> None:
    path.write_text(json.dumps(obj, indent=2, default=_json_default))


def _sha_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _rel_rms(a, b) -> float:
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    d = a - b
    den = math.sqrt(float(np.mean(b*b))) if b.size else 0.0
    num = math.sqrt(float(np.mean(d*d))) if d.size else 0.0
    return num / max(den, 1.0e-300)


def local_matter_source(T_up: np.ndarray, dG_dchi: np.ndarray) -> np.ndarray:
    """Compute J_A = -1/2 T^{mu nu} dG_{mu nu}/dchi^A.

    Parameters
    ----------
    T_up:
        Symmetric rank-2 stress-energy tensor, shape (4,4).
    dG_dchi:
        Metric-medium response derivatives, shape (nchi,4,4).

    This routine is algebra only.  It does not assert a particular PBUF
    dG/dchi.  Calling it with a synthetic kernel is a deterministic identity
    test, not a physical calculation.
    """
    T = np.asarray(T_up, dtype=np.float64)
    R = np.asarray(dG_dchi, dtype=np.float64)
    if T.shape != (4, 4):
        raise ValueError(f"T must have shape (4,4), got {T.shape}")
    if R.ndim != 3 or R.shape[1:] != (4, 4):
        raise ValueError(f"dG_dchi must have shape (nchi,4,4), got {R.shape}")
    if not np.allclose(T, T.T, rtol=0.0, atol=ALG_TOL):
        raise ValueError("T must be symmetric")
    if not np.allclose(R, np.swapaxes(R, 1, 2), rtol=0.0, atol=ALG_TOL):
        raise ValueError("each metric response tensor must be symmetric")
    return -0.5 * np.einsum("mn,amn->a", T, R, optimize=True)


def _algebra_tests() -> dict:
    # Deterministic synthetic tensors only.  Values have no physical meaning.
    T = np.array([
        [7.0, 0.2, -0.1, 0.0],
        [0.2, 2.0, 0.3, -0.2],
        [-0.1, 0.3, 3.0, 0.4],
        [0.0, -0.2, 0.4, 4.0],
    ], dtype=np.float64)
    R = np.array([
        [[1.0, .1, 0, 0], [.1, .3, 0, 0], [0, 0, -.2, 0], [0, 0, 0, .4]],
        [[-.5, 0, .2, 0], [0, .1, 0, 0], [.2, 0, .7, .1], [0, 0, .1, -.3]],
        [[.2, 0, 0, .1], [0, -.4, .2, 0], [0, .2, .2, 0], [.1, 0, 0, .9]],
    ], dtype=np.float64)

    J = local_matter_source(T, R)
    explicit = np.array([
        -0.5 * sum(T[m, n] * R[a, m, n] for m in range(4) for n in range(4))
        for a in range(R.shape[0])
    ])

    a, b = 2.75, -0.6
    T2 = np.array([
        [1.2, 0, 0, 0], [0, -.2, .1, 0], [0, .1, .8, 0], [0, 0, 0, .5]
    ], dtype=np.float64)
    left_T = local_matter_source(a*T + b*T2, R)
    right_T = a*local_matter_source(T, R) + b*local_matter_source(T2, R)

    R2 = 0.3 * R[::-1].copy()
    left_R = local_matter_source(T, a*R + b*R2)
    right_R = a*local_matter_source(T, R) + b*local_matter_source(T, R2)

    zero_T = local_matter_source(np.zeros((4,4)), R)
    zero_R = local_matter_source(T, np.zeros_like(R))

    return {
        "synthetic_only_not_physical_model": True,
        "explicit_contraction_relative_rms_error": _rel_rms(J, explicit),
        "linearity_in_stress_energy_relative_rms_error": _rel_rms(left_T, right_T),
        "linearity_in_metric_response_relative_rms_error": _rel_rms(left_R, right_R),
        "zero_stress_energy_source_norm": float(np.linalg.norm(zero_T)),
        "zero_metric_response_source_norm": float(np.linalg.norm(zero_R)),
        "explicit_contraction_pass": bool(_rel_rms(J, explicit) <= ALG_TOL),
        "linearity_in_stress_energy_pass": bool(_rel_rms(left_T, right_T) <= ALG_TOL),
        "linearity_in_metric_response_pass": bool(_rel_rms(left_R, right_R) <= ALG_TOL),
        "zero_source_pass": bool(np.linalg.norm(zero_T) <= ALG_TOL and np.linalg.norm(zero_R) <= ALG_TOL),
    }


def _iter_repo_text_files() -> Iterable[Path]:
    skip = {".git", "runs", "__pycache__", ".venv", "venv"}
    for path in ROOT.rglob("*"):
        if not path.is_file():
            continue
        if any(part in skip for part in path.parts):
            continue
        if path.suffix.lower() not in {".py", ".md", ".txt"}:
            continue
        yield path


def _marker_inventory() -> dict:
    hits = {k: [] for k in CLOSURE_MARKERS}
    for path in _iter_repo_text_files():
        try:
            text = path.read_text(errors="ignore")
        except Exception:
            continue
        rel = str(path.relative_to(ROOT))
        for key, markers in CLOSURE_MARKERS.items():
            found = [m for m in markers if m in text]
            if found:
                hits[key].append({"path": rel, "markers": found})
    return hits


def _legacy_inventory() -> dict:
    rows = []
    for path in AUDIT_FILES:
        item = {
            "path": str(path.relative_to(ROOT)),
            "exists": path.exists(),
            "sha256": _sha_file(path) if path.exists() else None,
        }
        if path.exists():
            text = path.read_text(errors="ignore")
            item.update({
                "contains_legacy_strength_0p18": ("0.18" in text),
                "contains_normalized_matter": ("_normalized(matter)" in text or "rho/rho_max" in text),
                "contains_STRENGTH_times_rho3": ("STRENGTH * rho3" in text),
                "contains_physical_mass_calibration": any(x in text for x in ("mass_density_si", "kg/m^3", "surface_density_kg_m2")),
            })
        rows.append(item)
    return {"files": rows}


def _requirements(marker_hits: dict) -> list[BridgeRequirement]:
    # Repository marker hits are an inventory only.  For closure we require a
    # normalized implementation that is explicitly identified as physical.
    # None of the current foundation artifacts supplies that closed chain.
    return [
        BridgeRequirement(
            "physical_stress_energy_tensor",
            "Independent physical T^{mu nu}(x), including energy density and where relevant pressure/momentum/anisotropic stress.",
            True, False, False,
            "Must come from independent observations/physics in SI or a fully documented unit system; never from kappa.",
        ),
        BridgeRequirement(
            "metric_medium_response_kernel",
            "R_{mu nu A}=delta G_{mu nu}/delta chi^A, or an equivalent normalized metric-medium map.",
            True, False, False,
            "Must be derived from the PBUF definition of the medium-to-metric map; no free matter coupling coefficient.",
        ),
        BridgeRequirement(
            "medium_action_or_closed_dynamics",
            "S_sigma or equivalent equations that determine propagation, stiffness and normalization of chi.",
            True, False, False,
            "Must be derived from frozen PBUF microphysics/constitutive dynamics and pass stability/causality checks.",
        ),
        BridgeRequirement(
            "coarse_graining_and_mode_projection",
            "Normalized C and P operators mapping chi/J_A into the finite lab state u and source s.",
            True, False, False,
            "Must be derived/documented before reducing full T^{mu nu} to a scalar/vector loading field.",
        ),
        BridgeRequirement(
            "absolute_source_geometry",
            "Physical cell areas/volumes/depths and coordinate units needed to convert observed matter to T^{mu nu}(x).",
            True, False, False,
            "Must come from measured geometry/redshifts/calibration, not from lensing benchmark morphology.",
        ),
    ]


def _bridge_spec(reqs: list[BridgeRequirement]) -> dict:
    source = SourceFunctionalSpec(
        action_architecture="S = S_EH[g] + S_sigma[g,chi;PBUF microphysics] + S_m[g,Psi]",
        exact_nonlocal_source=(
            "J_A(x)=-(1/(2 sqrt(-g(x)))) integral d4y sqrt(-g(y)) "
            "T_m^{mu nu}(y) delta G_{mu nu}(y)/delta chi^A(x)"
        ),
        exact_local_source="J_A=-(1/2) T_m^{mu nu} partial G_{mu nu}/partial chi^A",
        stress_energy_role="physical measured matter source; not a fitted amplitude",
        metric_response_role="PBUF-derived response kernel fixing which stress-energy projection loads each medium mode",
        independent_matter_coupling_constant_allowed=False,
    )
    return {
        "source_functional": asdict(source),
        "requirements": [asdict(r) for r in reqs],
        "forbidden_shortcuts": [
            "replace_0p18_with_another_fitted_scalar",
            "derive_cluster_mass_from_F160W_using_assumed_or_tuned_mass_to_light_ratio",
            "use_observed_kappa_or_shear_upstream",
            "reduce_Tmunu_to_rho_only_before_metric_medium_projection_is_derived",
            "use_alpha_epsilon0_Rmax_or_ksat_as_local_mass_coupling_without_derivation",
            "inject_Newtonian_or_GR_lensing_force_as_the_PBUF_response",
        ],
        "required_next_derivation_order": [
            "derive_or_freeze_metric_medium_map_g_of_chi",
            "derive_R_munuA_deltaG_deltaChi_and_normalization",
            "derive_medium_dynamics_S_sigma_or_equivalent",
            "derive_coarse_graining_C_and_projection_P",
            "acquire_independent_physical_Tmunu_source_with_absolute_geometry",
            "compute_J_A_without_free_coupling_scalar",
            "map_J_A_through_closed_medium_dynamics_to_initial_deformation_response",
            "only_then_run_M10_G3D_observer_and_compare_with_lensing_at_the_end",
        ],
    }


def main() -> int:
    started = time.perf_counter()
    OUT.mkdir(parents=True, exist_ok=True)

    repo = _repo_state()
    if repo["branch"] != "main":
        raise RuntimeError(f"runner must execute on main, got {repo['branch']}")
    if repo["tracked_changes"] or repo["staged_changes"]:
        raise RuntimeError("tracked or staged repository changes present")

    algebra = _algebra_tests()
    if not all(algebra[k] for k in (
        "explicit_contraction_pass",
        "linearity_in_stress_energy_pass",
        "linearity_in_metric_response_pass",
        "zero_source_pass",
    )):
        raise RuntimeError("source-functional algebra gate failed")

    markers = _marker_inventory()
    legacy = _legacy_inventory()
    reqs = _requirements(markers)
    spec = _bridge_spec(reqs)

    unresolved = [r.key for r in reqs if r.required and not r.present_in_current_foundation]
    closure_complete = (len(unresolved) == 0)

    validation = {
        "lab_id": LAB_ID,
        "outcome": "Outcome B — PHYSICAL MATTER-LOADING BRIDGE NOT YET CLOSED; NO FITTED REPLACEMENT PERMITTED",
        "head_sha": repo["head_sha"],
        "benchmark_pixel_values_loaded": False,
        "lensing_pipeline_executed": False,
        "hst_mass_conversion_executed": False,
        "legacy_strength": LEGACY_STRENGTH,
        "legacy_strength_physical_derivation_present": False,
        "legacy_strength_replacement_selected": False,
        "fitted_numbers_introduced": False,
        "independent_matter_coupling_constant_allowed": False,
        "stress_energy_is_fixed_physical_source": True,
        "metric_medium_response_kernel_required": True,
        "medium_action_or_closed_dynamics_required": True,
        "coarse_graining_projection_required": True,
        "physical_stress_energy_dataset_present": False,
        "physical_bridge_complete": closure_complete,
        "unresolved_required_components": unresolved,
        "local_source_formula": "J_A=-(1/2) T_m^{mu nu} partial G_{mu nu}/partial chi^A",
        "full_nonlocal_source_formula": "J_A(x)=-(1/(2 sqrt(-g(x)))) integral d4y sqrt(-g(y)) T_m^{mu nu}(y) delta G_{mu nu}(y)/delta chi^A(x)",
        "source_reduction_to_rho_only_authorized": False,
        "replace_legacy_strength_now_authorized": False,
        "physical_cluster_lensing_run_authorized": False,
        "next_derivation_authorized": True,
        "next_derivation_target": "metric_medium_map_and_medium_action_closure_before_physical_cluster_mass_loading",
        "algebra": algebra,
        "duration_seconds": time.perf_counter() - started,
    }

    _write_json(OUT / "validation.json", validation)
    _write_json(OUT / "physical_bridge_spec.json", spec)
    _write_json(OUT / "repository_marker_inventory.json", markers)
    _write_json(OUT / "legacy_placeholder_inventory.json", legacy)

    print(validation["outcome"])
    print(json.dumps(validation, indent=2))
    print(f"output_directory={OUT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
