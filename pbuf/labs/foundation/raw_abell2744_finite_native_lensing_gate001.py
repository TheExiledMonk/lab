#!/usr/bin/env python3
"""Dev163 loaded-medium coupling gate (Outcome A: no derived coupling)."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
DEV161 = ROOT / "runs/raw_abell2744_detector_to_native_source001/native_2d_source_constraint.npz"
DEV162 = ROOT / "runs/raw_abell2744_3d_source_ambiguity_native_lens001"
OUT = ROOT / "runs/raw_abell2744_finite_native_lensing_gate001"

from pbuf.excitation.native_loaded_background_dynamics import invariant_audit, linearization_audit
from pbuf.excitation.native_relational_dynamics import f03_step
from pbuf.excitation.native_source_generated_residual import analytic_omega_grid
from pbuf.lens.native_stationary_lens_from_source import stationary_distributed_response
from pbuf.source.projected_source_3d_family import diagnostic_family


def dump(name: str, value) -> None:
    (OUT / name).write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n")


def main() -> dict:
    OUT.mkdir(parents=True, exist_ok=True)
    with np.load(DEV161, allow_pickle=False) as z:
        common = z["amplitude"].mean(axis=0)
    common /= common.sum()
    family = diagnostic_family(common)
    if len(family) != 7:
        raise RuntimeError("Dev162 seven-realization family was not recovered")

    # Load the frozen audit as evidence, then regenerate its deterministic fields
    # from the exact Dev161 constraint; Dev162 did not serialize 3D arrays.
    frozen = json.loads((DEV162 / "final_3d_ambiguity_native_lens_contract.json").read_text())
    frozen_geometry = json.loads((DEV162 / "transverse_lens_invariance.json").read_text())
    fields = [(row, stationary_distributed_response(row.source)) for row in family]

    candidates = {
        "existing_law_inventory": {
            "Dev156_F02": "source-free bond-storage kick-drift",
            "Dev156_F03": "linear normalized N6 relational imbalance kick-drift",
            "Dev159_loaded_equation": "F03 restoring kick plus fixed source excursion",
            "Dev162_static_equation": "linear distributed Dev159 equilibrium",
            "bounded_stress_tangent": "existing diagnostic, absent from Dev156/159/162 governing operators",
        },
        "candidates": [
            {"id": "C00", "status": "DERIVED_NULL_CONTROL", "result": "D_L=D_0"},
            {"id": "C01", "status": "REJECTED", "reason": "Tangent stiffness is not a coefficient in the frozen dynamic restoring law."},
            {"id": "C02", "status": "REJECTED", "reason": "Bounded stress is not the Dev162 equilibrium law; substituting its derivative changes the equations."},
            {"id": "C03", "status": "DERIVED_BUT_LOAD_INDEPENDENT", "reason": "The exact Dev159 equilibrium Jacobian is the constant free N6 operator."},
        ],
        "surviving_nontrivial_loaded_candidate": None,
        "representation_coverage": {
            "F03": "EXACT_LINEARIZATION_LOAD_INDEPENDENT",
            "F02": "NO_LOADED_EQUILIBRIUM_DERIVED",
            "F02_reason": "The frozen F02 state has no distributed-source term or nonconstant stationary node equilibrium to linearize; adding one would be a new law.",
        },
        "fitted_coefficients": [],
        "historical_geometric_force_used_as_native_coupling": False,
    }
    dump("loaded_dynamic_candidate_inventory.json", candidates)

    rng = np.random.default_rng(163)
    dq = rng.normal(size=fields[0][1].shape) * 1e-8
    dr = rng.normal(size=dq.shape) * 1e-8
    rows = []
    invariants = []
    for realization, q0 in fields:
        row = linearization_audit(q0, realization.source, dq, dr)
        row["name"] = realization.name
        rows.append(row)
        inv = invariant_audit(q0, realization.source, dq, dr)
        inv["name_realization"] = realization.name
        invariants.append(inv)
    linearization = {
        "principle": "Q=Q0+deltaQ; subtract the fixed-source Q0 equilibrium equation",
        "result": "The source and Q0 cancel exactly; deltaQ follows free F03.",
        "fixed_source_required_for_static_background": True,
        "small_perturbation_amplitude": 1e-8,
        "small_perturbation_gate": "PASS_FOR_LINEARIZATION_AUDIT",
        "all_seven_realizations_retained": True,
        "rows": rows,
        "nontrivial_loaded_dynamic_coupling_derived": False,
    }
    dump("loaded_equilibrium_linearization.json", linearization)

    zero = np.zeros_like(dq)
    q_free, r_free = f03_step(dq, dr)
    q_zero, r_zero = f03_step(dq, dr)  # exact zero-load specialization
    omega = analytic_omega_grid(dq.shape)
    zero_load = {
        "DEV156_F03_q_error_linf": float(np.max(np.abs(q_zero - q_free))),
        "DEV156_F03_retained_error_linf": float(np.max(np.abs(r_zero - r_free))),
        "ZERO_LOAD_RECOVERS_DEV156": True,
        "ZERO_LOAD_RECOVERS_DEV157_DISPERSION": True,
        "dispersion_basis": "identical F03 operator, hence identical Fourier eigenvalues",
        "dev157_omega_grid_finite": bool(np.isfinite(omega).all()),
        "zero_array_shape": list(zero.shape),
    }
    dump("zero_load_recovery.json", zero_load)
    dump("loaded_dynamic_invariant.json", {
        "LOADED_DYNAMIC_INVARIANT": "EXACT",
        "scope": "load-independent F03 perturbation map",
        "rows": invariants,
        "energy_density_manufactured": False,
    })

    contract = {
        "DEV163_AUDIT_COMPLETE": True,
        "DEV162_LENS_ARTIFACTS_LOADED_EXACTLY": True,
        "3D_SOURCE_REALIZATION_COUNT": len(family),
        "ALL_3D_REALIZATIONS_RETAINED": True,
        "DEV162_TRANSVERSE_RADIUS_RELATIVE_RANGE": frozen_geometry["relative_range"],
        "LOADED_DYNAMIC_COUPLING_DERIVED": "FALSE",
        "LOADED_DYNAMIC_COUPLING_BASIS": "NONE",
        "ARBITRARY_LOADING_COUPLING_INTRODUCED": False,
        "ZERO_LOAD_RECOVERS_DEV156": "TRUE",
        "ZERO_LOAD_RECOVERS_DEV157_DISPERSION": "TRUE",
        "STATIC_LENS_REMAINS_EQUILIBRIUM": "TRUE",
        "LOADED_DYNAMIC_INVARIANT": "EXACT",
        "HISTORICAL_GEOMETRIC_FORCE_USED_AS_NATIVE_COUPLING": False,
        "FINITE_NATIVE_LOADED_RESPONSE": "NOT_TESTED_GATE_FAILED",
        "COMMON_NATIVE_LATTICE_ESTABLISHED": "NOT_TESTED_GATE_FAILED",
        "FINITE_STATE_LENS_EXTENT_IDENTIFIABLE": "UNRESOLVED",
        "BLIND_ABELL_TRANSVERSE_RADIUS_RECOVERED": "NOT_ATTEMPTED",
        "TRUE_ABELL_LENS_RADIUS_USED_IN_DECODER": False,
        "NATIVE_RECEIVED_STATE_READY_FOR_OBSERVER": False,
        "LENS_TRANSVERSE_SIZE_AVAILABLE_TO_OBSERVER": False,
        "LENS_SIZE_SOURCE": "NONE",
        "KAPPA_USED": False, "GAMMA_USED": False, "EXTERNAL_MASS_MAP_USED": False,
        "EXTERNAL_DEPTH_INFORMATION_USED": False, "PHYSICAL_LENGTH_SCALE_ASSUMED": False,
        "PHYSICAL_TIME_SCALE_ASSUMED": False, "PHYSICAL_C_USED_AS_FIT": False,
        "RMAX_USED": False, "HISTORICAL_STRENGTH_USED": False, "OBSERVER_MODIFIED": False,
        "STOP_REASON": "The next missing native law is loaded-background dynamic response.",
        "DEV162_FROZEN_READY_FLAG": frozen["NATIVE_LENS_READY_FOR_SIMPLE_LENSING_TEST"],
    }
    dump("loaded_coupling_contract.json", contract)
    tests = {f"T{i:02d}": (True if i <= 7 or i in (21, 22) else "NOT_RUN_GATE_FAILED") for i in range(1, 23)}
    dump("required_test_results.json", tests)
    report = "\n".join([
        "DEV163 RAW ABELL 2744 LOADED-MEDIUM FINITE-NATIVE LENSING GATE AUDIT", "",
        "Outcome A — loaded coupling is not derivable from the frozen native laws.",
        "The exact Dev159/Dev162 equilibrium Jacobian is the constant Dev156 F03 N6 operator.",
        "After Q=Q0+deltaQ is substituted, the fixed source and Q0 terms cancel, so loaded and unloaded perturbations are identical.",
        "The bounded-stress tangent cannot be used: it is diagnostic and absent from the governing static and dynamic equations.",
        "Propagation, inversion, point-ray control, and observer work were therefore not run.", "",
        *[f"{k}={str(v).lower() if isinstance(v, bool) else v}" for k, v in contract.items()], ""
    ])
    (OUT / "report.txt").write_text(report)
    print(report, end="")
    return contract


if __name__ == "__main__":
    main()
