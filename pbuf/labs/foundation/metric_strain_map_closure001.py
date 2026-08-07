#!/usr/bin/env python3
"""PBUF FOUNDATION — METRIC-STRAIN MAP CLOSURE 001.

No-fit analytic closure step following matter-loading physical closure 001.

Purpose
-------
Close ONLY the metric-medium map and its matter-loading kernel using the most
conservative no-extra-coupling identification available when the PBUF medium is
interpreted as spacetime geometry itself:

    chi_{mu nu} := 1/2 (g_{mu nu} - gbar_{mu nu})
    g_{mu nu}   := gbar_{mu nu} + 2 chi_{mu nu}

This is the standard metric-strain normalization (half the metric change), not
an observational fit and not a new matter coupling constant. It is a variable
identification/normalization. With minimal matter coupling,

    J_{alpha beta} = -1/2 T^{mu nu} d g_{mu nu}/d chi_{alpha beta}
                   = -T^{alpha beta}

for symmetric tensors in the local algebraic map.

What this DOES NOT do
---------------------
- It does not derive the PBUF medium action S_sigma.
- It does not derive local stiffness, gradient stiffness, dispersion, memory,
  or nonlinear constitutive terms.
- It does not reduce T^{mu nu} to rho.
- It does not load HST data or lensing benchmarks.
- It does not authorize a cluster lensing run.
- It does not claim chi is a new microscopic degree of freedom. In this branch
  chi is simply the normalized metric strain variable, i.e. a no-extra-DOF
  effective identification.

The lab algebraically verifies invertibility, tensor symmetry, exact response
kernel, exact matter-source contraction, basis covariance under synthetic
Lorentz transformations, and the absence of a free matter coupling scalar.
It then emits the remaining next target: derive/close S_sigma (or equivalent
medium dynamics) in this normalization without fitted numbers.
"""
from __future__ import annotations

import json
import subprocess
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "runs" / "metric_strain_map_closure001"
LAB_ID = "PBUF-FOUNDATION-METRIC-STRAIN-MAP-CLOSURE-001"
TOL = 1.0e-12


def _git(*args: str) -> str:
    return subprocess.check_output(["git", *args], cwd=str(ROOT), text=True).strip()


def _write_json(path: Path, obj) -> None:
    path.write_text(json.dumps(obj, indent=2, sort_keys=False) + "\n")


def _rms(x) -> float:
    a = np.asarray(x, dtype=np.float64)
    return float(np.sqrt(np.mean(a * a))) if a.size else 0.0


def _relative_rms(a, b) -> float:
    aa = np.asarray(a, dtype=np.float64)
    bb = np.asarray(b, dtype=np.float64)
    return _rms(aa - bb) / max(_rms(bb), 1.0e-30)


def _sym(a: np.ndarray) -> np.ndarray:
    return 0.5 * (a + a.T)


def metric_from_strain(gbar: np.ndarray, chi: np.ndarray) -> np.ndarray:
    """Exact defining map g = gbar + 2 chi."""
    return np.asarray(gbar, dtype=np.float64) + 2.0 * _sym(np.asarray(chi, dtype=np.float64))


def strain_from_metric(gbar: np.ndarray, g: np.ndarray) -> np.ndarray:
    """Exact inverse chi = (g-gbar)/2."""
    return 0.5 * _sym(np.asarray(g, dtype=np.float64) - np.asarray(gbar, dtype=np.float64))


def response_kernel() -> np.ndarray:
    """R[mu,nu,a,b] = d g_mn / d chi_ab for symmetric chi coordinates.

    The symmetric-coordinate derivative of g_mn=gbar_mn+2 chi_mn is

      R_mn,ab = delta_ma delta_nb + delta_mb delta_na.

    Contracted with a symmetric variation dchi_ab this gives 2 dchi_mn.
    """
    d = 4
    R = np.zeros((d, d, d, d), dtype=np.float64)
    eye = np.eye(d)
    for m in range(d):
        for n in range(d):
            for a in range(d):
                for b in range(d):
                    R[m, n, a, b] = eye[m, a] * eye[n, b] + eye[m, b] * eye[n, a]
    return R


def matter_source(Tup: np.ndarray, R: np.ndarray) -> np.ndarray:
    """J_ab = -1/2 T^mn R_mn,ab."""
    return -0.5 * np.einsum("mn,mnab->ab", np.asarray(Tup, dtype=np.float64), R)


def _lorentz_boost_x(beta: float) -> np.ndarray:
    if not (abs(beta) < 1.0):
        raise ValueError("|beta| must be <1")
    gamma = 1.0 / np.sqrt(1.0 - beta * beta)
    L = np.eye(4)
    L[0, 0] = gamma
    L[0, 1] = -gamma * beta
    L[1, 0] = -gamma * beta
    L[1, 1] = gamma
    return L


def _covariant_transform_cov2(A: np.ndarray, L: np.ndarray) -> np.ndarray:
    # Synthetic passive-coordinate tensor test: A' = L^{-T} A L^{-1}.
    Li = np.linalg.inv(L)
    return Li.T @ A @ Li


def _contravariant_transform_rank2(A: np.ndarray, L: np.ndarray) -> np.ndarray:
    # A'^{mu nu}=L^mu_a L^nu_b A^{ab}.
    return L @ A @ L.T


def algebra_audit() -> dict:
    rng = np.random.default_rng(20260808)
    eta = np.diag([-1.0, 1.0, 1.0, 1.0])

    chi = _sym(rng.normal(scale=1.0e-3, size=(4, 4)))
    g = metric_from_strain(eta, chi)
    chi_back = strain_from_metric(eta, g)
    invertibility_error = _relative_rms(chi_back, chi)

    dchi = _sym(rng.normal(scale=1.0e-7, size=(4, 4)))
    eps = 1.0e-4
    dg_fd = (metric_from_strain(eta, chi + eps * dchi) - metric_from_strain(eta, chi - eps * dchi)) / (2.0 * eps)
    dg_exact = 2.0 * dchi
    finite_difference_error = _relative_rms(dg_fd, dg_exact)

    R = response_kernel()
    dg_kernel = np.einsum("mnab,ab->mn", R, dchi)
    kernel_action_error = _relative_rms(dg_kernel, dg_exact)

    T = _sym(rng.normal(size=(4, 4)))
    J = matter_source(T, R)
    source_identity_error = _relative_rms(J, -T)

    scale_T = 3.25
    scale_chi = -1.75
    source_T_linearity = _relative_rms(matter_source(scale_T * T, R), scale_T * J)
    g_chi_linearity = _relative_rms(
        metric_from_strain(eta, scale_chi * chi) - eta,
        scale_chi * (metric_from_strain(eta, chi) - eta),
    )

    # Check action variation identity delta S_m density = (1/2) T^mn delta g_mn
    # and -J_ab delta chi_ab are identical with the adopted sign convention.
    lhs = 0.5 * float(np.einsum("mn,mn->", T, dg_exact))
    rhs = -float(np.einsum("ab,ab->", J, dchi))
    variation_identity_relative_error = abs(lhs - rhs) / max(abs(lhs), abs(rhs), 1.0e-30)

    # Synthetic covariance check. A background metric is transformed together with g,
    # so chi remains a covariant rank-2 tensor. T is transformed contravariantly.
    L = _lorentz_boost_x(0.37)
    eta_p = _covariant_transform_cov2(eta, L)
    g_p = _covariant_transform_cov2(g, L)
    chi_p_from_def = strain_from_metric(eta_p, g_p)
    chi_p_tensor = _covariant_transform_cov2(chi, L)
    chi_covariance_error = _relative_rms(chi_p_from_def, chi_p_tensor)

    T_p = _contravariant_transform_rank2(T, L)
    J_p_direct = -T_p
    # Since J carries the same component placement as T in this contraction convention,
    # compare to transformed J as contravariant rank-2 components.
    J_p_tensor = _contravariant_transform_rank2(J, L)
    source_covariance_error = _relative_rms(J_p_direct, J_p_tensor)

    zero_source_norm = _rms(matter_source(np.zeros((4, 4)), R))
    symmetry_error_R = _relative_rms(R, np.swapaxes(R, 0, 1)) + _relative_rms(R, np.swapaxes(R, 2, 3))

    checks = {
        "metric_strain_inverse_pass": invertibility_error <= TOL,
        "finite_difference_response_pass": finite_difference_error <= TOL,
        "response_kernel_action_pass": kernel_action_error <= TOL,
        "matter_source_identity_J_equals_minus_T_pass": source_identity_error <= TOL,
        "linearity_in_stress_energy_pass": source_T_linearity <= TOL,
        "linearity_in_metric_strain_pass": g_chi_linearity <= TOL,
        "matter_action_variation_identity_pass": variation_identity_relative_error <= TOL,
        "metric_strain_tensor_covariance_pass": chi_covariance_error <= TOL,
        "matter_source_tensor_covariance_pass": source_covariance_error <= TOL,
        "zero_stress_energy_zero_source_pass": zero_source_norm <= TOL,
        "response_kernel_pair_symmetry_pass": symmetry_error_R <= TOL,
    }
    return {
        "synthetic_only_not_physical_source": True,
        "tolerance": TOL,
        "metrics": {
            "metric_strain_inverse_relative_rms_error": invertibility_error,
            "finite_difference_response_relative_rms_error": finite_difference_error,
            "response_kernel_action_relative_rms_error": kernel_action_error,
            "matter_source_J_equals_minus_T_relative_rms_error": source_identity_error,
            "linearity_in_stress_energy_relative_rms_error": source_T_linearity,
            "linearity_in_metric_strain_relative_rms_error": g_chi_linearity,
            "matter_action_variation_identity_relative_error": variation_identity_relative_error,
            "metric_strain_tensor_covariance_relative_rms_error": chi_covariance_error,
            "matter_source_tensor_covariance_relative_rms_error": source_covariance_error,
            "zero_stress_energy_source_rms": zero_source_norm,
            "response_kernel_pair_symmetry_error": symmetry_error_R,
        },
        "checks": checks,
        "all_checks_pass": bool(all(checks.values())),
    }


def main() -> int:
    started = time.perf_counter()
    OUT.mkdir(parents=True, exist_ok=True)

    algebra = algebra_audit()
    if not algebra["all_checks_pass"]:
        _write_json(OUT / "algebra_failure.json", algebra)
        raise RuntimeError("metric-strain map algebra gate failed")

    metric_map_spec = {
        "status": "closed_as_minimal_no_extra_DOF_effective_identification",
        "interpretive_scope": "effective_metric_strain_variable_not_claimed_microscopic_ontology",
        "background_structure": "gbar_mu_nu must be specified and transformed covariantly with g_mu_nu",
        "medium_variable": "chi_mu_nu = (g_mu_nu - gbar_mu_nu)/2",
        "inverse_map": "g_mu_nu = gbar_mu_nu + 2 chi_mu_nu",
        "chi_dimensions": "dimensionless_metric_strain",
        "normalization_origin": "standard_half_metric_change_strain_definition_not_fit",
        "free_metric_medium_coupling_constant": False,
        "local_response_kernel": "R_mu_nu,alpha_beta = delta_mu_alpha delta_nu_beta + delta_mu_beta delta_nu_alpha",
        "minimal_matter_source": "J_alpha_beta = -T_m^{alpha beta}",
        "source_normalization_fixed_by_action_chain_rule": True,
        "rho_only_reduction_authorized": False,
        "new_propagating_degree_of_freedom_introduced": False,
        "uniqueness_claim": "unique only after choosing metric strain itself as the effective medium coordinate; does not prove this is PBUF microscopic ontology",
        "nonlinear_extension_note": "finite-strain/nonlinear metric parametrizations may be required outside the perturbative/effective regime and must be derived rather than fit",
    }

    remaining = {
        "metric_medium_map_closed_for_effective_metric_strain_branch": True,
        "metric_response_kernel_closed_for_effective_metric_strain_branch": True,
        "matter_loading_normalization_closed_for_effective_metric_strain_branch": True,
        "medium_action_or_closed_dynamics_present": False,
        "coarse_graining_to_existing_finite_PBUF_state_closed": False,
        "independent_physical_Tmunu_dataset_present": False,
        "absolute_source_geometry_present": False,
        "physical_bridge_complete": False,
        "next_derivation_authorized": True,
        "next_derivation_target": "derive_or_close_S_sigma_medium_dynamics_in_metric_strain_normalization_without_free_or_fitted_coefficients",
    }

    next_action_spec = {
        "target": remaining["next_derivation_target"],
        "required_equation": "delta S_sigma[gbar,chi;fixed_PBUF_microphysics]/delta chi_mu_nu = sqrt(-g) J^{mu nu}",
        "matter_source_fixed_by_this_lab": "J^{mu nu} = -T_m^{mu nu}",
        "must_derive": [
            "quadratic/static restoring operator or demonstrate why a different operator is required",
            "gradient/propagation operator and causal dynamical completion",
            "normalization/units of the finite laboratory state relative to chi_mu_nu",
            "stability and causal/hyperbolic conditions",
            "long_wavelength matching condition consistent with PBUF retained effective spacetime dynamics",
        ],
        "forbidden": [
            "fit_coefficients_to_kappa_or_shear",
            "replace_0p18_with_another_scalar",
            "identify_alpha_epsilon0_Rmax_ksat_or_1_over_137_as_local_modulus_without_derivation",
            "assume_F160W_mass_to_light_ratio",
            "inject_Newtonian_or_GR_lensing_force_into_PBUF_medium_response",
            "run_cluster_lensing_before_medium_dynamics_and_source_geometry_close",
        ],
    }

    repo = {
        "repository": "TheExiledMonk/lab",
        "branch": _git("rev-parse", "--abbrev-ref", "HEAD"),
        "head_sha": _git("rev-parse", "HEAD"),
        "tracked_changes": _git("diff", "--name-only"),
        "staged_changes": _git("diff", "--name-only", "--cached"),
    }

    validation = {
        "lab_id": LAB_ID,
        "outcome": "Outcome A — MINIMAL EFFECTIVE METRIC-STRAIN MAP AND MATTER-LOADING KERNEL CLOSED; MEDIUM DYNAMICS STILL OPEN",
        "head_sha": repo["head_sha"],
        "benchmark_pixel_values_loaded": False,
        "lensing_pipeline_executed": False,
        "hst_mass_conversion_executed": False,
        "fitted_numbers_introduced": False,
        "legacy_strength_used": False,
        "replacement_strength_selected": False,
        "metric_medium_map_branch_selected": "effective_metric_strain_no_extra_DOF",
        "metric_medium_map_is_observational_fit": False,
        "metric_medium_map_is_microscopic_PBUF_derivation": False,
        "metric_medium_map_is_effective_variable_identification": True,
        "metric_medium_map_closed_for_selected_branch": True,
        "metric_response_kernel_closed_for_selected_branch": True,
        "matter_loading_kernel_closed_for_selected_branch": True,
        "independent_matter_coupling_constant_allowed": False,
        "minimal_matter_source": "J_alpha_beta=-T_m^{alpha beta}",
        "rho_only_source_authorized": False,
        "medium_action_or_closed_dynamics_present": False,
        "physical_bridge_complete": False,
        "physical_cluster_lensing_run_authorized": False,
        "next_derivation_authorized": True,
        "next_derivation_target": remaining["next_derivation_target"],
        "algebra": algebra,
        "duration_seconds": time.perf_counter() - started,
    }

    _write_json(OUT / "validation.json", validation)
    _write_json(OUT / "metric_medium_map_spec.json", metric_map_spec)
    _write_json(OUT / "remaining_bridge_status.json", remaining)
    _write_json(OUT / "next_medium_action_derivation_spec.json", next_action_spec)
    _write_json(OUT / "repository_state.json", repo)

    print(validation["outcome"])
    print(json.dumps(validation, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
