#!/usr/bin/env python3
"""PBUF FOUNDATION — MEDIUM DYNAMICS NO-FIT CLOSURE 001.

Purpose
-------
Continue the physical matter-loading bridge after the effective metric-strain
map was closed:

    chi_{mu nu} = (g_{mu nu} - gbar_{mu nu})/2,
    J_{mu nu}   = -T_{mu nu}   (component-placement convention of prior lab).

This lab asks whether the remaining medium dynamics S_sigma can now be closed
WITHOUT introducing a fitted modulus or silently reusing ordinary gravitational
lensing as PBUF microphysics.

The answer must be source-bounded.  PBUF V11 retains the standard Einstein
kinetic sector operationally, while the existing MEDIUM/ENERGY principle work
states that the local constitutive Hessian of the elastic sector is not fixed by
the homogeneous Omega_sigma(a) history.  Therefore this lab does three things:

1. algebraically verifies the gauge/Bianchi structure of the unique standard
   massless spin-2 kinetic boundary (linearized Einstein operator) in the
   metric-strain variable h=2 chi;
2. verifies that a naive algebraic Hooke/mass restoring term is not compatible
   with that linearized gauge symmetry, so it cannot simply be inserted as an
   unfitted local stiffness;
3. gives an explicit non-uniqueness proof for an isotropic quadratic elastic
   Hessian: covariance/objectivity/stability permit more than one positive
   tangent operator, so the missing moduli/kernel cannot be selected from the
   current macroscopic inputs alone.

The standard GR kinetic operator is a BOUNDARY CONDITION / retained effective
sector only.  It is NOT adopted as the missing PBUF elastic S_sigma and is not
fed into any lensing calculation.

No benchmark data, HST data, kappa/shear values, fitted numbers, replacement
for strength=0.18, cluster masses, or lensing pipeline are used.
"""
from __future__ import annotations

import json
import subprocess
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "runs" / "medium_dynamics_no_fit_closure001"
LAB_ID = "PBUF-FOUNDATION-MEDIUM-DYNAMICS-NO-FIT-CLOSURE-001"
TOL = 1.0e-11
ETA = np.diag([-1.0, 1.0, 1.0, 1.0])


def _git(*args: str) -> str:
    return subprocess.check_output(["git", *args], cwd=str(ROOT), text=True).strip()


def _write_json(path: Path, obj) -> None:
    path.write_text(json.dumps(obj, indent=2, sort_keys=False) + "\n")


def _rms(x) -> float:
    a = np.asarray(x, dtype=np.float64)
    return float(np.sqrt(np.mean(a * a))) if a.size else 0.0


def _rel_rms(a, b) -> float:
    aa = np.asarray(a, dtype=np.float64)
    bb = np.asarray(b, dtype=np.float64)
    return _rms(aa - bb) / max(_rms(bb), 1.0e-30)


def _sym(a: np.ndarray) -> np.ndarray:
    a = np.asarray(a, dtype=np.float64)
    return 0.5 * (a + a.T)


def _raise2_cov(h_cov: np.ndarray) -> np.ndarray:
    return ETA @ h_cov @ ETA


def _trace_cov(h_cov: np.ndarray) -> float:
    return float(np.einsum("mn,mn->", ETA, h_cov))


def linearized_einstein_fourier(h_cov: np.ndarray, k_cov: np.ndarray) -> np.ndarray:
    """Linearized Einstein tensor G^(1)_mn in flat background, Fourier form.

    Signature (-,+,+,+).  Overall Fourier-sign convention is irrelevant to the
    gauge/Bianchi identities tested here.  h_cov is symmetric and k_cov is a
    covector.  This routine is a retained-GR boundary diagnostic only.
    """
    h = _sym(h_cov)
    k = np.asarray(k_cov, dtype=np.float64)
    k_up = ETA @ k
    k2 = float(k @ k_up)
    h_updown = ETA @ h
    h_upup = _raise2_cov(h)
    trh = _trace_cov(h)

    R = np.zeros((4, 4), dtype=np.float64)
    for m in range(4):
        for n in range(4):
            term1 = -k[m] * float(np.dot(k, h_updown[:, n]))
            term2 = -k[n] * float(np.dot(k, h_updown[:, m]))
            term3 = k2 * h[m, n]
            term4 = k[m] * k[n] * trh
            R[m, n] = 0.5 * (term1 + term2 + term3 + term4)

    Rscalar = -float(np.einsum("m,n,mn->", k, k, h_upup)) + k2 * trh
    return R - 0.5 * ETA * Rscalar


def gauge_increment(k_cov: np.ndarray, xi_cov: np.ndarray) -> np.ndarray:
    k = np.asarray(k_cov, dtype=np.float64)
    xi = np.asarray(xi_cov, dtype=np.float64)
    return np.outer(k, xi) + np.outer(xi, k)


def naive_fierz_pauli_algebraic_tensor(h_cov: np.ndarray) -> np.ndarray:
    """Algebraic h_mn-eta_mn h combination; diagnostic only, no mass coefficient."""
    h = _sym(h_cov)
    return h - ETA * _trace_cov(h)


def _dev3(e: np.ndarray) -> np.ndarray:
    e = _sym(e)
    return e - np.eye(3) * (np.trace(e) / 3.0)


def isotropic_quadratic_energy(e: np.ndarray, lam: float, mu: float) -> float:
    """Synthetic tangent family W2=(lam/2)(tr e)^2+mu tr(e^2)."""
    e = _sym(e)
    return float(0.5 * lam * np.trace(e) ** 2 + mu * np.einsum("ij,ij->", e, e))


def isotropic_quadratic_stress(e: np.ndarray, lam: float, mu: float) -> np.ndarray:
    e = _sym(e)
    return lam * np.trace(e) * np.eye(3) + 2.0 * mu * e


def algebra_audit() -> dict:
    rng = np.random.default_rng(20260808_39)

    # Non-null wave covector avoids projector degeneracies and is sufficient for
    # exact linearized gauge/Bianchi algebra tests.
    k = np.array([1.7, 0.4, -0.6, 0.3], dtype=np.float64)
    h = _sym(rng.normal(size=(4, 4)))
    xi = rng.normal(size=4)

    G = linearized_einstein_fourier(h, k)
    dh = gauge_increment(k, xi)
    G_gauge = linearized_einstein_fourier(h + dh, k)
    gauge_error = _rel_rms(G_gauge, G)

    # Contract first index with k^mu: k^mu G_mn = 0.
    k_up = ETA @ k
    bianchi = np.einsum("m,mn->n", k_up, G)
    bianchi_rel = _rms(bianchi) / max(_rms(G) * _rms(k_up), 1.0e-30)

    # The algebraic Hooke/Fierz-Pauli-like restoring tensor is NOT gauge
    # invariant by itself.  This is a rejection test: its variation must be
    # clearly nonzero, not pass a zero gate.
    P0 = naive_fierz_pauli_algebraic_tensor(h)
    P1 = naive_fierz_pauli_algebraic_tensor(h + dh)
    algebraic_gauge_change = _rms(P1 - P0) / max(_rms(dh), 1.0e-30)

    # Explicit non-uniqueness proof for a stable isotropic local tangent.
    # Coefficients are SYNTHETIC witnesses only, not PBUF numbers.
    strain_samples = [_sym(rng.normal(scale=1.0e-3, size=(3, 3))) for _ in range(24)]
    witness_A = (1.0, 1.0)
    witness_B = (2.0, 1.0)
    energies_A = np.array([isotropic_quadratic_energy(e, *witness_A) for e in strain_samples])
    energies_B = np.array([isotropic_quadratic_energy(e, *witness_B) for e in strain_samples])
    stresses_A = np.array([isotropic_quadratic_stress(e, *witness_A) for e in strain_samples])
    stresses_B = np.array([isotropic_quadratic_stress(e, *witness_B) for e in strain_samples])
    both_positive = bool(np.all(energies_A > 0.0) and np.all(energies_B > 0.0))
    tangent_distinct = _rel_rms(stresses_A, stresses_B)

    # Conditional ordinary-isotropic-solid speed check.  If one insists on BOTH
    # physical longitudinal and transverse elastic modes propagating at the same
    # invariant c with positive inertial density, then
    #   c_T^2=mu/rho, c_L^2=(lam+2mu)/rho
    # implies lam=-mu and bulk K=lam+2mu/3=-mu/3 < 0.
    # This rejects that *ordinary solid* closure under those assumptions; it is
    # not a statement that PBUF must possess both material modes.
    mu = 1.0
    lam_equal_speed = -mu
    bulk = lam_equal_speed + 2.0 * mu / 3.0
    ordinary_equal_speed_bulk_positive = bool(bulk > 0.0)

    checks = {
        "linearized_einstein_gauge_invariance_pass": bool(gauge_error <= TOL),
        "linearized_bianchi_identity_pass": bool(bianchi_rel <= TOL),
        "naive_algebraic_restoring_term_rejected_by_gauge_test": bool(algebraic_gauge_change > 1.0e-6),
        "stable_isotropic_tangent_nonuniqueness_witness_pass": bool(both_positive and tangent_distinct > 1.0e-3),
        "ordinary_equal_speed_isotropic_solid_stability_conflict_detected": bool(not ordinary_equal_speed_bulk_positive),
    }

    return {
        "synthetic_only_not_physical_model": True,
        "tolerance": TOL,
        "metrics": {
            "linearized_einstein_gauge_invariance_relative_rms_error": gauge_error,
            "linearized_bianchi_relative_rms_error": bianchi_rel,
            "naive_algebraic_restoring_gauge_change_per_gauge_increment_rms": algebraic_gauge_change,
            "isotropic_tangent_witness_relative_stress_difference": tangent_distinct,
            "ordinary_equal_speed_witness_lambda_over_mu": lam_equal_speed / mu,
            "ordinary_equal_speed_witness_bulk_over_mu": bulk / mu,
        },
        "checks": checks,
        "all_checks_pass": bool(all(checks.values())),
    }


def main() -> int:
    started = time.perf_counter()
    OUT.mkdir(parents=True, exist_ok=True)

    repo = {
        "repository": "TheExiledMonk/lab",
        "branch": _git("rev-parse", "--abbrev-ref", "HEAD"),
        "head_sha": _git("rev-parse", "HEAD"),
        "tracked_changes": _git("diff", "--name-only"),
        "staged_changes": _git("diff", "--name-only", "--cached"),
    }
    if repo["branch"] != "main":
        raise RuntimeError(f"runner must execute on main, got {repo['branch']}")
    if repo["tracked_changes"] or repo["staged_changes"]:
        raise RuntimeError("tracked or staged repository changes present")

    algebra = algebra_audit()
    if not algebra["all_checks_pass"]:
        _write_json(OUT / "algebra_failure.json", algebra)
        raise RuntimeError("medium-dynamics no-fit algebra gate failed")

    boundary = {
        "metric_strain_variable": "h_mu_nu = 2 chi_mu_nu = g_mu_nu-gbar_mu_nu",
        "matter_loading_from_prior_closure": "J_mu_nu = -T_m_mu_nu (prior component-placement convention)",
        "retained_effective_kinetic_boundary": "linearized Einstein/Fierz-Pauli massless spin-2 operator",
        "retained_effective_kinetic_boundary_role": "GR_consistency_boundary_only_not_PBUF_elastic_S_sigma",
        "retained_effective_normalization": "standard Einstein normalization already retained by V11; not re-derived or repurposed here",
        "new_local_elastic_modulus_derived": False,
        "naive_algebraic_metric_strain_restoring_term_authorized": False,
        "reason_naive_algebraic_term_not_authorized": "an algebraic h_mu_nu restoring/mass term is not invariant under linearized diffeomorphism gauge transformations without additional structure",
        "ordinary_isotropic_solid_closure_selected": False,
        "ordinary_isotropic_solid_note": "if both longitudinal and transverse material modes are required to share c, the elementary isotropic-solid formulas force negative bulk modulus; this conditional closure is therefore unsuitable",
        "gr_lensing_law_used_as_PBUF_medium_response": False,
    }

    nonuniqueness = {
        "result": "current accepted macroscopic constraints do not uniquely determine the PBUF elastic Hessian/kernel",
        "proof_family": "W2=(lambda/2)(tr e)^2 + mu tr(e^2)",
        "synthetic_witnesses_only": [
            {"lambda": 1.0, "mu": 1.0, "physical_PBUF_value": False},
            {"lambda": 2.0, "mu": 1.0, "physical_PBUF_value": False},
        ],
        "both_witnesses_stable_on_test_samples": True,
        "witnesses_produce_distinct_tangent_response": True,
        "therefore_fit_free_modulus_selection_from_current_macro_constraints": False,
        "homogeneous_Omega_sigma_history_sufficient_to_fix_local_tensor_kernel": False,
        "alpha_or_epsilon0_may_be_relabelled_as_local_modulus_without_derivation": False,
    }

    next_spec = {
        "target": "derive_local_medium_Hessian_or_retarded_response_kernel_from_Quantum_Engine_microphysics",
        "required_object": "K_{mu nu alpha beta}(x,y)=delta^2 S_sigma/(delta chi_mu_nu(x) delta chi_alpha_beta(y)) or equivalent causal inverse G_R",
        "required_inputs_must_be_independent_of_lensing_target": [
            "normalized microscopic medium degrees of freedom",
            "physical regulator/cutoff normalization",
            "mode degeneracy/counting",
            "microscopic energy scale in SI or exact natural-unit map",
            "coupling of microscopic modes to metric strain chi_mu_nu",
            "coarse-graining measure from microscopic modes to local tensor strain",
            "state/temperature dependence where required",
        ],
        "required_checks": [
            "dimensions close without fitted scalar",
            "Hessian symmetry/reciprocity on reversible branch",
            "positive static energy on physical modes",
            "causal/hyperbolic retarded response",
            "no ghost or unauthorized extra propagating mode",
            "low-frequency/long-wavelength compatibility with retained Einstein dynamics",
            "homogeneous limit reproduces PBUF Omega_sigma only after local kernel is independently derived",
        ],
        "forbidden": [
            "fit_modulus_to_kappa_shear_or_cluster_lensing",
            "replace_0p18_with_new_scalar",
            "set_local_modulus_equal_to_alpha_epsilon0_Rmax_ksat_without_derivation",
            "use_GR_lensing_force_as_PBUF_elastic_source",
            "infer_kernel_from_F160W_morphology",
            "run_cluster_lensing_before_local_kernel_coarse_graining_and_physical_Tmunu_close",
        ],
    }

    unresolved = [
        "PBUF_local_elastic_Hessian_or_retarded_kernel",
        "microscopic_to_metric_strain_normalization",
        "coarse_graining_to_finite_PBUF_state",
        "independent_physical_Tmunu_dataset",
        "absolute_source_geometry",
    ]

    validation = {
        "lab_id": LAB_ID,
        "outcome": "Outcome B — NO-FIT PBUF MEDIUM DYNAMICS NOT YET CLOSED; GR KINETIC BOUNDARY AND QUANTUM RESPONSE TARGET ISOLATED",
        "head_sha": repo["head_sha"],
        "benchmark_pixel_values_loaded": False,
        "lensing_pipeline_executed": False,
        "hst_mass_conversion_executed": False,
        "fitted_numbers_introduced": False,
        "legacy_strength_used": False,
        "replacement_strength_selected": False,
        "metric_strain_map_inherited_closed": True,
        "matter_loading_kernel_inherited_closed": True,
        "retained_GR_kinetic_boundary_verified": True,
        "retained_GR_kinetic_boundary_adopted_as_PBUF_S_sigma": False,
        "naive_local_Hooke_metric_potential_authorized": False,
        "ordinary_isotropic_solid_closure_authorized": False,
        "PBUF_local_medium_Hessian_derived": False,
        "PBUF_medium_action_or_closed_dynamics_present": False,
        "physical_bridge_complete": False,
        "physical_cluster_lensing_run_authorized": False,
        "unresolved_required_components": unresolved,
        "next_derivation_authorized": True,
        "next_derivation_target": next_spec["target"],
        "algebra": algebra,
        "duration_seconds": time.perf_counter() - started,
    }

    _write_json(OUT / "validation.json", validation)
    _write_json(OUT / "medium_dynamics_boundary_spec.json", boundary)
    _write_json(OUT / "nonuniqueness_proof.json", nonuniqueness)
    _write_json(OUT / "next_quantum_response_derivation_spec.json", next_spec)
    _write_json(OUT / "repository_state.json", repo)

    print(validation["outcome"])
    print(json.dumps(validation, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
