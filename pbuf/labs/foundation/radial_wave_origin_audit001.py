#!/usr/bin/env python3
"""PBUF FOUNDATION — RADIAL WAVE ORIGIN AUDIT 001.

Fact-finding only.

Purpose
-------
Test whether the 1/r-type spatial amplitude behavior that was numerically
interesting in ``spatial_coherence_kernel_audit001.py`` follows generically
from ordinary radial propagation in three spatial dimensions, rather than being
kept as an arbitrary candidate kernel.

This lab does NOT claim that the PBUF medium obeys the scalar wave equation,
that electromagnetism is derived from a deeper substrate, or that gravity has
been derived.  It asks a narrower mathematical question:

    If an underlying isotropic medium supports a source-free wave equation,
    what radial amplitude law follows from geometry alone?

For an isotropic d-dimensional wave, conservation of outward quadratic flux
through a shell gives

    shell_area_factor * amplitude^2 = const
    A_d(r) proportional to r^{-(d-1)/2}.

Thus in d=3, A_3(r) proportional to 1/r.  The same 1/r factor also appears in
source-free radial Helmholtz/Laplace solutions outside a compact source.

The lab verifies those statements numerically/analytically, then replaces the
previous arbitrary ``inverse_1 = 1/(1+r/L)`` fact-finding family with a bounded
wave-origin kernel

    W_wave(r;L) = 1              for r <= L
                  L/r            for r > L

whose 1/r exterior is fixed by 3D radial spreading while the bounded core is
only a regularization convention.  It compares this wave-origin kernel to the
previous inverse_1 asymptotically and propagates the same microscopic candidate
amplitudes through the same target-blind mass-radius cases.

Important circularity audit
---------------------------
The supplied conventional Planck length is carried forward from earlier labs.
That numerical length conventionally contains Newton's G.  This lab therefore
computes, as an identity only,

    G_embedded(l_P) = c^3 l_P^2 / hbar

without loading an external value of G.  Any appearance of the familiar
GM/(R c^2) scale after substituting conventional l_P is explicitly labelled an
embedded-constant identity, NOT an independent derivation of G or gravity.

Hard guardrails
---------------
- no candidate selection/ranking;
- no fit/tuning;
- no kappa/shear/HST/lens benchmark input;
- no GR/Newtonian force/potential/deflection law;
- no external numerical G loaded;
- no Quantum Engine input;
- existing PBUF baryon/alpha structure untouched;
- legacy 0.18 is not used;
- K_P* remains dimensional-candidate-only;
- stdout only; no run directory.
"""
from __future__ import annotations

import json
import math
import subprocess
from dataclasses import dataclass
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
LAB_ID = "PBUF-FOUNDATION-RADIAL-WAVE-ORIGIN-AUDIT-001"

C_M_S = 299_792_458.0
H_J_S = 6.626_070_15e-34
HBAR_J_S = H_J_S / (2.0 * math.pi)
PLANCK_LENGTH_M = 1.616_255e-35
PROTON_MASS_KG = 1.672_621_925_95e-27
PROTON_CHARGE_RADIUS_M = 0.8409e-15
EARTH_MASS_KG = 5.9722e24
EARTH_RADIUS_M = 6.371e6
SUN_MASS_KG = 1.98847e30
SUN_RADIUS_M = 6.957e8
KPC_M = 3.085_677_581_491_367e19
MPC_M = 1.0e3 * KPC_M


@dataclass(frozen=True)
class MicroCandidate:
    key: str
    length_m: float
    volume_m3: float
    eps1: float


@dataclass(frozen=True)
class MassCase:
    key: str
    mass_kg: float
    radius_m: float


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


def _cube(x: float) -> float:
    return x**3


def _sphere(x: float) -> float:
    return 4.0 * math.pi * x**3 / 3.0


def _micro_candidates() -> tuple[list[MicroCandidate], dict]:
    lp = PLANCK_LENGTH_M
    vp = lp**3
    ep_star = HBAR_J_S * C_M_S / lp
    kp_star = ep_star / vp
    eb = PROTON_MASS_KG * C_M_S**2
    reduced = HBAR_J_S / (PROTON_MASS_KG * C_M_S)
    compton = H_J_S / (PROTON_MASS_KG * C_M_S)
    rp = PROTON_CHARGE_RADIUS_M

    specs = [
        ("planck_cell_cube", lp, _cube(lp)),
        ("reduced_compton_cube", reduced, _cube(reduced)),
        ("reduced_compton_sphere", reduced, _sphere(reduced)),
        ("compton_cube", compton, _cube(compton)),
        ("compton_sphere", compton, _sphere(compton)),
        ("proton_charge_radius_cube", rp, _cube(rp)),
        ("proton_charge_radius_sphere", rp, _sphere(rp)),
    ]
    rows = []
    for key, length, volume in specs:
        eps1 = (eb / volume) / kp_star
        rows.append(MicroCandidate(key, length, volume, eps1))

    g_embedded = C_M_S**3 * lp**2 / HBAR_J_S
    return rows, {
        "Kp_star_J_m3": kp_star,
        "Ep_star_J": ep_star,
        "planck_length_m": lp,
        "G_embedded_from_conventional_lP_m3_kg_s2": g_embedded,
        "G_embedded_role": "algebraic_identity_from_supplied_conventional_planck_length_not_external_input",
    }


def _mass_cases() -> list[MassCase]:
    return [
        MassCase("one_proton", PROTON_MASS_KG, PROTON_CHARGE_RADIUS_M),
        MassCase("earth", EARTH_MASS_KG, EARTH_RADIUS_M),
        MassCase("sun", SUN_MASS_KG, SUN_RADIUS_M),
        MassCase("baryonic_galaxy_scale", 6.0e10 * SUN_MASS_KG, 15.0 * KPC_M),
        MassCase("baryonic_cluster_scale", 1.0e14 * SUN_MASS_KG, 1.0 * MPC_M),
    ]


def _flux_amplitude(d: int, x: float) -> float:
    """Dimensionless radial amplitude normalized A(1)=1 from shell spreading."""
    return x ** (-(d - 1.0) / 2.0)


def _shell_flux_factor(d: int, x: float) -> float:
    a = _flux_amplitude(d, x)
    return x ** (d - 1) * a * a


def _wave_kernel_moments(q: float) -> tuple[float, float]:
    """Uniform 3D sphere moments for W=1 inside L and L/r outside L."""
    if q <= 0.0 or not math.isfinite(q):
        raise ValueError(q)
    if q <= 1.0:
        return 1.0, 1.0
    mean_w = (3.0 * q * q - 1.0) / (2.0 * q**3)
    mean_w2 = (3.0 * q - 2.0) / q**3
    return mean_w, mean_w2


def _inverse1_moments(q: float) -> tuple[float, float]:
    i1 = 0.5 * q*q - q + math.log1p(q)
    i2 = q + 1.0 - 2.0 * math.log1p(q) - 1.0/(1.0+q)
    return 3.0*i1/q**3, 3.0*i2/q**3


def _helmholtz_dimensionless_residual(x: float, kappa: float) -> float:
    """Analytic radial 3D residual for psi=cos(kappa*x)/x.

    Checks psi'' + (2/x) psi' + kappa^2 psi = 0 for x>0.
    """
    c = math.cos(kappa*x)
    s = math.sin(kappa*x)
    psi = c/x
    d1 = -kappa*s/x - c/x**2
    d2 = -kappa*kappa*c/x + 2.0*kappa*s/x**2 + 2.0*c/x**3
    return d2 + 2.0*d1/x + kappa*kappa*psi


def _wave_geometry_checks() -> dict:
    xs = (1.0, 2.0, 10.0, 1.0e3, 1.0e9)
    dimensions = {}
    max_flux_err = 0.0
    for d in (1, 2, 3, 4):
        vals = []
        for x in xs:
            amp = _flux_amplitude(d, x)
            flux = _shell_flux_factor(d, x)
            max_flux_err = max(max_flux_err, abs(flux - 1.0))
            vals.append({"r_over_reference": x, "amplitude": amp, "shell_quadratic_flux_factor": flux})
        dimensions[str(d)] = {
            "amplitude_power": -(d - 1.0)/2.0,
            "samples": vals,
        }

    max_helmholtz_abs = 0.0
    residual_rows = []
    for kappa in (0.0, 0.1, 1.0, 10.0):
        for x in (1.1, 2.0, 10.0, 100.0):
            r = _helmholtz_dimensionless_residual(x, kappa)
            max_helmholtz_abs = max(max_helmholtz_abs, abs(r))
            residual_rows.append({"kappa": kappa, "x": x, "residual": r})

    return {
        "dimension_flux_spreading": dimensions,
        "d3_amplitude_law": "A proportional to 1/r",
        "max_shell_flux_conservation_abs_error": max_flux_err,
        "helmholtz_3d_psi_cos_kr_over_r_samples": residual_rows,
        "max_abs_analytic_helmholtz_residual": max_helmholtz_abs,
    }


def _evaluate_mass_bridge(micro: list[MicroCandidate], reference: dict) -> list[dict]:
    out = []
    g_embedded = reference["G_embedded_from_conventional_lP_m3_kg_s2"]
    for m in _mass_cases():
        nb = m.mass_kg / PROTON_MASS_KG
        mr = []
        for mc in micro:
            q = m.radius_m / mc.length_m
            mw, mw2 = _wave_kernel_moments(q)
            old_mw, _ = _inverse1_moments(q)
            coherent = mc.eps1 * nb * mw
            random_rms = mc.eps1 * math.sqrt(nb * mw2)
            asymptotic = mc.eps1 * nb * (1.5 / q)
            mr.append({
                "micro_key": mc.key,
                "L_m": mc.length_m,
                "eps1": mc.eps1,
                "q_R_over_L": q,
                "wave_origin_meanW": mw,
                "wave_origin_meanW2": mw2,
                "coherent_wave_origin_strain": coherent,
                "random_phase_wave_origin_rms_strain": random_rms,
                "large_q_asymptotic_3_over_2q_strain": asymptotic,
                "relative_error_to_large_q_asymptotic": abs(coherent-asymptotic)/max(abs(coherent), 1e-300),
                "previous_inverse1_meanW": old_mw,
                "wave_meanW_over_previous_inverse1_meanW": mw/old_mw,
            })

        # This identity is reported only for the Planck-length lane and only to
        # expose what is already embedded in conventional l_P.
        chi_planck_asym = 1.5 * m.mass_kg * C_M_S * PLANCK_LENGTH_M**2 / (HBAR_J_S * m.radius_m)
        chi_embedded_g = 1.5 * g_embedded * m.mass_kg / (m.radius_m * C_M_S**2)
        out.append({
            "mass_key": m.key,
            "mass_kg": m.mass_kg,
            "radius_m": m.radius_m,
            "baryon_count": nb,
            "micro_results": mr,
            "planck_lane_asymptotic_direct_hbar_c_lP_form": chi_planck_asym,
            "same_quantity_rewritten_using_G_embedded_from_lP": chi_embedded_g,
            "embedded_identity_relative_error": abs(chi_planck_asym-chi_embedded_g)/max(abs(chi_planck_asym), 1e-300),
            "embedded_identity_status": "tautological_rewrite_not_independent_gravity_derivation",
        })
    return out


def _checks(wave: dict, masses: list[dict]) -> dict:
    finite = True
    bounded = True
    asymptotic_planck = True
    embedded_identity = True
    inverse_asymptotic = True
    for m in masses:
        embedded_identity &= m["embedded_identity_relative_error"] < 1e-12
        for r in m["micro_results"]:
            vals = [r["q_R_over_L"], r["wave_origin_meanW"], r["wave_origin_meanW2"],
                    r["coherent_wave_origin_strain"], r["random_phase_wave_origin_rms_strain"]]
            finite &= all(math.isfinite(v) for v in vals)
            bounded &= 0.0 <= r["wave_origin_meanW2"] <= r["wave_origin_meanW"] <= 1.0
            if r["q_R_over_L"] > 1e6:
                inverse_asymptotic &= abs(r["wave_meanW_over_previous_inverse1_meanW"] - 1.0) < 1e-4
                if r["micro_key"] == "planck_cell_cube":
                    asymptotic_planck &= r["relative_error_to_large_q_asymptotic"] < 1e-12

    return {
        "all_reported_numbers_finite": finite,
        "wave_kernel_moments_bounded": bounded,
        "shell_flux_conservation_pass": wave["max_shell_flux_conservation_abs_error"] < 1e-12,
        "three_spatial_dimensions_give_inverse_r_amplitude": wave["dimension_flux_spreading"]["3"]["amplitude_power"] == -1.0,
        "analytic_3d_radial_helmholtz_inverse_r_solution_pass": wave["max_abs_analytic_helmholtz_residual"] < 1e-12,
        "wave_kernel_converges_to_previous_inverse1_large_q": inverse_asymptotic,
        "planck_lane_exact_wave_kernel_converges_to_3_over_2q_large_q": asymptotic_planck,
        "embedded_G_rewrite_identity_pass": embedded_identity,
        "external_numerical_G_not_loaded": True,
        "embedded_G_not_claimed_as_derivation": True,
        "no_candidate_selected_or_ranked": True,
        "no_fit_or_tuning": True,
        "target_blind_no_kappa_shear_HST_or_lens_data": True,
        "no_GR_Newtonian_force_potential_or_deflection_law": True,
        "quantum_engine_not_used": True,
        "existing_pbuf_baryon_alpha_structure_untouched": True,
        "Kp_star_not_promoted_to_PBUF_modulus": True,
        "stdout_only_no_run_directory_created": True,
    }


def main() -> int:
    state = _repo_state()
    micro, reference = _micro_candidates()
    wave = _wave_geometry_checks()
    masses = _evaluate_mass_bridge(micro, reference)
    checks = _checks(wave, masses)

    payload = {
        "lab_id": LAB_ID,
        "status": "FACT_FINDING_ONLY",
        "question": "does inverse-r amplitude follow from generic 3D radial wave geometry and reproduce the prior target-blind bridge behavior",
        "repo_state": state,
        "reference": reference,
        "wave_geometry": wave,
        "mass_bridge": masses,
        "checks": checks,
        "guardrails": {
            "scalar_wave_equation_role": "mathematical_candidate_medium_equation_only",
            "bounded_core_role": "regularization_convention_only",
            "conventional_planck_length_contains_G": True,
            "G_embedded_result_role": "identity_exposing_circularity_not_new_derivation",
            "no_external_G": True,
            "no_observational_target": True,
        },
    }

    print(LAB_ID)
    print("status=FACT_FINDING_ONLY")
    print(f"head_sha={state['head_sha']}")
    print("no_candidate_selected=true")
    print("no_fit=true")
    print("target_blind=true")
    print()
    print("WAVE_GEOMETRY")
    for d in (1, 2, 3, 4):
        p = wave["dimension_flux_spreading"][str(d)]["amplitude_power"]
        print(f"dimension={d} radial_amplitude_power={p:+.6f}")
    print(f"d3_law={wave['d3_amplitude_law']}")
    print(f"max_shell_flux_error={wave['max_shell_flux_conservation_abs_error']:.6e}")
    print(f"max_helmholtz_residual={wave['max_abs_analytic_helmholtz_residual']:.6e}")
    print()
    print("PLANCK_LENGTH_CIRCULARITY_AUDIT")
    print(f"planck_length_m={PLANCK_LENGTH_M:.17e}")
    print(f"G_embedded_from_lP={reference['G_embedded_from_conventional_lP_m3_kg_s2']:.17e}")
    print("G_embedded_status=identity_from_conventional_lP_not_external_input_not_derivation")
    print()
    print("MASS_BRIDGE_WAVE_ORIGIN")
    for m in masses:
        print(f"MASS {m['mass_key']} mass_kg={m['mass_kg']:.17e} radius_m={m['radius_m']:.17e}")
        for r in m["micro_results"]:
            print(
                f"  {r['micro_key']} q={r['q_R_over_L']:.6e} "
                f"meanW={r['wave_origin_meanW']:.6e} "
                f"coherent={r['coherent_wave_origin_strain']:.17e} "
                f"random_rms={r['random_phase_wave_origin_rms_strain']:.17e} "
                f"wave_over_old_inverse1={r['wave_meanW_over_previous_inverse1_meanW']:.9e}"
            )
        print(
            "  PLANCK_ASYMPTOTIC "
            f"direct_hbar_c_lP={m['planck_lane_asymptotic_direct_hbar_c_lP_form']:.17e} "
            f"embedded_G_rewrite={m['same_quantity_rewritten_using_G_embedded_from_lP']:.17e} "
            f"relerr={m['embedded_identity_relative_error']:.3e}"
        )
    print()
    print("CHECKS")
    for k, v in checks.items():
        print(f"{k}={str(v).lower() if isinstance(v, bool) else v}")
    print()
    print("INTERPRETATION_GUARDRAIL")
    print("3D inverse-r amplitude is a geometric wave-spreading result under the tested scalar-wave hypothesis")
    print("this does not establish that the PBUF substrate obeys that equation")
    print("the conventional numerical Planck length already embeds G, so the GM/Rc^2 rewrite is circular")
    print("no observational target was used to choose the wave law or microscopic scale")
    print("JSON=" + json.dumps(payload, separators=(",", ":"), allow_nan=False))

    return 0 if all(bool(v) for v in checks.values()) else 2


if __name__ == "__main__":
    raise SystemExit(main())
