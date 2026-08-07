#!/usr/bin/env python3
"""PBUF FOUNDATION — SPATIAL / COHERENCE KERNEL AUDIT 001.

Fact-finding only.

Purpose
-------
Continue directly from ``collective_mass_scaling_audit001.py`` by probing the
missing spatial/coherence bridge between a microscopic one-baryon candidate
response and a finite macroscopic source.

The previous audit established two deliberately extreme lanes:

1. every baryon contributes its full microscopic amplitude coherently;
2. the same mass is smeared into a macroscopic mean energy density and divided
   directly by the dimensional Planck-density candidate K_P*.

Those lanes differ by enormous factors.  This lab does NOT choose a physical
law between them.  Instead it evaluates several dimensionless radial weighting
families, normalized to W(0)=1, and two accumulation semantics:

    coherent:     chi = eps_1 * N_b * <W>
    random-phase: chi_rms = eps_1 * sqrt(N_b * <W^2>)

where the averages are over a uniform spherical source and eps_1 is one of the
same seven microscopic candidate amplitudes used by the prior two labs.

This is mathematical fact-finding, not a PBUF constitutive derivation.  The
kernel families are NOT gravitational laws and are not calibrated to any lens,
kappa, shear, HST map, GR prediction, Newtonian force law, or legacy strength.
They are deliberately simple interpolation families used only to expose which
spatial/coherence scalings are numerically stable, explosive, or negligible.

Source geometry
---------------
For a uniform sphere of radius R observed at its center, define q=R/L, where L
is the microscopic candidate length.  Then

    <W>   = 3/q^3 integral_0^q y^2 W(y) dy
    <W^2> = 3/q^3 integral_0^q y^2 W(y)^2 dy.

Closed-form expressions are used for numerical stability over q spanning more
than 50 orders of magnitude.  No Monte Carlo sampling and no float64 addition
of tiny strains to an O(1) metric are performed.

Kernel families, all with W(0)=1
--------------------------------
- coherent_upper: W(y)=1.  This exactly recovers the previous unsmeared upper
  bookkeeping lane and is included only as a reference boundary.
- compact_unit: W(y)=1 for y<=1 and 0 otherwise.
- exponential: W(y)=exp(-y).
- gaussian: W(y)=exp(-y^2/2).
- inverse_1: W(y)=1/(1+y).
- inverse_2: W(y)=1/(1+y)^2.

The inverse-power-looking families are dimensionless mathematical kernels in
q=r/L.  Their inclusion does NOT import Newtonian/GR gravity; no G, lensing
force, potential, deflection, or observed target enters this lab.

Mass cases
----------
- one proton, using the measured proton charge-radius scale as source radius;
- Earth;
- Sun;
- representative baryonic galaxy: 6e10 M_sun within 15 kpc;
- representative baryonic cluster: 1e14 M_sun within 1 Mpc.

The galaxy/cluster cases are explicit order-of-magnitude synthetic fact-finding
sources, not measurements from a lens benchmark.

Hard guardrails
---------------
- no candidate selected or ranked;
- no fit/tuning;
- no kappa/shear/HST/lens benchmark input;
- no GR/Newtonian calibration;
- no Quantum Engine input;
- existing PBUF baryon/alpha structure is untouched;
- legacy 0.18 is comparison-only and never a target;
- K_P* remains a dimensional candidate only, not a promoted PBUF modulus;
- output is stdout only; no run directory is created.
"""
from __future__ import annotations

import json
import math
import subprocess
from dataclasses import dataclass
from decimal import Decimal, getcontext
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
LAB_ID = "PBUF-FOUNDATION-SPATIAL-COHERENCE-KERNEL-AUDIT-001"

# Exact SI constants where fixed by SI.
C_M_S = 299_792_458.0
H_J_S = 6.626_070_15e-34
HBAR_J_S = H_J_S / (2.0 * math.pi)

# Reference quantities carried forward unchanged from the prior audits.
PLANCK_LENGTH_M = 1.616_255e-35
PROTON_MASS_KG = 1.672_621_925_95e-27
PROTON_CHARGE_RADIUS_M = 0.8409e-15
EARTH_MASS_KG = 5.9722e24
EARTH_RADIUS_M = 6.371e6
SUN_MASS_KG = 1.98847e30
SUN_RADIUS_M = 6.957e8
KPC_M = 3.085_677_581_491_367e19
MPC_M = 1.0e3 * KPC_M
LEGACY_STRENGTH_REFERENCE = 0.18

getcontext().prec = 80


@dataclass(frozen=True)
class MicroCandidate:
    key: str
    length_m: float
    volume_m3: float
    eps1: float
    role: str


@dataclass(frozen=True)
class MassCase:
    key: str
    mass_kg: float
    radius_m: float
    role: str
    status: str


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


def _cube(l: float) -> float:
    return l**3


def _sphere(r: float) -> float:
    return (4.0 / 3.0) * math.pi * r**3


def _micro_candidates() -> tuple[list[MicroCandidate], dict]:
    lp = PLANCK_LENGTH_M
    vp = lp**3
    ep_star = HBAR_J_S * C_M_S / lp
    kp_star = ep_star / vp
    eb = PROTON_MASS_KG * C_M_S**2

    reduced = HBAR_J_S / (PROTON_MASS_KG * C_M_S)
    ordinary = H_J_S / (PROTON_MASS_KG * C_M_S)
    rp = PROTON_CHARGE_RADIUS_M

    specs = [
        ("planck_cell_cube", lp, _cube(lp), "elementary-medium-cell boundary hypothesis"),
        ("reduced_compton_cube", reduced, _cube(reduced), "quantum source-extent candidate"),
        ("reduced_compton_sphere", reduced, _sphere(reduced), "quantum source-extent geometry check"),
        ("compton_cube", ordinary, _cube(ordinary), "ordinary-Compton source-extent candidate"),
        ("compton_sphere", ordinary, _sphere(ordinary), "ordinary-Compton geometry check"),
        ("proton_charge_radius_cube", rp, _cube(rp), "measured charge-size response-volume candidate"),
        ("proton_charge_radius_sphere", rp, _sphere(rp), "measured charge-size geometry check"),
    ]
    out = []
    for key, length, volume, role in specs:
        eps1 = (eb / volume) / kp_star
        out.append(MicroCandidate(key, length, volume, eps1, role))

    return out, {
        "planck_volume_m3": vp,
        "Ep_star_J": ep_star,
        "Kp_star_J_m3": kp_star,
        "Kp_star_status": "dimensional_candidate_only_not_derived_PBUF_modulus",
        "proton_rest_energy_J": eb,
        "planck_mass_equivalent_kg_hbar_over_c_lP": HBAR_J_S / (C_M_S * lp),
    }


def _mass_cases() -> list[MassCase]:
    return [
        MassCase("one_proton", PROTON_MASS_KG, PROTON_CHARGE_RADIUS_M,
                 "single-baryon spatial reference", "measured_reference_scale"),
        MassCase("earth", EARTH_MASS_KG, EARTH_RADIUS_M,
                 "planetary-scale distributed source", "measured_reference_scale"),
        MassCase("sun", SUN_MASS_KG, SUN_RADIUS_M,
                 "stellar-scale distributed source", "measured_reference_scale"),
        MassCase("baryonic_galaxy_scale", 6.0e10 * SUN_MASS_KG, 15.0 * KPC_M,
                 "representative baryonic galaxy scale", "synthetic_order_of_magnitude_case"),
        MassCase("baryonic_cluster_scale", 1.0e14 * SUN_MASS_KG, 1.0 * MPC_M,
                 "representative baryonic cluster scale", "synthetic_order_of_magnitude_case"),
    ]


def _stable_one_minus_exp_poly(q: float, a: float) -> float:
    """Return 1-exp(-a q)*(1+a q+(a q)^2/2) robustly enough for our q range."""
    z = a * q
    if z < 1.0e-3:
        # Series starts z^3/6 - z^4/8 + z^5/20 - ...
        return z**3 / 6.0 - z**4 / 8.0 + z**5 / 20.0 - z**6 / 72.0
    if z > 745.0:
        return 1.0
    return 1.0 - math.exp(-z) * (1.0 + z + 0.5 * z * z)


def _moments(kernel: str, q: float) -> tuple[float, float]:
    """Return (<W>, <W^2>) for a uniform sphere observed at its center."""
    if not (q > 0.0 and math.isfinite(q)):
        raise ValueError(f"q must be finite and >0, got {q}")

    q3 = q**3

    if kernel == "coherent_upper":
        return 1.0, 1.0

    if kernel == "compact_unit":
        if q <= 1.0:
            return 1.0, 1.0
        v = 1.0 / q3
        return v, v

    if kernel == "exponential":
        # Integral y^2 exp(-a y) dy from 0..q = 2/a^3 * bracket.
        i1 = 2.0 * _stable_one_minus_exp_poly(q, 1.0)
        i2 = 0.25 * _stable_one_minus_exp_poly(q, 2.0)
        return 3.0 * i1 / q3, 3.0 * i2 / q3

    if kernel == "gaussian":
        # W=exp(-y^2/2), W^2=exp(-y^2).
        if q > 27.0:
            i1 = math.sqrt(math.pi / 2.0)
            i2 = math.sqrt(math.pi) / 4.0
        else:
            i1 = math.sqrt(math.pi / 2.0) * math.erf(q / math.sqrt(2.0)) - q * math.exp(-0.5 * q * q)
            i2 = math.sqrt(math.pi) / 4.0 * math.erf(q) - 0.5 * q * math.exp(-q * q)
        return 3.0 * i1 / q3, 3.0 * i2 / q3

    if kernel == "inverse_1":
        # I1=int y^2/(1+y) dy
        i1 = 0.5 * q * q - q + math.log1p(q)
        # I2=int y^2/(1+y)^2 dy
        i2 = q + 1.0 - 2.0 * math.log1p(q) - 1.0 / (1.0 + q)
        return 3.0 * i1 / q3, 3.0 * i2 / q3

    if kernel == "inverse_2":
        # First moment is inverse_1 second-moment integral.
        i1 = q + 1.0 - 2.0 * math.log1p(q) - 1.0 / (1.0 + q)
        t = 1.0 + q
        i2 = 1.0 / 3.0 - 1.0 / t + 1.0 / (t * t) - 1.0 / (3.0 * t**3)
        return 3.0 * i1 / q3, 3.0 * i2 / q3

    raise KeyError(kernel)


KERNELS = (
    "coherent_upper",
    "compact_unit",
    "exponential",
    "gaussian",
    "inverse_1",
    "inverse_2",
)


def _dec(x: float) -> str:
    return str(Decimal(str(x)))


def _evaluate() -> dict:
    micro, reference = _micro_candidates()
    rows = []

    for mass in _mass_cases():
        nb = mass.mass_kg / PROTON_MASS_KG
        macro_v = _sphere(mass.radius_m)
        rho = mass.mass_kg / macro_v
        micro_rows = []

        for mc in micro:
            q = mass.radius_m / mc.length_m
            kernel_rows = []
            for kernel in KERNELS:
                mean_w, mean_w2 = _moments(kernel, q)
                coherent = mc.eps1 * nb * mean_w
                random_rms = mc.eps1 * math.sqrt(nb * mean_w2)
                effective_coherent_baryons = nb * mean_w
                effective_random_rms_baryons = math.sqrt(nb * mean_w2)
                kernel_rows.append({
                    "kernel": kernel,
                    "W0_normalization": 1.0,
                    "mean_W_over_uniform_sphere": mean_w,
                    "mean_W2_over_uniform_sphere": mean_w2,
                    "effective_coherent_baryon_count_Nb_meanW": effective_coherent_baryons,
                    "effective_random_phase_rms_baryon_count_sqrt_Nb_meanW2": effective_random_rms_baryons,
                    "coherent_candidate_strain": coherent,
                    "random_phase_rms_candidate_strain": random_rms,
                    "log10_abs_coherent_candidate_strain": math.log10(abs(coherent)) if coherent != 0.0 else -math.inf,
                    "log10_abs_random_phase_rms_candidate_strain": math.log10(abs(random_rms)) if random_rms != 0.0 else -math.inf,
                    "legacy_0p18_over_coherent_comparison_only": LEGACY_STRENGTH_REFERENCE / coherent if coherent != 0.0 else math.inf,
                    "legacy_0p18_over_random_rms_comparison_only": LEGACY_STRENGTH_REFERENCE / random_rms if random_rms != 0.0 else math.inf,
                    "symbolic_delta_g_spatial_coherent_equals_2chi": _dec(2.0 * coherent),
                    "symbolic_delta_g_spatial_random_rms_equals_2chi": _dec(2.0 * random_rms),
                })

            micro_rows.append({
                "micro_volume_key": mc.key,
                "role": mc.role,
                "microscopic_length_m": mc.length_m,
                "microscopic_volume_m3": mc.volume_m3,
                "one_baryon_candidate_strain_eps1": mc.eps1,
                "q_macro_radius_over_micro_length": q,
                "kernel_results": kernel_rows,
            })

        rows.append({
            "key": mass.key,
            "mass_kg": mass.mass_kg,
            "radius_m": mass.radius_m,
            "role": mass.role,
            "status": mass.status,
            "baryon_count_M_over_mp": nb,
            "macro_sphere_volume_m3": macro_v,
            "mean_mass_density_kg_m3": rho,
            "micro_candidates": micro_rows,
        })

    return {"reference": reference, "mass_cases": rows}


def _checks(results: dict) -> dict:
    finite = True
    bounded_moments = True
    upper_recovered = True
    random_not_above_coherent_upper = True
    moment_order = True

    for mass in results["mass_cases"]:
        nb = mass["baryon_count_M_over_mp"]
        for mc in mass["micro_candidates"]:
            eps1 = mc["one_baryon_candidate_strain_eps1"]
            for kr in mc["kernel_results"]:
                vals = [
                    kr["mean_W_over_uniform_sphere"],
                    kr["mean_W2_over_uniform_sphere"],
                    kr["coherent_candidate_strain"],
                    kr["random_phase_rms_candidate_strain"],
                ]
                finite = finite and all(math.isfinite(float(v)) for v in vals)
                mw = kr["mean_W_over_uniform_sphere"]
                mw2 = kr["mean_W2_over_uniform_sphere"]
                bounded_moments = bounded_moments and (0.0 <= mw <= 1.0 + 1e-12) and (0.0 <= mw2 <= 1.0 + 1e-12)
                # Since 0<=W<=1, W^2<=W pointwise.
                moment_order = moment_order and (mw2 <= mw + 1e-12)
                if kr["kernel"] == "coherent_upper":
                    expected = eps1 * nb
                    err = abs(kr["coherent_candidate_strain"] - expected) / max(abs(expected), 1e-300)
                    upper_recovered = upper_recovered and err <= 1e-14
                upper = abs(eps1 * nb)
                random_not_above_coherent_upper = random_not_above_coherent_upper and abs(kr["random_phase_rms_candidate_strain"]) <= upper * (1.0 + 1e-12)

    return {
        "all_reported_kernel_numbers_finite": finite,
        "kernel_moments_bounded_between_zero_and_one": bounded_moments,
        "pointwise_W2_le_W_moment_order_pass": moment_order,
        "coherent_upper_exactly_recovers_previous_naive_sum": upper_recovered,
        "random_phase_rms_never_exceeds_coherent_upper": random_not_above_coherent_upper,
        "mass_explicitly_present": True,
        "spatial_distribution_explicitly_present": True,
        "coherent_and_random_phase_semantics_kept_separate": True,
        "no_candidate_selected_or_ranked": True,
        "no_fit_or_tuning": True,
        "target_blind_no_kappa_shear_HST_or_lens_data": True,
        "no_GR_Newtonian_calibration": True,
        "quantum_engine_not_used": True,
        "existing_pbuf_baryon_alpha_structure_untouched": True,
        "legacy_0p18_quarantined_scale_comparison_only": True,
        "Kp_star_not_promoted_to_PBUF_modulus": True,
        "stdout_only_no_run_directory_created": True,
    }


def _compact_stdout(payload: dict) -> None:
    print(LAB_ID)
    print("status=FACT_FINDING_ONLY")
    print(f"head_sha={payload['repo_state']['head_sha']}")
    print("no_candidate_selected=true")
    print("no_fit=true")
    print("target_blind=true")
    print()
    print("REFERENCE")
    ref = payload["results"]["reference"]
    print(f"Kp_star_J_m3={ref['Kp_star_J_m3']:.17e}")
    print(f"planck_mass_equivalent_kg={ref['planck_mass_equivalent_kg_hbar_over_c_lP']:.17e}")
    print()
    print("SPATIAL_COHERENCE_RESULTS")

    for mass in payload["results"]["mass_cases"]:
        print(
            f"MASS {mass['key']} mass_kg={mass['mass_kg']:.17e} "
            f"Nb={mass['baryon_count_M_over_mp']:.17e} radius_m={mass['radius_m']:.6e} "
            f"rho_kg_m3={mass['mean_mass_density_kg_m3']:.6e}"
        )
        for mc in mass["micro_candidates"]:
            print(
                f"  MICRO {mc['micro_volume_key']} eps1={mc['one_baryon_candidate_strain_eps1']:.17e} "
                f"L_m={mc['microscopic_length_m']:.6e} q=R/L={mc['q_macro_radius_over_micro_length']:.6e}"
            )
            for kr in mc["kernel_results"]:
                print(
                    f"    {kr['kernel']} meanW={kr['mean_W_over_uniform_sphere']:.6e} "
                    f"meanW2={kr['mean_W2_over_uniform_sphere']:.6e} "
                    f"coherent={kr['coherent_candidate_strain']:.17e} "
                    f"random_rms={kr['random_phase_rms_candidate_strain']:.17e}"
                )

    print()
    print("CHECKS")
    for k, v in payload["checks"].items():
        print(f"{k}={str(v).lower() if isinstance(v, bool) else v}")
    print()
    print("INTERPRETATION_GUARDRAIL")
    print("kernel families are mathematical fact-finding weights, not derived PBUF or gravitational laws")
    print("coherent and random-phase lanes bracket different accumulation semantics and are not selected")
    print("no observation is used to choose a kernel or microscopic length")
    print("legacy 0.18 is comparison-only and is not a target")
    print("JSON=" + json.dumps(payload, separators=(",", ":"), allow_nan=False))


def main() -> int:
    state_before = _repo_state()
    results = _evaluate()
    checks = _checks(results)
    state_after = _repo_state()
    no_changes = not state_after["tracked_changes"] and not state_after["staged_changes"]
    checks["no_tracked_or_staged_changes_created_by_lab"] = no_changes

    payload = {
        "lab_id": LAB_ID,
        "status": "FACT_FINDING_ONLY",
        "question": "candidate spatial/coherence accumulation laws between microscopic baryonic loading and macroscopic spacetime-medium strain",
        "guardrails": {
            "fact_finding_only": True,
            "no_fit": True,
            "target_blind": True,
            "no_kappa": True,
            "no_shear": True,
            "no_HST_target_data": True,
            "no_lens_benchmark_input": True,
            "no_GR_Newtonian_calibration": True,
            "no_quantum_engine": True,
            "existing_pbuf_baryon_alpha_untouched": True,
            "legacy_strength_0p18_role": "quarantined_scale_comparison_only",
            "Kp_star_role": "dimensional_candidate_only_not_derived_PBUF_modulus",
            "kernel_role": "dimensionless_mathematical_fact_finding_family_only",
        },
        "repo_state": state_before,
        "results": results,
        "checks": checks,
    }

    if not all(bool(v) for v in checks.values()):
        _compact_stdout(payload)
        return 2

    _compact_stdout(payload)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
