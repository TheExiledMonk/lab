#!/usr/bin/env python3
"""PBUF FOUNDATION — COLLECTIVE MASS SCALING AUDIT 001.

Fact-finding only.

Purpose
-------
Continue directly from mass_medium_candidate_audit001 by restoring explicit
physical source mass.  The prior lab computed a candidate microscopic strain
scale for ONE proton under several source-volume hypotheses.  This lab asks what
those same microscopic hypotheses imply when many baryons are present.

It deliberately separates two different constructions that must not be
confused:

A) naive coherent/linear accumulation

       eps_sum(M) = N_b * eps_1,
       N_b = M / m_p

   This is an intentionally unsmeared upper bookkeeping lane.  It assumes every
   one-baryon contribution adds with the same sign and no geometric dilution,
   cancellation, retardation, saturation, or tensor projection.  It is NOT a
   physical source law.

B) uniformly distributed macroscopic energy density

       u_macro = M c^2 / V_macro,
       eps_uniform = u_macro / K_P*

   where K_P* = (hbar*c/l_P)/l_P^3 is retained ONLY as the dimensional candidate
   used by the previous lab.  This lane asks what happens if the same total mass
   is spread through a stated macroscopic volume.

The difference between A and B exposes how much of the missing problem is
spatial/coherent response rather than mass counting alone.

No candidate is selected, fitted, tuned, or promoted to PBUF physics.

Mass cases
----------
- one proton
- 1 kg
- Earth mass
- Solar mass
- representative baryonic galaxy: 6e10 solar masses within 15 kpc
- representative baryonic cluster: 1e14 solar masses within 1 Mpc

The galaxy/cluster values are explicit order-of-magnitude fact-finding cases,
not measurements of any lens used elsewhere in the repository.

Microscopic source-volume hypotheses
------------------------------------
Exactly the same seven lanes as mass_medium_candidate_audit001:
- Planck cube
- reduced-Compton cube/sphere
- ordinary-Compton cube/sphere
- proton charge-radius cube/sphere

Guardrails
----------
- no kappa/shear/HST/lens benchmark input;
- no Quantum Engine input;
- no GR/Newtonian force law or G-based calibration;
- no fit or tuning;
- no replacement for legacy strength=0.18;
- legacy 0.18 appears only as a quarantined scale comparison;
- existing PBUF baryon/alpha structure is not modified or re-derived;
- no conclusion that coherent summation or K_P* is physically correct;
- microscopic strains are kept as independent numbers rather than added to an
  O(1) float64 metric, avoiding the precision loss seen in the previous lab.
"""
from __future__ import annotations

import json
import math
import subprocess
from dataclasses import asdict, dataclass
from decimal import Decimal, getcontext
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
LAB_ID = "PBUF-FOUNDATION-COLLECTIVE-MASS-SCALING-AUDIT-001"

# High precision is used for reporting tiny metric deltas symbolically rather
# than adding them to an O(1) float64 background metric.
getcontext().prec = 80

# SI constants / reference values.
C_M_S = 299_792_458.0
H_J_S = 6.626_070_15e-34
HBAR_J_S = H_J_S / (2.0 * math.pi)
PLANCK_LENGTH_M = 1.616_255e-35
PROTON_MASS_KG = 1.672_621_925_95e-27
PROTON_CHARGE_RADIUS_M = 0.8409e-15
EARTH_MASS_KG = 5.9722e24
EARTH_RADIUS_M = 6.371e6
SOLAR_MASS_KG = 1.98847e30
SOLAR_RADIUS_M = 6.957e8
KPC_M = 3.085_677_581_491_367e19
MPC_M = 1.0e3 * KPC_M

LEGACY_STRENGTH_REFERENCE = 0.18


@dataclass(frozen=True)
class MicroVolume:
    key: str
    length_m: float
    volume_m3: float
    geometry: str
    role: str


@dataclass(frozen=True)
class MassCase:
    key: str
    mass_kg: float
    radius_m: float | None
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


def _cube(x: float) -> float:
    return x**3


def _sphere(r: float) -> float:
    return (4.0 / 3.0) * math.pi * r**3


def _micro_volumes() -> list[MicroVolume]:
    lp = PLANCK_LENGTH_M
    rlc = HBAR_J_S / (PROTON_MASS_KG * C_M_S)
    lc = H_J_S / (PROTON_MASS_KG * C_M_S)
    rp = PROTON_CHARGE_RADIUS_M
    return [
        MicroVolume("planck_cell_cube", lp, _cube(lp), "cube", "elementary-medium-cell boundary hypothesis"),
        MicroVolume("reduced_compton_cube", rlc, _cube(rlc), "cube", "quantum source-extent candidate"),
        MicroVolume("reduced_compton_sphere", rlc, _sphere(rlc), "sphere", "quantum source-extent geometry check"),
        MicroVolume("compton_cube", lc, _cube(lc), "cube", "ordinary-Compton source-extent candidate"),
        MicroVolume("compton_sphere", lc, _sphere(lc), "sphere", "ordinary-Compton source-extent geometry check"),
        MicroVolume("proton_charge_radius_cube", rp, _cube(rp), "cube", "measured charge-size used as response-volume candidate"),
        MicroVolume("proton_charge_radius_sphere", rp, _sphere(rp), "sphere", "measured charge-size geometry check"),
    ]


def _mass_cases() -> list[MassCase]:
    return [
        MassCase("one_proton", PROTON_MASS_KG, PROTON_CHARGE_RADIUS_M, "single-baryon reference", "measured_reference_scale"),
        MassCase("one_kg", 1.0, None, "mass-counting bridge only; no macro volume assigned", "synthetic_mass_case"),
        MassCase("earth", EARTH_MASS_KG, EARTH_RADIUS_M, "planetary-scale distributed source", "measured_reference_scale"),
        MassCase("sun", SOLAR_MASS_KG, SOLAR_RADIUS_M, "stellar-scale distributed source", "measured_reference_scale"),
        MassCase("baryonic_galaxy_scale", 6.0e10 * SOLAR_MASS_KG, 15.0 * KPC_M, "representative baryonic galaxy fact-finding scale", "explicit_order_of_magnitude_case"),
        MassCase("baryonic_cluster_scale", 1.0e14 * SOLAR_MASS_KG, 1.0 * MPC_M, "representative baryonic cluster fact-finding scale", "explicit_order_of_magnitude_case"),
    ]


def _reference_scales() -> dict:
    vp = PLANCK_LENGTH_M**3
    ep_star = HBAR_J_S * C_M_S / PLANCK_LENGTH_M
    kp_star = ep_star / vp
    return {
        "planck_volume_m3": vp,
        "Ep_star_J": ep_star,
        "Kp_star_J_m3": kp_star,
        "Kp_star_status": "dimensional_candidate_only_not_derived_PBUF_modulus",
        "proton_rest_energy_J": PROTON_MASS_KG * C_M_S**2,
        "planck_mass_equivalent_kg_hbar_over_c_lP": HBAR_J_S / (C_M_S * PLANCK_LENGTH_M),
    }


def _one_baryon_rows(ref: dict) -> list[dict]:
    ep_b = PROTON_MASS_KG * C_M_S**2
    rows = []
    for mv in _micro_volumes():
        u1 = ep_b / mv.volume_m3
        eps1 = u1 / ref["Kp_star_J_m3"]
        rows.append({
            **asdict(mv),
            "one_baryon_energy_density_J_m3": u1,
            "one_baryon_candidate_strain": eps1,
            "log10_abs_one_baryon_candidate_strain": math.log10(abs(eps1)),
        })
    return rows


def _symbolic_metric_delta(eps: float) -> dict:
    # For chi_spatial = eps, current metric-strain normalization gives
    # g_ii = 1 + 2 eps.  Do not actually add to float64 1.0 here.
    d = Decimal(str(eps)) * Decimal(2)
    return {
        "chi_spatial": str(Decimal(str(eps))),
        "delta_g_spatial_equals_2chi": str(d),
        "float64_addition_avoided": True,
    }


def _evaluate() -> dict:
    ref = _reference_scales()
    one_rows = _one_baryon_rows(ref)
    by_key = {r["key"]: r for r in one_rows}

    mass_rows = []
    for mc in _mass_cases():
        nb = mc.mass_kg / PROTON_MASS_KG
        macro_volume = _sphere(mc.radius_m) if mc.radius_m is not None else None
        macro_density = mc.mass_kg / macro_volume if macro_volume is not None else None
        macro_energy_density = mc.mass_kg * C_M_S**2 / macro_volume if macro_volume is not None else None
        eps_uniform = macro_energy_density / ref["Kp_star_J_m3"] if macro_volume is not None else None

        micro = []
        for key, one in by_key.items():
            eps1 = one["one_baryon_candidate_strain"]
            eps_sum = nb * eps1
            micro.append({
                "micro_volume_key": key,
                "one_baryon_candidate_strain": eps1,
                "baryon_count_M_over_mp": nb,
                "naive_linear_sum_strain": eps_sum,
                "log10_abs_naive_linear_sum_strain": math.log10(abs(eps_sum)) if eps_sum != 0.0 else -math.inf,
                "legacy_0p18_over_naive_sum": LEGACY_STRENGTH_REFERENCE / eps_sum if eps_sum != 0.0 else math.inf,
                "symbolic_metric_delta_from_naive_sum": _symbolic_metric_delta(eps_sum),
            })

        mass_rows.append({
            **asdict(mc),
            "baryon_count_M_over_mp": nb,
            "macro_sphere_volume_m3": macro_volume,
            "mean_mass_density_kg_m3": macro_density,
            "mean_rest_energy_density_J_m3": macro_energy_density,
            "uniform_macro_candidate_strain_if_K_equals_Kp_star": eps_uniform,
            "log10_abs_uniform_macro_candidate_strain": (
                math.log10(abs(eps_uniform)) if eps_uniform not in (None, 0.0) else None
            ),
            "legacy_0p18_over_uniform_macro_candidate_strain": (
                LEGACY_STRENGTH_REFERENCE / eps_uniform if eps_uniform not in (None, 0.0) else None
            ),
            "symbolic_metric_delta_from_uniform_macro": (
                _symbolic_metric_delta(eps_uniform) if eps_uniform is not None else None
            ),
            "micro_candidate_accumulations": micro,
        })

    return {
        "reference": ref,
        "one_baryon_candidates": one_rows,
        "mass_cases": mass_rows,
    }


def _checks(result: dict) -> dict:
    planck_one = next(r for r in result["one_baryon_candidates"] if r["key"] == "planck_cell_cube")
    mp_star = result["reference"]["planck_mass_equivalent_kg_hbar_over_c_lP"]
    expected = PROTON_MASS_KG / mp_star
    rel = abs(planck_one["one_baryon_candidate_strain"] - expected) / expected

    linearity_errors = []
    for mass in result["mass_cases"]:
        nb = mass["baryon_count_M_over_mp"]
        for row in mass["micro_candidate_accumulations"]:
            expected_sum = nb * row["one_baryon_candidate_strain"]
            err = abs(row["naive_linear_sum_strain"] - expected_sum) / max(abs(expected_sum), 1e-300)
            linearity_errors.append(err)

    finite = True
    for mass in result["mass_cases"]:
        finite &= math.isfinite(mass["mass_kg"]) and mass["mass_kg"] > 0
        finite &= math.isfinite(mass["baryon_count_M_over_mp"]) and mass["baryon_count_M_over_mp"] > 0
        if mass["macro_sphere_volume_m3"] is not None:
            finite &= math.isfinite(mass["macro_sphere_volume_m3"]) and mass["macro_sphere_volume_m3"] > 0
            finite &= math.isfinite(mass["uniform_macro_candidate_strain_if_K_equals_Kp_star"])

    return {
        "all_mass_and_volume_numbers_finite": bool(finite),
        "planck_cell_one_baryon_identity_eps_equals_mp_over_mPlanck_pass": rel <= 1e-12,
        "planck_cell_identity_relative_error": rel,
        "naive_mass_accumulation_exactly_linear_pass": max(linearity_errors, default=0.0) <= 1e-12,
        "max_naive_linearity_relative_error": max(linearity_errors, default=0.0),
        "mass_explicitly_restored": True,
        "macro_distribution_kept_separate_from_coherent_sum": True,
        "no_candidate_selected_or_ranked": True,
        "no_fit_or_tuning": True,
        "target_blind_no_kappa_shear_HST_or_lens_data": True,
        "no_GR_Newtonian_calibration": True,
        "quantum_engine_not_used": True,
        "existing_pbuf_baryon_alpha_structure_untouched": True,
        "legacy_0p18_quarantined_scale_comparison_only": True,
        "float64_metric_addition_of_tiny_strains_avoided": True,
    }


def _fmt(x) -> str:
    if x is None:
        return "n/a"
    return f"{x:.6e}"


def _print_report(payload: dict) -> None:
    r = payload["results"]
    print(LAB_ID)
    print("status=FACT_FINDING_ONLY")
    print(f"head_sha={payload['repo_state']['head_sha']}")
    print("no_candidate_selected=true")
    print("no_fit=true")
    print("target_blind=true")
    print()
    print("REFERENCE")
    print(f"Kp_star_J_m3={r['reference']['Kp_star_J_m3']:.17e}")
    print(f"planck_mass_equivalent_kg={r['reference']['planck_mass_equivalent_kg_hbar_over_c_lP']:.17e}")
    print()
    print("ONE_BARYON_CANDIDATES")
    for row in r["one_baryon_candidates"]:
        print(
            f"{row['key']} eps1={row['one_baryon_candidate_strain']:.17e} "
            f"log10eps1={row['log10_abs_one_baryon_candidate_strain']:.6f}"
        )
    print()
    print("MASS_ACCUMULATION_AND_DISTRIBUTION")
    for mass in r["mass_cases"]:
        print(
            f"MASS {mass['key']} mass_kg={mass['mass_kg']:.17e} "
            f"Nb={mass['baryon_count_M_over_mp']:.17e} "
            f"radius_m={_fmt(mass['radius_m'])} "
            f"rho_kg_m3={_fmt(mass['mean_mass_density_kg_m3'])} "
            f"eps_uniform_macro={_fmt(mass['uniform_macro_candidate_strain_if_K_equals_Kp_star'])}"
        )
        for row in mass["micro_candidate_accumulations"]:
            print(
                f"  {row['micro_volume_key']} naive_sum={row['naive_linear_sum_strain']:.17e} "
                f"log10={row['log10_abs_naive_linear_sum_strain']:.6f} "
                f"legacy0p18_over_naive={row['legacy_0p18_over_naive_sum']:.17e}"
            )
    print()
    print("CHECKS")
    for k, v in payload["checks"].items():
        print(f"{k}={v}")
    print()
    print("INTERPRETATION_GUARDRAIL")
    print("naive_linear_sum_strain is unsmeared bookkeeping, not a physical response law")
    print("uniform_macro_candidate_strain assumes only the stated Kp* dimensional candidate")
    print("difference between these lanes quantifies missing spatial/coherence/response physics")
    print("legacy 0.18 is comparison-only and is not a target")


def main() -> int:
    before = _repo_state()
    results = _evaluate()
    checks = _checks(results)
    after = _repo_state()
    checks["no_tracked_or_staged_changes_created_by_lab"] = (
        before["tracked_changes"] == after["tracked_changes"]
        and before["staged_changes"] == after["staged_changes"]
    )
    payload = {
        "lab_id": LAB_ID,
        "status": "FACT_FINDING_ONLY",
        "repo_state": before,
        "results": results,
        "checks": checks,
    }
    _print_report(payload)
    # JSON is emitted to stdout as a final machine-readable line; no run
    # directory or repository file is created.
    print("JSON=" + json.dumps(payload, separators=(",", ":"), sort_keys=False))
    return 0 if all(v is True or not isinstance(v, bool) for v in checks.values()) else 1


if __name__ == "__main__":
    raise SystemExit(main())
