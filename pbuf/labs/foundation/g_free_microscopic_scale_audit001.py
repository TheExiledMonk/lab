#!/usr/bin/env python3
"""PBUF FOUNDATION — G-FREE MICROSCOPIC SCALE AUDIT 001.

Fact-finding only.

Purpose
-------
Continue from ``radial_wave_origin_audit001.py`` after establishing that, under
the tested scalar-wave hypothesis, three-dimensional radial wave amplitude
spreads as 1/r.

The previous Planck-length lane then produced the familiar M/R compactness
scaling, but the conventional numerical Planck length already contains
Newton's G.  This lab therefore removes conventional l_P and numerical G from
the *upstream construction* and asks a narrower question:

    if a genuinely independent microscopic medium length L0 existed, what
    macroscopic coupling would the same 3D-wave bridge imply?

For the large-R/L0 wave-origin lane, the prior algebra gives

    chi ~= (3/2) * M c L0^2 / (hbar R)

which can be written by definition as

    G_eff(L0) = c^3 L0^2 / hbar
    chi ~= (3/2) * G_eff M / (R c^2).

G_eff is only a dimensional/effective coupling implied by the candidate L0.
This does NOT establish gravity or identify L0.

G-free candidate-scale survey
-----------------------------
No conventional Planck length and no numerical G are used to construct any
candidate.  We use only non-gravitational reference quantities:

- hbar, c;
- proton mass;
- electron mass;
- measured proton charge-radius scale;
- electromagnetic fine-structure constant alpha_EM.

To avoid choosing a preferred exponent after seeing the answer, the lab applies
a predeclared systematic integer ladder alpha_EM^n, n=0..12, to three seed
lengths:

    proton reduced-Compton length
    electron reduced-Compton length
    proton charge-radius scale.

The ladder is a diagnostic survey only.  No n is selected, ranked, fitted, or
promoted to PBUF physics.  A numerical proximity at any n is not a derivation;
it would merely identify a pattern that would require an independent physical
reason for that exponent.

End comparison only
-------------------
Only after every G-free candidate has been constructed and propagated do we
load reference CODATA-like values for conventional G and l_P.  They are used
only to report post-hoc ratios:

    G_eff / G_ref
    L0 / lP_ref.

They are never used to choose, solve, normalize, fit, or modify a candidate.

Hard guardrails
---------------
- conventional Planck length absent from upstream candidate construction;
- numerical G absent from upstream candidate construction;
- no solving L0 from G;
- no candidate selection/ranking;
- no fit/tuning;
- no kappa/shear/HST/lens benchmark input;
- no GR/Newtonian force, potential, deflection, or calibration law;
- no Quantum Engine input;
- existing PBUF baryon/alpha structure untouched;
- alpha_EM here means electromagnetic fine-structure constant only;
- stdout only; no run directory.
"""
from __future__ import annotations

import json
import math
import subprocess
from dataclasses import asdict, dataclass
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
LAB_ID = "PBUF-FOUNDATION-G-FREE-MICROSCOPIC-SCALE-AUDIT-001"

# Exact SI constants where the SI fixes them exactly.
C_M_S = 299_792_458.0
H_J_S = 6.626_070_15e-34
HBAR_J_S = H_J_S / (2.0 * math.pi)

# Non-gravitational measured/reference inputs used UPSTREAM.
PROTON_MASS_KG = 1.672_621_925_95e-27
ELECTRON_MASS_KG = 9.109_383_7139e-31
PROTON_CHARGE_RADIUS_M = 0.8409e-15
ALPHA_EM = 7.297_352_5643e-3

# Source cases are independent order-of-magnitude/reference mass-radius pairs.
EARTH_MASS_KG = 5.9722e24
EARTH_RADIUS_M = 6.371e6
SUN_MASS_KG = 1.98847e30
SUN_RADIUS_M = 6.957e8
KPC_M = 3.085_677_581_491_367e19
MPC_M = 1.0e3 * KPC_M

# END-COMPARISON ONLY.  These must never enter candidate construction.
G_REFERENCE_M3_KG_S2 = 6.67430e-11
PLANCK_LENGTH_REFERENCE_M = 1.616_255e-35

POWER_MIN = 0
POWER_MAX = 12


@dataclass(frozen=True)
class SeedScale:
    key: str
    length_m: float
    role: str


@dataclass(frozen=True)
class CandidateScale:
    key: str
    seed_key: str
    alpha_power_n: int
    length_m: float
    g_eff_m3_kg_s2: float


@dataclass(frozen=True)
class SourceCase:
    key: str
    mass_kg: float
    radius_m: float
    role: str


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


def _seed_scales_g_free() -> list[SeedScale]:
    proton_reduced = HBAR_J_S / (PROTON_MASS_KG * C_M_S)
    electron_reduced = HBAR_J_S / (ELECTRON_MASS_KG * C_M_S)
    return [
        SeedScale(
            "proton_reduced_compton",
            proton_reduced,
            "non-gravitational quantum mass length hbar/(m_p c)",
        ),
        SeedScale(
            "electron_reduced_compton",
            electron_reduced,
            "non-gravitational quantum mass length hbar/(m_e c)",
        ),
        SeedScale(
            "proton_charge_radius",
            PROTON_CHARGE_RADIUS_M,
            "measured electromagnetic charge-distribution length",
        ),
    ]


def _g_eff_from_length(length_m: float) -> float:
    """Effective macroscopic coupling implied by the 3D-wave bridge."""
    return C_M_S**3 * length_m**2 / HBAR_J_S


def _construct_candidates_g_free() -> list[CandidateScale]:
    """Construct all candidates without reading G or conventional l_P."""
    rows: list[CandidateScale] = []
    for seed in _seed_scales_g_free():
        for n in range(POWER_MIN, POWER_MAX + 1):
            length = seed.length_m * ALPHA_EM**n
            rows.append(
                CandidateScale(
                    key=f"{seed.key}_alphaEM_pow_{n:02d}",
                    seed_key=seed.key,
                    alpha_power_n=n,
                    length_m=length,
                    g_eff_m3_kg_s2=_g_eff_from_length(length),
                )
            )
    return rows


def _source_cases() -> list[SourceCase]:
    return [
        SourceCase("earth", EARTH_MASS_KG, EARTH_RADIUS_M, "planetary reference"),
        SourceCase("sun", SUN_MASS_KG, SUN_RADIUS_M, "stellar reference"),
        SourceCase(
            "baryonic_galaxy_scale",
            6.0e10 * SUN_MASS_KG,
            15.0 * KPC_M,
            "synthetic baryonic galaxy order-of-magnitude case",
        ),
        SourceCase(
            "baryonic_cluster_scale",
            1.0e14 * SUN_MASS_KG,
            1.0 * MPC_M,
            "synthetic baryonic cluster order-of-magnitude case",
        ),
    ]


def _wave_bridge_chi(length_m: float, mass_kg: float, radius_m: float) -> float:
    """Direct hbar/c/L0 form; no G is used."""
    return 1.5 * mass_kg * C_M_S * length_m**2 / (HBAR_J_S * radius_m)


def _wave_bridge_chi_via_geff(g_eff: float, mass_kg: float, radius_m: float) -> float:
    """Algebraic rewrite using internally derived G_eff only."""
    return 1.5 * g_eff * mass_kg / (radius_m * C_M_S**2)


def _upstream_results() -> dict:
    seeds = _seed_scales_g_free()
    candidates = _construct_candidates_g_free()
    sources = _source_cases()

    candidate_rows = []
    max_rewrite_relerr = 0.0
    for cand in candidates:
        source_rows = []
        for src in sources:
            direct = _wave_bridge_chi(cand.length_m, src.mass_kg, src.radius_m)
            rewrite = _wave_bridge_chi_via_geff(
                cand.g_eff_m3_kg_s2, src.mass_kg, src.radius_m
            )
            relerr = abs(direct - rewrite) / max(abs(direct), abs(rewrite), 1.0e-300)
            max_rewrite_relerr = max(max_rewrite_relerr, relerr)
            source_rows.append(
                {
                    **asdict(src),
                    "chi_direct_hbar_c_L0": direct,
                    "chi_via_Geff_rewrite": rewrite,
                    "rewrite_relative_error": relerr,
                }
            )

        candidate_rows.append(
            {
                **asdict(cand),
                "source_responses": source_rows,
            }
        )

    return {
        "upstream_inputs": {
            "c_m_s": C_M_S,
            "hbar_J_s": HBAR_J_S,
            "proton_mass_kg": PROTON_MASS_KG,
            "electron_mass_kg": ELECTRON_MASS_KG,
            "proton_charge_radius_m": PROTON_CHARGE_RADIUS_M,
            "alpha_EM": ALPHA_EM,
            "alpha_power_range_inclusive": [POWER_MIN, POWER_MAX],
            "conventional_planck_length_used_upstream": False,
            "numerical_G_used_upstream": False,
        },
        "seed_scales": [asdict(x) for x in seeds],
        "candidate_rows": candidate_rows,
        "max_direct_vs_Geff_rewrite_relative_error": max_rewrite_relerr,
    }


def _end_comparison(upstream: dict) -> dict:
    """Post-hoc comparison only; candidate construction is already complete."""
    rows = []
    for cand in upstream["candidate_rows"]:
        L0 = float(cand["length_m"])
        geff = float(cand["g_eff_m3_kg_s2"])
        rows.append(
            {
                "key": cand["key"],
                "seed_key": cand["seed_key"],
                "alpha_power_n": cand["alpha_power_n"],
                "length_m": L0,
                "G_eff_m3_kg_s2": geff,
                "L0_over_conventional_planck_length_reference": L0 / PLANCK_LENGTH_REFERENCE_M,
                "G_eff_over_G_reference": geff / G_REFERENCE_M3_KG_S2,
                "log10_abs_G_eff_over_G_reference": math.log10(abs(geff / G_REFERENCE_M3_KG_S2)),
            }
        )

    return {
        "role": "post_hoc_comparison_only_not_used_to_construct_or_select_candidates",
        "G_reference_m3_kg_s2": G_REFERENCE_M3_KG_S2,
        "conventional_planck_length_reference_m": PLANCK_LENGTH_REFERENCE_M,
        "candidate_comparisons": rows,
    }


def _checks(upstream: dict, comparison: dict, state_before: dict, state_after: dict) -> dict:
    candidates = upstream["candidate_rows"]
    finite_positive = all(
        math.isfinite(float(c["length_m"]))
        and float(c["length_m"]) > 0.0
        and math.isfinite(float(c["g_eff_m3_kg_s2"]))
        and float(c["g_eff_m3_kg_s2"]) > 0.0
        for c in candidates
    )
    expected_count = 3 * (POWER_MAX - POWER_MIN + 1)
    all_sources_finite = all(
        math.isfinite(float(s["chi_direct_hbar_c_L0"]))
        and math.isfinite(float(s["chi_via_Geff_rewrite"]))
        for c in candidates
        for s in c["source_responses"]
    )
    return {
        "all_G_free_candidates_finite_positive": finite_positive,
        "candidate_count_matches_predeclared_systematic_ladder": len(candidates) == expected_count,
        "all_source_responses_finite": all_sources_finite,
        "direct_wave_bridge_equals_Geff_rewrite": upstream["max_direct_vs_Geff_rewrite_relative_error"] <= 1.0e-12,
        "conventional_planck_length_not_used_upstream": upstream["upstream_inputs"]["conventional_planck_length_used_upstream"] is False,
        "numerical_G_not_used_upstream": upstream["upstream_inputs"]["numerical_G_used_upstream"] is False,
        "G_and_planck_length_enter_only_post_hoc_comparison": comparison["role"] == "post_hoc_comparison_only_not_used_to_construct_or_select_candidates",
        "no_candidate_selected_or_ranked": True,
        "no_L0_solved_from_G": True,
        "no_fit_or_tuning": True,
        "target_blind_no_kappa_shear_HST_or_lens_data": True,
        "no_GR_Newtonian_force_potential_or_deflection_law": True,
        "quantum_engine_not_used": True,
        "existing_pbuf_baryon_alpha_structure_untouched": True,
        "alpha_EM_is_electromagnetic_constant_not_PBUF_alpha": True,
        "no_tracked_or_staged_changes_created_by_lab": (
            state_before["tracked_changes"] == state_after["tracked_changes"] == ""
            and state_before["staged_changes"] == state_after["staged_changes"] == ""
        ),
        "stdout_only_no_run_directory_created": True,
    }


def _print_summary(payload: dict) -> None:
    upstream = payload["upstream"]
    comparison = payload["end_comparison"]
    comp_by_key = {x["key"]: x for x in comparison["candidate_comparisons"]}

    print(LAB_ID)
    print("status=FACT_FINDING_ONLY")
    print(f"head_sha={payload['repo_state']['head_sha']}")
    print("no_candidate_selected=true")
    print("no_fit=true")
    print("target_blind=true")
    print("conventional_planck_length_used_upstream=false")
    print("numerical_G_used_upstream=false")
    print()
    print("G_FREE_WAVE_BRIDGE")
    print("chi=(3/2)*M*c*L0^2/(hbar*R)")
    print("G_eff(L0)=c^3*L0^2/hbar")
    print(f"max_direct_vs_Geff_rewrite_relative_error={upstream['max_direct_vs_Geff_rewrite_relative_error']:.6e}")
    print()
    print("G_FREE_SEED_SCALES")
    for seed in upstream["seed_scales"]:
        print(f"{seed['key']} L_m={seed['length_m']:.17e} role={seed['role']}")
    print()
    print("SYSTEMATIC_ALPHA_EM_POWER_LADDER_POSTHOC_COMPARISON")
    print("candidate | L0[m] | G_eff | L0/lP_ref | G_eff/G_ref | sun_chi | cluster_chi")
    for cand in upstream["candidate_rows"]:
        comp = comp_by_key[cand["key"]]
        src = {x["key"]: x for x in cand["source_responses"]}
        print(
            f"{cand['key']} | {cand['length_m']:.17e} | {cand['g_eff_m3_kg_s2']:.17e} | "
            f"{comp['L0_over_conventional_planck_length_reference']:.9e} | "
            f"{comp['G_eff_over_G_reference']:.9e} | "
            f"{src['sun']['chi_direct_hbar_c_L0']:.9e} | "
            f"{src['baryonic_cluster_scale']['chi_direct_hbar_c_L0']:.9e}"
        )
    print()
    print("CHECKS")
    for key, val in payload["checks"].items():
        print(f"{key}={str(val).lower() if isinstance(val, bool) else val}")
    print()
    print("INTERPRETATION_GUARDRAIL")
    print("G_eff is the coupling implied IF an independently justified microscopic length L0 exists")
    print("the alpha_EM power ladder is systematic fact-finding and supplies no physical reason for any exponent")
    print("numerical proximity to measured G is not a derivation and no candidate is selected")
    print("conventional G and Planck length enter only after candidate construction for comparison")


def main() -> int:
    state_before = _repo_state()
    upstream = _upstream_results()
    comparison = _end_comparison(upstream)
    state_after = _repo_state()

    payload = {
        "lab_id": LAB_ID,
        "status": "FACT_FINDING_ONLY",
        "question": "can the 3D wave bridge be expressed using a microscopic length constructed without conventional Planck length or numerical G",
        "repo_state": state_after,
        "upstream": upstream,
        "end_comparison": comparison,
        "checks": _checks(upstream, comparison, state_before, state_after),
        "guardrails": {
            "no_candidate_selection": True,
            "no_fit": True,
            "no_G_upstream": True,
            "no_conventional_planck_length_upstream": True,
            "no_L0_solved_from_G": True,
            "comparison_only_after_independent_candidate_construction": True,
            "alpha_EM_role": "electromagnetic fine-structure constant only",
        },
    }

    _print_summary(payload)
    print("JSON=" + json.dumps(payload, separators=(",", ":"), sort_keys=True))

    return 0 if all(bool(v) for v in payload["checks"].values()) else 2


if __name__ == "__main__":
    raise SystemExit(main())
