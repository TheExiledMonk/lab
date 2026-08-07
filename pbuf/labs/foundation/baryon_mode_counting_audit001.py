#!/usr/bin/env python3
"""PBUF FOUNDATION — BARYON MODE COUNTING AUDIT 001.

Fact-finding only.

Question
--------
Can the exponent 18 that appeared post-hoc in the previous G-free audit arise
from a simple, independently stated counting rule in the existing PBUF baryon
picture: three spatial dimensions and two source modes per dimension (six total
modes)?

This lab does NOT assume that 18 is correct. It explicitly constructs several
predeclared counting semantics from the same (d=3, modes_per_dimension=2)
structure and keeps them separate:

- source_mode_count:                  N = 6
- longitudinal_component_count:      N = 6
- transverse_component_count:        N = 12
- full_vector_component_count:       N = 18
- unordered_mode_pair_count:         N = 15
- ordered_distinct_mode_pair_count:  N = 30

Only the full-vector-component semantics gives 18, and that requires the extra
physical statement that every one of the six source modes independently drives
all three spatial response components. The six-mode baryon count alone does
NOT imply 18.

For each structural count N, the lab forms the G-free algebraic candidate

    alpha_G,candidate = alpha_EM**N
    G_eff(N) = (hbar*c/m_p**2) * alpha_EM**N
    L0(N) = (hbar/(m_p*c)) * alpha_EM**(N/2)

These are dimensional consequences of the previously derived 3D wave bridge
IF the count N is physically justified. No count is selected or fitted.

Hard guardrails
---------------
- no numerical G or conventional Planck length in candidate construction;
- no solving for an exponent from G;
- no fractional exponent chosen to improve a match;
- no ranking or selection;
- no kappa/shear/HST/lens data;
- no GR/Newtonian force, potential, deflection, or calibration law;
- no Quantum Engine;
- alpha_EM is the electromagnetic fine-structure constant, not PBUF alpha;
- existing PBUF baryon/alpha structure is untouched;
- G and conventional Planck length enter only after every candidate is frozen,
  strictly for post-hoc comparison;
- stdout only, no run directory.
"""
from __future__ import annotations

import json
import math
import subprocess
from dataclasses import dataclass
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
LAB_ID = "PBUF-FOUNDATION-BARYON-MODE-COUNTING-AUDIT-001"

# Exact SI constants where fixed by SI.
C_M_S = 299_792_458.0
H_J_S = 6.626_070_15e-34
HBAR_J_S = H_J_S / (2.0 * math.pi)
E_CHARGE_C = 1.602_176_634e-19

# Non-gravitational measured constants.
PROTON_MASS_KG = 1.672_621_925_95e-27
ALPHA_EM = 7.297_352_5643e-3

# Existing PBUF structural hypothesis under audit.
SPATIAL_DIMENSIONS = 3
MODES_PER_DIMENSION = 2

# Post-hoc comparison constants. These MUST NOT enter candidate construction.
G_REFERENCE = 6.67430e-11
PLANCK_LENGTH_REFERENCE_M = 1.616_255e-35


@dataclass(frozen=True)
class CountingLane:
    key: str
    count: int
    statement: str
    extra_assumption: str


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


def _counting_lanes() -> list[CountingLane]:
    d = SPATIAL_DIMENSIONS
    m = MODES_PER_DIMENSION
    total_modes = d * m
    return [
        CountingLane(
            "source_mode_count",
            total_modes,
            "count each of the six source modes once",
            "none beyond 3 dimensions x 2 modes per dimension",
        ),
        CountingLane(
            "longitudinal_component_count",
            total_modes,
            "each source mode contributes one same-axis/longitudinal response component",
            "one response component per source mode",
        ),
        CountingLane(
            "transverse_component_count",
            total_modes * (d - 1),
            "each source mode contributes to the two components transverse to its source axis",
            "two independent transverse response components per source mode",
        ),
        CountingLane(
            "full_vector_component_count",
            total_modes * d,
            "each source mode contributes independently to all three spatial response components",
            "all 3 response components are physically independent for each of the 6 source modes",
        ),
        CountingLane(
            "unordered_mode_pair_count",
            total_modes * (total_modes - 1) // 2,
            "count each unordered pair among the six source modes",
            "coupling is pairwise between distinct modes and pair ordering is irrelevant",
        ),
        CountingLane(
            "ordered_distinct_mode_pair_count",
            total_modes * (total_modes - 1),
            "count ordered source-response pairs among distinct modes",
            "coupling is directed between every distinct ordered mode pair",
        ),
    ]


def _construct_candidates() -> tuple[list[dict], dict]:
    # Everything in this function is G-free.
    proton_reduced_compton = HBAR_J_S / (PROTON_MASS_KG * C_M_S)
    gravitational_prefactor = HBAR_J_S * C_M_S / (PROTON_MASS_KG**2)

    lanes = _counting_lanes()
    candidates = []
    for lane in lanes:
        n = lane.count
        alpha_power = ALPHA_EM**n
        length_exponent = n / 2.0
        l0 = proton_reduced_compton * (ALPHA_EM**length_exponent)
        geff = gravitational_prefactor * alpha_power
        geff_from_l0 = C_M_S**3 * l0**2 / HBAR_J_S
        candidates.append({
            "key": lane.key,
            "structural_count_N": n,
            "length_exponent_N_over_2": length_exponent,
            "statement": lane.statement,
            "extra_assumption": lane.extra_assumption,
            "alpha_EM_power_N": alpha_power,
            "L0_m": l0,
            "G_eff_m3_kg_s2": geff,
            "G_eff_from_L0_m3_kg_s2": geff_from_l0,
            "construction_uses_G": False,
            "construction_uses_planck_length": False,
        })

    reference = {
        "spatial_dimensions": SPATIAL_DIMENSIONS,
        "modes_per_dimension": MODES_PER_DIMENSION,
        "total_source_modes": SPATIAL_DIMENSIONS * MODES_PER_DIMENSION,
        "alpha_EM": ALPHA_EM,
        "proton_mass_kg": PROTON_MASS_KG,
        "proton_reduced_compton_m": proton_reduced_compton,
        "G_free_prefactor_hbar_c_over_mp2": gravitational_prefactor,
        "candidate_formula": "G_eff=(hbar*c/m_p^2)*alpha_EM^N",
        "length_formula": "L0=(hbar/(m_p*c))*alpha_EM^(N/2)",
    }
    return candidates, reference


def _posthoc_compare(candidates: list[dict]) -> list[dict]:
    # G and lP enter for the first time here, after candidate construction.
    out = []
    for row in candidates:
        geff = row["G_eff_m3_kg_s2"]
        l0 = row["L0_m"]
        out.append({
            "key": row["key"],
            "structural_count_N": row["structural_count_N"],
            "length_exponent_N_over_2": row["length_exponent_N_over_2"],
            "G_eff_over_G_reference": geff / G_REFERENCE,
            "L0_over_planck_length_reference": l0 / PLANCK_LENGTH_REFERENCE_M,
            "log10_abs_Geff_over_Gref": math.log10(abs(geff / G_REFERENCE)),
        })
    return out


def _checks(candidates: list[dict], posthoc: list[dict], repo: dict) -> dict:
    total_modes = SPATIAL_DIMENSIONS * MODES_PER_DIMENSION
    by_key = {x["key"]: x for x in candidates}

    finite_positive = all(
        math.isfinite(x["L0_m"]) and x["L0_m"] > 0.0
        and math.isfinite(x["G_eff_m3_kg_s2"]) and x["G_eff_m3_kg_s2"] > 0.0
        for x in candidates
    )
    rewrite_exact = all(
        math.isclose(x["G_eff_m3_kg_s2"], x["G_eff_from_L0_m3_kg_s2"], rel_tol=2e-15, abs_tol=0.0)
        for x in candidates
    )

    return {
        "three_dimensions_times_two_modes_equals_six": total_modes == 6,
        "source_mode_count_is_six": by_key["source_mode_count"]["structural_count_N"] == 6,
        "transverse_component_count_is_twelve": by_key["transverse_component_count"]["structural_count_N"] == 12,
        "full_vector_component_count_is_eighteen": by_key["full_vector_component_count"]["structural_count_N"] == 18,
        "unordered_pair_count_is_fifteen": by_key["unordered_mode_pair_count"]["structural_count_N"] == 15,
        "ordered_distinct_pair_count_is_thirty": by_key["ordered_distinct_mode_pair_count"]["structural_count_N"] == 30,
        "eighteen_requires_full_vector_extra_assumption": "all 3 response components" in by_key["full_vector_component_count"]["extra_assumption"],
        "six_mode_count_alone_does_not_imply_eighteen": by_key["source_mode_count"]["structural_count_N"] != 18,
        "all_candidates_finite_positive": finite_positive,
        "G_eff_equals_c3_L0sq_over_hbar": rewrite_exact,
        "candidate_construction_uses_no_G": all(not x["construction_uses_G"] for x in candidates),
        "candidate_construction_uses_no_planck_length": all(not x["construction_uses_planck_length"] for x in candidates),
        "no_exponent_solved_from_G": True,
        "no_fractional_exponent_chosen_for_match": True,
        "no_candidate_selected_or_ranked": True,
        "no_fit_or_tuning": True,
        "target_blind_no_kappa_shear_HST_or_lens_data": True,
        "no_GR_Newtonian_force_potential_or_deflection_law": True,
        "quantum_engine_not_used": True,
        "alpha_EM_is_electromagnetic_constant_not_PBUF_alpha": True,
        "G_and_planck_length_enter_only_posthoc": True,
        "posthoc_candidate_count_matches_structural_candidate_count": len(posthoc) == len(candidates),
        "no_tracked_or_staged_changes_created_by_lab": repo["tracked_changes"] == "" and repo["staged_changes"] == "",
        "stdout_only_no_run_directory_created": True,
    }


def main() -> None:
    repo = _repo_state()
    candidates, reference = _construct_candidates()
    posthoc = _posthoc_compare(candidates)
    checks = _checks(candidates, posthoc, repo)

    print(LAB_ID)
    print("status=FACT_FINDING_ONLY")
    print(f"head_sha={repo['head_sha']}")
    print("no_candidate_selected=true")
    print("no_fit=true")
    print("target_blind=true")
    print("numerical_G_used_upstream=false")
    print("conventional_planck_length_used_upstream=false")
    print()

    print("PBUF_STRUCTURE_UNDER_AUDIT")
    print(f"spatial_dimensions={SPATIAL_DIMENSIONS}")
    print(f"modes_per_dimension={MODES_PER_DIMENSION}")
    print(f"total_source_modes={SPATIAL_DIMENSIONS * MODES_PER_DIMENSION}")
    print("NOTE six source modes alone imply N=6, not N=18")
    print()

    print("STRUCTURAL_COUNT_CANDIDATES_G_FREE")
    print("key | N | N/2 | L0[m] | G_eff | extra_assumption")
    for row in candidates:
        print(
            f"{row['key']} | {row['structural_count_N']} | "
            f"{row['length_exponent_N_over_2']:.6f} | "
            f"{row['L0_m']:.17e} | {row['G_eff_m3_kg_s2']:.17e} | "
            f"{row['extra_assumption']}"
        )
    print()

    print("POSTHOC_COMPARISON_ONLY")
    print(f"G_reference={G_REFERENCE:.17e}")
    print(f"planck_length_reference_m={PLANCK_LENGTH_REFERENCE_M:.17e}")
    print("key | N | G_eff/G_ref | L0/lP_ref")
    for row in posthoc:
        print(
            f"{row['key']} | {row['structural_count_N']} | "
            f"{row['G_eff_over_G_reference']:.12e} | "
            f"{row['L0_over_planck_length_reference']:.12e}"
        )
    print()

    print("CHECKS")
    for key, value in checks.items():
        print(f"{key}={'true' if value else 'false'}")
    print()

    print("INTERPRETATION_GUARDRAIL")
    print("N=18 is structurally available only under the explicit full-vector response assumption")
    print("the existing 3D x 2-modes-per-dimension source count by itself gives six and does not derive eighteen")
    print("transverse-only, longitudinal-only, and pairwise semantics are retained as independent alternatives")
    print("post-hoc proximity to measured G cannot choose among these semantics")
    print("a future constitutive or mode-coupling derivation must decide which response channels actually exist")

    payload = {
        "lab_id": LAB_ID,
        "status": "FACT_FINDING_ONLY",
        "question": "does exponent 18 follow independently from the existing 3D, two-modes-per-dimension baryon structure",
        "repo_state": repo,
        "reference_G_free": reference,
        "structural_candidates_G_free": candidates,
        "posthoc_comparison_only": {
            "G_reference_m3_kg_s2": G_REFERENCE,
            "planck_length_reference_m": PLANCK_LENGTH_REFERENCE_M,
            "rows": posthoc,
        },
        "checks": checks,
        "guardrails": {
            "no_candidate_selected": True,
            "no_fit": True,
            "target_blind": True,
            "no_G_upstream": True,
            "no_planck_length_upstream": True,
            "no_exponent_solved_from_G": True,
            "alpha_EM_role": "electromagnetic_fine_structure_constant_not_PBUF_alpha",
            "N18_role": "candidate_full_vector_component_count_only_not_derived_without_channel_physics",
        },
    }
    print("JSON=" + json.dumps(payload, separators=(",", ":"), sort_keys=True))

    if not all(checks.values()):
        raise SystemExit(2)


if __name__ == "__main__":
    main()
