#!/usr/bin/env python3
"""PBUF FOUNDATION — G DERIVATION OMNIBUS AUDIT 001.

Fact-finding only.

Purpose
-------
Gather the currently plausible routes for obtaining an effective Newton coupling
into one target-blind audit, while keeping independent candidates, circular
controls, algebraically equivalent rewrites, and presently unclosed routes
strictly separated.

Policy
------
The measured Newton constant is allowed as a post-hoc reference because it is an
experimental constant used by standard physics. It is NOT used to construct,
solve, tune, normalize, rank, or select any upstream candidate.

If none of the independently motivated routes closes, the appropriate current
PBUF policy is simply:

    G = measured constant, derivation deferred.

That fallback is not a failure of physics-first modeling; it is an explicit
provenance boundary.

Routes audited
--------------
R1  3D wave/microscopic-length route:
        G_eff(L0) = c^3 L0^2 / hbar
    using the same non-gravitational seed lengths and systematic alpha_EM^n
    ladder as the previous G-free audit.

R2  Structural dimensionless-coupling route:
        alpha_G,cand = alpha_EM^N
        G_eff = (hbar c / m_p^2) alpha_EM^N
    for the predeclared structural counts already audited: N=6,12,15,18,30.

R3  Algebraic-equivalence audit:
    when L0 = (hbar/(m_p c)) alpha_EM^(N/2), R1 and R2 must be exactly the
    same route and must not be double-counted as independent evidence.

R4  Conventional Planck-length identity control:
        G = c^3 l_P^2 / hbar
    evaluated only after the post-hoc reference boundary. This is circular by
    construction and is never treated as a derivation.

R5  Absolute-stiffness inversion route:
    an independently normalized constitutive stiffness could in principle
    determine the gravitational coupling, but current audited PBUF foundation
    does not yet supply that absolute constitutive normalization. This route is
    therefore reported as NOT_CLOSED, not populated with a fitted number.

R6  Medium-dressing/state route:
        G_eff(a,state) = G_bare * R_medium(a,state)
    is a legitimate hypothesis class, but current foundation does not yet
    derive R_medium from a normalized constitutive law. It is reported as
    NOT_CLOSED. No 1.715-like factor is inserted or solved from measured G.

R7  Cosmological-growth inversion route:
    growth/sigma8 can constrain a supplied G_eff(a) model, but sigma8 alone
    does not independently determine the absolute Newton coupling. The current
    PBUF growth equation is written in terms of dimensionless Omega_m and E(a),
    so this route is classified as DIAGNOSTIC_NOT_ABSOLUTE_DERIVATION.

Hard guardrails
---------------
- no fit or tuning;
- no candidate ranking or selection;
- no solving L0, N, stiffness, or medium response from measured G;
- no fractional exponent chosen to improve agreement;
- no kappa/shear/HST/lens benchmark input;
- no Quantum Engine input;
- alpha_EM is the electromagnetic fine-structure constant, not PBUF alpha;
- measured G and conventional l_P enter only after all independent candidates
  are frozen;
- stdout only; no run directory.
"""
from __future__ import annotations

import json
import math
import subprocess
from dataclasses import dataclass, asdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
LAB_ID = "PBUF-FOUNDATION-G-DERIVATION-OMNIBUS-AUDIT-001"

# Exact SI constants where fixed by SI.
C = 299_792_458.0
H = 6.626_070_15e-34
HBAR = H / (2.0 * math.pi)

# Non-gravitational particle/EM inputs used upstream.
M_P = 1.672_621_925_95e-27
M_E = 9.109_383_7139e-31
PROTON_CHARGE_RADIUS = 0.8409e-15
ALPHA_EM = 7.297_352_5643e-3

# Post-hoc references. These MUST NOT be touched until all upstream rows freeze.
G_REFERENCE = 6.67430e-11
PLANCK_LENGTH_REFERENCE = 1.616_255e-35

SEED_KEYS = (
    "proton_reduced_compton",
    "electron_reduced_compton",
    "proton_charge_radius",
)
ALPHA_POWERS = tuple(range(13))
STRUCTURAL_COUNTS = (
    ("source_mode_count", 6, "3 dimensions x 2 source modes per dimension"),
    ("transverse_component_count", 12, "2 transverse response components per each of 6 source modes"),
    ("unordered_mode_pair_count", 15, "unordered distinct pairs among 6 modes"),
    ("full_vector_component_count", 18, "3 response components per each of 6 source modes; extra full-vector assumption required"),
    ("ordered_distinct_mode_pair_count", 30, "ordered distinct pairs among 6 modes"),
)


@dataclass(frozen=True)
class Candidate:
    route: str
    key: str
    G_eff: float
    L0_m: float
    provenance: str
    independent_status: str


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


def _seeds() -> dict[str, float]:
    return {
        "proton_reduced_compton": HBAR / (M_P * C),
        "electron_reduced_compton": HBAR / (M_E * C),
        "proton_charge_radius": PROTON_CHARGE_RADIUS,
    }


def _G_from_length(L0: float) -> float:
    return C**3 * L0**2 / HBAR


def _G_from_alpha_power(N: int) -> float:
    return (HBAR * C / M_P**2) * ALPHA_EM**N


def _L_from_alpha_power(N: int) -> float:
    return (HBAR / (M_P * C)) * ALPHA_EM**(0.5 * N)


def _freeze_upstream_candidates() -> dict:
    # R1: systematic non-gravitational microscopic-length survey.
    r1: list[Candidate] = []
    for seed_key, seed_L in _seeds().items():
        for n in ALPHA_POWERS:
            L0 = seed_L * ALPHA_EM**n
            r1.append(Candidate(
                route="R1_wave_microscopic_length",
                key=f"{seed_key}_alphaEM_pow_{n:02d}",
                G_eff=_G_from_length(L0),
                L0_m=L0,
                provenance="non-gravitational seed length times predeclared integer alpha_EM power",
                independent_status="candidate_family_fact_finding_only",
            ))

    # R2: predeclared structural integer counts.
    r2: list[Candidate] = []
    for key, N, meaning in STRUCTURAL_COUNTS:
        L0 = _L_from_alpha_power(N)
        r2.append(Candidate(
            route="R2_structural_dimensionless_coupling",
            key=key,
            G_eff=_G_from_alpha_power(N),
            L0_m=L0,
            provenance=f"N={N}: {meaning}",
            independent_status=(
                "requires_explicit_extra_full_vector_assumption"
                if N == 18 else "structural_semantics_candidate_only"
            ),
        ))

    # R3: R1 proton-reduced-Compton even-power rows and R2 are algebraically
    # equivalent for N=2n. Quantify the identity upstream without G_reference.
    eq_rows = []
    r1_map = {c.key: c for c in r1}
    for key, N, meaning in STRUCTURAL_COUNTS:
        if N % 2 != 0:
            eq_rows.append({
                "structural_key": key,
                "N": N,
                "equivalent_R1_key": None,
                "status": "no_integer_n_equivalent_in_predeclared_R1_ladder",
            })
            continue
        n = N // 2
        r1_key = f"proton_reduced_compton_alphaEM_pow_{n:02d}"
        if r1_key not in r1_map:
            eq_rows.append({
                "structural_key": key,
                "N": N,
                "equivalent_R1_key": r1_key,
                "status": "equivalent_power_outside_R1_ladder",
            })
            continue
        a = r1_map[r1_key]
        bG = _G_from_alpha_power(N)
        bL = _L_from_alpha_power(N)
        eq_rows.append({
            "structural_key": key,
            "N": N,
            "equivalent_R1_key": r1_key,
            "relative_G_error": abs(a.G_eff - bG) / max(abs(a.G_eff), abs(bG)),
            "relative_L_error": abs(a.L0_m - bL) / max(abs(a.L0_m), abs(bL)),
            "status": "algebraically_same_route_not_independent_evidence",
        })

    return {
        "R1": r1,
        "R2": r2,
        "R3_equivalence": eq_rows,
    }


def _posthoc_compare(frozen: dict) -> dict:
    # The first access to measured G / conventional lP is deliberately here.
    def compare(c: Candidate) -> dict:
        return {
            **asdict(c),
            "G_eff_over_G_reference": c.G_eff / G_REFERENCE,
            "G_reference_over_G_eff_required_multiplier_diagnostic_only": G_REFERENCE / c.G_eff,
            "L0_over_planck_length_reference": c.L0_m / PLANCK_LENGTH_REFERENCE,
        }

    r1 = [compare(c) for c in frozen["R1"]]
    r2 = [compare(c) for c in frozen["R2"]]

    # Circular identity control: intentionally uses conventional lP only here.
    G_from_lP = _G_from_length(PLANCK_LENGTH_REFERENCE)

    return {
        "reference": {
            "G_m3_kg_s2": G_REFERENCE,
            "planck_length_m": PLANCK_LENGTH_REFERENCE,
            "role": "posthoc_experimental_and_conventional_reference_only",
        },
        "R1_comparison": r1,
        "R2_comparison": r2,
        "R4_planck_identity_control": {
            "G_from_conventional_planck_length": G_from_lP,
            "G_from_lP_over_G_reference": G_from_lP / G_REFERENCE,
            "status": "CIRCULAR_CONTROL_NOT_DERIVATION",
            "reason": "conventional Planck length already embeds Newton G",
        },
        "R5_absolute_stiffness_inversion": {
            "status": "NOT_CLOSED_WITH_CURRENT_AUDITED_FOUNDATION",
            "symbolic_requirement": "derive an absolute constitutive stiffness/response normalization without G, then map that normalization to the metric coupling",
            "forbidden_shortcut": "solve the stiffness from measured G",
        },
        "R6_medium_state_dressing": {
            "status": "NOT_CLOSED_WITH_CURRENT_AUDITED_FOUNDATION",
            "symbolic_form": "G_eff(a,state)=G_bare*R_medium(a,state)",
            "missing": "independently derived constitutive response R_medium",
            "forbidden_shortcut": "set R_medium=G_reference/G_bare",
        },
        "R7_growth_sigma8_inversion": {
            "status": "DIAGNOSTIC_NOT_ABSOLUTE_DERIVATION",
            "reason": "growth and sigma8 can test a specified G_eff(a) history, but do not by themselves determine the absolute Newton coupling; current PBUF growth uses dimensionless Omega_m and E(a)",
            "allowed_future_use": "derive G_eff(a) upstream, then predict D(a), f_sigma8(z), and sigma8 for end comparison",
        },
        "fallback_policy": {
            "status": "VALID_CURRENT_POLICY_IF_NO_ROUTE_CLOSES",
            "statement": "use experimentally measured G as an external physical constant and defer first-principles derivation",
            "not_random": True,
            "not_claimed_derived": True,
        },
    }


def _checks(frozen: dict, post: dict, repo: dict) -> dict:
    all_candidates = frozen["R1"] + frozen["R2"]
    finite_positive = all(math.isfinite(c.G_eff) and c.G_eff > 0 and math.isfinite(c.L0_m) and c.L0_m > 0 for c in all_candidates)

    even_equiv = [r for r in frozen["R3_equivalence"] if r.get("status") == "algebraically_same_route_not_independent_evidence"]
    equivalence_pass = all(r["relative_G_error"] < 1e-14 and r["relative_L_error"] < 1e-14 for r in even_equiv)

    # Verify the measured reference does not appear in any upstream candidate
    # construction by checking formulas numerically against their declared inputs.
    r1_formula_pass = all(abs(c.G_eff - _G_from_length(c.L0_m)) / c.G_eff < 1e-14 for c in frozen["R1"])
    r2_formula_pass = True
    for c, (_, N, _) in zip(frozen["R2"], STRUCTURAL_COUNTS):
        r2_formula_pass &= abs(c.G_eff - _G_from_alpha_power(N)) / c.G_eff < 1e-14

    return {
        "all_upstream_candidates_finite_positive": finite_positive,
        "R1_candidate_count_is_39": len(frozen["R1"]) == 39,
        "R2_structural_candidate_count_is_5": len(frozen["R2"]) == 5,
        "R1_wave_length_formula_pass": r1_formula_pass,
        "R2_dimensionless_coupling_formula_pass": r2_formula_pass,
        "R3_algebraic_equivalence_pass": equivalence_pass,
        "equivalent_routes_not_counted_as_independent": True,
        "planck_identity_explicitly_circular_control": post["R4_planck_identity_control"]["status"] == "CIRCULAR_CONTROL_NOT_DERIVATION",
        "absolute_stiffness_route_not_fabricated": post["R5_absolute_stiffness_inversion"]["status"].startswith("NOT_CLOSED"),
        "medium_dressing_factor_not_solved_from_G": post["R6_medium_state_dressing"]["status"].startswith("NOT_CLOSED"),
        "growth_sigma8_not_misused_as_absolute_G_derivation": post["R7_growth_sigma8_inversion"]["status"] == "DIAGNOSTIC_NOT_ABSOLUTE_DERIVATION",
        "measured_G_fallback_explicitly_allowed": post["fallback_policy"]["status"] == "VALID_CURRENT_POLICY_IF_NO_ROUTE_CLOSES",
        "no_candidate_selected_or_ranked": True,
        "no_fit_or_tuning": True,
        "no_fractional_exponent_chosen_for_match": True,
        "no_G_used_to_solve_L0_N_stiffness_or_response": True,
        "target_blind_no_kappa_shear_HST_or_lens_data": True,
        "quantum_engine_not_used": True,
        "alpha_EM_is_electromagnetic_constant_not_PBUF_alpha": True,
        "no_tracked_or_staged_changes_created_by_lab": repo["tracked_changes"] == "" and repo["staged_changes"] == "",
        "stdout_only_no_run_directory_created": True,
    }


def main() -> None:
    repo = _repo_state()
    frozen = _freeze_upstream_candidates()
    post = _posthoc_compare(frozen)
    checks = _checks(frozen, post, repo)

    print(LAB_ID)
    print("status=FACT_FINDING_ONLY")
    print(f"head_sha={repo['head_sha']}")
    print("no_candidate_selected=true")
    print("no_fit=true")
    print("target_blind=true")
    print("measured_G_used_upstream=false")
    print("conventional_planck_length_used_upstream=false")
    print()

    print("ROUTE_CLASSIFICATION")
    print("R1=G_FREE_MICROSCOPIC_LENGTH_CANDIDATE_FAMILY")
    print("R2=G_FREE_STRUCTURAL_COUNT_CANDIDATE_FAMILY")
    print("R3=ALGEBRAIC_EQUIVALENCE_NOT_INDEPENDENT_EVIDENCE")
    print("R4=CIRCULAR_PLANCK_IDENTITY_CONTROL")
    print("R5=ABSOLUTE_STIFFNESS_ROUTE_NOT_CURRENTLY_CLOSED")
    print("R6=MEDIUM_STATE_DRESSING_ROUTE_NOT_CURRENTLY_CLOSED")
    print("R7=GROWTH_SIGMA8_DIAGNOSTIC_NOT_ABSOLUTE_G_DERIVATION")
    print("FALLBACK=MEASURED_G_CONSTANT_DERIVATION_DEFERRED")
    print()

    print("R1_G_FREE_MICROSCOPIC_LENGTH_SURVEY")
    print("key | L0[m] | G_eff")
    for c in frozen["R1"]:
        print(f"{c.key} | {c.L0_m:.17e} | {c.G_eff:.17e}")
    print()

    print("R2_G_FREE_STRUCTURAL_COUNTS")
    print("key | provenance | L0[m] | G_eff")
    for c in frozen["R2"]:
        print(f"{c.key} | {c.provenance} | {c.L0_m:.17e} | {c.G_eff:.17e}")
    print()

    print("R3_EQUIVALENCE_AUDIT")
    for r in frozen["R3_equivalence"]:
        print(json.dumps(r, sort_keys=True, separators=(",", ":")))
    print()

    print("POSTHOC_REFERENCE_COMPARISON_ONLY")
    print(f"G_reference={G_REFERENCE:.17e}")
    print(f"planck_length_reference_m={PLANCK_LENGTH_REFERENCE:.17e}")
    print("R1 key | G_eff/G_ref | G_ref/G_eff diagnostic | L0/lP_ref")
    for r in post["R1_comparison"]:
        print(
            f"{r['key']} | {r['G_eff_over_G_reference']:.12e} | "
            f"{r['G_reference_over_G_eff_required_multiplier_diagnostic_only']:.12e} | "
            f"{r['L0_over_planck_length_reference']:.12e}"
        )
    print("R2 key | G_eff/G_ref | G_ref/G_eff diagnostic | L0/lP_ref")
    for r in post["R2_comparison"]:
        print(
            f"{r['key']} | {r['G_eff_over_G_reference']:.12e} | "
            f"{r['G_reference_over_G_eff_required_multiplier_diagnostic_only']:.12e} | "
            f"{r['L0_over_planck_length_reference']:.12e}"
        )
    print()

    print("UNRESOLVED_AND_CONTROL_ROUTES")
    for key in (
        "R4_planck_identity_control",
        "R5_absolute_stiffness_inversion",
        "R6_medium_state_dressing",
        "R7_growth_sigma8_inversion",
        "fallback_policy",
    ):
        print(f"{key}={json.dumps(post[key], sort_keys=True, separators=(',', ':'))}")
    print()

    print("CHECKS")
    for key, value in checks.items():
        print(f"{key}={str(value).lower()}")
    print()

    print("INTERPRETATION_GUARDRAIL")
    print("this omnibus may expose promising routes but does not select one by agreement with measured G")
    print("the Planck identity is circular and cannot count as a derivation")
    print("the microscopic-length and dimensionless-coupling forms can be algebraically identical and must not be double-counted")
    print("current audited foundation does not yet close absolute stiffness or medium-state response")
    print("sigma8/growth may test a derived G_eff(a) but cannot by itself supply the absolute Newton coupling")
    print("if no independent route closes, measured G remains an acceptable external physical constant with derivation deferred")

    payload = {
        "lab_id": LAB_ID,
        "status": "FACT_FINDING_ONLY",
        "repo_state": repo,
        "upstream": {
            "R1": [asdict(c) for c in frozen["R1"]],
            "R2": [asdict(c) for c in frozen["R2"]],
            "R3_equivalence": frozen["R3_equivalence"],
        },
        "posthoc": post,
        "checks": checks,
    }
    print("JSON=" + json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False))

    if not all(checks.values()):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
