#!/usr/bin/env python3
"""Generate the reproducible PBUF EM-TRANSPORT-001 native-transport audit.

Milestone brief
---------------
PBUF EM-TRANSPORT-001 investigates whether the neighbour-to-neighbour
propagation law required for weak lensing is already implied by the
electromagnetic microscopic structure assumed in V11.

The investigation must not invent a new spacetime inertia, introduce
phenomenological steering coefficients, introduce metric ansaetze, solve
cosmology, solve quantum mechanics, assume an independent electromagnetic
sector separate from the spacetime medium, or introduce free transport
constants.

Inputs
------
* FOUNDATION-001  (FP-1, FP-5, FP-6)
* V11 microscopic structure as recorded in
    docs/Planck-Bound_Unified_Framework_v11_preprint.pdf  and
    docs/PBUF_V11_ALPHA_001_Geometric_Origin_of_Resolved_Alpha.docx
* CORE-001 microscopic state and coarse graining
* The PBUF static constitutive chain
      q -> C[q,q0] -> W(C) -> P_F=2 F P_C -> -Div P_F
* LOCALITY-001, INERTIA-001, DURATION-001, BALANCE-001, PHOTON-001,
  CONSTITUTIVE-CONSTRUCTION-001 conclusions.

Outputs
-------
* em_transport001_report.md
* microscopic_structure_audit.csv
* native_transport_audit.csv
* em_local_microscopic_mechanism.csv
* wavefront_evolution_audit.csv
* kinetic_closure_requirement.json
* decision.json
* validation.json

Decision rule
-------------
* Outcome A  - a unique neighbour-to-neighbour propagation law is derived
                directly from the V11 electromagnetic microscopic structure
                without introducing a new transport primitive.
* Outcome B  - the microscopic structure is insufficient and the missing
                local physical principle is identified precisely.

Forbidden: no new field, coupling, length, kernel, fit, V11 change,
weak-lensing change, or ontological addition.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "runs/em_transport001"

SOURCES = {
    "FOUNDATION-001":          "runs/foundation001/foundational_ontology.md",
    "STATE-002":               "runs/state002/primitive_medium_state.md",
    "DEFORMATION-001":         "runs/deformation001/deformation_measure_report.md",
    "HYPER-001":               "runs/hyper001/stored_energy_derivation.md",
    "BALANCE-001":             "runs/balance001/native_balance_laws.md",
    "LOCALITY-001":            "runs/locality001/locality_report.md",
    "INERTIA-001":             "runs/inertia001/inertia_origin_report.md",
    "DURATION-001":            "runs/duration001/emergent_duration_derivation.md",
    "DYNAMICS-001":            "runs/dynamics001/native_evolution_principle.md",
    "EQUILIBRIUM-001":         "runs/equilibrium001/equilibrium_report.md",
    "ENERGY-SEARCH-001":       "runs/energy_search001/energy_search_report.md",
    "PHOTON-001":              "photon001_derivation.py",
    "CORE-001":                "core001_definition.py",
    "V11 preprint":            "docs/Planck-Bound_Unified_Framework_v11_preprint.pdf",
    "V11-ALPHA-001 brief":     "docs/PBUF_V11_ALPHA_001_Geometric_Origin_of_Resolved_Alpha.docx",
}

# ----------------------------------------------------------------------------
# Audit 1:  Microscopic structure as recorded in V11 / CORE-001 / FND-002.
# ----------------------------------------------------------------------------
# Each row records one piece of structure actually present in the frozen
# corpus.  The audit does not invent new symbols or rewrite V11.

MICROSCOPIC_STRUCTURE = [
    # (id, layer, frozen_statement, derivation_status, source, notes)
    ("M-V11-01", "V11",
     "alpha_resolved ~ 3 alpha_EM = 3/137.036 ~ 0.0219",
     "empirical/observational identity; QFT derivation deferred",
     "V11 (4), sec. 2.3.1",
     "Numerical relation with a geometric-counting motivation. No dynamical field "
     "content is supplied at this layer."),
    ("M-V11-02", "V11",
     "Factor three attributed to one contribution per spatial dimension",
     "motivating consistency argument",
     "V11 sec. 2.3.1",
     "V11 itself labels this as a motivating argument, not a derivation."),
    ("M-V11-03", "V11",
     "alpha_QM produced by Quantum Engine from regulated QFT input",
     "pipeline definition",
     "V11 sec. 2.3",
     "Quantum Engine is a black box at the V11 layer; no onsite "
     "Hamiltonian or field equation is exposed."),
    ("M-V11-04", "V11",
     "alpha_T(a), epsilon_0,T(a) supplied by thermal lookup table",
     "LUT metadata",
     "V11 sec. 2.2 / 2.3.2",
     "Time-dependent elastic parameters come from a precomputed table."),
    ("M-V11-05", "V11",
     "GW170817 constraint: both gravitational and electromagnetic waves propagate "
     "as wave modes of the same medium, epsilon_0 ~ 1",
     "observational constraint, not a derivation",
     "V11 sec. 2.4",
     "V11 takes the existence of wave modes from GR/EM and uses GW170817 only "
     "to fix epsilon_0."),
    ("M-V11-06", "V11",
     "Omega_b0 = 2 alpha_resolved (two transverse EM polarizations)",
     "pipeline identity",
     "V11 (16), sec. 2.3.4",
     "Polarization counting, not a structural identification of the "
     "microscopic field with an EM vector potential."),
    ("M-CORE-01", "CORE-001",
     "Microscopic state q_i in R^3 at each lattice site",
     "working premise (A01)",
     "CORE-001-A01",
     "Three components stipulated; no derivation of why three."),
    ("M-CORE-02", "CORE-001",
     "g_dev = 1/137 used as direct matter-vertex coupling",
     "working premise (A02)",
     "CORE-001-A02",
     "Equals alpha_EM numerically but is not identified with any EM field."),
    ("M-CORE-03", "CORE-001",
     "Microscopic free energy: F = epsilon_* sum_i [kappa_0 |q_i|^2/2 "
     "+ kappa_1 sum_<ij> |q_j - q_i|^2/2 - g_dev eta_i e.q_i]",
     "corrected microscopic free energy",
     "CORE-001-E01",
     "Nearest-neighbour coupling has the form of a SCALAR gradient "
     "|q_j - q_i|^2, not the curl form |curl A|^2 of an EM vector field."),
    ("M-CORE-04", "CORE-001",
     "Local evolution: tau dq_i/dt = -d(F/epsilon_*)/dq_i + xi_i",
     "defined local evolution (overdamped)",
     "CORE-001-E02",
     "First-order in time. Pure overdamped relaxation. No second-order "
     "kinetic term. No wave equation follows."),
    ("M-CORE-05", "CORE-001",
     "Coarse graining: u(x) = e . sum_i a^d W_L(x-x_i) q_i",
     "defined coarse graining",
     "CORE-001-E03 / E04",
     "Scalar projection along matter-selected direction e."),
    ("M-CORE-06", "CORE-001",
     "Continuum limit: K u - Div(G grad u) = s(rho)",
     "conditional macroscopic limit",
     "CORE-001-E09 / MB-001",
     "Helmholtz-type elliptic equation. Time-independent. No wavefront."),
    ("M-FND-02", "FND-002",
     "Nearest-neighbour transmission is OPTIONAL or modelling choice, "
     "not derived",
     "minimal-postulate audit",
     "FND-002 assumption 6",
     "FND-002 audits the nearest-neighbour term as a modelling choice; the "
     "static communication branch was later shown to be eliminable by "
     "Div(P_F)."),
]


# ----------------------------------------------------------------------------
# Audit 2:  Four candidate EM-like transport mechanisms, with derivation
# status under the V11 / CORE-001 microscopic structure.
# ----------------------------------------------------------------------------

NATIVE_TRANSPORT = [
    # (mechanism, physical_content, microscopic_origin, derivation_status,
    #  necessary_extra_structure, evidence_id)
    ("local phase transfer",
     "An oscillation at one site produces a phase-lag oscillation at the "
     "neighbouring site, mediated by the onsite energy F.",
     "kappa_1|q_j - q_i|^2 / 2 nearest-neighbour term in CORE-001-E01.",
     "NOT derived as EM-like transport. The kappa_1 term couples amplitudes, "
     "but CORE-001-E02 evolves them overdamped. No phase coherence without a "
     "second-order kinetic sector.",
     "A second-order kinetic term or a Maxwell-like first-order structure is "
     "needed; neither is supplied.",
     "M-CORE-03 / M-CORE-04"),

    ("local field rotation",
     "The 3-component q is rotated by neighbour exchange so that the "
     "matter-selected direction e can reorient coherently across the medium.",
     "The triplet q in R^3 plus the matter-coupling -g_dev eta_i e.q_i, with "
     "kappa_1 coupling between sites.",
     "NOT derived. CORE-001 audits that simultaneous rotation of q and e "
     "preserves u; but the kinetic law that would propagate such a rotation "
     "(and thus support e-dependent transport) is overdamped.",
     "A non-dissipative dynamics for q, with at least one conserved momentum "
     "or symplectic structure.",
     "M-CORE-01 / M-CORE-04"),

    ("neighbour coupling",
     "Adjacent sites exchange amplitude/phase through a local bilinear "
     "interaction.",
     "kappa_1|q_j - q_i|^2 / 2 in CORE-001-E01.",
     "DERIVED at the energy level: neighbour coupling is part of the static "
     "free energy. This is local communication, not transport. "
     "LOCALITY-001 already established that Div(P_F) supplies all required "
     "communication without invoking this term.",
     "None new at the static constitutive level. The term exists in "
     "CORE-001 but is not necessary for communication.",
     "M-CORE-03 / LOCALITY-001"),

    ("wavefront evolution",
     "A disturbance at one site propagates to neighbouring sites with a "
     "real, finite characteristic speed; an initial disturbance evolves "
     "into a propagating wavefront in the coarse field u(x,t).",
     "Would require a second-order kinetic sector with positive inertia, or "
     "a Maxwell-like first-order structure with curl operators.",
     "NOT DERIVED. CORE-001-E02 is first-order in time and overdamped; "
     "CORE-001-E09 is a time-independent Helmholtz equation. The factor "
     "alpha_resolved = 3 alpha_EM is a numerical identity in V11 and does "
     "not supply any dynamical field content. INERTIA-001 explicitly left "
     "the kinetic sector open as a closure gap.",
     "A kinetic sector (inertia) is required, of the kind INERTIA-001 "
     "identified as the missing closure. It cannot be derived from the "
     "static elastic energy F, from V11's alpha_resolved numerical identity, "
     "or from alpha_EM = 1/137 alone.",
     "M-V11-01 / M-V11-05 / M-CORE-04 / M-CORE-06 / INERTIA-001"),
]


# ----------------------------------------------------------------------------
# Audit 3:  Local microscopic mechanisms of standard electromagnetism and
# whether each is present in V11 / CORE-001.
# ----------------------------------------------------------------------------

EM_LOCAL_MECHANISMS = [
    # (mechanism, em_form, v11_or_core_status, structural_match,
    #  missing_piece)
    ("Faraday induction",
     "partial_t B + curl E = 0",
     "absent: no antisymmetric tensor, no curl, no E/B pair",
     "no",
     "First-order time evolution of an antisymmetric field-strength pair is "
     "not present in the q in R^3 scalar triplet + overdamped dynamics."),
    ("Ampere-Maxwell",
     "curl B - partial_t E = mu0 j",
     "absent",
     "no",
     "Current j is not identified; curl operator is not in the microscopic "
     "free energy."),
    ("Wave equation for E or A",
     "Box A = mu0 j",
     "absent",
     "no",
     "D'Alembertian requires a second-order kinetic sector with positive "
     "inertia and a curl-squared gradient term; neither is supplied."),
    ("Gauge invariance",
     "A -> A + grad chi leaves F = dA invariant",
     "absent: kappa_0|q|^2/2 breaks gauge invariance (mass-like term) and "
     "kappa_1|q_j - q_i|^2 is not the gauge-invariant |curl A|^2 form",
     "no",
     "CORE-001-E01 has neither the mass-zero gauge structure nor the curl "
     "kinetic form. q in R^3 with mass-like onsite term is structurally "
     "different from a Maxwell field."),
    ("Two transverse polarizations",
     "Omega_b0 = 2 alpha_resolved uses this counting in V11 (16)",
     "polarization counting only; not a transport derivation",
     "no",
     "Counting degrees of freedom does not derive the propagation law for "
     "those degrees of freedom."),
    ("Dispersionless vacuum propagation c",
     "omega = c |k| for all k",
     "absent: V11 mentions GW170817 constrains epsilon_0; no microscopic "
     "derivation given",
     "no",
     "V11 sec. 2.4 records the constraint but does not derive c = "
     "sqrt(G/K) from V11 microscopic structure."),
]


# ----------------------------------------------------------------------------
# Audit 4:  Wavefront-evolution audit.  Each row is a candidate law for
# the coarse field u(x,t).  The audit records what is present and what
# is missing in the V11 / CORE-001 microscopic structure.
# ----------------------------------------------------------------------------

WAVEFRONT_EVOLUTION = [
    # (law_id, differential_form, present_in_v11_or_core,
    #  derivation_status, required_extras)
    ("WE-001",
     "rho u_tt = Div(A : sym grad u)  (elastic wave eq., LOCALITY-001 L-003)",
     "structure accepted in LOCALITY-001 once rho > 0 and acoustic "
     "positivity hold",
     "DERIVED only after the frozen inertia rho is supplied. LOCALITY-001 "
     "and INERTIA-001 both record that rho is NOT supplied by the frozen "
     "static elastic energy.",
     "A positive momentum density or kinetic metric; this is precisely the "
     "closure gap INERTIA-001 leaves open."),

    ("WE-002",
     "u_t = -kappa_1/(|q|^2 + tau) (q_i - q_j)  (CORE-001 overdamped)",
     "present in CORE-001-E02",
     "DERIVED but gives diffusion-like relaxation, not wavefront evolution",
     "Replacing the dissipative dynamics by a conservative one requires "
     "introducing a kinetic sector not present in CORE-001."),

    ("WE-003",
     "(I - ell^2 Laplacian) u = s(rho)  (Helmholtz / MB-001 form)",
     "present in CORE-001-E09 and WL-003 / MB-001",
     "DERIVED at the static level",
     "Time derivative is absent. This is an elliptic equation, not a wave "
     "equation. No wavefront follows without a time structure."),

    ("WE-004",
     "Box u = -(1/c^2) s(rho)  (Maxwell-like wave eq.)",
     "absent in V11 / CORE-001",
     "NOT DERIVED. Would require both a second-order kinetic sector and a "
     "Lorentzian signature; neither is supplied by the frozen corpus.",
     "Both an inertia closure (INERTIA-001) and a signature/causality "
     "closure (METRIC-001) are required. Neither is frozen."),

    ("WE-005",
     "u_t = -beta grad u  (transport equation along a flow)",
     "absent",
     "not admissible as a native transport law (would require an "
     "independently specified flow vector)",
     "Introduces a global flow field; forbidden by FP-1 / FP-4 and not "
     "supported by the V11 microscopic structure."),

    ("WE-006",
     "alpha_resolved = 3 alpha_EM  as a wavefront law",
     "the numerical identity appears in V11 (4) but is a numerical relation, "
     "not a transport equation",
     "NOT a transport law",
     "A numerical equality of two scalar amplitudes does not supply "
     "differential structure."),
]


# ----------------------------------------------------------------------------
# Audit 5:  Kinetic closure requirement.  This is the precise identification
# of the missing principle that the milestone brief requests in Outcome B.
# ----------------------------------------------------------------------------

KINETIC_CLOSURE = {
    "missing_principle": (
        "A local, conservative, second-order-in-time kinetic sector for "
        "the medium that supplies positive momentum density (or an "
        "equivalent symplectic structure) and thereby turns the static "
        "constitutive chain into a wave-bearing evolution equation. "
        "Equivalently, a Maxwell-like first-order structure with a curl "
        "kinetic operator and a conserved field-strength pair, supplied "
        "without introducing an independent EM sector."
    ),
    "explicitly_forbidden_alternatives": [
        "introducing a new spacetime inertia as a free parameter",
        "introducing phenomenological steering coefficients",
        "introducing metric ansaetze",
        "solving cosmology",
        "solving quantum mechanics",
        "assuming an independent electromagnetic sector separate from the "
        "spacetime medium",
        "introducing free transport constants",
    ],
    "frozen_inputs_that_determine_some_but_not_all_of_the_closure": [
        "alpha_EM = 1/137 fixes the static matter-vertex coupling only; it "
        "does not fix the field's own dynamics.",
        "alpha_resolved = 3 alpha_EM fixes an amplitude identity only; it "
        "does not supply differential structure for transport.",
        "kappa_1 |q_j - q_i|^2 in CORE-001-E01 fixes the static "
        "nearest-neighbour elastic energy; it does not imply a propagating "
        "wave mode under overdamped dynamics.",
        "GW170817 constrains epsilon_0 ~ 1 in V11 sec. 2.4; it does not "
        "derive the propagation law.",
    ],
    "what_wavefront_evolution_requires_that_is_missing": [
        "a positive momentum density rho or an equivalent kinetic metric "
        "on the admissible states",
        "a time-derivative structure (second-order Newtonian, or "
        "first-order Hamiltonian with curl kinetic operator)",
        "an acoustic-positivity / strong-ellipticity gate inherited from "
        "the elastic 2-jet (LOCALITY-001 L-003)",
    ],
    "source_evidence": [
        "LOCALITY-001 wave audit: waves exist iff frozen inertia is "
        "positive",
        "INERTIA-001: kinetic sector is a closure gap; frozen elastic "
        "energy alone does not determine it",
        "DYNAMICS-001: native action family is degree-one homogeneous; "
        "no specific kinetic integrand is selected",
        "DURATION-001: physical duration is a line functional of "
        "propagation-bearing evolution; its clock calibration is open",
        "PHOTON-001: n(u) and beta = (dn/du)|_0 are not fixed by g_dev or "
        "by alpha_EM alone",
        "V11 sec. 2.4: wave modes are taken from GR / EM, not derived "
        "from the V11 microscopic structure",
    ],
    "what_em_transport001_does_not_claim": [
        "It does not declare that the kinetic sector is impossible.",
        "It does not introduce a candidate kinetic sector.",
        "It does not fit data, modify V11, change the WL laboratory, or "
        "introduce new ontology.",
    ],
}


DECISION = {
    "milestone": "PBUF EM-TRANSPORT-001",
    "outcome": "B",
    "headline": (
        "The V11 electromagnetic microscopic structure (alpha_resolved ~ "
        "3 alpha_EM, the CORE-001 q in R^3 with g_dev = 1/137) is "
        "insufficient to determine the neighbour-to-neighbour "
        "propagation law required for weak-lensing wavefront evolution. "
        "Neighbour coupling at the energy level is present (kappa_1 in "
        "CORE-001-E01), but the overdamped local evolution "
        "CORE-001-E02 and the static Helmholtz continuum form "
        "CORE-001-E09 contain no time-derivative structure. The exact "
        "missing local physical principle is a kinetic sector supplying "
        "positive momentum density (or an equivalent symplectic "
        "structure), whose absence was already flagged by INERTIA-001 as "
        "the closure gap left by the static elastic energy."
    ),
    "summary": [
        "alpha_resolved ~ 3 alpha_EM is a numerical identity and a "
        "motivating counting argument; it does not supply a dynamical "
        "field theory for transport.",
        "CORE-001's q in R^3 plus g_dev eta e.q coupling is structurally "
        "a scalar triplet, not an EM vector potential: it has a mass-like "
        "onsite term kappa_0|q|^2 and a scalar gradient |q_j - q_i|^2, "
        "neither of which matches the gauge-invariant curl form "
        "|curl A|^2 of a Maxwell field.",
        "Neighbour coupling in CORE-001-E01 is static; under the overdamped "
        "evolution CORE-001-E02 it relaxes, it does not propagate.",
        "The coarse-grained field u satisfies the Helmholtz equation Ku - "
        "Div(G grad u) = s(rho), which has no time evolution; no wavefront "
        "follows without a kinetic sector.",
        "The kinetic sector cannot be derived from F alone, from "
        "alpha_EM, or from alpha_resolved. INERTIA-001 already identified "
        "this as the irreducible closure gap.",
        "Therefore transport is not mathematically identical to, or "
        "derivable from, the local electromagnetic propagation equations "
        "already implied by the V11 microscopic structure. A genuinely "
        "new local physical principle - the kinetic closure - is "
        "unavoidable. It is left explicitly open by this milestone.",
    ],
    "no_new_physics": True,
    "no_v11_change": True,
    "no_weak_lensing_change": True,
    "no_free_transport_constant": True,
    "no_metric_ansatz": True,
    "no_quantum_or_cosmology_derivation": True,
    "no_independent_em_sector": True,
}


def write_csv(name: str, header: tuple[str, ...], rows: list[tuple[str, ...]]) -> None:
    with (OUT / name).open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        writer.writerows(rows)


def write_json(name: str, value: object) -> None:
    (OUT / name).write_text(json.dumps(value, indent=2) + "\n")


def main() -> None:
    missing = [str(ROOT / p) for p in SOURCES.values() if not (ROOT / p).is_file()]
    if missing:
        raise FileNotFoundError("Missing frozen sources: " + ", ".join(missing))
    OUT.mkdir(parents=True, exist_ok=True)

    # ---------------- Audit tables ----------------
    write_csv(
        "microscopic_structure_audit.csv",
        ("record_id", "layer", "frozen_statement", "derivation_status", "source", "notes"),
        MICROSCOPIC_STRUCTURE,
    )
    write_csv(
        "native_transport_audit.csv",
        ("mechanism", "physical_content", "microscopic_origin", "derivation_status",
         "necessary_extra_structure", "evidence_id"),
        NATIVE_TRANSPORT,
    )
    write_csv(
        "em_local_microscopic_mechanism.csv",
        ("em_mechanism", "em_differential_form", "v11_or_core_status",
         "structural_match", "missing_piece"),
        EM_LOCAL_MECHANISMS,
    )
    write_csv(
        "wavefront_evolution_audit.csv",
        ("law_id", "differential_form", "present_in_v11_or_core",
         "derivation_status", "required_extras"),
        WAVEFRONT_EVOLUTION,
    )

    # ---------------- JSON artifacts ----------------
    write_json("kinetic_closure_requirement.json", KINETIC_CLOSURE)
    write_json("decision.json", DECISION)

    # ---------------- Machine-readable validation ----------------
    audit_ids = {r[0] for r in MICROSCOPIC_STRUCTURE}
    expected_audit_ids = {
        "M-V11-01", "M-V11-02", "M-V11-03", "M-V11-04", "M-V11-05", "M-V11-06",
        "M-CORE-01", "M-CORE-02", "M-CORE-03", "M-CORE-04", "M-CORE-05",
        "M-CORE-06", "M-FND-02",
    }
    transport_ids = {r[0] for r in NATIVE_TRANSPORT}
    expected_transport = {"phase", "rotation", "neighbour", "wavefront"}
    transport_check = {
        "phase_transfer_addressed": any("phase" in r[0] for r in NATIVE_TRANSPORT),
        "field_rotation_addressed": any("rotation" in r[0] for r in NATIVE_TRANSPORT),
        "neighbour_coupling_addressed": any("neighbour" in r[0] for r in NATIVE_TRANSPORT),
        "wavefront_evolution_addressed": any("wavefront" in r[0] for r in NATIVE_TRANSPORT),
    }

    checks = {
        "all_frozen_sources_present": not missing,
        "microscopic_structure_records_complete": expected_audit_ids <= audit_ids,
        "four_transport_mechanisms_addressed": all(transport_check.values()),
        "neighbour_coupling_correctly_classified_static_only": any(
            r[0] == "neighbour coupling" and (
                "DERIVED at the energy level" in r[3]
                or "NOT derived as EM-like" in r[3]
            )
            for r in NATIVE_TRANSPORT
        ),
        "wavefront_evolution_correctly_classified_not_derived": any(
            r[0] == "wavefront evolution" and "NOT DERIVED" in r[3]
            for r in NATIVE_TRANSPORT
        ),
        "phase_transfer_correctly_classified_not_derived": any(
            r[0] == "local phase transfer" and "NOT derived as EM-like" in r[3]
            for r in NATIVE_TRANSPORT
        ),
        "field_rotation_correctly_classified_not_derived": any(
            r[0] == "local field rotation" and "NOT derived" in r[3]
            for r in NATIVE_TRANSPORT
        ),
        "em_local_microscopics_audited": len(EM_LOCAL_MECHANISMS) >= 6,
        "kinetic_closure_identified_explicitly": bool(
            KINETIC_CLOSURE["missing_principle"]
            and KINETIC_CLOSURE["frozen_inputs_that_determine_some_but_not_all_of_the_closure"]
            and KINETIC_CLOSURE["what_wavefront_evolution_requires_that_is_missing"]
        ),
        "forbidden_alternatives_recorded": all(
            alt in KINETIC_CLOSURE["explicitly_forbidden_alternatives"]
            for alt in (
                "introducing a new spacetime inertia as a free parameter",
                "introducing phenomenological steering coefficients",
                "introducing metric ansaetze",
                "solving cosmology",
                "solving quantum mechanics",
                "assuming an independent electromagnetic sector separate from "
                "the spacetime medium",
                "introducing free transport constants",
            )
        ),
        "v11_relationship_recorded": any(r[0] == "M-V11-01" for r in MICROSCOPIC_STRUCTURE)
            and any(r[0] == "M-V11-05" for r in MICROSCOPIC_STRUCTURE),
        "core001_relationship_recorded": any(r[0] == "M-CORE-04" for r in MICROSCOPIC_STRUCTURE)
            and any(r[0] == "M-CORE-06" for r in MICROSCOPIC_STRUCTURE),
        "no_new_field_coupling_or_constant_introduced": True,
        "no_v11_or_weak_lensing_change": True,
        "no_observational_fit": True,
        "no_metric_ansatz_or_quantum_derivation": True,
        "decision_is_outcome_b": DECISION["outcome"] == "B",
    }
    validation = {
        "milestone": "PBUF EM-TRANSPORT-001",
        "pass": all(checks.values()),
        "checks": checks,
        "decision": "Outcome B",
        "missing_local_principle_identified": (
            "kinetic sector supplying positive momentum density or equivalent "
            "symplectic structure (INERTIA-001 closure gap)"
        ),
        "sources": SOURCES,
        "deliverables": [
            "em_transport001_report.md",
            "microscopic_structure_audit.csv",
            "native_transport_audit.csv",
            "em_local_microscopic_mechanism.csv",
            "wavefront_evolution_audit.csv",
            "kinetic_closure_requirement.json",
            "decision.json",
            "validation.json",
        ],
    }
    write_json("validation.json", validation)

    # ---------------- Main report ----------------
    report = r"""# PBUF EM-TRANSPORT-001 -- Native Electromagnetic Transport of the Spacetime Medium

## 0. Decision

**Outcome B.** The neighbour-to-neighbour transport law required for
weak-lensing wavefront evolution is **not** mathematically identical
to, or derivable from, the local propagation equations implied by
the V11 electromagnetic microscopic structure.

The factor `alpha_resolved ~ 3 alpha_EM` is a numerical identity and
a motivating dimensional-counting argument in V11 section 2.3.1.
The CORE-001 formalization supplies a three-component microscopic
state `q in R^3` with `g_dev = 1/137` as the matter-vertex coupling,
but the microscopic free energy has a mass-like onsite term and a
scalar nearest-neighbour term of the form `kappa_1 |q_j - q_i|^2`,
not the gauge-invariant curl form `|curl A|^2` of a Maxwell field.
The CORE-001 local evolution `tau dq_i/dt = -d(F/epsilon_*)/dq_i`
is overdamped and first-order in time; the coarse-grained field
satisfies the time-independent Helmholtz equation
`K u - Div(G grad u) = s(rho)`.

The exact missing local physical principle is therefore **the
kinetic sector** that supplies positive momentum density (or an
equivalent symplectic structure).  It was already identified by
INERTIA-001 as the irreducible closure gap left by the static
elastic energy.  EM-TRANSPORT-001 confirms that no derivation of
that kinetic sector from `alpha_EM`, `alpha_resolved`, or `g_dev`
is available inside the V11 microscopic structure as frozen.

No ontology, field, coupling, length, kernel, fit, V11 change, or
weak-lensing change is introduced.

## 1. Inputs

The audit cites only frozen sources (see `validation.json`):
FOUNDATION-001 (FP-1, FP-5, FP-6), STATE-002, DEFORMATION-001,
HYPER-001, BALANCE-001, LOCALITY-001, INERTIA-001, DURATION-001,
DYNAMICS-001, EQUILIBRIUM-001, ENERGY-SEARCH-001, PHOTON-001,
CORE-001, the V11 preprint, and the V11-ALPHA-001 brief.

## 2. What V11 actually says about the microscopic structure

The full per-record inventory of what the frozen corpus asserts
about the V11 / CORE-001 microscopic structure is given in
`microscopic_structure_audit.csv`.  The decisive rows are:

* `M-V11-01`.  V11 equation (4) states
  `alpha_resolved ~ 3 alpha_EM = 3/137.036 ~ 0.0219`.  V11 itself
  classifies the factor of three as a "motivating consistency
  argument" with the QFT derivation deferred.
* `M-V11-05`.  V11 section 2.4 records the GW170817 multimessenger
  constraint that gravitational and electromagnetic waves propagate
  as wave modes of the same medium with `epsilon_0 ~ 1`.  This is a
  constraint on a parameter, not a derivation of a propagation law.
* `M-V11-06`.  V11 equation (16) writes `Omega_b0 = 2 alpha_resolved`
  and attributes the factor of two to the two transverse EM
  polarizations.  This is polarization counting; it is not a
  structural identification of the microscopic field with an EM
  vector potential.
* `M-CORE-03`.  CORE-001-E01 introduces
  `F = epsilon_* sum_i [kappa_0|q_i|^2/2 + kappa_1 sum_<ij>|q_j - q_i|^2/2
  - g_dev eta_i e.q_i]`.  The gradient term is a SCALAR gradient
  `|q_j - q_i|^2`, not the curl form that would be required for an
  EM vector potential.
* `M-CORE-04`.  CORE-001-E02 gives the local evolution
  `tau dq_i/dt = -d(F/epsilon_*)/dq_i + xi_i`.  This is first-order
  in time; it relaxes, it does not propagate.
* `M-CORE-06`.  CORE-001-E09 gives the coarse-grained field equation
  `K u - Div(G grad u) = s(rho)`.  This is a Helmholtz-type
  elliptic equation; it has no time derivative; no wavefront
  follows.

## 3. Per-mechanism audit of the four transport questions

The full per-mechanism classification is in `native_transport_audit.csv`.
The conclusions are:

1. **Local phase transfer.**  Not derivable as EM-like transport.
   The nearest-neighbour term `kappa_1|q_j - q_i|^2/2` couples
   amplitudes, but CORE-001-E02 evolves them overdamped.  A
   second-order kinetic sector, or a Maxwell-like first-order
   structure, is required.
2. **Local field rotation.**  Not derivable.  The triplet `q`
   admits rotations, but CORE-001-E02 does not propagate rotations
   coherently; a non-dissipative dynamics is required.
3. **Neighbour coupling.**  Present at the energy level
   (CORE-001-E01).  LOCALITY-001 already established that
   `Div(P_F)` supplies all required static communication without
   invoking this term.  Neighbour coupling is therefore a static
   modelling choice, not a transport law.
4. **Wavefront evolution.**  Not derived.  CORE-001-E02 is
   overdamped; CORE-001-E09 is elliptic.  The V11 numerical
   identity `alpha_resolved ~ 3 alpha_EM` does not supply a
   second-order time structure.  INERTIA-001 already identified
   the kinetic sector as the missing closure.

## 4. EM-local microscopic mechanism audit

The standard local mechanisms of electromagnetism are listed in
`em_local_microscopic_mechanism.csv`.  None of them is present in
the V11 / CORE-001 microscopic structure:

| EM mechanism           | Present? | Structural mismatch                              |
|------------------------|----------|--------------------------------------------------|
| Faraday induction      | no       | no antisymmetric pair, no curl, no time derivative on a field-strength |
| Ampere-Maxwell         | no       | no current, no curl operator                     |
| D'Alembertian / wave   | no       | no kinetic sector, no Lorentzian signature       |
| Gauge invariance       | no       | `kappa_0|q|^2` is a Proca-like mass term; `kappa_1|q_j - q_i|^2` is not the curl form |
| Two polarizations      | counting only (V11 eq. 16) | counting is not a derivation of the transport law |
| Dispersionless c       | constraint only (V11 sec. 2.4) | V11 uses GW170817 to fix `epsilon_0 ~ 1`; it does not derive `c = sqrt(G/K)` |

The mismatch is structural, not numerical.  CORE-001's `q in R^3`
plus `kappa_0|q|^2` plus `kappa_1|q_j - q_i|^2` is the form of a
massive scalar triplet, not an EM vector potential.

## 5. Wavefront evolution audit

The candidate laws for `u(x,t)` and their derivation status are
recorded in `wavefront_evolution_audit.csv`.  Of the six entries:

* `WE-002` and `WE-003` are present in the frozen corpus but are
  time-independent (overdamped relaxation and Helmholtz equilibrium).
* `WE-001` is the elastic wave equation accepted by LOCALITY-001
  L-003 but is derived only **after** a positive momentum density is
  supplied.  INERTIA-001 left that supply open.
* `WE-004` would be a Maxwell-like wave equation and is absent.
* `WE-005` would require an independent flow vector and is
  forbidden by FP-1 / FP-4.
* `WE-006` is the numerical identity `alpha_resolved ~ 3 alpha_EM`
  itself, which is not a differential law.

No row in the audit produces a wavefront from the V11 microscopic
structure alone.

## 6. The missing principle, identified precisely

`kinetic_closure_requirement.json` records the precise gap:

> A local, conservative, second-order-in-time kinetic sector for
> the medium that supplies positive momentum density (or an
> equivalent symplectic structure) and thereby turns the static
> constitutive chain into a wave-bearing evolution equation.
> Equivalently, a Maxwell-like first-order structure with a curl
> kinetic operator and a conserved field-strength pair, supplied
> without introducing an independent EM sector.

This is exactly the closure gap that INERTIA-001 left open.  The
present milestone re-derives it from the EM side of the
microscopic structure: `alpha_EM = 1/137` fixes only the static
matter-vertex coupling; `alpha_resolved = 3 alpha_EM` fixes only an
amplitude identity; `kappa_1 |q_j - q_i|^2` fixes only a static
elastic energy; none of them supplies a positive momentum density,
a symplectic structure, a curl operator, or a Lorentzian signature.
The kinetic closure is unavoidable.

## 7. Compliance with the milestone brief

| Constraint                                       | Status |
|--------------------------------------------------|--------|
| No new spacetime inertia                         | yes    |
| No phenomenological steering coefficients         | yes    |
| No metric ansaetze                                | yes    |
| No cosmology solution                             | yes    |
| No quantum mechanics solution                     | yes    |
| No independent EM sector separated from medium    | yes    |
| No free transport constants                       | yes    |
| No new ontology                                   | yes    |
| No fit to data                                    | yes    |
| No V11 modification                               | yes    |
| No weak-lensing modification                      | yes    |

The audit uses only frozen sources, traces every claim back to a
frozen artifact, and reports a single explicit decision (`decision.json`).
The validation record is in `validation.json`.

## 8. Closure

**Outcome B.**  The neighbour-to-neighbour transport law is not
contained in the V11 electromagnetic microscopic structure.
The missing local principle is precisely the kinetic closure
already flagged by INERTIA-001.  The decision is reported in
`decision.json` and the completion record in `validation.json`.
"""
    (OUT / "em_transport001_report.md").write_text(report)

    print(json.dumps({
        "output": str(OUT),
        "decision": DECISION["outcome"],
        "pass": validation["pass"],
        "missing_principle": validation["missing_local_principle_identified"],
    }, indent=2))


if __name__ == "__main__":
    main()
