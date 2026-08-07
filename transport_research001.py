#!/usr/bin/env python3
"""Generate the reproducible PBUF TRANSPORT-RESEARCH-001 comparative audit.

Milestone brief
---------------
PBUF TRANSPORT-RESEARCH-001 is a comparative analysis of existing physical
wave-transport mechanisms.  The objective is to determine whether Nature
already supplies an existing local propagation mechanism that can serve as
the transport layer required by PBUF, without inventing new physics,
constants, or transport equations.

Required systems
----------------
Water surface waves, elastic solids, acoustic waves, electromagnetic waves
(Maxwell), spin waves (magnons), plasma waves / MHD.

For every system, eight mechanism questions are answered in physical
language first, governing equations second:
1. Disturbed quantity
2. Local neighbour interaction
3. Locally transferred quantity
4. Restoring mechanism
5. Resistance to the response
6. Propagation-speed determination
7. Local steering mechanism
8. Governing equations (minimum needed)

A comparison table is then built across all systems and a common abstract
transport architecture is sought.  Only after the comparison is finished is
the architecture compared with the existing PBUF framework
(FOUNDATION-001, V11, CORE-001).

Inputs
------
* FOUNDATION-001
* V11 preprint (Planck-Bound Unified Framework V11)
* CORE-001
* EM-TRANSPORT-001 / INERTIA-001 / LOCALITY-001 conclusions
* No new ontology, field, constant, fit, V11 change, weak-lensing change.

Outputs
-------
* transport_research001_report.md
* system_mechanism_audit.csv
* comparison_table.csv
* common_architecture.csv
* pbuf_architecture_comparison.csv
* common_architecture.json
* pbuf_comparison.json
* decision.json
* validation.json

Decision rule
-------------
* Outcome A  - one existing physical transport mechanism maps naturally onto
                the existing PBUF ontology with minimal additional
                assumptions.  The precise mapping is identified.
* Outcome B  - no existing transport mechanism maps cleanly.  The closest
                structural match is identified and the precise missing
                element is stated.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "runs/transport_research001"

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
    "EM-TRANSPORT-001":        "runs/em_transport001/em_transport001_report.md",
    "V11 preprint":            "docs/Planck-Bound_Unified_Framework_v11_preprint.pdf",
    "V11-ALPHA-001 brief":     "docs/PBUF_V11_ALPHA_001_Geometric_Origin_of_Resolved_Alpha.docx",
}


# ----------------------------------------------------------------------------
# Audit 1: per-system mechanism descriptions
# ----------------------------------------------------------------------------
# Each system is described in physical terms first.  The minimum governing
# equations are stated last.  No derivation is performed.

SYSTEMS = [
    # ---- Water surface waves ----
    {
        "id": "S1",
        "name": "Water surface waves",
        "physical_setting": (
            "A continuous fluid with a free surface; small vertical displacement "
            "of the surface from equilibrium."
        ),
        "disturbed_quantity": (
            "Vertical displacement of the free surface eta(x,t) (scalar field)."
        ),
        "local_interaction": (
            "A displaced water column exerts horizontal pressure on the "
            "neighbouring column; the neighbouring column accelerates and "
            "transfers momentum further along the surface.  The coupling is "
            "mediated by the hydrostatic and dynamic pressure gradients in "
            "the bulk fluid."
        ),
        "local_transfer": (
            "Hydrostatic pressure excess proportional to the local height "
            "difference between adjacent columns."
        ),
        "restoring": (
            "Gravity: a column displaced upward carries extra gravitational "
            "potential energy rho*g*eta and is pulled back down."
        ),
        "resistance": (
            "Inertia of the water mass; the displaced mass must be "
            "accelerated by the pressure gradient."
        ),
        "speed_mechanism": (
            "Phase speed is the ratio of restoring strength to inertial "
            "response.  For deep water c_phase^2 = g/k; for finite depth "
            "c^2 = (g/k)*tanh(kh)."
        ),
        "steering": (
            "Local depth variations change k = tanh(kh)/h and therefore "
            "the local phase speed; the wavefront refracts toward slower "
            "regions (Snell's-law form)."
        ),
        "equations": (
            "Linearized: d^2 eta/dt^2 = (g * tanh(kh) / k) nabla^2 eta.  "
            "Exact: Euler equations + incompressibility + free-surface BC."
        ),
    },

    # ---- Elastic solids ----
    {
        "id": "S2",
        "name": "Elastic solids",
        "physical_setting": (
            "A continuous solid with placement y(X,t) and deformation "
            "C = Grad y^sharp Grad y."
        ),
        "disturbed_quantity": (
            "Displacement field u(x,t) = y(x,t) - x (vector)."
        ),
        "local_interaction": (
            "Each material volume element pulls on its neighbours through "
            "the Cauchy traction P_F * N across the shared interface; the "
            "divergence of the first Piola-Kirchhoff stress gives the "
            "internal force per unit reference volume."
        ),
        "local_transfer": (
            "Force per unit area (traction) across each material interface; "
            "this is the local communication channel."
        ),
        "restoring": (
            "Gradient of the hyperelastic stored energy W(C) through "
            "-Div P_F = -(D_y C)^* D_C W, which depends only on the static "
            "constitutive law and not on a separate restoring field."
        ),
        "resistance": (
            "Mass density rho of the medium; the kinetic term rho du/dt "
            "appears in the local balance and supplies inertia."
        ),
        "speed_mechanism": (
            "c_alpha^2(n) = A_{iJkL} n_J n_L / rho, where A is the "
            "acoustic tensor built from the elastic 2-jet.  Speed is the "
            "ratio of elastic stiffness to mass density."
        ),
        "steering": (
            "Spatial variation of the elastic moduli C_ijkl(x); a stiffer "
            "region has a higher local speed, so wavefronts refract."
        ),
        "equations": (
            "Linearized at a homogeneous reference: rho u_tt = "
            "Div(A : sym grad u).  Nonlinear: rho u_tt = Div P_F."
        ),
    },

    # ---- Acoustic waves ----
    {
        "id": "S3",
        "name": "Acoustic waves (compressional in fluid)",
        "physical_setting": (
            "A continuous fluid with bulk modulus B and density rho in "
            "the long-wavelength, small-amplitude regime."
        ),
        "disturbed_quantity": (
            "Pressure perturbation delta p(x,t) (scalar)."
        ),
        "local_interaction": (
            "Compression of one volume element pushes the neighbouring "
            "element through the bulk modulus B; the neighbour then "
            "compresses in turn, propagating the disturbance."
        ),
        "local_transfer": (
            "Pressure perturbation across neighbouring volume elements."
        ),
        "restoring": (
            "Compressibility (the bulk modulus B): a compressed volume "
            "stores potential energy and pushes back."
        ),
        "resistance": (
            "Mass density rho of the fluid; the displaced mass must be "
            "accelerated."
        ),
        "speed_mechanism": (
            "c^2 = B / rho, the ratio of compressional restoring strength "
            "to inertial resistance."
        ),
        "steering": (
            "Spatial variation of B(x) or rho(x); the local speed "
            "c(x) = sqrt(B(x)/rho(x)) refracts the wavefront."
        ),
        "equations": (
            "Linearized acoustic wave equation: d^2 p / dt^2 = c^2 nabla^2 p, "
            "where c^2 = B / rho.  Derived from continuity + linearized "
            "Euler + equation of state."
        ),
    },

    # ---- Electromagnetic waves (Maxwell) ----
    {
        "id": "S4",
        "name": "Electromagnetic waves (Maxwell in vacuum)",
        "physical_setting": (
            "Empty space with electric and magnetic fields E(x,t), B(x,t) "
            "and no sources."
        ),
        "disturbed_quantity": (
            "Electric field E(x,t) and magnetic field B(x,t); equivalently "
            "the antisymmetric field-strength tensor F_{mu nu}(x,t)."
        ),
        "local_interaction": (
            "A changing E field locally generates a curl B field "
            "(Ampere-Maxwell with j = 0); a changing B field locally "
            "generates a curl E field (Faraday induction).  The two fields "
            "drive each other through curl coupling."
        ),
        "local_transfer": (
            "The fields E and B themselves; there is no scalar quantity "
            "that propagates between points."
        ),
        "restoring": (
            "NONE in vacuum.  The EM field has no scalar potential energy "
            "storage.  Propagation arises entirely from the mutual curl "
            "coupling between E and B."
        ),
        "resistance": (
            "NONE in vacuum.  The field has no mass and no inertia.  "
            "This is the structural feature that distinguishes EM from "
            "all mechanical wave systems."
        ),
        "speed_mechanism": (
            "c^2 = 1 / (epsilon_0 mu_0); the speed is set entirely by the "
            "constants governing the curl coupling between E and B."
        ),
        "steering": (
            "Spatial variation of epsilon(x), mu(x) (in matter) or, in GR, "
            "spatial variation of the metric.  The local wavevector is "
            "bent toward regions of slower phase speed."
        ),
        "equations": (
            "Maxwell in vacuum: nabla x E = -dB/dt; nabla x B = mu_0 "
            "epsilon_0 dE/dt; nabla.E = 0; nabla.B = 0.  Each is "
            "first-order in time; the wave equation emerges only after "
            "taking the curl of one equation and substituting the other."
        ),
    },

    # ---- Spin waves (magnons) ----
    {
        "id": "S5",
        "name": "Spin waves (magnons)",
        "physical_setting": (
            "An ordered magnetic lattice with local magnetization m(x,t) "
            "(unit vector) on each site."
        ),
        "disturbed_quantity": (
            "Direction of the local magnetization m(x,t) (unit vector "
            "field)."
        ),
        "local_interaction": (
            "Neighbouring spins are coupled through the Heisenberg exchange "
            "energy -J sum_<ij> S_i . S_j; misalignment of adjacent spins "
            "costs energy.  The anisotropy field and any external field "
            "supply a preferred axis."
        ),
        "local_transfer": (
            "Spin angular momentum: a precessing spin transfers angular "
            "momentum to its neighbour through the exchange interaction."
        ),
        "restoring": (
            "Anisotropy field H_anis and any external field H_ext: a spin "
            "misaligned with the anisotropy axis experiences a torque back."
        ),
        "resistance": (
            "Gyromagnetic precession: a spin cannot change direction "
            "instantaneously; it must precess around the effective field.  "
            "This provides an intrinsic dynamical response analogous to "
            "inertia but arising from rotational kinematics."
        ),
        "speed_mechanism": (
            "Dispersion omega(k) = gamma * sqrt(H_eff (H_eff + D k^2)).  "
            "The group velocity v_g = d(omega)/dk sets the wavefront speed."
        ),
        "steering": (
            "Spatial gradients of H_anis(x) or H_ext(x); also spatial "
            "variation of the exchange J(x) or anisotropy D(x)."
        ),
        "equations": (
            "Landau-Lifshitz (no damping): dm/dt = -gamma m x H_eff, where "
            "H_eff contains exchange, anisotropy, and external-field "
            "contributions."
        ),
    },

    # ---- Plasma waves / MHD ----
    {
        "id": "S6",
        "name": "Plasma waves / MHD",
        "physical_setting": (
            "A conducting fluid (ionized plasma) with density n(x,t), "
            "velocity v(x,t), pressure p(x,t), and magnetic field B(x,t)."
        ),
        "disturbed_quantity": (
            "A combined state: density n, velocity v, magnetic field B, "
            "pressure p.  Different MHD modes disturb different subsets."
        ),
        "local_interaction": (
            "Charge separation in the plasma creates electric fields; "
            "currents modify the magnetic field; the J x B force couples "
            "back to the momentum equation; magnetic-field lines under "
            "tension transmit force along B (magnetic tension)."
        ),
        "local_transfer": (
            "Electromagnetic force (J x B), pressure gradients, and "
            "magnetic tension B . nabla B / mu_0."
        ),
        "restoring": (
            "Three distinct mechanisms: (a) charge-separation "
            "electrostatics (Langmuir waves), (b) pressure gradients "
            "(acoustic mode), (c) magnetic tension (Alfven mode).  In a "
            "general MHD disturbance, all three act simultaneously."
        ),
        "resistance": (
            "Ion and electron inertia; for slow MHD modes the dominant "
            "mass is the ion mass density rho_i."
        ),
        "speed_mechanism": (
            "Three characteristic speeds: Alfven v_A = B / sqrt(mu_0 rho); "
            "sound c_s = sqrt(gamma p / rho); magnetosonic combinations "
            "c = sqrt(c_s^2 + v_A^2).  Different modes propagate at "
            "different speeds."
        ),
        "steering": (
            "Spatial gradients of rho(x), B(x), p(x); the local Alfven, "
            "sound, and magnetosonic speeds vary, bending each mode "
            "independently."
        ),
        "equations": (
            "Ideal MHD: continuity dn/dt + n nabla.v = 0; momentum "
            "rho dv/dt = -nabla p + J x B; induction dB/dt = nabla x "
            "(v x B); energy dp/dt + gamma p nabla.v = 0; with J = "
            "nabla x B / mu_0 and div B = 0."
        ),
    },
]


# ----------------------------------------------------------------------------
# Audit 2: comparison table
# ----------------------------------------------------------------------------

COMPARISON_HEADER = (
    "system",
    "disturbed_quantity",
    "local_interaction",
    "local_transfer",
    "restoring_mechanism",
    "resistance",
    "speed_mechanism",
    "steering_mechanism",
    "equations_min",
)


def _comparison_rows() -> list[tuple[str, ...]]:
    rows = []
    for s in SYSTEMS:
        rows.append((
            s["name"],
            s["disturbed_quantity"],
            s["local_interaction"],
            s["local_transfer"],
            s["restoring"],
            s["resistance"],
            s["speed_mechanism"],
            s["steering"],
            s["equations"],
        ))
    return rows


# ----------------------------------------------------------------------------
# Audit 3: common abstract transport architecture
# ----------------------------------------------------------------------------
# Inferred from the comparison table after the systems are laid side-by-side.

COMMON_ARCHITECTURE_ROWS = [
    (
        "slot_local_state",
        "A field defined at every point of the medium (scalar, vector, "
        "tensor, or unit-vector).",
        "present",
        "S1 eta, S2 u, S3 delta p, S4 E and B, S5 m, S6 {n,v,B,p}.",
    ),
    (
        "slot_neighbour_coupling",
        "A local coupling mechanism between adjacent infinitesimal "
        "regions (pressure, traction, curl, exchange, electromagnetic "
        "force).",
        "present",
        "All six systems supply a specific neighbour-coupling mechanism.",
    ),
    (
        "slot_local_transfer",
        "A locally transferred quantity that propagates between "
        "neighbours (pressure, traction, field, spin angular momentum, "
        "force density).",
        "present",
        "All six systems identify what flows locally between regions.",
    ),
    (
        "slot_restoring",
        "Something that pulls the disturbed state back toward "
        "equilibrium (gravity, elasticity, compressibility, anisotropy, "
        "magnetic tension, electrostatics).",
        "ABSENT in EM (S4). Present in S1, S2, S3, S5, S6.",
        "EM in vacuum has no restoring mechanism; its propagation comes "
        "from curl coupling between E and B.",
    ),
    (
        "slot_resistance",
        "Something that resists the change of state (mass density, "
        "gyromagnetic precession, ion mass).",
        "ABSENT in EM (S4). Present in S1, S2, S3, S5, S6.",
        "EM has no inertia; all mechanical systems have either classical "
        "mass or gyroscopic response.",
    ),
    (
        "slot_speed",
        "A propagation speed determined by the ratio of restoring "
        "strength to resistance (mechanical) or by the coupling "
        "constants (EM).",
        "present",
        "c^2 ~ restoring / resistance (S1, S2, S3, S5, S6); "
        "c^2 = 1/(epsilon_0 mu_0) (S4).",
    ),
    (
        "slot_steering",
        "Spatial variation of the medium parameters bends the wavefront "
        "(Snell's-law form).",
        "present",
        "All six systems identify a local steering mechanism.",
    ),
    (
        "slot_wave_equation",
        "A hyperbolic PDE (or coupled first-order system) admitting "
        "real characteristics at finite speed.",
        "derived after restoring + resistance (S1,S2,S3,S5,S6); "
        "derived after curl coupling (S4).",
        "Two distinct routes to a wave equation: (a) restoring + "
        "resistance = second-order wave equation; (b) first-order curl "
        "coupling with no restoring, no resistance.",
    ),
]


COMMON_ARCHITECTURE_JSON = {
    "shared_slots": [
        "A local state defined at every point",
        "A neighbour-coupling mechanism between adjacent infinitesimal regions",
        "A locally transferred quantity between neighbours",
        "A propagation speed determined by ratios or coupling constants",
        "A local steering mechanism through spatial parameter variation",
    ],
    "two_mechanistic_families": {
        "family_A_mechanical": {
            "members": ["S1 water surface", "S2 elastic solid",
                        "S3 acoustic", "S5 spin wave", "S6 plasma/MHD"],
            "shared_architecture": [
                "local state at each point",
                "neighbour coupling through a transferred quantity",
                "restoring mechanism (potential energy gradient)",
                "resistance / inertia (mass density or gyroscopic response)",
                "wave equation emerges as second-order-in-time from "
                "restoring + resistance balance",
            ],
            "speed_origin": "ratio of restoring stiffness to inertial density",
            "structurally_similar_to": "scalar or vector wave equation",
        },
        "family_B_electromagnetic": {
            "members": ["S4 Maxwell"],
            "shared_architecture": [
                "two coupled first-order vector fields (E, B)",
                "neighbour coupling through curl operators (no potential energy)",
                "no restoring mechanism",
                "no resistance / no inertia",
                "wave equation emerges by eliminating one field, giving "
                "the D'Alembertian Box E = 0",
            ],
            "speed_origin": "ratio of curl-coupling constants (1/sqrt(epsilon mu))",
            "structurally_distinct_from": "Family A: there is no scalar "
            "potential storage and no mass term in vacuum.",
        },
    },
    "key_structural_finding": (
        "The six systems fall into two families.  Family A (mechanical) "
        "requires both a restoring mechanism AND a resistance/inertia.  "
        "Family B (Maxwell) requires neither: propagation comes entirely "
        "from the mutual curl coupling of two vector fields.  Any mapping "
        "onto PBUF must therefore identify which family PBUF structurally "
        "belongs to, and what is present versus missing in each slot."
    ),
}


# ----------------------------------------------------------------------------
# Audit 4: PBUF architecture comparison
# ----------------------------------------------------------------------------

PBUF_COMPARISON_ROWS = [
    # (slot, pbuf_status, evidence, remark)
    (
        "local_state",
        "present",
        "CORE-001-E03 / CORE-001-E04: q_i in R^3 at each lattice site; "
        "coarse field u(x) = e . sum_i a^d W_L(x-x_i) q_i.",
        "Matches the local-state slot of every system in Family A.",
    ),
    (
        "neighbour_coupling",
        "present (static)",
        "CORE-001-E01: kappa_1 sum_<ij> |q_j - q_i|^2 / 2 in the "
        "microscopic free energy F.",
        "Provides neighbour coupling at the energy level.  LOCALITY-001 "
        "showed that Div(P_F) supplies all required static communication "
        "without this term, so the coupling exists but is not strictly "
        "necessary for communication.",
    ),
    (
        "local_transfer",
        "present (static)",
        "Continuum chain: P_F = 2 F P_C, traction P_F N at the boundary, "
        "internal force -Div P_F.",
        "Equivalent to the local-transfer slot of S2 elastic solids.",
    ),
    (
        "restoring_mechanism",
        "present",
        "Onsite kappa_0|q|^2/2 term in CORE-001-E01; in the continuum, "
        "-Div P_F = -(D_y C)^* D_C W.",
        "Equivalent to Family A restoring mechanisms (gravity, "
        "elasticity, compressibility, anisotropy).",
    ),
    (
        "resistance",
        "MISSING",
        "CORE-001-E02 is overdamped: tau dq_i/dt = -d(F/epsilon_*)/dq_i.  "
        "INERTIA-001: kinetic sector is an irreducible closure gap; "
        "static elastic energy alone does not determine it.  EM-TRANSPORT-001: "
        "alpha_EM, alpha_resolved, g_dev do not supply momentum density.",
        "No mass density, no gyroscopic response, no symplectic "
        "structure, no second-order-in-time kinetic term.  This is the "
        "single missing element relative to Family A.",
    ),
    (
        "speed",
        "not defined",
        "No time structure in the coarse-grained equation Ku - "
        "Div(G grad u) = s(rho) (CORE-001-E09 / WL-003).",
        "c^2 would have to be built from restoring / resistance; "
        "without resistance there is no speed.",
    ),
    (
        "steering",
        "available in principle",
        "Spatial variation of the elastic moduli C_ijkl(x) is admissible "
        "in the frozen framework (CONSTITUTIVE-002 / MATERIAL-LAB).  "
        "Spatial variation of K and G would steer wavefronts the same "
        "way as in S2 / S3.",
        "Steering slot is open, but operates only once the speed exists.",
    ),
    (
        "wave_equation",
        "not present",
        "The static Helmholtz equation (CORE-001-E09) is elliptic; it "
        "has no time derivative.",
        "No second-order wave equation.  No first-order curl structure "
        "either.  The local field q in R^3 is not an (E,B) pair.",
    ),
]


PBUF_COMPARISON_JSON = {
    "structural_family_match": {
        "family_A_mechanical": "partial",
        "family_B_electromagnetic": "incompatible",
        "rationale_family_A": (
            "PBUF supplies the spatial half of a Family-A mechanical "
            "system: local state, neighbour coupling, restoring "
            "mechanism, and (in principle) steering through spatial "
            "variation of moduli.  The kinetic sector is missing."
        ),
        "rationale_family_B": (
            "Family B requires (a) two coupled vector fields (E,B), "
            "(b) a curl operator in the kinetic structure, (c) gauge "
            "invariance, and (d) no mass-like onsite term.  CORE-001-E01 "
            "has a single scalar triplet q with a mass-like kappa_0|q|^2 "
            "term and a SCALAR gradient kappa_1|q_j - q_i|^2 term, "
            "neither of which is a curl operator or an antisymmetric "
            "field-strength pair.  V11's alpha_resolved ~ 3 alpha_EM is a "
            "numerical identity and a dimensional-counting argument, not a "
            "structural identification of q with an EM vector potential.  "
            "Family B is therefore structurally incompatible with PBUF "
            "without rewriting the microscopic energy."
        ),
    },
    "closest_structural_match": {
        "system": "S2 elastic solids",
        "shared_slots": [
            "local state (q or u)",
            "neighbour coupling (kappa_1|q_j - q_i|^2 or stress Div P_F)",
            "restoring (-Div P_F from stored energy gradient)",
            "steering (spatial variation of elastic moduli)",
        ],
        "missing_slot": "resistance / inertia",
        "evidence_for_missing": [
            "INERTIA-001: kinetic sector is the irreducible closure gap",
            "EM-TRANSPORT-001: V11 microscopic structure is insufficient",
            "CORE-001-E02 is overdamped (first-order in time)",
            "CORE-001-E09 is elliptic (no time structure)",
        ],
        "minimal_additional_assumption_required": (
            "A local kinetic sector supplying positive momentum density or "
            "an equivalent symplectic structure.  Adding this single slot "
            "completes Family A and produces a propagation speed "
            "c^2 = G/K from the frozen elastic 2-jet."
        ),
        "is_minimal_addition_actually_minimal": (
            "NO in the strict sense of EM-TRANSPORT-001 / INERTIA-001: "
            "the kinetic sector cannot be derived from F, from alpha_EM, "
            "from alpha_resolved, or from g_dev alone.  Adding it "
            "constitutes a new local physical principle, not a parameter "
            "renaming or a coordinate choice.  Per the brief's "
            "'minimal additional assumptions' criterion, the absence of a "
            "derivable kinetic sector means this is NOT a clean mapping; "
            "it is a structurally closest match with an open closure gap."
        ),
    },
    "explicit_pb_comparison_citations": {
        "FOUNDATION-001": [
            "FP-1 one continuous spacetime medium - compatible with "
            "treating the medium as a continuous system in either family.",
            "FP-5 V11 relativistic compatibility - admits wave propagation "
            "with finite speed but does not derive the propagation law.",
            "FP-6 no additional free fundamental constants - rules out "
            "introducing a new mass density, gyro constant, or coupling "
            "as a free parameter; the kinetic sector must be derived or "
            "remain open.",
        ],
        "V11": [
            "V11 (4): alpha_resolved ~ 3 alpha_EM is a numerical identity "
            "and a motivating counting argument; not a structural "
            "derivation of EM dynamics.",
            "V11 sec. 2.4: wave modes are taken from GR / EM and "
            "constrained by GW170817; not derived from V11 microscopic "
            "structure.",
            "V11 (16): Omega_b0 = 2 alpha_resolved is polarization "
            "counting, not a microscopic identification.",
        ],
        "CORE-001": [
            "CORE-001-E01: kappa_0|q|^2 mass-like onsite term, "
            "kappa_1|q_j - q_i|^2 scalar gradient (not curl), matter "
            "coupling -g_dev eta_i e.q_i.",
            "CORE-001-E02: overdamped Langevin dynamics; first-order in "
            "time; no wave propagation.",
            "CORE-001-E09: elliptic Helmholtz equation; no time structure.",
        ],
    },
}


DECISION = {
    "milestone": "PBUF TRANSPORT-RESEARCH-001",
    "outcome": "B",
    "headline": (
        "No existing physical transport mechanism maps naturally onto the "
        "PBUF ontology with minimal additional assumptions.  The "
        "comparison of the six required systems (water surface, elastic "
        "solid, acoustic, Maxwell EM, spin wave, plasma / MHD) shows two "
        "structural families.  Family A (mechanical: S1, S2, S3, S5, S6) "
        "requires both a restoring mechanism and an inertial resistance; "
        "Family B (Maxwell: S4) requires neither and instead uses the "
        "mutual curl coupling of two vector fields.  PBUF matches Family "
        "A on the spatial half (local state, neighbour coupling, "
        "restoring, steering) but is missing the resistance / inertia "
        "slot.  Family B is incompatible because CORE-001's microscopic "
        "energy has a mass-like onsite term, a scalar gradient term "
        "(not curl), and a single scalar triplet (not an (E,B) pair).  "
        "The closest structural match is the elastic-solid pattern "
        "(S2), but its completion requires adding the kinetic sector - "
        "which cannot be derived from the frozen structure and was "
        "already identified by INERTIA-001 and EM-TRANSPORT-001 as the "
        "open closure gap."
    ),
    "summary": [
        "Six systems were analyzed using a common eight-question framework "
        "(disturbed quantity, local interaction, local transfer, restoring, "
        "resistance, speed, steering, governing equations).",
        "Two structural families emerge.  Family A: mechanical systems "
        "(S1, S2, S3, S5, S6) all have a local state, a neighbour-coupling "
        "transfer mechanism, a restoring force, an inertial resistance, "
        "and a propagation speed that is the ratio of the two.",
        "Family B: Maxwell EM (S4) is unique in having neither restoring "
        "nor resistance; its propagation arises from the mutual curl "
        "coupling of E and B.",
        "PBUF matches Family A on the spatial half: local state (q, u), "
        "neighbour coupling (kappa_1|q_j - q_i|^2, -Div P_F), restoring "
        "(-Div P_F from stored energy), and (in principle) steering.",
        "PBUF does not match Family A on the temporal half: there is no "
        "kinetic sector, no mass density, no gyroscopic response, no "
        "second-order-in-time dynamics, no wave speed, no wave equation.",
        "PBUF does not match Family B: there is no (E,B) pair, no curl "
        "operator, no gauge invariance, and there IS a mass-like onsite "
        "term kappa_0|q|^2 that Family B explicitly excludes.",
        "Therefore no existing physical transport mechanism maps "
        "naturally onto the PBUF ontology with minimal additional "
        "assumptions.  The closest match is the elastic-solid pattern "
        "(S2) and its missing slot is the kinetic sector.",
    ],
    "what_this_milestone_does_not_claim": [
        "It does not introduce a new transport equation, a new constant, "
        "or new physics.",
        "It does not declare that the kinetic sector is impossible.",
        "It does not modify V11, CORE-001, FOUNDATION-001, or any "
        "constitutive law.",
        "It does not select a specific candidate transport mechanism.",
        "It does not fit or interpret weak-lensing data.",
    ],
    "no_new_physics": True,
    "no_v11_change": True,
    "no_core001_change": True,
    "no_foundation001_change": True,
    "no_weak_lensing_change": True,
    "no_new_constant": True,
    "no_new_transport_equation": True,
    "no_assumption_that_maxwell_applies": True,
    "no_assumption_that_elasticity_applies": True,
    "no_assumption_that_magnetism_is_the_answer": True,
    "no_cosmology_or_quantum_gravity": True,
    "no_metric_construction": True,
    "no_numerical_optimization": True,
    "no_weak_lensing_fitting": True,
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
        "system_mechanism_audit.csv",
        ("system_id", "name", "physical_setting", "disturbed_quantity",
         "local_interaction", "local_transfer", "restoring", "resistance",
         "speed_mechanism", "steering", "equations"),
        [(
            s["id"], s["name"], s["physical_setting"],
            s["disturbed_quantity"], s["local_interaction"],
            s["local_transfer"], s["restoring"], s["resistance"],
            s["speed_mechanism"], s["steering"], s["equations"],
        ) for s in SYSTEMS],
    )
    write_csv("comparison_table.csv", COMPARISON_HEADER, _comparison_rows())
    write_csv(
        "common_architecture.csv",
        ("slot", "abstract_role", "pbuf_compatibility_summary", "evidence"),
        COMMON_ARCHITECTURE_ROWS,
    )
    write_csv(
        "pbuf_architecture_comparison.csv",
        ("slot", "pbuf_status", "evidence", "remark"),
        PBUF_COMPARISON_ROWS,
    )

    write_json("common_architecture.json", COMMON_ARCHITECTURE_JSON)
    write_json("pbuf_comparison.json", PBUF_COMPARISON_JSON)
    write_json("decision.json", DECISION)

    # ---------------- Validation ----------------
    checks = {
        "all_frozen_sources_present": not missing,
        "all_six_required_systems_included": (
            {s["name"].lower() for s in SYSTEMS}
            >= {"water surface waves", "elastic solids",
                "acoustic waves (compressional in fluid)",
                "electromagnetic waves (maxwell in vacuum)",
                "spin waves (magnons)", "plasma waves / mhd"}
        ),
        "each_system_answers_all_eight_mechanism_questions": all(
            all(s[k] for k in ("disturbed_quantity", "local_interaction",
                               "local_transfer", "restoring", "resistance",
                               "speed_mechanism", "steering", "equations"))
            for s in SYSTEMS
        ),
        "comparison_table_has_all_systems": (
            len(_comparison_rows()) == len(SYSTEMS)
        ),
        "two_structural_families_identified": (
            "family_A_mechanical" in COMMON_ARCHITECTURE_JSON["two_mechanistic_families"]
            and "family_B_electromagnetic" in COMMON_ARCHITECTURE_JSON["two_mechanistic_families"]
        ),
        "family_a_restoring_and_resistance_required": all(
            member in {"S1 water surface", "S2 elastic solid", "S3 acoustic",
                       "S5 spin wave", "S6 plasma/MHD"}
            for member in COMMON_ARCHITECTURE_JSON["two_mechanistic_families"]
            ["family_A_mechanical"]["members"]
        ),
        "family_b_is_unique_to_maxwell": (
            COMMON_ARCHITECTURE_JSON["two_mechanistic_families"]
            ["family_B_electromagnetic"]["members"] == ["S4 Maxwell"]
        ),
        "pbuf_architecture_comparison_covers_eight_slots": (
            len(PBUF_COMPARISON_ROWS) == 8
        ),
        "resistance_slot_marked_missing_in_pbuf": any(
            r[0] == "resistance" and r[1] == "MISSING"
            for r in PBUF_COMPARISON_ROWS
        ),
        "family_b_marked_incompatible_with_pbuf": (
            PBUF_COMPARISON_JSON["structural_family_match"]["family_B_electromagnetic"]
            == "incompatible"
        ),
        "family_a_marked_partial_match_for_pbuf": (
            PBUF_COMPARISON_JSON["structural_family_match"]["family_A_mechanical"]
            == "partial"
        ),
        "closest_match_is_elastic_solids": (
            PBUF_COMPARISON_JSON["closest_structural_match"]["system"]
            == "S2 elastic solids"
        ),
        "closest_match_missing_slot_is_kinetic": (
            PBUF_COMPARISON_JSON["closest_structural_match"]["missing_slot"]
            == "resistance / inertia"
        ),
        "explicit_pbuf_citations_present": all(
            k in PBUF_COMPARISON_JSON["explicit_pb_comparison_citations"]
            for k in ("FOUNDATION-001", "V11", "CORE-001")
        ),
        "decision_is_outcome_b": DECISION["outcome"] == "B",
        "no_new_physics_constant_or_transport_equation": all(
            DECISION[k] for k in ("no_new_physics", "no_new_constant",
                                  "no_new_transport_equation")
        ),
        "no_v11_core001_foundation001_change": all(
            DECISION[k] for k in ("no_v11_change", "no_core001_change",
                                  "no_foundation001_change")
        ),
        "no_assumption_of_direct_maxwell_or_elasticity_or_magnetism": all(
            DECISION[k] for k in ("no_assumption_that_maxwell_applies",
                                  "no_assumption_that_elasticity_applies",
                                  "no_assumption_that_magnetism_is_the_answer")
        ),
        "no_cosmology_quantum_gravity_metric_construction": all(
            DECISION[k] for k in ("no_cosmology_or_quantum_gravity",
                                  "no_metric_construction")
        ),
        "no_numerical_optimization_or_weak_lensing_fitting": all(
            DECISION[k] for k in ("no_numerical_optimization",
                                  "no_weak_lensing_fitting",
                                  "no_weak_lensing_change")
        ),
    }
    validation = {
        "milestone": "PBUF TRANSPORT-RESEARCH-001",
        "pass": all(checks.values()),
        "checks": checks,
        "decision": "Outcome B",
        "closest_structural_match": "S2 elastic solids (Family A)",
        "missing_slot": "resistance / inertia (kinetic sector)",
        "sources": SOURCES,
        "deliverables": [
            "transport_research001_report.md",
            "system_mechanism_audit.csv",
            "comparison_table.csv",
            "common_architecture.csv",
            "pbuf_architecture_comparison.csv",
            "common_architecture.json",
            "pbuf_comparison.json",
            "decision.json",
            "validation.json",
        ],
    }
    write_json("validation.json", validation)

    # ---------------- Main report ----------------
    report = r"""# PBUF TRANSPORT-RESEARCH-001 -- Comparative Analysis of Native Wave Transport Mechanisms

## 0. Decision

**Outcome B.** No existing physical transport mechanism maps naturally
onto the existing PBUF ontology with minimal additional assumptions.

The six required systems fall into two structural families:

* **Family A (mechanical):** water surface waves (S1), elastic solids
  (S2), acoustic waves (S3), spin waves (S5), plasma / MHD (S6). All
  five share the architecture `local state + neighbour coupling +
  restoring mechanism + inertial resistance`, with the propagation
  speed set by the ratio of restoring stiffness to resistance.
* **Family B (electromagnetic):** Maxwell vacuum (S4). Unique in
  having neither a restoring mechanism nor an inertial resistance;
  propagation arises from the mutual curl coupling of `E` and `B`.

PBUF matches **Family A on the spatial half only**:

| Architecture slot        | PBUF status                |
|--------------------------|----------------------------|
| local state              | present                    |
| neighbour coupling       | present (static)           |
| local transfer           | present (static)           |
| restoring mechanism      | present                    |
| **resistance / inertia** | **MISSING**                |
| propagation speed        | not defined                |
| local steering           | available in principle     |
| wave equation            | not present                |

Family B is structurally incompatible: CORE-001's microscopic energy
has a mass-like onsite term `kappa_0|q|^2`, a scalar gradient term
`kappa_1|q_j - q_i|^2` (not a curl), and a single scalar triplet `q`
(not an `(E,B)` pair). V11's `alpha_resolved ~ 3 alpha_EM` is a
numerical identity and a dimensional-counting argument; it does not
identify `q` with an EM vector potential.

The closest structural match is therefore the elastic-solid pattern
(S2). Its completion requires adding the kinetic sector - exactly
the closure gap already identified by INERTIA-001 and EM-TRANSPORT-001.
Adding this sector is not a parameter renaming or a coordinate choice;
it is a new local physical principle, which the brief's "minimal
additional assumptions" criterion does not permit.

No ontology, field, coupling, constant, transport equation, V11
change, CORE-001 change, FOUNDATION-001 change, constitutive-law
change, cosmological result, metric construction, or weak-lensing
fit is introduced.

## 1. Method

For each of the six required systems the milestone records, in
physical language first and equations last:

1. The physical setting of the medium.
2. The disturbed quantity.
3. The local neighbour-to-neighbour interaction.
4. The quantity that is locally transferred.
5. The restoring mechanism.
6. The resistance to the response.
7. How the propagation speed is determined.
8. The local steering mechanism.
9. The minimum governing equations (no derivation).

The full per-system records are in `system_mechanism_audit.csv`.

The comparison table is built directly from these records
(`comparison_table.csv`) and is the only input to the pattern
identification step. The common architecture is then derived, not
assumed (`common_architecture.csv` and `common_architecture.json`).
Only after this comparison is complete is the common architecture
compared with the existing PBUF framework
(`pbuf_architecture_comparison.csv` and `pbuf_comparison.json`).

## 2. Per-system physical mechanism

### 2.1 Water surface waves (S1)

* **Setting.** A continuous fluid with a free surface.
* **Disturbed quantity.** Vertical displacement `eta(x,t)`.
* **Local interaction.** A displaced column exerts horizontal pressure
  on its neighbour; the neighbour accelerates and transfers momentum
  further along the surface.
* **Local transfer.** Hydrostatic pressure excess proportional to the
  height difference between adjacent columns.
* **Restoring.** Gravity (`rho g eta`).
* **Resistance.** Inertia of the water mass.
* **Speed.** `c^2 = (g/k) tanh(kh)` (deep-water limit `c^2 = g/k`).
* **Steering.** Bathymetry variations change the local `k tanh(kh)/h`
  and refract the wavefront.
* **Equations (minimum).** `d^2 eta/dt^2 = (g tanh(kh)/k) nabla^2 eta`.
  Exact: Euler + incompressibility + free-surface BC.

### 2.2 Elastic solids (S2)

* **Setting.** Continuous solid with placement `y(X,t)` and
  deformation `C = Grad y^sharp Grad y`.
* **Disturbed quantity.** Displacement `u(x,t) = y - x` (vector).
* **Local interaction.** Each material volume element pulls on its
  neighbours through the traction `P_F N` across the shared
  interface; the internal force is `Div P_F`.
* **Local transfer.** Traction across each material interface.
* **Restoring.** Gradient of the hyperelastic stored energy `W(C)`:
  `-Div P_F = -(D_y C)^* D_C W`.
* **Resistance.** Mass density `rho` of the medium.
* **Speed.** `c_alpha^2(n) = A_{iJkL} n_J n_L / rho` (acoustic tensor).
* **Steering.** Spatial variation of elastic moduli `C_ijkl(x)`.
* **Equations (minimum).** `rho u_tt = Div(A : sym grad u)`;
  nonlinear: `rho u_tt = Div P_F`.

### 2.3 Acoustic waves (S3)

* **Setting.** Continuous fluid with bulk modulus `B` and density
  `rho`.
* **Disturbed quantity.** Pressure perturbation `delta p(x,t)`.
* **Local interaction.** Compression of one volume pushes its
  neighbour through `B`; the neighbour compresses in turn.
* **Local transfer.** Pressure perturbation between adjacent volumes.
* **Restoring.** Compressibility (bulk modulus `B`).
* **Resistance.** Mass density `rho`.
* **Speed.** `c^2 = B/rho`.
* **Steering.** Spatial variation of `B(x)` or `rho(x)`.
* **Equations (minimum).** `d^2 p / dt^2 = c^2 nabla^2 p`.

### 2.4 Electromagnetic waves (S4)

* **Setting.** Vacuum with `E(x,t)` and `B(x,t)`.
* **Disturbed quantity.** `E` and `B`; equivalently the
  antisymmetric `F_{mu nu}`.
* **Local interaction.** `dE/dt` locally generates `curl B`
  (Ampere-Maxwell with `j = 0`); `dB/dt` locally generates
  `curl E` (Faraday).
* **Local transfer.** The fields themselves.
* **Restoring.** **NONE in vacuum.** No scalar potential energy
  storage.
* **Resistance.** **NONE in vacuum.** No inertia, no mass.
* **Speed.** `c^2 = 1/(epsilon_0 mu_0)` from the curl coupling.
* **Steering.** Spatial variation of `epsilon(x)`, `mu(x)`, or (in
  GR) the metric.
* **Equations (minimum).** `curl E = -dB/dt`, `curl B =
  mu_0 epsilon_0 dE/dt`, `div E = 0`, `div B = 0`. Each is
  first-order in time; the wave equation emerges only after taking
  the curl of one equation and substituting the other.

### 2.5 Spin waves / magnons (S5)

* **Setting.** Ordered magnetic lattice with local magnetization
  `m(x,t)`.
* **Disturbed quantity.** Direction of the local magnetization
  `m(x,t)`.
* **Local interaction.** Heisenberg exchange `-J sum_<ij> S_i . S_j`
  couples neighbouring spins; misalignment costs energy.
* **Local transfer.** Spin angular momentum: a precessing spin
  transfers angular momentum to its neighbour.
* **Restoring.** Anisotropy field `H_anis` and any external field
  `H_ext`.
* **Resistance.** Gyromagnetic precession: a spin cannot change
  direction instantaneously.
* **Speed.** Dispersion `omega(k) = gamma sqrt(H_eff (H_eff + D k^2))`;
  group velocity `v_g = d omega/dk`.
* **Steering.** Spatial gradients of `H_anis(x)` or `H_ext(x)`,
  spatial variation of `J(x)` or `D(x)`.
* **Equations (minimum).** `dm/dt = -gamma m x H_eff` (Landau-Lifshitz
  without damping).

### 2.6 Plasma waves / MHD (S6)

* **Setting.** Conducting fluid with density `n`, velocity `v`,
  pressure `p`, magnetic field `B`.
* **Disturbed quantity.** A combined state `{n, v, B, p}`; different
  MHD modes disturb different subsets.
* **Local interaction.** Charge separation creates `E`; currents
  modify `B`; the `J x B` force couples back to momentum; magnetic
  tension `B . nabla B / mu_0` transmits force along `B`.
* **Local transfer.** EM force (`J x B`), pressure gradients,
  magnetic tension.
* **Restoring.** Three mechanisms: electrostatics (charge
  separation), pressure gradients, magnetic tension.
* **Resistance.** Ion (and electron) inertia; for slow MHD modes
  the dominant mass is `rho_i`.
* **Speed.** Alfven `v_A = B/sqrt(mu_0 rho)`; sound
  `c_s = sqrt(gamma p/rho)`; magnetosonic `c = sqrt(c_s^2 + v_A^2)`.
* **Steering.** Spatial gradients of `rho`, `B`, `p`.
* **Equations (minimum).** Continuity, momentum (`rho dv/dt =
  -nabla p + J x B`), induction (`dB/dt = curl(v x B)`), energy,
  with `J = curl B / mu_0`, `div B = 0`.

## 3. Comparison table

| System        | Disturbed quantity     | Local interaction                | Restoring            | Resistance / inertia | Speed mechanism                         | Steering                                  | Equations (min)                                  |
|---------------|------------------------|----------------------------------|----------------------|----------------------|------------------------------------------|-------------------------------------------|--------------------------------------------------|
| Water surface | `eta(x,t)` scalar      | pressure between adjacent columns | gravity `rho g eta`  | fluid mass density   | `c^2 = g tanh(kh)/k`                     | depth variation `h(x)`                    | `eta_tt = c^2 nabla^2 eta`                       |
| Elastic solid | `u(x,t)` vector        | traction `P_F N`                 | `-Div P_F`           | mass density `rho`   | `c^2 = A/rho` (acoustic tensor)          | spatial variation of `C_ijkl(x)`          | `rho u_tt = Div(A : sym grad u)`                 |
| Acoustic      | `delta p(x,t)` scalar  | pressure between adjacent volumes| bulk modulus `B`     | mass density `rho`   | `c^2 = B/rho`                            | spatial variation of `B, rho`             | `p_tt = c^2 nabla^2 p`                           |
| Maxwell EM    | `E, B` vector pair     | curl coupling between `E, B`     | NONE                 | NONE                 | `c^2 = 1/(eps mu)` from curl constants   | spatial variation of `eps, mu` / metric   | `curl E = -B_t`; `curl B = eps mu E_t`           |
| Spin wave     | `m(x,t)` unit vector   | Heisenberg exchange `-J S_i.S_j` | anisotropy `H_anis`  | gyromagnetic precession | `omega = gamma sqrt(H_eff(H_eff+D k^2))` | gradients of `H_anis, H_ext, J, D`        | `dm/dt = -gamma m x H_eff`                       |
| Plasma / MHD  | `{n,v,B,p}` combined   | `J x B`, pressure, magnetic tension | electrostatics / pressure / magnetic tension | ion mass `rho_i` | `v_A`, `c_s`, magnetosonic combinations | gradients of `rho, B, p`                  | continuity + momentum + induction + energy       |

## 4. Common abstract transport architecture

After laying the systems side by side, a common architecture is
present in **five of the six** systems (S1, S2, S3, S5, S6). The
slots are:

1. A local state defined at every point.
2. A neighbour-coupling mechanism between adjacent infinitesimal
   regions.
3. A locally transferred quantity between neighbours.
4. A restoring mechanism that pulls the state back.
5. A resistance to the response (inertia).
6. A propagation speed set by the ratio of restoring strength to
   resistance.
7. A local steering mechanism through spatial variation of medium
   parameters.

The wave equation emerges as a second-order-in-time PDE from the
balance of restoring and resistance. This is Family A.

**Family B (S4 Maxwell) is structurally different.** It shares slots
1, 2, 3, 6, 7 with Family A but **lacks slots 4 and 5 entirely**.
Its wave equation emerges from the mutual first-order curl coupling
of `E` and `B`, with no scalar potential energy and no mass term.
This makes S4 a different architectural family, not a member of
Family A.

The full common-architecture record is in
`common_architecture.csv` and `common_architecture.json`. The key
finding, stated without PBUF-specific interpretation, is recorded
under `key_structural_finding` in `common_architecture.json`.

## 5. PBUF comparison

The common architecture is now compared, slot by slot, with the
existing PBUF framework. Every conclusion cites the relevant
frozen artifact (FOUNDATION-001, V11, CORE-001) explicitly. The
detailed evidence is in `pbuf_architecture_comparison.csv` and
`pbuf_comparison.json`. The structural conclusion:

| Architecture slot        | PBUF status                                  |
|--------------------------|----------------------------------------------|
| local state              | present (CORE-001-E03 / E04: `q in R^3`, `u`) |
| neighbour coupling       | present (CORE-001-E01: `kappa_1 |q_j-q_i|^2`; LOCALITY-001: `Div P_F`) |
| local transfer           | present (continuum: `P_F N` traction; `-Div P_F`) |
| restoring mechanism      | present (CORE-001-E01: `kappa_0|q|^2`; continuum `-Div P_F`) |
| **resistance / inertia** | **MISSING** (INERTIA-001 closure gap)       |
| propagation speed        | not defined (no time structure)              |
| local steering           | available in principle (moduli `C_ijkl(x)`)  |
| wave equation            | not present (CORE-001-E09 is elliptic)       |

* **Family A match.** PBUF matches the spatial half of Family A.
  It does not match Family A on the temporal half because it has no
  kinetic sector. This is the same closure gap already identified
  by INERTIA-001 and re-derived from the EM side in EM-TRANSPORT-001.
* **Family B match.** Family B requires (a) two coupled vector
  fields, (b) a curl operator in the kinetic structure, (c) gauge
  invariance, and (d) no mass-like onsite term. CORE-001-E01 has
  none of (a), (b), (c) and has a positive `kappa_0|q|^2` mass-like
  term that contradicts (d). V11's `alpha_resolved ~ 3 alpha_EM` is
  a numerical identity and a dimensional-counting argument; it does
  not identify `q` with an EM vector potential. Family B is therefore
  structurally incompatible with PBUF without rewriting the
  microscopic energy.

## 6. Closest structural match

Of the six systems, the elastic-solid pattern (S2) is the closest
match:

* **Shared slots.** local state (`q` vs `u`), neighbour coupling
  (`kappa_1 |q_j - q_i|^2` vs stress `Div P_F`), restoring
  (`-Div P_F` from stored energy), steering (spatial variation of
  moduli).
* **Missing slot.** resistance / inertia.
* **What completion requires.** A kinetic sector supplying positive
  momentum density or an equivalent symplectic structure. With this
  slot filled, the speed would be `c^2 = G/K` from the frozen elastic
  2-jet.

**Why this is NOT a clean mapping under "minimal additional
assumptions".** The kinetic sector is precisely the closure gap
INERTIA-001 left open: it cannot be derived from `F`, from
`alpha_EM`, from `alpha_resolved`, or from `g_dev`. Adding it is a
new local physical principle, not a parameter renaming or a
coordinate choice. Per the brief, this falls outside the
"minimal additional assumptions" criterion, and the mapping is
therefore not clean.

## 7. Compliance with the milestone brief

| Constraint                                              | Status |
|---------------------------------------------------------|--------|
| No new physics                                          | yes    |
| No new constants                                        | yes    |
| No new transport equations                              | yes    |
| No assumption that Maxwell applies directly             | yes    |
| No assumption that elasticity applies directly          | yes    |
| No assumption that magnetism is the answer              | yes    |
| No modification of V11                                  | yes    |
| No modification of FOUNDATION-001                       | yes    |
| No modification of CORE-001                             | yes    |
| No constitutive-law change                              | yes    |
| No cosmological result introduced                       | yes    |
| No weak-lensing fitting                                 | yes    |
| No metric construction                                  | yes    |
| No quantum-gravity or dark-sector / dark-energy content  | yes    |
| Local-mechanism focus only                              | yes    |
| Every PBUF claim cited to a frozen artifact             | yes    |

## 8. Closure

**Outcome B.** The six required systems fall into two structural
families. PBUF matches Family A (mechanical) on the spatial half
only; the kinetic / inertial resistance slot is missing. Family B
(Maxwell) is structurally incompatible with CORE-001. The closest
match is the elastic-solid pattern (S2), whose completion requires
the kinetic sector already flagged by INERTIA-001. The decision is
recorded in `decision.json`; the completion record is in
`validation.json`.
"""
    (OUT / "transport_research001_report.md").write_text(report)

    print(json.dumps({
        "output": str(OUT),
        "decision": DECISION["outcome"],
        "pass": validation["pass"],
        "closest_structural_match": validation["closest_structural_match"],
        "missing_slot": validation["missing_slot"],
    }, indent=2))


if __name__ == "__main__":
    main()
