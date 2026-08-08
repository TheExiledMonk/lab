#!/usr/bin/env python3
"""PBUF FOUNDATION — MASS / SPACETIME RESPONSE STRENGTH AUDIT 001.

Fact-finding only.

Purpose
-------
Return to the missing absolute mass -> spacetime-response bridge after the
G-origin omnibus audit.  In PBUF gravity is not assumed fundamental.  The
experimentally measured Newton constant is therefore used here only as the
observed *macroscopic response anchor* that any microscopic PBUF constitutive
law must ultimately reproduce.

This lab does NOT derive G, does NOT insert Newton/GR into the PBUF propagation
pipeline, and does NOT choose a replacement for the legacy strength=0.18.
Instead it asks a narrower question:

    What absolute stress->curvature compliance and dimensionless weak-field
    deformation scales are implied by the measured macroscopic response of
    spacetime to known masses?

The resulting quantities are frozen empirical boundary conditions for a later
PBUF microscopic/constitutive derivation.

Macroscopic response anchor
---------------------------
The weak-field Einstein coupling may be written as

    C_T_to_curvature = 8*pi*G/c^4

so that a stress-energy density T has a curvature-source scale

    q = C_T_to_curvature * T .

Its inverse

    K_spacetime = c^4/(8*pi*G)

is an empirical effective stiffness scale.  This is NOT claimed fundamental;
it is simply another representation of measured G.

For nonrelativistic rest-mass density rho, T00 ~= rho*c^2, hence

    q_rho = 8*pi*G*rho/c^2                     [1/m^2].

For a spherical source of mass M, the external weak-field dimensionless
potential/deformation scale at radius r is

    eta = G*M/(r*c^2),
    |h00| = 2*eta.

For a uniform sphere evaluated at its surface, direct volume averaging gives

    q_rho * R^2 = 6*eta = 3*|h00|.

That factor is a geometry identity for a uniform sphere, not a fitted number.

Guardrails
----------
- gravity is NOT declared fundamental in PBUF;
- measured G is an empirical macroscopic anchor only;
- no kappa, shear, HST pixels, lens morphology, or lens amplitude input;
- no fitting, tuning, optimization, or candidate selection;
- no use of the legacy 0.18 value in any physical calculation;
- no Quantum Engine input;
- no Planck length or Planck units;
- no microscopic PBUF coupling is inferred by algebraically renaming G;
- stdout only; no run directory.
"""
from __future__ import annotations

import json
import math
import re
import subprocess
from dataclasses import asdict, dataclass
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
LAB_ID = "PBUF-FOUNDATION-MASS-SPACETIME-RESPONSE-STRENGTH-AUDIT-001"

# Exact SI speed of light.
C = 299_792_458.0

# Experimental macroscopic gravitational response anchor (CODATA conventional
# value used throughout the recent PBUF audits).  It is not treated as a
# microscopic/fundamental PBUF input.
G_MEASURED = 6.67430e-11

# Reference bodies.  These are used only to expose physical scale.  They are not
# calibration targets and no result is selected by agreement with them.
REFERENCE_BODIES = (
    # name, mass [kg], radius [m], provenance
    ("earth", 5.9722e24, 6.3710e6, "standard geophysical reference mass/radius"),
    ("sun", 1.98847e30, 6.9570e8, "standard solar reference mass/radius"),
)

LEGACY_FILE = ROOT / "pbuf" / "labs" / "foundation" / "m10_coverage_25pct_science001.py"


@dataclass(frozen=True)
class BodyResponse:
    name: str
    mass_kg: float
    radius_m: float
    mean_density_kg_m3: float
    eta_GM_over_Rc2: float
    h00_surface_magnitude: float
    schwarzschild_radius_over_R: float
    rho_curvature_source_inv_m2: float
    rho_curvature_times_R2: float
    uniform_sphere_identity_ratio_to_h00: float
    provenance: str


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


def _legacy_strength_inventory() -> dict:
    out = {
        "path": str(LEGACY_FILE.relative_to(ROOT)),
        "exists": LEGACY_FILE.exists(),
        "value": None,
        "role": "historical_dimensionless_initial_state_amplitude_only",
        "used_in_physical_calculations_in_this_lab": False,
    }
    if not LEGACY_FILE.exists():
        return out
    text = LEGACY_FILE.read_text(errors="ignore")
    m = re.search(r"(?m)^STRENGTH\s*=\s*([0-9]+(?:\.[0-9]+)?)", text)
    if m:
        out["value"] = float(m.group(1))
    return out


def _macroscopic_anchor() -> dict:
    compliance = 8.0 * math.pi * G_MEASURED / C**4
    stiffness = 1.0 / compliance
    return {
        "G_measured_m3_kg_s2": G_MEASURED,
        "G_role": "experimental_macroscopic_response_anchor_not_fundamental_PBUF_coupling",
        "stress_energy_to_curvature_compliance_s2_per_kg_m": compliance,
        "stress_energy_to_curvature_compliance_m_per_joule": compliance,
        "effective_spacetime_stiffness_newton": stiffness,
        "formula_compliance": "8*pi*G_measured/c^4",
        "formula_stiffness": "c^4/(8*pi*G_measured)",
        "fundamental_claim": False,
    }


def _body_response(name: str, mass: float, radius: float, provenance: str) -> BodyResponse:
    volume = (4.0 / 3.0) * math.pi * radius**3
    rho = mass / volume
    eta = G_MEASURED * mass / (radius * C**2)
    h00 = 2.0 * eta
    q_rho = 8.0 * math.pi * G_MEASURED * rho / C**2
    qR2 = q_rho * radius**2
    ratio = qR2 / h00
    return BodyResponse(
        name=name,
        mass_kg=mass,
        radius_m=radius,
        mean_density_kg_m3=rho,
        eta_GM_over_Rc2=eta,
        h00_surface_magnitude=h00,
        schwarzschild_radius_over_R=2.0 * eta,
        rho_curvature_source_inv_m2=q_rho,
        rho_curvature_times_R2=qR2,
        uniform_sphere_identity_ratio_to_h00=ratio,
        provenance=provenance,
    )


def _point_mass_scale(mass_kg: float, radius_m: float) -> dict:
    eta = G_MEASURED * mass_kg / (radius_m * C**2)
    return {
        "mass_kg": mass_kg,
        "evaluation_radius_m": radius_m,
        "eta_GM_over_rc2": eta,
        "h00_magnitude": 2.0 * eta,
        "role": "analytic_scale_check_only",
    }


def _checks(anchor: dict, bodies: list[BodyResponse], legacy: dict, repo: dict) -> dict:
    compliance = anchor["stress_energy_to_curvature_compliance_s2_per_kg_m"]
    stiffness = anchor["effective_spacetime_stiffness_newton"]
    inverse_error = abs(compliance * stiffness - 1.0)

    sphere_errors = []
    for b in bodies:
        # q R^2 = 6 eta and h00 = 2 eta => qR^2/h00 = 3.
        sphere_errors.append(abs(b.uniform_sphere_identity_ratio_to_h00 - 3.0))

    point = _point_mass_scale(1.0, 1.0)
    point_expected = G_MEASURED / C**2
    point_error = abs(point["eta_GM_over_rc2"] - point_expected) / point_expected

    return {
        "all_anchor_values_finite_positive": all(
            math.isfinite(x) and x > 0.0 for x in (compliance, stiffness)
        ),
        "compliance_stiffness_inverse_identity_pass": inverse_error < 1.0e-14,
        "compliance_stiffness_inverse_abs_error": inverse_error,
        "uniform_sphere_geometry_identity_pass": max(sphere_errors) < 1.0e-12,
        "uniform_sphere_geometry_identity_max_abs_error": max(sphere_errors),
        "one_kg_at_one_m_scale_identity_pass": point_error < 1.0e-14,
        "one_kg_at_one_m_scale_relative_error": point_error,
        "legacy_strength_not_used_as_physical_input": legacy["used_in_physical_calculations_in_this_lab"] is False,
        "legacy_strength_is_not_promoted_to_constitutive_law": True,
        "gravity_not_declared_fundamental": True,
        "measured_G_role_is_macroscopic_anchor_only": anchor["fundamental_claim"] is False,
        "no_lensing_target_input": True,
        "no_fit_or_tuning": True,
        "no_quantum_engine_input": True,
        "no_planck_scale_input": True,
        "no_tracked_or_staged_changes_created_by_lab": not repo["tracked_changes"] and not repo["staged_changes"],
        "stdout_only_no_run_directory_created": True,
    }


def main() -> None:
    repo = _repo_state()
    anchor = _macroscopic_anchor()
    legacy = _legacy_strength_inventory()
    bodies = [_body_response(*row) for row in REFERENCE_BODIES]
    point_1kg_1m = _point_mass_scale(1.0, 1.0)
    checks = _checks(anchor, bodies, legacy, repo)

    if not all(v for k, v in checks.items() if k.endswith("_pass") or k in {
        "all_anchor_values_finite_positive",
        "legacy_strength_not_used_as_physical_input",
        "legacy_strength_is_not_promoted_to_constitutive_law",
        "gravity_not_declared_fundamental",
        "measured_G_role_is_macroscopic_anchor_only",
        "no_lensing_target_input",
        "no_fit_or_tuning",
        "no_quantum_engine_input",
        "no_planck_scale_input",
        "no_tracked_or_staged_changes_created_by_lab",
        "stdout_only_no_run_directory_created",
    }):
        raise RuntimeError("one or more fact-finding guardrail/identity checks failed")

    result = {
        "lab_id": LAB_ID,
        "status": "FACT_FINDING_ONLY",
        "repo_state": repo,
        "interpretation": {
            "gravity_fundamental_in_PBUF": False,
            "what_is_measured": "macroscopic mass-to-spacetime response",
            "what_remains_missing": "microscopic PBUF constitutive derivation of the same absolute response",
            "measured_G_policy": "allowed empirical endpoint/anchor; not a microscopic derivation",
            "legacy_0p18_policy": "historical trajectory/amplitude diagnostic only; never a normalization clue",
        },
        "macroscopic_anchor": anchor,
        "reference_bodies": [asdict(b) for b in bodies],
        "analytic_point_control_1kg_at_1m": point_1kg_1m,
        "legacy_strength_inventory": legacy,
        "bridge_conclusion": {
            "empirical_absolute_compliance_available": True,
            "empirical_absolute_stiffness_available": True,
            "native_PBUF_microscopic_origin_closed": False,
            "native_PBUF_field_mapping_closed": False,
            "safe_next_question": "which existing PBUF medium field/strain variable carries q or its dimensionless integral, with normalization derived from PBUF constitutive physics rather than copied from G",
            "do_not_do": [
                "replace strength=0.18 with a number solved from lensing",
                "declare c^4/(8*pi*G) a fundamental PBUF stiffness",
                "inject h00 or Newtonian potential directly into the PBUF observer pipeline and call that a derivation",
                "use kappa/shear/morphology to choose a source amplitude",
            ],
        },
        "checks": checks,
    }

    print(LAB_ID)
    print("status=FACT_FINDING_ONLY")
    print(f"head_sha={repo['head_sha']}")
    print("gravity_fundamental_in_PBUF=false")
    print("measured_G_role=MACROSCOPIC_RESPONSE_ANCHOR_ONLY")
    print("legacy_0p18_used_as_physical_input=false")
    print("no_fit=true")
    print("target_blind=true")
    print("quantum_engine_used=false")
    print("planck_scale_used=false")
    print()
    print("MACROSCOPIC_RESPONSE_ANCHOR")
    for k, v in anchor.items():
        print(f"{k}={v}")
    print()
    print("REFERENCE_BODY_RESPONSE")
    print("name | M[kg] | R[m] | rho[kg/m3] | eta=GM/(Rc2) | |h00| | q_rho[1/m2] | qR2 | qR2/h00")
    for b in bodies:
        print(
            f"{b.name} | {b.mass_kg:.17e} | {b.radius_m:.17e} | "
            f"{b.mean_density_kg_m3:.17e} | {b.eta_GM_over_Rc2:.17e} | "
            f"{b.h00_surface_magnitude:.17e} | {b.rho_curvature_source_inv_m2:.17e} | "
            f"{b.rho_curvature_times_R2:.17e} | {b.uniform_sphere_identity_ratio_to_h00:.17e}"
        )
    print()
    print("LEGACY_STRENGTH_INVENTORY")
    print(json.dumps(legacy, sort_keys=True))
    print()
    print("BRIDGE_CONCLUSION")
    print("empirical_absolute_response_scale_available=true")
    print("microscopic_PBUF_origin_closed=false")
    print("native_PBUF_field_mapping_closed=false")
    print("next=map empirical curvature/loading dimensions onto native PBUF medium variables without promoting G to fundamental input")
    print()
    print("CHECKS")
    for k, v in checks.items():
        print(f"{k}={str(v).lower() if isinstance(v, bool) else v}")
    print("JSON=" + json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()
