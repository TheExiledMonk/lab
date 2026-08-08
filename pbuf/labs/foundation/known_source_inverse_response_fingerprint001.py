#!/usr/bin/env python3
"""PBUF FOUNDATION — KNOWN SOURCE INVERSE RESPONSE FINGERPRINT 001.

Invert the mass->spacetime problem from the observable/effective side.

This lab does NOT claim to derive the microscopic substrate below the effective
spacetime/EM layer.  Instead it asks what response fingerprint a medium must
reproduce when the source mass and geometry are known independently by
construction.

Controlled uniform spheres are used so M, R and rho are inputs rather than
values inferred from gravitational motion.  The weak-field GR/Newton response
is used only as the empirical/macroscopic benchmark retained by PBUF, with G
explicitly labelled an observed effective coupling rather than a fundamental
microscopic constant.

No fitting, no cosmology, no Rmax, no lensing target and no microscopic claim.
"""
from __future__ import annotations

import json
import math
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
LAB_ID = "PBUF-FOUNDATION-KNOWN-SOURCE-INVERSE-RESPONSE-FINGERPRINT-001"
C = 299_792_458.0
G_MACRO = 6.67430e-11  # empirical macroscopic benchmark only


def git(*args: str) -> str:
    return subprocess.check_output(["git", *args], cwd=ROOT, text=True).strip()


def sphere(mass_kg: float, radius_m: float) -> dict:
    volume = 4.0 * math.pi * radius_m**3 / 3.0
    rho = mass_kg / volume
    # Effective weak-field response benchmarks.
    h00_surface = 2.0 * G_MACRO * mass_kg / (radius_m * C**2)
    q_rho = 8.0 * math.pi * G_MACRO * rho / C**2
    return {
        "M_kg": mass_kg,
        "R_m": radius_m,
        "rho_kg_m3": rho,
        "h00_surface": h00_surface,
        "q_rho_m-2": q_rho,
        "qR2": q_rho * radius_m**2,
        "identity_ratio_qR2_over_h00": (q_rho * radius_m**2) / h00_surface,
    }


def log_slope(x1: float, y1: float, x2: float, y2: float) -> float:
    return math.log(y2 / y1) / math.log(x2 / x1)


def main() -> int:
    branch = git("rev-parse", "--abbrev-ref", "HEAD")
    head = git("rev-parse", "HEAD")
    tracked = git("diff", "--name-only")
    staged = git("diff", "--name-only", "--cached")
    if tracked or staged:
        raise RuntimeError("tracked or staged repository changes present")

    # Base controlled source and orthogonal perturbations.
    base = sphere(1000.0, 0.50)
    mass2 = sphere(2000.0, 0.50)      # double M, fixed R
    radius2 = sphere(1000.0, 1.00)    # double R, fixed M
    density2 = sphere(8000.0, 1.00)   # same rho as base, doubled R

    # Independent probe radius outside every source.
    r_probe = 10.0
    def h00_far(src: dict) -> float:
        return 2.0 * G_MACRO * src["M_kg"] / (r_probe * C**2)

    fingerprints = {
        "surface_response_mass_exponent_fixed_R": log_slope(base["M_kg"], base["h00_surface"], mass2["M_kg"], mass2["h00_surface"]),
        "surface_response_radius_exponent_fixed_M": log_slope(base["R_m"], base["h00_surface"], radius2["R_m"], radius2["h00_surface"]),
        "local_curvature_density_exponent": log_slope(base["rho_kg_m3"], base["q_rho_m-2"], radius2["rho_kg_m3"], radius2["q_rho_m-2"]),
        "far_response_mass_exponent_fixed_probe_r": log_slope(base["M_kg"], h00_far(base), mass2["M_kg"], h00_far(mass2)),
        "same_density_radius_exponent_surface_response": log_slope(base["R_m"], base["h00_surface"], density2["R_m"], density2["h00_surface"]),
    }

    # Effective linear superposition benchmark at a field point. This is a
    # fingerprint to reproduce, not a microscopic interaction mechanism.
    M1, M2 = 1000.0, 3000.0
    r1, r2 = 10.0, 14.0
    h1 = 2.0 * G_MACRO * M1 / (r1 * C**2)
    h2 = 2.0 * G_MACRO * M2 / (r2 * C**2)
    hsum = h1 + h2
    superposition = {
        "h1": h1,
        "h2": h2,
        "linear_sum": hsum,
        "microscopic_reason_claimed": False,
        "role": "effective weak-field response fingerprint only",
    }

    checks = {
        "controlled_source_mass_not_gravity_inferred": True,
        "G_used_only_as_macroscopic_empirical_benchmark": True,
        "gravity_fundamental_in_PBUF": False,
        "microscopic_substrate_claimed": False,
        "fit_or_tuning_used": False,
        "Rmax_used": False,
        "replacement_free_parameter_introduced": False,
        "cosmology_used": False,
        "lensing_target_used": False,
        "legacy_0p18_used": False,
        "mass_exponent_is_one": abs(fingerprints["surface_response_mass_exponent_fixed_R"] - 1.0) < 1e-12,
        "radius_exponent_fixed_mass_is_minus_one": abs(fingerprints["surface_response_radius_exponent_fixed_M"] + 1.0) < 1e-12,
        "local_density_exponent_is_one": abs(fingerprints["local_curvature_density_exponent"] - 1.0) < 1e-12,
        "uniform_sphere_identity_qR2_equals_3h00": all(abs(s["identity_ratio_qR2_over_h00"] - 3.0) < 1e-12 for s in (base, mass2, radius2, density2)),
        "no_tracked_or_staged_changes": not tracked and not staged,
        "stdout_only_no_run_directory_created": True,
    }

    conclusion = {
        "status": "EFFECTIVE_RESPONSE_FINGERPRINT_DERIVED_MICROSCOPIC_CAUSE_SPECULATIVE",
        "effective_findings": [
            "for a controlled uniform source, local weak-field source curvature scales linearly with rho",
            "surface metric response scales as M/R",
            "far-field response at fixed probe radius scales linearly with total M",
            "at fixed density, surface response scales as R^2 because M scales as R^3",
            "weak-field benchmark is additive to first order",
        ],
        "inverse_constraint": "any PBUF medium law must reproduce these effective scalings in the weak-field regime before its microscopic interpretation is considered viable",
        "speculative_boundary": "the mechanism by which stress-energy produces the effective spacetime-medium response is not directly probed here and remains speculative",
        "safe_next": "compare the existing native PBUF local and accumulated fields against this frozen response fingerprint without rescaling their amplitude; identify which native variable reproduces rho-local loading, M/R surface accumulation, M/r far-field scaling, and linear weak-field superposition. Treat any microscopic explanation as speculative until independently constrained.",
    }

    payload = {
        "lab_id": LAB_ID,
        "status": conclusion["status"],
        "repo_state": {"branch": branch, "head_sha": head, "tracked_changes": tracked, "staged_changes": staged},
        "policy": {
            "epistemic_rule": "derive_effective_response_speculate_only_about_unprobeable_cause",
            "G_role": "measured_macroscopic_response_benchmark_not_fundamental_PBUF_input",
        },
        "controlled_sources": {"base": base, "double_mass_fixed_R": mass2, "double_radius_fixed_M": radius2, "same_density_double_radius": density2},
        "probe_radius_m": r_probe,
        "fingerprints": fingerprints,
        "superposition": superposition,
        "conclusion": conclusion,
        "checks": checks,
    }

    print(LAB_ID)
    print(f"status={conclusion['status']}")
    print(f"head_sha={head}")
    print("gravity_fundamental_in_PBUF=false")
    print("G_role=MACROSCOPIC_EMPIRICAL_BENCHMARK_ONLY")
    print("microscopic_substrate_claimed=false")
    print("fit_or_tuning_used=false")
    print()
    print("INVERSE_RESPONSE_FINGERPRINT")
    for k, v in fingerprints.items():
        print(f"{k}={v:.12g}")
    print("uniform_sphere_qR2_over_h00=3")
    print("weak_field_superposition=ADDITIVE_EFFECTIVE_BENCHMARK")
    print()
    print("CONCLUSION")
    print(f"status={conclusion['status']}")
    print(f"inverse_constraint={conclusion['inverse_constraint']}")
    print(f"speculative_boundary={conclusion['speculative_boundary']}")
    print(f"safe_next={conclusion['safe_next']}")
    print()
    print("CHECKS")
    for k, v in checks.items():
        print(f"{k}={str(v).lower()}")
    print(json.dumps(payload, sort_keys=True))

    if not all(checks.values()):
        raise RuntimeError("inverse-response fingerprint gate failed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
