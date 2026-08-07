#!/usr/bin/env python3
"""PBUF FOUNDATION — MASS / MEDIUM CANDIDATE AUDIT 001.

Fact-finding only.

Question under investigation
----------------------------
Given a known baryonic rest mass, what microscopic source length/volume and
medium modulus could turn that physical source into an absolute PBUF metric
strain without fitting a replacement for the historical dimensionless
``strength = 0.18``?

This lab DOES NOT attempt to decide that question.  It evaluates several
candidate dimensional bridges side-by-side and reports what each one implies.
No candidate is selected, ranked, fitted, or tuned.

Existing PBUF baryon/alpha structure is deliberately left untouched.  This lab
starts *after* that model structure and probes only the missing absolute
length/volume/modulus bridge.

Candidate source scales
-----------------------
1. Planck cell: one Planck length in each spatial dimension, V = l_P^3.
   This is an intentionally extreme boundary case if the full baryon rest
   energy is assigned to one cell; it is NOT asserted as physical.
2. Reduced Compton length: hbar/(m c).
3. Ordinary Compton wavelength: h/(m c).
4. Measured proton charge-radius scale.
5. Two-scale bookkeeping: source volumes above resolved into Planck cells.

Candidate medium modulus
------------------------
As a dimensional fact-finding reference only, construct

    E_P* = hbar c / l_P
    K_P* = E_P* / l_P^3 = hbar c / l_P^4

where l_P is taken as the supplied Planck-scale length.  K_P* has units J/m^3,
the same dimensions as pressure / bulk modulus / energy density.  Calling it a
PBUF spacetime modulus would require a physical derivation that this lab does
NOT provide.

For every candidate source volume V, the lab reports

    E_b      = m_b c^2
    u_b      = E_b / V
    eps_KP   = u_b / K_P*

and the number of Planck cubes N_P = V/l_P^3.  eps_KP is therefore a hypothesis
output under the single stated modulus candidate, not a prediction.

Current-foundation compatibility
--------------------------------
The resulting scalar eps_KP is also inserted, without rescaling, into the
already-adopted effective metric-strain normalization

    g_mu_nu = gbar_mu_nu + 2 chi_mu_nu

using a purely isotropic spatial diagnostic chi=diag(0,eps,eps,eps).  This tests
numerical/signature compatibility only; it is not a lensing calculation and it
does not claim that a baryonic source produces isotropic local chi.

Hard guardrails
---------------
- no observed kappa or shear;
- no HST/photometric target data;
- no GR/Newtonian force law used to calibrate a candidate;
- no Quantum Engine input;
- no fitted/tuned scalar;
- no candidate selection by agreement with the legacy 0.18 value;
- legacy 0.18 appears only as a quarantined scale comparison.
"""
from __future__ import annotations

import json
import math
import subprocess
import time
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "runs" / "mass_medium_candidate_audit001"
LAB_ID = "PBUF-FOUNDATION-MASS-MEDIUM-CANDIDATE-AUDIT-001"

# Exact SI constants where the SI fixes them exactly.
C_M_S = 299_792_458.0
H_J_S = 6.626_070_15e-34
HBAR_J_S = H_J_S / (2.0 * math.pi)

# Measured/reference quantities.  These are explicit inputs, not fit here.
# Values are suitable for order-of-magnitude candidate testing; this lab does
# not propagate their metrology uncertainties into a physics conclusion.
PLANCK_LENGTH_M = 1.616_255e-35
PROTON_MASS_KG = 1.672_621_925_95e-27
PROTON_CHARGE_RADIUS_M = 0.8409e-15

# Historical repository coefficient: diagnostic comparison ONLY.
LEGACY_STRENGTH_REFERENCE = 0.18


@dataclass(frozen=True)
class CandidateScale:
    key: str
    physical_role: str
    length_m: float
    volume_geometry: str
    volume_m3: float
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


def _cube_volume(length_m: float) -> float:
    return length_m**3


def _sphere_volume(radius_m: float) -> float:
    return (4.0 / 3.0) * math.pi * radius_m**3


def _candidate_scales(mass_kg: float) -> list[CandidateScale]:
    reduced_compton = HBAR_J_S / (mass_kg * C_M_S)
    compton = H_J_S / (mass_kg * C_M_S)
    rp = PROTON_CHARGE_RADIUS_M
    lp = PLANCK_LENGTH_M
    return [
        CandidateScale(
            key="planck_cell_cube",
            physical_role="elementary medium-cell boundary hypothesis",
            length_m=lp,
            volume_geometry="cube_l_cubed",
            volume_m3=_cube_volume(lp),
            status="speculative_boundary_case",
        ),
        CandidateScale(
            key="reduced_compton_cube",
            physical_role="quantum localization/source-extent candidate",
            length_m=reduced_compton,
            volume_geometry="cube_l_cubed",
            volume_m3=_cube_volume(reduced_compton),
            status="speculative_candidate",
        ),
        CandidateScale(
            key="reduced_compton_sphere",
            physical_role="same reduced-Compton scale with spherical geometry",
            length_m=reduced_compton,
            volume_geometry="sphere_radius_l",
            volume_m3=_sphere_volume(reduced_compton),
            status="speculative_candidate_geometry_check",
        ),
        CandidateScale(
            key="compton_cube",
            physical_role="ordinary de Broglie/Compton source-extent candidate",
            length_m=compton,
            volume_geometry="cube_l_cubed",
            volume_m3=_cube_volume(compton),
            status="speculative_candidate",
        ),
        CandidateScale(
            key="compton_sphere",
            physical_role="same ordinary Compton scale with spherical geometry",
            length_m=compton,
            volume_geometry="sphere_radius_l",
            volume_m3=_sphere_volume(compton),
            status="speculative_candidate_geometry_check",
        ),
        CandidateScale(
            key="proton_charge_radius_cube",
            physical_role="measured baryon charge-size scale used only as a length candidate",
            length_m=rp,
            volume_geometry="cube_l_cubed",
            volume_m3=_cube_volume(rp),
            status="measured_length_speculative_response_volume",
        ),
        CandidateScale(
            key="proton_charge_radius_sphere",
            physical_role="measured baryon charge-size scale with spherical geometry",
            length_m=rp,
            volume_geometry="sphere_radius_l",
            volume_m3=_sphere_volume(rp),
            status="measured_length_speculative_response_volume",
        ),
    ]


def _metric_strain_diagnostic(eps: float) -> dict:
    """Insert eps into current effective metric-strain normalization only.

    This is deliberately NOT a source law.  It checks what an isotropic spatial
    chi of the candidate amplitude would do to the local metric representation.
    """
    gbar = np.diag([-1.0, 1.0, 1.0, 1.0])
    chi = np.diag([0.0, eps, eps, eps])
    g = gbar + 2.0 * chi
    eig = np.linalg.eigvalsh(g)
    lorentzian_signature = bool(np.count_nonzero(eig < 0.0) == 1 and np.count_nonzero(eig > 0.0) == 3)
    return {
        "chi_isotropic_spatial_diag": [0.0, eps, eps, eps],
        "metric_diag": [float(g[i, i]) for i in range(4)],
        "metric_determinant": float(np.linalg.det(g)),
        "metric_eigenvalues": [float(x) for x in eig],
        "lorentzian_signature_preserved": lorentzian_signature,
        "max_abs_metric_delta": float(np.max(np.abs(g - gbar))),
    }


def _evaluate_candidates() -> dict:
    lp = PLANCK_LENGTH_M
    vp = lp**3
    baryon_energy = PROTON_MASS_KG * C_M_S**2

    # Dimensional Planck-cell modulus candidate.  The star is intentional: this
    # is a constructed reference scale, not a claimed PBUF constitutive modulus.
    ep_star = HBAR_J_S * C_M_S / lp
    kp_star = ep_star / vp

    rows = []
    for c in _candidate_scales(PROTON_MASS_KG):
        energy_density = baryon_energy / c.volume_m3
        eps_kp = energy_density / kp_star
        n_planck_cubes = c.volume_m3 / vp
        energy_per_planck_cube_uniform = baryon_energy / n_planck_cubes
        fraction_planck_energy_per_cell = energy_per_planck_cube_uniform / ep_star

        # Algebraic identity check: both routes are the same uniform-density
        # statement under K_P*=E_P*/V_P.
        identity_rel = abs(eps_kp - fraction_planck_energy_per_cell) / max(
            abs(eps_kp), abs(fraction_planck_energy_per_cell), 1.0e-300
        )

        rows.append({
            **asdict(c),
            "baryon_rest_energy_J": baryon_energy,
            "energy_density_J_m3": energy_density,
            "planck_cube_count_in_candidate_volume": n_planck_cubes,
            "uniform_energy_per_planck_cube_J": energy_per_planck_cube_uniform,
            "uniform_energy_per_cell_over_Ep_star": fraction_planck_energy_per_cell,
            "candidate_strain_if_K_equals_Kp_star": eps_kp,
            "log10_abs_candidate_strain": math.log10(abs(eps_kp)) if eps_kp != 0.0 else -math.inf,
            "legacy_0p18_over_candidate_strain": LEGACY_STRENGTH_REFERENCE / eps_kp if eps_kp != 0.0 else math.inf,
            "candidate_strain_over_legacy_0p18": eps_kp / LEGACY_STRENGTH_REFERENCE,
            "uniform_cell_identity_relative_error": identity_rel,
            "metric_strain_map_diagnostic": _metric_strain_diagnostic(eps_kp),
        })

    return {
        "reference_quantities": {
            "proton_mass_kg": PROTON_MASS_KG,
            "proton_rest_energy_J": baryon_energy,
            "planck_length_m": lp,
            "planck_cube_volume_m3": vp,
            "proton_reduced_compton_length_m": HBAR_J_S / (PROTON_MASS_KG * C_M_S),
            "proton_compton_wavelength_m": H_J_S / (PROTON_MASS_KG * C_M_S),
            "proton_charge_radius_m": PROTON_CHARGE_RADIUS_M,
            "Ep_star_J_equals_hbar_c_over_lP": ep_star,
            "Kp_star_J_m3_equals_Ep_star_over_lP3": kp_star,
            "Kp_star_status": "dimensional_candidate_only_not_derived_PBUF_modulus",
        },
        "candidate_rows": rows,
    }


def _checks(result: dict) -> dict:
    rows = result["candidate_rows"]
    finite = all(
        math.isfinite(float(r["volume_m3"]))
        and math.isfinite(float(r["energy_density_J_m3"]))
        and math.isfinite(float(r["candidate_strain_if_K_equals_Kp_star"]))
        and float(r["volume_m3"]) > 0.0
        for r in rows
    )
    identity = all(float(r["uniform_cell_identity_relative_error"]) <= 1.0e-12 for r in rows)
    signature = all(bool(r["metric_strain_map_diagnostic"]["lorentzian_signature_preserved"]) for r in rows)
    no_candidate_selected = True
    return {
        "all_candidate_numbers_finite_and_positive_volume": finite,
        "uniform_density_planck_cell_identity_pass": identity,
        "current_metric_strain_map_signature_preserved_for_all_candidates": signature,
        "no_candidate_selected_or_ranked": no_candidate_selected,
        "no_fit_performed": True,
        "target_blind_no_kappa_shear_or_HST_loaded": True,
        "quantum_engine_not_used": True,
        "existing_pbuf_baryon_alpha_structure_untouched": True,
        "legacy_strength_used_only_as_quarantined_scale_reference": True,
    }


def _write_text_report(payload: dict) -> None:
    ref = payload["fact_finding"]["reference_quantities"]
    rows = payload["fact_finding"]["candidate_rows"]
    lines = [
        LAB_ID,
        "FACT-FINDING ONLY — NO CANDIDATE SELECTED",
        "",
        f"head_sha={payload['repo_state']['head_sha']}",
        f"proton_rest_energy_J={ref['proton_rest_energy_J']:.17e}",
        f"planck_length_m={ref['planck_length_m']:.17e}",
        f"Ep_star_J={ref['Ep_star_J_equals_hbar_c_over_lP']:.17e}",
        f"Kp_star_J_m3={ref['Kp_star_J_m3_equals_Ep_star_over_lP3']:.17e}",
        "",
        "candidate | L[m] | V[m^3] | N_P | u[J/m^3] | eps(u/Kp*) | eps/0.18",
    ]
    for r in rows:
        lines.append(
            f"{r['key']} | {r['length_m']:.17e} | {r['volume_m3']:.17e} | "
            f"{r['planck_cube_count_in_candidate_volume']:.17e} | "
            f"{r['energy_density_J_m3']:.17e} | "
            f"{r['candidate_strain_if_K_equals_Kp_star']:.17e} | "
            f"{r['candidate_strain_over_legacy_0p18']:.17e}"
        )
    lines.extend([
        "",
        "Interpretation guardrail:",
        "eps(u/Kp*) is what follows IF the source volume and Kp* hypotheses are both adopted.",
        "It is not a PBUF prediction and is not calibrated to 0.18 or to lensing data.",
        "Existing PBUF baryon/alpha structure was not altered or re-derived here.",
        "",
        "Next use: inspect orders of magnitude, reject dimensionally/numerically implausible lanes,",
        "then promote only physically motivated survivors into a separate target-blind source-field test.",
    ])
    (OUT / "report.txt").write_text("\n".join(lines) + "\n")


def main() -> int:
    t0 = time.time()
    OUT.mkdir(parents=True, exist_ok=True)
    state_before = _repo_state()
    fact = _evaluate_candidates()
    checks = _checks(fact)

    payload = {
        "lab_id": LAB_ID,
        "status": "FACT_FINDING_ONLY",
        "question": "candidate absolute baryonic mass-energy to spacetime-medium strain bridge",
        "guardrails": {
            "fact_finding_only": True,
            "no_fit": True,
            "target_blind": True,
            "no_kappa": True,
            "no_shear": True,
            "no_HST_target_data": True,
            "no_GR_Newtonian_calibration": True,
            "no_quantum_engine": True,
            "existing_pbuf_baryon_alpha_untouched": True,
            "legacy_strength_0p18_role": "quarantined_scale_comparison_only",
        },
        "hypothesis_boundary": {
            "existing_PBUF_side": "baryon/alpha/multiplicity and current metric-strain framework are treated as existing inputs",
            "speculative_side": "choice of source response length/volume and Kp* dimensional modulus candidate",
            "candidate_selection_authorized": False,
        },
        "repo_state": state_before,
        "fact_finding": fact,
        "checks": checks,
        "runtime_seconds": time.time() - t0,
    }

    (OUT / "results.json").write_text(json.dumps(payload, indent=2, allow_nan=False) + "\n")
    _write_text_report(payload)

    state_after = _repo_state()
    # Runs/ is expected to be untracked/ignored; still refuse success if tracked
    # or staged repository content was modified by the lab itself.
    clean = not state_after["tracked_changes"] and not state_after["staged_changes"]
    payload["post_run_repo_state"] = state_after
    payload["checks"]["no_tracked_or_staged_changes_after_run"] = clean
    (OUT / "results.json").write_text(json.dumps(payload, indent=2, allow_nan=False) + "\n")

    print(json.dumps(payload, indent=2, allow_nan=False))
    return 0 if all(bool(v) for v in payload["checks"].values()) else 1


if __name__ == "__main__":
    raise SystemExit(main())
