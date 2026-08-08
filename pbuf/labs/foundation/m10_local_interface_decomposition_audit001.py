#!/usr/bin/env python3
"""PBUF FOUNDATION — M10 LOCAL INTERFACE DECOMPOSITION AUDIT 001.

Purpose
-------
Trace exactly where the historical unit-loading M10 propagation interface gets
its ~0.008*c_state scale, and compare that scale with local bond quantities in
the native accumulated bounded-strain medium.

This is an audit only. It does NOT construct a replacement propagation mapping.
It does NOT divide by 360, fit a coefficient, normalize amplitudes, rescale the
native field, or use observed lensing values.

Historical unit route decomposed here:
    rho3 -> unit A8 state
         -> positive-N6 slow/fast neighbour differences
         -> frozen T1 pair amplitudes
              A_ij = (dt*omega*K) * Delta u_fast
                   + (dt*tau_slow) * Delta u_slow
         -> PM1/PS2 pair response
         -> M10 midpoint rasterisation

Native accumulated medium inspected here:
    rho3 -> zero-flux raw c_state
         -> bounded-strain accumulated equilibrium u
         -> positive-N6 bond strain Delta u
         -> bounded-strain bond traction
              sigma = K0*Delta u / (1-(Delta u/epsilon_max)^2)

The native bond quantities are reported beside M10 only to identify structural
scale and locality. They are NOT fed into G3D and are NOT claimed to be the
correct propagation interface.
"""
from __future__ import annotations

import json
import math
import subprocess
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))

from a8_three_dimensional_projection_lab001 import construct_rho_3d
from pbuf.core import benchmark_data as BENCH
from pbuf.core import pair_enumeration as M05
from pbuf.core import pair_transfer as M08
from pbuf.core import midpoint_rasterization as M10
from pbuf.models import a8_state as A8
from pbuf.models import a8_pair_amplitude as PAIR
from pbuf.models import transverse_projector as PROJ
import pbuf.labs.foundation.m10_coverage_25pct_science001 as BASE
import pbuf.labs.foundation.native_accumulated_full_lensing_local_benchmark001 as LOCAL
import pbuf.labs.foundation.native_accumulated_full_lensing001 as FULL

LAB_ID = "PBUF-FOUNDATION-M10-LOCAL-INTERFACE-DECOMPOSITION-AUDIT-001"
EXPECTED_CLUSTER_IDS = (
    "Abell2744", "MACS0416", "MACS1149", "AbellS1063", "Abell370"
)


def git(*args: str) -> str:
    return subprocess.check_output(["git", *args], cwd=ROOT, text=True).strip()


def repo_state() -> dict:
    return {
        "repository": "TheExiledMonk/lab",
        "branch": git("rev-parse", "--abbrev-ref", "HEAD"),
        "head_sha": git("rev-parse", "HEAD"),
        "tracked_changes": git("diff", "--name-only"),
        "staged_changes": git("diff", "--name-only", "--cached"),
    }


def rms(a) -> float:
    x = np.asarray(a, dtype=np.float64)
    return float(np.sqrt(np.mean(x * x))) if x.size else float("nan")


def ratio(a: float, b: float) -> float:
    return float(a / max(abs(b), 1.0e-30))


def vector_rms(vector) -> float:
    comps = [np.asarray(v, dtype=np.float64) for v in vector]
    mag2 = np.zeros_like(comps[0])
    for c in comps:
        mag2 += c * c
    return float(np.sqrt(np.mean(mag2)))


def concat_positive_bonds(field: np.ndarray) -> np.ndarray:
    """Return all positive-N6 nearest-neighbour differences j-i."""
    f = np.asarray(field, dtype=np.float64)
    dx = (f[:, :, 1:] - f[:, :, :-1]).ravel()
    dy = (f[:, 1:, :] - f[:, :-1, :]).ravel()
    dz = (f[1:, :, :] - f[:-1, :, :]).ravel()
    return np.concatenate((dx, dy, dz))


def response_slot_vector_rms(response: dict, shape: tuple[int, int, int]) -> float:
    """RMS magnitude over the three positive-direction pair-response slots."""
    mags = []
    for keys, axis in (
        (("R_ij_xp", "R_ij_y_xp", "R_ij_z_xp"), "xp"),
        (("R_ij_yp", "R_ij_y_yp", "R_ij_z_yp"), "yp"),
        (("R_ij_zp", "R_ij_y_zp", "R_ij_z_zp"), "zp"),
    ):
        rx, ry, rz = (np.asarray(response[k], dtype=np.float64) for k in keys)
        if axis == "xp":
            sl = (slice(None), slice(None), slice(0, shape[2] - 1))
        elif axis == "yp":
            sl = (slice(None), slice(0, shape[1] - 1), slice(None))
        else:
            sl = (slice(0, shape[0] - 1), slice(None), slice(None))
        mags.append(np.sqrt(rx[sl] ** 2 + ry[sl] ** 2 + rz[sl] ** 2).ravel())
    return rms(np.concatenate(mags))


def local_rho3(cluster: dict) -> np.ndarray:
    kappa = BENCH.load_kappa(cluster)
    rho2 = BASE.construct_common_proxy(
        kappa, bins=BASE.OBS_BINS, extent=BASE.CFG["extent"]
    )
    return np.asarray(
        construct_rho_3d(rho2, BASE.NZ, profile=BASE.PROFILE),
        dtype=np.float64,
    )


def unit_decomposition(rho3: np.ndarray) -> dict:
    """Rebuild the frozen unit-loading A8 -> pair -> M10 chain explicitly."""
    rng = np.random.RandomState(BASE.SEED)
    eq = np.asarray(rho3, dtype=np.float64)
    noise = A8.A8_INIT_INJECTION_NOISE * rng.randn(*rho3.shape)
    initial = {"rho_3d": rho3.copy(), "u_slow0": eq.copy(), "u_fast0": eq + noise}
    state = BASE._evolve(initial)

    us = np.asarray(state["u_slow"], dtype=np.float64)
    uf = np.asarray(state["u_fast"], dtype=np.float64)
    c = np.asarray(state["c_state"], dtype=np.float64)
    shape = tuple(c.shape)

    pairs = M05.enumerate_internal_pairs(shape)
    ex, ey, ez, valid, gmag = PROJ.build_longitudinal_direction(c)
    projector = PROJ.build_transverse_projector(ex, ey, ez)
    amps = PAIR.compute_a8_pair_amplitudes(us, uf, c, pairs)
    response = M08.build_pair_responses(
        pairs, amps, projector,
        magnitude_formulation="PM1", pair_symmetrization="PS2",
    )
    interface = M10.rasterize_interface_field(response, shape)
    m10_vec = (
        np.asarray(interface["Rx_3d_interface"], dtype=np.float64),
        np.asarray(interface["Ry_3d_interface"], dtype=np.float64),
        np.asarray(interface["Rz_3d_interface"], dtype=np.float64),
    )

    slow_diff = concat_positive_bonds(us)
    fast_diff = concat_positive_bonds(uf)
    c_diff = concat_positive_bonds(c)

    coef_fast = A8.A8_INIT_DT * A8.A8_INIT_OMEGA * A8.A8_INIT_K
    coef_slow = A8.A8_INIT_DT * A8.A8_INIT_SLOW_TIMESCALE
    fast_contrib = coef_fast * fast_diff
    slow_contrib = coef_slow * slow_diff
    exact_pair_amp = fast_contrib + slow_contrib

    stored_pair_amp = np.concatenate((
        np.asarray(amps["A_xp"], dtype=np.float64)[:, :, :-1].ravel(),
        np.asarray(amps["A_yp"], dtype=np.float64)[:, :-1, :].ravel(),
        np.asarray(amps["A_zp"], dtype=np.float64)[:-1, :, :].ravel(),
    ))

    exact_relerr = rms(stored_pair_amp - exact_pair_amp) / max(rms(exact_pair_amp), 1.0e-30)
    pair_amp_rms = rms(stored_pair_amp)
    pair_response_rms = response_slot_vector_rms(response, shape)
    m10_rms = vector_rms(m10_vec)

    return {
        "c_state": c,
        "c_rms": rms(c),
        "c_bond_diff_rms": rms(c_diff),
        "slow_bond_diff_rms": rms(slow_diff),
        "fast_bond_diff_rms": rms(fast_diff),
        "coef_fast": float(coef_fast),
        "coef_slow": float(coef_slow),
        "fast_contrib_rms": rms(fast_contrib),
        "slow_contrib_rms": rms(slow_contrib),
        "pair_amp_rms": pair_amp_rms,
        "pair_amp_exact_formula_relative_rms_error": float(exact_relerr),
        "pair_response_slot_vector_rms": pair_response_rms,
        "m10_vector_rms": m10_rms,
        "pair_response_over_pair_amp": ratio(pair_response_rms, pair_amp_rms),
        "m10_over_pair_response": ratio(m10_rms, pair_response_rms),
        "pair_amp_over_c": ratio(pair_amp_rms, rms(c)),
        "m10_over_c": ratio(m10_rms, rms(c)),
        "c_bond_diff_over_c": ratio(rms(c_diff), rms(c)),
        "fast_contrib_over_c": ratio(rms(fast_contrib), rms(c)),
        "slow_contrib_over_c": ratio(rms(slow_contrib), rms(c)),
        "valid_longitudinal_count": int(np.count_nonzero(valid)),
        "gradient_rms": rms(gmag),
    }


def native_bond_decomposition(rho3: np.ndarray) -> dict:
    """Inspect local bonds of the frozen native accumulated equilibrium only."""
    build = LOCAL.native_accumulated_vector_zero_flux(rho3)
    c = np.asarray(build["c_state"], dtype=np.float64)
    u = np.asarray(build["accumulated"], dtype=np.float64)
    bonds = concat_positive_bonds(u)

    frac = np.abs(bonds) / FULL.EPSILON_MAX
    if np.any(frac >= 1.0):
        raise RuntimeError("native bond exceeded bounded-strain domain")
    traction = FULL.K0 * bonds / (1.0 - frac * frac)

    return {
        "c_rms": rms(c),
        "u_rms": rms(u),
        "bond_diff_rms": rms(bonds),
        "bond_traction_rms": rms(traction),
        "bond_diff_over_c": ratio(rms(bonds), rms(c)),
        "bond_traction_over_c": ratio(rms(traction), rms(c)),
        "traction_over_bond_diff": ratio(rms(traction), rms(bonds)),
        "max_abs_bond_fraction": float(np.max(frac)),
        "c_state_integral_relative_error": float(build["c_state_integral_relative_error"]),
        "accumulation_converged": bool(build["converged"]),
    }


def run_cluster(cluster: dict) -> dict:
    rho3 = local_rho3(cluster)
    unit = unit_decomposition(rho3)
    native = native_bond_decomposition(rho3)

    return {
        "cluster_id": cluster["id"],
        "source_rho3_rms": rms(rho3),
        "unit_c_rms": unit["c_rms"],
        "unit_c_bond_diff_rms": unit["c_bond_diff_rms"],
        "unit_c_bond_diff_over_c": unit["c_bond_diff_over_c"],
        "unit_fast_bond_diff_rms": unit["fast_bond_diff_rms"],
        "unit_slow_bond_diff_rms": unit["slow_bond_diff_rms"],
        "coef_fast": unit["coef_fast"],
        "coef_slow": unit["coef_slow"],
        "unit_fast_contrib_rms": unit["fast_contrib_rms"],
        "unit_slow_contrib_rms": unit["slow_contrib_rms"],
        "unit_fast_contrib_over_c": unit["fast_contrib_over_c"],
        "unit_slow_contrib_over_c": unit["slow_contrib_over_c"],
        "unit_pair_amp_rms": unit["pair_amp_rms"],
        "unit_pair_amp_over_c": unit["pair_amp_over_c"],
        "unit_pair_amp_exact_formula_relative_rms_error": unit["pair_amp_exact_formula_relative_rms_error"],
        "unit_pair_response_slot_vector_rms": unit["pair_response_slot_vector_rms"],
        "unit_pair_response_over_pair_amp": unit["pair_response_over_pair_amp"],
        "unit_m10_vector_rms": unit["m10_vector_rms"],
        "unit_m10_over_pair_response": unit["m10_over_pair_response"],
        "unit_m10_over_c": unit["m10_over_c"],
        "native_c_rms": native["c_rms"],
        "native_u_rms": native["u_rms"],
        "native_bond_diff_rms": native["bond_diff_rms"],
        "native_bond_diff_over_c": native["bond_diff_over_c"],
        "native_bond_traction_rms": native["bond_traction_rms"],
        "native_bond_traction_over_c": native["bond_traction_over_c"],
        "native_traction_over_bond_diff": native["traction_over_bond_diff"],
        "native_max_abs_bond_fraction": native["max_abs_bond_fraction"],
        "native_c_state_integral_relative_error": native["c_state_integral_relative_error"],
        "native_accumulation_converged": native["accumulation_converged"],
        "native_bond_diff_over_unit_m10": ratio(native["bond_diff_rms"], unit["m10_vector_rms"]),
        "native_bond_traction_over_unit_m10": ratio(native["bond_traction_rms"], unit["m10_vector_rms"]),
    }


def main() -> int:
    state = repo_state()
    clusters = list(BENCH.clusters())
    inventory = BENCH.inventory()
    ids = tuple(c["id"] for c in clusters)

    rows = []
    failures = []
    ready = bool(ids == EXPECTED_CLUSTER_IDS and len(inventory) == 5 and all(x["exists"] for x in inventory))
    if ready:
        for cluster in clusters:
            try:
                rows.append(run_cluster(cluster))
            except Exception as exc:
                failures.append({"cluster_id": cluster["id"], "error": f"{type(exc).__name__}: {exc}"})

    finite_keys = (
        "unit_c_rms", "unit_c_bond_diff_rms", "unit_pair_amp_rms",
        "unit_pair_response_slot_vector_rms", "unit_m10_vector_rms",
        "native_c_rms", "native_u_rms", "native_bond_diff_rms",
        "native_bond_traction_rms",
    )
    checks = {
        "canonical_five_cluster_inventory": ids == EXPECTED_CLUSTER_IDS,
        "all_five_local_fits_present": ready,
        "all_five_clusters_completed": bool(len(rows) == 5 and not failures),
        "all_measured_values_finite": bool(rows and all(all(math.isfinite(float(r[k])) for k in finite_keys) for r in rows)),
        "pair_amplitude_exact_formula_reproduced": bool(rows and all(r["unit_pair_amp_exact_formula_relative_rms_error"] <= 1.0e-12 for r in rows)),
        "native_c_state_integral_preserved": bool(rows and all(r["native_c_state_integral_relative_error"] <= 1.0e-12 for r in rows)),
        "native_accumulation_converged": bool(rows and all(r["native_accumulation_converged"] for r in rows)),
        "no_tracked_or_staged_changes": bool(not state["tracked_changes"] and not state["staged_changes"]),
    }
    execution_gate_pass = bool(all(checks.values()))

    if execution_gate_pass:
        status = "M10_LOCAL_INTERFACE_DECOMPOSITION_AUDIT_EXECUTED"
    elif rows:
        status = "M10_LOCAL_INTERFACE_DECOMPOSITION_AUDIT_PARTIAL_EXECUTION"
    else:
        status = "M10_LOCAL_INTERFACE_DECOMPOSITION_AUDIT_NOT_ESTABLISHED"

    result = {
        "lab_id": LAB_ID,
        "status": status,
        "repo_state": state,
        "model": {
            "unit_pair_amplitude_equation": "A_ij=(dt*omega*K)*Delta_u_fast+(dt*tau_slow)*Delta_u_slow",
            "unit_candidate": "PM1_PS2_M10",
            "native_bond_equation": "sigma=K0*Delta_u/(1-(Delta_u/epsilon_max)^2)",
            "observed_lensing_values_used": False,
            "network_access_used": False,
            "replacement_strength_scalar": None,
            "normalization_or_rescaling": False,
            "fit_or_tuning": False,
            "degrees_or_360_factor_used": False,
            "native_bonds_fed_to_G3D": False,
        },
        "rows": rows,
        "failures": failures,
        "checks": checks,
        "execution_gate_pass": execution_gate_pass,
        "interpretation_rule": "Determine whether the ~0.008 M10/c_state scale is generated primarily by frozen local T1 pair-amplitude coefficients, PM1/PS2 response construction, or M10 midpoint rasterisation. Native bond quantities are comparison diagnostics only.",
    }

    print(LAB_ID)
    print(f"status={status}")
    print(f"head_sha={state['head_sha']}")
    print("benchmark_loader=pbuf.core.benchmark_data")
    print("network_access_used=false")
    print("observed_lensing_values_used=false")
    print("degrees_or_360_factor_used=false")
    print("replacement_strength_scalar=none")
    print("native_response_rescaled=false")
    print("native_bonds_fed_to_G3D=false")
    print("fit_or_tuning=false")
    print(f"coef_fast={A8.A8_INIT_DT * A8.A8_INIT_OMEGA * A8.A8_INIT_K:.12g}")
    print(f"coef_slow={A8.A8_INIT_DT * A8.A8_INIT_SLOW_TIMESCALE:.12g}")
    print()
    print("CLUSTERS")
    for r in rows:
        print(
            f"cluster={r['cluster_id']} "
            f"c_bond_over_c={r['unit_c_bond_diff_over_c']:.12g} "
            f"pair_amp_over_c={r['unit_pair_amp_over_c']:.12g} "
            f"pair_response_over_amp={r['unit_pair_response_over_pair_amp']:.12g} "
            f"m10_over_pair_response={r['unit_m10_over_pair_response']:.12g} "
            f"m10_over_c={r['unit_m10_over_c']:.12g} "
            f"native_bond_over_c={r['native_bond_diff_over_c']:.12g} "
            f"native_traction_over_c={r['native_bond_traction_over_c']:.12g} "
            f"native_bond_over_m10={r['native_bond_diff_over_unit_m10']:.12g} "
            f"native_traction_over_m10={r['native_bond_traction_over_unit_m10']:.12g}"
        )
    for failure in failures:
        print(f"failure_cluster={failure['cluster_id']} error={failure['error']}")
    print()
    print("CHECKS")
    for key, value in checks.items():
        print(f"{key}={str(value).lower()}")
    print(f"execution_gate_pass={str(execution_gate_pass).lower()}")
    print(json.dumps(result, sort_keys=True, separators=(",", ":"), default=str))
    return 0 if execution_gate_pass else 1


if __name__ == "__main__":
    raise SystemExit(main())
