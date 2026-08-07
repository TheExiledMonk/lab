#!/usr/bin/env python3
"""PBUF FOUNDATION — ENDPOINT VS INTERFACE ROLE AUDIT 001.

Diagnostic-only audit of the representation used after the frozen
PL1_PM1_PS2 pair response. It holds the physical stochastic A8 initial
condition fixed under RC0..RC6, verifies the orientation parity of the
full pair response, then compares the reviewed M09 endpoint field and M10
interface field as candidate downstream physical representations.

No source modules are modified. No ray propagation/Jacobian/observational
fit is performed. LOS is used only after each 3D field has been inverse-
transformed to native RC0 coordinates, so the observer remains the native
z-axis in every comparison.
"""
from __future__ import annotations

import csv
import hashlib
import json
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
from astropy.io import fits

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))

from a8_three_dimensional_projection_lab001 import CLUSTERS, PRODUCTION, construct_common_proxy, construct_rho_3d
from pbuf.core import conventions as M01
from pbuf.core import coordinate_transforms as M02
from pbuf.core import vector_transforms as M03
from pbuf.core import pair_enumeration as M05
from pbuf.core import pair_transfer as M08
from pbuf.core import midpoint_rasterization as M10
from pbuf.core import los_projection as M14
from pbuf.models import a8_state as M06_state
from pbuf.models import a8_pair_amplitude as M06
from pbuf.models import transverse_projector as M07

LAB_ID = "PBUF-FOUNDATION-ENDPOINT-INTERFACE-ROLE-AUDIT-001"
OUT = ROOT / "runs" / "endpoint_interface_role_audit001"
BENCHMARK = ROOT / "PBUF_benchmark"
CLUSTER_ID = "MACS0416"
CANDIDATE_ID = "PL1_PM1_PS2"
NZ = 9
PROFILE = "gaussian"
STRENGTH = 0.18
SEED = 12345
PASS_TOL = 1e-8
CFG = dict(PRODUCTION)


def _git(*args):
    return subprocess.check_output(["git", *args], cwd=str(ROOT), text=True).strip()


def _json_default(obj):
    if isinstance(obj, np.floating): return float(obj)
    if isinstance(obj, np.integer): return int(obj)
    if isinstance(obj, np.bool_): return bool(obj)
    if isinstance(obj, tuple): return list(obj)
    return str(obj)


def _write_json(name, obj):
    (OUT / name).write_text(json.dumps(obj, indent=2, default=_json_default))


def _write_csv(name, rows):
    path = OUT / name
    if not rows:
        path.write_text("")
        return
    fields = []
    for row in rows:
        for key in row:
            if key not in fields: fields.append(key)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader(); writer.writerows(rows)


def _repo_state():
    return {
        "repository": "TheExiledMonk/lab",
        "branch": _git("rev-parse", "--abbrev-ref", "HEAD"),
        "head_sha": _git("rev-parse", "HEAD"),
        "tracked_changes": _git("diff", "--name-only"),
        "staged_changes": _git("diff", "--name-only", "--cached"),
    }


def _sha_file(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _norm_vector(v):
    x, y, z = (np.asarray(a, dtype=np.float64) for a in v)
    return float(np.sqrt(np.sum(x*x + y*y + z*z)))


def _err_vector(ref, test):
    rx, ry, rz = (np.asarray(a, dtype=np.float64) for a in ref)
    tx, ty, tz = (np.asarray(a, dtype=np.float64) for a in test)
    if rx.shape != tx.shape or ry.shape != ty.shape or rz.shape != tz.shape:
        return float("inf")
    num = float(np.sqrt(np.sum((tx-rx)**2 + (ty-ry)**2 + (tz-rz)**2)))
    return num / max(_norm_vector(ref), 1e-15)


def _err_pair2(ref1, ref2, tst1, tst2):
    r1, r2 = np.asarray(ref1), np.asarray(ref2)
    t1, t2 = np.asarray(tst1), np.asarray(tst2)
    if r1.shape != t1.shape or r2.shape != t2.shape: return float("inf")
    num = float(np.sqrt(np.sum((t1-r1)**2) + np.sum((t2-r2)**2)))
    den = float(np.sqrt(np.sum(r1*r1) + np.sum(r2*r2)))
    return num / max(den, 1e-15)


def _load_real_input():
    cluster = next(c for c in CLUSTERS if c["id"] == CLUSTER_ID)
    path = BENCHMARK / cluster["directory"] / f"hlsp_frontier_model_{cluster['slug']}_merten_v1_kappa.fits"
    if not path.exists(): raise FileNotFoundError(path)
    with fits.open(path) as hdul:
        kappa = np.asarray(hdul[0].data, dtype=np.float64)
    rho2 = construct_common_proxy(kappa, bins=CFG["bins"], extent=CFG["extent"])
    rho3 = construct_rho_3d(rho2, NZ, profile=PROFILE)
    return {"cluster": cluster, "fits_path": path, "fits_sha256": _sha_file(path), "rho3": rho3}


def _native_initial_state(rho3):
    rng = np.random.RandomState(SEED)
    eq = STRENGTH * rho3
    u_slow0 = eq.copy()
    noise = M06_state.A8_INIT_INJECTION_NOISE * STRENGTH * rng.randn(*rho3.shape)
    return {"rho_3d": rho3.copy(), "u_slow0": u_slow0, "u_fast0": eq + noise}


def _transform_initial(initial, rc):
    return {k: M02.transform_scalar_field(v, rc) for k, v in initial.items()}


def _evolve(initial):
    us, uf, history = M06_state.evolve_a8_transport_3d(
        initial["u_slow0"].copy(), initial["u_fast0"].copy(), stencil="N6", boundary="reflective"
    )
    return {"rho_3d": initial["rho_3d"].copy(), "u_slow": us, "u_fast": uf, "c_state": history[-1]}


def _candidate(state):
    shape = tuple(state["c_state"].shape)
    pairs = M05.enumerate_internal_pairs(shape)
    ex, ey, ez, valid, gmag = M07.build_longitudinal_direction(state["c_state"])
    P = M07.build_transverse_projector(ex, ey, ez)
    amp = M06.compute_a8_pair_amplitudes(state["u_slow"], state["u_fast"], state["c_state"], pairs)
    resp = M08.build_pair_responses(pairs, amp, P, magnitude_formulation="PM1", pair_symmetrization="PS2")
    end = M08.assemble_endpoint_field(resp, shape)
    iface = M10.rasterize_interface_field(resp, shape)
    return {"shape": shape, "pairs": pairs, "amp": amp, "resp": resp, "end": end, "iface": iface}


def _index_map(native_shape, rc):
    labels = np.arange(np.prod(native_shape), dtype=np.int64).reshape(native_shape)
    transformed = M02.transform_scalar_field(labels, rc)
    coords = np.indices(transformed.shape, dtype=np.int64).reshape(3, -1).T
    native_to_transformed = np.empty((labels.size, 3), dtype=np.int64)
    native_to_transformed[transformed.ravel().astype(np.int64)] = coords
    def map_index(idx):
        flat = np.ravel_multi_index(tuple(idx), native_shape)
        return tuple(int(x) for x in native_to_transformed[flat])
    return map_index


def _pair_lookup(pairs):
    return {frozenset((tuple(p.i_index), tuple(p.j_index))): p for p in pairs}


def _amp_at(amp, pair):
    return float(amp[f"A_{pair.axis}"][tuple(pair.i_index)])


def _response_at(resp, pair):
    i = tuple(pair.i_index)
    keys = {
        "xp": ("R_ij_xp", "R_ij_y_xp", "R_ij_z_xp"),
        "yp": ("R_ij_yp", "R_ij_y_yp", "R_ij_z_yp"),
        "zp": ("R_ij_zp", "R_ij_y_zp", "R_ij_z_zp"),
    }[pair.axis]
    return np.array([float(resp[k][i]) for k in keys], dtype=np.float64)


def _pair_parity(native, transformed, rc):
    map_index = _index_map(native["shape"], rc)
    lookup = _pair_lookup(transformed["pairs"])
    Q = np.asarray(M01.get_coordinate_matrix(rc, inverse=False), dtype=np.float64)
    a_ref, a_good, r_ref, r_even, r_odd = [], [], [], [], []
    reversed_count = missing = 0
    for p0 in native["pairs"]:
        it, jt = map_index(tuple(p0.i_index)), map_index(tuple(p0.j_index))
        pt = lookup.get(frozenset((it, jt)))
        if pt is None:
            missing += 1; continue
        same = tuple(pt.i_index) == it and tuple(pt.j_index) == jt
        sign = 1.0 if same else -1.0
        if sign < 0: reversed_count += 1
        a_ref.append(_amp_at(native["amp"], p0)); a_good.append(sign * _amp_at(transformed["amp"], pt))
        r0 = _response_at(native["resp"], p0); rb = Q.T @ _response_at(transformed["resp"], pt)
        r_ref.extend(r0.tolist()); r_even.extend(rb.tolist()); r_odd.extend((-rb).tolist())
    def rel(ref, test):
        ref, test = np.asarray(ref), np.asarray(test)
        return float(np.linalg.norm(test-ref) / max(np.linalg.norm(ref), 1e-15))
    return {
        "mapped_pair_count": len(a_ref), "missing_pair_count": missing,
        "orientation_reversing_pair_count": reversed_count,
        "E_pair_amplitude_oriented": rel(a_ref, a_good),
        "E_pair_response_even": rel(r_ref, r_even),
        "E_pair_response_odd": rel(r_ref, r_odd),
    }


def _endpoint_vec(cand):
    e = cand["end"]; return e["Rx_3d"], e["Ry_3d"], e["Rz_3d"]


def _interface_vec(cand):
    i = cand["iface"]; return i["Rx_3d_interface"], i["Ry_3d_interface"], i["Rz_3d_interface"]


def _back_vector(v, rc):
    return v if rc == "RC0" else M03.inverse_transform_vector_field(*v, rc)


def _los_native_z(v):
    out = M14.project_vector_to_image_plane(*v, los_axis="z")
    return out["comp_1"], out["comp_2"]


def main():
    started = time.perf_counter(); OUT.mkdir(parents=True, exist_ok=True)
    repo = _repo_state(); _write_json("repository_state.json", repo)
    if repo["branch"] != "main" or repo["tracked_changes"] or repo["staged_changes"]:
        v = {"lab_id": LAB_ID, "outcome": "REPOSITORY_GATE_FAILURE", "head_sha": repo["head_sha"]}
        _write_json("validation.json", v); print(json.dumps(v, indent=2)); return 2

    real = _load_real_input(); initial0 = _native_initial_state(real["rho3"])
    cand0 = _candidate(_evolve(initial0)); ref_end = _endpoint_vec(cand0); ref_iface = _interface_vec(cand0)
    ref_los_end, ref_los_iface = _los_native_z(ref_end), _los_native_z(ref_iface)
    _write_json("input_provenance.json", {
        "lab_id": LAB_ID, "cluster_id": CLUSTER_ID, "candidate_id": CANDIDATE_ID,
        "fits_path": str(real["fits_path"].relative_to(ROOT)), "fits_sha256": real["fits_sha256"],
        "native_shape": list(cand0["shape"]), "nz": NZ, "profile": PROFILE, "strength": STRENGTH,
        "seed": SEED, "observer_rule": "inverse-transform 3D field to RC0, then project native z LOS",
    })

    parity_rows=[]; role_rows=[]; endpoint_rows=[]; interface_rows=[]; los_rows=[]; wrong_rows=[]
    for rc in M01.RC_TRANSFORMS:
        cand_t = _candidate(_evolve(_transform_initial(initial0, rc)))
        parity = _pair_parity(cand0, cand_t, rc)
        end_back = _back_vector(_endpoint_vec(cand_t), rc); iface_back = _back_vector(_interface_vec(cand_t), rc)
        E_end, E_iface = _err_vector(ref_end, end_back), _err_vector(ref_iface, iface_back)
        los_end, los_iface = _los_native_z(end_back), _los_native_z(iface_back)
        E_los_end = _err_pair2(*ref_los_end, *los_end); E_los_iface = _err_pair2(*ref_los_iface, *los_iface)
        estats, istats = cand_t["end"]["statistics"], cand_t["iface"]["statistics"]
        expected = int(M08.expected_interface_pair_count(cand_t["shape"])); consumed = int(istats["consumed_pair_count_total"])
        pair_ok = expected == consumed
        parity_rows.append({"transform": rc, **parity})
        endpoint_rows.append({"transform": rc, "E_endpoint": E_end, "endpoint_energy": float(estats["endpoint_energy"]), "endpoint_global_sum_norm": float(estats["global_vector_sum_norm"]), "orientation_reversing_pair_count": parity["orientation_reversing_pair_count"]})
        interface_rows.append({"transform": rc, "E_interface": E_iface, "interface_energy": float(istats["interface_energy"]), "expected_pair_count": expected, "consumed_pair_count": consumed, "pair_count_ok": pair_ok, "orientation_reversing_pair_count": parity["orientation_reversing_pair_count"]})
        los_rows.append({"transform": rc, "E_LOS_endpoint_after_native_backtransform": E_los_end, "E_LOS_interface_after_native_backtransform": E_los_iface})
        role_rows.append({"transform": rc, "reversing_pair_count": parity["orientation_reversing_pair_count"], "E_pair_amplitude_oriented": parity["E_pair_amplitude_oriented"], "E_pair_response_even": parity["E_pair_response_even"], "E_pair_response_odd": parity["E_pair_response_odd"], "E_endpoint": E_end, "E_interface": E_iface, "E_LOS_endpoint": E_los_end, "E_LOS_interface": E_los_iface, "endpoint_closure": float(estats["global_vector_sum_norm"]), "endpoint_energy": float(estats["endpoint_energy"]), "interface_energy": float(istats["interface_energy"]), "interface_pair_count_ok": pair_ok})
        wrong_rows.append({"transform": rc, "WC_response_odd_parity": parity["E_pair_response_odd"], "WC_endpoint_production_orientation_dependence": E_end, "WC_interface_endpoint_swap_invariance_proxy": E_iface})

    _write_csv("pair_response_parity.csv", parity_rows); _write_csv("representation_role_audit.csv", role_rows)
    _write_csv("endpoint_covariance.csv", endpoint_rows); _write_csv("interface_covariance.csv", interface_rows)
    _write_csv("los_covariance.csv", los_rows); _write_csv("wrong_controls.csv", wrong_rows)

    pair_even_pass = all(r["E_pair_response_even"] <= PASS_TOL for r in parity_rows)
    endpoint_pass = all(r["E_endpoint"] <= PASS_TOL for r in role_rows)
    interface_pass = all(r["E_interface"] <= PASS_TOL for r in role_rows)
    los_interface_pass = all(r["E_LOS_interface"] <= PASS_TOL for r in role_rows)
    pair_counts_pass = all(bool(r["interface_pair_count_ok"]) for r in role_rows)
    reversed_rows = [r for r in role_rows if r["reversing_pair_count"] > 0]
    nonreversed_rows = [r for r in role_rows if r["reversing_pair_count"] == 0]
    endpoint_orientation_pattern = bool(reversed_rows and all(r["E_endpoint"] > PASS_TOL for r in reversed_rows) and all(r["E_endpoint"] <= PASS_TOL for r in nonreversed_rows))

    if not pair_even_pass: outcome = "Outcome D — PAIR RESPONSE REMAINS NONCOVARIANT"
    elif endpoint_pass and interface_pass: outcome = "Outcome B — BOTH REPRESENTATIONS COVARIANT"
    elif interface_pass and los_interface_pass and pair_counts_pass and endpoint_orientation_pattern: outcome = "Outcome A — INTERFACE FIELD IS COORDINATE-SAFE; ENDPOINT FIELD IS ORIENTATION-DEPENDENT"
    elif not interface_pass: outcome = "Outcome C — DOWNSTREAM PAIR ASSEMBLY STILL NONCOVARIANT"
    elif interface_pass and not los_interface_pass: outcome = "Outcome E — LOS PROJECTION ROLE/COVARIANCE FAILURE"
    else: outcome = "Outcome F — ROLE AUDIT INCONCLUSIVE"

    endpoint_bookkeeping = endpoint_orientation_pattern and all(r["endpoint_energy"] > 0 for r in role_rows)
    interface_physical = pair_even_pass and interface_pass and los_interface_pass and pair_counts_pass
    validation = {
        "lab_id": LAB_ID, "outcome": outcome, "head_sha": repo["head_sha"],
        "pair_response_even_covariant": pair_even_pass,
        "endpoint_representation_covariant": endpoint_pass,
        "endpoint_representation_orientation_dependent": endpoint_orientation_pattern,
        "interface_representation_coordinate_safe": interface_pass,
        "interface_native_z_projection_coordinate_safe": los_interface_pass,
        "interface_pair_counts_pass": pair_counts_pass,
        "endpoint_role_candidate": "conservation_bookkeeping" if endpoint_bookkeeping else None,
        "interface_role_candidate": "physical_pair_field" if interface_physical else None,
        "downstream_source_change_authorized": False, "ray_rerun_authorized": False,
        "duration_seconds": time.perf_counter()-started,
    }
    _write_json("validation.json", validation); _write_json("run.json", {"lab_id": LAB_ID, "head_sha": repo["head_sha"], "duration_seconds": validation["duration_seconds"]})

    lines=[f"# {LAB_ID}","",f"**Head:** `{repo['head_sha']}`","",f"**Outcome:** {outcome}","","| RC | reversed pairs | E response even | E response odd | E endpoint | E interface | E LOS endpoint | E LOS interface | endpoint closure | interface pairs ok |","|---|---:|---:|---:|---:|---:|---:|---:|---:|---|"]
    for r in role_rows:
        lines.append(f"| {r['transform']} | {r['reversing_pair_count']} | {r['E_pair_response_even']:.3e} | {r['E_pair_response_odd']:.3e} | {r['E_endpoint']:.3e} | {r['E_interface']:.3e} | {r['E_LOS_endpoint']:.3e} | {r['E_LOS_interface']:.3e} | {r['endpoint_closure']:.3e} | {r['interface_pair_count_ok']} |")
    lines += ["","## Role classification","",f"- pair response orientation-even covariance: `{pair_even_pass}`",f"- endpoint orientation dependence pattern: `{endpoint_orientation_pattern}`",f"- interface 3D covariance: `{interface_pass}`",f"- interface native-z LOS covariance after back-transform: `{los_interface_pass}`",f"- endpoint role candidate: `{validation['endpoint_role_candidate']}`",f"- interface role candidate: `{validation['interface_role_candidate']}`","","No source change or ray rerun is authorized by this lab."]
    (OUT / "report.md").write_text("\n".join(lines)+"\n")
    print(json.dumps(validation, indent=2)); return 0


if __name__ == "__main__":
    raise SystemExit(main())
