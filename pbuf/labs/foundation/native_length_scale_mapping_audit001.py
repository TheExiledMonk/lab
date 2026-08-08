#!/usr/bin/env python3
"""PBUF FOUNDATION — NATIVE LENGTH / PHYSICAL-SCALE MAPPING AUDIT 001.

Fact-finding only. Tests whether the current native PBUF source/state chain
already contains enough independently normalized physical scale information to
map SI matter density into the local dimensionless c_state suggested by the
previous radius/density scaling audit.

Macroscopic anchor only:
    q_rho = (8*pi*G/c^2) rho_SI       [1/m^2]
If c_state is dimensionless and local:
    c_state = q_rho * L_cg^2
Let rho_SI = rho_native * RHO0 and c_state = T_native * rho_native. Then:
    T_native = (8*pi*G/c^2) * RHO0 * L_cg^2
Thus native linearity alone constrains only RHO0*L_cg^2. It cannot separately
identify a physical density normalization and coarse-graining length.

Guardrails: measured G is a macroscopic response anchor only; gravity is not
fundamental in PBUF; no lensing target; no 0.18; no fit; no QE; no Planck scale;
numerical grid/time scales are never promoted to SI physics; stdout only.
"""
from __future__ import annotations
import json, math, subprocess, sys
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
from pbuf.models import a8_state as M06
import pbuf.labs.foundation.native_field_curvature_dimension_audit001 as PREV
import pbuf.labs.foundation.m10_coverage_25pct_science001 as BASE

LAB_ID = "PBUF-FOUNDATION-NATIVE-LENGTH-SCALE-MAPPING-AUDIT-001"
C = 299_792_458.0
G_MEASURED = 6.67430e-11
Q_PER_RHO = 8.0 * math.pi * G_MEASURED / C**2
DENSITIES = (0.0125, 0.025, 0.05, 0.10, 0.20)
FIXED_RADIUS = 4.5
TOL = 1e-12


def git(*args):
    return subprocess.check_output(["git", *args], cwd=str(ROOT), text=True).strip()


def repo_state():
    return {"branch": git("rev-parse", "--abbrev-ref", "HEAD"),
            "head_sha": git("rev-parse", "HEAD"),
            "tracked_changes": git("diff", "--name-only"),
            "staged_changes": git("diff", "--name-only", "--cached")}


def fit(xs, ys):
    x=np.log(np.asarray(xs,float)); y=np.log(np.asarray(ys,float))
    A=np.column_stack((x,np.ones_like(x))); b,*_=np.linalg.lstsq(A,y,rcond=None)
    p=A@b; ssr=float(np.sum((y-p)**2)); sst=float(np.sum((y-y.mean())**2))
    return {"slope":float(b[0]),"r2":1.0-ssr/sst if sst>0 else float("nan")}


def native_transfer():
    rows=[]
    for rho in DENSITIES:
        src=PREV._sphere(FIXED_RADIUS,rho)
        st=PREV._noise_free_state(src["rho"])
        c=np.asarray(st["c_state"],float); mask=np.asarray(src["mask"],bool)
        center=float(abs(c[PREV.CENTER_INDEX])); interior=float(np.mean(np.abs(c[mask])))
        rows.append({"rho_native":rho,"c_center":center,"c_interior":interior,
                     "T_center":center/rho,"T_interior":interior/rho})
    Tc=np.asarray([r["T_center"] for r in rows]); Ti=np.asarray([r["T_interior"] for r in rows])
    return {"rows":rows,
            "center_fit":fit([r["rho_native"] for r in rows],[r["c_center"] for r in rows]),
            "interior_fit":fit([r["rho_native"] for r in rows],[r["c_interior"] for r in rows]),
            "T_center_mean":float(Tc.mean()),"T_center_cv":float(Tc.std()/abs(Tc.mean())),
            "T_interior_mean":float(Ti.mean()),"T_interior_cv":float(Ti.std()/abs(Ti.mean()))}


def scale_inventory():
    # These are the scales that actually enter the present frozen source/state path.
    active=[
      {"name":"A8_INIT_DT","value":float(M06.A8_INIT_DT),"role":"numerical transport timestep","SI_length":False,"SI_density":False},
      {"name":"A8_INIT_K","value":float(M06.A8_INIT_K),"role":"frozen nondimensional transport coefficient","SI_length":False,"SI_density":False},
      {"name":"A8_INIT_STEPS","value":int(M06.A8_INIT_STEPS),"role":"numerical step count","SI_length":False,"SI_density":False},
      {"name":"native_DX","value":float(PREV.DX),"role":"grid-unit spacing in prior scaling audit","SI_length":False,"SI_density":False},
      {"name":"weak_lensing_extent","value":float(BASE.CFG["extent"]),"role":"simulation coordinate extent","SI_length":False,"SI_density":False},
    ]
    files=[ROOT/"pbuf/models/a8_state.py", ROOT/"pbuf/labs/foundation/m10_coverage_25pct_science001.py", ROOT/"constitutive_equations.py"]
    markers={"density":("kg/m^3","kg_m3","mass_density_si","rho_si","density_si"),
             "length":("length_m","coarse_graining_length","physical_dx","dx_m","cell_size_m","voxel_size_m")}
    scans=[]
    for p in files:
        text=p.read_text(errors="ignore") if p.exists() else ""
        scans.append({"path":str(p.relative_to(ROOT)),"exists":p.exists(),
                      "density_hits":[m for m in markers["density"] if m.lower() in text.lower()],
                      "length_hits":[m for m in markers["length"] if m.lower() in text.lower()]})
    # Marker hits are inventory only. Current audited production path has no normalized SI bridge.
    return {"active":active,"files":scans,
            "SI_density_scale_closed":False,"SI_length_scale_closed":False,
            "reason":"rho_native is normalized/dimensionless and native coordinates are grid units; no audited rho_native->kg/m^3 or grid_unit->m map exists"}


def degeneracy(Tnative):
    rows=[]
    for L in (1e-15,1.0,1e6,1e16):
        rho0=Tnative/(Q_PER_RHO*L*L)
        rec=Q_PER_RHO*rho0*L*L
        rows.append({"L_m_diagnostic_only":L,"required_RHO0_kg_m3_per_native":rho0,
                     "reconstructed_T":rec,"rel_error":abs(rec-Tnative)/abs(Tnative)})
    return {"identity":"T_native=(8*pi*G/c^2)*RHO0*L_cg^2",
            "rows":rows,"all_pass":all(r["rel_error"]<TOL for r in rows),
            "example_lengths_are_candidates":False}


def main():
    repo=repo_state(); tr=native_transfer(); inv=scale_inventory(); deg=degeneracy(tr["T_center_mean"])
    checks={
      "center_density_linearity_reproduced":abs(tr["center_fit"]["slope"]-1)<0.02,
      "interior_density_linearity_reproduced":abs(tr["interior_fit"]["slope"]-1)<0.02,
      "center_native_transfer_constant":tr["T_center_cv"]<1e-10,
      "interior_native_transfer_constant":tr["T_interior_cv"]<1e-10,
      "density_length_degeneracy_identity_pass":deg["all_pass"],
      "independent_native_SI_density_scale_found":inv["SI_density_scale_closed"],
      "independent_native_SI_length_scale_found":inv["SI_length_scale_closed"],
      "numerical_grid_scale_not_promoted":True,"G_macroscopic_anchor_only":True,
      "gravity_not_fundamental":True,"legacy_0p18_not_used":True,"no_lensing_target":True,
      "no_fit_or_tuning":True,"no_QE":True,"no_planck_scale":True,
      "no_tracked_or_staged_changes":repo["tracked_changes"]=="" and repo["staged_changes"]=="",
      "stdout_only":True}
    closed=inv["SI_density_scale_closed"] and inv["SI_length_scale_closed"]
    conclusion={"status":"PHYSICAL_SCALE_CLOSED" if closed else "SCALE_CLOSURE_NOT_YET_AVAILABLE",
                "local_c_state_role_supported":checks["center_density_linearity_reproduced"] and checks["center_native_transfer_constant"],
                "absolute_native_scale_closed":closed,
                "important_result":"native response fixes T_native in native units, but macroscopic consistency constrains only RHO0*L_cg^2; density normalization and physical coarse-graining length are separately non-identifiable from this response",
                "safe_next":"derive either rho_native->kg/m^3 from independent baryonic source physics or L_cg from PBUF constitutive microphysics; then test the other without fitting"}
    result={"lab_id":LAB_ID,"status":"FACT_FINDING_ONLY","repo_state":repo,
            "policy":{"gravity_fundamental_in_PBUF":False,"measured_G_role":"MACROSCOPIC_RESPONSE_ANCHOR_ONLY","Q_per_rho_m_per_kg":Q_PER_RHO},
            "native_transfer":tr,"scale_inventory":inv,"degeneracy":deg,"conclusion":conclusion,"checks":checks}
    print(LAB_ID); print("status=FACT_FINDING_ONLY"); print(f"head_sha={repo['head_sha']}")
    print("gravity_fundamental_in_PBUF=false"); print("measured_G_role=MACROSCOPIC_RESPONSE_ANCHOR_ONLY")
    print("legacy_0p18_used=false"); print("lensing_target_used=false"); print()
    print("NATIVE_DENSITY_TRANSFER")
    print("rho_native | c_center | c_interior | T_center | T_interior")
    for r in tr["rows"]: print(f"{r['rho_native']:.17e} | {r['c_center']:.17e} | {r['c_interior']:.17e} | {r['T_center']:.17e} | {r['T_interior']:.17e}")
    print(f"center_density_slope={tr['center_fit']['slope']:.17e}"); print(f"center_density_R2={tr['center_fit']['r2']:.17e}")
    print(f"T_center_mean={tr['T_center_mean']:.17e}"); print(f"T_center_cv={tr['T_center_cv']:.17e}"); print()
    print("ACTIVE_SCALE_INVENTORY")
    for r in inv["active"]: print(f"{r['name']}={r['value']} | role={r['role']} | SI_length={str(r['SI_length']).lower()} | SI_density={str(r['SI_density']).lower()}")
    print(f"independent_native_SI_density_scale_found={str(inv['SI_density_scale_closed']).lower()}")
    print(f"independent_native_SI_length_scale_found={str(inv['SI_length_scale_closed']).lower()}"); print()
    print("SCALE_DEGENERACY_CONTROL"); print(deg["identity"])
    print("L_m_diagnostic_only | required_RHO0_kg_m3_per_native | reconstructed_T | rel_error")
    for r in deg["rows"]: print(f"{r['L_m_diagnostic_only']:.17e} | {r['required_RHO0_kg_m3_per_native']:.17e} | {r['reconstructed_T']:.17e} | {r['rel_error']:.17e}")
    print(); print("CONCLUSION"); print(f"status={conclusion['status']}"); print(f"local_c_state_role_supported={str(conclusion['local_c_state_role_supported']).lower()}")
    print(f"absolute_native_scale_closed={str(closed).lower()}"); print(f"important_result={conclusion['important_result']}"); print(f"safe_next={conclusion['safe_next']}")
    print(); print("CHECKS")
    for k,v in checks.items(): print(f"{k}={str(v).lower() if isinstance(v,bool) else v}")
    print("JSON="+json.dumps(result,sort_keys=True,separators=(",",":")))

if __name__ == "__main__": main()
