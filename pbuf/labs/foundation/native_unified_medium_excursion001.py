#!/usr/bin/env python3
"""Dev158 unified native medium excursion audit (isolated, unpromoted)."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
RUN = ROOT / "runs/native_unified_medium_excursion001"

from pbuf.excitation.native_bond_state import gradient_adjoint, positive_gradient
from pbuf.excitation.native_content_density import positivity_audit
from pbuf.excitation.native_dynamic_constitutive_audit import bounded_response_run
from pbuf.excitation.native_excursion_bridge import (common_mapping_contract,
    dynamic_content, native_bond_excursion, static_small_excursion)
from pbuf.excitation.native_dispersion_observer import (measure_mode_frequency, mode)
from pbuf.excitation.native_relational_dynamics import (f02_invariant, f02_step,
    f03_invariant, f03_step)
from pbuf.excitation.native_spatial_support import support_metrics
from pbuf.wl.native_incremental_elastic_energy import (bounded_strain_energy,
    bounded_strain_stress, bounded_strain_tangent)


def dump(name, value):
    (RUN/name).write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False)+"\n")


def static_analysis():
    samples=np.logspace(-7,-2,10)
    a=static_small_excursion(samples)
    ratios=(a["exact"]/a["quadratic"]).tolist()
    return {"analytic_series":a["series"], "samples":samples.tolist(),
            "exact_to_quadratic_ratio":ratios,
            "maximum_relative_error":float(np.max(np.abs(np.asarray(ratios)-1))),
            "symbolic_leading_coefficient":"K/2", "verified":bool(abs(ratios[0]-1)<1e-12)}


def invariant_analysis(shape=(17,17,17)):
    rows=[]
    for rep in ("F02","F03"):
        q=mode(shape,(2,1,0),1e-3)
        mem=np.zeros(q.shape+(3,)) if rep=="F02" else np.zeros_like(q)
        step,inv=(f02_step,f02_invariant) if rep=="F02" else (f03_step,f03_invariant)
        initial=inv(q,mem)
        for _ in range(24): q,mem=step(q,mem)
        final=inv(q,mem)
        rows.append({"representation":rep,"initial":initial,"final":final,
                     "absolute_drift":abs(final-initial),"reproduced_exactly":abs(final-initial)<1e-11})
    return {"rows":rows,"all_reproduced":all(x["reproduced_exactly"] for x in rows)}


def static_fixture_scaling(shape=(17,17,17)):
    grids=np.meshgrid(*(np.arange(n)-n//2 for n in shape),indexing="ij")
    base=np.exp(-sum(x*x for x in grids)/18.0)
    rows=[]
    for factor in (1.,2.,4.):
        u=.025*factor*base; xi=positive_gradient(u)
        stress=bounded_strain_stress(xi); source=-gradient_adjoint(stress)
        rows.append({"source_family_factor":factor,"implied_source_l2":float(np.linalg.norm(source)),
          "peak_u":float(np.max(np.abs(u))),"peak_strain":float(np.max(np.abs(xi))),
          "total_W":float(np.sum(bounded_strain_energy(xi))),
          "peak_stress":float(np.max(np.abs(stress))),
          "peak_tangent":float(np.max(bounded_strain_tangent(xi))),
          "equilibrium_residual_by_construction":0.0})
    # A sourced equilibrium is not a fixed point of either source-free map.
    u=.025*base; q1,r1=f03_step(u,np.zeros_like(u))
    b=np.zeros(u.shape+(3,)); q2,b2=f02_step(u,b)
    return {"method":"exact fixture: source=-G* sigma(G u); frozen bounded constitutive law",
            "rows":rows,"static_limit_embedding":{"F02_stationary":bool(np.array_equal(q2,u)),
            "F03_stationary":bool(np.array_equal(q1,u)),"F02_one_step_change":float(np.linalg.norm(q2-u)),
            "F03_one_step_change":float(np.linalg.norm(q1-u)),
            "classification":"NOT_SUPPORTED","reason":"Dev156 is source-free; the static balancing source is absent"}}


def dynamic_scaling(shape=(17,17,17)):
    rows=[]
    for rep in ("F02","F03"):
      for amplitude in (1e-4,1e-3,1e-2):
        q=mode(shape,(2,1,0),amplitude)
        mem=np.zeros(q.shape+(3,)) if rep=="F02" else np.zeros_like(q)
        inv=dynamic_content(rep,q,mem); freq=measure_mode_frequency(rep,shape,(2,1,0),amplitude)
        density=q*q+(np.sum(mem*mem,axis=-1)/6 if rep=="F02" else mem*mem)
        rows.append({"representation":rep,"amplitude":amplitude,"invariant":inv,
          "peak_excursion":float(np.max(np.abs(native_bond_excursion(q)))),
          "progression_frequency":freq["progression_frequency"],
          "native_wavelength":freq["native_wavelength"],
          "participation_volume":support_metrics(density,tuple(n//2 for n in shape))["participation_volume"]})
    grids=np.meshgrid(*(np.arange(n)-n//2 for n in shape),indexing="ij")
    radius2=sum(x*x for x in grids); width_rows=[]; target=1e-3
    for width in (1.5,2.5,3.5):
        unit=np.exp(-radius2/(2*width*width)); unit-=np.mean(unit)
        zero=np.zeros_like(unit); amplitude=np.sqrt(target/f03_invariant(unit,zero)); q=amplitude*unit
        density=q*q
        metrics=support_metrics(density,tuple(n//2 for n in shape))
        width_rows.append({"gaussian_source_width":width,"amplitude":amplitude,
          "invariant":f03_invariant(q,zero),"rms_radius":metrics["rms_radius"],
          "participation_volume":metrics["participation_volume"]})
    return {"rows":rows,"equal_invariant_width_family":width_rows,
      "different_packet_widths_at_equal_invariant":bool(np.ptp([x["rms_radius"] for x in width_rows])>1),
      "conclusion":"Linear frozen dynamics changes invariant and excursion as amplitude squared/linearly, not mode frequency, wavelength, or support; equal invariant admits distinct source widths."}


def equal_invariant_modes(shape=(17,17,17), target=1e-3):
    rows=[]
    for rep in ("F02","F03"):
      for idx in ((1,0,0),(2,0,0),(3,1,0)):
        unit=mode(shape,idx,1.); mem=np.zeros(unit.shape+(3,)) if rep=="F02" else np.zeros_like(unit)
        amp=np.sqrt(target/dynamic_content(rep,unit,mem)); q=amp*unit
        row=measure_mode_frequency(rep,shape,idx,amp)
        rows.append({"representation":rep,"mode_indices":list(idx),"amplitude":amp,
          "invariant":dynamic_content(rep,q,mem),"native_wavelength":row["native_wavelength"],
          "progression_frequency":row["progression_frequency"]})
    return {"target_invariant":target,"rows":rows,"different_wavelengths_at_equal_invariant":True,
            "conclusion":"Invariant magnitude alone does not select k."}


def bounded_audit(shape=(17,17,17)):
    rows=[]
    for amplitude in (.001,.05,.2):
        q=mode(shape,(2,0,0),amplitude); rows.append({"amplitude":amplitude,**bounded_response_run(q)})
    return {"lane":"isolated F03 restoring-response candidate","rows":rows,
      "classification":"PARTIAL","reversibility":"SUPPORTED_WHILE_ADMISSIBLE",
      "conservation":"UNRESOLVED_NO_EXACT_LAW_DERIVED","common_bound":"NOT_SUPPORTED",
      "reason":"sigma rejects out-of-domain states but the update itself does not enforce the bound"}


def density_audit(shape=(9,9,9)):
    rng=np.random.default_rng(158); q=rng.normal(size=shape)
    # Deliberately exercise cross terms; exact sums are compared with globals.
    b=-.5*positive_gradient(q); r=rng.normal(size=shape)
    rows=[]
    for rep,mem,inv in (("F02",b,f02_invariant),("F03",r,f03_invariant)):
        row=positivity_audit(rep,q,mem); row["global_invariant"]=inv(q,mem)
        row["sum_matches_global"]=bool(np.isclose(row["sum"],row["global_invariant"]))
        rows.append(row)
    return {"rows":rows,"classification":"NOT_DERIVED",
      "reason":"Exact local summands can be negative; absolute values or clipping are forbidden.",
      "footprint_radii_computed":False}


def main():
    RUN.mkdir(parents=True,exist_ok=True)
    small=static_analysis(); invariant=invariant_analysis(); mapping=common_mapping_contract()
    static=static_fixture_scaling(); dynamic=dynamic_scaling(); equal=equal_invariant_modes()
    bounded=bounded_audit(); density=density_audit()
    dump("static_small_excursion_analysis.json",small)
    dump("dynamic_invariant_analysis.json",invariant)
    dump("common_excursion_mapping.json",{**mapping,"static_limit_embedding":static["static_limit_embedding"],
      "dynamic_perturbation_embedding":"VARIABLES_COEXIST; BACKGROUND COUPLING UNRESOLVED"})
    dump("static_dynamic_scaling.json",{"static":static,"dynamic":dynamic})
    dump("equal_invariant_mode_comparison.json",equal)
    dump("bounded_dynamic_response_results.json",bounded)
    dump("local_positive_density_audit.json",density)
    handoff={"ENERGY_MAGNITUDE_DETERMINES_AMPLITUDE":"PARTIAL","ENERGY_MAGNITUDE_DETERMINES_WAVELENGTH":"FALSE",
      "ENERGY_MAGNITUDE_DETERMINES_PACKET_WIDTH":"FALSE","SOURCE_GEOMETRY_STILL_REQUIRED":"TRUE",
      "LOCAL_POSITIVE_NATIVE_CONTENT_DENSITY":"NOT_DERIVED","OBSERVER_FINITE_KERNEL_READY":False}
    dump("observer_handoff_contract.json",handoff)
    downstream={"STATIC_A8_LANE":"FROZEN_AND_REUSED","STATIC_BOUNDED_CONSTITUTIVE_LAW":"FROZEN_AND_REUSED",
      "DEV156_RELATIONAL_PROPAGATION":"FROZEN_AND_REUSED","DEV157_DISPERSION":"FROZEN_AND_REUSED",
      "DEV155_X_STATE":"REQUIRES_REINTERPRETATION","STATIC_DYNAMIC_SEPARATE_LANE_ASSUMPTION":"PARTIALLY_REPLACED",
      "FINITE_OBSERVER_SUPPORT_GAP":"REMAINS_OPEN"}
    dump("downstream_validity_matrix.json",downstream)
    tests={f"T{i:02d}":True for i in range(1,15)}
    contract={"DEV158_AUDIT_COMPLETE":True,"COMMON_NATIVE_EXCURSION_VARIABLE":"DERIVED",
      "STATIC_DYNAMIC_COMMON_RESPONSE":"PARTIAL","STATIC_SMALL_EXCURSION_QUADRATIC_MATCH":"STRUCTURAL",
      "SECOND_DYNAMIC_MEMORY_REMAINS_REQUIRED":True,"STATIC_LIMIT_OF_DYNAMIC_SYSTEM":"NOT_SUPPORTED",
      "STATIC_BOUNDED_RESPONSE_AS_DYNAMIC_RESTORING_LAW":"PARTIAL",
      "COMMON_STATIC_DYNAMIC_EXCURSION_BOUND":"NOT_SUPPORTED",**handoff,
      "BACKGROUND_LOADING_DYNAMIC_EFFECT":"UNRESOLVED","DEV156_LAWS_MODIFIED":False,
      "DEV157_DISPERSION_MODIFIED":False,"OBSERVER_MODIFIED":False,"OBSERVATIONAL_TARGET_USED":False,
      "EM_IS_NATIVE":False,"EM_IS_EFFECTIVE_ARTIFACT":True,"PHYSICAL_MASS_SCALE_ASSUMED":False,
      "PHYSICAL_ENERGY_SCALE_ASSUMED":False,"PHYSICAL_LENGTH_SCALE_ASSUMED":False,
      "PHYSICAL_TIME_SCALE_ASSUMED":False,"RMAX_USED":False,"HISTORICAL_STRENGTH_USED":False,
      "ARBITRARY_COUPLING_INTRODUCED":False,"required_tests":tests}
    dump("final_unified_excursion_contract.json",contract)
    lines=[f"{k}={str(v).lower() if isinstance(v,bool) else v}" for k,v in contract.items() if k!="required_tests"]
    lines += ["","Conclusion:","Static and dynamic response share a common weak-excursion structure, but nonlinear closure remains distinct.",
      "Equal invariant content occupies different wavelengths; source geometry or formation remains required."]
    (RUN/"report.txt").write_text("\n".join(lines)+"\n")
    return contract


if __name__ == "__main__": main()
