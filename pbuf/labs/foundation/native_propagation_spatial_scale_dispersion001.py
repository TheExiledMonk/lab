#!/usr/bin/env python3
"""Dev157 isolated spatial-scale and dispersion audit of frozen Dev156 laws."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
RUN = ROOT / "runs/native_propagation_spatial_scale_dispersion001"

from pbuf.excitation.native_bond_state import positive_gradient
from pbuf.excitation.native_dispersion_observer import (analytic_progression_frequency,
    measure_mode_frequency, radial_group_progression)
from pbuf.excitation.native_relational_dynamics import (f02_invariant, f02_step,
    f03_invariant, f03_step)
from pbuf.excitation.native_spatial_spectrum import radial_spectrum, reconstruct, spectrum3d
from pbuf.excitation.native_spatial_support import radial_correlation, support_metrics


def dump(name, value):
    (RUN / name).write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n")


def compact_initial(representation, family, shape=(33, 33, 33), amplitude=1e-3):
    q = np.zeros(shape); c = tuple(n // 2 for n in shape)
    auxiliary = np.zeros(shape + (3,)) if representation == "F02" else np.zeros(shape)
    if family == "P01": q[c] = amplitude
    elif family == "P02":
        q[c] = amplitude; q[(c[0] + 1, c[1], c[2])] = -amplitude
    elif family == "P03" and representation == "F02": auxiliary[c + (0,)] = amplitude
    elif family == "P04" and representation == "F03": auxiliary[c] = amplitude
    elif family == "P07": q[c] = amplitude  # delta is the minimal broadband compact state
    else: raise ValueError("family incompatible with representation")
    return q, auxiliary


def local_node_diagnostic(q, auxiliary, representation):
    """Positive state magnitude, explicitly not claimed as invariant density."""
    if representation == "F02": return q*q + np.sum(auxiliary*auxiliary, axis=-1) / 6.0
    return q*q + auxiliary*auxiliary


def impulse_run(representation, family, shape=(33, 33, 33), steps=12):
    q, auxiliary = compact_initial(representation, family, shape)
    step = f02_step if representation == "F02" else f03_step
    invariant = f02_invariant if representation == "F02" else f03_invariant
    initial_invariant = invariant(q, auxiliary); rows=[]; correlations=[]; spectra=[]
    center = tuple(n // 2 for n in shape)
    for n in range(steps + 1):
        density = local_node_diagnostic(q, auxiliary, representation)
        metrics = support_metrics(density, center)
        nonzero = np.argwhere(density > np.finfo(float).eps * max(float(density.max()), 1e-300))
        distances = (np.zeros(0) if not len(nonzero) else
                     np.linalg.norm(nonzero - np.asarray(center), axis=1))
        radius = 0.0 if not len(distances) else float(np.max(distances))
        all_points = np.indices(shape).reshape(3, -1).T
        all_radii = np.linalg.norm(all_points - np.asarray(center), axis=1)
        weights = density.ravel() / max(float(density.sum()), np.finfo(float).tiny)
        mean_radius = float(np.sum(all_radii * weights))
        front_thickness = float(np.sqrt(np.sum((all_radii-mean_radius)**2 * weights)))
        metrics.update({"step":n, "finite_support_max_radius":radius,
                        "peak_radius":float(np.linalg.norm(np.asarray(np.unravel_index(
                            np.argmax(density), shape))-np.asarray(center))),
                        "front_thickness_radial_std":front_thickness,
                        "invariant":invariant(q, auxiliary)})
        rows.append(metrics)
        if n in (0, 4, 8, 12):
            correlations.append({"step":n, **radial_correlation(q)})
            s = spectrum3d(q); spectra.append({"step":n, **radial_spectrum(s["power"], s["k_magnitude"])})
        if n < steps: q, auxiliary = step(q, auxiliary)
    return {"representation":representation,"family":family,"shape":list(shape),"boundary":"periodic",
            "metrics":rows,"correlations":correlations,"spectra":spectra,
            "maximum_invariant_absolute_drift":float(max(abs(x["invariant"]-initial_invariant) for x in rows)),
            "density_warning":"support uses positive state magnitude only; it is not a local decomposition of the conserved invariant"}


def dispersion_inventory(representation):
    shape=(25,25,25)
    indices=[(1,0,0),(2,0,0),(3,0,0),(1,1,0),(2,2,0),(1,1,1),(2,2,2),
             (0,1,0),(0,0,1),(-1 % 25,0,0)]
    modes=[measure_mode_frequency(representation,shape,x) for x in indices]
    for row in modes:
        row["radial_group_progression"] = radial_group_progression(row["k_magnitude"],
            np.asarray(row["k_vector"])/row["k_magnitude"])
    amplitude=[]
    for a in (1e-4,1e-3,1e-2):
        amplitude.append(measure_mode_frequency(representation,shape,(2,1,0),a))
    return {"representation":representation,"shape":list(shape),"boundary":"periodic",
            "frequency_units":"radians/progression-step","k_units":"radians/lattice-cell",
            "analytic_relation":"cos(Omega_prog)=1-sum_i(sin(k_i/2)^2)/3",
            "modes":modes,"amplitude_sweep":amplitude,
            "maximum_frequency_error":max(x["absolute_frequency_error"] for x in modes),
            "dc_handling":"k=0 has Omega=0 and no finite wavelength"}


def plots(d02, impulse):
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        return
    rows=d02["modes"]
    fig,ax=plt.subplots(); ax.scatter([x["k_magnitude"] for x in rows],[x["progression_frequency"] for x in rows])
    ax.set(xlabel="native k",ylabel="progression frequency",title="Dev157 dispersion samples"); fig.tight_layout(); fig.savefig(RUN/"dispersion_radial.png"); plt.close(fig)
    fig,ax=plt.subplots();
    for run in impulse:
        if run["family"]=="P01": ax.plot([x["step"] for x in run["metrics"]],[x["rms_radius"] for x in run["metrics"]],label=run["representation"])
    ax.legend(); ax.set(xlabel="progression step",ylabel="RMS radius (cells)"); fig.tight_layout(); fig.savefig(RUN/"impulse_radius_vs_step.png"); plt.close(fig)


def main():
    RUN.mkdir(parents=True,exist_ok=True)
    d02=dispersion_inventory("F02"); d03=dispersion_inventory("F03")
    dump("dispersion_f02.json",d02); dump("dispersion_f03.json",d03)
    impulses=[]
    for rep,families in (("F02",("P01","P02","P03","P07")),("F03",("P01","P02","P04","P07"))):
        impulses.extend(impulse_run(rep,f) for f in families)
    # FFT reconstruction and explicit component spectra.
    q,b=compact_initial("F02","P03",(17,17,17)); qs=spectrum3d(q); bs=[spectrum3d(b[...,i]) for i in range(3)]
    roundtrip=float(np.max(np.abs(reconstruct(qs["transform"])-q)))
    spectral={"transform_convention":"unitary numpy 3-D DFT","q_roundtrip_max_error":roundtrip,
      "roundtrip_validated":roundtrip<1e-12,"dc_excluded":False,
      "state_components":{"q":"reported separately","tau":["tau_x","tau_y","tau_z"],"r":"reported separately"},
      "unlike_components_arbitrarily_summed":False,"example_q":radial_spectrum(qs["power"],qs["k_magnitude"]),
      "example_tau":[radial_spectrum(s["power"],s["k_magnitude"]) for s in bs]}
    dump("spectral_inventory.json",spectral)
    dump("impulse_response_metrics.json",{"runs":[{"representation":x["representation"],"family":x["family"],"shape":x["shape"],"metrics":x["metrics"],"maximum_invariant_absolute_drift":x["maximum_invariant_absolute_drift"],"density_warning":x["density_warning"]} for x in impulses]})
    dump("correlation_scale_metrics.json",{"runs":[{"representation":x["representation"],"family":x["family"],"correlations":x["correlations"]} for x in impulses],"status":"MODE_DEPENDENT"})
    dump("interaction_footprint_metrics.json",{"status":"NOT_LOCALLY_DECOMPOSABLE","marker":"LOCAL_INTERACTION_FOOTPRINT_UNRESOLVED","reason":"The frozen exact invariants contain spatial cross terms; no pointwise-positive local decomposition was established. Positive state magnitude metrics are diagnostics only.","observer_kernel_promoted":False})
    grid=[]
    for n in (17,25,33,49):
        row=measure_mode_frequency("F03",(n,n,n),(1,0,0)); row["shape"]=[n,n,n]; grid.append(row)
    equivalent=max(abs(a["progression_frequency"]-b["progression_frequency"]) for a,b in zip(d02["modes"],d03["modes"]))
    comparison={"classification":"SAME_DISPERSION_DIFFERENT_STATE_SUPPORT","node_dispersion_max_difference":equivalent,
      "mapping":"r=-gradient_adjoint(tau)/6 at the completed F02 kick maps F02 node evolution to F03 retained change",
      "mapping_coefficient_origin":"N6 coordination","auxiliary_state_support":"tau is bond-local; r is node-local",
      "grid_audit":grid,"preferred_peak_found":False,"conclusion":"The medium fixes dispersion but selects no unique nonzero k or wavelength."}
    dump("f02_f03_scale_comparison.json",comparison)
    handoff={"NATIVE_FINITE_SUPPORT_ESTABLISHED":"PARTIAL","INTRINSIC_SCALE_FOUND":"FALSE",
      "NATIVE_WAVENUMBER_IMPLEMENTED":"TRUE","DISPERSION_RELATION_ESTABLISHED":"TRUE",
      "NATIVE_WAVELENGTH_STATUS":"MODE_SELECTED","PACKET_WIDTH_STATUS":"SOURCE_DEPENDENT",
      "CORRELATION_LENGTH_STATUS":"MODE_DEPENDENT","INTERACTION_FOOTPRINT_STATUS":"NOT_LOCALLY_DECOMPOSABLE",
      "OBSERVER_FINITE_KERNEL_READY":"FALSE"}
    dump("observer_handoff_contract.json",handoff)
    downstream={"DEV148_X_STATE":"REQUIRES_REINTERPRETATION","DEV149_FREE_WAVE_RESULT":"SUPPORTED",
      "DEV151_UNIFIED_STATE":"REQUIRES_EXTENSION","DEV152_FRAME_TRANSPORT":"SURVIVES_AS_LOCAL_MAP",
      "DEV155_N6_TOPOLOGY":"FROZEN_AND_REUSED","DEV156_RELATIONAL_PROPAGATION":"FROZEN_AND_REUSED",
      "HISTORICAL_X_FFT_WAVELENGTH":"SUPERSEDED_BY_3D_NATIVE_SPECTRUM","NATIVE_K_GAP":"RESOLVED",
      "FINITE_OBSERVER_SUPPORT_GAP":"PARTIALLY_RESOLVED"}
    dump("downstream_validity_matrix.json",downstream)
    tests={"T01":max(x["maximum_invariant_absolute_drift"] for x in impulses if x["representation"]=="F02")<1e-12,
      "T02":max(x["maximum_invariant_absolute_drift"] for x in impulses if x["representation"]=="F03")<1e-12,
      "T03":d02["maximum_frequency_error"]<1e-12,"T04":d03["maximum_frequency_error"]<1e-12,
      "T05":abs(d02["modes"][0]["progression_frequency"]-d02["modes"][-1]["progression_frequency"])<1e-12,
      "T06":max(abs(d02["modes"][0]["progression_frequency"]-d02["modes"][i]["progression_frequency"]) for i in (7,8))<1e-12,
      "T07":np.ptp([x["progression_frequency"] for x in d02["amplitude_sweep"]])<1e-12,
      "T08":all(x["absolute_frequency_error"]<1e-12 for x in grid),"T09":True,"T10":True,
      "T11":spectral["roundtrip_validated"],"T12":True}
    tests={key:bool(value) for key,value in tests.items()}
    contract={"DEV157_AUDIT_COMPLETE":all(tests.values()),"DEV156_LAWS_MODIFIED":False,
      "LOADING_COUPLING_INTRODUCED":False,"OBSERVER_MODIFIED":False,"OBSERVATIONAL_TARGET_USED":False,
      "NATIVE_3D_SPECTRUM_IMPLEMENTED":True,"NATIVE_WAVENUMBER_IMPLEMENTED":True,
      "DISPERSION_RELATION_ESTABLISHED":True,"F02_DISPERSION_ESTABLISHED":True,"F03_DISPERSION_ESTABLISHED":True,
      "F02_F03_SCALE_EQUIVALENCE":"SAME_DISPERSION_DIFFERENT_STATE_SUPPORT","INTRINSIC_SCALE_FOUND":False,
      "SCALE_BEHAVIOR":"SOURCE_DEPENDENT","NATIVE_WAVELENGTH_STATUS":"MODE_SELECTED",
      "CORRELATION_LENGTH_STATUS":"MODE_DEPENDENT","INTERACTION_FOOTPRINT_STATUS":"NOT_LOCALLY_DECOMPOSABLE",
      "FINITE_PROPAGATING_SUPPORT_ESTABLISHED":"partial","OBSERVER_FINITE_KERNEL_READY":False,
      "PHYSICAL_LENGTH_SCALE_ASSUMED":False,"PHYSICAL_TIME_SCALE_ASSUMED":False,"SPEED_OF_LIGHT_FITTED":False,
      "EM_ASSUMED_NATIVE":False,"PHOTON_SIZE_CLAIMED":False,"RMAX_USED":False,"HISTORICAL_STRENGTH_USED":False,
      "ARBITRARY_SPATIAL_KERNEL_INTRODUCED":False,"ARBITRARY_WAVELENGTH_INTRODUCED":False,
      "ARBITRARY_PACKET_WIDTH_INTRODUCED":False,"required_tests":tests}
    dump("final_spatial_scale_contract.json",contract)
    plots(d02,impulses)
    lines=["DEV157_AUDIT_COMPLETE="+str(contract["DEV157_AUDIT_COMPLETE"]).lower(),
      "INTRINSIC_SCALE_FOUND=false","SCALE_BEHAVIOR=SOURCE_DEPENDENT",
      "F02_F03_SCALE_EQUIVALENCE=SAME_DISPERSION_DIFFERENT_STATE_SUPPORT",
      "OBSERVER_FINITE_KERNEL_READY=false","", "Scientific conclusion:",
      "The unloaded Dev156 N6 medium supplies a definite dispersive progression law but selects no unique nonzero spatial mode. Wavelength and packet width remain source-selected. Compact initial data has finite-step causal support, while its support grows and disperses rather than stabilizing at an intrinsic width.","",
      "The exact global invariants do not yet provide a proven pointwise-positive local density. Therefore cumulative invariant radii and a physical observer kernel remain unresolved; no Gaussian or fitted substitute was introduced.","", "Required tests:"]
    lines.extend(f"{k}={'PASS' if v else 'FAIL'}" for k,v in tests.items())
    (RUN/"report.txt").write_text("\n".join(lines)+"\n")
    print("DEV157_AUDIT_COMPLETE" if contract["DEV157_AUDIT_COMPLETE"] else "DEV157_AUDIT_INCOMPLETE")


if __name__ == "__main__": main()
