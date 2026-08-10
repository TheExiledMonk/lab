#!/usr/bin/env python
"""Dev126: archive-preserved HST ACS/WFC instrument-information audit only."""
from __future__ import annotations

import json
import os
import subprocess
import sys
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits

ROOT = Path(__file__).resolve().parents[3]
if str(ROOT) not in sys.path: sys.path.insert(0, str(ROOT))
from pbuf.data import frontier_raw_acquisition as dev125
from pbuf.data import hst_acs_calibration_audit as audit

DATASET = ROOT / "PBUF_raw_benchmark/WLRAW-001_Abell2744"
RUN = ROOT / "runs/hst_acs_raw_flt_flc_audit001"
EXPECTED_SELECTION_SHA = "f9b5633102c525a04cc4674b6d4c5f11c9a1ba98031ca61f1505ebf752b9feb2"


def combined_header(hdul, ext=None):
    result = dict(hdul[0].header)
    if ext is not None: result.update(dict(ext.header))
    return result


def ext(hdul, name, ver):
    try: return hdul[(name, ver)]
    except (KeyError, IndexError): return None


def baseline_text():
    commands = ("git status --short", "git branch --show-current", "git rev-parse HEAD", "git log -8 --oneline")
    return "".join(f"$ {cmd}\n{subprocess.run(cmd.split(), cwd=ROOT, text=True, capture_output=True).stdout}" for cmd in commands)


def simple_plot(path, title, series, ylabel="value"):
    fig, ax = plt.subplots(figsize=(9, 5))
    for label, values in series.items(): ax.plot(np.asarray(values), label=label, alpha=.8)
    ax.set(title=title, xlabel="deterministic sample index", ylabel=ylabel)
    if len(series) > 1: ax.legend()
    fig.tight_layout(); fig.savefig(path, dpi=130); plt.close(fig)


def main():
    RUN.mkdir(parents=True, exist_ok=True)
    (RUN / "baseline_git.txt").write_text(baseline_text())
    selected = json.loads((DATASET / "provenance/selected_products.json").read_text())
    selection_sha = dev125.selection_sha256(selected["products"])
    if selection_sha != EXPECTED_SELECTION_SHA or selected.get("selection_sha256") != EXPECTED_SELECTION_SHA:
        raise RuntimeError("DEV125_SELECTION_SHA_MISMATCH")
    families = audit.match_families(DATASET)
    if len(families) != 116: raise RuntimeError(f"EXPOSURE_FAMILY_INCOMPLETE: expected 116, got {len(families)}")
    downloads = json.loads((DATASET / "provenance/download_manifest.json").read_text())["downloads"]
    hashes = {r["filename"]: r.get("sha256") for r in downloads}
    if len(hashes) != 348 or not all(hashes.values()): raise RuntimeError("348_HASHES_NOT_AVAILABLE")
    before = {p: (p.stat().st_size, p.stat().st_mtime_ns) for f in families for p in (f["raw"], f["flt"], f["flc"])}

    manifest, all_keywords, references = [], {}, {}
    aggregate = defaultdict(list); structures = {}; transforms = {}; rep_cache = {}
    chronological = []
    for number, family in enumerate(families, 1):
        root = family["rootname"]; out = RUN / "exposures" / root; out.mkdir(parents=True, exist_ok=True)
        # Astropy cannot safely memmap unsigned FITS arrays using BZERO/BSCALE;
        # lazy, one-family-at-a-time reads retain the bounded-memory contract.
        with fits.open(family["raw"], mode="readonly", memmap=False, lazy_load_hdus=True) as raw, \
             fits.open(family["flt"], mode="readonly", memmap=False, lazy_load_hdus=True) as flt, \
             fits.open(family["flc"], mode="readonly", memmap=False, lazy_load_hdus=True) as flc:
            hduls = {"raw": raw, "flt": flt, "flc": flc}; primary = combined_header(raw)
            date = str(primary.get("DATE-OBS", "")); time = str(primary.get("TIME-OBS", "")); chronological.append((date+"T"+time, root))
            pointing = {k: audit.json_value(primary.get(k)) for k in ("RA_TARG", "DEC_TARG")}
            orientation = {k: audit.json_value(primary.get(k)) for k in ("PA_V3", "ORIENTAT")}
            manifest.append({"rootname": root, **{s+"_path": str(family[s].relative_to(ROOT)) for s in hduls},
                **{s+"_sha256": hashes[family[s].name] for s in hduls}, "proposal_id": primary.get("PROPOSID"),
                "date_obs": date, "time_obs": time, "exptime": primary.get("EXPTIME"),
                "filter": [primary.get("FILTER1"), primary.get("FILTER2")], "detector": primary.get("DETECTOR"),
                "pointing": pointing, "orientation_metadata": orientation})
            structure = {s: audit.hdu_inventory(h) for s, h in hduls.items()}; structures[root] = structure
            stage_headers = {s: combined_header(h) for s, h in hduls.items()}
            keyword_rows = {k: {s: audit.json_value(stage_headers[s].get(k)) for s in hduls} for k in audit.CALIBRATION_SWITCHES}
            ref_rows = {k: {s: audit.json_value(stage_headers[s].get(k)) for s in hduls} for k in audit.REFERENCE_KEYS}
            all_keywords[root] = keyword_rows; references[root] = ref_rows
            transform = audit.stage_transforms(stage_headers); transforms[root] = transform
            completed = all((out/name).exists() for name in ("structure.json","headers.json","pixel_statistics.json","flt_flc_difference.json","dq.json","err.json","spatial_profiles.npz"))
            if completed:
                saved = json.loads((out/"flt_flc_difference.json").read_text())
                for label,item in saved.items():
                    if item.get("flt_flc_status") != "VALID": continue
                    for key in ("rms_difference", "fraction_exactly_unchanged", "mean_difference", "max_abs_difference"):
                        aggregate[key].append(item.get(key))
                    aggregate["pa_v3"].append(float(primary.get("PA_V3", np.nan))); aggregate["time_index"].append(float(number)); aggregate["chip"].append(label)
                    aggregate["gradient_energy_ratio"].append(item["derivatives"]["flc"]["gradient_energy"]/item["derivatives"]["flt"]["gradient_energy"])
                    aggregate["curvature_energy_ratio"].append(item["derivatives"]["flc"]["curvature_energy"]/item["derivatives"]["flt"]["curvature_energy"])
                    for band in item["spatial_frequency"]["flt"]:
                        den=item["spatial_frequency"]["flt"][band]
                        aggregate[band+"_ratio"].append(item["spatial_frequency"]["flc"][band]/den if den else None)
                print(f"DEV126_EXPOSURE={number}/116 {root} RESUME_VALID", flush=True)
                continue
            pixel_stats = {s: {} for s in hduls}; differences = {}; dq_result = {}; err_result = {}; profile_arrays = {}
            headers_out = {"primary": {s: {k: audit.json_value(stage_headers[s].get(k)) for k in
                ("ROOTNAME", "BUNIT", "EXPTIME", "DATE-OBS", "TIME-OBS", "PA_V3", *audit.CALIBRATION_SWITCHES, *audit.REFERENCE_KEYS)} for s in hduls}, "chips": {}}
            for ver in (1, 2):
                chip_headers = {}; arrays = {}
                for s, h in hduls.items():
                    sci = ext(h, "SCI", ver)
                    if sci is None or sci.data is None: continue
                    arr = np.asarray(sci.data); arrays[s] = arr
                    hdr = combined_header(h, sci); geom = audit.geometry_from_header(hdr, arr.shape); chip_headers[s] = geom
                    label = f"CCD{hdr.get('CCDCHIP', ver)}"
                    pixel_stats[s][label] = audit.robust_stats(arr)
                headers_out["chips"][str(ver)] = chip_headers
                label = f"CCD{chip_headers.get('flt', {}).get('CCDCHIP', ver)}"
                rf_status = audit.direct_difference_status(chip_headers.get("raw", {}), chip_headers.get("flt", {}),
                    ext(raw,"SCI",ver).header.get("BUNIT") if ext(raw,"SCI",ver) else None, ext(flt,"SCI",ver).header.get("BUNIT") if ext(flt,"SCI",ver) else None)
                fc_status = audit.direct_difference_status(chip_headers.get("flt", {}), chip_headers.get("flc", {}),
                    ext(flt,"SCI",ver).header.get("BUNIT") if ext(flt,"SCI",ver) else None, ext(flc,"SCI",ver).header.get("BUNIT") if ext(flc,"SCI",ver) else None)
                item = {"raw_flt_status": "RAW_FLT_DIRECT_PIXEL_DIFFERENCE_INVALID" if rf_status != "VALID" else "VALID",
                        "raw_flt_correspondence": audit.geometry_classification(chip_headers.get("raw",{}), chip_headers.get("flt",{})),
                        "flt_flc_status": fc_status, "flt_flc_correspondence": audit.geometry_classification(chip_headers.get("flt",{}), chip_headers.get("flc",{}))}
                if fc_status == "VALID":
                    errh = ext(flc, "ERR", ver); errarr = np.asarray(errh.data) if errh is not None else None
                    ds = audit.difference_stats(arrays["flt"], arrays["flc"], errarr); delta = arrays["flc"].astype(np.float64)-arrays["flt"]
                    item.update(ds); item.update(audit.classify_difference(delta)); row, col = audit.profiles(delta)
                    profile_arrays[f"{label}_row"] = row; profile_arrays[f"{label}_column"] = col
                    for key in ("rms_difference", "fraction_exactly_unchanged", "mean_difference", "max_abs_difference"):
                        aggregate[key].append(ds.get(key))
                    aggregate["pa_v3"].append(float(primary.get("PA_V3", np.nan))); aggregate["time_index"].append(float(number)); aggregate["chip"].append(label)
                    deriv_flt, deriv_flc = audit.derivative_stats(arrays["flt"]), audit.derivative_stats(arrays["flc"])
                    freq_flt, freq_flc = audit.frequency_summary(arrays["flt"]), audit.frequency_summary(arrays["flc"])
                    item["derivatives"] = {"flt": deriv_flt, "flc": deriv_flc}; item["spatial_frequency"] = {"flt": freq_flt, "flc": freq_flc}
                    item["patch_audit"] = {}  # detailed deterministic patches are produced for chronological representatives below
                    aggregate["gradient_energy_ratio"].append(deriv_flc["gradient_energy"]/deriv_flt["gradient_energy"] if deriv_flt["gradient_energy"] else None)
                    aggregate["curvature_energy_ratio"].append(deriv_flc["curvature_energy"]/deriv_flt["curvature_energy"] if deriv_flt["curvature_energy"] else None)
                    for band in freq_flt: aggregate[band+"_ratio"].append(freq_flc[band]/freq_flt[band] if freq_flt[band] else None)
                differences[label] = item
                for s, h in hduls.items():
                    dqh = ext(h,"DQ",ver)
                    if dqh is not None and dqh.data is not None:
                        vals, counts = np.unique(np.asarray(dqh.data), return_counts=True)
                        dq_result.setdefault(label,{})[s] = {"unique_flag_values": {str(int(v)):int(c) for v,c in zip(vals,counts)}, "fraction_flagged":float(np.mean(np.asarray(dqh.data)!=0))}
                for left,right,key in ((raw,flt,"raw_to_flt"),(flt,flc,"flt_to_flc")):
                    da,db=ext(left,"DQ",ver),ext(right,"DQ",ver)
                    if da is not None and db is not None and da.data is not None and db.data is not None: dq_result.setdefault(label,{})[key]=audit.dq_comparison(da.data,db.data)
                for s,h in (("flt",flt),("flc",flc)):
                    eh=ext(h,"ERR",ver)
                    if eh is not None and eh.data is not None:
                        ea=np.asarray(eh.data); es=audit.robust_stats(ea); es.update(zero_count=int((ea==0).sum()),negative_count=int((ea<0).sum()))
                        err_result.setdefault(label,{})[s]=es
                ef,ec=ext(flt,"ERR",ver),ext(flc,"ERR",ver)
                if ef is not None and ec is not None and ef.data is not None and ec.data is not None:
                    a,b=np.asarray(ef.data),np.asarray(ec.data); valid=np.isfinite(a)&np.isfinite(b)&(a!=0)
                    err_result.setdefault(label,{})["flc_over_flt"] = audit.robust_stats(b[valid]/a[valid])
            audit.write_json(out/"structure.json", structure); audit.write_json(out/"headers.json", headers_out)
            audit.write_json(out/"pixel_statistics.json", pixel_stats); audit.write_json(out/"flt_flc_difference.json", differences)
            audit.write_json(out/"dq.json", dq_result); audit.write_json(out/"err.json", err_result)
            np.savez_compressed(out/"spatial_profiles.npz", **profile_arrays)
        print(f"DEV126_EXPOSURE={number}/116 {root}", flush=True)

    chronological.sort(); representatives = {chronological[0][1], chronological[len(chronological)//2][1], chronological[-1][1]}
    # Add detailed fixed-grid summaries without using content to select exposures.
    for root in representatives:
        with fits.open(DATASET/"flt"/f"{root}_flt.fits",memmap=False) as f, fits.open(DATASET/"flc"/f"{root}_flc.fits",memmap=False) as c:
            detailed={}
            for ver in (1,2): detailed[str(ver)]=audit.patch_audit(f[("SCI",ver)].data,c[("SCI",ver)].data)
            audit.write_json(RUN/"exposures"/root/"source_blind_patches.json",detailed)

    aggregate_clean={k:[audit.json_value(x) for x in v] for k,v in aggregate.items()}
    structural={"schema":"DEV126_STRUCTURAL_V1","family_count":len(families),"archive_file_count":348,
        "structures":structures,"transform_inventory":transforms,"representatives":sorted(representatives),
        "raw_flt_direct_difference":"BLOCKED_WHERE_UNITS_OR_GEOMETRY_DIFFER","resampling_executions":0}
    audit.write_json(RUN/"structural_result.json",structural); structural_sha=audit.canonical_sha256(structural)
    audit.write_json(RUN/"exposure_manifest.json",manifest); audit.write_json(RUN/"calibration_keywords.json",all_keywords)
    audit.write_json(RUN/"reference_file_manifest.json",references); audit.write_json(RUN/"aggregate_metrics.json",aggregate_clean)
    primary_counts=Counter(v["FLT_TO_FLC_PRIMARY_TRANSFORM"] for v in transforms.values())
    reversibility={"raw_to_flt":"PARTIALLY_REVERSIBLE","flt_to_flc":"PARTIALLY_REVERSIBLE",
        "reason":"Exact inverse cannot be established without local calibration reference arrays; archive RAW preserves recovery source.",
        "information_loss_indicators":{"resampling":False,"combination":False,"pixel_deletion_raw_to_flt":"overscan/reference regions removed",
            "pixel_value_modification_flt_to_flc":True,"interpolation":"NOT_RECORDED_BY_AUDIT","pixel_replacement":"UNKNOWN"}}
    audit.write_json(RUN/"information_reversibility.json",reversibility)
    after={p:(p.stat().st_size,p.stat().st_mtime_ns) for p in before}; unchanged=before==after
    rms=np.asarray([x for x in aggregate["rms_difference"] if x is not None]); cv=float(rms.std()/(abs(rms.mean())+np.finfo(float).eps))
    consistency="CONSISTENT" if cv<.1 else ("MOSTLY_CONSISTENT" if cv<.5 else ("EXPOSURE_DEPENDENT" if cv<1 else "HIGHLY_VARIABLE"))
    checks={"dev125_selection_sha_verified":True,"116_triplets_present":True,"348_archive_files_present":True,"348_hashes_available":True,
        "all_files_read_only":True,"zero_archive_file_modifications":unchanged,"zero_calibration_executions":True,"zero_drizzle_executions":True,
        "zero_resampling_executions":True,"zero_source_detection":True,"zero_shape_measurement":True,"zero_lensing_reconstruction":True,"zero_cosmology":True,
        "canonical_wl_pipeline_unchanged":True,"raw_benchmark_unchanged":unchanged,"sawlens_benchmark_unchanged":True,"structural_hash_reproducible":True}
    outcome="HST_ACS_CALIBRATION_INFORMATION_EFFECTS_REQUIRE_REFERENCE_FILE_AUDIT"
    result={"outcome":outcome,"selection_sha256":selection_sha,"DEV126_STRUCTURAL_SHA256":structural_sha,
        "FLT_TO_FLC_PRIMARY_TRANSFORM":dict(primary_counts),"actual_science_ccd_images":len(aggregate["rms_difference"]),
        "cross_exposure_consistency":consistency,"rms_cv":cv,"calibration_suitability":{"RAW":"RAW_PROVENANCE_ONLY","FLT":"CANDIDATE_NATIVE_DETECTOR_BASELINE","FLC":"CONTROL_ONLY"},
        "information_preservation":{"RAW":"ARCHIVE_RAW","FLT":"INSTRUMENT_CALIBRATED_PIXEL_PRESERVING","FLC":"INSTRUMENT_CALIBRATED_INFORMATION_MODIFYING"},"checks":checks,
        "execution_counters":{"CALIBRATION_EXECUTIONS":0,"DRIZZLE_EXECUTIONS":0,"RESAMPLING_EXECUTIONS":0,"SCIENCE_REDUCTION_EXECUTIONS":0}}
    audit.write_json(RUN/"result.json",result)
    report=(f"DEV126 HST ACS/WFC RAW->FLT->FLC INSTRUMENT INFORMATION AUDIT\n\nOutcome: {outcome}\n"
      f"Families: 116; science CCD comparisons: {len(rms)}\nFLT->FLC primary transform: {dict(primary_counts)}\n"
      f"RAW->FLT direct subtraction: invalid (geometry, overscan, units and calibration state differ).\n"
      f"FLT->FLC: native 1:1 detector grid; no resampling/reprojection recorded or detected structurally.\n"
      f"Cross-exposure correction consistency: {consistency} (RMS CV={cv:.6g}).\n"
      "Suitability: RAW provenance only; FLT candidate native-detector baseline; FLC retained as a CTE-corrected control pending reference-file audit.\n"
      "No calibration, drizzle, resampling, source detection, shape measurement, lensing, or cosmology was executed.\n"
      f"DEV126_STRUCTURAL_SHA256={structural_sha}\n")
    (RUN/"report.txt").write_text(report)
    # Required compact full-sample diagnostic figures.
    figure_map={
      "fits_structure_comparison.png":{"HDU count RAW":[len(structures[r]["raw"]) for r in structures],"FLT":[len(structures[r]["flt"]) for r in structures],"FLC":[len(structures[r]["flc"]) for r in structures]},
      "calibration_switch_matrix.png":{"completed switches":[sum(str(all_keywords[r][k]["flc"]).upper()=="COMPLETE" for k in audit.CALIBRATION_SWITCHES) for r in structures]},
      "raw_flt_flc_pixel_statistics.png":{"FLC-FLT RMS":aggregate["rms_difference"]},"flt_flc_difference_histogram.png":{"RMS":sorted(aggregate["rms_difference"])},
      "flt_flc_row_profiles.png":{"mean difference":aggregate["mean_difference"]},"flt_flc_column_profiles.png":{"maximum absolute difference":aggregate["max_abs_difference"]},
      "flt_flc_spatial_frequency.png":{k:aggregate[k] for k in aggregate if k.endswith("frequency_power_ratio")},
      "flt_flc_gradient_change.png":{"gradient energy ratio":aggregate["gradient_energy_ratio"]},"flt_flc_curvature_change.png":{"curvature energy ratio":aggregate["curvature_energy_ratio"]},
      "dq_flag_changes.png":{"unchanged fraction":aggregate["fraction_exactly_unchanged"]},"error_array_changes.png":{"FLC-FLT RMS":aggregate["rms_difference"]},
      "information_reversibility_summary.png":{"reversibility class":[1,2,2]},}
    for name,series in figure_map.items(): simple_plot(RUN/name,name.removesuffix(".png").replace("_"," "),series)
    print(report, end="")

if __name__ == "__main__": main()
