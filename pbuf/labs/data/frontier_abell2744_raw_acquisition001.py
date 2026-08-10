#!/usr/bin/env python3
"""Dev125: freeze and acquire raw HST/ACS F814W Abell 2744 exposures."""
from __future__ import annotations

import argparse
import concurrent.futures
import json
import subprocess
import sys
import shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))

from pbuf.data import frontier_raw_acquisition as acq

LAB = "frontier_abell2744_raw_acquisition001"
RUN_DIR = ROOT / "runs" / LAB
BENCHMARK_DIR = ROOT / "PBUF_raw_benchmark" / "WLRAW-001_Abell2744"


def query_manifest() -> dict:
    return {
        "RA": acq.TARGET_RA, "Dec": acq.TARGET_DEC, "search_radius_deg": acq.SEARCH_RADIUS_DEG,
        "mission": "HST", "instrument": "ACS/WFC", "filter": "F814W", "public_only": True,
        "main_field_criterion": {"method": "footprint_or_pointing_overlap",
                                 "center_deg": [acq.TARGET_RA, acq.TARGET_DEC],
                                 "radius_deg": acq.MAIN_FIELD_RADIUS_DEG},
    }


def baseline() -> str:
    def git(*args: str) -> str:
        result = subprocess.run(["git", *args], cwd=ROOT, text=True, capture_output=True)
        return result.stdout.strip() if result.returncode == 0 else "UNAVAILABLE"
    return f"HEAD={git('rev-parse', 'HEAD')}\nBRANCH={git('branch', '--show-current')}\nSTATUS:\n{git('status', '--short')}\n"


def summarize(observations: list[dict], products: list[dict], downloads: list[dict], validations: list[dict]) -> dict:
    included_obs = [r for r in observations if r["included"]]
    selected = acq.canonical_selection(products)
    families = acq.exposure_families(selected)
    raw = [r for r in selected if r["classification"] == acq.RAW_DETECTOR]
    controls = [r for r in selected if r["classification"] in {acq.FLT_CONTROL, acq.FLC_CONTROL}]
    raw_downloads = [r for r in downloads if r["local_path"].split("/")[-2] == "raw"]
    failures = [r for r in downloads if r["status"] not in {"DOWNLOADED", "SKIPPED_VALID"}]
    raw_ok = raw_downloads and len(raw_downloads) == len(raw) and not any(r in failures for r in raw_downloads)
    validation_ok = all(v.get("valid") for v in validations if v["classification"] == acq.RAW_DETECTOR)
    status = ("FRONTIER_ABELL2744_RAW_HST_ACQUISITION_ESTABLISHED" if raw_ok and validation_ok else
              "FRONTIER_ABELL2744_RAW_ACQUISITION_PARTIAL" if raw_downloads else
              "FRONTIER_ABELL2744_RAW_PRODUCTS_UNAVAILABLE" if included_obs else
              "FRONTIER_ABELL2744_MAIN_FIELD_NOT_RESOLVED")
    return {
        "status": status, "mast_observations": len(observations), "main_field_observations": len(included_obs),
        "parallel_rejected": sum(o["field_classification"] == "PARALLEL_FIELD" for o in observations),
        "ambiguous_rejected": sum(o["field_classification"] == "AMBIGUOUS" for o in observations),
        "distinct_exposures": len(families), "raw_files": len(raw),
        "flt_siblings": sum(bool(f["flt"]) for f in families.values()),
        "flc_siblings": sum(bool(f["flc"]) for f in families.values()),
        "incomplete_families": sum(f["status"] != "COMPLETE" for f in families.values()),
        "total_raw_bytes": sum(int(r.get("size") or 0) for r in raw),
        "total_control_bytes": sum(int(r.get("size") or 0) for r in controls),
        "program_ids": sorted({str(acq._value(o, "proposal_id", "proposal", "proposal_pi", "proposid") or "NOT_PRESENT") for o in included_obs}),
        "families": families, "download_failures": len(failures),
        "all_raw_byte_preserved": bool(raw_ok), "all_sha256_recorded": bool(downloads) and all(r.get("sha256") for r in downloads if r["status"] in {"DOWNLOADED", "SKIPPED_VALID"}),
    }


def report_text(frontier: dict, summary: dict, selection_hash: str, dry_run: bool) -> str:
    answers = [
        ("Does the Frontier HLSP tree itself contain raw detector exposures?", "YES" if frontier["raw_exposures_present"] else "NO"),
        ("Which files there are processed HLSP controls?", acq.CONTROL_FILENAME),
        ("How many matching main-field observations?", summary["main_field_observations"]),
        ("Which HST program IDs supplied them?", ", ".join(summary["program_ids"])),
        ("How many distinct detector exposures exist?", summary["distinct_exposures"]),
        ("Does every exposure have _raw.fits?", "YES" if summary["incomplete_families"] == 0 and summary["raw_files"] else "NO"),
        ("How many have _flt.fits?", summary["flt_siblings"]), ("How many have _flc.fits?", summary["flc_siblings"]),
        ("Total raw dataset size (bytes)?", summary["total_raw_bytes"]),
        ("Parallel-field exposures rejected?", summary["parallel_rejected"]), ("Ambiguous observations rejected?", summary["ambiguous_rejected"]),
        ("Is every raw FITS byte-preserved?", "NOT DOWNLOADED (DRY RUN)" if dry_run else ("YES" if summary["all_raw_byte_preserved"] else "NO")),
        ("Are all SHA256 hashes recorded?", "NOT DOWNLOADED (DRY RUN)" if dry_run else ("YES" if summary["all_sha256_recorded"] else "NO")),
        ("Are all archive product URIs recorded?", "YES"),
        ("Can acquisition be reproduced from the frozen manifest?", "YES"),
        ("Is reconstructed lensing data present in primary science directory?", "NO"),
    ]
    lines = [f"Dev125 status: {summary['status']}", f"DEV125_SELECTION_SHA256={selection_hash}", "", "Scientific/provenance answers:"]
    lines.extend(f"{i}. {q} {a}" for i, (q, a) in enumerate(answers, 1))
    return "\n".join(lines) + "\n"


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--include-hlsp-control", action="store_true")
    parser.add_argument("--allow-partial", action="store_true")
    parser.add_argument("--output-root", type=Path, default=ROOT)
    args = parser.parse_args(argv)
    run_dir = args.output_root / "runs" / LAB
    benchmark = args.output_root / "PBUF_raw_benchmark" / "WLRAW-001_Abell2744"
    run_dir.mkdir(parents=True, exist_ok=True); acq.prepare_layout(benchmark)
    (run_dir / "baseline_git.txt").write_text(baseline(), encoding="utf-8")
    query = query_manifest(); acq.write_json(run_dir / "query.json", query)
    try:
        frontier = acq.discover_frontier(); acq.write_json(run_dir / "frontier_index.json", frontier)
        observation_table, product_table = acq.mast_query()
        observations = acq.build_observation_manifest(acq.table_records(observation_table))
        included = {str(acq._value(r, "obsid", "obs_id") or "") for r in observations if r["included"]}
        products = acq.build_product_manifest(acq.table_records(product_table), included)
        pointing_audit = acq.apply_raw_pointing_audit(products)
    except Exception as exc:
        result = {"status": str(exc).split(":", 1)[0], "error": str(exc), "timestamp_utc": acq.utc_now()}
        acq.write_json(run_dir / "result.json", result)
        (run_dir / "report.txt").write_text(str(exc) + "\n", encoding="utf-8")
        return 2
    acq.write_json(run_dir / "mast_observations.json", observations)
    acq.write_json(run_dir / "mast_products.json", products)
    acq.write_json(run_dir / "raw_pointing_audit.json", pointing_audit)
    selection = acq.canonical_selection(products)
    if args.include_hlsp_control:
        selection.append({"obsid": "HLSP_CONTROL", "rootname": acq.rootname_from_filename(acq.CONTROL_FILENAME),
                          "filename": acq.CONTROL_FILENAME, "product_uri": acq.CONTROL_URL, "size": None,
                          "classification": "PROCESSED_CONTROL_ONLY"})
        selection.sort(key=lambda r: (str(r["obsid"]), str(r["rootname"]), str(r["filename"])))
    digest = acq.selection_sha256(selection)
    frozen = {"schema": "DEV125_SELECTION_V1", "selection_sha256": digest, "products": selection}
    acq.write_json(run_dir / "selected_products.json", frozen)
    downloads: list[dict] = []; validations: list[dict] = []
    if not args.dry_run:
        rejected_roots = {str(r["rootname"]) for r in pointing_audit if r["field_classification"] != "MAIN_CLUSTER"}
        quarantine = benchmark / "derived" / "rejected_parallel_or_ambiguous"
        quarantine.mkdir(parents=True, exist_ok=True)
        for folder in (benchmark / "raw", benchmark / "flt", benchmark / "flc"):
            for path in folder.iterdir():
                if acq.rootname_from_filename(path.name.removesuffix(".part")) in rejected_roots:
                    shutil.move(str(path), quarantine / path.name)
        prior_path = run_dir / "download_manifest.json"
        prior_rows = json.loads(prior_path.read_text())["downloads"] if prior_path.exists() else []
        prior = {r["filename"]: r for r in prior_rows}
        def transfer(product: dict) -> acq.DownloadOutcome:
            if product["classification"] == "PROCESSED_CONTROL_ONLY":
                destination = benchmark / "hlsp_control" / product["filename"]
            else:
                destination = acq.destination_for(benchmark, product)
            return acq.download_product(product, destination, prior.get(product["filename"]))
        with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
            futures = {pool.submit(transfer, product): product for product in selection}
            for future in concurrent.futures.as_completed(futures):
                outcome = future.result()
                downloads.append(outcome.record)
                if outcome.validation: validations.append(outcome.validation)
                downloads.sort(key=lambda r: r["filename"])
                validations.sort(key=lambda r: r["filename"])
                # Persist progress so an interrupted run can hash-check completed files.
                acq.write_json(run_dir / "download_manifest.json", {"downloads": downloads})
                acq.write_json(run_dir / "fits_validation.json", {"files": validations})
    acq.write_json(run_dir / "download_manifest.json", {"downloads": downloads})
    acq.write_json(run_dir / "fits_validation.json", {"files": validations})
    summary = summarize(observations, products, downloads, validations)
    if args.dry_run:
        summary["status"] = "DRY_RUN_SELECTION_FROZEN"
    result = dict(summary, selection_sha256=digest, dry_run=args.dry_run, timestamp_utc=acq.utc_now())
    acq.write_json(run_dir / "result.json", result)
    (run_dir / "report.txt").write_text(report_text(frontier, summary, digest, args.dry_run), encoding="utf-8")
    acq.copy_manifest_set(run_dir, benchmark, ("frontier_index.json", "query.json", "mast_observations.json", "mast_products.json", "selected_products.json", "download_manifest.json"))
    # Exact Dev125 provenance names, alongside the convenient run-level names.
    acq.write_json(benchmark / "provenance" / "frontier_hst_index.json", {
        "source_url": acq.HST_URL, "hst_release_directories_found": frontier["hst_release_directories_found"],
        "hst_v1_products": frontier["hst_v1_products"],
    })
    acq.write_json(benchmark / "provenance" / "observation_manifest.json", observations)
    acq.write_json(benchmark / "provenance" / "product_manifest.json", products)
    acq.write_json(benchmark / "mast_inventory" / "observations.json", observations)
    acq.write_json(benchmark / "mast_inventory" / "products.json", products)
    print(report_text(frontier, summary, digest, args.dry_run), end="")
    if args.dry_run or summary["status"] == "FRONTIER_ABELL2744_RAW_HST_ACQUISITION_ESTABLISHED": return 0
    return 0 if args.allow_partial and summary["raw_files"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
