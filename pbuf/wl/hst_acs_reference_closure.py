"""Official acquisition and offline verification of ACS calibration references."""
from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable
import hashlib
import urllib.request
from astropy.io import fits

CRDS_BASE = "https://hst-crds.stsci.edu/unchecked_get/references/hst"


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def _inventory(path: Path) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    # Headers and table/image shapes only. Reference calibration arrays are not HST science arrays.
    with fits.open(path, mode="readonly", memmap=True, lazy_load_hdus=True) as hdul:
        primary = hdul[0].header
        keys = ("FILETYPE", "PEDIGREE", "USEAFTER", "DESCRIP", "INSTRUME", "DETECTOR")
        metadata = {key: primary.get(key) for key in keys}
        extensions = []
        for index, hdu in enumerate(hdul):
            header = hdu.header
            extensions.append({"index": index, "extname": header.get("EXTNAME", "PRIMARY"),
                               "extver": header.get("EXTVER"), "naxis": header.get("NAXIS", 0),
                               "shape_from_header": [header.get(f"NAXIS{i}") for i in range(1, header.get("NAXIS", 0)+1)],
                               "ccdchip": header.get("CCDCHIP"), "filter1": header.get("FILTER1"),
                               "filter2": header.get("FILTER2"), "bunit": header.get("BUNIT")})
    return metadata, extensions


def acquire_references(logical_names: Iterable[str], cache: Path, *, offline: bool = False) -> dict[str, Any]:
    cache.mkdir(parents=True, exist_ok=True)
    records, missing = [], []
    for logical in logical_names:
        filename = logical.split("$", 1)[-1]
        target = cache / filename
        url = f"{CRDS_BASE}/{filename}"
        retrieved = None
        if not target.exists() and not offline:
            with urllib.request.urlopen(url, timeout=60) as response:
                payload = response.read()
            target.write_bytes(payload)
            retrieved = datetime.now(timezone.utc).isoformat()
        if not target.exists():
            missing.append(logical)
            continue
        try:
            metadata, extensions = _inventory(target)
        except Exception as exc:
            missing.append(logical)
            records.append({"logical_reference_name": logical, "resolved_filename": filename,
                            "status": "CORRUPT_FITS", "error": str(exc)})
            continue
        records.append({"logical_reference_name": logical, "resolved_filename": filename,
                        "source_url_service": url, "retrieval_utc": retrieved,
                        "byte_size": target.stat().st_size, "sha256": sha256_file(target),
                        "fits_primary_metadata": metadata, "reference_pedigree": metadata.get("PEDIGREE"),
                        "useafter_date_applicability": metadata.get("USEAFTER"), "extensions": extensions,
                        "status": "VERIFIED_LOCAL_OFFLINE" if offline else "ACQUIRED_OR_VERIFIED"})
    return {"reference_closure_complete": not missing, "records": records,
            "missing_references": missing, "official_service": CRDS_BASE, "offline": offline}
