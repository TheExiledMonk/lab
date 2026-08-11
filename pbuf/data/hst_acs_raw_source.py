"""Read-only inventory and calibrated-pixel access for HST ACS exposure families.

This module contains instrument semantics only.  In particular, it never turns
detector intensity into mass, convergence, shear, or a three-dimensional state.
"""
from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS


@dataclass(frozen=True)
class ExposureFamily:
    rootname: str
    raw: Path
    flt: Path
    flc: Path


def exposure_families(archive: Path) -> list[ExposureFamily]:
    groups = {}
    for kind in ("raw", "flt", "flc"):
        for path in sorted((archive / kind).glob(f"*_{kind}.fits")):
            root = path.name.removesuffix(f"_{kind}.fits")
            groups.setdefault(root, {})[kind] = path
    incomplete = {k: sorted(set(("raw", "flt", "flc")) - set(v)) for k, v in groups.items() if len(v) != 3}
    if incomplete:
        raise ValueError(f"incomplete RAW/FLT/FLC families: {incomplete}")
    return [ExposureFamily(k, **groups[k]) for k in sorted(groups)]


def _json(value):
    if isinstance(value, np.generic):
        return value.item()
    return value


def effective_filter(header) -> str:
    values = [str(header.get(k, "")).strip() for k in ("FILTER1", "FILTER2")]
    science = [x for x in values if x and x.upper() not in {"CLEAR", "CLEAR1L", "CLEAR2L", "N/A"}]
    return science[0] if science else "+".join(x for x in values if x) or "UNKNOWN"


def inventory_family(family: ExposureFamily) -> dict:
    products = {}
    for kind, path in (("raw", family.raw), ("flt", family.flt), ("flc", family.flc)):
        with fits.open(path, mode="readonly", memmap=False, lazy_load_hdus=True) as hdul:
            ph = hdul[0].header
            extensions = []
            for hdu in hdul:
                naxis = int(hdu.header.get("NAXIS", 0))
                shape = [int(hdu.header.get(f"NAXIS{i}", 0)) for i in range(naxis, 0, -1)] if naxis else None
                extensions.append({"name": hdu.name, "extver": _json(hdu.header.get("EXTVER")),
                    "shape": shape, "bunit": hdu.header.get("BUNIT"),
                    "ccdchip": _json(hdu.header.get("CCDCHIP"))})
            products[kind] = {"filename": path.name, "calibration": {k: ph.get(k) for k in
                ("CAL_VER", "BLEVCORR", "BIASCORR", "DARKCORR", "FLATCORR", "PCTECORR")},
                "extensions": extensions}
    with fits.open(family.raw, mode="readonly", memmap=False, lazy_load_hdus=True) as hdul:
        h = hdul[0].header
        metadata = {"rootname": family.rootname, "exposure_id": h.get("ROOTNAME", family.rootname),
            "detector": h.get("DETECTOR"), "filter": effective_filter(h), "exposure_time": _json(h.get("EXPTIME")),
            "pointing": {k: _json(h.get(k)) for k in ("RA_TARG", "DEC_TARG")},
            "orientation": {k: _json(h.get(k)) for k in ("PA_V3", "ORIENTAT")}}
    metadata["products"] = products
    return metadata


def archive_inventory(archive: Path) -> tuple[list[ExposureFamily], dict]:
    families = exposure_families(archive)
    rows = [inventory_family(f) for f in families]
    counts = {kind.upper(): sum(1 for _ in (archive / kind).glob(f"*_{kind}.fits")) for kind in ("raw", "flt", "flc")}
    return families, {"archive": str(archive), "family_count": len(families), "file_counts": counts,
        "filters": dict(Counter(r["filter"] for r in rows)), "exposures": rows}


def calibrated_chips(path: Path):
    """Yield SCI/ERR/DQ chip arrays, releasing each FITS family after use."""
    with fits.open(path, mode="readonly", memmap=False, lazy_load_hdus=True) as hdul:
        primary = hdul[0].header
        for ver in (1, 2):
            try:
                sci = hdul[("SCI", ver)]
            except KeyError:
                continue
            if sci.data is None:
                continue
            err = hdul[("ERR", ver)].data if ("ERR", ver) in hdul else None
            dq = hdul[("DQ", ver)].data if ("DQ", ver) in hdul else None
            # WCS is extension-local; passing the HDUList resolves ACS lookup-table
            # distortion without allowing the primary HDU's NAXIS=0 to shadow SCI.
            celestial_wcs = WCS(sci.header, hdul).celestial
            yield ver, np.asarray(sci.data, dtype=np.float32), (
                None if err is None else np.asarray(err, dtype=np.float32)), (
                None if dq is None else np.asarray(dq)), celestial_wcs
