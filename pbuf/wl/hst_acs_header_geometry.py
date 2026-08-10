"""Header-only ACS/WFC geometry inventory (Dev134).

This module deliberately never accesses ``HDU.data``.  Image dimensions are
read from FITS cards, so callers can prove that science arrays were not read.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Mapping

from astropy.io import fits

from pbuf.data.hst_acs_calibration_audit import canonical_sha256, match_families

SOURCE_CLASSES = ("LOCAL_FITS_HEADER", "LOCAL_REFERENCE_FILE",
                  "OFFICIAL_INSTRUMENT_DOCUMENT", "DERIVED_FROM_FROZEN_SOURCE")
GEOMETRY_KEYS = ("CCDCHIP", "NAXIS1", "NAXIS2", "LTV1", "LTV2", "LTM1_1", "LTM2_2",
                 "CRPIX1", "CRPIX2", "CRVAL1", "CRVAL2", "CD1_1", "CD1_2", "CD2_1", "CD2_2",
                 "CTYPE1", "CTYPE2", "ORIENTAT", "VAFACTOR", "IDCTAB", "NPOLFILE", "D2IMFILE")
PRIMARY_KEYS = ("ROOTNAME", "DETECTOR", "FILTER1", "FILTER2", "DATE-OBS", "TIME-OBS", "PA_V3",
                "IDCTAB", "NPOLFILE", "D2IMFILE")


@dataclass(frozen=True)
class GeometryValue:
    value: Any
    units: str | None
    source_class: str
    source_file: str
    source_location: str
    derivation: str = "direct FITS card"


@dataclass(frozen=True)
class ACSExposureGeometry:
    exposure_uid: str
    rootname: str
    raw_filename: str
    flt_filename: str
    flc_filename: str
    detector: str | None
    filter: str | None
    date_obs: str | None
    time_obs: str | None
    header_geometry_source: str
    chips: tuple[Mapping[str, GeometryValue], ...]
    orientation_metadata: Mapping[str, GeometryValue]
    reference_file_metadata: Mapping[str, GeometryValue]

    def manifest(self) -> dict[str, Any]:
        return asdict(self)


def _card(header: fits.Header, key: str, path: Path, location: str) -> GeometryValue:
    return GeometryValue(header.get(key), header.comments[key] if key in header else None,
                         "LOCAL_FITS_HEADER", str(path.resolve()), f"{location}:{key}")


def read_headers(path: Path) -> tuple[fits.Header, list[fits.Header]]:
    """Read only primary and SCI headers; no FITS array is opened or sampled."""
    primary = fits.getheader(path, 0)
    sci = []
    # ACS calibrated products conventionally have SCI,1 and SCI,2.  Discover
    # them through header reads and stop at the first missing extension.
    for extver in (1, 2):
        try:
            sci.append(fits.getheader(path, ("SCI", extver)))
        except (KeyError, IndexError, OSError):
            break
    return primary, sci


def load_exposure_geometry(row: Mapping[str, Path]) -> ACSExposureGeometry:
    path = Path(row["flc"])
    primary, sci = read_headers(path)
    if len(sci) != 2:
        raise RuntimeError(f"ACS_GEOMETRY_SOURCE_CONFLICT_UNRESOLVED: {path}: expected two SCI headers")
    chips = tuple({k: _card(h, k, path, f"SCI,{i}") for k in GEOMETRY_KEYS if k in h}
                  for i, h in enumerate(sci, 1))
    root = str(primary.get("ROOTNAME") or row["rootname"]).lower()
    filt = primary.get("FILTER2") if str(primary.get("FILTER2", "")).startswith("F") else primary.get("FILTER1")
    orient = {k: _card(primary, k, path, "PRIMARY") for k in ("PA_V3",) if k in primary}
    for i, h in enumerate(sci, 1):
        for k in ("ORIENTAT", "VAFACTOR"):
            if k in h: orient[f"CHIP{i}_{k}"] = _card(h, k, path, f"SCI,{i}")
    refs = {k: _card(primary, k, path, "PRIMARY") for k in ("IDCTAB", "NPOLFILE", "D2IMFILE") if k in primary}
    return ACSExposureGeometry(root, root, str(Path(row["raw"]).resolve()), str(Path(row["flt"]).resolve()),
        str(path.resolve()), primary.get("DETECTOR"), filt, primary.get("DATE-OBS"), primary.get("TIME-OBS"),
        "FLC primary plus SCI extension headers", chips, orient, refs)


def inventory(dataset: Path) -> tuple[ACSExposureGeometry, ...]:
    return tuple(load_exposure_geometry(row) for row in match_families(dataset))


def exposure_manifest_sha256(exposures: tuple[ACSExposureGeometry, ...]) -> str:
    return canonical_sha256([e.manifest() for e in exposures])

