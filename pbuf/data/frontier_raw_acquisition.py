"""Deterministic acquisition primitives for Frontier Fields detector exposures.

This module performs archive discovery and byte-preserving acquisition only.  It
contains no calibration, image processing, source extraction, lensing, or cosmology.
Network entry points are restricted to STScI's Frontier archive and MAST.
"""
from __future__ import annotations

import hashlib
import json
import math
import os
import re
import shutil
import concurrent.futures
import time
from dataclasses import dataclass
from datetime import datetime, timezone
from html.parser import HTMLParser
from pathlib import Path
from typing import Any, Callable, Iterable, Mapping, Sequence
from urllib.parse import urljoin, urlparse

import requests
from astropy.io import fits

FRONTIER_ROOT = "https://archive.stsci.edu/pub/hlsp/frontier/"
CLUSTER_URL = urljoin(FRONTIER_ROOT, "abell2744/")
HST_URL = urljoin(CLUSTER_URL, "images/hst/")
CONTROL_FILENAME = "hlsp_frontier_hst_acs-60mas_abell2744_f814w_v1.0_drz.fits"
CONTROL_URL = urljoin(HST_URL, "v1.0/" + CONTROL_FILENAME)
TARGET_RA = 3.58611
TARGET_DEC = -30.40024
SEARCH_RADIUS_DEG = 0.25
# The main ACS mosaic is a few arcminutes wide.  This frozen region rejects the
# HFF parallel field (~6 arcmin away) while accepting overlapping ACS pointings.
MAIN_FIELD_RADIUS_DEG = 0.065
ALLOWED_HOSTS = {"archive.stsci.edu", "mast.stsci.edu"}
RETRY_WAITS = (2, 4, 8, 16, 32)

RAW_DETECTOR = "RAW_DETECTOR"
FLT_CONTROL = "FLT_CALIBRATED_CONTROL"
FLC_CONTROL = "FLC_CTE_CALIBRATED_CONTROL"


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


def json_safe(value: Any) -> Any:
    """Convert numpy/astropy values and masked values to stable JSON values."""
    if value is None:
        return None
    try:
        if bool(getattr(value, "mask", False)):
            return None
    except (TypeError, ValueError):
        pass
    if isinstance(value, (str, int, bool)):
        return value
    if isinstance(value, float):
        return value if math.isfinite(value) else None
    if hasattr(value, "item"):
        try:
            return json_safe(value.item())
        except (ValueError, TypeError):
            pass
    if isinstance(value, Mapping):
        return {str(k): json_safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(v) for v in value]
    return str(value)


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    text = json.dumps(json_safe(value), indent=2, sort_keys=True, ensure_ascii=True) + "\n"
    path.write_text(text, encoding="utf-8")


class _Links(HTMLParser):
    def __init__(self) -> None:
        super().__init__()
        self.hrefs: list[str] = []

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        if tag.lower() == "a":
            href = dict(attrs).get("href")
            if href:
                self.hrefs.append(href)


def _allowed_url(url: str) -> None:
    parsed = urlparse(url)
    if parsed.scheme != "https" or parsed.hostname not in ALLOWED_HOSTS:
        raise ValueError(f"network URL outside frozen STScI allowlist: {url}")


def fetch_index(url: str, session: requests.Session | None = None) -> list[dict[str, Any]]:
    """Return a normalized Apache-index listing, never a rendered page snapshot."""
    _allowed_url(url)
    client = session or requests.Session()
    last_error: Exception | None = None
    for attempt, wait in enumerate(RETRY_WAITS, 1):
        try:
            response = client.get(url, timeout=60)
            response.raise_for_status()
            break
        except requests.RequestException as exc:
            last_error = exc
            if attempt < len(RETRY_WAITS):
                time.sleep(wait)
    else:
        raise RuntimeError(f"FRONTIER_DISCOVERY_FAILED: {last_error}") from last_error
    parser = _Links()
    parser.feed(response.text)
    rows: list[dict[str, Any]] = []
    for href in parser.hrefs:
        if href in ("../", "./", "/") or href.startswith(("?", "#")):
            continue
        absolute = urljoin(url, href)
        if urlparse(absolute).hostname != "archive.stsci.edu":
            continue
        name = href.rstrip("/").split("/")[-1]
        rows.append({"name": name, "url": absolute, "is_directory": href.endswith("/")})
    return sorted({r["url"]: r for r in rows}.values(), key=lambda r: (r["name"], r["url"]))


def classify_frontier_file(filename: str) -> str:
    name = filename.lower().split("?")[0]
    if Path(name).name == CONTROL_FILENAME:
        return "HLSP_DRIZZLED_CONTROL"
    if name.endswith(("_raw.fits", "_raw.fits.gz")):
        return RAW_DETECTOR
    if name.endswith(("_flt.fits", "_flt.fits.gz")):
        return FLT_CONTROL
    if name.endswith(("_flc.fits", "_flc.fits.gz")):
        return FLC_CONTROL
    if name.endswith(("_drz.fits", "_drc.fits", "_drz.fits.gz", "_drc.fits.gz")):
        return "HLSP_DRIZZLED"
    if name.endswith(("_wht.fits", "_wht.fits.gz")):
        return "HLSP_WEIGHT"
    if name.endswith(("_rms.fits", "_rms.fits.gz")):
        return "HLSP_RMS"
    if any(x in name for x in ("catalog", "catalogue", ".cat", ".csv")):
        return "CATALOG"
    if any(x in name for x in ("model", "kappa", "shear", "magnification", "massmap")):
        return "MODEL"
    if name.endswith((".txt", ".md", ".pdf", ".html")) or "readme" in name:
        return "DOCUMENTATION"
    return "UNKNOWN"


def classify_hst_product(filename: str) -> str:
    cls = classify_frontier_file(filename)
    return cls if cls in {RAW_DETECTOR, FLT_CONTROL, FLC_CONTROL} else (
        "HST_DRIZZLED" if cls in {"HLSP_DRIZZLED", "HLSP_DRIZZLED_CONTROL"} else "UNKNOWN_HST_PRODUCT"
    )


def rootname_from_filename(filename: str) -> str:
    base = Path(filename.lower().removesuffix(".gz")).name
    match = re.match(r"(.+?)_(?:raw|flt|flc|drz|drc|wht|rms)\.fits$", base)
    return match.group(1) if match else base.removesuffix(".fits")


def discover_frontier(session: requests.Session | None = None) -> dict[str, Any]:
    root = fetch_index(FRONTIER_ROOT, session)
    cluster = fetch_index(CLUSTER_URL, session)
    hst = fetch_index(HST_URL, session)
    releases = [r for r in hst if r["is_directory"] and re.match(r"v\d", r["name"])]
    v1_url = next((r["url"] for r in releases if r["name"] == "v1.0"), urljoin(HST_URL, "v1.0/"))
    v1 = fetch_index(v1_url, session)
    files = [dict(r, classification=classify_frontier_file(r["name"])) for r in v1 if not r["is_directory"]]
    raw_present = any(r["classification"] == RAW_DETECTOR for r in files)
    return {
        "retrieval_timestamp_utc": utc_now(),
        "source_url": FRONTIER_ROOT,
        "cluster_directories_found": [r["name"] for r in root if r["is_directory"]],
        "abell2744_subdirectories_found": [r["name"] for r in cluster if r["is_directory"]],
        "hst_release_directories_found": [r["name"] for r in releases],
        "hst_v1_products": files,
        "raw_exposures_present": raw_present,
        "raw_status": "PRESENT" if raw_present else "FRONTIER_HLSP_RAW_EXPOSURES_NOT_PRESENT",
    }


def mast_query() -> tuple[Any, Any]:
    """Query MAST using its recommended Astroquery interface."""
    try:
        from astroquery.mast import Observations
        from astropy.coordinates import SkyCoord
        import astropy.units as u
    except ImportError:
        return _mast_api_query()
    coord = SkyCoord(TARGET_RA, TARGET_DEC, unit="deg", frame="icrs")
    observations = Observations.query_criteria(
        coordinates=coord,
        radius=SEARCH_RADIUS_DEG * u.deg,
        obs_collection="HST",
        instrument_name="ACS/WFC",
        filters="F814W",
        dataproduct_type="image",
    )
    try:
        products = Observations.get_product_list(observations)
    except Exception as exc:
        raise RuntimeError(f"PRODUCT_ENUMERATION_FAILED: {exc}") from exc
    return observations, products


def _mast_invoke(service: str, params: Mapping[str, Any], pagesize: int = 5000) -> list[dict[str, Any]]:
    """Supported MAST JSON API fallback when optional Astroquery is unavailable."""
    url = "https://mast.stsci.edu/api/v0/invoke"
    _allowed_url(url)
    rows: list[dict[str, Any]] = []
    page = 1
    while True:
        request = {"service": service, "params": dict(params), "format": "json", "pagesize": pagesize, "page": page}
        last_error: Exception | None = None
        payload: dict[str, Any] = {}
        for attempt, wait in enumerate(RETRY_WAITS, 1):
            try:
                response = requests.post(url, data={"request": json.dumps(request)}, timeout=180)
                response.raise_for_status()
                payload = response.json()
                break
            except (requests.RequestException, ValueError) as exc:
                last_error = exc
                if attempt < len(RETRY_WAITS):
                    time.sleep(wait)
        else:
            raise RuntimeError(f"MAST_QUERY_FAILED: {last_error}") from last_error
        if payload.get("status") != "COMPLETE":
            raise RuntimeError(f"MAST_QUERY_FAILED: {payload.get('msg') or payload.get('status')}")
        rows.extend(dict(row) for row in payload.get("data", []))
        paging = payload.get("paging") or {}
        if page >= int(paging.get("pagesFiltered") or 1):
            return rows
        page += 1


def _mast_api_query() -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    observations = _mast_invoke("Mast.Caom.Filtered.Position", {
        "columns": "*",
        "filters": [
            {"paramName": "obs_collection", "values": ["HST"]},
            {"paramName": "instrument_name", "values": ["ACS/WFC"]},
            {"paramName": "filters", "values": ["F814W"]},
        ],
        "position": f"{TARGET_RA}, {TARGET_DEC}, {SEARCH_RADIUS_DEG}",
    })
    obsids = ",".join(str(row["obsid"]) for row in observations)
    products = _mast_invoke("Mast.Caom.Products", {"obsid": obsids}) if obsids else []
    return observations, products


def table_records(table: Any) -> list[dict[str, Any]]:
    if isinstance(table, (list, tuple)):
        return [dict(json_safe(row)) for row in table]
    names = list(getattr(table, "colnames", []))
    return [{name: json_safe(row[name]) for name in names} for row in table]


def _value(row: Mapping[str, Any], *names: str) -> Any:
    lower = {str(k).lower(): v for k, v in row.items()}
    return next((lower[n.lower()] for n in names if n.lower() in lower), None)


def angular_distance_deg(ra1: float, dec1: float, ra2: float, dec2: float) -> float:
    r1, d1, r2, d2 = map(math.radians, (ra1, dec1, ra2, dec2))
    x = math.sin(d1) * math.sin(d2) + math.cos(d1) * math.cos(d2) * math.cos(r1 - r2)
    return math.degrees(math.acos(max(-1.0, min(1.0, x))))


def _footprint_points(s_region: Any) -> list[tuple[float, float]]:
    if not isinstance(s_region, str):
        return []
    nums = [float(x) for x in re.findall(r"[-+]?\d+(?:\.\d*)?(?:[Ee][-+]?\d+)?", s_region)]
    # STC-S commonly starts "POLYGON ICRS"; numeric pairs then follow.
    return list(zip(nums[0::2], nums[1::2]))


def classify_field(row: Mapping[str, Any]) -> str:
    """Classify from footprint/pointing only; names deliberately play no role."""
    points = _footprint_points(_value(row, "s_region", "footprint"))
    ra = _value(row, "s_ra", "ra", "RA")
    dec = _value(row, "s_dec", "dec", "DEC")
    if ra is not None and dec is not None:
        points.append((float(ra), float(dec)))
    if not points:
        return "AMBIGUOUS"
    distances = [angular_distance_deg(TARGET_RA, TARGET_DEC, p[0], p[1]) for p in points]
    if min(distances) <= MAIN_FIELD_RADIUS_DEG:
        return "MAIN_CLUSTER"
    # Query candidates outside the frozen main region are parallel/off-field.
    if min(distances) <= SEARCH_RADIUS_DEG + 0.05:
        return "PARALLEL_FIELD"
    return "AMBIGUOUS"


def observation_inclusion(row: Mapping[str, Any]) -> tuple[bool, str, str]:
    collection = str(_value(row, "obs_collection") or "").upper()
    instrument = str(_value(row, "instrument_name", "instrument") or "").upper()
    filt = str(_value(row, "filters", "filter") or "").upper()
    rights = str(_value(row, "dataRights", "data_rights") or "PUBLIC").upper()
    intent = str(_value(row, "intentType", "intent_type") or "SCIENCE").upper()
    field = classify_field(row)
    if collection != "HST": return False, "NOT_HST", field
    if instrument not in {"ACS/WFC", "ACS_WFC"}: return False, "NOT_ACS_WFC", field
    if filt != "F814W": return False, "NOT_F814W", field
    if "SCIENCE" not in intent: return False, "NOT_SCIENCE", field
    if rights not in {"PUBLIC", ""}: return False, "PROTECTED_PRODUCT_UNEXPECTED", field
    if field != "MAIN_CLUSTER": return False, field, field
    return True, "MATCHES_FROZEN_TARGET", field


def build_observation_manifest(rows: Sequence[Mapping[str, Any]]) -> list[dict[str, Any]]:
    output = []
    for source in rows:
        row = dict(json_safe(source))
        included, reason, field = observation_inclusion(row)
        row.update({"included": included, "inclusion_reason": reason, "field_classification": field})
        output.append(row)
    return sorted(output, key=lambda r: str(_value(r, "obsid", "obs_id") or ""))


def build_product_manifest(
    rows: Sequence[Mapping[str, Any]], included_obsids: set[str]
) -> list[dict[str, Any]]:
    output = []
    for source in rows:
        row = dict(json_safe(source))
        filename = str(_value(row, "productFilename", "filename") or "")
        obsid = str(_value(row, "parent_obsid", "obsid", "obs_id") or "")
        classification = classify_hst_product(filename)
        selected_class = classification in {RAW_DETECTOR, FLT_CONTROL, FLC_CONTROL}
        included = obsid in included_obsids and selected_class
        reason = "SELECTED_RAW_OR_CONTROL" if included else (
            "OBSERVATION_EXCLUDED" if obsid not in included_obsids else "PRODUCT_CLASS_EXCLUDED"
        )
        row.update({
            "obsid": obsid,
            "rootname": rootname_from_filename(filename),
            "filename": filename,
            "product_uri": _value(row, "dataURI", "productURI", "product_uri"),
            "product_type": _value(row, "productType", "product_type"),
            "product_subgroup": _value(row, "productSubGroupDescription", "product_subgroup"),
            "size": _value(row, "size", "size_bytes"),
            "classification": classification,
            "included": included,
            "reason": reason,
        })
        output.append(row)
    return sorted(output, key=lambda r: (r["obsid"], r["rootname"], r["filename"]))


SELECTION_FIELDS = ("obsid", "rootname", "filename", "product_uri", "size", "classification")


def canonical_selection(products: Iterable[Mapping[str, Any]]) -> list[dict[str, Any]]:
    candidates = [row for row in products if row.get("included")]
    raw_roots = {str(row.get("rootname")) for row in candidates if row.get("classification") == RAW_DETECTOR}
    selected = [
        {k: json_safe(row.get(k)) for k in SELECTION_FIELDS}
        for row in candidates
        if row.get("classification") == RAW_DETECTOR or str(row.get("rootname")) in raw_roots
    ]
    # CAOM aggregate observation rows can return the same archive product under
    # multiple parents.  A URI/filename is one immutable archive byte stream and
    # must occur only once in the download selection.
    selected.sort(key=lambda r: (str(r["obsid"]), str(r["rootname"]), str(r["filename"])))
    unique: dict[tuple[str, str], dict[str, Any]] = {}
    for row in selected:
        unique.setdefault((str(row["product_uri"]), str(row["filename"])), row)
    return sorted(unique.values(), key=lambda r: (str(r["obsid"]), str(r["rootname"]), str(r["filename"])))


def selection_sha256(selection: Iterable[Mapping[str, Any]]) -> str:
    canonical = json.dumps(list(selection), sort_keys=True, separators=(",", ":"), ensure_ascii=True)
    return hashlib.sha256(canonical.encode("utf-8")).hexdigest()


def exposure_families(products: Sequence[Mapping[str, Any]]) -> dict[str, dict[str, Any]]:
    families: dict[str, dict[str, Any]] = {}
    for row in products:
        root = str(row["rootname"])
        family = families.setdefault(root, {"rootname": root, "raw": [], "flt": [], "flc": [], "other": []})
        key = {RAW_DETECTOR: "raw", FLT_CONTROL: "flt", FLC_CONTROL: "flc"}.get(str(row["classification"]), "other")
        family[key].append(row["filename"])
    for family in families.values():
        family["raw_exists"] = bool(family["raw"])
        family["status"] = "COMPLETE" if all(family[x] for x in ("raw", "flt", "flc")) else (
            "RAW_PRODUCT_MISSING" if not family["raw"] else "CALIBRATION_CONTROLS_INCOMPLETE"
        )
    return dict(sorted(families.items()))


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def resume_action(path: Path, expected_size: int | None, previous_sha256: str | None) -> str:
    part = path.with_name(path.name + ".part")
    if path.exists():
        size_ok = expected_size in (None, 0) or path.stat().st_size == expected_size
        hash_ok = previous_sha256 is None or sha256_file(path) == previous_sha256
        return "SKIP_VALID" if size_ok and hash_ok else "REDOWNLOAD_CORRUPT"
    if part.exists():
        return "REDOWNLOAD_PARTIAL"
    return "DOWNLOAD"


def product_download_url(uri: str) -> str:
    if uri.startswith("mast:"):
        from urllib.parse import quote
        return "https://mast.stsci.edu/api/v0.1/Download/file?uri=" + quote(uri, safe=":/")
    _allowed_url(uri)
    return uri


def resolve_raw_pointing(product: Mapping[str, Any], session: requests.Session | None = None) -> dict[str, Any]:
    """Resolve a raw exposure pointing from a range-read primary FITS header."""
    url = product_download_url(str(product.get("product_uri") or ""))
    client = session or requests.Session()
    last_error = ""
    for attempt, wait in enumerate(RETRY_WAITS, 1):
        try:
            response = client.get(url, headers={"Range": "bytes=0-65535"}, timeout=(10, 15))
            response.raise_for_status()
            data = response.content
            end = next((offset + 80 for offset in range(0, len(data) - 79, 80)
                        if data[offset:offset + 8] == b"END     "), None)
            if end is None:
                raise ValueError("FITS primary header END card not found in range response")
            header = fits.Header.fromstring(data[:end].decode("ascii", errors="strict"), sep="")
            ra, dec = header.get("RA_TARG"), header.get("DEC_TARG")
            field = classify_field({"s_ra": ra, "s_dec": dec})
            return {"rootname": product["rootname"], "filename": product["filename"],
                    "RA_TARG": json_safe(ra), "DEC_TARG": json_safe(dec),
                    "field_classification": field, "status": "HEADER_RANGE_RESOLVED"}
        except Exception as exc:
            last_error = str(exc)
            if attempt < len(RETRY_WAITS):
                time.sleep(wait)
    return {"rootname": product["rootname"], "filename": product["filename"],
            "field_classification": "AMBIGUOUS", "status": "HEADER_RANGE_FAILED", "error": last_error}


def apply_raw_pointing_audit(products: list[dict[str, Any]], workers: int = 2) -> list[dict[str, Any]]:
    """Range-resolve raw pointings and propagate decisions to exact-root siblings."""
    raw_by_root: dict[str, dict[str, Any]] = {}
    for row in products:
        if row.get("included") and row.get("classification") == RAW_DETECTOR:
            raw_by_root.setdefault(str(row["rootname"]), row)
    audits: list[dict[str, Any]] = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
        futures = [pool.submit(resolve_raw_pointing, row) for row in raw_by_root.values()]
        for future in concurrent.futures.as_completed(futures):
            audits.append(future.result())
    by_root = {str(row["rootname"]): row for row in audits}
    for row in products:
        audit = by_root.get(str(row["rootname"]))
        if not audit:
            continue
        row["raw_header_pointing"] = audit
        if audit["field_classification"] != "MAIN_CLUSTER":
            row["included"] = False
            row["reason"] = audit["field_classification"] + "_RAW_HEADER_REJECTED"
    return sorted(audits, key=lambda r: (str(r["rootname"]), str(r["filename"])))


def validate_fits(path: Path, classification: str) -> dict[str, Any]:
    keys = ("TELESCOP", "INSTRUME", "DETECTOR", "FILTER1", "FILTER2", "ROOTNAME", "PROPOSID",
            "DATE-OBS", "TIME-OBS", "EXPTIME", "RA_TARG", "DEC_TARG", "PA_V3")
    result: dict[str, Any] = {"filename": path.name, "classification": classification, "valid": False}
    try:
        with fits.open(path, mode="readonly", checksum=False, memmap=False) as hdus:
            # Force every header and data span to be reachable, detecting truncated HDUs.
            for hdu in hdus:
                _ = hdu.header
                if hdu.data is not None:
                    _ = hdu.data.shape
                    if hdu.data.size:
                        _ = hdu.data.reshape(-1)[-1]
            header = hdus[0].header
            metadata = {key: json_safe(header.get(key, "NOT_PRESENT")) for key in keys}
        filters = {str(metadata["FILTER1"]).upper(), str(metadata["FILTER2"]).upper()}
        metadata_ok = (
            str(metadata["TELESCOP"]).upper() == "HST"
            and str(metadata["INSTRUME"]).upper() == "ACS"
            and str(metadata["DETECTOR"]).upper() == "WFC"
            and "F814W" in filters
        )
        field = classify_field({"s_ra": metadata["RA_TARG"], "s_dec": metadata["DEC_TARG"]})
        result.update(metadata=metadata, metadata_valid=metadata_ok, field_classification=field,
                      valid=metadata_ok and field == "MAIN_CLUSTER")
        if not result["valid"]:
            result["status"] = "FITS_VALIDATION_FAILED"
    except Exception as exc:
        result.update(status="FITS_VALIDATION_FAILED", error=str(exc))
    return result


@dataclass
class DownloadOutcome:
    record: dict[str, Any]
    validation: dict[str, Any] | None


def download_product(
    product: Mapping[str, Any], destination: Path, previous: Mapping[str, Any] | None = None,
    session: requests.Session | None = None, waits: Sequence[int] = RETRY_WAITS,
    sleeper: Callable[[float], None] = time.sleep,
) -> DownloadOutcome:
    destination.parent.mkdir(parents=True, exist_ok=True)
    expected = product.get("size")
    expected_size = int(expected) if expected not in (None, "", 0) else None
    previous_sha = str(previous.get("sha256")) if previous and previous.get("sha256") else None
    action = resume_action(destination, expected_size, previous_sha)
    uri = str(product.get("product_uri") or "")
    base = {"filename": destination.name, "source_uri": uri, "local_path": str(destination), "action": action}
    if action == "SKIP_VALID":
        digest = sha256_file(destination)
        record = dict(base, size_bytes=destination.stat().st_size, sha256=digest,
                      download_timestamp_utc=utc_now(), status="SKIPPED_VALID")
        return DownloadOutcome(record, validate_fits(destination, str(product["classification"])))
    part = destination.with_name(destination.name + ".part")
    part.unlink(missing_ok=True)
    client = session or requests.Session()
    url = product_download_url(uri)
    last_error = ""
    for attempt, wait in enumerate(waits, start=1):
        try:
            with client.get(url, stream=True, timeout=(30, 180)) as response:
                if response.status_code in (401, 403):
                    return DownloadOutcome(dict(base, status="PROTECTED_PRODUCT_UNEXPECTED",
                                                http_status=response.status_code), None)
                response.raise_for_status()
                with part.open("wb") as stream:
                    for chunk in response.iter_content(1024 * 1024):
                        if chunk:
                            stream.write(chunk)
            if expected_size is not None and part.stat().st_size != expected_size:
                raise IOError(f"DOWNLOAD_PARTIAL: expected {expected_size}, got {part.stat().st_size}")
            os.replace(part, destination)
            digest = sha256_file(destination)
            validation = validate_fits(destination, str(product["classification"]))
            status = "DOWNLOADED" if validation["valid"] else "FITS_VALIDATION_FAILED"
            return DownloadOutcome(dict(base, size_bytes=destination.stat().st_size, sha256=digest,
                                        download_timestamp_utc=utc_now(), http_status=response.status_code,
                                        attempts=attempt, status=status), validation)
        except Exception as exc:
            last_error = str(exc)
            part.unlink(missing_ok=True)
            if attempt < len(waits):
                sleeper(wait)
    return DownloadOutcome(dict(base, status="DOWNLOAD_FAILED", attempts=len(waits), error=last_error), None)


def destination_for(base: Path, product: Mapping[str, Any]) -> Path:
    folder = {RAW_DETECTOR: "raw", FLT_CONTROL: "flt", FLC_CONTROL: "flc"}[str(product["classification"])]
    return base / folder / str(product["filename"])


def prepare_layout(base: Path) -> None:
    for name in ("provenance", "mast_inventory", "raw", "flt", "flc", "hlsp_control", "derived"):
        (base / name).mkdir(parents=True, exist_ok=True)


def copy_manifest_set(run_dir: Path, benchmark_dir: Path, names: Sequence[str]) -> None:
    for name in names:
        source = run_dir / name
        if source.exists():
            shutil.copy2(source, benchmark_dir / "provenance" / name)
