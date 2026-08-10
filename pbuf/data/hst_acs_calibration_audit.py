"""Read-only, target-blind ACS/WFC RAW--FLT--FLC audit primitives."""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
from typing import Any, Iterable, Mapping

import numpy as np

CALIBRATION_SWITCHES = (
    "DQICORR", "ATODCORR", "BLEVCORR", "BIASCORR", "FLSHCORR", "CRCORR",
    "SHADCORR", "DARKCORR", "FLATCORR", "PHOTCORR", "PCTECORR", "LFLGCORR",
    "GLINCORR", "NLINCORR", "SINKCORR",
)
REFERENCE_KEYS = (
    "BPIXTAB", "CCDTAB", "ATODTAB", "OSCNTAB", "BIASFILE", "DARKFILE",
    "PFLTFILE", "DFLTFILE", "LFLTFILE", "PCTETAB",
)
GEOMETRY_KEYS = (
    "CCDCHIP", "LTV1", "LTV2", "LTM1_1", "LTM1_2", "LTM2_1", "LTM2_2",
    "BINAXIS1", "BINAXIS2", "SUBARRAY", "SIZAXIS1", "SIZAXIS2", "NAXIS1",
    "NAXIS2", "CTYPE1", "CTYPE2", "CRPIX1", "CRPIX2", "CRVAL1", "CRVAL2",
    "CD1_1", "CD1_2", "CD2_1", "CD2_2", "ORIENTAT",
)


def json_value(value: Any) -> Any:
    if value is None or isinstance(value, (str, int, bool)):
        return value
    if isinstance(value, float):
        return value if math.isfinite(value) else None
    if isinstance(value, np.generic):
        return json_value(value.item())
    if isinstance(value, Mapping):
        return {str(k): json_value(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_value(v) for v in value]
    return str(value)


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(json_value(value), sort_keys=True, indent=2) + "\n")


def canonical_sha256(value: Any) -> str:
    payload = json.dumps(json_value(value), sort_keys=True, separators=(",", ":"), ensure_ascii=True)
    return hashlib.sha256(payload.encode()).hexdigest()


def sha256_file(path: Path, block: int = 8 << 20) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(block), b""):
            digest.update(chunk)
    return digest.hexdigest()


def match_families(dataset: Path) -> list[dict[str, Any]]:
    stages: dict[str, dict[str, Path]] = {}
    for stage in ("raw", "flt", "flc"):
        for path in sorted((dataset / stage).glob(f"*_{stage}.fits")):
            root = path.name[: -len(f"_{stage}.fits")].lower()
            stages.setdefault(root, {})[stage] = path
    incomplete = {root: sorted(set(("raw", "flt", "flc")) - set(rows)) for root, rows in stages.items() if len(rows) != 3}
    if incomplete:
        raise RuntimeError(f"EXPOSURE_FAMILY_INCOMPLETE: {incomplete}")
    return [{"rootname": root, **{s: rows[s] for s in ("raw", "flt", "flc")}}
            for root, rows in sorted(stages.items())]


def robust_stats(data: np.ndarray, percentiles: Iterable[float] = (.1, 1, 5, 25, 50, 75, 95, 99, 99.9)) -> dict[str, Any]:
    a = np.asarray(data)
    finite = np.isfinite(a)
    x = np.asarray(a[finite], dtype=np.float64)
    out: dict[str, Any] = {"count": int(a.size), "count_finite": int(x.size),
                           "count_nan": int(np.isnan(a).sum()), "count_inf": int(np.isinf(a).sum())}
    if not x.size:
        return {**out, **{k: None for k in ("minimum", "maximum", "mean", "median", "std", "mad", "percentiles")}}
    med = float(np.median(x))
    out.update(minimum=float(x.min()), maximum=float(x.max()), mean=float(x.mean()), median=med,
               std=float(x.std()), mad=float(np.median(np.abs(x - med))),
               percentiles={str(p): float(v) for p, v in zip(percentiles, np.percentile(x, list(percentiles)))})
    return out


def geometry_from_header(header: Mapping[str, Any], shape: tuple[int, ...] | None) -> dict[str, Any]:
    return {"shape": list(shape) if shape is not None else None,
            **{key: json_value(header.get(key)) for key in GEOMETRY_KEYS}}


def geometry_classification(a: Mapping[str, Any], b: Mapping[str, Any]) -> str:
    if a.get("shape") != b.get("shape"):
        return "PIXEL_REINDEXED" if a.get("CCDCHIP") == b.get("CCDCHIP") else "UNKNOWN"
    if a.get("CCDCHIP") != b.get("CCDCHIP"):
        return "UNKNOWN"
    transform_keys = ("LTV1", "LTV2", "LTM1_1", "LTM1_2", "LTM2_1", "LTM2_2", "BINAXIS1", "BINAXIS2")
    if all(a.get(k) == b.get(k) for k in transform_keys):
        return "PIXEL_PRESERVING"
    return "PIXEL_REINDEXED"


def units_compatible(a: str | None, b: str | None) -> bool:
    norm = lambda x: " ".join(str(x or "").upper().split())
    return bool(norm(a)) and norm(a) == norm(b)


def direct_difference_status(a_geometry: Mapping[str, Any], b_geometry: Mapping[str, Any],
                             a_unit: str | None, b_unit: str | None) -> str:
    if geometry_classification(a_geometry, b_geometry) != "PIXEL_PRESERVING" or not units_compatible(a_unit, b_unit):
        return "DIRECT_DIFFERENCE_NOT_VALID"
    return "VALID"


def difference_stats(a: np.ndarray, b: np.ndarray, err: np.ndarray | None = None) -> dict[str, Any]:
    d = np.asarray(b, dtype=np.float64) - np.asarray(a, dtype=np.float64)
    good = np.isfinite(d)
    x = d[good]
    if not x.size:
        return {"status": "NO_FINITE_PIXELS"}
    med = float(np.median(x)); result = {
        "status": "VALID", "mean_difference": float(x.mean()), "median_difference": med,
        "std_difference": float(x.std()), "mad_difference": float(np.median(np.abs(x-med))),
        "rms_difference": float(np.sqrt(np.mean(x*x))), "max_abs_difference": float(np.max(np.abs(x))),
        "fraction_exactly_unchanged": float(np.count_nonzero(x == 0) / x.size),
    }
    if err is not None and np.asarray(err).shape == d.shape:
        e = np.asarray(err, dtype=np.float64); valid = good & np.isfinite(e) & (e > 0)
        for sigma in (1, 3, 5):
            result[f"fraction_changed_gt_{sigma}_sigma"] = float(np.mean(np.abs(d[valid]) > sigma*e[valid])) if valid.any() else None
    return result


def dq_comparison(a: np.ndarray, b: np.ndarray) -> dict[str, Any]:
    aa = np.asarray(a, dtype=np.uint64); bb = np.asarray(b, dtype=np.uint64)
    if aa.shape != bb.shape:
        return {"status": "DIRECT_DIFFERENCE_NOT_VALID", "changed_pixel_count": None}
    values, counts = np.unique(bb, return_counts=True)
    changed = aa != bb
    newly = (aa == 0) & (bb != 0)
    return {"status": "VALID", "changed_pixel_count": int(changed.sum()),
            "newly_flagged_count": int(newly.sum()), "fraction_flagged": float(np.mean(bb != 0)),
            "new_value_histogram": {str(int(v)): int(c) for v, c in zip(values, counts)}}


def classify_difference(delta: np.ndarray) -> dict[str, Any]:
    d = np.asarray(delta, dtype=np.float64); finite = np.isfinite(d)
    if not finite.any(): return {"pattern": "UNKNOWN", "information": "UNKNOWN"}
    x = d[finite]; scale = max(float(np.max(np.abs(x))), 1.0)
    if np.all(np.abs(x - np.median(x)) <= np.finfo(float).eps * scale * 8):
        pattern = "UNCHANGED" if np.all(x == 0) else "UNIFORM"
    else:
        # Median profiles are saved separately; RMS profiles retain sparse trails.
        row = np.sqrt(np.nanmean(d*d, axis=1)); col = np.sqrt(np.nanmean(d*d, axis=0))
        rv, cv, av = np.nanstd(row), np.nanstd(col), np.nanstd(d)
        pattern = "ROW_DEPENDENT" if rv > .2*av and rv > cv*1.25 else ("COLUMN_DEPENDENT" if cv > .2*av and cv > rv*1.25 else "LOCALIZED")
    changed = np.count_nonzero(x != 0) / x.size
    return {"pattern": pattern, "changed_fraction": float(changed),
            "row_profile_activity": float(np.nanstd(np.sqrt(np.nanmean(d*d, axis=1)))),
            "column_profile_activity": float(np.nanstd(np.sqrt(np.nanmean(d*d, axis=0)))),
            "information": "INFORMATION_MODIFYING" if 0 < changed < 1 else "PIXEL_PRESERVING"}


def detect_resampling(source: np.ndarray, output: np.ndarray, metadata_pixel_preserving: bool = False) -> str:
    """Conservative diagnostic: metadata proof wins; fractional interpolation is resampling."""
    a, b = np.asarray(source), np.asarray(output)
    if metadata_pixel_preserving and a.shape == b.shape:
        return "PIXEL_PRESERVING"
    if a.shape != b.shape:
        return "RESAMPLED"
    source_values = np.unique(a)
    sample = b.ravel()[::max(1, b.size // 10000)]
    if np.mean(np.isin(sample, source_values)) < .95:
        return "RESAMPLED"
    return "UNKNOWN"


def profiles(delta: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    d = np.asarray(delta, dtype=np.float64)
    return np.nanmedian(d, axis=1), np.nanmedian(d, axis=0)


def derivative_stats(data: np.ndarray) -> dict[str, Any]:
    a = np.asarray(data, dtype=np.float64)
    dy, dx = np.gradient(a); dyy, dxy = np.gradient(dy); _, dxx = np.gradient(dx)
    return {"dx": robust_stats(dx), "dy": robust_stats(dy), "gradient_magnitude": robust_stats(np.hypot(dx, dy)),
            "dxx": robust_stats(dxx), "dxy": robust_stats(dxy), "dyy": robust_stats(dyy),
            "gradient_energy": float(np.nanmean(dx*dx + dy*dy)),
            "curvature_energy": float(np.nanmean(dxx*dxx + 2*dxy*dxy + dyy*dyy))}


def frequency_summary(data: np.ndarray) -> dict[str, float]:
    a = np.asarray(data, dtype=np.float64); a = np.nan_to_num(a - np.nanmedian(a))
    # Deterministic stride bounds FFT memory/time while retaining detector-wide coverage.
    stride = max(1, int(math.ceil(max(a.shape) / 1024)))
    z = a[::stride, ::stride]; power = np.abs(np.fft.rfft2(z)) ** 2
    fy = np.fft.fftfreq(z.shape[0])[:, None]; fx = np.fft.rfftfreq(z.shape[1])[None, :]
    radius = np.hypot(fx, fy)
    return {name: float(power[mask].mean()) if mask.any() else 0.0 for name, mask in {
        "low_frequency_power": radius < .05, "mid_frequency_power": (radius >= .05) & (radius < .2),
        "high_frequency_power": radius >= .2}.items()}


def patch_audit(a: np.ndarray, b: np.ndarray, sizes: Iterable[int] = (32, 64, 128)) -> dict[str, Any]:
    aa = np.asarray(a, dtype=np.float64); bb = np.asarray(b, dtype=np.float64); out = {}
    for size in sizes:
        rows = []
        for y in range(0, aa.shape[0] - size + 1, size):
            for x in range(0, aa.shape[1] - size + 1, size):
                p, q = aa[y:y+size, x:x+size], bb[y:y+size, x:x+size]
                if not (np.isfinite(p).all() and np.isfinite(q).all()): continue
                pdy,pdx=np.gradient(p); qdy,qdx=np.gradient(q)
                pdyy,pdxy=np.gradient(pdy); _,pdxx=np.gradient(pdx)
                qdyy,qdxy=np.gradient(qdy); _,qdxx=np.gradient(qdx)
                pge=float(np.mean(pdx*pdx+pdy*pdy)); qge=float(np.mean(qdx*qdx+qdy*qdy))
                pce=float(np.mean(pdxx*pdxx+2*pdxy*pdxy+pdyy*pdyy)); qce=float(np.mean(qdxx*qdxx+2*qdxy*qdxy+qdyy*qdyy))
                corr = np.corrcoef(p.ravel(), q.ravel())[0, 1] if p.std() and q.std() else None
                rows.append([float(np.mean(q-p)), float(np.var(q)-np.var(p)),
                             qge-pge, qce-pce, corr])
        arr = np.asarray([[np.nan if v is None else v for v in row] for row in rows])
        out[str(size)] = {"patch_count": len(rows), "median_metrics": np.nanmedian(arr, axis=0).tolist() if len(rows) else None,
                          "metric_order": ["mean_change", "variance_change", "gradient_energy_change", "curvature_energy_change", "cross_correlation"]}
    return out


def hdu_inventory(hdul: Any) -> list[dict[str, Any]]:
    result = []
    for index, hdu in enumerate(hdul):
        shape = tuple(hdu.shape) if getattr(hdu, "shape", None) else None
        bitpix = hdu.header.get("BITPIX")
        dtype = f"FITS_BITPIX_{bitpix}" if bitpix is not None else type(hdu).__name__
        result.append({"index": index, "EXTNAME": hdu.header.get("EXTNAME", "PRIMARY"), "EXTVER": hdu.header.get("EXTVER"),
                       "shape": list(shape) if shape else None, "dtype": dtype,
                       "BITPIX": hdu.header.get("BITPIX"), "NAXIS": hdu.header.get("NAXIS"),
                       "header_keyword_count": len(hdu.header), "data_present": bool(hdu.header.get("NAXIS", 0))})
    return result


def stage_transforms(headers: Mapping[str, Mapping[str, Any]]) -> dict[str, Any]:
    raw, flt, flc = headers["raw"], headers["flt"], headers["flc"]
    pctecorr = str(flc.get("PCTECORR", "")).upper()
    primary = "PIXEL_BASED_CTE_CORRECTION" if pctecorr == "COMPLETE" and str(flt.get("PCTECORR", "")).upper() != "COMPLETE" else "UNKNOWN"
    return {"raw_to_flt_categories": [name for key, name in (("BLEVCORR", "ELECTRONICS_CORRECTION"), ("BIASCORR", "BIAS_CORRECTION"),
             ("DARKCORR", "DARK_CORRECTION"), ("FLATCORR", "FLAT_FIELD_CORRECTION"), ("DQICORR", "DATA_QUALITY_PROPAGATION"),
             ("PHOTCORR", "PHOTOMETRIC_METADATA")) if str(flt.get(key, "")).upper() == "COMPLETE"],
            "FLT_TO_FLC_PRIMARY_TRANSFORM": primary}
