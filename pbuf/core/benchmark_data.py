"""Canonical local PBUF weak-lensing benchmark loader.

This module is the single repository-level entry point for the five frozen
Frontier Fields benchmark products already stored under ``PBUF_benchmark``.
It performs local filesystem I/O only: no URL discovery, downloads, network
fallback, normalization, source construction, or physics transformations.

Labs should import this module rather than rebuilding benchmark directory/file
names independently.
"""
from __future__ import annotations

from pathlib import Path
from typing import Mapping

import numpy as np
from astropy.io import fits

ROOT = Path(__file__).resolve().parents[2]
BENCHMARK_ROOT = ROOT / "PBUF_benchmark"

CLUSTERS = (
    {"id": "Abell2744", "label": "Abell 2744", "slug": "abell2744", "directory": "WL-001_Abell2744", "slug_dir": "abell_2744"},
    {"id": "MACS0416", "label": "MACS J0416", "slug": "macs0416", "directory": "WL-002_MACS0416", "slug_dir": "macs_j0416"},
    {"id": "MACS1149", "label": "MACS J1149", "slug": "macs1149", "directory": "WL-003_MACS1149", "slug_dir": "macs_j1149"},
    {"id": "AbellS1063", "label": "Abell S1063", "slug": "abells1063", "directory": "WL-004_AbellS1063", "slug_dir": "abell_s1063"},
    {"id": "Abell370", "label": "Abell 370", "slug": "abell370", "directory": "WL-005_Abell370", "slug_dir": "abell_370"},
)

PRODUCTS = ("kappa", "gamma", "gamma1", "gamma2")
_BY_ID = {row["id"]: row for row in CLUSTERS}


def clusters() -> tuple[dict, ...]:
    return tuple(dict(row) for row in CLUSTERS)


def resolve_cluster(cluster: str | Mapping[str, object]) -> dict:
    if isinstance(cluster, str):
        try:
            return dict(_BY_ID[cluster])
        except KeyError as exc:
            raise KeyError(f"unknown PBUF benchmark cluster id: {cluster}") from exc
    cid = str(cluster.get("id", ""))
    if cid in _BY_ID:
        return dict(_BY_ID[cid])
    raise KeyError(f"unknown PBUF benchmark cluster mapping id: {cid!r}")


def product_path(cluster: str | Mapping[str, object], product: str) -> Path:
    if product not in PRODUCTS:
        raise KeyError(f"unknown weak-lensing product: {product!r}; expected one of {PRODUCTS}")
    row = resolve_cluster(cluster)
    return BENCHMARK_ROOT / row["directory"] / f"hlsp_frontier_model_{row['slug']}_merten_v1_{product}.fits"


def require_product_path(cluster: str | Mapping[str, object], product: str) -> Path:
    path = product_path(cluster, product)
    if not path.is_file():
        raise FileNotFoundError(f"local PBUF benchmark FITS not found: {path}")
    return path


def load_product(cluster: str | Mapping[str, object], product: str) -> np.ndarray:
    path = require_product_path(cluster, product)
    with fits.open(path, memmap=True) as hdul:
        return np.asarray(hdul[0].data, dtype=np.float64).copy()


def kappa_path(cluster: str | Mapping[str, object]) -> Path:
    return product_path(cluster, "kappa")


def require_kappa_path(cluster: str | Mapping[str, object]) -> Path:
    return require_product_path(cluster, "kappa")


def load_kappa(cluster: str | Mapping[str, object]) -> np.ndarray:
    return load_product(cluster, "kappa")


def load_gamma(cluster: str | Mapping[str, object]) -> np.ndarray:
    return load_product(cluster, "gamma")


def load_gamma1(cluster: str | Mapping[str, object]) -> np.ndarray:
    return load_product(cluster, "gamma1")


def load_gamma2(cluster: str | Mapping[str, object]) -> np.ndarray:
    return load_product(cluster, "gamma2")


def load_header_shape(cluster: str | Mapping[str, object], product: str = "kappa"):
    path = require_product_path(cluster, product)
    hdr = fits.getheader(path, 0)
    shape = (int(hdr["NAXIS2"]), int(hdr["NAXIS1"]))
    return hdr, shape


def inventory() -> list[dict]:
    rows = []
    for cluster in CLUSTERS:
        products = {name: {"path": str(product_path(cluster, name)), "exists": product_path(cluster, name).is_file()} for name in PRODUCTS}
        rows.append({**dict(cluster), "products": products, "all_products_exist": all(v["exists"] for v in products.values())})
    return rows
