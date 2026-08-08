"""Canonical local PBUF weak-lensing benchmark loader.

This module is the single repository-level entry point for the five frozen
Frontier Fields benchmark FITS files already stored under ``PBUF_benchmark``.
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

# Frozen canonical five-cluster inventory.  These values are consolidated from
# the long-standing a8_three_dimensional_projection_lab001 benchmark definition
# and the existing PBUF_benchmark directory layout.
CLUSTERS = (
    {
        "id": "Abell2744",
        "label": "Abell 2744",
        "slug": "abell2744",
        "directory": "WL-001_Abell2744",
        "slug_dir": "abell_2744",
    },
    {
        "id": "MACS0416",
        "label": "MACS J0416",
        "slug": "macs0416",
        "directory": "WL-002_MACS0416",
        "slug_dir": "macs_j0416",
    },
    {
        "id": "MACS1149",
        "label": "MACS J1149",
        "slug": "macs1149",
        "directory": "WL-003_MACS1149",
        "slug_dir": "macs_j1149",
    },
    {
        "id": "AbellS1063",
        "label": "Abell S1063",
        "slug": "abells1063",
        "directory": "WL-004_AbellS1063",
        "slug_dir": "abell_s1063",
    },
    {
        "id": "Abell370",
        "label": "Abell 370",
        "slug": "abell370",
        "directory": "WL-005_Abell370",
        "slug_dir": "abell_370",
    },
)

_BY_ID = {row["id"]: row for row in CLUSTERS}


def clusters() -> tuple[dict, ...]:
    """Return copies of the frozen five-cluster records."""
    return tuple(dict(row) for row in CLUSTERS)


def resolve_cluster(cluster: str | Mapping[str, object]) -> dict:
    """Resolve a canonical cluster record from an id or compatible mapping."""
    if isinstance(cluster, str):
        try:
            return dict(_BY_ID[cluster])
        except KeyError as exc:
            raise KeyError(f"unknown PBUF benchmark cluster id: {cluster}") from exc

    cid = str(cluster.get("id", ""))
    if cid in _BY_ID:
        return dict(_BY_ID[cid])
    raise KeyError(f"unknown PBUF benchmark cluster mapping id: {cid!r}")


def kappa_path(cluster: str | Mapping[str, object]) -> Path:
    """Return the canonical local Merten v1 kappa FITS path."""
    row = resolve_cluster(cluster)
    return (
        BENCHMARK_ROOT
        / row["directory"]
        / f"hlsp_frontier_model_{row['slug']}_merten_v1_kappa.fits"
    )


def require_kappa_path(cluster: str | Mapping[str, object]) -> Path:
    """Return the local path, raising clearly if the repository data is missing."""
    path = kappa_path(cluster)
    if not path.is_file():
        raise FileNotFoundError(
            f"local PBUF benchmark FITS not found: {path}. "
            "This loader has no network fallback."
        )
    return path


def load_kappa(cluster: str | Mapping[str, object]) -> np.ndarray:
    """Load the canonical local kappa image as float64."""
    path = require_kappa_path(cluster)
    with fits.open(path, memmap=True) as hdul:
        return np.asarray(hdul[0].data, dtype=np.float64).copy()


def load_header_shape(cluster: str | Mapping[str, object]):
    """Load primary FITS header and (ny, nx) shape without loading image values."""
    path = require_kappa_path(cluster)
    hdr = fits.getheader(path, 0)
    shape = (int(hdr["NAXIS2"]), int(hdr["NAXIS1"]))
    return hdr, shape


def inventory() -> list[dict]:
    """Return local existence/path metadata for all five frozen benchmarks."""
    rows = []
    for cluster in CLUSTERS:
        path = kappa_path(cluster)
        rows.append(
            {
                **dict(cluster),
                "path": str(path),
                "exists": path.is_file(),
                "network_used": False,
            }
        )
    return rows
