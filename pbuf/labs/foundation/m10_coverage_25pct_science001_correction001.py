#!/usr/bin/env python3
"""Narrow correction for the M10 25% coverage science lab.

Only the broken diagnostic correlation helper is corrected: the reviewed M16
public API is safe_pearson / safe_spearman, not pearson / spearman. All physics,
coverage geometry, photon counts, propagation, Jacobian logic, and outputs stay
in the original lab unchanged.
"""
from __future__ import annotations

import sys
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))

from pbuf.core import observable_extraction as M16
from pbuf.labs.foundation import m10_coverage_25pct_science001 as lab


def corrected_corr(a, b):
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    mask = np.isfinite(a) & np.isfinite(b)
    n = int(mask.sum())
    if n < 2:
        return float("nan"), float("nan"), n
    return (
        float(M16.safe_pearson(a[mask], b[mask])),
        float(M16.safe_spearman(a[mask], b[mask])),
        n,
    )


if __name__ == "__main__":
    lab._corr = corrected_corr
    raise SystemExit(lab.main())
