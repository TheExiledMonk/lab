"""Target-blind, WCS-only coarse common-frame construction."""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np


@dataclass(frozen=True)
class CommonFrame:
    ra_min: float
    ra_max: float
    dec_min: float
    dec_max: float
    shape: tuple[int, int] = (128, 128)

    def bin_indices(self, ra, dec):
        # The archive is far from the RA wrap; preserve the explicit bounded rule.
        x = np.floor((ra - self.ra_min) / (self.ra_max - self.ra_min) * self.shape[1]).astype(int)
        y = np.floor((dec - self.dec_min) / (self.dec_max - self.dec_min) * self.shape[0]).astype(int)
        return y, x


def sampled_chip(sci, err, dq, celestial_wcs, samples=64):
    """Return deterministic detector samples and their celestial coordinates."""
    ny, nx = sci.shape
    ys = np.linspace(0, ny - 1, min(samples, ny), dtype=int)
    xs = np.linspace(0, nx - 1, min(samples, nx), dtype=int)
    yy, xx = np.meshgrid(ys, xs, indexing="ij")
    values = sci[yy, xx].astype(float)
    errors = np.ones_like(values) if err is None else err[yy, xx].astype(float)
    good = np.isfinite(values) & np.isfinite(errors) & (errors > 0)
    if dq is not None:
        good &= dq[yy, xx] == 0
    sky = celestial_wcs.pixel_to_world_values(xx, yy)
    return values, errors, good, np.asarray(sky[0]), np.asarray(sky[1])


def frame_from_bounds(bounds, shape=(128, 128)):
    ra = np.concatenate([x[0].ravel() for x in bounds]); dec = np.concatenate([x[1].ravel() for x in bounds])
    valid = np.isfinite(ra) & np.isfinite(dec)
    if not np.any(valid):
        raise ValueError("no finite celestial WCS samples")
    return CommonFrame(float(ra[valid].min()), float(ra[valid].max()),
        float(dec[valid].min()), float(dec[valid].max()), shape)


def combine_samples(samples, frame):
    """Inverse-variance mean with DQ rejection and a target-blind median background."""
    shape = frame.shape
    sum_w = np.zeros(shape); sum_wv = np.zeros(shape); coverage = np.zeros(shape, dtype=np.int32)
    backgrounds = []
    for values, errors, good, ra, dec in samples:
        finite = good & np.isfinite(ra) & np.isfinite(dec)
        background = float(np.median(values[finite])) if np.any(finite) else 0.0
        backgrounds.append(background)
        y, x = frame.bin_indices(ra[finite], dec[finite]); v = values[finite] - background
        w = 1.0 / np.square(errors[finite])
        inside = (x >= 0) & (x < shape[1]) & (y >= 0) & (y < shape[0])
        np.add.at(sum_w, (y[inside], x[inside]), w[inside])
        np.add.at(sum_wv, (y[inside], x[inside]), w[inside] * v[inside])
        np.add.at(coverage, (y[inside], x[inside]), 1)
    image = np.zeros(shape); uncertainty = np.full(shape, np.inf)
    occupied = sum_w > 0
    image[occupied] = sum_wv[occupied] / sum_w[occupied]
    uncertainty[occupied] = np.sqrt(1.0 / sum_w[occupied])
    return {"image": image, "uncertainty": uncertainty, "coverage": coverage,
        "backgrounds": backgrounds, "occupied": occupied}
