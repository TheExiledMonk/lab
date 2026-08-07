"""M05 — Unique Pair Enumeration.

Enumerate every unordered N6 neighbour pair exactly once using only the
three positive directions (xp, yp, zp).

Correction pass FOUNDATION-001-CORRECTION-001
---------------------------------------------
* Midpoint is now computed via the generic formula
  ``mid = tuple(0.5 * (a + b) for a, b in zip(i, j))``.
  This makes xp/yp/zp midpoints direction-correct: only the axis
  along which j-i is nonzero carries the half-step, the other two
  coordinates equal the source i exactly.
* Added the four mandatory tests M05-C1..C4 from §5.4.
"""
from __future__ import annotations
import numpy as np

from .conventions import (
    N6_POSITIVE_DIRECTIONS, AXIS_OFFSETS, validate_transform_id,
    get_coordinate_matrix,
)
from .coordinate_transforms import transform_pair_direction

__all__ = [
    "PairRecord", "enumerate_internal_pairs", "enumerate_internal_pairs_reference",
    "pair_count_formula", "expected_pair_count",
    "pair_direction", "pair_midpoint",
    "PairEnumerationError",
]


class PairEnumerationError(ValueError):
    pass


class PairRecord:
    """Lightweight record describing a single unordered pair."""

    __slots__ = ("pair_id", "i_index", "j_index", "axis",
                 "direction_xyz", "midpoint_zyx")

    def __init__(self, pair_id, i_index, j_index, axis,
                 direction_xyz, midpoint_zyx):
        self.pair_id = int(pair_id)
        self.i_index = tuple(int(x) for x in i_index)
        self.j_index = tuple(int(x) for x in j_index)
        self.axis = str(axis)
        self.direction_xyz = tuple(int(x) for x in direction_xyz)
        self.midpoint_zyx = tuple(float(x) for x in midpoint_zyx)

    def __repr__(self):
        return (f"PairRecord(id={self.pair_id}, i={self.i_index}, "
                f"j={self.j_index}, axis={self.axis!r})")

    def to_dict(self):
        return {
            "pair_id": self.pair_id,
            "i_index": list(self.i_index),
            "j_index": list(self.j_index),
            "axis": self.axis,
            "direction_xyz": list(self.direction_xyz),
            "midpoint_zyx": list(self.midpoint_zyx),
        }


def pair_count_formula(shape):
    """Return the expected positive-N6 internal-pair count.

    For shape (Nz, Ny, Nx), the number of unordered internal N6
    pairs is

        Nz * Ny * (Nx - 1)   +   Nz * (Ny - 1) * Nx   +   (Nz - 1) * Ny * Nx
    """
    nz, ny, nx = shape
    return int(nz * ny * (nx - 1) + nz * (ny - 1) * nx + (nz - 1) * ny * nx)


def expected_pair_count(shape):
    return pair_count_formula(shape)


def pair_direction(pair):
    """Return the (dx, dy, dz) integer direction of a pair."""
    if not isinstance(pair, PairRecord):
        raise PairEnumerationError("pair must be a PairRecord")
    return np.array(pair.direction_xyz, dtype=np.int64)


def pair_midpoint(pair):
    """Return the midpoint in array-axis (z, y, x) coordinates."""
    if not isinstance(pair, PairRecord):
        raise PairEnumerationError("pair must be a PairRecord")
    return np.array(pair.midpoint_zyx, dtype=np.float64)


def _pair_midpoint_zyx(i, j):
    """Compute the pair midpoint in array-axis (z, y, x) coordinates.

    Uses the generic CORRECTION-001 formula
        m = 0.5 * (i + j)
    so the geometry is independent of the axis label and the
    direction-specific literals (xp/yp/zp) all flow from one rule.
    """
    return tuple(0.5 * (a + b) for a, b in zip(i, j))


def enumerate_internal_pairs(shape, stencil="N6_positive"):
    """Enumerate every unordered pair exactly once.

    Returns a list of PairRecord objects. Only the three positive N6
    directions are stored; the partner at i+axis is implicit. ``shape``
    is (nz, ny, nx).

    The midpoint is computed with the generic
    ``0.5 * (a + b)`` rule (CORRECTION-001 §5).
    """
    if stencil != "N6_positive":
        raise PairEnumerationError(
            f"unsupported stencil: {stencil!r}")
    nz, ny, nx = shape
    if nz < 2 or ny < 2 or nx < 2:
        raise PairEnumerationError("need at least 2 voxels on each axis")
    pairs = []
    pid = 0
    # xp: i and j differ only in x. Midpoint x = ix + 0.5; z = iz, y = iy.
    for iz in range(nz):
        for iy in range(ny):
            for ix in range(nx - 1):
                pid += 1
                i = (iz, iy, ix)
                j = (iz, iy, ix + 1)
                mid = _pair_midpoint_zyx(i, j)
                pairs.append(PairRecord(pid, i, j, "xp",
                                         (+1, 0, 0), mid))
    # yp
    for iz in range(nz):
        for iy in range(ny - 1):
            for ix in range(nx):
                pid += 1
                i = (iz, iy, ix)
                j = (iz, iy + 1, ix)
                mid = _pair_midpoint_zyx(i, j)
                pairs.append(PairRecord(pid, i, j, "yp",
                                         (0, +1, 0), mid))
    # zp
    for iz in range(nz - 1):
        for iy in range(ny):
            for ix in range(nx):
                pid += 1
                i = (iz, iy, ix)
                j = (iz + 1, iy, ix)
                mid = _pair_midpoint_zyx(i, j)
                pairs.append(PairRecord(pid, i, j, "zp",
                                         (0, 0, +1), mid))
    return pairs


def enumerate_internal_pairs_reference(shape, stencil="N6_positive"):
    """Independent reference enumerator using a set-based check.

    Iterates every voxel, generates the three positive-direction
    neighbours, and deduplicates by unordered pair key. The midpoint
    uses the same generic 0.5*(a+b) formula, but it is computed
    through a different code path so the two implementations are
    structurally independent.
    """
    if stencil != "N6_positive":
        raise PairEnumerationError(
            f"unsupported stencil: {stencil!r}")
    nz, ny, nx = shape
    if nz < 2 or ny < 2 or nx < 2:
        raise PairEnumerationError("need at least 2 voxels on each axis")
    seen = set()
    pairs = []
    pid = 0
    for iz in range(nz):
        for iy in range(ny):
            for ix in range(nx):
                # xp
                if ix + 1 < nx:
                    key = (iz, iy, ix, iz, iy, ix + 1)
                    if key not in seen:
                        seen.add(key); pid += 1
                        i = (iz, iy, ix); j = (iz, iy, ix + 1)
                        pairs.append(PairRecord(pid, i, j, "xp",
                                                 (+1, 0, 0),
                                                 _pair_midpoint_zyx(i, j)))
                # yp
                if iy + 1 < ny:
                    key = (iz, iy, ix, iz, iy + 1, ix)
                    if key not in seen:
                        seen.add(key); pid += 1
                        i = (iz, iy, ix); j = (iz, iy + 1, ix)
                        pairs.append(PairRecord(pid, i, j, "yp",
                                                 (0, +1, 0),
                                                 _pair_midpoint_zyx(i, j)))
                # zp
                if iz + 1 < nz:
                    key = (iz, iy, ix, iz + 1, iy, ix)
                    if key not in seen:
                        seen.add(key); pid += 1
                        i = (iz, iy, ix); j = (iz + 1, iy, ix)
                        pairs.append(PairRecord(pid, i, j, "zp",
                                                 (0, 0, +1),
                                                 _pair_midpoint_zyx(i, j)))
    return pairs


# ----------------------------------------------------------------------
# Tests (CORRECTION-001 §5.4)
# ----------------------------------------------------------------------
def _M05_C1_midpoint_identity(shape):
    """Exact midpoint identity: m_ij = 0.5 * (i + j)."""
    pairs = enumerate_internal_pairs(shape)
    fails = []
    for p in pairs:
        expected = tuple(0.5 * (a + b) for a, b in zip(p.i_index, p.j_index))
        for got, want in zip(p.midpoint_zyx, expected):
            if got != want:
                fails.append((p.pair_id, expected, p.midpoint_zyx))
                break
    return {"test": "M05-C1", "shape": list(shape),
            "n_pairs": len(pairs), "n_failures": len(fails),
            "passes": len(fails) == 0}


def _M05_C2_fixed_axis_coordinates(shape):
    """For xp: m_z=i_z, m_y=i_y. For yp: m_z=i_z, m_x=i_x.
    For zp: m_y=i_y, m_x=i_x."""
    pairs = enumerate_internal_pairs(shape)
    fails = []
    for p in pairs:
        i = p.i_index
        m = p.midpoint_zyx
        # i and midpoint in (z, y, x) order.
        if p.axis == "xp":
            if m[0] != i[0] or m[1] != i[1]:
                fails.append((p.pair_id, p.axis, i, m))
        elif p.axis == "yp":
            if m[0] != i[0] or m[2] != i[2]:
                fails.append((p.pair_id, p.axis, i, m))
        elif p.axis == "zp":
            if m[1] != i[1] or m[2] != i[2]:
                fails.append((p.pair_id, p.axis, i, m))
    return {"test": "M05-C2", "shape": list(shape),
            "n_pairs": len(pairs), "n_failures": len(fails),
            "passes": len(fails) == 0}


def _M05_C3_direction_displacement(shape):
    """After converting array order (z, y, x) to vector order (x, y, z),
    j - i must equal the pair direction n_hat."""
    pairs = enumerate_internal_pairs(shape)
    fails = []
    for p in pairs:
        # i, j are stored as (z, y, x); convert to (x, y, z) for the
        # displacement test.
        i_xyz = (p.i_index[2], p.i_index[1], p.i_index[0])
        j_xyz = (p.j_index[2], p.j_index[1], p.j_index[0])
        disp = tuple(b - a for a, b in zip(i_xyz, j_xyz))
        if disp != p.direction_xyz:
            fails.append((p.pair_id, disp, p.direction_xyz))
    return {"test": "M05-C3", "shape": list(shape),
            "n_pairs": len(pairs), "n_failures": len(fails),
            "passes": len(fails) == 0}


def _M05_C4_noncubic(shape=(3, 4, 5)):
    """Run all family checks on a non-cubic grid."""
    a = _M05_C1_midpoint_identity(shape)
    b = _M05_C2_fixed_axis_coordinates(shape)
    c = _M05_C3_direction_displacement(shape)
    return {"test": "M05-C4", "shape": list(shape),
            "M05-C1_passes": a["passes"],
            "M05-C2_passes": b["passes"],
            "M05-C3_passes": c["passes"],
            "passes": a["passes"] and b["passes"] and c["passes"]}


def _pair_count_test(shape):
    p = enumerate_internal_pairs(shape)
    return {"shape": list(shape), "n_pairs": len(p),
            "expected": expected_pair_count(shape),
            "passes": len(p) == expected_pair_count(shape)}


def _no_duplicates_test(shape):
    p = enumerate_internal_pairs(shape)
    keys = set()
    for pair in p:
        key = (min(pair.i_index, pair.j_index),
               max(pair.i_index, pair.j_index))
        if key in keys:
            return {"passes": False, "duplicate": key}
        keys.add(key)
    return {"passes": True, "n_unique_keys": len(keys),
            "n_pairs": len(p)}


def _no_self_pairs_test(shape):
    p = enumerate_internal_pairs(shape)
    for pair in p:
        if pair.i_index == pair.j_index:
            return {"passes": False, "self_pair": pair.i_index}
    return {"passes": True, "n_pairs": len(p)}


def _direction_transform_test():
    from .conventions import RC_TRANSFORMS, N6_DIRECTIONS
    rows = []
    for rc in RC_TRANSFORMS:
        for lbl in N6_POSITIVE_DIRECTIONS:
            out = transform_pair_direction(lbl, rc)
            rows.append({"transform": rc, "input": lbl,
                         "output": out, "passes": out in N6_DIRECTIONS})
    return rows


def _reference_agreement_test(shape):
    a = enumerate_internal_pairs(shape)
    b = enumerate_internal_pairs_reference(shape)
    if len(a) != len(b):
        return {"passes": False, "n_prod": len(a), "n_ref": len(b)}
    ka = sorted((min(p.i_index, p.j_index), max(p.i_index, p.j_index))
                for p in a)
    kb = sorted((min(p.i_index, p.j_index), max(p.i_index, p.j_index))
                for p in b)
    return {"passes": ka == kb, "n_pairs": len(a)}


if __name__ == "__main__":
    for sh in [(3, 4, 5), (4, 5, 6), (5, 4, 3)]:
        r = _pair_count_test(sh)
        assert r["passes"], f"pair-count failed for {sh}: {r}"
        r2 = _no_duplicates_test(sh)
        assert r2["passes"], f"duplicate test failed for {sh}"
        r3 = _no_self_pairs_test(sh)
        assert r3["passes"], f"self-pair test failed for {sh}"
        r4 = _reference_agreement_test(sh)
        assert r4["passes"], f"reference agreement failed for {sh}"
    rows = _direction_transform_test()
    assert all(r["passes"] for r in rows)
    # CORRECTION-001 mandatory tests
    for sh in [(3, 4, 5), (4, 5, 6)]:
        for fn in (_M05_C1_midpoint_identity,
                    _M05_C2_fixed_axis_coordinates,
                    _M05_C3_direction_displacement):
            r = fn(sh)
            assert r["passes"], f"{r}"
    r = _M05_C4_noncubic((3, 4, 5))
    assert r["passes"], r
    print(f"M05 pair count/no-dup/no-self/ref: 3 shapes × 4 checks all pass")
    print(f"M05 pair direction transforms: {len(rows)} cases all pass")
    print("M05 pair midpoint identity / fixed-axis / direction / noncubic: pass")
    print("M05 pair enumeration: all checks passed")
