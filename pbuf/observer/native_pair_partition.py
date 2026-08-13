"""Geometry-only, periodic N6 Voronoi partition for a frozen pair of centers."""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np


def periodic_n6_distance(shape: tuple[int, int, int], center: np.ndarray) -> np.ndarray:
    """Canonical shortest-path distance on the periodic Cartesian N6 lattice."""
    points = np.indices(shape).reshape(3, -1).T
    sizes = np.asarray(shape, dtype=int)
    delta = np.abs((points - np.asarray(center, dtype=int) + sizes // 2) % sizes - sizes // 2)
    return np.sum(delta, axis=1).reshape(shape)


@dataclass(frozen=True)
class NativePairPartition:
    shape: tuple[int, int, int]
    center_a: tuple[int, int, int]
    center_b: tuple[int, int, int]
    omega_a: np.ndarray
    omega_b: np.ndarray
    omega_i: np.ndarray
    omega_d: np.ndarray
    distance_a: np.ndarray
    distance_b: np.ndarray


def derive_partition(shape: tuple[int, int, int], center_a, center_b) -> NativePairPartition:
    """Derive a disjoint three-region N6 Voronoi partition; ties are interface."""
    ca = tuple(np.asarray(center_a, dtype=int).tolist())
    cb = tuple(np.asarray(center_b, dtype=int).tolist())
    da, db = periodic_n6_distance(shape, ca), periodic_n6_distance(shape, cb)
    ai, bi = da < db, db < da
    interface = da == db
    domain = np.ones(shape, dtype=bool)
    return NativePairPartition(tuple(shape), ca, cb, ai, bi, interface, domain, da, db)


def translate_partition(partition: NativePairPartition, shift) -> NativePairPartition:
    """Translate centers and domain by an exact periodic integer translation."""
    s = np.asarray(shift, dtype=int)
    shape = np.asarray(partition.shape, dtype=int)
    return derive_partition(tuple(shape), (np.asarray(partition.center_a) + s) % shape,
                            (np.asarray(partition.center_b) + s) % shape)
