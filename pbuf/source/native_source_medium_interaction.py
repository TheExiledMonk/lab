"""Coefficient-free local source interaction on the frozen periodic N6 medium."""
from __future__ import annotations

import numpy as np

from pbuf.excitation.native_bond_state import relational_imbalance
from .native_source_state import NativeSourceState

N6_COORDINATION = 6


def source_imposed_excursion(shape: tuple[int, int, int], source: NativeSourceState) -> np.ndarray:
    """One occupied cell, with the constant zero mode removed as a gauge.

    The source contact itself is local.  The uniform subtraction is not a second
    interaction or a nonlocal force: a periodic relational Laplacian has no
    inverse on its constant mode, so this is the unique zero-mean compatibility
    projection.
    """
    source = source.wrapped(shape)
    out = np.zeros(shape, dtype=float)
    out[source.position] = source.amplitude
    out -= np.mean(out)
    return out


def source_medium_response(shape: tuple[int, int, int], source: NativeSourceState) -> np.ndarray:
    """F03 kick which translates the local equilibrium constraint."""
    return source_imposed_excursion(shape, source)


def stationary_response(shape: tuple[int, int, int], source: NativeSourceState) -> np.ndarray:
    """Solve imbalance(q)/6 + source_response == 0 in the zero-mean gauge."""
    forcing = source_medium_response(shape, source)
    axes = [2 * np.pi * np.fft.fftfreq(n) for n in shape]
    mesh = np.meshgrid(*axes, indexing="ij")
    # imbalance/6 has Fourier symbol -(2/3) sum sin^2(k_i/2).
    stiffness = (2.0 / 3.0) * sum(np.sin(k / 2.0) ** 2 for k in mesh)
    fhat = np.fft.fftn(forcing)
    qhat = np.zeros_like(fhat)
    mask = stiffness > 0
    qhat[mask] = fhat[mask] / stiffness[mask]
    q = np.fft.ifftn(qhat).real
    q -= np.mean(q)
    return q


def equilibrium_residual(q: np.ndarray, source: NativeSourceState) -> np.ndarray:
    return relational_imbalance(q) / N6_COORDINATION + source_medium_response(q.shape, source)


def medium_medium_response(q: np.ndarray) -> np.ndarray:
    """The frozen N6 restoring kick, named separately from source response."""
    return relational_imbalance(q) / N6_COORDINATION
