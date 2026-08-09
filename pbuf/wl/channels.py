"""Frozen 45-channel decoder bank with injectable spatial deposition."""

from contextlib import contextmanager
import threading
import os
import numpy as np

from pbuf.labs.foundation import native_full_received_state_information_retention001 as RET
from pbuf.labs.foundation import native_observable_extraction_method_sweep001 as EX
from .config import EXTENT, OBS_BINS
from .deposition import DepositionMethod, get_deposition_method
from .backends.vulkan_kde import make_kde_backend
from .observer_cache import ObserverPrimitiveCache, ObserverStateId


_HISTOGRAM_LOCK = threading.RLock()


@contextmanager
def _deposition_histogram(method: DepositionMethod):
    """Route frozen formula histogram calls through one selected depositor."""
    if method.name == "hard_bin_current":
        yield
        return
    original = np.histogram2d

    def injected(y, x, bins=10, range=None, density=None, weights=None):
        if density:
            raise ValueError("density-normalized histogram deposition is unsupported")
        if not isinstance(bins, (tuple, list)) or len(bins) != 2:
            return original(y, x, bins=bins, range=range, density=density, weights=weights)
        y_edges, x_edges = (np.asarray(edge, dtype=np.float64) for edge in bins)
        if (y_edges.ndim != 1 or x_edges.ndim != 1 or
                y_edges.size != x_edges.size or
                not np.array_equal(y_edges, x_edges) or
                not np.isclose(x_edges[0], -x_edges[-1])):
            raise ValueError("deposition audit supports square, symmetric detector grids only")
        count = method.deposit(x, y, weights, bins=x_edges.size - 1,
                               extent=float(x_edges[-1]))
        return count, x_edges, y_edges

    with _HISTOGRAM_LOCK:
        np.histogram2d = injected
        try:
            yield
        finally:
            np.histogram2d = original


def decode_full_channel_bank(
    screen: dict,
    received_state: dict,
    deposition_method: str | DepositionMethod | None = None,
    *, cache: ObserverPrimitiveCache | None = None,
    state_id: ObserverStateId | None = None,
    kde_backend=None,
) -> dict:
    method = get_deposition_method(deposition_method)
    backend = kde_backend
    owns_backend = False
    if backend is None:
        backend = make_kde_backend(os.environ.get("PBUF_OBSERVER_KDE_BACKEND", "cpu")); owns_backend = True
    original_diag = EX.OLD._diag_kde
    kde_ordinal = 0
    def cached_diag(data, bandwidth):
        nonlocal kde_ordinal
        role = "initial" if kde_ordinal == 0 else "final"
        kde_ordinal += 1
        reference = original_diag(data, bandwidth)
        def evaluate(points):
            # Only the all-ray self-query is O(N^2); preserve the frozen grid path.
            if points is not data and not np.array_equal(points, data): return reference(points)
            sid = state_id or ObserverStateId("anonymous_content_state", backend=backend.name)
            key = ObserverPrimitiveCache.key("pairwise_kde", sid,
                coordinates=(data[0], data[1]), parameters=(role,),
                translation_invariant=state_id is not None)
            def compute():
                if cache is not None and cache.profile is not None:
                    with cache.profile.time("pairwise_kde"):
                        return backend.evaluate(data[0], data[1], config=bandwidth)
                return backend.evaluate(data[0], data[1], config=bandwidth)
            return cache.get_or_compute(key, compute, "pairwise_kde") if cache else compute()
        return evaluate
    try:
        with _deposition_histogram(method):
            EX.OLD._diag_kde = cached_diag
            extracted = EX._extract_all(screen, EXTENT, OBS_BINS)
    finally:
        EX.OLD._diag_kde = original_diag
        if owns_backend and hasattr(backend, "close"): backend.close()
    bank, family = RET._decoded_bank(extracted, received_state)
    if len(bank) != 45:
        raise RuntimeError(f"expected exactly 45 decoded WL channels, got {len(bank)}")
    return {"bank": bank, "family": family, "deposition_method": method.name}
