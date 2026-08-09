"""Frozen 45-channel decoder bank with injectable spatial deposition."""

from contextlib import contextmanager
import threading
import numpy as np

from pbuf.labs.foundation import native_full_received_state_information_retention001 as RET
from pbuf.labs.foundation import native_observable_extraction_method_sweep001 as EX
from .config import EXTENT, OBS_BINS
from .deposition import DepositionMethod, get_deposition_method


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
) -> dict:
    method = get_deposition_method(deposition_method)
    with _deposition_histogram(method):
        extracted = EX._extract_all(screen, EXTENT, OBS_BINS)
    bank, family = RET._decoded_bank(extracted, received_state)
    if len(bank) != 45:
        raise RuntimeError(f"expected exactly 45 decoded WL channels, got {len(bank)}")
    return {"bank": bank, "family": family, "deposition_method": method.name}
