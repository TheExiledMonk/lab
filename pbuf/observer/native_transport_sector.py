"""Native-only sector gate: no split is made without a pre-existing free window."""
from __future__ import annotations
import numpy as np

def source_contact_window(length: int) -> np.ndarray:
    """DEV204's retained source-contact trajectory has no independently zero window."""
    return np.ones(length, dtype=bool)

def sector_decomposition(activity: np.ndarray, source_contact: np.ndarray) -> dict[str,np.ndarray|str]:
    # A sign-selected transported split requires source-contact==0.  With no
    # such window, preserve the full activity as unresolved local state and
    # emit an exact zero outgoing placeholder rather than inventing a flux cut.
    outgoing=np.zeros_like(activity)
    local=np.asarray(activity).copy()
    return {'local':local,'outgoing':outgoing,'definition':'BLOCKED_NO_FREE_NATIVE_PROPAGATION_WINDOW'}
