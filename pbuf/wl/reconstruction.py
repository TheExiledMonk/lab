"""Frozen PR #105 target-blind reconstruction inventory."""

import numpy as np
from pbuf.labs.foundation import native_full_state_2d_reconstruction_decoder_sweep001 as DEC


def build_reconstruction_candidates(bank: dict[str, np.ndarray], family: dict[str, str]) -> tuple[dict[str, np.ndarray], dict]:
    return DEC._build_candidates(bank, family)
