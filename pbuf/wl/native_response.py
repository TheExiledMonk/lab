"""Frozen native fast/slow response and pair-transfer preparation."""

import numpy as np

from pbuf.labs.foundation import current_native_five_cluster_observable_benchmark001 as CUR


def build_native_response(rho3: np.ndarray) -> dict:
    m10, raw = CUR.current_native_m10(rho3)
    return {
        "rho3": np.asarray(rho3),
        "u_fast": None,
        "u_slow": None,
        "c_state": None,
        "pair_fast_coefficient_from_A8": raw["pair_fast_coefficient_from_A8"],
        "pair_slow_coefficient_from_A8": raw["pair_slow_coefficient_from_A8"],
        "terminal_common_history_relative_rms_error": raw["terminal_common_history_relative_rms_error"],
        "m10_vector": m10,
        "raw": raw,
    }
