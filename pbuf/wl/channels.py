"""Frozen 45-channel decoder bank."""

from pbuf.labs.foundation import native_full_received_state_information_retention001 as RET
from pbuf.labs.foundation import native_observable_extraction_method_sweep001 as EX
from .config import EXTENT, OBS_BINS


def decode_full_channel_bank(screen: dict, received_state: dict) -> dict:
    extracted = EX._extract_all(screen, EXTENT, OBS_BINS)
    bank, family = RET._decoded_bank(extracted, received_state)
    if len(bank) != 45:
        raise RuntimeError(f"expected exactly 45 decoded WL channels, got {len(bank)}")
    return {"bank": bank, "family": family}
