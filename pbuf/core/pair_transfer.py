"""Compatibility entry point for the assistant-audited M08/M09/M10 core.

The scientific implementation lives in ``pair_transfer_verified``.
This wrapper preserves the historical import path used by existing labs
while ensuring all callers receive the second-review implementation.
"""
from .pair_transfer_verified import *  # noqa: F401,F403
