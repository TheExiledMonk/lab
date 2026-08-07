"""Compatibility entry point for the assistant-audited M13 Helmholtz core.

The scientific implementation lives in ``helmholtz_3d_verified``.
This wrapper preserves the historical import path used by existing labs
while ensuring callers receive the second-review implementation.
"""

from .helmholtz_3d_verified import *  # noqa: F401,F403
