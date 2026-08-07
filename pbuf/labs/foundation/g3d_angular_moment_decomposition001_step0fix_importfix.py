#!/usr/bin/env python3
"""Import-path correction runner for G3D angular moment decomposition step0fix.

This wrapper changes no science. It only ensures the repository root is on
sys.path before importing the existing step-0 fraction-gate correction module.
"""
from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))

import pbuf.labs.foundation.g3d_angular_moment_decomposition001_step0fix as FIX


if __name__ == "__main__":
    raise SystemExit(FIX.main())
