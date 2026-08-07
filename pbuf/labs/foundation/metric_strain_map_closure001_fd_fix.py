#!/usr/bin/env python3
"""Numerically conditioned rerun of METRIC-STRAIN MAP CLOSURE 001.

This wrapper changes no physics. The original closure lab failed only because
its finite-difference probe perturbed an O(1) background metric by O(1e-11),
so subtractive cancellation dominated a test of an exactly affine map.

The corrected probe uses an O(1e-3) symmetric strain increment and evaluates
an exact centered secant of the same map. All other algebra gates, equations,
interpretation, and generated bridge specifications come from the original
lab unchanged.
"""
from __future__ import annotations

from pathlib import Path
import numpy as np

from pbuf.labs.foundation import metric_strain_map_closure001 as base

ROOT = Path(__file__).resolve().parents[3]
base.OUT = ROOT / "runs" / "metric_strain_map_closure001_fd_fix"

_original_algebra_audit = base.algebra_audit


def _conditioned_algebra_audit() -> dict:
    audit = _original_algebra_audit()

    # Re-test only the numerically ill-conditioned finite-difference gate.
    # The map itself is exactly affine: g(gbar,chi)=gbar+2 chi.
    rng = np.random.default_rng(20260808 + 1)
    eta = np.diag([-1.0, 1.0, 1.0, 1.0])
    chi = base._sym(rng.normal(scale=1.0e-3, size=(4, 4)))
    dchi = base._sym(rng.normal(scale=1.0e-3, size=(4, 4)))

    # Centered secant with unit parameter displacement. This keeps the
    # metric perturbation ~1e-3 instead of ~1e-11 and therefore tests the
    # affine response rather than floating-point cancellation.
    dg_fd = (
        base.metric_from_strain(eta, chi + dchi)
        - base.metric_from_strain(eta, chi - dchi)
    ) / 2.0
    dg_exact = 2.0 * dchi
    err = base._relative_rms(dg_fd, dg_exact)

    audit["metrics"]["finite_difference_response_relative_rms_error"] = err
    audit["checks"]["finite_difference_response_pass"] = bool(err <= base.TOL)
    audit["all_checks_pass"] = bool(all(audit["checks"].values()))
    audit["finite_difference_conditioning_fix"] = {
        "physics_changed": False,
        "map_changed": False,
        "tolerance_changed": False,
        "original_failure_cause": "subtractive cancellation from O(1e-11) metric perturbation on O(1) background",
        "corrected_probe": "centered secant with O(1e-3) symmetric strain increment",
    }
    return audit


base.algebra_audit = _conditioned_algebra_audit


if __name__ == "__main__":
    raise SystemExit(base.main())
