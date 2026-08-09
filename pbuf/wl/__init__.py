"""Canonical, backend-ready weak-lensing prediction pipeline."""

from .pipeline import WLPipelineResult, compare_with_observations, run_wl_pipeline

__all__ = ["WLPipelineResult", "compare_with_observations", "run_wl_pipeline"]
