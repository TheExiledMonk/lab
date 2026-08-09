"""Canonical, backend-ready weak-lensing prediction pipeline."""

from .pipeline import WLPipelineResult, compare_with_observations, make_backend, run_wl_pipeline
from .backends import CpuReferenceBackend, VulkanBackend

__all__ = ["WLPipelineResult", "compare_with_observations", "make_backend", "run_wl_pipeline",
           "CpuReferenceBackend", "VulkanBackend"]
