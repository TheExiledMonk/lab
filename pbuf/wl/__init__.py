"""Canonical, backend-ready weak-lensing prediction pipeline."""

from .pipeline import WLPipelineResult, compare_with_observations, make_backend, run_wl_pipeline
from .backends import CpuReferenceBackend, VulkanBackend
from .deposition import DepositionMethod, METHODS, get_deposition_method

__all__ = ["WLPipelineResult", "compare_with_observations", "make_backend", "run_wl_pipeline",
           "CpuReferenceBackend", "VulkanBackend", "DepositionMethod", "METHODS",
           "get_deposition_method"]
