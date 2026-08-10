"""Propagation execution backends for the canonical WL pipeline."""

from .cpu import CpuReferenceBackend
from .vulkan import VulkanBackend
from .vulkan_runtime import (
    VulkanRuntime, VulkanUnavailableError, discover_vulkan_device,
    vulkan_available, vulkan_diagnostics,
)


def make_backend(name: str):
    if name == "cpu":
        return CpuReferenceBackend()
    if name == "vulkan":
        return VulkanBackend()
    raise ValueError(f"unsupported WL propagation backend: {name}")


__all__ = ["CpuReferenceBackend", "VulkanBackend", "VulkanRuntime",
           "VulkanUnavailableError", "discover_vulkan_device", "vulkan_available",
           "vulkan_diagnostics", "make_backend"]
