"""Minimal headless Vulkan runtime for canonical G3D compute propagation."""

from __future__ import annotations

import ctypes
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import threading

import numpy as np


_HERE = Path(__file__).resolve().parent
_BRIDGE_SOURCE = _HERE / "vulkan_bridge.c"
_SHADER_SOURCE = _HERE / "shaders" / "g3d_propagate.comp"
_BUILD_LOCK = threading.Lock()
_LIBRARY: ctypes.CDLL | None = None


class VulkanUnavailableError(RuntimeError):
    """Raised when the canonical float64 Vulkan compute path is unavailable."""


class _DeviceInfo(ctypes.Structure):
    _fields_ = [
        ("device_name", ctypes.c_char * 256),
        ("vendor_id", ctypes.c_uint32), ("device_id", ctypes.c_uint32),
        ("device_type", ctypes.c_uint32), ("api_version", ctypes.c_uint32),
        ("driver_version", ctypes.c_uint32), ("queue_family", ctypes.c_uint32),
        ("supports_float64", ctypes.c_uint32),
        ("max_invocations", ctypes.c_uint32), ("max_size_x", ctypes.c_uint32),
    ]


_DEVICE_TYPES = {
    0: "other", 1: "integrated_gpu", 2: "discrete_gpu", 3: "virtual_gpu", 4: "cpu",
}


def _version(value: int) -> str:
    return f"{value >> 22}.{(value >> 12) & 0x3ff}.{value & 0xfff}"


def _device_dict(info: _DeviceInfo) -> dict:
    return {
        "device_name": bytes(info.device_name).split(b"\0", 1)[0].decode(errors="replace"),
        "vendor_id": int(info.vendor_id), "device_id": int(info.device_id),
        "device_type": _DEVICE_TYPES.get(int(info.device_type), "unknown"),
        "api_version": _version(int(info.api_version)),
        "driver_version": _version(int(info.driver_version)),
        "compute_queue_family_index": int(info.queue_family),
        "supports_float64": bool(info.supports_float64),
        "max_compute_workgroup_invocations": int(info.max_invocations),
        "max_compute_workgroup_size_x": int(info.max_size_x),
    }


def _device_index() -> int:
    raw = os.environ.get("PBUF_VULKAN_DEVICE_INDEX")
    if raw is None:
        return -1
    try:
        value = int(raw)
    except ValueError as exc:
        raise VulkanUnavailableError(
            f"PBUF_VULKAN_DEVICE_INDEX must be an integer, got {raw!r}"
        ) from exc
    if value < 0:
        raise VulkanUnavailableError(f"PBUF_VULKAN_DEVICE_INDEX must be non-negative, got {value}")
    return value


def _build_paths() -> tuple[Path, Path]:
    digest = hashlib.sha256(_BRIDGE_SOURCE.read_bytes() + _SHADER_SOURCE.read_bytes()).hexdigest()[:16]
    root = Path(tempfile.gettempdir()) / "pbuf-vulkan" / digest
    return root / "libpbuf_vulkan.so", root / "g3d_propagate.spv"


def _run_build(command: list[str], dependency: str) -> None:
    try:
        result = subprocess.run(command, text=True, capture_output=True, check=False)
    except OSError as exc:
        raise VulkanUnavailableError(f"required Vulkan dependency {dependency} could not run: {exc}") from exc
    if result.returncode:
        detail = (result.stderr or result.stdout).strip()
        raise VulkanUnavailableError(f"{dependency} failed while preparing Vulkan backend: {detail}")


def _ensure_built() -> tuple[Path, Path]:
    library, spirv = _build_paths()
    with _BUILD_LOCK:
        library.parent.mkdir(parents=True, exist_ok=True)
        if not spirv.exists():
            compiler = shutil.which("glslangValidator")
            if compiler:
                _run_build([compiler, "-V", "-S", "comp", "-o", str(spirv), str(_SHADER_SOURCE)], "glslangValidator")
            else:
                compiler = shutil.which("glslc")
                if not compiler:
                    raise VulkanUnavailableError(
                        "Vulkan shader compiler missing; install glslangValidator or glslc (no downloads are performed at runtime)"
                    )
                _run_build([compiler, "-fshader-stage=compute", "-o", str(spirv), str(_SHADER_SOURCE)], "glslc")
        if not library.exists():
            cc = shutil.which("cc")
            if not cc:
                raise VulkanUnavailableError("C compiler missing; install cc and Vulkan development headers")
            _run_build([cc, "-O2", "-fPIC", "-shared", str(_BRIDGE_SOURCE), "-o", str(library), "-lvulkan"], "Vulkan C bridge build")
    return library, spirv


def _load_library() -> ctypes.CDLL:
    global _LIBRARY
    if _LIBRARY is None:
        library, _ = _ensure_built()
        try:
            lib = ctypes.CDLL(str(library))
        except OSError as exc:
            raise VulkanUnavailableError(f"Vulkan loader/bridge unavailable: {exc}") from exc
        lib.pbuf_vk_discover.argtypes = [ctypes.c_int, ctypes.POINTER(_DeviceInfo), ctypes.c_char_p]
        lib.pbuf_vk_discover.restype = ctypes.c_int
        lib.pbuf_vk_create.argtypes = [ctypes.c_char_p, ctypes.c_int, ctypes.c_uint32,
                                       ctypes.POINTER(_DeviceInfo), ctypes.c_char_p]
        lib.pbuf_vk_create.restype = ctypes.c_void_p
        lib.pbuf_vk_destroy.argtypes = [ctypes.c_void_p]
        lib.pbuf_vk_workgroup_size.argtypes = [ctypes.c_void_p]
        lib.pbuf_vk_workgroup_size.restype = ctypes.c_uint32
        _LIBRARY = lib
    return _LIBRARY


def discover_vulkan_device() -> dict:
    """Discover the deterministically selected compute+float64 Vulkan device."""
    lib = _load_library()
    info, error = _DeviceInfo(), ctypes.create_string_buffer(1024)
    if lib.pbuf_vk_discover(_device_index(), ctypes.byref(info), error):
        message = error.value.decode(errors="replace")
        if "float64" in message:
            raise RuntimeError("Vulkan backend requires shader float64 support for canonical parity")
        raise VulkanUnavailableError(message)
    return _device_dict(info)


def vulkan_diagnostics() -> dict:
    """Return availability diagnostics without raising for ordinary absence."""
    result = {"available": False, "shader_source": str(_SHADER_SOURCE)}
    try:
        library, spirv = _ensure_built()
        result.update({"bridge_library": str(library), "spirv": str(spirv),
                       "shader_compilation_pass": True, "device": discover_vulkan_device(),
                       "available": True})
    except Exception as exc:  # diagnostics API is deliberately non-raising
        result.update({"error": str(exc), "error_type": type(exc).__name__})
    return result


def vulkan_available() -> bool:
    return bool(vulkan_diagnostics()["available"])


class VulkanRuntime:
    """Persistent Vulkan instance/device/pipeline with per-run transfer buffers."""

    def __init__(self, workgroup_size: int = 256):
        lib = _load_library()
        _, spirv = _ensure_built()
        info, error = _DeviceInfo(), ctypes.create_string_buffer(1024)
        self._handle = lib.pbuf_vk_create(str(spirv).encode(), _device_index(), workgroup_size,
                                          ctypes.byref(info), error)
        if not self._handle:
            message = error.value.decode(errors="replace")
            if "float64" in message:
                raise RuntimeError("Vulkan backend requires shader float64 support for canonical parity")
            raise VulkanUnavailableError(message)
        self._lib = lib
        self.device = _device_dict(info)
        self.workgroup_size = int(lib.pbuf_vk_workgroup_size(self._handle))
        ptr = ctypes.POINTER(ctypes.c_double)
        lib.pbuf_vk_propagate.argtypes = [ctypes.c_void_p, ctypes.POINTER(ptr),
            ctypes.POINTER(ctypes.c_uint64), ctypes.POINTER(ptr), ctypes.c_uint32,
            ctypes.c_uint32, ctypes.c_uint32, ctypes.c_double, ctypes.c_uint32,
            ptr, ctypes.c_uint32, ctypes.c_char_p]
        lib.pbuf_vk_propagate.restype = ctypes.c_int

    def close(self) -> None:
        if getattr(self, "_handle", None):
            self._lib.pbuf_vk_destroy(self._handle)
            self._handle = None

    def __enter__(self): return self
    def __exit__(self, *_): self.close()
    def __del__(self): self.close()

    def propagate(self, arrays: list[np.ndarray], ray_count: int, nx: int, ny: int,
                  step: float, steps: int, checkpoints: tuple[int, ...]) -> list[np.ndarray]:
        arrays = [np.ascontiguousarray(a, dtype=np.float64) for a in arrays]
        outputs = [np.empty((len(checkpoints), ray_count), dtype=np.float64) for _ in range(6)]
        outputs.append(np.empty(ray_count, dtype=np.float64))
        cp = np.ascontiguousarray(checkpoints, dtype=np.float64)
        ptr = ctypes.POINTER(ctypes.c_double)
        in_ptrs = (ptr * len(arrays))(*(a.ctypes.data_as(ptr) for a in arrays))
        lengths = (ctypes.c_uint64 * len(arrays))(*(a.size for a in arrays))
        out_ptrs = (ptr * len(outputs))(*(a.ctypes.data_as(ptr) for a in outputs))
        error = ctypes.create_string_buffer(1024)
        rc = self._lib.pbuf_vk_propagate(
            self._handle, in_ptrs, lengths, out_ptrs, ray_count, nx, ny,
            step, steps, cp.ctypes.data_as(ptr), len(checkpoints), error,
        )
        if rc:
            raise RuntimeError(f"Vulkan propagation failed: {error.value.decode(errors='replace')}")
        return outputs
