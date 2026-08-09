"""Persistent Vulkan context for the exact pairwise KDE shader."""
import ctypes, hashlib, shutil, subprocess, tempfile, threading, time
from pathlib import Path
import numpy as np
from .vulkan_runtime import (_DeviceInfo, _device_dict, _device_index,
                             VulkanUnavailableError)

_HERE=Path(__file__).resolve().parent; _C=_HERE/"vulkan_bridge.c"; _S=_HERE/"shaders/kde_pairwise.comp"
_LOCK=threading.Lock()

def _built():
    digest=hashlib.sha256(_C.read_bytes()+_S.read_bytes()).hexdigest()[:16]
    root=Path(tempfile.gettempdir())/"pbuf-vulkan-kde"/digest; lib=root/"libpbuf_kde.so"; spv=root/"kde.spv"
    with _LOCK:
        root.mkdir(parents=True,exist_ok=True)
        compiler=shutil.which("glslangValidator")
        if not compiler: raise VulkanUnavailableError("Vulkan KDE requires glslangValidator")
        if not spv.exists():
            p=subprocess.run([compiler,"-V","-S","comp","-o",str(spv),str(_S)],capture_output=True,text=True)
            if p.returncode: raise VulkanUnavailableError((p.stderr or p.stdout).strip())
        if not lib.exists():
            p=subprocess.run(["cc","-O2","-fPIC","-shared",str(_C),"-o",str(lib),"-lvulkan"],capture_output=True,text=True)
            if p.returncode: raise VulkanUnavailableError(p.stderr.strip())
    return lib,spv

class VulkanKDERuntime:
    def __init__(self,workgroup_size=256):
        libpath,spv=_built();self._lib=ctypes.CDLL(str(libpath));info=_DeviceInfo();err=ctypes.create_string_buffer(1024)
        self._lib.pbuf_vk_create.argtypes=[ctypes.c_char_p,ctypes.c_int,ctypes.c_uint32,ctypes.POINTER(_DeviceInfo),ctypes.c_char_p];self._lib.pbuf_vk_create.restype=ctypes.c_void_p
        self._handle=self._lib.pbuf_vk_create(str(spv).encode(),_device_index(),workgroup_size,ctypes.byref(info),err)
        if not self._handle: raise VulkanUnavailableError(err.value.decode())
        self.device=_device_dict(info);self.workgroup_size=workgroup_size
        ptr=ctypes.POINTER(ctypes.c_double);self._lib.pbuf_vk_kde.argtypes=[ctypes.c_void_p,ptr,ptr,ptr,ctypes.c_uint32,ptr,ctypes.c_char_p];self._lib.pbuf_vk_kde.restype=ctypes.c_int
        self._lib.pbuf_vk_destroy.argtypes=[ctypes.c_void_p]
    def evaluate(self,u,v,h):
        out=np.empty(u.size,dtype=np.float64);h=np.ascontiguousarray(h,dtype=np.float64);ptr=ctypes.POINTER(ctypes.c_double);err=ctypes.create_string_buffer(1024);started=time.perf_counter()
        rc=self._lib.pbuf_vk_kde(self._handle,u.ctypes.data_as(ptr),v.ctypes.data_as(ptr),out.ctypes.data_as(ptr),u.size,h.ctypes.data_as(ptr),err)
        elapsed=time.perf_counter()-started
        if rc: raise RuntimeError(err.value.decode())
        return out,{"warm_total_seconds":elapsed,"input_buffer_bytes":int(u.nbytes+v.nbytes+h.nbytes),"output_buffer_bytes":int(out.nbytes),"temporary_buffer_bytes":0,"estimated_total_gpu_bytes":int(u.nbytes+v.nbytes+h.nbytes+out.nbytes),"no_n_squared_buffer_allocation":True}
    def close(self):
        if getattr(self,"_handle",None):self._lib.pbuf_vk_destroy(self._handle);self._handle=None
    def __del__(self):self.close()
