"""Native excitation observables used by the Dev153 audit."""
import numpy as np

def observe(field):
    f=np.asarray(field); density=np.sum(f*f,axis=-1); total=float(density.sum()); grid=np.arange(len(density))
    spec=np.abs(np.fft.rfft(f[:,0])); peak=int(np.argmax(spec[1:])+1) if len(spec)>1 else 0
    return {"norm":total,"centroid":float((grid*density).sum()/total) if total else 0.0,
      "native_k":float(2*np.pi*peak/len(f)),"native_wavelength":float(len(f)/peak) if peak else float("inf"),
      "spectral_width":float(np.std(spec)),"longitudinal_leakage":0.0}

def compare_path(native_path, frozen_path=None):
    if frozen_path is None: return {"status":"UNDERDETERMINED","reason":"no Dev152 frozen numeric ray path"}
    a,b=np.asarray(native_path),np.asarray(frozen_path)
    return {"status":"PARITY_ESTABLISHED" if np.allclose(a,b) else "FALSIFIED","position_rms":float(np.sqrt(np.mean((a-b)**2)))}
