"""Wave-state contracts and non-circular synthetic ratio controls for Dev140."""
from __future__ import annotations
from dataclasses import asdict, dataclass


@dataclass(frozen=True)
class NativeWaveState:
    wavelength_native: float | None = None
    frequency_native: float | None = None
    provenance: str = ""
    classification: str = "NATIVE_DIMENSIONLESS_WAVE_STATE"
    derived_from_L0: bool = False
    derived_from_T0: bool = False

    def to_dict(self): return asdict(self)


def redshift_from_wavelength(emitted, observed):
    if emitted <= 0 or observed <= 0: raise ValueError("wavelengths must be positive")
    return float(observed/emitted - 1.0)


def redshift_from_frequency(emitted, observed):
    if emitted <= 0 or observed <= 0: raise ValueError("frequencies must be positive")
    return float(emitted/observed - 1.0)


def frequency_closure(f_native, f_physical, *, derived_from_T0=False):
    if derived_from_T0: return {"status": "CIRCULAR", "T0": None}
    if f_native <= 0 or f_physical <= 0: raise ValueError("frequencies must be positive")
    return {"status": "ESTABLISHED", "T0": float(f_native/f_physical)}


def wavelength_closure(wavelength_native, wavelength_physical, *, derived_from_L0=False):
    if derived_from_L0: return {"status": "CIRCULAR", "L0": None}
    if wavelength_native <= 0 or wavelength_physical <= 0: raise ValueError("wavelengths must be positive")
    return {"status": "ESTABLISHED", "L0": float(wavelength_physical/wavelength_native)}


def triad_residual(speed_native, frequency_native, wavelength_native):
    return float(speed_native-frequency_native*wavelength_native)


def synthetic_wave_controls():
    rows=[]
    for z in (0,.01,.05,.1,.25,.5,1,2,4):
        for lam in (.25,.5,1,2,4):
            rows.append({"z":z,"wavelength_emit_native":lam,"wavelength_observed_native":lam*(1+z),
                         "frequency_emit_native":1/lam,"frequency_observed_native":1/(lam*(1+z))})
    return rows


def current_wave_inventory():
    return [
      {"name":"native_source_depth.D07 local variable phase","classification":"NUMERICAL_FREQUENCY_ONLY","usable":False,"reason":"covariance determinant label, not physical optical phase"},
      {"name":"trajectory path length","classification":"NOT_PHYSICALLY_USABLE","usable":False,"reason":"phase cannot be manufactured from path length"},
      {"name":"frequency/wavelength/period/wavenumber state","classification":"NOT_PHYSICALLY_USABLE","usable":False,"reason":"absent from frozen optical transport state"},
    ]
