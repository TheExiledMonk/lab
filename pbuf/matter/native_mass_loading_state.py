"""Dev145 native persistent-loading inventory and coefficient-free diagnostics."""
from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Any
import numpy as np

from pbuf.wl.native_incremental_elastic_energy import bounded_strain_energy, bounded_strain_stress

LOADING_NAMES = (
    "c_state", "source loading S", "accumulated response u", "|u|", "strain epsilon",
    "|epsilon|", "stress sigma(epsilon)", "bounded-strain energy W(epsilon)",
    "local response gradient", "local strain gradient", "integrated local deformation",
    "source-normalized loading", "response-normalized loading", "deformation energy above unloaded state",
    "central loading", "volume-integrated loading", "surface loading fingerprint",
    "far-response loading fingerprint", "compactness-like native ratio", "native loading/state composite",
)

_STATUS = (
    "NOT_LOADING", "SOURCE_DEPENDENT", "DERIVABLE_LOADING_PROXY", "DERIVABLE_LOADING_PROXY",
    "STRUCTURALLY_PLAUSIBLE", "DERIVABLE_LOADING_PROXY", "DERIVABLE_LOADING_PROXY",
    "DERIVABLE_LOADING_PROXY", "GEOMETRY_DEPENDENT", "GEOMETRY_DEPENDENT", "NONLOCAL",
    "DIMENSIONALLY_UNSUITABLE", "DIMENSIONALLY_UNSUITABLE", "DERIVABLE_LOADING_PROXY",
    "GEOMETRY_DEPENDENT", "DERIVABLE_LOADING_PROXY", "GEOMETRY_DEPENDENT", "NONLOCAL",
    "STRUCTURALLY_PLAUSIBLE", "UNRESOLVED",
)
_PERSISTENCE = (
    "PERSISTENT", "PERSISTENT", "PERSISTENT", "PERSISTENT", "PERSISTENT", "PERSISTENT",
    "PERSISTENT", "PERSISTENT", "MIXED", "MIXED", "PERSISTENT", "UNDETERMINED",
    "UNDETERMINED", "PERSISTENT", "PERSISTENT", "PERSISTENT", "PERSISTENT", "PERSISTENT",
    "PERSISTENT", "UNDETERMINED",
)

@dataclass(frozen=True)
class LoadingCandidate:
    id: str
    name: str
    status: str
    persistence: str
    local_or_integrated: str

def loading_inventory() -> list[dict[str, Any]]:
    rows=[]
    for i, (name, status, persistence) in enumerate(zip(LOADING_NAMES, _STATUS, _PERSISTENCE), 1):
        locality = "INTEGRATED" if i in (11, 16, 17, 18) else "LOCAL" if i not in (12, 13, 19, 20) else "COMPOSITE"
        rows.append(asdict(LoadingCandidate(f"L{i:02d}", name, status, persistence, locality)))
    return rows

def strain_loading_fraction(strain, epsilon_max: float = 1.0):
    """The sole immediate bounded, coefficient-free local loading proxy."""
    e=np.asarray(strain,dtype=float); em=float(epsilon_max)
    if not np.isfinite(em) or em <= 0 or np.any(~np.isfinite(e)) or np.any(np.abs(e)>=em):
        raise ValueError("requires finite abs(strain) < positive epsilon_max")
    return np.abs(e)/em

def loading_fingerprints(strain, epsilon_max: float = 1.0, K: float = 1.0, cell_volume=1.0):
    e=np.asarray(strain,dtype=float); ell=strain_loading_fraction(e,epsilon_max)
    vol=np.broadcast_to(np.asarray(cell_volume,dtype=float),e.shape)
    if np.any(vol <= 0): raise ValueError("cell_volume must be positive")
    return {"local_strain_fraction":ell, "peak_strain_fraction":float(ell.max(initial=0)),
            "integrated_abs_deformation":float(np.sum(np.abs(e)*vol)),
            "integrated_deformation_energy":float(np.sum(bounded_strain_energy(e,K,epsilon_max)*vol)),
            "stress":bounded_strain_stress(e,K,epsilon_max)}

def normalization_audit() -> dict[str, Any]:
    return {"strain_fraction":{"definition":"abs(epsilon)/epsilon_max","coefficient_free":True,"bounded":True},
            "stress":{"finite_native_maximum":False,"coefficient_free_fraction":False},
            "energy":{"finite_native_maximum":False,"coefficient_free_fraction":False,
                      "outcome":"NO_NATIVE_FINITE_W_NORMALIZATION"}}

def mass_loading_contract() -> dict[str, Any]:
    return {"contract":"PBUF_MASS_LOADING_V1","rest_loading_established":False,
            "rest_loading_definition":None,"persistent_loading":"SUPPORTED_PROXIES_NOT_MASS_IDENTIFICATION",
            "dimensionless_loading_available":True,"dimensionless_loading_definition":"abs(epsilon)/epsilon_max (local proxy)",
            "zero_loading_state_defined":True,"bounded_loading":True,
            "loading_geometry_dependence":"MIXED","loading_composite_behavior":"binding-dependent/geometry-dependent diagnostic"}
