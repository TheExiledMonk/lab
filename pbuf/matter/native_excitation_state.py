"""Dev146 native excitation candidates, without an energy identification."""
from __future__ import annotations
from dataclasses import asdict, dataclass, field
from typing import Any
import numpy as np

EXCITATION_NAMES=("Dev144 transported scalar q","signed excitation state","excitation magnitude + sign",
 "local excitation vector","transverse excitation vector","neighbor-link excitation","excitation density",
 "packet-integrated excitation","excitation norm","excitation-gradient state","fast-channel excitation",
 "slow-channel excitation","coupled fast/slow excitation","localized standing component","propagating component",
 "localized + propagating pair","conserved excitation norm","conserved pair-state excitation",
 "loading-bound excitation component","minimal additional excitation state required")

@dataclass(frozen=True)
class ExcitationCandidate:
    id: str
    name: str
    status: str
    source_initializable: bool
    conserved: str
    energy_like: bool = False

def excitation_registry() -> list[dict[str,Any]]:
    rows=[]
    for i,name in enumerate(EXCITATION_NAMES,1):
        if i==1: status,initial,conserved="ESTABLISHED_NEUTRAL_STATE",True,"IDENTITY_ONLY"
        elif i in (2,3,4,5,6,7,8,9,10,11,12,13): status,initial,conserved="STRUCTURALLY_DEFINABLE",True,"NOT_DERIVED"
        elif i in (14,15,16,19): status,initial,conserved="MISSING_DYNAMIC_DECOMPOSITION",False,"NOT_DERIVED"
        elif i in (17,18): status,initial,conserved="MISSING_NATIVE_NORM",False,"NOT_DERIVED"
        else: status,initial,conserved="REQUIRED_ADDITIONAL_DOF",False,"NOT_DERIVED"
        rows.append(asdict(ExcitationCandidate(f"X{i:02d}",name,status,initial,conserved)))
    return rows

@dataclass
class NativeExcitationState:
    q_source: float
    q_state: float | None = None
    spatial_index: int = 0
    history: list[dict[str,float]] = field(default_factory=list)
    def __post_init__(self):
        self.q_source=float(self.q_source); self.q_state=self.q_source if self.q_state is None else float(self.q_state)
        if not np.isfinite(self.q_source) or not np.isfinite(self.q_state): raise ValueError("excitation must be finite")
    def identity_progress(self, spatial_position: float):
        self.spatial_index += 1
        self.history.append({"spatial_index":self.spatial_index,"spatial_position":float(spatial_position),"q_state":self.q_state})

def excitation_requirements() -> dict[str,Any]:
    return {"SOURCE_EXCITATION_INITIAL_STATE":True,"ZERO_LOADING_PROPAGATION_AT_C":True,
            "EXCITATION_AMOUNT_CHANGES_C":False,"VACUUM_DISSIPATION":False,
            "TRAJECTORY_COUPLING_UNCHANGED":True,"Q_IDENTIFIED_AS_ENERGY":False}

def energy_like_classification() -> dict[str,Any]:
    gates={"source_supplied":True,"conserved_under_unloaded_identity_transport":True,"transferable":False,
           "additive_composable":False,"nonnegative":False,"consistent_exchange":False}
    return {"gates":gates,"EXCITATION_ENERGY_LIKE_STATE_ESTABLISHED":False,
            "classification":"PBUF_EXCITATION_PHYSICAL_DEFINITION_UNRESOLVED"}

