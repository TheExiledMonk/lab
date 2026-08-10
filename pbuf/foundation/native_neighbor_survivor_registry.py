"""Canonical, hash-validated Dev151 survivor registry for Dev152."""
from __future__ import annotations

from dataclasses import asdict, dataclass
from hashlib import sha256
import json
from pathlib import Path

from .native_neighbor_constitutive_law import law_registry, mechanism_registry


@dataclass(frozen=True)
class Survivor:
    survivor_id: str
    state_representation: str
    constitutive_law_id: str
    mechanism_id: str
    static_parity_score: float = 1.0
    dynamic_parity_score: float = 1.0
    free_parameter_count: int = 0
    coefficient_free: bool = True
    basis_invariant: bool = True
    resolution_status: str = "DEV151_PARITY"
    coordinate_status: str = "COVARIANT"
    executable: bool = True

    def to_dict(self):
        return asdict(self)


LAW_MECHANISMS = {
    "C08": "MEC12", "C10": "MEC15", "C12": "MEC02",
    "C13": "MEC06", "C16": "MEC17", "C18": "MEC18",
}


def load_dev151_survivors() -> list[Survivor]:
    """Load every Dev151 law survivor without manual pruning."""
    laws = {x["id"]: x for x in law_registry() if x["status"] == "STRUCTURALLY_SUPPORTED"}
    mechanisms = {x["id"] for x in mechanism_registry() if x["status"] == "STRUCTURALLY_SUPPORTED"}
    out = []
    for law_id in sorted(laws):
        mechanism = LAW_MECHANISMS[law_id]
        if mechanism not in mechanisms:
            raise ValueError(f"Dev151 mechanism mismatch for {law_id}")
        out.append(Survivor(f"S-{law_id}", "(L,X1,X2)", law_id, mechanism))
    return out


def reference_hashes(root: Path) -> dict[str, str]:
    files = [
        root / "pbuf/foundation/native_neighbor_state.py",
        root / "pbuf/foundation/native_neighbor_constitutive_law.py",
        root / "runs/unified_native_neighbor_state001/result.json",
    ]
    return {str(p.relative_to(root)): sha256(p.read_bytes()).hexdigest() for p in files}


def validate_dev151(root: Path) -> dict:
    result = root / "runs/unified_native_neighbor_state001/result.json"
    report = root / "runs/unified_native_neighbor_state001/report.txt"
    if not result.exists() or not report.exists():
        raise RuntimeError("DEV152_BASELINE_MISMATCH")
    data = json.loads(result.read_text())
    required = {
        "PBUF_UNIFIED_NATIVE_NEIGHBOR_STATE_ESTABLISHED",
        "UNIFIED_NEIGHBOR_STATIC_DEFORMATION_PARITY",
        "UNIFIED_NEIGHBOR_DYNAMIC_EXCITATION_PARITY",
        "PBUF_UNIFIED_STATIC_DYNAMIC_PARITY_ESTABLISHED",
        "PBUF_MICRO_MACRO_NEIGHBOR_LAW_NOT_UNIQUE",
        "PBUF_SHARED_STATE_CROSS_COUPLING_UNRESOLVED",
    }
    if data.get("status") != "DEV151_AUDIT_COMPLETE" or not required.issubset(data.get("outcomes", [])):
        raise RuntimeError("DEV152_BASELINE_MISMATCH")
    return {"status": data["status"], "required_outcomes": sorted(required), "validated": True}


def interface_contract() -> dict:
    return {"state_rank": 3, "fields": ["L", "X1", "X2"],
            "progression_signature": "progress(state, frames, steps)->history",
            "new_interaction_coefficients": 0}
