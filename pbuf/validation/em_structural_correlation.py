"""DEV199 Phase-B categorical comparison helpers (no numerical calibration)."""
from __future__ import annotations
import hashlib
from pathlib import Path


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def primary_references() -> list[dict[str, str]]:
    return [
        {"doi": "10.1103/PhysRevA.113.033732", "kind": "primary research",
         "claim": "Euler--Heisenberg nonlinear wave equation; phase modulation, birefringence, frequency mixing"},
        {"doi": "10.1103/PhysRevLett.119.250403", "kind": "primary research",
         "claim": "strong-field probe response depends on probe polarization"},
        {"doi": "10.1103/PhysRevD.113.085018", "kind": "primary research",
         "claim": "constant-background Heisenberg--Euler action is organized by electromagnetic Lorentz invariants"},
        {"doi": "10.1103/PhysRevD.100.036004", "kind": "primary research",
         "claim": "nonlinear interaction of electromagnetic waves in vacuum in the long-wavelength Euler--Heisenberg regime"},
    ]
