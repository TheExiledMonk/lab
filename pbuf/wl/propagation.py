"""Propagation contract and compatibility exports for execution backends."""

from dataclasses import dataclass
from typing import Protocol
from .launch import RayLaunch


@dataclass(frozen=True)
class PropagationConfig:
    step: float
    steps: int
    checkpoint: object


class PropagationBackend(Protocol):
    def propagate(self, field: dict, launch: RayLaunch, config: PropagationConfig,
                  step_observer=None) -> dict: ...


from .backends.cpu import CpuReferenceBackend  # noqa: E402 - compatibility export
