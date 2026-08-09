"""Canonical end-to-end weak-lensing pipeline."""

from dataclasses import dataclass
import numpy as np

from pbuf.labs.foundation import native_full_state_2d_reconstruction_decoder_sweep001 as DEC
from .channels import decode_full_channel_bank
from .config import CHECKPOINT, PROPAGATION_STEP, PROPAGATION_STEPS
from .interface import get_interface_vector
from .launch import RayLaunch, launch_25pct, launch_100pct
from .los import project_interface_to_los
from .native_response import build_native_response
from .propagation import CpuReferenceBackend, PropagationBackend, PropagationConfig
from .received_state import build_received_state
from .reconstruction import build_reconstruction_candidates
from .screen import build_detector_screen
from .source import load_cluster_source


@dataclass(frozen=True)
class WLPipelineResult:
    cluster_id: str
    coverage_label: str
    source: dict
    native_response: dict
    los: dict
    launch: RayLaunch
    propagation: dict
    screen: dict
    received_state: dict
    channel_bank: dict[str, np.ndarray]
    channel_family: dict[str, str]
    reconstruction_candidates: dict[str, np.ndarray]
    reconstruction_meta: dict


def run_wl_pipeline(cluster: dict, coverage: str, backend: PropagationBackend | None = None) -> WLPipelineResult:
    if coverage not in ("25pct", "100pct"):
        raise ValueError(f"unsupported WL coverage: {coverage}")
    source = load_cluster_source(cluster)
    native = build_native_response(source["rho3"])
    los = project_interface_to_los(get_interface_vector(native))
    launch = launch_25pct() if coverage == "25pct" else launch_100pct()
    propagation = (backend or CpuReferenceBackend()).propagate(
        los["field"], launch, PropagationConfig(PROPAGATION_STEP, PROPAGATION_STEPS, CHECKPOINT)
    )
    screen = build_detector_screen(launch, propagation)
    received = build_received_state(launch, propagation, screen)
    decoded = decode_full_channel_bank(screen, received)
    candidates, meta = build_reconstruction_candidates(decoded["bank"], decoded["family"])
    return WLPipelineResult(cluster["id"], launch.coverage_label, source, native, los, launch,
                            propagation, screen, received, decoded["bank"], decoded["family"], candidates, meta)


def compare_with_observations(result: WLPipelineResult, data: dict) -> dict:
    targets = DEC._targets_after_decoding(data)
    return DEC._compare_candidates(result.reconstruction_candidates, targets)


def make_backend(name: str) -> PropagationBackend:
    """Construct a propagation backend without changing pipeline physics."""
    from .backends import make_backend as _make_backend
    return _make_backend(name)
