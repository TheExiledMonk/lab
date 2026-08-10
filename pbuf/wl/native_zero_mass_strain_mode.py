"""Native strain packet extraction and scalar bookkeeping; no propagation law."""
from __future__ import annotations
from dataclasses import dataclass
import hashlib
import numpy as np
from pbuf.wl.native_incremental_elastic_energy import incremental_elastic_energy, integrate_packet


@dataclass(frozen=True)
class StrainPacket:
    packet_uid: str
    event_uid: str
    trajectory_uid: str
    mask: np.ndarray
    background: np.ndarray
    perturbation: np.ndarray
    method: str


def decompose_strain(total_strain, background_strain):
    total, bg = np.broadcast_arrays(np.asarray(total_strain, float), np.asarray(background_strain, float))
    if np.any(~np.isfinite(total)) or np.any(~np.isfinite(bg)):
        raise ValueError("strain states must be finite")
    return bg.copy(), total - bg


def extract_packet(perturbation, background=0.0, *, method="relative_threshold", threshold=0.05,
                   event_uid="event-0", trajectory_uid="trajectory-0"):
    de, bg = np.broadcast_arrays(np.asarray(perturbation, float), np.asarray(background, float))
    mag = np.abs(de)
    if method in ("relative_threshold", "PK01"):
        mask = mag >= float(threshold) * mag.max() if mag.max() else np.zeros(de.shape, bool)
    elif method in ("fixed_percentile", "PK03"):
        mask = mag >= np.percentile(mag, float(threshold))
    elif method in ("energy_fraction", "PK04"):
        fraction = float(threshold); order = np.argsort(mag.ravel())[::-1]
        weights = mag.ravel()[order] ** 2; n = np.searchsorted(np.cumsum(weights), fraction * weights.sum()) + 1
        mask = np.zeros(de.size, bool); mask[order[:n]] = True; mask = mask.reshape(de.shape)
    else:
        raise ValueError("unsupported deterministic packet method")
    digest = hashlib.sha256(np.packbits(mask).tobytes() + event_uid.encode() + trajectory_uid.encode()).hexdigest()[:16]
    return StrainPacket(f"packet-{digest}", event_uid, trajectory_uid, mask, bg.copy(), de.copy(), method)


def scalar_candidates(packet: StrainPacket, *, displacement=None, response_gradient=None,
                      cell_volume=1.0, K=1.0, epsilon_max=1.0):
    de, m = packet.perturbation, packet.mask; x = de[m]
    dw = incremental_elastic_energy(packet.background, de, K, epsilon_max)
    integ = integrate_packet(dw, m, cell_volume)
    vol = float(np.sum(np.broadcast_to(np.asarray(cell_volume, float), de.shape)[m]))
    out = {"A01": float(np.max(np.abs(x))) if x.size else 0.0,
           "A02": float(np.sqrt(np.mean(x*x))) if x.size else 0.0,
           "A03": float(np.sum(np.abs(x))), "A04": float(np.sqrt(np.sum(x*x))),
           "A08": float(np.sum(np.abs(x))), "A09": float(np.sum(x*x)),
           "A10": integ["signed"], "A11": integ["positive"], "A12": integ["signed"],
           "A13": vol, "A14": _covariance_volume(m),
           "A15": (float(np.max(np.abs(x))) if x.size else 0.0) * vol,
           "A16": (float(np.max(np.abs(x))) if x.size else 0.0) ** 2 * vol}
    u = de if displacement is None else np.asarray(displacement, float)
    out.update(A05=float(np.max(np.abs(u[m]))) if x.size else 0.0,
               A06=float(np.sqrt(np.mean(u[m]**2))) if x.size else 0.0)
    grad = np.gradient(de) if response_gradient is None else response_gradient
    if isinstance(grad, (list, tuple)): g2 = sum(np.asarray(g)**2 for g in grad)
    else: g2 = np.asarray(grad, float)**2
    out["A07"] = float(np.sqrt(np.sum(g2[m])))
    return out


def _covariance_volume(mask):
    coords = np.argwhere(mask)
    if len(coords) < 2: return float(len(coords))
    cov = np.atleast_2d(np.cov(coords, rowvar=False)); return float(np.sqrt(max(np.linalg.det(cov), 0.0)))


def packet_sample(packet, path_position_native, propagation_progression_index, candidate_scalar_state):
    return {"packet_uid": packet.packet_uid, "event_uid": packet.event_uid,
            "trajectory_uid": packet.trajectory_uid, "path_position_native": float(path_position_native),
            "propagation_progression_index": int(propagation_progression_index),
            "local_background_state": packet.background.copy(), "perturbation_state": packet.perturbation.copy(),
            "candidate_scalar_state": float(candidate_scalar_state)}
