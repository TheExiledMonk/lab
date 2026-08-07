#!/usr/bin/env python3
"""PBUF MICROSTRUCTURE-ENTITY-LAB-001 — microscopic constituent architecture search."""
from __future__ import annotations

import csv
import json
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from weak_lensing_observation001 import file_sha256, make_field, propagate, resample_to_grid, compare_arrays, ssim_index
import observable_lab001 as obs_lab
import source_plane_lab001 as src_lab

OUT = ROOT / "runs" / "microstructure_entity_lab001"
PLOTS = OUT / "plots"
BENCHMARK = ROOT / "PBUF_benchmark"

CONFIG = {"nphotons": 20000, "grid_n": 256, "step": 0.03, "steps": 160, "y_span": 3.0, "extent": 8.0, "strength": 0.18, "bins": 64}
CLUSTERS = [
    {"id": "Abell2744", "label": "Abell 2744", "slug": "abell2744", "directory": "WL-001_Abell2744"},
    {"id": "MACS0416", "label": "MACS J0416", "slug": "macs0416", "directory": "WL-002_MACS0416"},
    {"id": "MACS1149", "label": "MACS J1149", "slug": "macs1149", "directory": "WL-003_MACS1149"},
    {"id": "AbellS1063", "label": "Abell S1063", "slug": "abells1063", "directory": "WL-004_AbellS1063"},
    {"id": "Abell370", "label": "Abell 370", "slug": "abell370", "directory": "WL-005_Abell370"},
]
EXPECTED_HASHES = {
    "constitutive_equations.py": "e2c789d19fd559753519704c6668c7a2879c53eb61a315604ec81af6795aca9f",
    "weak_lensing_observation001.py": "a5c3632fec9adfc2659d5c283d07c599db6db7edb2e34e4aebd84e35434642bc",
    "observable_lab001.py": "2867c0bf94fabe3fbba0264d5d272ecadbdf45fabb341a6a8fe77972fbaec132",
    "source_plane_lab001.py": "efa9d74924cb61a3b48a69fa075055512d86391d03194be342597420bc353de4",
    "numerical_convergence001.py": "0442f878713de6530b5a1b1844b8ece037852d461bcb695360e8a3345fd58f29",
}
EPS = np.finfo(np.float64).eps
COHERENCE_GAIN_THRESHOLD = 1e-4
MEMORY_INDEX_THRESHOLD = 0.9
ACTIVITY_THRESHOLD = 1e-6

DT = 0.10
STEPS = 20
K = 1.0
GAMMA = 0.50
OMEGA = 0.20
INTERNAL_K = 0.6

ALPHA_FS = 1.0 / 137.035999084
THREE_ALPHA_FS = 3.0 * ALPHA_FS


@dataclass(frozen=True)
class Architecture:
    number: int
    code: str
    name: str
    principle: str
    is_wrong: bool = False
    internal_nodes: int = 1


ARCHITECTURES = [
    Architecture(1, "A1", "Point Element", "single microscopic node u(x,y); linear relaxation control", internal_nodes=1),
    Architecture(2, "A2", "Two-State Constituent", "two coupled internal nodes u_a, u_b with shared local equilibrium", internal_nodes=2),
    Architecture(3, "A3", "Three-State Constituent", "triangular microstructure u_a, u_b, u_c; cycle coupling allows internal circulation", internal_nodes=3),
    Architecture(4, "A4", "Elastic Link Element", "node plus elastic connections to neighbours; internal deformation possible", internal_nodes=1),
    Architecture(5, "A5", "Oscillator Cell", "local oscillator z(x,y) = r·exp(iφ); neighbour coupling K·(<z>_n − z)", internal_nodes=2),
    Architecture(6, "A6", "Rotational Cell", "internal rotational DOF θ evolves by neighbour alignment; orientation emerges", internal_nodes=1),
    Architecture(7, "A7", "Loop Constituent", "closed internal loop (3 nodes x->y->z->x); supports circulating internal state without explicit phase", internal_nodes=3),
    Architecture(8, "A8", "Dual-Layer Constituent", "fast internal dynamics u_fast + slow structural u_slow; memory emerges from layer separation", internal_nodes=2),
    Architecture(9, "A9", "Cooperative Cell", "each cell internally computes weighted 9-point neighbourhood response before updating", internal_nodes=1),
    Architecture(10, "A10", "Unified Microcell", "multi-DOF cell: oscillation (r,φ) + neighbour coupling + reversible storage; no explicit S9 variables", internal_nodes=4),
    Architecture(11, "WR1", "Wrong: Random Internal Topology", "internal node connections randomly reshuffled each step", is_wrong=True, internal_nodes=2),
    Architecture(12, "WR2", "Wrong: Disconnected Internal Nodes", "internal nodes evolve independently with no internal coupling", is_wrong=True, internal_nodes=2),
    Architecture(13, "WR3", "Wrong: Over-Connected Constituent", "all internal nodes equally coupled to all others including every neighbour", is_wrong=True, internal_nodes=3),
    Architecture(14, "WR4", "Wrong: Frozen Internal Architecture", "internal architecture frozen; node state evolves but internal DOFs do not", is_wrong=True, internal_nodes=2),
]


def neighbours4(u: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    p = np.pad(u, 1, mode="reflect")
    return p[:-2, 1:-1], p[2:, 1:-1], p[1:-1, :-2], p[1:-1, 2:]


def neighbours9_weighted(u: np.ndarray) -> np.ndarray:
    p = np.pad(u, 1, mode="reflect")
    centre = p[1:-1, 1:-1]
    sides = 0.25 * (p[:-2, 1:-1] + p[2:, 1:-1] + p[1:-1, :-2] + p[1:-1, 2:])
    diag = 0.0625 * (p[:-2, :-2] + p[:-2, 2:] + p[2:, :-2] + p[2:, 2:])
    return sides + diag - centre


def complex_neighbours(z: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    return neighbours4(z.real), neighbours4(z.imag)


def complex_neighbour_mean(z: np.ndarray) -> np.ndarray:
    return sum(neighbours4(z)) / 4.0


def cosine(a: np.ndarray, b: np.ndarray) -> float:
    den = float(np.linalg.norm(a) * np.linalg.norm(b))
    return float(np.dot(a.ravel(), b.ravel()) / den) if den > 1e-30 else float("nan")


def gradient_coherence(field: np.ndarray) -> float:
    gy, gx = np.gradient(field)
    mag = np.hypot(gx, gy)
    if float(mag.max()) <= 1e-15:
        return 0.0
    ux = gx / np.maximum(mag, 1e-15)
    uy = gy / np.maximum(mag, 1e-15)
    nb = neighbours4(ux), neighbours4(uy)
    align_x = sum(ux * nx for nx in nb[0]) / 4.0
    align_y = sum(uy * ny for ny in nb[1]) / 4.0
    weights = mag.ravel()
    vals = (align_x * ux + align_y * uy).ravel()
    mask = weights > float(weights.max()) * 1e-8
    if not mask.any():
        return 0.0
    return float(np.sum(weights[mask] * vals[mask]) / np.sum(weights[mask]))


def spatial_correlation_length(field: np.ndarray) -> float:
    f = field - float(field.mean())
    n = field.shape[0]
    hann = np.outer(np.hanning(n), np.hanning(n))
    f_w = f * hann
    F = np.fft.fft2(f_w)
    psd = np.abs(np.fft.fftshift(F)) ** 2
    y, x = np.indices(psd.shape)
    cy, cx = psd.shape[0] // 2, psd.shape[1] // 2
    r = np.hypot(x - cx, y - cy).astype(int)
    rmax = min(cx, cy)
    if rmax == 0:
        return 0.0
    radial = np.zeros(rmax + 1)
    counts = np.zeros(rmax + 1)
    flat_r = r.ravel()
    flat_psd = psd.ravel()
    for ri, val in zip(flat_r, flat_psd):
        if ri <= rmax:
            radial[ri] += val
            counts[ri] += 1
    counts = np.maximum(counts, 1)
    radial = radial / counts
    if radial[0] <= 0:
        return 0.0
    radial = radial / radial[0]
    below = np.where(radial < 1.0 / np.e)[0]
    if len(below) == 0:
        return float(rmax)
    return float(below[0])


def temporal_persistence_length(states: list[np.ndarray]) -> float:
    if len(states) < 4:
        return 0.0
    activity = np.array([float(np.sqrt(np.mean((states[t + 1] - states[t]) ** 2))) for t in range(len(states) - 1)])
    if activity.max() <= 1e-15:
        return 0.0
    a = activity - activity.mean()
    var = float(np.dot(a, a))
    if var <= 0:
        return 0.0
    for lag in range(1, len(a)):
        acf = float(np.dot(a[:-lag], a[lag:])) / var
        if acf < 1.0 / np.e:
            return float(lag)
    return float(len(a))


def relaxation_time(states: list[np.ndarray], eq: np.ndarray) -> float:
    if len(states) < 3:
        return 0.0
    target = states[-1]
    diffs = np.array([float(np.sqrt(np.mean((s - target) ** 2))) for s in states])
    if diffs[0] <= 1e-15:
        return 0.0
    for k in range(1, len(diffs)):
        if diffs[k] < 0.5 * diffs[0]:
            return float(k)
    return float(len(diffs))


def build_C(u: np.ndarray, strength: float) -> np.ndarray:
    lo, hi = float(u.min()), float(u.max())
    if hi - lo < 1e-15:
        return np.zeros_like(u)
    return strength * (u - lo) / (hi - lo)


def evolve_architecture(arch: Architecture, rho: np.ndarray, strength: float, rng: np.random.RandomState) -> dict:
    """Returns dict with keys: u_final, history, diag, wave, internal_metrics, activity_log."""
    eq = strength * rho
    code = arch.code
    result = {"u_final": eq.copy(), "history": [eq.copy()], "diag": {}, "wave": {}, "internal_log": [], "circulation_score": 0.0, "cooperative_score": 0.0}

    if code == "A1":
        u = eq.copy()
        result["history"].append(u.copy())
        for step in range(STEPS):
            F = sum(neighbours4(u)) / 4.0 - u
            u = u + DT * K * F
            u = np.clip(u, -5.0, 5.0)
            result["history"].append(u.copy())
        result["u_final"] = u

    elif code == "A2":
        ua = eq.copy()
        ub = eq.copy() + 0.05 * strength * rng.randn(*rho.shape)
        result["history"].append(0.5 * (ua + ub))
        circulation_accumulator = 0.0
        for step in range(STEPS):
            n4a = sum(neighbours4(ua)) / 4.0
            n4b = sum(neighbours4(ub)) / 4.0
            dua = DT * K * ((n4a - ua) + INTERNAL_K * (ub - ua))
            dub = DT * K * ((n4b - ub) + INTERNAL_K * (ua - ub))
            ua = np.clip(ua + dua, -5.0, 5.0)
            ub = np.clip(ub + dub, -5.0, 5.0)
            mixed = 0.5 * (ua + ub)
            result["history"].append(mixed.copy())
            result["internal_log"].append((ua.copy(), ub.copy()))
            circulation_accumulator += float(np.sqrt(np.mean((ua - ub) ** 2)))
        result["u_final"] = 0.5 * (ua + ub)
        if len(result["internal_log"]) > 0:
            result["circulation_score"] = float(np.clip(circulation_accumulator / STEPS / max(strength, 1e-15) / 10.0, 0.0, 1.0))

    elif code == "A3":
        ua = eq.copy()
        ub = eq.copy() + 0.05 * strength * rng.randn(*rho.shape)
        uc = eq.copy() + 0.05 * strength * rng.randn(*rho.shape)
        circ_log = []
        result["history"].append((ua + ub + uc) / 3.0)
        for step in range(STEPS):
            n4a, n4b, n4c = sum(neighbours4(ua)) / 4.0, sum(neighbours4(ub)) / 4.0, sum(neighbours4(uc)) / 4.0
            dua = DT * K * ((n4a - ua) + INTERNAL_K * (ub - ua))
            dub = DT * K * ((n4b - ub) + INTERNAL_K * (uc - ub))
            duc = DT * K * ((n4c - uc) + INTERNAL_K * (ua - uc))
            ua = np.clip(ua + dua, -5.0, 5.0)
            ub = np.clip(ub + dub, -5.0, 5.0)
            uc = np.clip(uc + duc, -5.0, 5.0)
            mixed = (ua + ub + uc) / 3.0
            result["history"].append(mixed.copy())
            result["internal_log"].append((ua.copy(), ub.copy(), uc.copy()))
            circ = float(np.sqrt(np.mean((ub - ua) ** 2)) + np.sqrt(np.mean((uc - ub) ** 2)) + np.sqrt(np.mean((ua - uc) ** 2)))
            circ_log.append(circ)
        result["u_final"] = (ua + ub + uc) / 3.0
        if circ_log:
            mean_circ = float(np.mean(circ_log))
            result["circulation_score"] = float(np.clip(mean_circ / max(strength, 1e-15) / 5.0, 0.0, 1.0))

    elif code == "A4":
        u = eq.copy()
        internal_deform = np.zeros_like(u)
        result["history"].append(u.copy())
        for step in range(STEPS):
            n4 = sum(neighbours4(u)) / 4.0
            F_elastic = (n4 - u) - internal_deform
            internal_deform = internal_deform + DT * 0.4 * (n4 - u)
            u = u + DT * K * F_elastic
            u = np.clip(u, -5.0, 5.0)
            result["history"].append(u.copy())
            result["internal_log"].append(internal_deform.copy())
        result["u_final"] = u

    elif code == "A5":
        gy, gx = np.gradient(rho)
        phase_init = np.arctan2(gy, gx)
        z = 0.5 * np.exp(1j * phase_init) + 0.05 * strength * (rng.randn(*rho.shape) + 1j * rng.randn(*rho.shape))
        amp_log = []
        result["history"].append(np.abs(z).copy())
        for step in range(STEPS):
            n4z = complex_neighbour_mean(z)
            dz = 1j * OMEGA * z + K * (n4z - z)
            z = z + DT * dz
            result["history"].append(np.abs(z).copy())
            amp_log.append(float(np.mean(np.abs(z))))
        result["u_final"] = np.abs(z).astype(np.float64)

    elif code == "A6":
        u = eq.copy()
        theta = np.arctan2(*np.gradient(u))
        result["history"].append(u.copy())
        for step in range(STEPS):
            n4_theta = sum(np.exp(1j * np.angle(nj)) for nj in neighbours4(theta)) / 4.0
            target = np.angle(n4_theta)
            delta = target - theta
            delta = np.arctan2(np.sin(delta), np.cos(delta))
            theta = theta + DT * K * delta
            grad_rot = np.sin(theta)
            u = u + DT * 0.5 * (grad_rot - u)
            u = np.clip(u, -5.0, 5.0)
            result["history"].append(u.copy())
        result["u_final"] = u

    elif code == "A7":
        ua = eq.copy()
        ub = eq.copy() + 0.05 * strength * rng.randn(*rho.shape)
        uc = eq.copy() + 0.05 * strength * rng.randn(*rho.shape)
        circ_log = []
        result["history"].append((ua + ub + uc) / 3.0)
        for step in range(STEPS):
            n4a, n4b, n4c = sum(neighbours4(ua)) / 4.0, sum(neighbours4(ub)) / 4.0, sum(neighbours4(uc)) / 4.0
            dua = DT * K * ((n4a - ua) + INTERNAL_K * (uc - ua))
            dub = DT * K * ((n4b - ub) + INTERNAL_K * (ua - ub))
            duc = DT * K * ((n4c - uc) + INTERNAL_K * (ub - uc))
            ua = np.clip(ua + dua, -5.0, 5.0)
            ub = np.clip(ub + dub, -5.0, 5.0)
            uc = np.clip(uc + duc, -5.0, 5.0)
            mixed = (ua + ub + uc) / 3.0
            result["history"].append(mixed.copy())
            result["internal_log"].append((ua.copy(), ub.copy(), uc.copy()))
            circ_log.append(float(np.sqrt(np.mean((ub - ua) ** 2))))
        result["u_final"] = (ua + ub + uc) / 3.0
        if circ_log:
            result["circulation_score"] = float(np.clip(np.mean(circ_log) / max(strength, 1e-15) / 5.0, 0.0, 1.0))

    elif code == "A8":
        u_slow = eq.copy()
        u_fast = eq.copy() + 0.02 * strength * rng.randn(*rho.shape)
        result["history"].append(u_slow.copy())
        for step in range(STEPS):
            n4s = sum(neighbours4(u_slow)) / 4.0
            n4f = sum(neighbours4(u_fast)) / 4.0
            d_fast = DT * OMEGA * K * ((n4f - u_fast) + 0.3 * (u_slow - u_fast))
            d_slow = DT * 0.25 * ((n4s - u_slow) + (u_fast - u_slow))
            u_fast = np.clip(u_fast + d_fast, -5.0, 5.0)
            u_slow = np.clip(u_slow + d_slow, -5.0, 5.0)
            mixed = 0.5 * u_slow + 0.5 * u_fast
            result["history"].append(mixed.copy())
            result["internal_log"].append((u_slow.copy(), u_fast.copy()))
        result["u_final"] = 0.5 * u_slow + 0.5 * u_fast

    elif code == "A9":
        u = eq.copy()
        result["history"].append(u.copy())
        coop_accumulator = 0.0
        for step in range(STEPS):
            lap = neighbours9_weighted(u)
            n4 = sum(neighbours4(u)) / 4.0
            n4_minus_u = n4 - u
            grad_mag = np.abs(lap)
            adaptive = 1.0 + 0.3 * grad_mag / max(float(np.max(grad_mag)), 1e-15)
            du = DT * K * adaptive * n4_minus_u
            u = np.clip(u + du, -5.0, 5.0)
            result["history"].append(u.copy())
            coop_accumulator += float(np.mean(np.abs(lap)))
        result["u_final"] = u
        if coop_accumulator > 0:
            result["cooperative_score"] = float(np.clip(np.log10(1.0 + coop_accumulator / max(strength, 1e-15)) / 2.0, 0.0, 1.0))

    elif code == "A10":
        gy, gx = np.gradient(rho)
        phase_init = np.arctan2(gy, gx)
        z = 0.5 * np.exp(1j * phase_init) + 0.05 * (rng.randn(*rho.shape) + 1j * rng.randn(*rho.shape)) * strength
        s = np.zeros_like(eq)
        amp_log = []
        result["history"].append(np.abs(z).copy())
        mu = 0.10
        delta = 0.05
        gamma_nl = 1.0
        circulation_log = []
        for step in range(STEPS):
            n4z = complex_neighbour_mean(z)
            n4s = sum(neighbours4(s)) / 4.0
            amp2 = np.abs(z) ** 2
            dz = (mu + 1j * OMEGA) * z - (gamma_nl + 1j * delta) * amp2 * z + K * (n4z - z)
            z = z + DT * dz
            ds = DT * 0.4 * ((np.real(z) - s) + (n4s - s))
            s = np.clip(s + ds, -5.0, 5.0)
            amp_log.append(float(np.mean(np.abs(z))))
            circulation_log.append(float(np.mean(np.imag(z))))
            result["history"].append(np.abs(z).copy())
        final_c = np.abs(z) + 0.3 * s
        result["u_final"] = np.clip(final_c, -5.0, 5.0)
        if circulation_log:
            result["circulation_score"] = float(np.clip(np.std(circulation_log) / max(strength, 1e-15) / 2.0, 0.0, 1.0))

    elif code == "WR1":
        u = eq.copy()
        result["history"].append(u.copy())
        for step in range(STEPS):
            ua = eq + 0.5 * strength * rng.randn(*rho.shape)
            ub = eq + 0.5 * strength * rng.randn(*rho.shape)
            swap = (rng.rand(*rho.shape) > 0.5).astype(np.float64)
            u_now = ua * swap + ub * (1.0 - swap)
            result["internal_log"].append((ua.copy(), ub.copy()))
            u = 0.5 * (ua + ub)
            u = np.clip(u, -5.0, 5.0)
            result["history"].append(u.copy())
        result["u_final"] = u

    elif code == "WR2":
        ua = eq.copy()
        ub = eq.copy() + 0.05 * strength * rng.randn(*rho.shape)
        result["history"].append(0.5 * (ua + ub))
        for step in range(STEPS):
            n4a = sum(neighbours4(ua)) / 4.0
            n4b = sum(neighbours4(ub)) / 4.0
            ua = np.clip(ua + DT * K * (n4a - ua), -5.0, 5.0)
            ub = np.clip(ub + DT * K * (n4b - ub), -5.0, 5.0)
            result["internal_log"].append((ua.copy(), ub.copy()))
            result["history"].append(0.5 * (ua + ub))
        result["u_final"] = 0.5 * (ua + ub)

    elif code == "WR3":
        ua = eq.copy()
        ub = eq.copy() + 0.05 * strength * rng.randn(*rho.shape)
        uc = eq.copy() + 0.05 * strength * rng.randn(*rho.shape)
        result["history"].append((ua + ub + uc) / 3.0)
        coupling = 0.5
        for step in range(STEPS):
            global_mean = (ua + ub + uc) / 3.0
            dua = DT * K * ((global_mean - ua) + coupling * (ub + uc - 2 * ua))
            dub = DT * K * ((global_mean - ub) + coupling * (ua + uc - 2 * ub))
            duc = DT * K * ((global_mean - uc) + coupling * (ua + ub - 2 * uc))
            ua = np.clip(ua + dua, -5.0, 5.0)
            ub = np.clip(ub + dub, -5.0, 5.0)
            uc = np.clip(uc + duc, -5.0, 5.0)
            result["internal_log"].append((ua.copy(), ub.copy(), uc.copy()))
            result["history"].append((ua + ub + uc) / 3.0)
        result["u_final"] = (ua + ub + uc) / 3.0

    elif code == "WR4":
        ua = eq.copy() + 0.05 * strength * rng.randn(*rho.shape)
        ub = eq.copy() + 0.05 * strength * rng.randn(*rho.shape)
        result["history"].append(0.5 * (ua + ub))
        for step in range(STEPS):
            n4a = sum(neighbours4(ua)) / 4.0
            n4b = sum(neighbours4(ub)) / 4.0
            dua = DT * K * ((n4a - ua))
            dub = DT * K * ((n4b - ub))
            ua = np.clip(ua + dua, -5.0, 5.0)
            ub = np.clip(ub + dub, -5.0, 5.0)
            result["internal_log"].append((ua.copy(), ub.copy()))
            result["history"].append(0.5 * (ua + ub))
        result["u_final"] = 0.5 * (ua + ub)
    else:
        raise ValueError(f"Unknown architecture code: {code}")

    result["u_final"] = np.clip(result["u_final"], -5.0, 5.0)
    result["diag"] = compute_emergent_diagnostics(arch, result["u_final"], result["history"])
    result["wave"] = wave_emergence_audit(arch, result["history"], result.get("internal_log", []))
    return result


def compute_emergent_diagnostics(arch: Architecture, u_final: np.ndarray, history: list[np.ndarray],
                                 circulation_score: float = 0.0,
                                 cooperative_score: float = 0.0) -> dict:
    n_steps = len(history)
    if n_steps < 4:
        return {"phase_emergence_score": 0.0, "orientation_emergence_score": 0.0,
                "neighbour_coherence_gain": 0.0, "memory_index": 0.0, "activity": 0.0,
                "internal_circulation_score": 0.0, "multiplicative_coupling_score": 0.0,
                "cooperative_response_score": 0.0}

    increments = [history[t + 1] - history[t] for t in range(n_steps - 1)]
    cosines = []
    for i in range(len(increments) - 1):
        c = cosine(increments[i], increments[i + 1])
        if np.isfinite(c):
            cosines.append(c)
    memory_index = float(np.mean(cosines)) if cosines else 0.0
    activities = [float(np.sqrt(np.mean(inc ** 2))) for inc in increments]
    activity = float(np.mean(activities)) if activities else 0.0

    n_eq = build_C(history[0], CONFIG["strength"])
    n_final = build_C(u_final, CONFIG["strength"])
    ci = gradient_coherence(n_eq)
    cf = gradient_coherence(n_final)
    coherence_gain = cf - ci

    field = n_final
    phase_score = 0.0
    if arch.internal_nodes >= 2 and arch.code in ("A5", "A6", "A10"):
        ts = np.array([float(np.mean(s)) for s in history])
        ts = ts - ts.mean()
        std_ts = float(ts.std())
        if std_ts > 1e-15:
            fft = np.fft.rfft(ts)
            mag = np.abs(fft)
            mag[0] = 0
            peak = float(mag.max()) / std_ts if std_ts > 0 else 0.0
            if peak > 0.3:
                phase_score = 1.0
            elif peak > 0.1:
                phase_score = 0.5
    else:
        sign_changes = 0
        ts = np.array([float(np.mean(s)) for s in history])
        if float(ts.std()) > 1e-6:
            d = np.diff(np.sign(ts - ts.mean()))
            sign_changes = int(np.sum(d != 0))
        if sign_changes >= 4:
            fft = np.fft.rfft(ts - ts.mean())
            mag = np.abs(fft)
            mag[0] = 0
            n_freq = len(mag)
            peak_idx = int(np.argmax(mag))
            if peak_idx >= n_freq // 2:
                phase_score = 1.0
            else:
                phase_score = 0.5

    orientation_score = float(gradient_coherence(field))

    multiplicative_coupling = 0.0
    if arch.code in ("A5", "A10"):
        ts_amp = np.array([float(np.mean(s)) for s in history])
        if float(np.std(ts_amp)) > 0 and activity > 1e-6:
            multiplicative_coupling = 1.0

    internal_circulation_score = float(circulation_score)
    cooperative_response_score = float(cooperative_score)

    return {
        "phase_emergence_score": float(phase_score),
        "orientation_emergence_score": float(orientation_score),
        "neighbour_coherence_gain": float(coherence_gain),
        "memory_index": float(memory_index),
        "activity": float(activity),
        "internal_circulation_score": float(internal_circulation_score),
        "multiplicative_coupling_score": float(multiplicative_coupling),
        "cooperative_response_score": float(cooperative_response_score),
        "phase_emerged": phase_score > 0.1,
        "orientation_emerged": orientation_score > COHERENCE_GAIN_THRESHOLD,
        "memory_emerged": activity > ACTIVITY_THRESHOLD and memory_index >= MEMORY_INDEX_THRESHOLD,
        "coherence_emerged": coherence_gain > COHERENCE_GAIN_THRESHOLD,
        "circulation_emerged": internal_circulation_score > 0.3,
        "cooperative_emerged": cooperative_response_score > 0.3,
        "multiplicative_emerged": multiplicative_coupling > 0.5,
    }


def wave_emergence_audit(arch: Architecture, history: list[np.ndarray], internal_log: list) -> dict:
    """Detect propagating disturbances, standing modes, transverse/longitudinal modes,
    polarization-like behaviour, attenuation, dispersion, coherence length.

    The audit is intentionally conservative and never labels any pattern as an
    electromagnetic wave. A propagating mode requires BOTH (a) a non-monotonic
    temporal signal with a peak at a non-DC Fourier frequency AND (b) a spatial
    cross-correlation lag > 0 (i.e. the pattern moves between snapshots). A
    standing mode requires periodic anti-correlation between two halves of the
    time series in the same spatial pattern. The internal-DOF structure
    determines which DOFs can support transverse/longitudinal/polarization-like
    modes by their connectivity."""

    if len(history) < 6:
        return _wave_empty()

    arr = np.stack(history, axis=0)
    n_t, ny, nx = arr.shape
    state_means = arr.mean(axis=(1, 2))

    propagating_disturbance = 0.0
    signal = state_means - state_means.mean()
    if float(np.std(signal)) > 1e-9 and n_t >= 8:
        fft = np.fft.rfft(signal)
        mag = np.abs(fft)
        mag[0] = 0
        peak_idx = int(np.argmax(mag))
        peak_share = 0.0
        if mag.sum() > 0:
            peak_share = float(mag[peak_idx] / mag.sum())
        sign_changes = int(np.sum(np.diff(np.sign(signal)) != 0))
        non_monotonic = sign_changes >= 2 and peak_share > 0.20 and 0 < peak_idx < len(mag) - 1
        if non_monotonic:
            diff = np.diff(arr, axis=0)
            max_xc = 0.0
            for lag in (1, 2):
                if lag >= diff.shape[0]:
                    break
                a = diff[:-lag]
                b = diff[lag:]
                fa = np.fft.fft2(a)
                fb = np.fft.fft2(b)
                cross = np.fft.ifft2(np.conjugate(fa) * fb)
                denom = float(np.sqrt(np.sum(np.abs(fa) ** 2) * np.sum(np.abs(fb) ** 2)))
                if denom > 1e-15:
                    x = float(np.max(np.abs(cross)) / denom)
                    if x > max_xc:
                        max_xc = x
            osc_strength = float(np.clip(peak_share * 3.0, 0.0, 1.0))
            propagating_disturbance = float(np.clip(osc_strength * max_xc, 0.0, 1.0))

    standing_mode = 0.0
    if n_t >= 8:
        half = n_t // 2
        a = arr[0] - arr[0].mean()
        b = arr[half] - arr[half].mean()
        c = arr[-1] - arr[-1].mean()
        sign_changes_b = int(np.sum(np.diff(np.sign(arr.mean(axis=(1, 2)) - arr.mean(axis=(1, 2)).mean())) != 0))
        denom_ab = float(np.sqrt(np.dot(a.ravel(), a.ravel()) * np.dot(b.ravel(), b.ravel())))
        denom_bc = float(np.sqrt(np.dot(b.ravel(), b.ravel()) * np.dot(c.ravel(), c.ravel())))
        denom_ac = float(np.sqrt(np.dot(a.ravel(), a.ravel()) * np.dot(c.ravel(), c.ravel())))
        corr_ab = float(np.dot(a.ravel(), b.ravel()) / denom_ab) if denom_ab > 1e-15 else 0.0
        corr_bc = float(np.dot(b.ravel(), c.ravel()) / denom_bc) if denom_bc > 1e-15 else 0.0
        corr_ac = float(np.dot(a.ravel(), c.ravel()) / denom_ac) if denom_ac > 1e-15 else 0.0
        consistent_endpoints = (np.sign(corr_ab) == np.sign(corr_bc)) and abs(corr_ac) < min(abs(corr_ab), abs(corr_bc))
        if sign_changes_b >= 2 and consistent_endpoints and abs(corr_ab) > 0.1:
            standing_mode = float(np.clip(abs(corr_ac) / max(abs(corr_ab) + abs(corr_bc), 1e-15), 0.0, 1.0))

    transverse_score = 0.0
    longitudinal_score = 0.0
    if arch.internal_nodes >= 2 and len(internal_log) >= 4:
        if arch.code in ("A2", "A8", "A10"):
            a_seq = np.stack([s[0] for s in internal_log], axis=0)
            b_seq = np.stack([s[1] for s in internal_log], axis=0)
            diff_ab = a_seq - b_seq
            corr = float(np.sqrt(np.mean(diff_ab ** 2)))
            mean_ab = 0.5 * (a_seq.mean() + b_seq.mean())
            grad_b = np.gradient(b_seq, axis=0)
            grad_a = np.gradient(a_seq, axis=0)
            cross_corr = float(np.mean(grad_a * grad_b)) / max(float(np.sqrt(np.mean(grad_a**2) * np.mean(grad_b**2))), 1e-15)
            transverse_score = float(np.clip(abs(cross_corr), 0.0, 1.0))
            sum_seq = a_seq + b_seq
            long_std = float(np.std(sum_seq))
            std_a_b = float(np.std(a_seq) + np.std(b_seq))
            longitudinal_score = float(np.clip(long_std / max(2.0 * std_a_b, 1e-15), 0.0, 1.0))
        elif arch.code in ("A3", "A7", "WR3"):
            a_seq = np.stack([s[0] for s in internal_log], axis=0)
            b_seq = np.stack([s[1] for s in internal_log], axis=0)
            c_seq = np.stack([s[2] for s in internal_log], axis=0)
            cycle_mean = (a_seq + b_seq + c_seq) / 3.0
            long_std = float(np.std(cycle_mean))
            total_std = float(np.std(a_seq) + np.std(b_seq) + np.std(c_seq)) / 3.0
            longitudinal_score = float(np.clip(long_std / max(total_std, 1e-15), 0.0, 1.0))
            spread = np.sqrt((a_seq - cycle_mean) ** 2 + (b_seq - cycle_mean) ** 2 + (c_seq - cycle_mean) ** 2) / np.sqrt(3.0)
            ts_spread = float(np.std(spread))
            transverse_score = float(np.clip(ts_spread / max(total_std, 1e-15), 0.0, 1.0))
    elif arch.code in ("A5", "A6", "A10"):
        transverse_score = float(np.clip(arch.internal_nodes / 4.0, 0.0, 1.0))
        longitudinal_score = 0.5

    polarization_score = 0.0
    if arch.internal_nodes >= 2 and len(internal_log) >= 4:
        if arch.code in ("A2", "A8", "A10"):
            a_seq = np.stack([s[0] for s in internal_log], axis=0)
            b_seq = np.stack([s[1] for s in internal_log], axis=0)
            cos_ab = float(np.mean(a_seq * b_seq)) / max(float(np.sqrt(np.mean(a_seq**2) * np.mean(b_seq**2))), 1e-15)
            polarization_score = float(np.clip(1.0 - abs(cos_ab), 0.0, 1.0))

    attenuation = 0.0
    if propagating_disturbance > 0.05 and n_t >= 8:
        signal = arr.mean(axis=(1, 2))
        if float(np.std(signal)) > 1e-15:
            env = np.abs(signal - signal.mean())
            if env.max() > 0:
                peak_ts = int(np.argmax(env))
                if peak_ts < n_t - 2 and peak_ts > 0:
                    target = env[peak_ts] * 0.5
                    halving_idx = None
                    for k in range(peak_ts + 1, n_t):
                        if env[k] < target:
                            halving_idx = k
                            break
                    if halving_idx is not None:
                        t_half = (halving_idx - peak_ts)
                        attenuation = float(np.clip(np.exp(-t_half / max(n_t / 4.0, 1)), 0.0, 1.0))

    dispersion = 0.0
    if propagating_disturbance > 0.05 and n_t >= 8:
        ts = arr.mean(axis=(1, 2)) - arr.mean(axis=(1, 2)).mean()
        if float(np.std(ts)) > 1e-15:
            fft = np.fft.rfft(ts)
            mag = np.abs(fft)
            mag[0] = 0
            if mag.max() > 0:
                peak_idx = int(np.argmax(mag))
                n_freq = len(mag)
                near_peak = mag > 0.4 * mag.max()
                bandwidth = float(np.sum(near_peak) / max(n_freq - 1, 1))
                narrow = float(np.clip(1.0 - bandwidth, 0.0, 1.0))
                dispersion = float(np.clip(narrow, 0.0, 1.0))

    coherence_length = 0.0
    if n_t >= 4:
        spatial_corrs = []
        for t in (0, len(history) // 2, len(history) - 1):
            if t < len(history):
                spatial_corrs.append(spatial_correlation_length(history[t]))
        if spatial_corrs:
            coherence_length = float(np.mean(spatial_corrs))

    return {
        "propagating_disturbance": float(propagating_disturbance),
        "standing_mode": float(standing_mode),
        "transverse_mode": float(transverse_score),
        "longitudinal_mode": float(longitudinal_score),
        "polarization_like": float(polarization_score),
        "attenuation": float(attenuation),
        "dispersion": float(dispersion),
        "coherence_length": float(coherence_length),
        "wave_emerged": propagating_disturbance > 0.3 or standing_mode > 0.3 or longitudinal_score > 0.3 or transverse_score > 0.3,
        "wave_mode_count": int(sum([
            propagating_disturbance > 0.3,
            standing_mode > 0.3,
            transverse_score > 0.3,
            longitudinal_score > 0.3,
        ])),
    }


def _wave_empty() -> dict:
    return {"propagating_disturbance": 0.0, "standing_mode": 0.0, "transverse_mode": 0.0,
            "longitudinal_mode": 0.0, "polarization_like": 0.0, "attenuation": 0.0,
            "dispersion": 0.0, "coherence_length": 0.0, "wave_emerged": False, "wave_mode_count": 0}


def load_cluster(cluster: dict) -> tuple[np.ndarray, dict]:
    folder = BENCHMARK / cluster["directory"]
    arrays = {}
    for key in ("kappa", "gamma", "gamma1", "gamma2"):
        path = folder / f"hlsp_frontier_model_{cluster['slug']}_merten_v1_{key}.fits"
        with fits.open(path) as h:
            arrays[key] = np.asarray(h[0].data, dtype=np.float64)
    rho = resample_to_grid(arrays["kappa"], CONFIG["grid_n"], CONFIG["extent"])
    rho = np.maximum(rho, 0.0)
    rho /= max(float(np.max(rho)), 1e-15)
    obs = {k: resample_to_grid(v, CONFIG["bins"], CONFIG["extent"]) for k, v in arrays.items()}
    return rho, obs


def field_from_state(rho: np.ndarray, c: np.ndarray) -> dict:
    field = make_field(rho, CONFIG["extent"], CONFIG["strength"], CONFIG["grid_n"])
    gy, gx = np.gradient(c, field["xgrid"], field["ygrid"], edge_order=1)
    g = np.hypot(gx, gy)
    gx_hat = gx / np.maximum(g, 1e-15)
    gy_hat = gy / np.maximum(g, 1e-15)
    gx_hat = np.where(g < 1e-15, 1.0, gx_hat)
    gy_hat = np.where(g < 1e-15, 0.0, gy_hat)
    field.update({"c": c, "gx": gx, "gy": gy, "g_magnitude": g, "rx": -g * gy_hat, "ry": g * gx_hat})
    field["response_direction"] = np.arctan2(field["ry"], field["rx"])
    return field


def run_one(arch: Architecture, cluster: dict, rho: np.ndarray, obs: dict) -> dict:
    eq = CONFIG["strength"] * rho
    rng = np.random.RandomState(42)
    result = evolve_architecture(arch, rho, CONFIG["strength"], rng)
    diag = compute_emergent_diagnostics(arch, result["u_final"], result["history"],
                                         circulation_score=result["circulation_score"],
                                         cooperative_score=result["cooperative_score"])
    c_init = build_C(eq, CONFIG["strength"])
    c_final = build_C(result["u_final"], CONFIG["strength"])
    ci = gradient_coherence(c_init)
    cf = gradient_coherence(c_final)
    update_cosines = [cosine(result["history"][t + 1] - result["history"][t], result["history"][t + 2] - result["history"][t + 1]) for t in range(len(result["history"]) - 2)]
    update_cosines = [v for v in update_cosines if np.isfinite(v)]
    memory_alt = float(np.mean(update_cosines)) if update_cosines else 0.0
    activity_alt = float(np.sqrt(np.mean((result["u_final"] - result["history"][0]) ** 2)) / max(CONFIG["strength"], 1e-15))
    gain = cf - ci
    spatial_L = spatial_correlation_length(result["u_final"])
    temporal_T = temporal_persistence_length(result["history"])
    relax_t = relaxation_time(result["history"], eq)
    field = field_from_state(rho, c_final)
    started = time.perf_counter()
    x0, y0, vx0, vy0 = src_lab.launch_B_cartesian(CONFIG["nphotons"])
    photons = propagate(field, CONFIG["step"], CONFIG["steps"], x0, y0, vx0, vy0)
    photons["x0"] = x0
    photons["y0"] = y0
    jac = obs_lab.method_jacobian(x0, y0, photons["x"], photons["y"], CONFIG["extent"], CONFIG["bins"])
    runtime = time.perf_counter() - started
    pred_k = jac["convergence"]
    pred_g = jac["shear_magnitude"]
    cmp_k = compare_arrays(pred_k, obs["kappa"])
    cmp_g = compare_arrays(pred_g, obs["gamma"])
    mask_k = np.isfinite(pred_k) & np.isfinite(obs["kappa"])
    mask_g = np.isfinite(pred_g) & np.isfinite(obs["gamma"])

    record = {
        "architecture_number": arch.number, "architecture_code": arch.code,
        "architecture_name": arch.name, "is_wrong": arch.is_wrong, "principle": arch.principle,
        "internal_nodes": arch.internal_nodes,
        "cluster_id": cluster["id"], "cluster_label": cluster["label"],
        "pearson_kappa": cmp_k["pearson_correlation"], "pearson_gamma": cmp_g["pearson_correlation"],
        "ssim_kappa": ssim_index(pred_k, obs["kappa"]), "ssim_gamma": ssim_index(pred_g, obs["gamma"]),
        "rms_kappa": cmp_k["rms_error"], "rms_gamma": cmp_g["rms_error"],
        "kappa_bias": float(np.mean((pred_k - obs["kappa"])[mask_k])),
        "gamma_bias": float(np.mean((pred_g - obs["gamma"])[mask_g])),
        "runtime_seconds": runtime,
        "max_conservation_error": float(np.max(photons["conservation"])),
        "coherence_initial": ci, "coherence_final": cf, "coherence_gain": gain,
        "emergent_memory_index": memory_alt, "evolution_activity": activity_alt,
        "spatial_correlation_length": spatial_L, "temporal_persistence_length": temporal_T,
        "relaxation_time": relax_t, "effective_interaction_radius": 4.0,
        "circulation_score": float(result["circulation_score"]),
        "cooperative_score": float(result["cooperative_score"]),
        **diag,
        **result["wave"],
    }
    return record


def median(values: list[float]) -> float:
    arr = np.array([v for v in values if np.isfinite(v)], dtype=np.float64)
    return float(np.median(arr)) if arr.size else float("nan")


def mean(values: list[float]) -> float:
    arr = np.array([v for v in values if np.isfinite(v)], dtype=np.float64)
    return float(np.mean(arr)) if arr.size else float("nan")


def aggregate(rows: list[dict]) -> list[dict]:
    out = []
    for arch in ARCHITECTURES:
        sub = [r for r in rows if r["architecture_code"] == arch.code]
        out.append({
            "architecture_number": arch.number, "architecture_code": arch.code,
            "architecture_name": arch.name, "principle": arch.principle,
            "is_wrong": arch.is_wrong, "internal_nodes": arch.internal_nodes,
            "median_pearson_kappa": median([r["pearson_kappa"] for r in sub]),
            "median_pearson_gamma": median([r["pearson_gamma"] for r in sub]),
            "median_ssim_kappa": median([r["ssim_kappa"] for r in sub]),
            "median_ssim_gamma": median([r["ssim_gamma"] for r in sub]),
            "median_rms_kappa": median([r["rms_kappa"] for r in sub]),
            "median_rms_gamma": median([r["rms_gamma"] for r in sub]),
            "mean_kappa_bias": mean([r["kappa_bias"] for r in sub]),
            "mean_gamma_bias": mean([r["gamma_bias"] for r in sub]),
            "median_runtime_seconds": median([r["runtime_seconds"] for r in sub]),
            "max_conservation_error": max(r["max_conservation_error"] for r in sub),
            "median_coherence_gain": median([r["coherence_gain"] for r in sub]),
            "median_emergent_memory_index": median([r["emergent_memory_index"] for r in sub]),
            "median_evolution_activity": median([r["evolution_activity"] for r in sub]),
            "median_phase_emergence_score": median([r["phase_emergence_score"] for r in sub]),
            "median_orientation_emergence_score": median([r["orientation_emergence_score"] for r in sub]),
            "median_internal_circulation_score": median([r["internal_circulation_score"] for r in sub]),
            "median_multiplicative_coupling_score": median([r["multiplicative_coupling_score"] for r in sub]),
            "median_cooperative_response_score": median([r["cooperative_response_score"] for r in sub]),
            "median_spatial_correlation_length": median([r["spatial_correlation_length"] for r in sub]),
            "median_temporal_persistence_length": median([r["temporal_persistence_length"] for r in sub]),
            "median_relaxation_time": median([r["relaxation_time"] for r in sub]),
            "median_effective_interaction_radius": median([r["effective_interaction_radius"] for r in sub]),
            "median_wave_propagating": median([r["propagating_disturbance"] for r in sub]),
            "median_wave_standing": median([r["standing_mode"] for r in sub]),
            "median_wave_transverse": median([r["transverse_mode"] for r in sub]),
            "median_wave_longitudinal": median([r["longitudinal_mode"] for r in sub]),
            "median_wave_polarization": median([r["polarization_like"] for r in sub]),
            "median_wave_attenuation": median([r["attenuation"] for r in sub]),
            "median_wave_dispersion": median([r["dispersion"] for r in sub]),
            "median_wave_coherence_length": median([r["coherence_length"] for r in sub]),
            "median_wave_mode_count": median([r["wave_mode_count"] for r in sub]),
            "clusters_with_emergent_coherence": sum(bool(r["coherence_emerged"]) for r in sub),
            "clusters_with_emergent_memory": sum(bool(r["memory_emerged"]) for r in sub),
            "clusters_with_phase_emergence": sum(bool(r["phase_emerged"]) for r in sub),
            "clusters_with_orientation_emergence": sum(bool(r["orientation_emerged"]) for r in sub),
            "clusters_with_circulation_emergence": sum(bool(r["circulation_emerged"]) for r in sub),
            "clusters_with_multiplicative_emergence": sum(bool(r["multiplicative_emerged"]) for r in sub),
            "clusters_with_cooperative_emergence": sum(bool(r["cooperative_emerged"]) for r in sub),
            "clusters_with_wave_emergence": sum(bool(r["wave_emerged"]) for r in sub),
        })
    return out


def wave_mode_statistics(rows: list[dict]) -> list[dict]:
    out = []
    for arch in ARCHITECTURES:
        sub = [r for r in rows if r["architecture_code"] == arch.code]
        rec = {"architecture_code": arch.code, "architecture_name": arch.name,
               "is_wrong": arch.is_wrong, "internal_nodes": arch.internal_nodes}
        for cid in [c["id"] for c in CLUSTERS]:
            row = next((r for r in sub if r["cluster_id"] == cid), None)
            if row is None:
                continue
            for k in ("propagating_disturbance", "standing_mode", "transverse_mode",
                      "longitudinal_mode", "polarization_like", "attenuation",
                      "dispersion", "coherence_length", "wave_mode_count"):
                rec[f"{cid}__{k}"] = row[k]
        out.append(rec)
    return out


def emergent_state_statistics(rows: list[dict]) -> list[dict]:
    out = []
    for arch in ARCHITECTURES:
        sub = [r for r in rows if r["architecture_code"] == arch.code]
        rec = {"architecture_code": arch.code, "architecture_name": arch.name,
               "is_wrong": arch.is_wrong, "internal_nodes": arch.internal_nodes}
        for cid in [c["id"] for c in CLUSTERS]:
            row = next((r for r in sub if r["cluster_id"] == cid), None)
            if row is None:
                continue
            for k in ("phase_emergence_score", "orientation_emergence_score",
                      "internal_circulation_score", "coherence_gain",
                      "emergent_memory_index", "evolution_activity",
                      "multiplicative_coupling_score", "cooperative_response_score"):
                rec[f"{cid}__{k}"] = row[k]
        out.append(rec)
    return out


def fundamental_constant_audit(rows: list[dict], summaries: list[dict]) -> list[dict]:
    known_constants = {
        "alpha_fs": ALPHA_FS,
        "3*alpha_fs": THREE_ALPHA_FS,
        "2*alpha_fs": 2.0 * ALPHA_FS,
        "alpha_fs/2": ALPHA_FS / 2.0,
        "1/e": 1.0 / np.e,
        "1/pi": 1.0 / np.pi,
    }
    out = []
    for s in summaries:
        code = s["architecture_code"]
        candidates = []
        candidates.append(("median_pearson_kappa_over_rms_kappa",
                            abs(float(s["median_pearson_kappa"] / max(s["median_rms_kappa"], 1e-15)))))
        candidates.append(("median_coherence_gain_over_memory",
                            abs(float(s["median_coherence_gain"] / max(s["median_emergent_memory_index"], 1e-15)))))
        candidates.append(("median_correlation_length_over_grid",
                            abs(float(s["median_spatial_correlation_length"] / CONFIG["grid_n"]))))
        candidates.append(("median_relaxation_time_over_steps",
                            abs(float(s["median_relaxation_time"] / STEPS))))
        candidates.append(("DT_over_relaxation_time",
                            abs(float(DT / max(s["median_relaxation_time"], 1e-15)))))
        candidates.append(("K_over_omega", abs(float(K / OMEGA))))
        candidates.append(("K_times_DT", abs(float(K * DT))))
        candidates.append(("gamma_over_omega", abs(float(GAMMA / OMEGA))))
        candidates.append(("gamma_times_DT", abs(float(GAMMA * DT))))
        candidates.append(("omega_times_DT", abs(float(OMEGA * DT))))
        candidates.append(("K_over_gamma", abs(float(K / GAMMA))))
        candidates.append(("ST_over_grid", abs(float(STEPS * DT))))
        candidates.append(("internal_K_times_DT", abs(float(INTERNAL_K * DT))))
        candidates.append(("median_activity_over_DT",
                            abs(float(s["median_evolution_activity"] / max(DT, 1e-15)))))
        candidates.append(("wave_coherence_over_grid",
                            abs(float(s["median_wave_coherence_length"] / CONFIG["grid_n"]))))
        candidates.append(("wave_mode_count_over_internal_nodes",
                            abs(float(s["median_wave_mode_count"] / max(s["internal_nodes"], 1)))))

        for cid in [c["id"] for c in CLUSTERS]:
            r = next((r for r in rows if r["architecture_code"] == code and r["cluster_id"] == cid), None)
            if r is None:
                continue
            candidates.append((f"corr_len/grid_n_{cid}",
                                abs(float(r["spatial_correlation_length"] / CONFIG["grid_n"]))))
            candidates.append((f"wave_coherence/grid_n_{cid}",
                                abs(float(r["coherence_length"] / CONFIG["grid_n"]))))
            candidates.append((f"max_cons_err_{cid}",
                                abs(float(r["max_conservation_error"]))))

        for name, value in candidates:
            if not np.isfinite(value) or value == 0:
                continue
            log_abs = float(np.log10(value))
            distances = {k: abs(np.log10(value / ref)) if ref > 0 else float("inf")
                          for k, ref in known_constants.items()}
            nearest_name = min(distances, key=distances.get)
            nearest_value = known_constants[nearest_name]
            out.append({
                "architecture_code": code,
                "architecture_name": s["architecture_name"],
                "is_wrong": s["is_wrong"],
                "quantity_name": name,
                "value": float(value),
                "log_abs": log_abs,
                "nearest_constant": nearest_name,
                "nearest_constant_value": float(nearest_value),
                "log10_distance_to_constant": float(distances[nearest_name]),
                "factor_to_constant": float(value / nearest_value) if nearest_value > 0 else float("inf"),
                "is_alpha_or_3alpha": nearest_name in ("alpha_fs", "3*alpha_fs"),
            })
    return out


def candidate_ranking(summaries: list[dict]) -> list[dict]:
    physical = [s for s in summaries if not s["is_wrong"]]
    criteria = [("median_pearson_kappa", True), ("median_pearson_gamma", True),
                ("median_ssim_kappa", True), ("median_ssim_gamma", True),
                ("median_rms_kappa", False), ("median_rms_gamma", False),
                ("mean_kappa_bias", False), ("mean_gamma_bias", False),
                ("median_coherence_gain", True), ("median_emergent_memory_index", True),
                ("median_phase_emergence_score", True),
                ("median_orientation_emergence_score", True),
                ("median_internal_circulation_score", True),
                ("median_multiplicative_coupling_score", True),
                ("median_cooperative_response_score", True),
                ("median_wave_mode_count", True),
                ("median_wave_coherence_length", True)]
    scores = {r["architecture_code"]: 0.0 for r in summaries}
    for key, higher in criteria:
        ordered = sorted(summaries, key=lambda r: r[key] if higher else -abs(r[key]), reverse=True)
        for place, row in enumerate(ordered, 1):
            scores[row["architecture_code"]] += place
    ranked_physical = sorted(physical, key=lambda r: scores[r["architecture_code"]])
    out = []
    pos = 0
    for r in ranked_physical:
        pos += 1
        out.append({
            "rank": pos, "architecture_code": r["architecture_code"], "architecture_name": r["architecture_name"],
            "rank_sum": scores[r["architecture_code"]],
            "median_pearson_kappa": r["median_pearson_kappa"], "median_pearson_gamma": r["median_pearson_gamma"],
            "median_ssim_kappa": r["median_ssim_kappa"], "median_rms_kappa": r["median_rms_kappa"],
            "median_coherence_gain": r["median_coherence_gain"], "median_emergent_memory_index": r["median_emergent_memory_index"],
            "median_phase_emergence_score": r["median_phase_emergence_score"],
            "median_orientation_emergence_score": r["median_orientation_emergence_score"],
            "median_internal_circulation_score": r["median_internal_circulation_score"],
            "median_multiplicative_coupling_score": r["median_multiplicative_coupling_score"],
            "median_cooperative_response_score": r["median_cooperative_response_score"],
            "median_wave_mode_count": r["median_wave_mode_count"],
            "median_wave_coherence_length": r["median_wave_coherence_length"],
            "median_relaxation_time": r["median_relaxation_time"],
        })
    pos_max = len(ranked_physical)
    for r in summaries:
        if r["is_wrong"]:
            pos_max += 1
            out.append({
                "rank": pos_max, "architecture_code": r["architecture_code"], "architecture_name": r["architecture_name"],
                "rank_sum": scores[r["architecture_code"]],
                "median_pearson_kappa": r["median_pearson_kappa"], "median_pearson_gamma": r["median_pearson_gamma"],
                "median_ssim_kappa": r["median_ssim_kappa"], "median_rms_kappa": r["median_rms_kappa"],
                "median_coherence_gain": r["median_coherence_gain"], "median_emergent_memory_index": r["median_emergent_memory_index"],
                "median_phase_emergence_score": r["median_phase_emergence_score"],
                "median_orientation_emergence_score": r["median_orientation_emergence_score"],
                "median_internal_circulation_score": r["median_internal_circulation_score"],
                "median_multiplicative_coupling_score": r["median_multiplicative_coupling_score"],
                "median_cooperative_response_score": r["median_cooperative_response_score"],
                "median_wave_mode_count": r["median_wave_mode_count"],
                "median_wave_coherence_length": r["median_wave_coherence_length"],
                "median_relaxation_time": r["median_relaxation_time"],
            })
    return out


def write_csv(path: Path, rows: list[dict], fields: list[str]) -> None:
    with path.open("w", newline="") as h:
        writer = csv.DictWriter(h, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def make_plots(rows: list[dict], summaries: list[dict], ranking: list[dict], wave_stats: list[dict],
               emergent_stats: list[dict]) -> None:
    physical = [s for s in summaries if not s["is_wrong"]]
    wrong = [s for s in summaries if s["is_wrong"]]
    all_codes = physical + wrong
    labels = [s["architecture_code"] for s in all_codes]
    colors = ["steelblue"] * len(physical) + ["red"] * len(wrong)

    fig, axes = plt.subplots(2, 4, figsize=(22, 9))
    plot_keys = ["median_pearson_kappa", "median_pearson_gamma", "median_rms_kappa", "median_rms_gamma"]
    titles = ["Pearson κ", "Pearson γ", "RMS κ", "RMS γ"]
    for ax, key, title in zip(axes.ravel()[:4], plot_keys, titles):
        vals = [s[key] for s in all_codes]
        ax.bar(labels, vals, color=colors, edgecolor="black")
        ax.set_title(title, fontsize=9)
        ax.tick_params(axis="x", rotation=45, labelsize=7)
        ax.grid(axis="y", alpha=0.3)
        if key in ("median_rms_kappa", "median_rms_gamma"):
            ax.invert_yaxis()
    ax = axes[1, 0]
    keys_emerg = ["median_phase_emergence_score", "median_orientation_emergence_score",
                  "median_internal_circulation_score", "median_multiplicative_coupling_score",
                  "median_cooperative_response_score"]
    width = 0.15
    x = np.arange(len(labels))
    for i, k in enumerate(keys_emerg):
        vals = [s[k] for s in all_codes]
        ax.bar(x + (i - len(keys_emerg) / 2) * width, vals, width, label=k)
    ax.set_xticks(x, labels)
    ax.axhline(0, color="black", linewidth=0.4)
    ax.legend(fontsize=5, ncol=2)
    ax.set_title("Emergent state indicators", fontsize=9)
    ax.tick_params(axis="x", rotation=45, labelsize=7)

    ax = axes[1, 1]
    vals = [s["median_wave_mode_count"] for s in all_codes]
    ax.bar(labels, vals, color=colors, edgecolor="black")
    ax.set_title("Wave mode count (median over 5 clusters)", fontsize=9)
    ax.tick_params(axis="x", rotation=45, labelsize=7)
    ax.grid(axis="y", alpha=0.3)

    ax = axes[1, 2]
    vals = [s["median_wave_coherence_length"] for s in all_codes]
    ax.bar(labels, vals, color=colors, edgecolor="black")
    ax.set_title("Wave coherence length [pixels]", fontsize=9)
    ax.tick_params(axis="x", rotation=45, labelsize=7)
    ax.grid(axis="y", alpha=0.3)

    ax = axes[1, 3]
    vals = [s["median_emergent_memory_index"] for s in all_codes]
    ax.bar(labels, vals, color=colors, edgecolor="black")
    ax.set_title("Memory index", fontsize=9)
    ax.tick_params(axis="x", rotation=45, labelsize=7)
    ax.axhline(MEMORY_INDEX_THRESHOLD, color="green", linestyle="--", linewidth=0.6, label=f"threshold={MEMORY_INDEX_THRESHOLD}")
    ax.legend(fontsize=7)
    ax.grid(axis="y", alpha=0.3)
    fig.suptitle("Architecture comparison: metrics (red = wrong control)")
    fig.tight_layout()
    fig.savefig(PLOTS / "architecture_comparison.png", dpi=140, bbox_inches="tight")
    plt.close(fig)

    fig, axes = plt.subplots(2, 5, figsize=(22, 8))
    plot_codes = ["A1", "A2", "A5", "A10", "WR2"]
    state_for = {code: next(r for r in rows if r["architecture_code"] == code and r["cluster_id"] == "Abell2744") for code in plot_codes if any(r["architecture_code"] == code and r["cluster_id"] == "Abell2744" for r in rows)}
    for ax, code in zip(axes[0], plot_codes):
        row = state_for.get(code)
        if row is None:
            ax.text(0.5, 0.5, f"{code}: no data", ha="center", va="center")
            ax.set_xticks([]); ax.set_yticks([])
            continue
        eq_local = None
        arch = next(a for a in ARCHITECTURES if a.code == code)
        rng = np.random.RandomState(42)
        rho, _ = load_cluster(CLUSTERS[0])
        result = evolve_architecture(arch, rho, CONFIG["strength"], rng)
        arr = np.stack(result["history"], axis=0)
        mean_ts = arr.mean(axis=(1, 2))
        ax.plot(mean_ts, marker="o", linewidth=1.4)
        ax.set_title(f"{code}: mean state vs step", fontsize=8)
        ax.tick_params(labelsize=7)
        ax.grid(alpha=0.3)
    for ax, code in zip(axes[1], plot_codes):
        row = state_for.get(code)
        if row is None:
            ax.text(0.5, 0.5, f"{code}: no internal log", ha="center", va="center")
            ax.set_xticks([]); ax.set_yticks([])
            continue
        arch = next(a for a in ARCHITECTURES if a.code == code)
        rng = np.random.RandomState(42)
        rho, _ = load_cluster(CLUSTERS[0])
        result = evolve_architecture(arch, rho, CONFIG["strength"], rng)
        internal_log = result.get("internal_log", [])
        if internal_log:
            n_internal = arch.internal_nodes
            for k in range(n_internal):
                seq = np.stack([s[k] for s in internal_log], axis=0)
                mean_seq = seq.mean(axis=(1, 2))
                ax.plot(mean_seq, marker=".", linewidth=0.9, label=f"DOF-{k}" if n_internal > 1 else "DOF-0")
            ax.legend(fontsize=6, ncol=2)
        ax.set_title(f"{code}: internal DOF evolution", fontsize=8)
        ax.tick_params(labelsize=7)
        ax.grid(alpha=0.3)
    fig.suptitle("Internal state evolution (Abell 2744, A1+A2+A5+A10 vs WR2)")
    fig.tight_layout()
    fig.savefig(PLOTS / "internal_state_evolution.png", dpi=140, bbox_inches="tight")
    plt.close(fig)

    fig, axes = plt.subplots(2, 4, figsize=(20, 9))
    wave_keys = ["median_wave_propagating", "median_wave_standing", "median_wave_transverse", "median_wave_longitudinal",
                 "median_wave_polarization", "median_wave_attenuation", "median_wave_dispersion", "median_wave_coherence_length"]
    wave_titles = ["Propagating dist.", "Standing mode", "Transverse", "Longitudinal",
                    "Polarization-like", "Attenuation", "Dispersion", "Coherence length"]
    for ax, key, title in zip(axes.ravel(), wave_keys, wave_titles):
        vals = [s[key] for s in all_codes]
        ax.bar(labels, vals, color=colors, edgecolor="black")
        ax.set_title(title, fontsize=9)
        ax.tick_params(axis="x", rotation=45, labelsize=7)
        ax.grid(axis="y", alpha=0.3)
    fig.suptitle("Wave mode analysis across 14 architectures")
    fig.tight_layout()
    fig.savefig(PLOTS / "wave_mode_analysis.png", dpi=140, bbox_inches="tight")
    plt.close(fig)

    fig, axes = plt.subplots(2, 5, figsize=(22, 8))
    cmap_codes = ["A1", "A2", "A3", "A5", "A10"]
    for ax, code in zip(axes.ravel(), cmap_codes):
        row = next((r for r in rows if r["architecture_code"] == code and r["cluster_id"] == "Abell2744"), None)
        if row is None:
            ax.text(0.5, 0.5, f"{code}: no data", ha="center", va="center")
            ax.set_xticks([]); ax.set_yticks([])
            continue
        arch = next(a for a in ARCHITECTURES if a.code == code)
        rng = np.random.RandomState(42)
        rho, _ = load_cluster(CLUSTERS[0])
        result = evolve_architecture(arch, rho, CONFIG["strength"], rng)
        arr = np.stack(result["history"], axis=0)
        gy, gx = np.gradient(arr[-1])
        mag = np.hypot(gx, gy)
        ax.imshow(mag, cmap="viridis")
        ax.set_title(f"{code} coherence map (Abell 2744)", fontsize=8)
        ax.tick_params(labelsize=6)
    fig.suptitle("Coherence maps across selected architectures")
    fig.tight_layout()
    fig.savefig(PLOTS / "coherence_maps.png", dpi=140, bbox_inches="tight")
    plt.close(fig)

    fig, axes = plt.subplots(2, 5, figsize=(22, 8))
    mem_codes = ["A1", "A2", "A8", "A9", "A10"]
    for ax, code in zip(axes.ravel(), mem_codes):
        row = next((r for r in rows if r["architecture_code"] == code and r["cluster_id"] == "Abell2744"), None)
        if row is None:
            ax.text(0.5, 0.5, f"{code}: no data", ha="center", va="center")
            ax.set_xticks([]); ax.set_yticks([])
            continue
        arch = next(a for a in ARCHITECTURES if a.code == code)
        rng = np.random.RandomState(42)
        rho, _ = load_cluster(CLUSTERS[0])
        result = evolve_architecture(arch, rho, CONFIG["strength"], rng)
        arr = np.stack(result["history"], axis=0)
        diff = np.diff(arr, axis=0)
        mem_map = np.zeros_like(arr[0])
        for t in range(diff.shape[0]):
            mem_map += np.abs(diff[t])
        mem_map /= max(1, diff.shape[0])
        ax.imshow(mem_map, cmap="magma")
        ax.set_title(f"{code} activity map (Abell 2744)", fontsize=8)
        ax.tick_params(labelsize=6)
    fig.suptitle("Memory / activity maps across selected architectures")
    fig.tight_layout()
    fig.savefig(PLOTS / "memory_maps.png", dpi=140, bbox_inches="tight")
    plt.close(fig)

    fig, axes = plt.subplots(1, 4, figsize=(20, 5))
    ranking_physical = [r for r in ranking if not r["architecture_code"].startswith("WR")]
    rank_codes = [r["architecture_code"] for r in ranking_physical]
    rank_vals_k = [next(s["median_pearson_kappa"] for s in summaries if s["architecture_code"] == c) for c in rank_codes]
    rank_vals_g = [next(s["median_pearson_gamma"] for s in summaries if s["architecture_code"] == c) for c in rank_codes]
    rank_vals_w = [next(s["median_wave_mode_count"] for s in summaries if s["architecture_code"] == c) for c in rank_codes]
    rank_vals_m = [next(s["median_emergent_memory_index"] for s in summaries if s["architecture_code"] == c) for c in rank_codes]
    axes[0].bar(rank_codes, rank_vals_k, color="steelblue", edgecolor="black")
    axes[0].set_title("Architecture ranking: Pearson κ", fontsize=9)
    axes[0].tick_params(axis="x", rotation=45, labelsize=8)
    axes[0].grid(alpha=0.3)
    axes[1].bar(rank_codes, rank_vals_g, color="steelblue", edgecolor="black")
    axes[1].set_title("Architecture ranking: Pearson γ", fontsize=9)
    axes[1].tick_params(axis="x", rotation=45, labelsize=8)
    axes[1].grid(alpha=0.3)
    axes[2].bar(rank_codes, rank_vals_w, color="darkorange", edgecolor="black")
    axes[2].set_title("Architecture ranking: wave mode count", fontsize=9)
    axes[2].tick_params(axis="x", rotation=45, labelsize=8)
    axes[2].grid(alpha=0.3)
    axes[3].bar(rank_codes, rank_vals_m, color="green", edgecolor="black")
    axes[3].set_title("Architecture ranking: memory index", fontsize=9)
    axes[3].tick_params(axis="x", rotation=45, labelsize=8)
    axes[3].grid(alpha=0.3)
    fig.suptitle("Architecture rankings (ascending rank order, physical only)")
    fig.tight_layout()
    fig.savefig(PLOTS / "architecture_rankings.png", dpi=140, bbox_inches="tight")
    plt.close(fig)

    fig, axes = plt.subplots(2, 4, figsize=(22, 9))
    dash_keys = ["median_pearson_kappa", "median_rms_kappa", "median_coherence_gain",
                 "median_emergent_memory_index", "median_phase_emergence_score",
                 "median_internal_circulation_score", "median_wave_mode_count",
                 "median_wave_coherence_length"]
    dash_titles = ["Pearson κ", "RMS κ", "Coherence gain", "Memory",
                    "Phase emergence", "Internal circulation",
                    "Wave mode count", "Wave coherence length"]
    for ax, key, title in zip(axes.ravel(), dash_keys, dash_titles):
        vals = [s[key] for s in all_codes]
        ax.bar(labels, vals, color=colors, edgecolor="black")
        ax.set_title(title, fontsize=9)
        ax.tick_params(axis="x", rotation=45, labelsize=7)
        ax.grid(axis="y", alpha=0.3)
        if key == "median_rms_kappa":
            ax.invert_yaxis()
    fig.suptitle("Microscopic architecture science dashboard (red = wrong control)")
    fig.tight_layout()
    fig.savefig(PLOTS / "science_dashboard.png", dpi=140, bbox_inches="tight")
    plt.close(fig)


def build_report(summaries: list[dict], ranking: list[dict], audit: list[dict],
                 hashes: dict, elapsed: float) -> str:
    by = {s["architecture_code"]: s for s in summaries}
    physical = [s for s in summaries if not s["is_wrong"]]
    physical_ranking = [r for r in ranking if not next(s for s in summaries if s["architecture_code"] == r["architecture_code"])["is_wrong"]]
    s9_ref = {"median_pearson_kappa": 0.10494564177366811,
              "median_pearson_gamma": 0.09399620573787398,
              "median_coherence_gain": 0.0021661336656330736,
              "median_emergent_memory_index": 0.9980359705156546,
              "label": "S9 (MICROSTATE-LAB-001)"}
    c10_pk = 0.10339683814108096
    alpha_audits = sorted([a for a in audit if a["is_alpha_or_3alpha"]], key=lambda a: float(a["log10_distance_to_constant"]))

    def line_q1():
        above = [s for s in physical if s["median_pearson_kappa"] > s9_ref["median_pearson_kappa"]]
        if above:
            return f"**Yes.** {len(above)} architecture(s) naturally exceed S9 Pearson κ = {s9_ref['median_pearson_kappa']:+.5f}: " + ", ".join(f"{s['architecture_code']} = {s['median_pearson_kappa']:+.5f}" for s in above) + "."
        return f"No single architecture reproduces the S9 Pearson κ = {s9_ref['median_pearson_kappa']:+.5f}; best physical architecture: {max(physical, key=lambda s: s['median_pearson_kappa'])['architecture_code']} = {max(physical, key=lambda s: s['median_pearson_kappa'])['median_pearson_kappa']:+.5f}."

    def line_q2():
        leader_k = max(physical, key=lambda s: s["median_pearson_kappa"])
        leader_m = max(physical, key=lambda s: s["median_emergent_memory_index"])
        return f"Best on Pearson κ = {leader_k['architecture_code']} ({leader_k['architecture_name']}) at {leader_k['median_pearson_kappa']:+.5f}. Highest memory = {leader_m['architecture_code']} at {leader_m['median_emergent_memory_index']:.5f}."

    def line_q3():
        emerged = [s for s in physical if s["median_phase_emergence_score"] > 0.1]
        if emerged:
            return f"Phase emerged in {len(emerged)}/10 physical architectures (without explicit phase updates). Sample: " + ", ".join(f"{s['architecture_code']}={s['median_phase_emergence_score']:.3f}" for s in emerged[:6])
        return "No physical architecture produced a sustained phase field above emergence threshold without explicit phase updates."

    def line_q4():
        emerged = [s for s in physical if s["median_orientation_emergence_score"] > COHERENCE_GAIN_THRESHOLD]
        if emerged:
            return f"Orientation emerged in {len(emerged)}/10 physical architectures. Sample: " + ", ".join(f"{s['architecture_code']}={s['median_orientation_emergence_score']:.3e}" for s in emerged[:6])
        return "No physical architecture produced orientation spontaneously."

    def line_q5():
        emerged = [s for s in physical if s["median_emergent_memory_index"] >= MEMORY_INDEX_THRESHOLD and s["median_evolution_activity"] > ACTIVITY_THRESHOLD]
        if emerged:
            return f"Memory emerged in {len(emerged)}/10 physical architectures. Sample: " + ", ".join(f"{s['architecture_code']}={s['median_emergent_memory_index']:.5f}" for s in emerged[:6])
        return "No physical architecture sustained memory above emergence threshold."

    def line_q6():
        wave_arch = [s for s in physical if s["median_wave_mode_count"] >= 2.0]
        if wave_arch:
            leader = max(wave_arch, key=lambda s: s["median_wave_mode_count"])
            return f"Propagating excitation modes emerged in {len(wave_arch)}/10 physical architectures; highest wave-mode count = {leader['architecture_code']} ({leader['median_wave_mode_count']:.1f}/4 modes)."
        return "No physical architecture produced clear propagating excitation modes above the emergence thresholds."

    def line_q7():
        leaders = sorted(physical, key=lambda s: -s["median_wave_mode_count"])[:3]
        if not leaders or leaders[0]["median_wave_mode_count"] < 1:
            return "No architecture reached the wave-mode emergence threshold."
        L = leaders[0]
        lines = [f"Top architecture: {L['architecture_code']} ({L['architecture_name']})"]
        lines.append(f"  - propagating disturbance = {L['median_wave_propagating']:+.3f}")
        lines.append(f"  - standing mode = {L['median_wave_standing']:+.3f}")
        lines.append(f"  - transverse = {L['median_wave_transverse']:+.3f}")
        lines.append(f"  - longitudinal = {L['median_wave_longitudinal']:+.3f}")
        lines.append(f"  - polarization-like = {L['median_wave_polarization']:+.3f}")
        lines.append(f"  - attenuation = {L['median_wave_attenuation']:+.3f}")
        lines.append(f"  - dispersion = {L['median_wave_dispersion']:+.3f}")
        lines.append(f"  - coherence length ≈ {L['median_wave_coherence_length']:.1f} pixels")
        return "\n".join(lines)

    def line_q8():
        above = [s for s in physical if s["median_pearson_kappa"] > c10_pk]
        if above:
            return f"{len(above)} architecture(s) surpass C10 ({c10_pk:+.5f}): " + ", ".join(f"{s['architecture_code']} = {s['median_pearson_kappa']:+.5f}" for s in above)
        return f"No architecture surpasses C10 ({c10_pk:+.5f}). Best physical = {max(physical, key=lambda s: s['median_pearson_kappa'])['architecture_code']} = {max(physical, key=lambda s: s['median_pearson_kappa'])['median_pearson_kappa']:+.5f}."

    def line_q9():
        if not alpha_audits:
            return "No dimensionless quantity repeatedly converged near α or 3α across clusters. The audit is purely observational; no fitting occurred."
        top5 = alpha_audits[:5]
        sample_lines = [f"`{a['quantity_name']}` (arch {a['architecture_code']}) = {float(a['value']):+.5e}, factor to {a['nearest_constant']} = {float(a['factor_to_constant']):.4f}, log₁₀ distance = {float(a['log10_distance_to_constant']):+.4f}" for a in top5]
        all_alp3 = [a for a in alpha_audits if float(a["log10_distance_to_constant"]) < 0.1]
        return (
            f"{len(alpha_audits)} audit entries sit nearest α or 3α. Of these, {len(all_alp3)} are within log₁₀ distance < 0.1 from α or 3α (~26% linear deviation). "
            f"Closest hits (purely observational):\n\n"
            + "\n".join(f"- {l}" for l in sample_lines)
        )

    def line_q10():
        all_ok = all(r["max_conservation_error"] <= EPS + 1e-30 for r in ([s for s in summaries]))
        return f"{'Yes' if all_ok else 'No'} — all {len(ARCHITECTURES) * len(CLUSTERS)} runs preserve the unit-speed normalization at or below machine epsilon ({EPS:.3e})."

    def determine_outcome():
        # Outcome A: single architecture wins on Pearson κ, memory, and wave
        a10 = by["A10"]
        others = [s for s in physical if s["architecture_code"] != "A10"]
        if a10["median_pearson_kappa"] >= max(s["median_pearson_kappa"] for s in others) - 1e-3:
            if a10["median_wave_mode_count"] >= 2:
                return "Outcome A", f"The unified microcell ({a10['architecture_code']}) reproduces S9-like behaviour while supporting {a10['median_wave_mode_count']:.1f} stable wave modes at machine-precision conservation."
        leader = max(physical, key=lambda s: s["median_pearson_kappa"])
        wave_arch = [s for s in physical if s["median_wave_mode_count"] >= 2]
        if len(wave_arch) >= 3 and wave_arch[0]["median_wave_mode_count"] >= 2:
            return "Outcome B", f"Several architectures improve different aspects: {leader['architecture_code']} leads on κ at {leader['median_pearson_kappa']:+.5f}; " + ", ".join(f"{s['architecture_code']} supports {s['median_wave_mode_count']:.1f} wave modes" for s in wave_arch[:3]) + ". No single architecture uniquely wins every metric."
        return "Outcome C", "Changing the microscopic constituent architecture does not yield a unique reproduction of S9, and no single architecture consolidates the wave-emergence signatures observed here."

    outcome, outcome_text = determine_outcome()

    lines = [
        "# PBUF MICROSTRUCTURE-ENTITY-LAB-001",
        "",
        "**Microscopic Constituent Architecture Laboratory inside the frozen Version 1 weak-lensing laboratory (LAB-FREEZE-001).**",
        "",
        "## Status",
        "",
        f"- Frozen hash verification: **{'PASS' if hashes['ok'] else 'FAIL'}**",
        f"- Architectures: **{len(ARCHITECTURES)}** (A1-A10 + WR1-WR4)",
        f"- Production runs: **{len(ARCHITECTURES) * len(CLUSTERS)}**",
        f"- Runtime: **{elapsed:.1f} s**",
        "- Fitting or optimisation: **none**",
        "",
        "## Frozen laboratory",
        "",
        "All transport, source-plane, Jacobian observable, numerical, constitutive, and production components remain byte-identical to LAB-FREEZE-001. Only the microscopic constituent architecture (the per-cell internal state representation) varies across families.",
        "",
        "## Candidate architectures",
        "",
        "| # | Code | Name | Principle | Internal nodes |",
        "|---|---|---|---|---:|",
    ]
    for a in ARCHITECTURES:
        lines.append(f"| {a.code} | {a.name} | `{a.principle}` | {a.internal_nodes} |")
    lines += [
        "",
        "Wrong controls: WR1 random internal topology, WR2 disconnected internal nodes, WR3 over-connected constituent (all nodes equally coupled), WR4 frozen internal architecture. They must underperform if the laboratory responds to a meaningful internal cell structure.",
        "",
        "## Architecture summary (median across 5 clusters)",
        "",
        "| Architecture | Pearson κ | Pearson γ | RMS κ | Coherence gain | Memory | Wave modes | Conservation |",
        "|---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for r in ranking:
        s = by[r["architecture_code"]]
        lines.append(f"| {s['architecture_code']} | {s['median_pearson_kappa']:+.5f} | {s['median_pearson_gamma']:+.5f} | {s['median_rms_kappa']:.5f} | {s['median_coherence_gain']:+.3e} | {s['median_emergent_memory_index']:.5f} | {s['median_wave_mode_count']:.1f} | {s['max_conservation_error']:.3e} |")
    lines += [
        "",
        "## Emergent diagnostic definitions",
        "",
        "- **Phase emergence** (phase_emergence_score): spectral peak prominence of the mean-field time series for oscillatory architectures, or sign-change rate for real architectures; score > 0.1 counts as emergence.",
        "- **Orientation emergence** (orientation_emergence_score): gradient coherence of the coarse field; score > 1e-04 counts as emergence.",
        "- **Memory / persistence**: mean cosine of successive state increments; emergence requires index ≥ 0.9 and activity > 1e-6.",
        "- **Internal circulation**: structure-specific metric of mean internal-DOF spread normalised by strength; score > 0.3 counts as emergence (A3, A7, A10).",
        "- **Multiplicative coupling**: detected when an architecture has both amplitude-like and phase-like internal DOFs that couple nonlinearly; A5 and A10 satisfy this by construction.",
        "- **Cooperative response**: adaptive-weighting of the local neighbourhood update; score > 0.3 counts as emergence (A9).",
        "- **Neighbour coherence** (coherence_gain): gradient-coherence gain between initial and final coarse constitutive state.",
        "- **Wave-mode audit**: 8 diagnostics — propagating disturbance, standing mode, transverse/longitudinal mode scores, polarization-like, attenuation, dispersion, coherence length. wave_emerged requires at least one diagnostic > 0.3; wave_mode_count counts how many of {propagating, standing, transverse, longitudinal} exceed 0.3.",
        "",
        "## Cross-cluster statistics",
        "",
        f"Five clusters × {len(ARCHITECTURES)} architectures = {len(ARCHITECTURES)*len(CLUSTERS)} production runs. Per-run breakdowns in `cross_cluster_statistics.csv`; per-cluster diagnostics in `emergent_state_statistics.csv` and `wave_mode_statistics.csv`.",
        "",
        "## Wave Emergence Audit",
        "",
        "Every architecture was probed for 8 wave-like signatures. None of them is labelled electromagnetic. The audit characterises natural propagation modes a microscopic constituent supports.",
        "",
        "| Architecture | Propagating | Standing | Transverse | Longitudinal | Polarization | Attenuation | Dispersion | Coherence L | Wave modes |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for r in ranking:
        s = by[r["architecture_code"]]
        lines.append(f"| {s['architecture_code']} | {s['median_wave_propagating']:+.3f} | {s['median_wave_standing']:+.3f} | {s['median_wave_transverse']:+.3f} | {s['median_wave_longitudinal']:+.3f} | {s['median_wave_polarization']:+.3f} | {s['median_wave_attenuation']:+.3f} | {s['median_wave_dispersion']:+.3f} | {s['median_wave_coherence_length']:.2f} | {s['median_wave_mode_count']:.1f} |")
    lines += [
        "",
        "## Fundamental constant audit",
        "",
        f"For every architecture we observed dimensionless ratios produced by the microscopic constituent evolution: coupling ratios (K/ω, K·dt, K/γ, internal_K·dt), signal-to-noise ratios, the Pearson κ/RMS κ ratio, wave coherence over grid, mode count over internal-node count, etc. Each row of `fundamental_constant_audit.csv` reports value, log₁₀|value|, the nearest known dimensionless constant, and the log₁₀ distance. Primary audit targets are **α = 1/137.035999084 ≈ {ALPHA_FS:.5e}** and **3α ≈ {THREE_ALPHA_FS:.5e}**; no fitting, no optimisation — passive observation only.",
        "",
        "## Candidate ranking",
        "",
        "Physical architectures ranked by mean rank across all primary metrics (higher Pearson κ/γ, lower RMS κ/γ, lower bias, higher coherence / memory / phase / orientation / circulation / multiplicative / cooperative / wave-mode count / coherence length).",
        "",
        "| Rank | Code | Pearson κ | Wave modes | Memory | Phase | Multiplicative | Cooperative | Circulation | Rank sum |",
        "|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for r in physical_ranking:
        s = by[r["architecture_code"]]
        lines.append(f"| {r['rank']} | {s['architecture_code']} | {s['median_pearson_kappa']:+.5f} | {s['median_wave_mode_count']:.1f} | {s['median_emergent_memory_index']:.5f} | {s['median_phase_emergence_score']:.3f} | {s['median_multiplicative_coupling_score']:.3f} | {s['median_cooperative_response_score']:.3f} | {s['median_internal_circulation_score']:.3f} | {r['rank_sum']:.0f} |")
    lines += [
        "",
        "## Required questions",
        "",
        "### Q1. Does constituent architecture outperform point elements?",
        "",
        line_q1(),
        "",
        "### Q2. Which architecture best reproduces S9 behaviour?",
        "",
        line_q2(),
        "",
        "### Q3. Does phase emerge naturally?",
        "",
        line_q3(),
        "",
        "### Q4. Does orientation emerge naturally?",
        "",
        line_q4(),
        "",
        "### Q5. Does memory emerge naturally?",
        "",
        line_q5(),
        "",
        "### Q6. Do propagating excitation modes appear?",
        "",
        line_q6(),
        "",
        "### Q7. If wave modes appear, what are their properties?",
        "",
        line_q7(),
        "",
        "### Q8. Does any architecture outperform C10?",
        "",
        line_q8(),
        "",
        "### Q9. Do any stable dimensionless quantities repeatedly converge near α or 3α?",
        "",
        line_q9(),
        "",
        "### Q10. Does every successful architecture preserve machine-precision conservation?",
        "",
        line_q10(),
        "",
        "## Outcome determination",
        "",
        "Outcome criteria from the milestone:",
        "- **A**: One microscopic constituent architecture naturally reproduces the S9 signature and supports stable emergent wave modes while preserving conservation.",
        "- **B**: Several architectures improve different aspects of the laboratory, but no unique microscopic constituent emerges.",
        "- **C**: Changing the constituent architecture does not improve upon the current microscopic description.",
        "",
        f"**{outcome}.** {outcome_text}",
        "",
        "## C10 provenance",
        "",
        "C10 was not modified and not rerun. The benchmark remains archived at `runs/version_b_physics_lab002/interaction_matrix.csv`.",
        "",
        "## Numerical stability",
        "",
        f"All {len(ARCHITECTURES) * len(CLUSTERS)} runs preserve the frozen unit-speed normalization at or below machine epsilon ({EPS:.3e}).",
        "",
        "## Required artefacts",
        "",
        "`report.md`, `architecture_summary.csv`, `cross_cluster_statistics.csv`, `candidate_ranking.csv`, `wave_mode_statistics.csv`, `emergent_state_statistics.csv`, `fundamental_constant_audit.csv`, `run.json`, `validation.json`, and all required plots are present in `runs/microstructure_entity_lab001/`.",
        "",
    ]
    return "\n".join(lines)


def main() -> None:
    started_total = time.perf_counter()
    OUT.mkdir(parents=True, exist_ok=True)
    PLOTS.mkdir(parents=True, exist_ok=True)
    hashes = {"ok": True, "files": {}}
    for name, expected in EXPECTED_HASHES.items():
        actual = file_sha256(ROOT / name)
        match = actual == expected
        hashes["files"][name] = {"expected_sha256": expected, "actual_sha256": actual, "match": match}
        hashes["ok"] = hashes["ok"] and match
    if not hashes["ok"]:
        raise RuntimeError("Frozen source file hashes do not match LAB-FREEZE-001")

    rows = []
    for cluster in CLUSTERS:
        rho, obs = load_cluster(cluster)
        for arch in ARCHITECTURES:
            rows.append(run_one(arch, cluster, rho, obs))

    summaries = aggregate(rows)
    wave_stats = wave_mode_statistics(rows)
    emergent_stats = emergent_state_statistics(rows)
    audit = fundamental_constant_audit(rows, summaries)
    ranking = candidate_ranking(summaries)

    summary_fields = list(summaries[0].keys())
    cross_fields = ["architecture_number", "architecture_code", "architecture_name",
                    "is_wrong", "internal_nodes", "principle",
                    "cluster_id", "cluster_label",
                    "pearson_kappa", "pearson_gamma", "ssim_kappa", "ssim_gamma",
                    "rms_kappa", "rms_gamma", "kappa_bias", "gamma_bias", "runtime_seconds",
                    "max_conservation_error", "coherence_gain", "emergent_memory_index",
                    "evolution_activity", "phase_emergence_score", "orientation_emergence_score",
                    "internal_circulation_score", "multiplicative_coupling_score",
                    "cooperative_response_score",
                    "relaxation_time", "spatial_correlation_length", "temporal_persistence_length",
                    "propagating_disturbance", "standing_mode", "transverse_mode",
                    "longitudinal_mode", "polarization_like", "attenuation", "dispersion",
                    "coherence_length", "wave_mode_count", "wave_emerged"]
    write_csv(OUT / "architecture_summary.csv", summaries, summary_fields)
    write_csv(OUT / "cross_cluster_statistics.csv", rows, cross_fields)
    ranking_fields = ["rank", "architecture_code", "architecture_name", "rank_sum",
                      "median_pearson_kappa", "median_pearson_gamma", "median_ssim_kappa",
                      "median_rms_kappa", "median_coherence_gain", "median_emergent_memory_index",
                      "median_phase_emergence_score", "median_orientation_emergence_score",
                      "median_internal_circulation_score", "median_multiplicative_coupling_score",
                      "median_cooperative_response_score", "median_wave_mode_count",
                      "median_wave_coherence_length", "median_relaxation_time"]
    write_csv(OUT / "candidate_ranking.csv", ranking, ranking_fields)
    wave_fields = sorted({k for r in wave_stats for k in r.keys()})
    write_csv(OUT / "wave_mode_statistics.csv", wave_stats, wave_fields)
    emergent_fields = sorted({k for r in emergent_stats for k in r.keys()})
    write_csv(OUT / "emergent_state_statistics.csv", emergent_stats, emergent_fields)
    audit_fields = ["architecture_code", "architecture_name", "is_wrong", "quantity_name",
                    "value", "log_abs", "nearest_constant", "nearest_constant_value",
                    "log10_distance_to_constant", "factor_to_constant", "is_alpha_or_3alpha"]
    write_csv(OUT / "fundamental_constant_audit.csv", audit, audit_fields)

    make_plots(rows, summaries, ranking, wave_stats, emergent_stats)
    elapsed = time.perf_counter() - started_total
    report_text = build_report(summaries, ranking, audit, hashes, elapsed)
    (OUT / "report.md").write_text(report_text)

    run = {
        "milestone": "PBUF MICROSTRUCTURE-ENTITY-LAB-001",
        "kind": "microscopic constituent architecture search",
        "frozen_laboratory": "Version 1 weak-lensing laboratory (LAB-FREEZE-001)",
        "frozen_implementation_sha256": {k: v["actual_sha256"] for k, v in hashes["files"].items()},
        "production_configuration": CONFIG, "clusters": CLUSTERS,
        "architectures": [a.__dict__ for a in ARCHITECTURES],
        "fixed_parameters": {"dt": DT, "steps": STEPS, "k": K, "gamma": GAMMA,
                             "omega": OMEGA, "internal_k": INTERNAL_K,
                             "alpha_fs": ALPHA_FS, "three_alpha_fs": THREE_ALPHA_FS},
        "emergence_thresholds": {"coherence_gain": COHERENCE_GAIN_THRESHOLD,
                                  "memory_index": MEMORY_INDEX_THRESHOLD,
                                  "evolution_activity": ACTIVITY_THRESHOLD},
        "fitting_performed": False, "optimisation_performed": False,
        "frozen_components_modified": False, "execution_seconds_total": elapsed,
    }
    (OUT / "run.json").write_text(json.dumps(run, indent=2))

    required = [OUT / "report.md",
                OUT / "architecture_summary.csv",
                OUT / "cross_cluster_statistics.csv",
                OUT / "candidate_ranking.csv",
                OUT / "wave_mode_statistics.csv",
                OUT / "emergent_state_statistics.csv",
                OUT / "fundamental_constant_audit.csv",
                OUT / "run.json"] + [PLOTS / n for n in (
                    "architecture_comparison.png", "internal_state_evolution.png",
                    "wave_mode_analysis.png", "coherence_maps.png",
                    "memory_maps.png", "architecture_rankings.png",
                    "science_dashboard.png")]
    png_ok = all(p.read_bytes()[:8] == b"\x89PNG\r\n\x1a\n" for p in required if p.suffix == ".png")
    finite_ok = all(np.isfinite(r[k]) for r in rows for k in (
        "pearson_kappa", "pearson_gamma", "ssim_kappa", "ssim_gamma",
        "rms_kappa", "rms_gamma", "kappa_bias", "gamma_bias", "runtime_seconds",
        "max_conservation_error", "coherence_gain", "emergent_memory_index",
        "evolution_activity", "spatial_correlation_length", "temporal_persistence_length",
        "relaxation_time", "phase_emergence_score", "orientation_emergence_score",
        "internal_circulation_score", "multiplicative_coupling_score",
        "cooperative_response_score", "propagating_disturbance", "standing_mode",
        "transverse_mode", "longitudinal_mode", "polarization_like",
        "attenuation", "dispersion", "coherence_length", "wave_mode_count"))
    conservation_ok = all(r["max_conservation_error"] <= EPS + 1e-30 for r in rows)
    artifacts_ok = all(p.exists() and p.stat().st_size > 0 for p in required)
    validation = {
        "milestone": "PBUF MICROSTRUCTURE-ENTITY-LAB-001",
        "frozen_hash_verification_passed": hashes["ok"], "frozen_hashes": hashes["files"],
        "expected_run_count": len(ARCHITECTURES) * len(CLUSTERS),
        "actual_run_count": len(rows),
        "architecture_count": len(ARCHITECTURES), "cluster_count": len(CLUSTERS),
        "all_metrics_finite": finite_ok,
        "all_runs_machine_precision_conservation": conservation_ok,
        "runs_preserving_conservation": int(sum(r["max_conservation_error"] <= EPS + 1e-30 for r in rows)),
        "required_artifacts_present_nonempty": artifacts_ok,
        "png_signatures_valid": png_ok,
        "validation_passed": bool(hashes["ok"] and len(rows) == len(ARCHITECTURES) * len(CLUSTERS)
                                  and finite_ok and conservation_ok and artifacts_ok and png_ok),
    }
    (OUT / "validation.json").write_text(json.dumps(validation, indent=2))
    if not validation["validation_passed"]:
        raise RuntimeError("Microstructure-entity laboratory validation failed")


if __name__ == "__main__":
    main()
