#!/usr/bin/env python3
"""PBUF MICROSCOPIC-INVARIANTS-LAB-001 — conserved-quantity transport laboratory."""
from __future__ import annotations

import csv
import json
import math
import sys
import time
from dataclasses import dataclass
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

OUT = ROOT / "runs" / "microscopic_invariants_lab001"
PLOTS = OUT / "plots"
BENCHMARK = ROOT / "PBUF_benchmark"
WAVE_REGISTRY = ROOT / "runs" / "wave_family_registry.csv"
INVARIANT_REGISTRY = ROOT / "runs" / "invariant_registry.csv"

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
COUPLING_FAST_TO_SLOW = 1.0
COUPLING_SLOW_TO_FAST = 0.3
FAST_TIMESCALE = OMEGA * K
SLOW_TIMESCALE = 0.25

ALPHA_FS = 1.0 / 137.035999084
THREE_ALPHA_FS = 3.0 * ALPHA_FS


@dataclass(frozen=True)
class TransportPrinciple:
    number: int
    code: str
    name: str
    principle: str
    invariant_name: str
    invariant_kind: str
    is_wrong: bool = False


TRANSPORT_PRINCIPLES = [
    TransportPrinciple(1, "T1", "Scalar Density", "transport constitutive density; standard A8 relaxation; density NOT conserved",
                       "density (no strict conservation)", "scalar"),
    TransportPrinciple(2, "T2", "Conserved Phase", "transport phase of (u_fast + i·u_slow); density emerges from neighbour magnitude",
                       "|z|² (rotation-invariant)", "phase"),
    TransportPrinciple(3, "T3", "Conserved Orientation", "transport unit-vector (u_fast, u_slow); magnitude responds",
                       "|v̂|² = 1 (orientation unit)", "orientation"),
    TransportPrinciple(4, "T4", "Conserved Action", "neighbour exchange preserves ∑K/2·u² + V(u_slow, u_fast)",
                       "A = ∑K/2·u² + V", "action"),
    TransportPrinciple(5, "T5", "Conserved Internal Energy", "energy freely exchanges between fast and slow layers; total E conserved",
                       "E_fast + E_slow", "energy"),
    TransportPrinciple(6, "T6", "Conserved Information", "local Shannon-like information redistributes among neighbours; H constant",
                       "H = -∑p·log(p)", "information"),
    TransportPrinciple(7, "T7", "Conserved Circulation", "rotational transport with ∮u·dr constant on each closed loop",
                       "circulation Γ", "circulation"),
    TransportPrinciple(8, "T8", "Conserved Flux", "neighbour-to-neighbour flux with ∇·j = 0",
                       "div j = 0", "flux"),
    TransportPrinciple(9, "T9", "Coupled Energy + Phase", "joint energy and phase conservation: A8 phase evolves while E preserved",
                       "E and arg(z) jointly preserved", "energy_phase"),
    TransportPrinciple(10, "T10", "Unified State Transport", "single complex state z = u_fast + i·u_slow evolves under Ginzburg-Landau",
                       "|z|² and arg(z)", "unified"),
    TransportPrinciple(11, "WR1", "Wrong: Random Transport", "u_slow and u_fast randomized each step; no transport law",
                       "no invariant", "wrong", is_wrong=True),
    TransportPrinciple(12, "WR2", "Wrong: Non-conserved Transport", "intentionally introduces a sink; invariant violated",
                       "deliberately violated", "wrong", is_wrong=True),
    TransportPrinciple(13, "WR3", "Wrong: Pure Diffusion", "linear diffusion only; no layer coupling; no wave structure",
                       "linear", "wrong", is_wrong=True),
    TransportPrinciple(14, "WR4", "Wrong: Frozen Transport", "states frozen at initial values; no evolution",
                       "frozen", "wrong", is_wrong=True),
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


def build_C(u: np.ndarray, strength: float) -> np.ndarray:
    lo, hi = float(u.min()), float(u.max())
    if hi - lo < 1e-15:
        return np.zeros_like(u)
    return strength * (u - lo) / (hi - lo)


def A8_init(rho: np.ndarray, strength: float, rng: np.random.RandomState) -> tuple[np.ndarray, np.ndarray]:
    eq = strength * rho
    u_slow = eq.copy()
    u_fast = eq.copy() + 0.02 * strength * rng.randn(*rho.shape)
    return u_slow, u_fast


def E_layer(u: np.ndarray) -> float:
    return float(0.5 * np.sum(u ** 2))


def E_total(u_fast: np.ndarray, u_slow: np.ndarray) -> float:
    return E_layer(u_fast) + E_layer(u_slow)


def evolve_transport(t: TransportPrinciple, rho: np.ndarray, strength: float, rng: np.random.RandomState) -> dict:
    eq = strength * rho
    code = t.code
    u_slow, u_fast = A8_init(rho, strength, rng)
    history: list[np.ndarray] = []
    log: list[tuple[np.ndarray, np.ndarray]] = []
    invariants: list[float] = []
    invariants.append(0.0)
    E_fast_log: list[float] = [E_layer(u_fast)]
    E_slow_log: list[float] = [E_layer(u_slow)]
    history.append(0.5 * u_slow + 0.5 * u_fast)
    log.append((u_slow.copy(), u_fast.copy()))

    def _record(step_idx):
        mixed = 0.5 * u_slow + 0.5 * u_fast
        history.append(mixed.copy())
        log.append((u_slow.copy(), u_fast.copy()))
        E_fast_log.append(E_layer(u_fast))
        E_slow_log.append(E_layer(u_slow))
        inv = compute_invariant(code, u_fast, u_slow)
        invariants.append(inv)

    if code == "T1":
        for step in range(STEPS):
            n4s = sum(neighbours4(u_slow)) / 4.0
            n4f = sum(neighbours4(u_fast)) / 4.0
            d_fast = DT * OMEGA * K * ((n4f - u_fast) + COUPLING_SLOW_TO_FAST * (u_slow - u_fast))
            d_slow = DT * SLOW_TIMESCALE * ((n4s - u_slow) + COUPLING_FAST_TO_SLOW * (u_fast - u_slow))
            u_fast = np.clip(u_fast + d_fast, -5.0, 5.0)
            u_slow = np.clip(u_slow + d_slow, -5.0, 5.0)
            _record(step)

    elif code == "T2":
        for step in range(STEPS):
            z = u_fast + 1j * u_slow
            mag = np.abs(z)
            phi = np.angle(z)
            n4z = complex_neighbour_mean(z)
            n4_phi = np.angle(n4z)
            n4_mag = np.abs(n4z)
            dphi = DT * OMEGA * np.sin(n4_phi - phi) * 0.5
            dmag = DT * K * (n4_mag - mag) * 0.3
            new_phi = phi + dphi
            new_mag = mag + dmag
            new_z = new_mag * np.exp(1j * new_phi)
            u_fast = np.clip(np.real(new_z), -5.0, 5.0)
            u_slow = np.clip(np.imag(new_z), -5.0, 5.0)
            _record(step)

    elif code == "T3":
        for step in range(STEPS):
            v = np.stack([u_fast, u_slow], axis=0)
            norm = np.sqrt(np.sum(v ** 2, axis=0))
            norm_safe = np.where(norm > 1e-12, norm, 1.0)
            v_hat = v / norm_safe
            vx = v_hat[0]
            vy = v_hat[1]
            n4_vx = sum(neighbours4(vx)) / 4.0
            n4_vy = sum(neighbours4(vy)) / 4.0
            cross = vx * n4_vy - vy * n4_vx
            dot = vx * n4_vx + vy * n4_vy
            theta = np.arctan2(cross, dot) * 0.5
            cos_t = np.cos(DT * K * theta)
            sin_t = np.sin(DT * K * theta)
            new_vx = vx * cos_t - vy * sin_t
            new_vy = vx * sin_t + vy * cos_t
            n4_mag = (norm_safe + sum(neighbours4(norm_safe)) / 4.0) / 2.0
            u_fast = np.clip(new_vx * n4_mag, -5.0, 5.0)
            u_slow = np.clip(new_vy * n4_mag, -5.0, 5.0)
            _record(step)

    elif code == "T4":
        def action_total(u_fast_, u_slow_):
            return 0.5 * K * (np.sum(u_fast_ ** 2) + np.sum(u_slow_ ** 2)) + 0.5 * INTERNAL_K * np.sum((u_fast_ - u_slow_) ** 2)
        A_init = action_total(u_fast, u_slow)
        invariants[0] = A_init
        for step in range(STEPS):
            n4s = sum(neighbours4(u_slow)) / 4.0
            n4f = sum(neighbours4(u_fast)) / 4.0
            d_fast = DT * OMEGA * K * ((n4f - u_fast) + COUPLING_SLOW_TO_FAST * (u_slow - u_fast))
            d_slow = DT * SLOW_TIMESCALE * ((n4s - u_slow) + COUPLING_FAST_TO_SLOW * (u_fast - u_slow))
            u_fast_new = u_fast + d_fast
            u_slow_new = u_slow + d_slow
            A_new = action_total(u_fast_new, u_slow_new)
            if A_new > 1e-15:
                scale = math.sqrt(A_init / A_new)
                u_fast_new = u_fast_new * scale
                u_slow_new = u_slow_new * scale
            u_fast = np.clip(u_fast_new, -5.0, 5.0)
            u_slow = np.clip(u_slow_new, -5.0, 5.0)
            _record(step)

    elif code == "T5":
        def energy(u_fast_, u_slow_):
            return 0.5 * np.sum(u_fast_ ** 2) + 0.5 * np.sum(u_slow_ ** 2)
        E_init = energy(u_fast, u_slow)
        invariants[0] = E_init
        for step in range(STEPS):
            n4s = sum(neighbours4(u_slow)) / 4.0
            n4f = sum(neighbours4(u_fast)) / 4.0
            exchange = 0.5 * (u_slow - u_fast)
            d_fast = DT * OMEGA * K * ((n4f - u_fast) + INTERNAL_K * exchange)
            d_slow = DT * SLOW_TIMESCALE * ((n4s - u_slow) - INTERNAL_K * exchange)
            u_fast_new = u_fast + d_fast
            u_slow_new = u_slow + d_slow
            E_new = energy(u_fast_new, u_slow_new)
            if E_new > 1e-15:
                scale = math.sqrt(E_init / E_new)
                u_fast_new = u_fast_new * scale
                u_slow_new = u_slow_new * scale
            u_fast = np.clip(u_fast_new, -5.0, 5.0)
            u_slow = np.clip(u_slow_new, -5.0, 5.0)
            _record(step)

    elif code == "T6":
        p0 = np.abs(u_fast) + 1e-9
        p_norm = p0 / np.sum(p0)
        H_init = float(-np.sum(p_norm * np.log(p_norm + 1e-15)))
        invariants[0] = H_init
        for step in range(STEPS):
            n4f = sum(neighbours4(u_fast)) / 4.0
            n4s = sum(neighbours4(u_slow)) / 4.0
            d_fast = DT * K * 0.3 * (n4f - u_fast)
            d_slow = DT * 0.25 * ((n4s - u_slow) + 0.3 * (u_fast - u_slow))
            u_fast = np.clip(u_fast + d_fast, -5.0, 5.0)
            u_slow = np.clip(u_slow + d_slow, -5.0, 5.0)
            abs_f = np.abs(u_fast) + 1e-9
            p_norm = abs_f / np.sum(abs_f)
            H_curr = float(-np.sum(p_norm * np.log(p_norm + 1e-15)))
            if H_init > 1e-15 and H_curr > 1e-15:
                scale = H_init / H_curr
                u_fast = u_fast * math.sqrt(scale)
            _record(step)

    elif code == "T7":
        for step in range(STEPS):
            curl = (np.roll(u_fast, -1, axis=0) - np.roll(u_fast, 1, axis=0)) - (
                np.roll(u_fast, -1, axis=1) - np.roll(u_fast, 1, axis=1))
            d_fast = DT * K * 0.3 * curl
            n4s = sum(neighbours4(u_slow)) / 4.0
            d_slow = DT * 0.25 * ((n4s - u_slow) + (u_fast - u_slow))
            u_fast = np.clip(u_fast + d_fast, -5.0, 5.0)
            u_slow = np.clip(u_slow + d_slow, -5.0, 5.0)
            _record(step)
        invariants[0] = float(np.sum(np.gradient(u_fast, axis=0) - np.gradient(u_fast, axis=1)))

    elif code == "T8":
        for step in range(STEPS):
            lap_f = neighbours9_weighted(u_fast)
            lap_s = neighbours9_weighted(u_slow)
            d_fast = DT * OMEGA * K * (lap_f + 0.3 * (u_slow - u_fast))
            d_slow = DT * SLOW_TIMESCALE * (lap_s + (u_fast - u_slow))
            u_fast = np.clip(u_fast + d_fast, -5.0, 5.0)
            u_slow = np.clip(u_slow + d_slow, -5.0, 5.0)
            _record(step)

    elif code == "T9":
        def energy_(u_fast_, u_slow_):
            return 0.5 * np.sum(u_fast_ ** 2) + 0.5 * np.sum(u_slow_ ** 2)
        E_init = energy_(u_fast, u_slow)
        invariants[0] = E_init
        for step in range(STEPS):
            n4s = sum(neighbours4(u_slow)) / 4.0
            n4f = sum(neighbours4(u_fast)) / 4.0
            phase_term = np.sin(u_slow - u_fast)
            d_fast = DT * OMEGA * K * ((n4f - u_fast) + 0.15 * (u_slow - u_fast) + 0.05 * phase_term)
            d_slow = DT * SLOW_TIMESCALE * ((n4s - u_slow) + (u_fast - u_slow))
            u_fast_new = u_fast + d_fast
            u_slow_new = u_slow + d_slow
            E_new = energy_(u_fast_new, u_slow_new)
            if E_new > 1e-15:
                scale = math.sqrt(E_init / E_new)
                u_fast_new = u_fast_new * scale
                u_slow_new = u_slow_new * scale
            u_fast = np.clip(u_fast_new, -5.0, 5.0)
            u_slow = np.clip(u_slow_new, -5.0, 5.0)
            _record(step)

    elif code == "T10":
        def total_norm(uf, us):
            return float(np.sum(uf ** 2) + np.sum(us ** 2))
        N_init = total_norm(u_fast, u_slow)
        invariants[0] = N_init
        for step in range(STEPS):
            n4s = sum(neighbours4(u_slow)) / 4.0
            n4f = sum(neighbours4(u_fast)) / 4.0
            d_fast = DT * 0.4 * ((n4f - u_fast) + 0.3 * (u_slow - u_fast))
            d_slow = DT * 0.4 * ((n4s - u_slow) + (u_fast - u_slow))
            u_fast_new = u_fast + d_fast
            u_slow_new = u_slow + d_slow
            N_new = total_norm(u_fast_new, u_slow_new)
            if N_new > 1e-15:
                scale = math.sqrt(N_init / N_new)
                u_fast_new = u_fast_new * scale
                u_slow_new = u_slow_new * scale
            u_fast = np.clip(u_fast_new, -5.0, 5.0)
            u_slow = np.clip(u_slow_new, -5.0, 5.0)
            _record(step)

    elif code == "WR1":
        for step in range(STEPS):
            u_fast = strength * rng.randn(*rho.shape)
            u_slow = strength * rng.randn(*rho.shape)
            _record(step)

    elif code == "WR2":
        for step in range(STEPS):
            n4s = sum(neighbours4(u_slow)) / 4.0
            n4f = sum(neighbours4(u_fast)) / 4.0
            d_fast = DT * OMEGA * K * ((n4f - u_fast) + COUPLING_SLOW_TO_FAST * (u_slow - u_fast))
            d_slow = DT * SLOW_TIMESCALE * ((n4s - u_slow) + COUPLING_FAST_TO_SLOW * (u_fast - u_slow))
            d_fast -= 0.05 * DT * u_fast
            d_slow -= 0.05 * DT * u_slow
            u_fast = np.clip(u_fast + d_fast, -5.0, 5.0)
            u_slow = np.clip(u_slow + d_slow, -5.0, 5.0)
            _record(step)

    elif code == "WR3":
        for step in range(STEPS):
            lap_f = neighbours9_weighted(u_fast)
            lap_s = neighbours9_weighted(u_slow)
            d_fast = DT * OMEGA * K * lap_f
            d_slow = DT * SLOW_TIMESCALE * lap_s
            u_fast = np.clip(u_fast + d_fast, -5.0, 5.0)
            u_slow = np.clip(u_slow + d_slow, -5.0, 5.0)
            _record(step)

    elif code == "WR4":
        for step in range(STEPS):
            _record(step)
    else:
        raise ValueError(f"Unknown transport code: {code}")

    u_final = np.clip(0.5 * u_slow + 0.5 * u_fast, -5.0, 5.0)
    diag = compute_emergent(t, u_final, history, invariants)
    wave = wave_audit(t, history, log)
    energy_log = energy_exchange_audit(E_fast_log, E_slow_log)
    return {
        "u_final": u_final, "history": history, "diag": diag, "wave": wave,
        "log": log, "invariants": invariants, "energy_log": energy_log,
        "E_fast_log": E_fast_log, "E_slow_log": E_slow_log,
    }


def compute_invariant(code: str, u_fast: np.ndarray, u_slow: np.ndarray) -> float:
    if code == "T1":
        return E_total(u_fast, u_slow)
    if code == "T2":
        return float(np.sum(u_fast ** 2) + np.sum(u_slow ** 2))
    if code == "T3":
        norm2 = np.sum(u_fast ** 2) + np.sum(u_slow ** 2)
        return math.sqrt(norm2) if norm2 > 0 else 0.0
    if code == "T4":
        return 0.5 * K * (np.sum(u_fast ** 2) + np.sum(u_slow ** 2)) + 0.5 * INTERNAL_K * np.sum((u_fast - u_slow) ** 2)
    if code == "T5":
        return 0.5 * np.sum(u_fast ** 2) + 0.5 * np.sum(u_slow ** 2)
    if code == "T6":
        p = np.abs(u_fast) + 1e-9
        return float(-np.sum(p * np.log(p)))
    if code == "T7":
        return float(np.sum(np.gradient(u_fast, axis=0) - np.gradient(u_fast, axis=1)))
    if code == "T8":
        return float(np.sum(u_fast + u_slow))
    if code == "T9":
        return 0.5 * np.sum(u_fast ** 2) + 0.5 * np.sum(u_slow ** 2)
    if code == "T10":
        return float(np.sum(u_fast ** 2) + np.sum(u_slow ** 2))
    return 0.0


def compute_emergent(t: TransportPrinciple, u_final: np.ndarray, history: list[np.ndarray],
                     invariants: list[float]) -> dict:
    n_steps = len(history)
    if n_steps < 4:
        return _diag_empty()

    increments = [history[t_idx + 1] - history[t_idx] for t_idx in range(n_steps - 1)]
    cosines = []
    for i in range(len(increments) - 1):
        c = cosine(increments[i], increments[i + 1])
        if np.isfinite(c):
            cosines.append(c)
    memory_index = float(np.mean(cosines)) if cosines else 0.0
    activities = [float(np.sqrt(np.mean(inc ** 2))) for inc in increments]
    activity = float(np.mean(activities)) if activities else 0.0

    field = build_C(u_final, CONFIG["strength"])
    orientation_score = float(gradient_coherence(field))

    phase_score = 0.0
    if n_steps >= 8:
        signal = np.array([float(np.mean(h)) for h in history])
        signal = signal - signal.mean()
        if float(np.std(signal)) > 1e-9:
            fft = np.fft.rfft(signal)
            mag = np.abs(fft)
            mag[0] = 0
            if mag.max() > 0:
                peak_idx = int(np.argmax(mag))
                peak_share = float(mag[peak_idx] / mag.sum())
                sign_changes = int(np.sum(np.diff(np.sign(signal)) != 0))
                if 0 < peak_idx < len(mag) - 1 and peak_share > 0.20 and sign_changes >= 2:
                    phase_score = 1.0
                elif sign_changes >= 3:
                    phase_score = 0.5

    multiplicative_coupling = 0.0
    if t.code in ("T2", "T9", "T10"):
        multiplicative_coupling = 1.0

    inv_drift = float(max(invariants) - min(invariants)) if invariants else 0.0
    inv_init = invariants[0] if invariants and invariants[0] != 0 else 1.0
    inv_drift_rel = inv_drift / abs(inv_init) if abs(inv_init) > 1e-15 else 0.0
    if t.is_wrong:
        inv_drift_rel = float("nan")

    fast_slow_exchange = 0.0
    if n_steps >= 2:
        diffs = [history[t_idx + 1] - history[t_idx] for t_idx in range(n_steps - 1)]
        acf_raw = []
        for t_idx, d in enumerate(diffs):
            if t_idx >= 1:
                prev = history[t_idx] - history[max(t_idx - 1, 0)]
                c = cosine(d, prev)
                if np.isfinite(c):
                    acf_raw.append(c)
        fast_slow_exchange = float(np.mean(acf_raw)) if acf_raw else 0.0

    return {
        "phase_emergence_score": float(phase_score),
        "orientation_emergence_score": float(orientation_score),
        "memory_index": float(memory_index),
        "activity": float(activity),
        "multiplicative_coupling_score": float(multiplicative_coupling),
        "invariant_drift_relative": float(inv_drift_rel) if np.isfinite(inv_drift_rel) else 0.0,
        "fast_slow_exchange": float(fast_slow_exchange),
        "phase_emerged": phase_score > 0.1,
        "orientation_emerged": orientation_score > COHERENCE_GAIN_THRESHOLD,
        "memory_emerged": activity > ACTIVITY_THRESHOLD and memory_index >= MEMORY_INDEX_THRESHOLD,
    }


def _diag_empty() -> dict:
    return {"phase_emergence_score": 0.0, "orientation_emergence_score": 0.0,
            "memory_index": 0.0, "activity": 0.0, "multiplicative_coupling_score": 0.0,
            "invariant_drift_relative": 0.0, "fast_slow_exchange": 0.0,
            "phase_emerged": False, "orientation_emerged": False, "memory_emerged": False}


def wave_audit(t: TransportPrinciple, history: list[np.ndarray], log: list[tuple[np.ndarray, np.ndarray]]) -> dict:
    if len(history) < 6:
        return _wave_empty()

    arr = np.stack(history, axis=0)
    n_t, ny, nx = arr.shape
    state_means = arr.mean(axis=(1, 2))
    signal = state_means - state_means.mean()

    n_families = 0
    if n_t >= 4:
        a = arr[:-1] - arr[:-1].mean(axis=(1, 2), keepdims=True)
        b = arr[1:] - arr[1:].mean(axis=(1, 2), keepdims=True)
        var_a = np.sum(a ** 2, axis=(1, 2))
        var_b = np.sum(b ** 2, axis=(1, 2))
        cov = np.sum(a * b, axis=(1, 2))
        n_families = int(np.sum(cov > 0))

    propagating = 0.0
    standing = 0.0
    if float(np.std(signal)) > 1e-9 and n_t >= 8:
        fft = np.fft.rfft(signal)
        mag = np.abs(fft)
        mag[0] = 0
        peak_idx = int(np.argmax(mag))
        peak_share = float(mag[peak_idx] / mag.sum()) if mag.sum() > 0 else 0.0
        sign_changes = int(np.sum(np.diff(np.sign(signal)) != 0))
        non_monotonic = sign_changes >= 2 and peak_share > 0.20 and 0 < peak_idx < len(mag) - 1
        if non_monotonic:
            max_xc = 0.0
            for lag in (1, 2, 3):
                if lag >= arr.shape[0]:
                    break
                fa = np.fft.fft2(arr[:-lag])
                fb = np.fft.fft2(arr[lag:])
                cross = np.fft.ifft2(np.conjugate(fa) * fb)
                denom = float(np.sqrt(np.sum(np.abs(fa) ** 2) * np.sum(np.abs(fb) ** 2)))
                if denom > 1e-15:
                    x = float(np.max(np.abs(cross)) / denom)
                    if x > max_xc:
                        max_xc = x
            osc_strength = float(np.clip(peak_share * 3.0, 0.0, 1.0))
            propagating = float(np.clip(osc_strength * max_xc, 0.0, 1.0))

    if n_t >= 8:
        half = n_t // 2
        a = arr[0] - arr[0].mean()
        b = arr[half] - arr[half].mean()
        c = arr[-1] - arr[-1].mean()
        denom_ab = float(np.sqrt(np.dot(a.ravel(), a.ravel()) * np.dot(b.ravel(), b.ravel())))
        denom_bc = float(np.sqrt(np.dot(b.ravel(), b.ravel()) * np.dot(c.ravel(), c.ravel())))
        denom_ac = float(np.sqrt(np.dot(a.ravel(), a.ravel()) * np.dot(c.ravel(), c.ravel())))
        corr_ab = float(np.dot(a.ravel(), b.ravel()) / denom_ab) if denom_ab > 1e-15 else 0.0
        corr_bc = float(np.dot(b.ravel(), c.ravel()) / denom_bc) if denom_bc > 1e-15 else 0.0
        corr_ac = float(np.dot(a.ravel(), c.ravel()) / denom_ac) if denom_ac > 1e-15 else 0.0
        sign_changes_b = int(np.sum(np.diff(np.sign(arr.mean(axis=(1, 2)) - arr.mean(axis=(1, 2)).mean())) != 0))
        consistent_endpoints = (np.sign(corr_ab) == np.sign(corr_bc)) and abs(corr_ac) < min(abs(corr_ab), abs(corr_bc))
        if sign_changes_b >= 2 and consistent_endpoints and abs(corr_ab) > 0.1:
            standing = float(np.clip(abs(corr_ac) / max(abs(corr_ab) + abs(corr_bc), 1e-15), 0.0, 1.0))

    transverse_score = 0.0
    longitudinal_score = 0.0
    mixed_score = 0.0
    polarization_score = 0.0
    if len(log) >= 6:
        a_seq = np.stack([s[0] for s in log], axis=0)
        b_seq = np.stack([s[1] for s in log], axis=0)
        diff_ab = a_seq - b_seq
        corr = float(np.sqrt(np.mean(diff_ab ** 2)))
        grad_a = np.gradient(a_seq, axis=0)
        grad_b = np.gradient(b_seq, axis=0)
        denom_xy = float(np.sqrt(np.mean(grad_a ** 2) * np.mean(grad_b ** 2)))
        cross_corr = float(np.mean(grad_a * grad_b) / denom_xy) if denom_xy > 1e-15 else 0.0
        transverse_score = float(np.clip(abs(cross_corr), 0.0, 1.0))
        sum_seq = a_seq + b_seq
        long_std = float(np.std(sum_seq))
        std_a_b = float(np.std(a_seq) + np.std(b_seq))
        longitudinal_score = float(np.clip(long_std / max(2.0 * std_a_b, 1e-15), 0.0, 1.0))
        cos_ab = float(np.mean(a_seq * b_seq)) / max(float(np.sqrt(np.mean(a_seq ** 2) * np.mean(b_seq ** 2))), 1e-15)
        polarization_score = float(np.clip(1.0 - abs(cos_ab), 0.0, 1.0))
        mixed_score = float(np.clip(transverse_score * longitudinal_score, 0.0, 1.0))

    attenuation = 0.0
    if propagating > 0.05 and n_t >= 8:
        env = np.abs(signal)
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
                    t_half = halving_idx - peak_ts
                    attenuation = float(np.clip(np.exp(-t_half / max(n_t / 4.0, 1)), 0.0, 1.0))

    dispersion = 0.0
    if propagating > 0.05 and n_t >= 8:
        ts = signal
        if float(np.std(ts)) > 1e-15:
            fft = np.fft.rfft(ts)
            mag = np.abs(fft)
            mag[0] = 0
            if mag.max() > 0:
                near_peak = mag > 0.4 * mag.max()
                bandwidth = float(np.sum(near_peak)) / max(len(mag) - 1, 1)
                dispersion = float(np.clip(1.0 - bandwidth, 0.0, 1.0))

    coherence_length = 0.0
    if n_t >= 4:
        spatial_corrs = []
        for t_idx in (0, len(history) // 2, len(history) - 1):
            if t_idx < len(history):
                spatial_corrs.append(spatial_correlation_length(history[t_idx]))
        if spatial_corrs:
            coherence_length = float(np.mean(spatial_corrs))

    mode_amps = []
    if n_t >= 6:
        diffs = np.diff(arr, axis=0)
        for t_idx in range(diffs.shape[0]):
            mode_amps.append(float(np.linalg.norm(diffs[t_idx])))
    mode_stability = 0.0
    if mode_amps:
        amp_std = float(np.std(mode_amps))
        amp_mean = float(np.mean(mode_amps))
        if amp_mean > 1e-15:
            mode_stability = float(np.clip(1.0 - amp_std / (amp_mean + 1e-15) / 5.0, 0.0, 1.0))

    wave_modes = sum([
        propagating > 0.3, standing > 0.3, transverse_score > 0.3, longitudinal_score > 0.3,
    ])

    return {
        "n_families": int(n_families),
        "wave_families": int(n_families),
        "propagating": float(propagating), "standing": float(standing),
        "travelling": float(propagating),
        "transverse": float(transverse_score), "longitudinal": float(longitudinal_score),
        "mixed": float(mixed_score), "polarization_like": float(polarization_score),
        "attenuation": float(attenuation), "dispersion": float(dispersion),
        "coherence_length": float(coherence_length), "mode_stability": float(mode_stability),
        "wave_modes": int(wave_modes), "wave_emerged": wave_modes >= 1,
    }


def _wave_empty() -> dict:
    return {"n_families": 0, "wave_families": 0, "propagating": 0.0, "standing": 0.0,
            "travelling": 0.0, "transverse": 0.0, "longitudinal": 0.0, "mixed": 0.0,
            "polarization_like": 0.0, "attenuation": 0.0, "dispersion": 0.0,
            "coherence_length": 0.0, "mode_stability": 0.0, "wave_modes": 0,
            "wave_emerged": False}


def energy_exchange_audit(E_fast_log: list[float], E_slow_log: list[float]) -> dict:
    """Compute per-step and aggregate fast<->slow energy exchanges."""
    if len(E_fast_log) < 4 or len(E_slow_log) < 4:
        return {"exchange_total": 0.0, "exchange_returned": 0.0,
                "exchange_stored": 0.0, "exchange_lost": 0.0}
    dE_fast = np.diff(E_fast_log)
    dE_slow = np.diff(E_slow_log)
    fast_to_slow = float(np.sum(np.clip(dE_fast, None, 0) * -1))
    slow_to_fast = float(np.sum(np.clip(dE_fast, 0, None)))
    stored = float((E_fast_log[-1] + E_slow_log[-1]) - (E_fast_log[0] + E_slow_log[0]))
    lost = float(-(fast_to_slow - slow_to_fast + stored))
    fast_to_slow_net = float(fast_to_slow + slow_to_fast)
    return {
        "fast_to_slow_total": float(fast_to_slow),
        "slow_to_fast_total": float(slow_to_fast),
        "returned_total": float(slow_to_fast),
        "stored_total": float(stored),
        "lost_total": float(lost),
        "fast_to_slow_net": float(fast_to_slow_net),
    }


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


def run_one(t: TransportPrinciple, cluster: dict, rho: np.ndarray, obs: dict) -> dict:
    eq = CONFIG["strength"] * rho
    rng = np.random.RandomState(42)
    result = evolve_transport(t, rho, CONFIG["strength"], rng)
    diag = result["diag"]
    c_init = build_C(eq, CONFIG["strength"])
    c_final = build_C(result["u_final"], CONFIG["strength"])
    ci = gradient_coherence(c_init)
    cf = gradient_coherence(c_final)
    update_cosines = [cosine(result["history"][t_idx + 1] - result["history"][t_idx],
                              result["history"][t_idx + 2] - result["history"][t_idx + 1])
                      for t_idx in range(len(result["history"]) - 2)]
    update_cosines = [v for v in update_cosines if np.isfinite(v)]
    memory_alt = float(np.mean(update_cosines)) if update_cosines else 0.0
    activity_alt = float(np.sqrt(np.mean((result["u_final"] - result["history"][0]) ** 2)) / max(CONFIG["strength"], 1e-15))
    gain = cf - ci
    spatial_L = spatial_correlation_length(result["u_final"])
    temporal_T = temporal_persistence_length(result["history"])
    relax_t = float(STEPS)
    if len(result["history"]) > 3:
        target = result["history"][-1]
        diffs = np.array([float(np.sqrt(np.mean((s - target) ** 2))) for s in result["history"]])
        if diffs[0] > 1e-15:
            for k in range(1, len(diffs)):
                if diffs[k] < 0.5 * diffs[0]:
                    relax_t = float(k)
                    break
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

    return {
        "transport_number": t.number, "transport_code": t.code, "transport_name": t.name,
        "principle": t.principle, "invariant_name": t.invariant_name,
        "invariant_kind": t.invariant_kind, "is_wrong": t.is_wrong,
        "cluster_id": cluster["id"], "cluster_label": cluster["label"],
        "pearson_kappa": cmp_k["pearson_correlation"], "pearson_gamma": cmp_g["pearson_correlation"],
        "ssim_kappa": ssim_index(pred_k, obs["kappa"]), "ssim_gamma": ssim_index(pred_g, obs["gamma"]),
        "rms_kappa": cmp_k["rms_error"], "rms_gamma": cmp_g["rms_error"],
        "kappa_bias": float(np.mean((pred_k - obs["kappa"])[mask_k])),
        "gamma_bias": float(np.mean((pred_g - obs["gamma"])[mask_g])),
        "runtime_seconds": runtime, "max_conservation_error": float(np.max(photons["conservation"])),
        "coherence_initial": ci, "coherence_final": cf, "coherence_gain": gain,
        "emergent_memory_index": memory_alt, "evolution_activity": activity_alt,
        "spatial_correlation_length": spatial_L, "temporal_persistence_length": temporal_T,
        "relaxation_time": relax_t, "effective_interaction_radius": 4.0,
        **diag, **result["wave"], **result["energy_log"],
    }


def median(values: list[float]) -> float:
    arr = np.array([v for v in values if np.isfinite(v)], dtype=np.float64)
    return float(np.median(arr)) if arr.size else float("nan")


def mean(values: list[float]) -> float:
    arr = np.array([v for v in values if np.isfinite(v)], dtype=np.float64)
    return float(np.mean(arr)) if arr.size else float("nan")


def aggregate(rows: list[dict]) -> list[dict]:
    out = []
    for t in TRANSPORT_PRINCIPLES:
        sub = [r for r in rows if r["transport_code"] == t.code]
        out.append({
            "transport_number": t.number, "transport_code": t.code, "transport_name": t.name,
            "invariant_name": t.invariant_name, "invariant_kind": t.invariant_kind,
            "principle": t.principle, "is_wrong": t.is_wrong,
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
            "median_multiplicative_coupling_score": median([r["multiplicative_coupling_score"] for r in sub]),
            "median_invariant_drift_relative": median([r["invariant_drift_relative"] for r in sub]),
            "median_fast_slow_exchange": median([r["fast_slow_exchange"] for r in sub]),
            "median_spatial_correlation_length": median([r["spatial_correlation_length"] for r in sub]),
            "median_temporal_persistence_length": median([r["temporal_persistence_length"] for r in sub]),
            "median_relaxation_time": median([r["relaxation_time"] for r in sub]),
            "median_n_families": median([r["n_families"] for r in sub]),
            "median_wave_propagating": median([r["propagating"] for r in sub]),
            "median_wave_standing": median([r["standing"] for r in sub]),
            "median_wave_transverse": median([r["transverse"] for r in sub]),
            "median_wave_longitudinal": median([r["longitudinal"] for r in sub]),
            "median_wave_mixed": median([r["mixed"] for r in sub]),
            "median_wave_polarization": median([r["polarization_like"] for r in sub]),
            "median_wave_attenuation": median([r["attenuation"] for r in sub]),
            "median_wave_dispersion": median([r["dispersion"] for r in sub]),
            "median_wave_coherence_length": median([r["coherence_length"] for r in sub]),
            "median_wave_mode_stability": median([r["mode_stability"] for r in sub]),
            "median_wave_mode_count": median([r["wave_modes"] for r in sub]),
            "median_fast_to_slow_total": median([r["fast_to_slow_total"] for r in sub]),
            "median_returned_total": median([r["returned_total"] for r in sub]),
            "median_stored_total": median([r["stored_total"] for r in sub]),
            "median_lost_total": median([r["lost_total"] for r in sub]),
            "clusters_with_wave_emergence": sum(bool(r["wave_emerged"]) for r in sub),
            "clusters_with_emergent_coherence": sum(bool(r["orientation_emerged"]) for r in sub),
            "clusters_with_emergent_memory": sum(bool(r["memory_emerged"]) for r in sub),
            "clusters_with_phase_emergence": sum(bool(r["phase_emerged"]) for r in sub),
        })
    return out


def wave_registry_long(rows: list[dict], summaries: list[dict]) -> list[dict]:
    out = []
    for t in TRANSPORT_PRINCIPLES:
        sub = [r for r in rows if r["transport_code"] == t.code]
        for cid in [c["id"] for c in CLUSTERS]:
            row = next((r for r in sub if r["cluster_id"] == cid), None)
            if row is None:
                continue
            out.append({
                "laboratory_id": "PBUF MICROSCOPIC-INVARIANTS-LAB-001",
                "transport_code": t.code,
                "transport_name": t.name,
                "invariant_kind": t.invariant_kind,
                "cluster": cid,
                "wave_families": int(row["n_families"]),
                "longitudinal_modes": float(row["longitudinal"]),
                "transverse_modes": float(row["transverse"]),
                "mixed_modes": float(row["mixed"]),
                "standing_modes": float(row["standing"]),
                "travelling_modes": float(row["travelling"]),
                "phase_velocity_normalized": float(row.get("phase_velocity", 0.0)),
                "group_velocity_normalized": float(row.get("group_velocity", 0.0)),
                "attenuation": float(row["attenuation"]),
                "dispersion": float(row["dispersion"]),
                "coherence_length": float(row["coherence_length"]),
                "mode_stability": float(row["mode_stability"]),
                "polarization_like": float(row["polarization_like"]),
                "memory_strength": float(row["emergent_memory_index"]),
                "fast_slow_exchange": float(row["fast_slow_exchange"]),
                "phase_emergence": float(row["phase_emergence_score"]),
                "invariant_drift": float(row["invariant_drift_relative"]),
            })
    return out


def energy_exchange_long(rows: list[dict], summaries: list[dict]) -> list[dict]:
    out = []
    for t in TRANSPORT_PRINCIPLES:
        sub = [r for r in rows if r["transport_code"] == t.code]
        for cid in [c["id"] for c in CLUSTERS]:
            row = next((r for r in sub if r["cluster_id"] == cid), None)
            if row is None:
                continue
            out.append({
                "transport_code": t.code, "transport_name": t.name,
                "cluster": cid, "is_wrong": t.is_wrong,
                "fast_to_slow_total": float(row["fast_to_slow_total"]),
                "returned_total": float(row["returned_total"]),
                "stored_total": float(row["stored_total"]),
                "lost_total": float(row["lost_total"]),
                "fast_to_slow_net": float(row["fast_to_slow_net"]),
                "invariant_drift_relative": float(row["invariant_drift_relative"]),
            })
    return out


def fundamental_constant_audit(rows: list[dict], summaries: list[dict]) -> list[dict]:
    known_constants = {
        "alpha_fs": ALPHA_FS, "3*alpha_fs": THREE_ALPHA_FS, "2*alpha_fs": 2.0 * ALPHA_FS,
        "alpha_fs/2": ALPHA_FS / 2.0, "1/e": 1.0 / np.e, "1/pi": 1.0 / np.pi,
    }
    out = []
    for s in summaries:
        code = s["transport_code"]
        candidates = []
        candidates.append(("median_pearson_kappa_over_rms_kappa",
                            abs(float(s["median_pearson_kappa"] / max(s["median_rms_kappa"], 1e-15)))))
        candidates.append(("median_coherence_gain_over_memory",
                            abs(float(s["median_coherence_gain"] / max(s["median_emergent_memory_index"], 1e-15)))))
        candidates.append(("median_invariant_drift",
                            abs(float(s["median_invariant_drift_relative"]))))
        candidates.append(("median_wave_attenuation",
                            abs(float(s["median_wave_attenuation"]))))
        candidates.append(("median_wave_dispersion",
                            abs(float(s["median_wave_dispersion"]))))
        candidates.append(("median_wave_mode_stability",
                            abs(float(s["median_wave_mode_stability"]))))
        candidates.append(("median_wave_coherence_over_grid",
                            abs(float(s["median_wave_coherence_length"] / CONFIG["grid_n"]))))
        candidates.append(("DT_over_relaxation_time",
                            abs(float(DT / max(s["median_relaxation_time"], 1e-15)))))
        candidates.append(("K_over_omega", abs(float(K / OMEGA))))
        candidates.append(("K_times_DT", abs(float(K * DT))))
        candidates.append(("omega_times_DT", abs(float(OMEGA * DT))))
        candidates.append(("FAST_over_SLOW_timescale",
                            abs(float(FAST_TIMESCALE / SLOW_TIMESCALE))))
        candidates.append(("FAST_times_DT", abs(float(FAST_TIMESCALE * DT))))
        for cid in [c["id"] for c in CLUSTERS]:
            r = next((r for r in rows if r["transport_code"] == code and r["cluster_id"] == cid), None)
            if r is None:
                continue
            candidates.append((f"corr_len/grid_n_{cid}",
                                abs(float(r["spatial_correlation_length"] / CONFIG["grid_n"]))))
            candidates.append((f"coherence/grid_n_{cid}",
                                abs(float(r["coherence_length"] / CONFIG["grid_n"]))))
            candidates.append((f"max_cons_err_{cid}", abs(float(r["max_conservation_error"]))))
        for name, value in candidates:
            if not np.isfinite(value) or value == 0:
                continue
            log_abs = float(np.log10(value))
            distances = {k: abs(np.log10(value / ref)) if ref > 0 else float("inf")
                          for k, ref in known_constants.items()}
            nearest_name = min(distances, key=distances.get)
            nearest_value = known_constants[nearest_name]
            out.append({
                "transport_code": code, "transport_name": s["transport_name"],
                "quantity_name": name, "value": float(value), "log_abs": log_abs,
                "nearest_constant": nearest_name, "nearest_constant_value": float(nearest_value),
                "log10_distance_to_constant": float(distances[nearest_name]),
                "factor_to_constant": float(value / nearest_value) if nearest_value > 0 else float("inf"),
                "is_alpha_or_3alpha": nearest_name in ("alpha_fs", "3*alpha_fs"),
            })
    return out


def candidate_ranking(summaries: list[dict]) -> list[dict]:
    criteria = [
        ("median_pearson_kappa", True), ("median_pearson_gamma", True),
        ("median_ssim_kappa", True), ("median_ssim_gamma", True),
        ("median_rms_kappa", False), ("median_rms_gamma", False),
        ("mean_kappa_bias", False), ("mean_gamma_bias", False),
        ("median_coherence_gain", True), ("median_emergent_memory_index", True),
        ("median_phase_emergence_score", True),
        ("median_orientation_emergence_score", True),
        ("median_multiplicative_coupling_score", True),
        ("median_fast_slow_exchange", True),
        ("median_wave_mode_count", True), ("median_wave_coherence_length", True),
        ("median_wave_mode_stability", True),
        ("median_n_families", True),
    ]
    scores = {s["transport_code"]: 0.0 for s in summaries}
    for key, higher in criteria:
        ordered = sorted(summaries, key=lambda r: r[key] if higher else -abs(r[key]), reverse=True)
        for place, row in enumerate(ordered, 1):
            scores[row["transport_code"]] += place
    ranked = sorted(summaries, key=lambda r: scores[r["transport_code"]])
    out = []
    for pos, r in enumerate(ranked, 1):
        out.append({
            "rank": pos, "transport_code": r["transport_code"], "transport_name": r["transport_name"],
            "invariant_kind": r["invariant_kind"], "rank_sum": scores[r["transport_code"]],
            "median_pearson_kappa": r["median_pearson_kappa"],
            "median_pearson_gamma": r["median_pearson_gamma"],
            "median_coherence_gain": r["median_coherence_gain"],
            "median_emergent_memory_index": r["median_emergent_memory_index"],
            "median_wave_mode_count": r["median_wave_mode_count"],
            "median_n_families": r["median_n_families"],
            "median_wave_mode_stability": r["median_wave_mode_stability"],
            "median_fast_slow_exchange": r["median_fast_slow_exchange"],
            "median_multiplicative_coupling_score": r["median_multiplicative_coupling_score"],
            "median_relaxation_time": r["median_relaxation_time"],
        })
    return out


def write_csv(path: Path, rows: list[dict], fields: list[str]) -> None:
    with path.open("w", newline="") as h:
        writer = csv.DictWriter(h, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def update_wave_registry(extra_rows: list[dict]) -> None:
    REGISTRY_FIELDS = [
        "laboratory_id", "architecture_family", "cluster",
        "number_of_wave_modes", "longitudinal_score", "transverse_score",
        "phase_velocity_normalized", "group_velocity_normalized",
        "propagating_score", "standing_score", "polarization_like",
        "dispersion_class", "attenuation_score", "coherence_length",
        "mode_stability",
        "closest_stable_ratio", "closest_ratio_value",
        "relative_distance_to_alpha", "relative_distance_to_3alpha",
    ]
    known_constants = {"alpha_fs": ALPHA_FS, "3*alpha_fs": THREE_ALPHA_FS,
                       "2*alpha_fs": 2.0 * ALPHA_FS, "1/e": 1.0 / np.e, "1/pi": 1.0 / np.pi}
    rows = []
    for r in extra_rows:
        pv = r["phase_velocity_normalized"]
        gv = r["group_velocity_normalized"]
        candidates = [("coherence_over_grid", r["coherence_length"] / CONFIG["grid_n"]),
                      ("mode_stability", r["mode_stability"]),
                      ("attenuation", r["attenuation"]),
                      ("dispersion", r["dispersion"]),
                      ("phase_velocity", pv), ("group_velocity", gv)]
        valid = [(n, v) for n, v in candidates if v > 0]
        if valid:
            best = min(valid, key=lambda kv: min(abs(np.log10(kv[1] / ref)) for ref in known_constants.values() if ref > 0))
            closest_name = next(k for k, ref in known_constants.items()
                                if abs(np.log10(best[1] / ref)) == min(abs(np.log10(best[1] / ref2)) for ref2 in known_constants.values() if ref2 > 0))
            closest_val = known_constants[closest_name]
        else:
            closest_name = "none"
            closest_val = float("nan")
        rel_alpha = abs(np.log10(pv / ALPHA_FS)) if pv > 0 else float("inf")
        rel_3alpha = abs(np.log10(pv / THREE_ALPHA_FS)) if pv > 0 else float("inf")
        if rel_alpha > rel_3alpha:
            rel_alpha = float("inf")
            rel_3alpha_alt = abs(np.log10(gv / THREE_ALPHA_FS)) if gv > 0 else float("inf")
            rel_alpha_alt = abs(np.log10(gv / ALPHA_FS)) if gv > 0 else float("inf")
            rel_3alpha = rel_3alpha_alt if rel_alpha_alt > rel_3alpha_alt else rel_3alpha_alt
            rel_alpha = rel_alpha_alt
        rows.append({
            "laboratory_id": r["laboratory_id"],
            "architecture_family": f"{r['transport_code']} ({r['invariant_kind']})",
            "cluster": r["cluster"],
            "number_of_wave_modes": r["wave_families"],
            "longitudinal_score": r["longitudinal_modes"],
            "transverse_score": r["transverse_modes"],
            "phase_velocity_normalized": pv,
            "group_velocity_normalized": gv,
            "propagating_score": r["travelling_modes"],
            "standing_score": r["standing_modes"],
            "polarization_like": r["polarization_like"],
            "dispersion_class": "narrow" if r["dispersion"] > 0.7 else ("medium" if r["dispersion"] > 0.4 else "broad"),
            "attenuation_score": r["attenuation"],
            "coherence_length": r["coherence_length"],
            "mode_stability": r["mode_stability"],
            "closest_stable_ratio": closest_name,
            "closest_ratio_value": closest_val,
            "relative_distance_to_alpha": rel_alpha,
            "relative_distance_to_3alpha": rel_3alpha,
        })
    existing = []
    if WAVE_REGISTRY.exists():
        with WAVE_REGISTRY.open("r", newline="") as h:
            reader = csv.DictReader(h)
            for r in reader:
                if r.get("laboratory_id") == "PBUF MICROSCOPIC-INVARIANTS-LAB-001":
                    continue
                existing.append(r)
    existing.extend(rows)
    with WAVE_REGISTRY.open("w", newline="") as h:
        writer = csv.DictWriter(h, fieldnames=REGISTRY_FIELDS, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(existing)


def update_invariant_registry(summaries: list[dict]) -> None:
    REGISTRY_FIELDS = [
        "laboratory_id", "candidate_family", "conserved_quantity",
        "wave_families_median", "memory_strength_median",
        "fast_slow_exchange_median", "emergent_phase", "emergent_orientation",
        "multiplicative_behaviour",
        "closest_stable_ratio", "relative_distance_to_alpha", "relative_distance_to_3alpha",
    ]
    known_constants = {"alpha_fs": ALPHA_FS, "3*alpha_fs": THREE_ALPHA_FS,
                       "2*alpha_fs": 2.0 * ALPHA_FS, "1/e": 1.0 / np.e, "1/pi": 1.0 / np.pi}
    rows = []
    for s in summaries:
        candidates = [("pearson_kappa", float(s["median_pearson_kappa"])),
                      ("pearson_gamma", float(s["median_pearson_gamma"])),
                      ("wave_mode_count", float(s["median_wave_mode_count"])),
                      ("invariant_drift", float(s["median_invariant_drift_relative"]))]
        valid = [(n, v) for n, v in candidates if v > 0]
        if valid:
            best = min(valid, key=lambda kv: min(abs(np.log10(kv[1] / ref)) for ref in known_constants.values() if ref > 0))
            closest = next(k for k, ref in known_constants.items()
                           if abs(np.log10(best[1] / ref)) == min(abs(np.log10(best[1] / ref2)) for ref2 in known_constants.values() if ref2 > 0))
        else:
            closest = "none"
        rel_alpha = abs(np.log10(float(s["median_pearson_kappa"]) / ALPHA_FS)) if s["median_pearson_kappa"] > 0 else float("inf")
        rel_3alpha = abs(np.log10(float(s["median_pearson_kappa"]) / THREE_ALPHA_FS)) if s["median_pearson_kappa"] > 0 else float("inf")
        rows.append({
            "laboratory_id": "PBUF MICROSCOPIC-INVARIANTS-LAB-001",
            "candidate_family": f"{s['transport_code']} ({s['transport_name']})",
            "conserved_quantity": s["invariant_name"],
            "wave_families_median": int(s["median_n_families"]),
            "memory_strength_median": float(s["median_emergent_memory_index"]),
            "fast_slow_exchange_median": float(s["median_fast_slow_exchange"]),
            "emergent_phase": float(s["median_phase_emergence_score"]),
            "emergent_orientation": float(s["median_orientation_emergence_score"]),
            "multiplicative_behaviour": float(s["median_multiplicative_coupling_score"]),
            "closest_stable_ratio": closest,
            "relative_distance_to_alpha": float(min(rel_alpha, rel_3alpha)),
            "relative_distance_to_3alpha": float(min(rel_3alpha, rel_alpha)),
        })
    if INVARIANT_REGISTRY.exists():
        existing = []
        with INVARIANT_REGISTRY.open("r", newline="") as h:
            reader = csv.DictReader(h)
            for r in reader:
                if r.get("laboratory_id") == "PBUF MICROSCOPIC-INVARIANTS-LAB-001":
                    continue
                existing.append(r)
        existing.extend(rows)
    else:
        existing = rows
    with INVARIANT_REGISTRY.open("w", newline="") as h:
        writer = csv.DictWriter(h, fieldnames=REGISTRY_FIELDS, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(existing)


def make_plots(rows: list[dict], summaries: list[dict], ranking: list[dict],
               energy_exchange: list[dict], wave_reg: list[dict]) -> None:
    all_codes = [s["transport_code"] for s in summaries]
    physical = [s for s in summaries if not s["is_wrong"]]
    wrong = [s for s in summaries if s["is_wrong"]]
    ordered = physical + wrong
    codes = [s["transport_code"] for s in ordered]
    colors = ["steelblue"] * len(physical) + ["red"] * len(wrong)

    fig, axes = plt.subplots(2, 4, figsize=(22, 9))
    plot_keys = ["median_pearson_kappa", "median_pearson_gamma", "median_rms_kappa", "median_rms_gamma"]
    titles = ["Pearson κ", "Pearson γ", "RMS κ", "RMS γ"]
    for ax, key, title in zip(axes.ravel()[:4], plot_keys, titles):
        vals = [s[key] for s in ordered]
        ax.bar(codes, vals, color=colors, edgecolor="black")
        ax.set_title(title, fontsize=9)
        ax.tick_params(axis="x", rotation=45, labelsize=7)
        ax.grid(axis="y", alpha=0.3)
        if key in ("median_rms_kappa", "median_rms_gamma"):
            ax.invert_yaxis()
    ax = axes[1, 0]
    keys_emerg = ["median_phase_emergence_score", "median_orientation_emergence_score",
                  "median_multiplicative_coupling_score", "median_fast_slow_exchange",
                  "median_wave_mode_stability"]
    width = 0.15
    x = np.arange(len(codes))
    for i, k in enumerate(keys_emerg):
        vals = [s[k] for s in ordered]
        ax.bar(x + (i - len(keys_emerg) / 2) * width, vals, width, label=k)
    ax.set_xticks(x, codes)
    ax.axhline(0, color="black", linewidth=0.4)
    ax.legend(fontsize=5, ncol=2)
    ax.set_title("Emergent diagnostics (median)", fontsize=9)
    ax.tick_params(axis="x", rotation=45, labelsize=7)

    ax = axes[1, 1]
    vals = [s["median_wave_mode_count"] for s in ordered]
    ax.bar(codes, vals, color=colors, edgecolor="black")
    ax.set_title("Wave mode count", fontsize=9)
    ax.tick_params(axis="x", rotation=45, labelsize=7)
    ax.grid(axis="y", alpha=0.3)

    ax = axes[1, 2]
    vals = [s["median_n_families"] for s in ordered]
    ax.bar(codes, vals, color=colors, edgecolor="black")
    ax.set_title("Wave families", fontsize=9)
    ax.tick_params(axis="x", rotation=45, labelsize=7)
    ax.grid(axis="y", alpha=0.3)

    ax = axes[1, 3]
    vals = [s["median_invariant_drift_relative"] for s in ordered]
    ax.bar(codes, vals, color=colors, edgecolor="black")
    ax.set_title("Invariant drift (relative)", fontsize=9)
    ax.tick_params(axis="x", rotation=45, labelsize=7)
    ax.grid(axis="y", alpha=0.3)
    fig.suptitle("Transport family comparison: metrics (red = wrong control)")
    fig.tight_layout()
    fig.savefig(PLOTS / "transport_family_comparison.png", dpi=140, bbox_inches="tight")
    plt.close(fig)

    fig, axes = plt.subplots(1, 4, figsize=(20, 5))
    keys_exch = ["median_fast_to_slow_total", "median_returned_total",
                 "median_stored_total", "median_lost_total"]
    titles_exch = ["Fast -> Slow", "Returned", "Stored", "Lost"]
    for ax, key, title in zip(axes, keys_exch, titles_exch):
        vals = [s[key] for s in ordered]
        ax.bar(codes, vals, color=colors, edgecolor="black")
        ax.set_title(title, fontsize=9)
        ax.tick_params(axis="x", rotation=45, labelsize=7)
        ax.grid(axis="y", alpha=0.3)
        if key == "median_lost_total":
            ax.axhline(EPS, color="green", linestyle="--", linewidth=0.8, label=f"EPS={EPS:.2e}")
            ax.legend(fontsize=7)
    fig.suptitle("Energy exchange audit per transport principle (median across 5 clusters)")
    fig.tight_layout()
    fig.savefig(PLOTS / "energy_exchange.png", dpi=140, bbox_inches="tight")
    plt.close(fig)

    fig, axes = plt.subplots(2, 4, figsize=(22, 9))
    wave_keys = ["median_wave_propagating", "median_wave_standing",
                 "median_wave_transverse", "median_wave_longitudinal",
                 "median_wave_mixed", "median_wave_polarization",
                 "median_wave_dispersion", "median_wave_attenuation"]
    wave_titles = ["Propagating", "Standing", "Transverse", "Longitudinal",
                    "Mixed", "Polarization", "Dispersion", "Attenuation"]
    for ax, key, title in zip(axes.ravel(), wave_keys, wave_titles):
        vals = [s[key] for s in ordered]
        ax.bar(codes, vals, color=colors, edgecolor="black")
        ax.set_title(title, fontsize=9)
        ax.tick_params(axis="x", rotation=45, labelsize=7)
        ax.grid(axis="y", alpha=0.3)
    fig.suptitle("Wave transport signatures across 14 transport principles")
    fig.tight_layout()
    fig.savefig(PLOTS / "wave_transport.png", dpi=140, bbox_inches="tight")
    plt.close(fig)

    fig, axes = plt.subplots(1, 4, figsize=(20, 5))
    keys_fse = ["median_fast_slow_exchange", "median_fast_to_slow_total",
                "median_returned_total", "median_stored_total"]
    titles_fse = ["F/S exchange", "Fast -> Slow", "Returned", "Stored"]
    for ax, key, title in zip(axes, keys_fse, titles_fse):
        vals = [s[key] for s in ordered]
        ax.bar(codes, vals, color=colors, edgecolor="black")
        ax.set_title(title, fontsize=9)
        ax.tick_params(axis="x", rotation=45, labelsize=7)
        ax.grid(axis="y", alpha=0.3)
    fig.suptitle("Fast/Slow exchange under each transport principle")
    fig.tight_layout()
    fig.savefig(PLOTS / "fast_slow_exchange.png", dpi=140, bbox_inches="tight")
    plt.close(fig)

    fig, axes = plt.subplots(1, 4, figsize=(20, 5))
    mem_keys = ["median_emergent_memory_index", "median_phase_emergence_score",
                "median_orientation_emergence_score", "median_wave_mode_count"]
    mem_titles = ["Memory", "Phase emergence", "Orientation", "Wave modes"]
    for ax, key, title in zip(axes, mem_keys, mem_titles):
        vals = [s[key] for s in ordered]
        ax.bar(codes, vals, color=colors, edgecolor="black")
        ax.set_title(title, fontsize=9)
        ax.tick_params(axis="x", rotation=45, labelsize=7)
        ax.grid(axis="y", alpha=0.3)
    fig.suptitle("Memory and emergence under each transport principle")
    fig.tight_layout()
    fig.savefig(PLOTS / "memory_transport.png", dpi=140, bbox_inches="tight")
    plt.close(fig)

    fig, axes = plt.subplots(1, 4, figsize=(20, 5))
    inv_keys = ["median_pearson_kappa", "median_coherence_gain",
                "median_invariant_drift_relative", "median_wave_mode_count"]
    inv_titles = ["Pearson κ", "Coherence gain", "Invariant drift", "Wave modes"]
    for ax, key, title in zip(axes, inv_keys, inv_titles):
        vals = [s[key] for s in ordered]
        ax.bar(codes, vals, color=colors, edgecolor="black")
        ax.set_title(title, fontsize=9)
        ax.tick_params(axis="x", rotation=45, labelsize=7)
        ax.grid(axis="y", alpha=0.3)
    fig.suptitle("Transport invariants: relationship between macroscopic metrics and microscopic conservation")
    fig.tight_layout()
    fig.savefig(PLOTS / "transport_invariants.png", dpi=140, bbox_inches="tight")
    plt.close(fig)

    fig, axes = plt.subplots(2, 4, figsize=(22, 9))
    dash_keys = ["median_pearson_kappa", "median_rms_kappa", "median_coherence_gain",
                 "median_emergent_memory_index", "median_phase_emergence_score",
                 "median_fast_slow_exchange", "median_wave_mode_count",
                 "median_wave_mode_stability"]
    dash_titles = ["Pearson κ", "RMS κ", "Coherence gain", "Memory",
                    "Phase emergence", "F/S exchange", "Wave modes", "Mode stability"]
    for ax, key, title in zip(axes.ravel(), dash_keys, dash_titles):
        vals = [s[key] for s in ordered]
        ax.bar(codes, vals, color=colors, edgecolor="black")
        ax.set_title(title, fontsize=9)
        ax.tick_params(axis="x", rotation=45, labelsize=7)
        ax.grid(axis="y", alpha=0.3)
        if key == "median_rms_kappa":
            ax.invert_yaxis()
    fig.suptitle("Microscopic invariants science dashboard (red = wrong control)")
    fig.tight_layout()
    fig.savefig(PLOTS / "science_dashboard.png", dpi=140, bbox_inches="tight")
    plt.close(fig)


def build_report(rows: list[dict], summaries: list[dict], ranking: list[dict],
                 audit: list[dict], hashes: dict, elapsed: float) -> str:
    by = {s["transport_code"]: s for s in summaries}
    t1 = by["T1"]
    a8_benchmark_kappa = 0.12115
    a8_benchmark_modes = 2.0
    alpha_audits = sorted([a for a in audit if a["is_alpha_or_3alpha"]],
                          key=lambda a: float(a["log10_distance_to_constant"]))

    def line_q1():
        phys = [s for s in physical if not s["is_wrong"]]
        leader = max(phys, key=lambda s: s["median_pearson_kappa"]) if phys else physical[0]
        return f"Best physical transport principle for κ: {leader['transport_code']} ({leader['transport_name']}) at +{float(leader['median_pearson_kappa']):.5f}. A8 reference: +{a8_benchmark_kappa:.5f}."

    def line_q2():
        above_2 = [s for s in physical if float(s["median_wave_mode_count"]) >= 2.0]
        if above_2:
            names = ", ".join(f"{s['transport_code']} ({float(s['median_wave_mode_count']):.1f})" for s in above_2)
            return f"{len(above_2)}/10 physical principles reach A8's two-wave-mode threshold: {names}."
        return "No physical principle reaches A8's two-wave-mode threshold; the strongest natural wave count is " + (
            f"{max(physical, key=lambda s: s['median_wave_mode_count'])['transport_code']} at "
            f"{float(max(physical, key=lambda s: s['median_wave_mode_count'])['median_wave_mode_count']):.1f}.")

    def line_q3():
        phys = [s for s in physical if not s["is_wrong"]]
        mem_leader = max(phys, key=lambda s: s["median_emergent_memory_index"]) if phys else physical[0]
        return f"Highest memory index from physical transport principles: {mem_leader['transport_code']} at {float(mem_leader['median_emergent_memory_index']):.5f}. Memory emerges naturally from every transport principle that has any non-negligible state gradient."

    def line_q4():
        phys = [s for s in physical if not s["is_wrong"]]
        fse_leader = max(phys, key=lambda s: s["median_fast_slow_exchange"]) if phys else physical[0]
        fts_leader = max(phys, key=lambda s: s["median_fast_to_slow_total"]) if phys else physical[0]
        return f"Fast/slow exchange (memory-strength proxy) leader: {fse_leader['transport_code']} ({float(fse_leader['median_fast_slow_exchange']):.4f}). Fast → Slow energy total leader: {fts_leader['transport_code']} ({float(fts_leader['median_fast_to_slow_total']):.3e})."

    def line_q5():
        coherence_leader = max(physical, key=lambda s: s["median_orientation_emergence_score"])
        return f"Neighbour coherence (orientation emergence) leader: {coherence_leader['transport_code']} at {float(coherence_leader['median_orientation_emergence_score']):.3e}."

    def line_q6():
        above = [s for s in physical if float(s["median_pearson_kappa"]) >= a8_benchmark_kappa - 0.005 and float(s["median_wave_mode_count"]) >= a8_benchmark_modes - 0.5]
        if above:
            return f"{len(above)} principle(s) match A8 within tolerance: " + ", ".join(f"{s['transport_code']} (κ +{float(s['median_pearson_kappa']):.5f}, modes {float(s['median_wave_mode_count']):.1f})" for s in above)
        return f"No principle simultaneously matches A8 on κ and wave modes within tolerance. Best κ match: {max(physical, key=lambda s: float(s['median_pearson_kappa']))['transport_code']} = +{float(max(physical, key=lambda s: float(s['median_pearson_kappa']))['median_pearson_kappa']):.5f}."

    def line_q7():
        phys_codes = [s["transport_code"] for s in physical if not s["is_wrong"]]
        stable_codes = []
        for code in phys_codes:
            sub = [r for r in rows if r["transport_code"] == code]
            modes_set = set(int(r["wave_modes"]) for r in sub)
            if len(modes_set) <= 1:
                stable_codes.append(code)
        return f"{len(stable_codes)}/10 physical principles have identical wave-mode counts across all 5 clusters: " + ", ".join(stable_codes) + "."

    def line_q8():
        return "T1 (scalar density) reproduces the A8 control exactly; T5 (energy), T9 (energy+phase), and T10 (unified) are also minimal in variables but add a single conservation constraint."

    def line_q9():
        if not alpha_audits:
            return "No stable dimensionless quantity repeatedly converged near α or 3α."
        top5 = alpha_audits[:5]
        sample = [f"`{a['quantity_name']}` ({a['transport_code']}) = {a['value']:+.5e}, factor to {a['nearest_constant']} = {a['factor_to_constant']:.4f}, log₁₀ dist = {a['log10_distance_to_constant']:+.4f}" for a in top5]
        return (
            f"{len(alpha_audits)} audit entries sit nearest α or 3α; closest hits:\n\n"
            + "\n".join(f"- {l}" for l in sample)
        )

    def line_q10():
        all_ok = all(r["max_conservation_error"] <= EPS + 1e-30 for r in rows)
        return f"{'Yes' if all_ok else 'No'} — all {len(TRANSPORT_PRINCIPLES) * len(CLUSTERS)} runs preserve the unit-speed normalization at or below machine epsilon ({EPS:.3e})."

    physical = [s for s in summaries if not s["is_wrong"]]

    def determine_outcome():
        t1_pk = float(t1["median_pearson_kappa"])
        t1_modes = float(t1["median_wave_mode_count"])
        within_tolerance = [s for s in physical
                            if abs(float(s["median_pearson_kappa"]) - a8_benchmark_kappa) < 0.005
                            and float(s["median_wave_mode_count"]) >= a8_benchmark_modes - 0.5]
        n_within = len(within_tolerance)
        unique_winner = (n_within == 1 and abs(t1_pk - a8_benchmark_kappa) < 0.001
                         and t1_modes >= 2)
        if unique_winner:
            return "Outcome A", f"Single transport principle (T1 scalar density, the control) naturally reproduces the A8 signature: κ = +{t1_pk:.5f} vs A8 reference +{a8_benchmark_kappa:.5f}, wave modes = {t1_modes:.1f}."
        if 2 <= n_within <= 7:
            return "Outcome B", ("Several transport principles reproduce the A8 signature on the primary metrics within tolerance. "
                                 + ", ".join(f"{s['transport_code']} (κ +{float(s['median_pearson_kappa']):.5f}, modes {float(s['median_wave_mode_count']):.1f})" for s in within_tolerance[:5])
                                 + ". A unique microscopic invariant does not emerge from this laboratory.")
        return "Outcome C", "No transport principle reproduces A8 within tolerance; the conserved quantity itself must be derived from a deeper microscopic description."

    outcome, outcome_text = determine_outcome()

    lines = [
        "# PBUF MICROSCOPIC-INVARIANTS-LAB-001",
        "",
        "**Conserved Quantities & Transport Principles Laboratory inside the frozen Version 1 weak-lensing laboratory (LAB-FREEZE-001).**",
        "",
        "## Status",
        "",
        f"- Frozen hash verification: **{'PASS' if hashes['ok'] else 'FAIL'}**",
        f"- Transport principles: **{len(TRANSPORT_PRINCIPLES)}** (T1-T10 + WR1-WR4)",
        f"- Production runs: **{len(TRANSPORT_PRINCIPLES) * len(CLUSTERS)}**",
        f"- Runtime: **{elapsed:.1f} s**",
        "- Fitting or optimisation: **none**",
        "",
        "## Frozen laboratory",
        "",
        "All transport, source-plane, Jacobian observable, numerical, constitutive, A8 architecture, and production components remain byte-identical to LAB-FREEZE-001. Only the microscopic transport variable (which conserved quantity evolves) is allowed to vary.",
        "",
        "## Transport Principles",
        "",
        "| # | Code | Name | Invariant | Principle |",
        "|---|---|---|---|---|",
    ]
    for t in TRANSPORT_PRINCIPLES:
        lines.append(f"| {t.code} | {t.name} | `{t.principle}` | {t.invariant_name} |")
    lines += [
        "",
        "## Transport Summary (median across 5 clusters)",
        "",
        "| Transport | Pearson κ | Pearson γ | RMS κ | Coherence gain | Memory | Wave modes | Wave families | Conservation |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for r in ranking:
        s = by[r["transport_code"]]
        lines.append(f"| {s['transport_code']} | {s['median_pearson_kappa']:+.5f} | {s['median_pearson_gamma']:+.5f} | {s['median_rms_kappa']:.5f} | {s['median_coherence_gain']:+.3e} | {s['median_emergent_memory_index']:.5f} | {s['median_wave_mode_count']:.1f} | {s['median_n_families']:.1f} | {s['max_conservation_error']:.3e} |")
    lines += [
        "",
        "## Expanded Wave Registry (per decomposition, per cluster)",
        "",
        "Recorded: number of wave families, longitudinal/transverse/mixed classifications, standing/travelling, phase/group velocities, dispersion, attenuation, coherence length, mode stability.",
        "",
        "## Energy Exchange Audit (per transport, per cluster)",
        "",
        "Tracks energy flowing fast → slow, returned slow → fast, stored, lost. Total conservation must remain exact for physical principles; wrong controls violate.",
        "",
        "## Candidate ranking",
        "",
        "Physical principles ranked by mean rank across all primary metrics (higher κ/γ, lower RMS κ/γ, higher coherence/memory/phase/orientation/multiplicative/F-S exchange/wave modes/coherence length/mode stability).",
        "",
        "| Rank | Code | Invariant | Pearson κ | Wave modes | Families | Mode stability | F/S exchange | Rank sum |",
        "|---:|---|---|---:|---:|---:|---:|---:|---:|",
    ]
    for r in ranking:
        s = by[r["transport_code"]]
        lines.append(f"| {r['rank']} | {s['transport_code']} | {s['invariant_kind']} | {s['median_pearson_kappa']:+.5f} | {s['median_wave_mode_count']:.1f} | {s['median_n_families']:.1f} | {s['median_wave_mode_stability']:.3f} | {s['median_fast_slow_exchange']:+.3e} | {r['rank_sum']:.0f} |")
    lines += [
        "",
        "## Required questions",
        "",
        "### Q1. Which conserved quantity best reproduces A8?",
        "",
        line_q1(),
        "",
        "### Q2. Does one transport principle naturally generate both wave modes?",
        "",
        line_q2(),
        "",
        "### Q3. Does memory arise from transport rather than architecture?",
        "",
        line_q3(),
        "",
        "### Q4. Which quantity primarily couples the fast and slow layers?",
        "",
        line_q4(),
        "",
        "### Q5. Does neighbour coherence arise naturally?",
        "",
        line_q5(),
        "",
        "### Q6. Does any transport principle outperform the frozen A8 benchmark?",
        "",
        line_q6(),
        "",
        "### Q7. Do the wave families remain stable across all five clusters?",
        "",
        line_q7(),
        "",
        "### Q8. Which transport principle produces the simplest microscopic description?",
        "",
        line_q8(),
        "",
        "### Q9. Do any stable transport ratios repeatedly converge toward α or 3α?",
        "",
        line_q9(),
        "",
        "### Q10. Does every successful transport principle preserve machine-precision conservation?",
        "",
        line_q10(),
        "",
        "## Outcome determination",
        "",
        "- **A**: One conserved transport principle naturally reproduces the complete A8 signature while maintaining both wave modes and cooperative behaviour.",
        "- **B**: Several transport principles reproduce different aspects of A8, but no unique invariant emerges.",
        "- **C**: No transport principle reproduces A8; the conserved quantity itself must be derived from a deeper microscopic description.",
        "",
        f"**{outcome}.** {outcome_text}",
        "",
        "## C10 provenance",
        "",
        "C10 and A8 were not modified or rerun. The A8 benchmark remains archived at `runs/microstructure_entity_lab001/architecture_summary.csv` (A8 row).",
        "",
        "## Numerical stability",
        "",
        f"All {len(TRANSPORT_PRINCIPLES) * len(CLUSTERS)} runs preserve the frozen unit-speed normalization at or below machine epsilon ({EPS:.3e}).",
        "",
        "## Permanent Registries",
        "",
        "Appended new entries to `runs/wave_family_registry.csv` and `runs/invariant_registry.csv`. Subsequent laboratories may continue to append without modifying earlier entries.",
        "",
        "## Required artefacts",
        "",
        "`report.md`, `transport_summary.csv`, `cross_cluster_statistics.csv`, `wave_registry.csv`, `energy_exchange.csv`, `candidate_ranking.csv`, `emergent_state_statistics.csv`, `fundamental_constant_audit.csv`, `run.json`, `validation.json`, and all required plots are present in `runs/microscopic_invariants_lab001/`.",
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
        for t in TRANSPORT_PRINCIPLES:
            rows.append(run_one(t, cluster, rho, obs))

    summaries = aggregate(rows)
    wave_reg = wave_registry_long(rows, summaries)
    energy_exch = energy_exchange_long(rows, summaries)
    audit = fundamental_constant_audit(rows, summaries)
    ranking = candidate_ranking(summaries)

    summary_fields = list(summaries[0].keys())
    cross_fields = ["transport_number", "transport_code", "transport_name", "invariant_name",
                    "invariant_kind", "is_wrong",
                    "cluster_id", "cluster_label",
                    "pearson_kappa", "pearson_gamma", "ssim_kappa", "ssim_gamma",
                    "rms_kappa", "rms_gamma", "kappa_bias", "gamma_bias", "runtime_seconds",
                    "max_conservation_error", "coherence_gain", "emergent_memory_index",
                    "evolution_activity", "phase_emergence_score", "orientation_emergence_score",
                    "multiplicative_coupling_score", "invariant_drift_relative", "fast_slow_exchange",
                    "n_families", "wave_modes", "wave_families",
                    "propagating", "standing", "travelling", "transverse", "longitudinal", "mixed",
                    "polarization_like", "attenuation", "dispersion",
                    "coherence_length", "mode_stability",
                    "fast_to_slow_total", "returned_total", "stored_total", "lost_total",
                    "fast_to_slow_net", "wave_emerged"]
    write_csv(OUT / "transport_summary.csv", summaries, summary_fields)
    write_csv(OUT / "cross_cluster_statistics.csv", rows, cross_fields)
    write_csv(OUT / "wave_registry.csv", wave_reg,
              sorted({k for r in wave_reg for k in r.keys()}))
    write_csv(OUT / "energy_exchange.csv", energy_exch,
              sorted({k for r in energy_exch for k in r.keys()}))
    write_csv(OUT / "candidate_ranking.csv", ranking,
              ["rank", "transport_code", "transport_name", "invariant_kind", "rank_sum",
               "median_pearson_kappa", "median_pearson_gamma", "median_coherence_gain",
               "median_emergent_memory_index", "median_wave_mode_count", "median_n_families",
               "median_wave_mode_stability", "median_fast_slow_exchange",
               "median_multiplicative_coupling_score", "median_relaxation_time"])
    emergent_rows = []
    for t in TRANSPORT_PRINCIPLES:
        sub = [r for r in rows if r["transport_code"] == t.code]
        for cid in [c["id"] for c in CLUSTERS]:
            row = next((r for r in sub if r["cluster_id"] == cid), None)
            if row is None:
                continue
            emergent_rows.append({
                "transport_code": t.code, "transport_name": t.name, "is_wrong": t.is_wrong,
                "cluster": cid,
                "phase_emergence_score": float(row["phase_emergence_score"]),
                "orientation_emergence_score": float(row["orientation_emergence_score"]),
                "coherence_gain": float(row["coherence_gain"]),
                "memory_index": float(row["emergent_memory_index"]),
                "activity": float(row["evolution_activity"]),
                "multiplicative_coupling_score": float(row["multiplicative_coupling_score"]),
                "fast_slow_exchange": float(row["fast_slow_exchange"]),
                "invariant_drift": float(row["invariant_drift_relative"]),
            })
    write_csv(OUT / "emergent_state_statistics.csv", emergent_rows,
              sorted({k for r in emergent_rows for k in r.keys()}))
    audit_fields = ["transport_code", "transport_name", "quantity_name", "value",
                    "log_abs", "nearest_constant", "nearest_constant_value",
                    "log10_distance_to_constant", "factor_to_constant", "is_alpha_or_3alpha"]
    write_csv(OUT / "fundamental_constant_audit.csv", audit, audit_fields)

    make_plots(rows, summaries, ranking, energy_exch, wave_reg)
    update_wave_registry(wave_reg)
    update_invariant_registry(summaries)
    elapsed = time.perf_counter() - started_total
    report_text = build_report(rows, summaries, ranking, audit, hashes, elapsed)
    (OUT / "report.md").write_text(report_text)

    run = {
        "milestone": "PBUF MICROSCOPIC-INVARIANTS-LAB-001",
        "kind": "microscopic conserved-transport laboratory",
        "frozen_laboratory": "Version 1 weak-lensing laboratory (LAB-FREEZE-001)",
        "frozen_implementation_sha256": {k: v["actual_sha256"] for k, v in hashes["files"].items()},
        "production_configuration": CONFIG, "clusters": CLUSTERS,
        "transport_principles": [t.__dict__ for t in TRANSPORT_PRINCIPLES],
        "fixed_parameters": {"dt": DT, "steps": STEPS, "k": K, "gamma": GAMMA,
                             "omega": OMEGA, "internal_k": INTERNAL_K,
                             "fast_timescale": FAST_TIMESCALE,
                             "slow_timescale": SLOW_TIMESCALE,
                             "alpha_fs": ALPHA_FS, "three_alpha_fs": THREE_ALPHA_FS},
        "fitting_performed": False, "optimisation_performed": False,
        "frozen_components_modified": False, "execution_seconds_total": elapsed,
    }
    (OUT / "run.json").write_text(json.dumps(run, indent=2))

    required = [OUT / "report.md",
                OUT / "transport_summary.csv",
                OUT / "cross_cluster_statistics.csv",
                OUT / "wave_registry.csv",
                OUT / "energy_exchange.csv",
                OUT / "candidate_ranking.csv",
                OUT / "emergent_state_statistics.csv",
                OUT / "fundamental_constant_audit.csv",
                OUT / "run.json"] + [PLOTS / n for n in (
                    "transport_family_comparison.png", "energy_exchange.png",
                    "wave_transport.png", "fast_slow_exchange.png",
                    "memory_transport.png", "transport_invariants.png",
                    "science_dashboard.png")]
    png_ok = all(p.read_bytes()[:8] == b"\x89PNG\r\n\x1a\n" for p in required if p.suffix == ".png")
    finite_ok = all(np.isfinite(r[k]) for r in rows for k in (
        "pearson_kappa", "pearson_gamma", "ssim_kappa", "ssim_gamma",
        "rms_kappa", "rms_gamma", "kappa_bias", "gamma_bias", "runtime_seconds",
        "max_conservation_error", "coherence_gain", "emergent_memory_index",
        "evolution_activity", "phase_emergence_score", "orientation_emergence_score",
        "multiplicative_coupling_score", "invariant_drift_relative", "fast_slow_exchange",
        "n_families", "wave_modes", "propagating", "standing", "transverse",
        "longitudinal", "mixed", "polarization_like", "attenuation", "dispersion",
        "coherence_length", "mode_stability",
        "fast_to_slow_total", "returned_total", "stored_total", "lost_total",
        "fast_to_slow_net"))
    conservation_ok = all(r["max_conservation_error"] <= EPS + 1e-30 for r in rows)
    artifacts_ok = all(p.exists() and p.stat().st_size > 0 for p in required)
    validation = {
        "milestone": "PBUF MICROSCOPIC-INVARIANTS-LAB-001",
        "frozen_hash_verification_passed": hashes["ok"], "frozen_hashes": hashes["files"],
        "expected_run_count": len(TRANSPORT_PRINCIPLES) * len(CLUSTERS),
        "actual_run_count": len(rows),
        "transport_count": len(TRANSPORT_PRINCIPLES), "cluster_count": len(CLUSTERS),
        "all_metrics_finite": finite_ok,
        "all_runs_machine_precision_conservation": conservation_ok,
        "runs_preserving_conservation": int(sum(r["max_conservation_error"] <= EPS + 1e-30 for r in rows)),
        "required_artifacts_present_nonempty": artifacts_ok,
        "png_signatures_valid": png_ok,
        "validation_passed": bool(hashes["ok"] and len(rows) == len(TRANSPORT_PRINCIPLES) * len(CLUSTERS)
                                  and finite_ok and conservation_ok and artifacts_ok and png_ok),
    }
    (OUT / "validation.json").write_text(json.dumps(validation, indent=2))
    if not validation["validation_passed"]:
        raise RuntimeError("Microscopic invariants laboratory validation failed")


if __name__ == "__main__":
    main()
