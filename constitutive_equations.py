"""Interchangeable PBUF constitutive equations for the WL-002 laboratory."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import numpy as np


@dataclass(frozen=True)
class Equation:
    version: str
    description: str
    response: str
    stiffness: str
    solve: Callable[[np.ndarray, object], np.ndarray]
    formula: str = ""
    motivation: str = ""
    assumptions: tuple[str, ...] = ()
    strengths: tuple[str, ...] = ()
    weaknesses: tuple[str, ...] = ()

    @property
    def equation_id(self) -> str:
        return f"PBUF-VERSION-{self.version}"


def _normalized(matter: np.ndarray) -> np.ndarray:
    return matter / max(float(matter.max()), 1e-15)


def version_a(matter: np.ndarray, c: object) -> np.ndarray:
    """WL-001 baseline: local linear loading."""
    return c.deformation_strength * _normalized(matter)


def version_b(matter: np.ndarray, c: object) -> np.ndarray:
    """Local quadratic loading: dilute matter deforms the medium less strongly."""
    loading = _normalized(matter)
    return c.deformation_strength * loading**2


def version_c(matter: np.ndarray, c: object) -> np.ndarray:
    """Local nonlinear response with rigidity increasing with matter loading."""
    loading = _normalized(matter)
    # K/K0=(1+loading)/2 keeps the documented peak response unchanged.
    return c.deformation_strength * 2.0 * loading / (1.0 + loading)


def version_d(matter: np.ndarray, c: object) -> np.ndarray:
    """Nonlocal elastic response to quadratic loading via a Helmholtz kernel."""
    source = version_b(matter, c)
    spacing = 2.0 * c.extent / (matter.shape[0] - 1)
    ky = 2.0 * np.pi * np.fft.fftfreq(matter.shape[0], d=spacing)
    kx = 2.0 * np.pi * np.fft.fftfreq(matter.shape[1], d=spacing)
    KX, KY = np.meshgrid(kx, ky, indexing="xy")
    # The propagation length is the observed baryonic width, not a fitted scale.
    transfer = 1.0 / (1.0 + c.mass_sigma**2 * (KX**2 + KY**2))
    return np.fft.ifft2(np.fft.fft2(source) * transfer).real


EQUATIONS = {
    "A": Equation("A", "deformation = strength * normalized baryonic matter", "local linear", "constant", version_a,
        "u = u0 rho/rho_max", "Minimal scalar, isotropic local response and the WL-001 limiting law.",
        ("deformation is scalar", "response is instantaneous and local"),
        ("fewest assumptions", "exact WL-001 baseline"), ("cannot propagate deformation",)),
    "B": Equation("B", "quadratic local matter loading", "local nonlinear", "constant", version_b,
        "u = u0 (rho/rho_max)^2", "Tests nonlinear loading without changing the medium model.",
        ("deformation is scalar", "quadratic loading is physically selected"),
        ("stable", "suppresses dilute loading"), ("quadratic exponent is not derived from PBUF", "no propagation")),
    "C": Equation("C", "matter-dependent rigidity K/K0=(1+loading)/2", "local nonlinear", "increases with loading", version_c,
        "u = 2 u0 q/(1+q), q=rho/rho_max", "A local compliance law tests a medium that stiffens with loading.",
        ("local equilibrium", "K/K0=(1+q)/2"),
        ("bounded", "recovers the baseline at zero and peak loading"), ("rigidity interpolation is postulated", "no propagation")),
    "D": Equation("D", "Helmholtz-propagated quadratic loading; propagation length = observed mass width", "nonlocal linear propagation of nonlinear loading", "constant propagation rigidity", version_d,
        "(1 - sigma_rho^2 Laplacian) u = u0 (rho/rho_max)^2",
        "Distributed elastic recovery balances local loading; the observed baryonic width supplies the only length scale.",
        ("scalar isotropic medium", "periodic numerical boundary", "recovery is distributed"),
        ("propagates deformation", "stable positive spectral response", "no fitted length"),
        ("scalar proxy cannot represent shear stress", "periodic boundary is a laboratory approximation")),
}


def get_equation(version: str) -> Equation:
    try:
        return EQUATIONS[version.upper()]
    except KeyError as exc:
        raise ValueError(f"unknown equation version {version!r}; choose from {', '.join(EQUATIONS)}") from exc
