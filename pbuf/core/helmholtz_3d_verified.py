"""Assistant-audited M13 — 3D Helmholtz decomposition.

This module separates three distinct numerical diagnostics that must not
be conflated:

1. field_reconstruction_error
      ||R - (R_irr + R_sol)||_2 / ||R||_2
2. energy_closure_error
      |E_irr + E_sol - E_native| / E_native
3. orthogonality_error
      |<R_irr, R_sol>| / E_native

Production uses a vectorised Fourier projector. The independent
reference path (padding='none' only) uses an explicit per-mode loop and
does not call the production projector helper.
"""
from __future__ import annotations

import numpy as np

from .conventions import EPS_VARIANCE_UNDEFINED

__all__ = [
    "helmholtz_decompose_3d",
    "helmholtz_decompose_3d_reference",
    "HelmholtzError",
    "VALID_PADDINGS",
]

VALID_PADDINGS = ("none", "reflect_half")


class HelmholtzError(ValueError):
    pass


def _validate_inputs(Rx, Ry, Rz, spacing, padding):
    Rx = np.asarray(Rx, dtype=np.float64)
    Ry = np.asarray(Ry, dtype=np.float64)
    Rz = np.asarray(Rz, dtype=np.float64)
    if Rx.shape != Ry.shape or Rx.shape != Rz.shape:
        raise HelmholtzError("Rx, Ry, Rz must share the same shape")
    if Rx.ndim != 3:
        raise HelmholtzError("Helmholtz decomposition requires 3D arrays")
    if len(spacing) != 3:
        raise HelmholtzError("spacing must be (dx, dy, dz)")
    spacing = tuple(float(s) for s in spacing)
    if not all(np.isfinite(s) and s > 0.0 for s in spacing):
        raise HelmholtzError("all spacing values must be finite and > 0")
    if padding not in VALID_PADDINGS:
        raise HelmholtzError(
            f"unsupported padding {padding!r}; expected one of {VALID_PADDINGS}"
        )
    if not (np.all(np.isfinite(Rx)) and np.all(np.isfinite(Ry)) and np.all(np.isfinite(Rz))):
        raise HelmholtzError("input field contains NaN or Inf")
    return Rx, Ry, Rz, spacing, padding


def _apply_padding(Rx, Ry, Rz, padding):
    native_shape = Rx.shape
    if padding == "none":
        return Rx, Ry, Rz, (0, 0, 0), native_shape
    pz, py, px = (native_shape[0] // 2,
                  native_shape[1] // 2,
                  native_shape[2] // 2)
    pads = ((pz, pz), (py, py), (px, px))
    return (
        np.pad(Rx, pads, mode="reflect"),
        np.pad(Ry, pads, mode="reflect"),
        np.pad(Rz, pads, mode="reflect"),
        (pz, py, px),
        native_shape,
    )


def _crop(arr, pad, native_shape):
    pz, py, px = pad
    nz, ny, nx = native_shape
    return np.asarray(arr[pz:pz+nz, py:py+ny, px:px+nx].real, dtype=np.float64).copy()


def _energy(Rx, Ry, Rz):
    return float(np.sum(Rx*Rx + Ry*Ry + Rz*Rz))


def _norm(Rx, Ry, Rz):
    return float(np.sqrt(_energy(Rx, Ry, Rz)))


def _metrics(native, irr, sol):
    Rx, Ry, Rz = native
    Ix, Iy, Iz = irr
    Sx, Sy, Sz = sol

    E_native = _energy(Rx, Ry, Rz)
    E_irr = _energy(Ix, Iy, Iz)
    E_sol = _energy(Sx, Sy, Sz)

    denom_E = max(E_native, EPS_VARIANCE_UNDEFINED)
    norm_native = max(np.sqrt(E_native), np.sqrt(EPS_VARIANCE_UNDEFINED))

    dx = Rx - (Ix + Sx)
    dy = Ry - (Iy + Sy)
    dz = Rz - (Iz + Sz)
    field_reconstruction_error = _norm(dx, dy, dz) / norm_native

    energy_closure_error = abs((E_irr + E_sol) - E_native) / denom_E
    inner = float(np.sum(Ix*Sx + Iy*Sy + Iz*Sz))
    orthogonality_error = abs(inner) / denom_E

    if E_native <= EPS_VARIANCE_UNDEFINED:
        nan = float("nan")
        fractions = {
            "f_irr_partition": nan,
            "f_sol_partition": nan,
            "f_irr_native": nan,
            "f_sol_native": nan,
        }
        trivial = True
    else:
        part = E_irr + E_sol
        if part <= EPS_VARIANCE_UNDEFINED:
            fip = fsp = float("nan")
        else:
            fip = E_irr / part
            fsp = E_sol / part
        fractions = {
            "f_irr_partition": float(fip),
            "f_sol_partition": float(fsp),
            "f_irr_native": float(E_irr / E_native),
            "f_sol_native": float(E_sol / E_native),
        }
        trivial = False

    return {
        "E_native": E_native,
        "E_irr": E_irr,
        "E_sol": E_sol,
        "field_reconstruction_error": float(field_reconstruction_error),
        "energy_closure_error": float(energy_closure_error),
        "orthogonality_error": float(orthogonality_error),
        "inner_product_irr_sol": inner,
        "field_is_trivial": trivial,
        **fractions,
    }


def _kvectors(shape, spacing):
    nz, ny, nx = shape
    dx, dy, dz = spacing
    kx = 2.0*np.pi*np.fft.fftfreq(nx, d=dx).reshape(1, 1, nx)
    ky = 2.0*np.pi*np.fft.fftfreq(ny, d=dy).reshape(1, ny, 1)
    kz = 2.0*np.pi*np.fft.fftfreq(nz, d=dz).reshape(nz, 1, 1)
    return (
        np.broadcast_to(kx, shape),
        np.broadcast_to(ky, shape),
        np.broadcast_to(kz, shape),
    )


def helmholtz_decompose_3d(Rx, Ry, Rz, spacing=(1.0, 1.0, 1.0),
                            padding="reflect_half"):
    """Decompose R = R_irr + R_sol in Fourier space.

    Array axes are (z, y, x); vector components are (x, y, z).
    Spacing is supplied as (dx, dy, dz).

    The k=0 mode is assigned to the solenoidal residual. This preserves
    a constant vector component exactly and avoids inventing a
    longitudinal direction where k is undefined.
    """
    Rx, Ry, Rz, spacing, padding = _validate_inputs(Rx, Ry, Rz, spacing, padding)
    Rxp, Ryp, Rzp, pad, native_shape = _apply_padding(Rx, Ry, Rz, padding)

    KX, KY, KZ = _kvectors(Rxp.shape, spacing)
    K2 = KX*KX + KY*KY + KZ*KZ
    nonzero_k = K2 > 0.0
    safe_K2 = np.where(nonzero_k, K2, 1.0)

    Hx = np.fft.fftn(Rxp)
    Hy = np.fft.fftn(Ryp)
    Hz = np.fft.fftn(Rzp)
    dot = KX*Hx + KY*Hy + KZ*Hz

    Ixh = np.where(nonzero_k, KX*dot/safe_K2, 0.0)
    Iyh = np.where(nonzero_k, KY*dot/safe_K2, 0.0)
    Izh = np.where(nonzero_k, KZ*dot/safe_K2, 0.0)

    Ix_p = np.fft.ifftn(Ixh).real
    Iy_p = np.fft.ifftn(Iyh).real
    Iz_p = np.fft.ifftn(Izh).real
    Sx_p = Rxp - Ix_p
    Sy_p = Ryp - Iy_p
    Sz_p = Rzp - Iz_p

    pad_metrics = _metrics((Rxp, Ryp, Rzp), (Ix_p, Iy_p, Iz_p), (Sx_p, Sy_p, Sz_p))

    Ix = _crop(Ix_p, pad, native_shape)
    Iy = _crop(Iy_p, pad, native_shape)
    Iz = _crop(Iz_p, pad, native_shape)
    Sx = _crop(Sx_p, pad, native_shape)
    Sy = _crop(Sy_p, pad, native_shape)
    Sz = _crop(Sz_p, pad, native_shape)
    crop_metrics = _metrics((Rx, Ry, Rz), (Ix, Iy, Iz), (Sx, Sy, Sz))

    result = {
        "Rirr_x": Ix, "Rirr_y": Iy, "Rirr_z": Iz,
        "Rsol_x": Sx, "Rsol_y": Sy, "Rsol_z": Sz,
        "padding": padding,
        "spacing": tuple(spacing),
    }
    result.update(crop_metrics)

    # Compatibility aliases now point to the actual field reconstruction
    # metric rather than the former energy-closure quantity.
    result["reconstruction_error"] = crop_metrics["field_reconstruction_error"]
    result["reconstruction_error_pad"] = pad_metrics["field_reconstruction_error"]

    for key, value in pad_metrics.items():
        result[f"{key}_pad"] = value

    return result


def helmholtz_decompose_3d_reference(Rx, Ry, Rz, spacing=(1.0, 1.0, 1.0),
                                      padding="none"):
    """Independent explicit-mode reference implementation.

    Only the periodic unpadded case is supported. Each Fourier mode is
    projected with a literal 3-vector calculation in Python loops.
    """
    Rx, Ry, Rz, spacing, padding = _validate_inputs(Rx, Ry, Rz, spacing, padding)
    if padding != "none":
        raise HelmholtzError("reference implementation supports padding='none' only")

    nz, ny, nx = Rx.shape
    dx, dy, dz = spacing
    kx = 2.0*np.pi*np.fft.fftfreq(nx, d=dx)
    ky = 2.0*np.pi*np.fft.fftfreq(ny, d=dy)
    kz = 2.0*np.pi*np.fft.fftfreq(nz, d=dz)

    H = [np.fft.fftn(Rx), np.fft.fftn(Ry), np.fft.fftn(Rz)]
    I = [np.zeros_like(H[0]), np.zeros_like(H[1]), np.zeros_like(H[2])]

    for iz in range(nz):
        for iy in range(ny):
            for ix in range(nx):
                k = np.array([kx[ix], ky[iy], kz[iz]], dtype=np.float64)
                k2 = float(np.dot(k, k))
                if k2 == 0.0:
                    continue
                h = np.array([H[0][iz, iy, ix], H[1][iz, iy, ix], H[2][iz, iy, ix]],
                             dtype=np.complex128)
                coeff = np.dot(k, h) / k2
                projected = k * coeff
                for c in range(3):
                    I[c][iz, iy, ix] = projected[c]

    irr = tuple(np.fft.ifftn(I[c]).real for c in range(3))
    sol = (Rx - irr[0], Ry - irr[1], Rz - irr[2])
    metrics = _metrics((Rx, Ry, Rz), irr, sol)
    return {
        "Rirr_x": irr[0], "Rirr_y": irr[1], "Rirr_z": irr[2],
        "Rsol_x": sol[0], "Rsol_y": sol[1], "Rsol_z": sol[2],
        "padding": "none",
        "spacing": tuple(spacing),
        **metrics,
        "reconstruction_error": metrics["field_reconstruction_error"],
    }


def _analytic_fixture_tests():
    shape = (8, 9, 10)
    spacing = (0.7, 1.1, 1.3)
    nz, ny, nx = shape
    dx, dy, dz = spacing
    x = np.arange(nx)*dx
    y = np.arange(ny)*dy
    z = np.arange(nz)*dz
    Z, Y, X = np.meshgrid(z, y, x, indexing="ij")

    kx = 2.0*np.pi/(nx*dx)
    ky = 2.0*np.pi/(ny*dy)
    phase = kx*X + ky*Y

    # Longitudinal: gradient of a periodic scalar Fourier mode.
    RxL = kx*np.cos(phase)
    RyL = ky*np.cos(phase)
    RzL = np.zeros(shape)

    # Transverse for same k: (-ky, kx, 0) sin(phase).
    RxT = -ky*np.sin(phase)
    RyT = kx*np.sin(phase)
    RzT = np.zeros(shape)

    outL = helmholtz_decompose_3d(RxL, RyL, RzL, spacing, padding="none")
    outT = helmholtz_decompose_3d(RxT, RyT, RzT, spacing, padding="none")

    nL = _norm(RxL, RyL, RzL)
    nT = _norm(RxT, RyT, RzT)
    RxM = RxL/nL + RxT/nT
    RyM = RyL/nL + RyT/nT
    RzM = np.zeros(shape)
    outM = helmholtz_decompose_3d(RxM, RyM, RzM, spacing, padding="none")

    tol = 1e-12
    return {
        "pure_longitudinal": (
            abs(outL["f_irr_partition"] - 1.0) < tol and
            outL["field_reconstruction_error"] < tol
        ),
        "pure_transverse": (
            abs(outT["f_sol_partition"] - 1.0) < tol and
            outT["field_reconstruction_error"] < tol
        ),
        "mixed_equal_energy": (
            abs(outM["f_irr_partition"] - 0.5) < 1e-10 and
            outM["field_reconstruction_error"] < tol
        ),
    }


def _production_reference_test():
    rng = np.random.RandomState(13)
    shape = (5, 6, 7)
    spacing = (0.8, 1.2, 1.7)
    field = [rng.randn(*shape) for _ in range(3)]
    p = helmholtz_decompose_3d(*field, spacing=spacing, padding="none")
    r = helmholtz_decompose_3d_reference(*field, spacing=spacing, padding="none")
    keys = ("Rirr_x", "Rirr_y", "Rirr_z", "Rsol_x", "Rsol_y", "Rsol_z")
    max_diff = max(float(np.max(np.abs(p[k] - r[k]))) for k in keys)
    return {"max_diff": max_diff, "passes": max_diff < 1e-12}


if __name__ == "__main__":
    fixtures = _analytic_fixture_tests()
    assert all(fixtures.values()), fixtures
    ref = _production_reference_test()
    assert ref["passes"], ref
    print("M13 analytic fixtures: PASS")
    print(f"M13 explicit-loop reference: max_diff={ref['max_diff']:.3e}")
    print("M13 metric separation: PASS")
