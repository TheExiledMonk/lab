#!/usr/bin/env python3
"""Narrow runner correction for G3D angular moment decomposition 001.

This wrapper changes no science. It fixes only the step-0 gate semantics for
centroid/spread fractions when the total angular second-moment energy is exactly
zero. At that checkpoint the fractions 0/0 are mathematically undefined and are
therefore represented as NaN by the underlying lab. The original gate treated
"no defined fractions" as a failure even when centroid, spread, and total energy
were all zero.

The corrected rule is:
- where total energy is non-zero, fractions must be finite, sum to one, and lie
  in [0,1] within the existing tolerance;
- where total energy is zero, both fractions must remain undefined and the
  centroid/spread energies must also be zero within the existing threshold.

All other identities, physics, propagation, observer coordinates, tolerances,
outputs, and authorization flags remain exactly those of the underlying 001 lab.
"""
from __future__ import annotations

import numpy as np

import pbuf.labs.foundation.g3d_angular_moment_decomposition001 as LAB


def _corrected_decomposition_gates(base: dict, f: dict) -> dict:
    total = f["total_angular_second_moment_energy"]
    energy_identity = LAB._safe_rel_rms(
        f["centroid_energy"] + f["spread_energy"] - total,
        total,
    )
    quadrature_identity = LAB._safe_rel_rms(
        f["quadrature_total_angle_mag"] - base["angular_rms_angle_mag"],
        base["angular_rms_angle_mag"],
    )

    centroid_fraction = f["centroid_energy_fraction"]
    spread_fraction = f["spread_energy_fraction"]
    fraction_sum = centroid_fraction + spread_fraction
    finite_frac = np.isfinite(fraction_sum)
    fraction_defined_count = int(np.sum(finite_frac))

    # Fractions are defined only where total energy is non-zero.
    fraction_sum_err = (
        float(np.max(np.abs(fraction_sum[finite_frac] - 1.0)))
        if fraction_defined_count
        else 0.0
    )

    frac_values = np.concatenate([
        centroid_fraction[np.isfinite(centroid_fraction)],
        spread_fraction[np.isfinite(spread_fraction)],
    ])
    fraction_min = float(np.min(frac_values)) if frac_values.size else float("nan")
    fraction_max = float(np.max(frac_values)) if frac_values.size else float("nan")

    finite_total = np.isfinite(total)
    nonzero_total = finite_total & (np.abs(total) > 1e-30)
    zero_total = finite_total & ~nonzero_total

    # Every non-zero-energy bin must have both fractions defined.
    nonzero_fraction_defined_pass = bool(
        np.all(np.isfinite(centroid_fraction[nonzero_total]))
        and np.all(np.isfinite(spread_fraction[nonzero_total]))
    )

    # At exact zero energy, 0/0 fractions are undefined by construction.
    # That is valid only if both physical decomposition energies are also zero.
    zero_fraction_undefined_pass = bool(
        np.any(zero_total)
        and np.all(~np.isfinite(centroid_fraction[zero_total]))
        and np.all(~np.isfinite(spread_fraction[zero_total]))
        and np.all(np.abs(f["centroid_energy"][zero_total]) <= 1e-30)
        and np.all(np.abs(f["spread_energy"][zero_total]) <= 1e-30)
        and np.all(np.abs(total[zero_total]) <= 1e-30)
    ) if np.any(zero_total) else True

    if frac_values.size:
        numerical_bounds_pass = bool(
            fraction_min >= -LAB.FRACTION_TOL
            and fraction_max <= 1.0 + LAB.FRACTION_TOL
        )
    else:
        # Step 0 is expected to land here: all supported bins have zero angular
        # energy, hence both fractions are intentionally undefined.
        numerical_bounds_pass = bool(np.any(zero_total))

    fraction_bounds_pass = bool(
        numerical_bounds_pass
        and nonzero_fraction_defined_pass
        and zero_fraction_undefined_pass
    )

    mxx = base["angular_second_moment_xx"]
    mxy = base["angular_second_moment_xy"]
    myy = base["angular_second_moment_yy"]
    rxx = f["centroid_outer_xx"] + base["angular_cov_xx"]
    rxy = f["centroid_outer_xy"] + base["angular_cov_xy"]
    ryy = f["centroid_outer_yy"] + base["angular_cov_yy"]
    tensor_diff = np.sqrt((rxx-mxx)**2 + 2.0*(rxy-mxy)**2 + (ryy-myy)**2)
    tensor_ref = f["full_second_moment_tensor_frobenius_mag"]
    tensor_identity = LAB._safe_rel_rms(tensor_diff, tensor_ref)

    fro_sq_rhs = (
        f["centroid_tensor_frobenius_mag"]**2
        + f["spread_tensor_frobenius_mag"]**2
        + 2.0*f["centroid_spread_tensor_inner_product"]
    )
    fro_sq_lhs = f["full_second_moment_tensor_frobenius_mag"]**2
    fro_identity = LAB._safe_rel_rms(fro_sq_lhs - fro_sq_rhs, fro_sq_lhs)

    return {
        "centroid_plus_spread_energy_identity_relative_rms_error": energy_identity,
        "quadrature_total_angle_vs_prior_rms_angle_relative_rms_error": quadrature_identity,
        "centroid_plus_spread_fraction_max_abs_error": fraction_sum_err,
        "energy_fraction_min": fraction_min,
        "energy_fraction_max": fraction_max,
        "energy_fraction_bounds_pass": fraction_bounds_pass,
        "energy_fraction_defined_count": fraction_defined_count,
        "zero_total_fraction_undefined_pass": zero_fraction_undefined_pass,
        "nonzero_total_fraction_defined_pass": nonzero_fraction_defined_pass,
        "second_moment_tensor_decomposition_relative_rms_error": tensor_identity,
        "second_moment_tensor_frobenius_identity_relative_rms_error": fro_identity,
    }


def main() -> int:
    LAB._decomposition_gates = _corrected_decomposition_gates
    return LAB.main()


if __name__ == "__main__":
    raise SystemExit(main())
