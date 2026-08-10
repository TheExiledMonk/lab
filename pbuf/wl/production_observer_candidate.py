"""Configuration-only Dev Doc 113 production observer candidate.

Generated selections are written by the compatibility lab to its result file;
this module intentionally does not alter the canonical observer default.
"""

PRODUCTION_OBSERVER_CANDIDATE = {
    "status": "WL_CHANNEL_SPECIFIC_OBSERVER_PARTIALLY_ESTABLISHED",
    "evidence": "runs/wl_channel_specific_compatibility001/result.json",
    "kappa": {
        "status": "WEAK_PREFERENCE",
        "deposition": "gaussian_sigma_half_cell",
        "channel_family": "all_except_depth_3d",
        "reconstruction": "nodepth_l1",
    },
    "gamma1": {"status": "TIED_SURVIVORS", "candidates": (
        ("tsc_3x3", "displacement_2d", "family_displacement_2d_signed"),
        ("gaussian_sigma_half_cell", "displacement_2d", "family_displacement_2d_signed"),
        ("bilinear_cic", "displacement_2d", "family_displacement_2d_signed"),
    )},
    "gamma2": {"status": "TIED_SURVIVORS", "candidates": (
        ("tsc_3x3", "density", "family_density_signed"),
        ("gaussian_sigma_half_cell", "density", "family_density_signed"),
        ("bilinear_cic", "density", "family_density_signed"),
    )},
}
