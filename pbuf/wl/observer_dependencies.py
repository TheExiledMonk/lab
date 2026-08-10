"""Audited dependency inventory for the frozen 45-channel observer bank."""

from pbuf.labs.foundation import native_full_received_state_information_retention001 as RET


def describe_observer_dependencies() -> dict:
    out = {}
    density = {"histogram_density": ("screen_coordinate", "deposition"),
               "kernel_density": ("screen_coordinate", "pairwise_kde", "deposition"),
               "knn_density": ("screen_coordinate", "pairwise_kde", "deposition")}
    for method in RET.EXTRACTION_METHODS:
        if method in density: deps = density[method]
        elif method in ("jacobian_affine", "polar_jacobian", "displacement_divergence"):
            deps = ("screen_coordinate", "local_differential_jacobian")
        else: deps = ("screen_coordinate", "covariance")
        for field in RET.EXTRACTION_FIELDS:
            name = f"{method}__{field}"
            out[name] = {"dependencies": deps + ("channel_assembly",),
                "uses_absolute_u": True, "uses_absolute_v": True,
                "uses_pairwise_distance": method in ("kernel_density", "knn_density"),
                "uses_deposition": method in density,
                "uses_values": ("u0", "v0", "uf", "vf")}
    for name in RET.PRIMARY_3D_BIN_CHANNELS:
        dep = "local_differential_jacobian" if name.startswith("j3_") else (
              "covariance" if name.startswith(("std_", "cov_")) else "received_state")
        out[name] = {"dependencies": ("received_state", "screen_coordinate", dep,
                                      "channel_assembly"),
            "uses_absolute_u": True, "uses_absolute_v": True,
            "uses_pairwise_distance": False, "uses_deposition": False,
            "uses_values": (name,)}
    if len(out) != 45: raise RuntimeError(f"expected 45 observer channels, got {len(out)}")
    return out


KDE_UNIFORM_TRANSLATION_INVARIANT = True
