"""PBUF core modules — verified numerical core (FOUNDATION-001).

Imports are deferred to attribute access to avoid partial-init failures
during test bootstrap.
"""
__all__ = [
    "conventions", "coordinate_transforms", "vector_transforms",
    "tensor_transforms", "pair_enumeration", "field_diagnostics",
    "differential_operators", "helmholtz_3d", "los_projection",
    "observable_extraction", "pair_transfer", "midpoint_rasterization",
    "ray_interface",
]


def __getattr__(name):
    import importlib
    if name in __all__:
        return importlib.import_module("." + name, __name__)
    raise AttributeError(name)