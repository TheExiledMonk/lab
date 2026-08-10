"""In-process cache for mathematically identical observer primitives."""

from dataclasses import dataclass
import hashlib
import numpy as np


@dataclass(frozen=True)
class ObserverStateId:
    base_state: str
    coordinate_transform: str = "identity"
    backend: str = "cpu"


def array_identity(value) -> tuple:
    """Content identity: safe even when a caller later mutates its array."""
    a = np.ascontiguousarray(value)
    return (a.shape, a.dtype.str, hashlib.sha256(a.view(np.uint8)).hexdigest())


class ObserverPrimitiveCache:
    def __init__(self, profile=None):
        self._values = {}; self.profile = profile

    @staticmethod
    def key(primitive, state_id, *, coordinates=(), values=None, parameters=(),
            translation_invariant=False):
        sid = state_id
        if translation_invariant:
            sid = ObserverStateId(state_id.base_state, "translation_invariant", state_id.backend)
        coordinate_key = () if translation_invariant else tuple(array_identity(x) for x in coordinates)
        return (primitive, sid, coordinate_key,
                None if values is None else array_identity(values), tuple(parameters))

    def get_or_compute(self, key, compute, category="other"):
        hit = key in self._values
        if self.profile: self.profile.cache(category, hit)
        if not hit: self._values[key] = compute()
        return self._values[key]

    def clear(self): self._values.clear()
