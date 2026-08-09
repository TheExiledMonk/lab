"""Lightweight, run-local observer timing and cache instrumentation."""

from collections import defaultdict
from contextlib import contextmanager
import time

CATEGORIES = ("pairwise_kde", "jacobian", "covariance", "deposition",
              "channel_assembly", "reconstruction", "other")


class ObserverProfile:
    def __init__(self):
        self._rows = defaultdict(lambda: {"call_count": 0, "cache_hit_count": 0,
                                         "cache_miss_count": 0, "total_seconds": 0.0})
        for name in CATEGORIES:
            self._rows[name]

    @contextmanager
    def time(self, category):
        if category not in CATEGORIES:
            raise ValueError(f"unsupported observer timing category: {category}")
        started = time.perf_counter(); self._rows[category]["call_count"] += 1
        try: yield
        finally: self._rows[category]["total_seconds"] += time.perf_counter() - started

    def cache(self, category, hit):
        self._rows[category]["cache_hit_count" if hit else "cache_miss_count"] += 1

    def describe(self):
        return {name: dict(self._rows[name]) for name in CATEGORIES}
