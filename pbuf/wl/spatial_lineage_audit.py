"""Coordinate-transform provenance primitives for the Dev136 spatial audit."""
from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Iterable


@dataclass(frozen=True)
class SpatialTransform:
    stage: str
    input_coordinate: str
    input_units: str
    operation: str
    output_coordinate: str
    output_units: str
    invertible: bool
    scale_metadata_retained: bool
    physical_provenance_retained: bool


def normalization_ledger(transforms: Iterable[SpatialTransform]) -> list[dict]:
    return [asdict(item) for item in transforms]


def coordinate_lineage(nodes: Iterable[str], transforms: Iterable[SpatialTransform],
                       start: str = "source_grid", end: str = "receiver_coordinates") -> dict:
    """Return a deterministic graph and whether *end* is reachable from *start*."""
    node_list = list(dict.fromkeys(nodes))
    edges = [asdict(item) for item in transforms]
    reached = {start}
    changed = True
    while changed:
        changed = False
        for edge in edges:
            if edge["input_coordinate"] in reached and edge["output_coordinate"] not in reached:
                reached.add(edge["output_coordinate"]); changed = True
    complete = end in reached
    return {"nodes": [{"id": n} for n in node_list], "edges": edges,
            "start": start, "end": end, "complete": complete,
            "broken_spatial_lineage_edges": 0 if complete else 1,
            "reachable_nodes": sorted(reached)}
