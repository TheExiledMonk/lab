"""Lossless additive packet/launch provenance for :class:`NativeReceivedState`."""
from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
import numpy as np

from pbuf.excitation.native_finite_receipt import NativeReceivedState


def _json_native(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, dict):
        return {key: _json_native(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [_json_native(item) for item in value]
    return value


@dataclass(frozen=True)
class PacketAwareReceiptCollection:
    """One or more unmodified DEV168 receipts plus compact aligned lineage indices."""
    native_received_state: NativeReceivedState
    receipt_launch_index: np.ndarray
    receipt_packet_index: np.ndarray
    receipt_realization_index: np.ndarray
    launch_manifest: tuple[dict, ...]
    packet_manifest: tuple[dict, ...]
    realization_manifest: tuple[dict, ...]

    def __post_init__(self):
        n = len(self.native_received_state.weights)
        arrays = ("receipt_launch_index", "receipt_packet_index", "receipt_realization_index")
        for name in arrays:
            value = np.asarray(getattr(self, name), dtype=np.int32)
            if value.shape != (n,):
                raise ValueError(f"{name} must have shape (N,)")
            object.__setattr__(self, name, value)
        for name, manifest, index in (("launch", self.launch_manifest, self.receipt_launch_index), ("packet", self.packet_manifest, self.receipt_packet_index), ("realization", self.realization_manifest, self.receipt_realization_index)):
            if not manifest and n:
                raise ValueError(f"{name} manifest is required for nonempty receipts")
            if np.any(index < 0) or np.any(index >= len(manifest)):
                raise ValueError(f"unknown {name} index")
        if len({row["launch_id"] for row in self.launch_manifest}) != len(self.launch_manifest):
            raise ValueError("launch IDs must be unique")
        if len({row["packet_id"] for row in self.packet_manifest}) != len(self.packet_manifest):
            raise ValueError("packet IDs must be unique")
        if len({row["realization_id"] for row in self.realization_manifest}) != len(self.realization_manifest):
            raise ValueError("realization IDs must be unique")

    def arrays(self) -> dict[str, np.ndarray]:
        return {**self.native_received_state.arrays(), "receipt_launch_index": self.receipt_launch_index,
                "receipt_packet_index": self.receipt_packet_index, "receipt_realization_index": self.receipt_realization_index}

    def write(self, npz_path: Path, manifest_path: Path) -> None:
        np.savez_compressed(npz_path, **self.arrays())
        metadata = {"schema": "NativePacketReceiptCollection/v1", "representation": self.native_received_state.representation,
                    "launch_manifest": self.launch_manifest, "packet_manifest": self.packet_manifest,
                    "realization_manifest": self.realization_manifest}
        manifest_path.write_text(json.dumps(_json_native(metadata), sort_keys=True, indent=2) + "\n")

    @classmethod
    def read(cls, npz_path: Path, manifest_path: Path) -> "PacketAwareReceiptCollection":
        meta = json.loads(manifest_path.read_text())
        if meta.get("schema") != "NativePacketReceiptCollection/v1":
            raise ValueError("unknown packet-aware receipt schema")
        with np.load(npz_path, allow_pickle=False) as z:
            missing = set(NativeReceivedState.__dataclass_fields__) - {"representation"} - set(z.files)
            required = {"receipt_launch_index", "receipt_packet_index", "receipt_realization_index"}
            if missing or not required <= set(z.files):
                raise ValueError("packet-aware receipt arrays are incomplete")
            state = NativeReceivedState(**{name: z[name].copy() for name in NativeReceivedState.__dataclass_fields__ if name != "representation"}, representation=meta["representation"])
            return cls(state, z["receipt_launch_index"].copy(), z["receipt_packet_index"].copy(), z["receipt_realization_index"].copy(),
                       tuple(meta["launch_manifest"]), tuple(meta["packet_manifest"]), tuple(meta["realization_manifest"]))
