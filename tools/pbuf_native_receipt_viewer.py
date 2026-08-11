"""Read-only diagnostic viewer for frozen DEV177/DEV178 native receipt NPZ files."""
from __future__ import annotations
import argparse
import sys
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path: sys.path.insert(0, str(ROOT))
from pbuf.labs.foundation.native_received_j3_dev177 import fit_j3

def load_receipt(path: Path) -> dict[str, np.ndarray]:
    """Load arrays without mutation; undefined values remain undefined."""
    with np.load(path, allow_pickle=False) as data:
        return {name: data[name].copy() for name in data.files}

def main():
    parser = argparse.ArgumentParser(description="Diagnostic-only native receipt viewer; RGB has no physical semantics.")
    parser.add_argument("--lane", choices=("baseline", "25pct"), default="baseline")
    parser.add_argument("--realization", type=int, default=0, choices=range(8))
    parser.add_argument("--channels", default="position,direction", help="comma-separated: position,displacement,direction,momentum,flux,weight,W01,W02,W03,W04,progression,lineage")
    parser.add_argument("--summary", action="store_true")
    parser.add_argument("--j3", action="store_true", help="report the read-only J3/G3 diagnostic; undefined is never shown as zero")
    args = parser.parse_args()
    root = ROOT / ("runs/dev177_full_native_received_state" if args.lane == "baseline" else "runs/dev178_high_density_native_vulkan")
    values = load_receipt(root / f"receipt_realization_{args.realization:02d}.npz")
    if args.summary or args.j3:
        report = {"lane": args.lane, "realization": args.realization, "records": len(values["weights"]), "channels": args.channels.split(","), "read_only": True, "undefined_preserved": True}
        if args.j3: report["J3_G3"] = fit_j3(values["source_positions"], values["received_positions"])
        print(report)
        return
    import matplotlib.pyplot as plt
    channels = set(args.channels.split(",")); p = values["received_positions"]
    fig, ax = plt.subplots(); ax.scatter(p[:, 2], p[:, 1], s=2, c=values["weights"], cmap="viridis", label="received")
    if "position" in channels:
        s = values["source_positions"]; ax.scatter(s[:, 2], s[:, 1], s=2, c="tab:red", alpha=.25, label="source")
    mapping = {"displacement": "local_displacement", "direction": "directions", "momentum": "local_momentum", "flux": "local_flux"}
    for name, key in mapping.items():
        if name in channels:
            v = values[key]; ax.quiver(p[:, 2], p[:, 1], v[:, 2], v[:, 1], angles="xy", scale_units="xy", scale=1, alpha=.35, label=name)
    ax.set(xlabel="native z", ylabel="native y", title="Diagnostic native receipt view — no observer/physical RGB semantics"); ax.legend(); plt.show()

if __name__ == "__main__": main()
