"""Plot generation for the macro-micro response-bridge diagnostic lab.

This module is imported by the main driver.  It produces every required
plot under ``runs/macro_micro_response_bridge_diagnostic_lab001/plots/``.
"""
from __future__ import annotations

import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def _safe_scalar(x):
    return x if (x is not None and math.isfinite(x)) else np.nan


def _four_panel(out_path, panels, title, cmap="viridis", symmetric=False):
    n = len(panels)
    fig, axes = plt.subplots(1, n, figsize=(5 * n, 5))
    if n == 1:
        axes = [axes]
    for ax, (lbl, f) in zip(axes, panels):
        if f is None:
            ax.set_title(f"{lbl} - unavailable")
            ax.axis("off")
            continue
        finite = f[np.isfinite(f)]
        if symmetric:
            vmax_abs = float(np.max(np.abs(finite))) if finite.size else 1.0
            im = ax.imshow(f, origin="lower", cmap=cmap,
                            vmin=-vmax_abs, vmax=vmax_abs)
        else:
            vmax = float(np.max(finite)) if finite.size else 1.0
            vmin = float(np.min(finite)) if finite.size else 0.0
            if vmin == vmax:
                vmax = vmin + 1e-12
            im = ax.imshow(f, origin="lower", cmap=cmap, vmin=vmin, vmax=vmax)
        ax.set_title(lbl, fontsize=10)
        ax.set_xlabel("x"); ax.set_ylabel("y")
        plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    fig.suptitle(title)
    fig.tight_layout()
    fig.savefig(out_path, dpi=110)
    plt.close(fig)


def _grid_panel(out_path, panels, title, cmap="viridis", ncols=5,
                 symmetric=False):
    n = len(panels)
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(4 * ncols, 4 * nrows))
    axes = axes.flatten() if hasattr(axes, "flatten") else [axes]
    for k in range(len(axes)):
        if k >= n:
            axes[k].axis("off"); continue
        lbl, f = panels[k]
        ax = axes[k]
        if f is None:
            ax.set_title(f"{lbl} - unavailable"); ax.axis("off"); continue
        finite = f[np.isfinite(f)]
        if symmetric:
            vmax = float(np.max(np.abs(finite))) if finite.size else 1.0
            im = ax.imshow(f, origin="lower", cmap=cmap,
                            vmin=-vmax, vmax=vmax)
        else:
            vmax = float(np.max(finite)) if finite.size else 1.0
            vmin = float(np.min(finite)) if finite.size else 0.0
            if vmin == vmax:
                vmax = vmin + 1e-12
            im = ax.imshow(f, origin="lower", cmap=cmap, vmin=vmin, vmax=vmax)
        ax.set_title(lbl, fontsize=9)
        ax.set_xlabel("x"); ax.set_ylabel("y")
        plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    fig.suptitle(title)
    fig.tight_layout()
    fig.savefig(out_path, dpi=110)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Main entry
# ---------------------------------------------------------------------------
def generate_all_plots(OUT, PLOTS, CLUSTERS, cluster_data, cluster_run_data,
                        all_stage_metrics, all_stage_vs_prev_metrics,
                        all_stage_geometric, all_stage_lag,
                        all_stage_longtrans, all_stage_comparison,
                        all_jacobian_verification, all_time_evolution,
                        all_radial, all_multipole, all_power, all_peak,
                        all_wrong_control, all_alpha,
                        first_divergence, loss_rows, STAGE_REGISTRY):
    PLOTS.mkdir(parents=True, exist_ok=True)
    cids = [c["id"] for c in CLUSTERS]

    # ====================================================================
    # 1. pipeline_stage_overview.png
    # ====================================================================
    first_cid = cids[0]
    c10 = cluster_run_data[first_cid]["c10_stages"]
    a8 = cluster_run_data[first_cid]["a8_stages"]
    overview_panels = [
        ("rho proxy (S00)", cluster_data[first_cid]["rho"]),
        ("GR kappa", cluster_run_data[first_cid]["gr_kappa"]),
        ("C10 c (S01)", c10["S01"]),
        ("C10 grad (S02)", c10["S02"]),
        ("A8 c_init (S01)", a8["S01"]),
        ("A8 grad (S02)", a8["S02"]),
        ("C10 R_mag (S13)", c10.get("S13_P0") if c10.get("S13_P0") is not None else np.zeros_like(c10["c"])),
        ("A8 R_mag (S13)", a8.get("S13_P0") if a8.get("S13_P0") is not None else np.zeros_like(a8["c_final"])),
        ("C10 kappa", cluster_run_data[first_cid]["c10_jac"]["convergence"]),
        ("A8 kappa", cluster_run_data[first_cid]["a8_jac"]["convergence"]),
    ]
    _grid_panel(PLOTS / "pipeline_stage_overview.png", overview_panels,
                 f"Pipeline stage overview - {first_cid}", ncols=5,
                 cmap="viridis")

    # ====================================================================
    # 2. stage_correlation_vs_gr.png
    # ====================================================================
    # Aggregate per stage (S0, N0) Pearson vs GR across clusters
    c10_field_ids = ["S00", "S01", "S02", "C10-S03", "C10-S04", "C10-S05",
                     "C10-S06", "S13_P0", "S13_P3", "S13_P4",
                     "S15_Dx", "S15_Dy", "S15_Dmag", "S16_mapping",
                     "S18_trace", "S18_det", "S19_kappa", "S20_gamma1",
                     "S20_gamma2", "S20_gamma_mag"]
    a8_field_ids = ["S00", "S01", "S02", "A8-S03", "A8-S04", "A8-S05",
                    "A8-S06", "A8-S07", "A8-S08", "A8-S09", "A8-S10",
                    "A8-S11", "A8-S12", "A8-P5", "A8-P6",
                    "S13_P0", "S13_P3", "S13_P4",
                    "S15_Dx", "S15_Dy", "S15_Dmag", "S16_mapping",
                    "S18_trace", "S18_det", "S19_kappa", "S20_gamma1",
                    "S20_gamma2", "S20_gamma_mag"]

    def _aggregate_per_stage(metrics, model, field_ids):
        out = []
        for sid in field_ids:
            vals = [m["pearson"] for m in metrics
                    if m["model"] == model and m["stage_id"] == sid
                    and m["smoothing"] == "S0" and m["norm_mode"] == "N0"]
            vals = [v for v in vals if math.isfinite(v)]
            if vals:
                out.append((sid, float(np.median(vals)), float(np.min(vals)),
                             float(np.max(vals))))
            else:
                out.append((sid, float("nan"), float("nan"), float("nan")))
        return out

    c10_agg = _aggregate_per_stage(all_stage_metrics, "C10", c10_field_ids)
    a8_agg = _aggregate_per_stage(all_stage_metrics, "A8", a8_field_ids)

    fig, axes = plt.subplots(2, 1, figsize=(14, 8), sharex=True)
    for ax, agg, label, color in [(axes[0], c10_agg, "C10", "tab:blue"),
                                    (axes[1], a8_agg, "A8/T1", "tab:red")]:
        ys = [a[1] for a in agg]
        ymins = [a[2] for a in agg]
        ymaxs = [a[3] for a in agg]
        x = np.arange(len(agg))
        ax.errorbar(x, ys, yerr=[np.array(ys) - np.array(ymins),
                                    np.array(ymaxs) - np.array(ys)],
                     fmt="o-", color=color, capsize=3)
        ax.set_xticks(x)
        ax.set_xticklabels([a[0] for a in agg], rotation=45, ha="right",
                            fontsize=8)
        ax.set_ylabel(f"Median Pearson(kappa) vs GR\n(S0, N0)")
        ax.set_title(f"{label} - stage correlation vs GR (median across 5 clusters, "
                      "min/max band)")
        ax.axhline(0, color="k", linewidth=0.5)
        ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(PLOTS / "stage_correlation_vs_gr.png", dpi=110)
    plt.close(fig)

    # ====================================================================
    # 3. stage_correlation_vs_previous.png
    # ====================================================================
    fig, axes = plt.subplots(2, 1, figsize=(14, 8), sharex=True)
    for ax, model_name, label, color in [
        (axes[0], "C10", "C10", "tab:blue"),
        (axes[1], "A8", "A8/T1", "tab:red"),
    ]:
        fids = c10_field_ids if model_name == "C10" else a8_field_ids
        fids = [f for f in fids if f in ["S00", "S01", "S02", "C10-S03",
                                            "C10-S04", "C10-S05", "C10-S06"]
                if model_name == "C10"] or [f for f in fids if f in
                ["S00", "S01", "S02", "A8-S03", "A8-S04", "A8-S05",
                 "A8-S06", "A8-S07", "A8-S08", "A8-S09", "A8-S10",
                 "A8-S11", "A8-S12"]]
        vals = []
        for sid in fids:
            rows = [r for r in all_stage_vs_prev_metrics
                    if r["model"] == model_name and r["stage_id"] == sid]
            drs = [r["pearson_vs_previous"] for r in rows
                   if math.isfinite(r["pearson_vs_previous"])]
            vals.append(float(np.median(drs)) if drs else float("nan"))
        x = np.arange(len(fids))
        ax.bar(x, vals, color=color, alpha=0.7)
        ax.set_xticks(x)
        ax.set_xticklabels(fids, rotation=45, ha="right", fontsize=8)
        ax.set_ylabel("Pearson vs previous stage (median, 5 clusters)")
        ax.set_title(f"{label} - stage-to-previous correlation")
        ax.axhline(0, color="k", linewidth=0.5)
        ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(PLOTS / "stage_correlation_vs_previous.png", dpi=110)
    plt.close(fig)

    # ====================================================================
    # 4. stage_amplitude_ratio.png
    # ====================================================================
    fig, axes = plt.subplots(2, 1, figsize=(14, 8), sharex=True)
    for ax, agg, label, color in [(axes[0], c10_agg, "C10", "tab:blue"),
                                    (axes[1], a8_agg, "A8/T1", "tab:red")]:
        amps = []
        for sid, _, _, _ in agg:
            rows = [m for m in all_stage_metrics
                    if m["model"] == label.replace("/T1", "")
                    and m["stage_id"] == sid
                    and m["smoothing"] == "S0" and m["norm_mode"] == "N0"]
            vs = [r["amplitude_ratio"] for r in rows
                  if math.isfinite(r["amplitude_ratio"])]
            amps.append(float(np.median(vs)) if vs else float("nan"))
        ax.bar(np.arange(len(agg)), amps, color=color, alpha=0.7)
        ax.set_xticks(np.arange(len(agg)))
        ax.set_xticklabels([a[0] for a in agg], rotation=45, ha="right",
                            fontsize=8)
        ax.set_ylabel("Median RMS(GR) / RMS(stage)")
        ax.set_title(f"{label} - amplitude ratio (S0, N0)")
        ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(PLOTS / "stage_amplitude_ratio.png", dpi=110)
    plt.close(fig)

    # ====================================================================
    # 5. stage_loss_score.png
    # ====================================================================
    loss_by_model = {}
    for r in loss_rows:
        loss_by_model.setdefault(r["model"], {}).setdefault(
            r["stage_id"], []).append(r["stage_loss_score"])
    fig, axes = plt.subplots(1, 2, figsize=(16, 6), sharey=True)
    for ax, model_name, label, color in [
        (axes[0], "C10", "C10", "tab:blue"),
        (axes[1], "A8", "A8/T1", "tab:red"),
    ]:
        d = loss_by_model.get(model_name, {})
        if d:
            keys = list(d.keys())
            medians = [float(np.median([v for v in d[k] if math.isfinite(v)]))
                        for k in keys]
        else:
            keys, medians = [], []
        ax.bar(np.arange(len(keys)), medians, color=color, alpha=0.7)
        ax.set_xticks(np.arange(len(keys)))
        ax.set_xticklabels(keys, rotation=45, ha="right", fontsize=8)
        ax.set_ylabel("Median L_i = -Δr + 0.25ΔD_NRMS")
        ax.set_title(f"{label} - stage loss score (5 clusters)")
        ax.axhline(0, color="k", linewidth=0.5)
        ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(PLOTS / "stage_loss_score.png", dpi=110)
    plt.close(fig)

    # ====================================================================
    # 6. first_divergence_by_cluster.png
    # ====================================================================
    fig, ax = plt.subplots(figsize=(10, 5))
    fd_c10 = [first_divergence["C10"].get(c, "none") for c in cids]
    fd_a8 = [first_divergence["A8"].get(c, "none") for c in cids]
    x = np.arange(len(cids))
    width = 0.4
    ax.bar(x - width / 2, [hash(f) % 20 for f in fd_c10], width,
            label="C10", color="tab:blue")
    ax.bar(x + width / 2, [hash(f) % 20 for f in fd_a8], width,
            label="A8/T1", color="tab:red")
    for i, c in enumerate(cids):
        ax.text(i - width / 2, hash(fd_c10[i]) % 20 + 0.5, fd_c10[i],
                ha="center", fontsize=7)
        ax.text(i + width / 2, hash(fd_a8[i]) % 20 + 0.5, fd_a8[i],
                ha="center", fontsize=7)
    ax.set_xticks(x)
    ax.set_xticklabels(cids, rotation=30)
    ax.set_ylabel("Stage index (proxy hash)")
    ax.set_title("First material divergence stage per cluster")
    ax.legend()
    fig.tight_layout()
    fig.savefig(PLOTS / "first_divergence_by_cluster.png", dpi=110)
    plt.close(fig)

    # ====================================================================
    # 7. cumulative_divergence.png
    # ====================================================================
    fig, ax = plt.subplots(figsize=(10, 5))
    # Use the correct Δr from loss_rows (r_GR(S_i) - r_GR(S_{i-1}))
    for k, model_name in enumerate(["C10", "A8"]):
        fids = ["S00", "S01", "S02", "C10-S03", "C10-S04", "C10-S05",
                "C10-S06"] if model_name == "C10" else [
            "S00", "S01", "S02", "A8-S03", "A8-S04", "A8-S05", "A8-S06",
            "A8-S07", "A8-S08", "A8-S09", "A8-S10", "A8-S11", "A8-S12"]
        cum = []
        running = 0.0
        for sid in fids[1:]:
            rows = [r for r in loss_rows
                    if r["model"] == model_name and r["stage_id"] == sid]
            drs = [r["delta_r_vs_gr"] for r in rows
                   if math.isfinite(r["delta_r_vs_gr"])]
            if drs:
                running += float(np.median(drs))
            cum.append(running)
        ax.plot(range(len(cum)), cum, marker="o",
                 label=f"{model_name} cumulative Δr",
                 color="tab:blue" if model_name == "C10" else "tab:red")
    ax.set_xticks(range(len(cum)))
    ax.set_xticklabels(fids[1:], rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("Cumulative Δr (median across clusters)")
    ax.axhline(0, color="k", linewidth=0.5)
    ax.set_title("Cumulative divergence trajectory (Δr = r_GR(S_i) - r_GR(S_{i-1}))")
    ax.legend(); ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(PLOTS / "cumulative_divergence.png", dpi=110)
    plt.close(fig)

    # ====================================================================
    # 8. c10_stage_maps.png  (first cluster)
    # ====================================================================
    panels = [
        ("S00 rho", cluster_data[first_cid]["rho"]),
        ("S01 c", c10["S01"]),
        ("S02 |grad|", c10["S02"]),
        ("C10-S03 coh", c10["C10-S03"]),
        ("C10-S04 mem", c10["C10-S04"]),
        ("C10-S05 int", c10["C10-S05"]),
        ("C10-S06 comb", c10["C10-S06"]),
        ("S13_P0 |R|", c10["S13_P0"]),
        ("kappa final", cluster_run_data[first_cid]["c10_jac"]["convergence"]),
    ]
    _grid_panel(PLOTS / "c10_stage_maps.png", panels,
                 f"C10 stage maps - {first_cid}", ncols=3,
                 cmap="viridis")

    # ====================================================================
    # 9. a8_stage_maps.png
    # ====================================================================
    panels = [
        ("S00 rho", cluster_data[first_cid]["rho"]),
        ("S01 c_init", a8["S01"]),
        ("S02 |grad|", a8["S02"]),
        ("A8-S03 fast_pre", a8["A8-S03"]),
        ("A8-S04 fast_post", a8["A8-S04"]),
        ("A8-S05 slow_pre", a8["A8-S05"]),
        ("A8-S06 slow_post", a8["A8-S06"]),
        ("A8-S07 J_FS", a8["A8-S07"]),
        ("A8-S08 J_SF", a8["A8-S08"]),
        ("A8-S09 J_net", a8["A8-S09"]),
        ("A8-S10 mean", a8["A8-S10"]),
        ("A8-S11 memory", a8["A8-S11"]),
        ("A8-S12 nbr", a8["A8-S12"]),
        ("S13_P0 |R|", a8["S13_P0"]),
        ("kappa final", cluster_run_data[first_cid]["a8_jac"]["convergence"]),
    ]
    _grid_panel(PLOTS / "a8_stage_maps.png", panels,
                 f"A8/T1 stage maps - {first_cid}", ncols=4,
                 cmap="viridis")

    # ====================================================================
    # 10. fast_slow_state_evolution.png (first cluster)
    # ====================================================================
    # Use one shot: show RMS vs snapshot index for fast_pre, fast_post, slow_pre, slow_post
    fig, ax = plt.subplots(figsize=(10, 5))
    colors = ["tab:blue", "tab:cyan", "tab:red", "tab:orange"]
    for ck, cid in enumerate(cids):
        rows_fp = [r for r in all_time_evolution if r["cluster_id"] == cid
                    and r["stage_name"] == "A8-fast_pre"]
        rows_fpo = [r for r in all_time_evolution if r["cluster_id"] == cid
                     and r["stage_name"] == "A8-fast_post"]
        rows_sp = [r for r in all_time_evolution if r["cluster_id"] == cid
                    and r["stage_name"] == "A8-slow_pre"]
        rows_spo = [r for r in all_time_evolution if r["cluster_id"] == cid
                     and r["stage_name"] == "A8-slow_post"]
        # We plot only the first cluster to keep it readable
        if ck == 0:
            ax.plot([r["snapshot_index"] for r in rows_fp],
                     [r["rms"] for r in rows_fp], "o-", color="tab:blue",
                     label="fast_pre")
            ax.plot([r["snapshot_index"] for r in rows_fpo],
                     [r["rms"] for r in rows_fpo], "s-", color="tab:cyan",
                     label="fast_post")
            ax.plot([r["snapshot_index"] for r in rows_sp],
                     [r["rms"] for r in rows_sp], "^-", color="tab:red",
                     label="slow_pre")
            ax.plot([r["snapshot_index"] for r in rows_spo],
                     [r["rms"] for r in rows_spo], "d-", color="tab:orange",
                     label="slow_post")
    ax.set_xlabel("Timestep index (21 uniformly spaced snapshots)")
    ax.set_ylabel("RMS")
    ax.set_title(f"Fast/slow state evolution (cluster {first_cid})")
    ax.legend(); ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(PLOTS / "fast_slow_state_evolution.png", dpi=110)
    plt.close(fig)

    # ====================================================================
    # 11. fast_slow_exchange_evolution.png
    # ====================================================================
    fig, ax = plt.subplots(figsize=(10, 5))
    rows_fs = [r for r in all_time_evolution if r["cluster_id"] == first_cid
                and r["stage_name"] == "A8-J_FS"]
    rows_sf = [r for r in all_time_evolution if r["cluster_id"] == first_cid
                and r["stage_name"] == "A8-J_SF"]
    rows_net = [r for r in all_time_evolution if r["cluster_id"] == first_cid
                 and r["stage_name"] == "A8-J_net"]
    if rows_fs:
        ax.plot([r["snapshot_index"] for r in rows_fs],
                 [r["mean"] for r in rows_fs], "o-", label="J_FS")
    if rows_sf:
        ax.plot([r["snapshot_index"] for r in rows_sf],
                 [r["mean"] for r in rows_sf], "s-", label="J_SF")
    if rows_net:
        ax.plot([r["snapshot_index"] for r in rows_net],
                 [r["mean"] for r in rows_net], "^-", label="J_net")
    ax.set_xlabel("Timestep index")
    ax.set_ylabel("Mean exchange")
    ax.set_title(f"Fast/slow exchange evolution (cluster {first_cid})")
    ax.legend(); ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(PLOTS / "fast_slow_exchange_evolution.png", dpi=110)
    plt.close(fig)

    # ====================================================================
    # 12. longitudinal_transverse_breakdown.png
    # ====================================================================
    fig, axes = plt.subplots(1, 2, figsize=(14, 5), sharey=True)
    for ax, model_name, label in [(axes[0], "C10", "C10"),
                                    (axes[1], "A8", "A8/T1")]:
        rows = [r for r in all_stage_longtrans if r["model"] == model_name
                and r["reference_observable"] == "kappa"]
        cids_l = sorted({r["cluster_id"] for r in rows})
        for obs in ("gamma1", "gamma2", "gamma_mag", "kappa"):
            vals = [next((r["r_longitudinal"] for r in rows
                            if r["cluster_id"] == c and r["reference_observable"] == obs),
                          float("nan")) for c in cids_l]
            ax.plot(cids_l, vals, marker="o", label=f"long vs {obs}")
        ax.set_xticklabels(cids_l, rotation=30)
        ax.set_ylabel("Pearson(kappa) - longitudinal vs kappa")
        ax.set_title(f"{label} - longitudinal S13 vs GR kappa")
        ax.axhline(0, color="k", linewidth=0.5)
        ax.legend(fontsize=7); ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(PLOTS / "longitudinal_transverse_breakdown.png", dpi=110)
    plt.close(fig)

    # ====================================================================
    # 13. ray_displacement_evolution.png
    # ====================================================================
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    for ax, model_name, label, phot in [
        (axes[0], "C10", "C10", cluster_run_data[first_cid]["c10_phot"]),
        (axes[1], "A8", "A8/T1", cluster_run_data[first_cid]["a8_phot"]),
    ]:
        # Per-step accumulated displacement magnitude
        acc_dx = phot["x"] - cluster_run_data[first_cid]["x0"]
        acc_dy = phot["y"] - cluster_run_data[first_cid]["y0"]
        mags = np.hypot(acc_dx, acc_dy)
        ax.hist(mags, bins=50, alpha=0.7, label=label)
        ax.set_xlabel("|Δ|")
        ax.set_ylabel("Photon count")
        ax.set_title(f"{label} - final ray displacement magnitudes ({first_cid})")
        ax.legend(); ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(PLOTS / "ray_displacement_evolution.png", dpi=110)
    plt.close(fig)

    # ====================================================================
    # 14. jacobian_component_maps.png
    # ====================================================================
    panels = [
        ("GR kappa", cluster_run_data[first_cid]["gr_kappa"]),
        ("C10 kappa", cluster_run_data[first_cid]["c10_jac"]["convergence"]),
        ("A8 kappa", cluster_run_data[first_cid]["a8_jac"]["convergence"]),
        ("C10 gamma1", cluster_run_data[first_cid]["c10_jac"]["shear_g1"]),
        ("A8 gamma1", cluster_run_data[first_cid]["a8_jac"]["shear_g1"]),
        ("C10 gamma2", cluster_run_data[first_cid]["c10_jac"]["shear_g2"]),
        ("A8 gamma2", cluster_run_data[first_cid]["a8_jac"]["shear_g2"]),
    ]
    _grid_panel(PLOTS / "jacobian_component_maps.png", panels,
                 f"Jacobian components - {first_cid}", ncols=4,
                 cmap="viridis")

    # ====================================================================
    # 15. kappa_extraction_breakdown.png
    # ====================================================================
    fig, ax = plt.subplots(figsize=(10, 5))
    for model_name, color in [("C10", "tab:blue"), ("A8", "tab:red")]:
        medians = []
        for cid in cids:
            rows = [m for m in all_stage_metrics
                    if m["cluster"] == cid and m["model"] == model_name
                    and m["stage_id"] == "S00"
                    and m["smoothing"] == "S0" and m["norm_mode"] == "N0"]
            if rows:
                medians.append(rows[0]["pearson"])
            else:
                medians.append(float("nan"))
        ax.plot(cids, medians, marker="o", color=color, label=f"{model_name} S00")
    for model_name, color in [("C10", "tab:blue"), ("A8", "tab:red")]:
        medians = []
        for cid in cids:
            rows = [m for m in all_stage_metrics
                    if m["cluster"] == cid and m["model"] == model_name
                    and m["stage_id"] == "S19_kappa"
                    and m["smoothing"] == "S0" and m["norm_mode"] == "N0"]
            if rows:
                medians.append(rows[0]["pearson"])
            else:
                medians.append(float("nan"))
        ax.plot(cids, medians, marker="s", color=color, linestyle="--",
                 label=f"{model_name} S19")
    ax.set_xticklabels(cids, rotation=30)
    ax.set_ylabel("Pearson vs GR_kappa")
    ax.set_title("Kappa extraction breakdown: S00 vs S19")
    ax.axhline(0, color="k", linewidth=0.5)
    ax.legend(); ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(PLOTS / "kappa_extraction_breakdown.png", dpi=110)
    plt.close(fig)

    # ====================================================================
    # 16. geometric_transform_dashboard.png
    # ====================================================================
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    for ax, model_name in [(axes[0], "C10"), (axes[1], "A8")]:
        rows = [r for r in all_stage_geometric if r["model"] == model_name
                and r["stage_id"] == "S19"]
        if not rows:
            ax.set_title(f"{model_name} - no geometric audit data"); continue
        labels = ["identity", "G1_sign", "G2_rot90", "G3_rot180",
                   "G4_rot270", "G5_hflip", "G6_vflip", "G7_diag",
                   "G8_anti_diag"]
        keys = ["identity", "G1_sign_reversal", "G2_rotation_90",
                "G3_rotation_180", "G4_rotation_270",
                "G5_horizontal_reflection", "G6_vertical_reflection",
                "G7_main_diagonal_transpose", "G8_anti_diagonal_transpose"]
        medians = []
        for k in keys:
            vals = [r[k] for r in rows if math.isfinite(r.get(k, float("nan")))]
            medians.append(float(np.median(vals)) if vals else float("nan"))
        ax.bar(np.arange(len(labels)), medians)
        ax.set_xticks(np.arange(len(labels)))
        ax.set_xticklabels(labels, rotation=30, ha="right", fontsize=8)
        ax.set_ylabel("Median Pearson vs GR")
        ax.set_title(f"{model_name} - geometric transform audit (S19)")
        ax.axhline(0, color="k", linewidth=0.5)
        ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(PLOTS / "geometric_transform_dashboard.png", dpi=110)
    plt.close(fig)

    # ====================================================================
    # 17. spatial_lag_dashboard.png
    # ====================================================================
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    for ax, model_name in [(axes[0], "C10"), (axes[1], "A8")]:
        rows = [r for r in all_stage_lag if r["model"] == model_name
                and r["stage_id"] == "S19"]
        if not rows:
            ax.set_title(f"{model_name} - no spatial lag data"); continue
        cids_l = sorted({r["cluster_id"] for r in rows})
        zeros = [next((r["zero_lag_correlation"] for r in rows
                        if r["cluster_id"] == c), float("nan")) for c in cids_l]
        bests = [next((r["best_lag_correlation"] for r in rows
                        if r["cluster_id"] == c), float("nan")) for c in cids_l]
        x = np.arange(len(cids_l))
        ax.bar(x - 0.2, zeros, 0.4, label="zero lag")
        ax.bar(x + 0.2, bests, 0.4, label="best fixed lag")
        ax.set_xticks(x)
        ax.set_xticklabels(cids_l, rotation=30)
        ax.set_ylabel("Pearson vs GR")
        ax.set_title(f"{model_name} - spatial lag audit (S19)")
        ax.legend(); ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(PLOTS / "spatial_lag_dashboard.png", dpi=110)
    plt.close(fig)

    # ====================================================================
    # 18. radial_stage_evolution.png
    # ====================================================================
    fig, axes = plt.subplots(2, 1, figsize=(14, 8), sharex=True)
    for ax, model_name, label in [(axes[0], "C10", "C10"),
                                    (axes[1], "A8", "A8/T1")]:
        fids = c10_field_ids if model_name == "C10" else a8_field_ids
        fids = [f for f in ["S00", "S01", "S02", "C10-S03", "C10-S04",
                            "C10-S05", "C10-S06", "S13_P0"]
                if model_name == "C10"] or [f for f in
                ["S00", "S01", "S02", "A8-S03", "A8-S04", "A8-S05",
                 "A8-S06", "A8-S07", "A8-S08", "A8-S09", "A8-S10",
                 "A8-S11", "A8-S12", "S13_P0"]]
        for sid in fids:
            vals = [r["integrated_abs_radial_difference"] for r in all_radial
                    if r["model"] == model_name and r["stage_id"] == sid]
            vals = [v for v in vals if math.isfinite(v)]
            if vals:
                ax.plot([sid], [float(np.median(vals))], "o", label=sid)
        ax.set_ylabel("Median integrated |Δr|")
        ax.set_title(f"{label} - radial stage evolution")
        ax.axhline(0, color="k", linewidth=0.5)
        ax.legend(fontsize=7, ncol=4); ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(PLOTS / "radial_stage_evolution.png", dpi=110)
    plt.close(fig)

    # ====================================================================
    # 19. multipole_stage_evolution.png
    # ====================================================================
    fig, axes = plt.subplots(2, 1, figsize=(14, 8), sharex=True)
    for ax, model_name, label in [(axes[0], "C10", "C10"),
                                    (axes[1], "A8", "A8/T1")]:
        fids = c10_field_ids if model_name == "C10" else a8_field_ids
        fids = [f for f in ["S00", "S01", "S02", "C10-S03", "C10-S04",
                            "C10-S05", "C10-S06"]
                if model_name == "C10"] or [f for f in
                ["S00", "S01", "S02", "A8-S03", "A8-S04", "A8-S05",
                 "A8-S06", "A8-S07", "A8-S08", "A8-S09", "A8-S10",
                 "A8-S11", "A8-S12"]]
        for sid in fids:
            vals = [r["D_Q"] for r in all_multipole
                    if r["model"] == model_name and r["stage_id"] == sid
                    and r["m"] == 1]
            vals = [v for v in vals if math.isfinite(v)]
            if vals:
                ax.plot([sid], [float(np.median(vals))], "o", label=sid)
        ax.set_ylabel("Median D_Q")
        ax.set_title(f"{label} - multipole stage evolution")
        ax.legend(fontsize=7, ncol=4); ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(PLOTS / "multipole_stage_evolution.png", dpi=110)
    plt.close(fig)

    # ====================================================================
    # 20. power_spectrum_stage_evolution.png
    # ====================================================================
    fig, axes = plt.subplots(2, 1, figsize=(14, 8), sharex=True)
    for ax, model_name, label in [(axes[0], "C10", "C10"),
                                    (axes[1], "A8", "A8/T1")]:
        fids = ["S00", "S01", "S02", "C10-S03", "C10-S04", "C10-S05",
                "C10-S06"] if model_name == "C10" else [
            "S00", "S01", "S02", "A8-S03", "A8-S04", "A8-S05", "A8-S06",
            "A8-S07", "A8-S08", "A8-S09", "A8-S10", "A8-S11", "A8-S12"]
        for sid in fids:
            vals = [r["D_P"] for r in all_power
                    if r["model"] == model_name and r["stage_id"] == sid
                    and r["bin_index"] == 0]
            vals = [v for v in vals if math.isfinite(v)]
            if vals:
                ax.plot([sid], [float(np.median(vals))], "o", label=sid)
        ax.set_ylabel("Median D_P")
        ax.set_title(f"{label} - power spectrum stage evolution")
        ax.legend(fontsize=7, ncol=4); ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(PLOTS / "power_spectrum_stage_evolution.png", dpi=110)
    plt.close(fig)

    # ====================================================================
    # 21. wrong_control_dashboard.png
    # ====================================================================
    fig, ax = plt.subplots(figsize=(10, 5))
    wc_names = ["WR1_stage_label_shuffle", "WR2_time_reversal_a8",
                "WR3_cell_shuffle_c10_s13",
                "WR4_component_swap_c10_RxRy",
                "WR5_jacobian_swap_A12_A21",
                "WR6_final_map_substitution_preJacobian"]
    for k, name in enumerate(wc_names):
        vals = [r["pearson_vs_GR_kappa"] for r in all_wrong_control
                if r["wrong_control"] == name]
        vals = [v for v in vals if math.isfinite(v)]
        if vals:
            ax.bar(k, float(np.median(vals)), width=0.7)
        ax.text(k, 0.02, name, ha="center", va="bottom", fontsize=7,
                rotation=20)
    ax.set_xticks([]); ax.set_ylabel("Median Pearson vs GR_kappa")
    ax.set_title("Wrong controls - median across 5 clusters")
    ax.axhline(0, color="k", linewidth=0.5)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(PLOTS / "wrong_control_dashboard.png", dpi=110)
    plt.close(fig)

    # ====================================================================
    # 22. science_dashboard.png
    # ====================================================================
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    # (a) Pearson per cluster per model S19
    ax = axes[0, 0]
    x = np.arange(len(cids)); width = 0.35
    for i, model_name in enumerate(["C10", "A8"]):
        vals = []
        for cid in cids:
            rows = [m for m in all_stage_metrics
                    if m["cluster"] == cid and m["model"] == model_name
                    and m["stage_id"] == "S19_kappa"
                    and m["smoothing"] == "S0" and m["norm_mode"] == "N0"]
            vals.append(rows[0]["pearson"] if rows else float("nan"))
        ax.bar(x + i * width, vals, width, label=model_name)
    ax.set_xticks(x + width / 2); ax.set_xticklabels(cids, rotation=30)
    ax.axhline(0, color="k", linewidth=0.5)
    ax.set_ylabel("Pearson(kappa) vs GR S19")
    ax.set_title("Final kappa vs GR")
    ax.legend()
    # (b) NRMSE per cluster
    ax = axes[0, 1]
    for i, model_name in enumerate(["C10", "A8"]):
        vals = []
        for cid in cids:
            rows = [m for m in all_stage_metrics
                    if m["cluster"] == cid and m["model"] == model_name
                    and m["stage_id"] == "S19_kappa"
                    and m["smoothing"] == "S0" and m["norm_mode"] == "N0"]
            vals.append(rows[0]["normalized_rms_difference"] if rows else float("nan"))
        ax.bar(x + i * width, vals, width, label=model_name)
    ax.set_xticks(x + width / 2); ax.set_xticklabels(cids, rotation=30)
    ax.set_ylabel("D_NRMS")
    ax.set_title("Final NRMSE")
    ax.legend()
    # (c) D_Q / D_P
    ax = axes[1, 0]
    for i, model_name in enumerate(["C10", "A8"]):
        dqs = []
        dps = []
        for cid in cids:
            rows_q = [r for r in all_multipole if r["cluster_id"] == cid
                       and r["model"] == model_name and r["stage_id"] == "S19"
                       and r["m"] == 1]
            dqs.append(float(np.median([r["D_Q"] for r in rows_q]))
                        if rows_q else float("nan"))
            rows_p = [r for r in all_power if r["cluster_id"] == cid
                       and r["model"] == model_name and r["stage_id"] == "S19"
                       and r["bin_index"] == 0]
            dps.append(float(np.median([r["D_P"] for r in rows_p]))
                        if rows_p else float("nan"))
        ax.bar(x + i * width - 0.2, dqs, 0.4, label=f"{model_name} D_Q")
        ax.bar(x + i * width + 0.2, dps, 0.4, label=f"{model_name} D_P")
    ax.set_xticks(x); ax.set_xticklabels(cids, rotation=30)
    ax.set_ylabel("D_Q / D_P")
    ax.set_title("Multipole / power spectrum distances")
    ax.legend(fontsize=8)
    # (d) First divergence per model
    ax = axes[1, 1]
    codes = {sid: k for k, sid in enumerate(
        ["S00", "S01", "S02", "C10-S03", "C10-S04", "C10-S05", "C10-S06",
         "A8-S03", "A8-S04", "A8-S05", "A8-S06", "A8-S07", "A8-S08",
         "A8-S09", "A8-S10", "A8-S11", "A8-S12", "S19"])}
    for i, model_name in enumerate(["C10", "A8"]):
        ys = [codes.get(first_divergence[model_name].get(c, "none"), 99)
              for c in cids]
        ax.bar(x + i * width, ys, width, label=model_name)
    ax.set_xticks(x + width / 2); ax.set_xticklabels(cids, rotation=30)
    ax.set_ylabel("First-divergence stage index")
    ax.set_title("First material divergence per cluster")
    ax.legend()
    fig.suptitle("Science dashboard - macro-micro response bridge diagnostic")
    fig.tight_layout()
    fig.savefig(PLOTS / "science_dashboard.png", dpi=110)
    plt.close(fig)

    # ====================================================================
    # 23. Per-cluster dashboards
    # ====================================================================
    for cluster in CLUSTERS:
        cid = cluster["id"]
        rd = cluster_run_data[cid]
        panels = [
            ("rho proxy", rd["rho"]),
            ("GR kappa", rd["gr_kappa"]),
            ("C10 c", rd["c10_stages"]["S01"]),
            ("C10 neigh (P0)", rd["c10_stages"]["C10-S03"]),
            ("A8 fast_pre", rd["a8_stages"]["A8-S03"]),
            ("A8 slow_pre", rd["a8_stages"]["A8-S05"]),
            ("A8 combined", rd["a8_stages"]["A8-S10"]),
            ("C10 R_mag", rd["c10_stages"]["S13_P0"]),
            ("C10 R_div", rd["c10_stages"]["S13_P1"]),
            ("C10 R_curl", rd["c10_stages"]["S13_P2"]),
            ("A8 R_mag", rd["a8_stages"]["S13_P0"]),
            ("|Δ| accumulated", np.hypot(rf := rd["c10_phot"]["x"] - rd["x0"],
                                          rd["c10_phot"]["y"] - rd["y0"]) if False else
             np.full((bins := 64, bins), float("nan"))),
            ("tr(A)", np.full((64, 64), float("nan"))),
            ("C10 kappa", rd["c10_jac"]["convergence"]),
            ("A8 kappa", rd["a8_jac"]["convergence"]),
            ("|C10 - GR|", rd["c10_jac"]["convergence"] - rd["gr_kappa"]),
            ("|A8 - GR|", rd["a8_jac"]["convergence"] - rd["gr_kappa"]),
        ]
        # Compute simple aggregated fields for the displacement-related panels
        x0, y0 = rd["x0"], rd["y0"]
        extent = 8.0; bins = 64
        for k, phot in [(11, rd["c10_phot"])]:
            d_mag = np.full((bins, bins), float("nan"))
            x_edges = np.linspace(-extent, extent, bins + 1)
            for i in range(bins):
                for j in range(bins):
                    in_bin = ((x0 >= x_edges[j]) & (x0 < x_edges[j + 1])
                              & (y0 >= x_edges[i]) & (y0 < x_edges[i + 1]))
                    if in_bin.sum() > 0:
                        d_mag[i, j] = float(np.mean(
                            np.hypot(phot["x"][in_bin] - x0[in_bin],
                                      phot["y"][in_bin] - y0[in_bin])))
            panels[k] = ("|Δ| accumulated", d_mag)
        # Compute |Δ| for A8 too
        d_a8 = np.full((bins, bins), float("nan"))
        for i in range(bins):
            for j in range(bins):
                in_bin = ((x0 >= x_edges[j]) & (x0 < x_edges[j + 1])
                          & (y0 >= x_edges[i]) & (y0 < x_edges[i + 1]))
                if in_bin.sum() > 0:
                    d_a8[i, j] = float(np.mean(
                        np.hypot(rd["a8_phot"]["x"][in_bin] - x0[in_bin],
                                  rd["a8_phot"]["y"][in_bin] - y0[in_bin])))
        # tr(A) for A8
        tr_a8 = np.full((bins, bins), float("nan"))
        for i in range(bins):
            for j in range(bins):
                in_bin = ((x0 >= x_edges[j]) & (x0 < x_edges[j + 1])
                          & (y0 >= x_edges[i]) & (y0 < x_edges[i + 1]))
                if in_bin.sum() >= 6:
                    x0c = x0[in_bin] - x0[in_bin].mean()
                    y0c = y0[in_bin] - y0[in_bin].mean()
                    Jx = np.linalg.lstsq(
                        np.column_stack([x0c, y0c]),
                        rd["a8_phot"]["x"][in_bin] - rd["a8_phot"]["x"][in_bin].mean(),
                        rcond=None)[0]
                    Jy = np.linalg.lstsq(
                        np.column_stack([x0c, y0c]),
                        rd["a8_phot"]["y"][in_bin] - rd["a8_phot"]["y"][in_bin].mean(),
                        rcond=None)[0]
                    tr_a8[i, j] = Jx[0] + Jy[1]
        _grid_panel(PLOTS / f"stage_dashboard_{cid}.png", panels,
                      f"Stage dashboard - {cluster['label']}", ncols=5,
                      cmap="viridis")
