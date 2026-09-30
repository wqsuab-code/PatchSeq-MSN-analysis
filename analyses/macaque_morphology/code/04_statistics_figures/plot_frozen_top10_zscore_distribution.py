#!/usr/bin/env python3
"""Distribution audit for the ten frozen Macaque-M radar Z-score features."""
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde


BASE = Path(__file__).resolve().parent
RUN = BASE / "m18_tempfreeze_NPC5_HCK4_res2.3"
ASSIGN = RUN / "01_temp_frozen_assignments_126.csv"
ZSCORE = BASE / "m18_adaptive_pca126" / "02_transformed_z_117.csv"
SELECTED = RUN / "panels_A_I" / "radar_group_top10_feature_selection_117_consensus.csv"
LABELS = BASE / "morphology_feature_abbreviations.csv"
OUT = RUN / "panels_A_I" / "Zscore_top10_radar"


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    features = pd.read_csv(SELECTED).sort_values("axis")["feature"].tolist()
    label_map = dict(pd.read_csv(LABELS).values)
    labels = [label_map[f] for f in features]
    assign = pd.read_csv(ASSIGN, dtype={"cell_label": str})
    keep = assign.loc[assign["concordant"].eq(True), "cell_label"]
    z = pd.read_csv(ZSCORE, dtype={"cell_label": str}).set_index("cell_label").loc[keep, features]
    if z.shape != (117, 10) or z.isna().any().any():
        raise RuntimeError(f"Unexpected frozen matrix: {z.shape}, missing={int(z.isna().sum().sum())}")

    pooled = z.to_numpy(float).ravel()
    q = np.quantile(pooled, [0, .005, .01, .025, .25, .5, .75, .975, .99, .995, 1])
    outside = {
        "n_total": pooled.size,
        "n_abs_gt_2": int(np.sum(np.abs(pooled) > 2)),
        "n_abs_gt_2p5": int(np.sum(np.abs(pooled) > 2.5)),
        "n_abs_gt_3": int(np.sum(np.abs(pooled) > 3)),
    }

    mpl.rcParams.update({
        "font.family": "Arial", "font.size": 8, "axes.titlesize": 9,
        "axes.labelsize": 8, "xtick.labelsize": 7, "ytick.labelsize": 7,
        "pdf.fonttype": 42, "ps.fonttype": 42,
    })
    fig = plt.figure(figsize=(7.2, 4.2), facecolor="white")
    gs = fig.add_gridspec(1, 2, width_ratios=[1.65, 1], wspace=.32,
                          left=.075, right=.98, top=.90, bottom=.13)
    ax = fig.add_subplot(gs[0, 0])
    bins = np.linspace(np.floor(pooled.min() * 4) / 4, np.ceil(pooled.max() * 4) / 4, 42)
    ax.hist(pooled, bins=bins, density=True, color="#6BAED6", alpha=.72,
            edgecolor="white", linewidth=.45)
    grid = np.linspace(bins[0], bins[-1], 600)
    ax.plot(grid, gaussian_kde(pooled)(grid), color="#1F4E79", lw=1.4)
    for x in (-3, 3):
        ax.axvline(x, color="#C43C39", lw=1.0, ls="--")
    ax.axvline(0, color="#444444", lw=.75, ls=":")
    ax.set_title("Pooled distribution of 10 frozen morphology Z-scores")
    ax.set_xlabel("Frozen cohort-wide Z-score")
    ax.set_ylabel("Density")
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", color="#D9D9D9", lw=.5, alpha=.55)
    ax.text(.02, .97,
            f"117 cells × 10 features = {pooled.size:,} values\n"
            f"median = {np.median(pooled):.2f}; 2.5–97.5% = {q[3]:.2f} to {q[7]:.2f}\n"
            f"|Z| > 3: {outside['n_abs_gt_3']} ({100*outside['n_abs_gt_3']/pooled.size:.2f}%)",
            transform=ax.transAxes, ha="left", va="top", fontsize=7,
            bbox={"facecolor": "white", "edgecolor": "#C8C8C8", "linewidth": .5, "alpha": .92})

    bx = fig.add_subplot(gs[0, 1])
    values = [z[f].to_numpy(float) for f in features]
    parts = bx.violinplot(values, positions=np.arange(1, 11), vert=False,
                          showmeans=False, showmedians=False, showextrema=False, widths=.72)
    for body in parts["bodies"]:
        body.set_facecolor("#B9D7EA")
        body.set_edgecolor("#4F81A8")
        body.set_linewidth(.55)
        body.set_alpha(.72)
    bx.boxplot(values, positions=np.arange(1, 11), vert=False, widths=.24,
               showfliers=False, patch_artist=True,
               boxprops={"facecolor": "white", "edgecolor": "#333333", "linewidth": .55},
               medianprops={"color": "#202020", "linewidth": .85},
               whiskerprops={"color": "#555555", "linewidth": .5},
               capprops={"color": "#555555", "linewidth": .5})
    for x in (-3, 3):
        bx.axvline(x, color="#C43C39", lw=.8, ls="--")
    bx.axvline(0, color="#444444", lw=.65, ls=":")
    bx.set_yticks(np.arange(1, 11), labels)
    bx.invert_yaxis()
    bx.set_xlabel("Frozen Z-score")
    bx.set_title("Feature-wise distributions")
    bx.spines[["top", "right", "left"]].set_visible(False)
    bx.tick_params(axis="y", length=0)
    bx.grid(axis="x", color="#D9D9D9", lw=.5, alpha=.55)

    fig.suptitle("Macaque M1–M4 radar input: frozen Z-score distribution", fontsize=10, y=.975)
    stem = OUT / "Macaque_M_top10_frozen_Zscore_distribution_n117"
    fig.savefig(stem.with_suffix(".png"), dpi=600, facecolor="white")
    fig.savefig(stem.with_suffix(".pdf"), facecolor="white")
    plt.close(fig)

    summary = pd.DataFrame({
        "statistic": ["min", "q0.5", "q1", "q2.5", "q25", "median", "q75", "q97.5", "q99", "q99.5", "max",
                      "n_abs_gt_2", "n_abs_gt_2.5", "n_abs_gt_3"],
        "value": list(q) + [outside["n_abs_gt_2"], outside["n_abs_gt_2p5"], outside["n_abs_gt_3"]],
    })
    summary.to_csv(OUT / "Macaque_M_top10_frozen_Zscore_distribution_summary.csv", index=False)
    print(stem.with_suffix(".png"))
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
