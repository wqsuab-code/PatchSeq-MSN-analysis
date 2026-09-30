from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(r"C:\Users\53461\OneDrive\Documentos\Patch-seq_Mouse_Acb_MSN_T-type_Visualization")
SRC = ROOT / "outputs/morph_qc/morph187_NPC3-4_HCK4-5_GC_allres_raw"
OUT = ROOT / "outputs/morph_qc/final_morph187_NPC3_HCK4_GCres0p50_figures/30_PCA_PC1-10_explained_and_PC1-3_loadings"

DISPLAY = {
    "M_soma_circularity_index": "Soma circularity",
    "M_soma_aspect_ratio": "Soma aspect ratio",
    "M_cell_max_radial_dist": "Max radial distance",
    "M_total_number_of_neurites": "Primary neurite number",
    "M_basal_dendrite_avg_tortuosity": "Mean tortuosity",
    "M_Total_neurite_length_(sections)": "Total neurite length",
    "M_Number_of_bifurcation_points": "Bifurcation points",
    "M_Maximum_branch_order": "Maximum branch order",
    "M_trunk_angle_min": "Minimum trunk angle",
    "M_trunk_angle_max": "Maximum trunk angle",
}


def configure() -> None:
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "DejaVu Sans"],
        "font.size": 7,
        "axes.titlesize": 7,
        "axes.labelsize": 7,
        "xtick.labelsize": 6.5,
        "ytick.labelsize": 6.5,
        "axes.linewidth": .75,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
    })


def save(fig: plt.Figure, stem: str) -> None:
    fig.savefig(OUT / f"{stem}.png", dpi=900, bbox_inches="tight", pad_inches=.03, facecolor="white")
    fig.savefig(OUT / f"{stem}.pdf", bbox_inches="tight", pad_inches=.03, facecolor="white")
    fig.savefig(OUT / f"{stem}.svg", bbox_inches="tight", pad_inches=.03, facecolor="white")
    plt.close(fig)


def plot_explained(variance: pd.DataFrame) -> None:
    pc = variance["PC"].astype(int).to_numpy()
    explained = variance["Explained_variance_percent"].astype(float).to_numpy()
    cumulative = variance["Cumulative_variance_percent"].astype(float).to_numpy()

    fig, ax = plt.subplots(figsize=(3.5, 2.45))
    ax.bar(pc, explained, color="#9E9E9E", width=.72, edgecolor="none", zorder=2)
    ax.plot(pc, cumulative, color="#ED0000", marker="o", ms=2.8, lw=.8, zorder=3)
    for x, value in zip(pc, explained):
        ax.text(x, value + 1.4, f"{value:.1f}", ha="center", va="bottom", fontsize=5.5)
    ax.set_xlim(.35, 10.65)
    ax.set_ylim(0, 105)
    ax.set_xticks(pc, [f"PC{i}" for i in pc], rotation=45, ha="right")
    ax.set_yticks(np.arange(0, 101, 20))
    ax.set_xlabel("Principal component")
    ax.set_ylabel("Explained variance (%)")
    ax.grid(axis="y", color="#E2E2E2", lw=.4)
    ax.set_axisbelow(True)
    ax.spines[["top", "right"]].set_visible(False)
    ax.text(.98, .42, "Bars: individual", transform=ax.transAxes, ha="right", color="#666666", fontsize=6)
    ax.text(.98, .34, "Red line: cumulative", transform=ax.transAxes, ha="right", color="#ED0000", fontsize=6)
    save(fig, "01_PC1-10_individual_and_cumulative_explained_variance")


def plot_loadings(loadings: pd.DataFrame) -> pd.DataFrame:
    work = loadings[["Feature", "PC1", "PC2", "PC3"]].copy()
    work["Display_name"] = work["Feature"].map(DISPLAY)
    feature_order = list(DISPLAY)
    work["Feature_order"] = work["Feature"].map({feature: i for i, feature in enumerate(feature_order)})
    work = work.sort_values("Feature_order")
    y = np.arange(len(work))

    fig, axes = plt.subplots(1, 3, figsize=(7.2, 3.15), sharey=True)
    for axis_index, (ax, pc_name) in enumerate(zip(axes, ["PC1", "PC2", "PC3"])):
        values = work[pc_name].astype(float).to_numpy()
        colors = np.where(values >= 0, "#ED0000", "#00468B")
        ax.barh(y, values, color=colors, height=.62, edgecolor="none")
        ax.axvline(0, color="#222222", lw=.65)
        ax.set_xlim(-.62, .62)
        ax.set_xticks([-.6, -.3, 0, .3, .6])
        ax.set_xlabel("PCA loading")
        ax.set_title(pc_name, pad=3)
        ax.grid(axis="x", color="#E2E2E2", lw=.4)
        ax.set_axisbelow(True)
        ax.spines[["top", "right", "left"]].set_visible(False)
        ax.tick_params(axis="y", length=0)
        if axis_index == 0:
            ax.set_yticks(y, work["Display_name"], fontsize=6.5)
        ax.invert_yaxis()
    fig.subplots_adjust(left=.19, right=.995, bottom=.15, top=.93, wspace=.24)
    save(fig, "02_PC1-PC3_signed_feature_loadings")

    rows = []
    for pc_name in ["PC1", "PC2", "PC3"]:
        ranked = work.assign(
            PC=pc_name,
            Loading=work[pc_name].astype(float),
            Absolute_loading=work[pc_name].astype(float).abs(),
        ).sort_values("Absolute_loading", ascending=False)
        ranked["Rank_by_absolute_loading"] = np.arange(1, len(ranked) + 1)
        rows.append(ranked[[
            "PC", "Rank_by_absolute_loading", "Feature", "Display_name",
            "Loading", "Absolute_loading",
        ]])
    return pd.concat(rows, ignore_index=True)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    configure()
    variance = pd.read_csv(SRC / "02_PCA_variance.csv")
    loadings = pd.read_csv(SRC / "03_PCA_loadings.csv")
    if len(variance) != 10 or set(DISPLAY) != set(loadings["Feature"]):
        raise RuntimeError("The frozen final-187 PCA inputs do not match the expected ten-feature model")
    plot_explained(variance)
    ranked = plot_loadings(loadings)
    variance.to_csv(OUT / "01a_PC1-10_explained_variance_values.csv", index=False)
    ranked.to_csv(OUT / "02a_PC1-PC3_loading_ranking.csv", index=False)


if __name__ == "__main__":
    main()
