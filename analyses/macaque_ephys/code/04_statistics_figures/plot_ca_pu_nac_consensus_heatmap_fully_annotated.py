from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, ListedColormap
from matplotlib.patches import Rectangle
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "outputs/R3_panels_v3"
ZFILE = ROOT / "outputs/dSTR_dSTRvSTR_E_QC/Ca_Pu_NAC_final19_consensus_scan/Ca_Pu_NAC_final19_transformed_Zscore_matrix_n390.csv"
CELL_ORDER = BASE / "E_consensus368_final19_diagonal_heatmap_cell_order.csv"
FEATURE_ORDER = BASE / "E_consensus368_final19_diagonal_feature_order_and_class_means.csv"
OUT = ROOT / "outputs/R3_heatmap_annotated"

E_ORDER = ["C1", "C2", "C3", "C4"]
E_COLORS = {"C1": "#F8766D", "C2": "#7CAE00", "C3": "#00BFC4", "C4": "#C77CFF"}
T_ORDER = ["D1", "D2", "Hybrid"]
T_COLORS = {"D1": "#D95F02", "D2": "#008F7A", "Hybrid": "#6B7280"}
ROI_ORDER = ["Ca", "Pu", "NAC"]
ROI_COLORS = {"Ca": "#246BB2", "Pu": "#3FA0E8", "NAC": "#7B1FA2"}
HEAT = LinearSegmentedColormap.from_list("magenta_black_yellow", ["#C51B7D", "#000000", "#FFD92F"])

# Manuscript-wide electrophysiology abbreviations.  These labels intentionally
# match the Macaque T-E RRR figure exactly.
RRR_LABELS = {
    "Epsy_width_adapt_ratio": "Width adapt",
    "Epsy_width_rheo": "Width rheo",
    "Epsy_ahp_delay_5spike": "AHP delay",
    "Epsy_ahp_delay_ratio_5spike": "AHP ratio",
    "Epsy_upstroke_downstroke_ratio_rheo": "Up/down ratio",
    "Epsy_peak_deltav_rheo": "Peak ΔV",
    "Epsy_peak_v_rheo": "Peak V",
    "Epsy_downstroke_adapt_ratio": "Downstroke adapt",
    "Epsy_upstroke_adapt_ratio": "Upstroke adapt",
    "Epsy_threshold_v_adapt_ratio": "Threshold adapt",
    "Epsy_postap_slope_rheo": "Post-AP slope",
    "Epsy_postap_slope_hero": "Post-AP hero",
    "Epsy_threshold_v_rheo": "Threshold V",
    "Epsy_peak_v_adapt_ratio": "Peak V adapt",
    "Epsy_trough_t_hero": "Trough t hero",
    "Epsy_trough_v_rheo": "Trough V",
    "Epsy_trough_t_rheo": "Trough t",
    "Epsy_fast_trough_v_rheo": "Fast trough",
    "Epsy_threshold_v_short_square": "Threshold SS",
}


def swatch_legend(fig, x, y, title, order, colors, dy=.055, first_offset=.055):
    fig.text(x, y, title, fontsize=4, weight="bold", ha="left", va="top", color="#222222")
    for i, item in enumerate(order):
        yy = y - first_offset - i * dy
        fig.add_artist(Rectangle((x, yy - .012), .012, .024, transform=fig.transFigure,
                                 facecolor=colors[item], edgecolor="none"))
        fig.text(x + .017, yy, item, fontsize=4, ha="left", va="center", color="#222222")


def main():
    mpl.rcParams.update({"font.family": "Arial", "font.size": 4,
                         "pdf.fonttype": 42, "ps.fonttype": 42})
    OUT.mkdir(parents=True, exist_ok=True)
    z = pd.read_csv(ZFILE, dtype={"cell_label": str})
    cells = pd.read_csv(CELL_ORDER, dtype={"cell_label": str})
    features = pd.read_csv(FEATURE_ORDER)["feature"].tolist()
    data = cells.merge(z, on="cell_label", validate="one_to_one")

    matrix_parts, e_parts, t_parts, roi_parts = [], [], [], []
    for i, group in enumerate(E_ORDER):
        block = data.loc[data.HC_class.eq(group)]
        matrix_parts.append(block[features].to_numpy(float).T)
        e_parts.append(np.full((1, len(block)), i, float))
        t_parts.append(np.array([[T_ORDER.index(x) for x in block.T_class]], float))
        roi_parts.append(np.array([[ROI_ORDER.index(x) for x in block.Lib_region_of_interest_label]], float))
        if i < 3:
            matrix_parts.append(np.full((len(features), 1), np.nan))
            e_parts.append(np.full((1, 1), np.nan)); t_parts.append(np.full((1, 1), np.nan)); roi_parts.append(np.full((1, 1), np.nan))
    matrix = np.concatenate(matrix_parts, axis=1)
    bars = [np.concatenate(parts, axis=1) for parts in [e_parts, t_parts, roi_parts]]

    fig = plt.figure(figsize=(8.2, 1.5), facecolor="white")
    gs = fig.add_gridspec(4, 1, height_ratios=[.045, .045, .045, 1],
                          left=.055, right=.620, bottom=.055, top=.975, hspace=.055)
    bar_specs = [
        (bars[0], ListedColormap([E_COLORS[x] for x in E_ORDER]), "E class"),
        (bars[1], ListedColormap([T_COLORS[x] for x in T_ORDER]), "T class"),
        (bars[2], ListedColormap([ROI_COLORS[x] for x in ROI_ORDER]), "ROI"),
    ]
    for i, (bar, cmap, label) in enumerate(bar_specs):
        cmap.set_bad("white")
        ax = fig.add_subplot(gs[i]); ax.imshow(bar, aspect="auto", interpolation="nearest", cmap=cmap)
        ax.axis("off"); ax.text(-.008, .5, label, transform=ax.transAxes, ha="right", va="center", fontsize=4)
    heat = HEAT.copy(); heat.set_bad("white")
    ax = fig.add_subplot(gs[3]); ax.imshow(matrix, aspect="auto", interpolation="nearest", cmap=heat, vmin=-2, vmax=2)
    ax.set_xticks([]); ax.set_yticks(np.arange(len(features))); ax.set_yticklabels(
        [RRR_LABELS.get(x, x.removeprefix("Epsy_").replace("_", " ")) for x in features], fontsize=4)
    ax.yaxis.tick_right(); ax.tick_params(axis="y", length=0, pad=1.3)
    boundaries = np.cumsum([2, 5, 6])
    for boundary in boundaries:
        ax.axhline(boundary - .5, color="white", lw=.35, zorder=5)
    for spine in ax.spines.values(): spine.set_visible(False)

    swatch_legend(fig, .830, .955, "E class", E_ORDER, E_COLORS, .050)
    swatch_legend(fig, .900, .955, "T class", T_ORDER, T_COLORS, .050)
    swatch_legend(fig, .960, .955, "ROI", ROI_ORDER, ROI_COLORS, .050)
    cax = fig.add_axes([.830, .105, .145, .045])
    gradient = np.linspace(-2, 2, 256)[None, :]
    cax.imshow(gradient, aspect="auto", cmap=HEAT, vmin=-2, vmax=2)
    cax.set_yticks([]); cax.set_xticks([0, 127.5, 255], ["−2", "0", "2"], fontsize=4)
    cax.tick_params(axis="x", length=1.3, width=.4, pad=1)
    cax.set_xlabel("Z-score", fontsize=4, labelpad=.5)
    for spine in cax.spines.values(): spine.set_linewidth(.35)

    stem = OUT / "E_consensus368_final19_heatmap_fully_annotated"
    fig.savefig(stem.with_suffix(".png"), dpi=900, facecolor="white")
    fig.savefig(stem.with_suffix(".pdf"), facecolor="white")
    plt.close(fig)
    print(stem.with_suffix(".png"))

    legend_fig = plt.figure(figsize=(1.8, 1.0), facecolor="white")
    swatch_legend(legend_fig, .055, .94, "E class", E_ORDER, E_COLORS, .105, .105)
    swatch_legend(legend_fig, .385, .94, "T class", T_ORDER, T_COLORS, .105, .105)
    swatch_legend(legend_fig, .700, .94, "ROI", ROI_ORDER, ROI_COLORS, .105, .105)
    legend_cax = legend_fig.add_axes([.07, .18, .84, .075])
    legend_cax.imshow(gradient, aspect="auto", cmap=HEAT, vmin=-2, vmax=2)
    legend_cax.set_yticks([]); legend_cax.set_xticks([0, 127.5, 255], ["−2", "0", "2"], fontsize=4)
    legend_cax.tick_params(axis="x", length=1.3, width=.4, pad=1)
    legend_cax.set_xlabel("Z-score", fontsize=4, labelpad=.5)
    for spine in legend_cax.spines.values(): spine.set_linewidth(.35)
    legend_stem = OUT / "E_consensus_heatmap_annotation_legend_only_W1p8_H1"
    legend_fig.savefig(legend_stem.with_suffix(".png"), dpi=1200, facecolor="white")
    legend_fig.savefig(legend_stem.with_suffix(".pdf"), facecolor="white")
    plt.close(legend_fig)
    print(legend_stem.with_suffix(".png"))


if __name__ == "__main__":
    main()
