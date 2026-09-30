from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, ListedColormap
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs/R3_panels_v3"
ASSIGN = OUT / "00_tuned_p80_ee12_random_seed777_coordinates.csv"
ZFILE = ROOT / "outputs/dSTR_dSTRvSTR_E_QC/Ca_Pu_NAC_final19_consensus_scan/Ca_Pu_NAC_final19_transformed_Zscore_matrix_n390.csv"
SCORES = ROOT / "outputs/dSTR_dSTRvSTR_E_QC/Ca_Pu_NAC_final19_consensus_scan/Ca_Pu_NAC_final19_PCA_scores_n390.csv"
ORDER = ["C1", "C2", "C3", "C4"]
CLASS_COLORS = {"C1": "#F8766D", "C2": "#7CAE00", "C3": "#00BFC4", "C4": "#C77CFF"}
T_ORDER = ["D1", "D2", "Hybrid"]
T_COLORS = {"D1": "#E66101", "D2": "#008F7A", "Hybrid": "#6B7280"}
ROI_ORDER = ["Ca", "Pu", "NAC"]
ROI_COLORS = {"Ca": "#246BB2", "Pu": "#3FA0E8", "NAC": "#7B1FA2"}
HEAT = LinearSegmentedColormap.from_list("magenta_black_yellow", ["#C51B7D", "#000000", "#FFD92F"])


def main():
    mpl.rcParams.update({"font.family": "Arial", "font.size": 4, "pdf.fonttype": 42, "ps.fonttype": 42})
    assign = pd.read_csv(ASSIGN, dtype={"cell_label": str})
    z = pd.read_csv(ZFILE, dtype={"cell_label": str})
    scores = pd.read_csv(SCORES, dtype={"cell_label": str}, usecols=["cell_label", "Lib_region_of_interest_label"])
    features = [x for x in z.columns if x != "cell_label"]
    data = assign.loc[assign["Consensus"].eq(True), ["cell_label", "HC_class", "T_class"]].merge(
        scores, on="cell_label", how="inner", validate="one_to_one"
    ).merge(z, on="cell_label", how="inner", validate="one_to_one")
    if len(data) != 368 or data[features].isna().any().any():
        raise RuntimeError(f"Expected 368 complete consensus cells; got {len(data)}")

    # Arrange rows by the class with the highest mean Z-score. Within each
    # class block, put the most class-selective feature first. Values and signs
    # are not altered; only the row order changes.
    class_means = data.groupby("HC_class", observed=True)[features].mean().T[ORDER]
    row_info = pd.DataFrame(index=features)
    row_info["peak_class"] = class_means.idxmax(axis=1)
    row_info["specificity_range"] = class_means.max(axis=1) - class_means.min(axis=1)
    for group in ORDER:
        row_info[f"mean_{group}"] = class_means[group]
    row_info["peak_order"] = row_info["peak_class"].map({x: i for i, x in enumerate(ORDER)})
    row_info = row_info.sort_values(["peak_order", "specificity_range"], ascending=[True, False])
    ordered_features = row_info.index.tolist()
    row_info.drop(columns="peak_order").reset_index(names="feature").to_csv(
        OUT / "E_consensus368_final19_diagonal_feature_order_and_class_means.csv", index=False
    )

    rng = np.random.default_rng(777)
    matrix_parts, class_parts, t_parts, roi_parts, order_rows = [], [], [], [], []
    for i, group in enumerate(ORDER):
        block = data.loc[data["HC_class"].eq(group)].copy()
        block = block.iloc[rng.permutation(len(block))]
        matrix_parts.append(block[ordered_features].to_numpy(float).T)
        class_parts.append(np.full((1, len(block)), i, float))
        t_parts.append(np.array([[T_ORDER.index(x) for x in block["T_class"]]], float))
        roi_parts.append(np.array([[ROI_ORDER.index(x) for x in block["Lib_region_of_interest_label"]]], float))
        order_rows.append(block[["cell_label", "HC_class", "T_class", "Lib_region_of_interest_label"]])
        if i < 3:
            matrix_parts.append(np.full((len(ordered_features), 1), np.nan))
            class_parts.append(np.full((1, 1), np.nan))
            t_parts.append(np.full((1, 1), np.nan))
            roi_parts.append(np.full((1, 1), np.nan))
    matrix = np.concatenate(matrix_parts, axis=1)
    class_bar = np.concatenate(class_parts, axis=1)
    t_bar = np.concatenate(t_parts, axis=1)
    roi_bar = np.concatenate(roi_parts, axis=1)
    pd.concat(order_rows, ignore_index=True).to_csv(
        OUT / "E_consensus368_final19_diagonal_heatmap_cell_order.csv", index=False
    )

    def render(display_matrix, stem_name, horizontal_boundaries=None, include_roi=False,
               row_labels=False):
        figure_size = (6.2, 1.35) if row_labels else (4.64, 1.0)
        heatmap_right = 0.755 if row_labels else 0.996
        fig = plt.figure(figsize=figure_size, facecolor="white")
        bars = [(class_bar, ListedColormap([CLASS_COLORS[x] for x in ORDER])),
                (t_bar, ListedColormap([T_COLORS[x] for x in T_ORDER]))]
        if include_roi:
            bars.append((roi_bar, ListedColormap([ROI_COLORS[x] for x in ROI_ORDER])))
        gs = fig.add_gridspec(len(bars) + 1, 1, height_ratios=[0.042] * len(bars) + [1],
                              left=0.004, right=heatmap_right, bottom=0.025, top=0.985, hspace=0.055)
        heat_cmap = HEAT.copy(); heat_cmap.set_bad("white")
        for i, (bar, cmap) in enumerate(bars):
            cmap.set_bad("white")
            ax = fig.add_subplot(gs[i]); ax.imshow(bar, aspect="auto", interpolation="nearest", cmap=cmap); ax.axis("off")
        ax = fig.add_subplot(gs[len(bars)])
        ax.imshow(display_matrix, aspect="auto", interpolation="nearest", cmap=heat_cmap, vmin=-2, vmax=2)
        if horizontal_boundaries:
            for boundary in horizontal_boundaries:
                ax.axhline(boundary - 0.5, color="white", lw=0.35, solid_capstyle="butt", zorder=5)
        if row_labels:
            ax.set_xticks([])
            ax.set_yticks(np.arange(len(ordered_features)))
            ax.set_yticklabels(ordered_features, fontsize=4, color="#333333")
            ax.yaxis.tick_right()
            ax.tick_params(axis="y", length=0, pad=1.2)
            for spine in ax.spines.values():
                spine.set_visible(False)
        else:
            ax.axis("off")
        stem = OUT / stem_name
        fig.savefig(stem.with_suffix(".png"), dpi=900, facecolor="white")
        fig.savefig(stem.with_suffix(".pdf"), facecolor="white")
        plt.close(fig)
        print(stem.with_suffix(".png"))

    render(matrix, "E_consensus368_final19_Z2_heatmap_diagonal_MagentaBlackYellow_W4p64_H1")
    # Minimal-gap version: thin vector rules mark feature blocks without
    # consuming a full raster row.
    peak_counts = row_info["peak_class"].value_counts().reindex(ORDER).fillna(0).astype(int)
    boundaries = np.cumsum(peak_counts.to_numpy())[:-1].tolist()
    render(matrix, "E_consensus368_final19_Z2_heatmap_diagonal_minGap_MagentaBlackYellow_W4p64_H1",
           horizontal_boundaries=boundaries)
    render(matrix,
           "E_consensus368_final19_Z2_heatmap_diagonal_minGap_MagentaBlackYellow_3annotations_W4p64_H1",
           horizontal_boundaries=boundaries, include_roi=True)
    render(matrix,
           "E_consensus368_final19_Z2_heatmap_diagonal_minGap_MagentaBlackYellow_3annotations_rowLabels_W6p2_H1p35",
           horizontal_boundaries=boundaries, include_roi=True, row_labels=True)
    print(row_info["peak_class"].value_counts().reindex(ORDER).to_string())


if __name__ == "__main__":
    main()
