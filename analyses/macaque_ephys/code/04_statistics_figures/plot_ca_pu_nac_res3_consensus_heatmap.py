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
IMPORTANCE = OUT / "11_radar_feature_ranking_eta_squared.csv"

ORDER = ["C1", "C2", "C3", "C4"]
CLASS_COLORS = {"C1": "#F8766D", "C2": "#7CAE00", "C3": "#00BFC4", "C4": "#C77CFF"}
T_ORDER = ["D1", "D2", "Hybrid"]
T_COLORS = {"D1": "#E66101", "D2": "#008F7A", "Hybrid": "#6B7280"}
HEAT = LinearSegmentedColormap.from_list("purple_black_yellow", ["#C77CFF", "#000000", "#FFD92F"])


def main():
    mpl.rcParams.update({"font.family": "Arial", "font.size": 4, "pdf.fonttype": 42, "ps.fonttype": 42})
    assign = pd.read_csv(ASSIGN, dtype={"cell_label": str})
    z = pd.read_csv(ZFILE, dtype={"cell_label": str})
    importance = pd.read_csv(IMPORTANCE)
    features = importance["feature"].tolist()
    data = assign.loc[assign["Consensus"].eq(True), ["cell_label", "HC_class", "T_class"]].merge(
        z[["cell_label"] + features], on="cell_label", how="inner", validate="one_to_one"
    )
    if len(data) != 368 or data[features].isna().any().any():
        raise RuntimeError(f"Expected 368 complete consensus cells; got {len(data)}")

    rng = np.random.default_rng(777)
    matrix_parts, class_parts, t_parts, order_rows = [], [], [], []
    for i, group in enumerate(ORDER):
        block = data.loc[data["HC_class"].eq(group)].copy()
        block = block.iloc[rng.permutation(len(block))]
        matrix_parts.append(block[features].to_numpy(float).T)
        class_parts.append(np.full((1, len(block)), i, float))
        t_parts.append(np.array([[T_ORDER.index(x) for x in block["T_class"]]], float))
        order_rows.append(block[["cell_label", "HC_class", "T_class"]])
        if i < len(ORDER) - 1:
            matrix_parts.append(np.full((len(features), 2), np.nan))
            class_parts.append(np.full((1, 2), np.nan))
            t_parts.append(np.full((1, 2), np.nan))

    matrix = np.concatenate(matrix_parts, axis=1)
    class_bar = np.concatenate(class_parts, axis=1)
    t_bar = np.concatenate(t_parts, axis=1)
    pd.concat(order_rows, ignore_index=True).to_csv(
        OUT / "E_consensus368_final19_heatmap_cell_order.csv", index=False
    )
    pd.DataFrame({"row": np.arange(1, len(features) + 1), "feature": features}).to_csv(
        OUT / "E_consensus368_final19_heatmap_feature_order.csv", index=False
    )

    fig = plt.figure(figsize=(4.64, 1.0), facecolor="white")
    gs = fig.add_gridspec(
        3, 1, height_ratios=[0.042, 0.042, 1.0],
        left=0.004, right=0.996, bottom=0.025, top=0.985, hspace=0.055
    )
    class_cmap = ListedColormap([CLASS_COLORS[x] for x in ORDER]); class_cmap.set_bad("white")
    t_cmap = ListedColormap([T_COLORS[x] for x in T_ORDER]); t_cmap.set_bad("white")
    heat_cmap = HEAT.copy(); heat_cmap.set_bad("white")
    for row, cmap, position in [(class_bar, class_cmap, 0), (t_bar, t_cmap, 1)]:
        ax = fig.add_subplot(gs[position])
        ax.imshow(row, aspect="auto", interpolation="nearest", cmap=cmap)
        ax.set_axis_off()
    ax = fig.add_subplot(gs[2])
    ax.imshow(matrix, aspect="auto", interpolation="nearest", cmap=heat_cmap, vmin=-2, vmax=2)
    ax.set_axis_off()

    stem = OUT / "E_consensus368_final19_Z2_heatmap_2annotations_W4p64_H1"
    fig.savefig(stem.with_suffix(".png"), dpi=900, facecolor="white")
    fig.savefig(stem.with_suffix(".pdf"), facecolor="white")
    plt.close(fig)
    print(stem.with_suffix(".png"))
    print(stem.with_suffix(".pdf"))
    print(data["HC_class"].value_counts().sort_index().to_string())


if __name__ == "__main__":
    main()
