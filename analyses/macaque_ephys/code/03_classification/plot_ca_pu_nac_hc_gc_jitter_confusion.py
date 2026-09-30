from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "outputs/R3_panels_v3/00_tuned_p80_ee12_random_seed777_coordinates.csv"
OUT = ROOT / "outputs/R3_confusion"
ORDER = ["C1", "C2", "C3", "C4"]
COLORS = {"C1": "#F8766D", "C2": "#7CAE00", "C3": "#00BFC4", "C4": "#C77CFF"}


def main():
    mpl.rcParams.update({"font.family": "Arial", "font.size": 4,
                         "pdf.fonttype": 42, "ps.fonttype": 42})
    OUT.mkdir(parents=True, exist_ok=True)
    data = pd.read_csv(INPUT, dtype={"cell_label": str})
    table = pd.crosstab(data.HC_class, data.GC_merged).reindex(index=ORDER, columns=ORDER, fill_value=0)
    agreement = int(np.trace(table.to_numpy()))
    n = int(table.to_numpy().sum())
    rng = np.random.default_rng(777)

    fig = plt.figure(figsize=(1.22, 1.0), facecolor="white")
    ax = fig.add_axes([.25, .18, .66, .66])
    for row_i, hc in enumerate(ORDER):
        y_cell = 3 - row_i
        for col_i, gc in enumerate(ORDER):
            ax.add_patch(Rectangle((col_i, y_cell), 1, 1, facecolor="white",
                                   edgecolor="#555555", linewidth=.48, zorder=0))
            count = int(table.loc[hc, gc])
            if count:
                sd = .105 if hc == gc else .145
                x = np.clip(rng.normal(col_i + .5, sd, count), col_i + .12, col_i + .88)
                y = np.clip(rng.normal(y_cell + .5, sd, count), y_cell + .12, y_cell + .88)
                ax.scatter(x, y, s=1.7, color="#202020", alpha=.50,
                           edgecolors="none", linewidths=0, zorder=2)

    for row_i, group in enumerate(ORDER):
        y_cell = 3 - row_i
        ax.add_patch(Rectangle((-.16, y_cell), .095, 1, facecolor=COLORS[group], edgecolor="none"))
        ax.text(-.27, y_cell + .5, f"HC{row_i + 1}", ha="right", va="center",
                fontsize=4, color=COLORS[group])
    for col_i, group in enumerate(ORDER):
        ax.add_patch(Rectangle((col_i, -.16), 1, .095, facecolor=COLORS[group], edgecolor="none"))
        ax.text(col_i + .5, -.28, f"GC{col_i + 1}", ha="center", va="top",
                fontsize=4, color=COLORS[group])

    ax.set_xlim(-.60, 4.02); ax.set_ylim(-.48, 4.03)
    ax.set_aspect("equal"); ax.axis("off")
    fig.text(.58, .965, f"HC × GC agreement:\n({agreement}/{n}, {agreement/n*100:.2f}%)",
             ha="center", va="top", fontsize=4.1, color="#111111")
    stem = OUT / "HC4_vs_mergedGC4_res3_jitter_confusion_H1in"
    fig.savefig(stem.with_suffix(".png"), dpi=1200, facecolor="white")
    fig.savefig(stem.with_suffix(".pdf"), facecolor="white")
    plt.close(fig)

    table.to_csv(OUT / "HC4_vs_mergedGC4_res3_counts.csv")
    table.div(table.sum(axis=1), axis=0).mul(100).to_csv(
        OUT / "HC4_vs_mergedGC4_res3_row_percent.csv")
    data[["cell_label", "HC_class", "GC_merged", "Consensus"]].to_csv(
        OUT / "HC4_vs_mergedGC4_res3_cell_assignments.csv", index=False)
    print(stem.with_suffix(".png"))
    print(table.to_string())


if __name__ == "__main__":
    main()
