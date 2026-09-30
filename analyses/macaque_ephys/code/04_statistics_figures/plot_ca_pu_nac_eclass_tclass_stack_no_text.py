from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "outputs/R3_panels_v3/00_tuned_p80_ee12_random_seed777_coordinates.csv"
OUT = ROOT / "outputs/R3_composition"
E_ORDER = ["C1", "C2", "C3", "C4"]
T_ORDER = ["D1", "D2", "Hybrid"]
COLORS = {"D1": "#D95F02", "D2": "#008F7A", "Hybrid": "#6B7280"}


def main():
    mpl.rcParams.update({"font.family": "Arial", "pdf.fonttype": 42, "ps.fonttype": 42})
    OUT.mkdir(parents=True, exist_ok=True)
    data = pd.read_csv(INPUT)
    data = data.loc[data["Consensus"].eq(True)]
    counts = pd.crosstab(data.HC_class, data.T_class).reindex(
        index=E_ORDER, columns=T_ORDER, fill_value=0)
    percent = counts.div(counts.sum(axis=1), axis=0).mul(100)

    fig, ax = plt.subplots(figsize=(1.25, 1.0), facecolor="white")
    bottom = np.zeros(len(E_ORDER))
    x = np.arange(len(E_ORDER))
    for group in T_ORDER:
        values = percent[group].to_numpy(float)
        ax.bar(x, values, bottom=bottom, width=.70, color=COLORS[group],
               edgecolor="none", linewidth=0)
        bottom += values

    ax.set_xlim(-.55, 3.55); ax.set_ylim(0, 100)
    ax.set_xticks(x); ax.set_xticklabels([])
    ax.set_yticks([0, 50, 100]); ax.set_yticklabels([])
    ax.tick_params(axis="both", which="major", direction="out", length=2.0,
                   width=.55, color="#111111", pad=0)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color("#111111"); ax.spines["bottom"].set_color("#111111")
    ax.spines["left"].set_linewidth(.55); ax.spines["bottom"].set_linewidth(.55)
    ax.set_xlabel(""); ax.set_ylabel(""); ax.set_title("")
    ax.set_position([.16, .14, .80, .82])

    stem = OUT / "Eclass_Tclass_100pct_noText_axesTicks_W1p25_H1"
    fig.savefig(stem.with_suffix(".png"), dpi=1200, facecolor="white")
    fig.savefig(stem.with_suffix(".pdf"), facecolor="white")
    plt.close(fig)
    counts.to_csv(OUT / "Eclass_Tclass_counts.csv")
    percent.to_csv(OUT / "Eclass_Tclass_percent.csv")
    print(stem.with_suffix(".png"))
    print(counts.to_string())


if __name__ == "__main__":
    main()
