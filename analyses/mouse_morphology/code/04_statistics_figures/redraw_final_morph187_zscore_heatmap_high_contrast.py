#!/usr/bin/env python
"""Redraw the frozen Morph heatmap with a high-contrast banded colour scale."""
from pathlib import Path
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm, ListedColormap
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs/morph_qc/final_morph187_NPC3_HCK4_GCres0p50_figures"
ZFILE = OUT / "08_heatmap_ordered_zscore_matrix.csv"
ASSIGN = OUT / "00_final_morph187_cell_assignments.csv"
LEVELS = ["M1", "M2", "M3", "M4"]
COLORS = {"M1": "#00468B", "M2": "#42B540", "M3": "#ED0000", "M4": "#0099B4"}
DISPLAY = [
    "Soma circularity", "Soma aspect ratio", "Max radial distance",
    "Primary neurite number", "Mean tortuosity", "Total neurite length",
    "Bifurcation points", "Maximum branch order", "Minimum trunk angle",
    "Maximum trunk angle",
]

def main():
    z = pd.read_csv(ZFILE, index_col=0)
    a = pd.read_csv(ASSIGN).set_index("MSN_unique_ID")
    classes = a.loc[z.index, "M_class"]
    mat = np.clip(z.to_numpy(float).T, -2, 2)
    bounds = np.array([-2, -1.5, -1.0, -0.5, -0.2, 0.2, 0.5, 1.0, 1.5, 2.0001])
    cmap = ListedColormap([
        "#6F00FF", "#8F00D8", "#69108D", "#351044", "#050505",
        "#625A00", "#A99B00", "#E1CD00", "#FFF200",
    ])
    norm = BoundaryNorm(bounds, cmap.N, clip=True)
    mpl.rcParams.update({"font.family":"Arial", "font.size":4, "axes.linewidth":.6,
                         "pdf.fonttype":42, "ps.fonttype":42, "svg.fonttype":"none"})
    fig = plt.figure(figsize=(7.5, 1.55))
    ax = fig.add_axes([.015, .08, .755, .82])
    cax = fig.add_axes([.955, .08, .012, .82])
    im = ax.imshow(mat, aspect="auto", cmap=cmap, norm=norm, interpolation="nearest")
    sizes = [int((classes == m).sum()) for m in LEVELS]
    boundaries = np.cumsum(sizes)
    for b in boundaries[:-1]: ax.axvline(b-.5, color="white", lw=2.4)
    for y in np.arange(.5, len(DISPLAY), 1): ax.axhline(y, color="white", lw=.15, alpha=.28)
    start = 0
    for m, n in zip(LEVELS, sizes):
        ax.add_patch(plt.Rectangle((start-.5, -.9), n, .35, color=COLORS[m], ec="none", clip_on=False))
        start += n
    ax.set_yticks(range(len(DISPLAY)), DISPLAY, fontsize=4)
    ax.yaxis.tick_right(); ax.tick_params(axis="y", length=0, pad=2); ax.set_xticks([])
    for s in ax.spines.values(): s.set_visible(False)
    cb = fig.colorbar(im, cax=cax, boundaries=bounds,
                      ticks=[-2,-1,-.5,0,.5,1,2], spacing="proportional")
    cb.set_label("Z score", fontsize=4); cb.ax.tick_params(labelsize=4, length=1.5, width=.45)
    fig.savefig(OUT/"08_final_M1-M4_10feature_Zscore_heatmap.png", dpi=900, bbox_inches="tight", pad_inches=.02)
    fig.savefig(OUT/"08_final_M1-M4_10feature_Zscore_heatmap.pdf", bbox_inches="tight", pad_inches=.02)
    fig.savefig(OUT/"08_final_M1-M4_10feature_Zscore_heatmap.svg", bbox_inches="tight", pad_inches=.02)
    plt.close(fig)

if __name__ == "__main__": main()
