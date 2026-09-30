#!/usr/bin/env python3
"""Core-18 GC-HC consensus heatmap with a class-mean summary panel."""

from pathlib import Path
import base64

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Rectangle
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
FEATURE_FILE = ROOT / "data/02_frozen_input/NPC3_HC5_Raw_Final18_Features.csv"
LABEL_FILE = ROOT / "data/04_frozen_classification/HC_GC_only_cell_assignments.csv"
OUT = ROOT / "figures/generated/core18_heatmap"
OUT.mkdir(parents=True, exist_ok=True)

VIS_DIR = ROOT / "interactive/generated"
VIS_DIR.mkdir(parents=True, exist_ok=True)

SEED = 777
LEVELS = ["E1", "E2", "E3", "E4", "E5"]
COLORS = {
    "E1": "#4E79A7", "E2": "#F28E2B", "E3": "#59A14F",
    "E4": "#2AA6B8", "E5": "#B07AA1",
}
FEATURES = [
    ("E_Rheobase..pA.", "Rheobase"),
    ("E_Latency....20pA.current..ms.", "Latency (-20 pA)"),
    ("E_Latency..ms.", "Latency"),
    ("E_AP.amplitude.adaptation.index", "AP amplitude adaptation"),
    ("E_Input.resistance..MOhm.", "Input resistance"),
    ("E_Membrane.time.constant..ms.", "Membrane tau"),
    ("E_Holding.MP..mV.", "Holding MP"),
    ("E_AP.threshold..mV.", "AP threshold"),
    ("E_Max.number.of.APs", "Max APs"),
    ("E_AP.amplitude..mV.", "AP amplitude"),
    ("E_Upstroke.to.downstroke.ratio", "Upstroke/downstroke"),
    ("E_ISI.adaptation.index", "ISI adaptation"),
    ("E_ISI.coefficient.of.variation", "ISI CV"),
    ("E_Sag.ratio", "Sag ratio"),
    ("E_Sag.time..s.", "Sag time"),
    ("E_Afterhyperpolarization..mV.", "AHP"),
    ("E_AP.width..ms.", "AP width"),
    ("E_AP.coefficient.of.variation", "AP CV"),
]
CMAP = LinearSegmentedColormap.from_list(
    "ephys_magenta_black_yellow", ["#D600A9", "#180D19", "#FFFF00"], N=256
)


def setup_style():
    mpl.rcParams.update({
        "font.family": "Arial", "font.size": 5,
        "pdf.fonttype": 42, "ps.fonttype": 42,
    })


def load_data():
    raw = pd.read_csv(FEATURE_FILE)
    labels = pd.read_csv(LABEL_FILE)
    cols = [x[0] for x in FEATURES]
    missing = sorted(set(cols) - set(raw.columns))
    if missing:
        raise ValueError(f"Missing core features: {missing}")
    keep = labels["HC_GC_consensus"].astype(str).str.lower().eq("true")
    data = labels.loc[keep, ["MSN_unique_ID", "HC_GC_consensus_E"]].merge(
        raw[["MSN_unique_ID", *cols]], on="MSN_unique_ID",
        how="left", validate="one_to_one",
    )
    if len(data) != 450 or data[cols].isna().any().any():
        raise ValueError("Expected 450 complete GC-HC consensus cells.")

    x = data[cols].to_numpy(float)
    z_raw = (x - x.mean(axis=0)) / x.std(axis=0, ddof=1)
    classes = data["HC_GC_consensus_E"].astype(str).to_numpy()
    means = np.vstack([z_raw[classes == e].mean(axis=0) for e in LEVELS]).T

    rng = np.random.default_rng(SEED)
    order, blocks = [], []
    start = 0
    for e in LEVELS:
        idx = rng.permutation(np.flatnonzero(classes == e))
        order.extend(idx.tolist())
        blocks.append((start, len(idx), e))
        start += len(idx)
    ordered = data.iloc[order].reset_index(drop=True)
    z = z_raw[np.asarray(order)]
    return ordered, z, means, blocks


def top_bar(ax, blocks=None, class_cells=False):
    if class_cells:
        for j, e in enumerate(LEVELS):
            ax.add_patch(Rectangle((j - 0.5, -1.02), 1, 0.18,
                                   color=COLORS[e], ec="none", clip_on=False))
            ax.text(j, -1.22, e, ha="center", va="bottom", fontsize=4.5)
    else:
        for start, n, e in blocks:
            ax.add_patch(Rectangle((start - 0.5, -1.02), n, 0.18,
                                   color=COLORS[e], ec="none", clip_on=False))
            ax.text(start + (n - 1) / 2, -1.22, e,
                    ha="center", va="bottom", fontsize=4.5)
            if start:
                ax.axvline(start - 0.5, color="white", lw=0.55)


def save_figure(ordered, z, means, blocks, z_limit):
    fig = plt.figure(figsize=(7.0, 2.35), facecolor="white")
    ax = fig.add_axes([0.145, 0.12, 0.675, 0.78])
    axm = fig.add_axes([0.842, 0.12, 0.075, 0.78], sharey=ax)
    cax = fig.add_axes([0.942, 0.25, 0.010, 0.48])

    im = ax.imshow(np.clip(z, -z_limit, z_limit).T,
                   aspect="auto", interpolation="nearest",
                   cmap=CMAP, vmin=-z_limit, vmax=z_limit)
    axm.imshow(np.clip(means, -z_limit, z_limit),
               aspect="auto", interpolation="nearest",
               cmap=CMAP, vmin=-z_limit, vmax=z_limit)
    top_bar(ax, blocks=blocks)
    top_bar(axm, class_cells=True)

    ax.set_yticks(range(len(FEATURES)), [x[1] for x in FEATURES], fontsize=4.5)
    ax.set_xticks([])
    axm.set_xticks([])
    axm.tick_params(left=False, labelleft=False)
    ax.tick_params(length=0)
    for x in np.arange(0.5, 4.5, 1):
        axm.axvline(x, color="white", lw=0.4)
    for a in (ax, axm):
        for s in a.spines.values():
            s.set_visible(False)

    ax.text(0.5, 1.075, "Single cells (n=450)", transform=ax.transAxes,
            ha="center", va="bottom", fontsize=5)
    axm.text(0.5, 1.075, "Class mean", transform=axm.transAxes,
             ha="center", va="bottom", fontsize=5)
    cb = fig.colorbar(im, cax=cax)
    cb.set_ticks([-z_limit, 0, z_limit])
    cb.set_label("Feature Z-score", fontsize=4.5)
    cb.ax.tick_params(labelsize=4, length=1, width=0.4)
    cb.outline.set_linewidth(0.4)

    stem = f"HC_GC_consensus450_core18_heatmap_with_class_means_Zm{z_limit}_p{z_limit}"
    fig.savefig(OUT / f"{stem}.png", dpi=900, facecolor="white")
    fig.savefig(OUT / f"{stem}.pdf", facecolor="white")
    fig.savefig(OUT / f"{stem}.tif", dpi=900, facecolor="white",
                pil_kwargs={"compression": "tiff_lzw"})
    preview = OUT / f"{stem}_preview.png"
    fig.savefig(preview, dpi=180, facecolor="white")
    plt.close(fig)
    return preview, stem


def write_inline_preview(preview):
    b64 = base64.b64encode(preview.read_bytes()).decode("ascii")
    html = f'''<div id="core18-ephys-heatmap" style="width:100%;background:transparent;color:var(--foreground)">
  <h3 style="margin:0 0 8px 0;font-weight:500">Core electrophysiological signatures of five E-types</h3>
  <img src="data:image/png;base64,{b64}" alt="Heatmap of 18 electrophysiological features across 450 GC-HC consensus cells, with E1 through E5 class means shown at right." style="display:block;width:100%;height:auto" />
</div>
'''
    path = VIS_DIR / "core18-ephys-heatmap.html"
    path.write_text(html, encoding="utf-8")
    return path


def main():
    setup_style()
    ordered, z, means, blocks = load_data()
    previews = []
    for z_limit in (2, 3, 4):
        previews.append(save_figure(ordered, z, means, blocks, z_limit))
    preview, stem = previews[0]
    ordered[["MSN_unique_ID", "HC_GC_consensus_E"]].assign(
        Heatmap_column=np.arange(1, len(ordered) + 1),
        Within_group_order="seed777_random_shuffle",
    ).to_csv(OUT / "HC_GC_consensus450_core18_cell_order.csv", index=False)
    pd.DataFrame(means, index=[x[0] for x in FEATURES], columns=LEVELS).to_csv(
        OUT / "HC_GC_consensus450_core18_class_mean_zscores.csv"
    )
    html_path = write_inline_preview(preview)
    (OUT / "run_log.txt").write_text("\n".join([
        "Core-18 GC-HC consensus heatmap",
        "N=450 consensus cells; 43 nonconsensus cells excluded",
        "Features=18 frozen clustering features",
        "Columns=E1-E5; within-class random shuffle seed=777",
        "Display z-score ranges generated: -2 to 2, -3 to 3, and -4 to 4",
        "Right panel=class mean z-score for each feature",
        f"Inline preview={html_path}",
    ]), encoding="utf-8")


if __name__ == "__main__":
    main()
