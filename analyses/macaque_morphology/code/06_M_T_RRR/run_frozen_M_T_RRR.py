#!/usr/bin/env python3
"""Macaque transcriptomic-to-frozen-morphology RRR four-panel display."""
from pathlib import Path
import sys

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import plot_ca_pu_nac_rrr_fourpanel as ref
from analyze_nac_d1d2_transcriptomic_ephys_rrr import (
    correlation_loadings, fit_rrr, orient_axes, read_h5ad_subset, transcriptomic_pca,
)

RUN = ROOT / "macaque_m" / "m18_tempfreeze_NPC5_HCK4_res2.3"
ASSIGN = RUN / "01_temp_frozen_assignments_126.csv"
ZFILE = ROOT / "macaque_m" / "m18_adaptive_pca126" / "02_transformed_z_117.csv"
OUT = RUN / "RRR_T_M"
PAIR12 = (1, 0)
PAIR13 = (2, 0)
T_ORDER = ["D1", "D2", "Hybrid"]
T_COLORS = {"D1": "#D95F02", "D2": "#008F7A", "Hybrid": "#6E7180"}
ELLIPSE_T_ORDER = ["D1", "D2"]
M_ORDER = ["M1", "M2", "M3", "M4"]
M_COLORS = {"M1": "#1F77B4", "M2": "#D9A400", "M3": "#8C564B", "M4": "#E377C2"}
M_LABELS = dict(pd.read_csv(ROOT / "macaque_m" / "morphology_feature_abbreviations.csv").values)


def main():
    mpl.rcParams.update({"font.family": "Arial", "font.size": 4,
                         "pdf.fonttype": 42, "ps.fonttype": 42})
    OUT.mkdir(parents=True, exist_ok=True)
    assign = pd.read_csv(ASSIGN, dtype={"cell_label": str})
    z = pd.read_csv(ZFILE, dtype={"cell_label": str})
    features = [x for x in z.columns if x != "cell_label"]
    data = assign.loc[assign.concordant.eq(True),
                      ["cell_label", "HC_K4", "Subclass"]].copy()
    data["T_class"] = data.Subclass.map({"STR D1 MSN": "D1", "STR D2 MSN": "D2",
                                          "STR Hybrid MSN": "Hybrid"})
    data["M_class"] = "M" + data.HC_K4.astype(str)
    data = data[["cell_label", "T_class", "M_class"]].merge(
        z, on="cell_label", validate="one_to_one")
    if len(data) != 117 or len(features) != 18 or data.T_class.isna().any():
        raise RuntimeError(f"Expected n=117, p=18; got n={len(data)}, p={len(features)}")

    expression, all_genes = read_h5ad_subset(data.cell_label.tolist())
    t_scores, _, genes, variance, gene_matrix = transcriptomic_pca(expression, all_genes)
    m = data[features].to_numpy(float)
    xm, ym, _, v_raw = fit_rrr(t_scores, m, rank=3)
    v = orient_axes(v_raw, features)
    beta = np.linalg.pinv(t_scores - xm) @ (m - ym)
    t_projection = (t_scores - xm) @ beta @ v
    m_projection = (m - ym) @ v
    t_loadings = correlation_loadings(gene_matrix, t_projection)
    m_loadings = correlation_loadings(m, m_projection)
    m_names = np.array([M_LABELS.get(x, x) for x in features], dtype=object)
    gene_names = np.array([
        x if len(x) <= 8 else f"{x[:4]}..{x[-2:]}" for x in genes
    ], dtype=object)
    t_groups = data.T_class.to_numpy(str)
    m_groups = data.M_class.to_numpy(str)
    m_display = np.where(t_groups == "Hybrid", "Hybrid", m_groups)
    m_display_colors = {**M_COLORS, "Hybrid": T_COLORS["Hybrid"]}

    ref.T_ORDER, ref.T_COLORS = ELLIPSE_T_ORDER, T_COLORS
    fig, axes = plt.subplots(1, 4, figsize=(3.36, 1.24),
                             gridspec_kw={"wspace": .34})
    ref.panel(axes[0], t_projection, t_loadings, gene_names, t_groups, T_ORDER, T_COLORS,
              t_groups, 4, PAIR12, vector_fontsize=4, vector_gap=.27,
              vector_force_side=-1, x_extent=1.45)
    ref.panel(axes[1], t_projection, t_loadings, gene_names, t_groups, T_ORDER, T_COLORS,
              t_groups, 4, PAIR13, vector_fontsize=4, vector_gap=.27,
              vector_force_side=-1, x_extent=1.45)
    ref.panel(axes[2], m_projection, m_loadings, m_names, m_display, M_ORDER + ["Hybrid"], m_display_colors,
              t_groups, 4, PAIR12, vector_fontsize=4, vector_gap=.27,
              vector_force_side=1, x_extent=1.45)
    ref.panel(axes[3], m_projection, m_loadings, m_names, m_display, M_ORDER + ["Hybrid"], m_display_colors,
              t_groups, 4, PAIR13, vector_fontsize=4, vector_gap=.27,
              vector_force_side=1, x_extent=1.45)
    fig.subplots_adjust(left=.055, right=.93, bottom=.02, top=.98)
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    label_boxes = [text.get_window_extent(renderer) for ax in axes for text in ax.texts
                   if text.get_text().strip()]
    canvas = fig.bbox
    clipped = [i for i, box in enumerate(label_boxes)
               if box.x0 < canvas.x0 or box.y0 < canvas.y0
               or box.x1 > canvas.x1 or box.y1 > canvas.y1]
    overlaps = [(i, j) for i, box in enumerate(label_boxes)
                for j in range(i + 1, len(label_boxes))
                if box.overlaps(label_boxes[j])]
    if clipped or overlaps:
        raise RuntimeError(f"Label QA failed: clipped={clipped}, overlaps={overlaps}")
    stem = OUT / "Macaque_CaPuNAC_RRR_T12_T13_M12_M13_frozen117_W3p36_H1p24_font4"
    fig.savefig(stem.with_suffix(".png"), dpi=1200, facecolor="white")
    fig.savefig(stem.with_suffix(".pdf"), facecolor="white")
    plt.close(fig)

    pd.DataFrame({"cell_label": data.cell_label, "T_class": data.T_class,
                  "M_class": data.M_class,
                  "T_Component1": t_projection[:, 0], "T_Component2": t_projection[:, 1],
                  "T_Component3": t_projection[:, 2], "M_Component1": m_projection[:, 0],
                  "M_Component2": m_projection[:, 1], "M_Component3": m_projection[:, 2]}).to_csv(
        OUT / "RRR_cell_scores_n117.csv", index=False)
    pd.DataFrame(m_loadings, index=features,
                 columns=["Component1", "Component2", "Component3"]).to_csv(
        OUT / "RRR_Mfeature_correlation_loadings.csv")
    pd.DataFrame(t_loadings, index=genes,
                 columns=["Component1", "Component2", "Component3"]).to_csv(
        OUT / "RRR_gene_correlation_loadings.csv")
    pd.DataFrame({"PC": np.arange(1, len(variance)+1),
                  "T_explained_variance_ratio": variance}).to_csv(
        OUT / "transcriptomic_PCA_variance.csv", index=False)
    pd.crosstab(data.T_class, data.M_class).to_csv(OUT / "Tclass_by_Mclass_counts.csv")
    (OUT / "README.txt").write_text(
        "Macaque Ca+Pu+NAC transcriptomic-to-morphology RRR display.\n"
        "Cells: 117 HC-GC consensus MSN cells (D1=47, D2=62, Hybrid=8).\n"
        "Response: 18 frozen transformed/Z-scored morphology features; 3_Sholl_PC1 excluded.\n"
        "Model: pooled T-to-M RRR, rank=3 for display; T PCA uses 1,000 eligible variable genes and 20 PCs.\n"
        "Panels: T Comp2/1, T Comp3/1, M Comp2/1, M Comp3/1.\n"
        "Hybrid cells are retained as gray points in all panels. Only D1/D2 receive class colors and 90% bivariate-normal ellipses.\n"
        "All cells are scaled by the 99th-percentile radius.\n"
        "The frozen M classification is not refit or changed. Higher RRR components are exploratory.\n",
        encoding="utf-8")
    (OUT / "compact_W3p36_H1p24_font4_QA.txt").write_text(
        "Canvas: 3.36 x 1.24 inches\n"
        "Font: Arial 4 pt for all vector labels\n"
        "Panels: four equal circular panels in one row\n"
        "Displayed labels: four strongest absolute loading vectors per panel\n"
        "Label clipping: 0\nLabel-label overlaps: 0\n"
        "RRR model, scores, cell set, ellipses, and saved full loadings: unchanged\n",
        encoding="utf-8")
    print(stem.with_suffix(".png"))
    print(pd.crosstab(data.T_class, data.M_class).to_string())


if __name__ == "__main__":
    main()
