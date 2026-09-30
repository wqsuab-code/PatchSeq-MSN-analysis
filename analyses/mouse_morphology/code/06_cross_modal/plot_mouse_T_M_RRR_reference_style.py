#!/usr/bin/env python3
"""Mouse NAc transcriptomic-to-morphology RRR, reference-style 2x2 panel."""
from pathlib import Path
import sys

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import plot_ca_pu_nac_rrr_fourpanel as ref
from analyze_nac_d1d2_transcriptomic_ephys_rrr import (
    correlation_loadings, fit_rrr, orient_axes,
)
from build_final_morph187_figure_suite import preprocess_train_test

ASSIGN = ROOT / "outputs/morph_qc/final_morph187_NPC3_HCK4_GCres0p50_figures/00_final_morph187_cell_assignments.csv"
RNA = ROOT / "outputs/morph_qc/morph_taxonomy_round5_biological_validation/43_transcriptomic_log1p_HVG2000_cell_by_gene.csv"
OUT = ROOT / "outputs/morph_qc/mouse_T_M_RRR_reference_style"

FEATURES = [
    "M_soma_circularity_index", "M_soma_aspect_ratio", "M_cell_max_radial_dist",
    "M_total_number_of_neurites", "M_basal_dendrite_avg_tortuosity",
    "M_Total_neurite_length_(sections)", "M_Number_of_bifurcation_points",
    "M_Maximum_branch_order", "M_trunk_angle_min", "M_trunk_angle_max",
]
M_NAMES = np.array([
    "Soma circularity", "Soma aspect ratio", "Max radial distance",
    "Primary neurites", "Mean tortuosity", "Total neurite length",
    "Bifurcations", "Max branch order", "Min trunk angle", "Max trunk angle",
])
T_ORDER = ["D1", "D2"]
T_COLORS = {"D1": "#ED0000", "D2": "#3C5488"}
M_ORDER = ["M1", "M2", "M3", "M4"]
M_COLORS = {"M1": "#00468B", "M2": "#42B540", "M3": "#ED0000", "M4": "#0099B4"}


def main() -> None:
    mpl.rcParams.update({"font.family": "Arial", "font.size": 4,
                         "pdf.fonttype": 42, "ps.fonttype": 42})
    OUT.mkdir(parents=True, exist_ok=True)

    a = pd.read_csv(ASSIGN)
    a = a[a.HC_GC_consensus.astype(str).str.lower().eq("true") & a.D1_D2.isin(T_ORDER)]
    rna = pd.read_csv(RNA)
    data = a.merge(rna, on="MSN_unique_ID", validate="one_to_one")
    if len(data) != 180:
        raise RuntimeError(f"Expected 180 Mouse consensus D1/D2 cells, found {len(data)}")

    # Re-rank the current label-blind HVG2000 within this 180-cell cohort and
    # retain 1,000 genes; scaling and PCA are independently refit here.
    candidates = [c for c in rna.columns if c != "MSN_unique_ID"]
    raw_t = data[candidates].to_numpy(float)
    variance = raw_t.var(axis=0, ddof=1)
    keep = np.argsort(variance)[::-1][:1000]
    genes = np.asarray(candidates, object)[keep]
    gene_matrix = raw_t[:, keep]
    gene_matrix = (gene_matrix - gene_matrix.mean(0)) / gene_matrix.std(0, ddof=0)
    gene_matrix = np.nan_to_num(gene_matrix)
    pca = PCA(n_components=20, svd_solver="full")
    t_scores = pca.fit_transform(gene_matrix)

    m_raw = data[FEATURES].to_numpy(float)
    m, _ = preprocess_train_test(m_raw, m_raw)
    xm, ym, _, v_raw = fit_rrr(t_scores, m, rank=3)
    v = orient_axes(v_raw, FEATURES)
    beta = np.linalg.pinv(t_scores - xm) @ (m - ym)
    t_projection = (t_scores - xm) @ beta @ v
    m_projection = (m - ym) @ v
    t_loadings = correlation_loadings(gene_matrix, t_projection)
    m_loadings = correlation_loadings(m, m_projection)

    t_group = data.D1_D2.to_numpy(str)
    m_group = data.M_class.to_numpy(str)
    ref.T_ORDER, ref.T_COLORS = T_ORDER, T_COLORS

    fig, axes = plt.subplots(2, 2, figsize=(3.36, 2.35),
                             gridspec_kw={"wspace": .34, "hspace": .16})
    pairs = [(0, 1), (0, 2)]  # Component 1 on x; Components 2/3 on y.
    for row, pair in enumerate(pairs):
        ref.panel(axes[row, 0], t_projection, t_loadings, genes,
                  t_group, T_ORDER, T_COLORS, t_group, 4, pair,
                  vector_fontsize=3.6, vector_gap=.24,
                  vector_force_side=-1, x_extent=1.32)
        ref.panel(axes[row, 1], m_projection, m_loadings, M_NAMES,
                  m_group, M_ORDER, M_COLORS, t_group, 4, pair,
                  vector_fontsize=3.6, vector_gap=.24,
                  vector_force_side=1, x_extent=1.32)

    axes[0, 0].set_title("Transcriptomic space", fontsize=5, pad=2)
    axes[0, 1].set_title("Morphological space", fontsize=5, pad=2)
    for ax in axes[1]:
        ax.set_xlabel("Component 1", fontsize=4, labelpad=1)
    axes[0, 0].set_ylabel("Component 2", fontsize=4, labelpad=1)
    axes[1, 0].set_ylabel("Component 3", fontsize=4, labelpad=1)
    fig.subplots_adjust(left=.12, right=.91, bottom=.08, top=.95)

    stem = OUT / "Mouse_NAc_T_M_RRR_rank3_T12_T13_M12_M13_n180"
    fig.savefig(stem.with_suffix(".png"), dpi=1200, facecolor="white")
    fig.savefig(stem.with_suffix(".pdf"), facecolor="white")
    fig.savefig(stem.with_suffix(".svg"), facecolor="white")
    plt.close(fig)

    pd.DataFrame({
        "MSN_unique_ID": data.MSN_unique_ID, "D1_D2": t_group, "M_class": m_group,
        **{f"T_Component{i+1}": t_projection[:, i] for i in range(3)},
        **{f"M_Component{i+1}": m_projection[:, i] for i in range(3)},
    }).to_csv(OUT / "RRR_cell_scores_n180.csv", index=False)
    pd.DataFrame(t_loadings, index=genes, columns=["Component1", "Component2", "Component3"]).to_csv(
        OUT / "RRR_gene_correlation_loadings.csv")
    pd.DataFrame(m_loadings, index=FEATURES, columns=["Component1", "Component2", "Component3"]).to_csv(
        OUT / "RRR_Mfeature_correlation_loadings.csv")
    pd.DataFrame({"PC": np.arange(1, 21), "variance_fraction": pca.explained_variance_ratio_}).to_csv(
        OUT / "transcriptomic_PCA_variance.csv", index=False)
    pd.crosstab(data.D1_D2, data.M_class).to_csv(OUT / "D1D2_by_Mclass_counts.csv")
    (OUT / "README.txt").write_text(
        "Mouse NAc T-to-M reduced-rank regression.\n"
        "Cohort: 180 HC-GC morphology-consensus cells with stable D1/D2 identity (D1=82, D2=98).\n"
        "T: top 1,000 variable genes from the current label-blind HVG2000 matrix; gene z-score; 20 PCs.\n"
        "M: 10 frozen core morphology features, transformed and z-scored within this cohort.\n"
        "RRR rank=3. Panels show Components 1-2 and 1-3 in independently fitted T and M spaces.\n"
        "Ellipses: 90% bivariate-normal D1/D2 regions. Display scaled by each panel's 99th-percentile radius.\n",
        encoding="utf-8")
    print(stem)
    print(pd.crosstab(data.D1_D2, data.M_class).to_string())


if __name__ == "__main__":
    main()
