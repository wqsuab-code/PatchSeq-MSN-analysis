#!/usr/bin/env python
"""PCA and classical three-standard-deviation outlier screen for Morph n=193."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA


ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "outputs/morph_qc/morph_taxonomy_round6_blind_de_novo/01_blinded_193cell_morphology_input.csv"
PRIORITY = ROOT / "outputs/morph_qc/morph193_remaining_extreme_influence_audit/04_priority_review_cells_joint_extreme_and_influential.csv"
OUT = ROOT / "outputs/morph_qc/morph193_pca_3sd_outlier_screen"
ID = "MSN_unique_ID"
FEATURES = [
    "M_soma_circularity_index",
    "M_soma_aspect_ratio",
    "M_cell_max_radial_dist",
    "M_total_number_of_neurites",
    "M_basal_dendrite_avg_tortuosity",
    "M_Total_neurite_length_(sections)",
    "M_Number_of_bifurcation_points",
    "M_Maximum_branch_order",
    "M_trunk_angle_min",
    "M_trunk_angle_max",
]
SHORT = [
    "Soma circularity", "Soma aspect ratio", "Max radial distance",
    "Primary neurites", "Mean tortuosity", "Total neurite length",
    "Bifurcation points", "Maximum branch order", "Min trunk angle",
    "Max trunk angle",
]


def preprocess(x: np.ndarray) -> np.ndarray:
    shifted = x - np.min(x, axis=0)
    total = shifted.sum(axis=0)
    if np.any(total <= 0):
        raise RuntimeError("Non-positive shifted feature sum")
    logged = np.log1p(shifted / total * 10000.0)
    return (logged - logged.mean(axis=0)) / logged.std(axis=0, ddof=1)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    dat = pd.read_csv(INPUT)
    if len(dat) != 193 or dat[ID].duplicated().any():
        raise RuntimeError("Expected 193 unique cells")
    x = dat[FEATURES].to_numpy(float)
    if not np.isfinite(x).all():
        raise RuntimeError("Non-finite values in input")

    z = preprocess(x)
    pca = PCA(n_components=len(FEATURES), svd_solver="full").fit(z)
    scores = pca.transform(z)
    score_sd = scores.std(axis=0, ddof=1)
    score_z = scores / score_sd
    flags = np.abs(score_z) > 3.0

    pc_names = [f"PC{i}" for i in range(1, len(FEATURES) + 1)]
    score_df = pd.DataFrame(scores, columns=[f"{p}_score" for p in pc_names])
    z_df = pd.DataFrame(score_z, columns=[f"{p}_SD" for p in pc_names])
    flag_df = pd.DataFrame(flags, columns=[f"{p}_outside_3SD" for p in pc_names])
    out = pd.concat([dat[[ID]].reset_index(drop=True), score_df, z_df, flag_df], axis=1)
    out["Any_PC_outside_3SD"] = flags.any(axis=1)
    out["N_PCs_outside_3SD"] = flags.sum(axis=1)
    out["Max_abs_PC_SD"] = np.abs(score_z).max(axis=1)
    out["Most_extreme_PC"] = np.array(pc_names)[np.abs(score_z).argmax(axis=1)]

    priority_ids: set[str] = set()
    if PRIORITY.exists():
        priority_ids = set(pd.read_csv(PRIORITY)[ID].astype(str))
    out["Previous_joint_robust_priority"] = out[ID].astype(str).isin(priority_ids)
    out["Overlap_3SD_and_previous_priority"] = out["Any_PC_outside_3SD"] & out["Previous_joint_robust_priority"]
    out.to_csv(OUT / "02_pc_scores_and_3SD_flags_all193.csv", index=False)

    variance = pd.DataFrame({
        "PC": pc_names,
        "Explained_variance_percent": pca.explained_variance_ratio_ * 100,
        "Cumulative_variance_percent": np.cumsum(pca.explained_variance_ratio_) * 100,
        "PC_score_SD": score_sd,
    })
    variance.to_csv(OUT / "01_pca_explained_variance.csv", index=False)

    loadings = pd.DataFrame(
        pca.components_.T,
        index=SHORT,
        columns=pc_names,
    )
    loadings.index.name = "Feature"
    loadings.reset_index().to_csv(OUT / "03_pca_loadings.csv", index=False)

    flagged = out[out["Any_PC_outside_3SD"]].copy().sort_values(
        ["Max_abs_PC_SD", ID], ascending=[False, True]
    )
    flagged.to_csv(OUT / "04_cells_outside_3SD_on_any_PC.csv", index=False)

    mpl.rcParams.update({
        "font.family": "Arial", "font.size": 6, "axes.linewidth": 0.8,
        "pdf.fonttype": 42, "ps.fonttype": 42,
    })
    fig, axes = plt.subplots(2, 2, figsize=(7.0, 6.0))
    ax = axes[0, 0]
    ax.bar(np.arange(1, 11), variance["Explained_variance_percent"], color="#4C78A8", width=0.72)
    ax2 = ax.twinx()
    ax2.plot(np.arange(1, 11), variance["Cumulative_variance_percent"], color="#E45756", marker="o", ms=2.5, lw=0.8)
    ax.set_xlabel("Principal component")
    ax.set_ylabel("Explained variance (%)")
    ax2.set_ylabel("Cumulative variance (%)")
    ax.set_xticks(np.arange(1, 11))
    ax.spines["top"].set_visible(False)
    ax2.spines["top"].set_visible(False)

    pairs = [(0, 1), (0, 2), (1, 2)]
    for ax, (i, j) in zip(axes.flat[1:], pairs):
        pair_flag = flags[:, i] | flags[:, j]
        colours = np.where(pair_flag, "#ED0000", "#BDBDBD")
        ax.scatter(score_z[:, i], score_z[:, j], c=colours, s=np.where(pair_flag, 18, 8), edgecolors="none")
        for v in (-3, 3):
            ax.axvline(v, color="#555555", ls="--", lw=0.55)
            ax.axhline(v, color="#555555", ls="--", lw=0.55)
        for row_idx in np.where(pair_flag)[0]:
            ax.annotate(dat.iloc[row_idx][ID], (score_z[row_idx, i], score_z[row_idx, j]),
                        xytext=(2, 2), textcoords="offset points", fontsize=4.5)
        ax.set_xlabel(f"PC{i+1} ({pca.explained_variance_ratio_[i]*100:.1f}%) [SD units]")
        ax.set_ylabel(f"PC{j+1} ({pca.explained_variance_ratio_[j]*100:.1f}%) [SD units]")
        ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(OUT / "05_PCA_and_3SD_outlier_overview.png", dpi=600, bbox_inches="tight")
    fig.savefig(OUT / "05_PCA_and_3SD_outlier_overview.pdf", bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(4.5, 3.1))
    im = ax.imshow(loadings.iloc[:, :6].to_numpy(), aspect="auto", cmap="RdBu_r", vmin=-0.7, vmax=0.7)
    ax.set_xticks(range(6), pc_names[:6])
    ax.set_yticks(range(len(SHORT)), SHORT)
    cbar = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.025)
    cbar.set_label("PCA loading")
    fig.tight_layout()
    fig.savefig(OUT / "06_PC1-PC6_loading_heatmap.png", dpi=600, bbox_inches="tight")
    fig.savefig(OUT / "06_PC1-PC6_loading_heatmap.pdf", bbox_inches="tight")
    plt.close(fig)

    summary = {
        "cells": len(dat),
        "features": FEATURES,
        "definition": "A cell is flagged when abs(PC score / sample SD of that PC) > 3 on any PC1-PC10.",
        "any_PC1_PC10_3SD_n": int(flagged.shape[0]),
        "any_PC1_PC3_3SD_n": int(flags[:, :3].any(axis=1).sum()),
        "flagged_cells": flagged[ID].tolist(),
        "overlap_with_previous_11_n": int(flagged["Overlap_3SD_and_previous_priority"].sum()),
        "overlap_with_previous_11_cells": flagged.loc[flagged["Overlap_3SD_and_previous_priority"], ID].tolist(),
        "automatic_deletion_n": 0,
    }
    (OUT / "00_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))
    print("\nFlagged details:")
    print(flagged[[ID, "Most_extreme_PC", "Max_abs_PC_SD", "N_PCs_outside_3SD", "Previous_joint_robust_priority"]].to_string(index=False))
    print("\nVariance:")
    print(variance.head(6).to_string(index=False))


if __name__ == "__main__":
    main()
