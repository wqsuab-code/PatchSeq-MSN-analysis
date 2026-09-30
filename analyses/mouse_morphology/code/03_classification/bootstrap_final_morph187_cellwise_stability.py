#!/usr/bin/env python
"""Cell-wise bootstrap recovery for the frozen Morph n=187 taxonomy.

Each replicate samples 80% of cells without replacement and reruns the full
feature transform, PCA (NPC=3), and Ward.D2 HC (K=4). Bootstrap labels are
matched to the frozen M1--M4 labels by maximum overlap. Pairwise co-clustering
probabilities use only replicates in which both cells were sampled.
"""

from __future__ import annotations

from pathlib import Path
import os

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, ListedColormap
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage, leaves_list, fcluster
from scipy.optimize import linear_sum_assignment
from scipy.spatial.distance import pdist, squareform
from sklearn.decomposition import PCA


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs/morph_qc/final_morph187_NPC3_HCK4_GCres0p50_figures"
RAW = ROOT / "outputs/morph_qc/morph187_after_incomplete_reconstruction_quarantine/01_morph187_discovery_input.csv"
ASSIGN = OUT / "00_final_morph187_cell_assignments.csv"
ID = "MSN_unique_ID"
FEATURES = [
    "M_soma_circularity_index", "M_soma_aspect_ratio", "M_cell_max_radial_dist",
    "M_total_number_of_neurites", "M_basal_dendrite_avg_tortuosity",
    "M_Total_neurite_length_(sections)", "M_Number_of_bifurcation_points",
    "M_Maximum_branch_order", "M_trunk_angle_min", "M_trunk_angle_max",
]
LEVELS = ["M1", "M2", "M3", "M4"]
COLORS = {"M1": "#00468B", "M2": "#42B540", "M3": "#ED0000", "M4": "#0099B4"}
LABEL_COLUMN = "M_class"
FILE_PREFIX = "19"
N_BOOT = int(os.environ.get("MORPH_N_BOOT", "1000"))
FRACTION = 0.80
SEED = 20260826


def transform(x: np.ndarray) -> np.ndarray:
    shifted = np.maximum(x - np.nanmin(x, axis=0), 0.0)
    totals = np.nansum(shifted, axis=0)
    totals[totals <= 0] = 1.0
    logged = np.log1p(shifted / totals * 10000.0)
    means = np.nanmean(logged, axis=0)
    sds = np.nanstd(logged, axis=0, ddof=1)
    sds[sds <= 0] = 1.0
    return (logged - means) / sds


def matched_labels(labels: np.ndarray, truth: np.ndarray) -> np.ndarray:
    table = np.zeros((4, 4), dtype=int)
    for i in range(4):
        for j in range(4):
            table[i, j] = np.sum((labels == i) & (truth == j))
    rows, cols = linear_sum_assignment(-table)
    mapping = dict(zip(rows, cols))
    return np.array([mapping[v] for v in labels], dtype=np.int8)


def main() -> None:
    raw = pd.read_csv(RAW)
    assign = pd.read_csv(ASSIGN)[[ID, LABEL_COLUMN]]
    data = assign.merge(raw[[ID] + FEATURES], on=ID, validate="one_to_one")
    if len(data) != 187 or data[FEATURES].isna().any().any():
        raise RuntimeError("Expected complete frozen n=187 Morph cohort")

    x = data[FEATURES].to_numpy(float)
    truth = data[LABEL_COLUMN].map({m: i for i, m in enumerate(LEVELS)}).to_numpy(np.int8)
    if np.any(pd.isna(truth)):
        raise RuntimeError(f"Unexpected labels in {LABEL_COLUMN}; expected {LEVELS}")
    n = len(data)
    sample_n = int(round(FRACTION * n))
    rng = np.random.default_rng(SEED)
    sampled_n = np.zeros(n, dtype=np.uint16)
    recovered_n = np.zeros(n, dtype=np.uint16)
    pair_seen = np.zeros((n, n), dtype=np.uint16)
    pair_together = np.zeros((n, n), dtype=np.uint16)

    for rep in range(N_BOOT):
        idx = np.sort(rng.choice(n, sample_n, replace=False))
        z = transform(x[idx])
        pc = PCA(n_components=3, svd_solver="full").fit_transform(z)
        lab = fcluster(linkage(pc, method="ward", metric="euclidean"), t=4, criterion="maxclust") - 1
        mapped = matched_labels(lab.astype(np.int8), truth[idx])
        sampled_n[idx] += 1
        recovered_n[idx] += mapped == truth[idx]
        pair_seen[np.ix_(idx, idx)] += 1
        same = mapped[:, None] == mapped[None, :]
        pair_together[np.ix_(idx, idx)] += same.astype(np.uint16)
        if (rep + 1) % 100 == 0:
            print(f"bootstrap {rep + 1}/{N_BOOT}", flush=True)

    co = np.divide(pair_together, pair_seen, out=np.zeros_like(pair_together, dtype=float), where=pair_seen > 0)
    recovery = recovered_n / sampled_n
    within = np.empty(n)
    max_other = np.empty(n)
    for i in range(n):
        same = (truth == truth[i]) & (np.arange(n) != i)
        within[i] = np.mean(co[i, same])
        max_other[i] = max(np.mean(co[i, truth == j]) for j in range(4) if j != truth[i])

    # Retain final-class blocks, then order cells by their co-clustering structure.
    order_parts = []
    for j in range(4):
        idx = np.where(truth == j)[0]
        if len(idx) > 2:
            dist = np.clip(1.0 - co[np.ix_(idx, idx)], 0, 1)
            np.fill_diagonal(dist, 0)
            idx = idx[leaves_list(linkage(squareform(dist, checks=False), method="average"))]
        order_parts.append(idx)
    order = np.concatenate(order_parts)
    ordered_truth = truth[order]
    coo = co[np.ix_(order, order)]

    out = data[[ID, LABEL_COLUMN]].copy()
    out["Bootstrap_sampled_n"] = sampled_n
    out["Bootstrap_label_recovery"] = recovery
    out["Mean_within_M_coassignment"] = within
    out["Maximum_other_M_coassignment"] = max_other
    out["Cellwise_stability_margin"] = within - max_other
    out["Stability_status"] = np.select(
        [recovery >= 0.90, recovery >= 0.80],
        ["high", "intermediate"], default="low"
    )
    out["Heatmap_order"] = pd.Series(np.arange(1, n + 1), index=order).sort_index().to_numpy()
    out.sort_values([LABEL_COLUMN, "Bootstrap_label_recovery"], ascending=[True, False]).to_csv(
        OUT / f"{FILE_PREFIX}_cellwise_bootstrap_stability_{N_BOOT}x.csv", index=False
    )
    pd.DataFrame(co, index=data[ID], columns=data[ID]).to_csv(
        OUT / f"{FILE_PREFIX}_cell_cell_coclustering_probability_{N_BOOT}x.csv"
    )

    summary = out.groupby(LABEL_COLUMN, sort=False).agg(
        N=(ID, "size"),
        Median_recovery=("Bootstrap_label_recovery", "median"),
        Mean_recovery=("Bootstrap_label_recovery", "mean"),
        Minimum_recovery=("Bootstrap_label_recovery", "min"),
        N_high=("Stability_status", lambda s: int((s == "high").sum())),
        N_intermediate=("Stability_status", lambda s: int((s == "intermediate").sum())),
        N_low=("Stability_status", lambda s: int((s == "low").sum())),
    ).reset_index()
    summary.to_csv(
        OUT / f"{FILE_PREFIX}_cellwise_bootstrap_stability_by_class_summary_{N_BOOT}x.csv",
        index=False,
    )

    mpl.rcParams.update({
        "font.family": "Arial", "font.size": 4, "axes.linewidth": 0.55,
        "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
        "savefig.facecolor": "white", "figure.facecolor": "white",
    })
    cmap = LinearSegmentedColormap.from_list(
        "coassignment", ["#050509", "#3B0F70", "#C9A227", "#FFF7B2"], N=256
    )
    class_cmap = ListedColormap([COLORS[m] for m in LEVELS])
    sizes = [np.sum(ordered_truth == j) for j in range(4)]
    bounds = np.cumsum(sizes)[:-1]
    starts = np.r_[0, np.cumsum(sizes)[:-1]]
    centers = starts + (np.array(sizes) - 1) / 2

    fig = plt.figure(figsize=(2.35, 2.05))
    ax_left = fig.add_axes([0.12, 0.15, 0.035, 0.76])
    ax_top = fig.add_axes([0.17, 0.925, 0.64, 0.035])
    ax = fig.add_axes([0.17, 0.15, 0.64, 0.76])
    ax_stab = fig.add_axes([0.825, 0.15, 0.025, 0.76])
    cax = fig.add_axes([0.895, 0.15, 0.025, 0.76])
    ax_left.imshow(ordered_truth[:, None], aspect="auto", interpolation="nearest",
                   cmap=class_cmap, vmin=-0.5, vmax=3.5)
    ax_top.imshow(ordered_truth[None, :], aspect="auto", interpolation="nearest",
                  cmap=class_cmap, vmin=-0.5, vmax=3.5)
    im = ax.imshow(coo, cmap=cmap, vmin=0, vmax=1, interpolation="nearest",
                   aspect="equal", rasterized=True)
    stab_map = LinearSegmentedColormap.from_list("stability", ["#4D4D4D", "#FFFFFF"], N=256)
    ax_stab.imshow(recovery[order, None], cmap=stab_map, vmin=0, vmax=1,
                   interpolation="nearest", aspect="auto", rasterized=True)
    for b in bounds:
        ax.axhline(b - 0.5, color="#F2F2F2", lw=0.45)
        ax.axvline(b - 0.5, color="#F2F2F2", lw=0.45)
    ax.set_xticks(centers, LEVELS)
    ax.set_yticks([])
    ax.tick_params(length=0, pad=1)
    ax.set_xlabel("Cells ordered within final M class", labelpad=2)
    for a in (ax_left, ax_top, ax_stab):
        a.set_xticks([]); a.set_yticks([])
        for s in a.spines.values(): s.set_visible(False)
    ax_stab.set_title("Cell\nstability", fontsize=4, pad=2)
    ax_left.set_yticks(centers, LEVELS)
    ax_left.tick_params(length=0, pad=1)
    cb = fig.colorbar(im, cax=cax, ticks=[0, .5, 1])
    cb.set_label("Co-clustering probability", fontsize=4, labelpad=2)
    cb.ax.tick_params(labelsize=4, length=1, width=.4, pad=1)
    stem = OUT / f"{FILE_PREFIX}_cellwise_bootstrap_coclustering_heatmap_{N_BOOT}x"
    fig.savefig(stem.with_suffix(".png"), dpi=900, bbox_inches="tight", pad_inches=.02)
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight", pad_inches=.02)
    fig.savefig(stem.with_suffix(".svg"), bbox_inches="tight", pad_inches=.02)
    plt.close(fig)

    print(summary.to_string(index=False))
    print(out["Stability_status"].value_counts().to_string())
    print(stem)


if __name__ == "__main__":
    main()
