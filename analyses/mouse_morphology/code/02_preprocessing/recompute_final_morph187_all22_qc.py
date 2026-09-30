#!/usr/bin/env python
"""Recompute all-22 morphology redundancy and skewness QC in the final n=187 cohort."""

from __future__ import annotations

import json
import textwrap
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse, FancyBboxPatch, Rectangle
from PIL import Image
from scipy.stats import skew

import rerun_morph_all22_redundancy_skew_after_quarantine22 as audit


ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "outputs/morph_qc/all22_shift_sum10000_log1p_zscore_215cells/01_all22_raw_215cells.csv"
FINAL = ROOT / "outputs/morph_qc/morph187_after_incomplete_reconstruction_quarantine/01_morph187_discovery_input.csv"
OUT = ROOT / "outputs/morph_qc/final_morph187_all22_recomputed_redundancy_skew"
ID = "MSN_unique_ID"
FINAL10 = [
    "M_soma_circularity_index", "M_soma_aspect_ratio", "M_cell_max_radial_dist",
    "M_total_number_of_neurites", "M_basal_dendrite_avg_tortuosity",
    "M_Total_neurite_length_(sections)", "M_Number_of_bifurcation_points",
    "M_Maximum_branch_order", "M_trunk_angle_min", "M_trunk_angle_max",
]

FEATURE_UNITS = {
    "M_soma_average_radius": "µm",
    "M_soma_maximal_radius": "µm",
    "M_soma_minimal_radius": "µm",
    "M_soma_max_pairwise_dist": "µm",
    "M_soma_perimeter": "µm",
    "M_soma_area": "µm²",
    "M_soma_circularity_index": "Index",
    "M_soma_shape_factor": "Index",
    "M_soma_aspect_ratio": "Ratio",
    "M_cell_max_radial_dist": "µm",
    "M_total_number_of_neurites": "Count",
    "M_basal_dendrite_longthest_path": "µm",
    "M_basal_dendrite_Nseg": "Count",
    "M_basal_dendrite_avg_tortuosity": "Index",
    "M_basal_dendrite_max_tortuosity": "Index",
    "M_Total_neurite_length_(sections)": "µm",
    "M_Total_neurite_volume": "µm³",
    "M_Total_neurite_area": "µm²",
    "M_Number_of_bifurcation_points": "Count",
    "M_Maximum_branch_order": "Order",
    "M_trunk_angle_min": "rad",
    "M_trunk_angle_max": "rad",
}


def combine_png(left: Path, right: Path, out: Path) -> None:
    a = Image.open(left).convert("RGB")
    b = Image.open(right).convert("RGB")
    height = max(a.height, b.height)
    if a.height != height:
        a = a.resize((round(a.width * height / a.height), height), Image.Resampling.LANCZOS)
    if b.height != height:
        b = b.resize((round(b.width * height / b.height), height), Image.Resampling.LANCZOS)
    gap = round(height * 0.025)
    canvas = Image.new("RGB", (a.width + gap + b.width, height), "white")
    canvas.paste(a, (0, 0)); canvas.paste(b, (a.width + gap, 0))
    canvas.save(out, dpi=(600, 600))


def plot_aligned_redundancy_and_skewness(
    corr: pd.DataFrame,
    skew_df: pd.DataFrame,
    heatmap_order: list[str],
) -> None:
    """Draw the two panels on one shared row coordinate system."""
    ordered = corr.loc[heatmap_order, heatmap_order]
    data = skew_df.set_index("Feature").loc[heatmap_order].reset_index()
    n = len(heatmap_order)

    # Final two-column publication width.  Keeping 7-pt labels on a 7.2-in
    # canvas prevents them from being visually reduced when the full figure is
    # placed at journal width.
    fig = plt.figure(figsize=(7.2, 4.05), facecolor="white")
    # Both axes use exactly the same bottom and height.  The generous central
    # gap accommodates the single set of feature labels on the skewness panel.
    ax_corr = fig.add_axes([0.015, 0.045, 0.535, 0.925])
    ax_skew = fig.add_axes([0.810, 0.045, 0.170, 0.925])

    cmap = mpl.colors.LinearSegmentedColormap.from_list(
        "corr_aligned",
        [(0, "#0571B0"), (.35, "#74ADD1"), (.5, "#EEEEEE"),
         (.65, "#F28E6B"), (1, "#CA0020")],
    )
    norm = mpl.colors.Normalize(-1, 1)
    for row in range(n):
        for col in range(row, n):
            rho = float(ordered.iat[row, col])
            ax_corr.add_patch(Rectangle(
                (col, row), 1, 1, facecolor="white",
                edgecolor="#111111", linewidth=.38,
            ))
            ax_corr.add_patch(Ellipse(
                (col + .5, row + .5),
                width=.86,
                height=max(.045, .86 * (1 - abs(rho))),
                angle=45 if rho >= 0 else -45,
                facecolor=cmap(norm(rho)),
                edgecolor="none",
            ))
            ax_corr.text(
                col + .5, row + .5, f"{rho:.2f}",
                ha="center", va="center", fontsize=3.5,
            )

    # Reconstruct the exact |rho| >= 0.80 connected components and outline
    # them in the already component-contiguous heatmap order.
    parent = {feature: feature for feature in heatmap_order}

    def find(feature: str) -> str:
        while parent[feature] != feature:
            parent[feature] = parent[parent[feature]]
            feature = parent[feature]
        return feature

    def union(first: str, second: str) -> None:
        root_first, root_second = find(first), find(second)
        if root_first != root_second:
            parent[root_second] = root_first

    for i, first in enumerate(heatmap_order):
        for second in heatmap_order[i + 1:]:
            if abs(float(corr.loc[first, second])) >= audit.RHO_THRESHOLD:
                union(first, second)

    start = 0
    while start < n:
        root = find(heatmap_order[start])
        end = start + 1
        while end < n and find(heatmap_order[end]) == root:
            end += 1
        size = end - start
        if size >= 2:
            ax_corr.add_patch(FancyBboxPatch(
                (start - .08, start - .08), size + .16, size + .16,
                boxstyle="round,pad=0,rounding_size=0.18",
                fill=False, edgecolor="#111111", linewidth=1.15,
                clip_on=False, zorder=6,
            ))
        start = end

    ax_corr.set_xlim(-.05, n + .05)
    ax_corr.set_ylim(n + .05, -.05)
    ax_corr.set_aspect("equal")
    # No repeated feature names on the correlation panel.
    ax_corr.set_xticks([])
    ax_corr.set_yticks([])
    for spine in ax_corr.spines.values():
        spine.set_visible(False)

    cax = ax_corr.inset_axes([.018, .31, .018, .28])
    colorbar = fig.colorbar(
        mpl.cm.ScalarMappable(norm=norm, cmap=cmap),
        cax=cax, ticks=[-1, -.5, 0, .5, 1],
    )
    colorbar.ax.tick_params(labelsize=5, length=2, width=.5)
    colorbar.outline.set_linewidth(.6)

    y = np.arange(n) + .5
    colors = np.where(data["Signed_skewness"].lt(0), "#00468B", "#ED0000")
    ax_skew.hlines(y, 0, data["Absolute_skewness"], color=colors, lw=1.15)
    ax_skew.scatter(data["Absolute_skewness"], y, c=colors, s=18, zorder=3)
    ax_skew.axvline(2, color="#ED0000", ls="--", lw=.85)
    ax_skew.set_yticks(y, [audit.SHORT[f] for f in heatmap_order], fontsize=7.0)
    ax_skew.tick_params(axis="y", length=0, pad=2)
    # Center every feature name within one fixed label column to the left of
    # the skewness axis instead of ragged right-aligning the text at x=0.
    for tick_label in ax_skew.get_yticklabels():
        tick_label.set_horizontalalignment("center")
        tick_label.set_x(-0.85)
    ax_skew.set_ylim(n + .05, -.05)
    ax_skew.set_xlabel("Absolute skewness after transformation", fontsize=7)
    upper = int(np.ceil(max(2.1, data["Absolute_skewness"].max())))
    ax_skew.set_xlim(-.03, upper + .08)
    ax_skew.set_xticks(np.arange(0, upper + 1, 1))
    ax_skew.grid(axis="x", color="#E6E6E6", lw=.45)
    ax_skew.set_axisbelow(True)
    ax_skew.spines[["top", "right"]].set_visible(False)

    # Matplotlib may slightly shrink the square correlation axis to preserve
    # equal cell geometry.  Use that final, rendered box as the vertical
    # reference so the right panel is centered on—and exactly as tall as—the
    # actual left plotting area rather than its provisional axes allocation.
    fig.canvas.draw()
    corr_box = ax_corr.get_position()
    skew_box = ax_skew.get_position()
    ax_skew.set_position([skew_box.x0, corr_box.y0, skew_box.width, corr_box.height])

    stem = OUT / "05_final187_recomputed_redundancy_and_skewness_combined"
    fig.savefig(stem.with_suffix(".png"), dpi=600, facecolor="white")
    fig.savefig(stem.with_suffix(".pdf"), facecolor="white")
    fig.savefig(stem.with_suffix(".svg"), facecolor="white")
    plt.close(fig)


def plot_all22_raw_histograms_before_redundancy(
    cohort: pd.DataFrame,
    heatmap_order: list[str],
) -> None:
    """Plot all candidate features before redundancy filtering on raw scales."""
    fig, axes = plt.subplots(4, 6, figsize=(7.2, 5.2), facecolor="white")
    for panel_index, (ax, feature) in enumerate(zip(axes.flat, heatmap_order)):
        values = pd.to_numeric(cohort[feature], errors="coerce").dropna().to_numpy(float)
        ax.hist(
            values,
            bins=20,
            color="#9E9E9E",
            alpha=1.0,
            edgecolor="white",
            linewidth=.25,
        )
        max_count = max((bar.get_height() for bar in ax.patches), default=1.0)
        tick_step = max(1, int(np.ceil((max_count * 1.06) / 4.0)))
        y_top = tick_step * 4
        ax.set_ylim(0, y_top)
        ax.set_yticks(np.arange(0, y_top + tick_step, tick_step))
        ax.axvline(
            float(np.median(values)), color="#ED0000", linewidth=.55,
            linestyle=(0, (3, 2)), zorder=3,
        )
        feature_title = textwrap.fill(audit.SHORT[feature], width=23)
        title = f"{feature_title}\n({FEATURE_UNITS[feature]})"
        ax.set_title(
            title,
            loc="center",
            fontsize=7,
            fontweight="normal",
            y=1.22,
            pad=0,
            verticalalignment="top",
            linespacing=.95,
        )
        ax.set_xlabel("")
        ax.set_ylabel("Cells" if panel_index % 6 == 0 else "", fontsize=7, labelpad=1.0)
        ax.tick_params(axis="both", which="major", labelsize=6.2, width=.5, length=1.8, pad=.8)
        ax.spines[["top", "right"]].set_visible(False)
        ax.spines["left"].set_linewidth(.55)
        ax.spines["bottom"].set_linewidth(.55)
        ax.grid(axis="y", color="#E1E1E1", linewidth=.3, alpha=.8)
        ax.set_axisbelow(True)

    for ax in axes.flat[len(heatmap_order):]:
        ax.axis("off")

    fig.subplots_adjust(
        left=.055, right=.995, bottom=.055, top=.925,
        wspace=.52, hspace=.52,
    )
    stem = OUT / "06_final187_all22_raw_histograms_before_redundancy_gray20bin"
    fig.savefig(stem.with_suffix(".png"), dpi=900, facecolor="white")
    fig.savefig(stem.with_suffix(".pdf"), facecolor="white")
    fig.savefig(stem.with_suffix(".svg"), facecolor="white")
    plt.close(fig)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    audit.OUT = OUT
    audit.prior.setup_plotting()

    raw = pd.read_csv(RAW)
    final_ids = pd.read_csv(FINAL, usecols=[ID])[ID]
    cohort = final_ids.to_frame().merge(raw[[ID] + audit.ALL22], on=ID, how="left", validate="one_to_one")
    if len(cohort) != 187 or cohort[audit.ALL22].isna().any().any():
        raise RuntimeError("Expected complete all-22 measurements for the final 187 unique cells")

    logged, z_all, transform = audit.prior.transcriptomic_like_transform(cohort, audit.ALL22)
    corr = logged.corr(method="spearman")
    decisions, retained = audit.greedy_direct_filter(corr)
    pairs = []
    for i, first in enumerate(audit.ALL22):
        for second in audit.ALL22[i + 1:]:
            rho = float(corr.loc[first, second])
            if abs(rho) >= audit.RHO_THRESHOLD:
                pairs.append({"Feature_1": first, "Feature_2": second,
                              "Spearman_rho": rho, "Absolute_rho": abs(rho)})

    corr.to_csv(OUT / "01a_final187_all22_spearman_correlation_matrix.csv")
    pd.DataFrame(pairs).sort_values("Absolute_rho", ascending=False).to_csv(
        OUT / "01b_final187_redundant_pairs_abs_rho_ge0p80.csv", index=False)
    decisions["In_frozen_final10"] = decisions["Feature"].isin(FINAL10)
    decisions.to_csv(OUT / "01c_final187_redundancy_decisions.csv", index=False)
    heatmap_order = audit.plot_redundancy(corr, labeled=False)
    audit.plot_redundancy(corr, labeled=True)

    rows = []
    for feature in audit.ALL22:
        signed = float(skew(logged[feature].to_numpy(float), bias=False))
        rows.append({
            "Feature": feature, "Short_name": audit.SHORT[feature],
            "Signed_skewness": signed, "Absolute_skewness": abs(signed),
            "Status": "highly_skewed" if abs(signed) >= 2 else
                      ("moderately_skewed" if abs(signed) >= 1 else "acceptable"),
            "Retained_after_abs_rho_filter": feature in retained,
            "In_frozen_final10": feature in FINAL10,
        })
    skew_all = pd.DataFrame(rows)
    skew_all["Heatmap_order"] = skew_all["Feature"].map(
        {feature: i + 1 for i, feature in enumerate(heatmap_order)})
    skew_all = skew_all.sort_values("Heatmap_order")
    skew_all.to_csv(OUT / "02a_final187_all22_skewness.csv", index=False)
    audit.plot_abs_skew(skew_all, heatmap_order, "02_final187_all22_absolute_skewness", 6.0)

    skew_retained = skew_all[skew_all["Retained_after_abs_rho_filter"]].copy()
    audit.plot_abs_skew(skew_retained, heatmap_order,
                        "03_final187_redundancy_filtered_absolute_skewness", 4.15)
    skew_retained.sort_values("Absolute_skewness", ascending=False).to_csv(
        OUT / "03a_final187_redundancy_filtered_skewness.csv", index=False)

    processed = z_all[retained].copy()
    processed.insert(0, ID, cohort[ID].to_numpy())
    processed.to_csv(OUT / "04_final187_redundancy_filtered_global_z_matrix.csv", index=False)
    transform["Retained_after_abs_rho_filter"] = transform["Feature"].isin(retained)
    transform["In_frozen_final10"] = transform["Feature"].isin(FINAL10)
    transform.to_csv(OUT / "04a_final187_preprocessing_parameters.csv", index=False)

    plot_aligned_redundancy_and_skewness(corr, skew_all, heatmap_order)
    plot_all22_raw_histograms_before_redundancy(cohort, heatmap_order)

    summary = {
        "cohort_n": 187,
        "candidate_features_n": 22,
        "preprocessing": "minimum translation -> feature-sum normalization x10000 -> log1p -> global feature Z score",
        "redundancy_metric": "Spearman rho on log1p values",
        "redundancy_threshold_abs_rho": 0.80,
        "redundant_pairs_n": len(pairs),
        "retained_after_direct_filter_n": len(retained),
        "retained_after_direct_filter": retained,
        "features_abs_skew_ge_2": skew_all.loc[skew_all.Absolute_skewness.ge(2), "Feature"].tolist(),
        "frozen_final10_all_retained_by_recomputed_filter": set(FINAL10).issubset(retained),
        "frozen_final10": FINAL10,
        "batch_correction": False,
    }
    (OUT / "00_final187_QC_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps(summary, indent=2, ensure_ascii=False))
    print(OUT)


if __name__ == "__main__":
    main()
