#!/usr/bin/env python3
"""Add Dunn-Holm pairwise significance brackets to all-25 paired plots."""

from __future__ import annotations

import json
from math import erfc, exp, sqrt
from pathlib import Path
import textwrap

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import Normalize
from matplotlib.ticker import MaxNLocator
import numpy as np
import pandas as pd

from plot_all25_violin_plus_feature_tsne import (
    E_COLORS,
    E_LEVELS,
    OUT_DIR as BASE_OUT_DIR,
    ZSCORE_CMAP,
    ZSCORE_LIMIT,
    add_tsne_panel,
    configure_style,
    load_data,
)


OUT_DIR = BASE_OUT_DIR / "posthoc_pairwise_significance"


def average_ranks(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    order = np.argsort(values, kind="mergesort")
    sorted_values = values[order]
    sorted_ranks = np.empty(len(values), dtype=float)
    start = 0
    while start < len(values):
        end = start + 1
        while end < len(values) and sorted_values[end] == sorted_values[start]:
            end += 1
        sorted_ranks[start:end] = ((start + 1) + end) / 2.0
        start = end
    ranks = np.empty(len(values), dtype=float)
    ranks[order] = sorted_ranks
    return ranks


def holm_adjust(p_values: np.ndarray) -> np.ndarray:
    p_values = np.asarray(p_values, dtype=float)
    order = np.argsort(p_values)
    adjusted_sorted = np.empty(len(p_values), dtype=float)
    running_max = 0.0
    total = len(p_values)
    for sorted_index, original_index in enumerate(order):
        candidate = min(1.0, (total - sorted_index) * p_values[original_index])
        running_max = max(running_max, candidate)
        adjusted_sorted[sorted_index] = running_max
    adjusted = np.empty(len(p_values), dtype=float)
    adjusted[order] = adjusted_sorted
    return adjusted


def bh_adjust(p_values: np.ndarray) -> np.ndarray:
    p_values = np.asarray(p_values, dtype=float)
    order = np.argsort(p_values)[::-1]
    total = len(p_values)
    adjusted = np.empty(total, dtype=float)
    running_min = 1.0
    for reverse_rank, original_index in enumerate(order, start=1):
        rank = total - reverse_rank + 1
        candidate = min(1.0, p_values[original_index] * total / rank)
        running_min = min(running_min, candidate)
        adjusted[original_index] = running_min
    return adjusted


def kruskal_wallis(values: np.ndarray, groups: np.ndarray) -> tuple[float, float]:
    values = np.asarray(values, dtype=float)
    groups = np.asarray(groups, dtype=str)
    ranks = average_ranks(values)
    total = len(values)
    rank_term = 0.0
    for group in E_LEVELS:
        mask = groups == group
        rank_term += ranks[mask].sum() ** 2 / mask.sum()
    statistic = 12.0 * rank_term / (total * (total + 1.0)) - 3.0 * (total + 1.0)

    _, tie_counts = np.unique(values, return_counts=True)
    correction = 1.0 - np.sum(tie_counts**3 - tie_counts) / (total**3 - total)
    if correction > 0:
        statistic /= correction

    # Five groups give 4 degrees of freedom. For df=4, the chi-square survival
    # function has the closed form exp(-x/2) * (1 + x/2).
    p_value = exp(-statistic / 2.0) * (1.0 + statistic / 2.0)
    return float(statistic), float(min(max(p_value, 0.0), 1.0))


def dunn_holm(values: np.ndarray, groups: np.ndarray) -> pd.DataFrame:
    values = np.asarray(values, dtype=float)
    groups = np.asarray(groups, dtype=str)
    ranks = average_ranks(values)
    total = len(values)
    _, tie_counts = np.unique(values, return_counts=True)
    tie_term = np.sum(tie_counts**3 - tie_counts)
    variance_constant = (
        total * (total + 1.0) / 12.0
        - tie_term / (12.0 * (total - 1.0))
    )

    group_rank_means = {
        group: float(ranks[groups == group].mean()) for group in E_LEVELS
    }
    group_sizes = {group: int(np.sum(groups == group)) for group in E_LEVELS}
    rows: list[dict[str, float | int | str]] = []
    for left_index in range(len(E_LEVELS) - 1):
        for right_index in range(left_index + 1, len(E_LEVELS)):
            group_1 = E_LEVELS[left_index]
            group_2 = E_LEVELS[right_index]
            standard_error = sqrt(
                variance_constant
                * (1.0 / group_sizes[group_1] + 1.0 / group_sizes[group_2])
            )
            z_value = (
                group_rank_means[group_1] - group_rank_means[group_2]
            ) / standard_error
            p_raw = erfc(abs(z_value) / sqrt(2.0))
            rows.append(
                {
                    "Group_1": group_1,
                    "Group_2": group_2,
                    "N_1": group_sizes[group_1],
                    "N_2": group_sizes[group_2],
                    "Dunn_Z": z_value,
                    "Dunn_r": abs(z_value) / sqrt(total),
                    "P_raw": p_raw,
                }
            )
    result = pd.DataFrame(rows)
    result["P_adj_Holm_within_feature"] = holm_adjust(result["P_raw"].to_numpy())
    return result


def significance_label(p_value: float) -> str:
    if p_value < 0.0001:
        return "****"
    if p_value < 0.001:
        return "***"
    if p_value < 0.01:
        return "**"
    if p_value < 0.05:
        return "*"
    return "ns"


def assign_bracket_levels(pairs: pd.DataFrame) -> pd.DataFrame:
    if pairs.empty:
        result = pairs.copy()
        result["Level"] = pd.Series(dtype=int)
        return result
    result = pairs.copy()
    result["X_1"] = result["Group_1"].map({name: index + 1 for index, name in enumerate(E_LEVELS)})
    result["X_2"] = result["Group_2"].map({name: index + 1 for index, name in enumerate(E_LEVELS)})
    result["Span"] = result["X_2"] - result["X_1"]
    processing_order = result.sort_values(
        ["Span", "X_1", "X_2", "P_adj_Holm_within_feature"]
    ).index
    intervals_by_level: list[list[tuple[int, int]]] = []
    levels = pd.Series(index=result.index, dtype=int)
    for row_index in processing_order:
        left = int(result.at[row_index, "X_1"])
        right = int(result.at[row_index, "X_2"])
        candidate = 0
        while True:
            if candidate == len(intervals_by_level):
                intervals_by_level.append([])
            overlaps = any(
                left <= existing_right and right >= existing_left
                for existing_left, existing_right in intervals_by_level[candidate]
            )
            if not overlaps:
                intervals_by_level[candidate].append((left, right))
                levels.at[row_index] = candidate + 1
                break
            candidate += 1
    result["Level"] = levels.astype(int)
    return result


def compact_title(title: str, unit: str) -> str:
    display_unit = {"MOhm": "MΩ", "AP count": "count", "mV s": "mV·s"}.get(unit, unit)
    if display_unit in {"Ratio", "Index", "CV", "Fano factor"}:
        return title
    return f"{title} ({display_unit})"


def add_violin_with_posthoc(
    ax: mpl.axes.Axes,
    data: pd.DataFrame,
    feature_row: pd.Series,
    pairwise: pd.DataFrame,
    display_threshold: float,
    seed: int,
) -> int:
    feature = feature_row["Feature"]
    grouped = [
        data.loc[data["E_type"] == e_type, feature].to_numpy(dtype=float)
        for e_type in E_LEVELS
    ]
    positions = np.arange(1, 6)
    violins = ax.violinplot(
        grouped,
        positions=positions,
        widths=0.78,
        showmeans=False,
        showmedians=False,
        showextrema=False,
        bw_method="scott",
    )
    for body, e_type in zip(violins["bodies"], E_LEVELS):
        body.set_facecolor(E_COLORS[e_type])
        body.set_edgecolor("#333333")
        body.set_linewidth(0.35)
        body.set_alpha(0.82)

    rng = np.random.default_rng(seed)
    for x_pos, values in zip(positions, grouped):
        jitter = rng.uniform(-0.105, 0.105, size=len(values))
        ax.scatter(
            np.full(len(values), x_pos) + jitter,
            values,
            s=0.75,
            c="#111111",
            alpha=0.28,
            linewidths=0,
            zorder=3,
        )
        q1, median, q3 = np.quantile(values, [0.25, 0.50, 0.75])
        ax.vlines(x_pos, q1, q3, color="#222222", linewidth=0.75, zorder=4)
        ax.scatter(
            [x_pos], [median], s=4.0, facecolor="white", edgecolor="#222222",
            linewidth=0.45, zorder=5
        )

    all_values = data[feature].to_numpy(dtype=float)
    data_min = float(all_values.min())
    data_max = float(all_values.max())
    data_span = data_max - data_min
    if not np.isfinite(data_span) or data_span <= 0:
        data_span = max(abs(data_min), abs(data_max), 1.0)

    visible_pairs = pairwise.loc[
        pairwise["P_adj_Holm_within_feature"] < display_threshold
    ].copy()
    visible_pairs = assign_bracket_levels(visible_pairs)
    lower_limit = data_min - 0.055 * data_span
    if visible_pairs.empty:
        upper_limit = data_max + 0.065 * data_span
    else:
        bracket_base = data_max + 0.045 * data_span
        bracket_step = 0.075 * data_span
        bracket_tip = 0.018 * data_span
        for row in visible_pairs.itertuples(index=False):
            bracket_y = bracket_base + (row.Level - 1) * bracket_step
            ax.plot(
                [row.X_1, row.X_1, row.X_2, row.X_2],
                [bracket_y - bracket_tip, bracket_y, bracket_y, bracket_y - bracket_tip],
                color="black",
                linewidth=0.35,
                solid_capstyle="butt",
                clip_on=False,
                zorder=6,
            )
            ax.text(
                (row.X_1 + row.X_2) / 2.0,
                bracket_y + 0.006 * data_span,
                significance_label(row.P_adj_Holm_within_feature),
                ha="center",
                va="bottom",
                fontsize=3.8,
                color="black",
                clip_on=False,
                zorder=7,
            )
        upper_limit = (
            bracket_base
            + (int(visible_pairs["Level"].max()) - 1) * bracket_step
            + 0.075 * data_span
        )

    ax.set_xlim(0.55, 5.45)
    ax.set_ylim(lower_limit, upper_limit)
    ax.set_xticks(positions, E_LEVELS)
    ax.set_ylabel(None)
    ax.set_title(
        textwrap.fill(
            compact_title(feature_row["Short_title"], feature_row["Y_label"]),
            width=30,
        ),
        fontweight="bold",
        fontsize=5.5,
        pad=1.5,
    )
    ax.yaxis.set_major_locator(MaxNLocator(nbins=3))
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_linewidth(0.5)
    ax.spines["bottom"].set_linewidth(0.5)
    ax.tick_params(direction="out", length=2.0, width=0.5, pad=1.2)
    ax.set_box_aspect(1.25)
    return len(visible_pairs)


def make_figure(
    data: pd.DataFrame,
    feature_row: pd.Series,
    pairwise: pd.DataFrame,
    threshold: float,
    feature_index: int,
) -> tuple[mpl.figure.Figure, int]:
    fig = plt.figure(figsize=(2.4, 1.4))
    grid = fig.add_gridspec(
        1, 2, width_ratios=[0.82, 1.18], left=0.105, right=0.985,
        bottom=0.16, top=0.94, wspace=0.20
    )
    violin_ax = fig.add_subplot(grid[0, 0])
    tsne_ax = fig.add_subplot(grid[0, 1])
    count = add_violin_with_posthoc(
        violin_ax,
        data,
        feature_row,
        pairwise,
        display_threshold=threshold,
        seed=777 + feature_index,
    )
    add_tsne_panel(tsne_ax, data, feature_row["Feature"])
    return fig, count


def calculate_statistics(
    data: pd.DataFrame, fmap: pd.DataFrame
) -> tuple[pd.DataFrame, pd.DataFrame]:
    omnibus_rows: list[dict[str, object]] = []
    pairwise_frames: list[pd.DataFrame] = []
    groups = data["E_type"].astype(str).to_numpy()
    for feature_index, feature_row in fmap.iterrows():
        feature = feature_row["Feature"]
        values = data[feature].to_numpy(dtype=float)
        statistic, p_value = kruskal_wallis(values, groups)
        omnibus_rows.append(
            {
                "Order": feature_index + 1,
                "Feature": feature,
                "Short_title": feature_row["Short_title"],
                "N": len(values),
                "Kruskal_Wallis_H": statistic,
                "df": 4,
                "P_raw": p_value,
            }
        )
        feature_pairs = dunn_holm(values, groups)
        feature_pairs.insert(0, "Short_title", feature_row["Short_title"])
        feature_pairs.insert(0, "Feature", feature)
        feature_pairs.insert(0, "Order", feature_index + 1)
        pairwise_frames.append(feature_pairs)

    omnibus = pd.DataFrame(omnibus_rows)
    omnibus["P_adj_BH_across_25_features"] = bh_adjust(
        omnibus["P_raw"].to_numpy()
    )
    pairwise = pd.concat(pairwise_frames, ignore_index=True)
    pairwise["P_adj_BH_global_250_tests"] = bh_adjust(
        pairwise["P_raw"].to_numpy()
    )
    pairwise["Significance_Holm"] = pairwise[
        "P_adj_Holm_within_feature"
    ].map(significance_label)
    return omnibus, pairwise


def generate_version(
    data: pd.DataFrame,
    fmap: pd.DataFrame,
    omnibus: pd.DataFrame,
    pairwise: pd.DataFrame,
    threshold: float,
) -> dict[str, object]:
    threshold_tag = f"p{threshold:g}".replace(".", "p")
    version_dir = OUT_DIR / f"DunnHolm_padj_lt_{threshold_tag}"
    png_dir = version_dir / "individual_png"
    pdf_dir = version_dir / "individual_pdf"
    png_dir.mkdir(parents=True, exist_ok=True)
    pdf_dir.mkdir(parents=True, exist_ok=True)
    multipage_pdf = version_dir / (
        f"All25_Violin_FeaturetSNE_DunnHolm_padj_lt_{threshold_tag}_25pages.pdf"
    )
    manifest_rows: list[dict[str, object]] = []

    with PdfPages(multipage_pdf) as pdf_pages:
        for feature_index, feature_row in fmap.iterrows():
            feature_pairs = pairwise.loc[
                pairwise["Feature"] == feature_row["Feature"]
            ]
            fig, visible_count = make_figure(
                data, feature_row, feature_pairs, threshold, feature_index + 1
            )
            stem = (
                f"{feature_index + 1:02d}_{feature_row['File_stub']}_"
                f"DunnHolm_padj_lt_{threshold_tag}_W2p4_H1p4"
            )
            png_file = png_dir / f"{stem}.png"
            pdf_file = pdf_dir / f"{stem}.pdf"
            fig.savefig(png_file, dpi=600, facecolor="white")
            fig.savefig(pdf_file, facecolor="white")
            pdf_pages.savefig(fig, facecolor="white")
            plt.close(fig)
            manifest_rows.append(
                {
                    "Order": feature_index + 1,
                    "Feature": feature_row["Feature"],
                    "Short_title": feature_row["Short_title"],
                    "Display_threshold": threshold,
                    "Significant_pairs_displayed": visible_count,
                    "PNG": png_file.as_posix(),
                    "PDF": pdf_file.as_posix(),
                }
            )

    overview = plt.figure(figsize=(7.2, 11.2))
    overview_grid = overview.add_gridspec(
        9, 3, left=0.035, right=0.995, bottom=0.02, top=0.995,
        wspace=0.08, hspace=0.13
    )
    for feature_index, feature_row in fmap.iterrows():
        outer = overview_grid[feature_index // 3, feature_index % 3].subgridspec(
            1, 2, width_ratios=[0.82, 1.18], wspace=0.16
        )
        left_ax = overview.add_subplot(outer[0, 0])
        right_ax = overview.add_subplot(outer[0, 1])
        feature_pairs = pairwise.loc[
            pairwise["Feature"] == feature_row["Feature"]
        ]
        add_violin_with_posthoc(
            left_ax,
            data,
            feature_row,
            feature_pairs,
            display_threshold=threshold,
            seed=777 + feature_index + 1,
        )
        add_tsne_panel(right_ax, data, feature_row["Feature"])
    for blank_index in range(len(fmap), 27):
        blank_ax = overview.add_subplot(
            overview_grid[blank_index // 3, blank_index % 3]
        )
        blank_ax.set_axis_off()

    overview_png = version_dir / (
        f"All25_Overview_DunnHolm_padj_lt_{threshold_tag}_W7p2_H11p2.png"
    )
    overview_pdf = version_dir / (
        f"All25_Overview_DunnHolm_padj_lt_{threshold_tag}_W7p2_H11p2.pdf"
    )
    overview.savefig(overview_png, dpi=400, facecolor="white")
    overview.savefig(overview_pdf, facecolor="white")
    plt.close(overview)

    manifest = pd.DataFrame(manifest_rows)
    manifest.to_csv(version_dir / "manifest.csv", index=False)
    return {
        "threshold": threshold,
        "version_dir": version_dir.as_posix(),
        "multipage_pdf": multipage_pdf.as_posix(),
        "overview_png": overview_png.as_posix(),
        "overview_pdf": overview_pdf.as_posix(),
        "displayed_pair_total": int(manifest["Significant_pairs_displayed"].sum()),
    }


def make_colorbar() -> None:
    fig = plt.figure(figsize=(1.2, 0.28))
    ax = fig.add_axes([0.12, 0.48, 0.76, 0.20])
    scalar_mappable = mpl.cm.ScalarMappable(
        norm=Normalize(-ZSCORE_LIMIT, ZSCORE_LIMIT), cmap=ZSCORE_CMAP
    )
    colorbar = fig.colorbar(scalar_mappable, cax=ax, orientation="horizontal")
    colorbar.set_ticks([-ZSCORE_LIMIT, 0, ZSCORE_LIMIT])
    colorbar.set_ticklabels([f"−{ZSCORE_LIMIT:g}", "0", f"{ZSCORE_LIMIT:g}"])
    colorbar.ax.tick_params(length=1.6, width=0.4, labelsize=5, pad=1)
    colorbar.outline.set_linewidth(0.4)
    fig.text(0.5, 0.94, "Feature Z-score", ha="center", va="top", fontsize=5.5)
    fig.savefig(OUT_DIR / "Shared_Feature_Zscore_Colorbar.png", dpi=600, transparent=True)
    fig.savefig(OUT_DIR / "Shared_Feature_Zscore_Colorbar.pdf", transparent=True)
    plt.close(fig)


def main() -> None:
    configure_style()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    data, fmap = load_data()
    omnibus, pairwise = calculate_statistics(data, fmap)
    omnibus.to_csv(OUT_DIR / "KruskalWallis_BH25_Omnibus.csv", index=False)
    pairwise.to_csv(OUT_DIR / "Dunn_Pairwise_Holm_within_feature.csv", index=False)
    make_colorbar()

    versions = [
        generate_version(data, fmap, omnibus, pairwise, threshold=0.05),
        generate_version(data, fmap, omnibus, pairwise, threshold=0.01),
    ]
    run_log = {
        "cells": int(len(data)),
        "cohort": "active GC×HC×EC triple-consensus cells",
        "features": int(len(fmap)),
        "omnibus_test": "two-sided Kruskal-Wallis; BH correction across 25 features",
        "posthoc_test": "two-sided Dunn test; Holm correction across 10 pairs within each feature",
        "star_key": {
            "*": "Holm-adjusted P < 0.05",
            "**": "Holm-adjusted P < 0.01",
            "***": "Holm-adjusted P < 0.001",
            "****": "Holm-adjusted P < 0.0001",
        },
        "versions": versions,
    }
    (OUT_DIR / "run_log.json").write_text(
        json.dumps(run_log, indent=2), encoding="utf-8"
    )
    print(OUT_DIR)
    for version in versions:
        print(
            f"threshold={version['threshold']}: "
            f"{version['displayed_pair_total']} displayed feature-pairs"
        )


if __name__ == "__main__":
    main()
