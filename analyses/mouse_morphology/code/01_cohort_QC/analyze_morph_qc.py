#!/usr/bin/env python3
"""Morphology missingness, cell QC, skewness, and redundancy audit."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import leaves_list, linkage
from scipy.spatial.distance import squareform


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / ".codex-work" / "e_type_qc" / "stage1_filtered.json"
EXCLUSIONS = ROOT / "config" / "e_type_global_exclusions.txt"
OUT = ROOT / "outputs" / "morph_qc" / "initial_qc_skew_redundancy"

ID_FIELDS = ["M_ID", "M_Neuron_id"]
COORD_FIELDS = ["M_center_X", "M_center_Y", "M_center_Z"]
NOT_APPLICABLE_OR_EMPTY = [
    "M_apical_dendrite_longthest_path",
    "M_axon_longthest_path",
]
SPARSE_FIELDS = ["M_trunk_angle_dispersion_index"]
ANALYSIS_EXCLUSIONS = set(ID_FIELDS + COORD_FIELDS + NOT_APPLICABLE_OR_EMPTY + SPARSE_FIELDS)


def load_source() -> tuple[pd.DataFrame, list[str]]:
    payload = json.loads(SOURCE.read_text(encoding="utf-8"))
    df = pd.DataFrame(payload["analysis_rows"], columns=payload["analysis_headers"])
    excluded = [
        line.strip()
        for line in EXCLUSIONS.read_text(encoding="utf-8").splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ]
    df = df.loc[~df["MSN_unique_ID"].isin(excluded)].copy()
    return df, excluded


def numeric(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series.replace({"None": np.nan, "NA": np.nan, "": np.nan}), errors="coerce")


def sample_skew(x: pd.Series) -> float:
    a = x.dropna().to_numpy(float)
    n = len(a)
    if n < 3:
        return np.nan
    sd = a.std(ddof=1)
    if sd == 0:
        return 0.0
    return float(n / ((n - 1) * (n - 2)) * np.sum(((a - a.mean()) / sd) ** 3))


def skew_status(value: float) -> str:
    a = abs(value)
    if a > 1:
        return "高度偏态"
    if a >= 0.5:
        return "中度偏态"
    return "近似对称"


def robust_outlier_flags(data: pd.DataFrame) -> pd.DataFrame:
    transformed = data.copy()
    for col in transformed:
        x = transformed[col]
        if x.min() >= 0 and sample_skew(x) > 1:
            transformed[col] = np.log1p(x)
    med = transformed.median()
    mad = (transformed - med).abs().median()
    robust_z = 0.67448975 * (transformed - med) / mad.replace(0, np.nan)
    return robust_z.abs() > 5


def plot_missingness(feature_df: pd.DataFrame, cell_df: pd.DataFrame, eligible: pd.DataFrame) -> None:
    mpl.rcParams.update({"font.family": "Arial", "font.size": 7, "axes.linewidth": 0.6})
    fig = plt.figure(figsize=(8.2, 8.4))
    gs = fig.add_gridspec(2, 2, height_ratios=[3.2, 1], width_ratios=[1.45, 1], hspace=0.32, wspace=0.35)

    ax = fig.add_subplot(gs[0, :])
    order = feature_df.sort_values("Missing_rate_among_morph_cells", ascending=True)
    colors = np.where(order["Recommended_role"].eq("Use"), "#4C78A8", "#B9BEC5")
    ax.barh(np.arange(len(order)), order["Missing_rate_among_morph_cells"] * 100, color=colors, height=0.72)
    ax.set_yticks(np.arange(len(order)), order["Feature"], fontsize=5.8)
    ax.set_xlabel("Missing among morphology-positive cells (%)")
    ax.set_xlim(0, 105)
    ax.axvline(20, color="#D62728", linestyle="--", linewidth=0.7)
    ax.grid(axis="x", color="#E6E6E6", linewidth=0.5)
    ax.set_axisbelow(True)
    ax.set_title("Morphology feature missingness", loc="left", fontweight="bold")

    ax2 = fig.add_subplot(gs[1, 0])
    vals = cell_df["Missing_rate_candidate22"].to_numpy() * 100
    ax2.hist(vals, bins=np.linspace(0, 100, 21), color="#4C78A8", edgecolor="white", linewidth=0.4)
    ax2.set_xlabel("Per-cell missingness across candidate 22 features (%)")
    ax2.set_ylabel("Cells")
    ax2.set_title("Cell-level missingness", loc="left", fontweight="bold")

    ax3 = fig.add_subplot(gs[1, 1])
    batch = eligible["T_Batch"].value_counts().sort_index()
    ax3.bar(batch.index.astype(str), batch.values, color="#72B7B2")
    ax3.tick_params(axis="x", rotation=45)
    ax3.set_ylabel("Morphology-positive cells")
    ax3.set_title("Eligible cells by batch", loc="left", fontweight="bold")
    fig.savefig(OUT / "Morph_missingness_QC_overview.png", dpi=400, bbox_inches="tight", facecolor="white")
    fig.savefig(OUT / "Morph_missingness_QC_overview.pdf", bbox_inches="tight", facecolor="white")
    plt.close(fig)


def plot_skewness(skew_df: pd.DataFrame) -> None:
    ordered = skew_df.sort_values("Skewness")
    colors = ["#B2182B" if abs(v) > 1 else "#E08214" if abs(v) >= 0.5 else "#7F7F7F" for v in ordered["Skewness"]]
    fig, ax = plt.subplots(figsize=(7.2, 6.4))
    ax.barh(np.arange(len(ordered)), ordered["Skewness"], color=colors, height=0.7)
    ax.set_yticks(np.arange(len(ordered)), ordered["Feature"], fontsize=6.2)
    ax.axvline(-1, color="#B2182B", linestyle="--", linewidth=0.7)
    ax.axvline(1, color="#B2182B", linestyle="--", linewidth=0.7)
    ax.axvline(-0.5, color="#E08214", linestyle=":", linewidth=0.7)
    ax.axvline(0.5, color="#E08214", linestyle=":", linewidth=0.7)
    ax.axvline(0, color="#444444", linewidth=0.6)
    ax.set_xlabel("Bias-corrected Fisher-Pearson skewness")
    ax.set_title("Raw morphology feature skewness", loc="left", fontweight="bold")
    ax.grid(axis="x", color="#E6E6E6", linewidth=0.5)
    ax.set_axisbelow(True)
    fig.savefig(OUT / "Morph_raw_feature_skewness.png", dpi=400, bbox_inches="tight", facecolor="white")
    fig.savefig(OUT / "Morph_raw_feature_skewness.pdf", bbox_inches="tight", facecolor="white")
    plt.close(fig)


def plot_correlation(corr: pd.DataFrame, pairs: pd.DataFrame) -> None:
    dist = 1 - corr.abs()
    np.fill_diagonal(dist.values, 0)
    order = leaves_list(linkage(squareform(dist.values, checks=False), method="average"))
    corr = corr.iloc[order, order]
    labels = [c.replace("M_", "").replace("_", " ") for c in corr.columns]
    n = len(labels)
    counts = pd.concat([pairs["Var1"], pairs["Var2"]]).value_counts().reindex(corr.columns, fill_value=0)

    fig = plt.figure(figsize=(10.2, 10.8))
    gs = fig.add_gridspec(2, 1, height_ratios=[4.5, 1.25], hspace=0.12)
    ax = fig.add_subplot(gs[0])
    cmap = mpl.colors.LinearSegmentedColormap.from_list("corr", ["#2166AC", "#F7F7F7", "#B2182B"])
    norm = mpl.colors.Normalize(-1, 1)
    for i in range(n):
        for j in range(n):
            if i < j:
                continue
            r = corr.iat[i, j]
            ax.add_patch(plt.Rectangle((j, i), 1, 1, facecolor="white", edgecolor="#D0D0D0", linewidth=0.35))
            if i == j:
                ax.text(j + 0.5, i + 0.5, "1.00", ha="center", va="center", fontsize=5.2)
            else:
                angle = 45 if r >= 0 else -45
                e = Ellipse((j + 0.5, i + 0.5), width=0.88, height=max(0.08, 0.88 * (1 - abs(r))), angle=angle,
                            facecolor=cmap(norm(r)), edgecolor="none", alpha=0.92)
                ax.add_patch(e)
                ax.text(j + 0.5, i + 0.5, f"{r:.2f}", ha="center", va="center", fontsize=4.4, color="#202020")
    ax.set_xlim(0, n); ax.set_ylim(n, 0); ax.set_aspect("equal")
    ax.set_xticks(np.arange(n) + 0.5, labels, rotation=60, ha="right", fontsize=5.5)
    ax.set_yticks(np.arange(n) + 0.5, labels, fontsize=5.5)
    ax.tick_params(length=0)
    ax.set_title("Morphology redundancy: Spearman correlation (|rho| > 0.85 flagged)", loc="left", fontweight="bold")
    sm = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
    fig.colorbar(sm, ax=ax, fraction=0.025, pad=0.02, label="Spearman rho")

    ax2 = fig.add_subplot(gs[1])
    bar_colors = ["#B2182B" if x else "#B9BEC5" for x in counts]
    ax2.bar(np.arange(n), counts.values, color=bar_colors, width=0.72)
    ax2.set_xticks(np.arange(n), labels, rotation=60, ha="right", fontsize=5.5)
    ax2.set_ylabel("Redundant pairs\n(|rho| > 0.85)")
    ax2.grid(axis="y", color="#E6E6E6", linewidth=0.5)
    ax2.set_axisbelow(True)
    fig.savefig(OUT / "Morph_Spearman_redundancy_matrix.png", dpi=400, bbox_inches="tight", facecolor="white")
    fig.savefig(OUT / "Morph_Spearman_redundancy_matrix.pdf", bbox_inches="tight", facecolor="white")
    plt.close(fig)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    df, exclusions = load_source()
    morph_cols = [c for c in df.columns if c.startswith("M_")]
    morph = df[morph_cols].copy()
    for col in morph_cols:
        if col != "M_Neuron_id":
            morph[col] = numeric(morph[col])
    has_m = df["has_M"].astype(str).str.lower().eq("true")
    morph_cells = df.loc[has_m].copy()
    candidate = [c for c in morph_cols if c not in ANALYSIS_EXCLUSIONS]
    x = morph.loc[has_m, candidate].apply(numeric)

    feature_rows = []
    for col in morph_cols:
        all_missing = morph[col].isna().mean()
        within_missing = morph.loc[has_m, col].isna().mean()
        role = "Use" if col in candidate else "Exclude from clustering"
        reason = "core quantitative morphology feature"
        if col in ID_FIELDS: reason = "identifier, not a phenotype"
        elif col in COORD_FIELDS: reason = "acquisition/location coordinate; M_center_Z is constant"
        elif col in NOT_APPLICABLE_OR_EMPTY: reason = "entirely missing; axon/apical field is not usable here"
        elif col in SPARSE_FIELDS: reason = "only six observed values; 97.3% missing among morphology cells"
        feature_rows.append({
            "Feature": col,
            "Nonmissing_all_MSN": int(morph[col].notna().sum()),
            "Missing_rate_all_MSN": float(all_missing),
            "Nonmissing_among_morph_cells": int(morph.loc[has_m, col].notna().sum()),
            "Missing_rate_among_morph_cells": float(within_missing),
            "Unique_nonmissing": int(morph[col].nunique(dropna=True)),
            "Recommended_role": role,
            "Reason": reason,
        })
    feature_df = pd.DataFrame(feature_rows)

    cell_df = df[["MSN_unique_ID", "METL_source_row", "T_Batch", "completeness", "has_M"]].copy()
    cell_df["Missing_n_candidate22"] = morph[candidate].isna().sum(axis=1)
    cell_df["Missing_rate_candidate22"] = cell_df["Missing_n_candidate22"] / len(candidate)
    cell_df["Morph_eligible_core22"] = has_m & cell_df["Missing_n_candidate22"].eq(0)

    skew_df = pd.DataFrame({
        "Feature": candidate,
        "N": [int(x[c].notna().sum()) for c in candidate],
        "Skewness": [sample_skew(x[c]) for c in candidate],
    })
    skew_df["Abs_skewness"] = skew_df["Skewness"].abs()
    skew_df["Status"] = skew_df["Skewness"].map(skew_status)
    skew_df["Suggested_transform"] = np.where(skew_df["Skewness"] > 1, "log1p then z-score", "z-score; Yeo-Johnson sensitivity check")
    skew_df = skew_df.sort_values("Abs_skewness", ascending=False).reset_index(drop=True)

    corr = x.corr(method="spearman")
    pair_rows = []
    for i, a in enumerate(candidate):
        for b in candidate[:i]:
            r = float(corr.loc[a, b])
            if abs(r) > 0.85:
                pair_rows.append({"Var1": a, "Var2": b, "Correlation": r, "Abs_correlation": abs(r)})
    pairs = pd.DataFrame(pair_rows).sort_values("Abs_correlation", ascending=False).reset_index(drop=True)

    logical = pd.DataFrame(index=x.index)
    logical["radius_order_invalid"] = ~((x["M_soma_minimal_radius"] <= x["M_soma_average_radius"]) & (x["M_soma_average_radius"] <= x["M_soma_maximal_radius"]))
    logical["trunk_angle_order_invalid"] = x["M_trunk_angle_min"] > x["M_trunk_angle_max"]
    positive = [c for c in candidate if c not in {"M_trunk_angle_min", "M_trunk_angle_max"}]
    logical["nonpositive_size_or_count"] = (x[positive] <= 0).any(axis=1)
    logical["circularity_outside_0_1"] = ~x["M_soma_circularity_index"].between(0, 1, inclusive="both")
    outlier_matrix = robust_outlier_flags(x)
    cell_qc = cell_df.loc[has_m, ["MSN_unique_ID", "METL_source_row", "T_Batch", "completeness"]].copy()
    cell_qc = cell_qc.join(logical)
    cell_qc["Robust_outlier_feature_count_abs_z_gt5"] = outlier_matrix.sum(axis=1)
    cell_qc["Hard_QC_fail"] = logical.any(axis=1)
    cell_qc["Review_flag"] = cell_qc["Hard_QC_fail"] | cell_qc["Robust_outlier_feature_count_abs_z_gt5"].gt(0)
    cell_qc["Recommended_action"] = np.where(cell_qc["Hard_QC_fail"], "Exclude after source-image verification", np.where(cell_qc["Review_flag"], "Manual ASC/image review; retain unless technical defect confirmed", "Retain"))

    plot_missingness(feature_df, cell_df, morph_cells)
    plot_skewness(skew_df)
    plot_correlation(corr, pairs)

    results = {
        "source_rows_after_global_exclusions": len(df),
        "global_exclusions": exclusions,
        "morphology_positive_cells": int(has_m.sum()),
        "candidate_features": candidate,
        "candidate_feature_count": len(candidate),
        "complete_candidate_cells": int(cell_df["Morph_eligible_core22"].sum()),
        "hard_qc_fail_cells": int(cell_qc["Hard_QC_fail"].sum()),
        "review_flag_cells": int(cell_qc["Review_flag"].sum()),
        "highly_skewed_features_abs_gt1": int((skew_df["Abs_skewness"] > 1).sum()),
        "redundant_pairs_abs_spearman_gt0_85": len(pairs),
        "method_notes": {
            "skewness": "bias-corrected Fisher-Pearson sample skewness on raw values",
            "redundancy": "pairwise Spearman correlation; flagged when absolute rho > 0.85",
            "outlier": "log1p for nonnegative features with skewness > 1, then median/MAD robust z; abs(z)>5 is review-only",
        },
    }
    (OUT / "Morph_QC_results.json").write_text(json.dumps(results, indent=2), encoding="utf-8")
    feature_df.to_json(OUT / "feature_missingness.json", orient="records", indent=2)
    cell_df.to_json(OUT / "cell_missingness.json", orient="records", indent=2)
    skew_df.to_json(OUT / "skewness.json", orient="records", indent=2)
    pairs.to_json(OUT / "redundant_pairs.json", orient="records", indent=2)
    corr.reset_index(names="Feature").to_json(OUT / "spearman_correlation.json", orient="records", indent=2)
    cell_qc.to_json(OUT / "cell_qc_flags.json", orient="records", indent=2)
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
