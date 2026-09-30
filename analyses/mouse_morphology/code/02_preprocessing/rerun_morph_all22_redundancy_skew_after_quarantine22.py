#!/usr/bin/env python
"""All-22-feature redundancy and skewness audit after 22-cell quarantine."""

from __future__ import annotations

from pathlib import Path
import json

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse, FancyBboxPatch, Rectangle
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, leaves_list, linkage
from scipy.spatial.distance import squareform
from scipy.stats import skew

import rerun_morph_pca_manual_pyramidal5_redundancy_first as prior
from rerun_morph_redundancy_skew_after_quarantine22 import flatten_quarantine


ROOT = Path(__file__).resolve().parents[1]
RAW_FILE = ROOT / "outputs/morph_qc/all22_shift_sum10000_log1p_zscore_215cells/01_all22_raw_215cells.csv"
QUARANTINE_FILE = ROOT / "config/morph_current_reconstruction_spe_local_pyramid_quarantine.json"
OUT = ROOT / "outputs/morph_qc/quarantine22_all22_including_soma_redundancy_skew_193cells"
ID = "MSN_unique_ID"
RHO_THRESHOLD = 0.80

ALL22 = [
    "M_soma_average_radius", "M_soma_maximal_radius", "M_soma_minimal_radius",
    "M_soma_max_pairwise_dist", "M_soma_perimeter", "M_soma_area",
    "M_soma_circularity_index", "M_soma_shape_factor", "M_soma_aspect_ratio",
    "M_cell_max_radial_dist", "M_total_number_of_neurites",
    "M_basal_dendrite_longthest_path", "M_basal_dendrite_Nseg",
    "M_basal_dendrite_avg_tortuosity", "M_basal_dendrite_max_tortuosity",
    "M_Total_neurite_length_(sections)", "M_Total_neurite_volume",
    "M_Total_neurite_area", "M_Number_of_bifurcation_points",
    "M_Maximum_branch_order", "M_trunk_angle_min", "M_trunk_angle_max",
]

SHORT = {
    "M_soma_average_radius": "Soma avg radius",
    "M_soma_maximal_radius": "Soma max radius",
    "M_soma_minimal_radius": "Soma min radius",
    "M_soma_max_pairwise_dist": "Soma max pairwise dist",
    "M_soma_perimeter": "Soma perimeter",
    "M_soma_area": "Soma area",
    "M_soma_circularity_index": "Soma circularity",
    "M_soma_shape_factor": "Soma shape factor",
    "M_soma_aspect_ratio": "Soma aspect ratio",
    "M_cell_max_radial_dist": "Max radial distance",
    "M_total_number_of_neurites": "Primary neurite number",
    "M_basal_dendrite_longthest_path": "Longest dendrite path",
    "M_basal_dendrite_Nseg": "Basal dendrite segments",
    "M_basal_dendrite_avg_tortuosity": "Mean tortuosity",
    "M_basal_dendrite_max_tortuosity": "Max tortuosity",
    "M_Total_neurite_length_(sections)": "Total neurite length",
    "M_Total_neurite_volume": "Total neurite volume",
    "M_Total_neurite_area": "Total neurite area",
    "M_Number_of_bifurcation_points": "Bifurcation points",
    "M_Maximum_branch_order": "Maximum branch order",
    "M_trunk_angle_min": "Minimum trunk angle",
    "M_trunk_angle_max": "Maximum trunk angle",
}

ABBR = {
    "M_soma_average_radius": "SomaAvgR", "M_soma_maximal_radius": "SomaMaxR",
    "M_soma_minimal_radius": "SomaMinR", "M_soma_max_pairwise_dist": "SomaMaxPD",
    "M_soma_perimeter": "SomaPerim", "M_soma_area": "SomaArea",
    "M_soma_circularity_index": "SomaCirc", "M_soma_shape_factor": "SomaShapeF",
    "M_soma_aspect_ratio": "SomaAR", "M_cell_max_radial_dist": "RadialMax",
    "M_total_number_of_neurites": "PrimaryN", "M_basal_dendrite_longthest_path": "PathMax",
    "M_basal_dendrite_Nseg": "DendSegN", "M_basal_dendrite_avg_tortuosity": "TortAvg",
    "M_basal_dendrite_max_tortuosity": "TortMax", "M_Total_neurite_length_(sections)": "NeuriteLen",
    "M_Total_neurite_volume": "NeuriteVol", "M_Total_neurite_area": "NeuriteArea",
    "M_Number_of_bifurcation_points": "BifurcN", "M_Maximum_branch_order": "BranchOrdMax",
    "M_trunk_angle_min": "TrunkAngMin", "M_trunk_angle_max": "TrunkAngMax",
}

# Prespecified order determines which member of a directly redundant pair is
# retained.  It preserves complementary soma minimum-radius and maximum-span
# measures because their direct correlation is below the 0.80 threshold.
PRIORITY = [
    "M_soma_minimal_radius", "M_soma_max_pairwise_dist", "M_soma_aspect_ratio",
    "M_soma_circularity_index", "M_soma_average_radius", "M_soma_maximal_radius",
    "M_soma_perimeter", "M_soma_area", "M_cell_max_radial_dist",
    "M_total_number_of_neurites", "M_basal_dendrite_avg_tortuosity",
    "M_Total_neurite_length_(sections)", "M_Total_neurite_volume",
    "M_Number_of_bifurcation_points", "M_Maximum_branch_order",
    "M_trunk_angle_min", "M_trunk_angle_max", "M_basal_dendrite_longthest_path",
    "M_basal_dendrite_max_tortuosity", "M_basal_dendrite_Nseg",
    "M_soma_shape_factor", "M_Total_neurite_area",
]


def greedy_direct_filter(corr: pd.DataFrame):
    retained: list[str] = []
    records = []
    for feature in PRIORITY:
        blockers = [r for r in retained if abs(float(corr.loc[feature, r])) >= RHO_THRESHOLD]
        if blockers:
            representative = max(blockers, key=lambda r: abs(float(corr.loc[feature, r])))
            records.append({
                "Feature": feature, "Retained": False, "Representative": representative,
                "Spearman_rho_with_representative": float(corr.loc[feature, representative]),
                "Decision": f"removed_direct_abs_rho_ge_{RHO_THRESHOLD:.2f}",
            })
        else:
            retained.append(feature)
            records.append({
                "Feature": feature, "Retained": True, "Representative": feature,
                "Spearman_rho_with_representative": 1.0, "Decision": "retained",
            })
    retained = [f for f in ALL22 if f in retained]
    return pd.DataFrame(records), retained


def plot_redundancy(corr: pd.DataFrame, labeled: bool):
    # Build exact threshold-connected components first.  Hierarchical-average
    # cutting can split a component even when one of its member pairs exceeds
    # the prespecified threshold, which would visually miss a required box.
    parent = {feature: feature for feature in corr.columns}

    def find(feature):
        while parent[feature] != feature:
            parent[feature] = parent[parent[feature]]
            feature = parent[feature]
        return feature

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra

    features = list(corr.columns)
    for i, a in enumerate(features):
        for b in features[i + 1:]:
            if abs(float(corr.loc[a, b])) >= RHO_THRESHOLD:
                union(a, b)

    distance = 1 - corr.abs()
    np.fill_diagonal(distance.values, 0)
    link = linkage(squareform(distance.values, checks=False), method="average")
    base_order = [features[i] for i in leaves_list(link)]
    base_position = {feature: i for i, feature in enumerate(base_order)}
    component_members: dict[str, list[str]] = {}
    for feature in features:
        component_members.setdefault(find(feature), []).append(feature)
    components = [sorted(members, key=base_position.get) for members in component_members.values()]
    components.sort(key=lambda members: min(base_position[f] for f in members))
    ordered_features = [feature for members in components for feature in members]
    ordered = corr.loc[ordered_features, ordered_features]
    n = len(ordered)
    fig_size = 7.2 if labeled else 6.2
    fig, ax = plt.subplots(figsize=(fig_size, fig_size))
    cmap = mpl.colors.LinearSegmentedColormap.from_list(
        "corr", [(0, "#0571B0"), (.35, "#74ADD1"), (.5, "#EEEEEE"), (.65, "#F28E6B"), (1, "#CA0020")]
    )
    norm = mpl.colors.Normalize(-1, 1)
    for row in range(n):
        for col in range(row, n):
            rho = float(ordered.iat[row, col])
            ax.add_patch(Rectangle((col, row), 1, 1, facecolor="white", edgecolor="#111111", linewidth=.38))
            ax.add_patch(Ellipse(
                (col + .5, row + .5), width=.86, height=max(.045, .86 * (1 - abs(rho))),
                angle=45 if rho >= 0 else -45, facecolor=cmap(norm(rho)), edgecolor="none",
            ))
            ax.text(col + .5, row + .5, f"{rho:.2f}", ha="center", va="center", fontsize=3.5)
    start = 0
    component_rows = []
    for component_id, members in enumerate(components, start=1):
        size = len(members)
        if size >= 2:
            ax.add_patch(FancyBboxPatch(
                (start - .08, start - .08), size + .16, size + .16,
                boxstyle="round,pad=0,rounding_size=0.18",
                fill=False, edgecolor="#111111", linewidth=1.15,
                clip_on=False, zorder=6,
            ))
        for member in members:
            component_rows.append({
                "Threshold_component": component_id,
                "Component_size": size,
                "Contains_abs_rho_ge_0p80_pair": size >= 2,
                "Feature": member,
            })
        start += size
    ax.set_xlim(-.05, n + .05); ax.set_ylim(n + .05, -.05); ax.set_aspect("equal")
    if labeled:
        labels = [ABBR[f] for f in ordered.columns]
        ax.set_xticks(np.arange(n) + .5, labels, rotation=90, fontsize=4)
        ax.xaxis.tick_top()
        ax.set_yticks(np.arange(n) + .5, labels, fontsize=4)
        ax.tick_params(length=0, pad=1)
        for spine in ax.spines.values(): spine.set_visible(False)
    else:
        ax.axis("off")
    cax = ax.inset_axes([.018, .31, .018, .28])
    cb = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax, ticks=[-1, -.5, 0, .5, 1])
    cb.ax.tick_params(labelsize=5, length=2, width=.5); cb.outline.set_linewidth(.6)
    suffix = "labeled" if labeled else "clean"
    fig.savefig(OUT / f"01_upper_triangle_all22_spearman_redundancy_{suffix}.png", dpi=600, bbox_inches="tight", facecolor="white")
    fig.savefig(OUT / f"01_upper_triangle_all22_spearman_redundancy_{suffix}.pdf", bbox_inches="tight", facecolor="white")
    plt.close(fig)
    component_table = pd.DataFrame(component_rows)
    component_table["Matrix_order"] = np.arange(1, n + 1)
    component_table["Abbreviation"] = component_table["Feature"].map(ABBR)
    component_table.to_csv(OUT / "01d_heatmap_feature_order_and_threshold_components.csv", index=False)
    return ordered_features


def plot_abs_skew(
    skew_df: pd.DataFrame,
    heatmap_order: list[str],
    output_stem: str,
    figure_height: float,
):
    # Preserve the left-to-right/top-to-bottom order of the redundancy map,
    # restricted to the features retained after direct redundancy filtering.
    retained_set = set(skew_df["Feature"])
    plot_order = [feature for feature in heatmap_order if feature in retained_set]
    data = skew_df.set_index("Feature").loc[plot_order].reset_index()
    fig, ax = plt.subplots(figsize=(3.35, figure_height))
    y = np.arange(len(data))
    colors = np.where(data["Signed_skewness"].lt(0), "#00468B", "#ED0000")
    ax.hlines(y, 0, data["Absolute_skewness"], color=colors, lw=1.15)
    ax.scatter(data["Absolute_skewness"], y, c=colors, s=18, zorder=3)
    ax.axvline(2, color="#ED0000", ls="--", lw=.85)
    ax.set_yticks(y, [SHORT[f] for f in data["Feature"]], fontsize=5.5)
    ax.invert_yaxis()
    # Explicit vertical padding prevents the first and last lollipops from
    # appearing clipped when the publication export uses a tight bounding box.
    ax.set_ylim(len(data) - .30, -.80)
    ax.set_xlabel("Absolute skewness after transformation", fontsize=7)
    upper = int(np.ceil(max(2.1, data["Absolute_skewness"].max())))
    ax.set_xlim(-.03, upper + .08); ax.set_xticks(np.arange(0, upper + 1, 1))
    ax.grid(axis="x", color="#E6E6E6", lw=.45); ax.set_axisbelow(True)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(
        OUT / f"{output_stem}.png", dpi=600, bbox_inches="tight",
        pad_inches=.14, facecolor="white",
    )
    fig.savefig(
        OUT / f"{output_stem}.pdf", bbox_inches="tight",
        pad_inches=.14, facecolor="white",
    )
    plt.close(fig)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    prior.setup_plotting()
    raw = pd.read_csv(RAW_FILE)
    config = json.loads(QUARANTINE_FILE.read_text(encoding="utf-8"))
    categories, quarantine = flatten_quarantine(config)
    cohort = raw.loc[~raw[ID].isin(quarantine)].copy().reset_index(drop=True)
    if len(cohort) != 193 or list(raw.columns[1:]) != ALL22:
        raise RuntimeError("Expected the frozen all-22 table and 193-cell post-quarantine cohort")

    logged, z_all, transform = prior.transcriptomic_like_transform(cohort, ALL22)
    corr = logged.corr(method="spearman")
    pairs = []
    for i, a in enumerate(ALL22):
        for b in ALL22[i + 1:]:
            rho = float(corr.loc[a, b])
            if abs(rho) >= RHO_THRESHOLD:
                pairs.append({"Feature_1": a, "Feature_2": b, "Spearman_rho": rho, "Absolute_rho": abs(rho)})
    decisions, retained = greedy_direct_filter(corr)
    corr.to_csv(OUT / "01a_all22_spearman_correlation_matrix.csv")
    pd.DataFrame(pairs).sort_values("Absolute_rho", ascending=False).to_csv(
        OUT / "01b_all22_redundant_pairs_abs_rho_ge0p80.csv", index=False
    )
    decisions.to_csv(OUT / "01c_all22_redundancy_decisions.csv", index=False)
    heatmap_order = plot_redundancy(corr, labeled=False)
    plot_redundancy(corr, labeled=True)

    skew_rows = []
    for feature in retained:
        signed = float(skew(logged[feature].to_numpy(float), bias=False))
        skew_rows.append({
            "Feature": feature, "Short_name": SHORT[feature], "Signed_skewness": signed,
            "Absolute_skewness": abs(signed),
            "Status": "highly_skewed" if abs(signed) >= 2 else (
                "moderately_skewed" if abs(signed) >= 1 else "acceptable"
            ),
            "Known_unit_discontinuity": feature in {"M_Total_neurite_volume", "M_Total_neurite_area"},
        })
    skew_df = pd.DataFrame(skew_rows).sort_values("Absolute_skewness", ascending=False)
    skew_df.to_csv(OUT / "02a_absolute_skewness_after_redundancy_removal.csv", index=False)
    plot_abs_skew(
        skew_df, heatmap_order,
        "02_absolute_skewness_all22_after_redundancy_removal", 4.15,
    )

    all22_skew_rows = []
    for feature in ALL22:
        signed = float(skew(logged[feature].to_numpy(float), bias=False))
        all22_skew_rows.append({
            "Feature": feature, "Short_name": SHORT[feature],
            "Signed_skewness": signed, "Absolute_skewness": abs(signed),
            "Status": "highly_skewed" if abs(signed) >= 2 else (
                "moderately_skewed" if abs(signed) >= 1 else "acceptable"
            ),
            "Retained_after_direct_redundancy_filter": feature in retained,
            "Known_unit_discontinuity": feature in {"M_Total_neurite_volume", "M_Total_neurite_area"},
        })
    all22_skew_df = pd.DataFrame(all22_skew_rows)
    all22_skew_df["Heatmap_order"] = all22_skew_df["Feature"].map(
        {feature: i + 1 for i, feature in enumerate(heatmap_order)}
    )
    all22_skew_df = all22_skew_df.sort_values("Heatmap_order")
    all22_skew_df.to_csv(OUT / "02b_absolute_skewness_all22_in_heatmap_order.csv", index=False)
    plot_abs_skew(
        all22_skew_df, heatmap_order,
        "02c_absolute_skewness_all22_in_heatmap_order", 6.0,
    )

    processed = z_all[retained].copy(); processed.insert(0, ID, cohort[ID].to_numpy())
    processed.to_csv(OUT / "03_all22_redundancy_filtered_z_matrix_193cells.csv", index=False)
    transform["Retained_after_direct_redundancy_filter"] = transform["Feature"].isin(retained)
    transform["Known_unit_discontinuity"] = transform["Feature"].isin({"M_Total_neurite_volume", "M_Total_neurite_area"})
    transform.to_csv(OUT / "03a_all22_preprocessing_parameters.csv", index=False)

    recommended = [f for f in retained if f not in {"M_Total_neurite_volume", "M_Total_neurite_area"}]
    summary = {
        "input_features_n": len(ALL22), "soma_features_included_n": 9,
        "frozen_input_cohort_n": 215, "temporary_quarantine_n": len(quarantine),
        "analysis_cohort_n": len(cohort), "redundancy_threshold_abs_spearman": RHO_THRESHOLD,
        "redundant_pairs_n": len(pairs), "retained_after_direct_redundancy_filter_n": len(retained),
        "retained_after_direct_redundancy_filter": retained,
        "recommended_for_downstream_after_known_unit_QC": recommended,
        "features_abs_skew_ge_2": skew_df.loc[skew_df["Absolute_skewness"].ge(2), "Feature"].tolist(),
        "batch_used": False, "permanent_global_exclusions_reapplied": False,
    }
    (OUT / "04_analysis_summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps(summary, indent=2, ensure_ascii=False))
    print("\nSkewness")
    print(skew_df.to_string(index=False))
    print("\nDecisions")
    print(decisions.to_string(index=False))


if __name__ == "__main__":
    main()
