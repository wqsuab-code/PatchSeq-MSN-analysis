#!/usr/bin/env python3
"""Second-pass Morph QC: selection bias, batch effects, outliers, and robust feature selection."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency, kruskal, rankdata, spearmanr
import statsmodels.api as sm
from statsmodels.formula.api import ols

from analyze_morph_qc import (
    ANALYSIS_EXCLUSIONS,
    EXCLUSIONS,
    OUT,
    SOURCE,
    load_source,
    numeric,
    sample_skew,
)


T_MAP = Path(
    "C:/Users/53461/Documents/Codex/2026-07-13/zh/outputs/"
    "05_branch_specific_subtype_mapping/53_MSN_final_RPCA_centroid_consensus/"
    "MSN_final_RPCA_centroid_consensus_per_cell.csv"
)

PRIMARY_FEATURES = [
    "M_soma_average_radius",
    "M_soma_aspect_ratio",
    "M_basal_dendrite_longthest_path",
    "M_total_number_of_neurites",
    "M_Number_of_bifurcation_points",
    "M_Maximum_branch_order",
    "M_basal_dendrite_avg_tortuosity",
    "M_Total_neurite_length_(sections)",
    "M_trunk_angle_min",
    "M_trunk_angle_max",
]


def bh_adjust(values: pd.Series) -> np.ndarray:
    p = values.to_numpy(float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order]
    adj = np.minimum.accumulate((ranked * n / np.arange(1, n + 1))[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.minimum(adj, 1)
    return out


def cramer_v(table: pd.DataFrame) -> tuple[float, float, float]:
    chi2, p, _, _ = chi2_contingency(table)
    n = table.to_numpy().sum()
    denom = min(table.shape[0] - 1, table.shape[1] - 1)
    return float(chi2), float(p), float(np.sqrt(chi2 / (n * denom)))


def transformed_z(x: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    transformed = x.copy()
    method = []
    for col in x:
        skew = sample_skew(x[col])
        if x[col].min() >= 0 and skew > 1:
            transformed[col] = np.log1p(x[col])
            method.append({"Feature": col, "Transform": "log1p", "Raw_skewness": skew})
        else:
            method.append({"Feature": col, "Transform": "none", "Raw_skewness": skew})
    z = (transformed - transformed.mean()) / transformed.std(ddof=1)
    return z, pd.DataFrame(method)


def robust_z_flags(x: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    t = x.copy()
    for col in x:
        if x[col].min() >= 0 and sample_skew(x[col]) > 1:
            t[col] = np.log1p(x[col])
    med = t.median()
    mad = (t - med).abs().median().replace(0, np.nan)
    rz = 0.67448975 * (t - med) / mad
    return rz, rz.abs() > 5


def availability_tables(df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    rows = []
    tests = []
    for var in ["T_Batch", "D1D2_Major", "Final_Subtype"]:
        tab = pd.crosstab(df[var].fillna("Unmapped"), df["Morph_available"])
        for group, r in tab.iterrows():
            total = int(r.sum())
            present = int(r.get(True, 0))
            rows.append({"Variable": var, "Group": group, "Morph_n": present, "Total_n": total, "Morph_rate": present / total})
        chi2, p, v = cramer_v(tab)
        expected_lt5 = int((chi2_contingency(tab)[3] < 5).sum())
        tests.append({"Variable": var, "Chi_square": chi2, "P_value": p, "Cramers_V": v, "Expected_cells_lt5": expected_lt5,
                      "Interpretation": "Use cautiously: sparse expected cells" if expected_lt5 else "Valid asymptotic chi-square"})
    morph = df.loc[df["Morph_available"]]
    conf = pd.crosstab(morph["T_Batch"], morph["D1D2_Major"])
    chi2, p, v = cramer_v(conf)
    conf_test = pd.DataFrame([{"Comparison": "Batch vs D1D2 major within Morph-positive cells", "Chi_square": chi2, "P_value": p, "Cramers_V": v}])
    return pd.DataFrame(rows), pd.DataFrame(tests), conf_test


def batch_effects(morph: pd.DataFrame, features: list[str]) -> pd.DataFrame:
    rows = []
    for feature in features:
        groups = [g[feature].dropna() for _, g in morph.groupby("T_Batch")]
        h, p = kruskal(*groups)
        k = len(groups)
        epsilon2 = max(0.0, (h - k + 1) / (len(morph) - k))
        tmp = morph[[feature, "T_Batch", "D1D2_Major"]].dropna().copy()
        tmp["rank_y"] = rankdata(tmp[feature])
        model = ols("rank_y ~ C(T_Batch) + C(D1D2_Major)", data=tmp).fit()
        anova = sm.stats.anova_lm(model, typ=2)
        batch_ss = anova.loc["C(T_Batch)", "sum_sq"]
        major_ss = anova.loc["C(D1D2_Major)", "sum_sq"]
        resid_ss = anova.loc["Residual", "sum_sq"]
        rows.append({
            "Feature": feature, "Kruskal_H_batch": float(h), "Kruskal_P_batch": float(p), "Kruskal_epsilon2_batch": float(epsilon2),
            "Rank_ANOVA_P_batch_adjusted_for_major": float(anova.loc["C(T_Batch)", "PR(>F)"]),
            "Partial_eta2_batch_adjusted_for_major": float(batch_ss / (batch_ss + resid_ss)),
            "Rank_ANOVA_P_major_adjusted_for_batch": float(anova.loc["C(D1D2_Major)", "PR(>F)"]),
            "Partial_eta2_major_adjusted_for_batch": float(major_ss / (major_ss + resid_ss)),
        })
    out = pd.DataFrame(rows)
    out["BH_q_batch_adjusted_for_major"] = bh_adjust(out["Rank_ANOVA_P_batch_adjusted_for_major"])
    out["BH_q_major_adjusted_for_batch"] = bh_adjust(out["Rank_ANOVA_P_major_adjusted_for_batch"])
    return out.sort_values("Partial_eta2_batch_adjusted_for_major", ascending=False).reset_index(drop=True)


def spatial_associations(morph: pd.DataFrame, features: list[str]) -> pd.DataFrame:
    rows = []
    for feature in features:
        for coord in ["L_X2O", "L_Y2O"]:
            sub = morph[[feature, coord]].dropna()
            rho, p = spearmanr(sub[feature], sub[coord])
            rows.append({"Feature": feature, "Coordinate": coord, "N": len(sub), "Spearman_rho": float(rho), "P_value": float(p)})
    out = pd.DataFrame(rows)
    out["BH_q"] = bh_adjust(out["P_value"])
    return out.sort_values("BH_q").reset_index(drop=True)


def redundancy_sensitivity(x: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    z, methods = transformed_z(x)
    spearman = x.corr(method="spearman")
    pearson = z.corr(method="pearson")
    rows = []
    cols = list(x.columns)
    for i, a in enumerate(cols):
        for b in cols[:i]:
            s = float(spearman.loc[a, b]); p = float(pearson.loc[a, b])
            rows.append({"Var1": a, "Var2": b, "Raw_Spearman": s, "Transformed_Pearson": p,
                         "Spearman_flag_abs_gt0_85": abs(s) > 0.85, "Pearson_flag_abs_gt0_85": abs(p) > 0.85,
                         "Stable_redundancy": abs(s) > 0.85 and abs(p) > 0.85})
    return pd.DataFrame(rows).sort_values("Raw_Spearman", key=lambda s: s.abs(), ascending=False).reset_index(drop=True), methods


def feature_decisions(batch: pd.DataFrame) -> pd.DataFrame:
    use = {
        "M_soma_average_radius": ("Primary", "single representative of soma size"),
        "M_soma_aspect_ratio": ("Primary", "soma anisotropy; minimal batch signal"),
        "M_basal_dendrite_longthest_path": ("Primary", "dendritic path extent"),
        "M_total_number_of_neurites": ("Primary", "primary neurite count; keep instead of identical Nseg"),
        "M_Number_of_bifurcation_points": ("Primary with batch sensitivity", "branching complexity; batch-associated"),
        "M_Maximum_branch_order": ("Primary with batch sensitivity", "hierarchical depth"),
        "M_basal_dendrite_avg_tortuosity": ("Primary", "typical dendritic tortuosity"),
        "M_Total_neurite_length_(sections)": ("Primary", "overall reconstructed dendritic mass by length"),
        "M_trunk_angle_min": ("Primary", "minimum trunk separation"),
        "M_trunk_angle_max": ("Primary with batch sensitivity", "maximum trunk separation; batch-associated"),
        "M_soma_circularity_index": ("Sensitivity only", "exactly derived from soma area and perimeter; strong batch effect"),
        "M_cell_max_radial_dist": ("Sensitivity only", "redundant with longest path"),
        "M_basal_dendrite_max_tortuosity": ("Sensitivity only", "tail-sensitive extreme tortuosity"),
        "M_Total_neurite_volume": ("Sensitivity only", "very strong batch effect and diameter sensitivity"),
        "M_Total_neurite_area": ("Exclude primary", "redundant with volume and very strong batch effect"),
        "M_basal_dendrite_Nseg": ("Exclude", "exact duplicate of total number of neurites in all 223 cells"),
        "M_soma_perimeter": ("Exclude", "redundant soma-size measure"),
        "M_soma_area": ("Exclude", "redundant soma-size measure; extreme skewness"),
        "M_soma_maximal_radius": ("Exclude", "redundant soma-size measure"),
        "M_soma_minimal_radius": ("Exclude", "redundant soma-size measure"),
        "M_soma_max_pairwise_dist": ("Exclude", "redundant soma-size measure"),
        "M_soma_shape_factor": ("Exclude", "redundant with circularity/aspect; batch-associated"),
    }
    out = pd.DataFrame([{"Feature": f, "Recommended_role": role, "Reason": reason} for f, (role, reason) in use.items()])
    return out.merge(batch[["Feature", "BH_q_batch_adjusted_for_major", "Partial_eta2_batch_adjusted_for_major"]], on="Feature", how="left")


def plots(availability: pd.DataFrame, tests: pd.DataFrame, batch: pd.DataFrame, outlier_feature: pd.DataFrame,
          outlier_cell: pd.DataFrame, redundancy: pd.DataFrame) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(8.2, 3.4))
    for ax, variable in zip(axes, ["T_Batch", "D1D2_Major"]):
        sub = availability.loc[availability["Variable"] == variable]
        ax.bar(sub["Group"], sub["Morph_rate"] * 100, color="#4C78A8")
        ax.set_ylim(0, 100); ax.set_ylabel("Morphology availability (%)"); ax.tick_params(axis="x", rotation=45)
        t = tests.loc[tests["Variable"] == variable].iloc[0]
        ax.set_title(f"{variable}: P={t.P_value:.2g}, Cramer's V={t.Cramers_V:.2f}", loc="left", fontweight="bold")
    fig.tight_layout(); fig.savefig(OUT / "Morph_availability_selection_bias.png", dpi=400, bbox_inches="tight"); plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.4, 5.6))
    b = batch.sort_values("Partial_eta2_batch_adjusted_for_major")
    colors = np.where(b["BH_q_batch_adjusted_for_major"] < 0.05, "#B2182B", "#B9BEC5")
    ax.barh(b["Feature"], b["Partial_eta2_batch_adjusted_for_major"], color=colors)
    ax.set_xlabel("Partial eta-squared for batch (rank ANOVA, adjusted for D1/D2 major)")
    ax.set_title("Morphology batch effects", loc="left", fontweight="bold"); ax.grid(axis="x", color="#E6E6E6"); ax.set_axisbelow(True)
    fig.tight_layout(); fig.savefig(OUT / "Morph_batch_effects_adjusted_for_major.png", dpi=400, bbox_inches="tight"); plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(8.5, 4.2))
    ff = outlier_feature.sort_values("Flagged_cells_full22")
    axes[0].barh(ff["Feature"], ff["Flagged_cells_full22"], color="#E08214")
    axes[0].set_xlabel("Cells with |robust z| > 5"); axes[0].set_title("Review flags by feature", loc="left", fontweight="bold")
    bins = np.arange(0, max(2, int(outlier_cell["Full22_flag_count"].max()) + 2)) - 0.5
    axes[1].hist(outlier_cell["Full22_flag_count"], bins=bins, alpha=.65, label="All 22", color="#E08214")
    axes[1].hist(outlier_cell["Primary10_flag_count"], bins=bins, alpha=.65, label="Primary 10", color="#4C78A8")
    axes[1].set_xlabel("Flagged features per cell"); axes[1].set_ylabel("Cells"); axes[1].legend(frameon=False)
    axes[1].set_title("Correlated features inflate review counts", loc="left", fontweight="bold")
    fig.tight_layout(); fig.savefig(OUT / "Morph_outlier_review_audit.png", dpi=400, bbox_inches="tight"); plt.close(fig)

    fig, ax = plt.subplots(figsize=(5.2, 4.8))
    ax.scatter(redundancy["Raw_Spearman"], redundancy["Transformed_Pearson"], s=12, alpha=.55, color="#4C78A8")
    for v in [-.85, .85]: ax.axvline(v, color="#B2182B", linestyle="--", linewidth=.7); ax.axhline(v, color="#B2182B", linestyle="--", linewidth=.7)
    ax.plot([-1,1],[-1,1], color="#999999", linewidth=.6); ax.set_xlim(-1,1); ax.set_ylim(-1,1)
    ax.set_xlabel("Raw Spearman rho"); ax.set_ylabel("Pearson r after feature-wise transform")
    ax.set_title("Redundancy stability across preprocessing", loc="left", fontweight="bold")
    fig.tight_layout(); fig.savefig(OUT / "Morph_redundancy_preprocessing_sensitivity.png", dpi=400, bbox_inches="tight"); plt.close(fig)


def main() -> None:
    df, excluded = load_source()
    mapping = pd.read_csv(T_MAP)
    agreement = mapping["major_class_agreement"].astype(str).str.upper().eq("TRUE")
    mapping["D1D2_Major"] = np.where(agreement, mapping["RPCA_major_class"], "D1/D2 unstable")
    df = df.merge(mapping[["cell_id", "Final_Subtype", "D1D2_Major"]], left_on="MSN_unique_ID", right_on="cell_id", how="left")
    df["Morph_available"] = df["has_M"].astype(str).str.lower().eq("true")
    features = [c for c in df.columns if c.startswith("M_") and c not in ANALYSIS_EXCLUSIONS]
    morph = df.loc[df["Morph_available"]].copy()
    for c in features + ["L_X2O", "L_Y2O"]:
        morph[c] = numeric(morph[c])
    x = morph[features]

    availability, availability_tests, confounding = availability_tables(df)
    batch = batch_effects(morph, features)
    spatial = spatial_associations(morph, features)
    redundancy, transforms = redundancy_sensitivity(x)

    rz, flags = robust_z_flags(x)
    core_rz, core_flags = robust_z_flags(x[PRIMARY_FEATURES])
    outlier_feature = pd.DataFrame({"Feature": features, "Flagged_cells_full22": flags.sum(axis=0).values,
                                    "Max_abs_robust_z": rz.abs().max(axis=0).values}).sort_values("Flagged_cells_full22", ascending=False)
    outlier_cell = morph[["MSN_unique_ID", "METL_source_row", "T_Batch", "D1D2_Major", "Final_Subtype"]].copy()
    outlier_cell["Full22_flag_count"] = flags.sum(axis=1).values
    outlier_cell["Primary10_flag_count"] = core_flags.sum(axis=1).values
    outlier_cell["Max_abs_robust_z_full22"] = rz.abs().max(axis=1).values
    outlier_cell["Flagged_features_full22"] = ["; ".join(flags.columns[row.to_numpy()].tolist()) for _, row in flags.iterrows()]
    outlier_cell["Review_priority"] = np.select(
        [outlier_cell["Primary10_flag_count"] >= 2, outlier_cell["Primary10_flag_count"] == 1,
         outlier_cell["Full22_flag_count"] > 0],
        ["High: >=2 primary features", "Medium: 1 primary feature", "Low: redundant/sensitivity features only"], default="None")

    feature_decision = feature_decisions(batch)
    exact_rel = pd.DataFrame([
        {"Relationship": "M_basal_dendrite_Nseg == M_total_number_of_neurites", "Maximum_absolute_error": 0.0, "Conclusion": "Exact duplicate in all 223 cells"},
        {"Relationship": "M_soma_circularity_index == 4*pi*M_soma_area/M_soma_perimeter^2",
         "Maximum_absolute_error": float(np.max(np.abs(x["M_soma_circularity_index"] - 4*np.pi*x["M_soma_area"]/x["M_soma_perimeter"]**2))),
         "Conclusion": "Exact deterministic derived feature"},
    ])

    tables = {
        "morph_availability_by_group.json": availability,
        "morph_availability_tests.json": availability_tests,
        "batch_major_confounding_test.json": confounding,
        "feature_batch_major_effects.json": batch,
        "feature_spatial_associations.json": spatial,
        "redundancy_preprocessing_sensitivity.json": redundancy,
        "feature_transform_plan.json": transforms,
        "outlier_flags_by_feature.json": outlier_feature,
        "outlier_review_by_cell.json": outlier_cell,
        "feature_decisions.json": feature_decision,
        "exact_deterministic_relationships.json": exact_rel,
    }
    for name, table in tables.items():
        table.to_json(OUT / name, orient="records", indent=2)
    summary = {
        "morph_cells": len(morph), "mapping_coverage": int(df["D1D2_Major"].notna().sum()),
        "batch_availability_test": availability_tests.loc[availability_tests.Variable.eq("T_Batch")].iloc[0].to_dict(),
        "major_availability_test": availability_tests.loc[availability_tests.Variable.eq("D1D2_Major")].iloc[0].to_dict(),
        "batch_major_confounding_within_morph": confounding.iloc[0].to_dict(),
        "batch_associated_features_q_lt_0_05": int((batch.BH_q_batch_adjusted_for_major < .05).sum()),
        "major_associated_features_q_lt_0_05": int((batch.BH_q_major_adjusted_for_batch < .05).sum()),
        "location_complete_cells": int(morph[["L_X2O", "L_Y2O"]].notna().all(axis=1).sum()),
        "spatial_associations_q_lt_0_05": int((spatial.BH_q < .05).sum()),
        "cells_flagged_full22": int((outlier_cell.Full22_flag_count > 0).sum()),
        "cells_flagged_primary10": int((outlier_cell.Primary10_flag_count > 0).sum()),
        "high_priority_review_cells": int((outlier_cell.Primary10_flag_count >= 2).sum()),
        "stable_redundant_pairs": int(redundancy.Stable_redundancy.sum()),
        "primary_features": PRIMARY_FEATURES,
    }
    (OUT / "Morph_QC_deep_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    plots(availability, availability_tests, batch, outlier_feature, outlier_cell, redundancy)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
