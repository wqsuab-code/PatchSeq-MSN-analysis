#!/usr/bin/env python3
"""Export an auditable statistical analysis package for the 25 E-features."""

from __future__ import annotations

import json
import math
from pathlib import Path
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from plot_all25_violin_plus_feature_tsne import E_LEVELS, load_data
from plot_all25_violin_plus_feature_tsne_posthoc import (
    bh_adjust,
    calculate_statistics,
    holm_adjust,
    significance_label,
)


ROOT = Path(__file__).resolve().parents[1]
OUT = (
    ROOT / "outputs" / "e_type_qc" / "fixed_npc3_res1.5_merged_k5_tsne"
    / "p45_i1500_seed777_four_label_views" / "all25_violin_plus_feature_tsne_triple438"
    / "statistical_analysis_audit"
)


def _betacf(a: float, b: float, x: float) -> float:
    """Continued fraction for the incomplete beta function."""
    max_iter, eps, fpmin = 300, 3e-14, 1e-300
    qab, qap, qam = a + b, a + 1.0, a - 1.0
    c = 1.0
    d = 1.0 - qab * x / qap
    if abs(d) < fpmin:
        d = fpmin
    d = 1.0 / d
    h = d
    for m in range(1, max_iter + 1):
        m2 = 2 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1.0 + aa * d
        if abs(d) < fpmin:
            d = fpmin
        c = 1.0 + aa / c
        if abs(c) < fpmin:
            c = fpmin
        d = 1.0 / d
        h *= d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1.0 + aa * d
        if abs(d) < fpmin:
            d = fpmin
        c = 1.0 + aa / c
        if abs(c) < fpmin:
            c = fpmin
        d = 1.0 / d
        delta = d * c
        h *= delta
        if abs(delta - 1.0) < eps:
            break
    return h


def regularized_beta(x: float, a: float, b: float) -> float:
    if x <= 0:
        return 0.0
    if x >= 1:
        return 1.0
    log_bt = (
        math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b)
        + a * math.log(x) + b * math.log1p(-x)
    )
    bt = math.exp(log_bt)
    if x < (a + 1.0) / (a + b + 2.0):
        return bt * _betacf(a, b, x) / a
    return 1.0 - bt * _betacf(b, a, 1.0 - x) / b


def f_survival(f_value: float, df1: int, df2: int) -> float:
    if not np.isfinite(f_value) or f_value < 0:
        return float("nan")
    # Complementary beta form avoids catastrophic cancellation for very small P.
    x = df2 / (df2 + df1 * f_value)
    return float(np.clip(regularized_beta(x, df2 / 2.0, df1 / 2.0), 0, 1))


def brown_forsythe(groups: list[np.ndarray]) -> tuple[float, int, int, float]:
    """Levene test centered on group medians (Brown-Forsythe variant)."""
    deviations = [np.abs(x - np.median(x)) for x in groups]
    k = len(deviations)
    n_total = sum(len(x) for x in deviations)
    grand = np.concatenate(deviations).mean()
    between = sum(len(x) * (x.mean() - grand) ** 2 for x in deviations)
    within = sum(np.sum((x - x.mean()) ** 2) for x in deviations)
    df1, df2 = k - 1, n_total - k
    statistic = (between / df1) / (within / df2) if within > 0 else float("inf")
    return float(statistic), df1, df2, f_survival(statistic, df1, df2)


def descriptive_rows(data: pd.DataFrame, fmap: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for order, feature_row in fmap.iterrows():
        feature = feature_row["Feature"]
        for group in E_LEVELS:
            x = data.loc[data["E_type"] == group, feature].to_numpy(dtype=float)
            n = len(x)
            mean = float(np.mean(x))
            sd = float(np.std(x, ddof=1))
            variance = sd**2
            sem = sd / math.sqrt(n)
            q1, median, q3 = np.quantile(x, [0.25, 0.50, 0.75])
            centered = x - mean
            m2 = float(np.mean(centered**2))
            skew = float(np.mean(centered**3) / m2**1.5) if m2 > 0 else 0.0
            kurt = float(np.mean(centered**4) / m2**2 - 3.0) if m2 > 0 else 0.0
            rows.append({
                "Order": order + 1,
                "Feature": feature,
                "Short_title": feature_row["Short_title"],
                "Unit": feature_row["Y_label"],
                "E_type": group,
                "N": n,
                "Missing_N": int(data.loc[data["E_type"] == group, feature].isna().sum()),
                "Mean": mean,
                "SD": sd,
                "Variance": variance,
                "SEM": sem,
                "Mean_CI95_low": mean - 1.96 * sem,
                "Mean_CI95_high": mean + 1.96 * sem,
                "Median": float(median),
                "Q1": float(q1),
                "Q3": float(q3),
                "IQR": float(q3 - q1),
                "Minimum": float(np.min(x)),
                "Maximum": float(np.max(x)),
                "Skewness_moment": skew,
                "Excess_kurtosis_moment": kurt,
                "CV_abs_percent": abs(sd / mean) * 100 if mean != 0 else float("nan"),
            })
    return pd.DataFrame(rows)


def variance_rows(data: pd.DataFrame, fmap: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for order, feature_row in fmap.iterrows():
        feature = feature_row["Feature"]
        groups = [data.loc[data["E_type"] == g, feature].to_numpy(float) for g in E_LEVELS]
        statistic, df1, df2, p = brown_forsythe(groups)
        variances = [np.var(x, ddof=1) for x in groups]
        rows.append({
            "Order": order + 1,
            "Feature": feature,
            "Short_title": feature_row["Short_title"],
            "Test": "Brown-Forsythe (median-centered Levene)",
            "Statistic_F": statistic,
            "df1": df1,
            "df2": df2,
            "P_raw": p,
            "Equal_variance_at_0.05": "Yes" if p >= 0.05 else "No",
            "Minimum_group_variance": float(np.min(variances)),
            "Maximum_group_variance": float(np.max(variances)),
            "Max_to_min_variance_ratio": float(np.max(variances) / np.min(variances)) if np.min(variances) > 0 else float("inf"),
        })
    result = pd.DataFrame(rows)
    result["P_adj_BH_across_25_features"] = bh_adjust(result["P_raw"].to_numpy())
    return result


def add_omnibus_effects(omnibus: pd.DataFrame) -> pd.DataFrame:
    out = omnibus.copy()
    out["Kruskal_epsilon_squared"] = np.maximum(
        0.0, (out["Kruskal_Wallis_H"] - 5 + 1) / (out["N"] - 5)
    )
    out["Significance_BH"] = out["P_adj_BH_across_25_features"].map(significance_label)
    out["Primary_inference"] = "Kruskal-Wallis, two-sided; BH adjusted across 25 features"
    return out


def add_pairwise_details(pairwise: pd.DataFrame, data: pd.DataFrame) -> pd.DataFrame:
    out = pairwise.copy()
    details = []
    for row in out.itertuples(index=False):
        x1 = data.loc[data["E_type"] == row.Group_1, row.Feature].to_numpy(float)
        x2 = data.loc[data["E_type"] == row.Group_2, row.Feature].to_numpy(float)
        details.append({
            "Mean_1": float(np.mean(x1)), "Mean_2": float(np.mean(x2)),
            "Mean_difference_1_minus_2": float(np.mean(x1) - np.mean(x2)),
            "Median_1": float(np.median(x1)), "Median_2": float(np.median(x2)),
            "Median_difference_1_minus_2": float(np.median(x1) - np.median(x2)),
            "Direction_by_median": row.Group_1 if np.median(x1) > np.median(x2) else row.Group_2,
        })
    detail_df = pd.DataFrame(details)
    insert_at = list(out.columns).index("Dunn_Z")
    out = pd.concat([out.iloc[:, :insert_at], detail_df, out.iloc[:, insert_at:]], axis=1)
    out["Reject_Holm_0.05"] = out["P_adj_Holm_within_feature"] < 0.05
    out["Reject_Holm_0.01"] = out["P_adj_Holm_within_feature"] < 0.01
    return out


def method_table() -> pd.DataFrame:
    return pd.DataFrame([
        ["Analysis cohort", "438 GC×HC×EC triple-consensus MSN cells; E1=141, E2=157, E3=44, E4=57, E5=39."],
        ["Features", "25 final electrophysiological metrics analyzed on their raw measurement scales."],
        ["Missing data", "No missing values in the active analysis matrix; all tests therefore use N=438."],
        ["Descriptive statistics", "N, mean, sample SD, sample variance, SEM, normal-approximation 95% CI of mean, median, Q1, Q3, IQR, range, moment skewness, excess kurtosis and absolute CV."],
        ["Variance analysis", "Brown-Forsythe test (median-centered Levene test), two-sided F reference distribution; BH correction across 25 features."],
        ["Primary omnibus comparison", "Two-sided Kruskal-Wallis rank test across E1-E5 with tie correction; df=4; BH correction across 25 features."],
        ["Omnibus effect size", "Epsilon-squared = max[0, (H-k+1)/(N-k)], where k=5."],
        ["Post hoc comparison", "Two-sided Dunn rank-sum comparisons for all 10 E-type pairs per feature, with tie correction."],
        ["Multiplicity: post hoc", "Holm family-wise correction across the 10 pairwise tests within each feature (primary pairwise inference); global BH values across all 250 tests are also reported."],
        ["Pairwise effect size", "Dunn r = |Z|/sqrt(N), using the total feature-level sample size."],
        ["Significance notation", "ns: P>=0.05; *: P<0.05; **: P<0.01; ***: P<0.001; ****: P<0.0001, based on Holm-adjusted within-feature P for pairwise results."],
        ["Software implementation", "Ranks, tie corrections, Kruskal-Wallis, Dunn, Holm and BH procedures were implemented in the project analysis scripts; Brown-Forsythe P values use the F distribution."],
        ["Interpretation", "The Kruskal-Wallis/Dunn workflow is the prespecified primary inference because several electrophysiological measures are skewed and group variances may differ."],
    ], columns=["Item", "Details"])


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    data, fmap = load_data()
    omnibus, pairwise = calculate_statistics(data, fmap)
    descriptive = descriptive_rows(data, fmap)
    variance = variance_rows(data, fmap)
    omnibus = add_omnibus_effects(omnibus)
    pairwise = add_pairwise_details(pairwise, data)
    methods = method_table()

    tables = {
        "methods": methods,
        "group_descriptives": descriptive,
        "variance_tests": variance,
        "omnibus_tests": omnibus,
        "pairwise_dunn_holm": pairwise,
    }
    for name, table in tables.items():
        table.to_csv(OUT / f"Ephys_Statistics_{name}.csv", index=False)
        (OUT / f"Ephys_Statistics_{name}.json").write_text(
            table.to_json(orient="records", double_precision=15), encoding="utf-8"
        )
    summary = {
        "N": len(data), "features": len(fmap), "groups": data["E_type"].value_counts().sort_index().to_dict(),
        "omnibus_BH_lt_0.05": int((omnibus["P_adj_BH_across_25_features"] < 0.05).sum()),
        "variance_BH_lt_0.05": int((variance["P_adj_BH_across_25_features"] < 0.05).sum()),
        "pairwise_Holm_lt_0.05": int((pairwise["P_adj_Holm_within_feature"] < 0.05).sum()),
        "pairwise_Holm_lt_0.01": int((pairwise["P_adj_Holm_within_feature"] < 0.01).sum()),
    }
    (OUT / "Ephys_Statistics_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))
    print(OUT)


if __name__ == "__main__":
    main()
