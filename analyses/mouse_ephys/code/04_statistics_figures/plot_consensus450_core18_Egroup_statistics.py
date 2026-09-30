#!/usr/bin/env python3
"""Raw-value E1-E5 comparisons for the frozen 18-feature E panel."""

from __future__ import annotations

import json
import math
import os
from pathlib import Path
import sys

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import MaxNLocator
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from plot_all25_violin_plus_feature_tsne_posthoc import (  # noqa: E402
    E_COLORS, E_LEVELS, assign_bracket_levels, bh_adjust, dunn_holm,
    kruskal_wallis,
)
from export_ephys_statistics_audit import brown_forsythe  # noqa: E402

FEATURE_FILE = ROOT / "outputs/e_type_qc/hc5_tree_distillation/NPC3_HC5_Raw_Final18_Features.csv"
LABEL_FILE = ROOT / "outputs/e_type_qc/HC_GC_only/HC_GC_only_cell_assignments.csv"
OUT = Path(os.environ.get(
    "MOUSE_E_STATS_OUT",
    ROOT / "outputs/e_type_qc/HC_GC_only/30_consensus450_core18_Egroup_statistics",
))

FEATURES = [
    ("E_Holding.MP..mV.", "Holding MP", "mV"),
    ("E_Input.resistance..MOhm.", "Input resistance", "MΩ"),
    ("E_Membrane.time.constant..ms.", "Membrane tau", "ms"),
    ("E_Rheobase..pA.", "Rheobase", "pA"),
    ("E_Sag.ratio", "Sag ratio", "ratio"),
    ("E_Sag.time..s.", "Sag time", "s"),
    ("E_AP.threshold..mV.", "AP threshold", "mV"),
    ("E_AP.amplitude..mV.", "AP amplitude", "mV"),
    ("E_AP.width..ms.", "AP width", "ms"),
    ("E_Upstroke.to.downstroke.ratio", "Up/down ratio", "ratio"),
    ("E_Afterhyperpolarization..mV.", "AHP", "mV"),
    ("E_Max.number.of.APs", "Max APs", "count"),
    ("E_Latency..ms.", "Latency", "ms"),
    ("E_Latency....20pA.current..ms.", "Latency @ +20 pA", "ms"),
    ("E_ISI.adaptation.index", "ISI adaptation", "index"),
    ("E_ISI.coefficient.of.variation", "ISI CV", "CV"),
    ("E_AP.amplitude.adaptation.index", "AP amplitude adaptation", "index"),
    ("E_AP.coefficient.of.variation", "AP CV", "CV"),
]


def load_data() -> pd.DataFrame:
    raw = pd.read_csv(FEATURE_FILE)
    labels = pd.read_csv(LABEL_FILE)
    keep = labels["HC_GC_consensus"].astype(str).str.lower().eq("true")
    data = labels.loc[keep, ["MSN_unique_ID", "HC_GC_consensus_E"]].merge(
        raw[["MSN_unique_ID", *[x[0] for x in FEATURES]]],
        on="MSN_unique_ID", how="inner", validate="one_to_one",
    ).rename(columns={"HC_GC_consensus_E": "E_type"})
    data["E_type"] = pd.Categorical(data["E_type"], E_LEVELS, ordered=True)
    if len(data) != 450 or data["E_type"].isna().any():
        raise ValueError(f"Expected 450 labeled consensus cells; found {len(data)}")
    if data[[x[0] for x in FEATURES]].isna().any().any():
        raise ValueError("The frozen core18 analysis matrix contains missing values")
    return data


def descriptive_statistics(data: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for order, (feature, title, unit) in enumerate(FEATURES, 1):
        for group in E_LEVELS:
            x = data.loc[data["E_type"] == group, feature].to_numpy(float)
            n, mean = len(x), float(x.mean())
            sd = float(x.std(ddof=1)); sem = sd / math.sqrt(n)
            q1, median, q3 = np.quantile(x, [.25, .5, .75])
            rows.append({
                "Order": order, "Feature": feature, "Display_name": title,
                "Unit": unit, "E_type": group, "N": n, "Mean": mean,
                "SD": sd, "Variance": sd**2, "SEM": sem,
                "Mean_CI95_low": mean - 1.96 * sem,
                "Mean_CI95_high": mean + 1.96 * sem,
                "Minimum": float(x.min()), "Q1": float(q1),
                "Median": float(median), "Q3": float(q3),
                "Maximum": float(x.max()), "IQR": float(q3-q1),
            })
    return pd.DataFrame(rows)


def inferential_statistics(data: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    groups = data["E_type"].astype(str).to_numpy()
    omnibus, pairs, variances = [], [], []
    for order, (feature, title, unit) in enumerate(FEATURES, 1):
        values = data[feature].to_numpy(float)
        h, p = kruskal_wallis(values, groups)
        omnibus.append({"Order": order, "Feature": feature, "Display_name": title,
                         "Unit": unit, "N": len(values), "Kruskal_Wallis_H": h,
                         "df": 4, "P_raw": p})
        ptab = dunn_holm(values, groups)
        ptab.insert(0, "Unit", unit); ptab.insert(0, "Display_name", title)
        ptab.insert(0, "Feature", feature); ptab.insert(0, "Order", order)
        medians = {g: float(np.median(data.loc[data.E_type == g, feature])) for g in E_LEVELS}
        ptab["Median_1"] = ptab.Group_1.map(medians)
        ptab["Median_2"] = ptab.Group_2.map(medians)
        ptab["Median_difference_1_minus_2"] = ptab.Median_1 - ptab.Median_2
        pairs.append(ptab)
        arrays = [data.loc[data.E_type == g, feature].to_numpy(float) for g in E_LEVELS]
        f, df1, df2, vp = brown_forsythe(arrays)
        variances.append({"Order": order, "Feature": feature, "Display_name": title,
                          "Unit": unit, "Test": "Brown-Forsythe (median-centered Levene)",
                          "Statistic_F": f, "df1": df1, "df2": df2, "P_raw": vp})
    omni = pd.DataFrame(omnibus)
    omni["P_adj_BH_across_18_features"] = bh_adjust(omni.P_raw.to_numpy())
    omni["Kruskal_epsilon_squared"] = np.maximum(0, (omni.Kruskal_Wallis_H - 4) / (omni.N - 5))
    pair = pd.concat(pairs, ignore_index=True)
    pair["P_adj_BH_global_180_tests"] = bh_adjust(pair.P_raw.to_numpy())
    pair["Reject_Holm_0.05"] = pair.P_adj_Holm_within_feature < .05
    variance = pd.DataFrame(variances)
    variance["P_adj_BH_across_18_features"] = bh_adjust(variance.P_raw.to_numpy())
    return omni, pair, variance


def p_text(value: float) -> str:
    return f"{value:.1e}".replace("e-0", "e−").replace("e-", "e−").replace("e+0", "e+")


def draw_panel(ax, data: pd.DataFrame, feature_info, pairwise: pd.DataFrame, seed: int) -> None:
    feature, title, unit = feature_info
    positions = np.arange(1, 6)
    grouped = [data.loc[data.E_type == g, feature].to_numpy(float) for g in E_LEVELS]
    rng = np.random.default_rng(seed)
    for pos, group, values in zip(positions, E_LEVELS, grouped):
        ax.scatter(pos + rng.uniform(-.23, .23, len(values)), values, s=4.0,
                   color=E_COLORS[group], alpha=.5, edgecolors="none", zorder=2)
    bp = ax.boxplot(grouped, positions=positions, widths=.48, notch=True,
                    patch_artist=True, showfliers=False, whis=1.5,
                    medianprops={"color": "#111111", "linewidth": .45},
                    boxprops={"facecolor": "white", "edgecolor": "#111111", "linewidth": .35},
                    whiskerprops={"color": "#111111", "linewidth": .35},
                    capprops={"color": "#111111", "linewidth": .35})
    for patch in bp["boxes"]:
        patch.set_alpha(.88)

    values = data[feature].to_numpy(float)
    dmin, dmax = float(values.min()), float(values.max())
    span = dmax - dmin if dmax > dmin else max(abs(dmin), 1.0)
    visible = assign_bracket_levels(pairwise.loc[pairwise.P_adj_Holm_within_feature < .05].copy())
    nlevels = int(visible.Level.max()) if not visible.empty else 0
    lower = dmin - .06 * span
    data_top = dmax + .08 * span
    upper = data_top + nlevels * .105 * span + (.06 * span if nlevels else 0)
    for row in visible.itertuples(index=False):
        y = data_top + (row.Level - .45) * .105 * span
        ax.plot([row.X_1, row.X_2], [y, y], color="#111111", lw=.32,
                solid_capstyle="butt", clip_on=False, zorder=4)
        ax.text((row.X_1 + row.X_2) / 2, y + .018 * span,
                p_text(row.P_adj_Holm_within_feature), ha="center", va="bottom",
                fontsize=3.3, color="#111111", clip_on=False, zorder=5)

    locator = MaxNLocator(nbins=3)
    ticks = locator.tick_values(dmin, dmax)
    ticks = ticks[(ticks >= lower) & (ticks <= data_top)]
    ax.set_yticks(ticks)
    ax.set_ylim(lower, upper)
    ax.set_xlim(.5, 5.5)
    ax.set_xticks(positions, E_LEVELS)
    ax.set_ylabel(f"{title} ({unit})" if unit not in {"ratio", "index", "CV"} else title,
                  fontsize=4.5, labelpad=3)
    ax.tick_params(axis="both", direction="out", length=1.4, width=.4, pad=1.6, labelsize=4)
    ax.tick_params(axis="y", labelrotation=90)
    ax.spines[["top", "right"]].set_visible(False)
    ax.spines[["left", "bottom"]].set_linewidth(.5)


def save_figures(data: pd.DataFrame, pairwise: pd.DataFrame) -> pd.DataFrame:
    mpl.rcParams.update({"font.family": "Arial", "pdf.fonttype": 42, "ps.fonttype": 42})
    png_dir, pdf_dir = OUT / "individual_png_900dpi", OUT / "individual_pdf"
    png_dir.mkdir(parents=True, exist_ok=True); pdf_dir.mkdir(parents=True, exist_ok=True)
    manifest = []
    multipage = OUT / "Core18_E1-E5_raw_distributions_Dunn-Holm_18pages.pdf"
    with PdfPages(multipage) as pages:
        for order, info in enumerate(FEATURES, 1):
            ptab = pairwise.loc[pairwise.Feature == info[0]]
            fig, ax = plt.subplots(figsize=(2.2, 2.25))
            fig.subplots_adjust(left=.22, right=.98, bottom=.14, top=.98)
            draw_panel(ax, data, info, ptab, 777 + order)
            stem = f"{order:02d}_{info[1].replace(' ', '_').replace('/', '-')}_E1-E5_Dunn-Holm"
            png, pdf = png_dir / f"{stem}.png", pdf_dir / f"{stem}.pdf"
            fig.savefig(png, dpi=900, facecolor="white")
            fig.savefig(pdf, facecolor="white"); pages.savefig(fig, facecolor="white")
            plt.close(fig)
            manifest.append({"Order": order, "Feature": info[0], "Display_name": info[1],
                             "Significant_Holm_pairs": int((ptab.P_adj_Holm_within_feature < .05).sum()),
                             "PNG": f"individual_png_900dpi/{png.name}",
                             "PDF": f"individual_pdf/{pdf.name}"})

    fig, axes = plt.subplots(6, 3, figsize=(7.2, 9.2))
    fig.subplots_adjust(left=.075, right=.99, bottom=.045, top=.99, wspace=.39, hspace=.25)
    for order, (ax, info) in enumerate(zip(axes.flat, FEATURES), 1):
        draw_panel(ax, data, info, pairwise.loc[pairwise.Feature == info[0]], 777 + order)
    fig.savefig(OUT / "Core18_E1-E5_raw_distributions_Dunn-Holm_overview_7p2x9p2in.png",
                dpi=900, facecolor="white")
    fig.savefig(OUT / "Core18_E1-E5_raw_distributions_Dunn-Holm_overview_7p2x9p2in.pdf",
                facecolor="white")
    plt.close(fig)
    return pd.DataFrame(manifest)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    data = load_data()
    descriptive = descriptive_statistics(data)
    omnibus, pairwise, variance = inferential_statistics(data)
    manifest = save_figures(data, pairwise)
    methods = pd.DataFrame([
        ["Cohort", "450 HC×GC consensus cells; E1=142, E2=160, E3=45, E4=61, E5=42; expert labels not used."],
        ["Data scale", "Raw values of the frozen 18-feature electrophysiology panel; no PCA values used for univariate plots or tests."],
        ["Omnibus test", "Two-sided Kruskal-Wallis test with tie correction; BH correction across 18 features."],
        ["Post hoc test", "All 10 two-sided Dunn comparisons with tie correction; Holm correction within each feature."],
        ["Pairwise effect size", "Dunn r = |Z|/sqrt(450)."],
        ["Variance test", "Brown-Forsythe median-centered Levene test; BH correction across 18 features."],
        ["Figure", "Raw observations, notched box plots, and exact Holm-adjusted P values for significant pairs only."],
    ], columns=["Item", "Details"])
    tables = {"Methods": methods, "Raw_data": data, "Group_descriptives": descriptive,
              "Omnibus_KW": omnibus, "Pairwise_Dunn_Holm": pairwise,
              "Variance_BrownForsythe": variance, "Figure_manifest": manifest}
    for name, table in tables.items():
        table.to_csv(OUT / f"Core18_Egroup_{name}.csv", index=False)
        (OUT / f"Core18_Egroup_{name}.json").write_text(
            table.to_json(orient="records", double_precision=15), encoding="utf-8"
        )
    summary = {"N": len(data), "group_counts": data.E_type.value_counts().sort_index().to_dict(),
               "features": len(FEATURES),
               "omnibus_BH_significant": int((omnibus.P_adj_BH_across_18_features < .05).sum()),
               "pairwise_Holm_significant": int(pairwise["Reject_Holm_0.05"].sum()),
               "variance_BH_significant": int((variance.P_adj_BH_across_18_features < .05).sum())}
    (OUT / "Core18_Egroup_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2)); print(OUT)


if __name__ == "__main__":
    main()
