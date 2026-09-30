#!/usr/bin/env python3
"""Export the four corrected Macaque E4 submission panels from cached tables."""

from pathlib import Path
import json

import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from PIL import Image

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "outputs" / "Macaque_E4_ML_submission_complete"
TABLES = BASE / "tables"
OUT = BASE / "individual_panels"
SUBSETS = ["All 19", "Remove top 2", "Remove top 5", "Bottom 3"]
HEIGHT_SCALE = 8 / 7
SUFFIX = "_Yx1p142857_H1p143in"


def style():
    mpl.rcParams.update({
        "font.family": "Arial", "font.size": 4, "axes.titlesize": 4,
        "axes.labelsize": 4, "xtick.labelsize": 4, "ytick.labelsize": 4,
        "axes.linewidth": .45, "xtick.major.width": .4, "ytick.major.width": .4,
        "xtick.major.size": 1.8, "ytick.major.size": 1.8,
        "pdf.fonttype": 42, "ps.fonttype": 42, "savefig.facecolor": "white",
    })


def save(fig, name):
    fig.savefig(OUT / f"{name}.png", dpi=900, facecolor="white")
    fig.savefig(OUT / f"{name}.pdf", facecolor="white")
    plt.close(fig)


def clean(ax):
    ax.spines[["top", "right"]].set_visible(False)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    style()
    summary = json.loads((BASE / "00_submission_completion_summary.json").read_text())
    perm = pd.read_csv(TABLES / "01_matched_ExtraTrees_within_donor_permutation_500.csv")
    importance = pd.read_csv(TABLES / "04_train_only_permutation_importance_summary.csv")
    ablation = pd.read_csv(TABLES / "05_nested_ablation_performance.csv")
    metadata = pd.read_csv(TABLES / "08_metadata_bias_corrected_permutation_tests.csv")

    fig, ax = plt.subplots(figsize=(1.05, HEIGHT_SCALE))
    sns.histplot(perm.mean_balanced_accuracy, bins=24, color="#BDBDBD",
                 edgecolor="white", linewidth=.2, ax=ax)
    ax.axvline(summary["observed_ET_grouped5fold_mean_BA"], color="#D55E00", lw=.7)
    ax.set(xlabel="Mean balanced accuracy", ylabel="Count", title="Matched ET permutation")
    clean(ax); fig.subplots_adjust(left=.29, right=.98, bottom=.25, top=.88)
    save(fig, "01_matched_ET_permutation" + SUFFIX)

    q = importance.sort_values("relative_importance").tail(10).copy()
    q["label"] = (q.feature.str.removeprefix("Epsy_").str.replace("upstroke_downstroke", "up/down")
                  .str.replace("adapt_ratio", "adapt.").str.replace("5spike", "5-spike")
                  .str.replace("_rheo", " (rheo)").str.replace("_", " "))
    fig, ax = plt.subplots(figsize=(1.62, HEIGHT_SCALE))
    ax.barh(q.label, q.relative_importance, color="#4C78A8", height=.66)
    ax.set(xlabel="Relative importance", title="Train-only importance")
    clean(ax); fig.subplots_adjust(left=.55, right=.98, bottom=.23, top=.88)
    save(fig, "02_train_only_permutation_importance" + SUFFIX)

    fig, ax = plt.subplots(figsize=(1.20, HEIGHT_SCALE))
    sns.boxplot(data=ablation, x="subset", y="balanced_accuracy", order=SUBSETS,
                color="#009E73", width=.55, showfliers=False, linewidth=.45, ax=ax)
    ax.axhline(.25, color="#777777", ls=":", lw=.45)
    ax.set(xlabel="", ylabel="Balanced accuracy", title="Nested feature ablation")
    ax.set_xticks(range(4), ["All 19", "−top 2", "−top 5", "Bottom 3"], rotation=25, ha="right")
    clean(ax); fig.subplots_adjust(left=.30, right=.98, bottom=.34, top=.88)
    save(fig, "03_nested_feature_ablation" + SUFFIX)

    q = metadata.copy()
    q["label"] = q.metadata.replace({"Lib_region_of_interest_label": "ROI",
                                      "T_class": "T class", "donor_label": "Donor"})
    q = q.sort_values("cramers_v_bias_corrected")
    fig, ax = plt.subplots(figsize=(1.40, HEIGHT_SCALE))
    ax.barh(q.label, q.cramers_v_bias_corrected, color="#8C8C8C", height=.62)
    for i, row in enumerate(q.itertuples()):
        ax.text(.425, i, f"P={row.empirical_p:.3g}",
                va="center", ha="right", fontsize=4)
    ax.set_xlim(0, .44)
    ax.set(xlabel="Bias-corrected Cramer's V", title="Metadata association")
    clean(ax); fig.subplots_adjust(left=.29, right=.98, bottom=.25, top=.88)
    save(fig, "04_metadata_association_corrected" + SUFFIX)

    panel_names = [
        "01_matched_ET_permutation" + SUFFIX + ".png",
        "02_train_only_permutation_importance" + SUFFIX + ".png",
        "03_nested_feature_ablation" + SUFFIX + ".png",
        "04_metadata_association_corrected" + SUFFIX + ".png",
    ]
    panels = [Image.open(OUT / name).convert("RGB") for name in panel_names]
    gap = 108  # 0.12 inch at 900 dpi
    canvas = Image.new("RGB", (sum(im.width for im in panels) + gap * 3,
                               round(900 * HEIGHT_SCALE)), "white")
    x = 0
    for im in panels:
        canvas.paste(im, (x, 0))
        x += im.width + gap
    canvas.save(OUT / ("00_four_corrected_panels_no_overlap" + SUFFIX + ".png"), dpi=(900, 900))
    canvas.save(OUT / ("00_four_corrected_panels_no_overlap" + SUFFIX + ".pdf"), resolution=900)
    for im in panels:
        im.close()

    print(OUT)


if __name__ == "__main__":
    main()
