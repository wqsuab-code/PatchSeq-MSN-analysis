#!/usr/bin/env python3
"""Frozen Macaque M1-M4 radar: existing ten features on the shared Z=-3..3 scale."""
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from PIL import Image
from scipy.interpolate import CubicSpline


BASE = Path(__file__).resolve().parent
RUN = BASE / "m18_tempfreeze_NPC5_HCK4_res2.3"
ASSIGN = RUN / "01_temp_frozen_assignments_126.csv"
ZSCORE = BASE / "m18_adaptive_pca126" / "02_transformed_z_117.csv"
SELECTION = RUN / "panels_A_I" / "radar_group_top10_feature_selection_117_consensus.csv"
LABELS = BASE / "morphology_feature_abbreviations.csv"
OUT = RUN / "panels_A_I" / "Zscore_top10_radar"

GROUPS = ["M1", "M2", "M3", "M4"]
EXPECTED = {"M1": 43, "M2": 42, "M3": 20, "M4": 12}
COLORS = {"M1": "#1F77B4", "M2": "#D9A400", "M3": "#8C564B", "M4": "#E377C2"}
Z_MIN, Z_MAX = -3.0, 3.0
CELL_LW, CELL_ALPHA = 0.35, 0.04
MEDIAN_LW, IQR_ALPHA = 0.64, 0.12
GRID_LW, GRID_ALPHA, FONT_SIZE = 0.28, 0.38, 4.0
FIGSIZE = (4.64, 1.00)


def display_radius(z):
    return np.clip((np.asarray(z, float) - Z_MIN) / (Z_MAX - Z_MIN), 0, 1)


def smooth(values, angles, smooth_angles):
    curve = CubicSpline(
        np.r_[angles, 2 * np.pi], np.r_[values, values[0]], bc_type="periodic"
    )(smooth_angles)
    return np.clip(curve, 0, 1)


def load_data():
    selected = pd.read_csv(SELECTION).sort_values("axis").copy()
    features = selected["feature"].tolist()
    if len(features) != 10 or len(set(features)) != 10:
        raise RuntimeError(f"Expected ten unique frozen radar features, got {features}")
    label_map = dict(pd.read_csv(LABELS).values)
    labels = [label_map[f] for f in features]

    assign = pd.read_csv(ASSIGN, dtype={"cell_label": str})
    assign = assign.loc[assign["concordant"].eq(True), ["cell_label", "HC_K4"]].copy()
    assign["M_class"] = "M" + assign["HC_K4"].astype(int).astype(str)
    z = pd.read_csv(ZSCORE, dtype={"cell_label": str})
    data = assign[["cell_label", "M_class"]].merge(
        z[["cell_label"] + features], on="cell_label", validate="one_to_one"
    )
    counts = data["M_class"].value_counts().reindex(GROUPS).to_dict()
    if len(data) != 117 or counts != EXPECTED:
        raise RuntimeError(f"Frozen cohort changed: n={len(data)}, counts={counts}")
    if data[features].isna().any().any():
        raise RuntimeError("Missing value in selected frozen morphology Z scores")
    return data, features, labels, selected


def draw_axis(ax, group, values, labels, angles, smooth_angles):
    color = COLORS[group]
    for row in display_radius(values):
        ax.plot(smooth_angles, smooth(row, angles, smooth_angles), color=color,
                lw=CELL_LW, alpha=CELL_ALPHA, zorder=1)

    q25, median, q75 = np.quantile(values, [0.25, 0.50, 0.75], axis=0)
    q25_curve = smooth(display_radius(q25), angles, smooth_angles)
    med_curve = smooth(display_radius(median), angles, smooth_angles)
    q75_curve = smooth(display_radius(q75), angles, smooth_angles)
    ax.fill_between(smooth_angles, q25_curve, q75_curve, color=color,
                    alpha=IQR_ALPHA, lw=0, zorder=2)
    ax.plot(smooth_angles, med_curve, color="#202020", lw=MEDIAN_LW, zorder=3)

    ax.set_theta_offset(np.pi / 2)
    ax.set_theta_direction(-1)
    ax.set_ylim(0, 1)
    ax.set_xticks(angles)
    ax.set_xticklabels(labels, fontsize=FONT_SIZE, color="#808080", alpha=0.70)
    ax.tick_params(axis="x", pad=-2.7, length=0)
    ax.set_yticks(np.linspace(0, 1, 7))
    ax.set_yticklabels([])
    ax.grid(color="#808080", lw=GRID_LW, alpha=GRID_ALPHA)
    ax.spines["polar"].set_color("#808080")
    ax.spines["polar"].set_alpha(0.45)
    ax.spines["polar"].set_linewidth(0.35)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    data, features, labels, selected = load_data()
    angles = np.linspace(0, 2 * np.pi, len(features), endpoint=False)
    smooth_angles = np.linspace(0, 2 * np.pi, 361)

    mpl.rcParams.update({
        "font.family": "Arial", "font.size": FONT_SIZE, "axes.titlesize": FONT_SIZE,
        "xtick.labelsize": FONT_SIZE, "ytick.labelsize": FONT_SIZE,
        "pdf.fonttype": 42, "ps.fonttype": 42, "savefig.pad_inches": 0,
    })
    fig = plt.figure(figsize=FIGSIZE, facecolor="white")
    centers = [0.125, 0.375, 0.625, 0.875]
    axes = [fig.add_axes([c - 0.071, 0.145, 0.142, 0.66], projection="polar") for c in centers]
    summaries = []

    for center, ax, group in zip(centers, axes, GROUPS):
        values = data.loc[data["M_class"].eq(group), features].to_numpy(float)
        draw_axis(ax, group, values, labels, angles, smooth_angles)
        q = np.quantile(values, [0.25, 0.50, 0.75], axis=0)
        for feature, q25, med, q75 in zip(features, q[0], q[1], q[2]):
            summaries.append({"M_class": group, "n": len(values), "feature": feature,
                              "q25_z": q25, "median_z": med, "q75_z": q75})
        fig.text(center, 0.975, f"{group} (n={len(values)})", ha="center", va="top",
                 fontsize=FONT_SIZE, color="#202020")

    stem = OUT / "M1-M4_top10_frozenZ_radar_Zminus3_to_plus3"
    png, pdf = stem.with_suffix(".png"), stem.with_suffix(".pdf")
    fig.savefig(png, dpi=900, facecolor="white")
    fig.savefig(pdf, facecolor="white")
    plt.close(fig)

    selected.insert(0, "display_order", np.arange(1, 11))
    selected["abbreviation"] = labels
    selected.to_csv(OUT / "top10_feature_selection_and_order.csv", index=False)
    data.to_csv(OUT / "top10_frozen_Zscore_cells_117.csv", index=False)
    pd.DataFrame(summaries).to_csv(OUT / "top10_class_median_IQR_Zscore.csv", index=False)

    with Image.open(png) as im:
        if im.size != (4176, 900):
            raise RuntimeError(f"Unexpected 900-dpi PNG size: {im.size}")
    lines = [
        "MACAQUE M1-M4 TOP-TEN FROZEN-Z RADAR",
        f"Cohort: 117 HC-GC consensus MSN cells; counts: {EXPECTED}.",
        f"Frozen Z-score source: {ZSCORE}",
        f"Frozen ten-feature selection/order source: {SELECTION}",
        "Feature order: " + "; ".join(f"{i+1}. {a} [{f}]" for i, (a, f) in enumerate(zip(labels, features))),
        "Shared display range: Z=-3 to +3; values outside this interval are clipped for display only.",
        "All panels use the same frozen cohort-wide Z scores; no within-class normalization is performed.",
        "Pale curves: individual cells; colored band: feature-wise IQR; black curve: feature-wise median.",
        "Periodic cubic splines produce closed curves. No median interior fill is used.",
        f"Canvas: {FIGSIZE[0]:.2f} x {FIGSIZE[1]:.2f} in; Arial 4 pt; PNG 900 dpi and vector PDF.",
        "No PCA, clustering, labels, feature selection, or Z scores were recalculated.",
    ]
    (OUT / "top10_frozenZ_radar_parameters.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(png)
    print(pdf)
    print(data["M_class"].value_counts().reindex(GROUPS).to_string())
    print("Features:", ", ".join(labels))


if __name__ == "__main__":
    main()
