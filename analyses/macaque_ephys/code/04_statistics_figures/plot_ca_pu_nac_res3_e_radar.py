from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.interpolate import CubicSpline


ROOT = Path(__file__).resolve().parents[1]
ASSIGN = ROOT / "outputs/R3_panels_v3/00_tuned_p80_ee12_random_seed777_coordinates.csv"
RAW = ROOT / "outputs/dSTR_dSTRvSTR_E_QC/primary10_sensitivity_cohorts/dSTR_plus_vSTR_A_all_stage1QC_raw_with_NA.csv"
OUT = ROOT / "outputs/R3_panels_v3"

FEATURES = [
    "fast_trough_v_rheo",
    "peak_v_rheo",
    "postap_slope_rheo",
    "threshold_v_rheo",
    "upstroke_downstroke_ratio_rheo",
    "upstroke_rheo",
    "width_rheo_ms",
    "avg_rate_rheo",
    "latency_rheo",
    "rheobase_i",
]
LABELS = [
    "Fast trough",
    "Peak V",
    "Post-AP slope",
    "Threshold V",
    "Up/down ratio",
    "Upstroke",
    "AP width",
    "Avg rate",
    "Latency",
    "Rheobase",
]
ORDER = ["C1", "C2", "C3", "C4"]
COLORS = {"C1": "#F8766D", "C2": "#7CAE00", "C3": "#00BFC4", "C4": "#C77CFF"}


def robust_scale(values: pd.DataFrame):
    scaled = values.copy().astype(float)
    rows = []
    for feature in FEATURES:
        x = values[feature].to_numpy(float)
        lo, hi = np.quantile(x, [0.025, 0.975])
        if not np.isfinite(lo) or not np.isfinite(hi) or hi <= lo:
            scaled[feature] = 0.5
        else:
            scaled[feature] = np.clip((x - lo) / (hi - lo), 0, 1)
        rows.append({"feature": feature, "q2.5": lo, "q97.5": hi})
    return scaled, pd.DataFrame(rows)


def draw(joined: pd.DataFrame, with_title: bool):
    scaled, limits = robust_scale(joined[FEATURES])
    angles = np.linspace(0, 2 * np.pi, len(FEATURES), endpoint=False)
    closed_angles = np.r_[angles, 2 * np.pi]
    smooth_angles = np.linspace(0, 2 * np.pi, 361)

    def smooth(row):
        curve = CubicSpline(closed_angles, np.r_[row, row[0]], bc_type="periodic")(smooth_angles)
        return np.clip(curve, 0, 1)

    mpl.rcParams.update({
        "font.family": "Arial",
        "font.size": 4,
        "axes.titlesize": 4,
        "xtick.labelsize": 4,
        "ytick.labelsize": 4,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })
    fig = plt.figure(figsize=(4.64, 1.0), facecolor="white")
    centers = [0.125, 0.375, 0.625, 0.875]
    bottom = 0.145 if not with_title else 0.125
    height = 0.66 if not with_title else 0.61
    axes = [fig.add_axes([c - 0.071, bottom, 0.142, height], projection="polar") for c in centers]
    means = []
    grey = "#808080"

    for ax, group in zip(axes, ORDER):
        mask = joined["E_class"].eq(group).to_numpy()
        matrix = scaled.loc[mask, FEATURES].to_numpy(float)
        trace_alpha = min(0.075, 4.0 / len(matrix))
        for row in matrix:
            ax.plot(smooth_angles, smooth(row), color=COLORS[group], lw=0.20, alpha=trace_alpha, zorder=1)
        mean = matrix.mean(axis=0)
        mean_curve = smooth(mean)
        ax.fill(smooth_angles, mean_curve, color=COLORS[group], alpha=0.14, zorder=2)
        ax.plot(smooth_angles, mean_curve, color="black", lw=0.64, zorder=3)
        means.extend(
            {"E_class": group, "n": len(matrix), "feature": f, "mean_scaled": v}
            for f, v in zip(FEATURES, mean)
        )

        ax.set_theta_offset(np.pi / 2)
        ax.set_theta_direction(-1)
        ax.set_ylim(0, 1.02)
        ax.set_xticks(angles)
        ax.set_xticklabels(LABELS, fontsize=4, color=grey, alpha=0.70)
        ax.tick_params(axis="x", pad=-2.7, length=0)
        ax.set_yticks([0.25, 0.50, 0.75, 1.00])
        ax.set_yticklabels([])
        ax.grid(color=grey, lw=0.28, alpha=0.38)
        ax.spines["polar"].set_color(grey)
        ax.spines["polar"].set_alpha(0.45)
        ax.spines["polar"].set_linewidth(0.35)

    if with_title:
        fig.text(0.5, 0.965, "Graph-based clustering (npcs=3, merged K=4)",
                 ha="center", va="top", fontsize=4, color="#333333")

    suffix = "title" if with_title else "clean"
    stem = f"E_GC_res3_mergedK4_radar10_consensus368_{suffix}_W4p64_H1"
    png = OUT / f"{stem}.png"
    pdf = OUT / f"{stem}.pdf"
    fig.savefig(png, dpi=900, facecolor="white")
    fig.savefig(pdf, facecolor="white")
    plt.close(fig)
    pd.DataFrame(means).to_csv(OUT / f"{stem}_class_means.csv", index=False)
    limits.to_csv(OUT / f"{stem}_robust_scaling.csv", index=False)
    return png, pdf


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    assign = pd.read_csv(ASSIGN, dtype={"cell_label": str})
    raw = pd.read_csv(RAW, dtype={"cell_label": str})
    joined = assign.loc[assign["Consensus"].eq(True), ["cell_label", "HC_class"]].rename(
        columns={"HC_class": "E_class"}
    ).merge(raw[["cell_label"] + FEATURES], on="cell_label", how="inner", validate="one_to_one")
    if len(joined) != 368:
        raise RuntimeError(f"Expected 368 HC-GC consensus cells, got {len(joined)}")
    if joined[FEATURES].isna().any().any():
        raise RuntimeError("Canonical radar features contain missing values")
    expected = {"C1": 59, "C2": 133, "C3": 57, "C4": 119}
    if joined["E_class"].value_counts().to_dict() != expected:
        raise RuntimeError(f"Unexpected class counts: {joined['E_class'].value_counts().to_dict()}")
    outputs = draw(joined, with_title=True) + draw(joined, with_title=False)
    for path in outputs:
        print(path)


if __name__ == "__main__":
    main()
