from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import skew


ROOT = Path(r"C:\Users\53461\OneDrive\Documentos\Patch-seq_Mouse_Acb_MSN_T-type_Visualization")
INPUT = ROOT / "outputs/morph_qc/final_morph187_NPC3_HCK4_GCres0p50_figures/00_final_morph187_cell_assignments.csv"
OUT = ROOT / "outputs/morph_qc/final_morph187_NPC3_HCK4_GCres0p50_figures/29_frozen10_histogram_distributions"


FEATURES = [
    ("M_soma_circularity_index", "Soma circularity", "Index", False),
    ("M_soma_aspect_ratio", "Soma aspect ratio", "Ratio", False),
    ("M_cell_max_radial_dist", "Max radial distance", "Distance (µm)", False),
    ("M_total_number_of_neurites", "Primary neurite number", "Count", True),
    ("M_basal_dendrite_avg_tortuosity", "Mean tortuosity", "Index", False),
    ("M_Total_neurite_length_(sections)", "Total neurite length", "Length (µm)", False),
    ("M_Number_of_bifurcation_points", "Bifurcation points", "Count", True),
    ("M_Maximum_branch_order", "Maximum branch order", "Order", True),
    ("M_trunk_angle_min", "Minimum trunk angle", "Angle (rad)", False),
    ("M_trunk_angle_max", "Maximum trunk angle", "Angle (rad)", False),
]


def fd_bins(values: np.ndarray) -> np.ndarray | int:
    values = np.asarray(values, dtype=float)
    q25, q75 = np.quantile(values, [0.25, 0.75])
    width = 2 * (q75 - q25) / np.cbrt(values.size)
    if not np.isfinite(width) or width <= 0:
        return 12
    n_bins = int(np.ceil((values.max() - values.min()) / width))
    return int(np.clip(n_bins, 8, 28))


def bins_for(values: np.ndarray, discrete: bool) -> np.ndarray | int:
    if not discrete:
        return fd_bins(values)
    lo = int(np.floor(values.min()))
    hi = int(np.ceil(values.max()))
    return np.arange(lo - 0.5, hi + 1.5, 1.0)


def style_axis(ax: plt.Axes, values: np.ndarray, title: str, xlabel: str, discrete: bool) -> None:
    ax.hist(
        values,
        bins=bins_for(values, discrete),
        color="#00468B",
        alpha=0.78,
        edgecolor="white",
        linewidth=0.35,
    )
    median = float(np.median(values))
    ax.axvline(median, color="#ED0000", linewidth=0.7, linestyle=(0, (3, 2)), zorder=3)
    ax.set_title(title, loc="left", pad=2.5, fontweight="normal")
    ax.set_xlabel(xlabel, labelpad=1.5)
    ax.set_ylabel("Cells", labelpad=1.5)
    ax.text(
        0.98,
        0.94,
        f"n={values.size}\nmedian={median:.2f}",
        transform=ax.transAxes,
        ha="right",
        va="top",
        color="#333333",
        fontsize=3.5,
    )
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_linewidth(0.65)
    ax.spines["bottom"].set_linewidth(0.65)
    ax.tick_params(axis="both", which="major", width=0.55, length=2.0, pad=1.2)
    ax.grid(axis="y", color="#D9D9D9", linewidth=0.35, alpha=0.7)
    ax.set_axisbelow(True)
    if discrete:
        unique = np.unique(values.astype(int))
        step = max(1, int(np.ceil(unique.size / 7)))
        ax.set_xticks(unique[::step])


def save_figure(fig: plt.Figure, stem: Path) -> None:
    fig.savefig(stem.with_suffix(".png"), dpi=600, bbox_inches="tight", facecolor="white")
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight", facecolor="white")
    fig.savefig(stem.with_suffix(".svg"), bbox_inches="tight", facecolor="white")


def style_compact_axis(
    ax: plt.Axes,
    values: np.ndarray,
    title: str,
    xlabel: str,
    show_ylabel: bool,
) -> None:
    ax.hist(
        values,
        bins=20,
        color="#9E9E9E",
        alpha=1.0,
        edgecolor="white",
        linewidth=0.2,
    )
    max_count = max((patch.get_height() for patch in ax.patches), default=1.0)
    tick_step = max(1, int(np.ceil((max_count * 1.06) / 4.0)))
    y_top = tick_step * 4
    ax.set_ylim(0, y_top)
    ax.set_yticks(np.arange(0, y_top + tick_step, tick_step))
    ax.axvline(
        float(np.median(values)),
        color="#ED0000",
        linewidth=0.45,
        linestyle=(0, (3, 2)),
        zorder=3,
    )
    ax.set_title(title, loc="center", pad=1.2, fontweight="normal", fontsize=3.2)
    ax.set_xlabel(xlabel, labelpad=0.8, fontsize=3.2)
    ax.set_ylabel("Cells" if show_ylabel else "", labelpad=0.8, fontsize=3.2)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_linewidth(0.45)
    ax.spines["bottom"].set_linewidth(0.45)
    ax.tick_params(axis="both", which="major", width=0.4, length=1.2, pad=0.5, labelsize=2.8)
    ax.grid(axis="y", color="#E1E1E1", linewidth=0.25, alpha=0.8)
    ax.set_axisbelow(True)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(INPUT)
    missing_columns = [name for name, *_ in FEATURES if name not in df.columns]
    if missing_columns:
        raise KeyError(f"Missing frozen feature columns: {missing_columns}")

    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "DejaVu Sans"],
            "font.size": 4.0,
            "axes.titlesize": 5.0,
            "axes.labelsize": 4.0,
            "xtick.labelsize": 3.6,
            "ytick.labelsize": 3.6,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
        }
    )

    summaries: list[dict[str, float | str | int]] = []
    for feature, label, unit, discrete in FEATURES:
        values = pd.to_numeric(df[feature], errors="coerce").dropna().to_numpy(float)
        summaries.append(
            {
                "Feature": feature,
                "Display_name": label,
                "Unit": unit,
                "N": values.size,
                "Missing": int(df.shape[0] - values.size),
                "Min": float(np.min(values)),
                "Q1": float(np.quantile(values, 0.25)),
                "Median": float(np.median(values)),
                "Mean": float(np.mean(values)),
                "Q3": float(np.quantile(values, 0.75)),
                "Max": float(np.max(values)),
                "Skewness_raw": float(skew(values, bias=False)),
            }
        )

        fig, ax = plt.subplots(figsize=(3.3, 2.2), constrained_layout=True)
        style_axis(ax, values, label, unit, discrete)
        save_figure(fig, OUT / f"individual_{FEATURES.index((feature, label, unit, discrete)) + 1:02d}_{feature}")
        plt.close(fig)

    fig, axes = plt.subplots(5, 2, figsize=(6.9, 8.2))
    for ax, (feature, label, unit, discrete) in zip(axes.flat, FEATURES):
        values = pd.to_numeric(df[feature], errors="coerce").dropna().to_numpy(float)
        style_axis(ax, values, label, unit, discrete)
    fig.subplots_adjust(left=0.09, right=0.99, bottom=0.055, top=0.985, hspace=0.62, wspace=0.30)
    save_figure(fig, OUT / "combined_frozen10_raw_histograms_final187")
    plt.close(fig)

    # Requested compact format: 5 panels per row, 2 rows; each panel's allotted
    # canvas is exactly 0.8 in wide by 1.0 in high (overall 4 x 2 in).
    fig, axes = plt.subplots(2, 5, figsize=(4.0, 2.0))
    for index, (ax, (feature, label, unit, _)) in enumerate(zip(axes.flat, FEATURES)):
        values = pd.to_numeric(df[feature], errors="coerce").dropna().to_numpy(float)
        style_compact_axis(ax, values, label, unit, show_ylabel=index % 5 == 0)
    fig.subplots_adjust(left=0.075, right=0.995, bottom=0.15, top=0.955, hspace=0.62, wspace=0.42)
    compact_stem = OUT / "combined_frozen10_gray20bin_2x5_panels0p8x1in_final187"
    fig.savefig(compact_stem.with_suffix(".png"), dpi=900, facecolor="white")
    fig.savefig(compact_stem.with_suffix(".pdf"), facecolor="white")
    fig.savefig(compact_stem.with_suffix(".svg"), facecolor="white")
    plt.close(fig)

    pd.DataFrame(summaries).to_csv(OUT / "frozen10_raw_distribution_summary_final187.csv", index=False)


if __name__ == "__main__":
    main()
