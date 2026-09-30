#!/usr/bin/env python
"""E-analysis-style Morph feature atlas for the frozen final n=187 taxonomy."""

from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.cm import ScalarMappable
from matplotlib.patches import Ellipse
import numpy as np
import pandas as pd
from scipy.stats import chi2, kruskal, mannwhitneyu, skew
from statsmodels.stats.multitest import multipletests


ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "outputs/morph_qc/final_morph187_NPC3_HCK4_GCres0p50_figures"
INPUT = BASE / "01_final_morph187_assignments_with_tSNE.csv"
RAW22 = ROOT / "outputs/morph_qc/all22_shift_sum10000_log1p_zscore_215cells/01_all22_raw_215cells.csv"
OUT = BASE / "33_all22_distribution_and_tSNE_atlas"
OUT.mkdir(parents=True, exist_ok=True)

FEATURES = [
    ("M_soma_average_radius", "Soma average radius", "µm"),
    ("M_soma_maximal_radius", "Soma maximum radius", "µm"),
    ("M_soma_minimal_radius", "Soma minimum radius", "µm"),
    ("M_soma_max_pairwise_dist", "Soma maximum pairwise distance", "µm"),
    ("M_soma_perimeter", "Soma perimeter", "µm"),
    ("M_soma_area", "Soma area", "µm²"),
    ("M_soma_circularity_index", "Soma circularity", "index"),
    ("M_soma_shape_factor", "Soma shape factor", "index"),
    ("M_soma_aspect_ratio", "Soma aspect ratio", "ratio"),
    ("M_cell_max_radial_dist", "Maximum radial distance", "µm"),
    ("M_total_number_of_neurites", "Primary neurite number", "count"),
    ("M_basal_dendrite_longthest_path", "Longest dendrite path", "µm"),
    ("M_basal_dendrite_Nseg", "Basal dendrite segments", "count"),
    ("M_basal_dendrite_avg_tortuosity", "Mean tortuosity", "ratio"),
    ("M_basal_dendrite_max_tortuosity", "Maximum tortuosity", "ratio"),
    ("M_Total_neurite_length_(sections)", "Total neurite length", "µm"),
    ("M_Total_neurite_volume", "Total neurite volume", "µm³"),
    ("M_Total_neurite_area", "Total neurite area", "µm²"),
    ("M_Number_of_bifurcation_points", "Bifurcation points", "count"),
    ("M_Maximum_branch_order", "Maximum branch order", "order"),
    ("M_trunk_angle_min", "Minimum trunk angle", "rad"),
    ("M_trunk_angle_max", "Maximum trunk angle", "rad"),
]
CORE10 = {
    "M_soma_circularity_index", "M_soma_aspect_ratio", "M_cell_max_radial_dist",
    "M_total_number_of_neurites", "M_basal_dendrite_avg_tortuosity",
    "M_Total_neurite_length_(sections)", "M_Number_of_bifurcation_points",
    "M_Maximum_branch_order", "M_trunk_angle_min", "M_trunk_angle_max",
}
M_LEVELS = ["M1", "M2", "M3", "M4"]
COLORS = {"M1": "#00468B", "M2": "#42B540", "M3": "#ED0000", "M4": "#0099B4"}
CMAP = LinearSegmentedColormap.from_list("feature_z", ["#2166AC", "#F7F7F7", "#B2182B"])
NORM = Normalize(-2.5, 2.5, clip=True)
SEED = 20260828


def configure() -> None:
    mpl.rcParams.update({
        "font.family": "Arial", "font.size": 5.5,
        "axes.linewidth": 0.8, "pdf.fonttype": 42, "ps.fonttype": 42,
        "svg.fonttype": "none", "savefig.facecolor": "white",
        "figure.facecolor": "white",
    })


def preprocess_global(x: np.ndarray) -> np.ndarray:
    """Frozen analysis transform: shift -> 10,000 scaling -> log1p -> global z."""
    minimum = x.min(axis=0)
    shifted = np.maximum(x - minimum, 0)
    total = shifted.sum(axis=0)
    total[total <= 0] = 1
    logged = np.log1p(shifted / total * 10000)
    sd = logged.std(axis=0, ddof=1)
    sd[sd <= 0] = 1
    return (logged - logged.mean(axis=0)) / sd


def format_p(p: float) -> str:
    if p < 1e-4:
        return f"{p:.1e}"
    if p < 0.001:
        return f"{p:.4f}"
    return f"{p:.3f}"


def covariance_ellipse(ax: plt.Axes, xy: np.ndarray, color: str, level: float = 0.85) -> None:
    if len(xy) < 3:
        return
    cov = np.cov(xy, rowvar=False)
    center = xy.mean(axis=0)
    vals, vecs = np.linalg.eigh(cov)
    order = vals.argsort()[::-1]
    vals, vecs = vals[order], vecs[:, order]
    width, height = 2 * np.sqrt(np.maximum(vals, 0) * chi2.ppf(level, 2))
    angle = np.degrees(np.arctan2(vecs[1, 0], vecs[0, 0]))
    ax.add_patch(Ellipse(center, width, height, angle=angle, fill=False,
                         edgecolor=color, linewidth=0.8, alpha=1.0, zorder=3))


def ellipse_bounds(xy: np.ndarray, level: float = 0.85) -> tuple[float, float, float, float]:
    cov = np.cov(xy, rowvar=False); center = xy.mean(axis=0)
    vals, vecs = np.linalg.eigh(cov); order = vals.argsort()[::-1]
    vals, vecs = vals[order], vecs[:, order]
    a, b = np.sqrt(np.maximum(vals, 0) * chi2.ppf(level, 2))
    angle = np.arctan2(vecs[1, 0], vecs[0, 0])
    xr = np.sqrt((a * np.cos(angle)) ** 2 + (b * np.sin(angle)) ** 2)
    yr = np.sqrt((a * np.sin(angle)) ** 2 + (b * np.cos(angle)) ** 2)
    return center[0] - xr, center[0] + xr, center[1] - yr, center[1] + yr


def distribution_axis(ax: plt.Axes, data: pd.DataFrame, feature: str, label: str, unit: str,
                      rng: np.random.Generator, marker: str) -> tuple[float, list[dict]]:
    arrays = [data.loc[data.M_class.eq(m), feature].dropna().to_numpy(float) for m in M_LEVELS]
    positions = np.arange(1, 5)
    vp = ax.violinplot(arrays, positions=positions, widths=0.68, showextrema=False)
    for body, m in zip(vp["bodies"], M_LEVELS):
        body.set_facecolor(COLORS[m]); body.set_edgecolor(COLORS[m])
        body.set_alpha(0.22); body.set_linewidth(0.65)
    bp = ax.boxplot(arrays, positions=positions, widths=0.28, patch_artist=True,
                    showfliers=False, medianprops={"color": "#222222", "linewidth": 0.75},
                    whiskerprops={"color": "#333333", "linewidth": 0.65},
                    capprops={"color": "#333333", "linewidth": 0.65})
    for box, m in zip(bp["boxes"], M_LEVELS):
        box.set_facecolor("white"); box.set_edgecolor(COLORS[m]); box.set_linewidth(0.75)
    for pos, vals, m in zip(positions, arrays, M_LEVELS):
        jitter = rng.uniform(-0.15, 0.15, len(vals))
        ax.scatter(pos + jitter, vals, s=20.0, color=COLORS[m], alpha=0.82,
                   edgecolors="none", rasterized=False, zorder=2)

    kw_p = float(kruskal(*arrays).pvalue)
    raw = []
    pairs = []
    for i in range(4):
        for j in range(i + 1, 4):
            p = float(mannwhitneyu(arrays[i], arrays[j], alternative="two-sided").pvalue)
            raw.append(p); pairs.append((i, j, p))
    qvals = multipletests(raw, method="fdr_bh")[1]
    records = []
    significant = []
    for (i, j, p), q in zip(pairs, qvals):
        records.append({"Feature": feature, "Comparison": f"{M_LEVELS[i]} vs {M_LEVELS[j]}",
                        "Mann_Whitney_p": p, "BH_FDR_q": float(q), "Significant_q_lt_0.05": q < .05})
        if q < .05:
            significant.append((i, j, float(q)))

    ymin = min(map(np.min, arrays)); ymax = max(map(np.max, arrays)); span = max(ymax - ymin, 1e-12)
    top = ymax + span * (0.15 + 0.10 * len(significant))
    ax.set_ylim(ymin - span * 0.06, top)
    base = ymax + span * 0.10
    step = span * 0.095
    for k, (i, j, q) in enumerate(significant):
        y = base + k * step
        ax.plot([i + 1, i + 1, j + 1, j + 1], [y - step*.12, y, y, y - step*.12],
                color="#222222", lw=0.55, clip_on=False)
        ax.text((i + j + 2) / 2, y + step*.05, format_p(q), ha="center", va="bottom", fontsize=4.0)

    ax.set_title(f"{label} {marker}\nKruskal–Wallis P={format_p(kw_p)}", fontsize=5.5, pad=2.0)
    ax.set_ylabel(f"{label} ({unit})", fontsize=5.0, labelpad=1.8)
    ax.set_xticks(positions, M_LEVELS, fontsize=5.0)
    ax.tick_params(axis="both", labelsize=4.8, width=.65, length=2.0, pad=1.2)
    ax.spines[["top", "right"]].set_visible(False)
    return kw_p, records


def embedding_axis(ax: plt.Axes, data: pd.DataFrame, z: np.ndarray, feature_index: int,
                   xlim: tuple[float, float], ylim: tuple[float, float]) -> None:
    xy = data[["tSNE1", "tSNE2"]].to_numpy(float)
    order = np.argsort(np.abs(z[:, feature_index]))
    ax.scatter(xy[order, 0], xy[order, 1], c=z[order, feature_index], cmap=CMAP, norm=NORM,
               s=20.0, alpha=0.95, edgecolors="none", rasterized=False, zorder=2)
    for m in M_LEVELS:
        covariance_ellipse(ax, xy[data.M_class.eq(m).to_numpy()], COLORS[m], level=0.85)
    ax.set_xticks([]); ax.set_yticks([])
    ax.spines[["top", "right", "left", "bottom"]].set_visible(False)
    ax.set_xlim(*xlim); ax.set_ylim(*ylim)
    ax.set_aspect("equal", adjustable="box")


def main() -> None:
    configure()
    data = pd.read_csv(INPUT)
    if len(data) != 187:
        raise ValueError(f"Expected frozen n=187, found n={len(data)}")
    raw22 = pd.read_csv(RAW22)
    raw_features = [f[0] for f in FEATURES]
    data = data.drop(columns=[c for c in raw_features if c in data.columns]).merge(
        raw22[["MSN_unique_ID"] + raw_features], on="MSN_unique_ID", how="left", validate="one_to_one")
    if data[raw_features].isna().any().any():
        raise ValueError("Missing all-22 morphology values after merging the frozen n=187 cohort")
    x = data[[f[0] for f in FEATURES]].to_numpy(float)
    z = preprocess_global(x)
    transformed_skew = skew(z, axis=0, bias=False, nan_policy="omit")
    rng = np.random.default_rng(SEED)

    # Match the frozen t-SNE coordinates used in the final M1-M4 figure.  A single
    # common square viewing window prevents panel-specific rescaling.
    xy = data[["tSNE1", "tSNE2"]].to_numpy(float)
    bounds = [ellipse_bounds(xy[data.M_class.eq(m).to_numpy()]) for m in M_LEVELS]
    xmin = min([xy[:, 0].min()] + [b[0] for b in bounds]); xmax = max([xy[:, 0].max()] + [b[1] for b in bounds])
    ymin = min([xy[:, 1].min()] + [b[2] for b in bounds]); ymax = max([xy[:, 1].max()] + [b[3] for b in bounds])
    cx, cy = (xmin + xmax) / 2, (ymin + ymax) / 2
    half = max(xmax - xmin, ymax - ymin) * .56
    xlim, ylim = (cx - half, cx + half), (cy - half, cy + half)

    # Four feature modules per row: 4-4-4-4-4-2. Derive a compact canvas
    # width from the four occupied columns and remove inter-column whitespace.
    n_columns = 4
    module_width_in = 4.0
    fig = plt.figure(figsize=(n_columns * module_width_in + .5, 16.5))
    outer = fig.add_gridspec(6, 4, left=.010, right=.995, bottom=.070, top=.985,
                             wspace=0.00, hspace=.43)
    global_records = []
    pairwise_records = []
    for idx, (feature, label, unit) in enumerate(FEATURES):
        row, col = (idx // 4, idx % 4)
        inner = outer[row, col].subgridspec(1, 2, width_ratios=[.92, 1.15], wspace=.035)
        ax_dist = fig.add_subplot(inner[0, 0])
        ax_emb = fig.add_subplot(inner[0, 1])
        marker = "†" if feature in CORE10 else "‡"
        if abs(float(transformed_skew[idx])) >= 2:
            marker += "§"
        kw_p, rec = distribution_axis(ax_dist, data, feature, label, unit, rng, marker)
        embedding_axis(ax_emb, data, z, idx, xlim, ylim)
        global_records.append({"Feature": feature, "Display": label, "Kruskal_Wallis_p": kw_p})
        pairwise_records.extend(rec)

    cax = fig.add_axes([.80, .024, .16, .010])
    cb = fig.colorbar(ScalarMappable(norm=NORM, cmap=CMAP), cax=cax, orientation="horizontal")
    cb.set_ticks([-2.5, 0, 2.5]); cb.ax.tick_params(labelsize=4.8, length=2.0, width=.6, pad=1)
    cb.set_label("Feature Z-score after frozen transformation", fontsize=5.0, labelpad=1)
    fig.text(.020, .022, "† frozen 10-feature core   ‡ excluded from the frozen core/redundant   § |skewness| ≥ 2 after transformation; pairwise labels are BH-FDR-adjusted q values.",
             fontsize=4.8, ha="left", va="center")

    stem = OUT / "all22_M1-M4_distributions_and_featureZ_tSNE_444442_final187"
    fig.savefig(stem.with_suffix(".png"), dpi=600, bbox_inches="tight", pad_inches=.03)
    fig.savefig(stem.with_name(stem.name + "_900dpi").with_suffix(".png"), dpi=900,
                bbox_inches="tight", pad_inches=.03)
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight", pad_inches=.03)
    fig.savefig(stem.with_suffix(".svg"), bbox_inches="tight", pad_inches=.03)
    plt.close(fig)

    pd.DataFrame(global_records).to_csv(OUT / "all22_global_Kruskal_Wallis_tests.csv", index=False)
    pd.DataFrame(pairwise_records).to_csv(OUT / "all22_pairwise_Mann_Whitney_BH_FDR.csv", index=False)
    zdf = pd.DataFrame(z, columns=[f"Z_{f[0]}" for f in FEATURES])
    pd.concat([data[["MSN_unique_ID", "M_class", "tSNE1", "tSNE2"]].reset_index(drop=True), zdf], axis=1).to_csv(
        OUT / "all22_transformed_Z_and_tSNE_coordinates.csv", index=False)


if __name__ == "__main__":
    main()
