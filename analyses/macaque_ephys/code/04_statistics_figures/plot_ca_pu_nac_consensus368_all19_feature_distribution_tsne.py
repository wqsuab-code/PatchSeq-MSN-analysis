#!/usr/bin/env python3
"""Draw all 19 frozen macaque E features as raw distributions plus feature t-SNE maps."""

from itertools import combinations
from math import erfc, sqrt
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.patches import Ellipse
from matplotlib.ticker import MaxNLocator
import numpy as np
import pandas as pd
from scipy.stats import chi2, kruskal, rankdata
from sklearn.covariance import MinCovDet


ROOT = Path(__file__).resolve().parents[1]
ASSIGN = ROOT / "outputs/R3_panels_v3/00_tuned_p80_ee12_random_seed777_coordinates.csv"
RAW = ROOT / "outputs/R3_E_QC_overview_redrawn/final19_raw_complete_n390.csv"
ZMAT = ROOT / "outputs/dSTR_dSTRvSTR_E_QC/Ca_Pu_NAC_final19_consensus_scan/Ca_Pu_NAC_final19_transformed_Zscore_matrix_n390.csv"
OUT = ROOT / "outputs/R3_E_feature_distribution_tsne_all19"

ORDER = ["C1", "C2", "C3", "C4"]
COLORS = {"C1": "#F8766D", "C2": "#7CAE00", "C3": "#00BFC4", "C4": "#C77CFF"}

# feature, display title, unit, raw-value display multiplier
FEATURES = [
    ("Epsy_width_rheo", "AP width", "ms", 1000.0),
    ("Epsy_fast_trough_v_rheo", "Fast trough", "mV", 1.0),
    ("Epsy_peak_deltav_rheo", "Peak delta V", "mV", 1.0),
    ("Epsy_peak_v_rheo", "Peak V", "mV", 1.0),
    ("Epsy_postap_slope_rheo", "Post-AP slope", "mV/ms", 1.0),
    ("Epsy_threshold_v_rheo", "Threshold V", "mV", 1.0),
    ("Epsy_trough_t_rheo", "Trough time", "s", 1.0),
    ("Epsy_trough_v_rheo", "Trough V", "mV", 1.0),
    ("Epsy_upstroke_downstroke_ratio_rheo", "Up/down ratio", "", 1.0),
    ("Epsy_ahp_delay_5spike", "AHP delay, 5-spike", "ms", 1000.0),
    ("Epsy_ahp_delay_ratio_5spike", "AHP delay ratio, 5-spike", "", 1.0),
    ("Epsy_postap_slope_hero", "Post-AP slope, hero", "mV/ms", 1.0),
    ("Epsy_trough_t_hero", "Trough time, hero", "s", 1.0),
    ("Epsy_downstroke_adapt_ratio", "Downstroke adapt ratio", "", 1.0),
    ("Epsy_peak_v_adapt_ratio", "Peak V adapt ratio", "", 1.0),
    ("Epsy_threshold_v_adapt_ratio", "Threshold adapt ratio", "", 1.0),
    ("Epsy_upstroke_adapt_ratio", "Upstroke adapt ratio", "", 1.0),
    ("Epsy_width_adapt_ratio", "Width adapt ratio", "", 1.0),
    ("Epsy_threshold_v_short_square", "Threshold V, short square", "mV", 1.0),
]

ZLIM = 2.5
TSNE_X_DISPLAY_SCALE = 0.75
MCD_SUPPORT_FRACTION = 0.75
MCD_COVERAGE = 0.80
RANDOM_STATE = 777
ZCMAP = LinearSegmentedColormap.from_list(
    "feature_zscore", ["#2166AC", "#C9DCEB", "#F3F3F3", "#E5B9BA", "#B2182B"], N=256
)


def holm_adjust(values):
    p = np.asarray(values, dtype=float)
    idx = np.argsort(p)
    out = np.empty_like(p)
    running = 0.0
    for j, i in enumerate(idx):
        running = max(running, min(1.0, (len(p) - j) * p[i]))
        out[i] = running
    return out


def bh_adjust(values):
    p = np.asarray(values, dtype=float)
    idx = np.argsort(p)[::-1]
    out = np.empty_like(p)
    running = 1.0
    m = len(p)
    for reverse_pos, i in enumerate(idx, start=1):
        rank = m - reverse_pos + 1
        running = min(running, min(1.0, p[i] * m / rank))
        out[i] = running
    return out


def dunn_holm(values, groups):
    values = np.asarray(values, dtype=float)
    groups = np.asarray(groups, dtype=str)
    ranks = rankdata(values, method="average")
    n = len(values)
    _, ties = np.unique(values, return_counts=True)
    variance = n * (n + 1) / 12 - np.sum(ties**3 - ties) / (12 * (n - 1))
    means = {g: ranks[groups == g].mean() for g in ORDER}
    sizes = {g: int(np.sum(groups == g)) for g in ORDER}
    rows = []
    for g1, g2 in combinations(ORDER, 2):
        se = sqrt(variance * (1 / sizes[g1] + 1 / sizes[g2]))
        z = (means[g1] - means[g2]) / se
        rows.append({"Group_1": g1, "Group_2": g2, "N_1": sizes[g1], "N_2": sizes[g2],
                     "Dunn_Z": z, "Dunn_r": abs(z) / sqrt(n), "P_raw": erfc(abs(z) / sqrt(2))})
    out = pd.DataFrame(rows)
    out["P_adj_Holm_within_feature"] = holm_adjust(out["P_raw"])
    return out


def density_halfwidth(values, maximum=0.27):
    values = np.asarray(values, dtype=float)
    n = max(len(values), 2)
    sd = np.std(values, ddof=1) if len(values) > 1 else 0.0
    span = np.ptp(values)
    bw = 1.06 * sd * n ** (-0.2)
    if not np.isfinite(bw) or bw <= 1e-12:
        bw = max(span * 0.08, 1e-6)
    scaled = (values[:, None] - values[None, :]) / bw
    density = np.exp(-0.5 * scaled**2).mean(axis=1)
    return np.clip(maximum * density / max(density.max(), 1e-12), maximum * 0.04, maximum)


def p_label(p):
    if p < 1e-3:
        return f"{p:.1e}"
    if p < 0.01:
        return f"{p:.3f}"
    return f"{p:.2f}"


def add_distribution(ax, data, feature, title, unit, multiplier, pairwise, seed):
    rng = np.random.default_rng(seed)
    plotted = data[feature] * multiplier
    for x, group in enumerate(ORDER, start=1):
        vals = plotted[data["E_class"].eq(group)].to_numpy(float)
        jitter = rng.uniform(-1, 1, len(vals)) * density_halfwidth(vals)
        ax.scatter(x + jitter, vals, s=1.75, color=COLORS[group], alpha=0.58,
                   linewidths=0, rasterized=True, zorder=3)
        ax.boxplot([vals], positions=[x], widths=0.42, whis=(2.5, 97.5), showfliers=False,
                   patch_artist=True, manage_ticks=False,
                   boxprops={"facecolor": "none", "edgecolor": "#333333", "linewidth": 0.32},
                   medianprops={"color": "#111111", "linewidth": 0.42},
                   whiskerprops={"color": "#555555", "linewidth": 0.32},
                   capprops={"color": "#555555", "linewidth": 0.32}, zorder=5)

    vals = plotted.to_numpy(float)
    ymin, ymax = np.nanmin(vals), np.nanmax(vals)
    span = max(ymax - ymin, abs(ymax) * 0.05, 1e-6)
    sig = pairwise[pairwise["P_adj_Holm_within_feature"] < 0.05].nsmallest(
        4, "P_adj_Holm_within_feature"
    )
    ax.set_ylim(ymin - 0.08 * span, ymax + span * (0.23 + 0.105 * len(sig)))
    for level, row in enumerate(sig.itertuples(index=False)):
        x1, x2 = ORDER.index(row.Group_1) + 1, ORDER.index(row.Group_2) + 1
        y = ymax + span * (0.13 + 0.10 * level)
        tick = span * 0.024
        ax.plot([x1, x1, x2, x2], [y - tick, y, y, y - tick], color="#333333", lw=0.32,
                clip_on=False)
        ax.text((x1 + x2) / 2, y + span * 0.012, p_label(row.P_adj_Holm_within_feature),
                ha="center", va="bottom", fontsize=2.85, color="#222222")

    ax.set_title(f"{title} ({unit})" if unit else title, fontsize=4, pad=1.2)
    ax.set_xticks(range(1, 5), ORDER)
    ax.yaxis.set_major_locator(MaxNLocator(5))
    ax.tick_params(axis="both", labelsize=3.0, length=1.35, width=0.32, pad=0.6)
    ax.spines[["top", "right"]].set_visible(False)
    ax.spines[["left", "bottom"]].set_linewidth(0.36)


def ellipse_specs(data):
    specs = []
    for group in ORDER:
        pts = data.loc[data["E_class"].eq(group), ["plot_tSNE1", "plot_tSNE2"]].to_numpy(float)
        mcd = MinCovDet(support_fraction=MCD_SUPPORT_FRACTION,
                        random_state=RANDOM_STATE).fit(pts)
        center = mcd.location_
        cov = mcd.covariance_
        values, vectors = np.linalg.eigh(cov)
        idx = values.argsort()[::-1]
        values, vectors = values[idx], vectors[:, idx]
        width, height = 2 * np.sqrt(np.maximum(values, 0) * chi2.ppf(MCD_COVERAGE, 2))
        angle = np.degrees(np.arctan2(vectors[1, 0], vectors[0, 0]))
        specs.append((group, center, width, height, angle))
    return specs


def shared_tsne_limits(data, specs):
    xmin, xmax = data.plot_tSNE1.min(), data.plot_tSNE1.max()
    ymin, ymax = data.plot_tSNE2.min(), data.plot_tSNE2.max()
    for _, c, w, h, a in specs:
        t = np.deg2rad(a)
        xr = 0.5 * np.sqrt((w * np.cos(t))**2 + (h * np.sin(t))**2)
        yr = 0.5 * np.sqrt((w * np.sin(t))**2 + (h * np.cos(t))**2)
        xmin, xmax = min(xmin, c[0] - xr), max(xmax, c[0] + xr)
        ymin, ymax = min(ymin, c[1] - yr), max(ymax, c[1] + yr)
    span = max(xmax - xmin, ymax - ymin)
    cx, cy = (xmin + xmax) / 2, (ymin + ymax) / 2
    half = span * 0.53
    return (cx - half, cx + half), (cy - half, cy + half)


def add_tsne(ax, data, feature, specs, xlim, ylim):
    ax.scatter(data.plot_tSNE1, data.plot_tSNE2, c=data[f"z__{feature}"], cmap=ZCMAP,
               norm=Normalize(-ZLIM, ZLIM), s=2.6, alpha=0.88,
               linewidths=0, rasterized=True, zorder=2)
    for group, center, width, height, angle in specs:
        ax.add_patch(Ellipse(center, width, height, angle=angle, fill=False,
                             edgecolor=COLORS[group], linewidth=0.56,
                             linestyle=(0, (2.2, 1.2)), zorder=4))
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_aspect("equal", adjustable="box")
    ax.axis("off")


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    mpl.rcParams.update({"font.family": "Arial", "font.size": 4,
                         "pdf.fonttype": 42, "ps.fonttype": 42})
    assign = pd.read_csv(ASSIGN, dtype={"cell_label": str})
    raw = pd.read_csv(RAW, dtype={"cell_label": str})
    z = pd.read_csv(ZMAT, dtype={"cell_label": str})
    names = [f[0] for f in FEATURES]
    data = assign.loc[assign.Consensus.eq(True),
                      ["cell_label", "HC_class", "T_class", "tSNE1", "tSNE2"]].rename(
                          columns={"HC_class": "E_class"})
    data = data.merge(raw[["cell_label", "Lib_region_of_interest_label"] + names],
                      on="cell_label", validate="one_to_one")
    data = data.merge(z[["cell_label"] + names].rename(columns={f: f"z__{f}" for f in names}),
                      on="cell_label", validate="one_to_one")
    if len(data) != 368 or data[names].isna().any().any():
        raise RuntimeError(f"Expected 368 complete consensus cells; got n={len(data)}")

    # Preserve the frozen numerical t-SNE coordinates and apply only the same
    # 75% horizontal display compression used in the R3 main consensus panel.
    data["plot_tSNE1"] = data["tSNE1"] * TSNE_X_DISPLAY_SCALE
    data["plot_tSNE2"] = data["tSNE2"]

    omnibus_rows, pairwise_frames = [], []
    for feature, title, unit, multiplier in FEATURES:
        arrays = [(data.loc[data.E_class.eq(g), feature] * multiplier).to_numpy() for g in ORDER]
        h, p = kruskal(*arrays)
        omnibus_rows.append({"Feature": feature, "Display": title, "Unit": unit,
                              "Kruskal_H": h, "df": 3, "P_raw": p})
        pw = dunn_holm((data[feature] * multiplier).to_numpy(), data.E_class.to_numpy())
        pw.insert(0, "Feature", feature)
        pw.insert(1, "Display", title)
        pairwise_frames.append(pw)
    omnibus = pd.DataFrame(omnibus_rows)
    omnibus["Q_BH_across_19_features"] = bh_adjust(omnibus.P_raw)
    pairwise = pd.concat(pairwise_frames, ignore_index=True)

    specs = ellipse_specs(data)
    xlim, ylim = shared_tsne_limits(data, specs)
    fig = plt.figure(figsize=(7.2, 6.35), facecolor="white")
    outer = fig.add_gridspec(4, 5, left=0.020, right=0.988, bottom=0.055, top=0.982,
                             wspace=0.055, hspace=0.30)
    for i, (feature, title, unit, multiplier) in enumerate(FEATURES):
        row, col = divmod(i, 5)
        inner = outer[row, col].subgridspec(1, 2, width_ratios=[0.94, 1.06], wspace=0.012)
        ax1 = fig.add_subplot(inner[0, 0])
        ax2 = fig.add_subplot(inner[0, 1])
        add_distribution(ax1, data, feature, title, unit, multiplier,
                         pairwise[pairwise.Feature.eq(feature)], 20260908 + i)
        add_tsne(ax2, data, feature, specs, xlim, ylim)

    # The twentieth slot is deliberately reserved for the common color key.
    holder = fig.add_subplot(outer[3, 4])
    holder.axis("off")
    cax = holder.inset_axes([0.10, 0.43, 0.80, 0.045])
    cb = mpl.colorbar.ColorbarBase(cax, cmap=ZCMAP, norm=Normalize(-ZLIM, ZLIM),
                                   orientation="horizontal", ticks=[-2.5, 0, 2.5])
    cb.ax.tick_params(labelsize=3.1, length=1.2, width=0.32, pad=0.5)
    cb.outline.set_linewidth(0.32)
    cb.set_label("Feature Z-score", fontsize=3.4, labelpad=0.2)

    stem = OUT / "Macaque_CaPuNAC_E4_all19_featureTSNE_mainParams_4x5"
    fig.savefig(stem.with_suffix(".png"), dpi=600, facecolor="white")
    fig.savefig(stem.with_suffix(".pdf"), facecolor="white")
    fig.savefig(stem.with_suffix(".svg"), facecolor="white")
    plt.close(fig)

    omnibus.to_csv(OUT / "all19_KruskalWallis_BH.csv", index=False)
    pairwise.to_csv(OUT / "all19_Dunn_Holm.csv", index=False)
    data.to_csv(OUT / "all19_plot_data_consensus368.csv", index=False)
    pd.DataFrame(FEATURES, columns=["Feature", "Display", "Unit", "Raw_display_multiplier"]).to_csv(
        OUT / "all19_feature_order_and_units.csv", index=False)
    (OUT / "main_tsne_display_parameters.txt").write_text(
        "Frozen t-SNE coordinates: PC1-PC3; perplexity=80; early_exaggeration=12; "
        "learning_rate=auto; init=random; max_iter=3000; random_state=777.\n"
        "Display: tSNE1 multiplied by 0.75; tSNE2 unchanged; identical shared limits "
        "for all 19 panels.\n"
        "Class envelopes: robust Minimum Covariance Determinant; support_fraction=0.75; "
        "random_state=777; chi-square coverage=80%.\n"
        "Cells: the same 368 frozen HC-GC consensus cells used in the main figure.\n",
        encoding="utf-8",
    )
    print(stem.with_suffix(".png"))
    print(omnibus[["Display", "Kruskal_H", "Q_BH_across_19_features"]].to_string(index=False))


if __name__ == "__main__":
    main()
