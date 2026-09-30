from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import numpy as np
import pandas as pd
from scipy.stats import chi2

from analyze_nac_d1d2_transcriptomic_ephys_rrr import (
    correlation_loadings, fit_rrr, orient_axes, read_h5ad_subset, transcriptomic_pca,
)


ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "outputs/dSTR_dSTRvSTR_E_QC/Ca_Pu_NAC_final19_consensus_scan"
ASSIGN = ROOT / "outputs/R3_panels_v3/00_tuned_p80_ee12_random_seed777_coordinates.csv"
ZFILE = BASE / "Ca_Pu_NAC_final19_transformed_Zscore_matrix_n390.csv"
OUT = ROOT / "outputs/R3_RRR_fourpanel"
PAIR12 = (1, 0)  # Component 2 on X, Component 1 on Y.
PAIR13 = (2, 0)  # Component 3 on X, Component 1 on Y.
T_ORDER = ["D1", "D2"]
T_COLORS = {"D1": "#D95F02", "D2": "#008F7A"}
E_ORDER = ["C1", "C2", "C3", "C4"]
E_COLORS = {"C1": "#F8766D", "C2": "#7CAE00", "C3": "#00BFC4", "C4": "#C77CFF"}
E_LABELS = {
    "Epsy_width_rheo": "Width rheo",
    "Epsy_fast_trough_v_rheo": "Fast trough",
    "Epsy_peak_deltav_rheo": "Peak ΔV",
    "Epsy_peak_v_rheo": "Peak V",
    "Epsy_postap_slope_rheo": "Post-AP slope",
    "Epsy_threshold_v_rheo": "Threshold V",
    "Epsy_trough_t_rheo": "Trough t",
    "Epsy_trough_v_rheo": "Trough V",
    "Epsy_upstroke_downstroke_ratio_rheo": "Up/down ratio",
    "Epsy_ahp_delay_5spike": "AHP delay",
    "Epsy_ahp_delay_ratio_5spike": "AHP ratio",
    "Epsy_postap_slope_hero": "Post-AP hero",
    "Epsy_trough_t_hero": "Trough t hero",
    "Epsy_downstroke_adapt_ratio": "Downstroke adapt",
    "Epsy_peak_v_adapt_ratio": "Peak V adapt",
    "Epsy_threshold_v_adapt_ratio": "Threshold adapt",
    "Epsy_upstroke_adapt_ratio": "Upstroke adapt",
    "Epsy_width_adapt_ratio": "Width adapt",
    "Epsy_threshold_v_short_square": "Threshold SS",
}


def short(feature):
    return E_LABELS.get(feature, feature.removeprefix("Epsy_").replace("_", " "))


def spread(values, gap=0.17):
    if len(values) < 2:
        return values.copy()
    order = np.argsort(values)
    placed = values[order].copy()
    for i in range(1, len(placed)):
        placed[i] = max(placed[i], placed[i - 1] + gap)
    if placed[-1] > .92:
        placed -= placed[-1] - .92
    for i in range(len(placed) - 2, -1, -1):
        placed[i] = min(placed[i], placed[i + 1] - gap)
    if placed[0] < -.92:
        placed += -.92 - placed[0]
    result = np.empty_like(placed); result[order] = placed
    return result


def vectors(ax, loadings, names, top_n, component_pair, *, fontsize=2.30,
            gap=0.17, label_x=.98, force_side=None):
    pair = loadings[:, list(component_pair)]
    selected = np.argsort(np.linalg.norm(pair, axis=1))[::-1][:top_n]
    vec = pair[selected]
    labels = np.asarray(names, dtype=object)[selected]
    for v in vec:
        ax.annotate("", xy=v, xytext=(0, 0),
                    arrowprops=dict(arrowstyle="-|>", color="black", lw=.42,
                                    shrinkA=0, shrinkB=0, mutation_scale=3.6), zorder=4)
    side = (np.full(len(vec), float(force_side)) if force_side is not None
            else np.where(vec[:, 0] >= 0, 1., -1.))
    label_y = np.zeros(len(vec))
    for s in [-1., 1.]:
        mask = side == s
        label_y[mask] = spread(np.clip(vec[mask, 1], -.90, .90), gap=gap)
    for v, label, s, ly in zip(vec, labels, side, label_y):
        lx = s * label_x
        ax.plot([v[0], lx], [v[1], ly], color="#555555", lw=.22, zorder=3)
        ax.text(lx, ly, str(label), ha="left" if s > 0 else "right", va="center",
                fontsize=fontsize, color="#202020",
                bbox=dict(boxstyle="round,pad=.08", facecolor="white",
                          edgecolor="#666666", linewidth=.22), clip_on=False, zorder=5)


def ellipse(ax, points, color, coverage=.90):
    if len(points) < 4:
        return
    cov = np.cov(points.T); vals, vecs = np.linalg.eigh(cov)
    idx = np.argsort(vals)[::-1]; vals, vecs = vals[idx], vecs[:, idx]
    angle = np.degrees(np.arctan2(vecs[1, 0], vecs[0, 0]))
    radius = np.sqrt(chi2.ppf(coverage, 2)); centre = points.mean(axis=0)
    width, height = 2 * radius * np.sqrt(np.maximum(vals, 0))
    ax.add_patch(Ellipse(centre, width, height, angle=angle, facecolor=color,
                         edgecolor="none", alpha=.07, zorder=0))
    ax.add_patch(Ellipse(centre, width, height, angle=angle, facecolor="none",
                         edgecolor=color, lw=.42, alpha=.95, zorder=1))


def scaled(projection, component_pair):
    xy = projection[:, list(component_pair)]
    scale = max(float(np.quantile(np.linalg.norm(xy, axis=1), .99)), 1e-12)
    return xy / scale * .72


def panel(ax, projection, loadings, names, point_groups, point_order, point_colors,
          t_groups, top_n, component_pair, *, vector_fontsize=2.30,
          vector_gap=.17, vector_label_x=.98, vector_force_side=None,
          x_extent=1.19):
    xy = scaled(projection, component_pair)
    theta = np.linspace(0, 2 * np.pi, 361)
    ax.plot(np.cos(theta), np.sin(theta), color="#4E4E4E", lw=.42, zorder=-1)
    ax.axhline(0, color="#AFAFAF", lw=.22, zorder=-1)
    ax.axvline(0, color="#AFAFAF", lw=.22, zorder=-1)
    for group in T_ORDER:
        ellipse(ax, xy[t_groups == group], T_COLORS[group], .90)
    for group in point_order:
        mask = point_groups == group
        ax.scatter(xy[mask, 0], xy[mask, 1], s=2.15, c=point_colors[group],
                   edgecolors="white", linewidths=.12, alpha=.88, zorder=2)
    vectors(ax, loadings, names, top_n, component_pair,
            fontsize=vector_fontsize, gap=vector_gap,
            label_x=vector_label_x, force_side=vector_force_side)
    ax.set_xlim(-x_extent, x_extent); ax.set_ylim(-1.12, 1.12)
    ax.set_aspect("equal", adjustable="box"); ax.set_xticks([]); ax.set_yticks([])
    ax.patch.set_alpha(0)
    for spine in ax.spines.values():
        spine.set_visible(False)


def main():
    mpl.rcParams.update({"font.family": "Arial", "font.size": 4,
                         "pdf.fonttype": 42, "ps.fonttype": 42})
    OUT.mkdir(parents=True, exist_ok=True)
    assign = pd.read_csv(ASSIGN, dtype={"cell_label": str})
    z = pd.read_csv(ZFILE, dtype={"cell_label": str})
    features = [x for x in z.columns if x != "cell_label"]
    data = assign.loc[assign.Consensus.eq(True) & assign.T_class.isin(T_ORDER),
                      ["cell_label", "T_class", "HC_class"]].merge(
        z, on="cell_label", validate="one_to_one")
    if len(data) != 346 or len(features) != 19:
        raise RuntimeError(f"Expected n=346, p=19; got n={len(data)}, p={len(features)}")

    expression, all_genes = read_h5ad_subset(data.cell_label.tolist())
    t_scores, _, genes, variance, gene_matrix = transcriptomic_pca(expression, all_genes)
    e = data[features].to_numpy(float)
    xm, ym, _, v_raw = fit_rrr(t_scores, e, rank=3)
    v = orient_axes(v_raw, features)
    beta = np.linalg.pinv(t_scores - xm) @ (e - ym)
    t_projection = (t_scores - xm) @ beta @ v
    e_projection = (e - ym) @ v
    t_loadings = correlation_loadings(gene_matrix, t_projection)
    e_loadings = correlation_loadings(e, e_projection)
    e_names = np.array([short(x) for x in features], dtype=object)
    t_groups = data.T_class.to_numpy(str); e_groups = data.HC_class.to_numpy(str)

    fig, axes = plt.subplots(1, 4, figsize=(5.0, 1.0), gridspec_kw={"wspace": .14})
    panel(axes[0], t_projection, t_loadings, genes, t_groups, T_ORDER, T_COLORS,
          t_groups, 9, PAIR12)
    panel(axes[1], t_projection, t_loadings, genes, t_groups, T_ORDER, T_COLORS,
          t_groups, 9, PAIR13)
    panel(axes[2], e_projection, e_loadings, e_names, e_groups, E_ORDER, E_COLORS,
          t_groups, 10, PAIR12)
    panel(axes[3], e_projection, e_loadings, e_names, e_groups, E_ORDER, E_COLORS,
          t_groups, 10, PAIR13)
    fig.subplots_adjust(left=.006, right=.925, bottom=.035, top=.965)
    stem = OUT / "Macaque_CaPuNAC_RRR_T12_T13_E12_E13_compact_W5_H1"
    fig.savefig(stem.with_suffix(".png"), dpi=1200, facecolor="white")
    fig.savefig(stem.with_suffix(".pdf"), facecolor="white")
    plt.close(fig)

    pd.DataFrame({"cell_label": data.cell_label, "T_class": data.T_class,
                  "E_class": data.HC_class,
                  "T_Component1": t_projection[:, 0], "T_Component2": t_projection[:, 1],
                  "T_Component3": t_projection[:, 2], "E_Component1": e_projection[:, 0],
                  "E_Component2": e_projection[:, 1], "E_Component3": e_projection[:, 2]}).to_csv(
        OUT / "RRR_cell_scores_n346.csv", index=False)
    pd.DataFrame(e_loadings, index=features,
                 columns=["Component1", "Component2", "Component3"]).to_csv(
        OUT / "RRR_Efeature_correlation_loadings.csv")
    pd.DataFrame(t_loadings, index=genes,
                 columns=["Component1", "Component2", "Component3"]).to_csv(
        OUT / "RRR_gene_correlation_loadings.csv")
    pd.DataFrame({"PC": np.arange(1, len(variance)+1), "T_explained_variance_ratio": variance}).to_csv(
        OUT / "transcriptomic_PCA_variance.csv", index=False)
    print(stem.with_suffix(".png"))
    print(pd.crosstab(data.T_class, data.HC_class).to_string())


if __name__ == "__main__":
    main()
