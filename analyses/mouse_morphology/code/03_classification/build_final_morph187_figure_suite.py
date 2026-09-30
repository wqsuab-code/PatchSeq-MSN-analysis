#!/usr/bin/env python
"""Freeze and plot the final Morph n=187 NPC3/HC4/GC-res0.50 taxonomy."""

from __future__ import annotations

import json
import argparse
from pathlib import Path

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm, LinearSegmentedColormap, ListedColormap
from matplotlib.patches import Ellipse
import neurom
import numpy as np
import pandas as pd
from scipy.stats import chi2
from sklearn.ensemble import RandomForestClassifier
from sklearn.inspection import permutation_importance
from sklearn.manifold import TSNE
from sklearn.metrics import balanced_accuracy_score
from sklearn.model_selection import StratifiedGroupKFold


ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "outputs/morph_qc/morph187_NPC3-4_HCK4-5_GC_allres_raw"
RAW = ROOT / "outputs/morph_qc/morph187_after_incomplete_reconstruction_quarantine/01_morph187_discovery_input.csv"
T_AUDIT = ROOT / "outputs/morph_qc/all228_Morph_E_T_QC_eligibility_audit/01_all228_Morph_E_T_QC_status.csv"
ASC_AUDIT = ROOT / "outputs/morph_qc/NPC3_res2p50_HCK4_G1-G4_ASC_matrix_cell_galleries/02_cell_to_ASC_year_metric_mapping_audit_193cells.csv"
ASC_ROOT = Path("E:/ASC files")
OUT = ROOT / "outputs/morph_qc/final_morph187_NPC3_HCK4_GCres0p50_figures"
ID = "MSN_unique_ID"
FEATURES = [
    "M_soma_circularity_index", "M_soma_aspect_ratio", "M_cell_max_radial_dist",
    "M_total_number_of_neurites", "M_basal_dendrite_avg_tortuosity",
    "M_Total_neurite_length_(sections)", "M_Number_of_bifurcation_points",
    "M_Maximum_branch_order", "M_trunk_angle_min", "M_trunk_angle_max",
]
DISPLAY = [
    "Soma circularity", "Soma aspect ratio", "Max radial distance",
    "Primary neurite number", "Mean tortuosity", "Total neurite length",
    "Bifurcation points", "Maximum branch order", "Minimum trunk angle",
    "Maximum trunk angle",
]
M_LEVELS = ["M1", "M2", "M3", "M4"]
COLORS = {"M1": "#00468B", "M2": "#42B540", "M3": "#ED0000", "M4": "#0099B4"}
GC_MAP = {"S3": "G1", "S2": "G2", "S0": "G3", "S1": "G4"}
G_TO_M = {"G1": "M1", "G2": "M2", "G3": "M3", "G4": "M4"}
SEED = 20260826


def preprocess_train_test(train: np.ndarray, test: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    minimum = train.min(axis=0)
    shifted_train = np.maximum(train - minimum, 0)
    shifted_test = np.maximum(test - minimum, 0)
    total = shifted_train.sum(axis=0)
    total[total <= 0] = 1
    log_train = np.log1p(shifted_train / total * 10000)
    log_test = np.log1p(shifted_test / total * 10000)
    mean = log_train.mean(axis=0)
    sd = log_train.std(axis=0, ddof=1); sd[sd <= 0] = 1
    return (log_train - mean) / sd, (log_test - mean) / sd


def configure() -> None:
    mpl.rcParams.update({
        "font.family": "Arial", "font.size": 4, "axes.linewidth": 0.75,
        "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
        "savefig.facecolor": "white", "figure.facecolor": "white",
    })


def save(fig: plt.Figure, stem: str, dpi: int = 900) -> None:
    fig.savefig(OUT / f"{stem}.png", dpi=dpi, bbox_inches="tight", pad_inches=0.02)
    fig.savefig(OUT / f"{stem}.pdf", bbox_inches="tight", pad_inches=0.02)
    fig.savefig(OUT / f"{stem}.svg", bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)


def load_final() -> tuple[pd.DataFrame, np.ndarray, np.ndarray, np.ndarray]:
    raw = pd.read_csv(RAW)
    scores = pd.read_csv(SRC / "01_PCA_scores_187cells.csv")
    hc = pd.read_csv(SRC / "04_HC_NPC3-4_K4-5_assignments.csv")
    hc = hc[(hc.NPC == 3) & (hc.HC_K == 4)][[ID, "HC_cluster"]]
    gc = pd.read_csv(SRC / "05_GC_NPC3-4_allres_raw_assignments.csv")
    gc = gc[(gc.NPC == 3) & np.isclose(gc.Resolution, 0.50)][[ID, "GC_raw_cluster"]]
    data = raw.merge(scores[[ID] + [f"PC{i}" for i in range(1, 11)]], on=ID, validate="one_to_one")
    data = data.merge(hc, on=ID, validate="one_to_one").merge(gc, on=ID, validate="one_to_one")
    data["M_class"] = data["HC_cluster"].str.replace("HC", "M", regex=False)
    data["GC_matched_class"] = data["GC_raw_cluster"].map(GC_MAP)
    data["GC_as_M_class"] = data["GC_matched_class"].map(G_TO_M)
    data["HC_GC_consensus"] = data["M_class"].eq(data["GC_as_M_class"])
    data["Final_taxonomy_parameter"] = "NPC3_HCK4_GCres0.50"
    t = pd.read_csv(T_AUDIT, low_memory=False)
    tcols = [ID, "Final_consensus_CellType", "Final_CellType_stability", "stable_class"]
    data = data.merge(t[tcols], on=ID, how="left", validate="one_to_one")
    data["D1_D2"] = np.where(data["Final_consensus_CellType"].isin(["D1", "D2"]), data["Final_consensus_CellType"], "Unstable")
    data.to_csv(OUT / "00_final_morph187_cell_assignments.csv", index=False)

    x = raw[FEATURES].to_numpy(float)
    z, _ = preprocess_train_test(x, x)
    pcs = data[["PC1", "PC2", "PC3"]].to_numpy(float)
    tsne = TSNE(
        n_components=2, perplexity=25, learning_rate=50, max_iter=3000,
        init="pca", method="barnes_hut", angle=0.0, random_state=SEED,
    ).fit_transform(pcs)
    data["tSNE1"], data["tSNE2"] = tsne[:, 0], tsne[:, 1]
    data.to_csv(OUT / "01_final_morph187_assignments_with_tSNE.csv", index=False)
    return data, z, pcs, tsne


def add_cov_ellipse(ax, xy: np.ndarray, color: str, level: float = 0.85, lw: float = 0.65) -> None:
    if len(xy) < 3:
        return
    cov = np.cov(xy, rowvar=False); mean = xy.mean(axis=0)
    vals, vecs = np.linalg.eigh(cov); order = vals.argsort()[::-1]
    vals, vecs = vals[order], vecs[:, order]
    width, height = 2 * np.sqrt(np.maximum(vals, 0) * chi2.ppf(level, 2))
    angle = np.degrees(np.arctan2(vecs[1, 0], vecs[0, 0]))
    ax.add_patch(Ellipse(mean, width, height, angle=angle, fill=False, edgecolor=color, lw=lw))


def style_embedding(ax, xlabel="t-SNE 1", ylabel="t-SNE 2") -> None:
    ax.set_xlabel(xlabel, fontsize=4); ax.set_ylabel(ylabel, fontsize=4)
    ax.set_xticks([]); ax.set_yticks([])
    ax.spines[["top", "right"]].set_visible(False)
    ax.spines[["left", "bottom"]].set_linewidth(0.75)
    ax.set_aspect("equal", adjustable="datalim")


def plot_confusions(data: pd.DataFrame) -> None:
    matrix = pd.crosstab(data["M_class"], data["GC_matched_class"]).reindex(index=M_LEVELS, columns=["G1", "G2", "G3", "G4"], fill_value=0)
    matrix.to_csv(OUT / "02_final_HC_rows_GCmatched_columns_counts.csv")
    rng = np.random.default_rng(SEED)
    fig, ax = plt.subplots(figsize=(1.25, 1.25))
    for i, m in enumerate(M_LEVELS):
        for j, g in enumerate(["G1", "G2", "G3", "G4"]):
            n = int(matrix.loc[m, g])
            if n:
                sigma = .11 if i == j else .15
                ax.scatter(j + np.clip(rng.normal(0, sigma, n), -.31, .31),
                           i + np.clip(rng.normal(0, sigma, n), -.31, .31),
                           s=1.0, c="#111111", edgecolors="none", alpha=0.4)
    for b in np.arange(-.5, 4.5, 1):
        ax.plot([-.5, 3.5], [b, b], color="#222222", lw=.5)
        ax.plot([b, b], [-.5, 3.5], color="#222222", lw=.5)
    for i, c in enumerate(COLORS.values()):
        ax.add_patch(plt.Rectangle((-.73, i-.5), .14, 1, color=c, ec="none", clip_on=False))
        ax.add_patch(plt.Rectangle((i-.5, 3.57), 1, .14, color=c, ec="none", clip_on=False))
    ax.set_xlim(-.77, 3.55); ax.set_ylim(3.73, -.55); ax.set_aspect("equal")
    ax.set_xticks(range(4), ["G1", "G2", "G3", "G4"], fontsize=4)
    ax.set_yticks(range(4), M_LEVELS, fontsize=4)
    ax.tick_params(length=0, pad=1); ax.set_xlabel("GC class"); ax.set_ylabel("HC-derived M class")
    for s in ax.spines.values(): s.set_visible(False)
    save(fig, "02_HC_GC_final_dot_confusion_matrix_1p25in")

    fig, ax = plt.subplots(figsize=(1.55, 1.45))
    cmap = LinearSegmentedColormap.from_list("blue", ["#FFFFFF", "#9ECAE1", "#08519C"])
    im = ax.imshow(matrix.to_numpy(), cmap=cmap, vmin=0, vmax=matrix.to_numpy().max())
    for i in range(4):
        for j in range(4):
            v = int(matrix.iloc[i, j]); ax.text(j, i, str(v), ha="center", va="center", fontsize=5, color="white" if v >= 35 else "black")
    ax.set_xticks(range(4), ["G1", "G2", "G3", "G4"]); ax.set_yticks(range(4), M_LEVELS)
    ax.set_xlabel("GC class"); ax.set_ylabel("HC-derived M class"); ax.tick_params(length=0, pad=1)
    ax.set_xticks(np.arange(-.5, 4, 1), minor=True); ax.set_yticks(np.arange(-.5, 4, 1), minor=True)
    ax.grid(which="minor", color="#222222", lw=.5); ax.tick_params(which="minor", bottom=False, left=False)
    for s in ax.spines.values(): s.set_visible(False)
    cb = fig.colorbar(im, ax=ax, fraction=.05, pad=.03); cb.set_label("Cell count", fontsize=4); cb.ax.tick_params(labelsize=4)
    save(fig, "03_HC_GC_final_count_confusion_heatmap")


def plot_embeddings(data: pd.DataFrame, tsne: np.ndarray) -> None:
    fig, ax = plt.subplots(figsize=(2.2, 2.2))
    for m in M_LEVELS:
        idx = data.M_class.eq(m).to_numpy(); ax.scatter(tsne[idx, 0], tsne[idx, 1], s=10, c=COLORS[m], edgecolors="none", label=m)
    style_embedding(ax); ax.legend(frameon=False, fontsize=4, ncol=2, loc="best", handletextpad=.2, columnspacing=.6)
    save(fig, "04_final_M1-M4_tSNE_p25_lr50")

    fig, ax = plt.subplots(figsize=(2.2, 2.2))
    discord = ~data.HC_GC_consensus.to_numpy()
    ax.scatter(tsne[discord, 0], tsne[discord, 1], s=9, c="#BDBDBD", edgecolors="none")
    for m in M_LEVELS:
        idx = data.M_class.eq(m).to_numpy() & ~discord
        ax.scatter(tsne[idx, 0], tsne[idx, 1], s=10, c=COLORS[m], edgecolors="none")
        add_cov_ellipse(ax, tsne[idx], COLORS[m])
    style_embedding(ax)
    save(fig, "05_HC_GC_consensus181_tSNE_discordant_gray")

    plot_d1d2_distributions(data, tsne)
    plot_gc_merged_distribution(data, tsne)


def plot_d1d2_distributions(data: pd.DataFrame, tsne: np.ndarray) -> None:
    """Plot D1/D2 highlights inside fixed all-consensus M-class ellipses."""
    consensus = data.HC_GC_consensus.astype(bool).to_numpy()
    audit = []
    for label, number in [("D1", "06"), ("D2", "07")]:
        fig, ax = plt.subplots(figsize=(2.2, 2.2))
        selected = consensus & data.D1_D2.eq(label).to_numpy()
        background = ~selected
        ax.scatter(tsne[background, 0], tsne[background, 1], s=8, c="#D0D0D0", edgecolors="none")
        for m in M_LEVELS:
            idx = selected & data.M_class.eq(m).to_numpy()
            consensus_m = consensus & data.M_class.eq(m).to_numpy()
            ax.scatter(tsne[idx, 0], tsne[idx, 1], s=10, c=COLORS[m], edgecolors="none")
            add_cov_ellipse(ax, tsne[consensus_m], COLORS[m])
            audit.append({"D1_D2": label, "M_consensus_class": m, "Cell_N": int(idx.sum())})
        style_embedding(ax); ax.set_title(label, fontsize=5, pad=2)
        save(fig, f"{number}_{label}_distribution_in_final_M_tSNE")
    pd.DataFrame(audit).to_csv(OUT / "06-07_D1D2_HC_GC_consensus_counts_by_M_class.csv", index=False)


def plot_gc_merged_distribution(data: pd.DataFrame, tsne: np.ndarray) -> None:
    """Plot all final cells by merged GC class with GC-derived ellipses."""
    g_levels = ["G1", "G2", "G3", "G4"]
    fig, ax = plt.subplots(figsize=(2.2, 2.2))
    audit = []
    for g, m in zip(g_levels, M_LEVELS):
        idx = data.GC_matched_class.eq(g).to_numpy()
        ax.scatter(tsne[idx, 0], tsne[idx, 1], s=10, c=COLORS[m], edgecolors="none", label=g)
        add_cov_ellipse(ax, tsne[idx], COLORS[m])
        audit.append({"GC_merged_class": g, "Cell_N": int(idx.sum())})
    style_embedding(ax)
    ax.set_title("GC merged", fontsize=5, pad=2)
    ax.legend(frameon=False, fontsize=4, ncol=2, loc="best", handletextpad=.2, columnspacing=.6)
    save(fig, "07b_GC_merged_distribution_in_final_M_tSNE")
    pd.DataFrame(audit).to_csv(OUT / "07b_GC_merged_counts.csv", index=False)


def plot_heatmap(data: pd.DataFrame, z: np.ndarray) -> None:
    rng = np.random.default_rng(SEED)
    order = np.concatenate([rng.permutation(np.where(data.M_class.eq(m))[0]) for m in M_LEVELS])
    mat = np.clip(z[order].T, -2, 2)
    # Nine visually distinct Z-score bands.  A narrow black interval around
    # zero keeps the requested purple--black--yellow direction while avoiding
    # the muddy near-black compression of the previous continuous gradient.
    bounds = np.array([-2, -1.5, -1.0, -0.5, -0.2, 0.2, 0.5, 1.0, 1.5, 2.0001])
    cmap = ListedColormap([
        "#6F00FF", "#8F00D8", "#69108D", "#351044", "#050505",
        "#625A00", "#A99B00", "#E1CD00", "#FFF200",
    ], name="purple_black_yellow_banded")
    norm = BoundaryNorm(bounds, cmap.N, clip=True)
    fig, ax = plt.subplots(figsize=(7.0, 1.55))
    im = ax.imshow(mat, aspect="auto", cmap=cmap, norm=norm, interpolation="nearest")
    sizes = [int(data.M_class.eq(m).sum()) for m in M_LEVELS]; boundaries = np.cumsum(sizes)
    for b in boundaries[:-1]: ax.axvline(b-.5, color="white", lw=2.4)
    for y in np.arange(.5, len(DISPLAY), 1):
        ax.axhline(y, color="#FFFFFF", lw=.15, alpha=.28)
    start = 0
    for m, n in zip(M_LEVELS, sizes):
        ax.add_patch(plt.Rectangle((start-.5, -.9), n, .35, color=COLORS[m], ec="none", clip_on=False)); start += n
    ax.set_yticks(range(len(DISPLAY)), DISPLAY, fontsize=4); ax.yaxis.tick_right(); ax.tick_params(axis="y", length=0, pad=2)
    ax.set_xticks([]); [s.set_visible(False) for s in ax.spines.values()]
    cb = fig.colorbar(
        im, ax=ax, fraction=.015, pad=.025,
        boundaries=bounds, ticks=[-2, -1, -.5, 0, .5, 1, 2], spacing="proportional"
    )
    cb.set_label("Z score", fontsize=4); cb.ax.tick_params(labelsize=4, length=1.5, width=.45)
    save(fig, "08_final_M1-M4_10feature_Zscore_heatmap")
    pd.DataFrame(z[order], index=data.iloc[order][ID], columns=FEATURES).to_csv(OUT / "08_heatmap_ordered_zscore_matrix.csv")


def plot_biplot(data: pd.DataFrame) -> None:
    load = pd.read_csv(SRC / "03_PCA_loadings.csv").set_index("Feature")
    stable = data.D1_D2.isin(["D1", "D2"]).to_numpy()
    dcolors = {"D1": "#E41A1C", "D2": "#377EB8"}
    for xpc, ypc, number in [("PC2", "PC1", "09"), ("PC3", "PC1", "10")]:
        xy = data.loc[stable, [xpc, ypc]].to_numpy(float)
        radial = np.sqrt((xy**2).sum(axis=1)); factor = .78 / np.quantile(radial, .985); xy *= factor
        fig, ax = plt.subplots(figsize=(2.4, 2.4))
        theta = np.linspace(0, 2*np.pi, 361); ax.plot(np.cos(theta), np.sin(theta), color="#666666", lw=.75)
        ax.axhline(0, color="#BBBBBB", lw=.35); ax.axvline(0, color="#BBBBBB", lw=.35)
        labels = data.loc[stable, "D1_D2"].to_numpy()
        for d in ["D1", "D2"]:
            idx = labels == d; ax.scatter(xy[idx, 0], xy[idx, 1], s=4.5, c=dcolors[d], edgecolors="none", alpha=.75)
            add_cov_ellipse(ax, xy[idx], dcolors[d], .85, .55)
        vectors = load.loc[FEATURES, [xpc, ypc]].to_numpy(float)
        vectors *= .72 / np.sqrt((vectors**2).sum(axis=1)).max()
        for k, (vx, vy) in enumerate(vectors):
            ax.arrow(0, 0, vx, vy, width=.002, head_width=.025, head_length=.035, length_includes_head=True, color="#111111")
            ax.text(vx*1.08, vy*1.08, DISPLAY[k], fontsize=3.2, ha="center", va="center",
                    bbox=dict(fc="white", ec="none", alpha=.6, pad=.1))
        ax.set_xlim(-1.15, 1.15); ax.set_ylim(-1.15, 1.15); ax.set_aspect("equal")
        ax.set_xlabel(f"Component {xpc[-1]}"); ax.set_ylabel(f"Component {ypc[-1]}")
        ax.set_xticks([]); ax.set_yticks([]); ax.spines[["top", "right"]].set_visible(False)
        ax.set_title("Morphological space", fontsize=4, pad=2)
        save(fig, f"{number}_morphological_{xpc}_vs_{ypc}_D1D2_85CI")


def plot_d1d2_stacked(data: pd.DataFrame) -> None:
    stable = data[data.HC_GC_consensus & data.D1_D2.isin(["D1", "D2"])]
    count = pd.crosstab(stable.M_class, stable.D1_D2).reindex(index=M_LEVELS, columns=["D1", "D2"], fill_value=0)
    pct = count.div(count.sum(axis=1), axis=0) * 100
    fig, ax = plt.subplots(figsize=(1.65, 1.55))
    x = np.arange(4); ax.bar(x, pct.D1, color="#FF379B", width=.72, edgecolor="none", label="D1")
    ax.bar(x, pct.D2, bottom=pct.D1, color="#00EEB3", width=.72, edgecolor="none", label="D2")
    ax.set_xticks(x, M_LEVELS); ax.set_ylim(0, 100); ax.set_ylabel("Cells (%)")
    ax.spines[["top", "right"]].set_visible(False); ax.legend(frameon=False, fontsize=4, ncol=2, loc="upper center", bbox_to_anchor=(.5, 1.13))
    save(fig, "11_D1D2_percentage_by_final_M_class")
    count.to_csv(OUT / "11_D1D2_counts_by_final_M_class.csv")


def plot_importance(data: pd.DataFrame) -> None:
    x = data[FEATURES].to_numpy(float); y = data.M_class.to_numpy(str); groups = data.Recording_day_proxy.to_numpy(str)
    rows = []; fold_metrics = []
    for repeat in range(10):
        cv = StratifiedGroupKFold(n_splits=5, shuffle=True, random_state=SEED + repeat)
        for fold, (tr, te) in enumerate(cv.split(x, y, groups), 1):
            ztr, zte = preprocess_train_test(x[tr], x[te])
            model = RandomForestClassifier(n_estimators=500, max_features="sqrt", min_samples_leaf=2,
                                           class_weight="balanced_subsample", random_state=SEED + repeat*10 + fold, n_jobs=1)
            model.fit(ztr, y[tr]); pred = model.predict(zte)
            baseline = balanced_accuracy_score(y[te], pred)
            pi = permutation_importance(model, zte, y[te], scoring="balanced_accuracy", n_repeats=15,
                                        random_state=SEED + repeat*100 + fold, n_jobs=1)
            for feature, value in zip(FEATURES, pi.importances_mean):
                rows.append({"Repeat": repeat+1, "Fold": fold, "Feature": feature, "Importance_raw": value})
            fold_metrics.append({"Repeat": repeat+1, "Fold": fold, "Balanced_accuracy": baseline,
                                 "N_train": len(tr), "N_test": len(te), "Train_days": len(np.unique(groups[tr])), "Test_days": len(np.unique(groups[te]))})
    imp = pd.DataFrame(rows); means = imp.groupby("Feature").Importance_raw.mean().clip(lower=0)
    denom = max(float(means.max()), 1e-12); imp["Relative_importance"] = np.clip(imp.Importance_raw / denom, 0, 1.0)
    summary = imp.groupby("Feature").Relative_importance.agg(["mean", "std", "median"]).sort_values("mean", ascending=False)
    order = summary.index.tolist(); display_map = dict(zip(FEATURES, DISPLAY))
    fig, ax = plt.subplots(figsize=(3.25, 3.15)); rng = np.random.default_rng(SEED)
    for yi, feature in enumerate(order):
        vals = imp.loc[imp.Feature.eq(feature), "Relative_importance"].to_numpy()
        ax.hlines(yi, 0, 1, color="#D9D9D9", lw=.55, zorder=0)
        ax.scatter(vals, yi + rng.normal(0, .055, len(vals)), s=6, c="#969696", edgecolors="none", alpha=.8)
        ax.scatter(summary.loc[feature, "mean"], yi, s=24, c="#4C78A8", edgecolors="none", zorder=3)
    ax.set_yticks(range(len(order)), [display_map[f] for f in order]); ax.invert_yaxis()
    ax.set_xlim(0, 1.02); ax.set_xticks([0, .25, .5, .75, 1]); ax.set_xlabel("Relative permutation importance")
    ax.spines[["top", "right", "left"]].set_visible(False); ax.tick_params(axis="y", length=0)
    save(fig, "12_groupedCV_random_forest_feature_importance")
    imp.to_csv(OUT / "12_feature_importance_all_folds.csv", index=False)
    summary.reset_index().to_csv(OUT / "12_feature_importance_summary.csv", index=False)
    pd.DataFrame(fold_metrics).to_csv(OUT / "12_groupedCV_fold_performance.csv", index=False)


def asc_segments(path: Path):
    neuron = neurom.load_morphology(str(path)); soma = np.asarray(neuron.soma.points, float)
    center = soma[:, :2].mean(axis=0) if soma.ndim == 2 and len(soma) else np.asarray(neuron.soma.center[:2], float)
    seg = [np.asarray(s.points[:, :2], float) - center for s in neuron.sections if len(s.points) >= 2]
    soma_xy = soma[:, :2] - center if soma.ndim == 2 and len(soma) else np.empty((0, 2))
    return seg, soma_xy


def plot_representatives(data: pd.DataFrame, pcs: np.ndarray) -> None:
    if not ASC_ROOT.is_dir():
        return
    audit = pd.read_csv(ASC_AUDIT); audit = audit.set_index(ID)
    candidates = []
    for i, row in data.iterrows():
        if row[ID] not in audit.index: continue
        raw_path = str(audit.loc[row[ID], "Selected_ASC_path"])
        path = Path(raw_path.replace("D:\\", "E:\\")) if raw_path and raw_path.lower() != "nan" else None
        if path and path.is_file(): candidates.append((i, row[ID], row.M_class, path))
    selected = []
    for m in M_LEVELS:
        idx = np.where(data.M_class.eq(m))[0]; center = pcs[idx].mean(axis=0)
        pool = [(i, cid, cl, p, float(np.linalg.norm(pcs[i] - center))) for i, cid, cl, p in candidates if cl == m]
        selected.extend(sorted(pool, key=lambda v: (v[4], v[1]))[:5])
    morphs = []
    for i, cid, m, path, dist in selected:
        seg, soma = asc_segments(path); morphs.append({"i":i,"id":cid,"m":m,"path":path,"dist":dist,"seg":seg,"soma":soma})
    points = [p for q in morphs for p in q["seg"]] + [q["soma"] for q in morphs if len(q["soma"])]
    allp = np.vstack(points); hx=max(170, max(abs(allp[:,0].min()),abs(allp[:,0].max()))*1.08); hy=max(150,max(abs(allp[:,1].min()),abs(allp[:,1].max()))*1.08)
    for number, m in enumerate(M_LEVELS, 13):
        subset=[q for q in morphs if q["m"]==m]; fig,axes=plt.subplots(1,5,figsize=(8.0,2.25))
        for ax,q in zip(axes,subset):
            for p in q["seg"]: ax.plot(p[:,0],p[:,1],color=COLORS[m],lw=.42)
            if len(q["soma"])>=3: ax.fill(q["soma"][:,0],q["soma"][:,1],color=COLORS[m],ec="none")
            bx=-hx+.07*(2*hx); by=-hy+.08*(2*hy); ax.plot([bx,bx+50],[by,by],color="black",lw=.65)
            ax.set_xlim(-hx,hx);ax.set_ylim(-hy,hy);ax.set_aspect("equal");ax.axis("off");ax.set_title(q["id"],fontsize=4,pad=1)
        fig.suptitle(f"{m}: centroid-nearest ASC reconstructions",fontsize=6,y=.99)
        fig.subplots_adjust(left=.01,right=.99,bottom=.01,top=.92,wspace=.02)
        save(fig,f"{number:02d}_{m}_five_centroid_nearest_ASC_reconstructions_common_scale")
    pd.DataFrame([{k:q[k] for k in ["id","m","path","dist"]} for q in morphs]).to_csv(OUT/"13-16_representative_ASC_manifest.csv",index=False)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--d1d2-only", action="store_true", help="Redraw only consensus-filtered D1/D2 t-SNE and percentage panels")
    args = parser.parse_args()
    OUT.mkdir(parents=True, exist_ok=True); configure()
    if args.d1d2_only:
        data = pd.read_csv(OUT / "01_final_morph187_assignments_with_tSNE.csv")
        tsne = data[["tSNE1", "tSNE2"]].to_numpy(float)
        plot_d1d2_distributions(data, tsne); plot_gc_merged_distribution(data, tsne); plot_d1d2_stacked(data)
        print(json.dumps({"HC_GC_consensus_D1D2_n": int((data.HC_GC_consensus & data.D1_D2.isin(["D1", "D2"])).sum())}, indent=2))
        return
    data, z, pcs, tsne = load_final()
    plot_confusions(data); plot_embeddings(data, tsne); plot_heatmap(data, z)
    plot_biplot(data); plot_d1d2_stacked(data); plot_importance(data); plot_representatives(data, pcs)
    summary = {
        "cells": len(data), "parameter": "NPC3 HC K4 GC resolution 0.50",
        "HC_GC_consensus_n": int(data.HC_GC_consensus.sum()),
        "HC_GC_consensus_percent": float(100*data.HC_GC_consensus.mean()),
        "M_sizes": data.M_class.value_counts().reindex(M_LEVELS).astype(int).to_dict(),
        "importance_target": "HC-derived M1-M4; descriptive recovery, not independent biological validation",
    }
    (OUT/"00_summary.json").write_text(json.dumps(summary,indent=2),encoding="utf-8")
    print(json.dumps(summary,indent=2))


if __name__ == "__main__":
    main()
