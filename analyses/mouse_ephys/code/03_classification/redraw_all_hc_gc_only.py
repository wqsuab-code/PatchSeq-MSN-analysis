#!/usr/bin/env python3
"""Regenerate the electrophysiology figure set using HC and GC labels only.

Frozen analysis choices are preserved:
  * graph-based clustering (GC): NPCs=3, Seurat resolution=1.5, merged K=5
  * hierarchical clustering (HC): NPCs=3, Ward.D2, K=5
  * primary E-type cohort: cells for which aligned HC and GC labels agree

No expert label is read, copied, used for filtering, or written by this script.
"""

from __future__ import annotations

import json
import math
from pathlib import Path
import re
import sys

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.patches import Rectangle
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import leaves_list, linkage
from scipy.stats import kruskal
from sklearn.decomposition import PCA
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from plot_all25_violin_plus_feature_tsne_posthoc import dunn_holm  # noqa: E402


LABEL_FILE = ROOT / "outputs/e_type_qc/active_cohort/Active_Frozen_triple_consensus_assignments.csv"
RAW25_FILE = ROOT / "outputs/e_type_qc/parallel_preprocessing/Ephys_QCpass_494_raw_25features.csv"
PCA18_FILE = ROOT / "outputs/e_type_qc/active_cohort/Active_Final18_PCA_input.csv"
SVM_OOF_FILE = ROOT / "outputs/e_type_qc/HC_GC_only/16_manuscript_ml_validation_figures/Figure2_RBF-SVM_grouped_OOF.csv"
D1D2_FILE = (
    ROOT
    / "outputs/e_type_qc/fixed_npc3_res1.5_merged_k5_tsne"
    / "D1_vs_D2_within_Eclass/D1_vs_D2_Within_Eclass_Analysis_Data_451.csv"
)
OUT = ROOT / "outputs/e_type_qc/HC_GC_only"

E_LEVELS = ["E1", "E2", "E3", "E4", "E5"]
E_COLORS = {
    "E1": "#4E79A7",
    "E2": "#F28E2B",
    "E3": "#59A14F",
    "E4": "#2AA6B8",
    "E5": "#B07AA1",
}
GREY = "#C7C7C7"
ZLIM = 2.5
ZCMAP = LinearSegmentedColormap.from_list(
    "feature_zscore_refined",
    [
        (0.00, "#2166AC"),
        (0.44, "#A9C7DD"),
        (0.485, "#D9D9D9"),
        (0.515, "#D9D9D9"),
        (0.56, "#D8A8AA"),
        (1.00, "#B2182B"),
    ],
    N=256,
)

FEATURES = [
    ("E_AP.amplitude..mV.", "AP amplitude", "mV"),
    ("E_Input.resistance..MOhm.", "Input resistance", "MOhm"),
    ("E_ISI.adaptation.index", "ISI adaptation", "index"),
    ("E_ISI.coefficient.of.variation", "ISI CV", "CV"),
    ("E_Rheobase..pA.", "Rheobase", "pA"),
    ("E_Latency....20pA.current..ms.", "Latency (−20 pA)", "ms"),
    ("E_AP.amplitude.adaptation.index", "AP amplitude adaptation", "index"),
    ("E_Membrane.time.constant..ms.", "Membrane tau", "ms"),
    ("E_AP.coefficient.of.variation", "AP CV", "CV"),
    ("E_AP.threshold..mV.", "AP threshold", "mV"),
    ("E_Sag.ratio", "Sag ratio", "ratio"),
    ("E_Afterhyperpolarization..mV.", "AHP", "mV"),
    ("E_Upstroke.to.downstroke.ratio", "Upstroke/downstroke", "ratio"),
    ("E_Max.number.of.APs", "Max APs", "count"),
    ("E_Sag.time..s.", "Sag time", "s"),
    ("E_Latency..ms.", "Latency", "ms"),
    ("E_Holding.MP..mV.", "Holding MP", "mV"),
    ("E_AP.width..ms.", "AP width", "ms"),
    ("E_AP.amplitude.average.adaptation.index", "AP amplitude adaptation avg", "index"),
    ("E_Fitted.MP..mV.", "Fitted Vm", "mV"),
    ("E_AP.Fano.factor", "AP Fano", "Fano"),
    ("E_Sag.area..mV.s.", "Sag area", "mV·s"),
    ("E_Rebound..mV.", "Rebound", "mV"),
    ("E_ISI.average.adaptation.index", "ISI adaptation avg", "index"),
    ("E_ISI.Fano.factor", "ISI Fano", "Fano"),
]


def style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 6,
            "axes.linewidth": 0.5,
            "axes.labelsize": 6,
            "axes.titlesize": 7,
            "xtick.labelsize": 5,
            "ytick.labelsize": 5,
            "xtick.major.size": 1.5,
            "ytick.major.size": 1.5,
            "xtick.major.width": 0.5,
            "ytick.major.width": 0.5,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def zscore(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    sd = x.std(ddof=1)
    return np.zeros_like(x) if not np.isfinite(sd) or sd == 0 else (x - x.mean()) / sd


def safe(text: str) -> str:
    return re.sub(r"[^A-Za-z0-9]+", "_", text).strip("_")


def save_all(fig: mpl.figure.Figure, stem: Path, dpi: int = 600, tif: bool = True) -> None:
    stem.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(stem.with_suffix(".png"), dpi=dpi, facecolor="white", bbox_inches="tight")
    fig.savefig(stem.with_suffix(".pdf"), facecolor="white", bbox_inches="tight")
    if tif:
        fig.savefig(
            stem.with_suffix(".tif"),
            dpi=dpi,
            facecolor="white",
            bbox_inches="tight",
            pil_kwargs={"compression": "tiff_lzw"},
        )


def load() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    labels = pd.read_csv(LABEL_FILE)
    required = {"MSN_unique_ID", "tSNE_1", "tSNE_2", "HC_E", "Seurat_E"}
    if not required.issubset(labels.columns):
        raise ValueError(f"Missing label columns: {sorted(required - set(labels.columns))}")
    labels = labels[list(required)].copy()
    labels = labels.rename(columns={"Seurat_E": "GC_E"})
    labels["HC_GC_consensus"] = labels["HC_E"].eq(labels["GC_E"])
    labels["HC_GC_consensus_E"] = np.where(labels["HC_GC_consensus"], labels["HC_E"], pd.NA)
    labels["HC_E"] = pd.Categorical(labels["HC_E"], E_LEVELS, ordered=True)
    labels["GC_E"] = pd.Categorical(labels["GC_E"], E_LEVELS, ordered=True)
    labels["HC_GC_consensus_E"] = pd.Categorical(labels["HC_GC_consensus_E"], E_LEVELS, ordered=True)
    svm = pd.read_csv(SVM_OOF_FILE, usecols=["MSN_unique_ID", "Predicted_E"])
    labels = labels.merge(svm, on="MSN_unique_ID", how="left", validate="one_to_one")
    labels["SVM_disagree"] = (
        labels["HC_GC_consensus"]
        & labels["Predicted_E"].notna()
        & labels["Predicted_E"].ne(labels["HC_GC_consensus_E"].astype("object"))
    )

    raw = pd.read_csv(RAW25_FILE)
    pca = pd.read_csv(PCA18_FILE)
    feature_names = [x[0] for x in FEATURES]
    missing = [x for x in feature_names if x not in raw.columns]
    if missing:
        raise ValueError(f"Missing raw features: {missing}")
    data = labels.merge(raw[["MSN_unique_ID", *feature_names]], on="MSN_unique_ID", how="left", validate="one_to_one")
    if data[feature_names].isna().any().any():
        bad = data[feature_names].isna().sum()
        raise ValueError(f"Missing raw values: {bad[bad > 0].to_dict()}")
    pca_data = labels.merge(pca, on="MSN_unique_ID", how="inner", validate="one_to_one")
    return labels, data, pca_data


def classification_metrics(labels: pd.DataFrame) -> dict[str, float | int]:
    h = labels["HC_E"].astype(str)
    g = labels["GC_E"].astype(str)
    agree = h.eq(g)
    return {
        "N_all": int(len(labels)),
        "N_HC_GC_consensus": int(agree.sum()),
        "N_discordant": int((~agree).sum()),
        "Agreement_fraction": float(agree.mean()),
        "Agreement_percent": float(100 * agree.mean()),
        "ARI": float(adjusted_rand_score(h, g)),
        "NMI": float(normalized_mutual_info_score(h, g)),
    }


def draw_tsne(ax, labels: pd.DataFrame, column: str, title: str, consensus_only: bool = False) -> None:
    if consensus_only:
        discord = ~labels["HC_GC_consensus"]
        ax.scatter(labels.loc[discord, "tSNE_1"], labels.loc[discord, "tSNE_2"], s=3.0, c=GREY, alpha=0.85, linewidths=0)
        use = labels.loc[~discord]
        key = "HC_GC_consensus_E"
    else:
        use = labels
        key = column
    for e in E_LEVELS:
        q = use[key].astype(str).eq(e)
        ax.scatter(use.loc[q, "tSNE_1"], use.loc[q, "tSNE_2"], s=3.0, c=E_COLORS[e], alpha=0.95, linewidths=0)
    # Mark consensus-labelled cells not recovered by grouped RBF-SVM OOF prediction.
    outline = labels["SVM_disagree"]
    ax.scatter(
        labels.loc[outline, "tSNE_1"], labels.loc[outline, "tSNE_2"],
        s=3.0, facecolors="none", edgecolors="black", linewidths=0.25,
        alpha=1.0, zorder=5,
    )
    ax.set_title(title, fontweight="bold", pad=2)
    ax.set_xlabel("t-SNE 1")
    ax.set_ylabel("t-SNE 2")
    ax.set_aspect("equal", adjustable="datalim")
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(length=1.5, width=0.5)


def plot_classification(labels: pd.DataFrame, metrics: dict) -> None:
    out = OUT / "01_classification"
    fig, axes = plt.subplots(1, 3, figsize=(5.1, 1.7), constrained_layout=True)
    draw_tsne(axes[0], labels, "GC_E", "GC K=5")
    draw_tsne(axes[1], labels, "HC_E", "Ward.D2 HC K=5")
    draw_tsne(
        axes[2], labels, "HC_GC_consensus_E",
        f"HC×GC consensus\n{metrics['N_HC_GC_consensus']} colored | {metrics['N_discordant']} gray",
        consensus_only=True,
    )
    save_all(fig, out / "HC_GC_tSNE_three_views_W5p1_H1p7")
    plt.close(fig)

    counts = pd.crosstab(labels["HC_E"], labels["GC_E"]).reindex(index=E_LEVELS, columns=E_LEVELS, fill_value=0)
    rowpct = counts.div(counts.sum(axis=1), axis=0) * 100
    counts.to_csv(out / "HC_GC_confusion_counts.csv")
    rowpct.to_csv(out / "HC_GC_confusion_row_percent.csv")
    fig, ax = plt.subplots(figsize=(2.2, 2.0))
    im = ax.imshow(rowpct.to_numpy(), cmap="Blues", vmin=0, vmax=100, aspect="equal")
    for i in range(5):
        for j in range(5):
            color = "white" if rowpct.iloc[i, j] >= 55 else "#222222"
            ax.text(j, i, f"{counts.iloc[i,j]}\n{rowpct.iloc[i,j]:.0f}%", ha="center", va="center", fontsize=5, color=color)
    ax.set_xticks(range(5), E_LEVELS)
    ax.set_yticks(range(5), E_LEVELS)
    ax.set_xlabel("GC classification")
    ax.set_ylabel("HC classification")
    ax.set_title(
        f"HC × GC agreement: {metrics['Agreement_percent']:.2f}%\nARI {metrics['ARI']:.3f} | NMI {metrics['NMI']:.3f}",
        fontsize=7,
    )
    for x in np.arange(-0.5, 5, 1):
        ax.axhline(x, color="#D0D0D0", lw=0.4)
        ax.axvline(x, color="#D0D0D0", lw=0.4)
    cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
    cb.set_label("Within-HC row (%)", fontsize=5)
    cb.ax.tick_params(labelsize=4, length=1)
    save_all(fig, out / "HC_GC_confusion_matrix_W2p2_H2p0")
    plt.close(fig)


def plot_pca(pca_data: pd.DataFrame) -> None:
    out = OUT / "02_PCA"
    out.mkdir(parents=True, exist_ok=True)
    cols = [c for c in pca_data.columns if c.startswith("E_")]
    X = pca_data[cols].to_numpy(float)
    pc = PCA(n_components=min(10, X.shape[1]), random_state=777).fit_transform(X)
    model = PCA(n_components=min(10, X.shape[1]), random_state=777).fit(X)
    score = pca_data[["MSN_unique_ID", "HC_E", "GC_E", "HC_GC_consensus", "HC_GC_consensus_E"]].copy()
    score["PC1"], score["PC2"] = pc[:, 0], pc[:, 1]
    score.to_csv(out / "HC_GC_PCA_scores.csv", index=False)
    pd.DataFrame(model.components_.T, index=cols, columns=[f"PC{i+1}" for i in range(model.n_components_)]).to_csv(out / "HC_GC_PCA_loadings.csv")
    pd.DataFrame({"PC": np.arange(1, model.n_components_+1), "Explained_fraction": model.explained_variance_ratio_, "Cumulative_fraction": np.cumsum(model.explained_variance_ratio_)}).to_csv(out / "HC_GC_PCA_variance.csv", index=False)
    fig, axes = plt.subplots(1, 3, figsize=(5.1, 1.7), constrained_layout=True)
    temp = pca_data.copy()
    temp["tSNE_1"], temp["tSNE_2"] = pc[:, 0], pc[:, 1]
    old_x, old_y = "tSNE_1", "tSNE_2"
    for ax, key, title, co in zip(axes, ["GC_E", "HC_E", "HC_GC_consensus_E"], ["GC K=5", "Ward.D2 HC K=5", "HC×GC consensus"], [False, False, True]):
        draw_tsne(ax, temp, key, title, consensus_only=co)
        ax.set_xlabel(f"PC1 ({model.explained_variance_ratio_[0]*100:.1f}%)")
        ax.set_ylabel(f"PC2 ({model.explained_variance_ratio_[1]*100:.1f}%)")
    save_all(fig, out / "HC_GC_PCA_PC1_PC2_three_views_W5p1_H1p7")
    plt.close(fig)


def feature_statistics(cons: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    desc_rows, omni_rows, pair_rows = [], [], []
    groups = cons["E_type"].astype(str).to_numpy()
    for order, (feature, title, unit) in enumerate(FEATURES, 1):
        for e in E_LEVELS:
            x = cons.loc[cons["E_type"].astype(str).eq(e), feature].to_numpy(float)
            desc_rows.append({"Order": order, "Feature": feature, "Title": title, "Unit": unit, "E_type": e, "N": len(x), "Mean": x.mean(), "SD": x.std(ddof=1), "SEM": x.std(ddof=1)/math.sqrt(len(x)), "Median": np.median(x), "Q1": np.quantile(x,.25), "Q3": np.quantile(x,.75), "Min": x.min(), "Max": x.max()})
        arrays = [cons.loc[cons["E_type"].astype(str).eq(e), feature].to_numpy(float) for e in E_LEVELS]
        kw = kruskal(*arrays)
        omni_rows.append({"Order": order, "Feature": feature, "Title": title, "Kruskal_H": kw.statistic, "df": 4, "P_value": kw.pvalue})
        q = dunn_holm(cons[feature].to_numpy(float), groups)
        q.insert(0, "Title", title)
        q.insert(0, "Feature", feature)
        q.insert(0, "Order", order)
        pair_rows.append(q)
    return pd.DataFrame(desc_rows), pd.DataFrame(omni_rows), pd.concat(pair_rows, ignore_index=True)


def density_halfwidth(values: np.ndarray, maximum: float = 0.32) -> np.ndarray:
    values = np.asarray(values, float)
    n = max(len(values), 2)
    sd = values.std(ddof=1) if len(values) > 1 else 0
    span = np.ptp(values)
    bw = 1.06 * sd * n ** (-.2)
    if not np.isfinite(bw) or bw <= 1e-12:
        bw = max(span * .08, 1e-3)
    scaled = (values[:, None] - values[None, :]) / bw
    den = np.exp(-.5 * scaled**2).mean(axis=1) / bw
    return maximum * den / max(den.max(), 1e-12)


def plot_distribution(ax, cons: pd.DataFrame, feature: str, title: str, unit: str, pairs: pd.DataFrame, seed: int) -> None:
    rng = np.random.default_rng(seed)
    for i, e in enumerate(E_LEVELS, 1):
        x = cons.loc[cons["E_type"].astype(str).eq(e), feature].to_numpy(float)
        jitter = rng.uniform(-.9, .9, len(x)) * density_halfwidth(x)
        ax.scatter(i+jitter, x, s=2.0, c=E_COLORS[e], alpha=.50, linewidths=0, zorder=2)
        np.random.seed(seed+i)
        ax.boxplot([x], positions=[i], widths=.46, notch=True, bootstrap=3000, whis=(2.5,97.5), showfliers=False, patch_artist=True, manage_ticks=False,
                   boxprops={"facecolor":"none","edgecolor":"#111111","linewidth":.25}, medianprops={"color":"#111111","linewidth":.25},
                   whiskerprops={"color":"#111111","linewidth":.25}, capprops={"color":"#111111","linewidth":.25})
    vals = cons[feature].to_numpy(float)
    lo, hi = vals.min(), vals.max()
    span = hi-lo if hi>lo else 1
    ax.set_ylim(lo-.06*span, hi+.24*span)
    ax.set_xlim(.5,5.5)
    ax.set_xticks(range(1,6), E_LEVELS)
    ax.set_ylabel(f"{title} ({unit})" if unit not in {"index","ratio","CV","Fano"} else title)
    ax.spines[["top","right"]].set_visible(False)
    ax.tick_params(length=1.5,width=.5,pad=1)
    visible = pairs[pairs.P_adj_Holm_within_feature < .01].copy()
    # Greedy compact interval rows; nonoverlapping comparisons share a row.
    occupied=[]; assigned=[]
    for row in visible.sort_values(["P_adj_Holm_within_feature"]).itertuples():
        x1=E_LEVELS.index(row.Group_1)+1; x2=E_LEVELS.index(row.Group_2)+1
        lev=0
        while lev < len(occupied) and any(not(x2+.15<a or b+.15<x1) for a,b in occupied[lev]): lev+=1
        if lev==len(occupied): occupied.append([])
        occupied[lev].append((x1,x2)); assigned.append((row,x1,x2,lev))
    if assigned:
        base=hi+.04*span; step=.18*span/max(len(occupied),1)
        for row,x1,x2,lev in assigned:
            y=base+lev*step
            ax.plot([x1,x2],[y,y],color="#111111",lw=.3,clip_on=False)
            ax.text((x1+x2)/2,y+.012*span,f"{row.P_adj_Holm_within_feature:.1e}".replace("e-","e−"),ha="center",va="bottom",fontsize=3.4)


def plot_feature_tsne(ax, cons: pd.DataFrame, feature: str) -> None:
    z = zscore(cons[feature].to_numpy(float)); order=np.argsort(np.abs(z))
    ax.scatter(cons["tSNE_1"].to_numpy()[order],cons["tSNE_2"].to_numpy()[order],c=z[order],cmap=ZCMAP,norm=Normalize(-ZLIM,ZLIM,clip=True),s=2.7,alpha=.9,linewidths=0)
    for e in E_LEVELS:
        p=cons.loc[cons.E_type.astype(str).eq(e),["tSNE_1","tSNE_2"]].to_numpy(float)
        # 95% covariance ellipse, used only as a stable visual boundary.
        mu=p.mean(axis=0); cov=np.cov(p.T); val,vec=np.linalg.eigh(cov); idx=np.argsort(val)[::-1]; val,vec=val[idx],vec[:,idx]
        theta=np.linspace(0,2*np.pi,240); ring=np.vstack([np.cos(theta),np.sin(theta)])
        boundary=mu[:,None]+vec@np.diag(np.sqrt(np.maximum(val,0)*5.991))@ring
        ax.plot(boundary[0],boundary[1],color=E_COLORS[e],lw=.55,ls=(0,(2.2,1.5)))
    ax.set_axis_off(); ax.set_aspect("equal",adjustable="datalim")


def plot_all25(data: pd.DataFrame) -> tuple[pd.DataFrame,pd.DataFrame,pd.DataFrame]:
    out=OUT/"03_all25_raw_stats_feature_tSNE"; (out/"individual").mkdir(parents=True,exist_ok=True)
    cons=data[data.HC_GC_consensus].copy(); cons["E_type"]=pd.Categorical(cons.HC_GC_consensus_E,E_LEVELS,ordered=True)
    desc,omni,pair=feature_statistics(cons)
    cons.to_csv(out/"HC_GC_consensus450_raw25_and_tSNE.csv",index=False)
    desc.to_csv(out/"HC_GC_consensus450_descriptive_statistics.csv",index=False)
    omni.to_csv(out/"HC_GC_consensus450_KruskalWallis.csv",index=False)
    pair.to_csv(out/"HC_GC_consensus450_Dunn_Holm_pairwise.csv",index=False)
    overview=plt.figure(figsize=(11.0,7.4),facecolor="white")
    gs=overview.add_gridspec(5,10,left=.035,right=.995,bottom=.045,top=.985,wspace=.40,hspace=.48)
    for k,(feature,title,unit) in enumerate(FEATURES):
        q=pair[pair.Feature.eq(feature)]
        fig=plt.figure(figsize=(2.4,1.55),facecolor="white")
        g=fig.add_gridspec(1,2,width_ratios=[.88,1.12],left=.12,right=.99,bottom=.14,top=.96,wspace=.18)
        a=fig.add_subplot(g[0,0]); b=fig.add_subplot(g[0,1]); plot_distribution(a,cons,feature,title,unit,q,777+k); plot_feature_tsne(b,cons,feature)
        save_all(fig,out/"individual"/f"{k+1:02d}_{safe(title)}_HC_GC_consensus450",dpi=600,tif=False); plt.close(fig)
        a=overview.add_subplot(gs[k//5,2*(k%5)]); b=overview.add_subplot(gs[k//5,2*(k%5)+1])
        plot_distribution(a,cons,feature,title,unit,q,777+k); plot_feature_tsne(b,cons,feature)
    cax=overview.add_axes([.43,.012,.14,.009]); cb=overview.colorbar(mpl.cm.ScalarMappable(norm=Normalize(-ZLIM,ZLIM),cmap=ZCMAP),cax=cax,orientation="horizontal"); cb.set_ticks([-2.5,0,2.5]); cb.ax.tick_params(labelsize=4,length=1); cb.set_label("Feature Z-score",fontsize=4,labelpad=1)
    save_all(overview,out/"HC_GC_consensus450_all25_overview",dpi=600); plt.close(overview)
    return desc,omni,pair


def within_group_order(matrix: np.ndarray, labels: np.ndarray) -> np.ndarray:
    order=[]
    for e in E_LEVELS:
        idx=np.where(labels==e)[0]
        if len(idx)>2:
            leaf=leaves_list(linkage(matrix[idx],method="ward",metric="euclidean")); idx=idx[leaf]
        order.extend(idx.tolist())
    return np.asarray(order,int)


def plot_heatmap(data: pd.DataFrame, label_col: str, stem_name: str, consensus: bool=False) -> None:
    out=OUT/"04_heatmaps"
    out.mkdir(parents=True, exist_ok=True)
    use=data[data.HC_GC_consensus].copy() if consensus else data.copy()
    labels=use[label_col].astype(str).to_numpy(); feats=[x[0] for x in FEATURES]
    X=use[feats].to_numpy(float); Z=(X-X.mean(axis=0))/X.std(axis=0,ddof=1); Z=np.clip(Z,-2,2)
    order=within_group_order(Z,labels); Z=Z[order]; labels=labels[order]; ids=use.MSN_unique_ID.to_numpy()[order]
    pd.DataFrame(Z,index=ids,columns=feats).to_csv(out/f"{stem_name}_zscore_matrix.csv")
    fig=plt.figure(figsize=(7.0,3.0),facecolor="white"); ax=fig.add_axes([.18,.16,.79,.74]); im=ax.imshow(Z.T,aspect="auto",interpolation="nearest",cmap="magma",vmin=-2,vmax=2)
    starts=[]; pos=0
    for e in E_LEVELS:
        n=int(np.sum(labels==e)); starts.append((pos,n,e)); pos+=n
    for pos,n,e in starts:
        ax.add_patch(Rectangle((pos-.5,-1.55),n,.32,facecolor=E_COLORS[e],edgecolor="none",clip_on=False)); ax.text(pos+(n-1)/2,-2.0,e,ha="center",va="bottom",fontsize=6)
        if pos>0: ax.axvline(pos-.5,color="white",lw=.8)
    ax.set_yticks(range(len(FEATURES)),[x[1] for x in FEATURES],fontsize=5); ax.set_xticks([]); ax.tick_params(length=0)
    cb=fig.colorbar(im,ax=ax,fraction=.018,pad=.012); cb.set_label("Feature Z-score",fontsize=5); cb.ax.tick_params(labelsize=4,length=1)
    save_all(fig,out/stem_name,dpi=600); plt.close(fig)


def plot_d1d2(labels: pd.DataFrame) -> None:
    if not D1D2_FILE.exists(): return
    out=OUT/"05_D1_D2"; old=pd.read_csv(D1D2_FILE)[["MSN_unique_ID","RPCA_major_class","Final_Subtype"]]
    out.mkdir(parents=True, exist_ok=True)
    d=labels.merge(old,on="MSN_unique_ID",how="left",validate="one_to_one"); d=d[d.HC_GC_consensus].copy(); d.to_csv(out/"HC_GC_consensus450_D1_D2_assignments.csv",index=False)
    fig,axes=plt.subplots(1,2,figsize=(3.2,1.6),constrained_layout=True)
    for ax,major in zip(axes,["D1","D2"]):
        bg=~d.RPCA_major_class.eq(major); ax.scatter(d.loc[bg,"tSNE_1"],d.loc[bg,"tSNE_2"],s=3,c=GREY,alpha=.55,linewidths=0)
        for e in E_LEVELS:
            q=d.RPCA_major_class.eq(major)&d.HC_GC_consensus_E.astype(str).eq(e); ax.scatter(d.loc[q,"tSNE_1"],d.loc[q,"tSNE_2"],s=3,c=E_COLORS[e],alpha=.95,linewidths=0)
        ax.set_title(f"{major} (n={d.RPCA_major_class.eq(major).sum()})",fontweight="bold"); ax.set_xlabel("t-SNE 1"); ax.set_ylabel("t-SNE 2"); ax.set_aspect("equal",adjustable="datalim"); ax.spines[["top","right"]].set_visible(False)
    save_all(fig,out/"HC_GC_consensus450_D1_D2_tSNE_W3p2_H1p6"); plt.close(fig)
    tab=pd.crosstab(d.HC_GC_consensus_E,d.RPCA_major_class).reindex(index=E_LEVELS,columns=["D1","D2"],fill_value=0); pct=tab.div(tab.sum(axis=1),axis=0)*100; tab.to_csv(out/"HC_GC_consensus450_D1_D2_counts.csv"); pct.to_csv(out/"HC_GC_consensus450_D1_D2_percent.csv")
    fig,ax=plt.subplots(figsize=(1.7,1.7)); bottom=np.zeros(5); colors={"D1":"#E31A1C","D2":"#377EB8"}
    for m in ["D1","D2"]: ax.bar(E_LEVELS,pct[m],bottom=bottom,color=colors[m],width=.72,label=m); bottom+=pct[m].to_numpy()
    for i,e in enumerate(E_LEVELS):
        y=0
        for m in ["D1","D2"]:
            v=pct.loc[e,m]; ax.text(i,y+v/2,f"{tab.loc[e,m]}\n{v:.0f}%",ha="center",va="center",fontsize=4,color="white" if v>20 else "#222"); y+=v
    ax.set_ylim(0,100); ax.set_ylabel("Cells within E-type (%)"); ax.legend(frameon=False,fontsize=5,handlelength=1); ax.spines[["top","right"]].set_visible(False)
    save_all(fig,out/"HC_GC_consensus450_D1_D2_100pct_stacked_W1p7_H1p7"); plt.close(fig)

    subtype_levels=sorted(d["Final_Subtype"].dropna().astype(str).unique(), key=lambda x:(x.split("_")[0],int(x.split("_")[1])))
    stab=pd.crosstab(d.HC_GC_consensus_E,d.Final_Subtype).reindex(index=E_LEVELS,columns=subtype_levels,fill_value=0)
    spct=stab.div(stab.sum(axis=1),axis=0)*100
    stab.to_csv(out/"HC_GC_consensus450_subtype_counts.csv"); spct.to_csv(out/"HC_GC_consensus450_subtype_percent.csv")
    subtype_colors=dict(zip(subtype_levels, mpl.colormaps["tab20"](np.linspace(0,1,len(subtype_levels)))))
    fig,ax=plt.subplots(figsize=(3.4,2.0)); bottom=np.zeros(5)
    for s in subtype_levels:
        ax.bar(E_LEVELS,spct[s],bottom=bottom,color=subtype_colors[s],width=.74,label=s,linewidth=0); bottom+=spct[s].to_numpy()
    ax.set_ylim(0,100); ax.set_ylabel("Cells within E-type (%)"); ax.spines[["top","right"]].set_visible(False); ax.legend(frameon=False,fontsize=3.6,ncol=2,bbox_to_anchor=(1.02,1),loc="upper left",handlelength=.8,columnspacing=.7)
    save_all(fig,out/"HC_GC_consensus450_16subtypes_100pct_stacked_W3p4_H2p0"); plt.close(fig)

    fig,axes=plt.subplots(4,4,figsize=(4.8,4.8),constrained_layout=True)
    for ax,s in zip(axes.flat,subtype_levels):
        q=d.Final_Subtype.astype(str).eq(s); ax.scatter(d.loc[~q,"tSNE_1"],d.loc[~q,"tSNE_2"],s=1.5,c=GREY,alpha=.35,linewidths=0)
        for e in E_LEVELS:
            z=q & d.HC_GC_consensus_E.astype(str).eq(e); ax.scatter(d.loc[z,"tSNE_1"],d.loc[z,"tSNE_2"],s=2.2,c=E_COLORS[e],alpha=.95,linewidths=0)
        ax.set_title(f"{s} (n={q.sum()})",fontsize=5); ax.set_axis_off(); ax.set_aspect("equal",adjustable="datalim")
    save_all(fig,out/"HC_GC_consensus450_16subtypes_tSNE_grid_W4p8_H4p8"); plt.close(fig)


def write_deprecation_manifest() -> None:
    rows=[
        ("Expert decision tree / phenotype strategy", "Deprecated", "Definition requires expert labels."),
        ("Expert-vs-consensus cross matrices", "Deprecated", "Expert classification removed."),
        ("Expert-labelled circular dendrogram outer ring", "Deprecated", "Replace with HC tree plus GC ring."),
        ("Expert phenotype names (LIR-HR, IE-RF, HIR-LA, HE-HA, SA-HT)", "Deprecated", "Use aligned E1-E5 only unless independently renamed later."),
        ("Representative raw traces selected by expert confidence", "Not regenerated", "Representative cells must be reselected as HC×GC medoids before plotting."),
        ("Expert-group radar plots", "Deprecated", "No expert grouping retained."),
    ]
    pd.DataFrame(rows,columns=["Figure_family","Status","Reason"]).to_csv(OUT/"Expert_dependent_figures_deprecation_manifest.csv",index=False)


def main() -> None:
    style(); OUT.mkdir(parents=True,exist_ok=True)
    labels,data,pca_data=load(); metrics=classification_metrics(labels)
    labels.to_csv(OUT/"HC_GC_only_cell_assignments.csv",index=False)
    (OUT/"HC_GC_only_run_summary.json").write_text(json.dumps(metrics,indent=2),encoding="utf-8")
    plot_classification(labels,metrics); plot_pca(pca_data); plot_all25(data)
    plot_heatmap(data,"HC_E","HC_K5_all493_raw25_zscore_heatmap")
    plot_heatmap(data,"GC_E","GC_K5_all493_raw25_zscore_heatmap")
    plot_heatmap(data,"HC_GC_consensus_E","HC_GC_consensus450_raw25_zscore_heatmap",consensus=True)
    plot_d1d2(labels); write_deprecation_manifest()
    manifest=sorted(str(p.relative_to(OUT)).replace("\\","/") for p in OUT.rglob("*") if p.is_file())
    (OUT/"output_manifest.txt").write_text("\n".join(manifest)+"\n",encoding="utf-8")
    print(json.dumps(metrics,indent=2)); print(f"Generated {len(manifest)} files in {OUT}")


if __name__ == "__main__":
    main()
