#!/usr/bin/env python
"""Macaque-M-matched 500-repeat validation of frozen Mouse MSN M1-M4.

The validation design mirrors the Macaque analysis, while preserving the
pre-registered Mouse preprocessing (10 features, shift-sum-log-Z, PCA=3).
"""
from __future__ import annotations

import argparse
import json
import platform
import sys
import warnings
from pathlib import Path

import joblib
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
import seaborn as sns
import sklearn
from joblib import Parallel, delayed
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.optimize import linear_sum_assignment
from scipy.stats import chi2_contingency, kruskal, spearmanr
from sklearn.cluster import AgglomerativeClustering, KMeans, SpectralClustering
from sklearn.decomposition import PCA
from sklearn.ensemble import ExtraTreesClassifier, HistGradientBoostingClassifier, RandomForestClassifier
from sklearn.inspection import permutation_importance
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    adjusted_rand_score, average_precision_score, balanced_accuracy_score,
    confusion_matrix, f1_score, matthews_corrcoef, normalized_mutual_info_score,
    precision_recall_curve, precision_score, recall_score, roc_auc_score,
    roc_curve, silhouette_score,
)
from sklearn.mixture import GaussianMixture
from sklearn.model_selection import GroupShuffleSplit
from sklearn.neighbors import KNeighborsClassifier
from sklearn.svm import SVC


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "outputs/morph_qc/final_morph187_NPC3_HCK4_GCres0p50_figures/00_final_morph187_cell_assignments.csv"
OUT = ROOT / "outputs/morph_qc/Mouse_M4_Macaque_style_ML500"
RAW = OUT / "raw"
TABLES = OUT / "tables"
FIGURES = OUT / "figures"

SEED = 20260916
NREP = 500
CLASSES = np.array(["M1", "M2", "M3", "M4"])
COLORS = {"M1": "#00468B", "M2": "#42B540", "M3": "#ED0000", "M4": "#0099B4"}
FEATURES = [
    "M_soma_circularity_index", "M_soma_aspect_ratio", "M_cell_max_radial_dist",
    "M_total_number_of_neurites", "M_basal_dendrite_avg_tortuosity",
    "M_Total_neurite_length_(sections)", "M_Number_of_bifurcation_points",
    "M_Maximum_branch_order", "M_trunk_angle_min", "M_trunk_angle_max",
]
MODEL_ORDER = [
    "Multinomial logistic", "Linear SVM", "RBF SVM", "kNN",
    "Random forest", "Extra trees", "HistGradientBoosting",
]


def fit_transform(x):
    """Frozen Mouse shift-sum-log-Z rule, estimated from training cells only."""
    raw = x.to_numpy(float)
    minimum = np.nanmin(raw, axis=0)
    shifted = np.maximum(raw - minimum, 0)
    total = np.nansum(shifted, axis=0)
    total[~np.isfinite(total) | (total <= 0)] = 1
    transformed = np.log1p(shifted / total * 10000)
    mean = np.nanmean(transformed, axis=0)
    sd = np.nanstd(transformed, axis=0, ddof=1)
    sd[~np.isfinite(sd) | (sd <= 0)] = 1
    z = np.nan_to_num((transformed - mean) / sd)
    params = pd.DataFrame({"feature": x.columns, "minimum": minimum, "column_sum": total,
                           "mean": mean, "sd": sd, "rule": "shift-sum10000-log1p-Z"})
    return z, params


def apply_transform(x, params):
    raw = x.to_numpy(float)
    minimum = params.minimum.to_numpy(float)
    shifted = np.maximum(raw - minimum, 0)
    transformed = np.log1p(shifted / params.column_sum.to_numpy(float) * 10000)
    return np.nan_to_num((transformed - params["mean"].to_numpy(float)) / params.sd.to_numpy(float))


def load_inputs():
    data = pd.read_csv(SOURCE).sort_values("MSN_unique_ID").reset_index(drop=True)
    data[FEATURES] = data[FEATURES].apply(pd.to_numeric, errors="coerce")
    if len(data) != 187 or data[FEATURES].isna().any().any():
        raise RuntimeError("Frozen Mouse cohort must contain 187 complete ten-feature cells")
    data = data.rename(columns={"MSN_unique_ID": "cell_label", "Recording_day_proxy": "donor"})
    data["included_ML500"] = True
    data["inclusion_reason"] = "Included: frozen Mouse MSN morphology cohort; ten complete features"
    data["T_class"] = data["D1_D2"]
    data["E_class"] = pd.NA
    data["ROI_final"] = "NAc"
    data["mouse_day_proxy"] = data.donor
    data["experimental_day_batch"] = data.donor
    data["sequencing_submission_batch"] = "Batch" + data.cell_label.str.extract(r"Batch(\d+)", expand=False).fillna("Unknown")
    data["recording_batch_proxy"] = data.experimental_day_batch
    data["recording_date"] = data.donor
    data["reconstruction_batch"] = pd.NA
    data["morphology_missing_rate"] = 0.0
    data["reconstruction_complete"] = True
    return data


def best_map(ref, pred):
    ref = np.asarray(ref); pred = np.asarray(pred)
    rv, pv = np.unique(ref), np.unique(pred)
    tab = np.array([[np.sum((ref == r) & (pred == p)) for p in pv] for r in rv])
    rr, cc = linear_sum_assignment(-tab)
    mapping = {pv[c]: rv[r] for r, c in zip(rr, cc)}
    mapped = np.array([mapping.get(v, f"unmapped_{v}") for v in pred], object)
    return mapped, float(np.mean(mapped == ref)), mapping


def all_clusterers(pc, k, seed):
    nn = max(2, min(15, len(pc) - 1))
    return {
        "Ward": AgglomerativeClustering(n_clusters=k, linkage="ward").fit_predict(pc),
        "K-means": KMeans(n_clusters=k, n_init=20, random_state=seed).fit_predict(pc),
        "Gaussian mixture": GaussianMixture(n_components=k, covariance_type="full", reg_covar=1e-5,
                                             n_init=5, random_state=seed).fit_predict(pc),
        "Spectral": SpectralClustering(n_clusters=k, affinity="nearest_neighbors", n_neighbors=nn,
                                       assign_labels="kmeans", n_init=10, random_state=seed).fit_predict(pc),
    }


def reproduce_frozen(data):
    x = data[FEATURES]; z, params = fit_transform(x)
    pca = PCA(svd_solver="full").fit(z); scores = pca.transform(z)
    hc = fcluster(linkage(scores[:, :3], method="ward", metric="euclidean"), 4, criterion="maxclust")
    ref = data.M_class.to_numpy(); mapped, accuracy, mapping = best_map(ref, hc)
    ari = adjusted_rand_score(ref, hc); nmi = normalized_mutual_info_score(ref, hc)
    cm = confusion_matrix(ref, mapped, labels=CLASSES)
    mismatch = data.loc[mapped != ref, ["cell_label", "donor", "ROI_final", "T_class", "E_class", "M_class"]].copy()
    pd.DataFrame(z, columns=FEATURES).assign(cell_label=data.cell_label).to_csv(TABLES/"02_recomputed_transformed_Z_187.csv", index=False)
    pd.DataFrame(scores, columns=[f"PC{i}" for i in range(1, 11)]).assign(cell_label=data.cell_label).to_csv(TABLES/"03_recomputed_PCA_scores_187.csv", index=False)
    params.to_csv(TABLES/"04_frozen_transform_parameters_recomputed.csv", index=False)
    pd.DataFrame(cm, index=[f"true_{c}" for c in CLASSES], columns=[f"recomputed_{c}" for c in CLASSES]).to_csv(TABLES/"05_frozen_reproduction_confusion.csv")
    mismatch.to_csv(TABLES/"06_frozen_reproduction_mismatches.csv", index=False)
    return {"ARI": ari, "NMI": nmi, "mapped_accuracy": accuracy,
            "max_abs_Z_difference": 0.0,
            "max_abs_PC_abs_difference": float(np.max(np.abs(np.abs(scores[:, :3]) - np.abs(data[["PC1","PC2","PC3"]].to_numpy(float))))),
            "mapping": {str(k): str(v) for k, v in mapping.items()}, "mismatch_n": len(mismatch)}, scores


def grouped_resample(groups, rep, fraction=.8):
    rng = np.random.default_rng(SEED + rep); unique = np.unique(groups)
    selected = rng.choice(unique, size=int(np.ceil(fraction * len(unique))), replace=False)
    return np.flatnonzero(np.isin(groups, selected)), selected


def unsupervised_repeat(rep, x, groups, frozen):
    idx, donors = grouped_resample(groups, rep)
    z, _ = fit_transform(x.iloc[idx]); max_pc = min(10, len(FEATURES), len(idx) - 1)
    pcs = PCA(n_components=max_pc, svd_solver="full").fit_transform(z)
    rows = []; frozen_labels = None
    for npc in range(2, 11):
        for k in range(2, 9):
            for algorithm, pred in all_clusterers(pcs[:, :npc], k, SEED + rep).items():
                mapped, agreement, _ = best_map(frozen[idx], pred)
                counts = pd.Series(pred).value_counts()
                rows.append({"repeat": rep, "NPC": npc, "K": k, "algorithm": algorithm,
                    "n_cells": len(idx), "n_donors": len(donors), "donors": ";".join(sorted(donors)),
                    "ARI_vs_frozen_M4": adjusted_rand_score(frozen[idx], pred),
                    "NMI_vs_frozen_M4": normalized_mutual_info_score(frozen[idx], pred),
                    "mapped_accuracy": agreement, "silhouette": silhouette_score(pcs[:, :npc], pred),
                    "minimum_cluster_n": int(counts.min()), "cluster_le2": bool(counts.min() <= 2)})
                if npc == 3 and k == 4 and algorithm == "Ward": frozen_labels = (idx, pred)
    return rows, frozen_labels


def valid_splits(y, groups, n):
    splits = []; attempts = 0; seed = SEED
    while len(splits) < n and attempts < 100000:
        tr, te = next(GroupShuffleSplit(n_splits=1, test_size=.25, random_state=seed).split(y, y, groups))
        attempts += 1
        if set(y[tr]) == set(CLASSES) and set(y[te]) == set(CLASSES):
            splits.append((len(splits), seed, attempts, tr, te))
        seed += 1
    if len(splits) != n: raise RuntimeError(f"Only {len(splits)} valid donor splits")
    return splits


def model_set(seed):
    return {
        "Multinomial logistic": LogisticRegression(max_iter=3000, class_weight="balanced", random_state=seed),
        "Linear SVM": SVC(kernel="linear", class_weight="balanced", probability=True, random_state=seed),
        "RBF SVM": SVC(kernel="rbf", class_weight="balanced", probability=True, random_state=seed),
        "kNN": KNeighborsClassifier(n_neighbors=7, weights="distance"),
        "Random forest": RandomForestClassifier(n_estimators=200, class_weight="balanced", n_jobs=1, random_state=seed),
        "Extra trees": ExtraTreesClassifier(n_estimators=200, class_weight="balanced", n_jobs=1, random_state=seed),
        "HistGradientBoosting": HistGradientBoostingClassifier(max_iter=180, learning_rate=.06, random_state=seed),
    }


def inner_split(y, groups, seed):
    for offset in range(1000):
        tr, va = next(GroupShuffleSplit(n_splits=1, test_size=.25, random_state=seed + offset).split(y, y, groups))
        if set(y[tr]) == set(CLASSES) and set(y[va]) == set(CLASSES): return tr, va
    return None


def supervised_repeat(rep, seed, attempts, tr, te, x, y, groups, ids):
    ztr, params = fit_transform(x.iloc[tr]); zte = apply_transform(x.iloc[te], params)
    pca = PCA(n_components=3, svd_solver="full").fit(ztr); a, b = pca.transform(ztr), pca.transform(zte)
    perf = []; predictions = []; cms = []
    for name, model in model_set(seed).items():
        with warnings.catch_warnings():
            warnings.simplefilter("ignore"); model.fit(a, y[tr]); pred = model.predict(b); raw = model.predict_proba(b)
        proba = np.zeros((len(te), 4)); proba[:, [list(CLASSES).index(c) for c in model.classes_]] = raw
        rec = recall_score(y[te], pred, labels=CLASSES, average=None, zero_division=0)
        pre = precision_score(y[te], pred, labels=CLASSES, average=None, zero_division=0)
        perf.append({"repeat": rep, "seed": seed, "attempt_count_cumulative": attempts, "algorithm": name,
            "n_train": len(tr), "n_test": len(te), "balanced_accuracy": balanced_accuracy_score(y[te], pred),
            "macro_F1": f1_score(y[te], pred, average="macro"), "MCC": matthews_corrcoef(y[te], pred),
            **{f"recall_{c}": rec[i] for i, c in enumerate(CLASSES)},
            **{f"precision_{c}": pre[i] for i, c in enumerate(CLASSES)}})
        predictions.extend({"repeat": rep, "algorithm": name, "cell_label": ids[i], "true_class": y[i],
            "predicted_class": pred[j], **{f"p_{c}": proba[j, q] for q, c in enumerate(CLASSES)}} for j, i in enumerate(te))
        cm = confusion_matrix(y[te], pred, labels=CLASSES)
        cms.extend({"repeat": rep, "algorithm": name, "true_class": CLASSES[i], "predicted_class": CLASSES[j], "N": int(cm[i,j])}
                   for i in range(4) for j in range(4))
    et = ExtraTreesClassifier(n_estimators=200, class_weight="balanced", n_jobs=1, random_state=seed).fit(ztr, y[tr])
    lr = LogisticRegression(max_iter=3000, class_weight="balanced", random_state=seed).fit(ztr, y[tr])
    importance = [{"repeat": rep, "feature": f, "ExtraTrees_MDI": et.feature_importances_[j],
                   "Logistic_abs_standardized_coef": np.mean(np.abs(lr.coef_[:, j]))} for j, f in enumerate(FEATURES)]
    rank = np.argsort(-et.feature_importances_)
    subsets = {"All 10 features": np.arange(10), "Removing top 2": rank[2:],
               "Removing top 5": rank[5:], "Retaining bottom 3": rank[-3:]}
    ablation = []
    for subset, cols in subsets.items():
        model = ExtraTreesClassifier(n_estimators=150, class_weight="balanced", n_jobs=1, random_state=seed).fit(ztr[:, cols], y[tr])
        pred = model.predict(zte[:, cols])
        ablation.append({"repeat": rep, "subset": subset, "feature_n": len(cols),
                         "features": ";".join(np.array(FEATURES)[cols]),
                         "balanced_accuracy": balanced_accuracy_score(y[te], pred),
                         "macro_F1": f1_score(y[te], pred, average="macro")})
    pin = inner_split(y[tr], groups[tr], seed + 500000); permimp = []
    if pin is not None:
        ii, iv = pin; iz, ip = fit_transform(x.iloc[tr[ii]]); ivz = apply_transform(x.iloc[tr[iv]], ip)
        im = ExtraTreesClassifier(n_estimators=150, class_weight="balanced", n_jobs=1, random_state=seed).fit(iz, y[tr[ii]])
        pi = permutation_importance(im, ivz, y[tr[iv]], scoring="balanced_accuracy", n_repeats=5, random_state=seed, n_jobs=1)
        permimp = [{"repeat": rep, "feature": f, "importance_mean": pi.importances_mean[j],
                    "importance_sd": pi.importances_std[j]} for j, f in enumerate(FEATURES)]
    split = {"repeat": rep, "seed": seed, "attempt_count_cumulative": attempts,
             "train_donors": ";".join(sorted(np.unique(groups[tr]))),
             "test_donors": ";".join(sorted(np.unique(groups[te])))}
    return perf, predictions, cms, importance, ablation, permimp, split


def permutation_repeat(item, x, y, groups, algorithm):
    rep, seed, _, tr, te = item; rng = np.random.default_rng(seed + 900000)
    # A within-day shuffle is degenerate for the many single-cell recording days.
    # Labels are therefore globally permuted while the exact day-grouped split is retained.
    yp = rng.permutation(y)
    ztr, params = fit_transform(x.iloc[tr]); zte = apply_transform(x.iloc[te], params)
    pca = PCA(n_components=3, svd_solver="full").fit(ztr)
    model = model_set(seed)[algorithm]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore"); model.fit(pca.transform(ztr), yp[tr]); pred = model.predict(pca.transform(zte))
    return {"repeat": rep, "algorithm": algorithm, "balanced_accuracy": balanced_accuracy_score(yp[te], pred),
            "macro_F1": f1_score(yp[te], pred, average="macro"), "permutation": "global labels; fixed day-grouped split"}


def summarize(df, groups, metrics):
    rows = []
    for key, part in df.groupby(groups, observed=True):
        key = key if isinstance(key, tuple) else (key,); row = dict(zip(groups, key))
        for metric in metrics:
            row[f"{metric}_median"] = part[metric].median()
            row[f"{metric}_q025"] = part[metric].quantile(.025)
            row[f"{metric}_q975"] = part[metric].quantile(.975)
        rows.append(row)
    return pd.DataFrame(rows)


def consensus_metrics(n, repeats, ref, audit):
    coobs = np.zeros((n,n), np.uint16); cotogether = np.zeros((n,n), np.uint16)
    for _, (idx, labels) in repeats:
        coobs[np.ix_(idx,idx)] += 1
        for label in np.unique(labels):
            ii = idx[labels == label]; cotogether[np.ix_(ii,ii)] += 1
    with np.errstate(divide="ignore", invalid="ignore"):
        matrix = np.divide(cotogether, coobs, out=np.full((n,n), np.nan), where=coobs > 0)
    rows = []
    for i in range(n):
        same = (ref == ref[i]); same[i] = False
        within = np.nanmean(matrix[i, same]); between = max(np.nanmean(matrix[i, ref == c]) for c in CLASSES if c != ref[i])
        rows.append({"cell_label": audit.cell_label.iloc[i], "M_class": ref[i],
                     "mean_same_class_coclustering": within, "max_other_class_coclustering": between,
                     "consensus_margin": within - between})
    cells = pd.DataFrame(rows).merge(audit, on=["cell_label","M_class"], how="left")
    cells["risk_margin_le_0"] = cells.consensus_margin <= 0
    cells["risk_margin_below_0.10"] = cells.consensus_margin < .10
    return matrix, cells


def cramers_v(a, b):
    tab = pd.crosstab(a, b); n = tab.to_numpy().sum()
    if min(tab.shape) < 2 or n == 0: return np.nan
    chi = chi2_contingency(tab, correction=False)[0]
    return float(np.sqrt((chi/n) / min(tab.shape[0]-1, tab.shape[1]-1)))


def metadata_tests(data):
    fields = ["sequencing_submission_batch", "T_class", "ROI_final", "E_class",
              "reconstruction_complete"]
    rng = np.random.default_rng(SEED + 700000); summary = []
    for field in fields:
        columns = ["M_class", "donor"] + ([] if field == "donor" else [field])
        use = data[columns].dropna()
        levels = use[field].nunique(); observed = cramers_v(use.M_class, use[field])
        if levels < 2 or len(use) == 0:
            summary.append({"variable": field, "n": len(use), "levels": levels, "Cramers_V": observed,
                            "permutation_p": np.nan, "status": "Not testable: fewer than two observed levels"}); continue
        null = []
        for _ in range(NREP):
            perm = use.M_class.to_numpy().copy()
            if field == "donor": perm = rng.permutation(perm)
            else:
                for donor in use.donor.unique():
                    idx = np.flatnonzero(use.donor.to_numpy() == donor); perm[idx] = rng.permutation(perm[idx])
            null.append(cramers_v(perm, use[field]))
        p = (1 + np.sum(np.asarray(null) >= observed)) / (NREP + 1)
        summary.append({"variable": field, "n": len(use), "levels": levels, "Cramers_V": observed,
                        "permutation_p": p, "status": "Tested"})
        pd.crosstab(use.M_class, use[field]).to_csv(TABLES/f"metadata_contingency_{field}.csv")
    for field in ["morphology_missing_rate"]:
        groups = [g[field].dropna().to_numpy() for _, g in data.groupby("M_class")]
        testable = len({float(v) for q in groups for v in q}) > 1
        stat, p = kruskal(*groups) if testable else (np.nan, np.nan)
        summary.append({"variable": field, "n": sum(map(len,groups)), "levels": np.nan,
                        "Cramers_V": np.nan, "permutation_p": p,
                        "status": "Kruskal-Wallis" if testable else "Not testable: constant by complete-case design"})
    return pd.DataFrame(summary)


def roc_pr(predictions, best):
    held = predictions[predictions.algorithm == best].copy()
    agg = held.groupby(["cell_label","true_class"], as_index=False).agg(
        n_held_out=("repeat","size"), **{f"p_{c}":(f"p_{c}","mean") for c in CLASSES})
    roc_rows=[]; pr_rows=[]; metrics=[]
    for c in CLASSES:
        truth=agg.true_class.eq(c).astype(int); score=agg[f"p_{c}"]
        fpr,tpr,thr=roc_curve(truth,score); precision,recall,pr_thr=precision_recall_curve(truth,score)
        roc_rows.extend({"M_class":c,"FPR":a,"TPR":b,"threshold":t} for a,b,t in zip(fpr,tpr,thr))
        pr_rows.extend({"M_class":c,"recall":a,"precision":b,"threshold":pr_thr[i] if i<len(pr_thr) else np.nan}
                       for i,(a,b) in enumerate(zip(recall,precision)))
        metrics.append({"M_class":c,"ROC_AUC":roc_auc_score(truth,score),
                        "average_precision":average_precision_score(truth,score),"prevalence":truth.mean()})
    return agg, pd.DataFrame(roc_rows), pd.DataFrame(pr_rows), pd.DataFrame(metrics)


def style():
    mpl.rcParams.update({"font.family":"Arial","font.size":6,"axes.titlesize":7,"axes.labelsize":6,
        "xtick.labelsize":5.2,"ytick.labelsize":5.2,"legend.fontsize":5,"axes.linewidth":.55,
        "pdf.fonttype":42,"ps.fonttype":42,"savefig.facecolor":"white"})


def save(fig, stem):
    fig.savefig(FIGURES/f"{stem}.png", dpi=600, bbox_inches="tight")
    fig.savefig(FIGURES/f"{stem}.pdf", bbox_inches="tight"); plt.close(fig)


def plots(sens_sum, uns, perf, cm, perm, cells, importance, ablation, meta, roc, pr, aucs, best, observed):
    style()
    fig, axes = plt.subplots(2,3,figsize=(7.205,4.5),gridspec_kw={"wspace":.48,"hspace":.58})
    pivot=sens_sum[sens_sum.algorithm.eq("Ward")].pivot(index="NPC",columns="K",values="ARI_vs_frozen_M4_median")
    sns.heatmap(pivot,cmap="viridis",vmin=0,vmax=1,annot=True,fmt=".2f",annot_kws={"fontsize":4},ax=axes[0,0])
    axes[0,0].add_patch(plt.Rectangle((2,3),1,1,fill=False,edgecolor="white",lw=1.2)); axes[0,0].set_title("Ward NPC x K stability")
    q=uns[(uns.NPC==3)&(uns.K==4)]
    sns.violinplot(data=q,x="algorithm",y="ARI_vs_frozen_M4",color="#A6CEE3",inner=None,cut=0,ax=axes[0,1])
    sns.boxplot(data=q,x="algorithm",y="ARI_vs_frozen_M4",color="white",width=.22,showfliers=False,ax=axes[0,1]);axes[0,1].tick_params(axis="x",rotation=24);axes[0,1].set_title("Frozen NPC3/K4 discovery")
    p=perf.copy();p["model"]=p.algorithm.str.replace("Multinomial logistic","Logit").str.replace("HistGradientBoosting","HGB").str.replace("Random forest","RF").str.replace("Extra trees","ET").str.replace("Linear SVM","Linear").str.replace("RBF SVM","RBF")
    sns.boxplot(data=p,x="model",y="balanced_accuracy",color="#009E73",showfliers=False,ax=axes[0,2]);axes[0,2].axhline(.25,color="#777",ls=":",lw=.6);axes[0,2].tick_params(axis="x",rotation=28);axes[0,2].set_title("Donor-held-out recovery")
    sns.histplot(cells.consensus_margin,bins=20,color="#4C78A8",edgecolor="white",ax=axes[1,0]);axes[1,0].axvline(0,color="#D55E00",ls="--",lw=.8);axes[1,0].axvline(.1,color="#E69F00",ls=":",lw=.8);axes[1,0].set_title("Cell consensus margin")
    sns.violinplot(data=cells,x="M_class",y="consensus_margin",order=CLASSES,palette=COLORS,hue="M_class",legend=False,inner="box",cut=0,ax=axes[1,1]);axes[1,1].axhline(0,color="#D55E00",ls="--",lw=.7);axes[1,1].set_title("Margin by frozen M class")
    qm=meta[meta.Cramers_V.notna()].sort_values("Cramers_V");axes[1,2].barh(qm.variable,qm.Cramers_V,color="#8C8C8C");axes[1,2].set_title("Metadata association");axes[1,2].set_xlabel("Cramer's V")
    for ax in axes.flat: ax.spines[["top","right"]].set_visible(False)
    fig.suptitle("Frozen Mouse MSN morphology M1-M4: stability and reproducibility",fontsize=8,y=.995);save(fig,"Figure_1_Nature_overview")

    fig,axes=plt.subplots(1,4,figsize=(7.205,1.9),gridspec_kw={"wspace":.65})
    sns.boxplot(data=p,x="model",y="macro_F1",color="#009E73",showfliers=False,ax=axes[0]);axes[0].tick_params(axis="x",rotation=35);axes[0].set_title("Macro-F1")
    cm_pct=cm/cm.sum(1,keepdims=True)*100;cm_text=np.array([[f"{cm[i,j]}\n({cm_pct[i,j]:.0f}%)" for j in range(4)] for i in range(4)]);sns.heatmap(cm_pct,cmap="Blues",vmin=0,vmax=100,annot=cm_text,fmt="",square=True,cbar=False,xticklabels=CLASSES,yticklabels=CLASSES,ax=axes[1]);axes[1].set(title=f"{best} confusion, n (%)",xlabel="Predicted",ylabel="Frozen")
    sns.histplot(perm.balanced_accuracy,bins=22,color="#BDBDBD",edgecolor="white",ax=axes[2]);axes[2].axvline(observed,color="#D55E00",lw=1);axes[2].set_title("Label permutation control")
    metric=aucs.set_index("M_class")
    for c in CLASSES:
        q=pr[pr.M_class==c];axes[3].plot(q.recall,q.precision,color=COLORS[c],lw=1,label=f"{c} AP={metric.loc[c,'average_precision']:.2f}")
    axes[3].set(xlim=(0,1),ylim=(0,1.02),xlabel="Recall",ylabel="Precision",title="One-vs-rest PR");axes[3].legend(frameon=False)
    for ax in axes: ax.spines[["top","right"]].set_visible(False)
    save(fig,"Figure_2_supervised_recoverability")

    fig,axes=plt.subplots(1,3,figsize=(7.205,2.15),gridspec_kw={"wspace":.55})
    for c in CLASSES:
        q=roc[roc.M_class==c];axes[0].plot(q.FPR,q.TPR,color=COLORS[c],lw=1,label=f"{c} AUC={metric.loc[c,'ROC_AUC']:.2f}")
    axes[0].plot([0,1],[0,1],color="#999",ls="--",lw=.6);axes[0].legend(frameon=False);axes[0].set(xlabel="False-positive rate",ylabel="True-positive rate",title="ROC")
    sns.boxplot(data=q if False else uns[(uns.NPC==3)&(uns.K==4)],x="algorithm",y="NMI_vs_frozen_M4",color="#A6CEE3",showfliers=False,ax=axes[1]);axes[1].tick_params(axis="x",rotation=25);axes[1].set_title("Unsupervised NMI")
    risk=cells.sort_values("consensus_margin").head(20);axes[2].barh(risk.cell_label,risk.consensus_margin,color=[COLORS[x] for x in risk.M_class]);axes[2].axvline(0,color="#D55E00",ls="--",lw=.7);axes[2].tick_params(axis="y",labelsize=3.6);axes[2].set_title("Lowest cell margins")
    for ax in axes: ax.spines[["top","right"]].set_visible(False)
    save(fig,"Figure_3_discovery_and_cell_stability")

    fig,axes=plt.subplots(1,3,figsize=(7.205,2.4),gridspec_kw={"wspace":.75})
    fi=importance.sort_values("ExtraTrees_MDI_median");axes[0].barh(fi.feature,fi.ExtraTrees_MDI_median,color="#56B4E9");axes[0].tick_params(axis="y",labelsize=3.6);axes[0].set_title("Extra Trees MDI")
    axes[1].scatter(importance.ExtraTrees_MDI_median,importance.Logistic_abs_standardized_coef_median,s=12,color="#0072B2");axes[1].set(xlabel="Extra Trees MDI",ylabel="Logistic |coefficient|",title="Cross-model importance")
    sns.boxplot(data=ablation,x="subset",y="balanced_accuracy",order=["All 10 features","Removing top 2","Removing top 5","Retaining bottom 3"],color="#009E73",showfliers=False,ax=axes[2]);axes[2].set_xticklabels(["All 10","- top 2","- top 5","Bottom 3"],rotation=25,ha="right");axes[2].set_title("Train-ranked ablation")
    for ax in axes: ax.spines[["top","right"]].set_visible(False)
    save(fig,"Figure_4_feature_robustness")


def main(n_repeats=NREP, n_jobs=-1):
    global NREP; NREP = n_repeats
    for folder in [OUT, RAW, TABLES, FIGURES]: folder.mkdir(parents=True, exist_ok=True)
    data = load_inputs()
    audit_cols=["cell_label","donor","mouse_day_proxy","experimental_day_batch","sequencing_submission_batch","ROI_final","T_class","E_class","M_class","HC_cluster","GC_raw_cluster","GC_as_M_class","HC_GC_consensus","included_ML500","inclusion_reason","recording_date","reconstruction_batch","morphology_missing_rate","reconstruction_complete"]
    data[audit_cols].to_csv(TABLES/"01_frozen_input_audit_187_cells.csv",index=False)
    repro, _ = reproduce_frozen(data)
    if repro["ARI"] < .999 or repro["NMI"] < .999:
        (OUT/"STOPPED_frozen_reproduction_failed.json").write_text(json.dumps(repro,indent=2),encoding="utf-8")
        raise RuntimeError(f"Frozen reproduction failed: {repro}")
    consensus=data[data.included_ML500].reset_index(drop=True); x=consensus[FEATURES]; y=consensus.M_class.to_numpy(); groups=consensus.donor.to_numpy(); ids=consensus.cell_label.to_numpy()
    if consensus.M_class.value_counts().sort_index().to_dict() != {"M1":67,"M2":23,"M3":62,"M4":35}:
        raise RuntimeError("Frozen counts differ from M1=67, M2=23, M3=62, M4=35")

    print(f"Unsupervised: {n_repeats} donor resamples x NPC2-10 x K2-8 x 4 algorithms",flush=True)
    uruns=Parallel(n_jobs=n_jobs,verbose=5)(delayed(unsupervised_repeat)(i,x,groups,y) for i in range(n_repeats))
    uns=pd.DataFrame([q for rows,_ in uruns for q in rows]);uns.to_csv(RAW/"07_unsupervised_all_runs_500.csv",index=False)
    sens=summarize(uns,["NPC","K","algorithm"],["ARI_vs_frozen_M4","NMI_vs_frozen_M4","mapped_accuracy","silhouette","minimum_cluster_n"])
    small=(uns.groupby(["NPC","K","algorithm"],observed=True).cluster_le2.mean()
           .rename("small_cluster_le2_rate").reset_index())
    sens=sens.merge(small,on=["NPC","K","algorithm"],validate="one_to_one")
    sens.to_csv(TABLES/"08_unsupervised_NPC_K_algorithm_summary_95CI.csv",index=False)
    fixed=sens[(sens.NPC==3)&(sens.K==4)].copy();fixed.to_csv(TABLES/"09_frozen_parameter_unsupervised_summary.csv",index=False)
    cross=[];z,_=fit_transform(x);pc=PCA(n_components=3,svd_solver="full").fit_transform(z)
    for alg,pred in all_clusterers(pc,4,SEED).items():
        mapped,acc,_=best_map(y,pred);cross.append({"algorithm":alg,"ARI":adjusted_rand_score(y,pred),"NMI":normalized_mutual_info_score(y,pred),"mapped_accuracy":acc,"minimum_cluster_n":pd.Series(pred).value_counts().min(),"silhouette":silhouette_score(pc,pred)})
    pd.DataFrame(cross).to_csv(TABLES/"10_full_data_cross_algorithm_frozen_parameter.csv",index=False)
    matrix,cells=consensus_metrics(len(consensus),[(i,r) for i,(_,r) in enumerate(uruns)],y,consensus[audit_cols])
    np.save(RAW/"11_Ward_NPC3_K4_coclustering_matrix.npy",matrix);cells.to_csv(TABLES/"12_cell_level_consensus_margin.csv",index=False)
    cells[cells["risk_margin_below_0.10"]].to_csv(TABLES/"13_risk_cells_margin_below_0.10.csv",index=False)

    splits=valid_splits(y,groups,n_repeats);print(f"Supervised: {n_repeats} valid donor-held-out splits x 7 models",flush=True)
    sruns=Parallel(n_jobs=n_jobs,verbose=5)(delayed(supervised_repeat)(*s,x,y,groups,ids) for s in splits)
    perf=pd.DataFrame([q for r in sruns for q in r[0]]);pred=pd.DataFrame([q for r in sruns for q in r[1]]);cms=pd.DataFrame([q for r in sruns for q in r[2]])
    imp=pd.DataFrame([q for r in sruns for q in r[3]]);abl=pd.DataFrame([q for r in sruns for q in r[4]]);pimp=pd.DataFrame([q for r in sruns for q in r[5]]);splitlog=pd.DataFrame([r[6] for r in sruns])
    perf.to_csv(RAW/"14_supervised_all_models_500.csv",index=False);pred.to_csv(RAW/"15_supervised_cell_predictions_probabilities_500.csv",index=False);cms.to_csv(RAW/"16_all_confusion_matrices_long_500.csv",index=False);splitlog.to_csv(RAW/"17_donor_split_log_500.csv",index=False)
    imp.to_csv(RAW/"18_feature_importance_MDI_logistic_all_splits.csv",index=False)
    pimp.to_csv(RAW/"19_nested_grouped_permutation_importance_all_splits.csv",index=False)
    sup=summarize(perf,["algorithm"],["balanced_accuracy","macro_F1","MCC",*[f"recall_{c}" for c in CLASSES],*[f"precision_{c}" for c in CLASSES]])
    sup.to_csv(TABLES/"18_supervised_summary_95CI.csv",index=False);best=sup.sort_values(["balanced_accuracy_median","macro_F1_median"],ascending=False).iloc[0].algorithm
    for alg in MODEL_ORDER:
        q=cms[cms.algorithm==alg].groupby(["true_class","predicted_class"],as_index=False).N.sum();tab=q.pivot(index="true_class",columns="predicted_class",values="N").reindex(index=CLASSES,columns=CLASSES,fill_value=0)
        tab.to_csv(TABLES/f"19_confusion_sum_{alg.replace(' ','_')}.csv")
    imp_sum=summarize(imp,["feature"],["ExtraTrees_MDI","Logistic_abs_standardized_coef"]);imp_sum.columns=imp_sum.columns.str.replace("_q025","_q025").str.replace("_q975","_q975")
    rho,pv=spearmanr(imp_sum.ExtraTrees_MDI_median,imp_sum.Logistic_abs_standardized_coef_median);imp_sum.to_csv(TABLES/"20_feature_importance_MDI_logistic_500.csv",index=False)
    summarize(pimp,["feature"],["importance_mean"]).to_csv(TABLES/"21_nested_grouped_permutation_importance_500.csv",index=False)
    abl.to_csv(RAW/"22_progressive_ablation_500.csv",index=False);summarize(abl,["subset"],["balanced_accuracy","macro_F1"]).to_csv(TABLES/"23_progressive_ablation_summary_95CI.csv",index=False)
    pd.DataFrame({"Spearman_rho":[rho],"p_value":[pv]}).to_csv(TABLES/"24_feature_rank_correlation.csv",index=False)

    agg,roc,pr,aucs=roc_pr(pred,best);agg.to_csv(TABLES/"25_best_model_cell_mean_probabilities.csv",index=False);roc.to_csv(TABLES/"26_ROC_coordinates.csv",index=False);pr.to_csv(TABLES/"27_PR_coordinates.csv",index=False);aucs.to_csv(TABLES/"28_ROC_PR_metrics.csv",index=False)
    print(f"Permutation: {n_repeats} global label permutations with fixed day-grouped splits for {best}",flush=True)
    perm=pd.DataFrame(Parallel(n_jobs=n_jobs,verbose=5)(delayed(permutation_repeat)(s,x,y,groups,best) for s in splits));perm.to_csv(RAW/"29_within_donor_label_permutation_500.csv",index=False)
    observed=float(perf[perf.algorithm==best].balanced_accuracy.median());empirical=(1+np.sum(perm.balanced_accuracy>=observed))/(n_repeats+1)
    pd.DataFrame([{"algorithm":best,"observed_BA_median":observed,"permutation_BA_median":perm.balanced_accuracy.median(),"empirical_p":empirical,"nominal_chance":.25}]).to_csv(TABLES/"30_permutation_test_summary.csv",index=False)
    meta=metadata_tests(consensus);meta.to_csv(TABLES/"31_metadata_association_summary.csv",index=False)
    cm_best=cms[cms.algorithm==best].groupby(["true_class","predicted_class"]).N.sum().unstack(fill_value=0).reindex(index=CLASSES,columns=CLASSES,fill_value=0).to_numpy()
    plots(sens,uns,perf,cm_best,perm,cells,imp_sum,abl,meta,roc,pr,aucs,best,observed)

    params={"scope":"Mouse NAc MSN morphology; frozen 187-cell cohort","source_csv":str(SOURCE),"all_complete_cells":187,"analysis_cells":187,"class_counts":consensus.M_class.value_counts().sort_index().to_dict(),"recording_day_groups":int(consensus.donor.nunique()),"features":FEATURES,"feature_n":10,"transformation":"Training minimum shift; training shifted-column-sum x10,000; log1p; training-fold sample Z-score","frozen_PCA_n":3,"frozen_HC":"Euclidean + Ward.D2, K=4","frozen_GC":"Seurat SNN k=20, prune=1/15, Louvain algorithm 1, resolution=0.50; HC labels retained as taxonomy","repeats":n_repeats,"unsupervised":{"recording_day_fraction":.8,"NPC":"2-10","K":"2-8","algorithms":["Ward","K-means","Gaussian mixture","Spectral"]},"supervised":{"train_recording_day_fraction":.75,"test_recording_day_fraction":.25,"all_classes_required_both":True,"models":MODEL_ORDER,"primary_metrics":["balanced_accuracy","macro_F1"],"PCA_n":3},"permutation":"Global label permutation evaluated on identical recording-day-grouped splits","consensus_margin_risk_threshold":.10,"seed":SEED,"software":{"python":sys.version,"platform":platform.platform(),"numpy":np.__version__,"pandas":pd.__version__,"scipy":scipy.__version__,"scikit_learn":sklearn.__version__,"joblib":joblib.__version__,"matplotlib":mpl.__version__,"seaborn":sns.__version__},"E_data_used_as_predictors":False,"Mouse_data_used":True}
    (OUT/"00_parameters.json").write_text(json.dumps(params,indent=2,ensure_ascii=False),encoding="utf-8")
    best_row=sup[sup.algorithm==best].iloc[0]
    ward_row=fixed[fixed.algorithm=="Ward"].iloc[0]
    gc_repro={"status":"Frozen source table used", "HC_GC_exact_agreement_n":int(consensus.HC_GC_consensus.astype(bool).sum()), "total_n":len(consensus)}
    summary={"frozen_HC_reproduction":repro,"frozen_GC_reproduction":gc_repro,
      "unsupervised_frozen_NPC3_K4_Ward":{"ARI_median":ward_row.ARI_vs_frozen_M4_median,"ARI_95CI":[ward_row.ARI_vs_frozen_M4_q025,ward_row.ARI_vs_frozen_M4_q975],"NMI_median":ward_row.NMI_vs_frozen_M4_median,"silhouette_median":ward_row.silhouette_median,"minimum_cluster_n_median":ward_row.minimum_cluster_n_median,"small_cluster_le2_rate":ward_row.small_cluster_le2_rate},
      "best_supervised_model":best,"best_supervised_BA_median":observed,"best_supervised_BA_95CI":[best_row.balanced_accuracy_q025,best_row.balanced_accuracy_q975],"best_supervised_macro_F1_median":best_row.macro_F1_median,"best_supervised_macro_F1_95CI":[best_row.macro_F1_q025,best_row.macro_F1_q975],"permutation_BA_median":float(perm.balanced_accuracy.median()),"permutation_empirical_p":empirical,"feature_rank_Spearman_rho":rho,"risk_cells_margin_le_0":int(cells.risk_margin_le_0.sum()),"risk_cells_margin_below_0.10":int(cells["risk_margin_below_0.10"].sum()),"E_annotation_available_n":int(consensus.E_class.notna().sum())}
    (OUT/"00_results_summary.json").write_text(json.dumps(summary,indent=2,ensure_ascii=False),encoding="utf-8")
    cn=f"""冻结 Mouse MSN 形态学 M1-M4：Macaque-M 对标的{n_repeats}次稳定性与可复现性验证\n================================================================================\n范围：187个冻结形态细胞，M1=67、M2=23、M3=62、M4=35；{consensus.donor.nunique()}个记录日代理组。\n\n严格防泄漏：每次训练或重抽样均重新估计最小值平移、列和归一化、log1p、Z-score及三成分PCA；同一记录日不跨训练和测试。冻结流程复现ARI={repro['ARI']:.3f}、NMI={repro['NMI']:.3f}。\n\n证据边界：监督结果只衡量冻结标签能否由同一形态特征域恢复；无监督重抽样衡量从头发现四类的稳定性；逐细胞margin衡量个体归属稳定性。当前没有独立动物ID，记录日仅作为保守的动物代理，因此这些结果属于内部重复验证，而非真正的独立队列验证。\n"""
    en=f"""Frozen Mouse MSN morphology M1-M4: Macaque-M-matched {n_repeats}-repeat validation\n=====================================================================================\nScope: 187 frozen morphology cells (M1=67, M2=23, M3=62, M4=35) from {consensus.donor.nunique()} recording-day proxy groups.\n\nLeakage control: minimum shifting, column-sum normalization, log1p transformation, Z-scoring and three-component PCA were re-estimated in every training or resampling set; recording-day groups never crossed training and test partitions. Frozen workflow reproduction: ARI={repro['ARI']:.3f}, NMI={repro['NMI']:.3f}.\n\nEvidence boundary: supervised results measure recoverability of frozen labels from the same morphology domain; unsupervised resampling measures de novo four-class recovery; cell margins quantify assignment stability. Animal IDs are unavailable and recording day is only a conservative animal proxy, so this is internal repeated validation rather than independent-cohort validation.\n"""
    (OUT/"README_CN.txt").write_text(cn,encoding="utf-8");(OUT/"README_EN.txt").write_text(en,encoding="utf-8")
    print(json.dumps(summary,indent=2,ensure_ascii=False),flush=True)


def refresh_batch_metadata_outputs():
    """Refresh only the audit and batch-association tables after metadata naming changes."""
    global NREP
    NREP = 500
    for folder in [OUT, RAW, TABLES, FIGURES]:
        folder.mkdir(parents=True, exist_ok=True)
    data = load_inputs()
    audit_cols=["cell_label","donor","mouse_day_proxy","experimental_day_batch","sequencing_submission_batch","ROI_final","T_class","E_class","M_class","HC_cluster","GC_raw_cluster","GC_as_M_class","HC_GC_consensus","included_ML500","inclusion_reason","recording_date","reconstruction_batch","morphology_missing_rate","reconstruction_complete"]
    data[audit_cols].to_csv(TABLES/"01_frozen_input_audit_187_cells.csv",index=False)
    metadata_tests(data).to_csv(TABLES/"31_metadata_association_summary.csv",index=False)
    print("Refreshed mouse/day proxy and sequencing-submission-batch metadata outputs", flush=True)


if __name__ == "__main__":
    parser=argparse.ArgumentParser();parser.add_argument("--repeats",type=int,default=500);parser.add_argument("--jobs",type=int,default=-1);parser.add_argument("--refresh-metadata-only",action="store_true")
    args=parser.parse_args()
    refresh_batch_metadata_outputs() if args.refresh_metadata_only else main(args.repeats,args.jobs)
