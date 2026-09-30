#!/usr/bin/env python
"""500-repeat validation of Macaque E taxonomy and exploratory M taxonomy.

Scope is fixed to MSN cells from Ca + Pu + NAC.  E uses the frozen E4 labels;
M is label-blind and scans K=2..8.  Every resample refits scaling and PCA.
"""

from __future__ import annotations

from pathlib import Path
import json
import warnings
from zipfile import ZipFile

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from joblib import Parallel, delayed
from scipy.optimize import linear_sum_assignment
from sklearn.cluster import AgglomerativeClustering, KMeans, SpectralClustering
from sklearn.decomposition import PCA
from sklearn.ensemble import ExtraTreesClassifier, HistGradientBoostingClassifier, RandomForestClassifier
from sklearn.metrics import (
    adjusted_rand_score, balanced_accuracy_score, f1_score,
    matthews_corrcoef, normalized_mutual_info_score, silhouette_score,
)
from sklearn.mixture import GaussianMixture
from sklearn.model_selection import GroupShuffleSplit
from sklearn.neighbors import KNeighborsClassifier
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import LabelEncoder, StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs" / "Macaque_EM_ML500"
E_FILE = ROOT / "outputs" / "R3_panels" / "00_res3_main_cell_assignments_and_tsne.csv"
M_FILE = ROOT / "outputs" / "Macaque_M_QC_MSN_CaPuNAc" / "01_MSN_CaPuNAc_cell_QC.csv"
ZIP_FILE = Path(r"C:\Users\53461\Downloads\Macaque-PatchSeq-BG.zip")
N_REPEATS = 500
SEED = 20260909

E19 = [
    "Epsy_width_rheo", "Epsy_fast_trough_v_rheo", "Epsy_peak_deltav_rheo",
    "Epsy_peak_v_rheo", "Epsy_postap_slope_rheo", "Epsy_threshold_v_rheo",
    "Epsy_trough_t_rheo", "Epsy_trough_v_rheo", "Epsy_upstroke_downstroke_ratio_rheo",
    "Epsy_ahp_delay_5spike", "Epsy_ahp_delay_ratio_5spike", "Epsy_postap_slope_hero",
    "Epsy_trough_t_hero", "Epsy_downstroke_adapt_ratio", "Epsy_peak_v_adapt_ratio",
    "Epsy_threshold_v_adapt_ratio", "Epsy_upstroke_adapt_ratio",
    "Epsy_width_adapt_ratio", "Epsy_threshold_v_short_square",
]


def fit_e_transform(x: pd.DataFrame):
    arr = np.empty(x.shape, float)
    params = []
    for j, c in enumerate(x.columns):
        v = x[c].to_numpy(float)
        if "ratio" in c.lower():
            if np.any(v <= 0):
                raise ValueError(f"Non-positive value in ratio feature {c}")
            z = np.log2(v)
            kind, a, b = "log2", 0.0, 1.0
        else:
            a = float(np.min(v)); b = float(np.sum(v - a))
            if b <= 0:
                raise ValueError(f"Constant E feature {c}")
            z = np.log1p((v - a) / b * 10000.0)
            kind = "shift_sum_log1p"
        mu, sd = float(z.mean()), float(z.std(ddof=0))
        arr[:, j] = (z - mu) / (sd if sd > 0 else 1.0)
        params.append((kind, a, b, mu, sd if sd > 0 else 1.0))
    return arr, params


def apply_e_transform(x: pd.DataFrame, params):
    arr = np.empty(x.shape, float)
    for j, (c, p) in enumerate(zip(x.columns, params)):
        kind, a, b, mu, sd = p
        v = x[c].to_numpy(float)
        if kind == "log2":
            z = np.log2(np.clip(v, np.finfo(float).tiny, None))
        else:
            # A held-out value may be below the training minimum; map it to zero.
            z = np.log1p(np.clip(v - a, 0, None) / b * 10000.0)
        arr[:, j] = (z - mu) / sd
    return arr


def preprocess(x: pd.DataFrame, kind: str, n_pc: int = 3):
    if kind == "E":
        z, _ = fit_e_transform(x)
    else:
        z = StandardScaler().fit_transform(x)
    return PCA(n_components=min(n_pc, z.shape[1]), svd_solver="full").fit_transform(z)


def best_map_accuracy(ref, pred):
    r = LabelEncoder().fit_transform(ref)
    p = LabelEncoder().fit_transform(pred)
    tab = pd.crosstab(r, p).to_numpy()
    rr, cc = linear_sum_assignment(-tab)
    return float(tab[rr, cc].sum() / len(r))


def cluster_jaccards(ref, pred):
    rvals, pvals = np.unique(ref), np.unique(pred)
    score = np.zeros((len(rvals), len(pvals)))
    for i, rv in enumerate(rvals):
        a = ref == rv
        for j, pv in enumerate(pvals):
            b = pred == pv
            score[i, j] = np.sum(a & b) / np.sum(a | b)
    rr, cc = linear_sum_assignment(-score)
    out = {str(rvals[i]): float(score[i, j]) for i, j in zip(rr, cc)}
    return out


def all_clusterers(pc, k, seed):
    nn = min(15, len(pc) - 1)
    return {
        "Ward": AgglomerativeClustering(n_clusters=k, linkage="ward").fit_predict(pc),
        "KMeans": KMeans(n_clusters=k, n_init=20, random_state=seed).fit_predict(pc),
        "GMM": GaussianMixture(n_components=k, covariance_type="full", reg_covar=1e-5,
                               n_init=5, random_state=seed).fit_predict(pc),
        "Spectral": SpectralClustering(n_clusters=k, affinity="nearest_neighbors",
                                       n_neighbors=nn, assign_labels="kmeans",
                                       n_init=10, random_state=seed).fit_predict(pc),
    }


def one_unsup_repeat(rep, x, groups, ref, kind, ks):
    rng = np.random.default_rng(SEED + rep)
    ug = np.unique(groups)
    chosen = rng.choice(ug, size=max(2, int(np.ceil(0.8 * len(ug)))), replace=False)
    idx = np.flatnonzero(np.isin(groups, chosen))
    pc = preprocess(x.iloc[idx], kind, 3)
    rows = []
    labels_for_consensus = {}
    for k in ks:
        labs = all_clusterers(pc, k, SEED + rep)
        labels_for_consensus[k] = (idx, labs["Ward"])
        ref_k = ref[k][idx]
        for alg, pred in labs.items():
            js = cluster_jaccards(ref_k, pred)
            rows.append({
                "repeat": rep, "K": k, "algorithm": alg, "n_cells": len(idx),
                "ARI_vs_reference": adjusted_rand_score(ref_k, pred),
                "NMI_vs_reference": normalized_mutual_info_score(ref_k, pred),
                "mapped_accuracy": best_map_accuracy(ref_k, pred),
                "min_cluster_Jaccard": min(js.values()),
                "mean_cluster_Jaccard": np.mean(list(js.values())),
            })
    return rows, labels_for_consensus


def consensus_from_repeats(n, repeats, k):
    coobs = np.zeros((n, n), dtype=np.uint16)
    coassign = np.zeros((n, n), dtype=np.uint16)
    for _, labdict in repeats:
        idx, lab = labdict[k]
        coobs[np.ix_(idx, idx)] += 1
        for c in np.unique(lab):
            ii = idx[lab == c]
            coassign[np.ix_(ii, ii)] += 1
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.divide(coassign, coobs, out=np.full_like(coassign, np.nan, dtype=float), where=coobs > 0)


def valid_group_splits(y, groups, n=N_REPEATS):
    out = []
    seed = SEED
    while len(out) < n and seed < SEED + 20000:
        tr, te = next(GroupShuffleSplit(n_splits=1, test_size=0.25, random_state=seed).split(y, y, groups))
        if len(np.unique(y[tr])) == len(np.unique(y)) and len(np.unique(y[te])) == len(np.unique(y)):
            out.append((len(out), seed, tr, te))
        seed += 1
    if len(out) < n:
        raise RuntimeError(f"Only {len(out)} valid donor-grouped splits")
    return out


def model_set(seed):
    return {
        "Multinomial logistic": LogisticRegression(max_iter=3000, class_weight="balanced", C=1.0),
        "Linear SVM": SVC(kernel="linear", class_weight="balanced", C=1.0),
        "RBF SVM": SVC(kernel="rbf", class_weight="balanced", C=1.0),
        "Random forest": RandomForestClassifier(n_estimators=150, class_weight="balanced", max_features="sqrt", n_jobs=1, random_state=seed),
        "Extra trees": ExtraTreesClassifier(n_estimators=150, class_weight="balanced", max_features="sqrt", n_jobs=1, random_state=seed),
        "kNN": KNeighborsClassifier(n_neighbors=7, weights="distance"),
        "HistGradientBoosting": HistGradientBoostingClassifier(max_iter=150, learning_rate=0.06, random_state=seed),
    }


def one_supervised(rep, seed, tr, te, x, y, kind):
    if kind == "E":
        ztr, par = fit_e_transform(x.iloc[tr]); zte = apply_e_transform(x.iloc[te], par)
    else:
        sc = StandardScaler().fit(x.iloc[tr]); ztr = sc.transform(x.iloc[tr]); zte = sc.transform(x.iloc[te])
    pca = PCA(n_components=3, svd_solver="full").fit(ztr)
    a, b = pca.transform(ztr), pca.transform(zte)
    rows = []
    for name, model in model_set(seed).items():
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            model.fit(a, y[tr]); pred = model.predict(b)
        rows.append({
            "repeat": rep, "seed": seed, "algorithm": name, "n_train": len(tr), "n_test": len(te),
            "balanced_accuracy": balanced_accuracy_score(y[te], pred),
            "macro_F1": f1_score(y[te], pred, average="macro"),
            "MCC": matthews_corrcoef(y[te], pred),
        })
    # Feature importance is collected without PCA so names remain interpretable.
    lr = LogisticRegression(max_iter=3000, class_weight="balanced", C=1.0).fit(ztr, y[tr])
    et = ExtraTreesClassifier(n_estimators=150, class_weight="balanced", n_jobs=1,
                              max_features="sqrt", random_state=seed).fit(ztr, y[tr])
    importance = (np.mean(np.abs(lr.coef_), axis=0), et.feature_importances_)
    return rows, importance


def one_permutation(rep, seed, tr, te, x, y, groups, kind):
    rng = np.random.default_rng(seed + 100000)
    # Permute donor-level label blocks by shuffling labels globally; split stays fixed.
    yp = rng.permutation(y)
    if kind == "E":
        ztr, par = fit_e_transform(x.iloc[tr]); zte = apply_e_transform(x.iloc[te], par)
    else:
        sc = StandardScaler().fit(x.iloc[tr]); ztr = sc.transform(x.iloc[tr]); zte = sc.transform(x.iloc[te])
    pca = PCA(n_components=3, svd_solver="full").fit(ztr)
    model = LogisticRegression(max_iter=3000, class_weight="balanced")
    model.fit(pca.transform(ztr), yp[tr])
    return {"repeat": rep, "balanced_accuracy": balanced_accuracy_score(yp[te], model.predict(pca.transform(zte)))}


def summarize_metric(df, group_cols, metrics):
    rows = []
    for key, g in df.groupby(group_cols, observed=True):
        key = (key,) if not isinstance(key, tuple) else key
        row = dict(zip(group_cols, key))
        for m in metrics:
            row[f"{m}_median"] = g[m].median()
            row[f"{m}_q025"] = g[m].quantile(.025)
            row[f"{m}_q975"] = g[m].quantile(.975)
        rows.append(row)
    return pd.DataFrame(rows)


def plot_consensus(mat, labels, path, title):
    order = np.argsort(labels.astype(str))
    fig, ax = plt.subplots(figsize=(3.1, 2.8))
    sns.heatmap(mat[np.ix_(order, order)], cmap="magma", vmin=0, vmax=1, ax=ax,
                xticklabels=False, yticklabels=False, cbar_kws={"label": "Co-clustering probability"})
    ax.set_title(title)
    fig.tight_layout(); fig.savefig(path, dpi=600); plt.close(fig)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    mpl.rcParams.update({"font.family": "Arial", "font.size": 7, "axes.linewidth": .7})

    e_labels = pd.read_csv(E_FILE).drop(columns=E19)
    with ZipFile(ZIP_FILE) as zf:
        e_raw = pd.read_csv(zf.open("Data/cell_metadata_AllCell.csv"), usecols=["cell_label", *E19])
    e = e_labels.merge(e_raw, on="cell_label", how="left", validate="one_to_one")
    e["Consensus"] = e["Consensus"].astype(str).str.lower().eq("true")
    ex = e[E19].apply(pd.to_numeric, errors="raise")
    ey = e["HC_class"].astype(str).to_numpy()
    eg = e["donor_label"].astype(str).to_numpy()
    epc = preprocess(ex, "E", 3)
    eward = AgglomerativeClustering(n_clusters=4, linkage="ward").fit_predict(epc)
    e_reproduction_ari = adjusted_rand_score(ey, eward)
    if e_reproduction_ari < .999:
        raise RuntimeError(f"Frozen E Ward reproduction failed: ARI={e_reproduction_ari}")
    eref = {4: ey}

    m = pd.read_csv(M_FILE)
    m = m[m["Complete_primary21"].astype(str).str.lower().eq("true")].reset_index(drop=True)
    primary21 = [c for c in m.columns if c.startswith("basal_dendrite_") or c.startswith("soma_") or c.startswith("3_Sholl_")]
    primary19 = [c for c in primary21 if c not in {"3_Sholl_PC0", "basal_dendrite_total_surface_area"}]
    direct18 = [c for c in primary19 if c != "3_Sholl_PC1"]
    panels = {"M19_primary": primary19, "M21_sensitivity": primary21, "M18_direct": direct18}
    mg = m["Lib_donor_label"].astype(str).to_numpy()

    # E: frozen K=4 donor-level resampling.
    print("E unsupervised 500 repeats", flush=True)
    eruns = Parallel(n_jobs=-1, verbose=5)(delayed(one_unsup_repeat)(i, ex, eg, eref, "E", [4]) for i in range(N_REPEATS))
    e_uns = pd.DataFrame([r for rows, _ in eruns for r in rows])
    e_uns.to_csv(OUT / "01_E4_donor_resampling_500.csv", index=False)

    # M: three feature panels, K=2..8, donor-level resampling.
    m_uns_all, m_runs_by_panel, mrefs = [], {}, {}
    for panel, feats in panels.items():
        print(f"{panel} unsupervised 500 x K2-8", flush=True)
        x = m[feats].apply(pd.to_numeric, errors="raise")
        pc = preprocess(x, "M", 3)
        refs = {k: AgglomerativeClustering(n_clusters=k, linkage="ward").fit_predict(pc) for k in range(2, 9)}
        mrefs[panel] = refs
        runs = Parallel(n_jobs=-1, verbose=5)(delayed(one_unsup_repeat)(i, x, mg, refs, "M", range(2, 9)) for i in range(N_REPEATS))
        d = pd.DataFrame([r for rows, _ in runs for r in rows]); d.insert(0, "panel", panel)
        m_uns_all.append(d); m_runs_by_panel[panel] = runs
    m_uns = pd.concat(m_uns_all, ignore_index=True)
    m_uns.to_csv(OUT / "02_M_K2-8_three_panels_donor_resampling_500.csv", index=False)

    e_sum = summarize_metric(e_uns, ["K", "algorithm"], ["ARI_vs_reference", "min_cluster_Jaccard", "mapped_accuracy"])
    m_sum = summarize_metric(m_uns, ["panel", "K", "algorithm"], ["ARI_vs_reference", "min_cluster_Jaccard", "mapped_accuracy"])
    e_sum.to_csv(OUT / "03_E4_stability_summary_95CI.csv", index=False)
    m_sum.to_csv(OUT / "04_M_K_selection_stability_summary_95CI.csv", index=False)

    # Automatic M choice: primary panel, Ward median ARI + min-Jaccard, with cluster-size guard.
    cand = m_sum[(m_sum.panel == "M19_primary") & (m_sum.algorithm == "Ward")].copy()
    sizes, sils = {}, {}
    mx = m[primary19]
    mpc = preprocess(mx, "M", 3)
    for k, lab in mrefs["M19_primary"].items():
        sizes[k] = int(pd.Series(lab).value_counts().min())
        sils[k] = silhouette_score(mpc, lab)
    cand["minimum_full_cluster_n"] = cand.K.map(sizes)
    cand["silhouette"] = cand.K.map(sils)
    cand["selection_score"] = .55 * cand.ARI_vs_reference_median + .35 * cand.min_cluster_Jaccard_median + .10 * ((cand.silhouette + 1) / 2)
    eligible = cand[cand.minimum_full_cluster_n >= 8]
    selected_k = int((eligible if len(eligible) else cand).sort_values("selection_score", ascending=False).iloc[0].K)
    cand.to_csv(OUT / "05_M_primary19_automatic_K_selection.csv", index=False)

    # Consensus matrices.
    econs = consensus_from_repeats(len(e), eruns, 4)
    mcons = consensus_from_repeats(len(m), m_runs_by_panel["M19_primary"], selected_k)
    np.save(OUT / "06_E4_consensus_matrix.npy", econs)
    np.save(OUT / "07_M_selectedK_consensus_matrix.npy", mcons)
    plot_consensus(econs, ey, OUT / "08_E4_consensus_heatmap.png", "E4: 500 donor-resamples")
    plot_consensus(mcons, mrefs["M19_primary"][selected_k], OUT / "09_M_consensus_heatmap.png", f"M19, K={selected_k}: 500 donor-resamples")

    # Donor-grouped supervised recoverability (E consensus cells only; M exploratory labels).
    sup_frames, imp_frames, perm_frames = [], [], []
    tasks = [
        ("E4_frozen_consensus", ex.loc[e.Consensus].reset_index(drop=True), ey[e.Consensus.to_numpy()], eg[e.Consensus.to_numpy()], "E", E19),
        (f"M19_exploratory_K{selected_k}", mx.reset_index(drop=True), mrefs["M19_primary"][selected_k], mg, "M", primary19),
    ]
    for name, x, y, groups, kind, features in tasks:
        print(f"{name} supervised + permutation 500 repeats", flush=True)
        splits = valid_group_splits(np.asarray(y), np.asarray(groups))
        sr = Parallel(n_jobs=-1, verbose=5)(delayed(one_supervised)(*s, x, np.asarray(y), kind) for s in splits)
        sdf = pd.DataFrame([r for rows, _ in sr for r in rows]); sdf.insert(0, "analysis", name); sup_frames.append(sdf)
        lr_imp = np.vstack([imp[0] for _, imp in sr]); et_imp = np.vstack([imp[1] for _, imp in sr])
        imp = pd.DataFrame({
            "analysis": name, "feature": features,
            "logistic_abs_coef_median": np.median(lr_imp, axis=0),
            "logistic_abs_coef_q025": np.quantile(lr_imp, .025, axis=0),
            "logistic_abs_coef_q975": np.quantile(lr_imp, .975, axis=0),
            "extra_trees_importance_median": np.median(et_imp, axis=0),
            "extra_trees_importance_q025": np.quantile(et_imp, .025, axis=0),
            "extra_trees_importance_q975": np.quantile(et_imp, .975, axis=0),
        }); imp_frames.append(imp)
        pr = Parallel(n_jobs=-1, verbose=5)(delayed(one_permutation)(*s, x, np.asarray(y), np.asarray(groups), kind) for s in splits)
        pdf = pd.DataFrame(pr); pdf.insert(0, "analysis", name); perm_frames.append(pdf)

    sup = pd.concat(sup_frames, ignore_index=True); sup.to_csv(OUT / "10_grouped_supervised_500_all_models.csv", index=False)
    sup_sum = summarize_metric(sup, ["analysis", "algorithm"], ["balanced_accuracy", "macro_F1", "MCC"])
    sup_sum.to_csv(OUT / "11_grouped_supervised_summary_95CI.csv", index=False)
    pd.concat(imp_frames, ignore_index=True).to_csv(OUT / "12_feature_importance_500_summary.csv", index=False)
    perm = pd.concat(perm_frames, ignore_index=True); perm.to_csv(OUT / "13_label_permutation_500.csv", index=False)

    # Compact publication/QC overview.
    fig, axes = plt.subplots(1, 3, figsize=(9.2, 2.7))
    sns.boxplot(data=e_uns, x="algorithm", y="ARI_vs_reference", color="#56B4E9", ax=axes[0], fliersize=1)
    axes[0].set_title("Frozen E4 stability"); axes[0].tick_params(axis="x", rotation=35)
    mw = m_uns[(m_uns.panel == "M19_primary") & (m_uns.algorithm == "Ward")]
    sns.boxplot(data=mw, x="K", y="ARI_vs_reference", color="#CC79A7", ax=axes[1], fliersize=1)
    axes[1].axvline(selected_k - 2, color="#D55E00", ls="--", lw=.8); axes[1].set_title("Exploratory M: Ward stability")
    sns.boxplot(data=sup, x="algorithm", y="balanced_accuracy", hue="analysis", ax=axes[2], fliersize=1)
    axes[2].set_title("Donor-held-out recoverability"); axes[2].tick_params(axis="x", rotation=50); axes[2].legend(fontsize=5)
    for ax in axes: ax.spines[["top", "right"]].set_visible(False); ax.set_ylim(-.02, 1.02)
    fig.tight_layout(); fig.savefig(OUT / "14_Macaque_EM_ML500_overview.png", dpi=600); fig.savefig(OUT / "14_Macaque_EM_ML500_overview.pdf"); plt.close(fig)

    summary = {
        "repeats": N_REPEATS,
        "scope": "Macaque MSN; Ca + Pu + NAC",
        "E": {"all_n": len(e), "consensus_n": int(e.Consensus.sum()), "donors": int(pd.Series(eg).nunique()),
              "features": 19, "frozen_K": 4, "full_reproduction_ARI": e_reproduction_ari},
        "M": {"complete_n": len(m), "donors": int(pd.Series(mg).nunique()), "primary_features": len(primary19),
              "sensitivity_features": len(primary21), "direct_features": len(direct18), "automatic_selected_K": selected_k,
              "classification_status": "exploratory; not frozen"},
        "resampling_unit": "donor (80% without replacement)",
        "supervised_test": "25% donor-group holdout; preprocessing and PCA refit within training only",
        "seed": SEED,
    }
    (OUT / "00_run_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    readme = f"""Macaque E/M 500-repeat stability and reproducibility analysis

Scope: MSN only, Ca + Pu + NAC only.
E: frozen E4 taxonomy; 390 complete cells, 368 consensus cells, 55 donors, 19 features.
M: 126 morphology-complete cells, 42 donors. M is de-novo/exploratory, not frozen.
Repeats: {N_REPEATS}; random seed: {SEED}.

Unsupervised stability:
- 80% of donors sampled without replacement in every repeat.
- Scaling/transformation and PCA are refit every time.
- E tests frozen K=4. M scans K=2..8 in primary19, sensitivity21, and direct18 panels.
- Algorithms: Ward, K-means, full-covariance GMM, spectral clustering.
- Metrics: ARI, NMI, mapped accuracy, matched per-cluster Jaccard.

Supervised recoverability:
- 500 valid 75/25 donor-group splits; all classes required in train and test.
- E uses only the 368 HC-GC consensus cells.
- M predicts the automatically selected full-data Ward pseudo-labels; this is not independent biological validation.
- Models: multinomial logistic, linear/RBF SVM, random forest, extra trees, kNN, histogram gradient boosting.
- Metrics: balanced accuracy, macro-F1, MCC; 2.5%-97.5% empirical intervals.
- A 500-repeat permuted-label logistic baseline is included.

Automatic exploratory M selection: K={selected_k}. Selection uses Ward median ARI, minimum matched-cluster Jaccard,
silhouette, and minimum full-cluster size >=8. It must not be described as a frozen Macaque M taxonomy without
external biological validation and a prespecified freeze decision.
"""
    (OUT / "README_methods_and_limits.txt").write_text(readme, encoding="utf-8")
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
