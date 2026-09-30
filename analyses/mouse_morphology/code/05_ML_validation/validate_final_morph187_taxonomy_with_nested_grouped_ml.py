#!/usr/bin/env python
"""Rigorous ML recovery test of the final n=187 M1-M4 morphology taxonomy.

This is a reproducibility/recoverability analysis of HC-derived labels, not an
independent biological validation. Recording day is used as the grouping proxy
because one mouse was recorded per day and animal IDs are unavailable.
"""

from __future__ import annotations

import json
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import bootstrap
from sklearn.base import clone
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.ensemble import ExtraTreesClassifier, RandomForestClassifier
from sklearn.inspection import permutation_importance
from sklearn.linear_model import LogisticRegression
from sklearn.manifold import TSNE
from sklearn.metrics import (
    accuracy_score, adjusted_rand_score, balanced_accuracy_score,
    confusion_matrix, f1_score, precision_recall_fscore_support,
    silhouette_score,
)
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.neighbors import KNeighborsClassifier
from sklearn.decomposition import PCA
from sklearn.svm import SVC


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "outputs" / "morph_qc" / "final_morph187_NPC3_HCK4_GCres0p50_figures" / "00_final_morph187_cell_assignments.csv"
OUT = ROOT / "outputs" / "morph_qc" / "final_morph187_NPC3_HCK4_GCres0p50_ML_validation"
ID = "MSN_unique_ID"
FEATURES = [
    "M_soma_circularity_index", "M_soma_aspect_ratio", "M_cell_max_radial_dist",
    "M_total_number_of_neurites", "M_basal_dendrite_avg_tortuosity",
    "M_Total_neurite_length_(sections)", "M_Number_of_bifurcation_points",
    "M_Maximum_branch_order", "M_trunk_angle_min", "M_trunk_angle_max",
]
DISPLAY = {
    "M_soma_circularity_index": "Soma circularity",
    "M_soma_aspect_ratio": "Soma aspect ratio",
    "M_cell_max_radial_dist": "Max radial distance",
    "M_total_number_of_neurites": "Primary neurite number",
    "M_basal_dendrite_avg_tortuosity": "Mean tortuosity",
    "M_Total_neurite_length_(sections)": "Total neurite length",
    "M_Number_of_bifurcation_points": "Bifurcation points",
    "M_Maximum_branch_order": "Maximum branch order",
    "M_trunk_angle_min": "Minimum trunk angle",
    "M_trunk_angle_max": "Maximum trunk angle",
}
CLASSES = ["M1", "M2", "M3", "M4"]
COLORS = {"M1": "#00468B", "M2": "#42B540", "M3": "#ED0000", "M4": "#0099B4"}
SEED = 20260826
OUTER_REPEATS = 5
OUTER_SPLITS = 5
INNER_SPLITS = 4


def configure() -> None:
    mpl.rcParams.update({
        "font.family": "Arial", "font.size": 4, "axes.linewidth": 0.75,
        "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
        "savefig.facecolor": "white", "figure.facecolor": "white",
    })


def save(fig: plt.Figure, stem: str, dpi: int = 900) -> None:
    fig.savefig(OUT / f"{stem}.png", dpi=dpi, bbox_inches="tight", pad_inches=0.025)
    fig.savefig(OUT / f"{stem}.pdf", bbox_inches="tight", pad_inches=0.025)
    fig.savefig(OUT / f"{stem}.svg", bbox_inches="tight", pad_inches=0.025)
    plt.close(fig)


def preprocess_train_test(train: np.ndarray, test: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    minimum = np.nanmin(train, axis=0)
    shifted_train = np.maximum(train - minimum, 0)
    shifted_test = np.maximum(test - minimum, 0)
    total = np.nansum(shifted_train, axis=0)
    total[~np.isfinite(total) | (total <= 0)] = 1
    log_train = np.log1p(shifted_train / total * 10000)
    log_test = np.log1p(shifted_test / total * 10000)
    mean = np.nanmean(log_train, axis=0)
    sd = np.nanstd(log_train, axis=0, ddof=1)
    sd[~np.isfinite(sd) | (sd <= 0)] = 1
    return np.nan_to_num((log_train - mean) / sd), np.nan_to_num((log_test - mean) / sd)


def model_specs() -> list[dict]:
    specs = []
    for c in [0.1, 1.0, 10.0]:
        specs.append({"family": "Regularized logistic", "name": f"C={c:g}", "model": LogisticRegression(
            C=c, penalty="l2", class_weight="balanced", solver="lbfgs", max_iter=5000, random_state=SEED)})
    for c in [0.5, 2.0, 8.0]:
        for gamma in ["scale", 0.1]:
            specs.append({"family": "RBF SVM", "name": f"C={c:g}, gamma={gamma}", "model": SVC(
                C=c, gamma=gamma, kernel="rbf", class_weight="balanced", probability=True, random_state=SEED)})
    for k in [3, 5, 9]:
        specs.append({"family": "kNN", "name": f"k={k}", "model": KNeighborsClassifier(
            n_neighbors=k, weights="distance", p=2)})
    for leaf in [1, 2, 4]:
        specs.append({"family": "Random forest", "name": f"leaf={leaf}", "model": RandomForestClassifier(
            n_estimators=300, max_features="sqrt", min_samples_leaf=leaf,
            class_weight="balanced_subsample", n_jobs=-1, random_state=SEED)})
        specs.append({"family": "Extra trees", "name": f"leaf={leaf}", "model": ExtraTreesClassifier(
            n_estimators=300, max_features="sqrt", min_samples_leaf=leaf,
            class_weight="balanced", n_jobs=-1, random_state=SEED)})
    specs.append({"family": "Shrinkage LDA", "name": "auto shrinkage", "model": LinearDiscriminantAnalysis(
        solver="lsqr", shrinkage="auto")})
    return specs


def aligned_proba(model, x: np.ndarray) -> np.ndarray:
    raw = model.predict_proba(x)
    out = np.zeros((len(x), len(CLASSES)), dtype=float)
    for j, cls in enumerate(model.classes_):
        out[:, CLASSES.index(str(cls))] = raw[:, j]
    return out


def inner_scores(specs, x, y, groups, train_idx, seed):
    scores = np.zeros(len(specs), dtype=float)
    cv = StratifiedGroupKFold(n_splits=INNER_SPLITS, shuffle=True, random_state=seed)
    for si, spec in enumerate(specs):
        fold_scores = []
        for itr, iva in cv.split(x[train_idx], y[train_idx], groups[train_idx]):
            tr, va = train_idx[itr], train_idx[iva]
            ztr, zva = preprocess_train_test(x[tr], x[va])
            model = clone(spec["model"])
            model.fit(ztr, y[tr])
            fold_scores.append(balanced_accuracy_score(y[va], model.predict(zva)))
        scores[si] = np.mean(fold_scores)
    return scores


def metrics(y_true, y_pred) -> dict:
    return {
        "accuracy": accuracy_score(y_true, y_pred),
        "balanced_accuracy": balanced_accuracy_score(y_true, y_pred),
        "macro_F1": f1_score(y_true, y_pred, average="macro"),
        "ARI": adjusted_rand_score(y_true, y_pred),
    }


def grouped_bootstrap_ci(y_true, y_pred, groups, iterations=3000):
    rng = np.random.default_rng(SEED)
    unique = np.unique(groups)
    values = defaultdict(list)
    for _ in range(iterations):
        chosen = rng.choice(unique, size=len(unique), replace=True)
        idx = np.concatenate([np.where(groups == g)[0] for g in chosen])
        if len(np.unique(y_true[idx])) < len(CLASSES):
            continue
        m = metrics(y_true[idx], y_pred[idx])
        for key, value in m.items():
            values[key].append(value)
    return {key: [float(np.quantile(v, .025)), float(np.quantile(v, .975))] for key, v in values.items()}


def benchmark_nested(x, y, groups, specs):
    family_rows, selection_rows = [], []
    probabilities = np.zeros((len(y), len(CLASSES)), dtype=float)
    prediction_votes = np.zeros((len(y), len(CLASSES)), dtype=int)
    test_counts = np.zeros(len(y), dtype=int)
    outer_splits = []
    families = list(dict.fromkeys(s["family"] for s in specs))

    for repeat in range(OUTER_REPEATS):
        outer = StratifiedGroupKFold(n_splits=OUTER_SPLITS, shuffle=True, random_state=SEED + repeat)
        for fold, (tr, te) in enumerate(outer.split(x, y, groups), 1):
            outer_splits.append((repeat, fold, tr.copy(), te.copy()))
            scores = inner_scores(specs, x, y, groups, tr, SEED + 1000 + repeat * 10 + fold)
            overall_best = int(np.argmax(scores))

            family_predictions = {}
            family_probabilities = {}
            for family in families:
                candidates = [i for i, s in enumerate(specs) if s["family"] == family]
                best = max(candidates, key=lambda i: (scores[i], -i))
                ztr, zte = preprocess_train_test(x[tr], x[te])
                model = clone(specs[best]["model"])
                model.fit(ztr, y[tr])
                pred = model.predict(zte)
                proba = aligned_proba(model, zte)
                family_predictions[family] = pred
                family_probabilities[family] = proba
                row = {"Repeat": repeat + 1, "Fold": fold, "Family": family,
                       "Selected_hyperparameters": specs[best]["name"], "Inner_balanced_accuracy": scores[best],
                       "N_train": len(tr), "N_test": len(te),
                       "Train_days": len(np.unique(groups[tr])), "Test_days": len(np.unique(groups[te]))}
                row.update(metrics(y[te], pred)); family_rows.append(row)

            chosen_family = specs[overall_best]["family"]
            chosen_pred = family_predictions[chosen_family]
            chosen_proba = family_probabilities[chosen_family]
            probabilities[te] += chosen_proba
            test_counts[te] += 1
            for local, idx in enumerate(te):
                prediction_votes[idx, CLASSES.index(str(chosen_pred[local]))] += 1
            selection_rows.append({
                "Repeat": repeat + 1, "Fold": fold, "Selected_family": chosen_family,
                "Selected_hyperparameters": specs[overall_best]["name"],
                "Inner_balanced_accuracy": scores[overall_best], **metrics(y[te], chosen_pred),
            })

    probabilities /= test_counts[:, None]
    aggregated_pred = np.array(CLASSES)[np.argmax(probabilities, axis=1)]
    vote_agreement = prediction_votes.max(axis=1) / test_counts
    return (pd.DataFrame(family_rows), pd.DataFrame(selection_rows), probabilities,
            aggregated_pred, vote_agreement, outer_splits)


def plot_benchmark(family_df):
    order = (family_df.groupby("Family")["balanced_accuracy"].mean().sort_values(ascending=False).index.tolist())
    rng = np.random.default_rng(SEED)
    fig, ax = plt.subplots(figsize=(3.35, 2.15))
    for i, family in enumerate(order):
        vals = family_df.loc[family_df.Family.eq(family), "balanced_accuracy"].to_numpy()
        ax.scatter(np.full(len(vals), i) + rng.normal(0, .055, len(vals)), vals,
                   s=6, c="#A0A0A0", edgecolors="none", alpha=.75)
        mean = vals.mean(); ci = np.quantile(vals, [.025, .975])
        ax.errorbar(i, mean, yerr=[[mean-ci[0]], [ci[1]-mean]], fmt="o", ms=3.7,
                    color="#4C78A8", ecolor="#4C78A8", elinewidth=.8, capsize=2)
    ax.set_xticks(range(len(order)), order, rotation=28, ha="right")
    ax.set_ylim(.45, 1.01); ax.set_ylabel("Outer-fold balanced accuracy")
    ax.spines[["top", "right"]].set_visible(False)
    save(fig, "01_nested_groupedCV_candidate_model_benchmark")


def plot_confusion(y, pred):
    cm = confusion_matrix(y, pred, labels=CLASSES)
    recall = cm / cm.sum(axis=1, keepdims=True)
    fig, ax = plt.subplots(figsize=(1.75, 1.55))
    im = ax.imshow(recall, cmap="Blues", vmin=0, vmax=1)
    for i in range(4):
        for j in range(4):
            ax.text(j, i, f"{cm[i,j]}\n{recall[i,j]*100:.0f}%", ha="center", va="center",
                    fontsize=4, color="white" if recall[i,j] > .55 else "black")
    ax.set_xticks(range(4), CLASSES); ax.set_yticks(range(4), CLASSES)
    ax.set_xlabel("ML OOF prediction"); ax.set_ylabel("HC-derived class")
    ax.tick_params(length=0, pad=1)
    for s in ax.spines.values(): s.set_visible(False)
    cb = fig.colorbar(im, ax=ax, fraction=.05, pad=.03, ticks=[0, .5, 1]); cb.set_label("Row fraction")
    save(fig, "02_nested_groupedCV_aggregated_OOF_confusion")
    pd.DataFrame(cm, index=CLASSES, columns=CLASSES).to_csv(OUT / "02_OOF_confusion_counts.csv")


def plot_stability(data):
    fig, axes = plt.subplots(1, 2, figsize=(3.5, 1.65), sharey=True)
    rng = np.random.default_rng(SEED)
    for ax, column, ylabel in [(axes[0], "ML_vote_agreement", "Repeated-fold agreement"),
                               (axes[1], "ML_mean_max_probability", "Mean predicted probability")]:
        for i, m in enumerate(CLASSES):
            vals = data.loc[data.M_class.eq(m), column].to_numpy()
            ax.scatter(i + rng.normal(0, .055, len(vals)), vals, s=5, c=COLORS[m], edgecolors="none", alpha=.75)
            ax.hlines(np.median(vals), i-.18, i+.18, color="black", lw=.8)
        ax.set_xticks(range(4), CLASSES); ax.set_ylim(.45, 1.015); ax.set_ylabel(ylabel)
        ax.spines[["top", "right"]].set_visible(False)
    save(fig, "03_per_cell_repeated_OOF_prediction_stability")


def tune_family_and_importance(x, y, groups, specs, winner, outer_splits):
    family_specs = [s for s in specs if s["family"] == winner]
    rows, selected_names = [], []
    for repeat, fold, tr, te in outer_splits:
        scores = inner_scores(family_specs, x, y, groups, tr, SEED + 4000 + repeat * 10 + fold)
        best = int(np.argmax(scores)); selected_names.append(family_specs[best]["name"])
        ztr, zte = preprocess_train_test(x[tr], x[te])
        model = clone(family_specs[best]["model"]); model.fit(ztr, y[tr])
        pi = permutation_importance(model, zte, y[te], scoring="balanced_accuracy", n_repeats=20,
                                    random_state=SEED + 5000 + repeat * 10 + fold, n_jobs=-1)
        for feature, value in zip(FEATURES, pi.importances_mean):
            rows.append({"Repeat": repeat + 1, "Fold": fold, "Feature": feature,
                         "Importance_raw": value, "Selected_hyperparameters": family_specs[best]["name"]})
    imp = pd.DataFrame(rows)
    means = imp.groupby("Feature").Importance_raw.mean().clip(lower=0)
    denom = max(float(means.max()), 1e-12)
    imp["Relative_importance"] = np.clip(imp.Importance_raw / denom, 0, 1)
    order = imp.groupby("Feature").Relative_importance.mean().sort_values(ascending=False).index.tolist()
    return imp, order, Counter(selected_names).most_common(1)[0][0]


def ablation_performance(x, y, groups, feature_order, winner_spec, outer_splits):
    subsets = [
        ("All 10 features", list(range(10))),
        ("Removing top 2", [i for i, f in enumerate(FEATURES) if f not in feature_order[:2]]),
        ("Removing top 5", [i for i, f in enumerate(FEATURES) if f not in feature_order[:5]]),
        ("Retaining bottom 3", [FEATURES.index(f) for f in feature_order[-3:]]),
    ]
    rows = []
    for label, cols in subsets:
        for repeat, fold, tr, te in outer_splits:
            ztr, zte = preprocess_train_test(x[tr][:, cols], x[te][:, cols])
            model = clone(winner_spec); model.fit(ztr, y[tr]); pred = model.predict(zte)
            rows.append({"Subset": label, "Repeat": repeat+1, "Fold": fold,
                         "N_features": len(cols), **metrics(y[te], pred)})
    return pd.DataFrame(rows), subsets


def plot_importance_tsne(x, y, ids, imp, order, subsets, ablation):
    summary = imp.groupby("Feature").Relative_importance.mean().sort_values(ascending=False)
    rng = np.random.default_rng(SEED)
    fig = plt.figure(figsize=(6.5, 5.6))
    gs = fig.add_gridspec(4, 2, width_ratios=[1.18, 1.0], wspace=.42, hspace=.34)
    ax = fig.add_subplot(gs[:, 0])
    for yi, feature in enumerate(order):
        vals = imp.loc[imp.Feature.eq(feature), "Relative_importance"].to_numpy()
        ax.hlines(yi, 0, 1, color="#D9D9D9", lw=.55, zorder=0)
        ax.scatter(vals, yi + rng.normal(0, .05, len(vals)), s=6, c="#969696", edgecolors="none", alpha=.8)
        ax.scatter(summary.loc[feature], yi, s=24, c="#4C78A8", edgecolors="none", zorder=3)
    ax.set_yticks(range(len(order)), [DISPLAY[f] for f in order]); ax.invert_yaxis()
    ax.set_xlim(0, 1.015); ax.set_xticks([0, .25, .5, .75, 1]); ax.set_xlabel("Relative permutation importance")
    ax.spines[["top", "right", "left"]].set_visible(False); ax.tick_params(axis="y", length=0)
    ax.text(-.20, 1.02, "a", transform=ax.transAxes, fontsize=8, weight="bold")

    # Panel b must be identical to the frozen main-figure embedding: the final
    # NPC3 (PC1-PC3) t-SNE coordinates, not a separate t-SNE of the 10-D matrix.
    frozen = pd.read_csv(ROOT / "outputs" / "morph_qc" /
                         "final_morph187_NPC3_HCK4_GCres0p50_figures" /
                         "01_final_morph187_assignments_with_tSNE.csv").set_index(ID)
    embeddings = [frozen.loc[list(ids), ["tSNE1", "tSNE2"]].to_numpy(float)]
    # For ablations, hold the dimensionality rule constant by refitting at most
    # three PCs before t-SNE after each requested feature removal.
    for label, cols in subsets[1:]:
        z, _ = preprocess_train_test(x[:, cols], x[:, cols])
        pc = PCA(n_components=min(3, len(cols)), random_state=SEED).fit_transform(z)
        embeddings.append(TSNE(n_components=2, perplexity=25, learning_rate=50, max_iter=3000,
                               init="pca", method="barnes_hut", angle=0.0,
                               random_state=SEED).fit_transform(pc))
    # Use one common, data-driven coordinate window for every panel.  X and Y
    # limits come from the global extrema across b-e (rather than a symmetric
    # artificial range), and equal data-unit scaling preserves original shape.
    x_min = min(float(e[:, 0].min()) for e in embeddings)
    x_max = max(float(e[:, 0].max()) for e in embeddings)
    y_min = min(float(e[:, 1].min()) for e in embeddings)
    y_max = max(float(e[:, 1].max()) for e in embeddings)
    x_pad = max((x_max - x_min) * 0.05, 1e-6)
    y_pad = max((y_max - y_min) * 0.05, 1e-6)
    shared_xlim = (x_min - x_pad, x_max + x_pad)
    shared_ylim = (y_min - y_pad, y_max + y_pad)
    letters = ["b", "c", "d", "e"]
    for rowi, ((label, cols), emb, letter) in enumerate(zip(subsets, embeddings, letters)):
        a = fig.add_subplot(gs[rowi, 1])
        for m in CLASSES:
            idx = y == m
            a.scatter(emb[idx, 0], emb[idx, 1], s=5.5, c=COLORS[m], edgecolors="none")
        a.set_xlim(*shared_xlim); a.set_ylim(*shared_ylim)
        a.set_aspect("equal", adjustable="box")
        a.set_xticks([]); a.set_yticks([]); a.spines[["top", "right"]].set_visible(False)
        a.set_xlabel("t-SNE 1" if rowi == 3 else ""); a.set_ylabel("t-SNE 2")
        ba = ablation.loc[ablation.Subset.eq(label), "balanced_accuracy"].mean()
        a.text(1.04, .58, f"{label}\n{len(cols)} features\nCV BA = {ba:.2f}", transform=a.transAxes,
               va="center", ha="left", fontsize=4)
        a.text(-.14, 1.02, letter, transform=a.transAxes, fontsize=8, weight="bold")
    save(fig, "04_bestML_feature_importance_with_progressive_ablation_tSNE")


def main():
    configure(); OUT.mkdir(parents=True, exist_ok=True)
    data = pd.read_csv(SOURCE)
    x = data[FEATURES].to_numpy(float); y = data.M_class.to_numpy(str)
    groups = data.Recording_day_proxy.fillna(data[ID]).to_numpy(str)
    specs = model_specs()
    family_df, selection_df, proba, pred, vote_agreement, outer_splits = benchmark_nested(x, y, groups, specs)

    family_summary = family_df.groupby("Family").agg(
        mean_balanced_accuracy=("balanced_accuracy", "mean"),
        sd_balanced_accuracy=("balanced_accuracy", "std"),
        mean_macro_F1=("macro_F1", "mean"), mean_ARI=("ARI", "mean")
    ).sort_values("mean_balanced_accuracy", ascending=False)
    winner = family_summary.index[0]
    winner_specs = [s for s in specs if s["family"] == winner]
    imp, importance_order, most_common_param = tune_family_and_importance(x, y, groups, specs, winner, outer_splits)
    winner_spec = next(s["model"] for s in winner_specs if s["name"] == most_common_param)
    ablation, subsets = ablation_performance(x, y, groups, importance_order, winner_spec, outer_splits)

    data["ML_OOF_predicted_class"] = pred
    data["ML_OOF_correct"] = data.M_class.eq(pred)
    data["ML_mean_max_probability"] = proba.max(axis=1)
    data["ML_vote_agreement"] = vote_agreement
    for j, m in enumerate(CLASSES): data[f"ML_probability_{m}"] = proba[:, j]

    aggregate = metrics(y, pred)
    aggregate_ci = grouped_bootstrap_ci(y, pred, groups)
    precision, recall, f1, support = precision_recall_fscore_support(y, pred, labels=CLASSES, zero_division=0)
    class_metrics = pd.DataFrame({"Class": CLASSES, "Support": support, "Precision": precision,
                                  "Recall": recall, "F1": f1})
    class_metrics["Median_vote_agreement"] = [data.loc[data.M_class.eq(m), "ML_vote_agreement"].median() for m in CLASSES]
    class_metrics["Median_probability"] = [data.loc[data.M_class.eq(m), "ML_mean_max_probability"].median() for m in CLASSES]

    family_df.to_csv(OUT / "01_outer_fold_candidate_family_performance.csv", index=False)
    family_summary.to_csv(OUT / "01_candidate_family_summary.csv")
    selection_df.to_csv(OUT / "01_nested_auto_selector_fold_results.csv", index=False)
    data.to_csv(OUT / "02_per_cell_aggregated_OOF_predictions_and_confidence.csv", index=False)
    class_metrics.to_csv(OUT / "03_per_class_ML_recovery_metrics.csv", index=False)
    imp.to_csv(OUT / "04_winning_model_permutation_importance_all_outer_folds.csv", index=False)
    imp.groupby("Feature").Relative_importance.agg(["mean", "std", "median"]).sort_values("mean", ascending=False).to_csv(
        OUT / "04_winning_model_permutation_importance_summary.csv")
    ablation.to_csv(OUT / "04_feature_ablation_outer_fold_performance.csv", index=False)

    plot_benchmark(family_df); plot_confusion(y, pred); plot_stability(data)
    plot_importance_tsne(x, y, data[ID].to_numpy(str), imp, importance_order, subsets, ablation)

    result = {
        "cohort_n": len(data), "target": "HC-derived final M1-M4 labels",
        "grouping": "Recording_day_proxy (one mouse per recording day)",
        "outer_validation": f"{OUTER_REPEATS} repeats x {OUTER_SPLITS}-fold StratifiedGroupKFold",
        "inner_selection": f"{INNER_SPLITS}-fold StratifiedGroupKFold",
        "candidate_families": family_summary.index.tolist(),
        "best_family_by_mean_outer_balanced_accuracy": winner,
        "most_common_selected_hyperparameters_within_best_family": most_common_param,
        "aggregated_repeated_OOF_metrics": aggregate,
        "recording_day_group_bootstrap_95CI": aggregate_ci,
        "correct_cells": int(np.sum(pred == y)),
        "total_cells": len(y),
        "mean_vote_agreement": float(vote_agreement.mean()),
        "median_vote_agreement": float(np.median(vote_agreement)),
        "interpretation_limit": "ML recovery tests reproducibility of labels derived from the same morphology features; it is not independent biological validation.",
    }
    (OUT / "00_ML_validation_summary.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(json.dumps(result, indent=2))
    print("\nFamily summary:\n", family_summary.to_string())
    print("\nPer-class metrics:\n", class_metrics.to_string(index=False))


if __name__ == "__main__":
    main()
