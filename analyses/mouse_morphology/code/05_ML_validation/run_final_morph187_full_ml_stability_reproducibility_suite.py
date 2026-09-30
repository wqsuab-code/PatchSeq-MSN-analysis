#!/usr/bin/env python
"""Full internal ML stability/reproducibility suite for final M1-M4 taxonomy.

The target labels are HC-derived from the same ten morphology measurements.
Accordingly, this suite quantifies recoverability, partition robustness, animal-
proxy/day portability, feature dependence and calibration. It is not external
biological validation of the taxonomy.
"""

from __future__ import annotations

import json
from itertools import combinations
from pathlib import Path

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import numpy as np
import pandas as pd
from scipy.stats import binomtest
from sklearn.base import clone
from sklearn.metrics import (
    accuracy_score, adjusted_rand_score, average_precision_score,
    balanced_accuracy_score, confusion_matrix, f1_score, log_loss,
    precision_recall_fscore_support, roc_auc_score,
)
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.svm import SVC

import validate_final_morph187_taxonomy_with_nested_grouped_ml as base


ROOT = Path(__file__).resolve().parents[1]
SOURCE = base.SOURCE
BASE_ML = base.OUT
OUT = BASE_ML / "07_full_ML_stability_reproducibility_suite"
OUT.mkdir(parents=True, exist_ok=True)

ID = base.ID
FEATURES = base.FEATURES
DISPLAY = base.DISPLAY
CLASSES = base.CLASSES
COLORS = base.COLORS
SEED = 20260909
N_SEEDS = 30
N_PERMUTATIONS = 500
N_LEARNING_REPEATS = 30


def configure() -> None:
    mpl.rcParams.update({
        "font.family": "Arial", "font.size": 6, "axes.linewidth": 0.8,
        "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
        "savefig.facecolor": "white", "figure.facecolor": "white",
    })


def save(fig: plt.Figure, stem: str, dpi: int = 600) -> None:
    fig.savefig(OUT / f"{stem}.png", dpi=dpi, bbox_inches="tight", pad_inches=.03)
    fig.savefig(OUT / f"{stem}.pdf", bbox_inches="tight", pad_inches=.03)
    fig.savefig(OUT / f"{stem}.svg", bbox_inches="tight", pad_inches=.03)
    plt.close(fig)


def metric_row(y: np.ndarray, pred: np.ndarray) -> dict[str, float]:
    return {
        "Accuracy": float(accuracy_score(y, pred)),
        "Balanced_accuracy": float(balanced_accuracy_score(y, pred)),
        "Macro_F1": float(f1_score(y, pred, average="macro")),
        "ARI": float(adjusted_rand_score(y, pred)),
    }


def rbf_model(c: float = .5, gamma: str | float = "scale", probability: bool = True,
              seed: int = SEED) -> SVC:
    return SVC(C=c, gamma=gamma, kernel="rbf", class_weight="balanced",
               probability=probability, random_state=seed)


def aligned_probability(model: SVC, x: np.ndarray) -> np.ndarray:
    raw = model.predict_proba(x)
    result = np.zeros((len(x), len(CLASSES)))
    for j, label in enumerate(model.classes_):
        result[:, CLASSES.index(str(label))] = raw[:, j]
    return result


def repeated_seed_grouped_oof(x, y, groups):
    prediction_matrix = np.empty((len(y), N_SEEDS), dtype=object)
    probability_cube = np.zeros((len(y), len(CLASSES), N_SEEDS))
    rows = []
    for si in range(N_SEEDS):
        seed = SEED + si
        pred = np.empty(len(y), dtype=object)
        proba = np.zeros((len(y), len(CLASSES)))
        splitter = StratifiedGroupKFold(n_splits=5, shuffle=True, random_state=seed)
        for fold, (tr, te) in enumerate(splitter.split(x, y, groups), 1):
            ztr, zte = base.preprocess_train_test(x[tr], x[te])
            model = rbf_model(probability=True, seed=seed + fold)
            model.fit(ztr, y[tr])
            pred[te] = model.predict(zte)
            proba[te] = aligned_probability(model, zte)
        prediction_matrix[:, si] = pred
        probability_cube[:, :, si] = proba
        rows.append({"Seed_index": si + 1, "Random_seed": seed, **metric_row(y, pred)})
    return prediction_matrix, probability_cube, pd.DataFrame(rows)


def nested_leave_one_day_out(x, y, groups):
    specs = [(c, gamma) for c in [.5, 2.0, 8.0] for gamma in ["scale", .1]]
    pred = np.empty(len(y), dtype=object)
    proba = np.zeros((len(y), len(CLASSES)))
    fold_rows = []
    unique_days = np.unique(groups)
    for di, day in enumerate(unique_days):
        te = np.where(groups == day)[0]
        tr = np.where(groups != day)[0]
        inner = StratifiedGroupKFold(n_splits=4, shuffle=True, random_state=SEED + di)
        scores = []
        for c, gamma in specs:
            values = []
            for itr, iva in inner.split(x[tr], y[tr], groups[tr]):
                a, b = tr[itr], tr[iva]
                ztr, zva = base.preprocess_train_test(x[a], x[b])
                model = rbf_model(c, gamma, probability=False, seed=SEED + di)
                model.fit(ztr, y[a])
                values.append(balanced_accuracy_score(y[b], model.predict(zva)))
            scores.append(float(np.mean(values)))
        best = int(np.argmax(scores)); c, gamma = specs[best]
        ztr, zte = base.preprocess_train_test(x[tr], x[te])
        model = rbf_model(c, gamma, probability=True, seed=SEED + di)
        model.fit(ztr, y[tr])
        pred[te] = model.predict(zte)
        proba[te] = aligned_probability(model, zte)
        fold_rows.append({
            "Recording_day": day, "N_test": len(te), "C": c, "Gamma": gamma,
            "Inner_BA": scores[best], "Day_accuracy": accuracy_score(y[te], pred[te]),
        })
    return pred, proba, pd.DataFrame(fold_rows)


def permutation_null(x, y, groups):
    splitter = list(StratifiedGroupKFold(n_splits=5, shuffle=True, random_state=SEED).split(x, y, groups))

    def evaluate(labels: np.ndarray) -> float:
        pred = np.empty(len(labels), dtype=object)
        for fold, (tr, te) in enumerate(splitter):
            ztr, zte = base.preprocess_train_test(x[tr], x[te])
            model = rbf_model(probability=False, seed=SEED + fold)
            model.fit(ztr, labels[tr])
            pred[te] = model.predict(zte)
        return float(balanced_accuracy_score(labels, pred))

    observed = evaluate(y)
    rng = np.random.default_rng(SEED)
    null = np.zeros(N_PERMUTATIONS)
    for i in range(N_PERMUTATIONS):
        null[i] = evaluate(rng.permutation(y))
    p = float((1 + np.sum(null >= observed)) / (N_PERMUTATIONS + 1))
    return observed, null, p


def subset_robustness(x, y, groups):
    subsets = [("All 10", list(range(10)), "Baseline")]
    subsets += [(f"Without {DISPLAY[f]}", [j for j in range(10) if j != i], "Leave-one-feature-out")
                for i, f in enumerate(FEATURES)]
    subsets += [(f"Only {DISPLAY[f]}", [i], "Single feature") for i, f in enumerate(FEATURES)]
    rows = []
    for repeat in range(5):
        splitter = StratifiedGroupKFold(n_splits=5, shuffle=True, random_state=SEED + repeat)
        for fold, (tr, te) in enumerate(splitter.split(x, y, groups), 1):
            for name, cols, kind in subsets:
                ztr, zte = base.preprocess_train_test(x[tr][:, cols], x[te][:, cols])
                model = rbf_model(probability=False, seed=SEED + repeat * 10 + fold)
                model.fit(ztr, y[tr]); p = model.predict(zte)
                rows.append({"Subset": name, "Analysis": kind, "N_features": len(cols),
                             "Repeat": repeat + 1, "Fold": fold, **metric_row(y[te], p)})
    return pd.DataFrame(rows)


def valid_group_split(y, groups, rng, test_fraction=.2):
    unique = np.unique(groups); n_test = max(1, round(len(unique) * test_fraction))
    for _ in range(2000):
        test_groups = rng.choice(unique, n_test, replace=False)
        te = np.isin(groups, test_groups); tr = ~te
        if len(np.unique(y[te])) == 4 and len(np.unique(y[tr])) == 4:
            return np.where(tr)[0], np.where(te)[0]
    raise RuntimeError("Unable to obtain a class-complete group split")


def learning_curve(x, y, groups):
    rng = np.random.default_rng(SEED)
    fractions = [.25, .40, .60, .80, 1.00]
    rows = []
    for repeat in range(N_LEARNING_REPEATS):
        pool, te = valid_group_split(y, groups, rng)
        pool_groups = np.unique(groups[pool])
        for fraction in fractions:
            target = max(4, round(len(pool_groups) * fraction))
            if fraction == 1:
                tr = pool
            else:
                tr = None
                for _ in range(1000):
                    chosen = rng.choice(pool_groups, target, replace=False)
                    candidate = pool[np.isin(groups[pool], chosen)]
                    if len(np.unique(y[candidate])) == 4:
                        tr = candidate; break
                if tr is None:
                    continue
            ztr, zte = base.preprocess_train_test(x[tr], x[te])
            model = rbf_model(probability=False, seed=SEED + repeat)
            model.fit(ztr, y[tr]); p = model.predict(zte)
            rows.append({"Repeat": repeat + 1, "Training_group_fraction": fraction,
                         "N_training_days": len(np.unique(groups[tr])), "N_training_cells": len(tr),
                         **metric_row(y[te], p)})
    return pd.DataFrame(rows)


def calibration_and_discrimination(y, proba):
    indices = np.array([CLASSES.index(v) for v in y])
    onehot = np.eye(len(CLASSES))[indices]
    pred_idx = np.argmax(proba, axis=1)
    confidence = np.max(proba, axis=1)
    correct = pred_idx == indices
    bins = np.linspace(0, 1, 11)
    bin_rows = []
    ece = 0.0
    for i in range(10):
        sel = (confidence >= bins[i]) & (confidence < bins[i+1] if i < 9 else confidence <= bins[i+1])
        if not np.any(sel):
            continue
        acc = float(correct[sel].mean()); conf = float(confidence[sel].mean())
        ece += sel.mean() * abs(acc - conf)
        bin_rows.append({"Bin_lower": bins[i], "Bin_upper": bins[i+1], "N": int(sel.sum()),
                         "Mean_confidence": conf, "Observed_accuracy": acc})
    class_rows = []
    for j, label in enumerate(CLASSES):
        truth = onehot[:, j]
        class_rows.append({"Class": label, "ROC_AUC_OVR": roc_auc_score(truth, proba[:, j]),
                           "Average_precision": average_precision_score(truth, proba[:, j])})
    summary = {
        "Multiclass_log_loss": float(log_loss(y, proba, labels=CLASSES)),
        "Multiclass_Brier_score": float(np.mean(np.sum((proba - onehot) ** 2, axis=1))),
        "Expected_calibration_error_10bins": float(ece),
        "Macro_OVR_ROC_AUC": float(roc_auc_score(onehot, proba, average="macro", multi_class="ovr")),
        "Macro_average_precision": float(np.mean([r["Average_precision"] for r in class_rows])),
    }
    return summary, pd.DataFrame(bin_rows), pd.DataFrame(class_rows)


def group_bootstrap_metrics(y, pred, groups, iterations=5000):
    rng = np.random.default_rng(SEED); unique = np.unique(groups); rows = []
    for _ in range(iterations):
        chosen = rng.choice(unique, len(unique), replace=True)
        idx = np.concatenate([np.where(groups == g)[0] for g in chosen])
        if len(np.unique(y[idx])) < 4:
            continue
        rows.append(metric_row(y[idx], pred[idx]))
    frame = pd.DataFrame(rows)
    return {c: [float(frame[c].quantile(.025)), float(frame[c].quantile(.975))] for c in frame.columns}


def dashboard(data, family_summary, nested_pred, nested_proba, lodo_pred, seed_metrics,
              seed_agreement, observed_perm, null, perm_p, calibration_bins, subset_df,
              learning_df):
    fig, axes = plt.subplots(3, 3, figsize=(10.5, 9.0))

    # a: candidate family benchmark
    ax = axes[0, 0]
    fam = family_summary.sort_values("mean_balanced_accuracy")
    ax.barh(fam.index, fam.mean_balanced_accuracy, color=["#00468B" if v == "RBF SVM" else "#AAAAAA" for v in fam.index])
    ax.set_xlim(.55, 1); ax.set_xlabel("Nested outer-fold BA"); ax.set_title("a  Model-family benchmark", loc="left", weight="bold")

    # b: nested repeated OOF confusion
    ax = axes[0, 1]; cm = confusion_matrix(data.M_class, nested_pred, labels=CLASSES)
    rec = cm / cm.sum(axis=1, keepdims=True); im = ax.imshow(rec, cmap="Blues", vmin=0, vmax=1)
    for i in range(4):
        for j in range(4):
            ax.text(j, i, f"{cm[i,j]}\n{rec[i,j]*100:.0f}%", ha="center", va="center",
                    fontsize=6, color="white" if rec[i,j] > .55 else "black")
    ax.set_xticks(range(4), CLASSES); ax.set_yticks(range(4), CLASSES)
    ax.set_xlabel("Repeated OOF prediction"); ax.set_ylabel("HC-derived class")
    ax.set_title("b  Leakage-safe recovery", loc="left", weight="bold")

    # c: nested versus leave-one-day-out class recall
    ax = axes[0, 2]; xloc = np.arange(4); width=.36
    nested_recall = np.diag(cm) / cm.sum(axis=1)
    lcm = confusion_matrix(data.M_class, lodo_pred, labels=CLASSES)
    lodo_recall = np.diag(lcm) / lcm.sum(axis=1)
    ax.bar(xloc-width/2, nested_recall, width, color="#00468B", label="Repeated grouped OOF")
    ax.bar(xloc+width/2, lodo_recall, width, color="#D4A72C", label="Nested leave-one-day-out")
    ax.set_xticks(xloc, CLASSES); ax.set_ylim(0, 1.05); ax.set_ylabel("Recall"); ax.legend(frameon=False, fontsize=5)
    ax.set_title("c  Recording-day portability", loc="left", weight="bold")

    # d: permutation null
    ax = axes[1, 0]; ax.hist(null, bins=24, color="#BDBDBD", edgecolor="white")
    ax.axvline(observed_perm, color="#ED0000", lw=1.4, label=f"Observed={observed_perm:.3f}\nP={perm_p:.4f}")
    ax.set_xlabel("Balanced accuracy"); ax.set_ylabel("Permutations"); ax.legend(frameon=False, fontsize=5)
    ax.set_title("d  Label-permutation test", loc="left", weight="bold")

    # e: calibration
    ax = axes[1, 1]; ax.plot([0,1], [0,1], "--", color="#777777", lw=.8)
    ax.plot(calibration_bins.Mean_confidence, calibration_bins.Observed_accuracy, "o-", color="#00468B", lw=1, ms=4)
    ax.set_xlim(0,1); ax.set_ylim(0,1); ax.set_xlabel("Mean confidence"); ax.set_ylabel("Observed accuracy")
    ax.set_title("e  OOF probability calibration", loc="left", weight="bold")

    # f: partition/seed stability
    ax = axes[1, 2]
    for i, label in enumerate(CLASSES):
        vals = seed_agreement[data.M_class.to_numpy() == label]
        rng = np.random.default_rng(SEED + i)
        ax.scatter(i + rng.normal(0,.055,len(vals)), vals, s=8, color=COLORS[label], alpha=.7, edgecolors="none")
        ax.hlines(np.median(vals), i-.18, i+.18, color="black", lw=1)
    ax.set_xticks(range(4), CLASSES); ax.set_ylim(0,1.02); ax.set_ylabel("Agreement across 30 partitions")
    ax.set_title("f  Cell-level reproducibility", loc="left", weight="bold")

    # g: leave-one-feature-out robustness
    ax = axes[2, 0]
    lofo = subset_df[subset_df.Analysis.eq("Leave-one-feature-out")].groupby("Subset").Balanced_accuracy.agg(["mean","std"])
    lofo = lofo.sort_values("mean")
    labels = [v.replace("Without ", "−") for v in lofo.index]
    ax.errorbar(lofo["mean"], np.arange(len(lofo)), xerr=lofo["std"], fmt="o", color="#00468B", ms=3, lw=.7)
    ax.set_yticks(np.arange(len(lofo)), labels, fontsize=4.5); ax.set_xlim(.65,1); ax.set_xlabel("Balanced accuracy")
    ax.set_title("g  Leave-one-feature-out", loc="left", weight="bold")

    # h: learning curve
    ax = axes[2, 1]
    lc = learning_df.groupby("Training_group_fraction").Balanced_accuracy.agg(["mean","std"]).reset_index()
    ax.errorbar(lc.Training_group_fraction*100, lc["mean"], yerr=lc["std"], marker="o", color="#00468B", capsize=2)
    ax.set_ylim(.45,1.02); ax.set_xlabel("Available training days (%)"); ax.set_ylabel("Held-day balanced accuracy")
    ax.set_title("h  Group-aware learning curve", loc="left", weight="bold")

    # i: confidence and errors on frozen t-SNE
    ax = axes[2, 2]; correct = data.M_class.to_numpy() == nested_pred
    ax.scatter(data.tSNE1[correct], data.tSNE2[correct], c=np.max(nested_proba[correct], axis=1), cmap="Blues",
               vmin=.25, vmax=1, s=15, edgecolors="none")
    ax.scatter(data.tSNE1[~correct], data.tSNE2[~correct], c="#ED0000", marker="x", s=23, lw=.8, label="OOF error")
    ax.set_xticks([]); ax.set_yticks([]); ax.spines[["top","right","left","bottom"]].set_visible(False)
    ax.set_title("i  Boundary-cell map", loc="left", weight="bold"); ax.legend(frameon=False, fontsize=5)

    for ax in axes.flat:
        if ax not in [axes[0,1], axes[2,2]]:
            ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout(w_pad=1.7, h_pad=1.8)
    save(fig, "01_full_ML_validation_dashboard")


def prediction_heatmap(data, prediction_matrix, seed_reference_agreement):
    true = data.M_class.to_numpy()
    order = np.concatenate([
        np.where(true == label)[0][np.argsort(seed_reference_agreement[true == label])] for label in CLASSES
    ])
    encoded = np.vectorize(CLASSES.index)(prediction_matrix[order])
    fig, ax = plt.subplots(figsize=(5.5, 4.2))
    ax.imshow(encoded, aspect="auto", interpolation="nearest",
              cmap=ListedColormap([COLORS[c] for c in CLASSES]), vmin=-.5, vmax=3.5)
    starts = np.cumsum([0] + [np.sum(true == c) for c in CLASSES])
    for s in starts[1:-1]: ax.axhline(s-.5, color="white", lw=1.2)
    ax.set_yticks([(starts[i]+starts[i+1]-1)/2 for i in range(4)], CLASSES)
    ax.set_xticks([0,9,19,29], [1,10,20,30]); ax.set_xlabel("Grouped-CV partition seed")
    ax.set_ylabel("Reference class; cells ordered by stability")
    ax.set_title("Per-cell prediction reproducibility across 30 grouped partitions")
    save(fig, "02_per_cell_30seed_prediction_stability_heatmap")


def main() -> None:
    configure()
    data = pd.read_csv(SOURCE)
    frozen_tsne = pd.read_csv(ROOT / "outputs/morph_qc/final_morph187_NPC3_HCK4_GCres0p50_figures/01_final_morph187_assignments_with_tSNE.csv")
    data = data.merge(frozen_tsne[[ID, "tSNE1", "tSNE2"]], on=ID, how="left", validate="one_to_one")
    x = data[FEATURES].to_numpy(float); y = data.M_class.to_numpy(str)
    groups = data.Recording_day_proxy.fillna(data[ID]).to_numpy(str)

    # Existing primary analysis: fully nested model-family benchmark and pure RBF repeated OOF.
    family_summary = pd.read_csv(BASE_ML / "01_candidate_family_summary.csv", index_col=0)
    nested = pd.read_csv(BASE_ML / "05_best_RBF_per_cell_aggregated_OOF_predictions.csv")
    if not np.array_equal(data[ID].to_numpy(str), nested[ID].to_numpy(str)):
        nested = data[[ID]].merge(nested, on=ID, validate="one_to_one")
    nested_pred = nested.RBF_OOF_predicted_class.to_numpy(str)
    nested_proba = nested[[f"RBF_probability_{c}" for c in CLASSES]].to_numpy(float)

    print("1/6 repeated grouped-CV partition reproducibility", flush=True)
    pred_matrix, proba_cube, seed_metrics = repeated_seed_grouped_oof(x, y, groups)
    seed_reference_agreement = np.mean(pred_matrix == y[:, None], axis=1)
    seed_modal = np.array([max(CLASSES, key=lambda c: np.sum(row == c)) for row in pred_matrix])
    seed_modal_agreement = np.array([np.mean(row == modal) for row, modal in zip(pred_matrix, seed_modal)])

    print("2/6 nested leave-one-recording-day-out", flush=True)
    lodo_pred, lodo_proba, lodo_folds = nested_leave_one_day_out(x, y, groups)

    print("3/6 label-permutation null", flush=True)
    observed_perm, null, perm_p = permutation_null(x, y, groups)

    print("4/6 feature-subset robustness", flush=True)
    subset_df = subset_robustness(x, y, groups)

    print("5/6 group-aware learning curve", flush=True)
    learning_df = learning_curve(x, y, groups)

    print("6/6 calibration, confidence and reporting", flush=True)
    calibration, calibration_bins, discrimination = calibration_and_discrimination(y, nested_proba)

    nested_metrics = metric_row(y, nested_pred)
    lodo_metrics = metric_row(y, lodo_pred)
    nested_ci = group_bootstrap_metrics(y, nested_pred, groups)
    lodo_ci = group_bootstrap_metrics(y, lodo_pred, groups)
    seed_pairwise_ari = [adjusted_rand_score(pred_matrix[:, i], pred_matrix[:, j])
                         for i, j in combinations(range(N_SEEDS), 2)]

    per_cell = data[[ID, "Recording_day_proxy", "M_class", "HC_GC_consensus", "tSNE1", "tSNE2"]].copy()
    per_cell["Nested_RBF_OOF_prediction"] = nested_pred
    per_cell["Nested_RBF_OOF_correct"] = nested_pred == y
    per_cell["Nested_RBF_OOF_confidence"] = nested_proba.max(axis=1)
    per_cell["LODO_prediction"] = lodo_pred
    per_cell["LODO_correct"] = lodo_pred == y
    per_cell["LODO_confidence"] = lodo_proba.max(axis=1)
    per_cell["Seed_modal_prediction"] = seed_modal
    per_cell["Seed_modal_agreement"] = seed_modal_agreement
    per_cell["Seed_reference_agreement"] = seed_reference_agreement
    per_cell["Reproducible_high_confidence"] = (
        per_cell.HC_GC_consensus.astype(bool) & per_cell.Nested_RBF_OOF_correct & per_cell.LODO_correct &
        (per_cell.Nested_RBF_OOF_confidence >= .70) & (per_cell.Seed_reference_agreement >= .80)
    )
    per_cell["Boundary_score"] = 1 - np.mean(np.column_stack([
        per_cell.Nested_RBF_OOF_confidence,
        per_cell.LODO_confidence,
        per_cell.Seed_reference_agreement,
        per_cell.HC_GC_consensus.astype(float),
    ]), axis=1)
    per_cell = per_cell.sort_values(["Boundary_score", "M_class"], ascending=[False, True])

    # Per-class metrics for both principal validation regimes.
    class_rows = []
    for regime, pred in [("Nested repeated grouped OOF", nested_pred), ("Nested leave-one-day-out", lodo_pred)]:
        precision, recall, f1, support = precision_recall_fscore_support(y, pred, labels=CLASSES, zero_division=0)
        for c, n, pr, re, fv in zip(CLASSES, support, precision, recall, f1):
            class_rows.append({"Regime": regime, "Class": c, "N": n, "Precision": pr, "Recall": re, "F1": fv})
    class_metrics = pd.DataFrame(class_rows)

    seed_metrics.to_csv(OUT / "01_30seed_grouped_OOF_metrics.csv", index=False)
    pd.DataFrame(pred_matrix, columns=[f"Seed_{i+1:02d}" for i in range(N_SEEDS)]).assign(**{ID: data[ID]}).to_csv(
        OUT / "02_per_cell_30seed_predictions.csv", index=False)
    lodo_folds.to_csv(OUT / "03_nested_leave_one_day_out_fold_results.csv", index=False)
    pd.DataFrame({ID: data[ID], "Reference": y, "LODO_prediction": lodo_pred,
                  "LODO_confidence": lodo_proba.max(axis=1)}).to_csv(OUT / "03_nested_LODO_per_cell_predictions.csv", index=False)
    pd.DataFrame({"Permutation": np.arange(1, N_PERMUTATIONS+1), "Null_balanced_accuracy": null}).to_csv(
        OUT / "04_label_permutation_null_500x.csv", index=False)
    subset_df.to_csv(OUT / "05_feature_subset_robustness.csv", index=False)
    learning_df.to_csv(OUT / "06_group_aware_learning_curve.csv", index=False)
    calibration_bins.to_csv(OUT / "07_probability_calibration_bins.csv", index=False)
    discrimination.to_csv(OUT / "07_classwise_ROC_AUC_and_average_precision.csv", index=False)
    class_metrics.to_csv(OUT / "08_per_class_recovery_metrics.csv", index=False)
    per_cell.to_csv(OUT / "09_per_cell_integrated_reproducibility_and_boundary_audit.csv", index=False)

    dashboard(data, family_summary, nested_pred, nested_proba, lodo_pred, seed_metrics,
              seed_reference_agreement, observed_perm, null, perm_p, calibration_bins,
              subset_df, learning_df)
    prediction_heatmap(data, pred_matrix, seed_reference_agreement)

    summary = {
        "scope": "Internal recoverability and reproducibility of HC-derived M1-M4 labels",
        "cohort_n": int(len(data)), "recording_day_groups": int(len(np.unique(groups))),
        "class_counts": data.M_class.value_counts().sort_index().astype(int).to_dict(),
        "frozen_features": FEATURES,
        "primary_nested_repeated_grouped_OOF": {"metrics": nested_metrics, "day_group_bootstrap_95CI": nested_ci},
        "nested_leave_one_recording_day_out": {"metrics": lodo_metrics, "day_group_bootstrap_95CI": lodo_ci},
        "partition_reproducibility_30_seeds": {
            "mean_balanced_accuracy": float(seed_metrics.Balanced_accuracy.mean()),
            "sd_balanced_accuracy": float(seed_metrics.Balanced_accuracy.std()),
            "median_cell_reference_agreement": float(np.median(seed_reference_agreement)),
            "mean_pairwise_ARI_between_seed_predictions": float(np.mean(seed_pairwise_ari)),
            "minimum_pairwise_ARI_between_seed_predictions": float(np.min(seed_pairwise_ari)),
        },
        "label_permutation_test": {
            "observed_balanced_accuracy": observed_perm,
            "null_mean": float(null.mean()), "null_95percentile": float(np.quantile(null, .95)),
            "empirical_P": perm_p, "permutations": N_PERMUTATIONS,
        },
        "probability_quality": calibration,
        "HC_GC_consensus_n": int(data.HC_GC_consensus.sum()),
        "integrated_reproducible_high_confidence_n": int(per_cell.Reproducible_high_confidence.sum()),
        "interpretation": (
            "Strong internal reproducibility supports recoverable M1-M4 morphology modules. "
            "Because HC labels and ML predictors use the same ten morphology features, these tests do not constitute "
            "independent biological or external-cohort validation."
        ),
    }
    (OUT / "00_full_ML_stability_reproducibility_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")

    report = f"""# Full ML stability and reproducibility audit of M1–M4

## Design

- Frozen cohort: {len(data)} cells; {len(np.unique(groups))} recording-day groups used as the animal proxy.
- Frozen predictors: 10 de-redundant morphology features.
- Primary target: HC-derived M1–M4 labels (M1={sum(y=='M1')}, M2={sum(y=='M2')}, M3={sum(y=='M3')}, M4={sum(y=='M4')}).
- Leakage control: every transformation is fitted on the training fold only; recording days never cross train/test folds.

## Results

- Nested repeated grouped OOF: accuracy {nested_metrics['Accuracy']:.3f}, balanced accuracy {nested_metrics['Balanced_accuracy']:.3f}, macro-F1 {nested_metrics['Macro_F1']:.3f}, ARI {nested_metrics['ARI']:.3f}.
- Nested leave-one-recording-day-out: accuracy {lodo_metrics['Accuracy']:.3f}, balanced accuracy {lodo_metrics['Balanced_accuracy']:.3f}, macro-F1 {lodo_metrics['Macro_F1']:.3f}, ARI {lodo_metrics['ARI']:.3f}.
- Thirty partition seeds: mean balanced accuracy {seed_metrics.Balanced_accuracy.mean():.3f} ± {seed_metrics.Balanced_accuracy.std():.3f}; mean pairwise prediction ARI {np.mean(seed_pairwise_ari):.3f}.
- Label permutation: observed balanced accuracy {observed_perm:.3f}; null mean {null.mean():.3f}; empirical P={perm_p:.4f} ({N_PERMUTATIONS} permutations).
- Probability quality: macro OVR ROC-AUC {calibration['Macro_OVR_ROC_AUC']:.3f}, macro average precision {calibration['Macro_average_precision']:.3f}, ECE {calibration['Expected_calibration_error_10bins']:.3f}.
- Integrated high-confidence reproducible cells: {int(per_cell.Reproducible_high_confidence.sum())}/{len(data)} under the prespecified HC–GC agreement, nested-OOF correctness, LODO correctness, OOF confidence ≥0.70 and 30-seed reference agreement ≥0.80 rule.

## Interpretation

These results test whether the four HC-derived morphology labels are recoverable in unseen recording-day groups and robust to model family, partition seed, feature removal and training-set size. They are internal reproducibility evidence. They do not independently establish that M1–M4 are discrete biological cell types because the labels and predictors originate from the same morphology measurements; external reconstruction cohorts, blinded manual annotation, spatial association or orthogonal transcriptomic/electrophysiological evidence remain independent validation layers.
"""
    (OUT / "00_full_ML_validation_report.md").write_text(report, encoding="utf-8")
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
