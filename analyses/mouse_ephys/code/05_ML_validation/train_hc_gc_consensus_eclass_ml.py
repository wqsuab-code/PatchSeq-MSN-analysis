#!/usr/bin/env python3
"""Audit learnability of the frozen HC x GC consensus E-class partition.

The labels were derived from the same electrophysiological feature space. Thus,
this is a reproducibility/distillation analysis, not independent biological
validation. Hyperparameters are selected within each outer training fold.
"""

from __future__ import annotations

import json
import os
import re
import warnings
from pathlib import Path

import joblib
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.base import BaseEstimator, TransformerMixin, clone
from sklearn.ensemble import RandomForestClassifier
from sklearn.inspection import permutation_importance
from sklearn.metrics import (
    accuracy_score,
    balanced_accuracy_score,
    confusion_matrix,
    f1_score,
    log_loss,
)
from sklearn.model_selection import GridSearchCV, LeaveOneGroupOut, StratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC


SEED = 777
ROOT = Path(__file__).resolve().parents[1]
FEATURE_FILE = ROOT / "outputs/e_type_qc/hc5_tree_distillation/NPC3_HC5_Raw_Final18_Features.csv"
ASSIGNMENT_FILE = ROOT / "outputs/e_type_qc/HC_GC_only/HC_GC_only_cell_assignments.csv"
OUT = Path(os.environ.get(
    "MOUSE_E_SVM_BASE_OUT",
    ROOT / "outputs/e_type_qc/HC_GC_only/15_supervised_classifier_audit",
))
OUT.mkdir(parents=True, exist_ok=True)

CLASSES = ["E1", "E2", "E3", "E4", "E5"]
COLORS = {
    "E1": "#4E79A7",
    "E2": "#F28E2B",
    "E3": "#59A14F",
    "E4": "#2AA6B8",
    "E5": "#B07AA1",
}
SHORT_FEATURES = {
    "E_Holding.MP..mV.": "Holding MP",
    "E_Input.resistance..MOhm.": "Input resistance",
    "E_Membrane.time.constant..ms.": "Membrane tau",
    "E_Rheobase..pA.": "Rheobase",
    "E_Sag.ratio": "Sag ratio",
    "E_Sag.time..s.": "Sag time",
    "E_AP.threshold..mV.": "AP threshold",
    "E_AP.amplitude..mV.": "AP amplitude",
    "E_AP.width..ms.": "AP width",
    "E_Upstroke.to.downstroke.ratio": "Up/down ratio",
    "E_Afterhyperpolarization..mV.": "AHP",
    "E_Max.number.of.APs": "Max APs",
    "E_Latency..ms.": "Latency",
    "E_Latency....20pA.current..ms.": "Latency (-20 pA)",
    "E_ISI.adaptation.index": "ISI adaptation",
    "E_ISI.coefficient.of.variation": "ISI CV",
    "E_AP.amplitude.adaptation.index": "AP amp adaptation",
    "E_AP.coefficient.of.variation": "AP CV",
}


def setup_style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 6,
            "axes.linewidth": 0.5,
            "xtick.major.width": 0.5,
            "ytick.major.width": 0.5,
            "xtick.major.size": 1.5,
            "ytick.major.size": 1.5,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def save(fig: mpl.figure.Figure, stem: str, width: float, height: float) -> None:
    fig.set_size_inches(width, height)
    fig.savefig(OUT / f"{stem}.png", dpi=900, facecolor="white", bbox_inches="tight", pad_inches=0.03)
    fig.savefig(OUT / f"{stem}.pdf", facecolor="white", bbox_inches="tight", pad_inches=0.03)
    plt.close(fig)


class TranscriptomeLikeTransformer(BaseEstimator, TransformerMixin):
    """Fold-fitted shift, column-sum normalization and log1p transform."""

    def fit(self, X: np.ndarray, y: np.ndarray | None = None) -> "TranscriptomeLikeTransformer":
        arr = np.asarray(X, dtype=float)
        self.minimum_ = np.nanmin(arr, axis=0)
        shifted = np.maximum(arr - self.minimum_, 0.0)
        self.column_sum_ = shifted.sum(axis=0)
        self.column_sum_[self.column_sum_ <= 0] = 1.0
        return self

    def transform(self, X: np.ndarray) -> np.ndarray:
        arr = np.asarray(X, dtype=float)
        # Values below a training-fold minimum are mapped to the translated floor.
        shifted = np.maximum(arr - self.minimum_, 0.0)
        normalized = shifted / self.column_sum_ * 10000.0
        return np.log1p(normalized)


def get_models() -> dict[str, tuple[object, dict[str, list[object]]]]:
    return {
        "Linear SVM": (
            Pipeline(
                [
                    ("transform", TranscriptomeLikeTransformer()),
                    ("scale", StandardScaler()),
                    ("model", SVC(kernel="linear", class_weight="balanced", probability=True, random_state=SEED)),
                ]
            ),
            {"model__C": [0.03, 0.1, 0.3, 1, 3, 10]},
        ),
        "RBF-SVM": (
            Pipeline(
                [
                    ("transform", TranscriptomeLikeTransformer()),
                    ("scale", StandardScaler()),
                    ("model", SVC(kernel="rbf", class_weight="balanced", probability=True, random_state=SEED)),
                ]
            ),
            {
                "model__C": [0.3, 1, 3, 10, 30],
                "model__gamma": ["scale", 0.03, 0.1, 0.3],
            },
        ),
        "Random forest": (
            Pipeline(
                [
                    ("transform", TranscriptomeLikeTransformer()),
                    (
                        "model",
                        RandomForestClassifier(
                            n_estimators=400,
                            class_weight="balanced_subsample",
                            random_state=SEED,
                            n_jobs=1,
                        ),
                    ),
                ]
            ),
            {
                "model__max_features": ["sqrt", 0.5],
                "model__min_samples_leaf": [1, 2, 4],
                "model__max_depth": [None, 8],
            },
        ),
    }


def metrics(y_true: np.ndarray, y_pred: np.ndarray, prob: np.ndarray) -> dict[str, float]:
    return {
        "Accuracy": accuracy_score(y_true, y_pred),
        "Balanced_accuracy": balanced_accuracy_score(y_true, y_pred),
        "Macro_F1": f1_score(y_true, y_pred, average="macro"),
        "Log_loss": log_loss(y_true, prob, labels=CLASSES),
    }


def batch_from_id(cell_id: str) -> str:
    match = re.search(r"_Batch(\d+)$", str(cell_id))
    if not match:
        raise ValueError(f"Cannot parse batch from {cell_id}")
    return f"Batch{match.group(1)}"


def nested_oof(
    X: np.ndarray,
    y: np.ndarray,
    ids: np.ndarray,
    model_name: str,
    estimator: object,
    grid: dict[str, list[object]],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    outer = StratifiedKFold(n_splits=5, shuffle=True, random_state=SEED)
    pred = np.empty(len(y), dtype=object)
    prob = np.full((len(y), len(CLASSES)), np.nan)
    fold_id = np.zeros(len(y), dtype=int)
    fold_rows: list[dict[str, object]] = []

    for fold, (tr, te) in enumerate(outer.split(X, y), start=1):
        inner = StratifiedKFold(n_splits=4, shuffle=True, random_state=SEED + fold)
        search = GridSearchCV(
            clone(estimator), grid, scoring="balanced_accuracy", cv=inner,
            n_jobs=-1, refit=True, return_train_score=False,
        )
        search.fit(X[tr], y[tr])
        p = search.best_estimator_.predict(X[te])
        pr0 = search.best_estimator_.predict_proba(X[te])
        order = list(search.best_estimator_.classes_)
        pr = np.column_stack([pr0[:, order.index(c)] for c in CLASSES])
        pred[te] = p
        prob[te] = pr
        fold_id[te] = fold
        row = {"Model": model_name, "Outer_fold": fold, "N_test": len(te), **metrics(y[te], p, pr)}
        row["Best_parameters"] = json.dumps(search.best_params_, sort_keys=True)
        fold_rows.append(row)

    oof = pd.DataFrame(
        {
            "MSN_unique_ID": ids,
            "True_E": y,
            "Predicted_E": pred,
            "Correct": pred == y,
            "Outer_fold": fold_id,
            "Prediction_confidence": np.nanmax(prob, axis=1),
            **{f"P_{c}": prob[:, i] for i, c in enumerate(CLASSES)},
        }
    )
    return oof, pd.DataFrame(fold_rows)


def summarize_models(oof_by_model: dict[str, pd.DataFrame], fold_table: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for model, oof in oof_by_model.items():
        y = oof["True_E"].to_numpy()
        p = oof["Predicted_E"].to_numpy()
        prob = oof[[f"P_{c}" for c in CLASSES]].to_numpy()
        row = {"Model": model, "N": len(oof), **metrics(y, p, prob)}
        fs = fold_table[fold_table["Model"].eq(model)]
        for col in ["Accuracy", "Balanced_accuracy", "Macro_F1", "Log_loss"]:
            row[f"Fold_mean_{col}"] = fs[col].mean()
            row[f"Fold_sd_{col}"] = fs[col].std(ddof=1)
        rows.append(row)
    return pd.DataFrame(rows).sort_values(["Balanced_accuracy", "Macro_F1"], ascending=False).reset_index(drop=True)


def plot_performance(summary: pd.DataFrame, fold_table: pd.DataFrame) -> None:
    metric_cols = ["Accuracy", "Balanced_accuracy", "Macro_F1"]
    names = summary["Model"].tolist()
    fig, axes = plt.subplots(1, 3, figsize=(5.4, 1.75), sharey=True)
    color_map = {"Linear SVM": "#7F7F7F", "RBF-SVM": "#4E79A7", "Random forest": "#F28E2B"}
    for ax, metric in zip(axes, metric_cols):
        for i, name in enumerate(names):
            vals = fold_table.loc[fold_table["Model"].eq(name), metric].to_numpy()
            ax.scatter(np.full(len(vals), i), vals, s=6, color=color_map[name], alpha=0.65, linewidths=0)
            ax.errorbar(i, vals.mean(), yerr=vals.std(ddof=1), fmt="o", ms=2.8,
                        color="black", ecolor="black", elinewidth=0.5, capsize=1.5, capthick=0.5)
        ax.set_xticks(range(len(names)), names, rotation=35, ha="right")
        ax.set_ylim(0.5, 1.01)
        ax.set_title(metric.replace("_", " "), fontsize=6)
        ax.spines[["top", "right"]].set_visible(False)
    axes[0].set_ylabel("Held-out score")
    fig.tight_layout(w_pad=1.2)
    save(fig, "01_NestedCV_Model_Performance", 5.4, 1.75)


def plot_confusions(oof_by_model: dict[str, pd.DataFrame], summary: pd.DataFrame) -> None:
    names = summary["Model"].tolist()
    fig, axes = plt.subplots(1, len(names), figsize=(5.7, 2.0))
    for ax, name in zip(axes, names):
        oof = oof_by_model[name]
        cm = confusion_matrix(oof["True_E"], oof["Predicted_E"], labels=CLASSES)
        pct = cm / cm.sum(axis=1, keepdims=True)
        ax.imshow(pct, cmap="Blues", vmin=0, vmax=1)
        for i in range(5):
            for j in range(5):
                color = "white" if pct[i, j] >= 0.55 else "black"
                ax.text(j, i, f"{cm[i,j]}\n{pct[i,j]*100:.0f}%", ha="center", va="center", fontsize=4, color=color)
        row = summary.loc[summary["Model"].eq(name)].iloc[0]
        ax.set_title(f"{name}\nBA={row['Balanced_accuracy']:.3f}", fontsize=6)
        ax.set_xticks(range(5), CLASSES)
        ax.set_yticks(range(5), CLASSES)
        ax.set_xlabel("Predicted")
        if ax is axes[0]:
            ax.set_ylabel("HC×GC consensus")
        for spine in ax.spines.values():
            spine.set_linewidth(0.5)
    fig.tight_layout(w_pad=1.2)
    save(fig, "02_NestedOOF_Confusion_Matrices", 5.7, 2.0)


def fit_final_model(
    X: np.ndarray,
    y: np.ndarray,
    estimator: object,
    grid: dict[str, list[object]],
) -> GridSearchCV:
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=SEED)
    search = GridSearchCV(clone(estimator), grid, scoring="balanced_accuracy", cv=cv, n_jobs=-1, refit=True)
    search.fit(X, y)
    return search


def grouped_batch_audit(
    X: np.ndarray,
    y: np.ndarray,
    ids: np.ndarray,
    groups: np.ndarray,
    estimator: object,
    grid: dict[str, list[object]],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    logo = LeaveOneGroupOut()
    rows = []
    cell_rows = []
    for fold, (tr, te) in enumerate(logo.split(X, y, groups), start=1):
        inner = StratifiedKFold(n_splits=4, shuffle=True, random_state=SEED + 100 + fold)
        search = GridSearchCV(clone(estimator), grid, scoring="balanced_accuracy", cv=inner, n_jobs=-1, refit=True)
        search.fit(X[tr], y[tr])
        p = search.best_estimator_.predict(X[te])
        pr0 = search.best_estimator_.predict_proba(X[te])
        order = list(search.best_estimator_.classes_)
        pr = np.column_stack([pr0[:, order.index(c)] for c in CLASSES])
        batch = groups[te][0]
        rows.append(
            {
                "Held_out_batch": batch,
                "N_test": len(te),
                **metrics(y[te], p, pr),
                "Best_parameters": json.dumps(search.best_params_, sort_keys=True),
            }
        )
        for local, idx in enumerate(te):
            cell_rows.append(
                {
                    "MSN_unique_ID": ids[idx],
                    "Batch": batch,
                    "True_E": y[idx],
                    "Predicted_E": p[local],
                    "Correct": p[local] == y[idx],
                    "Prediction_confidence": pr[local].max(),
                    **{f"P_{c}": pr[local, i] for i, c in enumerate(CLASSES)},
                }
            )
    return pd.DataFrame(rows), pd.DataFrame(cell_rows)


def best_model_permutation_importance(
    X: np.ndarray,
    y: np.ndarray,
    feature_names: list[str],
    estimator: object,
    grid: dict[str, list[object]],
) -> pd.DataFrame:
    outer = StratifiedKFold(n_splits=5, shuffle=True, random_state=SEED)
    records = []
    for fold, (tr, te) in enumerate(outer.split(X, y), start=1):
        inner = StratifiedKFold(n_splits=4, shuffle=True, random_state=SEED + fold)
        search = GridSearchCV(clone(estimator), grid, scoring="balanced_accuracy", cv=inner, n_jobs=-1, refit=True)
        search.fit(X[tr], y[tr])
        imp = permutation_importance(
            search.best_estimator_, X[te], y[te], scoring="balanced_accuracy",
            n_repeats=20, random_state=SEED + fold, n_jobs=-1,
        )
        for name, mean, sd in zip(feature_names, imp.importances_mean, imp.importances_std):
            records.append({"Outer_fold": fold, "Feature": name, "Importance": mean, "Within_fold_SD": sd})
    detail = pd.DataFrame(records)
    summary = (
        detail.groupby("Feature", as_index=False)
        .agg(Mean_importance=("Importance", "mean"), Fold_SD=("Importance", "std"))
        .sort_values("Mean_importance", ascending=False)
    )
    detail.to_csv(OUT / "BestModel_OOF_Permutation_Importance_byFold.csv", index=False)
    summary.to_csv(OUT / "BestModel_OOF_Permutation_Importance_Summary.csv", index=False)
    return summary


def plot_importance(table: pd.DataFrame) -> None:
    q = table.sort_values("Mean_importance", ascending=True)
    labels = [SHORT_FEATURES.get(x, x.removeprefix("E_")) for x in q["Feature"]]
    fig, ax = plt.subplots(figsize=(2.6, 3.2))
    ax.barh(range(len(q)), q["Mean_importance"], xerr=q["Fold_SD"],
            color="#4E79A7", alpha=0.85, ecolor="black", capsize=1,
            error_kw={"elinewidth": 0.45, "capthick": 0.45})
    ax.set_yticks(range(len(q)), labels)
    ax.set_xlabel("Decrease in held-out balanced accuracy")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    save(fig, "03_BestModel_OOF_Permutation_Importance", 2.6, 3.2)


def main() -> None:
    setup_style()
    features = pd.read_csv(FEATURE_FILE)
    assignments = pd.read_csv(ASSIGNMENT_FILE)
    assignments["HC_GC_consensus"] = assignments["HC_GC_consensus"].astype(str).str.lower().eq("true")
    data = features.merge(assignments, on="MSN_unique_ID", how="inner", validate="one_to_one")
    if len(data) != 493 or data["MSN_unique_ID"].nunique() != 493:
        raise ValueError("Expected the active 493-cell cohort.")
    feature_names = [c for c in features.columns if c != "MSN_unique_ID"]
    if len(feature_names) != 18:
        raise ValueError(f"Expected 18 final E-features, found {len(feature_names)}")
    if data[feature_names].isna().any().any() or not np.isfinite(data[feature_names].to_numpy()).all():
        raise ValueError("Feature matrix contains missing or non-finite values.")

    train = data[data["HC_GC_consensus"]].copy().reset_index(drop=True)
    boundary = data[~data["HC_GC_consensus"]].copy().reset_index(drop=True)
    if len(train) != 450 or len(boundary) != 43:
        raise ValueError("Expected 450 consensus and 43 nonconsensus cells.")

    X = train[feature_names].to_numpy(dtype=float)
    y = train["HC_GC_consensus_E"].astype(str).to_numpy()
    ids = train["MSN_unique_ID"].to_numpy()
    groups = np.array([batch_from_id(x) for x in ids])
    models = get_models()

    oof_by_model: dict[str, pd.DataFrame] = {}
    fold_tables = []
    for name, (estimator, grid) in models.items():
        print(f"Nested CV: {name}", flush=True)
        oof, folds = nested_oof(X, y, ids, name, estimator, grid)
        oof_by_model[name] = oof
        fold_tables.append(folds)
        oof.to_csv(OUT / f"NestedOOF_{name.replace(' ', '_')}_Predictions.csv", index=False)

    fold_table = pd.concat(fold_tables, ignore_index=True)
    summary = summarize_models(oof_by_model, fold_table)
    best_name = summary.iloc[0]["Model"]
    best_estimator, best_grid = models[best_name]

    fold_table.to_csv(OUT / "NestedCV_Fold_Performance.csv", index=False)
    summary.to_csv(OUT / "NestedCV_Model_Performance_Summary.csv", index=False)
    plot_performance(summary, fold_table)
    plot_confusions(oof_by_model, summary)

    print(f"Leave-one-batch-out: {best_name}", flush=True)
    batch_folds, batch_cells = grouped_batch_audit(X, y, ids, groups, best_estimator, best_grid)
    batch_folds.to_csv(OUT / "BestModel_LeaveOneBatchOut_Fold_Performance.csv", index=False)
    batch_cells.to_csv(OUT / "BestModel_LeaveOneBatchOut_Cell_Predictions.csv", index=False)
    batch_summary = pd.DataFrame([metrics(
        batch_cells["True_E"].to_numpy(),
        batch_cells["Predicted_E"].to_numpy(),
        batch_cells[[f"P_{c}" for c in CLASSES]].to_numpy(),
    )])
    batch_summary.insert(0, "Model", best_name)
    batch_summary.insert(1, "N", len(batch_cells))
    batch_summary.to_csv(OUT / "BestModel_LeaveOneBatchOut_Performance_Summary.csv", index=False)

    print(f"Final fit and boundary prediction: {best_name}", flush=True)
    final_search = fit_final_model(X, y, best_estimator, best_grid)
    joblib.dump(final_search.best_estimator_, OUT / "BestModel_Fitted_on_HC_GC_Consensus450.joblib")
    (OUT / "BestModel_Parameters.json").write_text(
        json.dumps({"Model": best_name, "Best_parameters": final_search.best_params_}, indent=2, sort_keys=True),
        encoding="utf-8",
    )

    Xb = boundary[feature_names].to_numpy(dtype=float)
    bp = final_search.best_estimator_.predict(Xb)
    bp0 = final_search.best_estimator_.predict_proba(Xb)
    order = list(final_search.best_estimator_.classes_)
    bprob = np.column_stack([bp0[:, order.index(c)] for c in CLASSES])
    boundary_out = boundary[["MSN_unique_ID", "HC_E", "GC_E", "tSNE_1", "tSNE_2"]].copy()
    boundary_out["ML_predicted_E"] = bp
    boundary_out["ML_confidence"] = bprob.max(axis=1)
    boundary_out["ML_agrees_HC"] = boundary_out["ML_predicted_E"].eq(boundary_out["HC_E"])
    boundary_out["ML_agrees_GC"] = boundary_out["ML_predicted_E"].eq(boundary_out["GC_E"])
    for i, c in enumerate(CLASSES):
        boundary_out[f"P_{c}"] = bprob[:, i]
    boundary_out.sort_values("ML_confidence", ascending=False).to_csv(
        OUT / "HC_GC_Nonconsensus43_BestModel_Predictions.csv", index=False
    )

    print(f"Permutation importance: {best_name}", flush=True)
    importance = best_model_permutation_importance(X, y, feature_names, best_estimator, best_grid)
    plot_importance(importance)

    class_counts = train["HC_GC_consensus_E"].value_counts().reindex(CLASSES)
    report = [
        "HC x GC consensus E-class supervised classifier audit",
        f"Active cohort: {len(data)} cells",
        f"Training labels: {len(train)} HC x GC consensus cells",
        f"Boundary deployment set: {len(boundary)} HC-GC nonconsensus cells",
        f"Features: {len(feature_names)} raw final E-features",
        "Fold-internal preprocessing: shift by training minimum, column-sum normalization x10,000, log1p, then z-score for SVM",
        f"Class counts: {class_counts.to_dict()}",
        "Primary evaluation: 5-fold nested stratified cross-validation; 4-fold inner tuning",
        f"Selected model: {best_name}",
        f"Nested OOF balanced accuracy: {summary.iloc[0]['Balanced_accuracy']:.6f}",
        f"Nested OOF macro-F1: {summary.iloc[0]['Macro_F1']:.6f}",
        f"Nested OOF accuracy: {summary.iloc[0]['Accuracy']:.6f}",
        f"Leave-one-batch-out balanced accuracy: {batch_summary.iloc[0]['Balanced_accuracy']:.6f}",
        f"Boundary predictions agreeing with HC: {boundary_out['ML_agrees_HC'].mean():.6f}",
        f"Boundary predictions agreeing with GC: {boundary_out['ML_agrees_GC'].mean():.6f}",
        "Interpretation: label learnability/reproducibility only; labels and predictors share the same E-feature space.",
    ]
    (OUT / "run_log.txt").write_text("\n".join(report), encoding="utf-8")
    print("\n".join(report), flush=True)


if __name__ == "__main__":
    warnings.filterwarnings("ignore", category=UserWarning)
    main()
