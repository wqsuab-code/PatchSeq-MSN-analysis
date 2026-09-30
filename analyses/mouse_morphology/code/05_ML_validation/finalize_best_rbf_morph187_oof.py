"""Generate pure nested-RBF-SVM repeated OOF predictions after family benchmarking."""
from __future__ import annotations

import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from sklearn.base import clone
from sklearn.metrics import precision_recall_fscore_support
from sklearn.model_selection import StratifiedGroupKFold

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import validate_final_morph187_taxonomy_with_nested_grouped_ml as m


def main():
    m.configure()
    out = m.OUT
    data = pd.read_csv(m.SOURCE)
    x = data[m.FEATURES].to_numpy(float)
    y = data.M_class.to_numpy(str)
    groups = data.Recording_day_proxy.fillna(data[m.ID]).to_numpy(str)
    specs = [s for s in m.model_specs() if s["family"] == "RBF SVM"]

    proba_sum = np.zeros((len(y), len(m.CLASSES)))
    vote = np.zeros((len(y), len(m.CLASSES)), dtype=int)
    counts = np.zeros(len(y), dtype=int)
    folds = []
    for repeat in range(m.OUTER_REPEATS):
        outer = StratifiedGroupKFold(n_splits=m.OUTER_SPLITS, shuffle=True, random_state=m.SEED + repeat)
        for fold, (tr, te) in enumerate(outer.split(x, y, groups), 1):
            scores = m.inner_scores(specs, x, y, groups, tr, m.SEED + 4000 + repeat * 10 + fold)
            best = int(np.argmax(scores))
            ztr, zte = m.preprocess_train_test(x[tr], x[te])
            model = clone(specs[best]["model"]); model.fit(ztr, y[tr])
            pred = model.predict(zte); proba = m.aligned_proba(model, zte)
            proba_sum[te] += proba; counts[te] += 1
            for local, idx in enumerate(te): vote[idx, m.CLASSES.index(str(pred[local]))] += 1
            folds.append({"Repeat": repeat+1, "Fold": fold, "Selected_hyperparameters": specs[best]["name"],
                          "Inner_balanced_accuracy": scores[best], **m.metrics(y[te], pred)})

    proba = proba_sum / counts[:, None]
    pred = np.asarray(m.CLASSES)[np.argmax(proba, axis=1)]
    agreement = vote.max(axis=1) / counts
    result = m.metrics(y, pred)
    ci = m.grouped_bootstrap_ci(y, pred, groups)

    data["RBF_OOF_predicted_class"] = pred
    data["RBF_OOF_correct"] = data.M_class.eq(pred)
    data["RBF_mean_max_probability"] = proba.max(axis=1)
    data["RBF_vote_agreement"] = agreement
    for j, cls in enumerate(m.CLASSES): data[f"RBF_probability_{cls}"] = proba[:, j]
    data.to_csv(out / "05_best_RBF_per_cell_aggregated_OOF_predictions.csv", index=False)
    pd.DataFrame(folds).to_csv(out / "05_best_RBF_outer_fold_results.csv", index=False)

    precision, recall, f1, support = precision_recall_fscore_support(y, pred, labels=m.CLASSES, zero_division=0)
    pd.DataFrame({"Class":m.CLASSES,"Support":support,"Precision":precision,"Recall":recall,"F1":f1,
                  "Median_vote_agreement":[np.median(agreement[y==c]) for c in m.CLASSES],
                  "Median_probability":[np.median(proba[y==c].max(axis=1)) for c in m.CLASSES]}).to_csv(
        out / "05_best_RBF_per_class_metrics.csv", index=False)
    m.plot_confusion(y, pred)
    # Rename the regenerated generic confusion to make its model identity explicit.
    for ext in ["png", "pdf", "svg"]:
        src = out / f"02_nested_groupedCV_aggregated_OOF_confusion.{ext}"
        src.replace(out / f"05_best_RBF_nested_groupedCV_aggregated_OOF_confusion.{ext}")
    (out / "02_OOF_confusion_counts.csv").replace(out / "05_best_RBF_OOF_confusion_counts.csv")

    summary_path = out / "00_ML_validation_summary.json"
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    summary["best_RBF_nested_repeated_OOF_metrics"] = result
    summary["best_RBF_recording_day_group_bootstrap_95CI"] = ci
    summary["best_RBF_correct_cells"] = int(np.sum(pred == y))
    summary["best_RBF_total_cells"] = len(y)
    summary["best_RBF_mean_vote_agreement"] = float(agreement.mean())
    summary["best_RBF_median_vote_agreement"] = float(np.median(agreement))
    summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps({"metrics":result,"CI":ci,"correct":int(np.sum(pred==y)),"n":len(y)}, indent=2))


if __name__ == "__main__":
    main()
