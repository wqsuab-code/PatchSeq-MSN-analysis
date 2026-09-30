"""Embed the completed Mouse M1–M4 ML validation tables in the static editor."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "outputs/morph_qc/Mouse_M4_Macaque_style_ML500"
TABLES = BASE / "tables"
RAW = BASE / "raw"
TARGET = ROOT / "interactive/mouse-morph-ml-studio/dist/ml-data.js"


def records(frame: pd.DataFrame) -> list[dict]:
    clean = frame.replace({np.nan: None, np.inf: None, -np.inf: None})
    return clean.to_dict(orient="records")


def main() -> None:
    summary = json.loads((BASE / "00_results_summary.json").read_text(encoding="utf-8"))
    best = summary["best_supervised_model"]

    supervised = pd.read_csv(RAW / "14_supervised_all_models_500.csv")
    supervised_summary = pd.read_csv(TABLES / "18_supervised_summary_95CI.csv")

    cm_name = f"19_confusion_sum_{best.replace(' ', '_')}.csv"
    cm = pd.read_csv(TABLES / cm_name, index_col=0).reindex(index=["M1", "M2", "M3", "M4"], columns=["M1", "M2", "M3", "M4"], fill_value=0)
    confusion = cm.reset_index(names="true_class").rename(columns={c: f"pred_{c}" for c in cm.columns})

    sensitivity = pd.read_csv(TABLES / "08_unsupervised_NPC_K_algorithm_summary_95CI.csv")
    sensitivity = sensitivity[sensitivity.algorithm.eq("Ward")].rename(columns={"ARI_vs_frozen_M4_median": "ARI_median"})
    discovery = pd.read_csv(RAW / "07_unsupervised_all_runs_500.csv")
    discovery = discovery[(discovery.NPC.eq(3)) & (discovery.K.eq(4))].rename(columns={"ARI_vs_frozen_M4": "ARI_vs_reference"})

    full = pd.read_csv(TABLES / "10_full_data_cross_algorithm_frozen_parameter.csv")
    agreement = pd.DataFrame({"algorithm_A": "Ward", "algorithm_B": full.algorithm, "ARI": full.ARI})

    consensus = pd.read_csv(TABLES / "12_cell_level_consensus_margin.csv").rename(columns={"M_class": "reference_class"})
    roc = pd.read_csv(TABLES / "26_ROC_coordinates.csv").rename(columns={"M_class": "Class"})
    pr = pd.read_csv(TABLES / "27_PR_coordinates.csv").rename(columns={"M_class": "Class", "recall": "Recall", "precision": "Precision"})
    roc_pr_metrics = pd.read_csv(TABLES / "28_ROC_PR_metrics.csv")

    permutation = pd.read_csv(RAW / "29_within_donor_label_permutation_500.csv").rename(columns={"balanced_accuracy": "mean_balanced_accuracy"})
    permutation_summary = pd.read_csv(TABLES / "30_permutation_test_summary.csv").rename(columns={"observed_BA_median": "observed_mean_balanced_accuracy"})

    metadata = pd.read_csv(TABLES / "31_metadata_association_summary.csv").rename(columns={"variable": "metadata", "Cramers_V": "cramers_v_bias_corrected"})

    nested = pd.read_csv(TABLES / "21_nested_grouped_permutation_importance_500.csv")
    maximum = nested.importance_mean_median.abs().max() or 1
    train_importance = pd.DataFrame({
        "feature": nested.feature,
        "relative_importance": nested.importance_mean_median / maximum,
        "std": (nested.importance_mean_q975 - nested.importance_mean_q025).abs() / (3.92 * maximum),
    })

    mdi = pd.read_csv(TABLES / "20_feature_importance_MDI_logistic_500.csv").rename(columns={
        "ExtraTrees_MDI_median": "ExtraTrees_median",
        "ExtraTrees_MDI_q025": "ExtraTrees_q025",
        "ExtraTrees_MDI_q975": "ExtraTrees_q975",
        "Logistic_abs_standardized_coef_median": "Logistic_abscoef_median",
    })
    ablation = pd.read_csv(RAW / "22_progressive_ablation_500.csv")

    data = {
        "supervised": records(supervised),
        "supervised_summary": records(supervised_summary),
        "confusion": records(confusion),
        "sensitivity": records(sensitivity),
        "discovery": records(discovery),
        "agreement": records(agreement),
        "consensus": records(consensus),
        "roc": records(roc),
        "pr": records(pr),
        "roc_pr_metrics": records(roc_pr_metrics),
        "permutation": records(permutation),
        "permSummary": records(permutation_summary),
        "metadata": records(metadata),
        "trainImp": records(train_importance),
        "mdi": records(mdi),
        "ablation": records(ablation),
        "summary": summary,
    }
    TARGET.write_text("window.ML_DATA=" + json.dumps(data, ensure_ascii=False, separators=(",", ":")) + ";\n", encoding="utf-8")
    print(f"Wrote {TARGET} ({TARGET.stat().st_size:,} bytes)")


if __name__ == "__main__":
    main()
