from pathlib import Path

import pandas as pd


ROOT = Path(r"C:\Users\53461\OneDrive\Documentos\Patch-seq_Mouse_Acb_MSN_T-type_Visualization")
BASE = ROOT / "outputs" / "dSTR_dSTRvSTR_E_QC"
OUT = BASE / "primary10_sensitivity_cohorts"
OUT.mkdir(parents=True, exist_ok=True)

FEATURES = [
    "fast_trough_v_rheo", "peak_v_rheo", "postap_slope_rheo", "threshold_v_rheo",
    "upstroke_downstroke_ratio_rheo", "upstroke_rheo", "width_rheo_ms",
    "avg_rate_rheo", "latency_rheo", "rheobase_i",
]
META = ["cell_label", "donor_label", "Lib_region_of_interest_label", "Subclass_name", "T_class"]

summary = []
imputation_audit = []
for dataset in ["dSTR", "dSTR_plus_vSTR"]:
    source = pd.read_csv(BASE / f"{dataset}_cell_missingness_QC.csv", low_memory=False)
    all_qc = source.loc[source.cell_qc_status.eq("Keep"), META + FEATURES].copy()
    complete = all_qc.loc[all_qc[FEATURES].notna().all(axis=1)].copy()

    # Preserve the unmodified all-cell matrix for full auditability.
    all_qc.to_csv(OUT / f"{dataset}_A_all_stage1QC_raw_with_NA.csv", index=False)
    complete.to_csv(OUT / f"{dataset}_B_primary10_complete_case.csv", index=False)

    imputed = all_qc.copy()
    missing_any = all_qc[FEATURES].isna().any(axis=1)
    imputed["any_imputation"] = missing_any
    imputed["imputed_feature_count"] = all_qc[FEATURES].isna().sum(axis=1)
    for feature in FEATURES:
        missing = all_qc[feature].isna()
        median = float(all_qc[feature].median())
        imputed[f"imputed__{feature}"] = missing
        imputed.loc[missing, feature] = median
        imputation_audit.append({
            "dataset": dataset, "feature": feature, "n_cells": len(all_qc),
            "n_missing_imputed": int(missing.sum()), "missing_rate": float(missing.mean()),
            "dataset_median_used": median,
        })
    ordered = META + FEATURES + ["any_imputation", "imputed_feature_count"] + [f"imputed__{f}" for f in FEATURES]
    imputed[ordered].to_csv(OUT / f"{dataset}_A_all_stage1QC_median_imputed.csv", index=False)

    summary.append({
        "dataset": dataset,
        "A_all_stage1QC_n": len(all_qc),
        "A_cells_with_any_imputation": int(missing_any.sum()),
        "A_cells_without_imputation": int((~missing_any).sum()),
        "B_complete_case_n": len(complete),
        "A_minus_B": len(all_qc) - len(complete),
        "imputation_method": "Feature-wise median estimated independently within dataset",
    })

pd.DataFrame(summary).to_csv(OUT / "primary10_two_cohort_summary.csv", index=False)
pd.DataFrame(imputation_audit).to_csv(OUT / "primary10_median_imputation_audit.csv", index=False)
(OUT / "README.txt").write_text(
    "Primary-10 sensitivity cohorts\n\n"
    "A: all cells passing stage-1 cell QC. Missing primary-10 values are replaced by the feature median estimated within that dataset.\n"
    "The original matrix with NA, per-cell any_imputation flag, feature count, and per-feature imputation flags are retained.\n"
    "B: strict complete cases; all ten primary features observed; no imputation.\n"
    "dSTR and dSTR+vSTR are processed independently. Future transformations must be fitted separately within each dataset and analysis fold.\n",
    encoding="utf-8",
)
print(pd.DataFrame(summary).to_string(index=False))
