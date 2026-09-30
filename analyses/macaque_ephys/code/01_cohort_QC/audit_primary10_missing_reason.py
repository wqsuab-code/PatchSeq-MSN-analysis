from pathlib import Path

import pandas as pd


ROOT = Path(r"C:\Users\53461\OneDrive\Documentos\Patch-seq_Mouse_Acb_MSN_T-type_Visualization")
DATA = ROOT / ".codex_tmp" / "macaque_patchseq_bg_audit" / "Data"
BASE = ROOT / "outputs" / "dSTR_dSTRvSTR_E_QC"
OUT = BASE / "missingness_overlap_audit"

ap = pd.read_csv(DATA / "ap_features.csv", low_memory=False)
sweep = pd.read_csv(DATA / "sweep_features.csv", low_memory=False)
cross = pd.read_csv(DATA / "cross_sweep_long_square_features.csv", low_memory=False)
meta = pd.read_csv(DATA / "cell_metadata_AllCell.csv", low_memory=False)

targets = ["postap_slope_rheo", "avg_rate_rheo", "latency_rheo"]
primary10 = [
    "fast_trough_v_rheo", "peak_v_rheo", "postap_slope_rheo", "threshold_v_rheo",
    "upstroke_downstroke_ratio_rheo", "upstroke_rheo", "width_rheo_ms",
    "avg_rate_rheo", "latency_rheo", "rheobase_i",
]
sources = {"ap_features": ap, "sweep_features": sweep, "cross_sweep_long_square_features": cross}
source_rows = []
for feature in targets:
    found = [name for name, df in sources.items() if feature in df.columns]
    source_rows.append({"feature": feature, "source_tables": "; ".join(found) or "not_found"})
pd.DataFrame(source_rows).to_csv(OUT / "primary10_missing_features_source_tables.csv", index=False)

detail_all = []
for dataset in ["dSTR", "dSTR_plus_vSTR"]:
    qc = pd.read_csv(BASE / f"{dataset}_cell_missingness_QC.csv", low_memory=False)
    subset = qc.loc[qc.cell_qc_status.eq("Keep") & qc[targets].isna().any(axis=1)].copy()
    for name, df in sources.items():
        idset = set(df.cell_label.astype(str))
        subset[f"present_in_{name}"] = subset.cell_label.astype(str).isin(idset)
        subset[f"nonmissing_columns_in_{name}"] = subset.cell_label.map(
            df.set_index("cell_label").notna().sum(axis=1) if df.cell_label.is_unique else pd.Series(dtype=float)
        )
    detail_all.append(subset[["dataset", "cell_label", "donor_label", "Lib_region_of_interest_label", "Subclass_name", "T_class",
                              *targets, *[c for c in subset.columns if c.startswith("present_in_") or c.startswith("nonmissing_columns_in_")]]])
detail = pd.concat(detail_all, ignore_index=True).drop_duplicates(["dataset", "cell_label"])
detail.to_csv(OUT / "primary10_partial_cells_source_presence.csv", index=False)

# Use dSTR+vSTR as the non-duplicated union cohort for distribution summaries.
union = detail.loc[detail.dataset.eq("dSTR_plus_vSTR")].copy()
by_donor = (union.groupby(["donor_label", "Lib_region_of_interest_label", "T_class"], dropna=False)
            .size().reset_index(name="n_partial_cells").sort_values("n_partial_cells", ascending=False))
by_donor.to_csv(OUT / "primary10_partial_cells_by_donor_ROI_Tclass.csv", index=False)

# Compare source availability with complete retained cells.
qc = pd.read_csv(BASE / "dSTR_plus_vSTR_cell_missingness_QC.csv", low_memory=False)
qc = qc.loc[qc.cell_qc_status.eq("Keep")].copy()
for name, df in sources.items():
    qc[f"present_in_{name}"] = qc.cell_label.astype(str).isin(set(df.cell_label.astype(str)))
availability = []
for status, mask in [("primary10_complete", qc[targets].notna().all(axis=1)), ("primary10_partial", qc[targets].isna().any(axis=1))]:
    for name in sources:
        availability.append({"group": status, "source_table": name, "n_cells": int(mask.sum()),
                             "n_present": int(qc.loc[mask, f"present_in_{name}"].sum()),
                             "present_rate": float(qc.loc[mask, f"present_in_{name}"].mean())})
pd.DataFrame(availability).to_csv(OUT / "primary10_source_table_availability_comparison.csv", index=False)

# Export strict complete-case matrices for the harmonized primary 10 features.
complete_summary = []
for dataset in ["dSTR", "dSTR_plus_vSTR"]:
    current = pd.read_csv(BASE / f"{dataset}_cell_missingness_QC.csv", low_memory=False)
    current = current.loc[current.cell_qc_status.eq("Keep")].copy()
    complete = current.loc[current[primary10].notna().all(axis=1),
                           ["cell_label", "donor_label", "Lib_region_of_interest_label", "Subclass_name", "T_class", *primary10]].copy()
    complete.to_csv(OUT / f"{dataset}_primary10_complete_case_matrix.csv", index=False)
    complete_summary.append({"dataset": dataset, "after_stage1_cell_QC": len(current),
                             "removed_for_primary10_incompleteness": len(current) - len(complete),
                             "primary10_complete_cells": len(complete),
                             "D1": int((complete.T_class == "D1").sum()),
                             "D2": int((complete.T_class == "D2").sum()),
                             "Hybrid": int((complete.T_class == "Hybrid").sum()),
                             "donors": int(complete.donor_label.nunique())})
pd.DataFrame(complete_summary).to_csv(OUT / "primary10_complete_case_summary.csv", index=False)

print("FEATURE SOURCES")
print(pd.DataFrame(source_rows).to_string(index=False))
print("\nPARTIAL CELL COUNTS")
print(detail.groupby("dataset").size().to_string())
print("\nSOURCE AVAILABILITY")
print(pd.DataFrame(availability).to_string(index=False))
print("\nTOP DONOR/ROI/T")
print(by_donor.head(30).to_string(index=False))
