from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(r"C:\Users\53461\OneDrive\Documentos\Patch-seq_Mouse_Acb_MSN_T-type_Visualization")
DATA = ROOT / ".codex_tmp" / "macaque_patchseq_bg_audit" / "Data"
OUT = ROOT / "outputs" / "dSTR_dSTRvSTR_E_QC"
OUT.mkdir(parents=True, exist_ok=True)

ROI_COL = "Lib_region_of_interest_label"
SUBCLASS_COL = "Subclass_name"
FEATURE_MISSING_THRESHOLD = 0.20
CELL_MISSING_THRESHOLD = 0.50
DATASETS = {
    "dSTR": ["Ca", "Pu", "DS"],
    "dSTR_plus_vSTR": ["Ca", "Pu", "DS", "NAC", "VS"],
}

meta = pd.read_csv(DATA / "cell_metadata_AllCell.csv", low_memory=False)
ap = pd.read_csv(DATA / "ap_features.csv", low_memory=False)
sweep = pd.read_csv(DATA / "sweep_features.csv", low_memory=False)
cross = pd.read_csv(DATA / "cross_sweep_long_square_features.csv", low_memory=False)
feature_order = json.loads(
    (ROOT / ".codex_tmp" / "ppt_missingness" / "plot_data_nac.json").read_text(encoding="utf-8")
)["feature_order"]

msn = meta.loc[
    meta[SUBCLASS_COL].fillna("").str.contains("MSN", case=False),
    ["cell_label", "donor_label", ROI_COL, SUBCLASS_COL],
].copy()
msn["T_class"] = msn[SUBCLASS_COL].map(
    {"STR D1 MSN": "D1", "STR D2 MSN": "D2", "STR Hybrid MSN": "Hybrid"}
).fillna("Other")

merged = msn.copy()
for frame in (ap, sweep, cross):
    available = [c for c in feature_order if c in frame.columns and c not in merged.columns]
    if available:
        merged = merged.merge(frame[["cell_label", *available]], on="cell_label", how="left", validate="one_to_one")

missing_features = [f for f in feature_order if f not in merged.columns]
if missing_features:
    raise RuntimeError(f"Missing rheobase columns: {missing_features}")

summary_rows = []
feature_tables = {}
cell_tables = {}
cohort_tables = {}

for dataset, rois in DATASETS.items():
    cohort = merged.loc[merged[ROI_COL].isin(rois)].copy()
    cohort["dataset"] = dataset
    feature_qc = pd.DataFrame({
        "feature_order": np.arange(1, len(feature_order) + 1),
        "feature": feature_order,
        "n_cells": len(cohort),
        "n_observed": [int(cohort[f].notna().sum()) for f in feature_order],
        "n_missing": [int(cohort[f].isna().sum()) for f in feature_order],
        "missing_rate": [float(cohort[f].isna().mean()) for f in feature_order],
    })
    feature_qc["feature_qc_status"] = np.where(
        feature_qc["missing_rate"] > FEATURE_MISSING_THRESHOLD, "Exclude", "Keep"
    )
    retained = feature_qc.loc[feature_qc["feature_qc_status"].eq("Keep"), "feature"].tolist()
    excluded = feature_qc.loc[feature_qc["feature_qc_status"].eq("Exclude"), "feature"].tolist()

    cohort["observed_n_20"] = cohort[feature_order].notna().sum(axis=1)
    cohort["missing_n_20"] = len(feature_order) - cohort["observed_n_20"]
    cohort["missing_rate_20"] = cohort["missing_n_20"] / len(feature_order)
    cohort["observed_n_retained"] = cohort[retained].notna().sum(axis=1)
    cohort["missing_n_retained"] = len(retained) - cohort["observed_n_retained"]
    cohort["missing_rate_retained"] = cohort["missing_n_retained"] / len(retained)
    cohort["cell_qc_status"] = np.where(
        cohort["missing_rate_retained"] > CELL_MISSING_THRESHOLD, "Exclude", "Keep"
    )

    sort_cols = [ROI_COL, "T_class", "donor_label", "missing_rate_retained", "cell_label"]
    cohort = cohort.sort_values(sort_cols, ascending=[True, True, True, False, True])
    feature_tables[dataset] = feature_qc
    cell_tables[dataset] = cohort[[
        "dataset", "cell_label", "donor_label", ROI_COL, SUBCLASS_COL, "T_class",
        "observed_n_20", "missing_n_20", "missing_rate_20",
        "observed_n_retained", "missing_n_retained", "missing_rate_retained", "cell_qc_status",
        *feature_order,
    ]]
    cohort_tables[dataset] = (
        cohort.groupby([ROI_COL, "T_class"], dropna=False)
        .agg(n_cells=("cell_label", "size"), n_donors=("donor_label", "nunique"),
             n_cell_qc_keep=("cell_qc_status", lambda x: int((x == "Keep").sum())),
             n_cell_qc_exclude=("cell_qc_status", lambda x: int((x == "Exclude").sum())))
        .reset_index()
    )
    summary_rows.append({
        "dataset": dataset,
        "ROI_definition": "+".join(rois),
        "n_cells": len(cohort),
        "n_donors": int(cohort["donor_label"].nunique()),
        "n_D1": int((cohort["T_class"] == "D1").sum()),
        "n_D2": int((cohort["T_class"] == "D2").sum()),
        "n_Hybrid": int((cohort["T_class"] == "Hybrid").sum()),
        "candidate_features": len(feature_order),
        "retained_features": len(retained),
        "excluded_features": len(excluded),
        "excluded_feature_names": "; ".join(excluded),
        "cell_qc_keep": int((cohort["cell_qc_status"] == "Keep").sum()),
        "cell_qc_exclude": int((cohort["cell_qc_status"] == "Exclude").sum()),
    })

    # Missingness heatmap: rows are cells; columns are the original 20 rheobase features.
    mat = cohort[feature_order].isna().astype(int).to_numpy()
    fig, ax = plt.subplots(figsize=(7.0, 3.4), dpi=300)
    ax.imshow(mat, aspect="auto", interpolation="nearest", cmap=plt.matplotlib.colors.ListedColormap(["#F2F3F5", "#D62728"]), vmin=0, vmax=1)
    ax.set_xticks(np.arange(len(feature_order)))
    ax.set_xticklabels(feature_order, rotation=55, ha="right", fontsize=5)
    ax.set_yticks([])
    ax.set_xlabel("20 rheobase features (red = missing)", fontsize=7)
    display_name = "dSTR+vSTR" if dataset == "dSTR_plus_vSTR" else "dSTR"
    ax.set_title(f"{display_name}: MSN cell-by-feature missingness (n={len(cohort)})", fontsize=9, loc="left")
    for idx in np.where(feature_qc["feature_qc_status"].eq("Exclude"))[0]:
        ax.add_patch(plt.Rectangle((idx - 0.5, -0.5), 1, len(cohort), fill=False, edgecolor="#7F0000", linewidth=0.8))
    fig.subplots_adjust(left=0.03, right=0.995, top=0.90, bottom=0.34)
    fig.savefig(OUT / f"{dataset}_rheobase20_missingness_heatmap.png", dpi=600, facecolor="white", bbox_inches="tight", pad_inches=0.04)
    plt.close(fig)

summary = pd.DataFrame(summary_rows)
summary.to_csv(OUT / "QC_summary.csv", index=False)
for dataset in DATASETS:
    feature_tables[dataset].to_csv(OUT / f"{dataset}_feature_missingness_QC.csv", index=False)
    cell_tables[dataset].to_csv(OUT / f"{dataset}_cell_missingness_QC.csv", index=False)
    cohort_tables[dataset].to_csv(OUT / f"{dataset}_ROI_Tclass_counts.csv", index=False)

payload = {
    "parameters": {
        "MSN_filter": f"{SUBCLASS_COL} contains MSN",
        "dSTR_ROIs": DATASETS["dSTR"],
        "dSTR_plus_vSTR_ROIs": DATASETS["dSTR_plus_vSTR"],
        "excluded_ROIs": ["ic"],
        "feature_missing_threshold": FEATURE_MISSING_THRESHOLD,
        "feature_rule": "Exclude when missing_rate > 0.20; exactly 0.20 is retained",
        "cell_missing_threshold": CELL_MISSING_THRESHOLD,
        "cell_rule": "Exclude only when retained-feature missing_rate > 0.50",
        "feature_order": feature_order,
    },
    "summary": summary.to_dict(orient="records"),
    "feature_qc": {k: json.loads(v.to_json(orient="records")) for k, v in feature_tables.items()},
    "cell_qc": {k: json.loads(v.to_json(orient="records")) for k, v in cell_tables.items()},
    "cohort_counts": {k: json.loads(v.to_json(orient="records")) for k, v in cohort_tables.items()},
}
(OUT / "QC_workbook_payload.json").write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
print(summary.to_string(index=False))
