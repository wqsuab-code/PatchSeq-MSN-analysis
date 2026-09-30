#!/usr/bin/env python
"""Audit E/T eligibility independently from morphology usability for all 228 Morph cells."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
STAGE1 = ROOT / ".codex-work/e_type_qc/stage1_filtered.json"
PROTOCOL = ROOT / ".codex-work/e_type_qc/protocol_qc_prepared.json"
T_LABELS = ROOT / (
    "outputs/20260821_strain_label_reaudit_v3_structured_dates/"
    "v3_strain_label_match_per_cell.csv"
)
MORPH_DISPOSITION = ROOT / (
    "outputs/morph_qc/all35_excluded_228_to_193_ASC_reason_galleries/"
    "01_all35_exclusion_reasons_and_ASC_mapping.csv"
)
OUT = ROOT / "outputs/morph_qc/all228_Morph_E_T_QC_eligibility_audit"

FEATURES13 = [
    "M_soma_minimal_radius",
    "M_soma_max_pairwise_dist",
    "M_soma_aspect_ratio",
    "M_cell_max_radial_dist",
    "M_Total_neurite_length_(sections)",
    "M_Total_neurite_volume",
    "M_Number_of_bifurcation_points",
    "M_Maximum_branch_order",
    "M_basal_dendrite_Nseg",
    "M_basal_dendrite_avg_tortuosity",
    "M_basal_dendrite_max_tortuosity",
    "M_trunk_angle_min",
    "M_trunk_angle_max",
]


def as_bool(series: pd.Series) -> pd.Series:
    return series.astype(str).str.strip().str.upper().isin(["TRUE", "T", "1"])


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    stage = json.loads(STAGE1.read_text(encoding="utf-8"))
    protocol = json.loads(PROTOCOL.read_text(encoding="utf-8"))
    raw = pd.DataFrame(stage["raw_rows"], columns=stage["raw_headers"])

    numeric_m = raw[FEATURES13].apply(pd.to_numeric, errors="coerce")
    morph = raw.loc[as_bool(raw["has_M"]) & numeric_m.notna().all(axis=1)].copy()
    morph = morph.sort_values("MSN_unique_ID").reset_index(drop=True)
    if len(morph) != 228 or morph["MSN_unique_ID"].duplicated().any():
        raise RuntimeError(f"Expected 228 unique complete Morph cells, found {len(morph)}")

    excluded = pd.DataFrame(stage["excluded_rows"], columns=stage["excluded_headers"])
    excluded = excluded.rename(
        columns={
            "n_missing_numeric_E": "E_numeric_missing_n_at_QC",
            "missing_numeric_E": "E_missing_metrics_at_QC",
            "exclusion_reason": "E_QC_reason",
        }
    )
    excluded_cols = [
        "MSN_unique_ID",
        "E_numeric_missing_n_at_QC",
        "E_missing_metrics_at_QC",
        "E_QC_reason",
    ]

    e_numeric = [row[0] for row in protocol["rec_rows"]]
    e_numeric = [name for name in e_numeric if name in morph.columns]
    e_values = morph[e_numeric].apply(pd.to_numeric, errors="coerce")
    morph["E_numeric_missing_n_all"] = e_values.isna().sum(axis=1)
    morph["E_missing_metrics_all"] = e_values.apply(
        lambda row: "; ".join(row.index[row.isna()].tolist()), axis=1
    )
    morph["E_data_present"] = as_bool(morph["has_E"])
    morph["T_data_present"] = as_bool(morph["has_T"])
    morph = morph.merge(excluded[excluded_cols], on="MSN_unique_ID", how="left")
    morph["E_QC_status"] = "Pass_frozen_stage1_E_QC"
    morph.loc[~morph["E_data_present"], "E_QC_status"] = "Missing_E_data"
    morph.loc[morph["E_QC_reason"].notna(), "E_QC_status"] = "Fail_frozen_stage1_E_QC"

    t = pd.read_csv(T_LABELS)
    t_cols = [
        "cell_id",
        "Final_consensus_CellType",
        "Final_CellType_stability",
        "Consensus_fraction",
        "Consensus_Root",
        "Root_stability",
        "stable_class",
        "analysis_disposition",
        "match_status_v3",
    ]
    morph = morph.merge(t[t_cols], left_on="MSN_unique_ID", right_on="cell_id", how="left")
    morph["T_QC_status"] = "Stable_MSN_D1_or_D2"
    morph.loc[~morph["T_data_present"], "T_QC_status"] = "Missing_T_data"
    morph.loc[morph["cell_id"].isna() & morph["T_data_present"], "T_QC_status"] = "T_label_unresolved"
    morph.loc[morph["stable_class"].eq("Unstable"), "T_QC_status"] = "Unstable_D1_D2"
    morph.loc[morph["stable_class"].eq("Stable_IN"), "T_QC_status"] = "Stable_non_MSN_IN"

    disposition = pd.read_csv(MORPH_DISPOSITION)
    disp_cols = [
        "MSN_unique_ID",
        "Exclusion_stage",
        "Primary_reason",
        "Status",
        "Selected_ASC_filename",
        "Selected_ASC_path",
        "ASC_parse_success",
        "ASC_drawable_for_audit",
        "Mapping_status",
    ]
    morph = morph.merge(disposition[disp_cols], on="MSN_unique_ID", how="left")
    morph["Morph_current_status"] = morph["Status"].fillna("retained in current 193-cell HC/GC cohort")
    morph["Morph_ASC_status"] = morph["Mapping_status"].fillna("not re-audited in all35 gallery")

    stable_msn = morph["T_QC_status"].eq("Stable_MSN_D1_or_D2")
    e_pass = morph["E_QC_status"].eq("Pass_frozen_stage1_E_QC")
    morph["Eligible_E_T_M_integrated"] = e_pass & stable_msn
    morph["Eligible_T_M_integrated"] = stable_msn
    morph["Morph_only_interpretation"] = "Morph eligibility must be adjudicated from Morph/ASC QC, independent of E/T status"
    morph.loc[
        morph["E_QC_status"].eq("Fail_frozen_stage1_E_QC")
        & morph["Mapping_status"].str.startswith("Verified", na=False),
        "Morph_only_interpretation",
    ] = "E-QC failure only; ASC verified—eligible for Morph-only sensitivity analysis after manual morphology review"

    output_cols = [
        "MSN_unique_ID",
        "M_ID",
        "M_Neuron_id",
        "M_Total_neurite_length_(sections)",
        "M_cell_max_radial_dist",
        "E_cellName",
        "E_data_present",
        "E_QC_status",
        "E_Max.number.of.APs",
        "E_numeric_missing_n_all",
        "E_missing_metrics_all",
        "E_QC_reason",
        "T_trans_ID",
        "cellID",
        "T_data_present",
        "T_QC_status",
        "Final_consensus_CellType",
        "Final_CellType_stability",
        "Consensus_fraction",
        "Consensus_Root",
        "Root_stability",
        "stable_class",
        "match_status_v3",
        "Morph_current_status",
        "Morph_ASC_status",
        "Selected_ASC_filename",
        "Selected_ASC_path",
        "ASC_parse_success",
        "ASC_drawable_for_audit",
        "Eligible_T_M_integrated",
        "Eligible_E_T_M_integrated",
        "Morph_only_interpretation",
    ]
    audit = morph[output_cols].copy()
    audit.to_csv(OUT / "01_all228_Morph_E_T_QC_status.csv", index=False, encoding="utf-8-sig")

    e_bad = audit.loc[audit["E_QC_status"].ne("Pass_frozen_stage1_E_QC")].copy()
    e_bad.to_csv(OUT / "02_Morph_cells_E_missing_or_failed_QC.csv", index=False, encoding="utf-8-sig")
    t_bad = audit.loc[audit["T_QC_status"].ne("Stable_MSN_D1_or_D2")].copy()
    t_bad.to_csv(OUT / "03_Morph_cells_T_missing_or_not_stable_MSN.csv", index=False, encoding="utf-8-sig")
    union_bad = audit.loc[
        audit["E_QC_status"].ne("Pass_frozen_stage1_E_QC")
        | audit["T_QC_status"].ne("Stable_MSN_D1_or_D2")
    ].copy()
    union_bad.to_csv(OUT / "04_Morph_cells_not_eligible_for_full_E_T_M_integration.csv", index=False, encoding="utf-8-sig")

    focus_ids = ["J5203_Batch6", "M1142_Batch5", "O8232_Batch2"]
    audit.loc[audit["MSN_unique_ID"].isin(focus_ids)].to_csv(
        OUT / "05_user_focus_three_cells_audit.csv", index=False, encoding="utf-8-sig"
    )
    current_193_t_unstable = audit.loc[
        audit["Morph_current_status"].eq("retained in current 193-cell HC/GC cohort")
        & audit["T_QC_status"].ne("Stable_MSN_D1_or_D2")
    ].copy()
    current_193_t_unstable.to_csv(
        OUT / "06_current193_cells_with_unstable_D1_D2.csv", index=False, encoding="utf-8-sig"
    )

    summary = pd.DataFrame(
        [
            ("Complete Morph cohort", len(audit)),
            ("E data present", int(audit["E_data_present"].sum())),
            ("E data missing", int((~audit["E_data_present"]).sum())),
            ("Explicit frozen E-QC failures", int(audit["E_QC_status"].eq("Fail_frozen_stage1_E_QC").sum())),
            ("T data present", int(audit["T_data_present"].sum())),
            ("T data missing", int((~audit["T_data_present"]).sum())),
            ("Stable D1", int(audit["stable_class"].eq("Stable_D1").sum())),
            ("Stable D2", int(audit["stable_class"].eq("Stable_D2").sum())),
            ("Unstable D1/D2", int(audit["T_QC_status"].eq("Unstable_D1_D2").sum())),
            ("Stable non-MSN IN", int(audit["T_QC_status"].eq("Stable_non_MSN_IN").sum())),
            ("Eligible T+M", int(audit["Eligible_T_M_integrated"].sum())),
            ("Eligible E+T+M", int(audit["Eligible_E_T_M_integrated"].sum())),
            ("E/T union requiring exclusion from full E+T+M", len(union_bad)),
            ("Current 193 with unstable D1/D2", len(current_193_t_unstable)),
        ],
        columns=["Metric", "N"],
    )
    summary.to_csv(OUT / "00_summary.csv", index=False, encoding="utf-8-sig")
    print(summary.to_string(index=False))
    print("\nE failures:")
    print(e_bad[["MSN_unique_ID", "E_Max.number.of.APs", "E_QC_reason", "stable_class"]].to_string(index=False))
    print("\nT not stable MSN:")
    print(t_bad[["MSN_unique_ID", "T_QC_status", "Final_CellType_stability", "Root_stability"]].to_string(index=False))


if __name__ == "__main__":
    main()
