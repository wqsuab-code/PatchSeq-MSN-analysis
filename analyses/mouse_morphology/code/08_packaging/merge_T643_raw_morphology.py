from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(r"C:\Users\53461\OneDrive\Documentos\Patch-seq_Mouse_Acb_MSN_T-type_Visualization")
T_META = Path(r"C:\Users\53461\Downloads\Cell_metadata_T643.txt")
M_SOURCE = Path(r"C:\Users\53461\Downloads\METL_20250126_with_cluster (2).xlsx")
E_MERGED = ROOT / "outputs/20260929_T643_with_raw_E/Cell_metadata_T643_with_raw_E.txt"
OUT = ROOT / "outputs/20260929_T643_with_raw_M"


def clean(value):
    if pd.isna(value):
        return ""
    return str(value).strip()


def choose_source_row(target, source):
    candidates = source[source["E_cellName"].eq(target["E_cellName"])]
    if candidates.empty:
        return None, "no_source_common_ID"

    by_batch = candidates[candidates["T_Batch"].eq(target["T_Batch"])]
    if not by_batch.empty:
        candidates = by_batch
        resolution = "common_ID_and_batch"
    else:
        resolution = "common_ID"

    target_m_id = target["M_ID"]
    target_neuron = target["M_Neuron_id"]
    if target_m_id not in {"", "NR"} and target_neuron not in {"", "NR"}:
        exact = candidates[
            candidates["M_ID"].eq(target_m_id)
            & candidates["M_Neuron_id"].eq(target_neuron)
        ]
        if not exact.empty:
            candidates = exact
            resolution = "existing_M_ID_and_neuron_ID"

    payload = candidates[M_COLS].drop_duplicates()
    if len(payload) > 1:
        return None, "ambiguous_source_M_payload"

    row = candidates.iloc[0]
    has_m = clean(row.get("has_M", "")).lower() == "true" and clean(row["M_ID"]) != ""
    if not has_m:
        return row, f"source_record_without_M__{resolution}"
    return row, f"matched_M__{resolution}"


def add_morphology(base, source):
    rows = []
    audits = []
    for idx, target in base.iterrows():
        src, status = choose_source_row(target, source)
        appended = {c: "NR" for c in M_MEASURE_COLS}
        source_row = "NR"
        source_has_m = "False"
        source_m_id = ""
        source_neuron = ""
        if src is not None:
            source_row = str(int(src["_source_row"]))
            source_has_m = clean(src.get("has_M", ""))
            source_m_id = clean(src.get("M_ID", ""))
            source_neuron = clean(src.get("M_Neuron_id", ""))
            if status.startswith("matched_M__"):
                appended = {c: clean(src[c]) or "NR" for c in M_MEASURE_COLS}

        out = target.to_dict()
        if out.get("M_ID", "") == "":
            out["M_ID"] = "NR"
        if out.get("M_Neuron_id", "") == "":
            out["M_Neuron_id"] = "NR"
        out.update(
            {
                "M_raw_match_status": status,
                "M_raw_source_workbook": M_SOURCE.name,
                "M_raw_source_sheet": "Sheet1",
                "M_raw_source_row": source_row,
                "M_raw_source_has_M": source_has_m,
            }
        )
        out.update(appended)
        rows.append(out)
        audits.append(
            {
                "metadata_row": idx + 2,
                "Sample_name": target["Sample_name"],
                "E_cellName": target["E_cellName"],
                "cellID": target["cellID"],
                "T_Batch": target["T_Batch"],
                "metadata_M_ID": target["M_ID"],
                "metadata_M_Neuron_id": target["M_Neuron_id"],
                "M_raw_match_status": status,
                "M_source_row": source_row,
                "M_source_has_M": source_has_m,
                "M_source_M_ID": source_m_id,
                "M_source_M_Neuron_id": source_neuron,
                "M_measurements_recorded": sum(v != "NR" for v in appended.values()),
            }
        )
    return pd.DataFrame(rows), pd.DataFrame(audits)


OUT.mkdir(parents=True, exist_ok=True)
meta = pd.read_csv(T_META, sep="\t", dtype=str, keep_default_na=False)
source = pd.read_excel(M_SOURCE, sheet_name="Sheet1", dtype=str).fillna("")
source["_source_row"] = np.arange(2, len(source) + 2)

M_COLS = [c for c in source.columns if c.startswith("M_")]
M_MEASURE_COLS = [c for c in M_COLS if c not in {"M_ID", "M_Neuron_id"}]

m_only, audit = add_morphology(meta, source)
m_only_path = OUT / "Cell_metadata_T643_with_raw_M.txt"
m_only.to_csv(m_only_path, sep="\t", index=False)
audit_path = OUT / "Cell_metadata_T643_M_merge_audit.csv"
audit.to_csv(audit_path, index=False)

combined_path = None
if E_MERGED.exists():
    # Preserve every existing E-field string exactly, including literal "nan"
    # tokens already present in the supplied E-enhanced metadata.
    e_base = pd.read_csv(E_MERGED, sep="\t", dtype=str, keep_default_na=False)
    combined, combined_audit = add_morphology(e_base, source)
    combined_path = OUT / "Cell_metadata_T643_with_raw_E_M.txt"
    combined.to_csv(combined_path, sep="\t", index=False)

counts = audit["M_raw_match_status"].value_counts().sort_index()
focus_ids = ["O875", "O8201", "O9273", "M377", "A2025711"]
focus = audit[audit["cellID"].str.contains("|".join(focus_ids), na=False)]

readme = [
    "T643 metadata merge with original morphology values",
    "",
    f"Input T metadata: {T_META}",
    f"Input morphology workbook: {M_SOURCE}",
    "Morphology source worksheet: Sheet1",
    "",
    "Matching hierarchy:",
    "1. Existing M_ID + M_Neuron_id when present;",
    "2. exact common ID (E_cellName) plus T_Batch;",
    "3. exact common ID only when unambiguous.",
    "",
    f"T metadata rows: {len(meta)}",
    f"Rows with verified original morphology records: {int(audit['M_raw_match_status'].str.startswith('matched_M__').sum())}",
    f"Rows without verified morphology records: {int((~audit['M_raw_match_status'].str.startswith('matched_M__')).sum())}",
    f"Raw morphology measurement columns appended: {len(M_MEASURE_COLS)}",
    "",
    "Status counts:",
]
readme.extend([f"- {k}: {v}" for k, v in counts.items()])
readme.extend(
    [
        "",
        "Important interpretation:",
        "The appended values are copied from the original METL morphology columns without",
        "transformation, normalization, imputation, QC filtering or recalculation.",
        "All missing M identifiers or measurements are represented explicitly as NR.",
        "NR means not recorded / no verified original morphology value available.",
        "A transcriptomic record is not treated as a morphology record unless the source row",
        "explicitly has has_M=True and contains an M_ID.",
        "",
        "The five restored T643 samples remain without verified M measurements:",
    ]
)
for _, row in focus.iterrows():
    readme.append(f"- {row['cellID'].split('_CKDL')[0]}_{row['T_Batch']}: {row['M_raw_match_status']}")
readme.extend(
    [
        "",
        f"M-only output: {m_only_path}",
        f"E+M combined output: {combined_path if combined_path else 'not generated'}",
        f"Row-level audit: {audit_path}",
    ]
)
(OUT / "README_M_raw_merge.txt").write_text("\n".join(readme), encoding="utf-8")

print("\n".join(readme))
