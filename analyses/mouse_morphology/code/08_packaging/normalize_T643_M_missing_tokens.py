"""Apply modality-aware missing tokens to the final T643 metadata table.

Rules for morphology columns:
- no morphology measurement anywhere in the row -> every M field is ``NR``;
- at least one morphology measurement exists -> missing individual M fields are
  ``NaN`` and recorded measurements are preserved verbatim.

All non-M columns and row order are preserved verbatim at the parsed-table level.
"""

from pathlib import Path

import pandas as pd


INPUT = Path(r"C:\Users\53461\Downloads\Cell_metadata_T643_final.txt")
OUTDIR = Path(
    r"C:\Users\53461\OneDrive\Documentos\Patch-seq_Mouse_Acb_MSN_T-type_Visualization"
    r"\outputs\20260929_T643_final_missing_tokens"
)
OUTPUT = OUTDIR / "Cell_metadata_T643_final_NR_NaN_corrected.txt"
AUDIT = OUTDIR / "Cell_metadata_T643_final_NR_NaN_audit.csv"

MISSING_TOKENS = {"", "nr", "na", "nan", "n/a", "none", "null"}


def is_missing(value: str) -> bool:
    return str(value).strip().lower() in MISSING_TOKENS


OUTDIR.mkdir(parents=True, exist_ok=True)
df = pd.read_csv(INPUT, sep="\t", dtype=str, keep_default_na=False)
m_cols = [c for c in df.columns if c.startswith("M_")]
if not m_cols:
    raise RuntimeError("No M_ columns were found in the input table")

non_m_before = df[[c for c in df.columns if c not in m_cols]].copy()
audit_rows = []

for idx in df.index:
    original = [df.at[idx, c] for c in m_cols]
    has_m_record = any(not is_missing(v) for v in original)
    if has_m_record:
        converted = ["NaN" if is_missing(v) else v for v in original]
        status = "M_record_present_partial_missing_as_NaN"
    else:
        converted = ["NR"] * len(m_cols)
        status = "entire_M_record_missing_as_NR"
    for c, value in zip(m_cols, converted):
        df.at[idx, c] = value
    audit_rows.append(
        {
            "metadata_row": idx + 2,
            "Sample_name": df.at[idx, "Sample_name"] if "Sample_name" in df else "",
            "T_Batch": df.at[idx, "T_Batch"] if "T_Batch" in df else "",
            "M_missingness_status": status,
            "M_recorded_field_count": sum(not is_missing(v) for v in original),
            "M_NaN_field_count": sum(v == "NaN" for v in converted),
            "M_NR_field_count": sum(v == "NR" for v in converted),
        }
    )

df.to_csv(OUTPUT, sep="\t", index=False)
audit = pd.DataFrame(audit_rows)
audit.to_csv(AUDIT, index=False)

non_m_after = df[non_m_before.columns]
assert non_m_before.equals(non_m_after), "A non-M field changed unexpectedly"
assert not df[m_cols].eq("").any().any(), "Blank M values remain"

print(f"rows={len(df)}")
print(f"M_columns={len(m_cols)}")
print(audit["M_missingness_status"].value_counts().to_string())
print(f"NR_cells={int(df[m_cols].eq('NR').sum().sum())}")
print(f"NaN_cells={int(df[m_cols].eq('NaN').sum().sum())}")
print(f"output={OUTPUT}")
print(f"audit={AUDIT}")
