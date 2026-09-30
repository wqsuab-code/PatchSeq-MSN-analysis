from pathlib import Path

import pandas as pd


ROOT = Path(r"C:\Users\53461\OneDrive\Documentos\Patch-seq_Mouse_Acb_MSN_T-type_Visualization")
BASE = ROOT / "outputs" / "dSTR_dSTRvSTR_E_QC"
OUT = BASE / "missingness_overlap_audit"
OUT.mkdir(parents=True, exist_ok=True)

excluded3 = ["first_isi_rheo", "first_isi_inv_rheo", "adp_v_last_rheo"]
features17 = [
    "fast_trough_v_last_rheo", "downstroke_rheo", "fast_trough_v_rheo",
    "peak_deltav_rheo", "peak_v_rheo", "postap_slope_rheo", "threshold_v_rheo",
    "trough_t_rheo", "trough_v_rheo", "upstroke_downstroke_ratio_rheo",
    "upstroke_rheo", "width_rheo", "width_rheo_ms", "width_suprathresh_rheo",
    "avg_rate_rheo", "latency_rheo", "rheobase_i",
]
primary10 = [
    "fast_trough_v_rheo", "peak_v_rheo", "postap_slope_rheo", "threshold_v_rheo",
    "upstroke_downstroke_ratio_rheo", "upstroke_rheo", "width_rheo_ms",
    "avg_rate_rheo", "latency_rheo", "rheobase_i",
]

summaries = []
feature_rows = []
pattern_rows = []
for dataset in ["dSTR", "dSTR_plus_vSTR"]:
    df = pd.read_csv(BASE / f"{dataset}_cell_missingness_QC.csv", low_memory=False)
    df["missing_excluded3_n"] = df[excluded3].isna().sum(axis=1)
    df["missing_retained17_n"] = df[features17].isna().sum(axis=1)
    df["missing_primary10_n"] = df[primary10].isna().sum(axis=1)
    df["missing_only_in_excluded3"] = (df.missing_excluded3_n > 0) & (df.missing_retained17_n == 0)
    df["missing_any_retained17"] = df.missing_retained17_n > 0
    df["missing_any_primary10"] = df.missing_primary10_n > 0
    df["missing_pattern_retained17"] = df.apply(
        lambda r: "; ".join(f for f in features17 if pd.isna(r[f])) or "Complete", axis=1
    )
    df[["dataset", "cell_label", "donor_label", "Lib_region_of_interest_label", "Subclass_name", "T_class",
        "cell_qc_status", "missing_excluded3_n", "missing_retained17_n", "missing_primary10_n",
        "missing_only_in_excluded3", "missing_any_retained17", "missing_any_primary10",
        "missing_pattern_retained17"]].to_csv(OUT / f"{dataset}_cellwise_missingness_overlap.csv", index=False)

    bad = df.cell_qc_status.eq("Exclude")
    summaries.append({
        "dataset": dataset, "n_cells": len(df),
        "cells_missing_only_excluded3": int(df.missing_only_in_excluded3.sum()),
        "cells_missing_any_retained17": int(df.missing_any_retained17.sum()),
        "cells_complete_retained17": int((df.missing_retained17_n == 0).sum()),
        "cell_qc_excluded": int(bad.sum()),
        "excluded_cells_missing_only_excluded3": int((bad & df.missing_only_in_excluded3).sum()),
        "excluded_cells_min_missing_retained17": int(df.loc[bad, "missing_retained17_n"].min()),
        "excluded_cells_max_missing_retained17": int(df.loc[bad, "missing_retained17_n"].max()),
        "excluded_cells_complete_primary10": int((bad & (df.missing_primary10_n == 0)).sum()),
        "excluded_cells_missing_primary10": int((bad & (df.missing_primary10_n > 0)).sum()),
    })
    for subset_name, mask in [("all", pd.Series(True, index=df.index)), ("cell_QC_excluded", bad), ("cell_QC_kept", ~bad)]:
        for feature in excluded3 + features17:
            feature_rows.append({"dataset": dataset, "subset": subset_name, "feature": feature,
                                 "feature_panel": "excluded3" if feature in excluded3 else ("primary10" if feature in primary10 else "retained17_extra"),
                                 "n_cells": int(mask.sum()), "n_missing": int(df.loc[mask, feature].isna().sum()),
                                 "missing_rate": float(df.loc[mask, feature].isna().mean())})
    pats = (df.loc[bad, "missing_pattern_retained17"].value_counts().rename_axis("missing_pattern_retained17").reset_index(name="n_cells"))
    pats.insert(0, "dataset", dataset)
    pattern_rows.extend(pats.to_dict("records"))

pd.DataFrame(summaries).to_csv(OUT / "missingness_overlap_summary.csv", index=False)
pd.DataFrame(feature_rows).to_csv(OUT / "feature_missingness_by_cell_QC_subset.csv", index=False)
pd.DataFrame(pattern_rows).to_csv(OUT / "excluded_cell_missingness_patterns.csv", index=False)
print(pd.DataFrame(summaries).to_string(index=False))
