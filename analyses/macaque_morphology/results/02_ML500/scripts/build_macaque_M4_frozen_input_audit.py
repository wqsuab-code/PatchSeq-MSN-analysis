#!/usr/bin/env python
"""Build the 486-cell frozen-input audit without assigning M labels to excluded cells."""
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
QC = ROOT / "outputs/Macaque_M_QC_MSN_CaPuNAc/01_MSN_CaPuNAc_cell_QC.csv"
FROZEN = ROOT / "macaque_m/m18_tempfreeze_NPC5_HCK4_res2.3/01_temp_frozen_assignments_126.csv"
E = ROOT / "outputs/R3_confusion/HC4_vs_mergedGC4_res3_cell_assignments.csv"
OUT = ROOT / "outputs/Macaque_M4_full_ML500/tables/00_frozen_input_audit_all486.csv"
FEATURES = [
    "basal_dendrite_bias_dorsal", "basal_dendrite_bias_medial",
    "basal_dendrite_calculate_number_of_stems", "basal_dendrite_extent_dorsal",
    "basal_dendrite_extent_medial", "basal_dendrite_max_branch_order",
    "basal_dendrite_max_euclidean_distance", "basal_dendrite_max_path_distance",
    "basal_dendrite_mean_contraction", "basal_dendrite_mean_diameter",
    "basal_dendrite_num_branches", "basal_dendrite_soma_percentile_dorsal",
    "basal_dendrite_soma_percentile_medial", "basal_dendrite_stem_exit_MedialLateral",
    "basal_dendrite_stem_exit_dorsal", "basal_dendrite_stem_exit_ventral",
    "basal_dendrite_total_length", "soma_surface_area",
]

qc = pd.read_csv(QC)
frozen = pd.read_csv(FROZEN)
frozen["M_class"] = "M" + frozen.HC_K4.astype(str)
e = pd.read_csv(E)
e["E_class"] = e.HC_class.astype(str).str.replace("C", "E", regex=False)
e.loc[~e.Consensus.astype(str).str.lower().eq("true"), "E_class"] = np.nan
d = qc.merge(frozen[["cell_label","M_class","HC_K4","GC_raw_K13","GC_merged_K4","concordant"]],
             on="cell_label", how="left", validate="one_to_one")
d = d.merge(e[["cell_label","E_class"]], on="cell_label", how="left", validate="one_to_one")
d["complete_frozen_18"] = d[FEATURES].notna().all(axis=1)
d["included_ML500"] = d.concordant.eq(True)
d["inclusion_reason"] = np.select(
    [~d.complete_frozen_18, d.complete_frozen_18 & ~d.concordant.eq(True), d.included_ML500],
    ["Excluded: incomplete morphology/reconstruction (one or more frozen 18 features missing)",
     "Excluded: morphology complete but HC-GC discordant",
     "Included: Macaque MSN; ROI Ca/Pu/NAC; complete frozen morphology; HC-GC consensus"],
    default="Excluded: not in frozen complete-case cohort")
d["frozen18_missing_n"] = d[FEATURES].isna().sum(axis=1)
d["frozen18_missing_rate"] = d.frozen18_missing_n / len(FEATURES)
columns = ["cell_label","Lib_donor_label","Lib_region_of_interest_label","Subclass_name","Group_name",
           "E_class","M_class","HC_K4","GC_raw_K13","GC_merged_K4","concordant","included_ML500",
           "inclusion_reason","complete_frozen_18","Complete_primary21","Observed_M_n","Missing_M_n",
           "frozen18_missing_n","frozen18_missing_rate","Any_extreme","Extreme_feature_n_abs_robust_z_ge3p5"]
d[columns].to_csv(OUT,index=False)
assert len(d)==486 and d.complete_frozen_18.sum()==126 and d.included_ML500.sum()==117
print(d.inclusion_reason.value_counts().to_string())
