"""Build the frozen-cell data file used by the interactive Macaque E radar editor."""

from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
ASSIGNMENTS = ROOT / "outputs/R3_panels_v3/00_tuned_p80_ee12_random_seed777_coordinates.csv"
RAW = ROOT / "outputs/dSTR_dSTRvSTR_E_QC/primary10_sensitivity_cohorts/dSTR_plus_vSTR_A_all_stage1QC_raw_with_NA.csv"
SCALING = ROOT / "outputs/R3_panels_v3/E_GC_res3_mergedK4_radar10_consensus368_clean_W4p64_H1_robust_scaling.csv"
DESTINATION = ROOT / "outputs/Macaque_E4_feature_explorer_site/dist/data/macaque_e_radar_cells.csv"

FEATURES = [
    "fast_trough_v_rheo",
    "peak_v_rheo",
    "postap_slope_rheo",
    "threshold_v_rheo",
    "upstroke_downstroke_ratio_rheo",
    "upstroke_rheo",
    "width_rheo_ms",
    "avg_rate_rheo",
    "latency_rheo",
    "rheobase_i",
]


def main() -> None:
    assignments = pd.read_csv(ASSIGNMENTS, dtype={"cell_label": str})
    raw = pd.read_csv(RAW, dtype={"cell_label": str})
    limits = pd.read_csv(SCALING).set_index("feature")

    cells = (
        assignments.loc[assignments["Consensus"].eq(True), ["cell_label", "HC_class"]]
        .rename(columns={"HC_class": "E_class"})
        .merge(raw[["cell_label", *FEATURES]], on="cell_label", how="inner", validate="one_to_one")
    )
    expected = {"C1": 59, "C2": 133, "C3": 57, "C4": 119}
    if len(cells) != 368 or cells["E_class"].value_counts().to_dict() != expected:
        raise RuntimeError("Frozen Macaque E consensus cohort no longer matches the expected 368 cells")
    if cells[FEATURES].isna().any().any():
        raise RuntimeError("Missing radar value in the frozen consensus cohort")

    cells["E_class"] = cells["E_class"].map({"C1": "E1", "C2": "E2", "C3": "E3", "C4": "E4"})
    for feature in FEATURES:
        values = cells[feature].to_numpy(float)
        z_score = (values - values.mean()) / values.std(ddof=0)
        cells[f"{feature}__z"] = z_score
        cells[f"{feature}__zscaled"] = np.clip((z_score + 3.0) / 6.0, 0, 1)
        lo = float(limits.loc[feature, "q2.5"])
        hi = float(limits.loc[feature, "q97.5"])
        cells[f"{feature}__scaled"] = np.clip((values - lo) / (hi - lo), 0, 1)

    DESTINATION.parent.mkdir(parents=True, exist_ok=True)
    cells.to_csv(DESTINATION, index=False, float_format="%.8g")
    print(DESTINATION)


if __name__ == "__main__":
    main()
