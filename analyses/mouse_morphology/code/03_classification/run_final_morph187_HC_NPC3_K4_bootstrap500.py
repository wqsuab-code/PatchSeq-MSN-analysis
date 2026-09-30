"""Run 500x cell-wise bootstrap stability for final HC NPC3, K=4."""

from pathlib import Path

import bootstrap_final_morph187_cellwise_stability as analysis


ROOT = Path(__file__).resolve().parents[1]

analysis.OUT = (
    ROOT
    / "outputs/morph_qc/final_morph187_NPC3_HCK4_GCres0p50_figures"
    / "31_HC_NPC3_K4_bootstrap500_stability"
)
analysis.OUT.mkdir(parents=True, exist_ok=True)
analysis.N_BOOT = 500
analysis.LABEL_COLUMN = "HC_cluster"
analysis.LEVELS = ["HC1", "HC2", "HC3", "HC4"]
analysis.COLORS = {
    "HC1": "#00468B",
    "HC2": "#42B540",
    "HC3": "#ED0000",
    "HC4": "#0099B4",
}
analysis.FILE_PREFIX = "01"


if __name__ == "__main__":
    analysis.main()
