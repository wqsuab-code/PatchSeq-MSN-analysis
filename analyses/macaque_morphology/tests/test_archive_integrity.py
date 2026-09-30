#!/usr/bin/env python3
import hashlib
import json
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
ZIP_SHA = "8b0aeaed726e27658066230fe467bdbb765d641d1ea1e3a115e079e4ff0c13a6"

def sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()

source = ROOT / "data/00_source/Macaque-PatchSeq-BG.zip"
if source.exists():
    assert sha256(source) == ZIP_SHA
else:
    source_manifest = json.loads((ROOT / "data/00_source/SOURCE_ARCHIVE_MANIFEST.json").read_text(encoding="utf-8"))
    assert source_manifest["sha256"] == ZIP_SHA
raw = pd.read_csv(ROOT / "data/02_frozen_input/01_raw_morphology_complete_126.csv")
z = pd.read_csv(ROOT / "data/03_transformed_and_PCA/02_transformed_z_all126.csv")
pca = pd.read_csv(ROOT / "data/03_transformed_and_PCA/04_pca_scores.csv")
labels = pd.read_csv(ROOT / "data/04_frozen_classification/01_temp_frozen_assignments_126.csv")
assert len(raw) == len(z) == len(pca) == len(labels) == 126
assert set(raw.cell_label) == set(z.cell_label) == set(pca.cell_label) == set(labels.cell_label)
assert len([c for c in z if c != "cell_label"]) == 18
consensus = labels.loc[labels.concordant.eq(True)].copy()
consensus["M"] = "M" + consensus.HC_K4.astype(int).astype(str)
assert len(consensus) == 117
assert consensus.M.value_counts().sort_index().to_dict() == {"M1": 43, "M2": 42, "M3": 20, "M4": 12}
rrr = pd.read_csv(ROOT / "results/03_M_T_RRR/RRR_cell_scores_n117.csv")
assert len(rrr) == 117 and set(rrr.cell_label) == set(consensus.cell_label)
params = json.loads((ROOT / "config/analysis_parameters.json").read_text(encoding="utf-8"))
assert params["feature_n"] == 18 and params["frozen_PCA_n"] == 5
print("Archive integrity checks passed: source hash, 126-cell input, 117-cell consensus, M1-M4 counts and RRR IDs.")
