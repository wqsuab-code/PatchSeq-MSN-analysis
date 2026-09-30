#!/usr/bin/env python3
"""Portable integrity checks for the frozen mouse E-type archive."""

from __future__ import annotations

import csv
import hashlib
import json
import math
import re
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def read_csv(relative: str) -> list[dict[str, str]]:
    with (ROOT / relative).open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def ids(rows: list[dict[str, str]], column: str = "MSN_unique_ID") -> set[str]:
    values = [row[column] for row in rows]
    assert all(values), f"Blank {column}"
    assert len(values) == len(set(values)), f"Duplicated {column}"
    return set(values)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


params = json.loads((ROOT / "config/analysis_parameters.json").read_text(encoding="utf-8"))
features = [line.strip() for line in (ROOT / "config/frozen_feature_order.txt").read_text().splitlines() if line.strip()]
assert len(features) == 18 and len(set(features)) == 18
assert all(name.startswith("E_") for name in features)
assert not any(re.search(r"morph|^M_", name, re.I) for name in features)

source_manifest = read_csv("data/00_source_manifest/source_files.csv")
assert all(re.fullmatch(r"[0-9a-f]{64}", row["sha256"]) for row in source_manifest)
for row in source_manifest:
    if row["repository_path"]:
        path = ROOT / row["repository_path"]
        assert path.is_file(), path
        assert path.stat().st_size == int(row["bytes"])
        assert sha256(path) == row["sha256"]

stage1 = json.loads((ROOT / "data/01_cohort_QC/stage1_filtered.json").read_text(encoding="utf-8"))
assert len(stage1["analysis_rows"]) == params["cohorts"]["stage1"] == 549

qc = read_csv("data/01_cohort_QC/Ephys_QCpass_494_raw_25features.csv")
active = read_csv("data/02_frozen_input/Active_Final18_PCA_input.csv")
pca = read_csv("data/03_transformed/NPC3_Frozen_PC_Scores_HC5.csv")
labels = read_csv("data/04_frozen_classification/HC_GC_only_cell_assignments.csv")
assert len(qc) == 494
assert len(active) == len(pca) == len(labels) == params["cohorts"]["active_taxonomy"] == 493
assert list(active[0].keys())[1:] == features
active_ids = ids(active)
assert ids(pca) == ids(labels) == active_ids
assert active_ids < ids(qc)
for row in active:
    for feature in features:
        assert row[feature] != ""
        assert math.isfinite(float(row[feature]))

consensus = [row for row in labels if row["HC_GC_consensus"].lower() == "true"]
assert len(consensus) == params["cohorts"]["gc_hc_consensus"] == 450
consensus_ids = ids(consensus)
assert Counter(row["HC_GC_consensus_E"] for row in consensus) == Counter(params["consensus_counts"])

strict_table = read_csv("data/04_frozen_classification/e_strict_t_identity.csv")
assert len(strict_table) == 450 and ids(strict_table) == consensus_ids
strict = [row for row in strict_table if row["major_class_agreement"].lower() == "true" and row["T_identity"] in {"D1", "D2"}]
assert len(strict) == params["cohorts"]["strict_t_identity"] == 441
strict_ids = ids(strict)
assert Counter(row["T_identity"] for row in strict) == Counter({"D1": 209, "D2": 232})
assert Counter(row["E_class"] for row in strict) == Counter({"E1": 140, "E2": 155, "E3": 45, "E4": 59, "E5": 42})

raw_stats = read_csv("results/01_statistics/Core18_Egroup_Raw_data.csv")
assert len(raw_stats) == 450 and ids(raw_stats) == consensus_ids
assert list(raw_stats[0].keys())[2:] == features
for row in raw_stats:
    assert all(row[feature] != "" and math.isfinite(float(row[feature])) for feature in features)

heatmap = read_csv("results/05_display_source/HC_GC_consensus450_core18_cell_order.csv")
radar = read_csv("results/05_display_source/E1-E5_core6_frozen_Zscore_cells.csv")
assert len(heatmap) == len(radar) == 450
assert ids(heatmap) == consensus_ids
assert ids(radar, "Cell_ID") == consensus_ids
assert len([column for column in radar[0] if column.startswith("E_") and column != "E_class"]) == 6

interactive = read_csv("interactive/data/final18_zscore_long.csv")
assert len(interactive) == 450 * 18
assert {row["Cell_ID"] for row in interactive} == consensus_ids
assert len({row["Feature_full"] for row in interactive}) == 18

confusion = read_csv("results/02_SVM_nested/01a_confusion_counts.csv")
assert sum(int(value) for row in confusion for key, value in row.items() if key) == 450
grouped_predictions = read_csv("results/03_ML_date_grouped/03_best_model_cell_predictions.csv")
assert len(grouped_predictions) == 450 and ids(grouped_predictions) == consensus_ids

rrr = read_csv("results/04_T_E_RRR/RRR_cell_scores_n441.csv")
assert len(rrr) == 441 and ids(rrr) == strict_ids
for row in rrr:
    for column in ["T_Component1", "T_Component2", "T_Component3", "E_Component1", "E_Component2", "E_Component3"]:
        assert row[column] != "" and math.isfinite(float(row[column]))

assert params["preprocessing"]["frozen_feature_n"] == 18
assert params["preprocessing"]["pca_components_for_clustering"] == 3
assert params["hierarchical_clustering"]["k"] == 5
assert params["graph_clustering"]["merged_k"] == 5
assert params["rrr"]["cohort_n"] == 441

print("PASS: source hash, 549/494/493/450/441 cohorts, 18-feature order, labels, statistics, ML, display data and RRR IDs.")
