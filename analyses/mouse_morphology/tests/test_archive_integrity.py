#!/usr/bin/env python3
"""Fast, dependency-free integrity checks for the frozen Mouse M archive."""

from __future__ import annotations

import csv
import hashlib
import json
import math
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CONFIG = json.loads((ROOT / "config" / "frozen_analysis.json").read_text(encoding="utf-8"))
FEATURES = CONFIG["features"]


def rows(path: Path, delimiter: str = ",") -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def as_bool(value: str) -> bool:
    return value.strip().lower() in {"true", "1", "yes"}


def test_source_manifest() -> None:
    manifest = rows(ROOT / "data" / "00_source_manifest" / "source_manifest.csv")
    assert len(manifest) == 10
    for item in manifest:
        path = ROOT / item["module_relative_path"]
        assert path.is_file(), path
        assert path.stat().st_size == int(item["bytes"])
        assert sha256(path) == item["sha256"]

    module_manifest = rows(ROOT / "MANIFEST.csv")
    manifest_paths = {item["module_relative_path"] for item in module_manifest}
    actual_paths = {
        path.relative_to(ROOT).as_posix()
        for path in ROOT.rglob("*")
        if path.is_file() and path.name != "MANIFEST.csv"
    }
    assert manifest_paths == actual_paths
    for item in module_manifest:
        path = ROOT / item["module_relative_path"]
        assert path.stat().st_size == int(item["bytes"])
        assert sha256(path) == item["sha256"]


def test_cohort_transitions_and_features() -> None:
    q228 = rows(ROOT / "data" / "02_frozen_input" / "cell_QC" / "01_all228_Morph_E_T_QC_status.csv")
    q193 = rows(ROOT / "data" / "02_frozen_input" / "morphology" / "01_blinded_193cell_morphology_input.csv")
    q187 = rows(ROOT / "data" / "02_frozen_input" / "morphology" / "01_morph187_discovery_input.csv")
    assert (len(q228), len(q193), len(q187)) == (228, 193, 187)
    assert list(q187[0])[2:] == FEATURES
    assert len({r["MSN_unique_ID"] for r in q187}) == 187
    for row in q187:
        for feature in FEATURES:
            assert row[feature] not in {"", "NA", "NaN", "NR"}
            assert math.isfinite(float(row[feature]))


def test_frozen_labels_pca_and_figure_source() -> None:
    frozen = rows(ROOT / "data" / "04_frozen_classification" / "00_final_morph187_cell_assignments.csv")
    pca = rows(ROOT / "data" / "03_transformed" / "01_PCA_scores_187cells.csv")
    heat = rows(ROOT / "data" / "04_frozen_classification" / "08_heatmap_ordered_zscore_matrix.csv")
    assert len(frozen) == len(pca) == len(heat) == 187
    ids = {r["MSN_unique_ID"] for r in frozen}
    assert ids == {r["MSN_unique_ID"] for r in pca} == {r["MSN_unique_ID"] for r in heat}
    assert Counter(r["M_class"] for r in frozen) == Counter({"M1": 67, "M2": 23, "M3": 62, "M4": 35})
    assert sum(as_bool(r["HC_GC_consensus"]) for r in frozen) == 181
    by_id = {r["MSN_unique_ID"]: r for r in frozen}
    for row in pca:
        target = by_id[row["MSN_unique_ID"]]
        for pc in ("PC1", "PC2", "PC3"):
            assert abs(float(row[pc]) - float(target[pc])) < 1e-12


def test_ml_and_rrr_cohorts() -> None:
    frozen = rows(ROOT / "data" / "04_frozen_classification" / "00_final_morph187_cell_assignments.csv")
    ids = {r["MSN_unique_ID"] for r in frozen}
    ml = rows(ROOT / "results" / "ml_nested_grouped" / "05_best_RBF_per_cell_aggregated_OOF_predictions.csv")
    assert len(ml) == 187
    assert {r["MSN_unique_ID"] for r in ml} == ids

    rrr = rows(ROOT / "results" / "t_m_rrr_strict168" / "RRR_cell_scores_strict168.csv")
    assert len(rrr) == 168
    assert Counter(r["T_identity"] for r in rrr) == Counter({"D1": 74, "D2": 94})
    assert Counter(r["M_class"] for r in rrr) == Counter({"M1": 63, "M2": 22, "M3": 55, "M4": 28})
    assert {r["MSN_unique_ID"] for r in rrr} <= {
        r["MSN_unique_ID"] for r in frozen if as_bool(r["HC_GC_consensus"])
    }

    broad = rows(ROOT / "audit" / "broad_D1D2_by_Mclass_n180.csv")
    assert sum(int(v) for row in broad for k, v in row.items() if k != "D1_D2") == 180
    strict = rows(ROOT / "audit" / "strict_RPCA_Pearson_D1D2_by_Mclass_n168.csv")
    assert sum(int(v) for row in strict for k, v in row.items() if k != "strict_T_identity") == 168

    stats = rows(ROOT / "results" / "feature_statistics" / "frozen10_global_Kruskal_Wallis_tests.csv")
    assert [r["Feature"] for r in stats] == FEATURES


def test_ml_parameters_and_module_boundary() -> None:
    params = json.loads((ROOT / "results" / "ml_500" / "00_parameters.json").read_text(encoding="utf-8"))
    assert params["features"] == FEATURES
    assert params["feature_n"] == 10
    assert params["all_complete_cells"] == 187
    assert params["recording_day_groups"] == 105
    assert params["repeats"] == 500
    assert params["E_data_used_as_predictors"] is False
    assert all(name.startswith("M_") for name in params["features"])


def test_no_unexpected_missing_values_in_key_tables() -> None:
    key_tables = [
        ROOT / "data" / "04_frozen_classification" / "00_final_morph187_cell_assignments.csv",
        ROOT / "data" / "03_transformed" / "01_PCA_scores_187cells.csv",
        ROOT / "results" / "t_m_rrr_strict168" / "RRR_cell_scores_strict168.csv",
        ROOT / "results" / "t_m_rrr_strict168" / "RRR_loadings_strict168.csv",
    ]
    for path in key_tables:
        for row in rows(path):
            assert all(value not in {"", "NA", "NaN", "NR"} for value in row.values()), path


def main() -> None:
    checks = [
        test_source_manifest,
        test_cohort_transitions_and_features,
        test_frozen_labels_pca_and_figure_source,
        test_ml_and_rrr_cohorts,
        test_ml_parameters_and_module_boundary,
        test_no_unexpected_missing_values_in_key_tables,
    ]
    for check in checks:
        check()
        print(f"PASS {check.__name__}")
    print("PASS all Mouse M archive integrity checks")


if __name__ == "__main__":
    main()
