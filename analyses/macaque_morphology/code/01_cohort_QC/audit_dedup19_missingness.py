#!/usr/bin/env python3
"""Audit missingness of the 19 correlation-deduplicated Macaque M features."""

from __future__ import annotations

import hashlib
import io
import json
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parent
OUT = ROOT / "out"
ZIP_PATH = Path(r"C:\Users\53461\Downloads\Macaque-PatchSeq-BG.zip")
EXPECTED_SHA = "8b0aeaed726e27658066230fe467bdbb765d641d1ea1e3a115e079e4ff0c13a6"
ROIS = {"Ca", "Pu", "NAC"}
MSN = {"STR D1 MSN", "STR D2 MSN", "STR Hybrid MSN"}
REMOVED_AS_REDUNDANT = {"3_Sholl_PC0", "basal_dendrite_total_surface_area"}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def broad_group(value: str) -> str:
    text = str(value)
    if "Hybrid" in text or "D1D2" in text:
        return "Hybrid"
    if "D1" in text:
        return "D1"
    if "D2" in text:
        return "D2"
    return "Unresolved"


def main() -> None:
    observed_sha = sha256(ZIP_PATH)
    if observed_sha != EXPECTED_SHA:
        raise RuntimeError("Unexpected source ZIP SHA-256")
    with zipfile.ZipFile(ZIP_PATH) as zf:
        meta = pd.read_csv(io.BytesIO(zf.read("Data/cell_metadata_AllCell.csv")))
        morph = pd.read_csv(io.BytesIO(zf.read("Data/morphology_features.csv")))
    data = meta.merge(morph, on="cell_label", how="inner", validate="one_to_one")
    target = data[
        data["Lib_region_of_interest_label"].isin(ROIS)
        & data["Subclass_name"].isin(MSN)
    ].copy()
    dendrite21 = [
        c for c in morph.columns
        if c.startswith("3_Sholl_") or c.startswith("basal_dendrite_") or c.startswith("soma_")
    ]
    dedup19 = [c for c in dendrite21 if c not in REMOVED_AS_REDUNDANT]
    if len(target) != 486 or len(dendrite21) != 21 or len(dedup19) != 19:
        raise RuntimeError("Unexpected cohort or feature count")

    values = target[dedup19].apply(pd.to_numeric, errors="coerce")
    finite = np.isfinite(values.to_numpy(dtype=float))
    feature_rows = []
    for j, feature in enumerate(dedup19):
        available_n = int(finite[:, j].sum())
        feature_rows.append({
            "Feature": feature,
            "MSN_n": len(target),
            "Available_n": available_n,
            "Missing_or_nonfinite_n": len(target) - available_n,
            "Available_percent": 100 * available_n / len(target),
            "Missing_percent": 100 * (len(target) - available_n) / len(target),
        })
    feature_table = pd.DataFrame(feature_rows)

    available_per_cell = finite.sum(axis=1)
    target["Dedup19_available_n"] = available_per_cell
    target["Dedup19_missing_n"] = len(dedup19) - available_per_cell
    target["Dedup19_status"] = np.select(
        [available_per_cell == len(dedup19), available_per_cell == 0],
        ["all_19_complete", "all_19_missing"],
        default="partial_missing",
    )
    target["broad_group"] = target["Group_name"].map(broad_group)

    pattern_table = (
        target.groupby("Dedup19_status", dropna=False)
        .size().rename("Cell_n").reset_index()
    )
    pattern_table["Percent_of_486"] = 100 * pattern_table["Cell_n"] / len(target)

    roi_table = (
        target.groupby(["Lib_region_of_interest_label", "Dedup19_status"], dropna=False)
        .size().rename("Cell_n").reset_index()
    )
    roi_totals = target.groupby("Lib_region_of_interest_label").size().rename("ROI_total_n")
    roi_table = roi_table.merge(roi_totals, on="Lib_region_of_interest_label")
    roi_table["Percent_within_ROI"] = 100 * roi_table["Cell_n"] / roi_table["ROI_total_n"]

    broad_table = (
        target.groupby(["broad_group", "Dedup19_status"], dropna=False)
        .size().rename("Cell_n").reset_index()
    )
    broad_totals = target.groupby("broad_group").size().rename("Broad_group_total_n")
    broad_table = broad_table.merge(broad_totals, on="broad_group")
    broad_table["Percent_within_broad_group"] = (
        100 * broad_table["Cell_n"] / broad_table["Broad_group_total_n"]
    )

    cell_cols = [
        "cell_label", "donor_label", "Lib_region_of_interest_label",
        "Subclass_name", "Group_name", "broad_group",
        "Dedup19_available_n", "Dedup19_missing_n", "Dedup19_status",
    ]
    feature_table.to_csv(OUT / "dedup19_feature_missingness_MSN486.csv", index=False, encoding="utf-8-sig")
    target[cell_cols].sort_values("cell_label").to_csv(
        OUT / "dedup19_cell_missingness_MSN486.csv", index=False, encoding="utf-8-sig"
    )
    pattern_table.to_csv(OUT / "dedup19_missingness_patterns.csv", index=False, encoding="utf-8-sig")
    roi_table.to_csv(OUT / "dedup19_missingness_by_ROI.csv", index=False, encoding="utf-8-sig")
    broad_table.to_csv(OUT / "dedup19_missingness_by_broad_group.csv", index=False, encoding="utf-8-sig")

    status_counts = target["Dedup19_status"].value_counts().to_dict()
    if status_counts != {"all_19_missing": 360, "all_19_complete": 126}:
        raise RuntimeError(f"Unexpected 19-feature missingness pattern: {status_counts}")
    result = {
        "verified": True,
        "source_zip_sha256": observed_sha,
        "target_MSN_n": len(target),
        "deduplicated_complete_feature_n": len(dedup19),
        "removed_as_high_correlation_redundancy": sorted(REMOVED_AS_REDUNDANT),
        "all_19_complete_n": status_counts["all_19_complete"],
        "all_19_missing_n": status_counts["all_19_missing"],
        "partial_missing_n": 0,
        "per_feature_available_n": sorted(feature_table["Available_n"].unique().tolist()),
        "per_feature_missing_n": sorted(feature_table["Missing_or_nonfinite_n"].unique().tolist()),
        "mouse_morphology_results_read": False,
    }
    (OUT / "dedup19_missingness_verification.json").write_text(
        json.dumps(result, indent=2, ensure_ascii=False), encoding="utf-8"
    )

    report = f"""# Macaque MSN去冗余19项形态指标缺失审计

- 数据直接来自 `{ZIP_PATH.as_posix()}`；SHA-256 `{observed_sha}`。
- 范围为Ca、Pu、NAC的486个MSN。
- 19项定义：21项完整树突/胞体候选中去掉高相关冗余 `3_Sholl_PC0` 与 `basal_dendrite_total_surface_area`。

每一项均为126/486可用（25.93%）、360/486缺失或非有限（74.07%）。缺失模式完全成块：126个细胞的19项全部完整，360个细胞的19项全部缺失，部分缺失细胞为0。

因此，去冗余不会把可用于形态分类的细胞数从126提高到更多；暂时隔离4个PC1-PC3三倍SD细胞后，当前分类队列为122个，且这122个细胞的19项仍全部完整。
"""
    (OUT / "dedup19_missingness_report_CN.md").write_text(report, encoding="utf-8")
    print(json.dumps(result, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
