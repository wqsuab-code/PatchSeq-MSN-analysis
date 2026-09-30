#!/usr/bin/env python3
"""Recompute all statistics used by the Macaque E1-E4 feature-comparison figure.

The script consumes the frozen 368-cell raw-value snapshot used by the
interactive figure. It produces descriptive summaries, omnibus two-sided
Kruskal-Wallis tests, tie-corrected two-sided Dunn contrasts, Holm adjustment
within each feature, and BH adjustment across all 114 pairwise tests.
"""

from __future__ import annotations

import argparse
import json
from itertools import combinations
from math import erfc, sqrt
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import kruskal, rankdata


GROUPS = ["E1", "E2", "E3", "E4"]
FEATURES = [
    ("Epsy_width_rheo", "AP width", "ms"),
    ("Epsy_fast_trough_v_rheo", "Fast trough", "mV"),
    ("Epsy_peak_deltav_rheo", "Peak delta V", "mV"),
    ("Epsy_peak_v_rheo", "Peak V", "mV"),
    ("Epsy_postap_slope_rheo", "Post-AP slope", "mV/ms"),
    ("Epsy_threshold_v_rheo", "Threshold V", "mV"),
    ("Epsy_trough_t_rheo", "Trough time", "s"),
    ("Epsy_trough_v_rheo", "Trough V", "mV"),
    ("Epsy_upstroke_downstroke_ratio_rheo", "Up/down ratio", "ratio"),
    ("Epsy_ahp_delay_5spike", "AHP delay, 5-spike", "ms"),
    ("Epsy_ahp_delay_ratio_5spike", "AHP delay ratio, 5-spike", "ratio"),
    ("Epsy_postap_slope_hero", "Post-AP slope, hero", "mV/ms"),
    ("Epsy_trough_t_hero", "Trough time, hero", "s"),
    ("Epsy_downstroke_adapt_ratio", "Downstroke adapt ratio", "ratio"),
    ("Epsy_peak_v_adapt_ratio", "Peak V adapt ratio", "ratio"),
    ("Epsy_threshold_v_adapt_ratio", "Threshold adapt ratio", "ratio"),
    ("Epsy_upstroke_adapt_ratio", "Upstroke adapt ratio", "ratio"),
    ("Epsy_width_adapt_ratio", "Width adapt ratio", "ratio"),
    ("Epsy_threshold_v_short_square", "Threshold V, short square", "mV"),
]


def holm_adjust(values: np.ndarray) -> np.ndarray:
    p = np.asarray(values, dtype=float)
    order = np.argsort(p)
    adjusted = np.empty_like(p)
    running = 0.0
    for rank, index in enumerate(order):
        running = max(running, min(1.0, (len(p) - rank) * p[index]))
        adjusted[index] = running
    return adjusted


def bh_adjust(values: np.ndarray) -> np.ndarray:
    p = np.asarray(values, dtype=float)
    order = np.argsort(p)[::-1]
    adjusted = np.empty_like(p)
    running = 1.0
    m = len(p)
    for reverse_position, index in enumerate(order, start=1):
        rank = m - reverse_position + 1
        running = min(running, min(1.0, p[index] * m / rank))
        adjusted[index] = running
    return adjusted


def dunn_table(values: np.ndarray, labels: np.ndarray) -> pd.DataFrame:
    values = np.asarray(values, dtype=float)
    labels = np.asarray(labels, dtype=str)
    ranks = rankdata(values, method="average")
    n = len(values)
    _, ties = np.unique(values, return_counts=True)
    variance = n * (n + 1) / 12 - np.sum(ties**3 - ties) / (12 * (n - 1))
    mean_ranks = {g: ranks[labels == g].mean() for g in GROUPS}
    sizes = {g: int(np.sum(labels == g)) for g in GROUPS}
    rows = []
    for g1, g2 in combinations(GROUPS, 2):
        se = sqrt(variance * (1 / sizes[g1] + 1 / sizes[g2]))
        z = (mean_ranks[g1] - mean_ranks[g2]) / se
        rows.append(
            {
                "Group_1": g1,
                "Group_2": g2,
                "N_1": sizes[g1],
                "N_2": sizes[g2],
                "Dunn_Z": z,
                "Dunn_r": abs(z) / sqrt(n),
                "P_raw": erfc(abs(z) / sqrt(2)),
            }
        )
    result = pd.DataFrame(rows)
    result["P_adj_Holm_within_feature"] = holm_adjust(result["P_raw"].to_numpy())
    return result


def load_raw(path: Path) -> pd.DataFrame:
    records = json.loads(path.read_text(encoding="utf-8"))
    frame = pd.DataFrame(records)
    expected = ["MSN_unique_ID", "E_type", *[f[0] for f in FEATURES]]
    missing = [c for c in expected if c not in frame.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")
    if len(frame) != 368 or frame["MSN_unique_ID"].nunique() != 368:
        raise ValueError("Expected 368 unique frozen consensus cells")
    if frame[[f[0] for f in FEATURES]].isna().any().any():
        raise ValueError("The frozen 19-feature input contains missing values")
    counts = frame["E_type"].value_counts().reindex(GROUPS).to_dict()
    if counts != {"E1": 59, "E2": 133, "E3": 57, "E4": 119}:
        raise ValueError(f"Frozen class counts changed: {counts}")
    return frame


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True, help="Frozen raw.json")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--reference-pairwise", type=Path)
    parser.add_argument("--reference-omnibus", type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    data = load_raw(args.input)
    labels = data["E_type"].astype(str).to_numpy()
    descriptions, omnibus_rows, pairwise_frames = [], [], []

    for feature, display, unit in FEATURES:
        for group in GROUPS:
            values = data.loc[data["E_type"].eq(group), feature].to_numpy(float)
            descriptions.append(
                {
                    "Feature": feature,
                    "Display_name": display,
                    "Unit": unit,
                    "E_type": group,
                    "N": len(values),
                    "Mean": np.mean(values),
                    "SD": np.std(values, ddof=1),
                    "Median": np.median(values),
                    "Q1": np.quantile(values, 0.25),
                    "Q3": np.quantile(values, 0.75),
                    "IQR": np.quantile(values, 0.75) - np.quantile(values, 0.25),
                    "Minimum": np.min(values),
                    "Maximum": np.max(values),
                }
            )
        arrays = [data.loc[data["E_type"].eq(g), feature].to_numpy(float) for g in GROUPS]
        statistic, p_value = kruskal(*arrays)
        omnibus_rows.append(
            {
                "Feature": feature,
                "Display_name": display,
                "Unit": unit,
                "N": len(data),
                "Kruskal_H": statistic,
                "df": 3,
                "P_raw": p_value,
                "Kruskal_epsilon_squared": max(0.0, (statistic - 4 + 1) / (len(data) - 4)),
            }
        )
        table = dunn_table(data[feature].to_numpy(float), labels)
        table.insert(0, "Feature", feature)
        table.insert(1, "Display_name", display)
        table.insert(2, "Unit", unit)
        pairwise_frames.append(table)

    descriptions = pd.DataFrame(descriptions)
    omnibus = pd.DataFrame(omnibus_rows)
    omnibus["P_adj_BH_across_19_features"] = bh_adjust(omnibus["P_raw"].to_numpy())
    omnibus["Reject_BH_0.05"] = omnibus["P_adj_BH_across_19_features"] < 0.05
    pairwise = pd.concat(pairwise_frames, ignore_index=True)
    pairwise["P_adj_BH_global_114_tests"] = bh_adjust(pairwise["P_raw"].to_numpy())
    pairwise["Reject_Holm_0.05"] = pairwise["P_adj_Holm_within_feature"] < 0.05
    pairwise["Reject_global_BH_0.05"] = pairwise["P_adj_BH_global_114_tests"] < 0.05

    descriptions.to_csv(args.output / "01_descriptive_statistics_by_E_class.csv", index=False)
    omnibus.to_csv(args.output / "02_omnibus_KruskalWallis_BH19.csv", index=False)
    pairwise.to_csv(args.output / "03_pairwise_Dunn_Holm_and_global_BH114.csv", index=False)
    data[["MSN_unique_ID", "E_type"]].to_csv(args.output / "04_frozen_cell_assignments.csv", index=False)
    pd.DataFrame(
        [{"E_type": g, "N": int(data["E_type"].eq(g).sum())} for g in GROUPS]
    ).to_csv(args.output / "05_class_counts.csv", index=False)

    audit = {
        "input_cells": len(data),
        "unique_cell_ids": int(data["MSN_unique_ID"].nunique()),
        "features": len(FEATURES),
        "class_counts": {g: int(data["E_type"].eq(g).sum()) for g in GROUPS},
        "omnibus_tests": len(omnibus),
        "pairwise_tests": len(pairwise),
        "omnibus_BH_lt_0.05": int(omnibus["Reject_BH_0.05"].sum()),
        "pairwise_Holm_lt_0.05": int(pairwise["Reject_Holm_0.05"].sum()),
        "pairwise_global_BH_lt_0.05": int(pairwise["Reject_global_BH_0.05"].sum()),
    }

    if args.reference_pairwise:
        ref = pd.DataFrame(json.loads(args.reference_pairwise.read_text(encoding="utf-8")))
        merged = pairwise.merge(
            ref,
            on=["Feature", "Group_1", "Group_2"],
            suffixes=("_new", "_reference"),
            validate="one_to_one",
        )
        fields = ["Dunn_Z", "Dunn_r", "P_raw", "P_adj_Holm_within_feature", "P_adj_BH_global_114_tests"]
        audit["pairwise_reference_max_abs_difference"] = {
            f: float(np.max(np.abs(merged[f + "_new"] - merged[f + "_reference"]))) for f in fields
        }
    if args.reference_omnibus:
        ref = pd.DataFrame(json.loads(args.reference_omnibus.read_text(encoding="utf-8")))
        merged = omnibus.merge(ref, on="Feature", suffixes=("_new", "_reference"), validate="one_to_one")
        audit["omnibus_reference_max_abs_difference"] = {
            "Kruskal_H": float(np.max(np.abs(merged["Kruskal_H_new"] - merged["Kruskal_H_reference"]))),
            "P_raw": float(np.max(np.abs(merged["P_raw_new"] - merged["P_raw_reference"]))),
            "P_adj_BH_across_19_features": float(
                np.max(
                    np.abs(
                        merged["P_adj_BH_across_19_features_new"]
                        - merged["P_adj_BH_across_19_features_reference"]
                    )
                )
            ),
        }

    (args.output / "06_recalculation_audit.json").write_text(
        json.dumps(audit, indent=2), encoding="utf-8"
    )
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
