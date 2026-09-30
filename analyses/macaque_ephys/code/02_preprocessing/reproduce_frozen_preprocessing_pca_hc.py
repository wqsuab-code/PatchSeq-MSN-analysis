#!/usr/bin/env python3
"""Reproduce the archived 390-cell transform, PCA and Ward partition.

This portable audit uses only analysis-ready files committed in this module.
It does not download or redistribute the third-party HMBA source archive.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment
from sklearn.cluster import AgglomerativeClustering
from sklearn.decomposition import PCA
from sklearn.metrics import adjusted_rand_score

MODULE = Path(__file__).resolve().parents[2]
PARAMS = json.loads((MODULE / "config/analysis_parameters.json").read_text(encoding="utf-8"))
FEATURES = PARAMS["features"]


def transform_frame(frame: pd.DataFrame) -> np.ndarray:
    """Apply the frozen full-data transform exactly as documented."""
    columns = []
    for feature in FEATURES:
        values = frame[feature].to_numpy(dtype=float)
        if "ratio" in feature.lower():
            transformed = np.log2(values)
        else:
            shifted = values - np.min(values)
            transformed = np.log1p(shifted / shifted.sum() * 10_000.0)
        sd = transformed.std(ddof=0)
        if not np.isfinite(sd) or sd == 0:
            raise ValueError(f"Invalid standard deviation for {feature}")
        columns.append((transformed - transformed.mean()) / sd)
    return np.column_stack(columns)


def best_label_mismatch(reference: np.ndarray, predicted: np.ndarray) -> int:
    ref_levels = sorted(pd.unique(reference))
    pred_levels = sorted(pd.unique(predicted))
    table = pd.crosstab(pd.Series(reference), pd.Series(predicted)).reindex(
        index=ref_levels, columns=pred_levels, fill_value=0
    )
    rows, cols = linear_sum_assignment(-table.to_numpy())
    return int(len(reference) - table.to_numpy()[rows, cols].sum())


def main() -> None:
    raw = pd.read_csv(MODULE / "data/02_frozen_input/final19_raw_complete_n390.csv")
    frozen_z = pd.read_csv(
        MODULE / "data/03_transformed/Ca_Pu_NAC_final19_transformed_Zscore_matrix_n390.csv"
    )
    frozen_pc = pd.read_csv(
        MODULE / "data/03_transformed/Ca_Pu_NAC_final19_PCA_scores_n390.csv"
    )
    assignments = pd.read_csv(
        MODULE / "data/04_frozen_classification/frozen_cell_assignments_and_tsne_n390.csv"
    )

    rebuilt_z = transform_frame(raw)
    archived_z = frozen_z.set_index("cell_label").loc[raw.cell_label, FEATURES].to_numpy(float)
    pca = PCA(n_components=len(FEATURES), svd_solver="full").fit(rebuilt_z)
    rebuilt_pc = pca.transform(rebuilt_z)
    pc_cols = [f"PC{i}" for i in range(1, len(FEATURES) + 1)]
    archived_pc = frozen_pc.set_index("cell_label").loc[raw.cell_label, pc_cols].to_numpy(float)
    signs = np.sign(np.sum(rebuilt_pc * archived_pc, axis=0))
    signs[signs == 0] = 1
    aligned_pc = rebuilt_pc * signs

    predicted = AgglomerativeClustering(n_clusters=4, linkage="ward").fit_predict(aligned_pc[:, :3])
    reference = assignments.set_index("cell_label").loc[raw.cell_label, "HC_class"].astype(str).to_numpy()
    report = {
        "cells": int(len(raw)),
        "features": int(len(FEATURES)),
        "max_abs_z_difference": float(np.max(np.abs(rebuilt_z - archived_z))),
        "max_abs_pc_difference_after_sign_alignment": float(np.max(np.abs(aligned_pc - archived_pc))),
        "pc1_pc3_cumulative_variance_percent": float(100 * pca.explained_variance_ratio_[:3].sum()),
        "ward_k4_ari": float(adjusted_rand_score(reference, predicted)),
        "ward_k4_minimum_label_mismatches": best_label_mismatch(reference, predicted),
    }
    print(json.dumps(report, indent=2))
    if report["max_abs_z_difference"] > 1e-10:
        raise SystemExit("Frozen Z-score matrix was not reproduced")
    if report["max_abs_pc_difference_after_sign_alignment"] > 1e-10:
        raise SystemExit("Frozen PCA scores were not reproduced")
    if report["ward_k4_ari"] < 1 - 1e-12 or report["ward_k4_minimum_label_mismatches"] != 0:
        raise SystemExit("Frozen Ward K=4 partition was not reproduced")


if __name__ == "__main__":
    main()
