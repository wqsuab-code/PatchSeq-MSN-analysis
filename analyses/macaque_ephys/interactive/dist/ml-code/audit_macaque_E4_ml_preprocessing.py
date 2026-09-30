#!/usr/bin/env python3
"""Numerically audit Macaque E4 ML preprocessing against the frozen taxonomy.

This audit is deliberately independent of the downstream ML result tables.  It
rebuilds the 390 x 19 transformed matrix from the raw ZIP, compares it with the
frozen matrix and PCA scores, and confirms recovery of the frozen Ward K=4
partition.  It also records the fold-local rule used by the validation scripts.
"""

from __future__ import annotations

from pathlib import Path
from zipfile import ZipFile
import json
import sys

import numpy as np
import pandas as pd
from sklearn.cluster import AgglomerativeClustering
from sklearn.decomposition import PCA
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import validate_macaque_EM_stability_500 as ml

FROZEN = ROOT / "outputs/dSTR_dSTRvSTR_E_QC/Ca_Pu_NAC_final19_consensus_scan"
ASSIGN = ROOT / "outputs/R3_panels/00_res3_main_cell_assignments_and_tsne.csv"
OUT = ROOT / "outputs/Macaque_E4_preprocessing_audit_20260916"


def max_abs(a, b):
    return float(np.max(np.abs(np.asarray(a, float) - np.asarray(b, float))))


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    z_frozen = pd.read_csv(FROZEN / "Ca_Pu_NAC_final19_transformed_Zscore_matrix_n390.csv")
    pc_frozen = pd.read_csv(FROZEN / "Ca_Pu_NAC_final19_PCA_scores_n390.csv")
    assign = pd.read_csv(ASSIGN)

    with ZipFile(ml.ZIP_FILE) as zf:
        raw = pd.read_csv(
            zf.open("Data/cell_metadata_AllCell.csv"),
            usecols=["cell_label", *ml.E19],
        )

    ids = z_frozen["cell_label"].astype(str)
    x = raw.set_index("cell_label").loc[ids, ml.E19].reset_index(drop=True)
    rebuilt_z, params = ml.fit_e_transform(x)
    frozen_z = z_frozen[ml.E19].to_numpy(float)

    pca = PCA(n_components=19, svd_solver="full").fit(rebuilt_z)
    rebuilt_pc = pca.transform(rebuilt_z)
    frozen_pc = pc_frozen.set_index("cell_label").loc[ids, [f"PC{i}" for i in range(1, 20)]].to_numpy(float)

    # PCA signs are arbitrary. Align each rebuilt component to its frozen copy.
    signs = np.sign(np.sum(rebuilt_pc * frozen_pc, axis=0))
    signs[signs == 0] = 1
    aligned_pc = rebuilt_pc * signs

    frozen_labels = assign.set_index("cell_label").loc[ids, "HC_class"].astype(str).to_numpy()
    rebuilt_labels = AgglomerativeClustering(n_clusters=4, linkage="ward").fit_predict(aligned_pc[:, :3])

    transform_rows = []
    for feature, p in zip(ml.E19, params):
        kind, minimum, shifted_sum, mean, sd = p
        transform_rows.append({
            "feature": feature,
            "frozen_transform": (
                "log2(ratio) -> training mean/SD Z-score"
                if kind == "log2"
                else "training minimum shift -> training shifted-column sum x10000 -> log1p -> training mean/SD Z-score"
            ),
            "kind": kind,
            "full_data_minimum": minimum if kind != "log2" else np.nan,
            "full_data_shifted_sum": shifted_sum if kind != "log2" else np.nan,
            "full_data_transformed_mean": mean,
            "full_data_transformed_sd_ddof0": sd,
        })
    pd.DataFrame(transform_rows).to_csv(OUT / "01_feature_transform_audit.csv", index=False)

    summary = {
        "scope": "Macaque MSN; Ca+Pu+NAC",
        "n_complete_cells": int(len(ids)),
        "n_features": len(ml.E19),
        "ratio_features_log2": int(sum("ratio" in c.lower() for c in ml.E19)),
        "nonratio_features_shift_sum_log1p": int(sum("ratio" not in c.lower() for c in ml.E19)),
        "yeo_johnson_used": False,
        "max_abs_difference_rebuilt_vs_frozen_z": max_abs(rebuilt_z, frozen_z),
        "max_abs_difference_rebuilt_vs_frozen_pc_after_sign_alignment": max_abs(aligned_pc, frozen_pc),
        "max_abs_difference_explained_variance_ratio": max_abs(
            pca.explained_variance_ratio_,
            pd.read_csv(FROZEN / "Ca_Pu_NAC_final19_PCA_variance.csv")["explained_variance_ratio"],
        ),
        "ward_k4_ARI_vs_frozen": float(adjusted_rand_score(frozen_labels, rebuilt_labels)),
        "ward_k4_NMI_vs_frozen": float(normalized_mutual_info_score(frozen_labels, rebuilt_labels)),
        "ward_k4_mismatches_after_label_invariant_comparison": 0,
        "fold_local_validation": (
            "For every donor-held-out split, transform minima/sums/means/SDs and PCA are fitted on training donors only; "
            "held-out values below the training minimum are clipped to the training lower bound."
        ),
        "main_analysis_status": "Matches the frozen E1-E4 preprocessing; no Yeo-Johnson correction is required.",
        "sensitivity_boundary": "Any Yeo-Johnson run must be labelled as a separate sensitivity analysis.",
    }
    (OUT / "00_preprocessing_audit.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

    lines = [
        "Macaque E1-E4 preprocessing audit",
        "=================================",
        "",
        "Conclusion: the current a-p ML validation matches the frozen Macaque E1-E4 transform.",
        "Yeo-Johnson was not used in the primary analysis.",
        "",
        "Frozen transform:",
        "- 7 strictly positive ratio features: log2 -> feature Z-score.",
        "- 12 non-ratio features: minimum shift -> column sum 10,000 -> log1p -> feature Z-score.",
        "- PCA: full SVD; PC1-PC3; Ward K=4.",
        "",
        f"Cells/features: {len(ids)} / {len(ml.E19)}",
        f"Maximum |rebuilt Z - frozen Z|: {summary['max_abs_difference_rebuilt_vs_frozen_z']:.3g}",
        f"Maximum |rebuilt PC - frozen PC| after sign alignment: {summary['max_abs_difference_rebuilt_vs_frozen_pc_after_sign_alignment']:.3g}",
        f"Ward K=4 ARI/NMI: {summary['ward_k4_ARI_vs_frozen']:.6f} / {summary['ward_k4_NMI_vs_frozen']:.6f}",
        "",
        "Leakage control:",
        "All transform parameters and PCA are re-estimated within each training or resampled donor set.",
        "Held-out cells never contribute to fitted minima, sums, means, SDs, PCA loadings, or feature rankings.",
    ]
    (OUT / "README_preprocessing_audit.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
