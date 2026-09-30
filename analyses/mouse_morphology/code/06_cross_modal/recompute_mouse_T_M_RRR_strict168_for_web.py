#!/usr/bin/env python3
"""Recompute Mouse T-M RRR for the strict RPCA-Pearson D1/D2 morph consensus cohort.

This script updates the interactive RRR editor payload:
  interactive/mouse-morph-feature-explorer/dist/rrr_scores.csv
  interactive/mouse-morph-feature-explorer/dist/rrr_loadings.csv
and writes an audit copy under outputs/morph_qc/mouse_T_M_RRR_strict168_web/.
"""

from __future__ import annotations

from pathlib import Path
import sys

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

DIST = ROOT / "interactive/mouse-morph-feature-explorer/dist"
ASSIGN = DIST / "data.csv"
RNA = ROOT / "outputs/morph_qc/morph_taxonomy_round5_biological_validation/43_transcriptomic_log1p_HVG2000_cell_by_gene.csv"
T_AUDIT = ROOT / "outputs/20260820_REFonly_D1D2_marker_classifier_scan/D1D2_marker_classifier_Query_predictions.csv"
OUT = ROOT / "outputs/morph_qc/mouse_T_M_RRR_strict168_web"

FEATURES = [
    "M_soma_circularity_index",
    "M_soma_aspect_ratio",
    "M_cell_max_radial_dist",
    "M_total_number_of_neurites",
    "M_basal_dendrite_avg_tortuosity",
    "M_Total_neurite_length_(sections)",
    "M_Number_of_bifurcation_points",
    "M_Maximum_branch_order",
    "M_trunk_angle_min",
    "M_trunk_angle_max",
]

M_NAMES = [
    "Soma circularity",
    "Soma aspect ratio",
    "Max radial distance",
    "Primary neurite number",
    "Mean tortuosity",
    "Total neurite length",
    "Bifurcation points",
    "Maximum branch order",
    "Minimum trunk angle",
    "Maximum trunk angle",
]


def preprocess_train_test(train: np.ndarray, test: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    minimum = train.min(axis=0)
    shifted_train = np.maximum(train - minimum, 0)
    shifted_test = np.maximum(test - minimum, 0)
    total = shifted_train.sum(axis=0)
    total[total <= 0] = 1
    log_train = np.log1p(shifted_train / total * 10000)
    log_test = np.log1p(shifted_test / total * 10000)
    mean = log_train.mean(axis=0)
    sd = log_train.std(axis=0, ddof=1)
    sd[sd <= 0] = 1
    return (log_train - mean) / sd, (log_test - mean) / sd


def fit_rrr(x: np.ndarray, y: np.ndarray, rank: int):
    x_mean = x.mean(axis=0)
    y_mean = y.mean(axis=0)
    xc = x - x_mean
    yc = y - y_mean
    beta = np.linalg.pinv(xc) @ yc
    fitted = xc @ beta
    _, _, vt = np.linalg.svd(fitted, full_matrices=False)
    v = vt.T[:, :rank]
    beta_rrr = beta @ v @ v.T
    return x_mean, y_mean, beta_rrr, v


def orient_axes(v: np.ndarray, _feature_names: list[str]) -> np.ndarray:
    v = v.copy()
    for axis in range(v.shape[1]):
        anchor = int(np.argmax(np.abs(v[:, axis])))
        if v[anchor, axis] < 0:
            v[:, axis] *= -1
    return v


def correlation_loadings(features: np.ndarray, scores: np.ndarray) -> np.ndarray:
    f = features - features.mean(axis=0, keepdims=True)
    s = scores - scores.mean(axis=0, keepdims=True)
    f_sd = np.sqrt(np.sum(f * f, axis=0, keepdims=True))
    s_sd = np.sqrt(np.sum(s * s, axis=0, keepdims=True))
    denom = f_sd.T @ s_sd
    return np.divide(f.T @ s, denom, out=np.zeros((f.shape[1], s.shape[1])), where=denom > 0)


def strict_boolean(series: pd.Series) -> pd.Series:
    return series.astype(str).str.upper().map({"TRUE": True, "FALSE": False}).fillna(False)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)

    assign = pd.read_csv(ASSIGN)
    audit = pd.read_csv(T_AUDIT).rename(columns={"cell_id": "MSN_unique_ID"})
    audit["strict_T_identity"] = audit["Final_CellType_stability"].str.replace("Stable_", "", regex=False)
    strict = audit[
        audit["Final_CellType_stability"].isin(["Stable_D1", "Stable_D2"])
        & strict_boolean(audit["agree_RPCA"])
        & strict_boolean(audit["agree_Pearson"])
    ][["MSN_unique_ID", "strict_T_identity", "Marker_major", "Pearson_major", "bootstrap_top1_fraction", "bootstrap_margin"]]

    data = assign.merge(strict, on="MSN_unique_ID", how="inner", validate="one_to_one")
    data = data[
        data["HC_GC_consensus"].astype(str).str.lower().eq("true")
        & data["M_class"].isin(["M1", "M2", "M3", "M4"])
        & data["strict_T_identity"].isin(["D1", "D2"])
    ].copy()
    data = data.sort_values("MSN_unique_ID").reset_index(drop=True)

    if len(data) != 168:
        raise RuntimeError(f"Expected 168 strict HC-GC consensus cells, found {len(data)}")

    rna = pd.read_csv(RNA)
    data = data.merge(rna, on="MSN_unique_ID", how="inner", validate="one_to_one")
    if len(data) != 168:
        raise RuntimeError(f"RNA merge changed cohort size to {len(data)}")

    # Keep the same Mouse T-M RRR processing policy as the previous web analysis:
    # rank the current label-blind HVG2000 within this cohort, retain 1,000 genes,
    # z-score genes, reduce transcriptomics to 20 PCs, and fit rank-3 RRR.
    candidates = [c for c in rna.columns if c != "MSN_unique_ID"]
    raw_t = data[candidates].to_numpy(float)
    variance = raw_t.var(axis=0, ddof=1)
    keep = np.argsort(variance)[::-1][:1000]
    genes = np.asarray(candidates, dtype=object)[keep]
    gene_matrix = raw_t[:, keep]
    gene_sd = gene_matrix.std(axis=0, ddof=0)
    gene_sd[gene_sd <= 0] = 1
    gene_matrix = (gene_matrix - gene_matrix.mean(axis=0)) / gene_sd
    gene_matrix = np.nan_to_num(gene_matrix, nan=0.0, posinf=0.0, neginf=0.0)

    pca = PCA(n_components=20, svd_solver="full")
    t_scores = pca.fit_transform(gene_matrix)

    m_raw = data[FEATURES].to_numpy(float)
    m_matrix, _ = preprocess_train_test(m_raw, m_raw)

    xm, ym, _, v_raw = fit_rrr(t_scores, m_matrix, rank=3)
    v = orient_axes(v_raw, FEATURES)
    beta = np.linalg.pinv(t_scores - xm) @ (m_matrix - ym)

    t_projection = (t_scores - xm) @ beta @ v
    m_projection = (m_matrix - ym) @ v
    t_loadings = correlation_loadings(gene_matrix, t_projection)
    m_loadings = correlation_loadings(m_matrix, m_projection)

    scores = pd.DataFrame(
        {
            "MSN_unique_ID": data["MSN_unique_ID"],
            "M_class": data["M_class"],
            "T_identity": data["strict_T_identity"],
            "T_RRR1": t_projection[:, 0],
            "T_RRR2": t_projection[:, 1],
            "T_RRR3": t_projection[:, 2],
            "M_RRR1": m_projection[:, 0],
            "M_RRR2": m_projection[:, 1],
            "M_RRR3": m_projection[:, 2],
        }
    )
    loadings = pd.concat(
        [
            pd.DataFrame(
                {
                    "Domain": "Transcriptomic",
                    "Feature": genes,
                    "RRR1": t_loadings[:, 0],
                    "RRR2": t_loadings[:, 1],
                    "RRR3": t_loadings[:, 2],
                }
            ),
            pd.DataFrame(
                {
                    "Domain": "Morphological",
                    "Feature": M_NAMES,
                    "RRR1": m_loadings[:, 0],
                    "RRR2": m_loadings[:, 1],
                    "RRR3": m_loadings[:, 2],
                }
            ),
        ],
        ignore_index=True,
    )

    scores.to_csv(DIST / "rrr_scores.csv", index=False)
    loadings.to_csv(DIST / "rrr_loadings.csv", index=False)
    scores.to_csv(OUT / "RRR_cell_scores_strict168.csv", index=False)
    loadings.to_csv(OUT / "RRR_loadings_strict168.csv", index=False)
    pd.DataFrame({"PC": np.arange(1, 21), "variance_fraction": pca.explained_variance_ratio_}).to_csv(
        OUT / "transcriptomic_PCA_variance_strict168.csv", index=False
    )
    pd.crosstab(data["strict_T_identity"], data["M_class"]).to_csv(OUT / "strict_D1D2_by_Mclass_counts.csv")
    (OUT / "README.txt").write_text(
        "Mouse T-M RRR strict primary web cohort.\n"
        "Cohort: n=168 HC-GC morphology-consensus cells with strict RPCA-Pearson concordant D1/D2 identity.\n"
        "Counts: D1=74, D2=94; M1=63, M2=22, M3=55, M4=28.\n"
        "T: top 1,000 variable genes re-ranked from the current label-blind HVG2000 matrix within the strict cohort; gene z-score; 20 PCs.\n"
        "M: 10 frozen morphology features using shift-sum-log1p-zscore preprocessing within the strict cohort.\n"
        "RRR: rank=3, refit from scratch after strict T filtering.\n",
        encoding="utf-8",
    )

    print("Updated", DIST / "rrr_scores.csv")
    print("Updated", DIST / "rrr_loadings.csv")
    print(f"n={len(scores)}")
    print(pd.crosstab(data["strict_T_identity"], data["M_class"]).to_string())


if __name__ == "__main__":
    main()
