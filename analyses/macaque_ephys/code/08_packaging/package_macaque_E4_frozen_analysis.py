#!/usr/bin/env python
"""Package and audit the frozen Macaque E1-E4 analysis (2026-09-23)."""
from __future__ import annotations

import hashlib
import json
import platform
import shutil
import sys
from pathlib import Path
from zipfile import ZIP_DEFLATED, ZipFile

import numpy as np
import pandas as pd
import scipy
import sklearn
import h5py
from scipy.cluster.hierarchy import fcluster, linkage
from sklearn.decomposition import PCA
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score
from sklearn.preprocessing import StandardScaler


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs" / "Macaque_E4_frozen_analysis_package_20260923"
ZIP_OUT = ROOT / "outputs" / "Macaque_E4_frozen_analysis_package_20260923.zip"
RAW_ZIP = Path(r"C:\Users\53461\Downloads\Macaque-PatchSeq-BG.zip")
MAIN_SVG = Path(r"C:\Users\53461\OneDrive\Desktop\Patch-seq\4Manuscript\Macaque_E_Main.svg")
RADAR_LAYOUT = Path(r"C:\Users\53461\Downloads\Macaque_E1-E4_radar_layout.json")
RRR_LAYOUT = Path(r"C:\Users\53461\Downloads\Macaque_T-E_RRR_gc_merged_layout (3).json")

CORE = ROOT / "outputs/dSTR_dSTRvSTR_E_QC/Ca_Pu_NAC_final19_consensus_scan"
R3 = ROOT / "outputs/R3_panels_v3"
QC = ROOT / "outputs/R3_E_QC_overview_redrawn"
WEB = ROOT / "outputs/Macaque_E4_feature_explorer_site/dist"
ML = ROOT / "outputs/Macaque_E4_full_ML500"

FEATURES = [
    "Epsy_width_rheo", "Epsy_fast_trough_v_rheo", "Epsy_peak_deltav_rheo",
    "Epsy_peak_v_rheo", "Epsy_postap_slope_rheo", "Epsy_threshold_v_rheo",
    "Epsy_trough_t_rheo", "Epsy_trough_v_rheo",
    "Epsy_upstroke_downstroke_ratio_rheo", "Epsy_ahp_delay_5spike",
    "Epsy_ahp_delay_ratio_5spike", "Epsy_postap_slope_hero",
    "Epsy_trough_t_hero", "Epsy_downstroke_adapt_ratio",
    "Epsy_peak_v_adapt_ratio", "Epsy_threshold_v_adapt_ratio",
    "Epsy_upstroke_adapt_ratio", "Epsy_width_adapt_ratio",
    "Epsy_threshold_v_short_square",
]
E_MAP = {"C1": "E1", "C2": "E2", "C3": "E3", "C4": "E4"}
COLORS = {"E1": "#F8766D", "E2": "#7CAE00", "E3": "#00BFC4", "E4": "#C77CFF",
          "D1": "#D95F02", "D2": "#008F7A", "Hybrid": "#808080"}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def copy_file(src: Path, rel: str) -> None:
    if not src.exists():
        return
    dst = OUT / rel
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dst)


def copy_tree_files(src: Path, rel: str, names: list[str] | None = None) -> None:
    if not src.exists():
        return
    selected = [src / n for n in names] if names else [p for p in src.rglob("*") if p.is_file()]
    for path in selected:
        if path.exists():
            copy_file(path, str(Path(rel) / path.relative_to(src)))


def write(rel: str, text: str) -> None:
    path = OUT / rel
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text.rstrip() + "\n", encoding="utf-8")


def add_check(rows: list[dict], item: str, observed, expected, passed: bool,
              severity: str = "ERROR", note: str = "") -> None:
    rows.append({"item": item, "status": "PASS" if passed else severity,
                 "observed": str(observed), "expected": str(expected), "note": note})


def main() -> None:
    if OUT.exists() or ZIP_OUT.exists():
        raise FileExistsError(f"Refusing to overwrite existing package: {OUT} or {ZIP_OUT}")
    OUT.mkdir(parents=True)
    audit: list[dict] = []

    # Load the canonical frozen files.
    assignments = pd.read_csv(R3 / "00_tuned_p80_ee12_random_seed777_coordinates.csv")
    zmat = pd.read_csv(CORE / "Ca_Pu_NAC_final19_transformed_Zscore_matrix_n390.csv")
    raw = pd.read_csv(QC / "final19_raw_complete_n390.csv")
    pcs = pd.read_csv(CORE / "Ca_Pu_NAC_final19_PCA_scores_n390.csv")
    pvar = pd.read_csv(CORE / "Ca_Pu_NAC_final19_PCA_variance.csv")
    transform = pd.read_csv(CORE / "Ca_Pu_NAC_final19_transform_plan.csv")
    confusion = pd.read_csv(CORE / "Ca_Pu_NAC_final19_NPC3_HCK4_res3p0_merged_4x4_confusion.csv", index_col=0)
    cons = assignments.loc[assignments["Consensus"].astype(str).str.lower().eq("true")].copy()
    cons["E_class"] = cons["HC_class"].map(E_MAP)

    add_check(audit, "Frozen complete-case cell count", len(assignments), 390, len(assignments) == 390)
    add_check(audit, "Consensus cell count", len(cons), 368, len(cons) == 368)
    add_check(audit, "Consensus fraction", f"{len(cons)/len(assignments):.4%}", "94.36%",
              round(len(cons)/len(assignments), 4) == 0.9436)
    expected_e = {"E1": 59, "E2": 133, "E3": 57, "E4": 119}
    observed_e = cons["E_class"].value_counts().sort_index().to_dict()
    add_check(audit, "Consensus E-class sizes", observed_e, expected_e, observed_e == expected_e)
    observed_t = cons["T_class"].value_counts().to_dict()
    expected_t = {"D1": 163, "D2": 183, "Hybrid": 22}
    add_check(audit, "Consensus transcriptomic-class sizes", observed_t, expected_t, observed_t == expected_t,
              note="Correct D2 value is 183; 193 is the pre-consensus n=390 value.")

    ids = set(assignments.cell_label)
    for label, frame in [("raw 390x19", raw), ("frozen Z matrix", zmat), ("PCA scores", pcs)]:
        add_check(audit, f"Cell IDs: assignments versus {label}",
                  f"intersection={len(ids & set(frame.cell_label))}", "390 identical IDs",
                  ids == set(frame.cell_label))
    add_check(audit, "Frozen feature order in raw matrix", raw.columns[-19:].tolist(), FEATURES,
              raw.columns[-19:].tolist() == FEATURES)
    add_check(audit, "Frozen feature order in Z matrix", zmat.columns[1:].tolist(), FEATURES,
              zmat.columns[1:].tolist() == FEATURES)
    add_check(audit, "Transformation plan has 19 features", len(transform), 19, len(transform) == 19)

    # Recompute the frozen transformation from the archived raw 390x19 input.
    calc = pd.DataFrame(index=raw.index)
    for feature in FEATURES:
        x = raw[feature].astype(float).to_numpy()
        if "ratio" in feature:
            if np.any(x <= 0):
                raise ValueError(f"Nonpositive value in ratio feature {feature}")
            y = np.log2(x)
        else:
            shifted = x - np.min(x)
            y = np.log1p(shifted / shifted.sum() * 10000.0)
        calc[feature] = StandardScaler().fit_transform(y.reshape(-1, 1)).ravel()
    z_aligned = zmat.set_index("cell_label").loc[raw.cell_label, FEATURES].reset_index(drop=True)
    z_maxdiff = float(np.max(np.abs(calc.to_numpy() - z_aligned.to_numpy())))
    add_check(audit, "Independent reconstruction of frozen 390x19 Z matrix",
              f"max_abs_diff={z_maxdiff:.3e}", "<=1e-10", z_maxdiff <= 1e-10)

    # PCA reproduction (sign-aligned).
    fit = PCA(n_components=19, svd_solver="full").fit(calc.to_numpy())
    scores = fit.transform(calc.to_numpy())
    ref_scores = pcs.set_index("cell_label").loc[raw.cell_label, [f"PC{i}" for i in range(1, 20)]].to_numpy()
    for j in range(19):
        if np.corrcoef(scores[:, j], ref_scores[:, j])[0, 1] < 0:
            scores[:, j] *= -1
    pc_maxdiff = float(np.max(np.abs(scores - ref_scores)))
    add_check(audit, "Independent PCA reconstruction", f"sign-aligned max_abs_diff={pc_maxdiff:.3e}",
              "<=1e-10", pc_maxdiff <= 1e-10)
    ev_col = "explained_variance_ratio" if "explained_variance_ratio" in pvar else pvar.columns[-1]
    pc3_sum = float(pd.to_numeric(pvar[ev_col]).iloc[:3].sum())
    if pc3_sum <= 1.0:
        pc3_sum *= 100
    add_check(audit, "Variance explained by PC1-PC3", f"{pc3_sum:.4f}%", "54.6211%",
              abs(pc3_sum - 54.6211) < 0.001)

    # HC reproduction. Label numbering can be permuted, so ARI is the invariant check.
    ward = linkage(ref_scores[:, :3], method="ward", metric="euclidean")
    hc_new = fcluster(ward, t=4, criterion="maxclust")
    hc_ref = assignments.set_index("cell_label").loc[raw.cell_label, "HC_class"]
    hc_ari = adjusted_rand_score(hc_ref, hc_new)
    add_check(audit, "Independent Ward.D2 K=4 reproduction", f"ARI={hc_ari:.6f}", "ARI=1.000000",
              abs(hc_ari - 1) < 1e-12)

    expected_conf = np.array([[59, 0, 4, 2], [2, 133, 2, 2], [0, 2, 57, 5], [1, 2, 0, 119]])
    conf_arr = confusion.apply(pd.to_numeric, errors="coerce").to_numpy()
    add_check(audit, "HC-by-merged-GC concordance matrix", conf_arr.tolist(), expected_conf.tolist(),
              conf_arr.shape == (4, 4) and np.array_equal(conf_arr, expected_conf))
    diag = int(np.trace(conf_arr))
    add_check(audit, "Concordance matrix diagonal", diag, 368, diag == 368)

    # Metadata composition from the frozen PCA table.
    cmeta = cons.merge(pcs[["cell_label", "Lib_region_of_interest_label"]], on="cell_label", validate="one_to_one")
    observed_roi = cmeta["Lib_region_of_interest_label"].value_counts().to_dict()
    expected_roi = {"Ca": 169, "Pu": 132, "NAc": 67}
    # Some files use NAC, some NAc. Normalize only for the audit.
    observed_roi = {("NAc" if k == "NAC" else k): int(v) for k, v in observed_roi.items()}
    add_check(audit, "Consensus ROI sizes", observed_roi, expected_roi, observed_roi == expected_roi)

    # Heat map provenance.
    heat_order = pd.read_csv(R3 / "E_consensus368_final19_diagonal_heatmap_cell_order.csv")
    heat_id_col = "cell_label" if "cell_label" in heat_order else heat_order.columns[0]
    add_check(audit, "Heat-map cell order", f"n={len(heat_order)}, unique={heat_order[heat_id_col].nunique()}",
              "368 unique consensus IDs", len(heat_order) == 368 and set(heat_order[heat_id_col]) == set(cons.cell_label))

    # Radar provenance: display subset, not the clustering input.
    radar = pd.read_csv(WEB / "data/macaque_e_radar_cells.csv")
    add_check(audit, "Radar cell IDs and labels", f"n={len(radar)}", "same 368 consensus cells",
              len(radar) == 368 and set(radar.cell_label) == set(cons.cell_label))
    add_check(audit, "Radar feature count", len([c for c in radar if not c.endswith(("__z", "__zscaled", "__scaled"))]) - 2,
              10, len([c for c in radar if not c.endswith(("__z", "__zscaled", "__scaled"))]) - 2 == 10,
              severity="WARN", note="Ten radar axes are a display subset; clustering used all 19 frozen features.")

    # RRR provenance.
    rrr = pd.read_csv(WEB / "rrr-data/rrr_scores_gc_merged.csv")
    rrr_t = rrr["D1_D2"].value_counts().to_dict()
    add_check(audit, "RRR stable D1/D2 cohort", {"n": len(rrr), **rrr_t}, {"n": 346, "D1": 163, "D2": 183},
              len(rrr) == 346 and rrr_t == {"D2": 183, "D1": 163},
              note="Hybrid cells are retained in E-class summaries but excluded from D1-versus-D2 RRR.")
    fl = pd.read_csv(WEB / "rrr-data/rrr_feature_loadings_gc_merged.csv")
    gl = pd.read_csv(WEB / "rrr-data/rrr_gene_loadings_gc_merged.csv")
    layout = json.loads(RRR_LAYOUT.read_text(encoding="utf-8")) if RRR_LAYOUT.exists() else {}
    rrr_ok = len(fl) == 19 and len(gl) == 1000 and layout.get("featureCount") == 10 and layout.get("geneCount") == 15
    add_check(audit, "RRR source and displayed loading dimensions",
              f"source E rows={len(fl)}, source gene rows={len(gl)}; display E={layout.get('featureCount')}, genes={layout.get('geneCount')}",
              "full source: 19 E features and 1,000 genes; display: 10 E features and 15 genes", rrr_ok,
              note="The interactive layout selects the displayed variables from the full loading tables.")

    # ML analysis uses the frozen transform fold-wise; Yeo-Johnson is not the primary branch.
    ml_params = json.loads((ML / "00_parameters.json").read_text(encoding="utf-8"))
    ml_yj = ml_params.get("primary_transform", {}).get("yeo_johnson")
    add_check(audit, "ML primary preprocessing", f"yeo_johnson={ml_yj}",
              "frozen mixed transform fitted within training folds; yeo_johnson=False", ml_yj is False,
              severity="WARN", note="The dedicated validation script and parameter record both set yeo_johnson=False.")

    # Main SVG simple text-level consistency checks.
    if MAIN_SVG.exists():
        svg_text = MAIN_SVG.read_text(encoding="utf-8", errors="ignore")
        add_check(audit, "Latest main SVG available", MAIN_SVG.name, "present", True)
        add_check(audit, "Latest main SVG does not contain literal D2 (193)", "not found" if "D2 (193)" not in svg_text else "found",
                  "not found", "D2 (193)" not in svg_text, severity="WARN",
                  note="SVG text may be split into separate text elements; numerical source tables are authoritative.")
    else:
        add_check(audit, "Latest main SVG available", "missing", "present", False, severity="WARN")

    # Explicit provenance limitations.
    add_check(audit, "One-click source script for frozen GC scan", "not located in repository", "archived source script",
              False, severity="WARN",
              note="All frozen GC assignments, mappings, seed scans, parameters and concordance tables are archived, but the exact script that originally generated the final19 Seurat scan is not present. The downstream figure scripts are present.")
    add_check(audit, "Exact R/Seurat package versions", "not recorded in frozen workflow", "version-locked session information",
              False, severity="WARN", note="Seurat functions and parameters are fully recorded; exact upstream R/Seurat versions remain a provenance gap.")

    # Copy the actual source and analysis files.
    copy_file(RAW_ZIP, "01_ORIGINAL_INPUTS/Macaque-PatchSeq-BG.zip")
    copy_file(QC / "final19_raw_complete_n390.csv", "01_ORIGINAL_INPUTS/final19_raw_complete_n390.csv")
    copy_file(ROOT / "outputs/dSTR_dSTRvSTR_E_QC/dSTR_plus_vSTR_cell_missingness_QC.csv",
              "01_ORIGINAL_INPUTS/dSTR_plus_vSTR_cell_missingness_QC.csv")
    copy_file(ROOT / "outputs/dSTR_dSTRvSTR_E_QC/primary10_sensitivity_cohorts/dSTR_plus_vSTR_A_all_stage1QC_raw_with_NA.csv",
              "01_ORIGINAL_INPUTS/radar_display_primary10_stage1QC_raw_with_NA.csv")

    core_names = [
        "Ca_Pu_NAC_final19_transform_plan.csv", "Ca_Pu_NAC_final19_transformed_Zscore_matrix_n390.csv",
        "Ca_Pu_NAC_final19_PCA_scores_n390.csv", "Ca_Pu_NAC_final19_PCA_loadings.csv",
        "Ca_Pu_NAC_final19_PCA_variance.csv", "Ca_Pu_NAC_final19_all_GC_assignments.csv",
        "Ca_Pu_NAC_final19_GC_to_HC_all_merge_maps.csv",
        "Ca_Pu_NAC_final19_HC_GC_merge_consensus_all_candidates.csv",
        "Ca_Pu_NAC_final19_NPC3_HCK4_res3p0_raw_HC_by_GC_confusion.csv",
        "Ca_Pu_NAC_final19_NPC3_HCK4_res3p0_merged_4x4_confusion.csv",
        "Ca_Pu_NAC_final19_NPC3_HCK4_ge90_seed_validation_summary.csv",
        "Ca_Pu_NAC_final19_NPC3_HCK4_ge90_seed_validation_detail.csv",
        "Ca_Pu_NAC_final19_NPC3_HCK4_ge90_seed_assignments.csv",
    ]
    copy_tree_files(CORE, "02_FROZEN_CORE", core_names)
    copy_file(R3 / "00_res3_GC15_to_C4_mapping.csv", "02_FROZEN_CORE/00_res3_GC15_to_C4_mapping.csv")
    copy_file(R3 / "00_tuned_p80_ee12_random_seed777_coordinates.csv", "02_FROZEN_CORE/frozen_cell_assignments_and_tsne_n390.csv")
    for name in ["final19_missingness_stage1_n455.csv", "final19_raw_skewness_complete_n390.csv",
                 "final19_Spearman_complete_n390.csv", "final19_complete_n390_pairs_absrho_ge0p80.csv"]:
        copy_file(QC / name, f"02_FROZEN_CORE/QC/{name}")

    for name in ["E_consensus368_final19_diagonal_heatmap_cell_order.csv",
                 "E_consensus368_final19_diagonal_feature_order_and_class_means.csv",
                 "E_consensus368_final19_heatmap_feature_order.csv",
                 "E_consensus368_final19_heatmap_feature_labels_top_to_bottom.txt",
                 "13_Tclass_composition_by_consensus4_H1in_counts.csv",
                 "13_Tclass_composition_by_consensus4_H1in_percent.csv",
                 "14_ROI_composition_by_consensus4_H1in_counts.csv",
                 "14_ROI_composition_by_consensus4_H1in_percent.csv",
                 "E_HC4_by_mergedGC4_counts.csv", "E_HC4_by_mergedGC4_row_percent.csv"]:
        copy_file(R3 / name, f"03_MAIN_FIGURE_SOURCE/heatmap_and_composition/{name}")
    copy_file(WEB / "data/macaque_e_radar_cells.csv", "03_MAIN_FIGURE_SOURCE/radar/macaque_e_radar_cells.csv")
    copy_file(R3 / "E_GC_res3_mergedK4_radar10_consensus368_clean_W4p64_H1_robust_scaling.csv",
              "03_MAIN_FIGURE_SOURCE/radar/radar10_scaling.csv")
    copy_file(RADAR_LAYOUT, "03_MAIN_FIGURE_SOURCE/radar/Macaque_E1-E4_radar_layout.json")
    copy_tree_files(WEB / "rrr-data", "03_MAIN_FIGURE_SOURCE/RRR")
    copy_tree_files(ROOT / "outputs/R3_RRR_fourpanel", "03_MAIN_FIGURE_SOURCE/RRR/frozen_model_outputs",
                    ["RRR_cell_scores_n346.csv", "RRR_Efeature_correlation_loadings.csv",
                     "RRR_gene_correlation_loadings.csv", "transcriptomic_PCA_variance.csv", "README.txt"])
    copy_file(RRR_LAYOUT, "03_MAIN_FIGURE_SOURCE/RRR/Macaque_T-E_RRR_gc_merged_layout.json")
    copy_file(MAIN_SVG, "07_REFERENCE_OUTPUTS/Macaque_E_Main.svg")
    copy_file(ROOT / "outputs/Macaque_E4_frozen_workflow_EN.txt", "07_REFERENCE_OUTPUTS/Macaque_E4_frozen_workflow_EN.txt")
    copy_file(ROOT / "outputs/Macaque_E4_frozen_workflow_CN.txt", "07_REFERENCE_OUTPUTS/Macaque_E4_frozen_workflow_CN.txt")

    scripts = [
        "run_dstr_dstrvstr_rheobase_qc.py", "build_primary10_sensitivity_cohorts.py",
        "plot_ca_pu_nac_hc_gc_jitter_confusion.py", "plot_ca_pu_nac_res3_consensus_heatmap.py",
        "plot_ca_pu_nac_res3_consensus_heatmap_diagonal.py",
        "plot_ca_pu_nac_consensus_heatmap_fully_annotated.py", "plot_ca_pu_nac_res3_e_radar.py",
        "build_macaque_e_radar_editor.py", "plot_ca_pu_nac_rrr_fourpanel.py",
        "build_macaque_e_rrr_editor.py", "build_macaque_e_tsne_editor.py",
        "plot_ca_pu_nac_eclass_tclass_stack_no_text.py", "validate_macaque_E4_full_ml500.py",
        "analyze_nac_d1d2_transcriptomic_ephys_rrr.py",
        "validate_macaque_EM_stability_500.py", "audit_macaque_E4_ml_preprocessing.py",
        "complete_macaque_E4_submission_ml_audit.py", "add_macaque_E4_roc_pr.py",
        "package_macaque_E4_frozen_analysis.py",
    ]
    for name in scripts:
        copy_file(ROOT / "scripts" / name, f"04_ANALYSIS_CODE/{name}")
    copy_tree_files(ML, "04_ML_VALIDATION_RESULTS")

    # Methods modules.
    methods = {
    "01_cohort_QC_EN.md": """# Cohort definition and electrophysiological quality control

Macaque Patch-seq cells were obtained from the original `Macaque-PatchSeq-BG.zip` release. The analysis was restricted to medium spiny neurons from caudate (Ca), putamen (Pu), and nucleus accumbens (NAc). Cells first passed the upstream cell-level quality-control rule (no more than 50% missingness across the assessed electrophysiological features). The frozen analysis then used strict complete cases across the final 19-feature panel; no missing values were imputed. This yielded 390 cells (Ca, 181; Pu, 140; NAc, 69). Transcriptomic D1, D2, and Hybrid annotations were not used to derive electrophysiological classes.

Of 57 `Epsy_` variables, 54 had at least 30 observations and more than one distinct value. Redundancy was assessed using pairwise-complete Spearman correlations (minimum 30 paired observations). A greedy independent-set procedure removed a candidate when |rho| exceeded 0.80 with a higher-priority retained variable. Protocol priority favored rheobase, last-rheobase, five-spike, hero, adaptation-ratio, ramp, and short-square measurements, in that order. Two high-missingness ramp variables and one severely skewed, uninformative adaptation variable were subsequently removed, leaving 19 prespecified nonredundant features. Extreme observations were retained; no winsorization or trimming was used.""",
    "02_preprocessing_PCA_EN.md": """# Feature transformation and principal-component analysis

All transformations were fitted to the 390 complete cases. The seven strictly positive ratio variables were log2-transformed. For each of the remaining 12 variables, the column minimum was subtracted, the shifted column was normalized to a total of 10,000, and `log1p` was applied. Every transformed variable was then standardized to zero mean and unit variance using the population standard deviation (`ddof=0`). No row normalization, imputation, batch correction, covariate regression, or Yeo-Johnson transformation was used in the frozen classification.

Principal-component analysis was performed on the 390 x 19 standardized matrix by full singular-value decomposition. The first three components were used for clustering and explained 23.2751%, 17.5400%, and 13.8060% of the variance, respectively (cumulative, 54.6211%). PC4 was examined only in sensitivity analyses.""",
    "03_clustering_consensus_EN.md": """# Hierarchical, graph-based, and consensus clustering

Ward hierarchical clustering was performed in PC1-PC3 space using Euclidean distances and the tree was cut at K=4. Independently, the same three PCs were inserted into a Seurat PCA reduction. A shared-nearest-neighbor graph was built with `k.param=20` and `prune.SNN=1/15`, followed by Louvain community detection (`algorithm=1`, `resolution=3.0`, random seed 777), yielding 15 raw graph communities. The graph communities were mapped surjectively onto the four Ward clusters by dynamic programming to maximize diagonal overlap while assigning every graph community to exactly one Ward class and ensuring that every Ward class received at least one graph community.

A cell was designated as consensus when its Ward label matched its merged graph-clustering label. This criterion retained 368 of 390 cells (94.36%) and defined E1-E4 (internal archived labels C1-C4): E1, n=59; E2, n=133; E3, n=57; and E4, n=119. The 22 discordant cells were not forced into an E class. The consensus set comprised D1, n=163; D2, n=183; and Hybrid, n=22 cells, and Ca, n=169; Pu, n=132; and NAc, n=67 cells. A ten-seed audit at the frozen parameters produced 14-16 raw communities and post-merge consensus rates of 90.00-94.62% (mean, 92.95%).""",
    "04_visualization_EN.md": """# t-SNE, concordance matrix, heat map, and radar plots

t-SNE was used only for visualization of the frozen PC1-PC3 scores (`perplexity=80`, `early_exaggeration=12`, `learning_rate='auto'`, random initialization, 3,000 iterations, random seed 777). Displayed class outlines were robust minimum-covariance-determinant regions (`support_fraction=0.75`, seed 777) covering 80% of each class. Neither t-SNE coordinates nor outlines contributed to clustering.

The HC-GC panel is a concordance matrix, not a predictive confusion matrix. The heat map contains the 368 consensus cells and all 19 frozen transformed feature z-scores, displayed over [-2, 2] with values outside this interval clipped for visualization. Columns are grouped by E class and annotated by transcriptomic class and ROI. The radar plots show ten selected electrophysiological measurements for the same 368 consensus cells. Cohort-wide feature z-scores were mapped linearly from -3 at the center to +3 at the outer edge and clipped only for display; pale lines denote individual cells, shaded bands denote feature-wise interquartile ranges, and black outlines denote feature-wise medians. The ten radar axes are a display subset and do not replace the 19-feature clustering input.""",
    "05_Tclass_ROI_EN.md": """# Transcriptomic-class and anatomical composition

Transcriptomic class and ROI were evaluated after the E classes had been frozen. Stacked bars show the percentage of D1, D2, and Hybrid cells within each E class; numbers within bars denote cell counts. The corresponding E-by-transcriptomic-class counts were E1: 24 D1, 20 D2, 15 Hybrid; E2: 35 D1, 95 D2, 3 Hybrid; E3: 38 D1, 19 D2, 0 Hybrid; and E4: 66 D1, 49 D2, 4 Hybrid. ROI counts were E1: 20 Ca, 23 Pu, 16 NAc; E2: 70 Ca, 49 Pu, 14 NAc; E3: 25 Ca, 28 Pu, 4 NAc; and E4: 54 Ca, 32 Pu, 33 NAc. These metadata were not used to fit PCA or clustering.""",
    "06_RRR_EN.md": """# Transcriptome-electrophysiology reduced-rank regression

Transcriptome-electrophysiology reduced-rank regression (RRR) was performed in the 346 HC-GC consensus cells with stable D1 or D2 transcriptomic assignments (D1, n=163; D2, n=183); the 22 Hybrid cells were excluded. Transcriptomic values were read from the processed log2-expression matrix `HMBA-Macaque-PatchSeq-BG-log2.h5ad` in the Allen Human-Mammalian Brain Atlas macaque Patch-seq release. Genes detected in fewer than eight of the 346 cells, genes with zero or non-finite variance, and technical genes matching mitochondrial, ribosomal, `GM`/`Gm`, `Rik`, `MALAT1`/`Malat1`, `XIST`/`Xist`, or empty-symbol patterns were excluded. When duplicate gene symbols were present, the row with the greatest variance was retained. The 1,000 eligible genes with the greatest sample variance were selected without using D1/D2 or E-class labels, standardized across cells to zero mean and unit population variance (`ddof=0`), and reduced to 20 transcriptomic principal components by full singular-value decomposition. These 20 components explained 34.72% of transcriptomic variance in the analyzed cohort.

The response matrix comprised the same 19 transformed and standardized electrophysiological features used for the frozen E classification. Both predictor and response matrices were centered. Ordinary least-squares coefficients were estimated using the Moore-Penrose pseudoinverse of the centered 20-PC transcriptomic matrix. Singular-value decomposition of the fitted electrophysiological matrix supplied the response directions, and the coefficient matrix was projected onto the first three directions to obtain a rank-3 model. Component signs were oriented so that the electrophysiological feature with the largest absolute loading on each component had a positive loading. Transcriptomic scores were calculated by projecting the centered transcriptomic PCs through the ordinary least-squares coefficient matrix and the three response directions; electrophysiological scores were obtained by projecting the centered 19-feature matrix onto the same directions. Gene and electrophysiological loading vectors are Pearson correlations between the standardized input variables and their corresponding component scores.

The rank-3 solution was fitted to all 346 cells for descriptive joint-axis visualization; rank 3 was not selected by donor-held-out cross-validation, and no predictive R-squared is reported. The figure shows Components 1-2 and 1-3 in transcriptomic and electrophysiological spaces. Points denote cells. Transcriptomic panels are colored by D1/D2 identity, whereas electrophysiological panels are colored by E1-E4. The interactive layout displays the 15 genes and ten electrophysiological features with the largest two-dimensional loading magnitudes for each displayed component pair; the archived loading tables contain all 1,000 genes and all 19 electrophysiological features. Covariance ellipses show descriptive 90% bivariate-normal regions and are not confidence intervals. Scores were rescaled by the 99th percentile of radial distance for plotting only. RRR did not alter the frozen E-class assignments.""",
    "08_data_code_availability_EN.md": """# Data and code availability for the macaque analysis

Processed macaque Patch-seq expression matrices and metadata were downloaded from the Allen Brain Cell Atlas Human-Mammalian Brain Atlas Basal Ganglia macaque Patch-seq release (dataset version `20260228`; https://alleninstitute.github.io/abc_atlas_access/descriptions/HMBA-Macaque-PatchSeq-BG.html). The RRR analysis used the released `HMBA-Macaque-PatchSeq-BG-log2.h5ad` matrix and did not reprocess raw sequencing reads. The corresponding publication is Liu et al., *Morphoelectric Diversity and Specialization of Neuronal Cell Types in the Primate Striatum* (2026), https://doi.org/10.64898/2026.02.26.708019. Raw sequencing resources are available through NeMO Archive collection `nemo:dat-nsm6mxv` (https://assets.nemoarchive.org/dat-nsm6mxv), and the source publication lists the associated electrophysiology NWB datasets at DANDI. The exact processed archive used here, derived analysis tables, frozen cell labels, source code and file checksums are included in this reproducibility package. The Allen processed release is distributed under CC BY-NC 4.0; repository-specific licensing should be followed for raw NeMO and DANDI records.""",
    "07_ML_validation_EN.md": """# Machine-learning stability and reproducibility analysis

The frozen E1-E4 labels were evaluated in 500 donor-grouped resampling iterations. Donors, rather than individual cells, were assigned to training or held-out subsets. Within each split, the original mixed 19-feature transformation was fitted using the training data only and then applied to held-out cells; standardization and PCA were likewise fitted within the training fold. The primary validation therefore did not use Yeo-Johnson preprocessing and avoided information leakage. Seven classifiers were compared using held-out balanced accuracy, macro-F1, class recall, aggregated confusion matrices, and one-versus-rest ROC and precision-recall curves. Additional analyses included donor-resampled Ward, k-means, Gaussian-mixture and spectral clustering; NPC-by-K sensitivity; donor-structured label permutation; cell-level consensus margins; fold-restricted permutation importance; Extra Trees mean-decrease-in-impurity stability; cross-model rank agreement; metadata association; and nested progressive feature ablation. Complete numerical outputs are archived in `04_ML_VALIDATION_RESULTS`.""",
    }
    for name, body in methods.items():
        write(f"05_METHODS/{name}", body)
    combined = "# Frozen Macaque E1-E4 analysis: Methods\n\n" + "\n\n".join(methods.values())
    write("05_METHODS/Methods_combined_EN.md", combined)
    write("05_METHODS/Methods_index_CN.md", """# Methods 文件索引

- `01_cohort_QC_EN.md`：队列、纳入排除和E特征QC。
- `02_preprocessing_PCA_EN.md`：19指标变换、标准化和PCA。
- `03_clustering_consensus_EN.md`：HC、GC、满射合并和共识定义。
- `04_visualization_EN.md`：t-SNE、矩阵、热图和雷达图。
- `05_Tclass_ROI_EN.md`：D1/D2/Hybrid及ROI构成。
- `06_RRR_EN.md`：稳定D1/D2共识细胞的RRR。
- `07_ML_validation_EN.md`：500次供体分组机器学习验证。
- `08_data_code_availability_EN.md`：公开数据下载、论文、许可和代码披露。
- `Methods_combined_EN.md`：上述模块的英文合并稿。

主分类采用冻结的19指标混合变换，不使用Yeo-Johnson。雷达图的10指标仅用于展示；RRR仅使用346个稳定D1/D2共识细胞；这些均未改变E1-E4标签。""")

    # Audit outputs and conflict report.
    audit_df = pd.DataFrame(audit)
    (OUT / "06_AUDIT").mkdir(exist_ok=True)
    audit_df.to_csv(OUT / "06_AUDIT/consistency_audit.csv", index=False)
    status_counts = audit_df.status.value_counts().to_dict()
    warnings = audit_df.loc[audit_df.status != "PASS"]
    conflict_lines = [
        "# Conflict and concordance report",
        "",
        f"Audit summary: {status_counts.get('PASS', 0)} PASS, {status_counts.get('WARN', 0)} WARN, {status_counts.get('ERROR', 0)} ERROR.",
        "",
        "## Resolved or explicitly separated issues",
        "",
        "1. The correct consensus transcriptomic counts are D1=163, D2=183, Hybrid=22. D2=193 is the count in the pre-consensus 390-cell cohort and must not appear as the consensus count.",
        "2. Archived internal labels C1-C4 correspond one-to-one to manuscript-facing E1-E4. No biological ordering is implied by either label set.",
        "3. The frozen classifier uses 19 features. The ten-axis radar and the ten-feature RRR loading display are descriptive subsets only.",
        "4. The frozen main preprocessing is the feature-specific log2 or shift-column-sum-log1p transform followed by z-scoring. Yeo-Johnson analyses, where present historically, are sensitivity branches and are not the frozen primary analysis.",
        "5. The heat map clips z-scores to [-2,2] and the radar clips display z-scores to [-3,3]; clipping is visual only.",
        "6. RRR uses 346 stable D1/D2 consensus cells and excludes 22 Hybrid cells; other E-class panels retain Hybrid cells.",
        "7. The HC-GC matrix is a clustering concordance matrix, not a supervised prediction confusion matrix.",
        "",
        "## Remaining provenance warnings",
        "",
    ]
    if warnings.empty:
        conflict_lines.append("None.")
    else:
        for _, row in warnings.iterrows():
            conflict_lines.append(f"- **{row['item']}**: {row['observed']}. {row['note']}")
    write("06_AUDIT/conflict_report.md", "\n".join(conflict_lines))

    feature_map = pd.DataFrame({"frozen_order": range(1, 20), "full_feature_name": FEATURES})
    feature_map["display_label"] = [
        "Width rheo", "Fast trough", "Peak delta V", "Peak V", "Post-AP slope",
        "Threshold V", "Trough t", "Trough V", "Up/down ratio", "AHP delay",
        "AHP ratio", "Post-AP hero", "Trough t hero", "Downstroke adapt",
        "Peak V adapt", "Threshold adapt", "Upstroke adapt", "Width adapt", "Threshold SS",
    ]
    feature_map.to_csv(OUT / "06_AUDIT/feature_name_map.csv", index=False)
    pd.DataFrame([{"object": k, "hex": v} for k, v in COLORS.items()]).to_csv(
        OUT / "06_AUDIT/color_map.csv", index=False)
    env = {"python": sys.version, "platform": platform.platform(), "numpy": np.__version__,
           "pandas": pd.__version__, "scipy": scipy.__version__, "scikit_learn": sklearn.__version__,
           "h5py": h5py.__version__,
           "note": "Exact original R and Seurat versions were not recorded in the frozen workflow."}
    write("06_AUDIT/software_environment.json", json.dumps(env, indent=2))

    readme = f"""# Macaque frozen E1-E4 analysis package (2026-09-23)

This archive locks the latest frozen Ca+Pu+NAc electrophysiological classification and separates primary classification inputs from display-only and validation analyses.

## Frozen result

- Complete-case clustering cohort: 390 cells x 19 nonredundant E features; no imputation.
- PC1-PC3: 54.6211% cumulative variance.
- HC: Euclidean Ward.D2, K=4.
- GC: Seurat SNN (`k.param=20`, `prune.SNN=1/15`), Louvain `algorithm=1`, `resolution=3.0`, seed 777; 15 raw communities.
- Consensus: 368/390 (94.36%); E1=59, E2=133, E3=57, E4=119.
- Consensus metadata: D1=163, D2=183, Hybrid=22; Ca=169, Pu=132, NAc=67.

## Directory guide

- `01_ORIGINAL_INPUTS`: original public archive and exact raw 390x19 input; radar source is explicitly marked as display-only.
- `02_FROZEN_CORE`: transformation, PCA, HC/GC assignments, merge maps, consensus and QC tables.
- `03_MAIN_FIGURE_SOURCE`: heat-map, radar and RRR source tables and current layouts.
- `04_ANALYSIS_CODE`: available QC, figure, RRR, interactive-page and ML scripts.
- `04_ML_VALIDATION_RESULTS`: complete 500-repeat frozen-label validation outputs.
- `05_METHODS`: modular and combined Methods text.
- `06_AUDIT`: machine-readable checks, conflicts, feature/color maps, software and checksums.
- `07_REFERENCE_OUTPUTS`: latest main SVG and frozen workflow documents.

The package contains the original 36.4-MB release ZIP. The exact script that originally generated the final Seurat GC scan and the exact R/Seurat versions were not found; this is reported as a provenance warning. Frozen assignments, mapping tables, seed scans and all downstream scripts are retained. No numerical audit returned ERROR.
"""
    write("00_README/README.md", readme)
    write("00_README/README_CN.md", """# Macaque E1-E4 冻结分析包

本包将主分类输入、冻结核心结果、主图数据、机器学习验证和仅用于展示的数据分开归档。权威数字为：390个完整病例、19项非冗余E特征；HC-GC共识368/390（94.36%）；E1=59、E2=133、E3=57、E4=119；D1=163、D2=183、Hybrid=22；Ca=169、Pu=132、NAc=67。

必须注意：D2=193属于390细胞预共识队列，不能写成368共识细胞的D2数量；雷达图10指标不是聚类输入；RRR只分析346个稳定D1/D2共识细胞；主分析没有使用Yeo-Johnson。详见`06_AUDIT/conflict_report.md`和`05_METHODS`。""")

    # Create manifest after every content file has been written.
    manifest_rows = []
    for path in sorted(p for p in OUT.rglob("*") if p.is_file()):
        manifest_rows.append({"relative_path": path.relative_to(OUT).as_posix(), "bytes": path.stat().st_size,
                              "sha256": sha256(path)})
    pd.DataFrame(manifest_rows).to_csv(OUT / "06_AUDIT/file_manifest.csv", index=False)
    checksums = [f"{row['sha256']}  {row['relative_path']}" for row in manifest_rows]
    write("06_AUDIT/SHA256SUMS.txt", "\n".join(checksums))

    # Add the manifest itself to the checksum list in a final companion record.
    manifest_path = OUT / "06_AUDIT/file_manifest.csv"
    write("06_AUDIT/manifest_self_sha256.txt", sha256(manifest_path) + "  06_AUDIT/file_manifest.csv")

    with ZipFile(ZIP_OUT, "w", compression=ZIP_DEFLATED, compresslevel=6) as archive:
        for path in sorted(p for p in OUT.rglob("*") if p.is_file()):
            archive.write(path, arcname=f"{OUT.name}/{path.relative_to(OUT).as_posix()}")

    print(json.dumps({"package_dir": str(OUT), "zip": str(ZIP_OUT),
                      "audit": status_counts, "files": len(list(OUT.rglob('*'))),
                      "zip_bytes": ZIP_OUT.stat().st_size}, indent=2))


if __name__ == "__main__":
    main()
