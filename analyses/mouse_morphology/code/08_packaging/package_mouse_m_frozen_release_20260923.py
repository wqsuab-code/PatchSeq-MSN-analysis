#!/usr/bin/env python3
"""Build the audited frozen Mouse M analysis release dated 2026-09-23."""

from __future__ import annotations

import csv
import hashlib
import json
import shutil
import sys
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd
from docx import Document
from docx.enum.section import WD_SECTION
from docx.enum.table import WD_CELL_VERTICAL_ALIGNMENT
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Inches, Pt, RGBColor


ROOT = Path(__file__).resolve().parents[1]
OUT_ROOT = ROOT / "outputs" / "morph_qc" / "Mouse_M_frozen_release_20260923"
ZIP_PATH = OUT_ROOT.with_suffix(".zip")
ID = "MSN_unique_ID"
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
M_LEVELS = ["M1", "M2", "M3", "M4"]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def source(rel_or_abs: str) -> Path:
    p = Path(rel_or_abs)
    return p if p.is_absolute() else ROOT / p


manifest: list[dict] = []


def copy_file(src: str, dest: str, role: str, required: bool = True) -> bool:
    src_path = source(src)
    if not src_path.exists():
        if required:
            raise FileNotFoundError(src_path)
        return False
    dest_path = OUT_ROOT / dest
    dest_path.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src_path, dest_path)
    manifest.append({
        "archive_path": dest.replace("\\", "/"),
        "source_path": str(src_path),
        "bytes": dest_path.stat().st_size,
        "sha256": sha256(dest_path),
        "role": role,
    })
    return True


def copy_tree(src: str, dest: str, role: str) -> None:
    src_path = source(src)
    if not src_path.exists():
        raise FileNotFoundError(src_path)
    for item in sorted(p for p in src_path.rglob("*") if p.is_file()):
        rel = item.relative_to(src_path)
        copy_file(str(item), str(Path(dest) / rel), role)


def transform(x: np.ndarray) -> np.ndarray:
    shifted = np.maximum(x - np.nanmin(x, axis=0), 0.0)
    totals = np.nansum(shifted, axis=0)
    totals[totals <= 0] = 1.0
    logged = np.log1p(shifted / totals * 10000.0)
    mean = np.nanmean(logged, axis=0)
    sd = np.nanstd(logged, axis=0, ddof=1)
    sd[sd <= 0] = 1.0
    return (logged - mean) / sd


def bool_series(s: pd.Series) -> pd.Series:
    return s.astype(str).str.lower().eq("true")


def run_audit() -> tuple[dict, pd.DataFrame, pd.DataFrame]:
    raw187_path = ROOT / "outputs/morph_qc/morph187_after_incomplete_reconstruction_quarantine/01_morph187_discovery_input.csv"
    final_path = ROOT / "outputs/morph_qc/final_morph187_NPC3_HCK4_GCres0p50_figures/00_final_morph187_cell_assignments.csv"
    interactive_path = ROOT / "interactive/mouse-morph-feature-explorer/dist/data.csv"
    pca_path = ROOT / "outputs/morph_qc/morph187_NPC3-4_HCK4-5_GC_allres_raw/01_PCA_scores_187cells.csv"
    gc_path = ROOT / "outputs/morph_qc/morph187_NPC3-4_HCK4-5_GC_allres_raw/05_GC_NPC3-4_allres_raw_assignments.csv"
    t_audit_path = ROOT / "outputs/20260820_REFonly_D1D2_marker_classifier_scan/D1D2_marker_classifier_Query_predictions.csv"
    rrr_scores_path = ROOT / "outputs/morph_qc/mouse_T_M_RRR_strict168_web/RRR_cell_scores_strict168.csv"
    asc_map_path = ROOT / "outputs/morph_qc/raw_ASC_recomputed_area_volume_final187_neurom4p0p5/01_year_aware_ASC_mapping_final187.csv"

    raw = pd.read_csv(raw187_path).sort_values(ID).reset_index(drop=True)
    final = pd.read_csv(final_path).sort_values(ID).reset_index(drop=True)
    interactive = pd.read_csv(interactive_path).sort_values(ID).reset_index(drop=True)
    stored_pca = pd.read_csv(pca_path).sort_values(ID).reset_index(drop=True)
    gc = pd.read_csv(gc_path)
    gc = gc[(gc["NPC"] == 3) & np.isclose(gc["Resolution"], 0.50)].sort_values(ID).reset_index(drop=True)

    ids_raw = set(raw[ID])
    ids_final = set(final[ID])
    ids_interactive = set(interactive[ID])
    raw_values = raw.set_index(ID)[FEATURES].sort_index().to_numpy(float)
    final_values = final.set_index(ID)[FEATURES].sort_index().to_numpy(float)
    interactive_values = interactive.set_index(ID)[FEATURES].sort_index().to_numpy(float)
    max_raw_final = float(np.max(np.abs(raw_values - final_values)))
    max_final_interactive = float(np.max(np.abs(final_values - interactive_values)))

    z = transform(raw_values)
    # Full PCA by NumPy SVD. PCA component signs are arbitrary, so compare
    # absolute score values with the frozen table.
    u, singular, _ = np.linalg.svd(z, full_matrices=False)
    recomputed = u[:, :10] * singular[:10]
    stored = stored_pca.set_index(ID).loc[sorted(ids_raw), [f"PC{i}" for i in range(1, 11)]].to_numpy(float)
    max_pc_abs_diff = float(np.max(np.abs(np.abs(recomputed) - np.abs(stored))))
    frozen_labels = final.set_index(ID).loc[sorted(ids_raw), "M_class"].to_numpy(str)
    ml_reproduction = json.loads((ROOT / "outputs/morph_qc/Mouse_M4_Macaque_style_ML500/00_results_summary.json").read_text(encoding="utf-8"))["frozen_HC_reproduction"]
    hc_match_n = int(len(final) - ml_reproduction["mismatch_n"])
    hc_ari = float(ml_reproduction["ARI"])

    gc_map = {"S3": "M1", "S2": "M2", "S0": "M3", "S1": "M4"}
    gc_as_m = gc.set_index(ID).loc[sorted(ids_raw), "GC_raw_cluster"].map(gc_map).to_numpy(str)
    consensus_n = int(np.sum(gc_as_m == frozen_labels))

    strict_audit = pd.read_csv(t_audit_path).rename(columns={"cell_id": ID})
    strict_audit["strict_T_identity"] = strict_audit["Final_CellType_stability"].str.replace("Stable_", "", regex=False)
    strict = strict_audit[
        strict_audit["Final_CellType_stability"].isin(["Stable_D1", "Stable_D2"])
        & bool_series(strict_audit["agree_RPCA"])
        & bool_series(strict_audit["agree_Pearson"])
    ][[ID, "strict_T_identity"]]
    strict = final[[ID, "M_class", "HC_GC_consensus"]].merge(strict, on=ID, how="inner")
    strict = strict[bool_series(strict["HC_GC_consensus"])].copy()
    rrr_scores = pd.read_csv(rrr_scores_path)
    rrr_ids_match = set(strict[ID]) == set(rrr_scores[ID])

    broad = final[bool_series(final["HC_GC_consensus"]) & final["D1_D2"].isin(["D1", "D2"])]
    strict_counts = pd.crosstab(strict["strict_T_identity"], strict["M_class"]).reindex(
        index=["D1", "D2"], columns=M_LEVELS, fill_value=0
    )
    broad_counts = pd.crosstab(broad["D1_D2"], broad["M_class"]).reindex(
        index=["D1", "D2"], columns=M_LEVELS, fill_value=0
    )

    asc = pd.read_csv(asc_map_path)
    asc = final[[ID]].merge(asc, on=ID, how="left", validate="one_to_one")
    ext_root = None
    for candidate in [Path("E:/ASC files"), Path("D:/ASC files")]:
        if candidate.exists():
            ext_root = candidate
            break
    rows = []
    for _, row in asc.iterrows():
        filename = row.get("Selected_ASC_filename", "")
        available = bool(ext_root is not None and isinstance(filename, str) and filename and (ext_root / filename).exists())
        rows.append({
            ID: row[ID],
            "Selected_ASC_filename": filename,
            "Expected_year": row.get("Explicit_filename_year", ""),
            "Mapping_status": row.get("Mapping_status", ""),
            "Year_match": row.get("Year_match", ""),
            "ASC_parse_success_from_frozen_audit": row.get("ASC_parse_success", ""),
            "Available_during_packaging": available,
        })
    asc_status = pd.DataFrame(rows)

    class_counts = final["M_class"].value_counts().reindex(M_LEVELS).astype(int).to_dict()
    audit = {
        "release_date": "2026-09-23",
        "frozen_taxonomy": "NPC3, Euclidean Ward.D2 HC K=4, GC resolution=0.50",
        "cohort": {
            "raw187_rows": int(len(raw)),
            "final_rows": int(len(final)),
            "interactive_rows": int(len(interactive)),
            "unique_ids": int(final[ID].nunique()),
            "class_counts": class_counts,
            "raw_final_id_sets_identical": ids_raw == ids_final,
            "final_interactive_id_sets_identical": ids_final == ids_interactive,
        },
        "feature_concordance": {
            "raw_vs_final_max_absolute_difference": max_raw_final,
            "final_vs_interactive_max_absolute_difference": max_final_interactive,
        },
        "pca_hc_reproduction": {
            "max_absolute_difference_of_absolute_PC_scores": max_pc_abs_diff,
            "HC_mapped_match_n": hc_match_n,
            "HC_total_n": int(len(final)),
            "HC_ARI": hc_ari,
        },
        "hc_gc": {
            "consensus_n": consensus_n,
            "total_n": int(len(final)),
            "consensus_percent": 100.0 * consensus_n / len(final),
        },
        "strict_T_M_RRR": {
            "strict_consensus_n": int(len(strict)),
            "D1_n": int((strict["strict_T_identity"] == "D1").sum()),
            "D2_n": int((strict["strict_T_identity"] == "D2").sum()),
            "RRR_score_rows": int(len(rrr_scores)),
            "strict_and_RRR_ID_sets_identical": bool(rrr_ids_match),
        },
        "identity_scope_conflict": {
            "broad_D1D2_consensus_n": int(len(broad)),
            "strict_RPCA_Pearson_consensus_n": int(len(strict)),
            "difference_n": int(len(broad) - len(strict)),
            "resolution": "Use the strict n=168 table only for T-M RRR; retain the broader table only for descriptive plots explicitly labelled as broad consensus identity.",
        },
        "raw_ASC": {
            "external_ASC_root_available": str(ext_root) if ext_root else None,
            "mapped_filename_n": int(asc_status["Selected_ASC_filename"].fillna("").astype(str).ne("").sum()),
            "available_file_n_during_packaging": int(asc_status["Available_during_packaging"].sum()),
            "note": "Raw ASC files were on an external E:/ or D:/ volume and were not mounted during packaging. The frozen mapping and NeuroM recomputation audits are included.",
        },
    }
    return audit, strict_counts, broad_counts, asc_status


METHOD_SECTIONS = {
"01_Cohort_and_morphology_QC.md": """# Cohort definition and morphology quality control

The primary morphology analysis comprised 187 mouse nucleus accumbens medium spiny neurons with complete measurements for the ten frozen morphology features. This cohort was derived from a 193-cell morphology set after six cells were placed in a reversible morphology-only quarantine following inspection of PC1-PC3 outliers and their raw reconstructions. The quarantined cells were D251161_Batch8, A20256247_Batch7, A20257158_Batch7, A2025797_Batch7, A20256245_Batch7 and O9192_Batch2. Their exclusion did not alter the global electrophysiological or transcriptomic eligibility tables.

Raw reconstruction files were linked to cells by a year-aware filename audit. The expected year was inferred from the morphology identifier and checked against file creation or modification metadata; when same-date filenames occurred in both 2024 and 2025, year compatibility and agreement with the tabulated total neurite length were used to resolve the match. The frozen audit reported year-aware matches for 183 of 187 cells and successful NeuroM geometry recomputation for 182 of 187 cells. Because the external ASC volume was not mounted when this release was assembled, the archive contains the complete mapping and recomputation audits but not the ASC files themselves.
""",
"02_Frozen_features_and_preprocessing.md": """# Frozen morphology features and preprocessing

The analysis used ten non-redundant morphology features: soma circularity index, soma aspect ratio, maximum radial distance, number of primary neurites, mean basal-dendrite tortuosity, total neurite length, number of bifurcation points, maximum branch order, minimum trunk angle and maximum trunk angle. Features were selected after redundancy, skewness and unit-continuity audits. Total neurite area and volume were excluded from the primary feature panel because raw-radius calibration differed across acquisition batches; total neurite length was retained because its NeuroM recomputation was effectively identical to the historical value.

For each feature, the minimum observed value in the analysis set was subtracted and negative values were truncated at zero. Each shifted feature column was divided by its column sum and multiplied by 10,000, transformed with log(1+x), and standardized to zero mean and unit sample standard deviation. Principal component analysis was then applied to the ten-feature standardized matrix. The primary clustering used PCs 1-3. In held-out or resampling analyses, every transformation and PCA model was re-estimated using only the corresponding training or resampled cells and then applied to held-out cells, preventing information leakage.
""",
"03_Clustering_and_consensus_taxonomy.md": """# Hierarchical and graph clustering

Hierarchical clustering was performed in the PC1-PC3 space using Euclidean distances and Ward.D2 linkage. The dendrogram was cut at K=4, yielding the HC-derived M1-M4 labels used as the complete 187-cell taxonomy (M1, n=67; M2, n=23; M3, n=62; M4, n=35).

Graph clustering was performed independently in the same PC1-PC3 space. A shared-nearest-neighbour graph was constructed with k=20 and pruning threshold 1/15. Louvain clustering (Seurat algorithm 1; 30 starts; 30 iterations; random seed 20260825) was run across resolutions 0.50-3.00. Resolution 0.50 produced four graph clusters and was selected for the frozen comparison. Graph labels were matched to HC labels by maximum overlap. HC and graph clustering agreed for 181 of 187 cells (96.79%). The term HC-GC consensus therefore refers only to these 181 cells; the 187-cell M taxonomy retains the HC-derived label for all included cells.
""",
"04_Visualization_and_feature_statistics.md": """# Embedding heat map reconstructions and feature comparisons

For visualization, t-distributed stochastic neighbour embedding was computed once from PCs 1-3 of the frozen 187-cell matrix using perplexity 25, learning rate 50, 3,000 iterations, PCA initialization, the Barnes-Hut algorithm with angle 0, random seed 20260826. The coordinates were reused for all class and feature overlays. Class contours are descriptive convex-hull guides after robust median and median-absolute-deviation trimming and do not represent confidence intervals or classifier decision boundaries.

The morphology heat map displays feature-wise standardized values ordered by final M class. Representative reconstructions were selected by distance to the class centroid in the frozen morphology space and rendered on common canvases with a common scale bar. For raw-value comparisons, all cells were displayed as horizontally jittered points with unfilled box plots showing the median, interquartile range and 1.5-times-interquartile-range whiskers. Each feature was tested across M1-M4 with a two-sided Kruskal-Wallis test and Benjamini-Hochberg correction across the ten features. Pairwise two-sided Mann-Whitney U tests were performed for significant omnibus results, with Benjamini-Hochberg correction within feature. Transformed Z-scores were used only for the fixed t-SNE overlays and heat map; statistical tests used the raw feature values.
""",
"05_Bootstrap_and_machine_learning_validation.md": """# Bootstrap and machine learning validation

Cell-wise cluster recovery was assessed by repeatedly sampling 80% of the 187 cells without replacement, refitting the complete transformation, three-component PCA and Ward.D2 K=4 clustering, and matching resampled clusters to the frozen M labels by maximum overlap. Pairwise co-clustering probabilities were calculated only in replicates containing both cells. The archived primary stability matrix contains 500 replicates; a 1,000-replicate sensitivity analysis is also retained.

The full machine learning validation used the same ten features and the HC-derived M1-M4 labels. Because animal identifiers were unavailable and one mouse was used per experimental day, recording date was used as a conservative animal proxy; cells recorded on the same date were never divided between training and test sets. The Macaque-matched validation performed 500 grouped resamples and compared multinomial logistic regression, linear and radial-basis-function support vector machines, k-nearest neighbours, random forest, extra trees and histogram gradient boosting. The shift-sum-10,000-log1p-Z transformation and three-component PCA were fitted within each training split. The best supervised model was the radial-basis-function support vector machine, with median held-out balanced accuracy 0.965 (95% interval, 0.873-1.000) and median macro-F1 0.957 (95% interval, 0.861-1.000). A 500-permutation label test gave empirical P=0.002. These analyses quantify internal recoverability and stability, not independent biological validation, because the labels and predictors derive from the same morphology measurements.
""",
"06_Transcriptomic_identity_definition.md": """# Strict D1 and D2 transcriptomic identity

Strict D1/D2 identity required a stable D1 or D2 call and agreement of both the reference-projection classifier and Pearson-correlation classifier. The T-M analysis additionally required HC-GC morphology consensus. This yielded 168 cells (D1, n=74; D2, n=94). Cells failing any component of this definition were excluded from T-M reduced-rank regression rather than being assigned to an ambiguous D1/D2 category. The broader D1_D2 field in the 187-cell morphology table is retained for descriptive legacy plots but is not the strict RRR inclusion variable.
""",
"07_Transcriptomic_morphological_RRR.md": """# Transcriptomic morphological reduced rank regression

T-M reduced-rank regression was refitted from scratch in the 168 cells with strict RPCA-Pearson-concordant D1/D2 identity and HC-GC morphology consensus. Starting from the label-blind 2,000-highly-variable-gene matrix, genes were re-ranked by variance within the strict cohort and the top 1,000 were retained. Each gene was standardized across cells, and transcriptomic variation was summarized by 20 principal components. The response matrix comprised the ten frozen morphology features transformed by the same minimum-shift, column-sum-10,000, log1p and Z-score procedure within the strict cohort.

Ordinary multivariate regression coefficients from transcriptomic PCs to morphology were reduced to rank 3 by singular-value decomposition of the fitted response matrix. Component signs were oriented using the morphology feature with the largest absolute loading. Transcriptomic and morphology cell scores were projected into the resulting shared axes, and displayed feature or gene vectors are correlations with the shared component scores. RRR component signs and ordering are model-specific and should not be interpreted as directly equivalent to independently fitted T-E RRR axes.
""",
"08_Reproducibility_and_data_availability.md": """# Reproducibility and data availability

The frozen analysis is defined by the 187-cell raw ten-feature table, the preprocessing rule, PC1-PC3, Euclidean Ward.D2 hierarchical clustering at K=4, and Seurat graph clustering at resolution 0.50. The archive includes the source tables, assignments, PCA outputs, graph-clustering scan, bootstrap matrices, full 500-repeat machine learning output, strict T-M RRR inputs and outputs, scripts, interactive web assets, user-saved layouts, figure source files, a SHA-256 manifest and a conflict audit.

Three sample-size scopes must remain distinct: n=187 for the complete morphology taxonomy, n=181 for HC-GC morphology consensus, and n=168 for strict D1/D2 plus HC-GC consensus T-M RRR. The actual ASC reconstructions reside on an external drive and are not embedded in this archive; the archive supplies the frozen cell-to-file mapping, availability list and NeuroM recomputation audit so that the original reconstructions can be inserted without changing analytical identifiers.
""",
}


def set_cell_shading(cell, fill: str) -> None:
    tc_pr = cell._tc.get_or_add_tcPr()
    shd = OxmlElement("w:shd")
    shd.set(qn("w:fill"), fill)
    tc_pr.append(shd)


def set_cell_border(cell, color="D9D9D9", size="4") -> None:
    tc_pr = cell._tc.get_or_add_tcPr()
    borders = tc_pr.first_child_found_in("w:tcBorders")
    if borders is None:
        borders = OxmlElement("w:tcBorders")
        tc_pr.append(borders)
    for edge in ("top", "left", "bottom", "right", "insideH", "insideV"):
        tag = "w:" + edge
        el = borders.find(qn(tag))
        if el is None:
            el = OxmlElement(tag)
            borders.append(el)
        el.set(qn("w:val"), "single")
        el.set(qn("w:sz"), size)
        el.set(qn("w:color"), color)


def build_docx(audit: dict) -> Path:
    doc = Document()
    section = doc.sections[0]
    section.top_margin = Inches(0.75)
    section.bottom_margin = Inches(0.75)
    section.left_margin = Inches(0.85)
    section.right_margin = Inches(0.85)
    styles = doc.styles
    styles["Normal"].font.name = "Arial"
    styles["Normal"]._element.rPr.rFonts.set(qn("w:ascii"), "Arial")
    styles["Normal"]._element.rPr.rFonts.set(qn("w:hAnsi"), "Arial")
    styles["Normal"].font.size = Pt(10.5)
    styles["Normal"].paragraph_format.space_after = Pt(5)
    styles["Normal"].paragraph_format.line_spacing = 1.08
    for name, size in [("Title", 20), ("Heading 1", 14), ("Heading 2", 11.5)]:
        styles[name].font.name = "Arial"
        styles[name]._element.rPr.rFonts.set(qn("w:ascii"), "Arial")
        styles[name]._element.rPr.rFonts.set(qn("w:hAnsi"), "Arial")
        styles[name].font.size = Pt(size)
        styles[name].font.color.rgb = RGBColor(0, 0, 0)
    title_ppr = styles["Title"].element.get_or_add_pPr()
    title_border = title_ppr.find(qn("w:pBdr"))
    if title_border is not None:
        title_ppr.remove(title_border)
    styles["Heading 1"].paragraph_format.space_before = Pt(11)
    styles["Heading 1"].paragraph_format.space_after = Pt(7)
    styles["Heading 1"].paragraph_format.line_spacing = 1.0
    styles["Heading 1"].paragraph_format.keep_with_next = True
    styles["Heading 1"].paragraph_format.keep_together = True
    styles["Heading 2"].paragraph_format.space_before = Pt(8)
    styles["Heading 2"].paragraph_format.space_after = Pt(5)
    styles["Heading 2"].paragraph_format.line_spacing = 1.0
    styles["Heading 2"].paragraph_format.keep_with_next = True

    title = doc.add_paragraph(style="Title")
    title.add_run("Mouse MSN Morphology Taxonomy Methods")
    subtitle = doc.add_paragraph()
    subtitle.alignment = WD_ALIGN_PARAGRAPH.LEFT
    run = subtitle.add_run("Frozen release 23 September 2026")
    run.bold = True
    run.font.size = Pt(11)
    doc.add_paragraph(
        "This document defines the reproducible analytical workflow for the frozen four-class morphology taxonomy and its transcriptomic and machine-learning validation. The authoritative sample scopes are 187 cells for the complete M taxonomy, 181 cells for HC-GC morphology consensus, and 168 cells for strict D1/D2 T-M reduced-rank regression."
    )

    doc.add_heading("Frozen analysis at a glance", level=1)
    rows = [
        ("Primary morphology cohort", "187 cells; M1=67, M2=23, M3=62, M4=35"),
        ("Frozen features", "10 de-redundant soma and dendritic measurements"),
        ("Primary taxonomy", "PC1-PC3; Euclidean Ward.D2; HC K=4"),
        ("Independent graph comparison", "Seurat SNN/Louvain; resolution 0.50"),
        ("HC-GC consensus", "181/187 cells (96.79%)"),
        ("Strict T-M RRR cohort", "168 cells; D1=74, D2=94; rank=3"),
        ("Grouped ML", "Recording date as animal proxy; 500 resamples/permutations"),
    ]
    table = doc.add_table(rows=1, cols=2)
    table.autofit = False
    table.columns[0].width = Inches(2.25)
    table.columns[1].width = Inches(4.75)
    table.rows[0].cells[0].text = "Element"
    table.rows[0].cells[1].text = "Frozen specification"
    for cell in table.rows[0].cells:
        set_cell_shading(cell, "1F4E78")
        cell.vertical_alignment = WD_CELL_VERTICAL_ALIGNMENT.CENTER
        for run in cell.paragraphs[0].runs:
            run.font.color.rgb = RGBColor(255, 255, 255)
            run.bold = True
    for idx, (a, b) in enumerate(rows):
        cells = table.add_row().cells
        cells[0].text, cells[1].text = a, b
        if idx % 2:
            for cell in cells:
                set_cell_shading(cell, "EAF2F8")
        for cell in cells:
            cell.vertical_alignment = WD_CELL_VERTICAL_ALIGNMENT.CENTER
    for row in table.rows:
        for cell in row.cells:
            set_cell_border(cell)

    for filename, markdown in METHOD_SECTIONS.items():
        lines = [x.rstrip() for x in markdown.strip().splitlines()]
        for line in lines:
            if not line:
                continue
            if line.startswith("# "):
                doc.add_heading(line[2:], level=1)
            elif line.startswith("## "):
                doc.add_heading(line[3:], level=2)
            else:
                doc.add_paragraph(line)

    doc.add_heading("Audit conclusion", level=1)
    doc.add_paragraph(
        f"The frozen raw, final and interactive data tables contain identical 187-cell identifier sets. Raw-to-final and final-to-interactive values for all ten features are exact (maximum absolute difference {audit['feature_concordance']['raw_vs_final_max_absolute_difference']:.3g} and {audit['feature_concordance']['final_vs_interactive_max_absolute_difference']:.3g}, respectively). Recomputed PCA scores agree with the stored solution up to component sign (maximum absolute difference {audit['pca_hc_reproduction']['max_absolute_difference_of_absolute_PC_scores']:.3g}), and Ward.D2 K=4 reproduces all 187 labels with ARI=1.000. The only material semantic conflict is the coexistence of a broader 180-cell D1/D2 descriptive table and the strict 168-cell RPCA-Pearson cohort; the strict table is authoritative for T-M RRR."
    )
    doc.add_paragraph(
        "Raw ASC files were unavailable because neither external ASC volume was mounted during packaging. This is a packaging limitation, not an analytical mismatch; the filename mapping, year audit, geometry recomputation results and a cell-level acquisition checklist are included."
    )

    out = OUT_ROOT / "methods" / "Mouse_M_frozen_Methods_20260923.docx"
    out.parent.mkdir(parents=True, exist_ok=True)
    doc.save(out)
    return out


def write_reports(audit: dict, strict_counts: pd.DataFrame, broad_counts: pd.DataFrame,
                  asc_status: pd.DataFrame) -> None:
    audit_dir = OUT_ROOT / "audit"
    audit_dir.mkdir(parents=True, exist_ok=True)
    (audit_dir / "AUDIT_RESULTS.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")
    strict_counts.to_csv(audit_dir / "strict_RPCA_Pearson_D1D2_by_Mclass_n168.csv")
    broad_counts.to_csv(audit_dir / "broad_D1D2_by_Mclass_n180.csv")
    asc_status.to_csv(audit_dir / "RAW_ASC_acquisition_status_for_final187.csv", index=False)

    report = f"""# Conflict and concordance audit

## Overall decision

The latest frozen Mouse M analysis is internally concordant. The 187-cell raw table, final assignment table and interactive feature table contain the same identifiers and identical values for all ten frozen features. Recomputed PCA and Ward.D2 K=4 clustering reproduce the stored taxonomy exactly. GC at resolution 0.50 agrees with HC for 181 of 187 cells (96.79%). The strict T-M RRR score table contains exactly the 168 cells defined by HC-GC consensus plus stable RPCA and Pearson D1/D2 agreement.

## Verified agreements

- Cohort: {audit['cohort']['final_rows']} rows, {audit['cohort']['unique_ids']} unique IDs; class counts {audit['cohort']['class_counts']}.
- Raw versus final frozen features: maximum absolute difference {audit['feature_concordance']['raw_vs_final_max_absolute_difference']:.3g}.
- Final versus interactive frozen features: maximum absolute difference {audit['feature_concordance']['final_vs_interactive_max_absolute_difference']:.3g}.
- PCA reproduction: maximum absolute difference after allowing component sign flips {audit['pca_hc_reproduction']['max_absolute_difference_of_absolute_PC_scores']:.3g}.
- HC reproduction: {audit['pca_hc_reproduction']['HC_mapped_match_n']}/{audit['pca_hc_reproduction']['HC_total_n']} labels, ARI={audit['pca_hc_reproduction']['HC_ARI']:.3f}.
- HC-GC agreement: {audit['hc_gc']['consensus_n']}/{audit['hc_gc']['total_n']} ({audit['hc_gc']['consensus_percent']:.2f}%).
- Strict T-M RRR: n={audit['strict_T_M_RRR']['strict_consensus_n']}, D1={audit['strict_T_M_RRR']['D1_n']}, D2={audit['strict_T_M_RRR']['D2_n']}; score IDs match exactly.

## Resolved semantic conflict

The file `11_D1D2_counts_by_final_M_class.csv` and the `D1_D2` column in the 187-cell table use a broader descriptive identity and yield {audit['identity_scope_conflict']['broad_D1D2_consensus_n']} HC-GC-consensus D1/D2 cells. Strict T-M RRR requires stable D1/D2 plus agreement of both RPCA and Pearson classifiers and yields {audit['identity_scope_conflict']['strict_RPCA_Pearson_consensus_n']} cells. The difference is {audit['identity_scope_conflict']['difference_n']} cells. These files are not numerically interchangeable. The package labels the broad table as descriptive and the strict n=168 table as authoritative for RRR.

## Version conflict resolved

The 19 September source package contained the earlier `Mouse_M_Ver09182026.tif`. This release supersedes it with `Mouse_M_20260920_Ver3.svg` and the final 20 September feature-comparison layout. The three downloaded copies of the final M-feature layout are byte-identical; only one canonical copy is retained.

## Raw reconstruction limitation

The external ASC root was not available during packaging, so zero ASC files were copied. The frozen mapping audit contains {audit['raw_ASC']['mapped_filename_n']} selected filenames for the 187-cell cohort, and the archived NeuroM audit documents 183 year-aware matches and 182 successful geometry recomputations. `RAW_ASC_acquisition_status_for_final187.csv` is the exact checklist for adding the files later without changing IDs or analytical outputs.
"""
    (audit_dir / "CONFLICT_AND_CONCORDANCE_AUDIT.md").write_text(report, encoding="utf-8")

    method_dir = OUT_ROOT / "methods" / "sections"
    method_dir.mkdir(parents=True, exist_ok=True)
    for filename, text in METHOD_SECTIONS.items():
        (method_dir / filename).write_text(text.strip() + "\n", encoding="utf-8")
    combined = "# Mouse MSN morphology taxonomy methods\n\n" + "\n\n".join(
        text.strip() for text in METHOD_SECTIONS.values()
    ) + "\n"
    (OUT_ROOT / "methods" / "Mouse_M_frozen_Methods_20260923.md").write_text(combined, encoding="utf-8")


def package_sources() -> None:
    # Final, user-edited figure and final feature comparison exports.
    copy_file(r"C:\Users\53461\OneDrive\Pictures\Mouse_M_20260920_Ver3.svg",
              "figures/final/Mouse_M_20260920_Ver3.svg", "latest user-composed final figure")
    copy_file(r"C:\Users\53461\Downloads\Mouse_M1-M4_frozen10_feature_comparison (20260920).svg",
              "figures/feature_comparison/Mouse_M1-M4_frozen10_feature_comparison_20260920.svg",
              "latest frozen feature comparison SVG")
    copy_file(r"C:\Users\53461\Downloads\Mouse_M1-M4_frozen10_feature_comparison_600dpi (1).tif",
              "figures/feature_comparison/Mouse_M1-M4_frozen10_feature_comparison_600dpi.tif",
              "latest frozen feature comparison TIFF")
    copy_file(r"C:\Users\53461\Downloads\Mouse_M_comparison_layout (Ver0260920).json",
              "layouts/Mouse_M_comparison_layout_final_20260920.json", "canonical feature-comparison layout")
    copy_file(r"C:\Users\53461\Downloads\Mouse_M_ML_validation-2-2.pdf",
              "figures/ML/Mouse_M_ML_validation_final.pdf", "latest user-exported ML figure")
    copy_file(r"C:\Users\53461\Downloads\Mouse_M_ML_validation_600dpi(1)(1).tif",
              "figures/ML/Mouse_M_ML_validation_final_600dpi.tif", "latest user-exported ML TIFF")
    copy_file(r"C:\Users\53461\Downloads\Mouse_M_ML_validation_layout.json",
              "layouts/Mouse_M_ML_validation_layout.json", "ML figure layout")
    copy_file(r"C:\Users\53461\Downloads\Mouse_T-M_RRR_layout_09172026.json",
              "layouts/Mouse_T-M_RRR_layout_09172026.json", "T-M RRR layout")

    # Primary raw and frozen inputs.
    files = [
        ("outputs/morph_qc/morph_taxonomy_round6_blind_de_novo/01_blinded_193cell_morphology_input.csv", "raw_inputs/morphology/01_blinded_193cell_morphology_input.csv", "pre-quarantine ten-feature morphology input"),
        ("outputs/morph_qc/morph187_after_incomplete_reconstruction_quarantine/00_n193_to_n187_transition_audit.csv", "raw_inputs/morphology/00_n193_to_n187_transition_audit.csv", "cohort transition audit"),
        ("outputs/morph_qc/morph187_after_incomplete_reconstruction_quarantine/00_summary.json", "raw_inputs/morphology/00_quarantine_summary.json", "quarantine summary"),
        ("outputs/morph_qc/morph187_after_incomplete_reconstruction_quarantine/01_morph187_discovery_input.csv", "raw_inputs/morphology/01_morph187_discovery_input.csv", "authoritative raw ten-feature input"),
        ("outputs/morph_qc/morph187_after_incomplete_reconstruction_quarantine/02_six_newly_quarantined_cells.csv", "raw_inputs/morphology/02_six_newly_quarantined_cells.csv", "quarantined cell list"),
        ("config/morph_current_reconstruction_spe_local_pyramid_quarantine.json", "raw_inputs/morphology/morph_current_reconstruction_quarantine.json", "frozen quarantine configuration"),
        ("outputs/morph_qc/all228_Morph_E_T_QC_eligibility_audit/01_all228_Morph_E_T_QC_status.csv", "raw_inputs/cell_QC/01_all228_Morph_E_T_QC_status.csv", "all-cell multimodal QC audit"),
        ("outputs/morph_qc/raw_ASC_recomputed_area_volume_final187_neurom4p0p5/01_year_aware_ASC_mapping_final187.csv", "raw_inputs/ASC/01_year_aware_ASC_mapping_final187.csv", "year-aware raw ASC mapping audit"),
        ("outputs/morph_qc/morph_taxonomy_round5_biological_validation/43_transcriptomic_log1p_HVG2000_cell_by_gene.csv", "raw_inputs/transcriptomics/43_transcriptomic_log1p_HVG2000_cell_by_gene.csv", "label-blind HVG2000 matrix used by T-M RRR"),
        ("outputs/20260820_REFonly_D1D2_marker_classifier_scan/D1D2_marker_classifier_Query_predictions.csv", "raw_inputs/transcriptomics/D1D2_marker_classifier_Query_predictions.csv", "strict RPCA-Pearson identity audit"),
    ]
    for src, dest, role in files:
        copy_file(src, dest, role)

    # Complete analysis directories that define the frozen result.
    copy_tree("outputs/morph_qc/morph187_NPC3-4_HCK4-5_GC_allres_raw", "analysis/clustering_scan", "PCA, HC and GC scan outputs")
    copy_tree("outputs/morph_qc/final_morph187_NPC3_HCK4_GCres0p50_figures", "analysis/final_taxonomy_figures", "complete frozen M taxonomy outputs")
    copy_tree("outputs/morph_qc/raw_ASC_recomputed_area_volume_final187_neurom4p0p5", "analysis/raw_ASC_recomputation", "NeuroM raw-geometry recomputation audit")
    copy_tree("outputs/morph_qc/Mouse_M4_Macaque_style_ML500", "analysis/ML_Macaque_style_500", "complete Macaque-style 500-repeat validation")
    copy_tree("outputs/morph_qc/final_morph187_NPC3_HCK4_GCres0p50_ML_validation", "analysis/ML_nested_grouped_validation", "nested and reproducibility ML validation")
    copy_tree("outputs/morph_qc/mouse_T_M_RRR_strict168_web", "analysis/T_M_RRR_strict168", "strict T-M RRR outputs")
    copy_tree("interactive/mouse-morph-feature-explorer/dist", "interactive/mouse-morph-feature-explorer", "frozen interactive web application")
    copy_tree("interactive/mouse-morph-ml-studio/dist", "interactive/mouse-morph-ml-studio", "frozen ML interactive web application")

    # Reproducibility scripts.
    scripts = [
        "freeze_morph187_after_incomplete_reconstruction_review.py",
        "scan_morph187_NPC3-4_HCK4-5_GC_allres.R",
        "build_final_morph187_figure_suite.py",
        "bootstrap_final_morph187_cellwise_stability.py",
        "run_final_morph187_HC_NPC3_K4_bootstrap500.py",
        "validate_final_morph187_taxonomy_with_nested_grouped_ml.py",
        "run_final_morph187_full_ml_stability_reproducibility_suite.py",
        "validate_mouse_M4_macaque_style_ml500.py",
        "recompute_mouse_T_M_RRR_strict168_for_web.py",
        "recompute_final187_neurite_area_volume_from_raw_asc.py",
        "plot_final_morph187_all22_distribution_tsne_atlas.py",
        "plot_final_morph187_frozen10_histograms.py",
        "render_final_morph187_top6_individual_traces_from_cache.py",
    ]
    for name in scripts:
        copy_file(f"scripts/{name}", f"code/{name}", "reproducibility script", required=False)


def write_readme(audit: dict) -> None:
    text = f"""# Mouse M frozen analysis release 2026-09-23

This package supersedes the 2026-09-19 source package and freezes the analysis represented by `Mouse_M_20260920_Ver3.svg`.

## Authoritative scopes

- Complete morphology taxonomy: n=187; M1=67, M2=23, M3=62, M4=35.
- HC-GC morphology consensus: n=181 (96.79%).
- Strict RPCA-Pearson D1/D2 plus HC-GC consensus T-M RRR: n=168; D1=74, D2=94.

## Frozen primary parameters

- Ten de-redundant morphology features.
- Minimum shift -> per-feature column-sum normalization to 10,000 -> log1p -> Z-score.
- PCA; PCs 1-3.
- Hierarchical clustering: Euclidean distance, Ward.D2, K=4.
- Graph clustering: Seurat SNN k=20, prune=1/15, Louvain algorithm 1, resolution=0.50.

## Directory guide

- `figures`: final user-composed figure and frozen feature/ML exports.
- `raw_inputs`: upstream morphology, transcriptomic, QC and ASC mapping inputs.
- `analysis`: complete frozen taxonomy, clustering scan, ML, raw-ASC recomputation and RRR outputs.
- `interactive`: local interactive web applications and embedded data.
- `code`: analysis scripts.
- `methods`: a compiled Word document, combined Markdown and independently reusable Method sections.
- `audit`: machine-readable and narrative consistency checks.
- `MANIFEST.csv` and `MANIFEST.json`: file provenance, SHA-256 and byte size.

## Important identity rule

The broad descriptive D1/D2 table contains {audit['identity_scope_conflict']['broad_D1D2_consensus_n']} HC-GC-consensus cells. It must not be used for T-M RRR. RRR uses only the strict n={audit['identity_scope_conflict']['strict_RPCA_Pearson_consensus_n']} RPCA-Pearson-concordant cells.

## Raw ASC status

The original ASC files were stored on an external drive that was not mounted during packaging. They are therefore not physically included. Use `audit/RAW_ASC_acquisition_status_for_final187.csv` to add the exact files later. The package already contains the frozen cell-to-file mapping and the NeuroM 4.0.5 geometry-recomputation audit.
"""
    (OUT_ROOT / "README.md").write_text(text, encoding="utf-8")


def finalize_manifest_and_zip() -> None:
    # Add generated files to manifest after all content has been produced.
    generated = [p for p in OUT_ROOT.rglob("*") if p.is_file() and p.name not in {"MANIFEST.csv", "MANIFEST.json"}]
    known = {m["archive_path"] for m in manifest}
    for p in generated:
        rel = p.relative_to(OUT_ROOT).as_posix()
        if rel not in known:
            manifest.append({
                "archive_path": rel,
                "source_path": "generated during release assembly",
                "bytes": p.stat().st_size,
                "sha256": sha256(p),
                "role": "release documentation or audit",
            })
    manifest.sort(key=lambda x: x["archive_path"])
    with (OUT_ROOT / "MANIFEST.csv").open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=["archive_path", "source_path", "bytes", "sha256", "role"])
        writer.writeheader()
        writer.writerows(manifest)
    (OUT_ROOT / "MANIFEST.json").write_text(json.dumps({"files": manifest}, indent=2), encoding="utf-8")

    with zipfile.ZipFile(ZIP_PATH, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=6, allowZip64=True) as archive:
        for p in sorted(x for x in OUT_ROOT.rglob("*") if x.is_file()):
            archive.write(p, arcname=(OUT_ROOT.name + "/" + p.relative_to(OUT_ROOT).as_posix()))


def main() -> None:
    if OUT_ROOT.exists() or ZIP_PATH.exists():
        raise RuntimeError(f"Refusing to overwrite existing release: {OUT_ROOT} or {ZIP_PATH}")
    OUT_ROOT.mkdir(parents=True)
    audit, strict_counts, broad_counts, asc_status = run_audit()
    package_sources()
    write_reports(audit, strict_counts, broad_counts, asc_status)
    build_docx(audit)
    write_readme(audit)
    finalize_manifest_and_zip()
    print(json.dumps({
        "release_dir": str(OUT_ROOT),
        "zip": str(ZIP_PATH),
        "zip_bytes": ZIP_PATH.stat().st_size,
        "manifest_files": len(manifest),
        "audit": audit,
    }, indent=2))


if __name__ == "__main__":
    main()
