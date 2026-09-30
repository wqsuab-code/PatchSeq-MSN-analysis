from __future__ import annotations

import csv
import hashlib
import json
import shutil
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, linkage
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score


ROOT = Path(r"C:\Users\53461\OneDrive\Documentos\Patch-seq_Mouse_Acb_MSN_T-type_Visualization")
OUT = ROOT / "outputs" / "Macaque_M_frozen_bundle_20260923_v1"
ZIP_OUT = ROOT / "outputs" / "Macaque_M_frozen_bundle_20260923_v1.zip"
SOURCE_ZIP = Path(r"C:\Users\53461\Downloads\Macaque-PatchSeq-BG.zip")
MAIN_SVG = Path(r"C:\Users\53461\OneDrive\Desktop\Patch-seq\4Manuscript\Macaque_M_Main.svg")
RRR_LAYOUT = Path(r"C:\Users\53461\Downloads\Macaque_T-M_RRR_layout (Macaque M1–M4 morphology radar).json")
RADAR_LAYOUT = Path(r"C:\Users\53461\Downloads\Macaque_M1-M4_radar_layout (Macaque M1–M4 morphology radar).json")

BASE = ROOT / "macaque_m"
PCA_DIR = BASE / "m18_adaptive_pca126"
FREEZE = BASE / "m18_tempfreeze_NPC5_HCK4_res2.3"
ML = ROOT / "outputs" / "Macaque_M4_full_ML500"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def cp(src: Path, rel: str) -> None:
    if not src.exists():
        raise FileNotFoundError(src)
    dst = OUT / rel
    dst.parent.mkdir(parents=True, exist_ok=True)
    if src.is_dir():
        shutil.copytree(src, dst)
    else:
        shutil.copy2(src, dst)


def write_text(rel: str, text: str) -> None:
    path = OUT / rel
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text.strip() + "\n", encoding="utf-8")


def audit() -> dict:
    manifest = json.loads((PCA_DIR / "00_manifest.json").read_text(encoding="utf-8-sig"))
    assignments = pd.read_csv(FREEZE / "01_temp_frozen_assignments_126.csv")
    raw = pd.read_csv(BASE / "m18" / "01_raw_126.csv")
    z = pd.read_csv(PCA_DIR / "02_transformed_z_117.csv")
    scores = pd.read_csv(PCA_DIR / "04_pca_scores.csv")
    merge_map = pd.read_csv(FREEZE / "02_GC13_to_GC4_merge_map.csv")
    confusion = pd.read_csv(FREEZE / "03_HC4_GC4_confusion.csv")
    tsne = pd.read_csv(FREEZE / "09_tSNE_optimized_coordinates.csv")
    rrr = pd.read_csv(FREEZE / "RRR_T_M" / "RRR_cell_scores_n117.csv")
    morph_load = pd.read_csv(FREEZE / "RRR_T_M" / "RRR_Mfeature_correlation_loadings.csv", index_col=0)
    radar = json.loads(RADAR_LAYOUT.read_text(encoding="utf-8-sig"))
    rrr_layout = json.loads(RRR_LAYOUT.read_text(encoding="utf-8-sig"))
    ml_params = json.loads((ML / "00_parameters.json").read_text(encoding="utf-8-sig"))

    features = manifest["features"]
    ids = set(assignments.cell_label)
    consensus = assignments[assignments.concordant.astype(str).str.lower().eq("true")].copy()
    consensus_ids = set(consensus.cell_label)
    map_dict = dict(zip(merge_map.GC_raw.astype(int), merge_map.GC_merged.astype(int)))
    mapped = assignments.GC_raw_K13.astype(int).map(map_dict)

    p = scores.set_index("cell_label").loc[assignments.cell_label]
    recalculated_hc = fcluster(
        linkage(p[[f"PC{i}" for i in range(1, 6)]].to_numpy(), method="ward", metric="euclidean"),
        4,
        criterion="maxclust",
    )
    hc_ari = adjusted_rand_score(assignments.HC_K4, recalculated_hc)
    hc_nmi = normalized_mutual_info_score(assignments.HC_K4, recalculated_hc)

    ct = pd.crosstab(assignments.HC_K4, assignments.GC_merged_K4).reindex(index=range(1, 5), columns=range(1, 5), fill_value=0)
    stored_ct = confusion.pivot(index="HC", columns="Merged_GC", values="N").reindex(index=range(1, 5), columns=range(1, 5), fill_value=0)

    svg = MAIN_SVG.read_text(encoding="utf-8-sig")
    colors = radar["colors"]
    class_counts = {f"M{k}": int((consensus.HC_K4 == k).sum()) for k in range(1, 5)}
    radar_counts = {k: int(v["title"].split("n=")[1].rstrip(")")) for k, v in radar["classes"].items()}

    results = {
        "authoritative_scope": "Macaque MSN only; Ca, Pu and NAC; 18 complete dendrite/soma features",
        "source_zip_sha256_observed": sha256(SOURCE_ZIP),
        "source_zip_sha256_expected": manifest["source_zip_sha256"],
        "source_zip_hash_match": sha256(SOURCE_ZIP).lower() == manifest["source_zip_sha256"].lower(),
        "n_complete": int(len(assignments)),
        "n_complete_donors": int(assignments.donor_label.nunique()),
        "n_consensus": int(len(consensus)),
        "n_consensus_donors": int(consensus.donor_label.nunique()),
        "class_counts": class_counts,
        "roi_counts": {str(k): int(v) for k, v in consensus.ROI.value_counts().sort_index().items()},
        "T_class_counts": {str(k): int(v) for k, v in consensus.Subclass.value_counts().sort_index().items()},
        "agreement": float(assignments.concordant.mean()),
        "adjusted_rand_HC_vs_merged_GC": float(adjusted_rand_score(assignments.HC_K4, assignments.GC_merged_K4)),
        "recomputed_HC_ARI": float(hc_ari),
        "recomputed_HC_NMI": float(hc_nmi),
        "raw_shape": list(raw.shape),
        "raw_missing_values": int(raw.isna().sum().sum()),
        "z_file_shape": list(z.shape),
        "z_missing_values": int(z.isna().sum().sum()),
        "PCA_score_shape": list(scores.shape),
        "feature_count": len(features),
        "feature_set_raw_match": set(features).issubset(raw.columns),
        "feature_set_z_match": set(features) == set(z.columns) - {"cell_label"},
        "feature_set_RRR_match": set(features) == set(morph_load.index.astype(str)),
        "IDs_raw_match_assignments": set(raw.cell_label) == ids,
        "IDs_z_match_assignments": set(z.cell_label) == ids,
        "IDs_PCA_match_assignments": set(scores.cell_label) == ids,
        "IDs_tSNE_match_assignments": set(tsne.cell_label) == ids,
        "IDs_RRR_match_consensus": set(rrr.cell_label) == consensus_ids,
        "stored_GC_merge_map_exact": bool((mapped.to_numpy() == assignments.GC_merged_K4.to_numpy()).all()),
        "stored_confusion_exact": bool(np.array_equal(ct.to_numpy(), stored_ct.to_numpy())),
        "radar_feature_count": len(radar["featureOrder"]),
        "radar_features_subset_of_18": set(radar["featureOrder"]).issubset(features),
        "radar_counts_match": radar_counts == class_counts,
        "radar_colors_present_in_main_SVG": {k: v.upper() in svg.upper() for k, v in colors.items()},
        "main_SVG_panel_order": {"c": "HC-GC confusion matrix", "d": "consensus clustering"},
        "RRR_rank_display": 3,
        "RRR_layout_radial_percentile": rrr_layout.get("radialPercentile"),
        "ML_class_counts_match": ml_params.get("class_counts") == class_counts,
        "ML_feature_count_match": ml_params.get("feature_n") == len(features),
        "Mouse_morphology_used": False,
        "Macaque_E_features_used": False,
    }
    return results


METHODS = {
"01_cohort_and_QC_EN.md": r'''# Macaque MSN cohort and morphology quality control

Analyses were restricted to macaque medium spiny neurons annotated as STR D1 MSN, STR D2 MSN or STR Hybrid MSN and sampled from the caudate nucleus, putamen or nucleus accumbens. Cell metadata and morphology measurements were read directly from `Macaque-PatchSeq-BG.zip`, whose SHA-256 checksum was verified before analysis. Mouse morphology results and macaque electrophysiology features were not used.

The source metadata contained 486 eligible MSN profiles. Morphology classification required complete observations across 18 nonredundant dendritic and somatic measurements. Axonal variables, axonal Sholl measurements, the circular variable `axon_exit_theta_coronal`, and the derived `3_Sholl_PC1` score were excluded. A total of 126 cells from 42 donors met the complete-case criterion. No imputation was performed and no cell was removed solely because of an extreme PCA score.
''',
"02_preprocessing_and_PCA_EN.md": r'''# Morphology transformation and principal component analysis

The 18 retained morphology features were transformed independently. Adjusted sample skewness was calculated on the 126-cell cohort. Features with an absolute adjusted skewness of at least 0.5 underwent maximum-likelihood Yeo–Johnson transformation, with the transformation parameter optimized over −5 to 5; the remaining variables were left on their original scale. Each resulting feature was centered and divided by its sample standard deviation. Principal component analysis was then performed by singular-value decomposition without further centering or scaling. The frozen classification used PC1–PC5.

The transformation and PCA were fitted to all 126 morphology-complete cells because they were part of the unsupervised discovery procedure. In donor-held-out machine-learning validation, transformation, standardization and PCA parameters were instead estimated within each training split and applied unchanged to the corresponding test cells.
''',
"03_consensus_clustering_EN.md": r'''# Hierarchical and graph-based consensus classification

Hierarchical clustering was performed in PC1–PC5 space using Euclidean distance and Ward's minimum-variance linkage (Ward.D2), and the dendrogram was cut at four groups. In parallel, a shared-nearest-neighbour graph was constructed with 20 nearest neighbours and an SNN pruning threshold of 1/15. Louvain community detection used algorithm 1, resolution 2.3 and random seed 777, yielding 13 graph communities. Each graph community was mapped to the hierarchical cluster containing the largest number of its cells, producing four merged graph classes.

Cells with identical hierarchical and merged graph assignments were defined as the frozen consensus set. Agreement was observed for 117 of 126 cells (92.86%; adjusted Rand index 0.825), yielding M1=43, M2=42, M3=20 and M4=12. Nine discordant cells were retained in displays where indicated but excluded from consensus-class summaries and supervised label-recovery analyses. T class, anatomical region and donor identity were not used to construct the M labels.
''',
"04_visualization_EN.md": r'''# Low-dimensional displays, dendrogram, radar plots and heat map

For visualization, exact t-SNE was applied to PC1–PC5. Perplexities of 10, 15, 20, 25 and 30 and random seeds 777–781 were compared using 10-nearest-neighbour preservation and hierarchical-class silhouette width. The selected embedding used perplexity 25, seed 777 and 2,000 iterations. For the displayed compact version only, within-class deviations from each hierarchical-cluster centroid were multiplied by 0.85; clustering, labels and inferential results were unchanged. Dashed t-SNE boundaries are 80% bivariate-normal regions fitted to consensus cells.

The circular dendrogram displays the complete Ward.D2 tree for all 126 cells. Leaf order was taken directly from the tree. The inner ring shows hierarchical labels and the outer ring shows merged graph labels.

Radar plots use ten prespecified features drawn from the frozen 18-feature standardized matrix. The final feature order and labels are stored in the accompanying radar-layout JSON. All classes share the same cohort-wide Z-score scale from −3 to 3; values outside this range are clipped only for display. Pale curves show individual cells and black closed curves show class medians. No within-class normalization or interior fill was used.

The heat map shows all 18 frozen feature Z-scores for the 117 consensus cells, clipped to −2.5 to 2.5 for display. Cells were grouped by M class and shuffled within class using a fixed seed. Features were assigned to the class with the largest mean Z-score and ordered within each module by Ward.D2 clustering. Annotation tracks show M class, T class and region.
''',
"05_T_class_and_RRR_EN.md": r'''# T-class composition and transcriptome-to-morphology reduced-rank regression

Transcriptomic D1, D2 and hybrid identities were treated as post hoc annotations and were not used to define morphology classes. The 117-cell consensus cohort contained 47 D1, 62 D2 and 8 hybrid MSNs. T-class composition was summarized as counts and within-M-class percentages.

Transcriptome-to-morphology coupling was examined by pooled reduced-rank regression (RRR) in the 117 consensus cells. Transcriptomic predictors were the first 20 principal components derived from 1,000 eligible variable genes. Responses were the 18 frozen transformed and standardized morphology features; `3_Sholl_PC1` was excluded. A rank-three model was used for display. Panels show components 2 versus 1 and 3 versus 1 in transcriptomic and morphological score space. D1 and D2 cells were assigned distinct colors and 90% bivariate-normal ellipses; hybrid cells were retained as grey points. The ten genes or morphology features with the largest correlation-loading magnitudes were shown in each component plane, with vectors uniformly rescaled for display.
''',
"06_machine_learning_validation_EN.md": r'''# Donor-grouped machine-learning validation

The frozen 117-cell M1–M4 consensus set was evaluated using 500 repeated donor-grouped analyses. Donors were kept intact so that cells from the same donor could not occur in both training and test sets. In every split or resample, feature transformation, standardization, PCA, feature selection and model fitting were estimated using training data only. Complete-case analysis was maintained without imputation.

Unsupervised discovery stability was assessed by repeatedly retaining approximately 80% of donors, refitting preprocessing and PCA, and comparing de novo Ward, k-means, Gaussian-mixture and spectral solutions with the frozen labels after optimal label matching. The scan covered 2–10 PCs and K=2–8. Supervised recoverability was evaluated on repeated 75%/25% donor-grouped partitions that contained all four classes in both subsets. Balanced accuracy and macro-F1 were the primary metrics. Label-permutation controls preserved donor structure. Cell-level stability was summarized from pairwise co-clustering probabilities, and feature robustness was evaluated using training-fold permutation importance, Extra Trees impurity importance, multinomial-logistic coefficients and progressive feature ablation. These analyses measure internal recoverability and resampling stability, not independent biological validation.
''',
"07_statistical_analysis_EN.md": r'''# Statistical analysis and interpretation

The main clustering figure is descriptive. HC–GC agreement was summarized by the proportion of identical assignments and the adjusted Rand index. T-class and regional compositions are reported as cell counts and percentages. Radar plots show cell-level profiles and class medians; heat-map values are standardized feature measurements. No null-hypothesis test was used to define M1–M4.

For the separate 18-feature class-comparison analysis, each feature was compared across M1–M4 using a two-sided Kruskal–Wallis test. The 18 omnibus P values were adjusted using the Benjamini–Hochberg procedure. Post hoc pairwise comparisons used Dunn tests with Holm correction. Cells were the displayed observational units; donor-grouped resampling was used for machine-learning validation to address within-donor dependence.
''',
}


METHODS_CN = {
"01_cohort_and_QC_CN.md": """# Macaque MSN队列与形态QC\n\n分析仅纳入Ca、Pu和NAc中的STR D1 MSN、STR D2 MSN及STR Hybrid MSN。数据直接读取自经SHA-256校验的Macaque-PatchSeq-BG.zip，未调用Mouse形态结果或Macaque电生理特征。486个元数据MSN中，126个细胞（42个供体）在18项非冗余树突/胞体指标上完整。轴突、轴突Sholl、圆周角变量axon_exit_theta_coronal及派生变量3_Sholl_PC1均未进入分类。采用完整病例分析，不插补；极端PCA分数不作为自动删除依据。""",
"02_preprocessing_and_PCA_CN.md": """# 数值变换与PCA\n\n在126细胞队列中逐指标计算校正样本偏度。绝对偏度不小于0.5的指标采用最大似然Yeo–Johnson变换（lambda在−5至5内优化），其余指标保留原尺度；随后逐指标中心化并除以样本标准差。对标准化矩阵执行PCA，冻结分类使用PC1–PC5。供体隔离的机器学习验证中，变换、标准化及PCA均仅在训练折内估计，再原样应用于测试折。""",
"03_consensus_clustering_CN.md": """# HC–GC共识分类\n\nHC在PC1–PC5空间使用欧氏距离和Ward.D2，并切为K=4。GC使用k=20、prune.SNN=1/15的SNN图及Louvain algorithm 1（resolution=2.3，seed=777），得到13个原始社区；每个社区按与HC的最大细胞重叠映射至4类。HC与合并GC标签一致的117/126个细胞构成冻结共识集（92.86%，ARI=0.825）：M1=43、M2=42、M3=20、M4=12。9个不一致细胞不进入共识类汇总和监督恢复性分析。T class、ROI和供体均未参与M标签构建。""",
"04_visualization_CN.md": """# 可视化\n\nt-SNE基于PC1–PC5，对perplexity 10、15、20、25、30及seed 777–781进行扫描，按10近邻保持率和HC silhouette选择perplexity=25、seed=777。主图紧凑显示仅将类内坐标相对类中心缩放为0.85，不改变聚类或标签。虚线为共识细胞拟合的80%二维正态区域。环形树保留完整Ward.D2拓扑及原始叶序，内环为HC、外环为合并GC。雷达图读取冻结Z-score矩阵中的10项指标，所有类别共享−3至3尺度，浅色线为单细胞、黑线为中位数，无组内归一化和内部填色。热图显示117个共识细胞的全部18项Z-score，显示范围为−2.5至2.5；列按M类分组并在组内固定种子打散，行按最高组均值形成模块。""",
"05_T_class_and_RRR_CN.md": """# T class组成与RRR\n\nD1、D2和Hybrid仅作为事后注释。117个共识细胞包括D1=47、D2=62、Hybrid=8。RRR使用1,000个合格变异基因得到的前20个转录组PC预测18项冻结形态特征，采用rank=3展示。3_Sholl_PC1未纳入。D1和D2分别配色并绘制90%二维正态椭圆，Hybrid保留为灰色点；每个分量平面显示相关载荷绝对值最高的10个基因或形态特征。""",
"06_machine_learning_validation_CN.md": """# 机器学习验证\n\n对117个冻结共识细胞进行500次供体分组验证，同一供体不跨训练集和测试集。每次均仅在训练数据内估计变换、标准化、PCA、特征选择和模型参数。无监督稳定性在约80%供体重抽样内重建Ward、k-means、GMM和spectral分类，并扫描NPC=2–10、K=2–8。监督恢复性采用约75%/25%供体训练/测试划分，以balanced accuracy和macro-F1为主要指标；置换检验保留供体结构。细胞共聚类概率用于计算个体稳定性，特征稳健性由训练折内置换重要性、Extra Trees MDI、逻辑回归系数及递进消融评估。这些结果反映内部可恢复性和重抽样稳定性，不等同于独立生物学验证。""",
"07_statistical_analysis_CN.md": """# 统计分析\n\n主分类图为描述性展示；HC–GC一致性以一致率和ARI汇总，T class及ROI以计数和百分比呈现。M1–M4并非由显著性检验定义。独立的18指标组间比较使用双侧Kruskal–Wallis检验，并对18个总体P值进行Benjamini–Hochberg校正；事后两两比较使用Dunn检验及Holm校正。机器学习采用供体分组重抽样，以处理同一供体内细胞的依赖性。""",
}


FIGURE_LEGEND = r'''# Figure legend

**Fig. 6 | Consensus morphological classes of macaque striatal medium spiny neurons.** Morphological classification was performed on 126 macaque medium spiny neurons (MSNs) from the caudate nucleus, putamen and nucleus accumbens with complete measurements for 18 nonredundant dendritic and somatic features. Four classes were defined by agreement between hierarchical clustering (HC) and graph-based clustering (GC).

**a,** t-SNE representation of all 126 cells, colored by merged GC assignment. Dashed contours denote 80% bivariate-normal regions estimated from HC–GC consensus cells. **b,** Circular HC dendrogram with HC and merged GC assignments shown by the inner and outer rings, respectively. Leaves retain the order determined by the dendrogram. **c,** HC-by-GC confusion matrix. The two methods agreed for 117 of 126 cells (92.86%; adjusted Rand index, 0.825). **d,** The 117 concordantly assigned cells formed four consensus morphology classes: M1 (n=43), M2 (n=42), M3 (n=20) and M4 (n=12). The nine discordant cells are shown in dark grey. Dashed contours denote 80% bivariate-normal regions.

**e,** Radar plots of ten representative morphology features. Values are cohort-wide Z-scores displayed from −3 to +3. Pale lines represent individual cells and black outlines indicate class medians; no within-class normalization or area fill was applied. **f,** Heat map of all 18 morphology feature Z-scores in the 117 consensus cells, displayed from −2.5 (magenta) through 0 (black) to +2.5 (yellow). Cells are grouped by M class and randomly ordered within each class; features are arranged into class-associated modules. Annotation tracks indicate M class, transcriptomic T class and anatomical region. The consensus cohort comprised 63 caudate, 41 putamen and 13 nucleus accumbens cells and included 47 D1, 62 D2 and 8 hybrid MSNs.

**g,** T-class composition of each M class. D1, D2 and hybrid counts were 18, 20 and 5 in M1; 18, 23 and 1 in M2; 6, 13 and 1 in M3; and 5, 6 and 1 in M4, respectively. Bar heights show percentages and numbers within bars show cell counts. **h,** Distribution of D1 (n=47) and D2 (n=62) cells in the fixed morphology t-SNE space. Cells are colored by M class; other consensus cells are light grey and HC–GC-discordant cells are dark grey. Dashed contours reproduce the consensus-class regions in **d**.

**i,** Reduced-rank regression (RRR) relating transcriptomic variation to morphology in the 117 consensus cells. Predictors were transcriptomic principal components derived from 1,000 variable genes, and responses were the 18 morphology features. Panels show components 2 versus 1 and 3 versus 1 for the transcriptomic and morphological RRR scores. D1 and D2 cells are colored separately, whereas hybrid cells are retained as grey points. Dashed ellipses indicate 90% bivariate-normal regions. Arrows show the ten genes or morphology features with the largest correlation loadings in each component plane.

All panels are descriptive and no null-hypothesis significance tests were performed. n denotes cells. Transcriptomic identity and anatomical region were not used to define M1–M4. t-SNE, t-distributed stochastic neighbour embedding; RRR, reduced-rank regression.
'''


def main() -> None:
    if OUT.exists() or ZIP_OUT.exists():
        raise RuntimeError(f"Refusing to overwrite existing bundle: {OUT} or {ZIP_OUT}")
    OUT.mkdir(parents=True)

    # Immutable source and final manuscript/layout files.
    cp(SOURCE_ZIP, "01_source/Macaque-PatchSeq-BG.zip")
    cp(MAIN_SVG, "03_main_figure/Macaque_M_Main.svg")
    cp(RADAR_LAYOUT, "03_main_figure/layouts/Macaque_M1-M4_radar_layout_FINAL.json")
    cp(RRR_LAYOUT, "03_main_figure/layouts/Macaque_T-M_RRR_layout_FINAL.json")

    # Frozen cell-level inputs and derived matrices.
    cp(BASE / "m18" / "01_raw_126.csv", "02_frozen_input/01_raw_morphology_complete_126.csv")
    for name in ["00_manifest.json", "01_transform_skewness_audit.csv", "02_transformed_z_117.csv",
                 "03_pca_variance.csv", "04_pca_scores.csv", "05_pca_loadings.csv",
                 "06_PC1_PC5_outlier_diagnostics.csv", "07_new_PC1_PC5_3SD_candidates.csv"]:
        target = "02_transformed_and_PCA/02_transformed_z_all126.csv" if name == "02_transformed_z_117.csv" else f"02_transformed_and_PCA/{name}"
        cp(PCA_DIR / name, target)
    for name in ["01_temp_frozen_assignments_126.csv", "02_GC13_to_GC4_merge_map.csv",
                 "03_HC4_GC4_confusion.csv", "08_tSNE_parameter_scan.csv",
                 "09_tSNE_optimized_coordinates.csv", "12_tSNE_display_contracted_coordinates.csv"]:
        cp(FREEZE / name, f"02_frozen_classification/{name}")

    # Main-figure source data and rendered components.
    cp(ROOT / "outputs" / "Macaque_M4_circular_dendrogram", "03_main_figure/circular_dendrogram")
    cp(FREEZE / "tSNE_consensus_GC_Tclass", "03_main_figure/tSNE_panels")
    cp(FREEZE / "panels_A_I" / "I_morphology_Zscore_heatmap.pdf", "03_main_figure/heatmap/I_morphology_Zscore_heatmap.pdf")
    cp(FREEZE / "panels_A_I" / "I_morphology_Zscore_heatmap.png", "03_main_figure/heatmap/I_morphology_Zscore_heatmap.png")
    cp(FREEZE / "panels_A_I" / "I_feature_peak_consensus_group.csv", "03_main_figure/heatmap/I_feature_peak_consensus_group.csv")
    cp(FREEZE / "panels_A_I" / "I_morphology_Zscore_heatmap_metric_names_only.txt", "03_main_figure/heatmap/feature_order.txt")
    cp(FREEZE / "panels_A_I" / "Zscore_top10_radar", "03_main_figure/radar_source")
    cp(FREEZE / "panels_A_I" / "K_consensusM4_Tclass_counts.csv", "03_main_figure/T_class_composition/Tclass_by_Mclass_counts.csv")
    for name in ["K_consensusM4_Tclass_stacked_percent_labeled.pdf", "K_consensusM4_Tclass_stacked_percent_labeled.png"]:
        cp(FREEZE / "panels_A_I" / name, f"03_main_figure/T_class_composition/{name}")
    cp(FREEZE / "RRR_T_M", "03_main_figure/RRR_T_M")

    # Feature-comparison statistics and full donor-grouped validation.
    for name in ["all18_plot_data_117.csv", "all18_KruskalWallis_BH.csv", "all18_Dunn_Holm.csv",
                 "Macaque_Morphology_ConsensusM4_all18_3perrow_W7p2_H8p76.pdf",
                 "Macaque_Morphology_ConsensusM4_all18_3perrow_W7p2_H8p76.png"]:
        cp(FREEZE / "M4_compare_all18" / name, f"04_statistics/feature_comparisons/{name}")
    cp(ML, "04_statistics/Macaque_M4_full_ML500")

    # Reproducible scripts directly used by this frozen chain.
    scripts = [
        BASE / "run_m18_adaptive_transform_pca.R", BASE / "freeze_npcs5_hck4_res23.R",
        BASE / "optimize_tempfreeze_k4_tsne.R", BASE / "plot_consensus_gc_tclass_tsne.R",
        BASE / "plot_consensusM4_Tclass_stacked.R", BASE / "make_panels_A_to_I.R",
        BASE / "plot_frozen_m_radar_top10_zscore.py", BASE / "plot_frozen_T_M_RRR_fourpanel.py",
        ROOT / "scripts" / "plot_macaque_M4_circular_dendrogram.R",
        ROOT / "scripts" / "validate_macaque_M4_full_ml500.py",
        ROOT / "scripts" / "reproduce_macaque_M4_frozen_gc.R",
        ROOT / "scripts" / "build_macaque_m_radar_editor.py",
        BASE / "build_interactive_macaque_m_explorer.py",
    ]
    for script in scripts:
        cp(script, f"05_scripts/{script.name}")
    cp(Path(__file__), "05_scripts/package_frozen_main_bundle.py")

    # Methods and figure legend.
    for name, text in METHODS.items():
        write_text(f"06_methods/EN/{name}", text)
    for name, text in METHODS_CN.items():
        write_text(f"06_methods/CN/{name}", text)
    all_en = "# Complete methods: frozen Macaque M1–M4 morphology analysis\n\n" + "\n\n".join(METHODS.values())
    all_cn = "# 完整方法：冻结Macaque M1–M4形态分析\n\n" + "\n\n".join(METHODS_CN.values())
    write_text("06_methods/Methods_ALL_EN.md", all_en)
    write_text("06_methods/Methods_ALL_CN.md", all_cn)
    write_text("06_methods/Figure_legend_EN.md", FIGURE_LEGEND)

    results = audit()
    write_text("07_audit/frozen_consistency_audit.json", json.dumps(results, ensure_ascii=False, indent=2))
    conflicts = [
        ["Legacy transformed matrix filename", "02_transformed_z_117.csv contains 126 rows", "Renamed in bundle to 02_transformed_z_all126.csv; consensus panels subset it by the 117 frozen IDs", "Resolved; content is correct"],
        ["Radar feature order", "Earlier top10 parameter text uses a different circular order from the final radar-layout JSON", "Final radar-layout JSON and current main SVG are authoritative", "Display-only; same ten features and same values"],
        ["Radar physical width", "Earlier static specification records 4.64 in; final layout requests 5.2 in", "Final layout JSON and manuscript SVG are authoritative", "Display-only"],
        ["RRR radial scaling", "Early README records 99th-percentile radius; final layout JSON records radialPercentile=100", "Final layout JSON is authoritative for the current main figure", "Display-only; RRR scores and loadings unchanged"],
        ["Legacy preprocessing branch", "Older m18 directory contains shift-sum-log1p outputs", "Frozen labels demonstrably use m18_adaptive_pca126: skew-adaptive Yeo-Johnson/direct Z-score and PC1-PC5", "Legacy branch excluded from bundle except raw126 source table"],
        ["Panel order", "Earlier drafts placed consensus t-SNE before the confusion matrix", "Current SVG is c=confusion matrix and d=consensus t-SNE", "Resolved in bundled legend"],
    ]
    conflict_path = OUT / "07_audit" / "conflict_resolution.tsv"
    conflict_path.parent.mkdir(parents=True, exist_ok=True)
    with conflict_path.open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["item", "observed_difference", "authoritative_resolution", "impact"])
        w.writerows(conflicts)

    readme = f'''# Frozen Macaque M1–M4 analysis bundle (2026-09-23)

This bundle contains the source archive, frozen 126-cell morphology input, 117-cell consensus labels, transformation/PCA outputs, HC and merged-GC assignments, visualization source data, RRR data, feature-comparison statistics, full 500-repeat donor-grouped machine-learning validation, final manuscript SVG, final layout JSON files, analysis scripts, split English/Chinese Methods, and consistency audits.

Authoritative frozen result:
- complete morphology cohort: 126 cells from 42 donors;
- consensus cohort: 117 cells from 41 donors;
- M1=43, M2=42, M3=20, M4=12;
- HC–GC agreement=92.86%, ARI=0.825;
- 18 dendrite/soma features, adaptive Yeo–Johnson/direct Z-score, PC1–PC5;
- HC: Euclidean Ward.D2, K=4;
- GC: SNN/Louvain resolution 2.3, raw K=13 merged to K=4.

The audit identifies display-only discrepancies among legacy plotting files and records which final file is authoritative. No Mouse morphology or Macaque E-class features are included as predictors.

Source ZIP SHA-256: {results['source_zip_sha256_observed']}
'''
    write_text("00_README.md", readme)

    # Hash manifest is generated last and covers every other file in the unpacked bundle.
    rows = []
    for path in sorted(p for p in OUT.rglob("*") if p.is_file()):
        rows.append([path.relative_to(OUT).as_posix(), path.stat().st_size, sha256(path)])
    manifest_path = OUT / "07_audit" / "file_manifest_sha256.csv"
    with manifest_path.open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow(["relative_path", "size_bytes", "sha256"])
        w.writerows(rows)

    # Create the final archive without deleting or changing any source file.
    with zipfile.ZipFile(ZIP_OUT, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=6) as zf:
        for path in sorted(p for p in OUT.rglob("*") if p.is_file()):
            zf.write(path, arcname=(OUT.name + "/" + path.relative_to(OUT).as_posix()))

    print(json.dumps({
        "bundle_dir": str(OUT),
        "zip": str(ZIP_OUT),
        "zip_size_bytes": ZIP_OUT.stat().st_size,
        "zip_sha256": sha256(ZIP_OUT),
        "audit": results,
    }, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
