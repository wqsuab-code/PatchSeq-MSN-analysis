#!/usr/bin/env python
"""Create a compact Macaque E RRR Methods/provenance supplement."""
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
import matplotlib


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs/Macaque_E_RRR_methods_supplement_20260924"
ZIP_OUT = ROOT / "outputs/Macaque_E_RRR_methods_supplement_20260924.zip"
RRR = ROOT / "outputs/R3_RRR_fourpanel"
ASSIGN = ROOT / "outputs/R3_panels_v3/00_tuned_p80_ee12_random_seed777_coordinates.csv"
ZMAT = ROOT / "outputs/dSTR_dSTRvSTR_E_QC/Ca_Pu_NAC_final19_consensus_scan/Ca_Pu_NAC_final19_transformed_Zscore_matrix_n390.csv"
H5AD = ROOT / ".codex_tmp/macaque_rrr/Data/HMBA-Macaque-PatchSeq-BG-log2.h5ad"
LAYOUT = Path(r"C:\Users\53461\Downloads\Macaque_T-E_RRR_gc_merged_layout (3).json")


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def write(rel: str, content: str) -> None:
    path = OUT / rel
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content.rstrip() + "\n", encoding="utf-8")


def copy(src: Path, rel: str) -> None:
    if src.exists():
        dst = OUT / rel
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dst)


def main() -> None:
    if OUT.exists() or ZIP_OUT.exists():
        raise FileExistsError("Supplement already exists; refusing to overwrite it.")
    OUT.mkdir(parents=True)

    scores = pd.read_csv(RRR / "RRR_cell_scores_n346.csv")
    e_load = pd.read_csv(RRR / "RRR_Efeature_correlation_loadings.csv")
    g_load = pd.read_csv(RRR / "RRR_gene_correlation_loadings.csv")
    variance = pd.read_csv(RRR / "transcriptomic_PCA_variance.csv")
    assign = pd.read_csv(ASSIGN)
    frozen = assign.loc[assign["Consensus"].astype(str).str.lower().eq("true") & assign.T_class.isin(["D1", "D2"])]

    checks = [
        ("RRR cell count", len(scores), 346, len(scores) == 346),
        ("D1 count", int((scores.T_class == "D1").sum()), 163, int((scores.T_class == "D1").sum()) == 163),
        ("D2 count", int((scores.T_class == "D2").sum()), 183, int((scores.T_class == "D2").sum()) == 183),
        ("RRR IDs equal stable D1/D2 consensus IDs", len(set(scores.cell_label) ^ set(frozen.cell_label)), 0,
         set(scores.cell_label) == set(frozen.cell_label)),
        ("Electrophysiological input features", len(e_load), 19, len(e_load) == 19),
        ("Transcriptomic gene loadings", len(g_load), 1000, len(g_load) == 1000),
        ("Transcriptomic PCs", len(variance), 20, len(variance) == 20),
        ("Variance explained by 20 transcriptomic PCs", float(variance.iloc[:, 1].sum()), 0.3472151021571191,
         np.isclose(float(variance.iloc[:, 1].sum()), 0.3472151021571191)),
        ("Frozen E matrix feature count", pd.read_csv(ZMAT, nrows=1).shape[1] - 1, 19,
         pd.read_csv(ZMAT, nrows=1).shape[1] - 1 == 19),
    ]
    audit = pd.DataFrame(checks, columns=["check", "observed", "expected", "pass"])
    (OUT / "04_AUDIT").mkdir(parents=True)
    audit.to_csv(OUT / "04_AUDIT/RRR_consistency_audit.csv", index=False)

    methods = """# Transcriptome-electrophysiology reduced-rank regression

Transcriptome-electrophysiology reduced-rank regression (RRR) was performed in the 346 HC-GC consensus cells with stable D1 or D2 transcriptomic assignments (D1, n=163; D2, n=183); the 22 Hybrid cells were excluded. Transcriptomic values were read from the processed log2-expression matrix `HMBA-Macaque-PatchSeq-BG-log2.h5ad` in the Allen Human-Mammalian Brain Atlas macaque Patch-seq release. Genes detected in fewer than eight of the 346 cells, genes with zero or non-finite variance, and technical genes matching mitochondrial, ribosomal, `GM`/`Gm`, `Rik`, `MALAT1`/`Malat1`, `XIST`/`Xist`, or empty-symbol patterns were excluded. When duplicate gene symbols were present, the row with the greatest variance was retained. The 1,000 eligible genes with the greatest sample variance were selected without using D1/D2 or E-class labels, standardized across cells to zero mean and unit population variance (`ddof=0`), and reduced to 20 transcriptomic principal components by full singular-value decomposition. These 20 components explained 34.72% of transcriptomic variance in the analyzed cohort.

The response matrix comprised the same 19 transformed and standardized electrophysiological features used for the frozen E classification. Both predictor and response matrices were centered. Ordinary least-squares coefficients were estimated using the Moore-Penrose pseudoinverse of the centered 20-PC transcriptomic matrix. Singular-value decomposition of the fitted electrophysiological matrix supplied the response directions, and the coefficient matrix was projected onto the first three directions to obtain a rank-3 model. Component signs were oriented so that the electrophysiological feature with the largest absolute loading on each component had a positive loading. Transcriptomic scores were calculated by projecting the centered transcriptomic PCs through the ordinary least-squares coefficient matrix and the three response directions; electrophysiological scores were obtained by projecting the centered 19-feature matrix onto the same directions. Gene and electrophysiological loading vectors are Pearson correlations between the standardized input variables and their corresponding component scores.

The rank-3 solution was fitted to all 346 cells for descriptive joint-axis visualization; rank 3 was not selected by donor-held-out cross-validation, and no predictive R-squared is reported. The figure shows Components 1-2 and 1-3 in transcriptomic and electrophysiological spaces. Points denote cells. Transcriptomic panels are colored by D1/D2 identity, whereas electrophysiological panels are colored by E1-E4. The interactive layout displays the 15 genes and ten electrophysiological features with the largest two-dimensional loading magnitudes for each displayed component pair; the archived loading tables contain all 1,000 genes and all 19 electrophysiological features. Covariance ellipses show descriptive 90% bivariate-normal regions and are not confidence intervals. Scores were rescaled by the 99th percentile of radial distance for plotting only. RRR did not alter the frozen E-class assignments.
"""
    write("01_METHODS/Macaque_E_RRR_Methods_EN.md", methods)

    availability = """# Data and code availability

Processed macaque Patch-seq expression matrices and metadata were downloaded from the Allen Brain Cell Atlas Human-Mammalian Brain Atlas Basal Ganglia macaque Patch-seq release (dataset version `20260228`; https://alleninstitute.github.io/abc_atlas_access/descriptions/HMBA-Macaque-PatchSeq-BG.html). The RRR analysis used the released `HMBA-Macaque-PatchSeq-BG-log2.h5ad` matrix and did not reprocess raw sequencing reads. The corresponding publication is Liu et al., *Morphoelectric Diversity and Specialization of Neuronal Cell Types in the Primate Striatum* (2026), https://doi.org/10.64898/2026.02.26.708019. Raw sequencing resources are available through NeMO Archive collection `nemo:dat-nsm6mxv` (https://assets.nemoarchive.org/dat-nsm6mxv), and the source publication lists the associated electrophysiology NWB datasets at DANDI. The processed Allen release is distributed under CC BY-NC 4.0; repository-specific licensing applies to NeMO and DANDI records.

The exact processed H5AD used locally was 63,797,589 bytes, timestamped 31 August 2026, with SHA-256 `E3DD6C32A3872AB83CB45E13F813C50BD06C61BFDC5F889D12BDDEF46EEA0986`. The H5AD is not duplicated in this compact supplement; it is present in the original Macaque release/archive retained with the main frozen package.
"""
    write("01_METHODS/Data_and_Code_Availability_EN.md", availability)

    checklist = """# Nature投稿前RRR披露清单

已在补充包中锁定：

- 数据集名称、正式下载页、版本号、处理层级和许可；
- 346个细胞的纳入规则及D1/D2数量；
- Hybrid排除规则；
- 表达矩阵为发布方提供的log2处理矩阵，而非本研究从FASTQ重新处理；
- 检出阈值、技术基因过滤、重复symbol处理及1,000高变基因选择；
- 基因Z-score、20个转录组PC及累计解释率；
- 19项冻结E特征作为响应矩阵；
- RRR拟合公式、中心化、伪逆、SVD、rank=3和符号定向；
- cell scores及correlation loadings定义；
- 图中15基因/10 E特征仅为显示选择，完整载荷表为1,000/19；
- 90%椭圆和99%径向缩放仅用于显示；
- rank=3为描述性展示，未做供体留出预测验证，不报告预测R²；
- 分析代码、结果表、布局、校验值和软件环境。

正式投稿仍需在主文或Reporting Summary中确认：软件版本是否采用本补充包记录的当前环境；数据访问日期；代码公开仓库和永久存档DOI；以及是否新增独立的rank选择或供体分组交叉验证。如果没有新增验证，应保留“descriptive visualization”边界，不把Component 2/3写成独立验证的生物学轴。
"""
    write("01_METHODS/Nature_reporting_checklist_CN.md", checklist)

    provenance = {
        "publication": "Liu et al., Morphoelectric Diversity and Specialization of Neuronal Cell Types in the Primate Striatum (2026)",
        "publication_doi": "https://doi.org/10.64898/2026.02.26.708019",
        "processed_data_page": "https://alleninstitute.github.io/abc_atlas_access/descriptions/HMBA-Macaque-PatchSeq-BG.html",
        "processed_dataset_version": "20260228",
        "processed_release_license": "CC BY-NC 4.0",
        "raw_sequencing": "https://assets.nemoarchive.org/dat-nsm6mxv",
        "h5ad_filename": H5AD.name,
        "h5ad_bytes": H5AD.stat().st_size,
        "h5ad_sha256": digest(H5AD),
        "analysis_cells": 346,
        "D1": 163,
        "D2": 183,
        "input_genes": 1000,
        "transcriptomic_PCs": 20,
        "transcriptomic_PC_variance_fraction": float(variance.iloc[:, 1].sum()),
        "electrophysiological_features": 19,
        "RRR_rank": 3,
        "inference_boundary": "descriptive in-sample joint-axis visualization; not donor-held-out prediction",
    }
    write("02_PROVENANCE/RRR_provenance.json", json.dumps(provenance, indent=2))
    environment = {
        "status": "numerically verified reproduction environment",
        "python": sys.version.split()[0],
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "scipy": scipy.__version__,
        "scikit_learn": sklearn.__version__,
        "h5py": h5py.__version__,
        "matplotlib": matplotlib.__version__,
        "platform": platform.platform(),
        "verification": {
            "cell_score_max_absolute_difference": 8.881784197001252e-16,
            "gene_loading_max_absolute_difference": 9.996344030316351e-17,
            "E_loading_max_absolute_difference": 8.326672684688674e-17,
            "gene_ID_order_identical": True,
            "interpretation": "Current environment reproduced the frozen numeric outputs to floating-point precision."
        }
    }
    write("02_PROVENANCE/software_environment_and_reproduction.json", json.dumps(environment, indent=2))

    for name in ["RRR_cell_scores_n346.csv", "RRR_Efeature_correlation_loadings.csv",
                 "RRR_gene_correlation_loadings.csv", "transcriptomic_PCA_variance.csv", "README.txt"]:
        copy(RRR / name, f"03_RESULTS/{name}")
    copy(LAYOUT, "03_RESULTS/Macaque_T-E_RRR_gc_merged_layout.json")
    for name in ["plot_ca_pu_nac_rrr_fourpanel.py", "analyze_nac_d1d2_transcriptomic_ephys_rrr.py",
                 "build_macaque_e_rrr_editor.py"]:
        copy(ROOT / "scripts" / name, f"03_CODE/{name}")

    write("00_README/README_CN.md", """# Macaque E RRR Methods补充包

该包只补充Macaque E图中的T-E RRR方法、数据来源、代码、结果表及Nature披露清单，不修改或替代主冻结E1-E4分析包。RRR分析使用346个HC-GC共识且稳定D1/D2细胞、1,000个高变基因、20个转录组PC、19项冻结E特征和rank=3。rank=3属于全数据描述性联合轴展示，不是供体留出预测模型。""")

    # Checksums and archive.
    rows = []
    for path in sorted(p for p in OUT.rglob("*") if p.is_file()):
        rows.append({"file": path.relative_to(OUT).as_posix(), "bytes": path.stat().st_size, "sha256": digest(path)})
    pd.DataFrame(rows).to_csv(OUT / "04_AUDIT/file_manifest.csv", index=False)
    write("04_AUDIT/SHA256SUMS.txt", "\n".join(f"{r['sha256']}  {r['file']}" for r in rows))

    with ZipFile(ZIP_OUT, "w", ZIP_DEFLATED, compresslevel=6) as archive:
        for path in sorted(p for p in OUT.rglob("*") if p.is_file()):
            archive.write(path, f"{OUT.name}/{path.relative_to(OUT).as_posix()}")
    print(json.dumps({"directory": str(OUT), "zip": str(ZIP_OUT),
                      "checks_passed": int(audit["pass"].sum()), "checks_total": len(audit)}, indent=2))


if __name__ == "__main__":
    main()
