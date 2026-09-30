from __future__ import annotations

import csv
import hashlib
import json
import shutil
import zipfile
from datetime import datetime
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
STAMP = "20260921"
BUNDLE_NAME = f"Macaque_M4_ML_writing_bundle_{STAMP}"
OUT_ROOT = ROOT / "outputs"
BUNDLE = OUT_ROOT / BUNDLE_NAME
ARCHIVE = OUT_ROOT / f"{BUNDLE_NAME}.zip"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def checked_remove(path: Path) -> None:
    resolved = path.resolve()
    output_resolved = OUT_ROOT.resolve()
    if resolved == output_resolved or output_resolved not in resolved.parents:
        raise RuntimeError(f"Refusing to remove path outside outputs: {resolved}")
    if path.is_dir():
        shutil.rmtree(path)
    elif path.exists():
        path.unlink()


def copy_file(source: Path, relative_destination: Path) -> None:
    if not source.is_file():
        raise FileNotFoundError(source)
    destination = BUNDLE / relative_destination
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, destination)


def copy_tree(source: Path, relative_destination: Path) -> None:
    if not source.is_dir():
        raise FileNotFoundError(source)
    destination = BUNDLE / relative_destination
    shutil.copytree(
        source,
        destination,
        ignore=shutil.ignore_patterns("__pycache__", "*.pyc", "*.tmp", "~$*"),
    )


def write_text(relative_path: str, text: str) -> None:
    path = BUNDLE / relative_path
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text.rstrip() + "\n", encoding="utf-8-sig")


def main() -> None:
    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    if BUNDLE.exists():
        checked_remove(BUNDLE)
    if ARCHIVE.exists():
        checked_remove(ARCHIVE)
    BUNDLE.mkdir(parents=True)

    # Current authoritative ML500 analysis: figures, tables, raw repetitions,
    # code snapshots, methods, parameters, conclusions, and bilingual READMEs.
    copy_tree(
        ROOT / "outputs" / "Macaque_M4_full_ML500",
        Path("01_authoritative_ML500"),
    )

    # Frozen M1-M4 inputs and their preprocessing/PCA provenance.
    frozen_dir = ROOT / "macaque_m" / "m18_tempfreeze_NPC5_HCK4_res2.3"
    for name in (
        "01_temp_frozen_assignments_126.csv",
        "02_GC13_to_GC4_merge_map.csv",
        "03_HC4_GC4_confusion.csv",
        "04_temp_frozen_crossplot.pdf",
        "04_temp_frozen_crossplot.png",
    ):
        copy_file(frozen_dir / name, Path("02_frozen_M4_provenance") / name)
    copy_tree(
        ROOT / "macaque_m" / "m18_adaptive_pca126",
        Path("02_frozen_M4_provenance") / "adaptive_transform_and_PCA",
    )
    copy_tree(
        ROOT / "macaque_m" / "m18_adaptive126_all_confusions",
        Path("02_frozen_M4_provenance") / "HC_GC_parameter_search",
    )
    for source_name in (
        "run_m18_adaptive_transform_pca.R",
        "freeze_npcs5_hck4_res23.R",
    ):
        copy_file(
            ROOT / "macaque_m" / source_name,
            Path("02_frozen_M4_provenance") / "source_scripts" / source_name,
        )

    # Manuscript-facing morphology panels that explain the frozen phenotype.
    for directory_name in (
        "HC_GC_jitter",
        "M4_compare_all18",
        "M4_compare_top10",
        "panels_A_I",
        "RRR_T_M",
        "tSNE_consensus_GC_Tclass",
    ):
        copy_tree(
            frozen_dir / directory_name,
            Path("03_supporting_M4_figures") / directory_name,
        )

    # Interactive reconstruction of the ML figure and its reproducible builder.
    copy_tree(
        ROOT / "interactive" / "macaque-morph-ml-studio",
        Path("04_interactive_figure_studio") / "macaque-morph-ml-studio",
    )
    copy_file(
        ROOT / "scripts" / "build_macaque_morph_ml_eplatform.py",
        Path("04_interactive_figure_studio") / "build_macaque_morph_ml_eplatform.py",
    )
    copy_file(
        Path(__file__).resolve(),
        Path("07_bundle_builder") / Path(__file__).name,
    )

    # These files are copied strictly as visual layout references. They are not
    # data inputs and must never be used as Macaque M predictors or preprocessing.
    layout_json = Path(
        r"C:\Users\53461\OneDrive\Desktop\Patch-seq\4Manuscript\Macaque_E4_ML_layout (Ver20260921).json"
    )
    layout_svg = Path(
        r"C:\Users\53461\OneDrive\Desktop\Patch-seq\4Manuscript\Macaque_E4_ML_interactive_export.svg"
    )
    copy_file(layout_json, Path("05_layout_reference_only_DO_NOT_USE_AS_DATA") / layout_json.name)
    copy_file(layout_svg, Path("05_layout_reference_only_DO_NOT_USE_AS_DATA") / layout_svg.name)
    write_text(
        "05_layout_reference_only_DO_NOT_USE_AS_DATA/README.txt",
        "These E4 files are included only to document the visual layout used by the interactive figure.\n"
        "They are not Macaque M data, predictors, labels, PCA scores, preprocessing parameters, or statistical results.",
    )

    raw_zip = Path(r"C:\Users\53461\Downloads\Macaque-PatchSeq-BG.zip")
    raw_zip_record = {
        "source_path": str(raw_zip),
        "exists_at_packaging_time": raw_zip.is_file(),
        "bytes": raw_zip.stat().st_size if raw_zip.is_file() else None,
        "sha256": sha256(raw_zip) if raw_zip.is_file() else None,
        "included_in_bundle": False,
        "reason": "External source archive is not duplicated; its exact path and checksum are recorded for provenance.",
    }
    write_text(
        "06_source_data_provenance/external_raw_archive.json",
        json.dumps(raw_zip_record, ensure_ascii=False, indent=2),
    )

    readme_cn = """
Macaque MSN 形态学 M1–M4：ML500 撰写与复现文件包

用途
本文件包汇总了撰写 Figure legend、Methods 和 Results 所需的当前权威文件，并保留冻结分类与交互图的可追溯来源。分析对象仅为 Macaque MSN；Mouse 数据和 E 分类特征未作为分析输入。

冻结样本与分类
- 完整形态重建病例：126 个。
- HC–GC 共识细胞：117 个；M1=43、M2=42、M3=20、M4=12。
- 冻结主流程：按偏态选择 Yeo–Johnson 或不变换，随后 Z-score、PCA（5 PCs）、Ward.D2 HC（K=4）；GC resolution=2.3 的 13 个原始群合并为 4 类。
- ML 验证中的所有预处理均在训练供体/重抽样内部重新估计，再原样应用于留出供体，以避免数据泄漏。

当前权威结果摘要
- 冻结流程复现：ARI=1.000，NMI=1.000，错配 0/117。
- 供体外监督恢复（RBF SVM）：balanced accuracy 中位数 0.909（95%经验区间 0.730–1.000）；macro-F1 0.899（0.713–1.000）；供体结构保留置换检验 P=0.001996。
- 从头无监督发现（冻结位置 Ward、NPC=5、K=4）：ARI 中位数 0.322（0.107–0.626），NMI 中位数约 0.437，silhouette 中位数约 0.155。
- 细胞稳定性：16 个细胞 consensus margin≤0；27 个细胞 margin<0.1。
- 跨模型特征排序：Spearman rho=0.674。
- 元数据关联：ROI Cramér’s V=0.288，P=0.00599；donor、T class 和 E class 在 0.05 水平未见显著关联。

目录导航
1. 01_authoritative_ML500：最终图、统计表、500 次原始结果、ROC/PR 坐标、参数、脚本及中英文结论。写 Methods/Results 时以此目录为权威来源。
2. 02_frozen_M4_provenance：冻结细胞名单、M1–M4/HC/GC 标签、变换审计、PCA、离群诊断及聚类参数搜索。
3. 03_supporting_M4_figures：解释冻结形态分类的主图/补图，包括 18 指标分布、热图、t-SNE、RRR 和 HC–GC 对照。
4. 04_interactive_figure_studio：当前交互式 a–p 机器学习图及其构建脚本。
5. 05_layout_reference_only_DO_NOT_USE_AS_DATA：仅用于证明版式来源的 E4 JSON/SVG，严禁作为 M 分析数据或参数。
6. 06_source_data_provenance：外部原始 ZIP 的路径、大小和 SHA-256；为避免重复，ZIP 本体未再次打包。

撰写边界
高监督准确率表示 M1–M4 在用于定义它们的形态特征中可恢复，不等同于独立生物学验证。监督可恢复性、无监督发现稳定性、细胞级稳定性和独立生物学支持必须分开表述。

完整逐文件清单与用途见 FILE_INDEX.csv；完整校验值见 SHA256SUMS.txt。
"""
    readme_en = """
Macaque MSN morphology M1–M4: ML500 writing and reproducibility bundle

Purpose
This bundle collects the current authoritative files needed to write the figure legend, Methods, and Results, together with the provenance of the frozen classification and interactive figure. The analysis is restricted to macaque MSNs; mouse data and E-class predictors were not used as analytical inputs.

Frozen cohort and classification
- Complete morphology reconstructions: n=126.
- HC–GC consensus cells: n=117; M1=43, M2=42, M3=20, and M4=12.
- Frozen main workflow: skewness-adaptive Yeo–Johnson transformation or no transformation, followed by z-scoring, PCA (five PCs), Ward.D2 HC (K=4), and merging of 13 raw GC communities at resolution 2.3 into four classes.
- During ML validation, every preprocessing parameter was re-estimated within the training donors/resample and then applied unchanged to held-out donors to prevent data leakage.

Authoritative result summary
- Exact frozen-pipeline reproduction: ARI=1.000, NMI=1.000, and 0/117 mismatches.
- Donor-held-out supervised recoverability (RBF SVM): median balanced accuracy 0.909 (95% empirical interval, 0.730–1.000), median macro-F1 0.899 (0.713–1.000), and donor-structure-preserving permutation P=0.001996.
- De novo unsupervised discovery at the frozen Ward/NPC=5/K=4 position: median ARI 0.322 (0.107–0.626), median NMI approximately 0.437, and median silhouette approximately 0.155.
- Cell-level stability: 16 cells had consensus margin <=0 and 27 had margin <0.1.
- Cross-model feature-rank agreement: Spearman rho=0.674.
- Metadata association: ROI Cramer's V=0.288, P=0.00599; donor, T class, and E class were not significant at alpha=0.05.

Directory guide
1. 01_authoritative_ML500: final figures, statistical tables, all 500-run raw outputs, ROC/PR coordinates, parameters, scripts, and bilingual conclusions. Treat this directory as authoritative when writing Methods and Results.
2. 02_frozen_M4_provenance: frozen cell list and M/HC/GC labels, transformation audit, PCA, outlier diagnostics, and clustering parameter search.
3. 03_supporting_M4_figures: manuscript-facing plots describing the frozen morphology classes, including all-18-feature distributions, heatmap, t-SNE, RRR, and HC–GC comparisons.
4. 04_interactive_figure_studio: the current interactive a–p ML figure and reproducible builder.
5. 05_layout_reference_only_DO_NOT_USE_AS_DATA: E4 JSON/SVG retained only to document visual layout; never use them as M data or parameters.
6. 06_source_data_provenance: path, size, and SHA-256 of the external raw ZIP. The raw ZIP itself is not duplicated.

Interpretive boundary
High supervised accuracy demonstrates recoverability of M1–M4 from the morphology features used to define them; it is not independent biological validation. Supervised recoverability, unsupervised discovery stability, cell-level stability, and independent biological support must be reported separately.

See FILE_INDEX.csv for the complete file-by-file guide and SHA256SUMS.txt for integrity checks.
"""
    write_text("README_CN.txt", readme_cn)
    write_text("README_EN.txt", readme_en)

    # Build a complete, auditable index after all content has been staged.
    category_purpose = {
        "01_authoritative_ML500": "Authoritative ML Methods, Results, figures, tables, raw repetitions, and code",
        "02_frozen_M4_provenance": "Frozen input, preprocessing, PCA, HC/GC labels, and clustering provenance",
        "03_supporting_M4_figures": "Supporting manuscript figures and their source tables",
        "04_interactive_figure_studio": "Interactive figure, embedded data, layout, and reproducible builder",
        "05_layout_reference_only_DO_NOT_USE_AS_DATA": "Visual layout reference only; never an analytical input",
        "06_source_data_provenance": "External raw-data archive identity and checksum",
        "07_bundle_builder": "Script that rebuilds this curated writing and reproducibility bundle",
    }
    rows = []
    for path in sorted(p for p in BUNDLE.rglob("*") if p.is_file()):
        rel = path.relative_to(BUNDLE).as_posix()
        top = rel.split("/", 1)[0]
        rows.append(
            {
                "category": top,
                "relative_path": rel,
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
                "suggested_use": category_purpose.get(top, "Bundle documentation"),
            }
        )
    index_path = BUNDLE / "FILE_INDEX.csv"
    with index_path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)

    checksum_lines = []
    for path in sorted(p for p in BUNDLE.rglob("*") if p.is_file()):
        if path.name == "SHA256SUMS.txt":
            continue
        checksum_lines.append(f"{sha256(path)}  {path.relative_to(BUNDLE).as_posix()}")
    write_text("SHA256SUMS.txt", "\n".join(checksum_lines))

    manifest = {
        "bundle": BUNDLE_NAME,
        "created_local": datetime.now().isoformat(timespec="seconds"),
        "scope": "Macaque MSN morphology M1-M4 ML500 writing and reproducibility",
        "file_count_excluding_manifest": len([p for p in BUNDLE.rglob("*") if p.is_file()]),
        "external_raw_archive": raw_zip_record,
    }
    write_text("BUNDLE_MANIFEST.json", json.dumps(manifest, ensure_ascii=False, indent=2))

    with zipfile.ZipFile(ARCHIVE, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=6) as archive:
        for path in sorted(p for p in BUNDLE.rglob("*") if p.is_file()):
            archive.write(path, Path(BUNDLE_NAME) / path.relative_to(BUNDLE))

    archive_record = {
        "archive": str(ARCHIVE),
        "bytes": ARCHIVE.stat().st_size,
        "sha256": sha256(ARCHIVE),
        "packaged_files": len([p for p in BUNDLE.rglob("*") if p.is_file()]),
    }
    print(json.dumps(archive_record, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
