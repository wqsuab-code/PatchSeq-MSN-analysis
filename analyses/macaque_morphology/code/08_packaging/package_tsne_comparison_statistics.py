from __future__ import annotations

import csv
import hashlib
import json
import shutil
import zipfile
from datetime import datetime
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
NAME = "Macaque_M_tsne_comparison_statistics_bundle_20260921"
OUTPUTS = ROOT / "outputs"
BUNDLE = OUTPUTS / NAME
ZIP_PATH = OUTPUTS / f"{NAME}.zip"
RUN = ROOT / "macaque_m" / "m18_tempfreeze_NPC5_HCK4_res2.3"
COMPARE = RUN / "M4_compare_all18"
EDITOR = ROOT / "interactive" / "macaque-morph-tsne-comparison-editor" / "dist"
_layout_dir = Path(r"C:\Users\53461\Downloads")
_layout_candidates = sorted(
    _layout_dir.glob("Macaque_M_morphology_tsne_comparison_layout*.json"),
    key=lambda path: path.stat().st_mtime,
    reverse=True,
)
LATEST_LAYOUT = _layout_candidates[0] if _layout_candidates else _layout_dir / "Macaque_M_morphology_tsne_comparison_layout.json"


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def reset_exact_output(path: Path) -> None:
    resolved = path.resolve()
    output_root = OUTPUTS.resolve()
    if resolved == output_root or output_root not in resolved.parents:
        raise RuntimeError(f"Unsafe output target: {resolved}")
    if path.is_dir():
        shutil.rmtree(path)
    elif path.exists():
        path.unlink()


def copy(source: Path, destination: str) -> None:
    if not source.is_file():
        raise FileNotFoundError(source)
    target = BUNDLE / destination
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, target)


def copy_tree(source: Path, destination: str) -> None:
    if not source.is_dir():
        raise FileNotFoundError(source)
    shutil.copytree(source, BUNDLE / destination, ignore=shutil.ignore_patterns("__pycache__", "*.pyc", "~$*"))


def write(relative: str, text: str) -> None:
    path = BUNDLE / relative
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text.rstrip() + "\n", encoding="utf-8-sig")


def main() -> None:
    OUTPUTS.mkdir(parents=True, exist_ok=True)
    if BUNDLE.exists():
        reset_exact_output(BUNDLE)
    if ZIP_PATH.exists():
        reset_exact_output(ZIP_PATH)
    BUNDLE.mkdir(parents=True)

    # Final figure-ready cell-level data and complete statistical outputs.
    for filename in (
        "all18_plot_data_117.csv",
        "all18_Dunn_Holm.csv",
        "all18_KruskalWallis_BH.csv",
        "Macaque_Morphology_ConsensusM4_all18_3perrow_W7p2_H8p76.pdf",
        "Macaque_Morphology_ConsensusM4_all18_3perrow_W7p2_H8p76.png",
        "Macaque_Morphology_ConsensusM4_all18_distribution_tSNE.png",
    ):
        copy(COMPARE / filename, f"01_final_figure_data_and_results/{filename}")
    for filename in (
        "STATISTICS_pairwise_Dunn_Holm_BH108.csv",
        "STATISTICS_omnibus_KruskalWallis_BH18.csv",
        "FEATURE_IMPORTANCE_PANEL_ORDER.csv",
        "DATA_AND_LAYOUT_AUDIT.json",
    ):
        copy(EDITOR / filename, f"01_final_figure_data_and_results/{filename}")

    # Frozen inputs required to regenerate the figure data without Mouse/E data.
    inputs = {
        ROOT / "macaque_m" / "m18" / "01_raw_126.csv": "02_frozen_inputs/01_raw_morphology_126.csv",
        RUN / "01_temp_frozen_assignments_126.csv": "02_frozen_inputs/02_frozen_HC_GC_assignments_126.csv",
        ROOT / "macaque_m" / "m18_adaptive_pca126" / "02_transformed_z_117.csv": "02_frozen_inputs/03_transformed_Z_complete126.csv",
        ROOT / "macaque_m" / "m18_adaptive_pca126" / "05_pca_loadings.csv": "02_frozen_inputs/04_frozen_PCA_loadings_18features.csv",
        RUN / "09_tSNE_optimized_coordinates.csv": "02_frozen_inputs/05_frozen_tSNE_coordinates.csv",
    }
    for source, destination in inputs.items():
        copy(source, destination)

    # Statistical and visualization code used by the full chain.
    code_files = {
        ROOT / "macaque_m" / "run_m18_adaptive_transform_pca.R": "03_code/01_transform_Zscore_PCA.R",
        ROOT / "macaque_m" / "freeze_npcs5_hck4_res23.R": "03_code/02_freeze_HC_GC_M4.R",
        ROOT / "macaque_m" / "optimize_tempfreeze_k4_tsne.R": "03_code/03_optimize_tSNE.R",
        ROOT / "macaque_m" / "plot_tsne_final_85pct.R": "03_code/04_original_tSNE_display.R",
        ROOT / "macaque_m" / "plot_consensusM4_top10_comparison.py": "03_code/05_KruskalWallis_Dunn_and_static_figure.py",
        ROOT / "macaque_m" / "build_macaque_m_tsne_comparison_editor.py": "03_code/06_build_interactive_figure_and_BH108.py",
        Path(__file__).resolve(): "03_code/07_build_this_statistics_bundle.py",
    }
    for source, destination in code_files.items():
        copy(source, destination)

    # ML importance sources used only to determine panel order.
    copy(
        ROOT / "outputs" / "Macaque_M4_full_ML500" / "tables" / "20_feature_importance_MDI_logistic_500.csv",
        "04_feature_importance_sources/ExtraTrees_MDI_and_logistic_500.csv",
    )
    copy(
        ROOT / "outputs" / "Macaque_M4_full_ML500" / "tables" / "21_nested_grouped_permutation_importance_500.csv",
        "04_feature_importance_sources/donor_grouped_permutation_importance_500.csv",
    )
    copy(
        ROOT / "outputs" / "Macaque_M4_full_ML500" / "tables" / "24_feature_rank_correlation.csv",
        "04_feature_importance_sources/cross_model_rank_correlation.csv",
    )

    # Latest user-edited layout and complete self-contained editor.
    copy(LATEST_LAYOUT, f"05_layout/{LATEST_LAYOUT.name}")
    copy_tree(EDITOR, "06_interactive_editor")

    methods = """
# Statistical methods and package guide

## Scope and cohort

The figure is restricted to the 117 frozen HC–GC consensus macaque MSN morphology cells (M1, n=43; M2, n=42; M3, n=20; M4, n=12). Only the 18 features present in the frozen PCA loading table are shown. No mouse observations, E-class predictors, or E preprocessing parameters enter the analysis.

## Raw-value group comparisons

Each box/scatter plot uses the untransformed raw value of the corresponding morphology feature. An omnibus two-sided Kruskal–Wallis test compares M1–M4. The 18 omnibus P values are adjusted by the Benjamini–Hochberg procedure. Pairwise contrasts use Dunn's rank-sum test with tie correction. For each feature, its six pairwise P values are adjusted by Holm's family-wise-error procedure. The editor additionally provides Benjamini–Hochberg adjustment across all 18 × 6 = 108 pairwise tests. Dunn effect size is reported as |Z|/sqrt(N), where N=117.

## t-SNE feature maps

All 18 panels reuse the same frozen t-SNE coordinates; only point color changes. Point color is the saved transformed feature Z score, displayed on a default range of -2.5 to 2.5. The current class boundary is a descriptive visualization, not a confidence region: within each M class, coordinates are centered by the coordinate-wise median and scaled by MAD; cells are ranked by robust radial distance; the selected fraction is enclosed by a convex hull, expanded about the robust center, and smoothed with Chaikin subdivision. The supplied layout currently sets coverage=0.90, expansion=1.00, and six smoothing passes.

## Feature-importance panel order

The 18 panels are ordered by the mean of three descending ranks calculated from the donor-grouped ML validation: Extra Trees median mean decrease in impurity, median absolute standardized multinomial-logistic coefficient, and median nested donor-grouped permutation importance. Ties in mean rank are resolved by higher Extra Trees MDI and then the feature name. This ranking changes panel order only; it does not alter the group-comparison tests.

## Display-only scaling

Total dendritic length tick labels are divided by 1,000 and the axis is labeled ×10^3 µm. Raw values, point positions, and statistical calculations remain unchanged.

## Directory guide

- 01_final_figure_data_and_results: cell-level plotting table, all omnibus and pairwise results, importance order, audit, and figure files.
- 02_frozen_inputs: raw morphology, frozen assignments, transformed Z scores, PCA whitelist/loadings, and t-SNE coordinates.
- 03_code: ordered end-to-end scripts for preprocessing, frozen classification, t-SNE, statistical testing, interactive rendering, and package construction.
- 04_feature_importance_sources: the three ML files used to derive panel order.
- 05_layout: latest user-exported presentation parameters.
- 06_interactive_editor: self-contained editor source and embedded data.

## Interpretive boundary

The tests describe differences among labels derived from the same morphology feature space. They do not constitute independent biological validation of M1–M4. The t-SNE boundaries are descriptive display aids and must not be interpreted as inferential confidence limits.
"""
    methods_cn = """
# 统计方法与文件说明

该图仅包含117个冻结HC–GC共识Macaque MSN细胞：M1=43、M2=42、M3=20、M4=12；仅使用冻结PCA载荷表中的18项形态指标。未使用Mouse数据、E特征或E预处理参数。

箱线/散点图使用原始形态数值。四组总体比较采用双侧Kruskal–Wallis检验，18项总体P值采用Benjamini–Hochberg校正。两两比较采用带并列值校正的Dunn秩和检验，每项指标内部6次比较采用Holm校正；互动页面另提供全部108次两两检验的全局BH校正。Dunn效应量定义为|Z|/sqrt(N)，N=117。

18个t-SNE面板共用相同冻结坐标，只有点颜色随指标变化；颜色为保存的变换后Z-score，默认显示范围−2.5至2.5。轮廓是描述性边界而非置信区间：每个M组以坐标中位数为中心、MAD为尺度，按稳健径向距离选取指定比例细胞，计算凸包，围绕中心扩张，再进行Chaikin平滑。当前布局为覆盖率0.90、扩张1.00、平滑6次。

panel顺序使用三种供体分组机器学习重要性的平均名次：Extra Trees MDI、多项逻辑回归绝对标准化系数、嵌套供体分组置换重要性。该排序只改变版面顺序，不改变统计结果。

Total length的Y轴刻度仅在显示时除以1000并标记为×10³ µm；原始数据、点的位置和统计计算均不改变。

注意：这些检验描述在同一形态特征空间中定义的M标签之间的差异，不能作为M1–M4的独立生物学验证。t-SNE轮廓也不能解释为统计置信范围。
"""
    write("README_METHODS_EN.md", methods)
    write("README_METHODS_CN.md", methods_cn)

    rows = []
    for path in sorted(item for item in BUNDLE.rglob("*") if item.is_file()):
        rows.append(
            {
                "relative_path": path.relative_to(BUNDLE).as_posix(),
                "bytes": path.stat().st_size,
                "sha256": digest(path),
            }
        )
    with (BUNDLE / "FILE_INDEX.csv").open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=["relative_path", "bytes", "sha256"])
        writer.writeheader()
        writer.writerows(rows)

    checksum_lines = []
    for path in sorted(item for item in BUNDLE.rglob("*") if item.is_file()):
        checksum_lines.append(f"{digest(path)}  {path.relative_to(BUNDLE).as_posix()}")
    write("SHA256SUMS.txt", "\n".join(checksum_lines))
    manifest = {
        "bundle": NAME,
        "created": datetime.now().isoformat(timespec="seconds"),
        "scope": "Macaque M 18-feature t-SNE plus raw-value group comparisons",
        "cells": 117,
        "class_counts": {"M1": 43, "M2": 42, "M3": 20, "M4": 12},
        "features": 18,
        "pairwise_tests": 108,
        "files": len([item for item in BUNDLE.rglob("*") if item.is_file()]),
    }
    write("MANIFEST.json", json.dumps(manifest, ensure_ascii=False, indent=2))

    with zipfile.ZipFile(ZIP_PATH, "w", zipfile.ZIP_DEFLATED, compresslevel=6) as archive:
        for path in sorted(item for item in BUNDLE.rglob("*") if item.is_file()):
            archive.write(path, Path(NAME) / path.relative_to(BUNDLE))
    print(
        json.dumps(
            {
                "zip": str(ZIP_PATH),
                "bytes": ZIP_PATH.stat().st_size,
                "sha256": digest(ZIP_PATH),
                "files": len([item for item in BUNDLE.rglob("*") if item.is_file()]),
            },
            ensure_ascii=False,
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
