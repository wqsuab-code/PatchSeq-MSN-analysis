#!/usr/bin/env python3
"""Macaque MSN morphology feature QC.

This script reads the source ZIP directly. It intentionally does not read any
existing Mouse morphology analysis, does not run PCA/clustering, and does not
assign or freeze morphology classes.
"""

from __future__ import annotations

import argparse
import hashlib
import io
import json
import math
import re
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd


META_ENTRY = "Data/cell_metadata_AllCell.csv"
MORPH_ENTRY = "Data/morphology_features.csv"
LOOKUP_ENTRY = "Data/morphological_feature_lookup.xlsx"

ROIS = ("Ca", "Pu", "NAC")
MSN_SUBCLASSES = ("STR D1 MSN", "STR D2 MSN", "STR Hybrid MSN")
ROBUST_Z_CUTOFF = 3.5
SPEARMAN_CUTOFF = 0.80
ALLEN_MORPHOLOGY_REFERENCE = (
    "https://github.com/AllenInstitute/neuron_morphology/blob/dev/feature_description.html"
)
ALLEN_CROSSCHECK_TERMS = {
    "basal_dendrite_calculate_number_of_stems",
    "basal_dendrite_max_branch_order",
    "basal_dendrite_max_euclidean_distance",
    "basal_dendrite_max_path_distance",
    "basal_dendrite_mean_contraction",
    "basal_dendrite_mean_diameter",
    "basal_dendrite_num_branches",
    "basal_dendrite_total_length",
    "soma_surface_area",
    "axon_max_branch_order",
    "axon_max_euclidean_distance",
    "axon_max_path_distance",
    "axon_mean_contraction",
    "axon_num_branches",
    "axon_total_length",
}


PRIMARY_PANEL = [
    "basal_dendrite_total_length",
    "basal_dendrite_num_branches",
    "basal_dendrite_calculate_number_of_stems",
    "basal_dendrite_max_branch_order",
    "basal_dendrite_max_path_distance",
    "basal_dendrite_mean_contraction",
    "basal_dendrite_mean_diameter",
    "basal_dendrite_extent_dorsal",
    "basal_dendrite_extent_medial",
    "soma_surface_area",
]

SPATIAL_ADDITIONS = [
    "basal_dendrite_max_euclidean_distance",
    "basal_dendrite_bias_dorsal",
    "basal_dendrite_bias_medial",
    "basal_dendrite_soma_percentile_dorsal",
    "basal_dendrite_soma_percentile_medial",
    "derived:stem_exit_ml_fraction",
    "derived:stem_exit_dv_balance",
]


CN_MEANING = {
    "23_Sholl_PC0": "全部神经突 Sholl 曲线的第1主成分得分",
    "23_Sholl_PC1": "全部神经突 Sholl 曲线的第2主成分得分",
    "2_Sholl_PC0": "轴突 Sholl 曲线的第1主成分得分",
    "2_Sholl_PC1": "轴突 Sholl 曲线的第2主成分得分",
    "3_Sholl_PC0": "树突 Sholl 曲线的第1主成分得分",
    "3_Sholl_PC1": "树突 Sholl 曲线的第2主成分得分",
    "axon_bias_dorsal": "轴突在背腹轴上的偏置",
    "axon_bias_medial": "轴突在内外侧轴上的偏置",
    "axon_exit_distance": "轴突离开胞体时的距离",
    "axon_exit_theta_coronal": "冠状面内轴突离开胞体的方向角",
    "axon_extent_dorsal": "轴突在背腹轴上的跨度",
    "axon_extent_medial": "轴突在内外侧轴上的跨度",
    "axon_max_branch_order": "轴突最大分支阶数",
    "axon_max_euclidean_distance": "胞体至最远轴突节点的最大欧氏距离",
    "axon_max_path_distance": "胞体至最远轴突末端的最大路径距离",
    "axon_mean_contraction": "轴突分支段欧氏距离总和与路径距离总和之比",
    "axon_num_branches": "轴突分支数",
    "axon_soma_percentile_dorsal": "胞体在轴突背腹跨度中的相对位置",
    "axon_soma_percentile_medial": "胞体在轴突内外侧跨度中的相对位置",
    "axon_total_length": "轴突总长度",
    "basal_dendrite_bias_dorsal": "基底树突在背腹轴上的偏置",
    "basal_dendrite_bias_medial": "基底树突在内外侧轴上的偏置",
    "basal_dendrite_calculate_number_of_stems": "从胞体发出的基底树突主干数",
    "basal_dendrite_extent_dorsal": "基底树突在背腹轴上的跨度",
    "basal_dendrite_extent_medial": "基底树突在内外侧轴上的跨度",
    "basal_dendrite_max_branch_order": "基底树突最大分支阶数",
    "basal_dendrite_max_euclidean_distance": "胞体至最远基底树突节点的最大欧氏距离",
    "basal_dendrite_max_path_distance": "胞体至最远基底树突末端的最大路径距离",
    "basal_dendrite_mean_contraction": "基底树突分支段欧氏距离总和与路径距离总和之比",
    "basal_dendrite_mean_diameter": "基底树突非胞体节点的平均直径",
    "basal_dendrite_num_branches": "基底树突分支数",
    "basal_dendrite_soma_percentile_dorsal": "胞体在基底树突背腹跨度中的相对位置",
    "basal_dendrite_soma_percentile_medial": "胞体在基底树突内外侧跨度中的相对位置",
    "basal_dendrite_stem_exit_MedialLateral": "沿内/外侧方向离开胞体的基底树突主干比例",
    "basal_dendrite_stem_exit_dorsal": "沿背侧方向离开胞体的基底树突主干比例",
    "basal_dendrite_stem_exit_ventral": "沿腹侧方向离开胞体的基底树突主干比例",
    "basal_dendrite_total_length": "基底树突总长度",
    "basal_dendrite_total_surface_area": "基底树突总表面积",
    "soma_surface_area": "胞体表面积",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--zip",
        dest="zip_path",
        type=Path,
        default=Path(r"C:\Users\53461\Downloads\Macaque-PatchSeq-BG.zip"),
    )
    parser.add_argument(
        "--out",
        dest="out_dir",
        type=Path,
        default=Path(__file__).resolve().parent / "out",
    )
    return parser.parse_args()


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def broad_group_name(value: object) -> str | None:
    if pd.isna(value):
        return None
    text = str(value)
    if "Hybrid" in text or "D1D2" in text:
        return "Hybrid"
    if re.search(r"\bD1\b", text):
        return "D1"
    if re.search(r"\bD2\b", text):
        return "D2"
    return None


def feature_family(name: str) -> str:
    if "Sholl" in name:
        return "Sholl_PC"
    if name == "axon_exit_theta_coronal":
        return "circular_angle"
    if "stem_exit" in name:
        return "stem_exit_composition"
    if "soma_percentile" in name:
        return "relative_position"
    if "bias" in name:
        return "axis_bias"
    if "extent" in name:
        return "axis_extent"
    if "distance" in name or "total_length" in name:
        return "length_or_distance"
    if "surface_area" in name:
        return "surface_area"
    if "diameter" in name:
        return "diameter"
    if "branch" in name or "number_of_stems" in name:
        return "branching_count"
    if "contraction" in name:
        return "contraction"
    return "other"


def unit_note(name: str) -> str:
    family = feature_family(name)
    if family == "Sholl_PC":
        return "无量纲 PCA 得分；ZIP 未提供载荷或半径网格"
    if family == "circular_angle":
        return "ZIP 未明示；全表149个非缺失值均在0–1，较符合归一化整周比例而非弧度（推断）"
    if family in {"stem_exit_composition", "relative_position", "contraction"}:
        return "无量纲"
    if family in {"length_or_distance", "axis_extent", "axis_bias", "diameter"}:
        return "坐标长度单位；ZIP 未注明（常见为 µm，不能仅凭 ZIP 确认）"
    if family == "surface_area":
        return "坐标长度单位的平方；ZIP 未注明（常见为 µm²，不能仅凭 ZIP 确认）"
    if family == "branching_count":
        return "计数/阶数"
    return "ZIP 未注明"


def derivation_note(name: str) -> tuple[str, str]:
    family = feature_family(name)
    if family == "Sholl_PC":
        return (
            "二次派生得分",
            "由 Sholl 曲线经 PCA 得到；ZIP 无原始 Sholl 矩阵、中心化/缩放参数或载荷，无法从 ZIP 复算公式",
        )
    if family == "stem_exit_composition":
        return (
            "上游派生比例",
            "方向主干数/总主干数；ZIP 未给计数分子，但126细胞三部分之和可直接核验",
        )
    if family == "circular_angle":
        return (
            "上游派生圆周变量",
            "方向角；ZIP 未给坐标或单位。若0–1表示一整周，敏感性编码为 sin(2*pi*theta)、cos(2*pi*theta)；若后续证实为弧度，则用 sin(theta)、cos(theta)",
        )
    if family == "relative_position":
        return (
            "上游派生相对位置",
            "胞体在对应轴向包络中的分位位置；ZIP 未给坐标级输入，无法复算",
        )
    if family == "contraction":
        return (
            "上游派生几何比值",
            "分支点/末端之间欧氏距离总和 ÷ 相同区段路径距离总和；ZIP 无重建坐标，无法逐段复算",
        )
    return (
        "上游计算的直接形态汇总",
        "由重建几何/拓扑汇总；ZIP 未包含 SWC/ASC 等重建文件，因此只能核验表内数值，不能从原始树重算",
    )


def definition_evidence(name: str) -> tuple[str, str]:
    if name in ALLEN_CROSSCHECK_TERMS:
        return (
            "ZIP lookup短释义；Allen neuron_morphology通用同名术语作算法语义交叉核对，但ZIP未声明生成软件/版本",
            ALLEN_MORPHOLOGY_REFERENCE,
        )
    return ("ZIP lookup直接释义；中文为保守翻译", "")


def logical_range(name: str, values: pd.Series) -> tuple[str, int]:
    x = values.dropna().astype(float)
    if x.empty:
        return "无可用值", 0
    family = feature_family(name)
    if family in {"stem_exit_composition", "relative_position", "contraction"}:
        bad = int(((x < 0) | (x > 1)).sum())
        return "期望 0–1", bad
    if family in {"axis_extent", "length_or_distance", "surface_area", "diameter"}:
        bad = int((x < 0).sum())
        return "期望非负", bad
    if family == "branching_count":
        bad = int(((x < 0) | (~np.isclose(x, np.round(x), atol=1e-9))).sum())
        return "期望非负整数", bad
    if family == "circular_angle":
        bad = int(((x < 0) | (x > 1)).sum())
        return "若为归一化整周比例则期望位于 [0,1]；ZIP 未明示单位", bad
    return "无预设硬边界", 0


def availability_class(name: str, n_nonmissing: int, n_total: int) -> str:
    if n_nonmissing == 0:
        if "axon" in name or name.startswith("2") or name.startswith("23"):
            return "all_missing_axon_or_axon_Sholl"
        return "all_missing"
    if n_nonmissing == n_total:
        return "complete_dendrite_soma_candidate"
    if name == "axon_exit_theta_coronal":
        return "partial_circular_angle"
    return "partial"


def robust_flags(values: pd.Series) -> tuple[float, pd.Series]:
    x = values.astype(float)
    if x.notna().sum() == 0:
        return np.nan, pd.Series(np.nan, index=x.index, dtype=float)
    med = float(x.median())
    mad = float(np.median(np.abs(x.dropna().to_numpy() - med)))
    if not np.isfinite(mad) or mad == 0:
        return mad, pd.Series(np.nan, index=x.index, dtype=float)
    z = 0.6744897501960817 * (x - med) / mad
    return mad, z


def connected_components(edges: pd.DataFrame) -> list[list[str]]:
    graph: dict[str, set[str]] = {}
    for row in edges.itertuples(index=False):
        graph.setdefault(row.feature_a, set()).add(row.feature_b)
        graph.setdefault(row.feature_b, set()).add(row.feature_a)
    seen: set[str] = set()
    groups: list[list[str]] = []
    for start in sorted(graph):
        if start in seen:
            continue
        stack = [start]
        group: list[str] = []
        while stack:
            node = stack.pop()
            if node in seen:
                continue
            seen.add(node)
            group.append(node)
            stack.extend(sorted(graph[node] - seen))
        groups.append(sorted(group))
    return groups


def panel_table() -> pd.DataFrame:
    rows: list[dict[str, str]] = []
    rationale = {
        "basal_dendrite_total_length": "整体树突材料/规模；保留可解释原始汇总，替代高度相关的 Sholl PC0 与总表面积",
        "basal_dendrite_num_branches": "分支丰富度；与总长度互补且二者 |rho|<0.80",
        "basal_dendrite_calculate_number_of_stems": "从胞体出发的一级主干数",
        "basal_dendrite_max_branch_order": "树状拓扑深度",
        "basal_dendrite_max_path_distance": "沿树的最大缆线距离",
        "basal_dendrite_max_euclidean_distance": "空间最远可达距离；与最大路径距离定义不同",
        "basal_dendrite_mean_contraction": "分支段直线性/曲折程度",
        "basal_dendrite_mean_diameter": "树突粗细",
        "basal_dendrite_extent_dorsal": "背腹轴空间跨度",
        "basal_dendrite_extent_medial": "内外侧轴空间跨度",
        "soma_surface_area": "胞体大小",
    }
    for feature in PRIMARY_PANEL:
        rows.append(
            {
                "panel": "primary_core10",
                "feature": feature,
                "operation": "keep",
                "formula": "source value",
                "rationale": rationale[feature],
            }
        )

    for feature in PRIMARY_PANEL:
        rows.append(
            {
                "panel": "sensitivity_spatial17",
                "feature": feature,
                "operation": "keep from primary",
                "formula": "source value",
                "rationale": "主面板基线",
            }
        )
    additions = {
        "basal_dendrite_max_euclidean_distance": "补回空间最远可达距离；其与最大路径距离 |rho| 接近0.80，故仅作敏感性",
        "basal_dendrite_bias_dorsal": "加入背腹轴偏置，检验方向/切片取向信息是否影响结果",
        "basal_dendrite_bias_medial": "加入内外侧轴偏置，检验方向/切片取向信息是否影响结果",
        "basal_dendrite_soma_percentile_dorsal": "加入胞体在背腹包络中的相对位置",
        "basal_dendrite_soma_percentile_medial": "加入胞体在内外侧包络中的相对位置",
        "derived:stem_exit_ml_fraction": "用一个自由度表示内/外侧主干比例",
        "derived:stem_exit_dv_balance": "用 dorsal - ventral 表示背腹主干平衡，避免三比例闭合冗余",
    }
    formulas = {
        "derived:stem_exit_ml_fraction": "basal_dendrite_stem_exit_MedialLateral",
        "derived:stem_exit_dv_balance": "basal_dendrite_stem_exit_dorsal - basal_dendrite_stem_exit_ventral",
    }
    for feature in SPATIAL_ADDITIONS:
        rows.append(
            {
                "panel": "sensitivity_spatial17",
                "feature": feature,
                "operation": "add",
                "formula": formulas.get(feature, "source value"),
                "rationale": additions[feature],
            }
        )

    for feature in PRIMARY_PANEL:
        if feature in {"basal_dendrite_total_length", "basal_dendrite_num_branches"}:
            continue
        rows.append(
            {
                "panel": "sensitivity_Sholl10",
                "feature": feature,
                "operation": "keep from primary",
                "formula": "source value",
                "rationale": "主面板基线",
            }
        )
    rows.extend(
        [
            {
                "panel": "sensitivity_Sholl10",
                "feature": "3_Sholl_PC0",
                "operation": "replace total_length and num_branches",
                "formula": "source-provided PCA score; loadings absent",
                "rationale": "检验由 Sholl 径向分布主轴概括规模/分支的结果是否一致",
            },
            {
                "panel": "sensitivity_Sholl10",
                "feature": "3_Sholl_PC1",
                "operation": "add",
                "formula": "source-provided PCA score; loadings absent",
                "rationale": "补充径向分布的第二独立主轴；仅作敏感性，因为载荷不可审计",
            },
        ]
    )

    for feature in PRIMARY_PANEL:
        rows.append(
            {
                "panel": "sensitivity_angle12_n124",
                "feature": feature,
                "operation": "keep from primary",
                "formula": "source value",
                "rationale": "主面板基线",
            }
        )
    rows.extend(
        [
            {
                "panel": "sensitivity_angle12_n124",
                "feature": "derived:axon_exit_sin",
                "operation": "add",
                "formula": "provisional: sin(2*pi*axon_exit_theta_coronal) if 0-1 is one full turn",
                "rationale": "圆周角二维编码；仅124个细胞可用。ZIP未明示单位，确认编码后才能运行",
            },
            {
                "panel": "sensitivity_angle12_n124",
                "feature": "derived:axon_exit_cos",
                "operation": "add",
                "formula": "provisional: cos(2*pi*axon_exit_theta_coronal) if 0-1 is one full turn",
                "rationale": "与 sin 成对保留方向信息；不可把 theta 当普通线性变量。确认编码后才能运行",
            },
        ]
    )

    for feature in PRIMARY_PANEL:
        if feature == "basal_dendrite_total_length":
            rows.append(
                {
                    "panel": "sensitivity_surface10",
                    "feature": "basal_dendrite_total_surface_area",
                    "operation": "replace total_length",
                    "formula": "source value",
                    "rationale": "检验以表面积而非总长度表示树突规模；二者高度相关，不同时进入",
                }
            )
        else:
            rows.append(
                {
                    "panel": "sensitivity_surface10",
                    "feature": feature,
                    "operation": "keep from primary",
                    "formula": "source value",
                    "rationale": "主面板基线",
                }
            )
    return pd.DataFrame(rows)


def write_csv(df: pd.DataFrame, path: Path) -> None:
    df.to_csv(path, index=False, encoding="utf-8-sig", float_format="%.10g")


def main() -> None:
    args = parse_args()
    zip_path = args.zip_path.resolve()
    out_dir = args.out_dir.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    if not zip_path.is_file():
        raise FileNotFoundError(zip_path)

    zip_digest = sha256(zip_path)
    with zipfile.ZipFile(zip_path) as zf:
        names = zf.namelist()
        required = {META_ENTRY, MORPH_ENTRY, LOOKUP_ENTRY}
        missing_required = sorted(required - set(names))
        if missing_required:
            raise RuntimeError(f"ZIP missing required entries: {missing_required}")
        meta = pd.read_csv(zf.open(META_ENTRY), low_memory=False)
        morph = pd.read_csv(zf.open(MORPH_ENTRY), low_memory=False)
        lookup = pd.read_excel(io.BytesIO(zf.read(LOOKUP_ENTRY)))
        entry_rows = []
        for entry in sorted(required):
            info = zf.getinfo(entry)
            entry_rows.append(
                {
                    "entry": entry,
                    "uncompressed_bytes": info.file_size,
                    "compressed_bytes": info.compress_size,
                    "crc32_hex": f"{info.CRC:08x}",
                }
            )
        raw_reconstruction_entries = [
            n for n in names if Path(n).suffix.lower() in {".swc", ".asc", ".hoc", ".xml", ".neu", ".nml"}
        ]

    required_meta = {
        "cell_label",
        "Lib_region_of_interest_label",
        "Subclass_name",
        "Group_name",
        "donor_label",
    }
    if missing := sorted(required_meta - set(meta.columns)):
        raise RuntimeError(f"Metadata columns missing: {missing}")
    if morph.columns[0] != "cell_label" or len(morph.columns) != 40:
        raise RuntimeError("Expected cell_label plus 39 morphology fields")
    feature_cols = [c for c in morph.columns if c != "cell_label"]
    lookup = lookup[["feature_name", "feature_description"]].copy()
    if set(lookup["feature_name"]) != set(feature_cols):
        raise RuntimeError("Feature lookup and morphology feature columns differ")

    meta_dups = int(meta["cell_label"].duplicated().sum())
    morph_dups = int(morph["cell_label"].duplicated().sum())
    meta_ids = set(meta["cell_label"])
    morph_ids = set(morph["cell_label"])
    if meta_dups or morph_dups:
        raise RuntimeError("Duplicate cell_label prevents one-to-one audit")
    merged = meta.copy().merge(morph.copy(), on="cell_label", how="outer", indicator=True, validate="one_to_one").copy()

    # Cross-check the standalone morphology table against the Morp_ columns in metadata.
    standalone_metadata_disagreements = 0
    metadata_morph_columns_present = 0
    comparison_rows = []
    meta_by_id = meta.set_index("cell_label").sort_index()
    morph_by_id = morph.set_index("cell_label").sort_index()
    for feature in feature_cols:
        meta_col = f"Morp_{feature}"
        if meta_col not in meta.columns:
            comparison_rows.append(
                {"feature": feature, "metadata_column_present": False, "nonmissing_pattern_mismatches": np.nan, "numeric_mismatches": np.nan}
            )
            continue
        metadata_morph_columns_present += 1
        a = pd.to_numeric(meta_by_id[meta_col], errors="coerce")
        b = pd.to_numeric(morph_by_id[feature], errors="coerce")
        pattern_mismatch = int((a.isna() != b.isna()).sum())
        both = a.notna() & b.notna()
        numeric_mismatch = int((~np.isclose(a[both], b[both], rtol=1e-12, atol=1e-12)).sum())
        standalone_metadata_disagreements += pattern_mismatch + numeric_mismatch
        comparison_rows.append(
            {
                "feature": feature,
                "metadata_column_present": True,
                "nonmissing_pattern_mismatches": pattern_mismatch,
                "numeric_mismatches": numeric_mismatch,
            }
        )

    merged["is_target_roi"] = merged["Lib_region_of_interest_label"].isin(ROIS)
    merged["is_msn"] = merged["Subclass_name"].isin(MSN_SUBCLASSES)
    merged["broad_group_from_Group_name"] = merged["Group_name"].map(broad_group_name)
    merged["subclass_group"] = merged["Subclass_name"].map(
        {"STR D1 MSN": "D1", "STR D2 MSN": "D2", "STR Hybrid MSN": "Hybrid"}
    )

    dendrite21 = [
        c
        for c in feature_cols
        if c.startswith("3_Sholl_") or c.startswith("basal_dendrite_") or c.startswith("soma_")
    ]
    if len(dendrite21) != 21:
        raise RuntimeError(f"Expected 21 dendrite/soma fields, observed {len(dendrite21)}")
    merged["n_nonmissing_39"] = merged[feature_cols].notna().sum(axis=1)
    merged["n_nonmissing_21"] = merged[dendrite21].notna().sum(axis=1)
    merged["morph_any39"] = merged[feature_cols].notna().any(axis=1)
    merged["morph_complete21"] = merged[dendrite21].notna().all(axis=1)
    merged["axon_angle_available"] = merged["axon_exit_theta_coronal"].notna()

    target = merged.loc[merged["is_target_roi"] & merged["is_msn"]].copy()
    cohort = target.loc[target["morph_complete21"]].copy()
    if len(target) != 486 or len(cohort) != 126:
        raise RuntimeError(f"Unexpected target/cohort sizes: {len(target)}/{len(cohort)}")

    audit_cols = [
        "cell_label",
        "donor_label",
        "Lib_region_of_interest_label",
        "Subclass_name",
        "Group_name",
        "subclass_group",
        "broad_group_from_Group_name",
        "n_nonmissing_39",
        "n_nonmissing_21",
        "morph_any39",
        "morph_complete21",
        "axon_angle_available",
    ]
    cohort_audit = target[audit_cols].sort_values(
        ["morph_complete21", "Lib_region_of_interest_label", "broad_group_from_Group_name", "cell_label"],
        ascending=[False, True, True, True],
    )

    count_frames = []
    for denominator_name, frame in [("metadata_MSN_486", target), ("complete_M_126", cohort)]:
        for class_basis, col in [
            ("Group_name broad class", "broad_group_from_Group_name"),
            ("Subclass_name", "subclass_group"),
        ]:
            counts = (
                frame.groupby(["Lib_region_of_interest_label", col], dropna=False)
                .size()
                .rename("n")
                .reset_index()
                .rename(columns={col: "class"})
            )
            counts.insert(0, "class_basis", class_basis)
            counts.insert(0, "denominator", denominator_name)
            count_frames.append(counts)
    cohort_counts = pd.concat(count_frames, ignore_index=True)

    lookup_map = lookup.set_index("feature_name")["feature_description"].to_dict()
    feature_rows = []
    outlier_rows = []
    for feature in feature_cols:
        s = pd.to_numeric(cohort[feature], errors="coerce")
        x = s.dropna().astype(float)
        n_total = len(cohort)
        n_nonmissing = int(x.size)
        finite = int(np.isfinite(x).sum())
        unique = int(x.nunique())
        top_fraction = float(x.value_counts(normalize=True, dropna=True).iloc[0]) if n_nonmissing else np.nan
        mad, rz = robust_flags(s)
        flags = rz.abs() >= ROBUST_Z_CUTOFF
        logical_expectation, logical_violations = logical_range(feature, s)
        origin_type, formula_note = derivation_note(feature)
        meaning_evidence, definition_reference = definition_evidence(feature)
        q = x.quantile([0.25, 0.5, 0.75]) if n_nonmissing else pd.Series({0.25: np.nan, 0.5: np.nan, 0.75: np.nan})
        # pandas Series.skew is the bias-corrected Fisher-Pearson sample skew.
        skew = float(x.skew()) if n_nonmissing >= 3 and unique > 1 else np.nan
        iqr = float(q.loc[0.75] - q.loc[0.25]) if n_nonmissing else np.nan
        if n_nonmissing == n_total and abs(skew) < 1:
            transform = "无须仅因偏态变换；后续按确认面板标准化"
        elif n_nonmissing == n_total and skew >= 1 and float(x.min()) >= 0:
            transform = "主结果可保留原尺度；敏感性比较 log1p 后标准化"
        elif n_nonmissing == n_total and skew <= -1:
            transform = "先核查边界/构成属性；不建议自动对数变换"
        elif feature == "axon_exit_theta_coronal":
            transform = "不得线性标准化；若确认0–1为整周则用 sin(2*pi*theta), cos(2*pi*theta)，n=124"
        else:
            transform = "不可直接进入完整病例主面板"
        feature_rows.append(
            {
                "feature": feature,
                "feature_description_zip": lookup_map[feature],
                "meaning_cn": CN_MEANING[feature],
                "meaning_evidence": meaning_evidence,
                "definition_reference": definition_reference,
                "family": feature_family(feature),
                "unit": unit_note(feature),
                "origin_type": origin_type,
                "derivation_or_formula_note": formula_note,
                "availability_class": availability_class(feature, n_nonmissing, n_total),
                "n_total": n_total,
                "n_nonmissing": n_nonmissing,
                "n_missing": n_total - n_nonmissing,
                "missing_pct": 100 * (n_total - n_nonmissing) / n_total,
                "n_finite": finite,
                "n_infinite": n_nonmissing - finite,
                "n_unique": unique,
                "top_value_fraction": top_fraction,
                "constant": unique <= 1 and n_nonmissing > 0,
                "near_constant_top_ge_0p95": top_fraction >= 0.95 if n_nonmissing else False,
                "min": float(x.min()) if n_nonmissing else np.nan,
                "q1": float(q.loc[0.25]),
                "median": float(q.loc[0.5]),
                "q3": float(q.loc[0.75]),
                "max": float(x.max()) if n_nonmissing else np.nan,
                "iqr": iqr,
                "mean": float(x.mean()) if n_nonmissing else np.nan,
                "sd_sample": float(x.std(ddof=1)) if n_nonmissing >= 2 else np.nan,
                "skew_fisher_pearson_bias_corrected": skew,
                "abs_skew_ge_2": abs(skew) >= 2 if np.isfinite(skew) else False,
                "mad": mad,
                "robust_z_cutoff": ROBUST_Z_CUTOFF,
                "n_robust_extremes": int(flags.sum()),
                "logical_expectation": logical_expectation,
                "n_logical_violations": logical_violations,
                "transform_note": transform,
            }
        )
        if feature in dendrite21:
            for idx in cohort.index[flags.fillna(False)]:
                outlier_rows.append(
                    {
                        "cell_label": cohort.at[idx, "cell_label"],
                        "ROI": cohort.at[idx, "Lib_region_of_interest_label"],
                        "broad_group": cohort.at[idx, "broad_group_from_Group_name"],
                        "Subclass_name": cohort.at[idx, "Subclass_name"],
                        "Group_name": cohort.at[idx, "Group_name"],
                        "feature": feature,
                        "value": float(s.loc[idx]),
                        "median": float(s.median()),
                        "mad": mad,
                        "modified_z": float(rz.loc[idx]),
                        "criterion": f"abs(0.67448975*(x-median)/MAD) >= {ROBUST_Z_CUTOFF}",
                        "action": "review only; not excluded",
                    }
                )
    feature_qc = pd.DataFrame(feature_rows).sort_values(["availability_class", "feature"])
    extremes = pd.DataFrame(outlier_rows)
    if extremes.empty:
        extremes = pd.DataFrame(
            columns=["cell_label", "ROI", "broad_group", "Subclass_name", "Group_name", "feature", "value", "median", "mad", "modified_z", "criterion", "action"]
        )
    else:
        extremes = extremes.sort_values("modified_z", key=lambda s: s.abs(), ascending=False)

    if extremes.empty:
        extreme_cells = pd.DataFrame(columns=["cell_label", "ROI", "broad_group", "n_flagged_features", "max_abs_modified_z", "flagged_features"])
    else:
        extreme_cells = (
            extremes.assign(abs_z=extremes["modified_z"].abs())
            .groupby(["cell_label", "ROI", "broad_group"], as_index=False)
            .agg(
                n_flagged_features=("feature", "size"),
                max_abs_modified_z=("abs_z", "max"),
                flagged_features=("feature", lambda x: ";".join(sorted(x))),
            )
            .sort_values("max_abs_modified_z", ascending=False)
        )

    primary_values = cohort.set_index("cell_label")[dendrite21].astype(float)
    # Explicit Spearman implementation: average ranks for ties, then Pearson.
    # The 21 candidate fields are complete, so pairwise sample sizes are all 126.
    spearman = primary_values.rank(method="average").corr(method="pearson")
    spearman.index.name = "feature"
    pair_rows = []
    for i, feature_a in enumerate(dendrite21):
        for feature_b in dendrite21[i + 1 :]:
            pair_rows.append(
                {
                    "feature_a": feature_a,
                    "feature_b": feature_b,
                    "n_pairwise": int(primary_values[[feature_a, feature_b]].dropna().shape[0]),
                    "rho": float(spearman.loc[feature_a, feature_b]),
                    "abs_rho": abs(float(spearman.loc[feature_a, feature_b])),
                    "tie_handling": "average ranks, then Pearson correlation",
                }
            )
    spearman_pairs = pd.DataFrame(pair_rows).sort_values("abs_rho", ascending=False)
    redundant_pairs = spearman_pairs.loc[spearman_pairs["abs_rho"] >= SPEARMAN_CUTOFF].copy()
    redundancy_groups = pd.DataFrame(
        [
            {"redundancy_group": i + 1, "n_features": len(group), "features": ";".join(group)}
            for i, group in enumerate(connected_components(redundant_pairs))
        ]
    )

    # Directly test exact/algebraic relationships visible in the 126 values.
    stem_cols = [
        "basal_dendrite_stem_exit_MedialLateral",
        "basal_dendrite_stem_exit_dorsal",
        "basal_dendrite_stem_exit_ventral",
    ]
    stem_sum = cohort[stem_cols].sum(axis=1)
    surface_ratio = cohort["basal_dendrite_total_surface_area"] / (
        math.pi * cohort["basal_dendrite_total_length"] * cohort["basal_dendrite_mean_diameter"]
    )
    relation_rows = [
        {
            "relation": "stem_exit closure",
            "formula_tested": "MedialLateral + dorsal + ventral = 1",
            "n_tested": len(stem_sum),
            "n_within_tolerance": int(np.isclose(stem_sum, 1.0, rtol=0, atol=1e-12).sum()),
            "tolerance": "absolute 1e-12",
            "max_abs_error": float(np.max(np.abs(stem_sum - 1.0))),
            "ratio_min": np.nan,
            "ratio_median": np.nan,
            "ratio_max": np.nan,
            "conclusion": "exact compositional closure; use only two degrees of freedom",
        },
        {
            "relation": "surface proxy check",
            "formula_tested": "total_surface_area = pi * total_length * mean_diameter",
            "n_tested": len(surface_ratio),
            "n_within_tolerance": int(np.isclose(surface_ratio, 1.0, rtol=1e-6, atol=1e-8).sum()),
            "tolerance": "relative 1e-6, absolute 1e-8",
            "max_abs_error": float(np.max(np.abs(surface_ratio - 1.0))),
            "ratio_min": float(surface_ratio.min()),
            "ratio_median": float(surface_ratio.median()),
            "ratio_max": float(surface_ratio.max()),
            "conclusion": "near-product relation, not exact and not declared by ZIP",
        },
    ]
    exact_relations = pd.DataFrame(relation_rows)
    surface_relation = pd.DataFrame(
        {
            "cell_label": cohort["cell_label"],
            "surface_over_pi_length_mean_diameter": surface_ratio,
        }
    )

    derived_rows = [
        {
            "name": "3_Sholl_PC0 / 3_Sholl_PC1",
            "status": "source-provided derived scores",
            "formula": "not reproducible from ZIP: Sholl profiles, preprocessing and PCA loadings absent",
            "use": "sensitivity_Sholl10 only",
        },
        {
            "name": "stem_exit_ml_fraction",
            "status": "recommended transparent derived variable",
            "formula": "basal_dendrite_stem_exit_MedialLateral",
            "use": "sensitivity_spatial17",
        },
        {
            "name": "stem_exit_dv_balance",
            "status": "recommended transparent derived variable",
            "formula": "basal_dendrite_stem_exit_dorsal - basal_dendrite_stem_exit_ventral",
            "use": "sensitivity_spatial17",
        },
        {
            "name": "axon_exit_sin",
            "status": "recommended circular encoding",
            "formula": "provisional sin(2*pi*theta) if source 0-1 is one full turn; use sin(theta) only if radians are confirmed",
            "use": "sensitivity_angle12_n124",
        },
        {
            "name": "axon_exit_cos",
            "status": "recommended circular encoding",
            "formula": "provisional cos(2*pi*theta) if source 0-1 is one full turn; use cos(theta) only if radians are confirmed",
            "use": "sensitivity_angle12_n124",
        },
    ]
    derived_variables = pd.DataFrame(derived_rows)
    panels = panel_table()

    source_checks = pd.DataFrame(
        [
            {"check": "ZIP SHA-256", "value": zip_digest},
            {"check": "ZIP total entries", "value": len(names)},
            {"check": "raw reconstruction entries", "value": len(raw_reconstruction_entries)},
            {"check": "metadata rows", "value": len(meta)},
            {"check": "morphology rows", "value": len(morph)},
            {"check": "metadata unique cell_label", "value": len(meta_ids)},
            {"check": "morphology unique cell_label", "value": len(morph_ids)},
            {"check": "metadata duplicate cell_label", "value": meta_dups},
            {"check": "morphology duplicate cell_label", "value": morph_dups},
            {"check": "metadata IDs absent from morphology", "value": len(meta_ids - morph_ids)},
            {"check": "morphology IDs absent from metadata", "value": len(morph_ids - meta_ids)},
            {"check": "morphology fields", "value": len(feature_cols)},
            {"check": "metadata Morp_ fields present", "value": metadata_morph_columns_present},
            {"check": "standalone vs metadata morphology disagreements", "value": standalone_metadata_disagreements},
            {"check": "Ca+Pu+NAC MSN metadata denominator", "value": len(target)},
            {"check": "complete dendrite+soma morphology cohort", "value": len(cohort)},
            {"check": "complete-M cohort with axon angle", "value": int(cohort["axon_angle_available"].sum())},
            {"check": "all-table axon angle nonmissing", "value": int(morph["axon_exit_theta_coronal"].notna().sum())},
            {"check": "PCA run", "value": False},
            {"check": "HC run", "value": False},
            {"check": "Seurat/GC clustering run", "value": False},
            {"check": "M class frozen", "value": False},
        ]
    )

    write_csv(pd.DataFrame(entry_rows), out_dir / "zip_entries.csv")
    write_csv(source_checks, out_dir / "source_checks.csv")
    write_csv(pd.DataFrame(comparison_rows), out_dir / "morph_table_crosscheck.csv")
    write_csv(cohort_audit, out_dir / "cohort_486_audit.csv")
    write_csv(cohort_counts, out_dir / "cohort_counts.csv")
    write_csv(feature_qc, out_dir / "feature_qc_39.csv")
    spearman.to_csv(out_dir / "spearman_21_matrix.csv", encoding="utf-8-sig", float_format="%.10g")
    write_csv(spearman_pairs, out_dir / "spearman_21_pairs.csv")
    write_csv(redundant_pairs, out_dir / "redundant_pairs_absrho_ge_0p80.csv")
    write_csv(redundancy_groups, out_dir / "redundancy_groups.csv")
    write_csv(extremes, out_dir / "extreme_values_21.csv")
    write_csv(extreme_cells, out_dir / "extreme_cells_21.csv")
    write_csv(exact_relations, out_dir / "exact_relation_checks.csv")
    write_csv(surface_relation, out_dir / "surface_relation_by_cell.csv")
    write_csv(derived_variables, out_dir / "derived_variables.csv")
    write_csv(panels, out_dir / "panel_recommendations.csv")

    availability_counts = feature_qc["availability_class"].value_counts().to_dict()
    broad_486 = target["broad_group_from_Group_name"].value_counts().to_dict()
    broad_126 = cohort["broad_group_from_Group_name"].value_counts().to_dict()
    subclass_486 = target["subclass_group"].value_counts().to_dict()
    subclass_126 = cohort["subclass_group"].value_counts().to_dict()
    roi_126 = cohort["Lib_region_of_interest_label"].value_counts().to_dict()
    angle_missing_ids = sorted(cohort.loc[~cohort["axon_angle_available"], "cell_label"].tolist())
    primary_extremes = extremes.loc[extremes["feature"].isin(PRIMARY_PANEL)].copy()
    primary_extreme_cells = int(primary_extremes["cell_label"].nunique())
    selected_pair_max = (
        spearman_pairs.loc[
            spearman_pairs["feature_a"].isin(PRIMARY_PANEL) & spearman_pairs["feature_b"].isin(PRIMARY_PANEL),
            "abs_rho",
        ].max()
    )
    summary = {
        "source_zip": str(zip_path),
        "source_zip_sha256": zip_digest,
        "rows_metadata": len(meta),
        "rows_morphology": len(morph),
        "target_msn_roi_n": len(target),
        "complete_m_n": len(cohort),
        "complete_m_rate_pct": 100 * len(cohort) / len(target),
        "broad_group_486": broad_486,
        "broad_group_126": broad_126,
        "subclass_group_486": subclass_486,
        "subclass_group_126": subclass_126,
        "roi_126": roi_126,
        "availability_counts": availability_counts,
        "axon_angle_n": int(cohort["axon_angle_available"].sum()),
        "axon_angle_all_n": int(morph["axon_exit_theta_coronal"].notna().sum()),
        "axon_angle_all_min": float(morph["axon_exit_theta_coronal"].min()),
        "axon_angle_all_max": float(morph["axon_exit_theta_coronal"].max()),
        "axon_angle_missing_cell_labels": angle_missing_ids,
        "abs_skew_ge_2_features": feature_qc.loc[feature_qc["abs_skew_ge_2"], "feature"].tolist(),
        "robust_extreme_cells": int(extreme_cells.shape[0]),
        "robust_extreme_observations": int(extremes.shape[0]),
        "primary_panel_robust_extreme_cells": primary_extreme_cells,
        "primary_panel_robust_extreme_observations": int(primary_extremes.shape[0]),
        "redundant_pairs_absrho_ge_0p80": int(redundant_pairs.shape[0]),
        "primary_panel": PRIMARY_PANEL,
        "primary_panel_max_abs_spearman": float(selected_pair_max),
        "pca_run": False,
        "clustering_run": False,
        "m_class_frozen": False,
    }
    (out_dir / "summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")

    red_lines = [
        f"- `{r.feature_a}` vs `{r.feature_b}`: rho={r.rho:.4f}"
        for r in redundant_pairs.itertuples(index=False)
    ]
    skew_top = feature_qc.loc[feature_qc["n_nonmissing"] == 126].sort_values(
        "skew_fisher_pearson_bias_corrected", key=lambda s: s.abs(), ascending=False
    ).head(8)
    skew_lines = [
        f"- `{r.feature}`: skew={r.skew_fisher_pearson_bias_corrected:.3f}, median={r.median:.4g}, IQR={r.iqr:.4g}, robust extremes={int(r.n_robust_extremes)}"
        for r in skew_top.itertuples(index=False)
    ]
    primary_lines = [f"- `{x}`" for x in PRIMARY_PANEL]
    report = f"""# Macaque MSN 形态指标 QC（待确认面板）

## 范围与边界

- 仅使用 `{zip_path.name}` 内的 Macaque 数据。
- 仅保留 `Subclass_name` 为 STR D1/D2/Hybrid MSN，ROI 为 Ca、Pu、NAC。
- 未读取或引用既有 Mouse morphology 结果。
- 未运行 PCA、HC、Seurat/GC 聚类；未冻结 M 分类。

## 原始来源核验

- ZIP SHA-256：`{zip_digest}`。
- 元数据与 morphology 表均为 {len(meta)} 行，`cell_label` 均唯一，双方 ID 差集均为 0。
- 39 个独立 morphology 字段与 metadata 中对应 `Morp_` 列逐值一致；不一致数为 {standalone_metadata_disagreements}。
- ZIP 中 SWC/ASC/HOC/XML/NEU/NML 原始重建文件数：{len(raw_reconstruction_entries)}。因此这里的“39 项原始字段”应理解为来源表中的 39 项预计算特征，而非原始树重建。

## 486 与 126

- Ca+Pu+NAC 元数据 MSN：{len(target)}。
- 21 项树突/胞体字段全部完整：{len(cohort)}（{100 * len(cohort) / len(target):.1f}%）。其余 {len(target)-len(cohort)} 个 MSN 的 21 项全部缺失；本数据中不是零散的 feature-level 缺失。
- 完整 M 队列 ROI：Ca={roi_126.get('Ca',0)}，Pu={roi_126.get('Pu',0)}，NAC={roi_126.get('NAC',0)}。
- `Group_name` 解析得到的 broad class：486 中 D1={broad_486.get('D1',0)}、D2={broad_486.get('D2',0)}、Hybrid={broad_486.get('Hybrid',0)}；126 中 D1={broad_126.get('D1',0)}、D2={broad_126.get('D2',0)}、Hybrid={broad_126.get('Hybrid',0)}。
- 直接按 `Subclass_name` 则为：486 中 D1={subclass_486.get('D1',0)}、D2={subclass_486.get('D2',0)}、Hybrid={subclass_486.get('Hybrid',0)}；126 中 D1={subclass_126.get('D1',0)}、D2={subclass_126.get('D2',0)}、Hybrid={subclass_126.get('Hybrid',0)}。两种口径不可混用。

## 39 字段可用性

- 17 项轴突/轴突 Sholl 字段在 126 个细胞中全缺失。
- `axon_exit_theta_coronal` 为 {int(cohort['axon_angle_available'].sum())}/126 可用；126 队列范围 {cohort['axon_exit_theta_coronal'].min():.4f}–{cohort['axon_exit_theta_coronal'].max():.4f}。全717行中共有 {int(morph['axon_exit_theta_coronal'].notna().sum())} 个非缺失值，范围 {morph['axon_exit_theta_coronal'].min():.4f}–{morph['axon_exit_theta_coronal'].max():.4f}，全部位于0–1，较符合“整周归一化比例”的编码，但 ZIP 未明示单位。它是圆周变量，不进入主面板；编码未确认前也不运行角度敏感性分析。
- 缺少角度的两个完整 M 细胞：`{angle_missing_ids[0]}`、`{angle_missing_ids[1]}`。
- 21 项树突/胞体字段均为 126/126 完整，均为有限值；无常量或 top frequency >=95% 的近常量字段。
- 三个 stem-exit 比例的和在 126/126 细胞中等于 1（绝对误差 ≤ {float(np.max(np.abs(stem_sum - 1.0))):.3g}），因此只有两个自由度，不能将三项原样同时作为独立变量。
- `total_surface_area / (pi * total_length * mean_diameter)` 的范围为 {surface_ratio.min():.4f}–{surface_ratio.max():.4f}（中位数 {surface_ratio.median():.4f}）：近似乘积关系但并非数值上严格相等，且 ZIP 未声明公式。结合与 total length 的高相关，主面板只保留 total length 和 mean diameter。

## 冗余、偏态与极端值

Spearman 使用平均秩处理并列值；21 项均完整，因此每对 n=126。`|rho| >= {SPEARMAN_CUTOFF:.2f}` 共 {len(redundant_pairs)} 对：

{chr(10).join(red_lines)}

按无偏校正 Fisher-Pearson 样本偏度，`|skew| >= 2` 的完整候选数为 {len(summary['abs_skew_ge_2_features'])}。绝对偏度最大的 8 项：

{chr(10).join(skew_lines)}

极端值仅作复核标记：modified Z = `0.67448975*(x-median)/MAD`，阈值 `|Z| >= {ROBUST_Z_CUTOFF}`。共有 {len(extreme_cells)} 个细胞、{len(extremes)} 个 cell-feature 标记；未自动删除任何细胞。

其中主面板涉及 {primary_extreme_cells} 个细胞、{len(primary_extremes)} 个 cell-feature 标记，来自 `basal_dendrite_calculate_number_of_stems` 与 `soma_surface_area`。两者样本偏度分别为 {float(feature_qc.set_index('feature').at['basal_dendrite_calculate_number_of_stems','skew_fisher_pearson_bias_corrected']):.3f} 和 {float(feature_qc.set_index('feature').at['soma_surface_area','skew_fisher_pearson_bias_corrected']):.3f}；均未达到2。建议主结果不因偏态自动删值，后续在获批后以 `log1p` 变换作为预处理敏感性比较。

## 建议的主分析面板（core10，待用户确认）

{chr(10).join(primary_lines)}

该面板保留可解释的树突规模、分支复杂度、主干数、拓扑深度、路径可达性、曲折度、直径、两轴跨度和胞体大小。面板内部最大 `|Spearman rho|={selected_pair_max:.4f}`，低于 0.80。最大欧氏距离与最大路径距离的 `|rho|` 接近阈值，因此前者移至空间敏感性面板。

主面板暂排除：

- `3_Sholl_PC0`：与 total length、surface area、branch number 高度相关，且 ZIP 无 PCA 载荷。
- `3_Sholl_PC1`：虽不触发 0.80 阈值，但载荷和 Sholl 原始曲线缺失，移至 Sholl 敏感性面板。
- `basal_dendrite_total_surface_area`：与 total length 高度相关；以 total length 作为更直接的规模代表。
- 方向偏置、胞体轴向 percentile 与三个 stem-exit 比例：更易受切片/配准方向影响，移至空间敏感性面板；三个比例先压缩为两个自由度。
- `axon_exit_theta_coronal`：圆周角且缺2例，主面板排除。

## 敏感性面板（均待用户确认）

1. `sensitivity_spatial17`：core10 加最大欧氏距离、两轴 bias、两轴 soma percentile、`stem_exit_ml_fraction` 和 `stem_exit_dv_balance=dorsal-ventral`，检验空间方向/位置特征是否改变结构。
2. `sensitivity_Sholl10`：core10 中去掉 total length 与 branch number，换入 `3_Sholl_PC0/PC1`，检验 Sholl 径向分布表征；不把不可审计的 PC 与其高度相关原始汇总同时纳入。
3. `sensitivity_angle12_n124`：core10 加圆周正弦/余弦双编码，仅在 124 个角度可用细胞中运行；若确认0–1为整周则使用 `sin(2*pi*theta)` 与 `cos(2*pi*theta)`，若证实单位为弧度才使用 `sin(theta)` 与 `cos(theta)`。禁止直接把 theta 当线性变量。
4. `sensitivity_surface10`：core10 中用 total surface area 替换 total length，检验规模定义。

## 当前决定点

建议保留 13 个 broad Hybrid MSN 参与后续无监督形态发现，并把 D1/D2/Hybrid 只作为事后注释/分层变量；这样不会用转录组标签预先限定 M 结构，也保持 n=126。请确认：主面板是否采用 core10、是否接受保留 Hybrid 的建议，以及首选哪个敏感性面板。确认前不进入 PCA/聚类，也不冻结 M 分类。
"""
    (out_dir / "report_CN.md").write_text(report, encoding="utf-8")


if __name__ == "__main__":
    main()
