#!/usr/bin/env python
"""Audit raw Neurolucida ASC reconstructions and draw class prototypes.

File matching is year-aware. An explicit year in the master-table filename is
a hard constraint. Otherwise LastWriteTime distinguishes 2024/2025 homonyms,
and agreement with frozen raw measurements resolves the dated version. Batch
number is not used as a proxy for acquisition year.
"""

from __future__ import annotations

import hashlib
import math
import os
import re
import shutil
from pathlib import Path

import matplotlib.pyplot as plt
import neurom
import numpy as np
import pandas as pd
from neurom import features as nf


ROOT = Path(__file__).resolve().parents[1]
ASC_ROOT = Path(r"D:\ASC files")
SOURCE = ROOT / "outputs" / "morph_qc" / "hc_seurat_shift_sum13_global8_review8_215cells"
CLASS_FILE = SOURCE / "13_morph_hc_seurat_classification_215cells.csv"
RAW_FILE = SOURCE / "02_raw_13features_215cells.csv"
EC_FILE = (
    ROOT / "outputs" / "morph_qc" / "npc4_hc5derived_merged4_pseudoexpert_3nn"
    / "02_pseudoexpert_9fold_oof_predictions.csv"
)
GC_FILE = (
    ROOT / "outputs" / "morph_qc" / "seurat_npc4_res2_tsne_separation_optimization"
    / "merged_S1_S4" / "03_Morph_Seurat_MergedK4_cell_assignments_coordinates.csv"
)
Z_FILE = (
    ROOT / "outputs" / "morph_qc" / "postprocessed_pearson_redundancy_then_skew_core11"
    / "05_core11_final_zscore_215cells.csv"
)
OUT = ROOT / "outputs" / "morph_qc" / "asc_reconstruction_qc_and_classic_examples"
INDIVIDUAL = OUT / "individual_examples"

ID = "MSN_unique_ID"
HC_TO_M = {"HC1": "Morph-M1", "HC2": "Morph-M1", "HC3": "Morph-M2", "HC4": "Morph-M3", "HC5": "Morph-M4"}
GC_TO_M = {"S-1": "Morph-M1", "S-2": "Morph-M3", "S-3": "Morph-M2", "S-4": "Morph-M4"}
COLORS = {
    "Morph-M1": "#00468B",
    "Morph-M2": "#42B540",
    "Morph-M3": "#ED0000",
    "Morph-M4": "#0099B4",
}


def date_key(text: str) -> tuple[int, int, int] | None:
    nums = [int(x) for x in re.findall(r"\d+", Path(str(text)).stem)]
    nums = [x for x in nums if x not in (2024, 2025)]
    if len(nums) < 3:
        return None
    return tuple(nums[:3])


def explicit_year(text: str) -> int | None:
    match = re.search(r"(?:^|[^0-9])(2024|2025)(?:[^0-9]|$)", str(text))
    return int(match.group(1)) if match else None


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def safe_float(value, default=np.nan) -> float:
    try:
        return float(value)
    except Exception:
        return float(default)


def branch_orders(neuron) -> list[int]:
    values: list[int] = []
    for neurite in neuron.neurites:
        stack = [(neurite.root_node, 1)]
        while stack:
            section, order = stack.pop()
            values.append(order)
            for child in section.children:
                stack.append((child, order + 1))
    return values


def parse_asc(path: Path) -> tuple[dict, object | None]:
    base = {
        "ASC_parse_success": False,
        "ASC_parse_error": "",
        "ASC_soma_present": False,
        "ASC_soma_points": 0,
        "ASC_neurites": np.nan,
        "ASC_sections": np.nan,
        "ASC_points": np.nan,
        "ASC_bifurcation_sections": np.nan,
        "ASC_terminal_sections": np.nan,
        "ASC_max_branch_order": np.nan,
        "ASC_total_length_um": np.nan,
        "ASC_max_radial_distance_um": np.nan,
        "ASC_width_um": np.nan,
        "ASC_height_um": np.nan,
        "ASC_depth_um": np.nan,
        "ASC_soma_radius_um": np.nan,
        "ASC_root_gap_max_um": np.nan,
        "ASC_zero_length_segments": np.nan,
        "ASC_nonfinite_values": np.nan,
        "ASC_nonpositive_radii": np.nan,
        "ASC_segment_median_um": np.nan,
        "ASC_segment_p99_um": np.nan,
        "ASC_segment_max_um": np.nan,
        "ASC_long_jump_segments": np.nan,
        "ASC_incomplete_markers": np.nan,
    }
    try:
        neuron = neurom.load_morphology(str(path))
        soma_points = np.asarray(neuron.soma.points, dtype=float)
        sections = list(neuron.sections)
        all_points = [np.asarray(s.points, dtype=float) for s in sections]
        seg_lengths: list[float] = []
        nonfinite = 0
        nonpositive_radii = 0
        for points in all_points:
            nonfinite += int((~np.isfinite(points)).sum())
            if points.shape[1] >= 4:
                nonpositive_radii += int(np.sum(points[:, 3] <= 0))
            if len(points) > 1:
                seg_lengths.extend(np.linalg.norm(np.diff(points[:, :3], axis=0), axis=1).tolist())
        seg = np.asarray(seg_lengths, dtype=float)
        positive = seg[seg > 1e-9]
        median_seg = float(np.median(positive)) if positive.size else np.nan
        long_threshold = max(30.0, 10.0 * median_seg) if np.isfinite(median_seg) else 30.0
        center = np.asarray(neuron.soma.center, dtype=float)
        soma_radius = safe_float(neuron.soma.radius)
        root_gaps = []
        for neurite in neuron.neurites:
            root = np.asarray(neurite.root_node.points[0, :3], dtype=float)
            root_gaps.append(max(0.0, float(np.linalg.norm(root - center)) - soma_radius))
        text = path.read_text(encoding="utf-8", errors="ignore")
        incomplete = len(re.findall(r"\b(?:incomplete|low|high)\b", text, flags=re.I))
        orders = branch_orders(neuron)

        base.update({
            "ASC_parse_success": True,
            "ASC_soma_present": bool(len(soma_points) >= 3),
            "ASC_soma_points": int(len(soma_points)),
            "ASC_neurites": int(len(neuron.neurites)),
            "ASC_sections": int(len(sections)),
            "ASC_points": int(sum(len(x) for x in all_points)),
            "ASC_bifurcation_sections": int(sum(len(s.children) > 0 for s in sections)),
            "ASC_terminal_sections": int(sum(len(s.children) == 0 for s in sections)),
            "ASC_max_branch_order": int(max(orders)) if orders else 0,
            "ASC_total_length_um": float(sum(nf.get("total_length_per_neurite", neuron))),
            "ASC_max_radial_distance_um": safe_float(nf.get("max_radial_distance", neuron)),
            "ASC_width_um": safe_float(nf.get("total_width", neuron)),
            "ASC_height_um": safe_float(nf.get("total_height", neuron)),
            "ASC_depth_um": safe_float(nf.get("total_depth", neuron)),
            "ASC_soma_radius_um": soma_radius,
            "ASC_root_gap_max_um": float(max(root_gaps)) if root_gaps else np.nan,
            "ASC_zero_length_segments": int(np.sum(seg <= 1e-9)),
            "ASC_nonfinite_values": nonfinite,
            "ASC_nonpositive_radii": nonpositive_radii,
            "ASC_segment_median_um": median_seg,
            "ASC_segment_p99_um": float(np.quantile(seg, 0.99)) if seg.size else np.nan,
            "ASC_segment_max_um": float(np.max(seg)) if seg.size else np.nan,
            "ASC_long_jump_segments": int(np.sum(seg > long_threshold)),
            "ASC_incomplete_markers": incomplete,
        })
        return base, neuron
    except Exception as exc:
        base["ASC_parse_error"] = f"{type(exc).__name__}: {exc}"
        return base, None


def candidate_score(parsed: dict, row: pd.Series) -> tuple[float, float, float]:
    expected_length = safe_float(row["M_Total_neurite_length_(sections)"])
    expected_radial = safe_float(row["M_cell_max_radial_dist"])
    if not parsed["ASC_parse_success"]:
        return math.inf, math.inf, math.inf
    rel_length = abs(parsed["ASC_total_length_um"] - expected_length) / max(abs(expected_length), 1.0)
    rel_radial = abs(parsed["ASC_max_radial_distance_um"] - expected_radial) / max(abs(expected_radial), 1.0)
    return 4.0 * rel_length + rel_radial, rel_length, rel_radial


def build_inventory() -> pd.DataFrame:
    records = []
    for path in sorted(ASC_ROOT.glob("*.asc"), key=lambda p: p.name.lower()):
        stat = path.stat()
        records.append({
            "ASC_filename": path.name,
            "ASC_path": str(path),
            "Date_key": "-".join(map(str, date_key(path.name))) if date_key(path.name) else "",
            "Modification_time": pd.Timestamp(stat.st_mtime, unit="s"),
            "Modification_year": pd.Timestamp(stat.st_mtime, unit="s").year,
            "File_size_bytes": stat.st_size,
            "Contains_trace": bool(re.search(r"trace", path.name, flags=re.I)),
            "Contains_at_symbol": "@" in path.name,
            "SHA256": sha256(path),
        })
    return pd.DataFrame(records)


def match_cells(cells: pd.DataFrame, inventory: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, object]]:
    parsed_cache: dict[str, tuple[dict, object | None]] = {}
    output = []
    for _, row in cells.iterrows():
        key = date_key(row["M_Neuron_id"])
        year = explicit_year(row["M_Neuron_id"])
        key_text = "-".join(map(str, key)) if key else ""
        same_key = inventory[inventory["Date_key"] == key_text]
        candidates = same_key if year is None else same_key[same_key["Modification_year"] == year]
        record = {
            ID: row[ID],
            "M_Neuron_id": row["M_Neuron_id"],
            "T_Batch": row["T_Batch"],
            "Explicit_filename_year": year,
            "Date_key": key_text,
            "All_same_key_files": " | ".join(same_key["ASC_filename"].tolist()),
            "Allowed_candidate_count": len(candidates),
            "Selected_ASC_filename": "",
            "Selected_ASC_path": "",
            "Selected_modification_time": pd.NaT,
            "Selected_modification_year": np.nan,
            "Year_match": False if year is not None else np.nan,
            "Mapping_status": "Missing_matching_ASC",
            "Candidate_metric_score": np.nan,
            "Length_relative_error": np.nan,
            "Radial_relative_error": np.nan,
        }
        if len(candidates) == 0:
            output.append(record)
            continue
        scored = []
        for _, candidate in candidates.iterrows():
            path = Path(candidate["ASC_path"])
            if str(path) not in parsed_cache:
                parsed_cache[str(path)] = parse_asc(path)
            parsed, _ = parsed_cache[str(path)]
            score, len_err, rad_err = candidate_score(parsed, row)
            scored.append((score, len_err, rad_err, candidate, parsed))
        scored.sort(key=lambda x: (x[0], x[3]["ASC_filename"].lower()))
        score, len_err, rad_err, chosen, parsed = scored[0]
        if not parsed["ASC_parse_success"]:
            status = "Parse_failure"
        elif len_err <= 0.005:
            status = "Verified_year_and_length"
        elif len_err <= 0.02:
            status = "Plausible_year_and_length"
        else:
            status = "Review_metric_mismatch"
        record.update({
            "Selected_ASC_filename": chosen["ASC_filename"],
            "Selected_ASC_path": chosen["ASC_path"],
            "Selected_modification_time": chosen["Modification_time"],
            "Selected_modification_year": chosen["Modification_year"],
            "Year_match": bool(chosen["Modification_year"] == year) if year is not None else np.nan,
            "Mapping_status": status,
            "Candidate_metric_score": score,
            "Length_relative_error": len_err,
            "Radial_relative_error": rad_err,
            **parsed,
        })
        output.append(record)
    return pd.DataFrame(output), parsed_cache


def add_qc_status(df: pd.DataFrame) -> pd.DataFrame:
    statuses = []
    flags_all = []
    for _, row in df.iterrows():
        hard = []
        review = []
        if row["Mapping_status"] == "Missing_matching_ASC":
            hard.append("missing_matching_ASC")
        elif row["Mapping_status"] == "Parse_failure" or not bool(row.get("ASC_parse_success", False)):
            hard.append("ASC_parse_failure")
        else:
            if not bool(row["ASC_soma_present"]):
                hard.append("soma_missing_or_lt3_points")
            if row["ASC_neurites"] < 1 or row["ASC_sections"] < 1:
                hard.append("no_neurite_tree")
            if row["ASC_nonfinite_values"] > 0:
                hard.append("nonfinite_coordinates")
            if row["Mapping_status"] == "Review_metric_mismatch":
                review.append("table_length_mismatch_gt2pct")
            if row["ASC_nonpositive_radii"] > 0:
                review.append("nonpositive_point_radius")
            if row["ASC_zero_length_segments"] > 0:
                review.append("zero_length_segment")
            if row["ASC_long_jump_segments"] > 0:
                review.append("long_coordinate_jump")
            if row["ASC_incomplete_markers"] > 0:
                review.append("incomplete_terminal_marker")
            root_limit = max(10.0, 2.0 * safe_float(row["ASC_soma_radius_um"], 0.0))
            if row["ASC_root_gap_max_um"] > root_limit:
                review.append("neurite_root_far_from_soma")
            if row["ASC_soma_points"] < 5:
                review.append("soma_contour_lt5_points")
        if hard:
            statuses.append("Unavailable_or_hard_fail")
        elif review:
            statuses.append("Review")
        else:
            statuses.append("Pass")
        flags_all.append(";".join(hard + review))
    result = df.copy()
    result["ASC_QC_status"] = statuses
    result["ASC_QC_flags"] = flags_all
    return result


def select_classics(df: pd.DataFrame, z: pd.DataFrame) -> pd.DataFrame:
    feature_cols = [c for c in z.columns if c != ID]
    z_renamed = z.rename(columns={c: f"Z__{c}" for c in feature_cols})
    z_cols = [f"Z__{c}" for c in feature_cols]
    merged = df.merge(z_renamed, on=ID, how="left", validate="one_to_one")
    selected = []
    for label in COLORS:
        group = merged[
            (merged["Triple_consensus"])
            & (merged["M_consensus_class"] == label)
            & (merged["ASC_QC_status"] == "Pass")
            & (merged["Mapping_status"] == "Verified_year_and_length")
        ].copy()
        matrix = group[z_cols].to_numpy(float)
        center = np.nanmedian(matrix, axis=0)
        # Robust distance to the multivariate class centre: low = archetypal.
        mad = np.nanmedian(np.abs(matrix - center), axis=0)
        mad[mad < 0.15] = 0.15
        group["Class_representativeness_distance"] = np.sqrt(np.sum(((matrix - center) / mad) ** 2, axis=1))
        # Penalize small residual table/reconstruction discrepancies without
        # allowing QC to replace biological representativeness.
        group["Classic_selection_score"] = (
            group["Class_representativeness_distance"]
            + 5.0 * group["Length_relative_error"].fillna(1.0)
            + 1.0 * group["Radial_relative_error"].fillna(1.0)
        )
        group = group.sort_values(["Classic_selection_score", ID]).head(5).copy()
        if len(group) < 5:
            raise RuntimeError(f"Only {len(group)} fully QC-passing classic candidates for {label}")
        group["Rank_within_class"] = np.arange(1, 6)
        selected.append(group)
    return pd.concat(selected, ignore_index=True)


def draw_neuron(ax, neuron, color: str, center: np.ndarray, half_span: float, title: str | None = None):
    for section in neuron.sections:
        points = np.asarray(section.points, dtype=float)
        ax.plot(points[:, 0] - center[0], points[:, 1] - center[1], color=color,
                lw=0.65, solid_capstyle="round", solid_joinstyle="round")
    soma = np.asarray(neuron.soma.points, dtype=float)
    if len(soma) >= 3:
        ax.fill(soma[:, 0] - center[0], soma[:, 1] - center[1], color=color, alpha=1.0, linewidth=0)
    ax.set_xlim(-half_span, half_span)
    ax.set_ylim(-half_span, half_span)
    ax.set_aspect("equal", adjustable="box")
    ax.axis("off")
    if title:
        ax.set_title(title, fontsize=6, pad=1.5)


def render_examples(selected: pd.DataFrame, cache: dict[str, tuple[dict, object | None]]):
    loaded = {}
    cell_extrema = {}
    for _, row in selected.iterrows():
        path = row["Selected_ASC_path"]
        neuron = cache[path][1] if path in cache else neurom.load_morphology(path)
        loaded[row[ID]] = neuron
        center = np.asarray(neuron.soma.center, dtype=float)
        extrema = []
        for section in neuron.sections:
            xy = np.asarray(section.points[:, :2], dtype=float) - center[:2]
            extrema.extend(np.abs(xy).ravel().tolist())
        cell_extrema[row[ID]] = max(extrema) if extrema else 50.0
    # Use the actual maximum, never a quantile, so no reconstruction is clipped.
    half_span = max(50.0, math.ceil(max(cell_extrema.values()) / 10.0) * 10.0) * 1.05

    # Publication-scale requirement: 50 um must occupy exactly 0.2 inch in the
    # final PNG/PDF. Because every axis spans 2 * half_span um, this fixes the
    # physical square plotting area. Axes are positioned explicitly so figure
    # margins and layout engines cannot change the micrometre-to-inch mapping.
    scale_bar_um = 50.0
    scale_bar_in = 0.2
    axis_size_in = (2.0 * half_span / scale_bar_um) * scale_bar_in
    panel_gap_in = 0.04
    left_margin_in = 0.04
    right_margin_in = 0.04
    bottom_margin_in = 0.04
    title_space_in = 0.50
    combined_width_in = (
        left_margin_in + 5.0 * axis_size_in + 4.0 * panel_gap_in + right_margin_in
    )
    combined_height_in = bottom_margin_in + axis_size_in + title_space_in

    for label, group in selected.groupby("M_consensus_class", sort=False):
        group = group.sort_values("Rank_within_class")
        fig = plt.figure(figsize=(combined_width_in, combined_height_in))
        axes = []
        for panel_index in range(5):
            x_in = left_margin_in + panel_index * (axis_size_in + panel_gap_in)
            axes.append(
                fig.add_axes(
                    [
                        x_in / combined_width_in,
                        bottom_margin_in / combined_height_in,
                        axis_size_in / combined_width_in,
                        axis_size_in / combined_height_in,
                    ]
                )
            )
        for ax, (_, row) in zip(axes, group.iterrows()):
            neuron = loaded[row[ID]]
            center = np.asarray(neuron.soma.center, dtype=float)
            draw_neuron(ax, neuron, COLORS[label], center, half_span, str(row[ID]).replace("_Batch", " B"))
        # One common 50-um scale bar; it is exactly 0.2 inch in the saved file.
        x0 = -half_span * 0.82
        y0 = -half_span * 0.84
        axes[0].plot(
            [x0, x0 + scale_bar_um], [y0, y0],
            color="black", lw=1.1, solid_capstyle="butt"
        )
        axes[0].text(
            x0 + scale_bar_um / 2.0, y0 + half_span * 0.06, "50 µm",
            ha="center", va="bottom", fontsize=5
        )
        fig.suptitle(
            f"{label}: five archetypal triple-consensus reconstructions",
            fontsize=8, color=COLORS[label], y=0.985
        )
        stem = OUT / f"07_{label.replace('Morph-', '')}_five_classic_ASC_reconstructions"
        fig.savefig(stem.with_suffix(".png"), dpi=600, facecolor="white")
        fig.savefig(stem.with_suffix(".pdf"), facecolor="white")
        plt.close(fig)

        for _, row in group.iterrows():
            fig = plt.figure(figsize=(axis_size_in, axis_size_in))
            ax = fig.add_axes([0, 0, 1, 1])
            neuron = loaded[row[ID]]
            center = np.asarray(neuron.soma.center, dtype=float)
            draw_neuron(ax, neuron, COLORS[label], center, half_span)
            x0 = -half_span * 0.82
            y0 = -half_span * 0.84
            ax.plot(
                [x0, x0 + scale_bar_um], [y0, y0],
                color="black", lw=1.1, solid_capstyle="butt"
            )
            stem = INDIVIDUAL / f"{label}_{int(row['Rank_within_class']):02d}_{row[ID]}"
            fig.savefig(stem.with_suffix(".png"), dpi=600, facecolor="white")
            fig.savefig(stem.with_suffix(".pdf"), facecolor="white")
            plt.close(fig)
    return half_span, axis_size_in, combined_width_in, combined_height_in


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    if INDIVIDUAL.exists():
        shutil.rmtree(INDIVIDUAL)
    INDIVIDUAL.mkdir(parents=True, exist_ok=True)
    labels = pd.read_csv(CLASS_FILE)
    raw = pd.read_csv(RAW_FILE)
    ec = pd.read_csv(EC_FILE)[[ID, "OOF_prediction"]]
    gc = pd.read_csv(GC_FILE)[[ID, "Merged_class"]]
    z = pd.read_csv(Z_FILE)
    cells = labels[[ID, "M_Neuron_id", "T_Batch", "D1D2_Major"]].merge(raw, on=ID, validate="one_to_one")
    cells = cells.merge(ec, on=ID, validate="one_to_one").merge(gc, on=ID, validate="one_to_one")
    cells["HC_class"] = labels.set_index(ID).loc[cells[ID], "Morph_HC4"].to_numpy()
    # Use the frozen NPC4 HC K5 teacher labels for the M1-M4 correspondence.
    hc5 = pd.read_csv(
        ROOT / "outputs" / "morph_qc" / "core11_npc4_hc4_5_6_7_8_9_seurat_res05_30"
        / "04_NPC4_HC_Seurat_cell_assignments.csv"
    )[[ID, "HC_NPC4_K5"]]
    cells = cells.merge(hc5, on=ID, validate="one_to_one")
    cells["M_consensus_class"] = cells["OOF_prediction"]
    cells["Triple_consensus"] = (
        cells["OOF_prediction"] == cells["HC_NPC4_K5"].map(HC_TO_M)
    ) & (
        cells["OOF_prediction"] == cells["Merged_class"].map(GC_TO_M)
    )
    if len(cells) != 215 or int(cells["Triple_consensus"].sum()) != 189:
        raise RuntimeError("Frozen 215-cell cohort or 189-cell triple consensus was not reproduced")

    inventory = build_inventory()
    mapping_qc, cache = match_cells(cells, inventory)
    combined = cells.merge(mapping_qc, on=[ID, "M_Neuron_id", "T_Batch"], validate="one_to_one")
    combined = add_qc_status(combined)

    inventory.to_csv(OUT / "01_ASC_file_inventory_with_modification_time.csv", index=False)
    mapping_qc.to_csv(OUT / "02_ASC_year_aware_cell_mapping_audit_215cells.csv", index=False)
    combined.to_csv(OUT / "03_ASC_reconstruction_QC_215cells.csv", index=False)

    summary_rows = []
    for field in ["Mapping_status", "ASC_QC_status"]:
        counts = combined[field].fillna("Missing").value_counts()
        for category, count in counts.items():
            summary_rows.append({"Summary_dimension": field, "Category": category, "Cell_count": int(count)})
    flag_counts = combined["ASC_QC_flags"].fillna("").str.split(";").explode()
    flag_counts = flag_counts[flag_counts != ""].value_counts()
    for category, count in flag_counts.items():
        summary_rows.append({"Summary_dimension": "ASC_QC_flag", "Category": category, "Cell_count": int(count)})
    pd.DataFrame(summary_rows).to_csv(OUT / "04_ASC_QC_summary_counts.csv", index=False)

    selected = select_classics(combined, z)
    selected_columns = [
        "M_consensus_class", "Rank_within_class", ID, "D1D2_Major", "M_Neuron_id",
        "Selected_ASC_filename", "Selected_modification_time", "Explicit_filename_year",
        "Mapping_status", "ASC_QC_status", "Class_representativeness_distance",
        "Classic_selection_score", "ASC_neurites", "ASC_sections", "ASC_total_length_um",
        "ASC_max_radial_distance_um", "ASC_bifurcation_sections", "ASC_terminal_sections",
        "ASC_max_branch_order", "ASC_root_gap_max_um", "Length_relative_error",
        "Radial_relative_error",
    ]
    selected[selected_columns].to_csv(OUT / "05_M1_M4_five_classic_cells_per_class.csv", index=False)
    half_span, axis_size_in, combined_width_in, combined_height_in = render_examples(selected, cache)

    methods = f"""ASC reconstruction QC and classic-cell selection

Raw directory: {ASC_ROOT}
Frozen analysis cohort: 215 cells
Triple consensus definition: EC OOF pseudoexpert = NPC4 HC K5-derived M class = merged Seurat GC class
Triple-consensus cells: {int(combined['Triple_consensus'].sum())}

Year-aware matching:
- An explicit 2024/2025 token in the table filename is a hard LastWriteTime-year constraint.
- If the table filename has no explicit year, all matching date keys are evaluated; the
  2024/2025 LastWriteTime separates homonyms and the frozen morphology values identify
  which dated version belongs to the cell. Batch number is not used to infer year.
- Date key = month-day-cell index after removing an explicit 2024/2025 token.
- Candidates are ranked by agreement with frozen total neurite length and maximum
  radial distance; total length receives fourfold weight.

ASC-derived QC fields include parser success, soma contour presence/point count,
neurite roots, sections, traced points, bifurcation/terminal sections, maximum branch
order, total length, radial distance, width/height/depth, soma radius, root-to-soma
gap, zero-length segments, non-finite values, non-positive radii, segment length
distribution, long coordinate jumps, and incomplete terminal markers.

Classic examples:
- restricted to EC-HC-GC triple-consensus cells
- restricted to year-and-length verified ASC files with ASC QC status Pass
- ranked by robust multivariate distance to the class median in the frozen core11
  postprocessed z-score space, with small penalties for reconstruction/table residuals
- five lowest-scoring cells retained per Morph-M1 through Morph-M4
- all example panels share one absolute XY scale: ±{half_span:.1f} µm about soma centre
- physical rendering scale: 50 µm = 0.2 inch in every PNG and PDF
- square plotting area per neuron: {axis_size_in:.3f} x {axis_size_in:.3f} inches
- five-cell plate size: {combined_width_in:.3f} x {combined_height_in:.3f} inches
- native orientation retained; no rotation, mirroring, or size normalization
"""
    (OUT / "06_ASC_QC_and_classic_selection_methods.txt").write_text(methods, encoding="utf-8")
    print(pd.crosstab(combined["Mapping_status"], combined["ASC_QC_status"], margins=True))
    print(selected[selected_columns].to_string(index=False))
    print(f"Common plot half-span: {half_span:.1f} um")
    print(f"Physical scale: 50 um = 0.2 inch; axis size = {axis_size_in:.3f} inch")
    print(f"Five-cell plate: {combined_width_in:.3f} x {combined_height_in:.3f} inch")
    print(OUT)


if __name__ == "__main__":
    main()
