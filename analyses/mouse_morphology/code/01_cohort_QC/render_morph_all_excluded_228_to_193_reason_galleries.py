#!/usr/bin/env python
"""Render all 35 cells excluded/quarantined between Morph n=228 and n=193.

The output is an audit, not an assertion that every cell is biologically invalid.
Five cells are frozen upstream exclusions, eight are user-confirmed numerical
review exclusions, and twenty-two are reversible morphology-only quarantines.
All ASC reconstructions retain native XY orientation and one absolute scale.
"""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from plot_all_frozen_morph_classes_one_page_each import (  # noqa: E402
    FOOTER_IN,
    HEADER_IN,
    ID,
    INCH_PER_UM,
    LABEL_HEIGHT_IN,
    PAGE_MARGIN_IN,
    SCALE_BAR_IN,
    SCALE_BAR_UM,
    draw_neuron,
    load_item,
    shortened_id,
)
from render_HC_G_matrix_cell_ASC_galleries import (  # noqa: E402
    balanced_pack,
    build_inventory,
    locate_asc_root,
    map_cells,
)


MASTER_JSON = ROOT / ".codex-work/e_type_qc/stage1_filtered.json"
EXCLUSION_AUDIT = ROOT / (
    "outputs/morph_qc/hc_seurat_shift_sum13_global8_review8_215cells/"
    "00_exclusions_applied.csv"
)
OUTLIER_REVIEW = ROOT / (
    "outputs/morph_qc/initial_qc_skew_redundancy/outlier_review_by_cell.json"
)
QUARANTINE_JSON = ROOT / "config/morph_current_reconstruction_spe_local_pyramid_quarantine.json"
OUT = ROOT / "outputs/morph_qc/all35_excluded_228_to_193_ASC_reason_galleries"
INDIVIDUAL = OUT / "individual_900dpi"

PAGE_PNG_DPI = 600
INDIVIDUAL_PNG_DPI = 900
LABEL_BLOCK_IN = 0.34
MAX_TILES_PER_PAGE = 6

STAGE_COLORS = {
    "Permanent upstream exclusion": "#6A3D9A",
    "Morph numerical-review exclusion": "#E31A1C",
    "Temporary reconstruction/ASC quarantine": "#FF7F00",
    "Temporary morphology-extreme quarantine": "#1F78B4",
}


def configure_plotting() -> None:
    mpl.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 4,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
            "savefig.facecolor": "white",
            "savefig.pad_inches": 0,
        }
    )


def load_master_raw() -> pd.DataFrame:
    payload = json.loads(MASTER_JSON.read_text(encoding="utf-8"))
    raw = pd.DataFrame(payload["raw_rows"], columns=payload["raw_headers"])
    required = [
        ID,
        "has_M",
        "M_ID",
        "M_Neuron_id",
        "M_Total_neurite_length_(sections)",
        "M_cell_max_radial_dist",
    ]
    missing = [column for column in required if column not in raw.columns]
    if missing:
        raise RuntimeError(f"Master raw table is missing: {missing}")
    return raw


def load_stage_lists() -> tuple[list[str], list[str], dict]:
    audit = pd.read_csv(EXCLUSION_AUDIT)
    permanent = audit.loc[
        audit["Reason"].eq("Frozen global exclusion")
        & audit["Present_in_raw_complete_Morph"].astype(bool),
        "Exclusion",
    ].tolist()
    manual = audit.loc[
        audit["Reason"].str.contains("Morph manual-review", na=False)
        & audit["Present_in_raw_complete_Morph"].astype(bool),
        "Exclusion",
    ].tolist()
    quarantine = json.loads(QUARANTINE_JSON.read_text(encoding="utf-8"))
    if len(permanent) != 5 or len(manual) != 8:
        raise RuntimeError(f"Expected 5 permanent and 8 manual exclusions, got {len(permanent)}, {len(manual)}")
    return permanent, manual, quarantine


def quarantine_reason_map(config: dict) -> dict[str, list[str]]:
    reasons: dict[str, list[str]] = {}

    def add(cell: str, reason: str) -> None:
        reasons.setdefault(cell, []).append(reason)

    suspected = config["suspected_reconstruction_failure"]
    for cell in suspected["possible_under_reconstruction"]:
        add(cell, "possible under-reconstruction")
    for cell in suspected["direct_ASC_technical_review"]:
        add(cell, "direct ASC technical review")
    for cell in suspected["ASC_parse_failure"]:
        add(cell, "ASC parse failure")
    for cell in config["SPE_extreme"]:
        add(cell, "S/P/E extreme morphology")
    for cell in config["strict_local_multivariate_outlier"]:
        add(cell, "strict local multivariate outlier")
    for cell in config["selected_pyramidal_like_candidates"]["cells"]:
        add(cell, "selected pyramidal-like candidate")
    return reasons


def build_manifest() -> pd.DataFrame:
    permanent, manual, quarantine = load_stage_lists()
    review_rows = json.loads(OUTLIER_REVIEW.read_text(encoding="utf-8"))
    review = {row[ID]: row for row in review_rows}
    q_reasons = quarantine_reason_map(quarantine)

    rows = []
    for cell in permanent:
        rows.append(
            {
                ID: cell,
                "Exclusion_stage": "Permanent upstream exclusion",
                "Primary_reason": "frozen pre-QC/global exclusion",
                "All_reasons": "frozen pre-QC/global exclusion; not morphology-shape adjudicated here",
                "Status": "permanently excluded",
            }
        )
    for cell in manual:
        source = review[cell]
        rows.append(
            {
                ID: cell,
                "Exclusion_stage": "Morph numerical-review exclusion",
                "Primary_reason": f"robust multivariate feature outlier ({source['Full22_flag_count']} flags)",
                "All_reasons": (
                    f"user-confirmed review exclusion; flagged features: "
                    f"{source['Flagged_features_full22']}; "
                    f"maximum |robust z|={float(source['Max_abs_robust_z_full22']):.2f}"
                ),
                "Status": "excluded from frozen 215-cell Morph cohort",
            }
        )

    technical_terms = {
        "possible under-reconstruction",
        "direct ASC technical review",
        "ASC parse failure",
    }
    for cell, reasons in q_reasons.items():
        if any(reason in technical_terms for reason in reasons):
            stage = "Temporary reconstruction/ASC quarantine"
        else:
            stage = "Temporary morphology-extreme quarantine"
        rows.append(
            {
                ID: cell,
                "Exclusion_stage": stage,
                "Primary_reason": reasons[0],
                "All_reasons": "; ".join(reasons),
                "Status": "temporary reversible quarantine",
            }
        )

    manifest = pd.DataFrame(rows).drop_duplicates(ID, keep="first")
    stage_order = list(STAGE_COLORS)
    manifest["Stage_order"] = manifest["Exclusion_stage"].map({value: i for i, value in enumerate(stage_order)})
    manifest = manifest.sort_values(["Stage_order", "Primary_reason", ID]).drop(columns="Stage_order").reset_index(drop=True)
    if len(manifest) != 35 or manifest[ID].nunique() != 35:
        raise RuntimeError(f"Expected 35 unique excluded/quarantined cells, found {len(manifest)}")
    return manifest


def map_all_cells(manifest: pd.DataFrame) -> pd.DataFrame:
    raw = load_master_raw()
    columns = [
        ID,
        "M_ID",
        "M_Neuron_id",
        "M_Total_neurite_length_(sections)",
        "M_cell_max_radial_dist",
    ]
    cells = manifest.merge(raw[columns], on=ID, how="left", validate="one_to_one")
    if cells[columns[1:]].isna().any().any():
        missing = cells.loc[cells[columns[1:]].isna().any(axis=1), ID].tolist()
        raise RuntimeError(f"Missing morphology metadata for: {missing}")

    # Compatibility columns used only by the year-aware mapping function.
    cells["HC_cluster"] = "not_used"
    cells["GC_cluster"] = "not_used"
    cells["PC1"] = 0.0
    cells["PC2"] = 0.0
    cells["PC3"] = 0.0

    asc_root = locate_asc_root()
    inventory = build_inventory(asc_root)
    mapped = map_cells(cells, inventory)
    mapped = cells.drop(columns=["HC_cluster", "GC_cluster", "PC1", "PC2", "PC3"]).merge(
        mapped.drop(columns=["HC_cluster", "GC_cluster", "PC1", "PC2", "PC3", "M_ID", "M_Neuron_id"]),
        on=ID,
        how="left",
        validate="one_to_one",
    )
    # For this visual QC audit, a same-date/year parseable file remains visible
    # even when its scalar metrics disagree; the mismatch is reported, not hidden.
    mapped["ASC_drawable_for_audit"] = (
        mapped["ASC_parse_success"].fillna(False).astype(bool)
        & mapped["Selected_ASC_path"].fillna("").map(lambda value: Path(str(value)).is_file())
    )
    return mapped


def reason_code(row: pd.Series) -> str:
    text = str(row["All_reasons"])
    replacements = [
        ("frozen pre-QC/global exclusion; not morphology-shape adjudicated here", "PERM: upstream frozen exclusion"),
        ("possible under-reconstruction", "RECON: possible under-reconstruction"),
        ("direct ASC technical review", "ASC: direct technical review"),
        ("ASC parse failure", "ASC: parse failure"),
        ("S/P/E extreme morphology", "SPE: extreme morphology"),
        ("strict local multivariate outlier", "LOCAL: multivariate outlier"),
        ("selected pyramidal-like candidate", "PYR: pyramidal-like candidate"),
        ("user-confirmed review exclusion; flagged features: ", "OUTLIER: "),
        ("; maximum |robust z|=", "; max|z|="),
    ]
    for old, new in replacements:
        text = text.replace(old, new)
    return text


def prepare_items(group: pd.DataFrame):
    items = []
    for index, (_, row) in enumerate(group.iterrows()):
        plot_row = row.copy()
        plot_row["ASC_parse_success"] = bool(row["ASC_drawable_for_audit"])
        item = load_item(index, plot_row)
        if item.drawable:
            item.height_in += LABEL_BLOCK_IN - LABEL_HEIGHT_IN
        else:
            item.width_in = max(item.width_in, 1.20)
            item.height_in = 0.78
        items.append(item)
    return balanced_pack(items)


def add_scale_bar(fig: plt.Figure, page_width: float, page_height: float) -> None:
    x0 = PAGE_MARGIN_IN + 0.03
    y = PAGE_MARGIN_IN + 0.10
    fig.add_artist(
        Line2D(
            [x0 / page_width, (x0 + SCALE_BAR_IN) / page_width],
            [y / page_height, y / page_height],
            transform=fig.transFigure,
            color="black",
            linewidth=1.0,
            solid_capstyle="butt",
        )
    )
    fig.text(
        (x0 + SCALE_BAR_IN / 2) / page_width,
        (y + 0.025) / page_height,
        "50 µm",
        ha="center",
        va="bottom",
        fontsize=4,
        family="Arial",
    )


def draw_page(group: pd.DataFrame, page_title: str, page_index: int, pdf: PdfPages) -> tuple[Path, list[dict]]:
    items, content_width, content_height = prepare_items(group)
    page_width = max(5.1, content_width + 2 * PAGE_MARGIN_IN)
    page_height = content_height + HEADER_IN + FOOTER_IN + 2 * PAGE_MARGIN_IN
    fig = plt.figure(figsize=(page_width, page_height), facecolor="white")
    by_index = {i: row for i, (_, row) in enumerate(group.iterrows())}
    x_offset = (page_width - content_width) / 2
    audits = []

    for item in items:
        row = by_index[item.row_index]
        x = x_offset + item.x_in
        y = PAGE_MARGIN_IN + FOOTER_IN + item.y_in
        stage_color = STAGE_COLORS[row["Exclusion_stage"]]
        morph_height = item.height_in - LABEL_BLOCK_IN if item.drawable else 0.38
        ax = fig.add_axes([x / page_width, y / page_height, item.width_in / page_width, morph_height / page_height])
        if item.drawable:
            draw_neuron(ax, item, "#222222")
        else:
            ax.set_xlim(0, 1)
            ax.set_ylim(0, 1)
            ax.axis("off")
            ax.text(0.5, 0.5, "ASC unavailable / unparseable", ha="center", va="center", fontsize=3.6, color="#777777")

        asc_name = str(row.get("Selected_ASC_filename", "")) or "unavailable"
        year = row.get("Selected_modification_year", np.nan)
        year_text = "" if pd.isna(year) else f" [{int(year)}]"
        fig.text(
            (x + item.width_in / 2) / page_width,
            (y + morph_height + 0.015) / page_height,
            f"{shortened_id(row[ID])}\nM: {row['M_Neuron_id']} | ASC: {asc_name}{year_text}\n{reason_code(row)}",
            ha="center",
            va="bottom",
            fontsize=3.4,
            linespacing=1.02,
            color=stage_color,
            family="Arial",
            wrap=True,
        )
        audits.append(
            {
                ID: row[ID],
                "Gallery_page": page_index,
                "Drawable_ASC": item.drawable,
                "Panel_width_in": item.width_in,
                "Panel_height_in": item.height_in,
                "Physical_scale_um_per_in": 1.0 / INCH_PER_UM,
            }
        )

    fig.text(0.5, 1 - (PAGE_MARGIN_IN + 0.01) / page_height, page_title, ha="center", va="top", fontsize=7, family="Arial")
    add_scale_bar(fig, page_width, page_height)
    pdf.savefig(fig, facecolor="white")
    png = OUT / f"{page_index:02d}_{page_title.replace(' ', '_').replace('/', '-')}.png"
    fig.savefig(png, dpi=PAGE_PNG_DPI, facecolor="white")
    plt.close(fig)
    return png, audits


def draw_individual(row: pd.Series, order: int) -> dict:
    plot_row = row.copy()
    plot_row["ASC_parse_success"] = bool(row["ASC_drawable_for_audit"])
    item = load_item(0, plot_row)
    stage_color = STAGE_COLORS[row["Exclusion_stage"]]
    if item.drawable:
        width = max(2.15, item.width_in + 0.28)
        height = item.height_in - LABEL_HEIGHT_IN + 0.72
    else:
        width, height = 2.4, 1.45
    fig = plt.figure(figsize=(width, height), facecolor="white")
    if item.drawable:
        morph_h = height - 0.66
        ax = fig.add_axes([0.07 / width, 0.30 / height, (width - 0.14) / width, morph_h / height])
        draw_neuron(ax, item, "#222222")
    else:
        ax = fig.add_axes([0.08, 0.30, 0.84, 0.50])
        ax.axis("off")
        ax.text(0.5, 0.5, "ASC unavailable / unparseable", ha="center", va="center", fontsize=5, color="#777777")
    asc_name = str(row.get("Selected_ASC_filename", "")) or "unavailable"
    fig.text(
        0.5,
        0.985,
        f"{row[ID]} | {row['Exclusion_stage']}",
        ha="center",
        va="top",
        fontsize=6,
        color=stage_color,
        family="Arial",
    )
    fig.text(
        0.5,
        0.17,
        f"M: {row['M_Neuron_id']} | ASC: {asc_name}\n{reason_code(row)}",
        ha="center",
        va="bottom",
        fontsize=4,
        color=stage_color,
        family="Arial",
        wrap=True,
    )
    add_scale_bar(fig, width, height)
    stem = INDIVIDUAL / f"{order:02d}_{row[ID]}"
    fig.savefig(stem.with_suffix(".png"), dpi=INDIVIDUAL_PNG_DPI, facecolor="white")
    fig.savefig(stem.with_suffix(".pdf"), facecolor="white")
    plt.close(fig)
    return {ID: row[ID], "Individual_PNG": str(stem.with_suffix('.png')), "Individual_PDF": str(stem.with_suffix('.pdf'))}


def page_groups(data: pd.DataFrame):
    for stage in STAGE_COLORS:
        stage_data = data[data["Exclusion_stage"] == stage].reset_index(drop=True)
        for start in range(0, len(stage_data), MAX_TILES_PER_PAGE):
            part = stage_data.iloc[start : start + MAX_TILES_PER_PAGE].reset_index(drop=True)
            suffix = f"{start + 1}-{start + len(part)} of {len(stage_data)}" if len(stage_data) > MAX_TILES_PER_PAGE else f"n={len(stage_data)}"
            yield part, f"{stage} ({suffix})"


def main() -> None:
    configure_plotting()
    OUT.mkdir(parents=True, exist_ok=True)
    INDIVIDUAL.mkdir(parents=True, exist_ok=True)
    manifest = build_manifest()
    data = map_all_cells(manifest)
    if len(data) != 35 or data[ID].nunique() != 35:
        raise RuntimeError("Mapped audit cohort is not 35 unique cells")

    data["Reason_code"] = data.apply(reason_code, axis=1)
    data.to_csv(OUT / "01_all35_exclusion_reasons_and_ASC_mapping.csv", index=False)

    pdf_path = OUT / "02_all35_excluded_or_quarantined_ASC_reconstructions_same_scale.pdf"
    page_audits = []
    page_manifest = []
    with PdfPages(pdf_path) as pdf:
        for page_index, (group, title) in enumerate(page_groups(data), start=1):
            png, audits = draw_page(group, title, page_index, pdf)
            page_audits.extend(audits)
            page_manifest.append(
                {
                    "Gallery_page": page_index,
                    "Page_title": title,
                    "N_cells": len(group),
                    "N_drawable_ASC": int(group["ASC_drawable_for_audit"].sum()),
                    "PNG": str(png),
                }
            )

    individual_rows = [draw_individual(row, i) for i, (_, row) in enumerate(data.iterrows(), start=1)]
    pd.DataFrame(page_audits).to_csv(OUT / "03_gallery_panel_scale_audit.csv", index=False)
    pd.DataFrame(page_manifest).to_csv(OUT / "04_gallery_page_manifest.csv", index=False)
    pd.DataFrame(individual_rows).to_csv(OUT / "05_individual_figure_manifest.csv", index=False)

    stage_summary = (
        data.groupby(["Exclusion_stage", "Status"], sort=False)
        .agg(N_cells=(ID, "size"), N_ASC_drawable=("ASC_drawable_for_audit", "sum"))
        .reset_index()
    )
    stage_summary.to_csv(OUT / "06_stage_summary.csv", index=False)
    (OUT / "07_methods_and_interpretation.txt").write_text(
        "\n".join(
            [
                "Audit of all 35 cells removed or quarantined between the original complete Morph cohort (n=228) and current HC/GC cohort (n=193).",
                "Five upstream frozen exclusions and eight user-confirmed Morph numerical-review exclusions are excluded; twenty-two cells are in reversible morphology-only quarantine.",
                "The gallery is not evidence that every displayed cell is biologically invalid.",
                "ASC homonyms are resolved first by date key and hard acquisition-year evidence from M_ID/M_Neuron_id, then by agreement with scalar total length and radial distance.",
                "Parseable same-date/year ASC files with scalar metric mismatch are intentionally shown for visual review and flagged in the mapping table.",
                "Every neuron retains native XY orientation, soma-centered display, and 7.5 um biological padding; no rotation, mirroring or size normalization is applied.",
                "Absolute rendering scale is 250 um per inch, so every 50 um scale bar is exactly 0.2 inch in page PNG, individual PNG and PDF outputs.",
                "Overview PNG resolution is 600 dpi; individual PNG resolution is 900 dpi; all PDFs are vector outputs.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    print(stage_summary.to_string(index=False))
    print(data[["Mapping_status", "ASC_drawable_for_audit"]].value_counts(dropna=False).to_string())
    print(f"Output: {OUT}")


if __name__ == "__main__":
    main()
