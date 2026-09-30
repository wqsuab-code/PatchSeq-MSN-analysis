#!/usr/bin/env python
"""Render the six PC1-PC3 3-SD Morph candidates at one absolute scale."""

from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import neurom
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
MAPPING = ROOT / "outputs/morph_qc/NPC3_res2p50_HCK4_G1-G4_ASC_matrix_cell_galleries/02_cell_to_ASC_year_metric_mapping_audit_193cells.csv"
PC = ROOT / "outputs/morph_qc/morph193_pca_3sd_outlier_screen/02_pc_scores_and_3SD_flags_all193.csv"
OUT = ROOT / "outputs/morph_qc/morph193_pca_3sd_six_ASC_review"
ASC_ROOT = Path("E:/ASC files")
IDS = [
    "D251161_Batch8",
    "A20256247_Batch7",
    "A20257158_Batch7",
    "A2025797_Batch7",
    "A20256245_Batch7",
    "O9192_Batch2",
]
FALLBACK_FILES = {
    "A20256247_Batch7": "6-24-7.asc",
    "A20257158_Batch7": "2025-7-15-8.asc",
    "A2025797_Batch7": "2025-7-9-7.asc",
    "A20256245_Batch7": "6-24-5.asc",
    "O9192_Batch2": "9-19-2.asc",
}
UM_PER_INCH = 250.0
SCALE_BAR_UM = 50.0
DPI = 900
LINE_WIDTH_PT = 0.42


def segments_and_soma(path: Path) -> tuple[list[np.ndarray], np.ndarray, np.ndarray]:
    neuron = neurom.load_morphology(str(path))
    soma = np.asarray(neuron.soma.points, dtype=float)
    if soma.ndim == 2 and soma.shape[0] and soma.shape[1] >= 2:
        center = soma[:, :2].mean(axis=0)
        soma_xy = soma[:, :2] - center
    else:
        center = np.asarray(neuron.soma.center[:2], dtype=float)
        soma_xy = np.empty((0, 2))
    segments = []
    for section in neuron.sections:
        points = np.asarray(section.points[:, :2], dtype=float)
        if len(points) >= 2 and np.isfinite(points).all():
            segments.append(points - center)
    return segments, soma_xy, center


def draw(ax, item: dict, xlim: tuple[float, float], ylim: tuple[float, float], title: bool = True) -> None:
    if item["path"] is None:
        ax.text(0.5, 0.57, "ASC reconstruction unavailable", transform=ax.transAxes,
                ha="center", va="center", fontsize=6, color="#ED0000")
        ax.text(0.5, 0.43, "Morphometric row retained; geometry cannot be adjudicated",
                transform=ax.transAxes, ha="center", va="center", fontsize=4.6, color="#555555")
    else:
        for points in item["segments"]:
            ax.plot(points[:, 0], points[:, 1], color="#202020", lw=LINE_WIDTH_PT,
                    solid_capstyle="round", solid_joinstyle="round")
        soma = item["soma"]
        if len(soma) >= 3:
            ax.fill(soma[:, 0], soma[:, 1], facecolor="#ED0000", edgecolor="none")
        # The exact same 50-um bar is drawn in every panel.
        bar_x = xlim[0] + 0.06 * (xlim[1] - xlim[0])
        bar_y = ylim[0] + 0.07 * (ylim[1] - ylim[0])
        ax.plot([bar_x, bar_x + SCALE_BAR_UM], [bar_y, bar_y], color="black", lw=0.65,
                solid_capstyle="butt")
        ax.text(bar_x + SCALE_BAR_UM / 2, bar_y + 0.025 * (ylim[1] - ylim[0]), "50 µm",
                ha="center", va="bottom", fontsize=4.5)
    if title:
        sd = item["sd"]
        filename = item["path"].name if item["path"] else "missing ASC"
        ax.text(
            0.5, 0.985, item["id"], transform=ax.transAxes,
            ha="center", va="top", fontsize=5.4,
        )
        ax.text(
            0.5, 0.945,
            f'{filename} | PC1 {sd[0]:+.2f}, PC2 {sd[1]:+.2f}, PC3 {sd[2]:+.2f} SD',
            transform=ax.transAxes, ha="center", va="top", fontsize=4.1,
        )
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.set_aspect("equal", adjustable="box")
    ax.axis("off")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    if not ASC_ROOT.is_dir():
        raise FileNotFoundError(f"ASC root is not mounted: {ASC_ROOT}")
    mapping = pd.read_csv(MAPPING)
    pc = pd.read_csv(PC).set_index("MSN_unique_ID")
    by_id = mapping.set_index("MSN_unique_ID")

    items = []
    for cell in IDS:
        candidate = None
        if cell in by_id.index:
            raw_path = str(by_id.loc[cell, "Selected_ASC_path"])
            if raw_path and raw_path.lower() != "nan":
                path = Path(raw_path.replace("D:\\", "E:\\"))
                if path.is_file():
                    candidate = path
        if candidate is None and cell in FALLBACK_FILES:
            path = ASC_ROOT / FALLBACK_FILES[cell]
            if path.is_file():
                candidate = path
        sd = pc.loc[cell, ["PC1_SD", "PC2_SD", "PC3_SD"]].astype(float).to_numpy()
        item = {"id": cell, "path": candidate, "sd": sd, "segments": [], "soma": np.empty((0, 2))}
        if candidate is not None:
            item["segments"], item["soma"], _ = segments_and_soma(candidate)
        items.append(item)

    arrays = [p for item in items for p in item["segments"]]
    arrays += [item["soma"] for item in items if len(item["soma"])]
    if not arrays:
        raise RuntimeError("No drawable morphology found")
    all_points = np.vstack(arrays)
    max_abs_x = max(abs(all_points[:, 0].min()), abs(all_points[:, 0].max()))
    max_abs_y = max(abs(all_points[:, 1].min()), abs(all_points[:, 1].max()))
    half_x = max(175.0, max_abs_x * 1.08)
    half_y = max(150.0, max_abs_y * 1.10)
    xlim, ylim = (-half_x, half_x), (-half_y, half_y)

    mpl.rcParams.update({
        "font.family": "Arial", "font.size": 5, "pdf.fonttype": 42,
        "ps.fonttype": 42, "savefig.facecolor": "white",
    })
    panel_w = (2 * half_x) / UM_PER_INCH
    panel_h = (2 * half_y) / UM_PER_INCH
    # A roomy comparison plate; identical coordinate limits ensure identical
    # morphology and scale-bar sizes across the six tiles.
    fig, axes = plt.subplots(2, 3, figsize=(9.0, 8.0), dpi=600)
    for ax, item in zip(axes.flat, items):
        draw(ax, item, xlim, ylim, title=True)
    fig.subplots_adjust(left=0.015, right=0.985, bottom=0.015, top=0.985, wspace=0.03, hspace=0.04)
    fig.savefig(OUT / "01_six_PC1-PC3_3SD_candidates_common_absolute_scale.png", dpi=900)
    fig.savefig(OUT / "01_six_PC1-PC3_3SD_candidates_common_absolute_scale.pdf")
    plt.close(fig)

    records = []
    for order, item in enumerate(items, 1):
        fig = plt.figure(figsize=(panel_w, panel_h), dpi=DPI)
        ax = fig.add_axes([0, 0, 1, 1])
        draw(ax, item, xlim, ylim, title=True)
        stem = OUT / f"individual_{order:02d}_{item['id']}_common_scale"
        fig.savefig(stem.with_suffix(".png"), dpi=DPI)
        fig.savefig(stem.with_suffix(".pdf"))
        plt.close(fig)
        records.append({
            "Order": order,
            "MSN_unique_ID": item["id"],
            "ASC_filename": item["path"].name if item["path"] else "",
            "ASC_path": str(item["path"]) if item["path"] else "",
            "Drawable": item["path"] is not None,
            "PC1_SD": item["sd"][0], "PC2_SD": item["sd"][1], "PC3_SD": item["sd"][2],
            "Common_x_field_um": 2 * half_x,
            "Common_y_field_um": 2 * half_y,
            "um_per_inch": UM_PER_INCH,
            "Scale_bar_um": SCALE_BAR_UM,
            "Scale_bar_inches": SCALE_BAR_UM / UM_PER_INCH,
            "Line_width_pt": LINE_WIDTH_PT,
        })
    pd.DataFrame(records).to_csv(OUT / "00_render_manifest.csv", index=False)
    print(pd.DataFrame(records).to_string(index=False))


if __name__ == "__main__":
    main()
