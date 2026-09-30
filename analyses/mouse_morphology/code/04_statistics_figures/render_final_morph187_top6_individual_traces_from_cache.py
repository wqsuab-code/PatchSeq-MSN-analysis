#!/usr/bin/env python
"""Render six centroid-proximal, consensus-supported traces per final M class.

The source gallery was rendered directly from ASC at 600 dpi and a frozen
physical scale of 250 um/in.  This script extracts morphology-only pixels
without geometric resizing, recolours them to the final Lancet M palette, and
places every trace on one common square canvas.  Thus the canvas dimensions,
physical morphology scale, and 50-um scale-bar length are identical for all
24 independent outputs.
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image, ImageDraw, ImageFont


ROOT = Path(__file__).resolve().parents[1]
FINAL_DIR = ROOT / "outputs" / "morph_qc" / "final_morph187_NPC3_HCK4_GCres0p50_figures"
ASSIGNMENTS = FINAL_DIR / "00_final_morph187_cell_assignments.csv"
GALLERY_DIR = ROOT / "outputs" / "morph_qc" / "NPC3_res2p50_HCK4_G1-G4_ASC_matrix_cell_galleries"
AUDIT = GALLERY_DIR / "03_all_gallery_tiles_and_source_files.csv"
PAGE_SUMMARY = GALLERY_DIR / "04_matrix_cell_page_summary.csv"
OUT = ROOT / "outputs" / "morph_qc" / "final_morph187_NPC3_HCK4_GCres0p50_representative_traces"

ID = "MSN_unique_ID"
M_ORDER = ["M1", "M2", "M3", "M4"]
M_COLORS = {
    "M1": "#00468B",
    "M2": "#42B540",
    "M3": "#ED0000",
    "M4": "#0099B4",
}
OLD_G_COLORS = {
    "G1": "#00468B",
    "G2": "#42B540",
    "G3": "#ED0000",
    "G4": "#0099B4",
}

SOURCE_SCALE_UM_PER_IN = 250.0
SOURCE_LABEL_HEIGHT_IN = 0.285
PAGE_MARGIN_IN = 0.10
FOOTER_IN = 0.28
OUTPUT_DPI = 900
SCALE_BAR_UM = 50.0
SCALE_BAR_IN = SCALE_BAR_UM / SOURCE_SCALE_UM_PER_IN


def rgb(hex_colour: str) -> np.ndarray:
    value = hex_colour.lstrip("#")
    return np.array([int(value[i:i + 2], 16) for i in (0, 2, 4)], dtype=float)


def arial(size: int, bold: bool = False) -> ImageFont.FreeTypeFont:
    name = "arialbd.ttf" if bold else "arial.ttf"
    path = Path(r"C:\Windows\Fonts") / name
    return ImageFont.truetype(str(path), size=size)


def page_png(page_number: int) -> Path:
    matches = sorted(GALLERY_DIR.glob(f"{page_number:02d}_*_ASC_gallery.png"))
    if len(matches) != 1:
        raise RuntimeError(f"Expected one cached gallery for page {page_number}, found {len(matches)}")
    return matches[0]


def extract_trace(row: pd.Series, pages: pd.DataFrame) -> Image.Image:
    summary = pages.loc[pages["PDF_page"] == int(row["PDF_page"])].iloc[0]
    page_w_in = float(summary["Page_width_in"])
    page_h_in = float(summary["Page_height_in"])
    same_page = audit.loc[audit["PDF_page"] == int(row["PDF_page"])]
    content_w_in = float((same_page["Panel_x_in"] + same_page["Panel_width_in"]).max())
    content_x_offset = (page_w_in - content_w_in) / 2.0

    image = Image.open(page_png(int(row["PDF_page"]))).convert("RGB")
    xdpi = image.width / page_w_in
    ydpi = image.height / page_h_in
    x0_in = content_x_offset + float(row["Panel_x_in"])
    y0_in = PAGE_MARGIN_IN + FOOTER_IN + float(row["Panel_y_in"])
    width_in = float(row["Panel_width_in"])
    height_in = float(row["Panel_height_in"]) - SOURCE_LABEL_HEIGHT_IN

    left = max(0, int(round(x0_in * xdpi)))
    right = min(image.width, int(round((x0_in + width_in) * xdpi)))
    top = max(0, int(round(image.height - (y0_in + height_in) * ydpi)))
    bottom = min(image.height, int(round(image.height - y0_in * ydpi)))
    crop = np.asarray(image.crop((left, top, right, bottom))).astype(float)

    base = rgb(OLD_G_COLORS[str(row["GC_cluster"])])
    denominators = 255.0 - base
    channels = denominators > 15
    alpha_stack = (255.0 - crop[:, :, channels]) / denominators[channels]
    alpha = np.clip(np.median(alpha_stack, axis=2), 0.0, 1.0)
    # Exclude neutral text/background residue; morphology is chromatic.
    chroma = crop.max(axis=2) - crop.min(axis=2)
    alpha[chroma < 5.0] = 0.0
    alpha[alpha < 0.025] = 0.0

    yy, xx = np.where(alpha > 0)
    if not len(xx):
        raise RuntimeError(f"No trace pixels detected for {row[ID]}")
    pad = 3
    x0 = max(0, int(xx.min()) - pad)
    x1 = min(alpha.shape[1], int(xx.max()) + pad + 1)
    y0 = max(0, int(yy.min()) - pad)
    y1 = min(alpha.shape[0], int(yy.max()) + pad + 1)
    alpha = alpha[y0:y1, x0:x1]

    target = rgb(M_COLORS[str(row["M_class"])])
    recoloured = np.full((alpha.shape[0], alpha.shape[1], 4), 255, dtype=np.uint8)
    recoloured[:, :, :3] = target.astype(np.uint8)
    recoloured[:, :, 3] = np.round(alpha * 255).astype(np.uint8)
    return Image.fromarray(recoloured, mode="RGBA")


def bool_value(value) -> bool:
    return str(value).strip().lower() in {"true", "1", "yes"}


assign = pd.read_csv(ASSIGNMENTS)
audit = pd.read_csv(AUDIT)
pages = pd.read_csv(PAGE_SUMMARY)
for col in ["PDF_page", "Panel_x_in", "Panel_y_in", "Panel_width_in", "Panel_height_in"]:
    audit[col] = pd.to_numeric(audit[col], errors="coerce")
pages["PDF_page"] = pd.to_numeric(pages["PDF_page"], errors="coerce")

centroids = assign.groupby("M_class")[["PC1", "PC2", "PC3"]].mean()
for pc in ["PC1", "PC2", "PC3"]:
    assign[pc] = pd.to_numeric(assign[pc], errors="coerce")
assign["distance_to_M_centroid_PC1-PC3"] = [
    float(np.linalg.norm(row[["PC1", "PC2", "PC3"]].to_numpy(float) - centroids.loc[row["M_class"]].to_numpy(float)))
    for _, row in assign.iterrows()
]

availability = audit[[
    ID, "PDF_page", "GC_cluster", "Selected_ASC_filename", "Selected_ASC_path",
    "ASC_renderable", "Panel_x_in", "Panel_y_in", "Panel_width_in", "Panel_height_in",
]].drop_duplicates(ID)
data = assign.merge(availability, on=ID, how="left", suffixes=("", "_gallery"))
data["ASC_renderable_bool"] = data["ASC_renderable"].map(bool_value)
data["consensus_bool"] = data["HC_GC_consensus"].map(bool_value)

selected_rows = []
for m in M_ORDER:
    group = data.loc[
        (data["M_class"] == m) & data["consensus_bool"] & data["ASC_renderable_bool"]
    ].sort_values(["distance_to_M_centroid_PC1-PC3", ID], kind="mergesort")
    if len(group) < 6:
        raise RuntimeError(f"Only {len(group)} consensus/renderable cells for {m}")
    chosen = group.head(6).copy()
    chosen["representative_rank"] = np.arange(1, 7)
    selected_rows.append(chosen)
selected = pd.concat(selected_rows, ignore_index=True)

# First extract all traces to determine one common canvas that fits every cell.
traces: dict[str, Image.Image] = {}
for _, row in selected.iterrows():
    traces[row[ID]] = extract_trace(row, pages)

max_w_in = max(im.width / 600.0 for im in traces.values())
max_h_in = max(im.height / 600.0 for im in traces.values())
required = max(max_w_in + 0.38, max_h_in + 0.68, 1.80)
canvas_in = math.ceil(required * 10.0) / 10.0
canvas_px = int(round(canvas_in * OUTPUT_DPI))

OUT.mkdir(parents=True, exist_ok=True)
manifest_rows = []
font_title = arial(round(6 / 72 * OUTPUT_DPI), bold=True)
font_sub = arial(round(4 / 72 * OUTPUT_DPI), bold=False)
font_scale = arial(round(4 / 72 * OUTPUT_DPI), bold=False)
top_reserved = int(round(0.28 * OUTPUT_DPI))
bottom_reserved = int(round(0.30 * OUTPUT_DPI))

for _, row in selected.iterrows():
    trace = traces[row[ID]]
    # Convert 600-dpi cached pixels to 900-dpi pixels while preserving inches.
    new_size = (int(round(trace.width * OUTPUT_DPI / 600.0)), int(round(trace.height * OUTPUT_DPI / 600.0)))
    trace = trace.resize(new_size, Image.Resampling.LANCZOS)
    canvas = Image.new("RGB", (canvas_px, canvas_px), "white")
    draw = ImageDraw.Draw(canvas)
    x = (canvas_px - trace.width) // 2
    available_h = canvas_px - top_reserved - bottom_reserved
    y = top_reserved + (available_h - trace.height) // 2
    canvas.paste(trace, (x, y), trace)

    draw.text((int(0.10 * OUTPUT_DPI), int(0.07 * OUTPUT_DPI)),
              f"{row['M_class']}  {row[ID]}", fill=M_COLORS[row["M_class"]], font=font_title)
    draw.text((int(0.10 * OUTPUT_DPI), int(0.16 * OUTPUT_DPI)),
              f"ASC: {row['Selected_ASC_filename']}", fill="#333333", font=font_sub)

    bar_px = int(round(SCALE_BAR_IN * OUTPUT_DPI))
    bar_x1 = canvas_px - int(round(0.12 * OUTPUT_DPI))
    bar_x0 = bar_x1 - bar_px
    bar_y = canvas_px - int(round(0.14 * OUTPUT_DPI))
    line_width = max(1, int(round(1.0 / 72.0 * OUTPUT_DPI)))
    draw.line((bar_x0, bar_y, bar_x1, bar_y), fill="black", width=line_width)
    label = "50 um"
    bbox = draw.textbbox((0, 0), label, font=font_scale)
    draw.text(((bar_x0 + bar_x1 - (bbox[2] - bbox[0])) / 2, bar_y - (bbox[3] - bbox[1]) - 5),
              label, fill="black", font=font_scale)

    stem = f"{row['M_class']}_{int(row['representative_rank']):02d}_{row[ID]}_common_canvas_scale"
    png = OUT / f"{stem}.png"
    pdf = OUT / f"{stem}.pdf"
    canvas.save(png, dpi=(OUTPUT_DPI, OUTPUT_DPI), compress_level=6)
    canvas.save(pdf, "PDF", resolution=OUTPUT_DPI)
    manifest_rows.append({
        "M_class": row["M_class"],
        "representative_rank": int(row["representative_rank"]),
        ID: row[ID],
        "HC_GC_consensus": bool(row["consensus_bool"]),
        "distance_to_M_centroid_PC1-PC3": row["distance_to_M_centroid_PC1-PC3"],
        "ASC_filename": row["Selected_ASC_filename"],
        "cached_ASC_source_path": row["Selected_ASC_path"],
        "PNG_output": str(png),
        "PDF_output": str(pdf),
        "canvas_width_in": canvas_in,
        "canvas_height_in": canvas_in,
        "output_DPI": OUTPUT_DPI,
        "physical_scale_um_per_in": SOURCE_SCALE_UM_PER_IN,
        "scale_bar_um": SCALE_BAR_UM,
        "scale_bar_length_in": SCALE_BAR_IN,
    })

manifest = pd.DataFrame(manifest_rows)
manifest.to_csv(OUT / "00_top6_representative_trace_manifest.csv", index=False)

# Record closer cells that were bypassed only because a reliable cached ASC was absent
# or because HC and GC did not agree.
bypass_rows = []
for m in M_ORDER:
    cutoff = float(manifest.loc[manifest["M_class"] == m, "distance_to_M_centroid_PC1-PC3"].max())
    closer = data.loc[(data["M_class"] == m) & (data["distance_to_M_centroid_PC1-PC3"] <= cutoff)].copy()
    closer = closer.loc[~(closer["consensus_bool"] & closer["ASC_renderable_bool"])]
    for _, row in closer.iterrows():
        bypass_rows.append({
            "M_class": m,
            ID: row[ID],
            "distance_to_M_centroid_PC1-PC3": row["distance_to_M_centroid_PC1-PC3"],
            "HC_GC_consensus": row["consensus_bool"],
            "ASC_renderable_from_cache": row["ASC_renderable_bool"],
            "reason_bypassed": "HC-GC discordant" if not row["consensus_bool"] else "reliable cached ASC unavailable",
        })
pd.DataFrame(bypass_rows).to_csv(OUT / "00_closer_cells_bypassed_audit.csv", index=False)

(OUT / "00_rendering_methods.txt").write_text(
    "Selection: within each final HC-derived M1-M4 class, Euclidean distance to the class centroid "
    "in PC1-PC3; restricted to HC-GC-consensus cells with a verified/renderable cached ASC trace.\n"
    "Rendering: morphology coordinates were extracted from prior direct-ASC 600-dpi pages without "
    "geometric normalization or per-cell fitting. Every trace retains 250 um/in. All output canvases "
    f"are {canvas_in:.1f} x {canvas_in:.1f} in at {OUTPUT_DPI} dpi. Every scale bar is 50 um = "
    f"{SCALE_BAR_IN:.1f} in. Colour is the final Lancet M palette.\n",
    encoding="utf-8",
)

print(manifest[["M_class", "representative_rank", ID, "ASC_filename", "distance_to_M_centroid_PC1-PC3"]].to_string(index=False))
print(f"Canvas: {canvas_in:.1f} x {canvas_in:.1f} in; {OUTPUT_DPI} dpi; outputs={len(manifest)}")
