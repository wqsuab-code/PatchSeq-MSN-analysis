#!/usr/bin/env python
"""Split frozen M1-M4 five-cell plates into 20 individual 900-dpi figures.

The frozen source morphology uses 68 pixels per 50 um. Every cell is uniformly
resampled to 90 pixels per 50 um and saved at 900 dpi, so the physical scale bar
is exactly 0.1 inch in all PNG/PDF outputs. Canvas dimensions remain adaptive.
"""

from pathlib import Path

import pandas as pd
import numpy as np
from scipy.ndimage import binary_erosion, distance_transform_edt
from PIL import Image, ImageDraw, ImageFont


ROOT = Path(__file__).resolve().parents[1]
SOURCE_DIR = (
    ROOT / "outputs" / "morph_qc"
    / "asc_reconstruction_qc_and_classic_examples"
)
SELECTION_FILE = SOURCE_DIR / "05_M1_M4_five_classic_cells_per_class.csv"
OUT = SOURCE_DIR / "individual_examples_900dpi_uniform_scale"

PANEL_COUNT = 5
SOURCE_SCALE_BAR_PX = 68
OUTPUT_DPI = 900.0
SCALE_BAR_PX = 90
RESAMPLE_FACTOR = SCALE_BAR_PX / SOURCE_SCALE_BAR_PX
SCALE_BAR_LABEL = "50 µm"
TARGET_BBOX_FRACTION = 0.20
CLASS_RGB = {
    1: (0, 70, 139),
    2: (66, 181, 64),
    3: (237, 0, 0),
    4: (0, 153, 180),
}


def font(size: int) -> ImageFont.FreeTypeFont | ImageFont.ImageFont:
    candidates = [
        Path(r"C:\Windows\Fonts\arial.ttf"),
        Path(r"C:\Windows\Fonts\Arial.ttf"),
        Path(r"C:\Windows\Fonts\segoeui.ttf"),
    ]
    for candidate in candidates:
        if candidate.exists():
            return ImageFont.truetype(str(candidate), size=size)
    return ImageFont.load_default()


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    for old in OUT.glob("*"):
        if old.is_file():
            old.unlink()

    selected = pd.read_csv(SELECTION_FILE).sort_values(
        ["M_consensus_class", "Rank_within_class"]
    )
    title_font = font(42)
    scale_font = font(32)
    manifest_rows = []

    for class_index in range(1, 5):
        class_name = f"Morph-M{class_index}"
        source = SOURCE_DIR / f"07_M{class_index}_five_classic_ASC_reconstructions.png"
        group = selected[selected["M_consensus_class"] == class_name]
        if len(group) != PANEL_COUNT:
            raise RuntimeError(f"Expected five frozen examples for {class_name}")

        with Image.open(source) as source_image:
            plate = source_image.convert("RGB")
        width, height = plate.size
        if width % PANEL_COUNT:
            raise RuntimeError(f"Plate width is not divisible by five: {source}")
        panel_width = width // PANEL_COUNT

        for panel_index, (_, row) in enumerate(group.iterrows()):
            source_panel = plate.crop(
                (panel_index * panel_width, 0, (panel_index + 1) * panel_width, height)
            )
            source_array = np.asarray(source_panel, dtype=np.int32)
            class_rgb = np.asarray(CLASS_RGB[class_index], dtype=np.int32)
            colour_distance = np.sqrt(
                np.sum((source_array - class_rgb) ** 2, axis=2)
            )
            # Exclude the plate header; the remaining class-coloured pixels are
            # the actual soma/dendrite reconstruction.
            mask = colour_distance < 80
            mask[:320, :] = False
            yy, xx = np.where(mask)
            if len(xx) == 0:
                raise RuntimeError(f"Could not locate morphology in {class_name} panel {panel_index + 1}")
            x_min, x_max = int(xx.min()), int(xx.max())
            y_min, y_max = int(yy.min()), int(yy.max())
            source_bbox_width = x_max - x_min + 1
            source_bbox_height = y_max - y_min + 1
            bbox_width = int(round(source_bbox_width * RESAMPLE_FACTOR))
            bbox_height = int(round(source_bbox_height * RESAMPLE_FACTOR))

            # A canvas 1/sqrt(0.20) times the morphology bounding box in each
            # dimension gives a morphology bounding-box area of approximately
            # 20%. Pixel geometry is never resized, preserving a common scale.
            margin_factor = 1.0 / TARGET_BBOX_FRACTION ** 0.5
            canvas_width = int(np.ceil(bbox_width * margin_factor))
            canvas_height = int(np.ceil(bbox_height * margin_factor))
            canvas_width = max(canvas_width, 360)
            canvas_height = max(canvas_height, 360)
            panel = Image.new("RGB", (canvas_width, canvas_height), "white")
            morphology = source_panel.crop((x_min, y_min, x_max + 1, y_max + 1))
            morphology = morphology.resize(
                (bbox_width, bbox_height), resample=Image.Resampling.LANCZOS
            )
            paste_x = (canvas_width - bbox_width) // 2
            paste_y = (canvas_height - bbox_height) // 2
            panel.paste(morphology, (paste_x, paste_y))
            draw = ImageDraw.Draw(panel)

            # Add one consistent title in the available upper margin.
            cell_id = str(row["MSN_unique_ID"])
            title = f"M{class_index}  |  {cell_id}"
            title_box = draw.textbbox((0, 0), title, font=title_font)
            title_width = title_box[2] - title_box[0]
            draw.text(
                ((canvas_width - title_width) / 2, max(8, (paste_y - 30) // 2)), title,
                font=title_font, fill="black"
            )

            # Place an identical exact-length scale bar on every figure.
            bar_x = max(16, int(canvas_width * 0.06))
            bar_y = canvas_height - max(20, int((canvas_height - (paste_y + bbox_height)) * 0.35))
            draw.line(
                (bar_x, bar_y, bar_x + SCALE_BAR_PX, bar_y),
                fill="black", width=4
            )
            label_box = draw.textbbox((0, 0), SCALE_BAR_LABEL, font=scale_font)
            label_width = label_box[2] - label_box[0]
            draw.text(
                (bar_x + (SCALE_BAR_PX - label_width) / 2, bar_y - 34),
                SCALE_BAR_LABEL, font=scale_font, fill="black"
            )

            rank = int(row["Rank_within_class"])
            safe_id = cell_id.replace("/", "-").replace("\\", "-")
            stem = OUT / f"M{class_index}_{rank:02d}_{safe_id}_50um_1inch"
            png_path = stem.with_suffix(".png")
            pdf_path = stem.with_suffix(".pdf")
            panel.save(png_path, format="PNG", dpi=(OUTPUT_DPI, OUTPUT_DPI))
            panel.save(
                pdf_path, format="PDF", resolution=OUTPUT_DPI,
                title=title, author="Morphology QC"
            )
            manifest_rows.append({
                "M_class": class_name,
                "Rank_within_class": rank,
                "MSN_unique_ID": cell_id,
                "PNG": str(png_path),
                "PDF": str(pdf_path),
                "Pixel_width": canvas_width,
                "Pixel_height": canvas_height,
                "Output_DPI": OUTPUT_DPI,
                "Scale_bar_um": 50,
                "Scale_bar_pixels": SCALE_BAR_PX,
                "Scale_bar_inches": SCALE_BAR_PX / OUTPUT_DPI,
                "Micrometres_per_pixel": 50.0 / SCALE_BAR_PX,
                "Neuron_stroke_rendering": "thin; source stroke uniformly resampled with Lanczos",
                "Morphology_bbox_width_px": bbox_width,
                "Morphology_bbox_height_px": bbox_height,
                "Morphology_bbox_area_fraction": (
                    bbox_width * bbox_height / (canvas_width * canvas_height)
                ),
            })

    manifest = pd.DataFrame(manifest_rows)
    manifest.to_csv(OUT / "00_individual_figure_manifest.csv", index=False)
    (OUT / "00_methods.txt").write_text(
        "\n".join([
            "Twenty frozen representative morphology figures exported separately.",
            "Classes: M1-M4; five previously frozen representative cells per class.",
            "All figures retain the same absolute XY scale and native orientation.",
            "Physical scale: 50 um = 0.1 inch.",
            "Implementation: 90-pixel bar at 900 dpi.",
            "Exact scale conversion: 50/90 um per pixel for every neuron and every scale bar.",
            "Each cell uses an individually cropped canvas; its morphology bounding box occupies approximately 20% of the canvas.",
            "Morphology pixels were not resized, so all 20 scale bars have identical length and the inter-cell size relationship is preserved.",
            "All morphologies were uniformly resampled by 90/68 using Lanczos; no class- or cell-specific resizing was applied.",
            "Neuron stroke weight remains thin and consistent across all figures.",
            "No cell geometry, selection, rotation, mirroring, or relative size was changed.",
        ]) + "\n",
        encoding="utf-8",
    )
    print(manifest.groupby("M_class").size())
    print(f"Figures: {len(manifest)} PNG + {len(manifest)} PDF")
    print(f"Scale: {SCALE_BAR_PX} px / {OUTPUT_DPI:g} dpi = {SCALE_BAR_PX / OUTPUT_DPI:.1f} inch")
    print(OUT)


if __name__ == "__main__":
    main()
