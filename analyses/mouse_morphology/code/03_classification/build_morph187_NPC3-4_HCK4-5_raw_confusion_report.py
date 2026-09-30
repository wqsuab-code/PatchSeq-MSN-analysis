#!/usr/bin/env python
"""Build the 204-page-ready raw HC-by-GC matrix report for Morph n=187."""

from __future__ import annotations

import importlib.util
from pathlib import Path

import pandas as pd
from reportlab.lib import colors
from reportlab.pdfgen import canvas


ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "scripts/build_morph193_NPC3-4_HCK4-5_raw_confusion_report.py"
SRC = ROOT / "outputs/morph_qc/morph187_NPC3-4_HCK4-5_GC_allres_raw"
OUT = ROOT / "output/pdf/morph187_NPC3-4_HCK4-5_all_raw_GC_confusion_matrices.pdf"


def load_base():
    spec = importlib.util.spec_from_file_location("raw_report_base", BASE)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def main() -> None:
    m = load_base()
    counts = pd.read_csv(SRC / "06_all_raw_cross_confusion_counts_long.csv")
    metrics = pd.read_csv(SRC / "07_all_NPC_HCK_resolution_raw_metrics.csv")
    if len(metrics) != 204:
        raise RuntimeError(f"Expected 204 comparisons, found {len(metrics)}")
    OUT.parent.mkdir(parents=True, exist_ok=True)

    def summary_page(c):
        m.header(
            c,
            "Morph n=187: raw HC x GC scan",
            "NPC=3-4; HC K=4-5; GC resolution=0.50-3.00. GC classes are never merged.",
        )
        y = m.PAGE_H - 75
        m.draw_text(c, m.MARGIN, y, "Best raw ARI at each PCA/HC setting", 9, True)
        y -= 22
        best = metrics.sort_values(
            ["NPC", "HC_K", "Raw_ARI", "Raw_NMI"],
            ascending=[True, True, False, False],
        ).groupby(["NPC", "HC_K"]).head(1)
        heads = ["NPC", "HC K", "Resolution", "Raw GC K", "ARI", "NMI", "GC-to-HC purity", "HC sizes", "GC sizes"]
        widths = [42, 42, 65, 58, 55, 55, 85, 95, 105]
        xx = m.MARGIN
        for h, w in zip(heads, widths):
            c.setFillColor(m.LIGHT); c.rect(xx, y - 4, w, 18, fill=1, stroke=0)
            m.draw_text(c, xx + w / 2, y + 2, h, 6.2, True, "center"); xx += w
        y -= 20
        for row in best.itertuples():
            values = [row.NPC, row.HC_K, f"{row.Resolution:.2f}", row.GC_raw_K,
                      f"{row.Raw_ARI:.3f}", f"{row.Raw_NMI:.3f}",
                      f"{100 * row.Raw_GC_to_HC_purity:.1f}%", row.HC_sizes, row.GC_raw_sizes]
            xx = m.MARGIN
            for value, w in zip(values, widths):
                m.draw_text(c, xx + w / 2, y + 2, value, 6.5, align="center")
                c.setStrokeColor(colors.HexColor("#CCCCCC")); c.line(xx, y - 4, xx + w, y - 4); xx += w
            y -= 18

        y -= 20
        m.draw_text(c, m.MARGIN, y, "Raw GC class count across resolution", 9, True); y -= 17
        for npc in [3, 4]:
            parts = []
            for start, end, k in m.raw_k_runs(metrics, npc):
                span = f"{start:.2f}" if abs(start - end) < 1e-9 else f"{start:.2f}-{end:.2f}"
                parts.append(f"{span}: K={k}")
            m.draw_text(c, m.MARGIN + 8, y, f"NPC={npc}: " + ";  ".join(parts), 6.3)
            y -= 20

        y -= 8
        m.draw_text(c, m.MARGIN, y, "Interpretation", 9, True)
        notes = [
            "Rows are Ward.D2 HC classes; columns are all unmerged graph-based clusters.",
            "ARI and NMI remain interpretable when HC K and raw GC K differ.",
            "Purity is directional and must not be described as classification accuracy.",
            "Each matrix contains actual counts for all 187 retained cells.",
            "The six newly quarantined incomplete reconstructions are absent from PCA, HC and GC fitting.",
        ]
        for note in notes:
            y -= 15; m.draw_text(c, m.MARGIN + 8, y, "- " + note, 7)
        m.footer(c, 1); c.showPage()

    c = canvas.Canvas(str(OUT), pagesize=(m.PAGE_W, m.PAGE_H), pageCompression=1)
    c.setTitle("Morph n=187 NPC3-4 HC K4-5 versus all unmerged GC resolutions")
    summary_page(c)
    pages = m.matrix_pages(c, counts, metrics)
    c.save()
    print(f"Created {OUT}")
    print(f"Pages={pages}; matrices={len(metrics)}")


if __name__ == "__main__":
    main()
