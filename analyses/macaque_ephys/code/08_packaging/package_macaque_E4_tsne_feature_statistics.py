#!/usr/bin/env python3
"""Build the complete statistics and reproducibility bundle for the Macaque E4 figure."""

from __future__ import annotations

import hashlib
import json
import platform
import shutil
import subprocess
import sys
from pathlib import Path

import numpy
import pandas
import scipy


ROOT = Path(__file__).resolve().parents[1]
NAME = "Macaque_E4_TSNE_feature_comparison_stats_bundle_20260921"
OUT = ROOT / "outputs" / NAME
ZIP = ROOT / "outputs" / f"{NAME}.zip"
SITE = ROOT / "outputs" / "Macaque_E4_feature_explorer_site" / "dist"
STATIC = ROOT / "outputs" / "R3_E_feature_distribution_tsne_all19"
REFERENCE_LAYOUT = Path(
    r"C:\Users\53461\OneDrive\Desktop\Patch-seq\4Manuscript\Macaque_E4_ML_layout (Ver20260921).json"
)


def copy(source: Path, destination: Path) -> None:
    if not source.exists():
        raise FileNotFoundError(source)
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, destination)


def write_readme() -> None:
    text = """Macaque E1-E4 t-SNE + feature-comparison statistics bundle
================================================================

Scope
-----
This archive contains the data snapshot, statistical code, complete numerical
results, interactive figure code, layout files and static exports used for the
19-feature Macaque E1-E4 comparison figure. The analyzed population is fixed at
368 HC-GC consensus MSN cells from Ca + Pu + NAc: E1=59, E2=133, E3=57,
E4=119. The main-analysis class colors are retained: E1 #F8766D, E2 #7CAE00,
E3 #00BFC4 and E4 #C77CFF.

Primary inference
-----------------
For every feature, raw displayed values are compared among E1-E4 using a
two-sided Kruskal-Wallis rank test (df=3). The 19 omnibus P values are adjusted
with the Benjamini-Hochberg procedure. All six pairwise contrasts per feature
use a two-sided Dunn rank-sum test with average ranks and tie correction.
Pairwise P values are Holm-adjusted within each feature. The archive also
reports BH adjustment across all 19 x 6 = 114 pairwise tests. Dunn effect size
is |Z|/sqrt(368). Omnibus epsilon-squared is max[0,(H-k+1)/(N-k)], with k=4.
No observations are removed, winsorized or imputed in these comparisons.

Plot definitions
----------------
Box plots show the median, interquartile range and whiskers. The interactive
figure defaults to 1.5 x IQR whiskers; the original static script used the
2.5th and 97.5th percentiles. Individual cells are displayed with jitter. The
t-SNE panels use the frozen coordinates and feature Z-scores clipped to
[-2.5, 2.5]. E1-E4 outlines and colors are display annotations and do not alter
the statistical analysis. The interactive figure displays Holm-adjusted
within-feature P values by default and can switch to global BH-adjusted P values.

Folder map
----------
01_code: exact recomputation script, original static-figure/statistics script,
         package builder, and browser rendering code.
02_frozen_input: exact 368-cell raw values, frozen t-SNE/class assignments,
                 19-feature Z-score table and feature definitions.
03_statistics: complete descriptive, omnibus and pairwise results plus a
               numerical audit against the values embedded in the webpage.
04_figure_and_layout: static PNG/PDF/SVG, interactive HTML/CSS/JS and both the
                      current comparison layout and supplied ML layout reference.
05_provenance: software versions and SHA-256 manifest.

Reproduction
------------
From the repository root, run:

python scripts/recompute_macaque_E4_feature_statistics.py \
  --input outputs/Macaque_E4_feature_explorer_site/dist/data/raw.json \
  --reference-pairwise outputs/Macaque_E4_feature_explorer_site/dist/data/pairwise.json \
  --reference-omnibus outputs/Macaque_E4_feature_explorer_site/dist/data/omnibus.json \
  --output outputs/Macaque_E4_recomputed_statistics

The recalculation audit should report only floating-point-level differences
from the webpage results. Run the package builder to reconstruct this archive:

python scripts/package_macaque_E4_tsne_feature_statistics.py

Important distinction
---------------------
The raw-value group comparisons and the continuous feature t-SNE overlays use
the same 368 frozen consensus cells, but they serve different purposes. The
rank tests operate on raw displayed values; t-SNE point colors use standardized
feature Z-scores. Neither the t-SNE embedding nor ellipse/outlining parameters
enter the Kruskal-Wallis or Dunn tests.
"""
    (OUT / "README.txt").write_text(text, encoding="utf-8")


def manifest() -> None:
    rows = []
    for path in sorted(p for p in OUT.rglob("*") if p.is_file()):
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        rows.append(f"{digest}  {path.relative_to(OUT).as_posix()}")
    (OUT / "05_provenance" / "SHA256SUMS.txt").write_text("\n".join(rows) + "\n", encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=False)
    for folder in ["01_code", "02_frozen_input", "03_statistics", "04_figure_and_layout", "05_provenance"]:
        (OUT / folder).mkdir()

    copy(ROOT / "scripts" / "recompute_macaque_E4_feature_statistics.py", OUT / "01_code" / "recompute_macaque_E4_feature_statistics.py")
    copy(ROOT / "scripts" / "plot_ca_pu_nac_consensus368_all19_feature_distribution_tsne.py", OUT / "01_code" / "plot_ca_pu_nac_consensus368_all19_feature_distribution_tsne.py")
    copy(ROOT / "scripts" / "package_macaque_E4_tsne_feature_statistics.py", OUT / "01_code" / "package_macaque_E4_tsne_feature_statistics.py")
    for name in ["comparison-editor.html", "comparison-editor.js", "comparison-ml-theme.css", "styles.css"]:
        copy(SITE / name, OUT / "01_code" / "interactive_figure" / name)

    for name in ["raw.json", "tsne_cells.csv", "final19_zscore_long.csv", "pairwise.json", "omnibus.json", "comparison_default_layout.json"]:
        copy(SITE / "data" / name, OUT / "02_frozen_input" / name)
    copy(STATIC / "all19_feature_order_and_units.csv", OUT / "02_frozen_input" / "all19_feature_order_and_units.csv")
    copy(STATIC / "main_tsne_display_parameters.txt", OUT / "02_frozen_input" / "main_tsne_display_parameters.txt")

    subprocess.run(
        [
            sys.executable,
            str(ROOT / "scripts" / "recompute_macaque_E4_feature_statistics.py"),
            "--input", str(SITE / "data" / "raw.json"),
            "--reference-pairwise", str(SITE / "data" / "pairwise.json"),
            "--reference-omnibus", str(SITE / "data" / "omnibus.json"),
            "--output", str(OUT / "03_statistics"),
        ],
        check=True,
    )
    copy(SITE / "data" / "pairwise.json", OUT / "03_statistics" / "07_webpage_pairwise_results.json")
    copy(SITE / "data" / "omnibus.json", OUT / "03_statistics" / "08_webpage_omnibus_results.json")

    stem = "Macaque_CaPuNAC_E4_all19_featureTSNE_mainParams_4x5"
    for suffix in [".png", ".pdf", ".svg"]:
        copy(STATIC / f"{stem}{suffix}", OUT / "04_figure_and_layout" / f"{stem}{suffix}")
    for name in ["comparison-editor.html", "comparison-editor.js", "comparison-ml-theme.css", "styles.css"]:
        copy(SITE / name, OUT / "04_figure_and_layout" / "interactive_webpage" / name)
    for name in ["raw.json", "tsne_cells.csv", "final19_zscore_long.csv", "pairwise.json", "omnibus.json", "comparison_default_layout.json"]:
        copy(SITE / "data" / name, OUT / "04_figure_and_layout" / "interactive_webpage" / "data" / name)
    copy(SITE / "data" / "comparison_default_layout.json", OUT / "04_figure_and_layout" / "comparison_default_layout.json")
    if REFERENCE_LAYOUT.exists():
        copy(REFERENCE_LAYOUT, OUT / "04_figure_and_layout" / REFERENCE_LAYOUT.name)

    write_readme()
    versions = {
        "python": platform.python_version(),
        "platform": platform.platform(),
        "numpy": numpy.__version__,
        "pandas": pandas.__version__,
        "scipy": scipy.__version__,
    }
    (OUT / "05_provenance" / "software_versions.json").write_text(json.dumps(versions, indent=2), encoding="utf-8")
    manifest()
    shutil.make_archive(str(ZIP.with_suffix("")), "zip", OUT.parent, OUT.name)
    print(OUT)
    print(ZIP)


if __name__ == "__main__":
    main()
