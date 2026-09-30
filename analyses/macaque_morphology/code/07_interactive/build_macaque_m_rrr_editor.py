from __future__ import annotations

import base64
import json
import re
import shutil
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE_UI = ROOT / "interactive" / "mouse-morph-feature-explorer" / "dist"
DEST = ROOT / "interactive" / "macaque-morph-feature-explorer" / "dist"
REFERENCE = Path(
    r"C:\Users\53461\Downloads\Macaque_T-E_RRR_gc_merged_layout (Ver20260922).json"
)

M_COLORS = {
    "D1": "#D95F02",
    "D2": "#008F7A",
    "M1": "#1F77B4",
    "M2": "#D9A400",
    "M3": "#8C564B",
    "M4": "#E377C2",
}


def encode(path: Path) -> str:
    return base64.b64encode(path.read_bytes()).decode("ascii")


def build_layout() -> dict:
    ref = json.loads(REFERENCE.read_text(encoding="utf-8-sig"))
    axis_display = {
        "T12": ref["axisDisplay"]["T12"],
        "T13": ref["axisDisplay"]["T13"],
        "M12": ref["axisDisplay"]["E12"],
        "M13": ref["axisDisplay"]["E13"],
    }
    pairs = {
        "T12": ref["pairs"]["T12"],
        "T13": ref["pairs"]["T13"],
        "M12": ref["pairs"]["E12"],
        "M13": ref["pairs"]["E13"],
    }
    # T-space labels can reuse their E-figure placements when the same gene is
    # selected. E-space label positions are intentionally not transferred to
    # unrelated morphology features.
    positions = {key: value for key, value in ref.get("positions", {}).items() if key.startswith("T")}
    return {
        "geneCount": ref["geneCount"],
        "featureCount": ref["featureCount"],
        "fontSize": ref["fontSize"],
        "labelAlign": ref["labelAlign"],
        "boxOpacity": ref["boxOpacity"],
        "cellSize": ref["cellSize"],
        "cellOpacity": ref["cellOpacity"],
        "panelLayout": ref["panelLayout"],
        "plotRadius": ref["plotRadius"],
        "panelWidth": ref["panelWidth"],
        "panelHeight": ref["panelHeight"],
        "columnGap": ref["columnGap"],
        "rowGap": ref["rowGap"],
        "marginLeft": ref["marginLeft"],
        "marginRight": ref["marginRight"],
        "marginTop": ref["marginTop"],
        "marginBottom": ref["marginBottom"],
        "radialPercentile": ref["radialPercentile"],
        "ellipseLevel": ref["ellipseLevel"],
        "loadingScale": ref["loadingScale"],
        "fontFamily": ref["fontFamily"],
        "arrowWidth": ref["arrowWidth"],
        "connectorWidth": ref["connectorWidth"],
        "ellipseWidth": ref["ellipseWidth"],
        "ellipseOpacity": ref["ellipseOpacity"],
        "circleWidth": ref["circleWidth"],
        "axisWidth": ref["axisWidth"],
        "titleSize": ref["titleSize"],
        "axisLabelSize": ref["axisLabelSize"],
        "axisColor": ref["axisColor"],
        "axisTextColor": ref["axisTextColor"],
        "showTicks": ref["showTicks"],
        "tickCount": ref["tickCount"],
        "tickDecimals": ref["tickDecimals"],
        "tickSize": ref["tickSize"],
        "titleT": ref["titleT"],
        "titleM": "Morphological space",
        "canvasColor": ref["canvasColor"],
        "exportWidth": ref["exportWidth"],
        "exportDpi": ref["exportDpi"],
        "showBoxes": ref["showBoxes"],
        "showConnectors": ref["showConnectors"],
        "showEllipses": ref["showEllipses"],
        "showLegend": ref["showLegend"],
        "showCircle": ref["showCircle"],
        "showPoints": ref["showPoints"],
        "axisDisplay": axis_display,
        "pairs": pairs,
        "positions": positions,
        "aliases": {},
        "colors": M_COLORS,
        "layoutReference": str(REFERENCE),
        "dataScope": "Frozen 117-cell Macaque M1-M4 HC-GC consensus cohort",
    }


def build_html() -> None:
    text = (SOURCE_UI / "rrr-editor.html").read_text(encoding="utf-8")
    text = text.replace("Mouse T–M RRR figure editor", "Macaque T–M RRR figure editor")
    text = text.replace(
        "Interactive editor for the Mouse transcriptomic-to-morphological reduced-rank regression figure.",
        "Interactive editor for the Macaque transcriptomic-to-morphological reduced-rank regression figure.",
    )
    text = re.sub(
        r"Rank 3 · n = 168 .*? Mouse_M data and colours",
        "Rank 3 · n = 117 HC–GC consensus MSN cells · Macaque-E reference geometry with frozen Macaque-M data and colours",
        text,
    )
    text = text.replace("1 × 4 · Mouse_E reference", "1 × 4 · Macaque-E reference")
    text = text.replace("Live Mouse_M preview", "Live Macaque-M preview")
    text = text.replace(
        "Interactive Mouse transcriptomic-to-morphological reduced-rank regression figure",
        "Interactive Macaque transcriptomic-to-morphological reduced-rank regression figure",
    )
    text = text.replace(
        '<a class="button" href="d1d2-distribution-editor.html">D1/D2 distribution editor</a><a class="button" href="percentage-editor.html">D1/D2 percentage editor</a><a class="button" href="index.html">Back to feature studio</a>',
        '<a class="button" href="index.html">Back to M feature studio</a>',
    )
    text = text.replace(
        '<label>Vertical gap <output id="rowGapValue"></output><input id="rowGap" type="range" min="0" max="180" step="5" value="10"></label>',
        '<label>Vertical gap <output id="rowGapValue"></output><input id="rowGap" type="range" min="0" max="180" step="5" value="0"></label>'
        '<label>Left page margin <output id="marginLeftValue"></output><input id="marginLeft" type="range" min="0" max="300" step="5" value="70"></label>'
        '<label>Right page margin <output id="marginRightValue"></output><input id="marginRight" type="range" min="0" max="300" step="5" value="60"></label>'
        '<label>Top page margin <output id="marginTopValue"></output><input id="marginTop" type="range" min="0" max="300" step="5" value="0"></label>'
        '<label>Bottom page margin <output id="marginBottomValue"></output><input id="marginBottom" type="range" min="0" max="300" step="5" value="0"></label>',
    )
    text = text.replace('src="rrr-data.js"', 'src="macaque-m-rrr-data.js?v=20260922"')
    text = text.replace('src="rrr-default-layout.js"', 'src="macaque-m-rrr-default-layout.js?v=20260922"')
    text = text.replace('src="rrr-editor.js"', 'src="macaque-m-rrr-editor.js?v=20260922"')
    (DEST / "rrr-editor.html").write_text(text, encoding="utf-8", newline="\n")


def build_js() -> None:
    text = (SOURCE_UI / "rrr-editor.js").read_text(encoding="utf-8")
    text = text.replace("window.MOUSE_T_M_RRR_DEFAULT_LAYOUT", "window.MACAQUE_T_M_RRR_DEFAULT_LAYOUT")
    text = text.replace("window.MOUSE_M_RRR_DATA", "window.MACAQUE_M_RRR_DATA")
    text = text.replace(
        '    "Maximum trunk angle": "TrunkAng_Max"\n  };',
        '    "Maximum trunk angle": "TrunkAng_Max",\n'
        '    "basal_dendrite_bias_dorsal": "BD-Bias-Dor",\n'
        '    "basal_dendrite_bias_medial": "BD-Bias-Med",\n'
        '    "basal_dendrite_calculate_number_of_stems": "BD-NStem",\n'
        '    "basal_dendrite_extent_dorsal": "BDE-Dor",\n'
        '    "basal_dendrite_extent_medial": "BDE-Med",\n'
        '    "basal_dendrite_max_branch_order": "BD-MaxBO",\n'
        '    "basal_dendrite_max_euclidean_distance": "BD-MaxED",\n'
        '    "basal_dendrite_max_path_distance": "BD-MaxPD",\n'
        '    "basal_dendrite_mean_contraction": "BD-MeanContr",\n'
        '    "basal_dendrite_mean_diameter": "BD-MeanDiam",\n'
        '    "basal_dendrite_num_branches": "BD-NBranch",\n'
        '    "basal_dendrite_soma_percentile_dorsal": "BD-SomaPctl-Dor",\n'
        '    "basal_dendrite_soma_percentile_medial": "BD-SomaPctl-Med",\n'
        '    "basal_dendrite_stem_exit_MedialLateral": "BD-StemExit-ML",\n'
        '    "basal_dendrite_stem_exit_dorsal": "BD-StemExit-Dor",\n'
        '    "basal_dendrite_stem_exit_ventral": "BD-StemExit-Ven",\n'
        '    "basal_dendrite_total_length": "BD-TotLen",\n'
        '    "soma_surface_area": "Soma-SA"\n'
        '  };',
    )
    text = text.replace("Mouse T–M", "Macaque T–M")
    text = text.replace("Mouse T-M", "Macaque T-M")
    text = text.replace("Mouse_T-M_RRR_", "Macaque_T-M_RRR_")
    text = text.replace("fitted Mouse T–M RRR model", "frozen Macaque T–M RRR model")
    text = text.replace("fitted Macaque T–M RRR model", "frozen Macaque T–M RRR model")
    text = text.replace('D1: "#ff379b", D2: "#00eeb3", M1: "#00468B", M2: "#42B540", M3: "#ED0000", M4: "#0099B4"',
                        'D1: "#D95F02", D2: "#008F7A", M1: "#1F77B4", M2: "#D9A400", M3: "#8C564B", M4: "#E377C2"')
    text = text.replace(
        "panelWidth: 440, panelHeight: 370, columnGap: 0, rowGap: 10,",
        "panelWidth: 440, panelHeight: 370, columnGap: 0, rowGap: 10,\n"
        "    marginLeft: 0, marginRight: 0, marginTop: 0, marginBottom: 0,",
    )
    text = text.replace(
        "const x = col * (state.panelWidth + state.columnGap);\n    const y = row * (state.panelHeight + state.rowGap);",
        "const x = state.marginLeft + col * (state.panelWidth + state.columnGap);\n"
        "    const y = state.marginTop + row * (state.panelHeight + state.rowGap);",
    )
    text = text.replace(
        "width: (grid ? 2 : 4) * state.panelWidth + (grid ? 1 : 3) * state.columnGap,\n      height: (grid ? 2 : 1) * state.panelHeight + (grid ? state.rowGap : 0)",
        "width: state.marginLeft + (grid ? 2 : 4) * state.panelWidth + (grid ? 1 : 3) * state.columnGap + state.marginRight,\n"
        "      height: state.marginTop + (grid ? 2 : 1) * state.panelHeight + (grid ? state.rowGap : 0) + state.marginBottom",
    )
    text = text.replace(
        '"panelWidth","panelHeight","columnGap","rowGap","radialPercentile"',
        '"panelWidth","panelHeight","columnGap","rowGap","marginLeft","marginRight","marginTop","marginBottom","radialPercentile"',
    )
    (DEST / "macaque-m-rrr-editor.js").write_text(text, encoding="utf-8", newline="\n")


def build_data() -> None:
    scores = encode(DEST / "rrr_scores.csv")
    loadings = encode(DEST / "rrr_loadings.csv")
    payload = (
        "window.MACAQUE_M_RRR_DATA={"
        f"scores:atob('{scores}'),"
        f"loadings:atob('{loadings}')"
        "};\n"
    )
    (DEST / "macaque-m-rrr-data.js").write_text(payload, encoding="utf-8", newline="\n")


def main() -> None:
    if not REFERENCE.exists():
        raise FileNotFoundError(REFERENCE)
    DEST.mkdir(parents=True, exist_ok=True)
    shutil.copy2(SOURCE_UI / "rrr-editor.css", DEST / "rrr-editor.css")
    layout = build_layout()
    (DEST / "Macaque_T-M_RRR_Macaque-E_reference_layout.json").write_text(
        json.dumps(layout, ensure_ascii=False, indent=2), encoding="utf-8", newline="\n"
    )
    (DEST / "macaque-m-rrr-default-layout.js").write_text(
        "window.MACAQUE_T_M_RRR_DEFAULT_LAYOUT = "
        + json.dumps(layout, ensure_ascii=False, indent=2)
        + ";\n",
        encoding="utf-8",
        newline="\n",
    )
    abbreviation_rows = [
        ("basal_dendrite_extent_dorsal", "BDE-Dor"),
        ("basal_dendrite_extent_medial", "BDE-Med"),
        ("basal_dendrite_max_euclidean_distance", "BD-MaxED"),
        ("basal_dendrite_max_path_distance", "BD-MaxPD"),
        ("basal_dendrite_calculate_number_of_stems", "BD-NStem"),
        ("basal_dendrite_soma_percentile_dorsal", "BD-SomaPctl-Dor"),
        ("basal_dendrite_max_branch_order", "BD-MaxBO"),
        ("basal_dendrite_num_branches", "BD-NBranch"),
        ("basal_dendrite_total_length", "BD-TotLen"),
        ("basal_dendrite_bias_dorsal", "BD-Bias-Dor"),
        ("basal_dendrite_mean_contraction", "BD-MeanContr"),
        ("basal_dendrite_mean_diameter", "BD-MeanDiam"),
        ("soma_surface_area", "Soma-SA"),
        ("basal_dendrite_stem_exit_dorsal", "BD-StemExit-Dor"),
        ("basal_dendrite_soma_percentile_medial", "BD-SomaPctl-Med"),
        ("basal_dendrite_stem_exit_ventral", "BD-StemExit-Ven"),
        ("basal_dendrite_bias_medial", "BD-Bias-Med"),
        ("basal_dendrite_stem_exit_MedialLateral", "BD-StemExit-ML"),
    ]
    (DEST / "Macaque_M_morphology_abbreviations.csv").write_text(
        "raw_feature,display_abbreviation\n"
        + "\n".join(f"{raw},{short}" for raw, short in abbreviation_rows)
        + "\n",
        encoding="utf-8-sig",
        newline="\n",
    )
    build_data()
    build_html()
    build_js()
    print(DEST / "rrr-editor.html")


if __name__ == "__main__":
    main()
