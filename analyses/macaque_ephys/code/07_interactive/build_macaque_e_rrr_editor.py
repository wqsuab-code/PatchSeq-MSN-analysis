from __future__ import annotations

import csv
import json
import re
import shutil
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MOUSE = ROOT / "sites" / "ephys-core18-explorer" / "dist"
SOURCE = ROOT / "outputs" / "R3_RRR_fourpanel"
DEST = ROOT / "outputs" / "Macaque_E4_feature_explorer_site" / "dist"
DATA = DEST / "rrr-data"
UPLOADED_LAYOUT = Path(r"C:\Users\53461\Downloads\Macaque_T-E_RRR_gc_merged_layout.json")

T_COLORS = {"D1": "#D95F02", "D2": "#008F7A"}
E_COLORS = {"E1": "#F8766D", "E2": "#7CAE00", "E3": "#00BFC4", "E4": "#C77CFF"}
E_NAMES = {
    "Epsy_width_rheo": "Width rheo",
    "Epsy_fast_trough_v_rheo": "Fast trough",
    "Epsy_peak_deltav_rheo": "Peak ΔV",
    "Epsy_peak_v_rheo": "Peak V",
    "Epsy_postap_slope_rheo": "Post-AP slope",
    "Epsy_threshold_v_rheo": "Threshold V",
    "Epsy_trough_t_rheo": "Trough t",
    "Epsy_trough_v_rheo": "Trough V",
    "Epsy_upstroke_downstroke_ratio_rheo": "Up/down ratio",
    "Epsy_ahp_delay_5spike": "AHP delay",
    "Epsy_ahp_delay_ratio_5spike": "AHP ratio",
    "Epsy_postap_slope_hero": "Post-AP hero",
    "Epsy_trough_t_hero": "Trough t hero",
    "Epsy_downstroke_adapt_ratio": "Downstroke adapt",
    "Epsy_peak_v_adapt_ratio": "Peak V adapt",
    "Epsy_threshold_v_adapt_ratio": "Threshold adapt",
    "Epsy_upstroke_adapt_ratio": "Upstroke adapt",
    "Epsy_width_adapt_ratio": "Width adapt",
    "Epsy_threshold_v_short_square": "Threshold SS",
}


def write_csv(path: Path, rows: list[dict], fields: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def build_data() -> None:
    DATA.mkdir(parents=True, exist_ok=True)
    with (SOURCE / "RRR_cell_scores_n346.csv").open(encoding="utf-8-sig", newline="") as handle:
        source_rows = list(csv.DictReader(handle))
    score_rows = []
    for row in source_rows:
        score_rows.append({
            "MSN_unique_ID": row["cell_label"],
            "D1_D2": row["T_class"],
            "E_class": row["E_class"].replace("C", "E", 1),
            **{f"T_Component{i}": row[f"T_Component{i}"] for i in range(1, 4)},
            **{f"E_Component{i}": row[f"E_Component{i}"] for i in range(1, 4)},
        })
    score_fields = ["MSN_unique_ID", "D1_D2", "E_class"] + [f"T_Component{i}" for i in range(1, 4)] + [f"E_Component{i}" for i in range(1, 4)]
    write_csv(DATA / "rrr_scores_gc_merged.csv", score_rows, score_fields)
    write_csv(DATA / "rrr_scores_gc_unmerged.csv", score_rows, score_fields)

    with (SOURCE / "RRR_gene_correlation_loadings.csv").open(encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.reader(handle))
    genes = [{"Gene": row[0], "Component1": row[1], "Component2": row[2], "Component3": row[3]} for row in rows[1:]]
    write_csv(DATA / "rrr_gene_loadings_gc_merged.csv", genes, ["Gene", "Component1", "Component2", "Component3"])
    write_csv(DATA / "rrr_gene_loadings_gc_unmerged.csv", genes, ["Gene", "Component1", "Component2", "Component3"])

    with (SOURCE / "RRR_Efeature_correlation_loadings.csv").open(encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.reader(handle))
    features = [{"Display": E_NAMES.get(row[0], row[0].removeprefix("Epsy_").replace("_", " ")),
                 "Component1": row[1], "Component2": row[2], "Component3": row[3]} for row in rows[1:]]
    write_csv(DATA / "rrr_feature_loadings_gc_merged.csv", features, ["Display", "Component1", "Component2", "Component3"])
    write_csv(DATA / "rrr_feature_loadings_gc_unmerged.csv", features, ["Display", "Component1", "Component2", "Component3"])
    write_csv(DATA / "rrr_gc_versions_manifest.csv", [
        {"version_key": "gc_merged", "agreement_n": "368", "agreement_denominator": "390"},
        {"version_key": "gc_unmerged", "agreement_n": "368", "agreement_denominator": "390"},
    ], ["version_key", "agreement_n", "agreement_denominator"])


def build_layout() -> None:
    layout = json.loads(UPLOADED_LAYOUT.read_text(encoding="utf-8-sig"))
    layout["datasetKey"] = "gc_merged"
    layout.setdefault("marginLeft", 0)
    layout.setdefault("marginRight", 0)
    layout.setdefault("marginTop", 0)
    layout.setdefault("marginBottom", 0)
    layout["colors"] = {**T_COLORS, **E_COLORS}
    # Preserve the uploaded geometry, typography, component order and label counts.
    (DATA / "Macaque_T-E_RRR_layout_reference.json").write_text(
        json.dumps(layout, indent=2, ensure_ascii=False), encoding="utf-8"
    )


def build_page() -> None:
    shutil.copy2(MOUSE / "rrr-editor.css", DEST / "rrr-editor.css")
    html = (MOUSE / "rrr-editor.html").read_text(encoding="utf-8")
    html = html.replace("Mouse T–E RRR figure editor", "Macaque T–E RRR figure editor")
    html = html.replace("Interactive editor for the Mouse transcriptomic-to-electrophysiological reduced-rank regression figure.",
                        "Interactive editor for the Macaque transcriptomic-to-electrophysiological reduced-rank regression figure.")
    html = html.replace("GC–HC consensus workflow · loading analysis version…", "Ca + Pu + NAc frozen E1–E4 workflow · loading…")
    html = html.replace('<a class="button" href="tsne-editor.html">t-SNE workflow</a><a class="button" href="stability-statistics.html">E stability statistics</a><a class="button" href="composition-editor.html">Strict T-identity composition</a><a class="button" href="index.html">Back to feature studio</a>',
                        '<a class="button" href="comparison-editor.html">t-SNE + feature comparisons</a><a class="button" href="ml-editor.html">Machine-learning editor</a><a class="button" href="index.html">Back to feature studio</a>')
    html = re.sub(r'<label>Clustering version<select id="datasetVersion">.*?</select></label>',
                  '<label>Frozen cohort<select id="datasetVersion"><option value="gc_merged">Macaque E1–E4 · Ca + Pu + NAc</option></select></label>', html)
    html = html.replace("Each option loads an independently selected GC–HC consensus cohort and an independently refitted rank-3 T–E RRR model.",
                        "Frozen HC–GC consensus E1–E4 cohort; Hybrid cells are excluded from the D1/D2 RRR display.")
    html = html.replace("1 × 4 · Mouse_E reference", "1 × 4 · uploaded RRR reference")
    html = html.replace(
        '<label>Vertical gap <output id="rowGapValue"></output><input id="rowGap" type="range" min="0" max="180" step="5" value="10"></label>',
        '<label>Vertical gap <output id="rowGapValue"></output><input id="rowGap" type="range" min="0" max="180" step="5" value="10"></label>'
        '<label>Left page margin <output id="marginLeftValue"></output><input id="marginLeft" type="range" min="0" max="300" step="5" value="0"></label>'
        '<label>Right page margin <output id="marginRightValue"></output><input id="marginRight" type="range" min="0" max="300" step="5" value="0"></label>'
        '<label>Top page margin <output id="marginTopValue"></output><input id="marginTop" type="range" min="0" max="300" step="5" value="0"></label>'
        '<label>Bottom page margin <output id="marginBottomValue"></output><input id="marginBottom" type="range" min="0" max="300" step="5" value="0"></label>'
    )
    html = html.replace("Save current parameters as default", "Set current settings as temporary default")
    html = html.replace("Restore saved default", "Restore temporary default")
    html = html.replace("Live Mouse_E preview", "Live Macaque-E preview")
    html = html.replace("Loading the frozen Mouse T–E RRR scores and loadings…", "Loading frozen Macaque T–E RRR scores and loadings…")
    html = html.replace("Interactive Mouse transcriptomic-to-electrophysiological reduced-rank regression figure", "Interactive Macaque transcriptomic-to-electrophysiological reduced-rank regression figure")
    html = html.replace('rrr-editor.js?v=layout-upload2-20260922', 'macaque-rrr-editor.js?v=macaque-e4-margins-20260922')
    (DEST / "rrr-editor.html").write_text(html, encoding="utf-8")

    js = (MOUSE / "rrr-editor.js").read_text(encoding="utf-8")
    js = js.replace('const STORAGE_KEY = "mouse-t-e-rrr-default-layout-upload2-v1";', 'const STORAGE_KEY = "macaque-e4-t-e-rrr-reference-layout-v2";')
    js = js.replace('const DATA_VERSION = "layout-upload2-20260922";', 'const DATA_VERSION = "macaque-e4-margins-20260922";')
    js = re.sub(r'const DATASETS = \{.*?\n  \};', 'const DATASETS = {\n    gc_merged: {label: "Macaque E1–E4 · Ca + Pu + NAc", stem: "macaque_e4", classLabel: "E class"}\n  };', js, count=1, flags=re.S)
    js = js.replace('const MERGED_GROUPS = ["E1", "E2", "E3", "E4", "E5"];', 'const MERGED_GROUPS = ["E1", "E2", "E3", "E4"];')
    js = js.replace('D1: "#E41A1C", D2: "#377EB8",\n    E1: "#4E79A7", E2: "#F28E2B", E3: "#59A14F", E4: "#2AA6B8", E5: "#B07AA1",',
                    'D1: "#D95F02", D2: "#008F7A",\n    E1: "#F8766D", E2: "#7CAE00", E3: "#00BFC4", E4: "#C77CFF",')
    js = js.replace(
        'panelWidth: 440, panelHeight: 370, columnGap: 0, rowGap: 10,',
        'panelWidth: 440, panelHeight: 370, columnGap: 0, rowGap: 10,\n    marginLeft: 0, marginRight: 0, marginTop: 0, marginBottom: 0,'
    )
    js = js.replace(
        'const x = col * (state.panelWidth + state.columnGap);\n    const y = row * (state.panelHeight + state.rowGap);',
        'const x = state.marginLeft + col * (state.panelWidth + state.columnGap);\n    const y = state.marginTop + row * (state.panelHeight + state.rowGap);'
    )
    js = js.replace(
        'width: (grid ? 2 : 4) * state.panelWidth + (grid ? 1 : 3) * state.columnGap,\n      height: (grid ? 2 : 1) * state.panelHeight + (grid ? state.rowGap : 0)',
        'width: state.marginLeft + (grid ? 2 : 4) * state.panelWidth + (grid ? 1 : 3) * state.columnGap + state.marginRight,\n      height: state.marginTop + (grid ? 2 : 1) * state.panelHeight + (grid ? state.rowGap : 0) + state.marginBottom'
    )
    js = js.replace(
        '"panelWidth","panelHeight","columnGap","rowGap","radialPercentile"',
        '"panelWidth","panelHeight","columnGap","rowGap","marginLeft","marginRight","marginTop","marginBottom","radialPercentile"'
    )
    js = js.replace('"E1–E5"', '"E1–E4"')
    js = js.replace('strict T-stable cells', 'frozen D1/D2 cells')
    js = js.replace('Primary analysis: GC merged K=5 matched to HC K=5; strict D1/D2 cells are shown as E1–E5.',
                    'Frozen Macaque E1–E4 analysis: Ca + Pu + NAc HC–GC consensus cells with D1/D2 identity; Hybrid cells are excluded.')
    js = js.replace('loading vectors remain fixed to this independently fitted Mouse T–E RRR model.',
                    'loading vectors remain fixed to the frozen Macaque T–E RRR model.')
    js = js.replace('Mouse_T-E_RRR_', 'Macaque_T-E_RRR_')
    js = js.replace('data/rrr_', 'rrr-data/rrr_')
    js = js.replace('data/Mouse_T-E_RRR_layout_09172026.json', 'rrr-data/Macaque_T-E_RRR_layout_reference.json')
    js = js.replace('data/Macaque_T-E_RRR_layout_09172026.json', 'rrr-data/Macaque_T-E_RRR_layout_reference.json')
    (DEST / "macaque-rrr-editor.js").write_text(js, encoding="utf-8")

    comparison = (DEST / "comparison-editor.html").read_text(encoding="utf-8")
    if 'href="rrr-editor.html"' not in comparison:
        comparison = comparison.replace('<a class="button" href="ml-editor.html">Machine-learning editor</a>',
                                        '<a class="button" href="rrr-editor.html">RRR editor</a><a class="button" href="ml-editor.html">Machine-learning editor</a>')
        (DEST / "comparison-editor.html").write_text(comparison, encoding="utf-8")


def main() -> None:
    DEST.mkdir(parents=True, exist_ok=True)
    build_data()
    build_layout()
    build_page()
    print(DEST / "rrr-editor.html")


if __name__ == "__main__":
    main()
