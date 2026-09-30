"""Build the interactive frozen Macaque M1-M4 ten-feature radar editor."""
from pathlib import Path
import re

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
E_SITE = ROOT / "outputs" / "Macaque_E4_feature_explorer_site" / "dist"
M_SITE = ROOT / "interactive" / "macaque-morph-feature-explorer" / "dist"
M_BASE = ROOT / "macaque_m"
RUN = M_BASE / "m18_tempfreeze_NPC5_HCK4_res2.3"
RAW = M_BASE / "m18" / "01_raw_126.csv"
ZSCORE = M_BASE / "m18_adaptive_pca126" / "02_transformed_z_117.csv"
ASSIGN = RUN / "01_temp_frozen_assignments_126.csv"
SELECTION = RUN / "panels_A_I" / "radar_group_top10_feature_selection_117_consensus.csv"
ABBREVIATIONS = M_BASE / "morphology_feature_abbreviations.csv"

EXPECTED = {"M1": 43, "M2": 42, "M3": 20, "M4": 12}
COLORS = {"M1": "#1F77B4", "M2": "#D9A400", "M3": "#8C564B", "M4": "#E377C2"}


def build_data():
    selected = pd.read_csv(SELECTION).sort_values("axis")
    features = selected["feature"].tolist()
    labels = dict(pd.read_csv(ABBREVIATIONS).values)
    assign = pd.read_csv(ASSIGN, dtype={"cell_label": str})
    assign = assign.loc[assign["concordant"].eq(True), ["cell_label", "HC_K4"]].copy()
    assign["M_class"] = "M" + assign["HC_K4"].astype(int).astype(str)
    raw = pd.read_csv(RAW, dtype={"cell_label": str})
    z = pd.read_csv(ZSCORE, dtype={"cell_label": str})
    data = assign[["cell_label", "M_class"]].merge(
        raw[["cell_label"] + features], on="cell_label", validate="one_to_one"
    ).merge(z[["cell_label"] + features], on="cell_label", suffixes=("", "__z"), validate="one_to_one")
    for feature in features:
        lo, hi = np.quantile(data[feature], [0.025, 0.975])
        data[feature + "__scaled"] = np.clip((data[feature] - lo) / (hi - lo), 0, 1)
    counts = data["M_class"].value_counts().reindex(EXPECTED).to_dict()
    if len(data) != 117 or counts != EXPECTED:
        raise RuntimeError(f"Frozen cohort mismatch: n={len(data)}, counts={counts}")
    order = ["cell_label", "M_class"]
    for feature in features:
        order += [feature, feature + "__z", feature + "__scaled"]
    data = data[order]
    data_dir = M_SITE / "data"
    data_dir.mkdir(parents=True, exist_ok=True)
    data.to_csv(data_dir / "macaque_m_radar_cells.csv", index=False)
    return features, [labels[f] for f in features]


def build_html():
    html = (E_SITE / "radar-editor.html").read_text(encoding="utf-8")
    html = html.replace("Macaque E1–E4 radar editor", "Macaque M1–M4 morphology radar editor")
    html = html.replace("frozen Macaque E1–E4 consensus cells", "frozen Macaque M1–M4 consensus cells")
    html = html.replace("Mouse-E reference style · frozen 368-cell consensus cohort", "Frozen 117-cell HC–GC consensus morphology cohort")
    html = html.replace("tsne-editor.html?v=macaque-e4-workflow-20260922", "index.html")
    html = html.replace(">t-SNE workflow<", ">Feature studio<")
    html = html.replace("comparison-editor.html", "index.html")
    html = html.replace(">Feature comparisons<", ">t-SNE comparisons<")
    html = html.replace("Mouse-E profile style", "Profile and shared Z-score space")
    html = html.replace("Mouse-compatible Z-score (−3 to +3)", "Frozen cohort-wide Z-score")
    html = html.replace("Original Macaque percentile scale (0–1)", "Robust raw-value percentile scale (0–1)")
    html = html.replace(
        "</select></label>\n          <label>Center statistic",
        "</select></label>\n"
        "          <label>Z-score lower bound <output id=\"zMinValue\"></output><input id=\"zMin\" type=\"range\" min=\"-6\" max=\"0\" step=\"0.25\"></label>\n"
        "          <label>Z-score upper bound <output id=\"zMaxValue\"></output><input id=\"zMax\" type=\"range\" min=\"0\" max=\"6\" step=\"0.25\"></label>\n"
        "          <label>Center statistic"
    )
    html = html.replace("Reset Mouse-E preset", "Reset frozen M preset")
    html = html.replace(
        "The original ten Macaque indicators are restored. Use the Z-score setting to match the Mouse-E radial rule or the percentile setting to reproduce the earlier Macaque scale. Neither option recalculates E1–E4.",
        "The ten frozen M features and their order are retained. The Z-score bounds alter display mapping only; frozen Z-scores, M1–M4 labels and statistics are never recalculated."
    )
    html = html.replace("Loading 368 consensus cells…", "Loading 117 consensus cells…")
    html = html.replace("Interactive Mouse-E-style radar profiles for Macaque E1 through E4", "Interactive radar profiles for frozen Macaque M1 through M4")
    html = html.replace("styles.css", "radar-editor.css")
    html = re.sub(r"macaque-radar-editor\.js\?v=[^\"]+", "macaque-m-radar-editor.js?v=zrange-20260923", html)
    (M_SITE / "radar-editor.html").write_text(html, encoding="utf-8")
    (M_SITE / "radar-editor.css").write_text((E_SITE / "styles.css").read_text(encoding="utf-8"), encoding="utf-8")


def build_js(features, labels):
    js = (E_SITE / "macaque-radar-editor.js").read_text(encoding="utf-8")
    feature_block = "const features = [\n" + ",\n".join(
        f"  {{id:'{feature}', label:'{label}'}}" for feature, label in zip(features, labels)
    ) + "\n];"
    js = re.sub(r"const features = \[.*?\n\];", feature_block, js, count=1, flags=re.S)
    js = js.replace("macaque-e4-radar-mouse-e-style-default-v5", "macaque-m4-radar-frozen-zrange-default-v1")
    js = js.replace("const classes = ['E1', 'E2', 'E3', 'E4'];", "const classes = ['M1', 'M2', 'M3', 'M4'];")
    js = re.sub(r"const expected = \{.*?\};", "const expected = {M1:43,M2:42,M3:20,M4:12};", js, count=1)
    js = js.replace("colors:{E1:'#F8766D',E2:'#7CAE00',E3:'#00BFC4',E4:'#C77CFF'},",
                    "colors:{M1:'#1F77B4',M2:'#D9A400',M3:'#8C564B',M4:'#E377C2'},")
    js = js.replace("scaleMode:'zscore', centerStatistic:'median'", "scaleMode:'zscore', zMin:-3, zMax:3, centerStatistic:'median'")
    js = js.replace("cellWidth:0.35, cellOpacity:0.08", "cellWidth:0.35, cellOpacity:0.04")
    js = js.replace("centerFillOpacity:0.18", "centerFillOpacity:0.0")
    js = js.replace(
        "function featureValues(row){const suffix=state.scaleMode==='percentile'?'__scaled':'__zscaled';return state.featureOrder.map(id=>+row[id+suffix])}",
        "function featureValue(row,id){if(state.scaleMode==='percentile')return +row[id+'__scaled'];const span=state.zMax-state.zMin;return Math.max(0,Math.min(1,(+row[id+'__z']-state.zMin)/span))}\n"
        "function featureValues(row){return state.featureOrder.map(id=>featureValue(row,id))}"
    )
    js = js.replace("'titleSize','scaleMode','centerStatistic'", "'titleSize','scaleMode','zMin','zMax','centerStatistic'")
    js = js.replace(
        "if(!rows.length)return;\n  sync();",
        "if(!rows.length)return;\n  sync();\n  if(state.scaleMode==='zscore' && !(state.zMin < state.zMax)){$('radarStatus').textContent='Z-score lower bound must be smaller than the upper bound.';return;}"
    )
    js = js.replace(
        "value=state.scaleMode==='zscore'?-3+6*fraction:fraction",
        "value=state.scaleMode==='zscore'?state.zMin+(state.zMax-state.zMin)*fraction:fraction"
    )
    js = js.replace(
        "const columns=state.featureOrder.map(id=>groupRows.map(row=>+row[id+'__scaled']));",
        "const columns=state.featureOrder.map(id=>groupRows.map(row=>featureValue(row,id)));"
    )
    js = js.replace("original ten Macaque indicators", "frozen ten morphology indicators")
    js = js.replace("'Z-score −3 to +3'", "`Z-score ${state.zMin} to ${state.zMax}`")
    js = js.replace("Macaque_E1-E4_radar_MouseE_style", "Macaque_M1-M4_top10_radar")
    js = js.replace("Macaque_E1-E4_radar_frozen368_cells.csv", "Macaque_M1-M4_radar_frozen117_cells.csv")
    js = js.replace("Macaque_E1-E4_radar_layout.json", "Macaque_M1-M4_radar_layout.json")
    js = js.replace("data/macaque_e_radar_cells.csv?v=mouse-e-style-20260922", "data/macaque_m_radar_cells.csv?v=zrange-20260923")
    js = js.replace("if(rows.length!==368)throw Error(`Expected 368 cells, found ${rows.length}`);", "if(rows.length!==117)throw Error(`Expected 117 cells, found ${rows.length}`);")
    js = js.replace("Unable to load Macaque radar data", "Unable to load Macaque M radar data")
    js = js.replace("Select at least one E class.", "Select at least one M class.")
    js = js.replace("r.E_class", "r.M_class")
    (M_SITE / "macaque-m-radar-editor.js").write_text(js, encoding="utf-8")


def add_navigation_links():
    rrr = M_SITE / "rrr-editor.html"
    text = rrr.read_text(encoding="utf-8")
    if "radar-editor.html" not in text:
        text = text.replace(
            '<a href="index.html">Back to M feature studio</a>',
            '<a href="index.html">Back to M feature studio</a> <a href="radar-editor.html">Radar editor</a>'
        )
        rrr.write_text(text, encoding="utf-8")


def main():
    M_SITE.mkdir(parents=True, exist_ok=True)
    features, labels = build_data()
    build_html()
    build_js(features, labels)
    add_navigation_links()
    print(M_SITE / "radar-editor.html")
    print("Features:", ", ".join(labels))


if __name__ == "__main__":
    main()
