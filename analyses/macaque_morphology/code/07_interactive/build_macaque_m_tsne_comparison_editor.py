"""Build the Macaque M morphology t-SNE + feature-comparison editor.

The UI and default layout are cloned from the final E-feature comparison editor,
while every analytical datum is replaced by the frozen Macaque M1-M4 cohort.
"""
from __future__ import annotations

import json
import math
import re
import shutil
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
REFERENCE = ROOT / "sites" / "ephys-core18-explorer" / "dist"
SOURCE = ROOT / "macaque_m" / "m18_tempfreeze_NPC5_HCK4_res2.3" / "M4_compare_all18"
PCA_LOADINGS = ROOT / "macaque_m" / "m18_adaptive_pca126" / "05_pca_loadings.csv"
_layout_dir = Path(r"C:\Users\53461\Downloads")
_layout_candidates = sorted(
    _layout_dir.glob("Macaque_M_morphology_tsne_comparison_layout*.json"),
    key=lambda path: path.stat().st_mtime,
    reverse=True,
)
USER_LAYOUT = _layout_candidates[0] if _layout_candidates else _layout_dir / "Macaque_M_morphology_tsne_comparison_layout.json"
IMPORTANCE_MAIN = ROOT / "outputs" / "Macaque_M4_full_ML500" / "tables" / "20_feature_importance_MDI_logistic_500.csv"
IMPORTANCE_PERM = ROOT / "outputs" / "Macaque_M4_full_ML500" / "tables" / "21_nested_grouped_permutation_importance_500.csv"
OUT = ROOT / "interactive" / "macaque-morph-tsne-comparison-editor" / "dist"

FEATURES = [
    ("basal_dendrite_total_length", "Total length", "µm"),
    ("basal_dendrite_num_branches", "Number of branches", "count"),
    ("basal_dendrite_calculate_number_of_stems", "Number of stems", "count"),
    ("basal_dendrite_max_branch_order", "Max branch order", "order"),
    ("basal_dendrite_max_path_distance", "Max path distance", "µm"),
    ("basal_dendrite_max_euclidean_distance", "Max Euclidean distance", "µm"),
    ("basal_dendrite_mean_contraction", "Mean contraction", "ratio"),
    ("basal_dendrite_mean_diameter", "Mean diameter", "µm"),
    ("basal_dendrite_extent_dorsal", "Dorsal extent", "µm"),
    ("basal_dendrite_extent_medial", "Medial extent", "µm"),
    ("basal_dendrite_bias_dorsal", "Dorsal bias", "a.u."),
    ("basal_dendrite_bias_medial", "Medial bias", "a.u."),
    ("basal_dendrite_soma_percentile_dorsal", "Dorsal soma percentile", ""),
    ("basal_dendrite_soma_percentile_medial", "Medial soma percentile", ""),
    ("basal_dendrite_stem_exit_dorsal", "Dorsal stem exit", "a.u."),
    ("basal_dendrite_stem_exit_ventral", "Ventral stem exit", "a.u."),
    ("basal_dendrite_stem_exit_MedialLateral", "ML stem exit", "a.u."),
    ("soma_surface_area", "Soma surface area", "µm²"),
]
GROUPS = ["M1", "M2", "M3", "M4"]
COLORS = {"M1": "#1F77B4", "M2": "#D9A400", "M3": "#8C564B", "M4": "#E377C2"}
CORE6 = [key for key, _, _ in FEATURES[:6]]


def bh_adjust(values: list[float]) -> list[float]:
    n = len(values)
    order = sorted(range(n), key=lambda i: values[i])
    adjusted = [1.0] * n
    running = 1.0
    for rank_index in range(n - 1, -1, -1):
        original_index = order[rank_index]
        rank = rank_index + 1
        running = min(running, values[original_index] * n / rank)
        adjusted[original_index] = min(1.0, running)
    return adjusted


def safe_number(value):
    if pd.isna(value):
        return None
    if isinstance(value, bool):
        return value
    if hasattr(value, "item"):
        value = value.item()
    return value


def records(frame: pd.DataFrame) -> list[dict]:
    return [{key: safe_number(value) for key, value in row.items()} for row in frame.to_dict("records")]


def build_payload() -> dict:
    source = pd.read_csv(SOURCE / "all18_plot_data_117.csv")
    feature_keys = [item[0] for item in FEATURES]
    pca_features = pd.read_csv(PCA_LOADINGS)["feature"].astype(str).tolist()
    if len(pca_features) != 18 or set(feature_keys) != set(pca_features):
        extra = sorted(set(feature_keys) - set(pca_features))
        missing = sorted(set(pca_features) - set(feature_keys))
        raise RuntimeError(
            "Comparison features no longer match the frozen PCA whitelist; "
            f"non-PCA={extra}, omitted-PCA={missing}"
        )

    importance = pd.read_csv(IMPORTANCE_MAIN).merge(
        pd.read_csv(IMPORTANCE_PERM)[["feature", "importance_mean_median"]],
        on="feature",
        how="inner",
        validate="one_to_one",
    )
    if set(importance["feature"]) != set(feature_keys):
        raise RuntimeError("ML importance tables do not match the frozen 18-feature PCA panel")
    importance["rank_ExtraTrees_MDI"] = importance["ExtraTrees_MDI_median"].rank(ascending=False, method="average")
    importance["rank_logistic_abs_coef"] = importance["Logistic_abs_standardized_coef_median"].rank(ascending=False, method="average")
    importance["rank_grouped_permutation"] = importance["importance_mean_median"].rank(ascending=False, method="average")
    importance["consensus_mean_rank"] = importance[
        ["rank_ExtraTrees_MDI", "rank_logistic_abs_coef", "rank_grouped_permutation"]
    ].mean(axis=1)
    importance = importance.sort_values(
        ["consensus_mean_rank", "ExtraTrees_MDI_median", "feature"],
        ascending=[True, False, True],
    ).reset_index(drop=True)
    importance.insert(1, "panel_rank", range(1, len(importance) + 1))
    importance_order = importance["feature"].tolist()
    required = ["cell_label", "M", "x", "y", *feature_keys, *[f"{key}__z" for key in feature_keys]]
    missing = [column for column in required if column not in source]
    if missing:
        raise RuntimeError(f"Missing frozen columns: {missing}")
    if len(source) != 117 or source["M"].value_counts().sort_index().to_dict() != {"M1": 43, "M2": 42, "M3": 20, "M4": 12}:
        raise RuntimeError("Frozen Macaque M cohort or class counts changed")

    raw = source[["cell_label", "M", *feature_keys]].rename(
        columns={"cell_label": "MSN_unique_ID", "M": "E_type"}
    )

    pairwise = pd.read_csv(SOURCE / "all18_Dunn_Holm.csv")
    raw_ps = pairwise["P_raw"].astype(float).tolist()
    pairwise["P_adj_BH_global_108_tests"] = bh_adjust(raw_ps)
    feature_meta = {key: (index + 1, label, unit) for index, (key, label, unit) in enumerate(FEATURES)}
    enriched = []
    for row in pairwise.to_dict("records"):
        key = row["Feature"]
        group_1, group_2 = row["Group_1"], row["Group_2"]
        values_1 = source.loc[source["M"] == group_1, key].dropna().astype(float)
        values_2 = source.loc[source["M"] == group_2, key].dropna().astype(float)
        order, label, unit = feature_meta[key]
        enriched.append(
            {
                "Order": order,
                "Feature": key,
                "Display_name": label,
                "Unit": unit,
                "Group_1": group_1,
                "Group_2": group_2,
                "N_1": len(values_1),
                "N_2": len(values_2),
                "Dunn_Z": float(row["Dunn_Z"]),
                "Dunn_r": float(row["Dunn_r"]),
                "P_raw": float(row["P_raw"]),
                "P_adj_Holm_within_feature": float(row["P_adj_Holm_within_feature"]),
                "P_adj_BH_global_108_tests": float(row["P_adj_BH_global_108_tests"]),
                "Median_1": float(values_1.median()),
                "Median_2": float(values_2.median()),
                "Median_difference_1_minus_2": float(values_1.median() - values_2.median()),
            }
        )

    omnibus_source = pd.read_csv(SOURCE / "all18_KruskalWallis_BH.csv")
    omnibus = []
    for row in omnibus_source.to_dict("records"):
        omnibus.append(
            {
                "Feature": row["Feature"],
                "Display_name": row["Display"],
                "Unit": row["Unit"] if not pd.isna(row["Unit"]) else "",
                "Kruskal_H": float(row["Kruskal_H"]),
                "df": int(row["df"]),
                "P_raw": float(row["P_raw"]),
                "P_adj_BH_across_18_features": float(row["Q_BH_across_18_features"]),
            }
        )

    tsne = source[["cell_label", "M", "x", "y"]].rename(
        columns={
            "cell_label": "MSN_unique_ID",
            "M": "HC_GC_consensus_E",
            "x": "tSNE_1",
            "y": "tSNE_2",
        }
    )
    tsne.insert(1, "HC_GC_consensus", "True")
    z_long = source.melt(
        id_vars=["cell_label"],
        value_vars=[f"{key}__z" for key in feature_keys],
        var_name="Feature_full",
        value_name="Z_score",
    )
    z_long["Feature_full"] = z_long["Feature_full"].str.replace(r"__z$", "", regex=True)
    z_long = z_long.rename(columns={"cell_label": "Cell_ID"})

    x_min, x_max = float(tsne["tSNE_1"].min()), float(tsne["tSNE_1"].max())
    y_min, y_max = float(tsne["tSNE_2"].min()), float(tsne["tSNE_2"].max())
    x_pad = (x_max - x_min) * 0.05
    y_pad = (y_max - y_min) * 0.05
    layout = {
        "pairLayout": "side-tsne-first",
        "pairGap": 60,
        "uniformYAxisHeight": True,
        "statsAspectRatio": 0,
        "columns": 4,
        "panelWidth": 800,
        "panelHeight": 620,
        "panelGapX": 0,
        "panelGapY": 19,
        "canvasColor": "#ffffff",
        "showPanelTitles": True,
        "panelTitleSize": 25,
        "showPoints": True,
        "pointSize": 3,
        "pointOpacity": 0.75,
        "jitterWidth": 0.4,
        "showBoxes": True,
        "boxWidth": 0.62,
        "boxLineWidth": 1,
        "boxFillOpacity": 0,
        "whiskerWidth": 0.25,
        "medianLineWidth": 2,
        "showP": True,
        "pSource": "P_adj_Holm_within_feature",
        "pCutoff": 0.05,
        "pFontSize": 25,
        "pLineWidth": 1,
        "pRowGap": 25,
        "pLabelGap": 5,
        "pNotation": "scientific",
        "showOmnibus": False,
        "showAxes": True,
        "showTicks": True,
        "showTickLabels": True,
        "showXLabels": True,
        "showYLabels": True,
        "axisLineWidth": 1,
        "tickLength": 7,
        "tickCount": 5,
        "axisFontSize": 25,
        "yLabelSize": 25,
        "yLabelOffset": 46,
        "totalLengthAxisDivisor": 1000,
        "fontFamily": "Arial",
        "dataPadding": 0.15,
        "showTsne": True,
        "tsneFraction": 0.4,
        "tsnePointSize": 8,
        "tsnePointOpacity": 1,
        "tsneOutline": True,
        "tsneOutlineWidth": 2.5,
        "tsneOutlineOpacity": 1,
        "tsneDash": "10 5",
        "tsneBoundaryCoverage": 0.80,
        "tsneBoundaryExpansion": 1.03,
        "tsneBoundarySmoothing": 3,
        "tsneShowLabels": True,
        "tsneShowAxes": False,
        "tsneShowColorbar": False,
        "tsneZMin": -2.5,
        "tsneZMax": 2.5,
        "tsneLowColor": "#2166ac",
        "tsneMidColor": "#d9d9d9",
        "tsneHighColor": "#b2182b",
        "tsneXMin": round(x_min - x_pad, 3),
        "tsneXMax": round(x_max + x_pad, 3),
        "tsneYMin": round(y_min - y_pad, 3),
        "tsneYMax": round(y_max + y_pad, 3),
        "tsneUnitRatio": 1.5,
        "colors": COLORS,
        "order": feature_keys,
        "visible": {key: True for key in feature_keys},
        "labels": {key: {"name": label, "unit": unit, "min": "", "max": ""} for key, label, unit in FEATURES},
        "exportWidthIn": 7.2,
        "exportDpi": 600,
        "exportBackground": "white",
    }
    if USER_LAYOUT.is_file():
        layout.update(json.loads(USER_LAYOUT.read_text(encoding="utf-8-sig")))
    layout.update(
        {
            "colors": COLORS,
            "order": importance_order,
            "visible": {key: True for key in feature_keys},
        }
    )
    layout.setdefault("tsneBoundaryCoverage", 0.80)
    layout.setdefault("tsneBoundaryExpansion", 1.03)
    layout.setdefault("tsneBoundarySmoothing", 3)
    layout.setdefault("totalLengthAxisDivisor", 1000)
    return {
        "rawData": records(raw),
        "pairData": enriched,
        "omnibusData": omnibus,
        "tsneText": tsne.to_csv(index=False, lineterminator="\n"),
        "zText": z_long.to_csv(index=False, lineterminator="\n"),
        "savedLayout": layout,
        "importanceTable": records(importance),
        "importanceOrder": importance_order,
    }


def build() -> Path:
    OUT.mkdir(parents=True, exist_ok=True)
    shutil.copy2(REFERENCE / "styles.css", OUT / "styles.css")
    payload = build_payload()
    pd.DataFrame(payload["importanceTable"]).to_csv(OUT / "FEATURE_IMPORTANCE_PANEL_ORDER.csv", index=False)
    pd.DataFrame(payload["pairData"]).to_csv(OUT / "STATISTICS_pairwise_Dunn_Holm_BH108.csv", index=False)
    pd.DataFrame(payload["omnibusData"]).to_csv(OUT / "STATISTICS_omnibus_KruskalWallis_BH18.csv", index=False)
    (OUT / "data.js").write_text(
        "window.MACAQUE_M_DATA=" + json.dumps(payload, ensure_ascii=False, separators=(",", ":")) + ";\n",
        encoding="utf-8",
    )

    feature_js = "const features=[\n" + ",\n".join(
        f" {{key:{json.dumps(key)},name:{json.dumps(label)},unit:{json.dumps(unit)}}}"
        for key, label, unit in FEATURES
    ) + "\n];"
    js = (REFERENCE / "comparison-editor.js").read_text(encoding="utf-8")
    js = re.sub(r"const NS='http://www.w3.org/2000/svg',groups=\[.*?\];", "const NS='http://www.w3.org/2000/svg',groups=['M1','M2','M3','M4'];", js, count=1)
    js = re.sub(r"const features=\[.*?\];\nconst core6=", lambda _: feature_js + "\nconst core6=", js, count=1, flags=re.S)
    js = re.sub(r"const core6=\[.*?\],counts=\{.*?\};", "const core6=" + json.dumps(payload["importanceOrder"][0:6]) + ",counts={M1:43,M2:42,M3:20,M4:12};", js, count=1)
    js = re.sub(r"colors:\{E1:'.*?'\s*,E2:'.*?'\s*,E3:'.*?'\s*,E4:'.*?'\s*,E5:'.*?'\}", "colors:{M1:'#1F77B4',M2:'#D9A400',M3:'#8C564B',M4:'#E377C2'}", js, count=1)
    js = js.replace(
        "tsneOutlineOpacity:1,tsneDash:'4 3'",
        "tsneOutlineOpacity:1,tsneDash:'4 3',tsneBoundaryCoverage:.8,tsneBoundaryExpansion:1.03,tsneBoundarySmoothing:3",
    )
    js = js.replace("yLabelOffset:46,fontFamily", "yLabelOffset:46,totalLengthAxisDivisor:1000,fontFamily")
    js = js.replace("(right-left)/5", "(right-left)/groups.length")
    js = js.replace("*(right-left)/5", "*(right-left)/groups.length")
    js = js.replace("${list.length}/18 features · 450 cells", "${list.length}/18 features · 117 cells")
    js = js.replace("P_adj_BH_global_180_tests", "P_adj_BH_global_108_tests")
    js = js.replace("Mouse_E1-E5_final_features_comparison", "Macaque_M1-M4_morphology_tsne_comparison")
    js = js.replace("Mouse_E_comparison_", "Macaque_M_morphology_tsne_comparison_")
    js = js.replace("Mouse_E_selected_features_pairwise_statistics", "Macaque_M_selected_features_pairwise_statistics")
    js = js.replace("Mouse_E_comparison_layout", "Macaque_M_tsne_comparison_layout")
    js = js.replace("[{Group_1:'E1',Group_2:'E2'},{Group_1:'E3',Group_2:'E4'}]", "[{Group_1:'M1',Group_2:'M2'},{Group_1:'M3',Group_2:'M4'}]")
    js = js.replace(
        "cut=quantile(scored.map(x=>x.d),.95),poly=hull(scored.filter(x=>x.d<=cut).map(x=>x.p));return{poly:chaikin(poly.map(p=>[mx+(p[0]-mx)*1.03,my+(p[1]-my)*1.03])),center:[mx,my]}",
        "cut=quantile(scored.map(x=>x.d),state.tsneBoundaryCoverage),poly=hull(scored.filter(x=>x.d<=cut).map(x=>x.p)),expanded=poly.map(p=>[mx+(p[0]-mx)*state.tsneBoundaryExpansion,my+(p[1]-my)*state.tsneBoundaryExpansion]);return{poly:chaikin(expanded,state.tsneBoundarySmoothing),center:[mx,my]}",
    )
    js = js.replace(
        "function zColor(v)",
        "function updateRegions(){regions={};for(const name of groups)regions[name]=makeRegion(tsne.filter(r=>r.E_type===name).map(r=>[r.x,r.y]))}\nfunction zColor(v)",
    )
    js = js.replace("function render(){if(!raw.length)return;sync();", "function render(){if(!raw.length)return;updateRegions();sync();")
    js = js.replace(
        ",g=el('g',{'data-feature':feature.key});g.append",
        ",displayDivisor=feature.key==='basal_dendrite_total_length'?state.totalLengthAxisDivisor:1,g=el('g',{'data-feature':feature.key});g.append",
    )
    js = js.replace(
        "'font-size':state.axisFontSize},numberText(t)))",
        "'font-size':state.axisFontSize},numberText(t/displayDivisor)))",
    )
    js = js.replace(
        "const cfg=state.labels[feature.key],label=cfg.unit&&cfg.unit!=='index'&&cfg.unit!=='CV'?`${cfg.name} (${cfg.unit})`:cfg.name,labelX=left-state.yLabelOffset;",
        "const cfg=state.labels[feature.key],scaledUnit=displayDivisor===1000?`×10³ ${cfg.unit}`:cfg.unit,label=scaledUnit&&scaledUnit!=='index'&&scaledUnit!=='CV'?`${cfg.name} (${scaledUnit})`:cfg.name,labelX=left-state.yLabelOffset;",
    )
    js = js.replace(
        "tsneOutlineOpacity:'num',tsneDash:'str'",
        "tsneOutlineOpacity:'num',tsneDash:'str',tsneBoundaryCoverage:'num',tsneBoundaryExpansion:'num',tsneBoundarySmoothing:'int'",
    )
    js = js.replace("yLabelOffset:'num',fontFamily", "yLabelOffset:'num',totalLengthAxisDivisor:'num',fontFamily")
    fetch_pattern = re.compile(
        r"const \[rawData,pairData,omnibusData,tsneText,zText,savedLayout\]=await Promise\.all\(\[fetch\('data/raw\.json'\).*?fetch\('data/comparison_default_layout\.json\?v=y-title-distance-20260920'\)\.then\(r=>r\.json\(\)\)\]\);",
        flags=re.S,
    )
    js, replaced = fetch_pattern.subn(
        "const {rawData,pairData,omnibusData,tsneText,zText,savedLayout}=window.MACAQUE_M_DATA;",
        js,
        count=1,
    )
    if replaced != 1:
        raise RuntimeError("Could not replace reference data loader")
    js = js.replace("if(raw.length!==450||zById.size!==450||tsne.length!==450", "if(raw.length!==117||zById.size!==117||tsne.length!==117")
    js = js.replace("Expected 450 matched cells with 18 frozen Z-scores", "Expected 117 matched cells with 18 frozen morphology Z-scores")
    (OUT / "comparison-editor.js").write_text(js, encoding="utf-8")

    html = (REFERENCE / "comparison-editor.html").read_text(encoding="utf-8")
    html = html.replace("E-feature group-comparison editor", "Macaque M morphology t-SNE + comparison editor")
    html = html.replace("Editable combined E1-E5 comparisons for the final mouse electrophysiology features.", "Editable comparisons of 18 frozen morphology features across Macaque M1-M4.")
    html = html.replace("450 GC–HC consensus cells · raw values · Holm-adjusted Dunn comparisons", "117 HC–GC consensus Macaque MSN cells · raw morphology values · Holm-adjusted Dunn comparisons")
    html = re.sub(r'<div class="header-actions">.*?</div>', '<div class="header-actions"><a class="button" href="../../macaque-morph-feature-explorer/dist/index.html">M morphology explorer</a><a class="button" href="../../macaque-morph-ml-studio/dist/index.html">ML validation</a></div>', html, count=1)
    html = html.replace("E1–E5 boundaries", "M1–M4 boundaries")
    html = html.replace("E labels", "M labels")
    html = html.replace("E1–E5 labels", "M1–M4 labels")
    html = html.replace("BH across 180 tests", "BH across 108 tests")
    html = html.replace("Core 6", "Top 6 importance")
    html = html.replace(
        '<input id="tsnePointSize" type="range" min=".2" max="8" step=".1" value="2.2">',
        '<input id="tsnePointSize" type="range" min=".2" max="20" step=".1" value="2.2">',
    )
    html = html.replace(
        '<label>Dash pattern<input id="tsneDash" type="text" value="4 3"></label>',
        '<label>Dash pattern<input id="tsneDash" type="text" value="4 3"></label>'
        '<label>Boundary coverage <output id="tsneBoundaryCoverageValue"></output><input id="tsneBoundaryCoverage" type="range" min=".50" max=".99" step=".01" value=".80"></label>'
        '<label>Boundary expansion <output id="tsneBoundaryExpansionValue"></output><input id="tsneBoundaryExpansion" type="range" min=".85" max="1.30" step=".01" value="1.03"></label>'
        '<label>Smoothing passes <output id="tsneBoundarySmoothingValue"></output><input id="tsneBoundarySmoothing" type="range" min="0" max="6" step="1" value="3"></label>',
    )
    html = html.replace(
        '<label>Y-title distance <output id="yLabelOffsetValue"></output><input id="yLabelOffset" type="range" min="0" max="140" step="1" value="46"></label>',
        '<label>Y-title distance <output id="yLabelOffsetValue"></output><input id="yLabelOffset" type="range" min="0" max="140" step="1" value="46"></label>'
        '<label>Total length Y-axis divisor<select id="totalLengthAxisDivisor"><option value="1">1</option><option value="1000" selected>1,000</option></select></label>',
    )
    html = html.replace("Editable combined E-feature comparisons", "Editable Macaque M morphology t-SNE and feature comparisons")
    html = html.replace("Loading frozen data and adjusted statistics…", "Loading frozen Macaque M morphology data and adjusted statistics…")
    html = html.replace('<script type="module" src="comparison-editor.js?v=markers-bottom-up-20260920"></script>', '<script src="data.js"></script><script type="module" src="comparison-editor.js"></script>')
    (OUT / "comparison-editor.html").write_text(html, encoding="utf-8")
    (OUT / "index.html").write_text(html, encoding="utf-8")

    audit = {
        "cohort": "Macaque MSN HC-GC consensus morphology cells",
        "n": 117,
        "counts": {"M1": 43, "M2": 42, "M3": 20, "M4": 12},
        "features": [key for key, _, _ in FEATURES],
        "feature_whitelist_source": str(PCA_LOADINGS),
        "feature_rule": "Only the 18 rows present in the frozen PCA loading table are included",
        "panel_order_rule": "Ascending consensus mean rank across Extra Trees MDI, absolute standardized multinomial-logistic coefficients, and nested donor-grouped permutation importance; each component is ranked from highest to lowest importance.",
        "panel_order_file": str(OUT / "FEATURE_IMPORTANCE_PANEL_ORDER.csv"),
        "importance_sources": [str(IMPORTANCE_MAIN), str(IMPORTANCE_PERM)],
        "class_colors": COLORS,
        "raw_value_source": str(SOURCE / "all18_plot_data_117.csv"),
        "pairwise_source": str(SOURCE / "all18_Dunn_Holm.csv"),
        "omnibus_source": str(SOURCE / "all18_KruskalWallis_BH.csv"),
        "tsne_coordinates": "Frozen x/y columns in all18_plot_data_117.csv",
        "tsne_feature_colors": "Frozen transformed feature Z scores clipped visually to [-2.5, 2.5]",
        "boundary": "80% robust-distance subset, convex hull, Chaikin smoothing",
        "boundary_definition": "Within each M class, t-SNE coordinates are centered by the coordinate-wise median and scaled by MAD. Cells are ranked by robust radial distance; the selected coverage fraction is enclosed by a convex hull, expanded about the robust center, and smoothed by the chosen number of Chaikin passes.",
        "boundary_controls": ["tsneBoundaryCoverage", "tsneBoundaryExpansion", "tsneBoundarySmoothing"],
        "total_length_y_axis": "Raw values and tests are unchanged; displayed ticks are divided by 1000 and the axis unit is ×10³ µm.",
        "layout_reference": str(REFERENCE / "data" / "comparison_default_layout.json"),
        "user_layout_imported": str(USER_LAYOUT) if USER_LAYOUT.is_file() else None,
        "layout_only_no_E_data_used": True,
    }
    (OUT / "DATA_AND_LAYOUT_AUDIT.json").write_text(json.dumps(audit, ensure_ascii=False, indent=2), encoding="utf-8")
    (OUT / "README.txt").write_text(
        "Macaque M morphology t-SNE + comparison editor\n"
        "The panel geometry and controls reproduce the reference E comparison editor.\n"
        "All plotted data are from the frozen 117-cell Macaque M1-M4 consensus cohort.\n"
        "M1=43, M2=42, M3=20, M4=12; 18 raw morphology features; Dunn-Holm pairwise tests.\n"
        "The feature set is restricted exactly to the 18 features in the frozen PCA loading table.\n"
        "Panels are ordered by the mean of descending ranks from Extra Trees MDI, absolute standardized logistic coefficients, and donor-grouped permutation importance.\n"
        "Frozen section colors are retained: M1 #1F77B4, M2 #D9A400, M3 #8C564B, M4 #E377C2.\n"
        "Feature t-SNE color uses the saved frozen Z score with a default range of -2.5 to 2.5.\n"
        "Class boundaries use robust median/MAD distance, an adjustable coverage fraction, convex hull, adjustable expansion, and adjustable Chaikin smoothing.\n"
        "The Total length Y-axis defaults to values divided by 1000 (×10³ µm); raw values and statistics are unchanged.\n"
        "The supplied Macaque_M_morphology_tsne_comparison_layout.json is used as the presentation default.\n"
        "The reference E data, E labels, E features, and E preprocessing are not used.\n",
        encoding="utf-8",
    )
    return OUT / "index.html"


if __name__ == "__main__":
    print(build())
