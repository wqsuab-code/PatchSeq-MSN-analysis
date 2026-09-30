from __future__ import annotations

import csv
import json
import re
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "sites" / "ephys-core18-explorer" / "dist"
ANALYSIS = ROOT / "outputs" / "morph_qc" / "Mouse_M4_Macaque_style_ML500"
TARGET = ROOT / "interactive" / "mouse-morph-ml-studio" / "dist"
LAYOUT_IN = Path(r"C:\Users\53461\Downloads\Mouse_E_ML_validation_layout (1).json")

M_COLORS = {
    "M1": "#00468B",
    "M2": "#42B540",
    "M3": "#ED0000",
    "M4": "#0099B4",
    "neutral": "#BDBDBD",
    "accent": "#E76F00",
    "heatLow": "#F7FBFF",
    "heatHigh": "#08306B",
    "error": "#111111",
}


def rows(path: Path):
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def f(value):
    return float(value) if value not in (None, "") else None


def clean_feature(name: str) -> str:
    names = {
        "M_soma_circularity_index": "Soma circularity",
        "M_soma_aspect_ratio": "Soma aspect ratio",
        "M_cell_max_radial_dist": "Max radial distance",
        "M_total_number_of_neurites": "Primary neurite number",
        "M_basal_dendrite_avg_tortuosity": "Mean tortuosity",
        "M_Total_neurite_length_(sections)": "Total neurite length",
        "M_Number_of_bifurcation_points": "Bifurcation points",
        "M_Maximum_branch_order": "Maximum branch order",
        "M_trunk_angle_min": "Minimum trunk angle",
        "M_trunk_angle_max": "Maximum trunk angle",
    }
    return names.get(name, name.removeprefix("M_").replace("_", " "))


def build_data():
    raw = ANALYSIS / "raw"
    tab = ANALYSIS / "tables"
    audit = rows(tab / "01_frozen_input_audit_187_cells.csv")
    classes = Counter(r["M_class"] for r in audit)

    model_map = {
        "Extra trees": "ET",
        "kNN": "kNN",
        "RBF SVM": "RBF",
        "Linear SVM": "Linear",
        "Multinomial logistic": "Logit",
        "Random forest": "RF",
        "HistGradientBoosting": "HGB",
    }
    model_folds = []
    for r in rows(raw / "14_supervised_all_models_500.csv"):
        model_folds.append({
            "fold": int(r["repeat"]) + 1,
            "model": model_map[r["algorithm"]],
            "balanced_accuracy": f(r["balanced_accuracy"]),
            "macro_F1": f(r["macro_F1"]),
        })

    confusion_rows = rows(tab / "19_confusion_sum_RBF_SVM.csv")
    class_names = ["M1", "M2", "M3", "M4"]
    counts = [[int(r[c]) for c in class_names] for r in confusion_rows]
    fractions = [[v / sum(line) for v in line] for line in counts]

    roc_rows, pr_rows = rows(tab / "26_ROC_coordinates.csv"), rows(tab / "27_PR_coordinates.csv")
    roc = {c: [] for c in class_names}
    pr = {c: [] for c in class_names}
    for r in roc_rows:
        roc[r["M_class"]].append([f(r["FPR"]), f(r["TPR"])])
    for r in pr_rows:
        pr[r["M_class"]].append([f(r["recall"]), f(r["precision"])])
    auc = [{"class": r["M_class"], "ROC_AUC": f(r["ROC_AUC"]), "PR_AP": f(r["average_precision"])}
           for r in rows(tab / "28_ROC_PR_metrics.csv")]

    npc_k = []
    for r in rows(tab / "08_unsupervised_NPC_K_algorithm_summary_95CI.csv"):
        if r["algorithm"] == "Ward":
            npc_k.append({"NPC": int(r["NPC"]), "K": int(r["K"]), "ARI": f(r["ARI_vs_frozen_M4_median"])})

    alg_map = {"Gaussian mixture": "GMM", "K-means": "K-means", "Ward": "Ward", "Spectral": "Spectral"}
    discovery = []
    for r in rows(raw / "07_unsupervised_all_runs_500.csv"):
        if int(r["NPC"]) == 3 and int(r["K"]) == 4:
            discovery.append({"algorithm": alg_map[r["algorithm"]], "ARI": f(r["ARI_vs_frozen_M4"])})
    full_agreement = [{"algorithm": alg_map[r["algorithm"]], "ARI": f(r["ARI"])}
                      for r in rows(tab / "10_full_data_cross_algorithm_frozen_parameter.csv")]
    margins = [{"class": r["M_class"], "margin": f(r["consensus_margin"])}
               for r in rows(tab / "12_cell_level_consensus_margin.csv")]
    associations = [{"metadata": r["variable"].replace("sequencing_submission_batch", "Sequencing batch").replace("T_class", "T class"),
                     "cramers_v": f(r["Cramers_V"]), "empirical_P": f(r["permutation_p"])}
                    for r in rows(tab / "31_metadata_association_summary.csv") if r["Cramers_V"]]

    importance = [{"feature": r["feature"], "name": clean_feature(r["feature"]), "importance": f(r["importance_mean"])}
                  for r in rows(raw / "19_nested_grouped_permutation_importance_all_splits.csv")]
    mdi = [{"feature": r["feature"], "name": clean_feature(r["feature"]), "MDI": f(r["ExtraTrees_MDI"])}
           for r in rows(raw / "18_feature_importance_MDI_logistic_all_splits.csv")]

    mdi_logit_rows = rows(tab / "20_feature_importance_MDI_logistic_500.csv")
    mdi_summary = {r["feature"]: f(r["ExtraTrees_MDI_median"]) for r in mdi_logit_rows}
    logit_summary = {r["feature"]: f(r["Logistic_abs_standardized_coef_median"]) for r in mdi_logit_rows}
    mdi_rank = {name: i + 1 for i, (name, _) in enumerate(sorted(mdi_summary.items(), key=lambda q: q[1], reverse=True))}
    logit_rank = {name: i + 1 for i, (name, _) in enumerate(sorted(logit_summary.items(), key=lambda q: q[1], reverse=True))}
    rank_agreement = [{"feature": clean_feature(name), "ET_rank": mdi_rank[name], "Logit_rank": logit_rank[name]}
                      for name in mdi_rank]

    subset_map = {
        "All 10 features": "All 10",
        "Removing top 2": "Remove top 2",
        "Removing top 5": "Remove top 5",
        "Retaining bottom 3": "Bottom 3",
    }
    ablation = [{"subset": subset_map.get(r["subset"], r["subset"]), "balanced_accuracy": f(r["balanced_accuracy"])}
                for r in rows(raw / "22_progressive_ablation_500.csv")]
    permutation = [f(r["balanced_accuracy"]) for r in rows(raw / "29_within_donor_label_permutation_500.csv")]
    perm_summary = rows(tab / "30_permutation_test_summary.csv")[0]
    rho = rows(tab / "24_feature_rank_correlation.csv")[0]

    experimental_days = len(set(r["recording_date"] for r in audit))
    return {
        "meta": {
            "n": len(audit), "features": 10, "classes": dict(classes),
            "experimentalDays": experimental_days, "sequencingBatches": 7,
            "splits": 500, "bestModel": "RBF SVM",
            "permutationP": f(perm_summary["empirical_p"]),
            "observedPermutationScore": f(perm_summary["observed_BA_median"]),
            "rankRho": f(rho["Spearman_rho"]), "seed": 20260916,
            "methods": "500 recording-day-grouped held-out splits; same-day cells kept together. Within each training split: feature-minimum shift, column-sum normalization to 10,000, log1p, z-scoring and PCA (3 components), with fitted parameters applied to held-out cells. M1–M4 colors follow the frozen morphology taxonomy.",
        },
        "colors": M_COLORS,
        "features": sorted(mdi_summary),
        "modelFolds": model_folds,
        "confusion": {"counts": counts, "fractions": fractions},
        "roc": roc, "pr": pr, "auc": auc,
        "permutation": permutation, "npcK": npc_k, "discovery": discovery,
        "fullAgreement": full_agreement, "margins": margins,
        "associations": associations, "importance": importance, "mdi": mdi,
        "rankAgreement": rank_agreement, "ablation": ablation,
    }


def adapt_layout():
    layout = json.loads(LAYOUT_IN.read_text(encoding="utf-8-sig"))
    # Retain the E_ML file exactly as the visual template. Only class colors
    # are exchanged because the Mouse morphology result has four classes.
    non_class_colors = {k: v for k, v in layout["colors"].items() if not re.fullmatch(r"E[1-5]", k)}
    layout["colors"] = {
        "M1": M_COLORS["M1"],
        "M2": M_COLORS["M2"],
        "M3": M_COLORS["M3"],
        "M4": M_COLORS["M4"],
        **non_class_colors,
    }
    layout["panels"]["h"]["yLabel"] = "ARI vs frozen M1–M4"
    return layout


def adapt_html(text: str) -> str:
    text = text.replace("Mouse E-class ML validation editor", "Mouse M-class ML validation editor")
    text = text.replace("Mouse E-class robustness editor", "Mouse M1–M4 robustness editor")
    text = text.replace("450 HC–GC consensus cells", "187 frozen morphology cells")
    text = text.replace("experimental-day-grouped held-out analysis", "recording-day-grouped held-out analysis")
    text = re.sub(r'<div class="header-actions">.*?</div>', '<div class="header-actions"><a class="button" href="../../mouse-morph-feature-explorer/dist/index.html">Morphology feature editor</a></div>', text)
    text = text.replace('<script type="module" src="ml-validation-editor.js?v=mouse-e-transform-1"></script>',
                        '<script src="data.js?v=mouse-m-eformat-2"></script><script type="module" src="ml-validation-editor.js?v=mouse-m-eformat-2"></script>')
    return text


def adapt_js(text: str) -> str:
    text = re.sub(r"colors:\{E1:.*?error:'#111111'\}",
                  "colors:{M1:'#00468B',M2:'#42B540',M3:'#ED0000',M4:'#0099B4',neutral:'#BDBDBD',accent:'#E76F00',heatLow:'#F7FBFF',heatHigh:'#08306B',error:'#111111'}",
                  text, count=1)
    text = text.replace("ARI vs frozen E5", "ARI vs frozen M1–M4")
    text = text.replace("const CLASSES=['E1','E2','E3','E4','E5'];", "const CLASSES=['M1','M2','M3','M4'];")
    text = text.replace("const b=chartBox(p),n=5,gap=p.heatGap", "const b=chartBox(p),n=CLASSES.length,gap=p.heatGap")
    text = text.replace("axis(g,p,b,[0,5],[0,5],[...CLASSES],[...CLASSES])", "axis(g,p,b,[0,n],[0,n],[...CLASSES],[...CLASSES])")
    text = text.replace("['All 18','Remove top 2','Remove top 5','Bottom 3']", "['All 10','Remove top 2','Remove top 5','Bottom 3']")
    text = text.replace("Mouse_E_ML_validation", "Mouse_M_ML_validation")
    text = text.replace("let S=structuredClone(defaults),D=null,selected='a',drag=null;",
                        "let S=structuredClone(defaults),LOADED_DEFAULT=null,D=null,selected='a',drag=null;")
    text = text.replace("$('reset').onclick=()=>{S=structuredClone(defaults);controls();panelControls();render()}",
                        "$('reset').onclick=()=>{S=structuredClone(LOADED_DEFAULT||defaults);controls();panelControls();render()}")
    load_pattern = re.compile(r"try\{const \[data,layout\]=await Promise\.all\(\[fetch\('data/mouse_ml_validation\.json'\),fetch\('data/Mouse_M_ML_validation_default_layout\.json'\)\]\);.*?\}\s*catch\(e\)\{\$\('status'\)\.textContent=e\.message\}\s*$", re.S)
    replacement = "try{D=window.MOUSE_ML_DATA;const x=window.MOUSE_ML_LAYOUT,base=structuredClone(defaults);if(!D||!x)throw Error('embedded analysis data not found');S={...base,...x,colors:{...base.colors,...x.colors},panels:{...base.panels}};for(const id of ids)S.panels[id]={...base.panels[id],...(x.panels?.[id]||{})};LOADED_DEFAULT=structuredClone(S);controls();render()}catch(e){$('status').textContent=e.message}"
    text, n = load_pattern.subn(replacement, text)
    if n != 1:
        raise RuntimeError("Could not replace the E-platform data loader")
    return text


def main():
    TARGET.mkdir(parents=True, exist_ok=True)
    (TARGET / "data").mkdir(exist_ok=True)
    data, layout = build_data(), adapt_layout()
    data_text = json.dumps(data, ensure_ascii=False, separators=(",", ":"))
    layout_text = json.dumps(layout, ensure_ascii=False, indent=2)
    (TARGET / "data" / "mouse_ml_validation.json").write_text(data_text, encoding="utf-8")
    (TARGET / "data" / "Mouse_M_ML_validation_default_layout.json").write_text(layout_text, encoding="utf-8")
    (TARGET / "data.js").write_text(f"window.MOUSE_ML_DATA={data_text};\nwindow.MOUSE_ML_LAYOUT={json.dumps(layout, ensure_ascii=False, separators=(',', ':'))};\n", encoding="utf-8")
    html = adapt_html((SOURCE / "ml-validation-editor.html").read_text(encoding="utf-8"))
    js = adapt_js((SOURCE / "ml-validation-editor.js").read_text(encoding="utf-8"))
    css = (SOURCE / "styles.css").read_text(encoding="utf-8")
    (TARGET / "index.html").write_text(html, encoding="utf-8")
    (TARGET / "ml-validation-editor.js").write_text(js, encoding="utf-8")
    (TARGET / "styles.css").write_text(css, encoding="utf-8")
    standalone = html.replace('<link rel="stylesheet" href="styles.css">', f"<style>{css}</style>")
    standalone = standalone.replace('<script src="data.js?v=mouse-m-eformat-2"></script><script type="module" src="ml-validation-editor.js?v=mouse-m-eformat-2"></script>',
                                    f"<script>{(TARGET / 'data.js').read_text(encoding='utf-8')}</script><script type=\"module\">{js}</script>")
    (TARGET.parent / "Mouse_M1-M4_ML_validation_studio.html").write_text(standalone, encoding="utf-8")
    print(json.dumps({"target": str(TARGET), "cells": data["meta"]["n"], "days": data["meta"]["experimentalDays"], "layout": [layout["canvasW"], layout["canvasH"]], "models": len(data["modelFolds"])}, indent=2))


if __name__ == "__main__":
    main()
