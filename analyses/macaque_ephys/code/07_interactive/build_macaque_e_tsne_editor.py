from __future__ import annotations

import csv
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE_SITE = ROOT / "sites" / "ephys-core18-explorer" / "dist"
DEST = ROOT / "outputs" / "Macaque_E4_feature_explorer_site" / "dist"
DATA_DIR = DEST / "data"
ASSIGN = ROOT / "outputs" / "R3_panels_v3" / "00_tuned_p80_ee12_random_seed777_coordinates.csv"


def build_data() -> None:
    with ASSIGN.open(encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle))
    out = []
    for row in rows:
        consensus = row["Consensus"].lower() == "true"
        t_identity = row["T_class"] if consensus else ""
        out.append({
            "MSN_unique_ID": row["cell_label"],
            "tSNE_1": row["tSNE1"],
            "tSNE_2": row["tSNE2"],
            "GC_E": row["GC_merged"].replace("C", "E", 1),
            "HC_E": row["HC_class"].replace("C", "E", 1),
            "HC_GC_consensus": str(consensus),
            "HC_GC_consensus_E": row["HC_class"].replace("C", "E", 1) if consensus else "",
            "T_identity": t_identity,
            "strict_T_stable": str(consensus and t_identity in {"D1", "D2"}),
        })
    fields = ["MSN_unique_ID", "tSNE_1", "tSNE_2", "GC_E", "HC_E", "HC_GC_consensus", "HC_GC_consensus_E", "T_identity", "strict_T_stable"]
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    with (DATA_DIR / "macaque_tsne_workflow_cells.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(out)


def build_page() -> None:
    html = (SOURCE_SITE / "tsne-editor.html").read_text(encoding="utf-8")
    html = html.replace("Mouse E-classification t-SNE workflow editor", "Macaque E-classification t-SNE workflow editor")
    html = html.replace("GC classification -> HC-GC E consensus -> strict T identity -> RRR input set",
                        "GC-merged → HC–GC E consensus → stable D1 → stable D2")
    html = html.replace("Fully editable GC classification, HC-GC consensus, strict T identity and RRR-input t-SNE panels.",
                        "Fully editable Macaque GC-merged, HC–GC consensus, stable D1 and stable D2 t-SNE panels.")
    html = re.sub(r'<div class="header-actions">.*?</div>',
                  '<div class="header-actions"><a class="button" href="comparison-editor.html">Feature comparisons</a><a class="button" href="rrr-editor.html">RRR editor</a><a class="button" href="ml-editor.html">Machine-learning editor</a><a class="button" href="index.html">Main editor</a></div>', html, count=1)
    html = html.replace("3 · D1-T stable only", "3 · Stable D1 consensus")
    html = html.replace("4 · D2-T stable only", "4 · Stable D2 consensus")
    html = html.replace("Loading 493 cells...", "Loading 390 Macaque cells...")
    html = html.replace("Editable GC classification, HC-GC consensus, strict T identity and RRR input t-SNE panels",
                        "Editable Macaque GC-merged, HC-GC consensus, stable D1 and stable D2 t-SNE panels")
    html = html.replace("Each page shows one t-SNE figure. Shared drawing parameters apply to all four figures. GC-merged uses GC-derived regions; the other three use the HC-GC consensus regions.",
                        "Each page shows one Macaque t-SNE figure. GC-merged uses all 390 complete E-data cells; HC–GC consensus uses 368 cells; stable D1 and D2 panels use the 346 consensus cells with unambiguous D1/D2 identity.")
    html = html.replace('tsne-editor.js?v=single-tsne-save-20260918', 'macaque-tsne-editor.js?v=macaque-e4-workflow-20260922')
    (DEST / "tsne-editor.html").write_text(html, encoding="utf-8")

    js = (SOURCE_SITE / "tsne-editor.js").read_text(encoding="utf-8")
    js = js.replace("const STORAGE_KEY = 'mouse-e-tsne-single-figure-default-v1';", "const STORAGE_KEY = 'macaque-e4-tsne-workflow-default-v1';")
    js = js.replace("const classes = ['E1','E2','E3','E4','E5'];", "const classes = ['E1','E2','E3','E4'];")
    js = js.replace("const tClasses = ['D1','D2','Ambiguous'];", "const tClasses = ['D1','D2','Hybrid'];")
    js = js.replace(
        "{id:'d1stable', title:'3  D1-T consensus stable cells only', active:r=>r.strictT && r.T_identity==='D1', color:r=>r.HC_GC_consensus_E, region:'consensus'},",
        "{id:'d1stable', title:'3  Stable D1 consensus cells', active:r=>r.consensus && r.T_identity==='D1', color:r=>r.HC_GC_consensus_E, region:'consensus'},"
    )
    js = js.replace(
        "{id:'d2stable', title:'4  D2-T consensus stable cells only', active:r=>r.strictT && r.T_identity==='D2', color:r=>r.HC_GC_consensus_E, region:'consensus'}",
        "{id:'d2stable', title:'4  Stable D2 consensus cells', active:r=>r.consensus && r.T_identity==='D2', color:r=>r.HC_GC_consensus_E, region:'consensus'}"
    )
    js = js.replace("colors:{E1:'#4E79A7',E2:'#F28E2B',E3:'#59A14F',E4:'#2AA6B8',E5:'#B07AA1',D1:'#ef3b9a',D2:'#00d7a7',Ambiguous:'#9aa0a6'}",
                    "colors:{E1:'#F8766D',E2:'#7CAE00',E3:'#00BFC4',E4:'#C77CFF',D1:'#D95F02',D2:'#008F7A',Hybrid:'#9aa0a6'}")
    js = js.replace("xMin:-12.55,xMax:9.68,yMin:-19.01,yMax:19.38,unitRatio:1.2", "xMin:-8.25,xMax:7.95,yMin:-6.95,yMax:8.85,unitRatio:1")
    js = js.replace("unitRatio:1.2", "unitRatio:1")
    js = js.replace("Mouse_E_", "Macaque_E_")
    js = js.replace("'Mouse_E_workflow_tSNE_cell_data.csv'", "'Macaque_E_workflow_tSNE_cell_data.csv'")
    js = js.replace("'Mouse_E_workflow_tSNE_layout.json'", "'Macaque_E_workflow_tSNE_layout.json'")
    js = js.replace("all 493", "all 390")
    js = js.replace("ambiguous ${ambiguous}", "Hybrid ${ambiguous}")
    init_start = js.index("async function init()")
    init_end = js.index("init().catch", init_start)
    replacement = """async function init(){controls();const text=await fetch('data/macaque_tsne_workflow_cells.csv?v=macaque-e4-workflow-20260922').then(r=>{if(!r.ok)throw Error(r.status);return r.text()});rows=parseCsv(text).map(r=>({...r,x:+r.tSNE_1,y:+r.tSNE_2,consensus:/true/i.test(r.HC_GC_consensus),T_identity:r.T_identity,strictT:/true/i.test(r.strict_T_stable)}));if(rows.length!==390)throw Error(`Expected 390 cells, found ${rows.length}`);sync();render()}\n"""
    js = js[:init_start] + replacement + js[init_end:]
    js = js.replace("Unable to load t-SNE data", "Unable to load Macaque t-SNE data")
    (DEST / "macaque-tsne-editor.js").write_text(js, encoding="utf-8")

    comparison = (DEST / "comparison-editor.html").read_text(encoding="utf-8")
    if 'href="tsne-editor.html"' not in comparison:
        comparison = comparison.replace('<a class="button" href="rrr-editor.html">RRR editor</a>',
                                        '<a class="button" href="tsne-editor.html">t-SNE workflow</a><a class="button" href="rrr-editor.html">RRR editor</a>')
        (DEST / "comparison-editor.html").write_text(comparison, encoding="utf-8")


def main() -> None:
    build_data()
    build_page()
    print(DEST / "tsne-editor.html")


if __name__ == "__main__":
    main()
