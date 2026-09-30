"""Build the Macaque-M interactive explorer from the frozen Mouse-M UI template."""
from __future__ import annotations

import base64
import re
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
TEMPLATE = ROOT / "interactive" / "mouse-morph-feature-explorer" / "dist" / "index.html"
SOURCE = ROOT / "macaque_m" / "m18_tempfreeze_NPC5_HCK4_res2.3"
COMPARE = SOURCE / "M4_compare_all18"
RRR = SOURCE / "RRR_T_M"
OUT = ROOT / "interactive" / "macaque-morph-feature-explorer" / "dist"
ABBREVIATIONS = dict(pd.read_csv(ROOT / "macaque_m" / "morphology_feature_abbreviations.csv").values)

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


def encoded(text: str) -> str:
    return base64.b64encode(text.encode("utf-8")).decode("ascii")


def build_data() -> tuple[str, str, str, str]:
    raw = pd.read_csv(COMPARE / "all18_plot_data_117.csv")
    keys = [x[0] for x in FEATURES]
    out = raw[["cell_label", "M", *keys, *[f"{x}__z" for x in keys], "x", "y"]].copy()
    out = out.rename(columns={"cell_label": "MSN_unique_ID", "M": "M_class", "x": "tSNE1", "y": "tSNE2"})
    out.insert(2, "HC_GC_consensus", "True")
    assert len(out) == 117 and out.M_class.value_counts().sort_index().tolist() == [43, 42, 20, 12]

    scores = pd.read_csv(RRR / "RRR_cell_scores_n117.csv").rename(columns={
        "cell_label": "MSN_unique_ID", "T_class": "T_identity",
        "T_Component1": "T_RRR1", "T_Component2": "T_RRR2", "T_Component3": "T_RRR3",
        "M_Component1": "M_RRR1", "M_Component2": "M_RRR2", "M_Component3": "M_RRR3",
    })
    scores = scores[["MSN_unique_ID", "M_class", "T_identity", "T_RRR1", "T_RRR2", "T_RRR3", "M_RRR1", "M_RRR2", "M_RRR3"]]
    assert len(scores) == 117

    genes = pd.read_csv(RRR / "RRR_gene_correlation_loadings.csv", index_col=0).reset_index()
    genes.columns = ["Feature", "RRR1", "RRR2", "RRR3"]
    genes.insert(0, "Domain", "Transcriptomic")
    morph = pd.read_csv(RRR / "RRR_Mfeature_correlation_loadings.csv", index_col=0).reset_index()
    morph.columns = ["Feature", "RRR1", "RRR2", "RRR3"]
    morph["Feature"] = morph["Feature"].map(ABBREVIATIONS).fillna(morph["Feature"])
    morph.insert(0, "Domain", "Morphological")
    loadings = pd.concat([genes, morph], ignore_index=True)

    pairs = pd.read_csv(COMPARE / "all18_Dunn_Holm.csv")
    pairs["Comparison"] = pairs["Group_1"] + " vs " + pairs["Group_2"]
    pairs["BH_FDR_q"] = pairs["P_adj_Holm_within_feature"]
    pairs["Mann_Whitney_p"] = pairs["P_raw"]
    pairs = pairs[["Feature", "Comparison", "Mann_Whitney_p", "BH_FDR_q"]]

    OUT.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT / "data.csv", index=False)
    scores.to_csv(OUT / "rrr_scores.csv", index=False)
    loadings.to_csv(OUT / "rrr_loadings.csv", index=False)
    pairs.to_csv(OUT / "pairwise.csv", index=False)
    return tuple(x.to_csv(index=False, lineterminator="\n") for x in (out, scores, loadings, pairs))


def build_html(csvs: tuple[str, str, str, str]) -> None:
    html = TEMPLATE.read_text(encoding="utf-8")
    data, scores, loadings, pairs = csvs
    embedded = (
        "<script>\nconst decode64=s=>new TextDecoder().decode(Uint8Array.from(atob(s),c=>c.charCodeAt(0)));\n"
        "window.EMBEDDED_CSV={"
        f"data:decode64('{encoded(data)}'),rrrScores:decode64('{encoded(scores)}'),"
        f"rrrLoadings:decode64('{encoded(loadings)}'),pairwise:decode64('{encoded(pairs)}')"
        "};\n</script>"
    )
    html, n = re.subn(r"<script>\s*const decode64=.*?</script>", embedded, html, count=1, flags=re.S)
    assert n == 1

    fblock = "const F=[\n" + ",\n".join(
        f"  {{key:'{key}',label:'{label}',unit:'{unit}'}}" for key, label, unit in FEATURES
    ) + "\n];\nconst G="
    html, n = re.subn(r"const F=\[.*?\];\s*const G=", fblock, html, count=1, flags=re.S)
    assert n == 1

    replacements = {
        "Mouse M-class feature studio": "Macaque M-class feature studio",
        "Interactive Mouse MSN morphology": "Interactive Macaque MSN morphology",
        "Mouse T→M RRR": "Macaque T→M RRR",
        "Mouse T–M RRR": "Macaque T–M RRR",
        "Mouse transcriptomic–morphological": "Macaque transcriptomic–morphological",
        "Mouse M-class morphology": "Macaque M-class morphology",
        "Mouse morphology": "Macaque morphology",
        "Final 187": "Frozen consensus 117",
        "HC–GC consensus 181": "HC–GC consensus 117",
        "HC–GC consensus cells (n = 181)": "HC–GC consensus cells (n = 117)",
        "181 consensus cells": "117 consensus cells",
        "frozen 10-feature morphology response": "frozen 18-feature morphology response",
        "Frozen 10-feature morphology panel": "Frozen 18-feature morphology panel",
        "All 10": "All 18",
        "Frozen 10": "Frozen 18",
        "Mann–Whitney · BH FDR": "Dunn · Holm within feature",
        "Two-sided Mann–Whitney U": "Two-sided Dunn test",
        "BH FDR within feature": "Holm within feature",
        "Mouse_M1-M4_frozen10": "Macaque_M1-M4_frozen18",
        "Mouse_M_selected": "Macaque_M_selected",
        "Mouse_M_comparison": "Macaque_M_comparison",
        "Mouse_T-M_RRR": "Macaque_T-M_RRR",
    }
    for old, new in replacements.items():
        html = html.replace(old, new)

    # Freeze the current Macaque M palette and retain the established Z-score scale.
    color_map = {"#00468B": "#1F77B4", "#42B540": "#D9A400", "#ED0000": "#8C564B", "#0099B4": "#E377C2"}
    for old, new in color_map.items():
        html = html.replace(old, new)
    html = html.replace("columns:4", "columns:3", 1)
    html = html.replace("<option>3</option><option selected>4</option>", "<option selected>3</option><option>4</option>")

    # The feature maps must use the frozen transformed Z values rather than recomputing Z from raw values.
    html = html.replace(
        "z:(+r[key]-mu)/sigma",
        "z:Number.isFinite(+r[key+'__z'])?+r[key+'__z']:(+r[key]-mu)/sigma",
    )
    html = html.replace("pSource:'BH_FDR_q'", "pSource:'BH_FDR_q'")
    html = html.replace("Mann_Whitney_p,p.BH_FDR_q", "Mann_Whitney_p,p.BH_FDR_q")
    (OUT / "index.html").write_text(html, encoding="utf-8")

    (OUT / "README.txt").write_text(
        "Macaque M interactive explorer cloned from the final Mouse-M user interface.\n"
        "Panels: frozen 18-feature comparison/t-SNE editor and four-panel T-M RRR editor.\n"
        "Cells: 117 HC-GC consensus MSN cells; M1=43, M2=42, M3=20, M4=12.\n"
        "Feature values and transformed Z scores come from the frozen Macaque M analysis.\n"
        "t-SNE Z range is fixed by default at -2.5 to 2.5; saved __z columns are used directly.\n"
        "M colors: M1 #1F77B4; M2 #D9A400; M3 #8C564B; M4 #E377C2.\n"
        "RRR radial display scaling remains the 99th-percentile rule used by the source panel.\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    build_html(build_data())
    print(OUT / "index.html")
