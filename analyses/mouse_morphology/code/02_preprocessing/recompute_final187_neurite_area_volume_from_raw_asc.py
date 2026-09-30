#!/usr/bin/env python3
"""Recompute neurite length, area, volume and radii from the raw ASC files.

This is a diagnostic remeasurement.  It reuses the established year-aware
file-matching logic, but computes every cell with one NeuroM version and one
definition.  Frozen taxonomy inputs and outputs are never overwritten.
"""

from __future__ import annotations

import importlib.util
from pathlib import Path

import matplotlib.pyplot as plt
import neurom
import numpy as np
import pandas as pd
import seaborn as sns
from neurom import features as nf
from scipy.stats import spearmanr
from sklearn.metrics import roc_auc_score


ROOT = Path(__file__).resolve().parents[1]
ASC_ROOT = Path(r"E:\ASC files")
MATCHER_FILE = ROOT / "scripts" / "analyze_asc_reconstruction_qc_and_classic_examples.py"
FINAL_ASSIGN = (
    ROOT / "outputs" / "morph_qc" / "final_morph187_NPC3_HCK4_GCres0p50_figures"
    / "00_final_morph187_cell_assignments.csv"
)
LABEL_FILE = (
    ROOT / "outputs" / "morph_qc" / "hc_seurat_shift_sum13_global8_review8_215cells"
    / "13_morph_hc_seurat_classification_215cells.csv"
)
RAW_FILE = (
    ROOT / "outputs" / "morph_qc" / "all22_shift_sum10000_log1p_zscore_215cells"
    / "01_all22_raw_215cells.csv"
)
OUT = (
    ROOT / "outputs" / "morph_qc"
    / "raw_ASC_recomputed_area_volume_final187_neurom4p0p5"
)
ID = "MSN_unique_ID"


def load_matcher():
    spec = importlib.util.spec_from_file_location("asc_matcher", MATCHER_FILE)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    module.ASC_ROOT = ASC_ROOT
    return module


def safe_auc(y: pd.Series, x: pd.Series) -> float:
    keep = y.notna() & x.notna()
    auc = roc_auc_score(y[keep].astype(int), x[keep].astype(float))
    return float(max(auc, 1.0 - auc))


def cell_geometry(path: Path) -> dict:
    neuron = neurom.load_morphology(str(path))
    radii = []
    segment_mid_radii = []
    segment_lengths = []
    for section in neuron.sections:
        pts = np.asarray(section.points, dtype=float)
        if pts.ndim != 2 or pts.shape[1] < 4:
            continue
        radii.extend(pts[:, 3].tolist())
        if len(pts) > 1:
            lengths = np.linalg.norm(np.diff(pts[:, :3], axis=0), axis=1)
            mids = 0.5 * (pts[:-1, 3] + pts[1:, 3])
            segment_lengths.extend(lengths.tolist())
            segment_mid_radii.extend(mids.tolist())
    radii = np.asarray(radii, dtype=float)
    seg_len = np.asarray(segment_lengths, dtype=float)
    seg_rad = np.asarray(segment_mid_radii, dtype=float)
    positive = np.isfinite(radii) & (radii > 0)
    seg_ok = np.isfinite(seg_len) & np.isfinite(seg_rad) & (seg_len > 0) & (seg_rad > 0)
    length = float(sum(nf.get("total_length_per_neurite", neuron)))
    area = float(sum(nf.get("total_area_per_neurite", neuron)))
    volume = float(sum(nf.get("total_volume_per_neurite", neuron)))
    return {
        "Recalc_total_length_um": length,
        "Recalc_total_area_um2": area,
        "Recalc_total_volume_um3": volume,
        "Recalc_effective_radius_2V_over_A_um": 2.0 * volume / area if area > 0 else np.nan,
        "ASC_point_radius_n": int(positive.sum()),
        "ASC_point_radius_min_um": float(np.min(radii[positive])) if positive.any() else np.nan,
        "ASC_point_radius_q25_um": float(np.quantile(radii[positive], 0.25)) if positive.any() else np.nan,
        "ASC_point_radius_median_um": float(np.median(radii[positive])) if positive.any() else np.nan,
        "ASC_point_radius_q75_um": float(np.quantile(radii[positive], 0.75)) if positive.any() else np.nan,
        "ASC_point_radius_max_um": float(np.max(radii[positive])) if positive.any() else np.nan,
        "ASC_length_weighted_segment_radius_um": (
            float(np.average(seg_rad[seg_ok], weights=seg_len[seg_ok])) if seg_ok.any() else np.nan
        ),
        "ASC_unique_positive_radius_n_6dp": int(len(np.unique(np.round(radii[positive], 6)))) if positive.any() else 0,
    }


def batch_summary(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    features = [
        "Recalc_total_length_um",
        "Recalc_total_area_um2",
        "Recalc_total_volume_um3",
        "Recalc_area_per_length_um",
        "Recalc_volume_per_length_um2",
        "Recalc_effective_radius_2V_over_A_um",
        "ASC_point_radius_median_um",
        "ASC_length_weighted_segment_radius_um",
    ]
    for feature in features:
        for batch, group in df.groupby("Batch", observed=True):
            rest = df[df["Batch"] != batch]
            med_in = float(group[feature].median())
            med_out = float(rest[feature].median())
            rows.append({
                "Feature": feature,
                "Batch": batch,
                "Batch_n": int(group[feature].notna().sum()),
                "Batch_median": med_in,
                "Other_batches_median": med_out,
                "Median_ratio": med_in / med_out if med_out != 0 else np.nan,
                "One_vs_rest_absolute_AUC": safe_auc(df["Batch"].eq(batch), df[feature]),
            })
    return pd.DataFrame(rows)


def comparison_summary(df: pd.DataFrame) -> pd.DataFrame:
    pairs = [
        ("Length", "M_Total_neurite_length_(sections)", "Recalc_total_length_um"),
        ("Area", "M_Total_neurite_area", "Recalc_total_area_um2"),
        ("Volume", "M_Total_neurite_volume", "Recalc_total_volume_um3"),
    ]
    rows = []
    for label, old, new in pairs:
        valid = df[[old, new]].dropna()
        rho, p = spearmanr(valid[old], valid[new])
        ratio = valid[new] / valid[old].replace(0, np.nan)
        rows.append({
            "Metric": label,
            "Matched_cell_n": len(valid),
            "Spearman_rho_old_vs_recalculated": float(rho),
            "Spearman_p": float(p),
            "Median_recalculated_over_historical": float(ratio.median()),
            "Median_absolute_relative_error": float(
                np.median(np.abs(valid[new] - valid[old]) / np.maximum(np.abs(valid[old]), 1e-12))
            ),
        })
    return pd.DataFrame(rows)


def diagnostic_figure(df: pd.DataFrame):
    sns.set_theme(style="ticks", font="Arial", font_scale=0.8)
    palette = {False: "#BDBDBD", True: "#ED0000"}
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 5.6))
    ax = axes[0, 0]
    for flag, group in df.groupby("Is_Batch8"):
        ax.scatter(group["M_Total_neurite_area"], group["Recalc_total_area_um2"],
                   s=11, alpha=0.75, c=palette[bool(flag)], linewidths=0,
                   label="Batch8" if flag else "Other batches")
    lo = max(1e-6, np.nanmin([df["M_Total_neurite_area"].min(), df["Recalc_total_area_um2"].min()]))
    hi = np.nanmax([df["M_Total_neurite_area"].max(), df["Recalc_total_area_um2"].max()])
    ax.plot([lo, hi], [lo, hi], ls="--", lw=0.7, color="black")
    ax.set(xscale="log", yscale="log", xlabel="Historical area (µm²)",
           ylabel="Recalculated area (µm²)", title="Historical versus ASC-recalculated area")
    ax.legend(frameon=False, fontsize=6)

    ax = axes[0, 1]
    sns.boxplot(data=df, x="Batch", y="ASC_point_radius_median_um", color="white",
                fliersize=0, linewidth=0.8, ax=ax)
    sns.stripplot(data=df, x="Batch", y="ASC_point_radius_median_um", hue="Is_Batch8",
                  palette=palette, size=2.2, alpha=0.75, linewidth=0, legend=False, ax=ax)
    ax.set(xlabel="", ylabel="Median ASC point radius (µm)", title="Raw radius encoded in ASC")
    ax.tick_params(axis="x", rotation=45)

    ax = axes[1, 0]
    sns.boxplot(data=df, x="Batch", y="Recalc_area_per_length_um", color="white",
                fliersize=0, linewidth=0.8, ax=ax)
    sns.stripplot(data=df, x="Batch", y="Recalc_area_per_length_um", hue="Is_Batch8",
                  palette=palette, size=2.2, alpha=0.75, linewidth=0, legend=False, ax=ax)
    ax.set(xlabel="", ylabel="Recalculated area / length (µm)", title="Surface scale after uniform recomputation")
    ax.tick_params(axis="x", rotation=45)

    ax = axes[1, 1]
    sns.boxplot(data=df, x="Batch", y="Recalc_effective_radius_2V_over_A_um", color="white",
                fliersize=0, linewidth=0.8, ax=ax)
    sns.stripplot(data=df, x="Batch", y="Recalc_effective_radius_2V_over_A_um", hue="Is_Batch8",
                  palette=palette, size=2.2, alpha=0.75, linewidth=0, legend=False, ax=ax)
    ax.set(xlabel="", ylabel="Effective radius, 2V/A (µm)", title="Area-volume implied radius")
    ax.tick_params(axis="x", rotation=45)
    sns.despine(fig)
    fig.tight_layout()
    fig.savefig(OUT / "06_raw_ASC_radius_area_volume_diagnostic.png", dpi=600, facecolor="white")
    fig.savefig(OUT / "06_raw_ASC_radius_area_volume_diagnostic.pdf", facecolor="white")
    plt.close(fig)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    if not ASC_ROOT.exists():
        raise FileNotFoundError(ASC_ROOT)
    matcher = load_matcher()
    final = pd.read_csv(FINAL_ASSIGN)[[ID, "M_class"]]
    labels = pd.read_csv(LABEL_FILE)[[ID, "M_Neuron_id", "T_Batch"]]
    raw = pd.read_csv(RAW_FILE)
    cells = final.merge(labels, on=ID, validate="one_to_one").merge(raw, on=ID, validate="one_to_one")

    inventory = matcher.build_inventory()
    mapping, _ = matcher.match_cells(cells, inventory)
    mapping.to_csv(OUT / "01_year_aware_ASC_mapping_final187.csv", index=False)

    geometry_rows = []
    for row in mapping.itertuples(index=False):
        record = {ID: getattr(row, ID), "Recalc_success": False, "Recalc_error": ""}
        path_text = getattr(row, "Selected_ASC_path")
        if not path_text:
            record["Recalc_error"] = "Missing_matching_ASC"
        else:
            try:
                record.update(cell_geometry(Path(path_text)))
                record["Recalc_success"] = True
            except Exception as exc:
                record["Recalc_error"] = f"{type(exc).__name__}: {exc}"
        geometry_rows.append(record)
    geometry = pd.DataFrame(geometry_rows)

    keep_map = [
        ID, "Selected_ASC_filename", "Selected_ASC_path", "Selected_modification_time",
        "Selected_modification_year", "Mapping_status", "Length_relative_error",
        "Radial_relative_error", "ASC_parse_success", "ASC_parse_error",
    ]
    result = cells.merge(mapping[keep_map], on=ID, validate="one_to_one").merge(
        geometry, on=ID, validate="one_to_one"
    )
    result["Batch"] = result[ID].str.extract(r"(Batch\d+)$", expand=False)
    result["Is_Batch8"] = result["Batch"].eq("Batch8")
    result["Recalc_area_per_length_um"] = result["Recalc_total_area_um2"] / result["Recalc_total_length_um"]
    result["Recalc_volume_per_length_um2"] = result["Recalc_total_volume_um3"] / result["Recalc_total_length_um"]
    result["Historical_area_per_length_um"] = result["M_Total_neurite_area"] / result["M_Total_neurite_length_(sections)"]
    result["Historical_volume_per_length_um2"] = result["M_Total_neurite_volume"] / result["M_Total_neurite_length_(sections)"]
    result["Area_recalc_over_historical"] = result["Recalc_total_area_um2"] / result["M_Total_neurite_area"]
    result["Volume_recalc_over_historical"] = result["Recalc_total_volume_um3"] / result["M_Total_neurite_volume"]
    result.to_csv(OUT / "02_recomputed_ASC_geometry_final187.csv", index=False)

    comparisons = comparison_summary(result[result["Recalc_success"]])
    comparisons.to_csv(OUT / "03_historical_vs_recomputed_summary.csv", index=False)
    batches = batch_summary(result[result["Recalc_success"]])
    batches.to_csv(OUT / "04_recomputed_batch_scale_audit.csv", index=False)

    status = pd.concat([
        mapping["Mapping_status"].fillna("Missing").value_counts().rename_axis("Category").reset_index(name="Cell_n").assign(Dimension="Mapping_status"),
        result["Recalc_success"].astype(str).value_counts().rename_axis("Category").reset_index(name="Cell_n").assign(Dimension="Recalc_success"),
    ], ignore_index=True)[["Dimension", "Category", "Cell_n"]]
    status.to_csv(OUT / "05_recomputation_status_summary.csv", index=False)
    diagnostic_figure(result[result["Recalc_success"]])

    b8 = batches[(batches["Batch"] == "Batch8")].copy()
    summary_lines = [
        "Raw ASC recomputation of final Morph n=187",
        f"ASC inventory: {len(inventory)} files from {ASC_ROOT}",
        f"Year-aware matches: {mapping['Selected_ASC_path'].astype(bool).sum()}/187",
        f"Successful geometry recomputation: {result['Recalc_success'].sum()}/187",
        "NeuroM version: 4.0.5",
        "",
        comparisons.to_string(index=False),
        "",
        "Batch8 one-versus-rest recomputed scale audit:",
        b8.to_string(index=False),
    ]
    (OUT / "07_run_summary.txt").write_text("\n".join(summary_lines), encoding="utf-8")
    print("\n".join(summary_lines))


if __name__ == "__main__":
    main()
