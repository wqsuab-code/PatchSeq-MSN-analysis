#!/usr/bin/env python
"""Focused Macaque morphology QC: MSN only, Ca + Pu + NAC only."""

from __future__ import annotations

from pathlib import Path
from zipfile import ZipFile
import json

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import skew, spearmanr


ROOT = Path(__file__).resolve().parents[1]
ZIP_PATH = Path(r"C:\Users\53461\Downloads\Macaque-PatchSeq-BG.zip")
OUT = ROOT / "outputs" / "Macaque_M_QC_MSN_CaPuNAc"
ID = "cell_label"
ROIS = ["Ca", "Pu", "NAC"]


def broad_t(group: object) -> str:
    s = "" if pd.isna(group) else str(group)
    x = s.lower()
    if "d1d2" in x or "hybrid" in x:
        return "Hybrid"
    if "d1" in x and "msn" in x:
        return "D1"
    if "d2" in x and "msn" in x:
        return "D2"
    return "Unresolved"


def robust_z(frame: pd.DataFrame) -> pd.DataFrame:
    med = frame.median(axis=0)
    scale = (frame.sub(med).abs().median(axis=0) * 1.4826).replace(0, np.nan)
    scale = scale.fillna(frame.std(axis=0, ddof=1)).replace(0, 1.0)
    return frame.sub(med).div(scale)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    with ZipFile(ZIP_PATH) as zf:
        morph = pd.read_csv(zf.open("Data/morphology_features.csv"))
        meta = pd.read_csv(zf.open("Data/cell_metadata_AllCell.csv"))

    features = [c for c in morph if c != ID]
    morph[features] = morph[features].apply(pd.to_numeric, errors="coerce")
    keep_meta = [ID, "Lib_region_of_interest_label", "Group_name", "Subclass_name", "Class_name", "Lib_donor_label"]
    dat = morph.merge(meta[keep_meta], on=ID, validate="one_to_one")
    dat["T_class_broad"] = dat["Group_name"].map(broad_t)
    dat["MSN_flag"] = dat["Group_name"].fillna("").str.contains("MSN", case=False)
    dat = dat.loc[dat["MSN_flag"] & dat["Lib_region_of_interest_label"].isin(ROIS)].copy()
    if set(dat["T_class_broad"]) - {"D1", "D2", "Hybrid"}:
        raise RuntimeError("Unexpected unresolved MSN group in Ca/Pu/NAC")

    primary21 = [
        c for c in features
        if c.startswith("basal_dendrite_") or c.startswith("soma_") or c.startswith("3_Sholl_")
    ]
    if len(primary21) != 21:
        raise RuntimeError(f"Expected 21 primary dendrite+soma features, found {len(primary21)}")

    dat["Observed_M_n"] = dat[features].notna().sum(axis=1)
    dat["Missing_M_n"] = dat[features].isna().sum(axis=1)
    dat["Complete_primary21"] = dat[primary21].notna().all(axis=1)
    dat["Complete_all39"] = dat[features].notna().all(axis=1)

    rows = []
    for feature in features:
        x = dat[feature]
        n = int(x.notna().sum())
        rows.append(
            {
                "Feature": feature,
                "Observed_n": n,
                "Missing_n": int(x.isna().sum()),
                "Missing_rate": float(x.isna().mean()),
                "Raw_skewness": float(skew(x.dropna(), bias=False)) if n >= 3 and x.nunique() > 1 else np.nan,
                "Minimum": x.min(),
                "Median": x.median(),
                "Maximum": x.max(),
            }
        )
    fq = pd.DataFrame(rows).sort_values(["Missing_rate", "Feature"], ascending=[False, True])

    corr_rows = []
    for i, a in enumerate(primary21):
        for b in primary21[i + 1 :]:
            pair = dat[[a, b]].dropna()
            rho, p = spearmanr(pair[a], pair[b])
            corr_rows.append({"Feature_A": a, "Feature_B": b, "N": len(pair), "Spearman_rho": rho, "Abs_rho": abs(rho), "P_value": p})
    corr = pd.DataFrame(corr_rows).sort_values("Abs_rho", ascending=False)

    rz = robust_z(dat[primary21])
    dat["Extreme_feature_n_abs_robust_z_ge3p5"] = (rz.abs() >= 3.5).sum(axis=1)
    dat["Maximum_abs_robust_z"] = rz.abs().max(axis=1)
    dat["Any_extreme"] = dat["Extreme_feature_n_abs_robust_z_ge3p5"] > 0

    counts = (
        dat.groupby(["Lib_region_of_interest_label", "T_class_broad"], observed=False)
        .agg(
            Cells=(ID, "size"),
            Donors=("Lib_donor_label", "nunique"),
            Complete_primary21=("Complete_primary21", "sum"),
            Complete_all39=("Complete_all39", "sum"),
            Extreme_screen_cells=("Any_extreme", "sum"),
        )
        .reset_index()
    )
    group_counts = dat.groupby("Group_name").size().rename("Cells").reset_index().sort_values("Cells", ascending=False)

    dat.to_csv(OUT / "01_MSN_CaPuNAc_cell_QC.csv", index=False)
    fq.to_csv(OUT / "02_MSN_CaPuNAc_feature_missingness_skewness.csv", index=False)
    counts.to_csv(OUT / "03_ROI_D1D2Hybrid_counts.csv", index=False)
    group_counts.to_csv(OUT / "04_Group_name_counts.csv", index=False)
    corr.to_csv(OUT / "05_primary21_spearman_pairs.csv", index=False)
    corr.loc[corr["Abs_rho"] >= 0.80].to_csv(OUT / "06_primary21_redundant_pairs_abs_rho_ge0p80.csv", index=False)
    dat.loc[dat["Any_extreme"], [ID, "Lib_region_of_interest_label", "T_class_broad", "Group_name", "Lib_donor_label", "Extreme_feature_n_abs_robust_z_ge3p5", "Maximum_abs_robust_z"]].sort_values(
        "Maximum_abs_robust_z", ascending=False
    ).to_csv(OUT / "07_primary21_extreme_screen_cells.csv", index=False)

    mpl.rcParams.update({"font.family": "Arial", "font.size": 7, "axes.linewidth": 0.7})
    fig, ax = plt.subplots(figsize=(4.5, 3.8))
    plot = fq.sort_values(["Missing_rate", "Feature"])
    y = np.arange(len(plot))
    ax.barh(y, plot["Missing_rate"] * 100, color="#3C5488", edgecolor="none")
    ax.set_yticks(y, plot["Feature"], fontsize=4.2)
    ax.axvline(20, ls="--", color="#F39B7F", lw=0.8)
    ax.axvline(50, ls="--", color="#B91C1C", lw=0.8)
    ax.set_xlabel("Missing cells (%)")
    ax.set_title("Macaque MSN morphology missingness: Ca + Pu + NAC")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(OUT / "08_MSN_CaPuNAc_feature_missingness.png", dpi=600, bbox_inches="tight")
    fig.savefig(OUT / "08_MSN_CaPuNAc_feature_missingness.pdf", bbox_inches="tight")
    plt.close(fig)

    cm = dat[primary21].corr(method="spearman")
    mask = np.triu(np.ones_like(cm, dtype=bool), 1)
    fig, ax = plt.subplots(figsize=(5.3, 4.8))
    im = ax.imshow(np.ma.array(cm.to_numpy(), mask=mask), cmap="RdBu_r", vmin=-1, vmax=1)
    labels = [x.replace("basal_dendrite_", "").replace("soma_", "soma ").replace("3_Sholl_", "dend Sholl ") for x in primary21]
    ax.set_xticks(np.arange(21), labels, rotation=60, ha="right", fontsize=4.2)
    ax.set_yticks(np.arange(21), labels, fontsize=4.2)
    ax.set_title("Macaque MSN primary21 Spearman correlation")
    cb = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.02)
    cb.set_label("Spearman rho")
    fig.tight_layout()
    fig.savefig(OUT / "09_MSN_CaPuNAc_primary21_spearman.png", dpi=600, bbox_inches="tight")
    fig.savefig(OUT / "09_MSN_CaPuNAc_primary21_spearman.pdf", bbox_inches="tight")
    plt.close(fig)

    skew_plot = fq.loc[fq["Feature"].isin(primary21)].copy()
    skew_plot = skew_plot.sort_values("Raw_skewness")
    fig, ax = plt.subplots(figsize=(4.6, 3.7))
    y = np.arange(len(skew_plot))
    ax.hlines(y, 0, skew_plot["Raw_skewness"], color="#BDBDBD", lw=0.7)
    ax.scatter(skew_plot["Raw_skewness"], y, color="#3C5488", s=14)
    ax.axvline(-2, color="#D55E00", ls="--", lw=0.7)
    ax.axvline(2, color="#D55E00", ls="--", lw=0.7)
    ax.set_yticks(y, [x.replace("basal_dendrite_", "").replace("soma_", "soma ").replace("3_Sholl_", "dend Sholl ") for x in skew_plot["Feature"]], fontsize=4.3)
    ax.set_xlabel("Raw skewness")
    ax.set_title("Macaque MSN primary21 skewness: Ca + Pu + NAC")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(OUT / "10_MSN_CaPuNAc_primary21_skewness.png", dpi=600, bbox_inches="tight")
    fig.savefig(OUT / "10_MSN_CaPuNAc_primary21_skewness.pdf", bbox_inches="tight")
    plt.close(fig)

    summary = {
        "scope": "Macaque MSN only; ROI Ca + Pu + NAC",
        "cells": len(dat),
        "donors": int(dat["Lib_donor_label"].nunique()),
        "D1": int((dat["T_class_broad"] == "D1").sum()),
        "D2": int((dat["T_class_broad"] == "D2").sum()),
        "Hybrid": int((dat["T_class_broad"] == "Hybrid").sum()),
        "primary21_complete_cells": int(dat["Complete_primary21"].sum()),
        "all39_complete_cells": int(dat["Complete_all39"].sum()),
        "axon_exit_theta_coronal_observed": int(dat["axon_exit_theta_coronal"].notna().sum()),
        "features_missing_gt50pct": fq.loc[fq["Missing_rate"] > 0.50, "Feature"].tolist(),
        "primary21_redundant_pairs_abs_rho_ge0p80": int((corr["Abs_rho"] >= 0.80).sum()),
        "primary21_abs_raw_skew_ge2": fq.loc[fq["Feature"].isin(primary21) & fq["Raw_skewness"].abs().ge(2), "Feature"].tolist(),
        "primary21_extreme_screen_cells": int(dat["Any_extreme"].sum()),
        "raw_reconstruction_files_in_zip": 0,
    }
    (OUT / "11_MSN_CaPuNAc_QC_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))
    print("\nCounts")
    print(counts.to_string(index=False))
    print("\nRedundant pairs")
    print(corr.loc[corr["Abs_rho"] >= 0.80].to_string(index=False))


if __name__ == "__main__":
    main()
