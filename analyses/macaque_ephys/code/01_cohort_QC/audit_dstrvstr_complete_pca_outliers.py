from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import chi2


ROOT = Path(r"C:\Users\53461\OneDrive\Documentos\Patch-seq_Mouse_Acb_MSN_T-type_Visualization")
PCA_DIR = ROOT / "outputs" / "dSTR_dSTRvSTR_E_QC" / "primary10_PCA_four_cohorts"
RAW_DIR = ROOT / "outputs" / "dSTR_dSTRvSTR_E_QC" / "primary10_sensitivity_cohorts"
OUT = PCA_DIR / "dSTR_plus_vSTR_B_complete_outlier_audit"
OUT.mkdir(parents=True, exist_ok=True)

FEATURES = [
    "fast_trough_v_rheo", "peak_v_rheo", "postap_slope_rheo", "threshold_v_rheo",
    "upstroke_downstroke_ratio_rheo", "upstroke_rheo", "width_rheo_ms",
    "avg_rate_rheo", "latency_rheo", "rheobase_i",
]

raw = pd.read_csv(RAW_DIR / "dSTR_plus_vSTR_B_primary10_complete_case.csv", low_memory=False)
zmat = pd.read_csv(PCA_DIR / "dSTR_plus_vSTR_B_complete_YeoJohnson_Zscore_matrix.csv")
scores = pd.read_csv(PCA_DIR / "dSTR_plus_vSTR_B_complete_PCA_scores.csv")
ev = pd.read_csv(PCA_DIR / "dSTR_plus_vSTR_B_complete_PCA_explained_variance.csv")
assert raw.cell_label.tolist() == zmat.cell_label.tolist() == scores.cell_label.tolist()

# Feature-wise standard z score after Yeo-Johnson transformation.  This is the
# actual matrix supplied to PCA.  A MAD score is unsuitable for avg_rate here
# because its discrete distribution has MAD=0.
Z = zmat[FEATURES].to_numpy(float)
feature_flag = np.abs(Z) > 3.5

# Tukey extreme fences on raw values (3 IQR; intentionally stricter than ordinary 1.5 IQR).
X = raw[FEATURES].to_numpy(float)
q1, q3 = np.quantile(X, [.25, .75], axis=0)
iqr = q3 - q1
raw_extreme = (X < (q1 - 3 * iqr)) | (X > (q3 + 3 * iqr))

# PCA-space squared Mahalanobis distance; PCA covariance is diagonal.
S = scores[[f"PC{i}" for i in range(1, 11)]].to_numpy(float)
eig = ev.eigenvalue.to_numpy(float)
md2_3 = np.sum(S[:, :3] ** 2 / eig[:3], axis=1)
md2_10 = np.sum(S ** 2 / eig, axis=1)
thr3_99, thr3_999 = chi2.ppf([.99, .999], 3)
thr10_99, thr10_999 = chi2.ppf([.99, .999], 10)

audit = raw[["cell_label", "donor_label", "Lib_region_of_interest_label", "Subclass_name", "T_class"]].copy()
audit["n_transformed_abs_z_gt3p5"] = feature_flag.sum(axis=1)
audit["n_raw_extreme_3IQR"] = raw_extreme.sum(axis=1)
audit["max_abs_transformed_z"] = np.max(np.abs(Z), axis=1)
audit["PC_MD2_first3"] = md2_3
audit["PC_MD2_all10"] = md2_10
audit["PC3D_above_chi2_99"] = md2_3 > thr3_99
audit["PC3D_above_chi2_99p9"] = md2_3 > thr3_999
audit["PC10D_above_chi2_99"] = md2_10 > thr10_99
audit["PC10D_above_chi2_99p9"] = md2_10 > thr10_999
audit["any_extreme_flag"] = feature_flag.any(axis=1) | raw_extreme.any(axis=1) | (md2_10 > thr10_999)
audit = audit.sort_values(["PC_MD2_all10", "max_abs_transformed_z"], ascending=False)
audit.to_csv(OUT / "cell_level_PCA_outlier_audit.csv", index=False)

detail = []
for i, row in raw.iterrows():
    for j, feature in enumerate(FEATURES):
        if feature_flag[i, j] or raw_extreme[i, j]:
            detail.append({
                "cell_label": row.cell_label, "ROI": row.Lib_region_of_interest_label,
                "T_class": row.T_class, "feature": feature, "raw_value": X[i, j],
                "transformed_z": Z[i, j], "transformed_abs_z_gt3p5": bool(feature_flag[i, j]),
                "raw_beyond_3IQR": bool(raw_extreme[i, j]),
            })
pd.DataFrame(detail).sort_values(["cell_label", "feature"]).to_csv(OUT / "flagged_cell_feature_details.csv", index=False)

summary = []
for j, feature in enumerate(FEATURES):
    summary.append({
        "feature": feature, "raw_min": X[:, j].min(), "raw_q1": q1[j], "raw_median": np.median(X[:, j]),
        "raw_q3": q3[j], "raw_max": X[:, j].max(), "n_raw_beyond_3IQR": int(raw_extreme[:, j].sum()),
        "n_abs_z_gt3p5_after_YJ": int(feature_flag[:, j].sum()),
        "max_abs_z_after_YJ": float(np.max(np.abs(Z[:, j]))),
    })
pd.DataFrame(summary).to_csv(OUT / "feature_level_extreme_summary.csv", index=False)

pd.DataFrame([{
    "n_cells": len(raw), "n_any_feature_abs_z_gt3p5": int(feature_flag.any(axis=1).sum()),
    "n_any_raw_beyond_3IQR": int(raw_extreme.any(axis=1).sum()),
    "n_PC3D_above_chi2_99": int((md2_3 > thr3_99).sum()),
    "n_PC3D_above_chi2_99p9": int((md2_3 > thr3_999).sum()),
    "n_PC10D_above_chi2_99": int((md2_10 > thr10_99).sum()),
    "n_PC10D_above_chi2_99p9": int((md2_10 > thr10_999).sum()),
    "chi2_3df_99": thr3_99, "chi2_3df_99p9": thr3_999,
    "chi2_10df_99": thr10_99, "chi2_10df_99p9": thr10_999,
}]).to_csv(OUT / "PCA_outlier_audit_summary.csv", index=False)

fig, axes = plt.subplots(1, 2, figsize=(5.0, 2.25), dpi=300)
flag = md2_10 > thr10_999
axes[0].scatter(S[~flag, 0], S[~flag, 1], s=7, c="#A9A9A9", alpha=.55, edgecolors="none")
axes[0].scatter(S[flag, 0], S[flag, 1], s=13, c="#D62728", alpha=.9, edgecolors="none")
for i in np.where(flag)[0]:
    axes[0].annotate(str(raw.cell_label.iloc[i]), (S[i, 0], S[i, 1]), fontsize=3.5, xytext=(2, 2), textcoords="offset points")
axes[0].set_xlabel("PC1", fontsize=5); axes[0].set_ylabel("PC2", fontsize=5)
order = np.argsort(md2_10)[::-1]
axes[1].scatter(np.arange(len(raw)), md2_10[order], s=6, c=np.where(flag[order], "#D62728", "#777777"), alpha=.7, edgecolors="none")
axes[1].axhline(thr10_99, color="#E69F00", lw=.7, ls="--", label="chi-square 99%")
axes[1].axhline(thr10_999, color="#D62728", lw=.7, ls="--", label="chi-square 99.9%")
axes[1].set_xlabel("Cells ranked by PCA distance", fontsize=5); axes[1].set_ylabel("Squared PCA distance (10 PCs)", fontsize=5)
for ax in axes:
    ax.tick_params(labelsize=4, length=2)
axes[1].legend(frameon=False, fontsize=4)
fig.suptitle("dSTR+vSTR complete-case PCA outlier audit (n=462)", x=.02, ha="left", fontsize=7, fontweight="bold")
fig.subplots_adjust(left=.10, right=.98, bottom=.18, top=.86, wspace=.35)
fig.savefig(OUT / "PCA_outlier_diagnostic.png", dpi=600, bbox_inches="tight", facecolor="white")
fig.savefig(OUT / "PCA_outlier_diagnostic.pdf", bbox_inches="tight", facecolor="white")
plt.close(fig)

print(pd.read_csv(OUT / "PCA_outlier_audit_summary.csv").to_string(index=False))
print("\nTop PCA-distance cells:")
print(audit.head(15).to_string(index=False))
print("\nFeature summary:")
print(pd.read_csv(OUT / "feature_level_extreme_summary.csv").to_string(index=False))
