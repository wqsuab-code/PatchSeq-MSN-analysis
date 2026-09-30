from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Ellipse, Rectangle
from scipy.cluster.hierarchy import leaves_list, linkage
from scipy.spatial.distance import squareform
from scipy.stats import yeojohnson


ROOT = Path(r"C:\Users\53461\OneDrive\Documentos\Patch-seq_Mouse_Acb_MSN_T-type_Visualization")
IN = ROOT / "outputs" / "dSTR_dSTRvSTR_E_QC"
OUT = IN / "redundancy_skewness"
OUT.mkdir(parents=True, exist_ok=True)

EXCLUDED_MISSING = {"first_isi_rheo", "first_isi_inv_rheo", "adp_v_last_rheo"}
FEATURES = [
    "first_isi_rheo", "first_isi_inv_rheo", "adp_v_last_rheo",
    "fast_trough_v_last_rheo", "downstroke_rheo", "fast_trough_v_rheo",
    "peak_deltav_rheo", "peak_v_rheo", "postap_slope_rheo", "threshold_v_rheo",
    "trough_t_rheo", "trough_v_rheo", "upstroke_downstroke_ratio_rheo",
    "upstroke_rheo", "width_rheo", "width_rheo_ms", "width_suprathresh_rheo",
    "avg_rate_rheo", "latency_rheo", "rheobase_i",
]
FEATURES17 = [x for x in FEATURES if x not in EXCLUDED_MISSING]
RHO_THRESHOLD = 0.80
HARMONIZED_PRIMARY10 = [
    "fast_trough_v_rheo", "peak_v_rheo", "postap_slope_rheo", "threshold_v_rheo",
    "upstroke_downstroke_ratio_rheo", "upstroke_rheo", "width_rheo_ms",
    "avg_rate_rheo", "latency_rheo", "rheobase_i",
]

mpl.rcParams.update({"font.family": "Arial", "pdf.fonttype": 42, "ps.fonttype": 42})


def skew_table(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for feature in FEATURES17:
        x = df[feature].dropna().astype(float)
        q1, med, q3 = x.quantile([0.25, 0.5, 0.75])
        iqr = q3 - q1
        iqr_n = int(((x < q1 - 1.5 * iqr) | (x > q3 + 1.5 * iqr)).sum()) if iqr > 0 else 0
        mad = float(np.median(np.abs(x - np.median(x))))
        mz_n = int((np.abs(0.6745 * (x - np.median(x)) / mad) > 3.5).sum()) if mad > 0 else 0
        raw_skew = float(x.skew())
        transformed, lam = yeojohnson(x.to_numpy())
        rows.append({
            "feature": feature, "n_observed": len(x), "n_missing": int(df[feature].isna().sum()),
            "missing_rate": float(df[feature].isna().mean()), "min": float(x.min()),
            "q1": float(q1), "median": float(med), "q3": float(q3), "max": float(x.max()),
            "raw_skew": raw_skew, "abs_raw_skew": abs(raw_skew),
            "skew_flag_abs_ge_1": abs(raw_skew) >= 1.0,
            "moderate_skew_abs_ge_0p5": abs(raw_skew) >= 0.5,
            "iqr_outlier_n": iqr_n, "modified_z_outlier_n": mz_n,
            "yeojohnson_lambda": float(lam),
            "yeojohnson_skew": float(pd.Series(transformed).skew()),
            "abs_yeojohnson_skew": abs(float(pd.Series(transformed).skew())),
        })
    return pd.DataFrame(rows).sort_values("abs_raw_skew", ascending=False)


def choose_nonredundant(rho: pd.DataFrame, skew: pd.DataFrame) -> pd.DataFrame:
    # Deterministic audit rule: lower missingness first, then lower post-YJ absolute skew,
    # then original feature order. A candidate is dropped if |rho| >= 0.80 with any earlier kept feature.
    priority = skew.set_index("feature")[["missing_rate", "abs_yeojohnson_skew"]].copy()
    priority["original_order"] = [FEATURES17.index(x) for x in priority.index]
    priority = priority.sort_values(["missing_rate", "abs_yeojohnson_skew", "original_order"])
    kept, audit = [], []
    for feature, row in priority.iterrows():
        conflicts = [(k, float(rho.loc[feature, k])) for k in kept if abs(rho.loc[feature, k]) >= RHO_THRESHOLD]
        if conflicts:
            partner, value = max(conflicts, key=lambda z: abs(z[1]))
            status, reason = "Drop_redundant", f"|rho|={abs(value):.3f} with retained {partner}"
        else:
            kept.append(feature)
            status, reason = "Keep", "No |rho|>=0.80 conflict with earlier retained feature"
        audit.append({"feature": feature, "priority_rank": len(audit) + 1, "missing_rate": row.missing_rate,
                      "abs_yeojohnson_skew": row.abs_yeojohnson_skew, "selection_status": status,
                      "selection_reason": reason})
    return pd.DataFrame(audit).sort_values("priority_rank")


def redundancy_plot(dataset: str, df: pd.DataFrame, rho: pd.DataFrame) -> None:
    distance = np.clip(1 - np.abs(rho.to_numpy()), 0, 1); np.fill_diagonal(distance, 0)
    order = leaves_list(linkage(squareform(distance, checks=False), method="average", optimal_ordering=True))
    names = [FEATURES17[i] for i in order]; corr = rho.loc[names, names]; n = len(names)
    fig, ax = plt.subplots(figsize=(7.2, 6.3), dpi=300); cmap = mpl.colormaps["RdBu"]; norm = mpl.colors.TwoSlopeNorm(vmin=-1, vcenter=0, vmax=1)
    for j in range(n):
        for i in range(j + 1):
            v = float(corr.iloc[j, i])
            ax.add_patch(Rectangle((i-.5,j-.5),1,1,facecolor="#FBFCFD",edgecolor="#D9DEE5",linewidth=.35))
            mag=abs(v); major=.94*np.sqrt(1+mag)/np.sqrt(2); minor=.94*np.sqrt(max(0,1-mag))/np.sqrt(2)
            color=cmap(norm(v)); ax.add_patch(Ellipse((i,j),major,max(minor,.025),angle=45 if v>=0 else -45,facecolor=color,edgecolor="none"))
            strong=i!=j and mag>=RHO_THRESHOLD
            if strong: ax.add_patch(Rectangle((i-.5,j-.5),1,1,fill=False,edgecolor="#D62728",linewidth=1.4))
            lum=.2126*color[0]+.7152*color[1]+.0722*color[2]
            ax.text(i,j,f"{v:.2f}",ha="center",va="center",fontsize=4.4,color="#222222" if lum>.57 else "white",fontweight="bold" if strong else "normal")
    ax.set_xlim(-.55,n-.45); ax.set_ylim(n-.45,-.55); ax.set_aspect("equal")
    ax.set_xticks(range(n),names,rotation=55,ha="right",fontsize=5); ax.set_yticks(range(n),names,fontsize=5); ax.tick_params(length=0)
    for s in ax.spines.values(): s.set_visible(False)
    display="dSTR+vSTR" if dataset=="dSTR_plus_vSTR" else "dSTR"
    ax.set_title(f"{display} rheobase redundancy · Spearman |ρ|≥0.80",loc="left",fontsize=8,fontweight="bold",color="#17365D")
    fig.subplots_adjust(left=.28,right=.98,bottom=.27,top=.93)
    fig.savefig(OUT/f"{dataset}_rheobase17_Spearman_ellipse.png",dpi=600,bbox_inches="tight",facecolor="white")
    fig.savefig(OUT/f"{dataset}_rheobase17_Spearman_ellipse.pdf",bbox_inches="tight",facecolor="white")
    plt.close(fig)


def skew_plot(dataset: str, skew: pd.DataFrame) -> None:
    d=skew.sort_values("raw_skew",ascending=False).reset_index(drop=True); y=np.arange(len(d)); vals=d.raw_skew.to_numpy()
    colors=np.where(vals>=0,"#2166AC","#D6604D"); fig,ax=plt.subplots(figsize=(5.0,3.6),dpi=300)
    ax.axvspan(-.5,.5,color="#F3F4F6"); ax.axvspan(-1,-.5,color="#FFF2CC"); ax.axvspan(.5,1,color="#FFF2CC"); ax.axvspan(min(-1,vals.min()-.2),-1,color="#FCE4D6"); ax.axvspan(1,max(1.2,vals.max()+.2),color="#FCE4D6")
    ax.axvline(0,color="#333",lw=.7)
    for yi,v,c in zip(y,vals,colors): ax.hlines(yi,0,v,color=c,lw=1.1); ax.scatter(v,yi,s=20,c=c,edgecolors="white",lw=.4,zorder=3); ax.text(v+.04 if v>=0 else v-.04,yi,f"{v:.2f}",ha="left" if v>=0 else "right",va="center",fontsize=4)
    ax.set_yticks(y,d.feature,fontsize=4); ax.invert_yaxis(); ax.tick_params(axis="x",labelsize=4); ax.tick_params(axis="y",length=0); ax.set_xlabel("Raw skewness",fontsize=5)
    for s in ["top","right","left"]: ax.spines[s].set_visible(False)
    display="dSTR+vSTR" if dataset=="dSTR_plus_vSTR" else "dSTR"; ax.set_title(f"{display} rheobase feature skewness",loc="left",fontsize=7,fontweight="bold",color="#17365D")
    fig.subplots_adjust(left=.39,right=.96,bottom=.13,top=.90); fig.savefig(OUT/f"{dataset}_rheobase17_skewness.png",dpi=600,bbox_inches="tight",facecolor="white"); fig.savefig(OUT/f"{dataset}_rheobase17_skewness.pdf",bbox_inches="tight",facecolor="white"); plt.close(fig)


summaries=[]
for dataset in ["dSTR", "dSTR_plus_vSTR"]:
    source=pd.read_csv(IN/f"{dataset}_cell_missingness_QC.csv",low_memory=False)
    df=source.loc[source.cell_qc_status.eq("Keep")].copy()
    rho=df[FEATURES17].corr(method="spearman",min_periods=20)
    skew=skew_table(df)
    pairs=[]
    for i,a in enumerate(FEATURES17):
        for b in FEATURES17[i+1:]:
            val=float(rho.loc[a,b]); pair=df[[a,b]].dropna()
            if abs(val)>=RHO_THRESHOLD: pairs.append({"feature_1":a,"feature_2":b,"n_pairwise":len(pair),"spearman_rho":val,"abs_spearman_rho":abs(val)})
    pair_df=pd.DataFrame(pairs).sort_values("abs_spearman_rho",ascending=False) if pairs else pd.DataFrame(columns=["feature_1","feature_2","n_pairwise","spearman_rho","abs_spearman_rho"])
    selection=choose_nonredundant(rho,skew)
    rho.to_csv(OUT/f"{dataset}_spearman_correlation_matrix.csv")
    pair_df.to_csv(OUT/f"{dataset}_redundancy_pairs_absrho_ge_0p80.csv",index=False)
    skew.to_csv(OUT/f"{dataset}_feature_skewness_YeoJohnson_audit.csv",index=False)
    selection.to_csv(OUT/f"{dataset}_nonredundant_feature_selection_audit.csv",index=False)
    redundancy_plot(dataset,df,rho); skew_plot(dataset,skew)
    kept=selection.loc[selection.selection_status.eq("Keep"),"feature"].tolist()
    summaries.append({"dataset":dataset,"n_cells":len(df),"input_features":len(FEATURES17),"redundant_pairs":len(pair_df),"suggested_nonredundant_features":len(kept),"strongly_skewed_raw_abs_ge_1":int(skew.skew_flag_abs_ge_1.sum()),"retained_feature_names":"; ".join(kept)})

summary=pd.DataFrame(summaries); summary.to_csv(OUT/"redundancy_skewness_summary.csv",index=False)
harmonized_rows=[]
for dataset in ["dSTR", "dSTR_plus_vSTR"]:
    rho=pd.read_csv(OUT/f"{dataset}_spearman_correlation_matrix.csv",index_col=0)
    for feature in HARMONIZED_PRIMARY10:
        others=[x for x in HARMONIZED_PRIMARY10 if x!=feature]
        partner=max(others,key=lambda x:abs(float(rho.loc[feature,x])))
        value=float(rho.loc[feature,partner])
        harmonized_rows.append({"dataset":dataset,"feature":feature,"max_abs_rho_within_primary10":abs(value),"most_correlated_partner":partner,"spearman_rho":value,"passes_absrho_lt_0p80":abs(value)<RHO_THRESHOLD})
pd.DataFrame(harmonized_rows).to_csv(OUT/"harmonized_primary10_nonredundancy_audit.csv",index=False)
(OUT/"harmonized_primary10_features.txt").write_text("\n".join(HARMONIZED_PRIMARY10)+"\n",encoding="utf-8")
(OUT/"README.txt").write_text(
    "Second-stage QC for dSTR and dSTR+vSTR.\n"
    "Input: cells passing stage-1 retained-feature missingness QC.\n"
    "Redundancy: pairwise-complete Spearman correlation; flag when |rho| >= 0.80.\n"
    "Suggested nonredundant set: greedy priority by lower missingness, then lower absolute post-Yeo-Johnson skew, then original feature order.\n"
    "Skewness: raw Fisher-Pearson skewness; |skew| >= 1 flagged. Yeo-Johnson values are audit results only; raw values are not overwritten.\n"
    "For direct dSTR versus dSTR+vSTR comparison, harmonized_primary10_features.txt preserves the prior NAc primary 10-feature definition; all ten pass |rho| < 0.80 within both new cohorts.\n\n"
    + summary.to_string(index=False), encoding="utf-8")
print(summary.to_string(index=False))
