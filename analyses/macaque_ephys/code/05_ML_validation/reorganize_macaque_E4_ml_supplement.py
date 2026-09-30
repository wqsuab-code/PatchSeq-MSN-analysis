#!/usr/bin/env python
"""Reorganize existing Macaque E4 ML results into three focused supplements."""
from pathlib import Path
import os
import shutil
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import spearmanr
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

ROOT = Path(__file__).resolve().parents[1]
ML = ROOT / "outputs/Macaque_E4_full_ML500"
ABL = ROOT / "outputs/R3_E_feature_ablation"
CORR = ROOT / "outputs/Macaque_E4_ML_submission_complete/tables"
NO_PANEL_LETTERS = os.environ.get("E4_NO_PANEL_LETTERS") == "1"
OUT = ROOT / ("outputs/Macaque_E4_ML_supplement_reorganized_corrected_Elabels_no_panel_letters"
              if NO_PANEL_LETTERS else
              "outputs/Macaque_E4_ML_supplement_reorganized_corrected_Elabels")
CLASSES = ["C1", "C2", "C3", "C4"]
CLASS_COLORS = {"C1":"#F8766D", "C2":"#7CAE00", "C3":"#00BFC4", "C4":"#C77CFF"}
DISPLAY = {c: c.replace("C", "E") for c in CLASSES}
DISPLAY_COLORS = {DISPLAY[c]: CLASS_COLORS[c] for c in CLASSES}
MODEL_ORDER = ["Extra trees", "kNN", "RBF SVM", "Linear SVM", "Multinomial logistic", "Random forest", "HistGradientBoosting"]
MODEL_SHORT = {"Extra trees":"ET", "kNN":"kNN", "RBF SVM":"RBF", "Linear SVM":"Linear",
               "Multinomial logistic":"Logit", "Random forest":"RF", "HistGradientBoosting":"HGB"}
ALG_ORDER = ["Ward", "KMeans", "GMM", "Spectral"]
PANEL_W, PANEL_H = 1.55, 1.45
HGAP, VGAP = .55, .68


def style():
    mpl.rcParams.update({"font.family":"Arial", "font.size":6, "axes.titlesize":7,
        "axes.labelsize":6, "xtick.labelsize":5.3, "ytick.labelsize":5.3,
        "legend.fontsize":5, "axes.linewidth":.55, "xtick.major.width":.5,
        "ytick.major.width":.5, "xtick.major.size":2.2, "ytick.major.size":2.2,
        "pdf.fonttype":42, "ps.fonttype":42, "savefig.facecolor":"white"})


def clean(ax):
    ax.spines[["top", "right"]].set_visible(False)


def letter(ax, s):
    if not NO_PANEL_LETTERS:
        ax.text(-.16, 1.08, s, transform=ax.transAxes, fontsize=8,
                weight="bold", va="top")


def fixed_grid(nrows, ncols, left=.65, right=.25, bottom=.55, top=.45,
               hgap=HGAP, vgap=VGAP):
    width = left + ncols * PANEL_W + (ncols - 1) * hgap + right
    height = bottom + nrows * PANEL_H + (nrows - 1) * vgap + top
    fig, axes = plt.subplots(nrows, ncols, figsize=(width, height), squeeze=False)
    fig.subplots_adjust(left=left / width, right=1 - right / width,
                        bottom=bottom / height, top=1 - top / height,
                        wspace=hgap / PANEL_W, hspace=vgap / PANEL_H)
    return fig, axes


def save(fig, stem):
    fig.savefig(OUT / f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(OUT / f"{stem}.png", dpi=600, bbox_inches="tight")
    plt.close(fig)


def feature_name(s):
    return (s.removeprefix("Epsy_").replace("upstroke_downstroke", "up/down")
            .replace("deltav", "delta V").replace("5spike", "5-spike")
            .replace("adapt_ratio", "adaptation ratio").replace("_rheo", " (rheo)")
            .replace("_", " "))


def supervised_figure(perf, summary, cm, perm, observed_ba, roc, pr, aucs):
    fig, grid = fixed_grid(2, 4, top=.36)
    axes = list(grid[0]) + list(grid[1, :3])
    grid[1, 3].axis("off")
    # Centre the three-panel second row under the four-panel first row.
    step = grid[0,1].get_position().x0 - grid[0,0].get_position().x0
    for ax in grid[1,:3]:
        p = ax.get_position()
        ax.set_position([p.x0 + step / 2, p.y0, p.width, p.height])
    p = perf.assign(model=perf.algorithm.map(MODEL_SHORT))
    order = [MODEL_SHORT[x] for x in MODEL_ORDER]
    for ax, metric, title, lab in [(axes[0], "balanced_accuracy", "Held-out balanced accuracy", "a"),
                                   (axes[1], "macro_F1", "Held-out macro-F1", "b")]:
        sns.boxplot(data=p, x="model", y=metric, order=order, color="#009E73", width=.58,
                    showfliers=False, linewidth=.55, ax=ax)
        ax.axhline(.25, color="#777777", ls=":", lw=.65)
        ax.set(xlabel="", ylabel=metric.replace("_", " "), ylim=(.2, 1.01), title=title)
        ax.tick_params(axis="x", rotation=30); clean(ax); letter(ax, lab)
    sns.heatmap(cm, cmap="Blues", vmin=0, vmax=100, annot=True, fmt=".0f",
                cbar=False, xticklabels=[DISPLAY[c] for c in CLASSES],
                yticklabels=[DISPLAY[c] for c in CLASSES], annot_kws={"fontsize":6}, ax=axes[2])
    axes[2].set(xlabel="Predicted class", ylabel="Frozen class", title="Confusion (%)")
    axes[2].tick_params(length=0); letter(axes[2], "c")

    et = perf[perf.algorithm.eq("Extra trees")]
    recall = et.melt(id_vars=["repeat"], value_vars=[f"recall_{c}" for c in CLASSES],
                     var_name="Class", value_name="Recall")
    recall["Class"] = recall.Class.str.removeprefix("recall_").map(DISPLAY)
    display_order = [DISPLAY[c] for c in CLASSES]
    sns.boxplot(data=recall, x="Class", y="Recall", order=display_order,
                palette=DISPLAY_COLORS, hue="Class", legend=False, width=.58, showfliers=False,
                linewidth=.55, ax=axes[3])
    axes[3].set(xlabel="",ylim=(.45,1.01), title="Class recall"); clean(axes[3]); letter(axes[3], "d")

    sns.histplot(perm.mean_balanced_accuracy, bins=24, stat="density", color="#BDBDBD", edgecolor="white", linewidth=.3, ax=axes[4])
    axes[4].axvline(perm.mean_balanced_accuracy.median(), color="#666666", lw=.8, label="Permuted")
    axes[4].axvline(observed_ba, color="#D55E00", lw=1, label="Observed")
    axes[4].set(xlabel="Balanced accuracy", ylabel="Density", title="Permutation control", xlim=(.15,1.02))
    axes[4].legend(frameon=False); clean(axes[4]); letter(axes[4], "e")

    metric=aucs.set_index("Class")
    for c in CLASSES:
        q=roc[roc.Class.eq(c)]; axes[5].plot(q.FPR,q.TPR,color=CLASS_COLORS[c],lw=1,label=f"{DISPLAY[c]} {metric.loc[c,'ROC_AUC']:.3f}")
        q=pr[pr.Class.eq(c)]; axes[6].plot(q.Recall,q.Precision,color=CLASS_COLORS[c],lw=1,label=f"{DISPLAY[c]} {metric.loc[c,'Average_precision']:.3f}")
    axes[5].plot([0,1],[0,1],color="#999",ls="--",lw=.6)
    axes[5].set(xlabel="False-positive rate",ylabel="True-positive rate",xlim=(0,1),ylim=(0,1.01),title="One-vs-rest ROC (AUC)")
    axes[6].set(xlabel="Recall",ylabel="Precision",xlim=(0,1),ylim=(0,1.01),title="One-vs-rest PR (AP)")
    for ax,lab in [(axes[5],"f"),(axes[6],"g")]: ax.legend(frameon=False,loc="lower right"); clean(ax); letter(ax,lab)
    fig.suptitle("Macaque E4: donor-independent supervised recoverability", y=.988, fontsize=8)
    save(fig, "Supplementary_ML1_supervised_recoverability")
    recall.to_csv(OUT/"SourceData_ML1d_class_recall.csv", index=False)


def discovery_figure(sens, uns, cross, cells, conf):
    # Extra horizontal room is reserved for the heat-map color bar and its label;
    # the data axes themselves retain the same fixed geometry as every panel.
    fig, axes = fixed_grid(2, 3, hgap=.82, top=.36)
    pv = sens.pivot(index="NPC", columns="K", values="ARI_median").loc[range(2,11),range(2,9)]
    cax = inset_axes(axes[0,0], width="4%", height="100%", loc="lower left",
                     bbox_to_anchor=(1.05, 0, 1, 1), bbox_transform=axes[0,0].transAxes,
                     borderpad=0)
    sns.heatmap(pv, cmap="viridis", vmin=0, vmax=.8, annot=True, fmt=".2f",
                annot_kws={"fontsize":4.1}, cbar_kws={},
                cbar_ax=cax, ax=axes[0,0])
    cax.set_title("ARI", fontsize=5, pad=3)
    axes[0,0].add_patch(plt.Rectangle((2,1),1,1,fill=False,edgecolor="white",lw=1.2))
    axes[0,0].set(title="NPC and K sensitivity"); letter(axes[0,0], "a")

    sns.violinplot(data=uns, x="algorithm", y="ARI_vs_reference", order=ALG_ORDER,
                   color="#A6CEE3", inner=None, cut=0, linewidth=.45, ax=axes[0,1])
    sns.boxplot(data=uns, x="algorithm", y="ARI_vs_reference", order=ALG_ORDER,
                color="white", width=.23, showfliers=False, linewidth=.5, ax=axes[0,1])
    axes[0,1].set(xlabel="", ylabel="ARI vs frozen E4", ylim=(0,1.01), title="Donor-resampled discovery")
    axes[0,1].tick_params(axis="x", rotation=25); clean(axes[0,1]); letter(axes[0,1], "b")

    q = cross[cross.algorithm_B.eq("Frozen_HC")].copy(); q["algorithm"] = q.algorithm_A
    axes[0,2].bar(q.algorithm, q.ARI, color="#4C78A8")
    axes[0,2].set(ylabel="ARI vs frozen HC", ylim=(0,1.05), title="Full-data algorithm agreement")
    axes[0,2].tick_params(axis="x", rotation=25); clean(axes[0,2]); letter(axes[0,2], "c")

    sns.histplot(cells.consensus_margin, bins=25, color="#4C78A8", edgecolor="white", linewidth=.3, ax=axes[1,0])
    axes[1,0].axvline(0,color="#D55E00",ls="--",lw=.8)
    axes[1,0].set(xlabel="Within-class minus max other-class\nco-clustering probability", ylabel="Cells", title="Cell-level consensus margin")
    clean(axes[1,0]); letter(axes[1,0], "d")

    sns.violinplot(data=cells, x="reference_class", y="consensus_margin", order=CLASSES,
                   palette=CLASS_COLORS, hue="reference_class", legend=False, inner="box", cut=0, linewidth=.45, ax=axes[1,1])
    axes[1,1].axhline(0,color="#D55E00",ls="--",lw=.7)
    axes[1,1].set(xlabel="Frozen class", ylabel="Consensus margin", title="Assignment stability by class")
    axes[1,1].set_xticks(range(4), [DISPLAY[c] for c in CLASSES])
    clean(axes[1,1]); letter(axes[1,1], "e")

    names={"Lib_region_of_interest_label":"ROI", "T_class":"T class", "donor_label":"Donor"}
    q=conf.assign(label=conf.metadata.map(names)).sort_values("cramers_v_bias_corrected")
    axes[1,2].barh(q.label,q.cramers_v_bias_corrected,color="#8C8C8C")
    for i,row in enumerate(q.itertuples()):
        axes[1,2].text(.425,i,f"P={row.empirical_p:.3g}",ha="right",va="center",fontsize=4.5)
    axes[1,2].set(xlabel="Bias-corrected Cramer's V", xlim=(0,.44), title="Association with metadata")
    clean(axes[1,2]); letter(axes[1,2], "f")
    fig.suptitle("Macaque E4: unsupervised discovery and assignment stability", y=.988, fontsize=8)
    save(fig, "Supplementary_ML2_discovery_stability")


def feature_figure(fi, perm_imp, abl):
    fig, axes = fixed_grid(2, 2, left=1.55, hgap=.78, top=.36)
    q = perm_imp.sort_values("mean").copy(); q["name"] = q.feature.map(feature_name)
    axes[0,0].barh(q.name,q["mean"],xerr=q["std"],color="#4C78A8",height=.65,error_kw={"lw":.45,"capsize":1})
    axes[0,0].set(xlabel="Held-out Δ balanced accuracy", title="Train-only permutation importance")
    axes[0,0].tick_params(axis="y",labelsize=4.1); clean(axes[0,0]); letter(axes[0,0],"a")

    top=fi.sort_values("ExtraTrees_median").copy(); top["name"]=top.feature.map(feature_name)
    y=np.arange(len(top)); med=top.ExtraTrees_median.to_numpy()
    axes[0,1].barh(y,med,color="#56B4E9",height=.65)
    axes[0,1].errorbar(med,y,xerr=[med-top.ExtraTrees_q025,top.ExtraTrees_q975-med],fmt="none",ecolor="#333",lw=.45,capsize=1)
    axes[0,1].set_yticks(y,top.name,fontsize=4.1); axes[0,1].set(xlabel="Extra Trees MDI",title="MDI stability (500 held-donor splits)")
    clean(axes[0,1]); letter(axes[0,1],"b")

    m=fi[["feature","ExtraTrees_median","Logistic_abscoef_median"]].copy()
    m["ET_rank"]=m.ExtraTrees_median.rank(ascending=False); m["LR_rank"]=m.Logistic_abscoef_median.rank(ascending=False)
    rho,p=spearmanr(m.ET_rank,m.LR_rank)
    axes[1,0].scatter(m.ET_rank,m.LR_rank,s=12,color="#0072B2",edgecolor="white",linewidth=.25)
    axes[1,0].plot([1,19],[1,19],color="#999",ls="--",lw=.6)
    axes[1,0].set(xlabel="Extra Trees rank",ylabel="Logistic coefficient rank",xlim=(.3,19.7),ylim=(19.7,.3),title=f"Cross-model rank agreement (rho={rho:.2f})")
    axes[1,0].invert_xaxis(); clean(axes[1,0]); letter(axes[1,0],"c")

    order=["All 19","Remove top 2","Remove top 5","Bottom 3"]
    sns.boxplot(data=abl,x="subset",y="balanced_accuracy",order=order,color="#009E73",width=.55,showfliers=False,linewidth=.55,ax=axes[1,1])
    sns.stripplot(data=abl,x="subset",y="balanced_accuracy",order=order,color="#333",size=1.3,alpha=.45,jitter=.18,ax=axes[1,1])
    axes[1,1].axhline(.25,color="#777",ls=":",lw=.65)
    axes[1,1].set(xlabel="",ylabel="Balanced accuracy",ylim=(.2,1.01),title="Progressive feature ablation")
    axes[1,1].set_xticklabels(["All 19","− top 2","− top 5","Bottom 3"],rotation=22,ha="right")
    clean(axes[1,1]); letter(axes[1,1],"d")
    fig.suptitle("Macaque E4: feature-level interpretability and perturbation robustness", y=.988, fontsize=8)
    save(fig,"Supplementary_ML3_feature_robustness")
    m.to_csv(OUT/"SourceData_ML3c_cross_model_feature_ranks.csv",index=False)


def main():
    OUT.mkdir(parents=True, exist_ok=True); style()
    perf=pd.read_csv(ML/"04_supervised_500_all_models_class_recall.csv")
    summary=pd.read_csv(ML/"05_supervised_summary_95CI.csv")
    cm=pd.read_csv(ML/"06_confusion_sum_Extra_trees.csv",index_col=0).to_numpy(float); cm=cm/cm.sum(1,keepdims=True)*100
    perm=pd.read_csv(CORR/"01_matched_ExtraTrees_within_donor_permutation_500.csv")
    observed_ba=pd.read_csv(CORR/"02_matched_ExtraTrees_permutation_summary.csv").observed_mean_balanced_accuracy.iloc[0]
    roc=pd.read_csv(ML/"16_ExtraTrees_ROC_curve.csv"); pr=pd.read_csv(ML/"17_ExtraTrees_PR_curve.csv"); aucs=pd.read_csv(ML/"18_ExtraTrees_ROC_PR_metrics.csv")
    sens=pd.read_csv(ML/"03_NPC_K_sensitivity_summary.csv")
    uns=pd.read_csv(ML/"08_cross_algorithm_E4_donor_resampling_500.csv")
    cross=pd.read_csv(ML/"09_full_data_cross_algorithm_ARI.csv")
    cells=pd.read_csv(ML/"11_cell_level_consensus_stability.csv")
    conf=pd.read_csv(CORR/"08_metadata_bias_corrected_permutation_tests.csv")
    fi=pd.read_csv(ML/"07_feature_importance_stability_500.csv")
    perm_imp=pd.read_csv(CORR/"04_train_only_permutation_importance_summary.csv")
    abl=pd.read_csv(CORR/"05_nested_ablation_performance.csv")
    supervised_figure(perf,summary,cm,perm,observed_ba,roc,pr,aucs)
    discovery_figure(sens,uns,cross,cells,conf)
    feature_figure(fi,perm_imp,abl)
    for p in [ML/"05_supervised_summary_95CI.csv", ML/"03_NPC_K_sensitivity_summary.csv",
              ML/"11_cell_level_consensus_stability.csv",
              ML/"15_ExtraTrees_cell_mean_heldout_probabilities.csv", ML/"16_ExtraTrees_ROC_curve.csv",
              ML/"17_ExtraTrees_PR_curve.csv", ML/"18_ExtraTrees_ROC_PR_metrics.csv",
              CORR/"01_matched_ExtraTrees_within_donor_permutation_500.csv",
              CORR/"02_matched_ExtraTrees_permutation_summary.csv",
              CORR/"04_train_only_permutation_importance_summary.csv",
              CORR/"05_nested_ablation_performance.csv",
              CORR/"08_metadata_bias_corrected_permutation_tests.csv"]:
        shutil.copy2(p,OUT/p.name)
    readme="""Macaque E4 machine-learning supplementary material — reorganized
================================================================
Scope: MSN cells from Ca + Pu + NAc; 390 complete E19 cells; frozen HC–GC consensus n=368 (C1=59, C2=133, C3=57, C4=119); 55 donors.

Supplementary ML1 — supervised recoverability
a,b: balanced accuracy and macro-F1 across 500 donor-held-out 75/25 group splits and seven model families.
c: row-normalized pooled held-out Extra Trees confusion matrix.
d: Extra Trees class-specific recall across the same 500 splits.
e: observed Extra Trees balanced accuracy versus 500 within-donor label-permutation controls using the matched Extra Trees pipeline.
f,g: one-vs-rest ROC and precision-recall curves from each cell's mean probability across all occasions on which that cell was held out. Each of the 368 cells contributes once to each curve.

Supplementary ML2 — discovery and assignment stability
a: 500 donor-resampled Ward analyses across NPC=2–10 and K=2–8; the white box marks frozen NPC=3, K=4.
b: donor-resampled agreement with frozen E4 for Ward, K-means, GMM and spectral clustering.
c: full-data agreement of alternative algorithms with frozen Ward labels.
d,e: cell-level co-clustering consensus margin overall and by frozen class.
f: bias-corrected Cramer's V with structured permutation P values for donor, ROI and T class; this is a confounding/biological-association audit, not a classifier score.

Supplementary ML3 — feature robustness
a: permutation importance calculated strictly within the held-out portions of 5 repeats x 5 donor-grouped folds.
b: Extra Trees MDI stability across 500 donor-held-out splits.
c: rank agreement between Extra Trees MDI and multinomial-logistic absolute coefficients.
d: nested donor-grouped feature ablation with feature ranking recomputed inside each training fold (all 19, remove top 2, remove top 5, retain bottom 3).

Interpretation boundary
These analyses demonstrate separability, reproducibility and perturbation robustness of the frozen electrophysiology-defined classes. Because labels and predictors derive from the same electrophysiological measurements, they are not independent biological validation. Orthogonal support must come from transcriptomic identity, morphology, anatomy or connectivity.

Not included
ROC/PR source probabilities were generated with the unchanged frozen preprocessing, train-only PCA and Extra Trees specification across the same 500 donor-held-out splits. Repeated test predictions were averaged per cell before curve construction, preventing repeated appearances of one cell from being treated as independent observations.
"""
    (OUT/"README_reorganized_supplement_CN.txt").write_text(readme,encoding="utf-8")
    rows=[]
    from PIL import Image
    for p in sorted(OUT.glob("Supplementary_*.png")):
        with Image.open(p) as im: rows.append({"file":p.name,"width_px":im.width,"height_px":im.height,"dpi":im.info.get("dpi",(None,None))[0]})
    pd.DataFrame(rows).to_csv(OUT/"figure_manifest.csv",index=False)
    print(f"Wrote {len(rows)} supplementary figures to {OUT}")


if __name__ == "__main__":
    main()
