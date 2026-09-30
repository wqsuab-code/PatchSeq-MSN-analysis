#!/usr/bin/env python
"""Export six Macaque E4 main panels separately, with/without text, H=1 inch."""
from pathlib import Path
import sys
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))
import plot_macaque_E4_nature_validation as base

SRC=ROOT/"outputs"/"Macaque_E4_full_ML500"
ASSIGN=ROOT/"outputs"/"R3_panels"/"00_res3_main_cell_assignments_and_tsne.csv"
COORD=ROOT/"outputs"/"R3_panels"/"04b_tuned_consensus368_p80_ee12_random_coordinates.csv"
OUT=ROOT/"outputs"/"Macaque_E4_Nature_validation"/"individual_H1in"
WIDTH={"a":1.05,"b":1.00,"c":1.12,"d":1.05,"e":1.00,"f":1.12}


def setup():
    mpl.rcParams.update({"font.family":"Arial","font.size":4,"axes.titlesize":4.8,"axes.labelsize":4.2,
        "xtick.labelsize":3.8,"ytick.labelsize":3.8,"legend.fontsize":3.5,"axes.linewidth":.45,
        "xtick.major.width":.4,"ytick.major.width":.4,"xtick.major.size":1.8,"ytick.major.size":1.8,
        "pdf.fonttype":42,"savefig.facecolor":"white"})


def strip_text(fig,ax):
    ax.set_title(""); ax.set_xlabel(""); ax.set_ylabel("")
    ax.set_xticklabels([]); ax.set_yticklabels([])
    if ax.get_legend() is not None: ax.get_legend().remove()
    for t in list(ax.texts): t.remove()
    for other in fig.axes:
        if other is ax: continue
        other.set_title(""); other.set_xlabel(""); other.set_ylabel("")
        other.set_xticklabels([]); other.set_yticklabels([])
        for t in list(other.texts): t.remove()


def save_both(letter,drawer):
    for labeled,folder in [(True,"with_text"),(False,"no_text")]:
        fig,ax=plt.subplots(figsize=(WIDTH[letter],1.0))
        fig.subplots_adjust(left=.18,right=.96,bottom=.20,top=.83)
        drawer(fig,ax,labeled)
        if not labeled:
            strip_text(fig,ax)
        target=OUT/folder; target.mkdir(parents=True,exist_ok=True)
        fig.savefig(target/f"panel_{letter}_H1in_{folder}.png",dpi=600,bbox_inches="tight",pad_inches=.01)
        plt.close(fig)


def main():
    setup()
    d=pd.read_csv(ASSIGN).drop(columns=["tSNE1","tSNE2"])
    d=d.merge(pd.read_csv(COORD,usecols=["cell_label","tSNE1","tSNE2"]),on="cell_label",validate="one_to_one")
    d=d[d.Consensus.astype(str).str.lower().eq("true")]
    uns=pd.read_csv(SRC/"08_cross_algorithm_E4_donor_resampling_500.csv")
    perf=pd.read_csv(SRC/"04_supervised_500_all_models_class_recall.csv")
    fi=pd.read_csv(SRC/"07_feature_importance_stability_500.csv")
    sens=pd.read_csv(SRC/"03_NPC_K_sensitivity_summary.csv")
    cm=pd.read_csv(SRC/"06_confusion_sum_Extra_trees.csv",index_col=0).to_numpy(float); cm=cm/cm.sum(1,keepdims=True)*100

    def draw_a(fig,ax,labeled):
        for c in base.COLORS:
            q=d[d.HC_class==c]; xy=q[["tSNE1","tSNE2"]].to_numpy(); base.ellipse80(ax,xy,base.COLORS[c])
            ax.scatter(xy[:,0],xy[:,1],s=2.6,c=base.COLORS[c],edgecolors="white",linewidths=.12,label=f"{c} n={len(q)}")
        ax.set_xticks([]); ax.set_yticks([]); base.freeze_tsne_axes(ax); base.clean(ax)
        if labeled: ax.set(xlabel="t-SNE 1",ylabel="t-SNE 2",title="Frozen E4 consensus"); ax.legend(frameon=False,loc="lower right",handletextpad=.2,borderpad=.1)
    save_both("a",draw_a)

    def draw_b(fig,ax,labeled):
        sns.violinplot(data=uns,x="algorithm",y="ARI_vs_reference",order=base.ALG_ORDER,color="#A6CEE3",inner=None,cut=0,linewidth=.35,ax=ax)
        sns.boxplot(data=uns,x="algorithm",y="ARI_vs_reference",order=base.ALG_ORDER,width=.22,showfliers=False,color="white",linewidth=.4,ax=ax)
        ax.axhline(.75,color="#D55E00",ls="--",lw=.45); ax.set_ylim(0,1.02); base.clean(ax)
        if labeled: ax.set(xlabel="",ylabel="ARI vs frozen E4",title="Donor-resampled discovery"); ax.tick_params(axis="x",rotation=25)
    save_both("b",draw_b)

    def draw_c(fig,ax,labeled):
        q=perf.copy(); q["model"]=q.algorithm.map(base.MODEL_SHORT); order=[base.MODEL_SHORT[x] for x in base.MODEL_ORDER]
        sns.boxplot(data=q,x="model",y="balanced_accuracy",order=order,color="#009E73",width=.55,showfliers=False,linewidth=.4,ax=ax)
        ax.axhline(.25,color="#777777",ls=":",lw=.45); ax.set_ylim(.2,1.02); base.clean(ax)
        if labeled: ax.set(xlabel="",ylabel="Balanced accuracy",title="Held-out donor recovery")
    save_both("c",draw_c)

    def draw_d(fig,ax,labeled):
        sns.heatmap(cm,cmap="Blues",vmin=0,vmax=100,annot=labeled,fmt=".0f",square=True,cbar=False,
                    xticklabels=["C1","C2","C3","C4"],yticklabels=["C1","C2","C3","C4"],annot_kws={"fontsize":3.8},ax=ax)
        ax.tick_params(length=0)
        if labeled: ax.set(xlabel="Predicted class",ylabel="Frozen class",title="Extra Trees confusion (%)")
    save_both("d",draw_d)

    def draw_e(fig,ax,labeled):
        top=fi.head(8).sort_values("ExtraTrees_median")
        # Keep electrophysiological display names identical to the RRR figure.
        short={"Epsy_ahp_delay_5spike":"AHP delay","Epsy_ahp_delay_ratio_5spike":"AHP ratio","Epsy_upstroke_adapt_ratio":"Upstroke adapt","Epsy_width_adapt_ratio":"Width adapt","Epsy_upstroke_downstroke_ratio_rheo":"Up/down ratio","Epsy_threshold_v_adapt_ratio":"Threshold adapt","Epsy_downstroke_adapt_ratio":"Downstroke adapt","Epsy_peak_deltav_rheo":"Peak ΔV"}
        ax.barh(range(len(top)),top.ExtraTrees_median,color="#4C78A8",height=.62)
        ax.errorbar(top.ExtraTrees_median,range(len(top)),xerr=[top.ExtraTrees_median-top.ExtraTrees_q025,top.ExtraTrees_q975-top.ExtraTrees_median],fmt="none",ecolor="#222",elinewidth=.4,capsize=1)
        ax.set_yticks(range(len(top)),top.feature.map(short),fontsize=3.2); base.clean(ax)
        if labeled: ax.set(xlabel="Extra Trees importance",title="Stable E4 features")
    save_both("e",draw_e)

    def draw_f(fig,ax,labeled):
        pv=sens.pivot(index="NPC",columns="K",values="ARI_median").loc[range(2,11),range(2,9)]
        sns.heatmap(pv,cmap="viridis",vmin=0,vmax=.8,annot=labeled,fmt=".2f",annot_kws={"fontsize":2.8},
                    cbar=True,cbar_kws={"label":"Median ARI","pad":.03},ax=ax)
        ax.add_patch(plt.Rectangle((4-2,3-2),1,1,fill=False,edgecolor="white",lw=.8))
        if labeled: ax.set(xlabel="K",ylabel="NPC",title="NPC and K sensitivity")
    save_both("f",draw_f)

    rows=[]
    for folder in ["with_text","no_text"]:
        from PIL import Image
        for p in sorted((OUT/folder).glob("*.png")):
            with Image.open(p) as im: rows.append({"file":str(p.relative_to(ROOT)),"width_px":im.width,"height_px":im.height,"dpi":im.info.get("dpi",(None,None))[0]})
    pd.DataFrame(rows).to_csv(OUT/"manifest.csv",index=False)
    assert all(abs(r["dpi"]-600)<1 for r in rows)
    print(f"Exported {len(rows)} panels with text-adaptive canvases at 600 dpi")


if __name__=="__main__": main()
