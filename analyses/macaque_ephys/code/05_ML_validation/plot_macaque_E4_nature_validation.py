#!/usr/bin/env python
"""Nature-style main and Extended Data figures for Macaque E4 ML500 results."""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import seaborn as sns
from sklearn.cluster import AgglomerativeClustering

ROOT=Path(__file__).resolve().parents[1]
SRC=ROOT/"outputs"/"Macaque_E4_full_ML500"
ASSIGN=ROOT/"outputs"/"R3_panels"/"00_res3_main_cell_assignments_and_tsne.csv"
TUNED_COORDS=ROOT/"outputs"/"R3_panels"/"04b_tuned_consensus368_p80_ee12_random_coordinates.csv"
OUT=ROOT/"outputs"/"Macaque_E4_Nature_validation"
COLORS={"C1":"#F8766D","C2":"#7CAE00","C3":"#00BFC4","C4":"#C77CFF"}
ALG_ORDER=["Ward","KMeans","GMM","Spectral"]
MODEL_ORDER=["Extra trees","kNN","RBF SVM","Linear SVM","Multinomial logistic","Random forest","HistGradientBoosting"]
MODEL_SHORT={"Extra trees":"ET","kNN":"kNN","RBF SVM":"RBF","Linear SVM":"Linear",
             "Multinomial logistic":"Logit","Random forest":"RF","HistGradientBoosting":"HGB"}
TSNE_XLIM=(-11.1203062245,10.5763037016)
TSNE_YLIM=(-9.4204480594,11.2429899654)


def style():
    mpl.rcParams.update({"font.family":"Arial","font.size":6,"axes.titlesize":7,"axes.labelsize":6,
        "xtick.labelsize":5.5,"ytick.labelsize":5.5,"legend.fontsize":5.5,"axes.linewidth":.55,
        "xtick.major.width":.5,"ytick.major.width":.5,"xtick.major.size":2.4,"ytick.major.size":2.4,
        "pdf.fonttype":42,"ps.fonttype":42,"savefig.facecolor":"white"})


def panel(ax,letter):
    ax.text(-.18,1.08,letter,transform=ax.transAxes,fontsize=8,fontweight="bold",va="top",ha="left")


def ellipse80(ax,xy,color):
    cov=np.cov(xy.T); vals,vecs=np.linalg.eigh(cov); order=vals.argsort()[::-1]; vals=vals[order]; vecs=vecs[:,order]
    ang=np.degrees(np.arctan2(vecs[1,0],vecs[0,0])); scale=np.sqrt(3.2188758248682006)
    e=Ellipse(xy.mean(0),2*scale*np.sqrt(vals[0]),2*scale*np.sqrt(vals[1]),angle=ang,
              facecolor=color,edgecolor=color,alpha=.10,lw=.75,ls="--",zorder=1)
    ax.add_patch(e)


def clean(ax):
    ax.spines[["top","right"]].set_visible(False)


def freeze_tsne_axes(ax):
    ax.set_xlim(*TSNE_XLIM); ax.set_ylim(*TSNE_YLIM); ax.set_aspect("equal",adjustable="box")


def main():
    OUT.mkdir(parents=True,exist_ok=True); style()
    d=pd.read_csv(ASSIGN).drop(columns=["tSNE1","tSNE2"])
    tuned=pd.read_csv(TUNED_COORDS,usecols=["cell_label","tSNE1","tSNE2"])
    d=d.merge(tuned,on="cell_label",validate="one_to_one")
    d=d[d.Consensus.astype(str).str.lower().eq("true")].copy()
    uns=pd.read_csv(SRC/"08_cross_algorithm_E4_donor_resampling_500.csv")
    perf=pd.read_csv(SRC/"04_supervised_500_all_models_class_recall.csv")
    summ=pd.read_csv(SRC/"05_supervised_summary_95CI.csv")
    fi=pd.read_csv(SRC/"07_feature_importance_stability_500.csv")
    sens=pd.read_csv(SRC/"03_NPC_K_sensitivity_summary.csv")
    sens_raw=pd.read_csv(SRC/"02_NPC2-10_K2-8_donor_resampling_500.csv")
    cell=pd.read_csv(SRC/"11_cell_level_consensus_stability.csv")
    cm=pd.read_csv(SRC/"06_confusion_sum_Extra_trees.csv",index_col=0).to_numpy(float)
    cm=cm/cm.sum(1,keepdims=True)*100

    # 183 mm wide, publication-facing 2x3 composition.
    fig=plt.figure(figsize=(7.205,4.95)); gs=fig.add_gridspec(2,3,width_ratios=[1.04,1,1.12],hspace=.52,wspace=.72)
    ax=fig.add_subplot(gs[0,0]); panel(ax,"a")
    for c in COLORS:
        q=d[d.HC_class==c]; xy=q[["tSNE1","tSNE2"]].to_numpy(); ellipse80(ax,xy,COLORS[c])
        ax.scatter(xy[:,0],xy[:,1],s=6,c=COLORS[c],edgecolors="white",linewidths=.18,label=f"{c}  n={len(q)}",zorder=2)
    ax.set(xlabel="t-SNE 1",ylabel="t-SNE 2",title="Frozen E4 consensus")
    ax.set_xticks([]); ax.set_yticks([]); freeze_tsne_axes(ax); clean(ax); ax.legend(frameon=False,loc="lower right",handletextpad=.25,borderpad=.1,fontsize=4.6)

    ax=fig.add_subplot(gs[0,1]); panel(ax,"b")
    sns.violinplot(data=uns,x="algorithm",y="ARI_vs_reference",order=ALG_ORDER,color="#A6CEE3",inner=None,cut=0,linewidth=.45,ax=ax)
    sns.boxplot(data=uns,x="algorithm",y="ARI_vs_reference",order=ALG_ORDER,width=.22,showfliers=False,color="white",linewidth=.55,ax=ax)
    ax.axhline(.75,color="#D55E00",ls="--",lw=.6); ax.set(xlabel="",ylabel="ARI versus frozen E4",title="Donor-resampled discovery")
    ax.yaxis.label.set_size(5); ax.yaxis.labelpad=1
    ax.set_ylim(0,1.02); ax.tick_params(axis="x",rotation=25); clean(ax)

    ax=fig.add_subplot(gs[0,2]); panel(ax,"c")
    tmp=perf.copy(); tmp["model"]=tmp.algorithm.map(MODEL_SHORT)
    short_order=[MODEL_SHORT[x] for x in MODEL_ORDER]
    sns.boxplot(data=tmp,x="model",y="balanced_accuracy",order=short_order,color="#009E73",width=.55,showfliers=False,linewidth=.55,ax=ax)
    ax.axhline(.25,color="#777777",ls=":",lw=.65); ax.set(xlabel="",ylabel="",title="Held-out donor balanced accuracy")
    ax.set_ylim(.2,1.02); ax.tick_params(axis="x",rotation=0); clean(ax)

    ax=fig.add_subplot(gs[1,0]); panel(ax,"d")
    sns.heatmap(cm,cmap="Blues",vmin=0,vmax=100,annot=True,fmt=".0f",square=True,cbar=False,
                xticklabels=["C1","C2","C3","C4"],yticklabels=["C1","C2","C3","C4"],annot_kws={"fontsize":5.5},ax=ax)
    ax.set(xlabel="Predicted class",ylabel="Frozen class",title="Extra Trees confusion (%)")
    ax.tick_params(length=0)

    ax=fig.add_subplot(gs[1,1]); panel(ax,"e")
    top=fi.head(8).sort_values("ExtraTrees_median")
    rename={
      "Epsy_ahp_delay_5spike":"AHP delay", "Epsy_ahp_delay_ratio_5spike":"AHP delay ratio",
      "Epsy_upstroke_adapt_ratio":"Upstroke adaptation", "Epsy_width_adapt_ratio":"AP-width adaptation",
      "Epsy_upstroke_downstroke_ratio_rheo":"Up/down ratio (rheo)",
      "Epsy_threshold_v_adapt_ratio":"Threshold-V adaptation",
      "Epsy_downstroke_adapt_ratio":"Downstroke adaptation", "Epsy_peak_deltav_rheo":"Peak ΔV (rheo)"}
    names=top.feature.map(rename)
    ax.barh(range(len(top)),top.ExtraTrees_median,color="#4C78A8",height=.62)
    ax.errorbar(top.ExtraTrees_median,range(len(top)),xerr=[top.ExtraTrees_median-top.ExtraTrees_q025,top.ExtraTrees_q975-top.ExtraTrees_median],fmt="none",ecolor="#222222",elinewidth=.55,capsize=1.5)
    ax.set_yticks(range(len(top)),names,fontsize=4.2); ax.set(xlabel="Extra Trees importance",title="Stable E4 features"); clean(ax)

    ax=fig.add_subplot(gs[1,2]); panel(ax,"f")
    npcs=list(range(2,11)); ks=list(range(2,9))
    ari=sens.pivot(index="NPC",columns="K",values="ARI_median").loc[npcs,ks]
    sil=sens.pivot(index="NPC",columns="K",values="silhouette_median").loc[npcs,ks]
    sf=(sens_raw.assign(small=lambda q:q.minimum_cluster_n<=2).groupby(["NPC","K"]).small.mean().unstack().loc[npcs,ks])
    full_min={}
    for npc in npcs:
        coords=d[[f"PC{i}" for i in range(1,npc+1)]].to_numpy()
        for k in ks:
            lab=AgglomerativeClustering(n_clusters=k,linkage="ward").fit_predict(coords)
            full_min[(npc,k)]=int(pd.Series(lab).value_counts().min())
    im=ax.imshow(ari,cmap="viridis",vmin=0,vmax=.8,aspect="auto")
    for i,npc in enumerate(npcs):
        for j,k in enumerate(ks):
            v=ari.loc[npc,k]; sv=sil.loc[npc,k]; freq=sf.loc[npc,k]
            ax.text(j,i-.06,f"{v:.2f}\n{sv:.2f}",ha="center",va="center",fontsize=3.7,
                    color="white" if v<.55 else "#222222",linespacing=.85)
            ax.add_patch(plt.Rectangle((j-.43,i+.31),.86*freq,.07,color="#D55E00",lw=0))
            if full_min[(npc,k)]<=2:
                ax.add_patch(plt.Rectangle((j-.49,i-.49),.98,.98,fill=False,edgecolor="#D55E00",lw=.7))
    # Frozen specification.
    ax.add_patch(plt.Rectangle((4-2-.49,3-2-.49),.98,.98,fill=False,edgecolor="white",lw=1.25))
    ax.set_xticks(range(len(ks)),ks); ax.set_yticks(range(len(npcs)),npcs)
    ax.set(xlabel="K",ylabel="NPC",title="ARI / silhouette sensitivity")
    cb=fig.colorbar(im,ax=ax,fraction=.045,pad=.025); cb.set_label("Median ARI",fontsize=5); cb.ax.tick_params(labelsize=5)
    ax.text(0,-.22,"red bar: P(min n≤2); red box: full min n≤2",transform=ax.transAxes,fontsize=4.1,ha="left")

    fig.text(.5,.006,"Full-data sizes: PC3/K4 64,65,122,139; PC3/K5 1,64,64,122,139; PC3/K6 1,16,64,64,106,139; PC4/K4 1,64,82,243.",ha="center",fontsize=4.3)
    fig.align_labels(); fig.savefig(OUT/"Macaque_E4_ML500_Nature_main.pdf",bbox_inches="tight")
    fig.savefig(OUT/"Macaque_E4_ML500_Nature_main.png",dpi=600,bbox_inches="tight")
    fig.savefig(OUT/"Macaque_E4_ML500_Nature_main.tiff",dpi=600,bbox_inches="tight",pil_kwargs={"compression":"tiff_lzw"}); plt.close(fig)

    # Extended Data 1: parameter sensitivity and cell-level stability.
    fig,ax=plt.subplots(1,2,figsize=(7.205,2.75),gridspec_kw={"wspace":.38})
    panel(ax[0],"a"); pv=sens.pivot(index="NPC",columns="K",values="ARI_median")
    sns.heatmap(pv,cmap="viridis",vmin=0,vmax=.8,annot=True,fmt=".2f",annot_kws={"fontsize":5},cbar_kws={"label":"Median ARI"},ax=ax[0])
    ax[0].add_patch(plt.Rectangle((4-2,3-2),1,1,fill=False,edgecolor="white",lw=1.2)); ax[0].set_title("NPC and K sensitivity")
    panel(ax[1],"b"); sns.histplot(cell.consensus_margin,bins=25,color="#4C78A8",edgecolor="white",linewidth=.3,ax=ax[1])
    ax[1].axvline(0,color="#D55E00",ls="--",lw=.7); ax[1].set(xlabel="Within-class minus maximum other-class\nco-clustering probability",ylabel="Cells",title="Cell-level consensus margin"); clean(ax[1])
    fig.savefig(OUT/"Macaque_E4_ML500_Nature_ExtendedData1.pdf",bbox_inches="tight"); fig.savefig(OUT/"Macaque_E4_ML500_Nature_ExtendedData1.png",dpi=600,bbox_inches="tight"); plt.close(fig)

    # Source-data tables used directly by the figures.
    uns.to_csv(OUT/"SourceData_main_b.csv",index=False); perf.to_csv(OUT/"SourceData_main_c.csv",index=False)
    pd.DataFrame(cm,index=["C1","C2","C3","C4"],columns=["C1","C2","C3","C4"]).to_csv(OUT/"SourceData_main_d.csv")
    fi.to_csv(OUT/"SourceData_main_e.csv",index=False)
    sens_display=sens.merge(sf.stack().rename("fraction_resamples_min_cluster_le2").reset_index(),on=["NPC","K"])
    sens_display["full_data_minimum_cluster_n"]=[full_min[(int(a),int(b))] for a,b in zip(sens_display.NPC,sens_display.K)]
    sens_display.to_csv(OUT/"SourceData_main_f.csv",index=False)
    summ.to_csv(OUT/"SourceData_class_recall.csv",index=False); sens.to_csv(OUT/"SourceData_ED1a.csv",index=False); cell.to_csv(OUT/"SourceData_ED1b.csv",index=False)
    (OUT/"FIGURE_SPEC.txt").write_text("""Nature-style specification
Main width: 183 mm (double column); Arial; panel letters 8 pt bold; other text 5-7 pt; line widths 0.25-1 pt.
RGB, colour-blind-friendly class palette. PDF retains vector text/lines; PNG/TIFF are 600 dpi.
Panel a: frozen consensus cells and 80% covariance ellipses. Panel b: 500 donor-resampled unsupervised fits.
Panel c-d: 500 donor-grouped train/test splits. Panel e: median feature importance and empirical 95% interval.
Panel f: NPC/K sensitivity; cell colour is ARI, text is ARI/silhouette, red bars show small-cluster frequency,
red boxes flag full-data solutions with minimum cluster n<=2, and the white box marks frozen PC3/K4.
Extended Data retains a larger sensitivity rendering and cell-level co-clustering margin.
""",encoding="utf-8")
    print(OUT)

if __name__=="__main__": main()
