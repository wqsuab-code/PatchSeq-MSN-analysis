#!/usr/bin/env python3
"""Consensus M1-M4 morphology comparison for 10 representative features."""
from itertools import combinations
from math import erfc, sqrt
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.patches import Ellipse
from matplotlib.ticker import MaxNLocator
import numpy as np
import pandas as pd
from scipy.stats import chi2, kruskal, rankdata

BASE = Path(r"C:\Users\53461\OneDrive\Documentos\Patch-seq_Mouse_Acb_MSN_T-type_Visualization\macaque_m")
RUN = BASE / "m18_tempfreeze_NPC5_HCK4_res2.3"
OUT = RUN / "M4_compare_all18"
ORDER = ["M1", "M2", "M3", "M4"]
COLORS = {"M1": "#1F77B4", "M2": "#D9A400", "M3": "#8C564B", "M4": "#E377C2"}
FEATURES = [
    ("basal_dendrite_total_length", "Total length", "um"),
    ("basal_dendrite_num_branches", "Number of branches", ""),
    ("basal_dendrite_calculate_number_of_stems", "Number of stems", ""),
    ("basal_dendrite_max_branch_order", "Max branch order", ""),
    ("basal_dendrite_max_path_distance", "Max path distance", "um"),
    ("basal_dendrite_max_euclidean_distance", "Max Euclidean distance", "um"),
    ("basal_dendrite_mean_contraction", "Mean contraction", "ratio"),
    ("basal_dendrite_mean_diameter", "Mean diameter", "um"),
    ("basal_dendrite_extent_dorsal", "Dorsal extent", "um"),
    ("basal_dendrite_extent_medial", "Medial extent", "um"),
    ("basal_dendrite_bias_dorsal", "Dorsal bias", "a.u."),
    ("basal_dendrite_bias_medial", "Medial bias", "a.u."),
    ("basal_dendrite_soma_percentile_dorsal", "Dorsal soma percentile", ""),
    ("basal_dendrite_soma_percentile_medial", "Medial soma percentile", ""),
    ("basal_dendrite_stem_exit_dorsal", "Dorsal stem exit", "a.u."),
    ("basal_dendrite_stem_exit_ventral", "Ventral stem exit", "a.u."),
    ("basal_dendrite_stem_exit_MedialLateral", "ML stem exit", "a.u."),
    ("soma_surface_area", "Soma surface area", "um2"),
]
CMAP = LinearSegmentedColormap.from_list("morph_z", ["#2166AC", "#C9DCEB", "#F3F3F3", "#E5B9BA", "#B2182B"])


def adjust(p, method):
    p = np.asarray(p, float); n = len(p); out = np.empty_like(p)
    if method == "holm":
        order = np.argsort(p); running = 0.0
        for j, i in enumerate(order):
            running = max(running, min(1.0, (n-j)*p[i])); out[i] = running
    else:
        order = np.argsort(p)[::-1]; running = 1.0
        for j, i in enumerate(order, 1):
            rank = n-j+1; running = min(running, min(1.0, p[i]*n/rank)); out[i] = running
    return out


def dunn(values, groups):
    values=np.asarray(values,float);groups=np.asarray(groups,str);ranks=rankdata(values);n=len(values)
    _,ties=np.unique(values,return_counts=True)
    variance=n*(n+1)/12-np.sum(ties**3-ties)/(12*(n-1))
    means={g:ranks[groups==g].mean() for g in ORDER};sizes={g:(groups==g).sum() for g in ORDER};rows=[]
    for g1,g2 in combinations(ORDER,2):
        se=sqrt(variance*(1/sizes[g1]+1/sizes[g2]));z=(means[g1]-means[g2])/se
        rows.append([g1,g2,z,abs(z)/sqrt(n),erfc(abs(z)/sqrt(2))])
    out=pd.DataFrame(rows,columns=["Group_1","Group_2","Dunn_Z","Dunn_r","P_raw"])
    out["P_adj_Holm_within_feature"]=adjust(out.P_raw,"holm")
    return out


def p_label(p):
    return f"{p:.1e}" if p < 1e-3 else (f"{p:.3f}" if p < .01 else f"{p:.2f}")


def distribution(ax,data,feature,title,unit,pairwise,seed):
    rng=np.random.default_rng(seed)
    for x,g in enumerate(ORDER,1):
        v=data.loc[data.M==g,feature].to_numpy(float)
        ax.scatter(x+rng.uniform(-.16,.16,len(v)),v,s=2.8,color=COLORS[g],alpha=.62,linewidths=0,rasterized=True)
        ax.boxplot([v],positions=[x],widths=.43,whis=(2.5,97.5),showfliers=False,patch_artist=True,manage_ticks=False,
                   boxprops={"facecolor":"none","edgecolor":"#333","linewidth":.4},medianprops={"color":"#111","linewidth":.55},
                   whiskerprops={"color":"#555","linewidth":.4},capprops={"color":"#555","linewidth":.4})
    allv=data[feature].to_numpy(float);lo,hi=allv.min(),allv.max();span=max(hi-lo,abs(hi)*.05,1e-6)
    sig=pairwise[pairwise.P_adj_Holm_within_feature<.05].sort_values("P_adj_Holm_within_feature").head(4)
    ax.set_ylim(lo-.08*span,hi+span*(.22+.11*len(sig)))
    for level,(_,r) in enumerate(sig.iterrows()):
        x1=ORDER.index(r.Group_1)+1;x2=ORDER.index(r.Group_2)+1;y=hi+span*(.12+.105*level);tick=.025*span
        ax.plot([x1,x1,x2,x2],[y-tick,y,y,y-tick],color="#333",lw=.38,clip_on=False)
        ax.text((x1+x2)/2,y+.012*span,p_label(r.P_adj_Holm_within_feature),ha="center",va="bottom",fontsize=3.1)
    ax.set_title(f"{title} ({unit})" if unit else title,fontsize=4.4,pad=1.5)
    ax.set_xticks(range(1,5),ORDER);ax.yaxis.set_major_locator(MaxNLocator(5))
    ax.tick_params(labelsize=3.2,length=1.4,width=.35,pad=.8);ax.spines[["top","right"]].set_visible(False)
    ax.spines[["left","bottom"]].set_linewidth(.4)


def tsne(ax,data,feature):
    ax.scatter(data.x,data.y,c=data[feature+"__z"],cmap=CMAP,norm=Normalize(-2.5,2.5),s=3.1,alpha=.8,linewidths=0,rasterized=True)
    for g in ORDER:
        p=data.loc[data.M==g,["x","y"]].to_numpy(float);center=p.mean(0);cov=np.cov(p,rowvar=False)
        vals,vecs=np.linalg.eigh(cov);idx=vals.argsort()[::-1];vals,vecs=vals[idx],vecs[:,idx]
        w,h=2*np.sqrt(np.maximum(vals,0)*chi2.ppf(.85,2));ang=np.degrees(np.arctan2(vecs[1,0],vecs[0,0]))
        ax.add_patch(Ellipse(center,w,h,angle=ang,fill=False,edgecolor=COLORS[g],linewidth=.65,linestyle=(0,(2.2,1.2))))
    ax.set_aspect("equal");ax.set_axis_off()


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    a=pd.read_csv(RUN/"01_temp_frozen_assignments_126.csv");a=a[a.concordant==True].copy();a["M"]="M"+a.HC_K4.astype(str)
    raw=pd.read_csv(BASE/"m18"/"01_raw_126.csv");z=pd.read_csv(BASE/"m18_adaptive_pca126"/"02_transformed_z_117.csv")
    xy=pd.read_csv(RUN/"09_tSNE_optimized_coordinates.csv")
    f=[x[0] for x in FEATURES]
    data=a[["cell_label","M"]].merge(raw[["cell_label"]+f],on="cell_label",validate="one_to_one")
    data=data.merge(z[["cell_label"]+f].rename(columns={x:x+"__z" for x in f}),on="cell_label",validate="one_to_one")
    data=data.merge(xy[["cell_label","tSNE1","tSNE2"]],on="cell_label",validate="one_to_one")
    centers=data.groupby("M")[["tSNE1","tSNE2"]].transform("mean");data["x"]=centers.tSNE1+.85*(data.tSNE1-centers.tSNE1);data["y"]=centers.tSNE2+.85*(data.tSNE2-centers.tSNE2)
    assert len(data)==117 and data.M.value_counts().sort_index().tolist()==[43,42,20,12]

    omnibus=[];pairs=[]
    for feature,title,unit in FEATURES:
        arrays=[data.loc[data.M==g,feature].to_numpy(float) for g in ORDER];h,p=kruskal(*arrays)
        omnibus.append([feature,title,unit,h,3,p]);pw=dunn(data[feature],data.M);pw.insert(0,"Feature",feature);pairs.append(pw)
    omnibus=pd.DataFrame(omnibus,columns=["Feature","Display","Unit","Kruskal_H","df","P_raw"])
    omnibus["Q_BH_across_18_features"]=adjust(omnibus.P_raw,"bh");pairs=pd.concat(pairs,ignore_index=True)

    mpl.rcParams.update({"font.family":"Arial","font.size":4,"pdf.fonttype":42,"ps.fonttype":42})
    # Publication-size overview: three feature groups per row on a fixed 7.2 x 8.76 inch canvas.
    fig=plt.figure(figsize=(7.2,8.76));outer=fig.add_gridspec(6,3,left=.045,right=.985,bottom=.055,top=.985,wspace=.10,hspace=.34)
    for i,(feature,title,unit) in enumerate(FEATURES):
        row,col=divmod(i,3);inner=outer[row,col].subgridspec(1,2,width_ratios=[.93,1.07],wspace=.025)
        ax1=fig.add_subplot(inner[0,0]);ax2=fig.add_subplot(inner[0,1]);distribution(ax1,data,feature,title,unit,pairs[pairs.Feature==feature],777+i);tsne(ax2,data,feature)
    cax=fig.add_axes([.73,.018,.245,.007]);cb=mpl.colorbar.ColorbarBase(cax,cmap=CMAP,norm=Normalize(-2.5,2.5),orientation="horizontal",ticks=[-2.5,0,2.5]);cb.ax.tick_params(labelsize=3.2,length=1.2,width=.35,pad=.5);cb.outline.set_linewidth(.35);cb.set_label("Feature Z-score",fontsize=3.5,labelpad=.15)
    stem="Macaque_Morphology_ConsensusM4_all18_3perrow_W7p2_H8p76"
    png=OUT/f"{stem}.png";pdf=OUT/f"{stem}.pdf"
    fig.savefig(png,dpi=900,facecolor="white")
    fig.savefig(pdf,facecolor="white")
    plt.close(fig)

    individual=OUT/"individual";individual.mkdir(exist_ok=True)
    for i,(feature,title,unit) in enumerate(FEATURES,1):
        one=plt.figure(figsize=(3.0,2.0));gs=one.add_gridspec(1,2,left=.08,right=.98,bottom=.18,top=.92,width_ratios=[.93,1.07],wspace=.03)
        ax1=one.add_subplot(gs[0,0]);ax2=one.add_subplot(gs[0,1])
        distribution(ax1,data,feature,title,unit,pairs[pairs.Feature==feature],777+i);tsne(ax2,data,feature)
        cax=one.add_axes([.64,.065,.30,.022]);cb=mpl.colorbar.ColorbarBase(cax,cmap=CMAP,norm=Normalize(-2.5,2.5),orientation="horizontal",ticks=[-2.5,0,2.5]);cb.ax.tick_params(labelsize=3.2,length=1,width=.35,pad=.4);cb.outline.set_linewidth(.35);cb.set_label("Feature Z-score",fontsize=3.4,labelpad=.1)
        slug=feature.replace("basal_dendrite_","")
        one.savefig(individual/f"{i:02d}_{slug}_M1-M4_distribution_tSNE.png",dpi=600,facecolor="white")
        plt.close(one)
    assert len(list(individual.glob("*.png")))==18

    omnibus.to_csv(OUT/"all18_KruskalWallis_BH.csv",index=False);pairs.to_csv(OUT/"all18_Dunn_Holm.csv",index=False);data.to_csv(OUT/"all18_plot_data_117.csv",index=False)
    print(png);print(pdf);print(omnibus[["Display","Kruskal_H","P_raw","Q_BH_across_18_features"]].to_string(index=False));print(data.M.value_counts().sort_index().to_string())


if __name__=="__main__": main()
