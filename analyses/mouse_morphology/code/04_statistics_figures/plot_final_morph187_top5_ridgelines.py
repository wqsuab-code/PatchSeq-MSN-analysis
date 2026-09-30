#!/usr/bin/env python
"""One top-five morphology ridgeline figure per final M class."""
from pathlib import Path
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.interpolate import PchipInterpolator

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "outputs/morph_qc/final_morph187_NPC3_HCK4_GCres0p50_figures"
ZFILE = BASE / "08_heatmap_ordered_zscore_matrix.csv"
ASSIGN = BASE / "00_final_morph187_cell_assignments.csv"
OUT = BASE / "28_top5_feature_ridgelines_global_zscore_minus6_plus3_by_M"
LEVELS = ["M1", "M2", "M3", "M4"]
COLORS = {"M1":"#00468B", "M2":"#42B540", "M3":"#ED0000", "M4":"#0099B4"}
METRICS = ["Soma AR", "Soma circ.", "Primary N", "Bifurcations", "Total length"]
FEATURES = [
    "M_soma_aspect_ratio", "M_soma_circularity_index", "M_total_number_of_neurites",
    "M_Number_of_bifurcation_points", "M_Total_neurite_length_(sections)"
]
BW = 0.35

def smooth_profile(values: np.ndarray, bases: np.ndarray, xmin: float, xmax: float) -> tuple[np.ndarray, np.ndarray]:
    """Shape-preserving smooth curve through the five feature points."""
    y = bases[::-1]
    x = np.asarray(values, dtype=float)[::-1]
    y_dense = np.linspace(float(y.min()), float(y.max()), 240)
    x_dense = PchipInterpolator(y, x)(y_dense)
    return np.clip(x_dense, xmin, xmax), y_dense

def fixed_kde(values: np.ndarray, grid: np.ndarray, bandwidth: float=BW) -> np.ndarray:
    """Fixed-bandwidth Gaussian KDE in the common global-Z coordinate system."""
    values = np.asarray(values, dtype=float)
    u = (grid[:, None] - values[None, :]) / bandwidth
    density = np.exp(-0.5*u*u).sum(axis=1) / (len(values)*bandwidth*np.sqrt(2*np.pi))
    maximum = density.max()
    return density/maximum if maximum > 0 else density

def main():
    OUT.mkdir(parents=True, exist_ok=True)
    z = pd.read_csv(ZFILE)[["MSN_unique_ID"] + FEATURES]
    a = pd.read_csv(ASSIGN)[["MSN_unique_ID", "M_class"]]
    wide = a.merge(z, on="MSN_unique_ID", validate="one_to_one")
    data = wide.melt(id_vars=["MSN_unique_ID", "M_class"], value_vars=FEATURES,
                     var_name="Feature", value_name="Scaled_value")
    label_map = dict(zip(FEATURES, METRICS))
    data["Metric"] = pd.Categorical(data["Feature"].map(label_map), categories=METRICS, ordered=True)
    data["Global_Z_raw"] = data["Scaled_value"]
    # Display winsorization only; the exported raw-Z column remains unchanged.
    data["Scaled_value"] = data["Global_Z_raw"].clip(-6, 3)
    xmin, xmax = -6.0, 3.0
    grid = np.linspace(xmin, xmax, 800)
    mpl.rcParams.update({
        "font.family":"Arial", "font.size":4, "axes.linewidth":.55,
        "pdf.fonttype":42, "ps.fonttype":42, "svg.fonttype":"none",
        "savefig.facecolor":"white", "figure.facecolor":"white",
    })
    manifest=[]
    for m in LEVELS:
        subset=data[data.M_class.eq(m)].copy()
        n=subset.MSN_unique_ID.nunique()
        fig,ax=plt.subplots(figsize=(2.45,2.25))
        bases=np.arange(len(METRICS))[::-1].astype(float)
        medians=np.array([
            np.median(subset.loc[subset.Metric.eq(metric),"Scaled_value"].to_numpy(float))
            for metric in METRICS
        ])
        # One thin profile per cell.  All paths share the same five axes and
        # sit beneath the semi-transparent ridge fills.
        cell_profiles=(subset.pivot(index="MSN_unique_ID",columns="Metric",values="Scaled_value")
                       .reindex(columns=METRICS))
        for profile in cell_profiles.to_numpy(float):
            xx,yy=smooth_profile(profile,bases,xmin,xmax)
            ax.plot(xx,yy,color="#8C8C8C",alpha=.30,lw=.25,ls="-",zorder=.35)
        # Shared reference thresholds across all four M-class panels.
        for xref in (-2, 0, 2):
            ax.axvline(xref, color="#ED0000", lw=.40, ls="--", alpha=.85, zorder=.25)
        # A darker grey median curve remains distinguishable from cell curves.
        xx,yy=smooth_profile(medians,bases,xmin,xmax)
        ax.plot(xx,yy,color="#4D4D4D",alpha=.90,lw=.25,ls="-",zorder=.65)
        for base,metric in zip(bases,METRICS):
            vals=subset.loc[subset.Metric.eq(metric),"Scaled_value"].to_numpy(float)
            curve=fixed_kde(vals,grid)*.82
            ax.fill_between(grid,base,base+curve,color=COLORS[m],alpha=.30,linewidth=0,zorder=1)
            ax.plot(grid,base+curve,color=COLORS[m],lw=.55,zorder=2)
            ax.plot(vals,np.full_like(vals,base-.035),"|",color=COLORS[m],ms=2.0,mew=.35,alpha=.65,zorder=3)
            ax.scatter(np.median(vals),base,s=5.0,c="black",edgecolors="none",zorder=4)
            manifest.append({"M_class":m,"N":n,"Metric":metric,"Median_scaled":np.median(vals),
                             "Mean_scaled":np.mean(vals),"SD_scaled":np.std(vals,ddof=1),
                             "KDE_bandwidth_Z":BW})
        ax.set_xlim(xmin,xmax); ax.set_ylim(-.14,4.95)
        ax.set_yticks(bases+.20,METRICS,fontsize=4)
        ticks = [-6, -4, -2, 0, 2, 3]
        ax.set_xticks(ticks, [f"{v:g}" for v in ticks])
        ax.set_xlabel("Global Z score",fontsize=4,labelpad=2)
        ax.set_title(f"{m} (n={n})",fontsize=6,fontweight="bold",pad=2)
        ax.spines[["top","right","left"]].set_visible(False)
        ax.spines["bottom"].set_linewidth(.55)
        ax.tick_params(axis="x",labelsize=4,width=.5,length=2,pad=1)
        ax.tick_params(axis="y",width=0,length=0,pad=2)
        fig.subplots_adjust(left=.31,right=.98,top=.93,bottom=.15)
        stem=OUT/f"28_{m}_Top5_Ridgeline_GlobalZminus6plus3_GrayCurves_CommonScale_W2p45_H2p25"
        fig.savefig(stem.with_suffix(".png"),dpi=900,bbox_inches="tight",pad_inches=.02)
        fig.savefig(stem.with_suffix(".pdf"),bbox_inches="tight",pad_inches=.02)
        fig.savefig(stem.with_suffix(".svg"),bbox_inches="tight",pad_inches=.02)
        plt.close(fig)
    pd.DataFrame(manifest).to_csv(OUT/"28_M1-M4_Top5_Ridgeline_GlobalZminus6plus3_distribution_summary.csv",index=False)
    data.to_csv(OUT/"28_M1-M4_Top5_GlobalZ_raw_and_display_clipped_values.csv",index=False)
    (OUT/"28_M1-M4_Top5_Ridgeline_GlobalZminus6plus3_run_log.txt").write_text(
        "One independent figure per M class.\n"
        "Input: feature-wise Z scores after the frozen shift/sum*10000/log1p preprocessing, using all 187 cells as the reference.\n"
        "All figures share x=-6 to +3, fixed Gaussian KDE bandwidth=0.35, and the same 2.45 x 2.25 inch canvas.\n"
        "Values outside [-6,+3] are clipped only for display; raw global Z scores are preserved in the exported CSV.\n"
        "Red vertical dashed reference lines are drawn at global Z=-2, 0, and +2.\n"
        "Ridge fill alpha=0.30 and is layered above the connecting paths; coloured baseline ticks are individual cells; black point is the median.\n"
        "Each cell is connected across the five metrics by a PCHIP curve: grey #8C8C8C, alpha=0.30, 0.25 pt.\n"
        "The median PCHIP curve is dark grey #4D4D4D, alpha=0.90, 0.25 pt; median points remain black.\n",
        encoding="utf-8")
    print(OUT)

if __name__=="__main__": main()
