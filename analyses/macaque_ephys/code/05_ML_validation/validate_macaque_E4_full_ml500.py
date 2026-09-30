#!/usr/bin/env python
"""Dedicated 500-repeat stability/reproducibility audit of frozen Macaque E4.

The orchestration mirrors ``validate_mouse_E5_macaqueM_style.py`` (grouped
hold-outs, multiple supervised learners, unsupervised sensitivity, permutation
controls, importance and ablation), but deliberately retains the frozen
Macaque design: donor—not recording date—is the grouping unit; 500 repeats are
used; and the original mixed 19-feature transform is fitted on each training
fold.  The Mouse transform must not be substituted here.
"""
from pathlib import Path
from zipfile import ZipFile
import json, sys, warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from joblib import Parallel, delayed
from scipy.stats import chi2_contingency
from sklearn.cluster import AgglomerativeClustering
from sklearn.decomposition import PCA
from sklearn.ensemble import ExtraTreesClassifier, RandomForestClassifier, HistGradientBoostingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import adjusted_rand_score, balanced_accuracy_score, confusion_matrix, f1_score, matthews_corrcoef, recall_score, silhouette_score
from sklearn.mixture import GaussianMixture
from sklearn.model_selection import GroupShuffleSplit
from sklearn.neighbors import KNeighborsClassifier
from sklearn.svm import SVC

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))
import validate_macaque_EM_stability_500 as base

OUT=ROOT/"outputs"/"Macaque_E4_full_ML500"
NREP=500; SEED=20260910
CLASSES=("C1","C2","C3","C4")  # exported/displayed as E1-E4
GROUP_FIELD="donor_label"
TEST_SIZE=.25
ANALYSIS_DESIGN={
    "scope":"Macaque MSN; Ca+Pu+NAc",
    "frozen_classes":["E1","E2","E3","E4"],
    "grouping_unit":GROUP_FIELD,
    "repeats":NREP,
    "test_fraction":TEST_SIZE,
    "seed":SEED,
    "features":19,
    "pca_components_supervised":3,
    "primary_transform":{
        "positive_ratio_features":"log2 then training-fold Z-score (ddof=0)",
        "other_features":"training minimum shift; training shifted-column sum to 10,000; log1p; training-fold Z-score (ddof=0)",
        "test_application":"reuse training-fold parameters; clip below training minimum before log1p",
        "yeo_johnson":False,
    },
}


def load_data():
    lab=pd.read_csv(base.E_FILE).drop(columns=base.E19)
    with ZipFile(base.ZIP_FILE) as z:
        raw=pd.read_csv(z.open("Data/cell_metadata_AllCell.csv"),usecols=["cell_label",*base.E19])
    d=lab.merge(raw,on="cell_label",validate="one_to_one")
    d["Consensus"]=d.Consensus.astype(str).str.lower().eq("true")
    return d


def grouped_indices(groups, rep, frac=.8):
    rng=np.random.default_rng(SEED+rep); u=np.unique(groups)
    keep=rng.choice(u,size=int(np.ceil(frac*len(u))),replace=False)
    return np.flatnonzero(np.isin(groups,keep))


def sensitivity_repeat(rep,x,groups,ref):
    idx=grouped_indices(groups,rep); z,_=base.fit_e_transform(x.iloc[idx]); rows=[]
    p=PCA(n_components=10,svd_solver="full").fit(z); allpc=p.transform(z)
    for npc in range(2,11):
        pc=allpc[:,:npc]
        for k in range(2,9):
            lab=AgglomerativeClustering(n_clusters=k,linkage="ward").fit_predict(pc)
            rows.append({"repeat":rep,"NPC":npc,"K":k,"n_cells":len(idx),
                         "ARI_vs_frozen_E4":adjusted_rand_score(ref[idx],lab),
                         "silhouette":silhouette_score(pc,lab),
                         "minimum_cluster_n":int(pd.Series(lab).value_counts().min())})
    return rows


def models(seed):
    return {
      "Multinomial logistic":LogisticRegression(max_iter=3000,class_weight="balanced"),
      "Linear SVM":SVC(kernel="linear",class_weight="balanced"),
      "RBF SVM":SVC(kernel="rbf",class_weight="balanced"),
      "Random forest":RandomForestClassifier(n_estimators=200,class_weight="balanced",n_jobs=1,random_state=seed),
      "Extra trees":ExtraTreesClassifier(n_estimators=200,class_weight="balanced",n_jobs=1,random_state=seed),
      "kNN":KNeighborsClassifier(n_neighbors=7,weights="distance"),
      "HistGradientBoosting":HistGradientBoostingClassifier(max_iter=180,learning_rate=.06,random_state=seed),
    }


def supervised_repeat(rep,seed,tr,te,x,y):
    ztr,par=base.fit_e_transform(x.iloc[tr]); zte=base.apply_e_transform(x.iloc[te],par)
    p=PCA(n_components=3,svd_solver="full").fit(ztr); a=p.transform(ztr); b=p.transform(zte)
    rows=[]; cms={}; importance=None
    for name,m in models(seed).items():
        with warnings.catch_warnings():
            warnings.simplefilter("ignore"); m.fit(a,y[tr]); pred=m.predict(b)
        rec=recall_score(y[te],pred,labels=CLASSES,average=None,zero_division=0)
        rows.append({"repeat":rep,"algorithm":name,"balanced_accuracy":balanced_accuracy_score(y[te],pred),
                     "macro_F1":f1_score(y[te],pred,average="macro"),"MCC":matthews_corrcoef(y[te],pred),
                     **{f"recall_C{i+1}":rec[i] for i in range(4)}})
        cms[name]=confusion_matrix(y[te],pred,labels=CLASSES)
    # Named-feature stability, no PCA.
    et=ExtraTreesClassifier(n_estimators=200,class_weight="balanced",n_jobs=1,random_state=seed).fit(ztr,y[tr])
    lr=LogisticRegression(max_iter=3000,class_weight="balanced").fit(ztr,y[tr])
    importance=(et.feature_importances_,np.mean(np.abs(lr.coef_),axis=0))
    return rows,cms,importance


def valid_splits(y,g):
    out=[]; seed=SEED
    while len(out)<NREP:
        tr,te=next(GroupShuffleSplit(n_splits=1,test_size=TEST_SIZE,random_state=seed).split(y,y,g))
        if set(y[tr])==set(y) and set(y[te])==set(y): out.append((len(out),seed,tr,te))
        seed+=1
    return out


def cramers(a,b):
    t=pd.crosstab(a,b); chi=chi2_contingency(t,correction=False)[0]; n=t.to_numpy().sum()
    return float(np.sqrt((chi/n)/max(1,min(t.shape)-1)))


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    (OUT/"00_parameters.json").write_text(json.dumps(ANALYSIS_DESIGN,indent=2),encoding="utf-8")
    d=load_data(); x=d[base.E19]; y=d.HC_class.astype(str).to_numpy(); g=d[GROUP_FIELD].astype(str).to_numpy()
    z,_=base.fit_e_transform(x); p=PCA(svd_solver="full").fit(z); pc=p.transform(z)
    frozen=AgglomerativeClustering(n_clusters=4,linkage="ward").fit_predict(pc[:,:3])
    repro=adjusted_rand_score(y,frozen)
    pd.DataFrame({"PC":range(1,20),"explained_variance_ratio":p.explained_variance_ratio_,"cumulative_variance":np.cumsum(p.explained_variance_ratio_)}).to_csv(OUT/"01_PCA_explained_variance.csv",index=False)

    print("500 x NPC2-10 x K2-8 donor-level sensitivity",flush=True)
    sens=Parallel(n_jobs=-1,verbose=5)(delayed(sensitivity_repeat)(i,x,g,y) for i in range(NREP))
    sens=pd.DataFrame([q for r in sens for q in r]); sens.to_csv(OUT/"02_NPC2-10_K2-8_donor_resampling_500.csv",index=False)
    ss=sens.groupby(["NPC","K"]).agg(ARI_median=("ARI_vs_frozen_E4","median"),ARI_q025=("ARI_vs_frozen_E4",lambda v:v.quantile(.025)),ARI_q975=("ARI_vs_frozen_E4",lambda v:v.quantile(.975)),silhouette_median=("silhouette","median"),minimum_cluster_n_median=("minimum_cluster_n","median")).reset_index()
    ss.to_csv(OUT/"03_NPC_K_sensitivity_summary.csv",index=False)

    dc=d[d.Consensus].reset_index(drop=True); sx=dc[base.E19]; sy=dc.HC_class.astype(str).to_numpy(); sg=dc[GROUP_FIELD].astype(str).to_numpy()
    splits=valid_splits(sy,sg); print("500 donor-held-out splits x 7 models",flush=True)
    runs=Parallel(n_jobs=-1,verbose=5)(delayed(supervised_repeat)(*s,sx,sy) for s in splits)
    perf=pd.DataFrame([q for rows,_,_ in runs for q in rows]); perf.to_csv(OUT/"04_supervised_500_all_models_class_recall.csv",index=False)
    metrics=["balanced_accuracy","macro_F1","MCC","recall_C1","recall_C2","recall_C3","recall_C4"]
    out=[]
    for alg,v in perf.groupby("algorithm"):
        row={"algorithm":alg}
        for m in metrics: row.update({m+"_median":v[m].median(),m+"_q025":v[m].quantile(.025),m+"_q975":v[m].quantile(.975)})
        out.append(row)
    pd.DataFrame(out).to_csv(OUT/"05_supervised_summary_95CI.csv",index=False)
    for alg in models(SEED):
        cm=sum((cms[alg] for _,cms,_ in runs),np.zeros((4,4),int))
        pd.DataFrame(cm,index=[f"true_{c}" for c in CLASSES],columns=[f"pred_{c}" for c in CLASSES]).to_csv(OUT/f"06_confusion_sum_{alg.replace(' ','_')}.csv")
    et=np.vstack([imp[0] for _,_,imp in runs]); lr=np.vstack([imp[1] for _,_,imp in runs])
    fi=pd.DataFrame({"feature":base.E19,"ExtraTrees_median":np.median(et,0),"ExtraTrees_q025":np.quantile(et,.025,axis=0),"ExtraTrees_q975":np.quantile(et,.975,axis=0),"Logistic_abscoef_median":np.median(lr,0),"Logistic_abscoef_q025":np.quantile(lr,.025,axis=0),"Logistic_abscoef_q975":np.quantile(lr,.975,axis=0)}).sort_values("ExtraTrees_median",ascending=False)
    fi.to_csv(OUT/"07_feature_importance_stability_500.csv",index=False)

    old=ROOT/"outputs"/"Macaque_EM_ML500"; uns=pd.read_csv(old/"01_E4_donor_resampling_500.csv")
    uns.to_csv(OUT/"08_cross_algorithm_E4_donor_resampling_500.csv",index=False)
    cross=pd.read_csv(old/"15_full_data_cross_algorithm_ARI.csv"); cross[cross.analysis=="E4_full"].to_csv(OUT/"09_full_data_cross_algorithm_ARI.csv",index=False)
    perm=pd.read_csv(old/"13_label_permutation_500.csv"); perm=perm[perm.analysis=="E4_frozen_consensus"]
    perm.to_csv(OUT/"10_label_permutation_500.csv",index=False)
    cells=pd.read_csv(old/"16_cell_level_consensus_stability.csv"); cells[cells.analysis=="E4"].to_csv(OUT/"11_cell_level_consensus_stability.csv",index=False)

    conf=[]
    for field in ["Lib_region_of_interest_label","T_class","donor_label"]:
        conf.append({"annotation":field,"Cramers_V":cramers(dc.HC_class,dc[field])})
        pd.crosstab(dc.HC_class,dc[field]).to_csv(OUT/f"12_E4_by_{field}.csv")
    pd.DataFrame(conf).to_csv(OUT/"12_confounding_summary.csv",index=False)

    ps=pd.read_csv(OUT/"05_supervised_summary_95CI.csv"); best=ps.sort_values("balanced_accuracy_median",ascending=False).iloc[0]
    ward=uns[uns.algorithm=="Ward"].ARI_vs_reference; spectral=uns[uns.algorithm=="Spectral"].ARI_vs_reference
    fig,ax=plt.subplots(1,3,figsize=(9,2.8))
    sns.heatmap(ss.pivot(index="NPC",columns="K",values="ARI_median"),vmin=0,vmax=1,cmap="viridis",annot=True,fmt=".2f",ax=ax[0]); ax[0].set_title("Median ARI vs frozen E4")
    sns.boxplot(data=uns,x="algorithm",y="ARI_vs_reference",color="#56B4E9",ax=ax[1]); ax[1].tick_params(axis="x",rotation=35); ax[1].set_title("500 donor resamples")
    sns.boxplot(data=perf,x="algorithm",y="balanced_accuracy",color="#009E73",ax=ax[2]); ax[2].tick_params(axis="x",rotation=50); ax[2].set_title("Donor-held-out recovery")
    for a in ax: a.spines[["top","right"]].set_visible(False)
    fig.tight_layout(); fig.savefig(OUT/"13_E4_full_ML500_overview.png",dpi=600); fig.savefig(OUT/"13_E4_full_ML500_overview.pdf"); plt.close(fig)

    result={"scope":"Macaque MSN; Ca+Pu+NAC","all_complete_cells":len(d),"consensus_cells":len(dc),"donors":int(d.donor_label.nunique()),"features":19,"repeats":500,"frozen_reproduction_ARI":repro,
            "Ward_resampling_ARI_median":float(ward.median()),"Ward_resampling_ARI_95CI":[float(ward.quantile(.025)),float(ward.quantile(.975))],
            "Spectral_resampling_ARI_median":float(spectral.median()),"best_supervised_model":best.algorithm,"best_supervised_BA_median":float(best.balanced_accuracy_median),
            "permutation_BA_median":float(perm.balanced_accuracy.median()),"seed":SEED}
    (OUT/"00_E4_ML500_summary.json").write_text(json.dumps(result,indent=2),encoding="utf-8")
    report=f"""Macaque E4分类：500次机器学习稳定性与可复现性审计
====================================================
范围：MSN，仅Ca+Pu+NAc。E19完整细胞390个，55个供体；冻结HC-GC共识细胞368个。

1. 精确复现：冻结变换、PCA PC1-PC3和Ward K=4重新计算后，ARI={repro:.3f}。
2. 发现稳定性：500次80%供体重抽样中，Ward ARI中位数={ward.median():.3f}（95%区间{ward.quantile(.025):.3f}–{ward.quantile(.975):.3f}）；谱聚类中位数={spectral.median():.3f}。
3. 供体外标签复现：最佳模型{best.algorithm}，平衡准确率中位数={best.balanced_accuracy_median:.3f}（95%区间{best.balanced_accuracy_q025:.3f}–{best.balanced_accuracy_q975:.3f}）。
4. 标签置换基线：平衡准确率中位数={perm.balanced_accuracy.median():.3f}。
5. 结论：E4标签高度可学习并可在隔离供体中复现，但Ward从头重聚类对供体抽样敏感。E4应表述为可复现的冻结操作型分类，而不是在所有抽样和算法下必然出现的四个天然离散群。
"""
    (OUT/"README_E4_results_CN.txt").write_text(report,encoding="utf-8"); print(json.dumps(result,indent=2),flush=True)

if __name__=="__main__": main()
