#!/usr/bin/env python
"""Add Extra Trees held-donor probabilities and one-vs-rest ROC/PR data."""
from pathlib import Path
import sys
import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from sklearn.decomposition import PCA
from sklearn.ensemble import ExtraTreesClassifier
from sklearn.metrics import average_precision_score, precision_recall_curve, roc_auc_score, roc_curve

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))
import validate_macaque_E4_full_ml500 as ml

OUT=ROOT/"outputs/Macaque_E4_full_ML500"; CLASSES=["C1","C2","C3","C4"]


def one_split(rep,seed,tr,te,x,y,ids):
    ztr,par=ml.base.fit_e_transform(x.iloc[tr]); zte=ml.base.apply_e_transform(x.iloc[te],par)
    p=PCA(n_components=3,svd_solver="full").fit(ztr)
    model=ExtraTreesClassifier(n_estimators=200,class_weight="balanced",n_jobs=1,random_state=seed)
    model.fit(p.transform(ztr),y[tr]); raw=model.predict_proba(p.transform(zte))
    proba=np.zeros((len(te),4)); proba[:,[CLASSES.index(c) for c in model.classes_]]=raw
    return pd.DataFrame({"repeat":rep,"cell_label":ids[te],"true_class":y[te],
                         **{f"p_{c}":proba[:,i] for i,c in enumerate(CLASSES)}})


def main():
    d=ml.load_data(); d=d[d.Consensus].reset_index(drop=True)
    x=d[ml.base.E19]; y=d.HC_class.astype(str).to_numpy(); g=d.donor_label.astype(str).to_numpy(); ids=d.cell_label.astype(str).to_numpy()
    rows=Parallel(n_jobs=-1,verbose=5)(delayed(one_split)(*s,x,y,ids) for s in ml.valid_splits(y,g))
    held=pd.concat(rows,ignore_index=True); held.to_csv(OUT/"14_ExtraTrees_heldout_probabilities_500.csv",index=False)
    agg=held.groupby(["cell_label","true_class"],as_index=False).agg(
        n_held_out=("repeat","size"),**{f"p_{c}":(f"p_{c}","mean") for c in CLASSES})
    agg.to_csv(OUT/"15_ExtraTrees_cell_mean_heldout_probabilities.csv",index=False)
    roc_rows=[]; pr_rows=[]; metrics=[]
    for c in CLASSES:
        truth=agg.true_class.eq(c).astype(int); score=agg[f"p_{c}"]
        fpr,tpr,thr=roc_curve(truth,score); precision,recall,pr_thr=precision_recall_curve(truth,score)
        roc_rows.extend({"Class":c,"FPR":a,"TPR":b,"Threshold":t} for a,b,t in zip(fpr,tpr,thr))
        pr_rows.extend({"Class":c,"Recall":a,"Precision":b,"Threshold":t if i<len(pr_thr) else np.nan}
                       for i,(a,b,t) in enumerate(zip(recall,precision,np.r_[pr_thr,np.nan])))
        metrics.append({"Class":c,"ROC_AUC":roc_auc_score(truth,score),"Average_precision":average_precision_score(truth,score),"Prevalence":truth.mean()})
    pd.DataFrame(roc_rows).to_csv(OUT/"16_ExtraTrees_ROC_curve.csv",index=False)
    pd.DataFrame(pr_rows).to_csv(OUT/"17_ExtraTrees_PR_curve.csv",index=False)
    pd.DataFrame(metrics).to_csv(OUT/"18_ExtraTrees_ROC_PR_metrics.csv",index=False)
    assert len(agg)==368 and agg.n_held_out.min()>0 and np.allclose(agg[[f"p_{c}" for c in CLASSES]].sum(1),1)
    print(pd.DataFrame(metrics).to_string(index=False)); print(f"Held-out rows={len(held)}; cells={len(agg)}")


if __name__=="__main__": main()
