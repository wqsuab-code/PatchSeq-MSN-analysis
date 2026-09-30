#!/usr/bin/env python3
"""Complete the leakage-safe, submission-facing Macaque E4 ML checks."""

from __future__ import annotations

from pathlib import Path
import json
import platform
import sys

import joblib
from joblib import Parallel, delayed
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
from scipy.stats import chi2_contingency
import seaborn as sns
import sklearn
from sklearn.decomposition import PCA
from sklearn.ensemble import ExtraTreesClassifier, RandomForestClassifier
from sklearn.inspection import permutation_importance
from sklearn.metrics import balanced_accuracy_score, f1_score
from sklearn.model_selection import StratifiedGroupKFold

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import validate_macaque_E4_full_ml500 as ml

OUT = ROOT / "outputs" / "Macaque_E4_ML_submission_complete"
TABLES = OUT / "tables"
SEED = 20260912
N_PERM = 500
N_META_PERM = 5000
CLASSES = ["C1", "C2", "C3", "C4"]
SUBSETS = ["All 19", "Remove top 2", "Remove top 5", "Bottom 3"]


def et(seed):
    return ExtraTreesClassifier(
        n_estimators=200, class_weight="balanced", n_jobs=1, random_state=seed
    )


def rf(seed):
    return RandomForestClassifier(
        n_estimators=500, max_features="sqrt", min_samples_leaf=1,
        class_weight="balanced_subsample", n_jobs=1, random_state=seed,
    )


def valid_sgkf(y, groups, n_splits, seed):
    all_classes = set(np.unique(y))
    for candidate in range(seed, seed + 1000):
        folds = list(StratifiedGroupKFold(
            n_splits=n_splits, shuffle=True, random_state=candidate
        ).split(np.zeros(len(y)), y, groups))
        if all(set(y[tr]) == all_classes and set(y[te]) == all_classes for tr, te in folds):
            return candidate, folds
    raise RuntimeError("Could not construct grouped folds containing all classes")


def transform_fold(x, train, test):
    ztr, params = ml.base.fit_e_transform(x.iloc[train])
    return ztr, ml.base.apply_e_transform(x.iloc[test], params)


def within_donor_permutation(y, groups, seed):
    rng = np.random.default_rng(seed)
    out = np.asarray(y).copy()
    for donor in np.unique(groups):
        idx = np.flatnonzero(groups == donor)
        out[idx] = rng.permutation(out[idx])
    return out


def et_cv_mean(x, y, folds, seed):
    scores = []
    for fold, (tr, te) in enumerate(folds):
        ztr, zte = transform_fold(x, tr, te)
        pca = PCA(n_components=3, svd_solver="full").fit(ztr)
        model = et(seed + fold).fit(pca.transform(ztr), y[tr])
        scores.append(balanced_accuracy_score(y[te], model.predict(pca.transform(zte))))
    return scores


def one_null(rep, x, y, groups, folds):
    yp = within_donor_permutation(y, groups, SEED + 10000 + rep)
    scores = et_cv_mean(x, yp, folds, SEED + 20000 + rep * 10)
    return {"permutation": rep + 1, "mean_balanced_accuracy": np.mean(scores),
            **{f"fold_{i + 1}": v for i, v in enumerate(scores)}}


def one_importance(task, x, y):
    repeat, fold, tr, te = task
    ztr, zte = transform_fold(x, tr, te)
    model = rf(SEED + repeat * 100 + fold).fit(ztr, y[tr])
    result = permutation_importance(
        model, zte, y[te], scoring="balanced_accuracy", n_repeats=20,
        random_state=SEED + 1000 * repeat + fold, n_jobs=1,
    )
    rows = []
    for j, feature in enumerate(x.columns):
        for permutation, value in enumerate(result.importances[j], 1):
            rows.append({"repeat": repeat, "fold": fold, "feature": feature,
                         "permutation": permutation, "importance": value})
    return rows


def one_nested_ablation(task, x, y, groups):
    repeat, fold, tr, te = task
    inner_seed, inner = valid_sgkf(y[tr], groups[tr], 3, SEED + 30000 + repeat * 100 + fold)
    inner_values = []
    for inner_fold, (itr0, iva0) in enumerate(inner):
        itr, iva = tr[itr0], tr[iva0]
        ztr, zva = transform_fold(x, itr, iva)
        model = rf(inner_seed + inner_fold).fit(ztr, y[itr])
        pi = permutation_importance(
            model, zva, y[iva], scoring="balanced_accuracy", n_repeats=10,
            random_state=inner_seed + 100 + inner_fold, n_jobs=1,
        )
        inner_values.append(pi.importances_mean)
    rank = np.argsort(np.mean(inner_values, axis=0))[::-1]
    definitions = {
        "All 19": np.arange(x.shape[1]),
        "Remove top 2": rank[2:],
        "Remove top 5": rank[5:],
        "Bottom 3": rank[-3:],
    }
    ztr, zte = transform_fold(x, tr, te)
    performance = []
    for i, label in enumerate(SUBSETS):
        cols = definitions[label]
        model = rf(SEED + 40000 + repeat * 1000 + fold * 10 + i).fit(ztr[:, cols], y[tr])
        pred = model.predict(zte[:, cols])
        performance.append({
            "repeat": repeat, "fold": fold, "subset": label,
            "n_features": len(cols),
            "balanced_accuracy": balanced_accuracy_score(y[te], pred),
            "macro_F1": f1_score(y[te], pred, average="macro"),
        })
    ranking = pd.DataFrame({
        "repeat": repeat, "fold": fold, "rank": np.arange(1, x.shape[1] + 1),
        "feature": np.asarray(x.columns)[rank], "inner_seed": inner_seed,
    })
    return performance, ranking


def cramers_v(a, b, corrected):
    table = pd.crosstab(a, b)
    chi2 = chi2_contingency(table, correction=False)[0]
    n = table.to_numpy().sum()
    phi2 = chi2 / n
    r, k = table.shape
    if not corrected:
        return float(np.sqrt(phi2 / max(1, min(r - 1, k - 1))))
    phi2 = max(0.0, phi2 - (k - 1) * (r - 1) / (n - 1))
    rc = r - (r - 1) ** 2 / (n - 1)
    kc = k - (k - 1) ** 2 / (n - 1)
    return float(np.sqrt(phi2 / max(np.finfo(float).eps, min(rc - 1, kc - 1))))


def metadata_test(labels, values, donors, name):
    observed = cramers_v(labels, values, corrected=True)
    rng = np.random.default_rng(SEED + 50000 + sum(map(ord, name)))
    null = np.empty(N_META_PERM)
    for i in range(N_META_PERM):
        if name == "donor_label":
            shuffled = rng.permutation(labels)
        else:
            shuffled = labels.copy()
            for donor in np.unique(donors):
                idx = np.flatnonzero(donors == donor)
                shuffled[idx] = rng.permutation(shuffled[idx])
        null[i] = cramers_v(shuffled, values, corrected=True)
    return {
        "metadata": name,
        "cramers_v_uncorrected": cramers_v(labels, values, corrected=False),
        "cramers_v_bias_corrected": observed,
        "permutation_scheme": "global E-label shuffle" if name == "donor_label" else "within-donor E-label shuffle",
        "n_permutations": N_META_PERM,
        "null_median": float(np.median(null)),
        "null_q025": float(np.quantile(null, .025)),
        "null_q975": float(np.quantile(null, .975)),
        "empirical_p": float((np.sum(null >= observed) + 1) / (N_META_PERM + 1)),
    }


def save_figure(permutation, observed, importance, ablation, metadata):
    mpl.rcParams.update({"font.family": "Arial", "font.size": 6, "pdf.fonttype": 42,
                         "axes.linewidth": .55, "savefig.facecolor": "white"})
    fig, axes = plt.subplots(1, 4, figsize=(7.205, 1.9), gridspec_kw={"wspace": .58})
    sns.histplot(permutation.mean_balanced_accuracy, bins=24, color="#BDBDBD",
                 edgecolor="white", linewidth=.3, ax=axes[0])
    axes[0].axvline(observed, color="#D55E00", lw=1)
    axes[0].set(xlabel="Mean balanced accuracy", title="Matched ET permutation")

    q = importance.sort_values("relative_importance").tail(10)
    axes[1].barh(q.feature.str.removeprefix("Epsy_"), q.relative_importance, color="#4C78A8")
    axes[1].set(xlabel="Relative importance", title="Train-only permutation importance")
    axes[1].tick_params(axis="y", labelsize=3.8)

    sns.boxplot(data=ablation, x="subset", y="balanced_accuracy", order=SUBSETS,
                color="#009E73", showfliers=False, linewidth=.55, ax=axes[2])
    axes[2].axhline(.25, color="#777", ls=":", lw=.6)
    axes[2].set(xlabel="", ylabel="Balanced accuracy", title="Nested feature ablation")
    axes[2].set_xticks(range(4), ["All 19", "−top 2", "−top 5", "Bottom 3"], rotation=25, ha="right")

    q = metadata.sort_values("cramers_v_bias_corrected")
    axes[3].barh(q.metadata.str.replace("Lib_region_of_interest_label", "ROI"),
                 q.cramers_v_bias_corrected, color="#8C8C8C")
    for i, row in enumerate(q.itertuples()):
        axes[3].text(row.cramers_v_bias_corrected + .01, i, f"P={row.empirical_p:.3g}", va="center", fontsize=5)
    axes[3].set(xlabel="Bias-corrected Cramer's V", title="Metadata association")
    for ax in axes:
        ax.spines[["top", "right"]].set_visible(False)
    fig.savefig(OUT / "Supplementary_ML4_submission_corrections.png", dpi=600, bbox_inches="tight")
    fig.savefig(OUT / "Supplementary_ML4_submission_corrections.pdf", bbox_inches="tight")
    plt.close(fig)


def main():
    TABLES.mkdir(parents=True, exist_ok=True)
    data = ml.load_data()
    data = data[data.Consensus].reset_index(drop=True)
    x = data[ml.base.E19].apply(pd.to_numeric, errors="raise")
    y = data.HC_class.astype(str).to_numpy()
    groups = data.donor_label.astype(str).to_numpy()
    assert len(data) == 368 and x.shape == (368, 19) and set(y) == set(CLASSES)

    fold_seed, folds = valid_sgkf(y, groups, 5, SEED)
    observed_folds = et_cv_mean(x, y, folds, SEED + 1000)
    observed = float(np.mean(observed_folds))
    permutation = pd.DataFrame(Parallel(n_jobs=-1, verbose=5)(
        delayed(one_null)(i, x, y, groups, folds) for i in range(N_PERM)
    ))
    permutation.to_csv(TABLES / "01_matched_ExtraTrees_within_donor_permutation_500.csv", index=False)
    p_perm = float((np.sum(permutation.mean_balanced_accuracy >= observed) + 1) / (N_PERM + 1))
    pd.DataFrame([{"observed_mean_balanced_accuracy": observed,
                   "observed_fold_scores": ";".join(f"{v:.8f}" for v in observed_folds),
                   "null_median": permutation.mean_balanced_accuracy.median(),
                   "null_q025": permutation.mean_balanced_accuracy.quantile(.025),
                   "null_q975": permutation.mean_balanced_accuracy.quantile(.975),
                   "empirical_p": p_perm, "n_permutations": N_PERM,
                   "fold_seed": fold_seed, "scheme": "within-donor E-label permutation"}]).to_csv(
        TABLES / "02_matched_ExtraTrees_permutation_summary.csv", index=False
    )

    tasks = []
    for repeat in range(5):
        _, repeat_folds = valid_sgkf(y, groups, 5, 777 + repeat)
        tasks.extend((repeat + 1, fold + 1, tr, te) for fold, (tr, te) in enumerate(repeat_folds))
    imp_rows = Parallel(n_jobs=-1, verbose=5)(delayed(one_importance)(task, x, y) for task in tasks)
    importance_raw = pd.DataFrame([row for block in imp_rows for row in block])
    importance_raw.to_csv(TABLES / "03_train_only_permutation_importance_raw.csv", index=False)
    importance = importance_raw.groupby("feature").importance.agg(["mean", "std", "median"]).reset_index()
    importance["positive_mean"] = importance["mean"].clip(lower=0)
    importance["relative_importance"] = importance.positive_mean / max(importance.positive_mean.max(), 1e-12)
    importance = importance.sort_values("relative_importance", ascending=False)
    importance.to_csv(TABLES / "04_train_only_permutation_importance_summary.csv", index=False)

    nested = Parallel(n_jobs=-1, verbose=5)(delayed(one_nested_ablation)(task, x, y, groups) for task in tasks)
    ablation = pd.DataFrame([row for perf, _ in nested for row in perf])
    rankings = pd.concat([rank for _, rank in nested], ignore_index=True)
    ablation.to_csv(TABLES / "05_nested_ablation_performance.csv", index=False)
    rankings.to_csv(TABLES / "06_nested_ablation_training_only_rankings.csv", index=False)
    ablation.groupby(["subset", "n_features"]).agg(
        balanced_accuracy_median=("balanced_accuracy", "median"),
        balanced_accuracy_q025=("balanced_accuracy", lambda v: v.quantile(.025)),
        balanced_accuracy_q975=("balanced_accuracy", lambda v: v.quantile(.975)),
        macro_F1_median=("macro_F1", "median"),
    ).reset_index().to_csv(TABLES / "07_nested_ablation_summary.csv", index=False)

    metadata = pd.DataFrame([
        metadata_test(y, data.donor_label.astype(str).to_numpy(), groups, "donor_label"),
        metadata_test(y, data.T_class.astype(str).to_numpy(), groups, "T_class"),
        metadata_test(y, data.Lib_region_of_interest_label.astype(str).to_numpy(), groups,
                      "Lib_region_of_interest_label"),
    ])
    metadata.to_csv(TABLES / "08_metadata_bias_corrected_permutation_tests.csv", index=False)
    for field in ["donor_label", "T_class", "Lib_region_of_interest_label"]:
        pd.crosstab(data.HC_class, data[field]).to_csv(TABLES / f"09_contingency_E4_by_{field}.csv")

    environment = {
        "python": sys.version, "platform": platform.platform(), "numpy": np.__version__,
        "pandas": pd.__version__, "scipy": scipy.__version__, "scikit_learn": sklearn.__version__,
        "joblib": joblib.__version__, "matplotlib": mpl.__version__, "seaborn": sns.__version__,
        "seed": SEED, "permutations": N_PERM, "metadata_permutations": N_META_PERM,
    }
    (OUT / "10_software_and_parameters.json").write_text(json.dumps(environment, indent=2), encoding="utf-8")
    save_figure(permutation, observed, importance, ablation, metadata)

    summary = {
        "observed_ET_grouped5fold_mean_BA": observed,
        "matched_within_donor_null_median": float(permutation.mean_balanced_accuracy.median()),
        "matched_permutation_empirical_p": p_perm,
        "top5_train_only_permutation_features": importance.head(5).feature.tolist(),
        "metadata": metadata.to_dict(orient="records"),
    }
    (OUT / "00_submission_completion_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    assert len(permutation) == N_PERM and len(ablation) == 25 * 4 and len(metadata) == 3
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
