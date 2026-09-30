#!/usr/bin/env python3
"""Run the final Mouse-E Python analyses in one recorded environment."""

from __future__ import annotations

import contextlib
import csv
import hashlib
import importlib.util
import json
import os
import platform
import shutil
import subprocess
import sys
import argparse
from pathlib import Path

import h5py
import joblib
import matplotlib
import numpy as np
import pandas as pd
import scipy
import sklearn


ROOT = Path(__file__).resolve().parents[1]
RELEASE = ROOT / "outputs" / "Mouse_E_final_performed_using_20260924"
FROZEN = ROOT / "outputs" / "Mouse_E_frozen_analysis_bundle_20260923" / "02_frozen_analysis"
PY_OUT = RELEASE / "03_python_analyses"
ENV_OUT = RELEASE / "00_environment"
AUDIT_OUT = RELEASE / "04_audit"


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot load {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def capture_environment() -> None:
    ENV_OUT.mkdir(parents=True, exist_ok=True)
    versions = {
        "analysis_date": "2026-09-24",
        "python": platform.python_version(),
        "python_executable": sys.executable,
        "platform": platform.platform(),
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "scipy": scipy.__version__,
        "scikit_learn": sklearn.__version__,
        "h5py": h5py.__version__,
        "matplotlib": matplotlib.__version__,
        "joblib": joblib.__version__,
    }
    (ENV_OUT / "python_environment.json").write_text(json.dumps(versions, indent=2), encoding="utf-8")
    freeze = subprocess.run(
        [sys.executable, "-m", "pip", "freeze"], check=True, text=True, capture_output=True
    ).stdout
    (ENV_OUT / "python_requirements_frozen.txt").write_text(freeze, encoding="utf-8")


def run_stats() -> None:
    target = PY_OUT / "statistics"
    target.mkdir(parents=True, exist_ok=True)
    env = os.environ.copy(); env["MOUSE_E_STATS_OUT"] = str(target)
    subprocess.run([sys.executable, str(ROOT / "scripts" / "plot_consensus450_core18_Egroup_statistics.py")], check=True, env=env)


def run_nested_svm() -> Path:
    target = PY_OUT / "SVM_nested_base"
    target.mkdir(parents=True, exist_ok=True)
    env = os.environ.copy(); env["MOUSE_E_SVM_BASE_OUT"] = str(target)
    subprocess.run([sys.executable, str(ROOT / "scripts" / "train_hc_gc_consensus_eclass_ml.py")], check=True, env=env)
    return target


def run_svm_support(base: Path) -> None:
    target = PY_OUT / "SVM_nested_support"
    target.mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    env["MOUSE_E_SVM_BASE_OUT"] = str(base)
    env["MOUSE_E_SVM_SUPPORT_OUT"] = str(target)
    subprocess.run([sys.executable, str(ROOT / "scripts" / "draw_svm_supplementary_10_individual_figures.py")], check=True, env=env)


def run_grouped_ml() -> None:
    target = PY_OUT / "ML_date_grouped"
    target.mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    env["MOUSE_E_GROUPED_ML_OUT"] = str(target)
    env["MOUSE_E_GROUPED_ML_SITE_OUT"] = str(target / "mouse_ml_validation.json")
    subprocess.run([sys.executable, str(ROOT / "scripts" / "validate_mouse_E5_macaqueM_style.py")], check=True, env=env)


def compare_csv(
    name: str,
    new_path: Path,
    frozen_path: Path,
    key_columns: list[str] | None = None,
    atol: float = 1e-10,
) -> dict:
    new = pd.read_csv(new_path)
    old = pd.read_csv(frozen_path)
    if key_columns:
        new = new.sort_values(key_columns).reset_index(drop=True)
        old = old.sort_values(key_columns).reset_index(drop=True)
    same_shape = new.shape == old.shape
    same_columns = list(new.columns) == list(old.columns)
    numeric_max = 0.0
    strings_equal = True
    if same_shape and same_columns:
        for column in new.columns:
            if pd.api.types.is_numeric_dtype(new[column]) and pd.api.types.is_numeric_dtype(old[column]):
                a = new[column].to_numpy(float)
                b = old[column].to_numpy(float)
                finite = np.isfinite(a) & np.isfinite(b)
                if finite.any():
                    numeric_max = max(numeric_max, float(np.max(np.abs(a[finite] - b[finite]))))
                if not np.array_equal(np.isnan(a), np.isnan(b)):
                    numeric_max = float("inf")
            else:
                strings_equal &= new[column].fillna("<NA>").astype(str).equals(
                    old[column].fillna("<NA>").astype(str)
                )
    passed = same_shape and same_columns and strings_equal and numeric_max <= atol
    return {
        "analysis": name,
        "new_file": str(new_path.relative_to(RELEASE)).replace("\\", "/"),
        "frozen_file": str(frozen_path.relative_to(ROOT)).replace("\\", "/"),
        "same_shape": same_shape,
        "same_columns": same_columns,
        "categorical_values_exact": strings_equal,
        "maximum_absolute_numeric_difference": numeric_max,
        "tolerance": atol,
        "pass": passed,
    }


def audit_outputs() -> None:
    rows = []
    rows.append(compare_csv(
        "E-group omnibus statistics",
        PY_OUT / "statistics" / "Core18_Egroup_Omnibus_KW.csv",
        FROZEN / "statistics" / "Core18_Egroup_Omnibus_KW.csv",
        ["Feature"], 1e-10,
    ))
    rows.append(compare_csv(
        "E-group pairwise statistics",
        PY_OUT / "statistics" / "Core18_Egroup_Pairwise_Dunn_Holm.csv",
        FROZEN / "statistics" / "Core18_Egroup_Pairwise_Dunn_Holm.csv",
        ["Feature", "Group_1", "Group_2"], 1e-10,
    ))
    rows.append(compare_csv(
        "Nested SVM confusion counts",
        PY_OUT / "SVM_nested_support" / "01a_confusion_counts.csv",
        FROZEN / "SVM_nested" / "01a_confusion_counts.csv",
        None, 0,
    ))
    rows.append(compare_csv(
        "Nested SVM class metrics",
        PY_OUT / "SVM_nested_support" / "04_classwise_metrics.csv",
        FROZEN / "SVM_nested" / "04_classwise_metrics.csv",
        ["Class"], 1e-10,
    ))
    rows.append(compare_csv(
        "Date-grouped model summary",
        PY_OUT / "ML_date_grouped" / "02_model_summary.csv",
        FROZEN / "ML_date_grouped" / "02_model_summary.csv",
        ["model"], 1e-10,
    ))
    rows.append(compare_csv(
        "Date-grouped confusion counts",
        PY_OUT / "ML_date_grouped" / "04_confusion_counts.csv",
        FROZEN / "ML_date_grouped" / "04_confusion_counts.csv",
        None, 0,
    ))

    # Compare strict441 RRR tables produced by the recorded R environment.
    rrr_new = RELEASE / "02_RRR_R"
    rrr_old = FROZEN / "RRR_strict441"
    rows.append(compare_csv(
        "Strict441 RRR cell scores", rrr_new / "RRR_cell_scores_n441.csv",
        rrr_old / "RRR_cell_scores_n441.csv", ["MSN_unique_ID"], 1e-10,
    ))
    rows.append(compare_csv(
        "Strict441 RRR E-feature loadings", rrr_new / "RRR_Efeature_correlation_loadings_n441.csv",
        rrr_old / "RRR_Efeature_correlation_loadings_n441.csv", ["Feature"], 1e-10,
    ))
    rows.append(compare_csv(
        "Strict441 RRR gene loadings", rrr_new / "RRR_gene_correlation_loadings_n441.csv",
        rrr_old / "RRR_gene_correlation_loadings_n441.csv", ["Gene"], 1e-10,
    ))

    core = pd.read_csv(RELEASE / "01_core_taxonomy_R" / "core_taxonomy_concordance.csv")
    for _, row in core.iterrows():
        rows.append({
            "analysis": f"Core taxonomy: {row['check']}",
            "new_file": "01_core_taxonomy_R/core_taxonomy_concordance.csv",
            "frozen_file": "Mouse_E_frozen_analysis_bundle_20260923",
            "same_shape": True,
            "same_columns": True,
            "categorical_values_exact": True,
            "maximum_absolute_numeric_difference": row["observed"],
            "tolerance": "as specified by core audit",
            "pass": bool(row["pass"]),
        })

    AUDIT_OUT.mkdir(parents=True, exist_ok=True)
    audit = pd.DataFrame(rows)
    audit.to_csv(AUDIT_OUT / "performed_using_concordance.csv", index=False)
    summary = {
        "checks": len(audit),
        "passed": int(audit["pass"].sum()),
        "failed": int((~audit["pass"]).sum()),
        "all_pass": bool(audit["pass"].all()),
    }
    (AUDIT_OUT / "performed_using_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    if not summary["all_pass"]:
        print(audit.loc[~audit["pass"]].to_string(index=False))
        raise RuntimeError("At least one performed-using concordance check failed")


def copy_scripts_and_manifest() -> None:
    code = RELEASE / "05_code"
    code.mkdir(parents=True, exist_ok=True)
    names = [
        "reproduce_mouse_E_final_environment_20260924.R",
        "reproduce_mouse_E_python_analyses_20260924.py",
        "plot_mouse_T_E_RRR_strict441.R",
        "plot_consensus450_core18_Egroup_statistics.py",
        "train_hc_gc_consensus_eclass_ml.py",
        "draw_svm_supplementary_10_individual_figures.py",
        "validate_mouse_E5_macaqueM_style.py",
        "plot_all25_violin_plus_feature_tsne_posthoc.py",
        "export_ephys_statistics_audit.py",
    ]
    for name in names:
        shutil.copy2(ROOT / "scripts" / name, code / name)

    records = []
    for path in sorted(p for p in RELEASE.rglob("*") if p.is_file()):
        records.append({
            "file": path.relative_to(RELEASE).as_posix(),
            "bytes": path.stat().st_size,
            "sha256": sha256(path),
        })
    with (AUDIT_OUT / "file_manifest.csv").open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["file", "bytes", "sha256"])
        writer.writeheader()
        writer.writerows(records)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--audit-only", action="store_true")
    args = parser.parse_args()
    PY_OUT.mkdir(parents=True, exist_ok=True)
    capture_environment()
    if not args.audit_only:
        run_stats()
        base = run_nested_svm()
        run_svm_support(base)
        run_grouped_ml()
    audit_outputs()
    copy_scripts_and_manifest()
    print((AUDIT_OUT / "performed_using_summary.json").read_text(encoding="utf-8"))


if __name__ == "__main__":
    main()
