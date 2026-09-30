from __future__ import annotations

import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
import numpy as np
import pandas as pd


ROOT = Path(r"C:\Users\53461\OneDrive\Documentos\Patch-seq_Mouse_Acb_MSN_T-type_Visualization")
ML = ROOT / "outputs/morph_qc/final_morph187_NPC3_HCK4_GCres0p50_ML_validation"
MORPH = ROOT / "outputs/morph_qc/final_morph187_NPC3_HCK4_GCres0p50_figures"
OUT = ML / "06_ML_workflow_and_multilevel_validation"

CLASSES = ["M1", "M2", "M3", "M4"]
COLORS = {"M1": "#00468B", "M2": "#42B540", "M3": "#ED0000", "M4": "#0099B4"}


def configure() -> None:
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "DejaVu Sans"],
        "font.size": 7,
        "axes.titlesize": 7,
        "axes.labelsize": 7,
        "xtick.labelsize": 6.5,
        "ytick.labelsize": 6.5,
        "axes.linewidth": .7,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
    })


def save(fig: plt.Figure, stem: str) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT / f"{stem}.png", dpi=900, bbox_inches="tight", pad_inches=.03, facecolor="white")
    fig.savefig(OUT / f"{stem}.pdf", bbox_inches="tight", pad_inches=.03, facecolor="white")
    fig.savefig(OUT / f"{stem}.svg", bbox_inches="tight", pad_inches=.03, facecolor="white")
    plt.close(fig)


def add_workflow(ax: plt.Axes) -> None:
    ax.set_axis_off()
    boxes = [
        ("Frozen input", "187 cells\n10 de-redundant metrics"),
        ("Reference labels", "PCA NPC3 → Ward.D2\nHC K=4 (M1–M4)"),
        ("Leakage-safe transform", "Fit in training fold only\nshift → ×10,000 → log1p → Z"),
        ("Nested grouped CV", "Outer: 5 repeats × 5 folds\nInner: 4 folds; group=recording day"),
        ("Model benchmark", "6 model families\nRBF SVM selected"),
        ("Repeated OOF output", "Each cell tested in every repeat\naggregate probabilities + votes"),
    ]
    positions = [(.01, .57), (.355, .57), (.70, .57), (.70, .10), (.355, .10), (.01, .10)]
    width, height = .285, .32
    for index, ((heading, body), (x, y)) in enumerate(zip(boxes, positions)):
        face = "#EEF4F8" if index not in (1, 4) else ("#F3F3F3" if index == 1 else "#E6F2FA")
        patch = FancyBboxPatch(
            (x, y), width, height,
            boxstyle="round,pad=0.012,rounding_size=.025",
            facecolor=face, edgecolor="#444444", linewidth=.65,
        )
        ax.add_patch(patch)
        ax.text(x + width / 2, y + .235, f"{index + 1}  {heading}", ha="center", va="center",
                fontsize=6.3, fontweight="bold")
        ax.text(x + width / 2, y + .105, body, ha="center", va="center", fontsize=5.35, linespacing=1.15)
    arrow_style = dict(arrowstyle="-|>", mutation_scale=7, linewidth=.7, color="#555555")
    ax.add_patch(FancyArrowPatch((.305, .73), (.345, .73), **arrow_style))
    ax.add_patch(FancyArrowPatch((.65, .73), (.69, .73), **arrow_style))
    ax.add_patch(FancyArrowPatch((.842, .56), (.842, .43), **arrow_style))
    ax.add_patch(FancyArrowPatch((.70, .26), (.65, .26), **arrow_style))
    ax.add_patch(FancyArrowPatch((.355, .26), (.305, .26), **arrow_style))
    ax.text(-.012, .96, "a", fontsize=9, fontweight="bold", ha="left", va="top")


def add_confusion(ax: plt.Axes, summary: dict) -> None:
    cm = pd.read_csv(ML / "05_best_RBF_OOF_confusion_counts.csv", index_col=0).loc[CLASSES, CLASSES].to_numpy(int)
    frac = cm / cm.sum(axis=1, keepdims=True)
    image = ax.imshow(frac, cmap="Blues", vmin=0, vmax=1, interpolation="nearest")
    for row in range(4):
        for col in range(4):
            color = "white" if frac[row, col] > .58 else "#222222"
            ax.text(col, row, f"{cm[row, col]}\n{100*frac[row, col]:.0f}%",
                    ha="center", va="center", fontsize=6, color=color)
    ax.set_xticks(range(4), CLASSES)
    ax.set_yticks(range(4), CLASSES)
    ax.set_xlabel("Repeated OOF prediction")
    ax.set_ylabel("HC-derived class")
    ax.tick_params(length=0, pad=1)
    ax.set_title("RBF-SVM grouped OOF confusion", pad=3)
    ax.text(-.18, 1.10, "b", transform=ax.transAxes, fontsize=9, fontweight="bold")
    metrics = summary["best_RBF_nested_repeated_OOF_metrics"]
    ci = summary["best_RBF_recording_day_group_bootstrap_95CI"]
    ax.text(
        .5, -0.30,
        f"Accuracy {metrics['accuracy']:.3f} ({summary['best_RBF_correct_cells']}/{summary['best_RBF_total_cells']})  |  "
        f"BA {metrics['balanced_accuracy']:.3f} [95% CI {ci['balanced_accuracy'][0]:.3f}–{ci['balanced_accuracy'][1]:.3f}]\n"
        f"Macro-F1 {metrics['macro_F1']:.3f}  |  ARI {metrics['ARI']:.3f}  |  mean vote agreement {summary['best_RBF_mean_vote_agreement']:.3f}",
        transform=ax.transAxes, ha="center", va="top", fontsize=5.7, linespacing=1.25,
    )
    return image


def add_benchmark(ax: plt.Axes) -> None:
    family = pd.read_csv(ML / "01_candidate_family_summary.csv").sort_values("mean_balanced_accuracy")
    y = np.arange(len(family))
    colors = np.where(family["Family"].eq("RBF SVM"), "#00468B", "#A7A7A7")
    ax.barh(y, family["mean_balanced_accuracy"], xerr=family["sd_balanced_accuracy"],
            color=colors, height=.58, edgecolor="none", error_kw={"lw": .55, "capsize": 1.5})
    for yi, value in zip(y, family["mean_balanced_accuracy"]):
        ax.text(value + .008, yi, f"{value:.3f}", va="center", fontsize=5.7)
    ax.set_yticks(y, [""] * len(y))
    ax.set_xlim(.60, 1.01)
    ax.set_xticks([.6, .7, .8, .9, 1.0])
    ax.set_xlabel("Outer-fold balanced accuracy")
    ax.set_title("Candidate-model benchmark", pad=3)
    ax.grid(axis="x", color="#E2E2E2", lw=.4)
    ax.set_axisbelow(True)
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.tick_params(axis="y", length=0)
    for yi, (name, value) in enumerate(zip(family["Family"], family["mean_balanced_accuracy"])):
        ax.text(.607, yi, name, va="center", ha="left", fontsize=5.6,
                color="white" if name == "RBF SVM" else "#222222")
    ax.text(-.12, 1.10, "c", transform=ax.transAxes, fontsize=9, fontweight="bold")


def add_multilevel_stability(ax: plt.Axes) -> None:
    per_class = pd.read_csv(ML / "05_best_RBF_per_class_metrics.csv").set_index("Class").loc[CLASSES]
    boot = pd.read_csv(MORPH / "19_cellwise_bootstrap_stability_by_class_summary_500x.csv").set_index("M_class").loc[CLASSES]
    assignments = pd.read_csv(MORPH / "00_final_morph187_cell_assignments.csv")
    hc_gc = assignments["HC_GC_consensus"].astype(str).str.lower().eq("true").mean()

    x = np.arange(4)
    width = .34
    ax.bar(x - width / 2, per_class["Recall"], width, color="#00468B", label="ML OOF recall")
    ax.bar(x + width / 2, boot["Median_recovery"], width, color="#C9A227", label="HC bootstrap median")
    ax.axhline(hc_gc, color="#ED0000", lw=.8, ls=(0, (3, 2)), label=f"HC–GC agreement ({hc_gc:.3f})")
    for xi, value in zip(x - width / 2, per_class["Recall"]):
        ax.text(xi, value + .025, f"{value:.2f}", ha="center", fontsize=5.5)
    for xi, value in zip(x + width / 2, boot["Median_recovery"]):
        ax.text(xi, value + .025, f"{value:.2f}", ha="center", fontsize=5.5)
    ax.set_xticks(x, CLASSES)
    ax.set_ylim(0, 1.10)
    ax.set_yticks(np.arange(0, 1.01, .2))
    ax.set_ylabel("Recovery / agreement")
    ax.set_title("Class recovery and resampling stability", pad=3)
    ax.grid(axis="y", color="#E2E2E2", lw=.4)
    ax.set_axisbelow(True)
    ax.spines[["top", "right"]].set_visible(False)
    legend = ax.legend(
        frameon=True,
        fancybox=False,
        framealpha=1.0,
        facecolor="white",
        edgecolor="#222222",
        fontsize=5.3,
        loc="lower right",
        handlelength=1.5,
        borderpad=.45,
    )
    legend.get_frame().set_linewidth(.55)
    ax.text(-.12, 1.10, "d", transform=ax.transAxes, fontsize=9, fontweight="bold")


def main() -> None:
    configure()
    summary = json.loads((ML / "00_ML_validation_summary.json").read_text(encoding="utf-8"))
    fig = plt.figure(figsize=(7.2, 6.0), facecolor="white")
    ax_flow = fig.add_axes([.03, .68, .94, .29])
    ax_cm = fig.add_axes([.055, .20, .28, .38])
    ax_bench = fig.add_axes([.405, .22, .245, .34])
    ax_stab = fig.add_axes([.705, .22, .275, .34])
    add_workflow(ax_flow)
    add_confusion(ax_cm, summary)
    add_benchmark(ax_bench)
    add_multilevel_stability(ax_stab)
    fig.text(
        .5, .020,
        "Interpretation: ML quantifies recoverability of HC-derived labels from the same morphology features; "
        "HC–GC agreement and cell-resampling add algorithmic and stability evidence, but external biological validation remains independent.",
        ha="center", va="bottom", fontsize=5.4,
    )
    save(fig, "ML_workflow_and_multilevel_validation_final187")


if __name__ == "__main__":
    main()
