"""Publication radar redraw from six predeclared features in the frozen 450-cell Z-score matrix."""
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from PIL import Image
from scipy.interpolate import CubicSpline


ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "interactive/data/final18_zscore_long.csv"
OUT = ROOT / "figures/generated/core6_radar"
GROUPS = ["S-1", "S-2", "S-3", "S-4", "S-5"]
DISPLAY = dict(zip(GROUPS, ["E1", "E2", "E3", "E4", "E5"]))
EXPECTED = {"S-1": 142, "S-2": 160, "S-3": 45, "S-4": 61, "S-5": 42}
COLORS = {"S-1": "#4E79A7", "S-2": "#F28E2B", "S-3": "#59A14F", "S-4": "#2AA6B8", "S-5": "#B07AA1"}
CORE_FEATURES = [
    "E_Input.resistance..MOhm.", "E_Membrane.time.constant..ms.", "E_Rheobase..pA.",
    "E_Max.number.of.APs", "E_Afterhyperpolarization..mV.", "E_AP.amplitude.adaptation.index",
]
SHORT = ["Rin", "Tau", "Rheo", "MaxAP", "AHP", "APad"]
Z_MIN, Z_MAX = -3.0, 3.0
CELL_LW, CELL_ALPHA = 0.25, 0.08
MEDIAN_LW, MEDIAN_FILL_ALPHA, IQR_ALPHA = 0.90, 0.18, 0.12
GRID_LW, GRID_ALPHA, FONT_SIZE = 0.35, 0.42, 4.0
INDIVIDUAL_SIZE, COMBINED_SIZE = (1.25, 1.00), (6.25, 1.00)


def load_frozen():
    long = pd.read_csv(SOURCE)
    order = long.sort_values("row_order").drop_duplicates("Feature_full").set_index("Feature_full")
    features = CORE_FEATURES
    labels = order.loc[features, "Feature_label"].tolist()
    meta = long.drop_duplicates("Cell_ID")[["Cell_ID", "E_consensus"]]
    counts = meta["E_consensus"].value_counts().reindex(GROUPS).to_dict()
    if len(meta) != 450 or counts != EXPECTED:
        raise RuntimeError(f"Frozen consensus cohort changed: n={len(meta)}, counts={counts}")
    matrix = long.pivot(index="Cell_ID", columns="Feature_full", values="Z_score").reindex(columns=features)
    if matrix.isna().any().any():
        raise RuntimeError("Missing value in frozen 18-feature Z-score matrix")
    matrix = meta.set_index("Cell_ID").join(matrix, validate="one_to_one")
    return matrix, features, labels


def display_radius(z):
    return np.clip((np.asarray(z, float) - Z_MIN) / (Z_MAX - Z_MIN), 0, 1)


def smooth(values, angles, smooth_angles):
    closed = np.r_[angles, 2 * np.pi]
    return np.clip(CubicSpline(closed, np.r_[values, values[0]], bc_type="periodic")(smooth_angles), 0, 1)


def draw_axis(ax, group, values, angles, smooth_angles):
    color = COLORS[group]
    for row in display_radius(values):
        ax.plot(smooth_angles, smooth(row, angles, smooth_angles), color=color, lw=CELL_LW, alpha=CELL_ALPHA, zorder=1)
    q25, med, q75 = np.nanquantile(values, [0.25, 0.50, 0.75], axis=0)
    q25, med, q75 = (smooth(display_radius(x), angles, smooth_angles) for x in (q25, med, q75))
    ax.fill_between(smooth_angles, q25, q75, color=color, alpha=IQR_ALPHA, lw=0, zorder=2)
    ax.fill(smooth_angles, med, color=color, alpha=MEDIAN_FILL_ALPHA, zorder=3)
    ax.plot(smooth_angles, med, color="#202020", lw=MEDIAN_LW, zorder=4)
    ax.set_theta_offset(np.pi / 2); ax.set_theta_direction(-1); ax.set_ylim(0, 1)
    ax.set_xticks(angles); ax.set_xticklabels([]); ax.tick_params(axis="x", length=0)
    for i, (label, theta, deg) in enumerate(zip(SHORT, angles, np.degrees(angles))):
        rotation = -deg
        if -270 < rotation < -90:
            rotation += 180
        ax.text(theta, 1.10, label, fontsize=FONT_SIZE, color="#303030",
                rotation=rotation, rotation_mode="anchor", ha="center", va="center", clip_on=False)
    ax.set_yticks([1 / 6, 2 / 6, 3 / 6, 4 / 6, 5 / 6, 1]); ax.set_yticklabels([])
    ax.grid(color="#A8A8A8", lw=GRID_LW, alpha=GRID_ALPHA)
    ax.spines["polar"].set(color="#A8A8A8", linewidth=GRID_LW, alpha=GRID_ALPHA)
    ax.set_title(f"{DISPLAY[group]} (n={len(values)})", fontsize=FONT_SIZE, color=color, y=1.30, pad=0)


def make_figure(data, features, groups, size):
    angles = np.linspace(0, 2 * np.pi, len(features), endpoint=False)
    smooth_angles = np.linspace(0, 2 * np.pi, 361)
    fig = plt.figure(figsize=size)
    panel_w = 1 / len(groups)
    axes = []
    for i, group in enumerate(groups):
        ax = fig.add_axes([i * panel_w + 0.075 * panel_w, 0.225, 0.850 * panel_w, 0.535], projection="polar")
        draw_axis(ax, group, data.loc[data["E_consensus"].eq(group), features].to_numpy(float), angles, smooth_angles)
        axes.append(ax)
    return fig, axes


def qa(fig, axes):
    fig.canvas.draw(); renderer = fig.canvas.get_renderer(); canvas = fig.bbox; issues = []
    for i, ax in enumerate(axes, 1):
        boxes = [x.get_window_extent(renderer).expanded(0.98, 0.98) for x in list(ax.texts) + [ax.title]]
        for j, box in enumerate(boxes):
            if box.x0 < canvas.x0 - .5 or box.y0 < canvas.y0 - .5 or box.x1 > canvas.x1 + .5 or box.y1 > canvas.y1 + .5:
                issues.append(f"panel {i}: label {j + 1} outside canvas")
        for a in range(len(boxes)):
            for b in range(a + 1, len(boxes)):
                if boxes[a].overlaps(boxes[b]): issues.append(f"panel {i}: labels {a + 1}/{b + 1} overlap")
    return issues


def save(fig, stem):
    paths = []
    for name, transparent in [("white", False), ("transparent", True)]:
        png = OUT / f"{stem}_{name}_600dpi.png"; pdf = OUT / f"{stem}_{name}_vector.pdf"
        fig.savefig(png, dpi=600, facecolor="none" if transparent else "white", transparent=transparent)
        fig.savefig(pdf, facecolor="none" if transparent else "white", transparent=transparent)
        paths += [png, pdf]
    return paths


def main():
    mpl.rcParams.update({"font.family": "Arial", "font.size": FONT_SIZE, "axes.titlesize": FONT_SIZE,
                         "pdf.fonttype": 42, "ps.fonttype": 42, "savefig.pad_inches": 0})
    OUT.mkdir(parents=True, exist_ok=True)
    data, features, full_labels = load_frozen(); paths = []; qa_log = []
    for group in GROUPS:
        fig, axes = make_figure(data, features, [group], INDIVIDUAL_SIZE)
        issues = qa(fig, axes); qa_log.append(f"{DISPLAY[group]}: " + ("PASS" if not issues else "FAIL; " + "; ".join(issues)))
        paths += save(fig, f"{DISPLAY[group]}_core6_radar_H1in"); plt.close(fig)
    fig, axes = make_figure(data, features, GROUPS, COMBINED_SIZE)
    issues = qa(fig, axes); qa_log.append("Combined: " + ("PASS" if not issues else "FAIL; " + "; ".join(issues)))
    paths += save(fig, "E1-E5_core6_radar_horizontal_H1in"); plt.close(fig)
    corner_alpha = Image.open(OUT / "E1-E5_core6_radar_horizontal_H1in_transparent_600dpi.png").convert("RGBA").getpixel((0, 0))[3]
    if issues or any("FAIL" in x for x in qa_log) or corner_alpha != 0:
        raise RuntimeError("QA failed: " + " | ".join(qa_log) + f" | alpha={corner_alpha}")
    z = data[features].to_numpy(float)
    lines = [
        "CURRENT FIVE-CLASS CORE-SIX E RADAR — FROZEN INPUTS", f"Source matrix: {SOURCE}",
        "Cohort: 450 frozen GC-HC consensus cells.",
        "Classes: E1 n=142; E2 n=160; E3 n=45; E4 n=61; E5 n=42.",
        "No PCA, clustering, class assignment, feature order, or feature Z-score was recalculated.", "",
        "Feature order (identical in all panels):",
        *[f"{i+1}. {f} [{s}; source label: {label}]" for i, (f, s, label) in enumerate(zip(features, SHORT, full_labels))], "",
        f"Frozen six-feature Z-score matrix observed range: {np.nanmin(z):.6f} to {np.nanmax(z):.6f}.",
        f"Shared radar display range: Z={Z_MIN:.1f} to {Z_MAX:.1f}; values outside are display-clipped only; stored Z-scores are unchanged.",
        "Zero angle: 12 o'clock; clockwise; identical radius and feature angles in every panel.",
        "Center: feature-wise median; interval: feature-wise 25th-75th percentiles.",
        f"Cells: {CELL_LW:.2f} pt, alpha={CELL_ALPHA:.2f}; median: {MEDIAN_LW:.2f} pt; median fill alpha={MEDIAN_FILL_ALPHA:.2f}; IQR alpha={IQR_ALPHA:.2f}.",
        f"Grid: {GRID_LW:.2f} pt, alpha={GRID_ALPHA:.2f}; Arial {FONT_SIZE:.1f} pt; no radial values; no legend.",
        "Colors: " + "; ".join(f"{DISPLAY[g]}={COLORS[g]}" for g in GROUPS) + ".",
        f"Individual canvas: {INDIVIDUAL_SIZE[0]:.2f} x {INDIVIDUAL_SIZE[1]:.2f} in; combined: {COMBINED_SIZE[0]:.2f} x {COMBINED_SIZE[1]:.2f} in.",
        "Outputs: white and true-transparent PNG at 600 dpi; white and transparent vector PDF.", "",
        "Pre-export QA:", *qa_log, "Transparent PNG alpha channel: PASS (corner alpha=0).",
    ]
    (OUT / "E1-E5_core6_radar_parameters_and_QA.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    web = data.reset_index()[["Cell_ID", "E_consensus"] + features].copy()
    web["E_class"] = web["E_consensus"].map(DISPLAY)
    web.drop(columns="E_consensus").to_csv(OUT / "E1-E5_core6_frozen_Zscore_cells.csv", index=False)
    print("\n".join(map(str, paths))); print(OUT / "E1-E5_core6_radar_parameters_and_QA.txt")


if __name__ == "__main__":
    main()
