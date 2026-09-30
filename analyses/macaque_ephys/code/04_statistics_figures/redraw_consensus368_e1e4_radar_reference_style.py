"""Redraw the frozen 368-cell E radar matrix without refitting any analysis."""
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from PIL import Image
from scipy.interpolate import CubicSpline


ROOT = Path(__file__).resolve().parents[1]
ASSIGN = ROOT / "outputs/R3_panels_v3/00_tuned_p80_ee12_random_seed777_coordinates.csv"
RAW = ROOT / "outputs/dSTR_dSTRvSTR_E_QC/primary10_sensitivity_cohorts/dSTR_plus_vSTR_A_all_stage1QC_raw_with_NA.csv"
SCALE = ROOT / "outputs/R3_panels_v3/E_GC_res3_mergedK4_radar10_consensus368_clean_W4p64_H1_robust_scaling.csv"
OUT = ROOT / "outputs/R3_panels_v3/E1-E4_radar_reference_style_frozen368"

FEATURES = [
    "fast_trough_v_rheo", "peak_v_rheo", "postap_slope_rheo", "threshold_v_rheo",
    "upstroke_downstroke_ratio_rheo", "upstroke_rheo", "width_rheo_ms",
    "avg_rate_rheo", "latency_rheo", "rheobase_i",
]
LABELS = ["F.tr.", "Peak", "PostAP", "Thr.", "U/D", "Upstr.", "Width", "Rate", "Lat.", "Rheo."]
SOURCE_CLASSES = ["C1", "C2", "C3", "C4"]
DISPLAY_CLASSES = {"C1": "E1", "C2": "E2", "C3": "E3", "C4": "E4"}
COLORS = {"C1": "#F8766D", "C2": "#7CAE00", "C3": "#00BFC4", "C4": "#C77CFF"}
EXPECTED = {"C1": 59, "C2": 133, "C3": 57, "C4": 119}

R_MIN, R_MAX = 0.0, 1.0
CELL_LW, CELL_ALPHA = 0.30, 0.10
MEDIAN_LW, MEDIAN_ALPHA = 0.90, 1.0
MEDIAN_FILL_ALPHA, IQR_ALPHA = 0.18, 0.12
GRID_LW, GRID_ALPHA = 0.35, 0.42
FONT_SIZE = 4.0
INDIVIDUAL_SIZE = (1.30, 1.00)
COMBINED_SIZE = (5.20, 1.00)


def load_frozen_matrix() -> pd.DataFrame:
    assign = pd.read_csv(ASSIGN, dtype={"cell_label": str})
    raw = pd.read_csv(RAW, dtype={"cell_label": str})
    limits = pd.read_csv(SCALE).set_index("feature").reindex(FEATURES)
    if limits[["q2.5", "q97.5"]].isna().any().any():
        raise RuntimeError("Saved radar scaling parameters are incomplete")
    joined = assign.loc[assign["Consensus"].eq(True), ["cell_label", "HC_class"]].merge(
        raw[["cell_label"] + FEATURES], on="cell_label", how="inner", validate="one_to_one"
    )
    counts = joined["HC_class"].value_counts().reindex(SOURCE_CLASSES).to_dict()
    if len(joined) != 368 or counts != EXPECTED or joined[FEATURES].isna().any().any():
        raise RuntimeError(f"Frozen cohort changed: n={len(joined)}, counts={counts}")
    for feature in FEATURES:
        lo, hi = limits.loc[feature, ["q2.5", "q97.5"]]
        joined[feature] = np.clip((joined[feature].astype(float) - lo) / (hi - lo), 0, 1)
    return joined


ANGLES = np.linspace(0, 2 * np.pi, len(FEATURES), endpoint=False)
SMOOTH_ANGLES = np.linspace(0, 2 * np.pi, 361)


def smooth(values: np.ndarray) -> np.ndarray:
    closed_angles = np.r_[ANGLES, 2 * np.pi]
    curve = CubicSpline(closed_angles, np.r_[values, values[0]], bc_type="periodic")(SMOOTH_ANGLES)
    return np.clip(curve, R_MIN, R_MAX)


def style_axis(ax, group: str, matrix: np.ndarray) -> None:
    color = COLORS[group]
    for row in matrix:
        ax.plot(SMOOTH_ANGLES, smooth(row), color=color, lw=CELL_LW, alpha=CELL_ALPHA, zorder=1)
    q25, median, q75 = np.nanquantile(matrix, [0.25, 0.50, 0.75], axis=0)
    q25s, meds, q75s = smooth(q25), smooth(median), smooth(q75)
    ax.fill_between(SMOOTH_ANGLES, q25s, q75s, color=color, alpha=IQR_ALPHA, lw=0, zorder=2)
    ax.fill(SMOOTH_ANGLES, meds, color=color, alpha=MEDIAN_FILL_ALPHA, zorder=3)
    ax.plot(SMOOTH_ANGLES, meds, color="#202020", lw=MEDIAN_LW, alpha=MEDIAN_ALPHA, zorder=4)

    ax.set_theta_offset(np.pi / 2)
    ax.set_theta_direction(-1)
    ax.set_ylim(R_MIN, R_MAX)
    ax.set_xticks(ANGLES)
    ax.set_xticklabels([])
    ax.tick_params(axis="x", length=0)
    for label, theta, angle in zip(LABELS, ANGLES, np.degrees(ANGLES)):
        rotation = -angle
        if -270 < rotation < -90:
            rotation += 180
        ax.text(theta, 1.085, label, fontsize=FONT_SIZE, color="#303030",
                rotation=rotation, rotation_mode="anchor", ha="center", va="center", clip_on=False)
    ax.set_yticks([0.25, 0.50, 0.75, 1.00])
    ax.set_yticklabels([])
    ax.xaxis.grid(True, color="#A8A8A8", lw=GRID_LW, alpha=GRID_ALPHA)
    ax.yaxis.grid(True, color="#A8A8A8", lw=GRID_LW, alpha=GRID_ALPHA)
    ax.spines["polar"].set_color("#A8A8A8")
    ax.spines["polar"].set_linewidth(GRID_LW)
    ax.spines["polar"].set_alpha(GRID_ALPHA)
    ax.set_title(f"{DISPLAY_CLASSES[group]} (n={len(matrix)})", fontsize=FONT_SIZE, color=color, y=1.22, pad=0)


def make_figure(data: pd.DataFrame, groups: list[str], size: tuple[float, float]):
    fig = plt.figure(figsize=size)
    panel_w = 1 / len(groups)
    axes = []
    for i, group in enumerate(groups):
        # Identical physical polar areas; only their horizontal origins differ.
        left = i * panel_w + 0.060 * panel_w
        ax = fig.add_axes([left, 0.175, 0.880 * panel_w, 0.625], projection="polar")
        matrix = data.loc[data["HC_class"].eq(group), FEATURES].to_numpy(float)
        style_axis(ax, group, matrix)
        axes.append(ax)
    return fig, axes


def save_pair(fig, stem: str) -> list[Path]:
    paths = []
    for background, transparent in [("white", False), ("transparent", True)]:
        png = OUT / f"{stem}_{background}_600dpi.png"
        pdf = OUT / f"{stem}_{background}_vector.pdf"
        fig.savefig(png, dpi=600, facecolor="white" if not transparent else "none", transparent=transparent)
        fig.savefig(pdf, facecolor="white" if not transparent else "none", transparent=transparent)
        paths += [png, pdf]
    return paths


def text_qa(fig, axes) -> tuple[bool, list[str]]:
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    canvas = fig.bbox
    issues = []
    for i, ax in enumerate(axes, 1):
        labels = list(ax.texts) + [ax.title]
        boxes = [x.get_window_extent(renderer).expanded(0.98, 0.98) for x in labels if x.get_visible()]
        for j, box in enumerate(boxes):
            if box.x0 < canvas.x0 - 0.5 or box.y0 < canvas.y0 - 0.5 or box.x1 > canvas.x1 + 0.5 or box.y1 > canvas.y1 + 0.5:
                issues.append(f"panel {i}: label {j + 1} outside canvas")
        for a in range(len(boxes)):
            for b in range(a + 1, len(boxes)):
                if boxes[a].overlaps(boxes[b]):
                    issues.append(f"panel {i}: labels {a + 1} and {b + 1} overlap")
    return not issues, issues


def main() -> None:
    mpl.rcParams.update({
        "font.family": "Arial", "font.size": FONT_SIZE, "axes.titlesize": FONT_SIZE,
        "xtick.labelsize": FONT_SIZE, "ytick.labelsize": FONT_SIZE,
        "pdf.fonttype": 42, "ps.fonttype": 42, "savefig.pad_inches": 0,
    })
    OUT.mkdir(parents=True, exist_ok=True)
    data = load_frozen_matrix()
    outputs, qa_lines = [], []
    for group in SOURCE_CLASSES:
        fig, axes = make_figure(data, [group], INDIVIDUAL_SIZE)
        ok, issues = text_qa(fig, axes)
        qa_lines.append(f"{DISPLAY_CLASSES[group]}: {'PASS' if ok else 'FAIL'}" + ("; " + "; ".join(issues) if issues else ""))
        outputs += save_pair(fig, f"{DISPLAY_CLASSES[group]}_radar_H1in")
        plt.close(fig)
    fig, axes = make_figure(data, SOURCE_CLASSES, COMBINED_SIZE)
    ok, issues = text_qa(fig, axes)
    qa_lines.append(f"Combined: {'PASS' if ok else 'FAIL'}" + ("; " + "; ".join(issues) if issues else ""))
    outputs += save_pair(fig, "E1-E4_radar_horizontal_H1in")
    plt.close(fig)

    # Verify the transparent export really contains a zero-alpha corner pixel.
    alpha = Image.open(OUT / "E1-E4_radar_horizontal_H1in_transparent_600dpi.png").convert("RGBA").getpixel((0, 0))[3]
    if alpha != 0:
        raise RuntimeError(f"Transparent PNG failed alpha check: corner alpha={alpha}")
    if any("FAIL" in line for line in qa_lines):
        raise RuntimeError("Layout QA failed: " + " | ".join(qa_lines))

    params = [
        "Frozen analysis inputs", f"Assignments: {ASSIGN}", f"Raw feature matrix: {RAW}",
        f"Saved scaling parameters: {SCALE}", "Cohort: 368 HC-GC consensus cells",
        "Displayed classes: C1->E1 (n=59), C2->E2 (n=133), C3->E3 (n=57), C4->E4 (n=119)",
        "No PCA, clustering, class assignment, feature selection, or scaling parameters were refitted.", "",
        "Feature order (identical in every panel):",
        *[f"{i + 1}. {feature} [{label}]" for i, (feature, label) in enumerate(zip(FEATURES, LABELS))], "",
        "Original radar transformation: saved feature-wise 2.5th-97.5th percentile robust min-max scaling, clipped to [0,1].",
        "Radial range: 0.00 to 1.00; zero angle at 12 o'clock; clockwise direction.",
        "Center statistic: feature-wise median; dispersion: feature-wise 25th-75th percentile band.",
        f"Single-cell curves: {CELL_LW:.2f} pt, alpha={CELL_ALPHA:.2f}, class color, no markers.",
        f"Median contour: {MEDIAN_LW:.2f} pt, #202020; median fill alpha={MEDIAN_FILL_ALPHA:.2f}.",
        f"IQR band: class color, alpha={IQR_ALPHA:.2f}.",
        f"Grid and outer spine: {GRID_LW:.2f} pt, #A8A8A8, alpha={GRID_ALPHA:.2f}; radial labels hidden.",
        f"Font: Arial {FONT_SIZE:.1f} pt; no legend.",
        f"Class colors: E1={COLORS['C1']}; E2={COLORS['C2']}; E3={COLORS['C3']}; E4={COLORS['C4']}.",
        f"Individual canvas: {INDIVIDUAL_SIZE[0]:.2f} x {INDIVIDUAL_SIZE[1]:.2f} in; combined: {COMBINED_SIZE[0]:.2f} x {COMBINED_SIZE[1]:.2f} in.",
        "Raster output: PNG 600 dpi; vector output: PDF; white and true transparent backgrounds.", "",
        "Pre-export QA:", *qa_lines, "Transparent PNG alpha channel: PASS (corner alpha=0).",
    ]
    (OUT / "E1-E4_radar_parameters_and_QA.txt").write_text("\n".join(params) + "\n", encoding="utf-8")
    print("\n".join(map(str, outputs)))
    print(OUT / "E1-E4_radar_parameters_and_QA.txt")


if __name__ == "__main__":
    main()
