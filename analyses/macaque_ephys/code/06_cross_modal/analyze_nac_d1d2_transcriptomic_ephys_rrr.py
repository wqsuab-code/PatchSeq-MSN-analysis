from __future__ import annotations

from pathlib import Path
import re

import h5py
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.decomposition import PCA
from sklearn.model_selection import RepeatedStratifiedKFold


ROOT = Path(r"C:\Users\53461\OneDrive\Documentos\Patch-seq_Mouse_Acb_MSN_T-type_Visualization")
H5AD = ROOT / ".codex_tmp" / "macaque_rrr" / "Data" / "HMBA-Macaque-PatchSeq-BG-log2.h5ad"
PCA_DIR = ROOT / "outputs" / "NAc_E_rheobase_QC" / "PCA_YeoJohnson_Zscore_redundancy10"
ASSIGN = PCA_DIR / "HC_Seurat_scan" / "final_npcs3_res1p5" / "consensus4" / "NAc_MSN_Consensus4_cell_assignments.csv"
EPHYS = PCA_DIR / "NAc_MSN_rheobase10_YeoJohnson_Zscore_matrix.csv"
OUT = ROOT / "outputs" / "NAc_T_E_RRR_D1D2"

E_FEATURES = [
    "fast_trough_v_rheo",
    "peak_v_rheo",
    "postap_slope_rheo",
    "threshold_v_rheo",
    "upstroke_downstroke_ratio_rheo",
    "upstroke_rheo",
    "width_rheo_ms",
    "avg_rate_rheo",
    "latency_rheo",
    "rheobase_i",
]
E_LABELS = {
    "fast_trough_v_rheo": "Fast trough",
    "peak_v_rheo": "Peak V",
    "postap_slope_rheo": "Post-AP slope",
    "threshold_v_rheo": "Threshold V",
    "upstroke_downstroke_ratio_rheo": "Up/down ratio",
    "upstroke_rheo": "Upstroke",
    "width_rheo_ms": "AP width",
    "avg_rate_rheo": "Avg rate",
    "latency_rheo": "Latency",
    "rheobase_i": "Rheobase",
}
CLASS_COLORS = {"D1": "#D95F02", "D2": "#008F7A"}


def decode(values: np.ndarray) -> list[str]:
    return [x.decode("utf-8") if isinstance(x, bytes) else str(x) for x in values]


def read_h5ad_subset(cell_ids: list[str]) -> tuple[np.ndarray, np.ndarray]:
    with h5py.File(H5AD, "r") as handle:
        obs_ids = decode(handle["obs"]["cell_label"][:])
        obs_lookup = {cell: i for i, cell in enumerate(obs_ids)}
        missing = [cell for cell in cell_ids if cell not in obs_lookup]
        if missing:
            raise ValueError(f"Transcriptomic matrix is missing {len(missing)} cells: {missing[:5]}")

        x_group = handle["X"]
        shape = tuple(int(x) for x in x_group.attrs["shape"])
        matrix = sparse.csr_matrix(
            (x_group["data"][:], x_group["indices"][:], x_group["indptr"][:]),
            shape=shape,
        )
        selected = matrix[[obs_lookup[cell] for cell in cell_ids], :].toarray().astype(np.float64)

        gene_group = handle["var"]["gene_symbol"]
        categories = decode(gene_group["categories"][:])
        codes = gene_group["codes"][:]
        genes = np.array([categories[int(code)] if int(code) >= 0 else "" for code in codes], dtype=object)
    return selected, genes


def technical_gene(gene: str) -> bool:
    return bool(
        re.match(r"^(MT-|mt-|RPL|RPS|Rpl|Rps|GM\d|Gm\d)", gene)
        or re.search(r"Rik$", gene)
        or gene in {"MALAT1", "Malat1", "XIST", "Xist", ""}
    )


def transcriptomic_pca(expr: np.ndarray, genes: np.ndarray, n_hvg: int = 1000, n_pc: int = 20):
    detection = np.sum(expr > 0, axis=0)
    variance = np.var(expr, axis=0, ddof=1)
    eligible = np.array(
        [
            i
            for i, gene in enumerate(genes)
            if detection[i] >= 8 and np.isfinite(variance[i]) and variance[i] > 0 and not technical_gene(str(gene))
        ],
        dtype=int,
    )
    eligible = eligible[np.argsort(variance[eligible])[::-1]]
    # Keep one representative per gene symbol, prioritizing the most variable row.
    selected = []
    seen = set()
    for i in eligible:
        gene = str(genes[i])
        if gene in seen:
            continue
        seen.add(gene)
        selected.append(int(i))
        if len(selected) >= n_hvg:
            break
    selected = np.array(selected, dtype=int)
    x = expr[:, selected]
    x = (x - x.mean(axis=0)) / x.std(axis=0, ddof=0)
    x = np.nan_to_num(x, nan=0.0, posinf=0.0, neginf=0.0)
    pca = PCA(n_components=min(n_pc, x.shape[0] - 1, x.shape[1]), svd_solver="full")
    scores = pca.fit_transform(x)
    return scores, pca.components_.T, genes[selected], pca.explained_variance_ratio_, x


def class_center(matrix: np.ndarray, labels: np.ndarray) -> np.ndarray:
    out = matrix.copy()
    for label in np.unique(labels):
        mask = labels == label
        out[mask] -= out[mask].mean(axis=0, keepdims=True)
    return out


def fit_rrr(x: np.ndarray, y: np.ndarray, rank: int):
    x_mean = x.mean(axis=0)
    y_mean = y.mean(axis=0)
    xc = x - x_mean
    yc = y - y_mean
    beta = np.linalg.pinv(xc) @ yc
    fitted = xc @ beta
    _, _, vt = np.linalg.svd(fitted, full_matrices=False)
    v = vt.T[:, :rank]
    beta_rrr = beta @ v @ v.T
    return x_mean, y_mean, beta_rrr, v


def cv_rank_scan(x: np.ndarray, y: np.ndarray, labels: np.ndarray, max_rank: int = 5) -> pd.DataFrame:
    splitter = RepeatedStratifiedKFold(n_splits=5, n_repeats=20, random_state=20260902)
    records = []
    for rank in range(1, max_rank + 1):
        repeat_r2 = []
        for repeat in range(20):
            pred = np.zeros_like(y)
            test_seen = np.zeros(y.shape[0], dtype=bool)
            fold_iter = list(splitter.split(x, labels))[repeat * 5 : (repeat + 1) * 5]
            for train, test in fold_iter:
                xm, ym, beta, _ = fit_rrr(x[train], y[train], rank)
                pred[test] = (x[test] - xm) @ beta + ym
                test_seen[test] = True
            if not test_seen.all():
                raise RuntimeError("Cross-validation did not predict every cell")
            sse = float(np.sum((y - pred) ** 2))
            sst = float(np.sum((y - y.mean(axis=0)) ** 2))
            repeat_r2.append(1.0 - sse / sst)
        records.append(
            {
                "rank": rank,
                "CV_R2_mean": float(np.mean(repeat_r2)),
                "CV_R2_sd": float(np.std(repeat_r2, ddof=1)),
                "CV_R2_min": float(np.min(repeat_r2)),
                "CV_R2_max": float(np.max(repeat_r2)),
            }
        )
    return pd.DataFrame(records)


def ellipse(ax, xy: np.ndarray, color: str) -> None:
    if len(xy) < 4:
        return
    cov = np.cov(xy.T)
    vals, vecs = np.linalg.eigh(cov)
    order = np.argsort(vals)[::-1]
    vals, vecs = vals[order], vecs[:, order]
    angle = np.degrees(np.arctan2(vecs[1, 0], vecs[0, 0]))
    radius80 = np.sqrt(3.2188758249)
    patch = Ellipse(
        xy.mean(axis=0),
        width=2 * radius80 * np.sqrt(max(vals[0], 0)),
        height=2 * radius80 * np.sqrt(max(vals[1], 0)),
        angle=angle,
        facecolor=color,
        edgecolor=color,
        alpha=0.10,
        linewidth=0.8,
        zorder=0,
    )
    ax.add_patch(patch)
    outline = Ellipse(
        xy.mean(axis=0),
        width=patch.width,
        height=patch.height,
        angle=angle,
        facecolor="none",
        edgecolor=color,
        alpha=0.9,
        linewidth=0.75,
        zorder=1,
    )
    ax.add_patch(outline)


def orient_axes(v: np.ndarray, e_names: list[str]) -> np.ndarray:
    v = v.copy()
    for axis in range(v.shape[1]):
        anchor = int(np.argmax(np.abs(v[:, axis])))
        if v[anchor, axis] < 0:
            v[:, axis] *= -1
    return v


def make_plot(
    x: np.ndarray,
    y: np.ndarray,
    labels: np.ndarray,
    gene_pca_loadings: np.ndarray,
    genes: np.ndarray,
    cv: pd.DataFrame,
    mode: str,
) -> tuple[Path, pd.DataFrame, pd.DataFrame]:
    display_rank = 2
    xm, ym, _, v_raw = fit_rrr(x, y, display_rank)
    v = orient_axes(v_raw, E_FEATURES)
    beta = np.linalg.pinv(x - xm) @ (y - ym)
    beta_rrr = beta @ v @ v.T
    predicted = (x - xm) @ beta_rrr
    scores = predicted @ v
    gene_coeff = gene_pca_loadings @ (beta @ v)

    e_mag = np.linalg.norm(v, axis=1)
    e_idx = np.argsort(e_mag)[::-1][:7]
    g_mag = np.linalg.norm(gene_coeff, axis=1)
    g_idx = np.argsort(g_mag)[::-1][:6]

    fig, (ax_score, ax_load) = plt.subplots(1, 2, figsize=(6.4, 2.8), gridspec_kw={"wspace": 0.43})
    for group in ["D1", "D2"]:
        mask = labels == group
        ellipse(ax_score, scores[mask], CLASS_COLORS[group])
        ax_score.scatter(
            scores[mask, 0],
            scores[mask, 1],
            s=11,
            color=CLASS_COLORS[group],
            edgecolor="white",
            linewidth=0.30,
            alpha=0.82,
            label=f"{group}  n={int(mask.sum())}",
            zorder=2,
        )
    ax_score.axhline(0, color="#B8B8B8", lw=0.45, zorder=-1)
    ax_score.axvline(0, color="#B8B8B8", lw=0.45, zorder=-1)
    ax_score.set_xlabel("RRR1")
    ax_score.set_ylabel("RRR2")
    ax_score.legend(frameon=False, loc="best")
    ax_score.set_title("Cell scores", fontweight="bold")

    ax_load.axhline(0, color="#B8B8B8", lw=0.45)
    ax_load.axvline(0, color="#B8B8B8", lw=0.45)
    e_scale = 0.87 / max(np.max(np.abs(v[e_idx])), 1e-12)
    g_scale = 0.87 / max(np.max(np.abs(gene_coeff[g_idx])), 1e-12)
    for i in e_idx:
        xx, yy = v[i, :2] * e_scale
        ax_load.annotate("", xy=(xx, yy), xytext=(0, 0), arrowprops=dict(arrowstyle="-|>", color="#E69F00", lw=0.70))
        ax_load.text(xx * 1.05, yy * 1.05, E_LABELS[E_FEATURES[i]], color="#8A5A00", fontsize=4, ha="center", va="center")
    for i in g_idx:
        xx, yy = gene_coeff[i, :2] * g_scale
        ax_load.annotate("", xy=(xx, yy), xytext=(0, 0), arrowprops=dict(arrowstyle="-|>", color="#7B3294", lw=0.60))
        ax_load.text(xx * 1.06, yy * 1.06, str(genes[i]), color="#6A1B7B", fontsize=4, ha="center", va="center")
    ax_load.set_xlim(-1.12, 1.12)
    ax_load.set_ylim(-1.12, 1.12)
    ax_load.set_aspect("equal", adjustable="box")
    ax_load.set_xlabel("RRR1 loading")
    ax_load.set_ylabel("RRR2 loading")
    ax_load.set_title("T-gene and E-feature directions", fontweight="bold")
    ax_load.text(0.02, 0.98, "T genes", color="#7B3294", transform=ax_load.transAxes, va="top", fontweight="bold")
    ax_load.text(0.02, 0.91, "E features", color="#E69F00", transform=ax_load.transAxes, va="top", fontweight="bold")

    best = cv.loc[cv["CV_R2_mean"].idxmax()]
    mode_title = "D1/D2 pooled T→E RRR" if mode == "pooled" else "D1/D2-adjusted T→E RRR"
    fig.suptitle(
        f"NAc {mode_title} | 2D display; CV-best rank={int(best['rank'])}, R²={best['CV_R2_mean']:.3f}",
        fontsize=6,
        fontweight="bold",
        y=0.985,
    )
    fig.subplots_adjust(left=0.09, right=0.98, bottom=0.15, top=0.88)
    path = OUT / f"NAc_D1D2_T_E_RRR_{mode}.png"
    fig.savefig(path, dpi=600, facecolor="white", bbox_inches="tight", pad_inches=0.035)
    plt.close(fig)

    score_df = pd.DataFrame({"T_class": labels, "RRR1": scores[:, 0], "RRR2": scores[:, 1]})
    loading_df = pd.concat(
        [
            pd.DataFrame(
                {
                    "domain": "E",
                    "feature": E_FEATURES,
                    "RRR1": v[:, 0],
                    "RRR2": v[:, 1],
                }
            ),
            pd.DataFrame(
                {
                    "domain": "T",
                    "feature": genes,
                    "RRR1": gene_coeff[:, 0],
                    "RRR2": gene_coeff[:, 1],
                }
            ),
        ],
        ignore_index=True,
    )
    return path, score_df, loading_df


def correlation_loadings(features: np.ndarray, scores: np.ndarray) -> np.ndarray:
    f = features - features.mean(axis=0, keepdims=True)
    s = scores - scores.mean(axis=0, keepdims=True)
    f_sd = np.sqrt(np.sum(f * f, axis=0, keepdims=True))
    s_sd = np.sqrt(np.sum(s * s, axis=0, keepdims=True))
    denom = f_sd.T @ s_sd
    return np.divide(f.T @ s, denom, out=np.zeros((f.shape[1], s.shape[1])), where=denom > 0)


def spread_label_y(y_values: np.ndarray, minimum_gap: float = 0.115) -> np.ndarray:
    if len(y_values) < 2:
        return y_values.copy()
    order = np.argsort(y_values)
    placed = y_values[order].copy()
    for i in range(1, len(placed)):
        placed[i] = max(placed[i], placed[i - 1] + minimum_gap)
    overflow = max(0.0, placed[-1] - 0.98)
    placed -= overflow
    for i in range(len(placed) - 2, -1, -1):
        placed[i] = min(placed[i], placed[i + 1] - minimum_gap)
    underflow = max(0.0, -0.98 - placed[0])
    placed += underflow
    result = np.empty_like(placed)
    result[order] = placed
    return result


def draw_correlation_vectors(
    ax: plt.Axes,
    loadings: np.ndarray,
    names: np.ndarray,
    pair: tuple[int, int],
    top_n: int,
) -> None:
    vectors = loadings[:, list(pair)]
    magnitude = np.linalg.norm(vectors, axis=1)
    selected = np.argsort(magnitude)[::-1][:top_n]
    vectors = vectors[selected]
    labels = np.asarray(names, dtype=object)[selected]

    for vector in vectors:
        ax.annotate(
            "",
            xy=(vector[0], vector[1]),
            xytext=(0, 0),
            arrowprops=dict(arrowstyle="-|>", color="black", lw=0.55, shrinkA=0, shrinkB=0, mutation_scale=5),
            zorder=3,
        )

    side = np.where(vectors[:, 0] >= 0, 1.0, -1.0)
    label_y = np.zeros(len(vectors))
    for side_value in [-1.0, 1.0]:
        mask = side == side_value
        label_y[mask] = spread_label_y(np.clip(vectors[mask, 1] * 1.06, -0.96, 0.96))
    label_x = side * 1.015

    for vector, label, lx, ly, side_value in zip(vectors, labels, label_x, label_y, side):
        ax.plot([vector[0], lx], [vector[1], ly], color="#505050", lw=0.32, zorder=2)
        ax.text(
            lx,
            ly,
            str(label),
            ha="left" if side_value > 0 else "right",
            va="center",
            fontsize=4,
            color="#202020",
            bbox=dict(boxstyle="round,pad=0.12", facecolor="white", edgecolor="#666666", linewidth=0.35),
            zorder=4,
            clip_on=False,
        )


def make_nature_style_plot(
    x: np.ndarray,
    y: np.ndarray,
    gene_matrix: np.ndarray,
    genes: np.ndarray,
    labels: np.ndarray,
    cv: pd.DataFrame,
    mode: str,
) -> Path:
    xm, ym, _, v_raw = fit_rrr(x, y, rank=3)
    v = orient_axes(v_raw, E_FEATURES)
    beta = np.linalg.pinv(x - xm) @ (y - ym)
    t_projection = (x - xm) @ beta @ v
    e_projection = (y - ym) @ v
    t_loadings = correlation_loadings(gene_matrix, t_projection)
    e_loadings = correlation_loadings(y, e_projection)

    projections = [t_projection, e_projection]
    loadings = [t_loadings, e_loadings]
    feature_names = [genes, np.array([E_LABELS[x] for x in E_FEATURES], dtype=object)]
    top_n = [10, 10]
    pairs = [(0, 1), (0, 2)]

    fig, axes = plt.subplots(2, 2, figsize=(5.8, 5.25), gridspec_kw={"wspace": 0.16, "hspace": 0.11})
    for row, pair in enumerate(pairs):
        for col in range(2):
            ax = axes[row, col]
            scores = projections[col][:, list(pair)]
            radius = np.sqrt(np.sum(scores * scores, axis=1))
            scale = max(float(np.quantile(radius, 0.985)), 1e-12)
            plotted = scores / scale * 0.92
            theta = np.linspace(0, 2 * np.pi, 361)
            ax.plot(np.cos(theta), np.sin(theta), color="#4A4A4A", lw=0.62, zorder=0)
            ax.axhline(0, color="#B0B0B0", lw=0.36, zorder=0)
            ax.axvline(0, color="#B0B0B0", lw=0.36, zorder=0)
            for group in ["D1", "D2"]:
                mask = labels == group
                ax.scatter(
                    plotted[mask, 0],
                    plotted[mask, 1],
                    s=9.5,
                    color=CLASS_COLORS[group],
                    alpha=0.80,
                    edgecolor="none",
                    zorder=1,
                    label=group,
                )
            draw_correlation_vectors(ax, loadings[col], feature_names[col], pair, top_n[col])
            ax.set_xlim(-1.18, 1.18)
            ax.set_ylim(-1.18, 1.18)
            ax.set_aspect("equal", adjustable="box")
            ax.set_xticks([])
            ax.set_yticks([])
            for spine in ax.spines.values():
                spine.set_visible(False)
            if row == 0:
                ax.set_title("Transcriptomic space" if col == 0 else "Electrophysiological space", fontsize=6, pad=4)
            if col == 0:
                ax.set_ylabel(f"Component {pair[1] + 1}", fontsize=5, labelpad=1)
            if row == 1:
                ax.set_xlabel("Component 1", fontsize=5, labelpad=1)

    axes[0, 0].text(-1.32, 1.13, "a", fontsize=7, fontweight="bold", ha="left", va="top", clip_on=False)
    axes[1, 0].text(-1.32, 1.13, "b", fontsize=7, fontweight="bold", ha="left", va="top", clip_on=False)
    axes[0, 1].legend(
        frameon=False,
        loc="upper right",
        bbox_to_anchor=(1.02, 1.03),
        handletextpad=0.25,
        borderaxespad=0,
        markerscale=0.9,
    )
    best = cv.loc[cv["CV_R2_mean"].idxmax()]
    mode_label = "pooled" if mode == "pooled" else "D1/D2-adjusted"
    fig.text(
        0.5,
        0.018,
        f"NAc D1/D2 {mode_label} PCA-RRR; n=80. Rank-3 display; repeated 5-fold CV favours rank {int(best['rank'])} (R²={best['CV_R2_mean']:.3f}).",
        ha="center",
        va="bottom",
        fontsize=4,
        color="#404040",
    )
    fig.subplots_adjust(left=0.065, right=0.985, bottom=0.065, top=0.965)
    path = OUT / f"NAc_D1D2_T_E_RRR_NatureFig2_style_{mode}.png"
    fig.savefig(path, dpi=600, facecolor="white", bbox_inches="tight", pad_inches=0.04)
    plt.close(fig)
    return path


def main() -> None:
    mpl.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 4,
            "axes.titlesize": 5,
            "axes.labelsize": 4,
            "xtick.labelsize": 4,
            "ytick.labelsize": 4,
            "legend.fontsize": 4,
            "axes.linewidth": 0.5,
        }
    )
    OUT.mkdir(parents=True, exist_ok=True)

    assignment = pd.read_csv(ASSIGN)
    assignment = assignment.loc[assignment["T_class"].isin(["D1", "D2"])].copy()
    ephys = pd.read_csv(EPHYS)
    data = assignment[["cell_label", "donor_label", "T_class", "Consensus4"]].merge(
        ephys[["cell_label"] + E_FEATURES], on="cell_label", how="inner", validate="one_to_one"
    )
    if len(data) != 80:
        raise ValueError(f"Expected 80 D1/D2 cells, found {len(data)}")

    expr, all_genes = read_h5ad_subset(data["cell_label"].tolist())
    t_scores, gene_loadings, genes, t_var, gene_matrix = transcriptomic_pca(expr, all_genes)
    y = data[E_FEATURES].to_numpy(float)
    labels = data["T_class"].to_numpy(str)

    all_cv = []
    outputs = []
    for mode in ["pooled", "class_adjusted"]:
        x_model = t_scores if mode == "pooled" else class_center(t_scores, labels)
        y_model = y if mode == "pooled" else class_center(y, labels)
        cv = cv_rank_scan(x_model, y_model, labels)
        cv.insert(0, "model", mode)
        all_cv.append(cv)
        path, scores, loadings = make_plot(x_model, y_model, labels, gene_loadings, genes, cv, mode)
        gene_model = gene_matrix if mode == "pooled" else class_center(gene_matrix, labels)
        nature_path = make_nature_style_plot(x_model, y_model, gene_model, genes, labels, cv, mode)
        scores.insert(0, "cell_label", data["cell_label"].to_numpy())
        scores.to_csv(OUT / f"NAc_D1D2_T_E_RRR_{mode}_cell_scores.csv", index=False)
        loadings.to_csv(OUT / f"NAc_D1D2_T_E_RRR_{mode}_loadings.csv", index=False)
        outputs.extend([path, nature_path])

    cv_all = pd.concat(all_cv, ignore_index=True)
    cv_all.to_csv(OUT / "NAc_D1D2_T_E_RRR_rank1-5_repeated5fold_CV.csv", index=False)
    pd.DataFrame(
        {
            "T_PC": [f"PC{i+1}" for i in range(len(t_var))],
            "explained_variance_ratio": t_var,
            "cumulative_variance": np.cumsum(t_var),
        }
    ).to_csv(OUT / "NAc_D1D2_transcriptomic_PCA_variance.csv", index=False)
    (OUT / "README.txt").write_text(
        "NAc D1/D2 T-to-E reduced-rank regression.\n"
        "Cohort: 80 cells (D1=44, D2=36); Hybrid cells excluded.\n"
        "T: log2 expression from HMBA-Macaque-PatchSeq-BG-log2.h5ad; top 1000 non-technical variable genes; gene z-score; top 20 PCs.\n"
        "E: 10 retained nonredundant rheobase features after Yeo-Johnson and z-score.\n"
        "Rank scan: ranks 1-5, repeated stratified 5-fold CV (20 repeats).\n"
        "Pooled model includes the D1/D2 group contrast. Class-adjusted model subtracts D1/D2-specific means from both T-PC and E matrices before RRR.\n",
        encoding="utf-8",
    )
    print("\n".join(str(x) for x in outputs))
    print(cv_all.to_string(index=False))


if __name__ == "__main__":
    main()
