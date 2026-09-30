#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
})

options(stringsAsFactors = FALSE)
set.seed(20260914)

root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
assign_file <- file.path(
  root, "outputs", "morph_qc", "final_morph187_NPC3_HCK4_GCres0p50_figures",
  "00_final_morph187_cell_assignments.csv"
)
morph_z_file <- file.path(
  root, "outputs", "morph_qc", "final_morph187_NPC3_HCK4_GCres0p50_figures",
  "08_heatmap_ordered_zscore_matrix.csv"
)
raw_rds <- file.path(
  "C:/Users/53461/Documents/Codex/2026-07-13/zh/outputs",
  "05_branch_specific_subtype_mapping", "01_MSN_input",
  "MSN_raw_branch_counts_and_metadata.rds"
)
out_dir <- file.path(
  root, "outputs", "morph_qc", "final_morph187_Mouse_T_M_RRR"
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

id_col <- "MSN_unique_ID"
morph_features <- c(
  "M_soma_circularity_index",
  "M_soma_aspect_ratio",
  "M_cell_max_radial_dist",
  "M_total_number_of_neurites",
  "M_basal_dendrite_avg_tortuosity",
  "M_Total_neurite_length_(sections)",
  "M_Number_of_bifurcation_points",
  "M_Maximum_branch_order",
  "M_trunk_angle_min",
  "M_trunk_angle_max"
)
morph_labels <- c(
  "Soma circularity",
  "Soma aspect ratio",
  "Max radial distance",
  "Primary neurite number",
  "Mean tortuosity",
  "Total neurite length",
  "Bifurcation points",
  "Maximum branch order",
  "Minimum trunk angle",
  "Maximum trunk angle"
)
names(morph_labels) <- morph_features
m_levels <- c("M1", "M2", "M3", "M4")
m_colours <- c(M1 = "#00468B", M2 = "#42B540", M3 = "#ED0000", M4 = "#0099B4")

assignments <- fread(assign_file)
morph_z <- fread(morph_z_file)
stopifnot(
  nrow(assignments) == 187L,
  uniqueN(assignments[[id_col]]) == 187L,
  all(morph_features %in% names(morph_z))
)

cohort <- assignments[HC_GC_consensus == TRUE & D1_D2 %chin% c("D1", "D2")]
cohort[, M_class := factor(M_class, levels = m_levels)]
cohort[, T_display := as.character(D1_D2)]
stopifnot(nrow(cohort) == 180L, uniqueN(cohort[[id_col]]) == 180L)

z_dt <- morph_z[, c(id_col, morph_features), with = FALSE]
cohort <- merge(cohort, z_dt, by = id_col, all.x = TRUE, sort = FALSE, suffixes = c("", ".z"))
# The assignment table contains raw morphology values with the same names.  After
# merge(), the frozen transformed values carry the .z suffix; rename them explicitly.
z_names <- paste0(morph_features, ".z")
stopifnot(all(z_names %in% names(cohort)), !anyNA(cohort[, ..z_names]))
y <- as.matrix(cohort[, ..z_names])
colnames(y) <- morph_features

raw <- readRDS(raw_rds)
counts_all <- raw$query_counts
stopifnot(inherits(counts_all, "Matrix"), all(cohort[[id_col]] %in% colnames(counts_all)))
counts <- counts_all[, cohort[[id_col]], drop = FALSE]
library_size <- Matrix::colSums(counts)
stopifnot(all(is.finite(library_size)), all(library_size > 0))

norm_counts <- counts %*% Diagonal(x = 10000 / library_size)
dimnames(norm_counts) <- dimnames(counts)
log_norm <- log1p(norm_counts)
gene_names <- rownames(log_norm)
gene_mean <- Matrix::rowMeans(log_norm)
gene_sq_mean <- Matrix::rowMeans(log_norm ^ 2)
gene_var <- pmax(gene_sq_mean - gene_mean ^ 2, 0)
gene_detected <- Matrix::rowSums(counts > 0)
technical <- grepl(
  "^(mt-|Mt-|MT-|Rpl|Rps|RPL|RPS|Gm[0-9]|GM[0-9]|[0-9].*Rik$|Malat1$|MALAT1$|Xist$|XIST$)",
  gene_names
)
eligible <- which(!technical & gene_detected >= 8L & is.finite(gene_var) & gene_var > 0)
eligible <- eligible[order(gene_var[eligible], decreasing = TRUE)]
# Keep one row per gene symbol, prioritizing the most variable occurrence.
eligible <- eligible[!duplicated(gene_names[eligible])]
hvg_idx <- eligible[seq_len(min(1000L, length(eligible)))]
hvg <- gene_names[hvg_idx]
rna_log <- as.matrix(log_norm[hvg_idx, , drop = FALSE])
rna_scaled <- t(scale(t(rna_log), center = TRUE, scale = TRUE))
finite_gene <- apply(rna_scaled, 1L, function(x) all(is.finite(x)))
rna_scaled <- rna_scaled[finite_gene, , drop = FALSE]
hvg <- hvg[finite_gene]
stopifnot(nrow(rna_scaled) == 1000L)

rna_pca <- prcomp(t(rna_scaled), center = FALSE, scale. = FALSE, rank. = 20L)
x <- rna_pca$x[, seq_len(20L), drop = FALSE]
rownames(x) <- cohort[[id_col]]

fit_rrr <- function(x, y, rank = 3L) {
  x_mean <- colMeans(x)
  y_mean <- colMeans(y)
  xc <- sweep(x, 2L, x_mean, "-")
  yc <- sweep(y, 2L, y_mean, "-")
  beta <- qr.solve(xc, yc, tol = 1e-10)
  fitted <- xc %*% beta
  sv <- svd(fitted, nu = 0L, nv = min(rank, ncol(y)))
  v <- sv$v[, seq_len(rank), drop = FALSE]
  # Fix display signs using the largest absolute morphology coefficient.
  for (j in seq_len(ncol(v))) {
    anchor <- which.max(abs(v[, j]))
    if (v[anchor, j] < 0) v[, j] <- -v[, j]
  }
  beta_rrr <- beta %*% v %*% t(v)
  list(
    x_mean = x_mean,
    y_mean = y_mean,
    beta = beta,
    beta_rrr = beta_rrr,
    v = v,
    t_scores = xc %*% beta %*% v,
    m_scores = yc %*% v
  )
}

rrr <- fit_rrr(x, y, rank = 3L)
colnames(rrr$t_scores) <- paste0("RRR", 1:3)
colnames(rrr$m_scores) <- paste0("RRR", 1:3)

cor_loadings <- function(features_by_cell, scores) {
  out <- cor(features_by_cell, scores, use = "pairwise.complete.obs")
  out[!is.finite(out)] <- 0
  out
}
t_loadings <- cor_loadings(t(rna_scaled), rrr$t_scores)
m_loadings <- cor_loadings(y, rrr$m_scores)
rownames(m_loadings) <- morph_labels[rownames(m_loadings)]

score_dt <- data.table(
  MSN_unique_ID = cohort[[id_col]],
  M_class = cohort$M_class,
  T_identity = cohort$T_display,
  T_RRR1 = rrr$t_scores[, 1], T_RRR2 = rrr$t_scores[, 2], T_RRR3 = rrr$t_scores[, 3],
  M_RRR1 = rrr$m_scores[, 1], M_RRR2 = rrr$m_scores[, 2], M_RRR3 = rrr$m_scores[, 3]
)
fwrite(score_dt, file.path(out_dir, "01_Mouse_HC-GC_consensus180_D1D2_T_M_RRR_scores.csv"))

loading_long <- rbindlist(list(
  data.table(Domain = "Transcriptomic", Feature = rownames(t_loadings),
             RRR1 = t_loadings[, 1], RRR2 = t_loadings[, 2], RRR3 = t_loadings[, 3]),
  data.table(Domain = "Morphological", Feature = rownames(m_loadings),
             RRR1 = m_loadings[, 1], RRR2 = m_loadings[, 2], RRR3 = m_loadings[, 3])
))
fwrite(loading_long, file.path(out_dir, "02_Mouse_T_M_RRR_correlation_loadings.csv"))
fwrite(data.table(Gene = hvg), file.path(out_dir, "03_Mouse_T_M_RRR_HVG1000.csv"))
fwrite(
  data.table(PC = seq_along(rna_pca$sdev),
             Variance_fraction = rna_pca$sdev^2 / sum(rna_pca$sdev^2),
             Cumulative_variance = cumsum(rna_pca$sdev^2 / sum(rna_pca$sdev^2))),
  file.path(out_dir, "04_Mouse_transcriptomic_PCA20_variance.csv")
)

circle <- data.table(theta = seq(0, 2 * pi, length.out = 361L))
circle[, `:=`(x = cos(theta), y = sin(theta))]

scale_scores_99 <- function(scores, pair) {
  xy <- scores[, pair, drop = FALSE]
  radius <- sqrt(rowSums(xy^2))
  sf <- 0.92 / max(as.numeric(quantile(radius, 0.99, na.rm = TRUE)), 1e-12)
  xy * sf
}

make_panel <- function(domain, scores, loadings, pair, title_text, top_n = 4L) {
  xy <- scale_scores_99(scores, pair)
  pts <- data.table(
    MSN_unique_ID = cohort[[id_col]],
    M_class = cohort$M_class,
    T_identity = cohort$T_display,
    x = xy[, 1], y = xy[, 2]
  )
  vectors <- data.table(
    Feature = rownames(loadings),
    x = loadings[, pair[1]],
    y = loadings[, pair[2]]
  )
  vectors[, Magnitude := sqrt(x^2 + y^2)]
  setorder(vectors, -Magnitude)
  vectors <- vectors[seq_len(min(top_n, .N))]

  p <- ggplot() +
    geom_path(data = circle, aes(x, y), colour = "grey35", linewidth = 1 / 2.8453) +
    geom_hline(yintercept = 0, colour = "grey75", linewidth = 0.25 / 2.8453) +
    geom_vline(xintercept = 0, colour = "grey75", linewidth = 0.25 / 2.8453) +
    stat_ellipse(
      data = pts[T_identity != "Hybrid/unstable"],
      aes(x, y, colour = M_class, group = M_class),
      type = "norm", level = 0.85, linewidth = 0.55 / 2.8453,
      show.legend = FALSE
    ) +
    geom_point(
      data = pts[T_identity != "Hybrid/unstable"],
      aes(x, y, colour = M_class), size = 0.72, alpha = 0.76,
      stroke = 0, show.legend = domain == "Transcriptomic"
    ) +
    geom_point(
      data = pts[T_identity == "Hybrid/unstable"],
      aes(x, y), colour = "grey55", size = 0.82, alpha = 0.90,
      stroke = 0, show.legend = FALSE
    ) +
    geom_segment(
      data = vectors,
      aes(x = 0, y = 0, xend = x, yend = y),
      arrow = arrow(length = unit(1.6, "pt"), type = "closed"),
      colour = "black", linewidth = 0.25 / 2.8453
    ) +
    geom_label_repel(
      data = vectors,
      aes(x, y, label = Feature),
      family = "Arial", size = 4 / 2.8453,
      colour = "grey20", fill = adjustcolor("white", alpha.f = 0.86),
      label.size = 0.10, label.padding = unit(0.035, "lines"),
      box.padding = unit(0.08, "lines"), point.padding = unit(0.02, "lines"),
      segment.colour = "grey55", segment.size = 0.12,
      min.segment.length = 0, max.overlaps = Inf, seed = 20260914,
      show.legend = FALSE
    ) +
    scale_colour_manual(values = m_colours, drop = FALSE) +
    coord_fixed(xlim = c(-1.45, 1.45), ylim = c(-1.16, 1.16), clip = "off") +
    labs(
      title = title_text,
      x = "Component 1",
      y = paste("Component", pair[2])
    ) +
    theme_void(base_family = "Arial", base_size = 4) +
    theme(
      text = element_text(family = "Arial", size = 4, colour = "black"),
      plot.title = element_text(size = 4, hjust = 0.5, margin = margin(b = 1.5, unit = "pt")),
      axis.title.x = element_text(size = 4, margin = margin(t = 1.5, unit = "pt")),
      axis.title.y = element_text(size = 4, angle = 90, margin = margin(r = 1.5, unit = "pt")),
      legend.position = "inside",
      legend.position.inside = c(0.79, 0.18),
      legend.title = element_blank(),
      legend.text = element_text(size = 3.7),
      legend.key.height = unit(5, "pt"),
      legend.key.width = unit(6, "pt"),
      legend.spacing.y = unit(0, "pt"),
      plot.margin = margin(3, 4, 3, 4, unit = "pt")
    )
  p
}

t12 <- make_panel("Transcriptomic", rrr$t_scores, t_loadings, c(1L, 2L), "Transcriptomic space")
m12 <- make_panel("Morphological", rrr$m_scores, m_loadings, c(1L, 2L), "Morphological space")
t13 <- make_panel("Transcriptomic", rrr$t_scores, t_loadings, c(1L, 3L), "Transcriptomic space")
m13 <- make_panel("Morphological", rrr$m_scores, m_loadings, c(1L, 3L), "Morphological space")

png_device <- if (requireNamespace("ragg", quietly = TRUE)) ragg::agg_png else "png"
save_plot <- function(stem, plot, width, height) {
  ggsave(file.path(out_dir, paste0(stem, ".png")), plot,
         width = width, height = height, units = "in", dpi = 900,
         device = png_device, bg = "white", limitsize = FALSE)
  ggsave(file.path(out_dir, paste0(stem, ".pdf")), plot,
         width = width, height = height, units = "in", device = cairo_pdf,
         bg = "white", limitsize = FALSE)
  ggsave(file.path(out_dir, paste0(stem, ".svg")), plot,
         width = width, height = height, units = "in", device = svglite::svglite,
         bg = "white", limitsize = FALSE)
}

save_plot("05_Mouse_T_M_RRR_Comp2_vs_Comp1", t12 | m12, 3.36, 1.24)
save_plot("06_Mouse_T_M_RRR_Comp3_vs_Comp1", t13 | m13, 3.36, 1.24)
save_plot("07_Mouse_T_M_RRR_4panel", (t12 | m12) / (t13 | m13), 3.36, 2.48)
save_plot("08a_Mouse_T_RRR_Comp2_vs_Comp1", t12, 1.68, 1.24)
save_plot("08b_Mouse_M_RRR_Comp2_vs_Comp1", m12, 1.68, 1.24)
save_plot("08c_Mouse_T_RRR_Comp3_vs_Comp1", t13, 1.68, 1.24)
save_plot("08d_Mouse_M_RRR_Comp3_vs_Comp1", m13, 1.68, 1.24)

t_counts <- cohort[, .N, by = T_display][order(T_display)]
m_counts <- cohort[, .N, by = M_class][order(M_class)]
summary_lines <- c(
  "Mouse T-M reduced-rank regression (RRR)",
  "=======================================",
  sprintf("Final morphology cohort: %d cells", nrow(assignments)),
  sprintf("HC-GC consensus cohort used for T-M RRR: %d cells", nrow(cohort)),
  sprintf("T identities: %s", paste(sprintf("%s=%d", t_counts$T_display, t_counts$N), collapse = ", ")),
  sprintf("M classes: %s", paste(sprintf("%s=%d", m_counts$M_class, m_counts$N), collapse = ", ")),
  "Morphology response: frozen 10-feature, redundancy-reduced transformed Z-score matrix used by the final Mouse M1-M4 taxonomy.",
  "Transcriptomics: library-size normalization to 10,000; log1p; technical-gene removal; detection in >=8 cells; top 1,000 variable genes; gene-wise Z-score; 20 PCs.",
  "RRR: pooled T->M model, rank=3. T panels show fitted transcriptomic shared scores; M panels show observed morphology response scores.",
  "Cohort restriction: only cells explicitly labelled D1 or D2 were included in model fitting; non-D1/D2 cells were excluded before preprocessing and RRR estimation.",
  "Loadings: Pearson correlations of HVG expression or morphology features with the displayed RRR scores; four largest bivariate loading magnitudes labelled per panel.",
  "Display: each score plane scaled independently so its 99th percentile radial distance equals 0.92; circles are therefore display guides, not confidence boundaries.",
  "Ellipses in the static script outputs: 85% normal-theory covariance ellipses by M class.",
  "Colours: M1 #00468B; M2 #42B540; M3 #ED0000; M4 #0099B4.",
  "Important: this is Mouse-only and was independently refitted; it does not reuse the Macaque 117-cell/18-feature RRR coordinates."
)
writeLines(summary_lines, file.path(out_dir, "00_Mouse_T_M_RRR_method_summary.txt"))

cat(paste(summary_lines, collapse = "\n"), "\n")
cat("Output:", normalizePath(out_dir, winslash = "/"), "\n")
