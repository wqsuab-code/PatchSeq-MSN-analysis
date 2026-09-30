#!/usr/bin/env Rscript

# Macaque-only 18-feature morphology PCA/outlier analysis.
# 3_Sholl_PC1 is excluded. Source values are read directly from the ZIP.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(jsonlite)
  library(digest)
  library(MASS)
  library(dbscan)
  library(ggrepel)
  library(pheatmap)
})

args <- commandArgs(trailingOnly = TRUE)
zip_path <- if (length(args) >= 1L) args[[1L]] else "C:/Users/53461/Downloads/Macaque-PatchSeq-BG.zip"
out <- if (length(args) >= 2L) args[[2L]] else file.path(getwd(), "macaque_m", "m18")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

expected_sha <- "8b0aeaed726e27658066230fe467bdbb765d641d1ea1e3a115e079e4ff0c13a6"
features <- c(
  "basal_dendrite_bias_dorsal",
  "basal_dendrite_bias_medial",
  "basal_dendrite_calculate_number_of_stems",
  "basal_dendrite_extent_dorsal",
  "basal_dendrite_extent_medial",
  "basal_dendrite_max_branch_order",
  "basal_dendrite_max_euclidean_distance",
  "basal_dendrite_max_path_distance",
  "basal_dendrite_mean_contraction",
  "basal_dendrite_mean_diameter",
  "basal_dendrite_num_branches",
  "basal_dendrite_soma_percentile_dorsal",
  "basal_dendrite_soma_percentile_medial",
  "basal_dendrite_stem_exit_MedialLateral",
  "basal_dendrite_stem_exit_dorsal",
  "basal_dendrite_stem_exit_ventral",
  "basal_dendrite_total_length",
  "soma_surface_area"
)
id <- "cell_label"
write_dt <- function(x, name) fwrite(as.data.table(x), file.path(out, name), bom = TRUE, na = "")

stopifnot(file.exists(zip_path))
zip_sha <- digest(file = zip_path, algo = "sha256", serialize = FALSE)
stopifnot(identical(tolower(zip_sha), expected_sha))
meta <- as.data.table(read.csv(unz(zip_path, "Data/cell_metadata_AllCell.csv"), check.names = FALSE))
morph <- as.data.table(read.csv(unz(zip_path, "Data/morphology_features.csv"), check.names = FALSE))
stopifnot(nrow(meta) == 717L, nrow(morph) == 717L, !anyDuplicated(meta[[id]]), !anyDuplicated(morph[[id]]))
dat <- merge(meta, morph, by = id, all = FALSE, sort = FALSE)
target <- dat[Lib_region_of_interest_label %chin% c("Ca", "Pu", "NAC") &
              Subclass_name %chin% c("STR D1 MSN", "STR D2 MSN", "STR Hybrid MSN")]
cohort <- target[complete.cases(target[, ..features])]
setorder(cohort, cell_label)
stopifnot(nrow(target) == 486L, nrow(cohort) == 126L, length(features) == 18L,
          !"3_Sholl_PC1" %chin% features)

raw_x <- as.matrix(cohort[, ..features]); rownames(raw_x) <- cohort[[id]]
mins <- apply(raw_x, 2L, min)
shifted <- sweep(raw_x, 2L, mins, "-")
shifted[shifted < 0 & shifted > -1e-12] <- 0
totals <- colSums(shifted); stopifnot(all(totals > 0))
sum10000 <- sweep(shifted, 2L, totals, "/") * 10000
logged <- log1p(sum10000)
z <- scale(logged, center = TRUE, scale = TRUE)
stopifnot(all(is.finite(z)), max(abs(colMeans(z))) < 1e-12,
          max(abs(apply(z, 2, sd) - 1)) < 1e-12)

pca <- prcomp(z, center = FALSE, scale. = FALSE)
for (j in seq_len(ncol(pca$rotation))) {
  anchor <- which.max(abs(pca$rotation[, j]))
  if (pca$rotation[anchor, j] < 0) {
    pca$rotation[, j] <- -pca$rotation[, j]
    pca$x[, j] <- -pca$x[, j]
  }
}
explained <- 100 * pca$sdev^2 / sum(pca$sdev^2)
cum_explained <- cumsum(explained)
npc80 <- which(cum_explained >= 80)[1L]
score_z <- scale(pca$x)

exceedances <- rbindlist(lapply(seq_len(ncol(pca$x)), function(j) {
  data.table(cell_label = cohort[[id]], PC = paste0("PC", j), PC_number = j,
             score = pca$x[, j], classical_z = score_z[, j], abs_z = abs(score_z[, j]))[abs_z > 3]
}))
cells3 <- exceedances[PC_number <= 3, .(
  exceeding_PCs = paste(PC, collapse = ";"), max_abs_z = max(abs_z)), by = cell_label]
cells80 <- exceedances[PC_number <= npc80, .(
  exceeding_PCs = paste(PC, collapse = ";"), max_abs_z = max(abs_z)), by = cell_label]

set.seed(777L)
mve <- cov.rob(pca$x[, seq_len(npc80), drop = FALSE], method = "mve")
robust_d2 <- mahalanobis(pca$x[, seq_len(npc80), drop = FALSE], mve$center, mve$cov)
d2_cutoff <- qchisq(1 - 0.05 / nrow(cohort), df = npc80)
lof5 <- lof(pca$x[, seq_len(npc80), drop = FALSE], minPts = 5L)

diag <- data.table(
  cell_label = cohort[[id]], donor_label = cohort$donor_label,
  ROI = cohort$Lib_region_of_interest_label, Subclass = cohort$Subclass_name,
  PC1 = pca$x[,1], PC2 = pca$x[,2], PC3 = pca$x[,3],
  any_3SD_PC1_PC3 = cohort[[id]] %chin% cells3$cell_label,
  any_3SD_PC1_PC80 = cohort[[id]] %chin% cells80$cell_label,
  max_abs_z_PC1_PC80 = apply(abs(score_z[, seq_len(npc80), drop = FALSE]), 1L, max),
  robust_D2_PC80 = robust_d2, robust_D2_cutoff = d2_cutoff,
  robust_D2_flag = robust_d2 > d2_cutoff,
  LOF5_PC80 = lof5, LOF5_gt2 = lof5 > 2
)
diag[, multimethod_flag := any_3SD_PC1_PC80 & (robust_D2_flag | LOF5_gt2)]
setorder(diag, -any_3SD_PC1_PC80, -max_abs_z_PC1_PC80, cell_label)

raw_out <- data.table(
  cell_label = cohort[[id]], donor_label = cohort$donor_label,
  ROI = cohort$Lib_region_of_interest_label, Subclass = cohort$Subclass_name,
  Group_name = cohort$Group_name, raw_x)
shift_out <- data.table(cell_label = cohort[[id]], shifted)
sum_out <- data.table(cell_label = cohort[[id]], sum10000)
log_out <- data.table(cell_label = cohort[[id]], logged)
z_out <- data.table(cell_label = cohort[[id]], z)
scores_out <- data.table(cell_label = cohort[[id]], pca$x)
variance_out <- data.table(
  PC = paste0("PC", seq_along(explained)), eigenvalue = pca$sdev^2,
  explained_variance_percent = explained,
  cumulative_variance_percent = cum_explained,
  included_in_80pct_space = seq_along(explained) <= npc80)
loadings_out <- data.table(feature = rownames(pca$rotation), pca$rotation)
transform_audit <- data.table(
  feature = features, minimum_before_shift = mins,
  shifted_minimum = apply(shifted, 2L, min), sum10000_column_sum = colSums(sum10000),
  z_mean = colMeans(z), z_sample_SD = apply(z, 2L, sd))

write_dt(raw_out, "01_raw_126.csv")
write_dt(shift_out, "02_shifted.csv")
write_dt(sum_out, "03_sum10000.csv")
write_dt(log_out, "04_log1p.csv")
write_dt(z_out, "05_z.csv")
write_dt(scores_out, "06_pca_scores.csv")
write_dt(variance_out, "07_pca_variance.csv")
write_dt(loadings_out, "08_pca_loadings.csv")
write_dt(transform_audit, "09_transform_audit.csv")
write_dt(exceedances, "10_pca_3sd_exceedances_all18PC.csv")
write_dt(diag, "11_pca_outlier_diagnostics_126.csv")
write_dt(cells80, sprintf("12_isolation_candidates_PC1_PC%d.csv", npc80))

# Authoritative diagnostics for downstream scans: only NPC = 3, 4, 5.
allowed_npcs <- 3:5
npc_diagnostics <- rbindlist(lapply(allowed_npcs, function(npc_i) {
  set.seed(777L + npc_i)
  mve_i <- cov.rob(pca$x[, seq_len(npc_i), drop = FALSE], method = "mve")
  d2_i <- mahalanobis(pca$x[, seq_len(npc_i), drop = FALSE], mve_i$center, mve_i$cov)
  cutoff_i <- qchisq(1 - 0.05 / nrow(cohort), df = npc_i)
  lof_i <- lof(pca$x[, seq_len(npc_i), drop = FALSE], minPts = 5L)
  max_z_i <- apply(abs(score_z[, seq_len(npc_i), drop = FALSE]), 1L, max)
  data.table(
    npcs = npc_i, cell_label = cohort[[id]], donor_label = cohort$donor_label,
    ROI = cohort$Lib_region_of_interest_label, Subclass = cohort$Subclass_name,
    max_abs_PC_z = max_z_i, classical_3SD_flag = max_z_i > 3,
    robust_D2 = d2_i, robust_D2_cutoff = cutoff_i, robust_D2_flag = d2_i > cutoff_i,
    LOF5 = lof_i, LOF5_gt2 = lof_i > 2,
    multimethod_flag = max_z_i > 3 & (d2_i > cutoff_i | lof_i > 2)
  )
}))
npc_summary <- npc_diagnostics[, .(
  explained_variance_percent = cum_explained[first(npcs)],
  classical_3SD_n = sum(classical_3SD_flag),
  robust_D2_n = sum(robust_D2_flag),
  LOF5_gt2_n = sum(LOF5_gt2),
  multimethod_n = sum(multimethod_flag),
  classical_3SD_ids = paste(cell_label[classical_3SD_flag], collapse = ";")
), by = npcs]
write_dt(npc_diagnostics, "17_NPC3to5_outlier_diagnostics.csv")
write_dt(npc_summary, "18_NPC3to5_outlier_summary.csv")
for (npc_i in allowed_npcs) {
  write_dt(npc_diagnostics[npcs == npc_i & classical_3SD_flag == TRUE,
                           .(cell_label, donor_label, ROI, Subclass, max_abs_PC_z,
                             robust_D2_flag, LOF5_gt2, multimethod_flag)],
           sprintf("19_isolation_candidates_NPC%d.csv", npc_i))
}

# Explicit NPC3 isolation package with source metadata and transcriptomic identity.
npc3_ids <- npc_diagnostics[npcs == 3 & classical_3SD_flag == TRUE, cell_label]
pc3_trigger <- exceedances[PC_number <= 3, .(
  exceeding_PCs = paste(PC, collapse = ";"),
  PC_z_values = paste(sprintf("%s:%+.4f", PC, classical_z), collapse = ";"),
  max_abs_PC_z = max(abs_z)
), by = cell_label]
tclass_cols <- c(
  "cell_label", "dataset_label", "library_label", "donor_label",
  "Lib_region_of_interest_label", "Neighborhood_name", "Class_name",
  "Subclass_name", "Group_name", "Cluster_label", "cluster_alias",
  "Group_bootstrapping_probability", "Group_aggregate_probability",
  "Cluster_bootstrapping_probability", "Cluster_aggregate_probability"
)
npc3_identity <- merge(cohort[cell_label %chin% npc3_ids, ..tclass_cols], pc3_trigger,
                       by = "cell_label", all.x = TRUE, sort = FALSE)
setorder(npc3_identity, -max_abs_PC_z, cell_label)
write_dt(npc3_identity, "21_NPC3_isolated_with_Tclass_source.csv")
write_dt(raw_out[cell_label %chin% npc3_ids], "22_NPC3_isolated_raw18.csv")
write_dt(raw_out[!cell_label %chin% npc3_ids], "23_NPC3_analysis_cells_after_isolation_119.csv")

# Compare to the previous all-19 QC without using any Mouse output.
m19_manifest_path <- file.path(dirname(out), "qc19", "00_manifest.json")
comparison <- data.table()
if (file.exists(m19_manifest_path)) {
  m19 <- fromJSON(m19_manifest_path)
  old <- as.character(m19$classical_3sd_pc1_pc80_ids)
  new <- cells80$cell_label
  comparison <- data.table(
    category = c("both", "18_only", "19_only"),
    n = c(length(intersect(new, old)), length(setdiff(new, old)), length(setdiff(old, new))),
    cell_labels = c(paste(intersect(new, old), collapse = ";"),
                    paste(setdiff(new, old), collapse = ";"),
                    paste(setdiff(old, new), collapse = ";")))
  write_dt(comparison, "13_outlier_comparison_vs_m19.csv")
}

p_scree <- ggplot(variance_out, aes(x = seq_along(explained), y = explained_variance_percent)) +
  geom_col(fill = "#4C78A8", width = .8) + geom_line(aes(y = cumulative_variance_percent), colour = "#D95F02") +
  geom_point(aes(y = cumulative_variance_percent), colour = "#D95F02", size = 1.4) +
  geom_hline(yintercept = 80, linetype = 2, colour = "#555555") +
  scale_x_continuous(breaks = seq_along(explained)) +
  labs(title = "Macaque 18-feature PCA variance", x = "PC", y = "Percent") + theme_classic(base_size = 9)
ggsave(file.path(out, "14_pca_scree.png"), p_scree, width = 7, height = 4, dpi = 600, bg = "white")
ggsave(file.path(out, "14_pca_scree.pdf"), p_scree, width = 7, height = 4, bg = "white")

plot_dt <- rbindlist(list(
  diag[, .(cell_label, X=PC1, Y=PC2, pair="PC1 vs PC2", any_3SD_PC1_PC80, LOF5_gt2)],
  diag[, .(cell_label, X=PC1, Y=PC3, pair="PC1 vs PC3", any_3SD_PC1_PC80, LOF5_gt2)]))
plot_dt[, status := fifelse(any_3SD_PC1_PC80 & LOF5_gt2, "3SD + LOF",
                            fifelse(any_3SD_PC1_PC80, "Any retained PC >3SD", "Not 3SD"))]
p_out <- ggplot(plot_dt, aes(X, Y, colour = status)) + geom_point(size = 1.7, alpha = .86) +
  geom_text_repel(data = plot_dt[any_3SD_PC1_PC80 == TRUE], aes(label = cell_label),
                  size = 1.8, seed = 777, max.overlaps = Inf, show.legend = FALSE) +
  facet_wrap(~pair, scales = "free", nrow = 1) +
  scale_colour_manual(values = c("Not 3SD"="#B8B8B8", "Any retained PC >3SD"="#E69F00", "3SD + LOF"="#C00000")) +
  labs(title = sprintf("18-feature PCA outliers; PC1-PC%d = %.2f%%", npc80, cum_explained[npc80]),
       x = NULL, y = NULL, colour = NULL) + theme_classic(base_size = 8) + theme(legend.position = "top")
ggsave(file.path(out, "15_pca_outlier_map.png"), p_out, width = 8, height = 3.8, dpi = 600, bg = "white")
ggsave(file.path(out, "15_pca_outlier_map.pdf"), p_out, width = 8, height = 3.8, bg = "white")

npc_plot <- merge(
  npc_diagnostics[, .(npcs, cell_label, classical_3SD_flag, LOF5_gt2)],
  data.table(cell_label = cohort[[id]], PC1 = pca$x[,1], PC2 = pca$x[,2]),
  by = "cell_label", all.x = TRUE, sort = FALSE)
npc_plot[, status := fifelse(classical_3SD_flag & LOF5_gt2, "3SD + LOF",
                             fifelse(classical_3SD_flag, "3SD", "Not 3SD"))]
p_npc <- ggplot(npc_plot, aes(PC1, PC2, colour = status)) +
  geom_point(size = 1.35, alpha = .82) +
  geom_text_repel(data = npc_plot[classical_3SD_flag == TRUE], aes(label = cell_label),
                  size = 1.55, seed = 777, max.overlaps = Inf, show.legend = FALSE) +
  facet_wrap(~npcs, nrow = 1, labeller = label_both) +
  scale_colour_manual(values = c("Not 3SD"="#B8B8B8", "3SD"="#E69F00", "3SD + LOF"="#C00000")) +
  labs(title = "18-feature PCA: NPC-specific outliers used for downstream scans",
       x = "PC1", y = "PC2", colour = NULL) +
  theme_classic(base_size = 8) + theme(legend.position = "top")
ggsave(file.path(out, "20_NPC3to5_outlier_map.png"), p_npc, width = 12, height = 4, dpi = 600, bg = "white")
ggsave(file.path(out, "20_NPC3to5_outlier_map.pdf"), p_npc, width = 12, height = 4, bg = "white")

cor_s <- cor(raw_x, method = "spearman")
pheatmap(cor_s, clustering_distance_rows = "correlation", clustering_distance_cols = "correlation",
         clustering_method = "average", color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(101),
         breaks = seq(-1, 1, length.out = 102), border_color = NA, fontsize = 6.5,
         filename = file.path(out, "16_spearman_heatmap.png"), width = 10, height = 9)

manifest <- list(
  analysis = "Macaque-only 18-feature PCA and outlier detection; 3_Sholl_PC1 excluded",
  source_zip = normalizePath(zip_path, winslash = "/"), source_zip_sha256 = zip_sha,
  mouse_morphology_results_read = FALSE, target_metadata_msn_n = nrow(target),
  complete_morphology_n = nrow(cohort), feature_n = length(features), excluded_feature = "3_Sholl_PC1",
  features = features,
  transformation = "per feature: shift by observed minimum; column-sum normalize to 10000; log1p; sample z-score",
  pca_80pct_npc = npc80, pca_80pct_cumulative_variance_percent = cum_explained[npc80],
  pc1_pc3_variance_percent = sum(explained[1:3]),
  classical_3sd_pc1_pc3_n = nrow(cells3), classical_3sd_pc1_pc80_n = nrow(cells80),
  classical_3sd_pc1_pc3_ids = cells3$cell_label,
  classical_3sd_pc1_pc80_ids = cells80$cell_label,
  robust_D2_flag_n = sum(diag$robust_D2_flag), LOF5_gt2_n = sum(diag$LOF5_gt2),
  multimethod_flag_n = sum(diag$multimethod_flag),
  downstream_allowed_npcs = allowed_npcs,
  downstream_npc_specific_outlier_summary = npc_summary,
  authoritative_downstream_outlier_file = "17_NPC3to5_outlier_diagnostics.csv",
  automatic_exclusion = FALSE)
write_json(manifest, file.path(out, "00_manifest.json"), pretty = TRUE, auto_unbox = TRUE)

cat(sprintf("n=%d p=%d NPC80=%d cumulative=%.3f%% PC1-PC3=%.3f%% 3SD3=%d 3SD80=%d D2=%d LOF=%d\n",
            nrow(cohort), length(features), npc80, cum_explained[npc80], sum(explained[1:3]),
            nrow(cells3), nrow(cells80), sum(diag$robust_D2_flag), sum(diag$LOF5_gt2)))
cat("3SD retained-PC cells:", paste(cells80$cell_label, collapse = "; "), "\n")
if (nrow(comparison)) print(comparison)
cat("Output:", normalizePath(out, winslash = "/"), "\n")
