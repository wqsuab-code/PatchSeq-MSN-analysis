#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(mclust)
  library(Seurat)
})

seed <- 20260825L
input_file <- file.path(
  "outputs", "morph_qc", "morph187_after_incomplete_reconstruction_quarantine",
  "01_morph187_discovery_input.csv"
)
out <- file.path(
  "outputs", "morph_qc", "morph187_NPC3-4_HCK4-5_GC_allres_raw"
)
dir.create(out, recursive = TRUE, showWarnings = FALSE)

id <- "MSN_unique_ID"
features <- c(
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
npc_grid <- 3:4
hc_grid <- 4:5
res_grid <- round(seq(0.50, 3.00, by = 0.05), 2)

dat <- fread(input_file)
stopifnot(nrow(dat) == 187L, uniqueN(dat[[id]]) == 187L, all(features %in% names(dat)))
x <- as.matrix(dat[, ..features])
stopifnot(all(is.finite(x)))

# Frozen transcriptomic-like Morph preprocessing, refit after quarantine.
shifted <- sweep(x, 2L, apply(x, 2L, min), FUN = "-")
totals <- colSums(shifted)
stopifnot(all(totals > 0))
logged <- log1p(sweep(shifted, 2L, totals, FUN = "/") * 10000)
z <- scale(logged, center = TRUE, scale = TRUE)
pca <- prcomp(z, center = FALSE, scale. = FALSE)
scores <- as.data.table(pca$x)
scores[, (id) := dat[[id]]]
setcolorder(scores, c(id, paste0("PC", seq_len(ncol(pca$x)))))

variance <- data.table(
  PC = seq_along(pca$sdev),
  Explained_variance_percent = 100 * pca$sdev^2 / sum(pca$sdev^2),
  Cumulative_variance_percent = 100 * cumsum(pca$sdev^2) / sum(pca$sdev^2)
)
loadings <- as.data.table(pca$rotation, keep.rownames = "Feature")
fwrite(scores, file.path(out, "01_PCA_scores_187cells.csv"))
fwrite(variance, file.path(out, "02_PCA_variance.csv"))
fwrite(loadings, file.path(out, "03_PCA_loadings.csv"))

nmi <- function(a, b) {
  tab <- table(as.integer(as.factor(a)), as.integer(as.factor(b)))
  pxy <- tab / sum(tab); px <- rowSums(pxy); py <- colSums(pxy)
  nz <- which(pxy > 0, arr.ind = TRUE)
  mi <- sum(pxy[nz] * log(pxy[nz] / (px[nz[, 1L]] * py[nz[, 2L]])))
  hx <- -sum(px[px > 0] * log(px[px > 0])); hy <- -sum(py[py > 0] * log(py[py > 0]))
  if (hx <= 0 || hy <= 0) return(NA_real_)
  mi / sqrt(hx * hy)
}

relabel_hc <- function(raw, tree) {
  leaf <- tree$labels[tree$order]
  map <- setNames(seq_along(unique(raw[leaf])), unique(raw[leaf]))
  result <- as.integer(map[as.character(raw)])
  names(result) <- names(raw)
  result
}

relabel_gc <- function(raw, pc1) {
  raw <- as.integer(as.factor(raw))
  centers <- tapply(pc1, raw, median)
  map <- setNames(seq_along(centers), names(sort(centers)))
  result <- as.integer(map[as.character(raw)])
  names(result) <- names(pc1)
  result
}

hc_assignments <- list(); gc_assignments <- list(); counts <- list(); metrics <- list()
hi <- gi <- ci <- mi <- 0L

for (npc in npc_grid) {
  coords <- as.matrix(scores[, paste0("PC", seq_len(npc)), with = FALSE])
  rownames(coords) <- scores[[id]]
  tree <- hclust(dist(coords), method = "ward.D2")
  hc_store <- list()
  for (k in hc_grid) {
    hc <- relabel_hc(cutree(tree, k = k), tree)
    hc_store[[as.character(k)]] <- hc
    hi <- hi + 1L
    hc_assignments[[hi]] <- data.table(
      MSN_unique_ID = names(hc), NPC = npc, HC_K = k,
      HC_cluster = paste0("HC", hc)
    )
  }

  snn <- FindNeighbors(
    coords, k.param = 20L, compute.SNN = TRUE,
    prune.SNN = 1 / 15, verbose = FALSE
  )[["snn"]]

  for (resolution in res_grid) {
    set.seed(seed)
    fit <- FindClusters(
      snn, algorithm = 1, resolution = resolution,
      random.seed = seed, n.start = 30, n.iter = 30, verbose = FALSE
    )
    raw <- as.integer(as.factor(fit[, ncol(fit)]))
    names(raw) <- rownames(coords)
    gc <- relabel_gc(raw, coords[, 1])
    raw_k <- uniqueN(gc)
    gi <- gi + 1L
    gc_assignments[[gi]] <- data.table(
      MSN_unique_ID = names(gc), NPC = npc, Resolution = resolution,
      GC_raw_K = raw_k, GC_raw_cluster = paste0("S", gc - 1L)
    )

    for (k in hc_grid) {
      hc <- hc_store[[as.character(k)]]
      tab <- table(
        HC = factor(paste0("HC", hc), levels = paste0("HC", seq_len(k))),
        GC_raw = factor(paste0("S", gc - 1L), levels = paste0("S", 0:(raw_k - 1L)))
      )
      long <- as.data.table(as.table(tab)); setnames(long, "N", "Cell_count")
      long[, `:=`(NPC = npc, HC_K = k, Resolution = resolution, GC_raw_K = raw_k)]
      ci <- ci + 1L; counts[[ci]] <- long
      mi <- mi + 1L
      metrics[[mi]] <- data.table(
        NPC = npc, HC_K = k, Resolution = resolution, GC_raw_K = raw_k,
        Raw_ARI = adjustedRandIndex(hc, gc),
        Raw_NMI = nmi(hc, gc),
        Raw_GC_to_HC_purity = sum(apply(tab, 2L, max)) / length(hc),
        Raw_HC_to_GC_purity = sum(apply(tab, 1L, max)) / length(hc),
        HC_sizes = paste(as.integer(table(hc)), collapse = ";"),
        GC_raw_sizes = paste(as.integer(table(gc)), collapse = ";")
      )
    }
  }
}

hc_dt <- rbindlist(hc_assignments)
gc_dt <- rbindlist(gc_assignments)
count_dt <- rbindlist(counts)
metric_dt <- rbindlist(metrics)

best_any <- metric_dt[order(NPC, HC_K, -Raw_ARI, -Raw_NMI), .SD[1], by = .(NPC, HC_K)]
best_equal <- metric_dt[GC_raw_K == HC_K][order(NPC, HC_K, -Raw_ARI, -Raw_NMI), .SD[1], by = .(NPC, HC_K)]
best_purity <- metric_dt[order(NPC, HC_K, -Raw_GC_to_HC_purity, -Raw_ARI), .SD[1], by = .(NPC, HC_K)]

fwrite(hc_dt, file.path(out, "04_HC_NPC3-4_K4-5_assignments.csv"))
fwrite(gc_dt, file.path(out, "05_GC_NPC3-4_allres_raw_assignments.csv"))
fwrite(count_dt, file.path(out, "06_all_raw_cross_confusion_counts_long.csv"))
fwrite(metric_dt, file.path(out, "07_all_NPC_HCK_resolution_raw_metrics.csv"))
fwrite(best_any, file.path(out, "08_best_raw_ARI_each_NPC_HCK.csv"))
fwrite(best_equal, file.path(out, "09_best_equal_K_each_NPC_HCK.csv"))
fwrite(best_purity, file.path(out, "10_best_directional_purity_each_NPC_HCK.csv"))

metric_plot <- melt(
  metric_dt,
  id.vars = c("NPC", "HC_K", "Resolution", "GC_raw_K"),
  measure.vars = c("Raw_ARI", "Raw_NMI", "Raw_GC_to_HC_purity"),
  variable.name = "Metric", value.name = "Value"
)
metric_plot[, Panel := sprintf("NPC=%d | HC K=%d", NPC, HC_K)]
p <- ggplot(metric_plot, aes(Resolution, Value, colour = Metric)) +
  geom_line(linewidth = 0.45) + geom_point(size = 0.5) +
  facet_wrap(~Panel, ncol = 2) +
  scale_colour_manual(
    values = c(Raw_ARI = "#00468B", Raw_NMI = "#42B540", Raw_GC_to_HC_purity = "#ED0000"),
    labels = c("Raw ARI", "Raw NMI", "GC-to-HC purity")
  ) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(x = "GC resolution", y = "Unmerged HC-GC concordance", colour = NULL,
       title = "Morph n=187: raw GC resolution scan") +
  theme_classic(base_family = "sans", base_size = 7) +
  theme(strip.background = element_blank(), strip.text = element_text(face = "bold"),
        legend.position = "top", plot.title = element_text(face = "bold"))
ggsave(file.path(out, "11_raw_concordance_across_all_resolutions.png"), p,
       width = 6.2, height = 4.7, units = "in", dpi = 600, bg = "white")
ggsave(file.path(out, "11_raw_concordance_across_all_resolutions.pdf"), p,
       width = 6.2, height = 4.7, units = "in", bg = "white")

cat("Best raw ARI at every NPC/HC K:\n")
print(best_any[, .(NPC, HC_K, Resolution, GC_raw_K, Raw_ARI, Raw_NMI,
                   Raw_GC_to_HC_purity, Raw_HC_to_GC_purity, HC_sizes, GC_raw_sizes)])
cat("\nBest raw ARI restricted to equal raw K:\n")
print(best_equal[, .(NPC, HC_K, Resolution, GC_raw_K, Raw_ARI, Raw_NMI,
                     Raw_GC_to_HC_purity, HC_sizes, GC_raw_sizes)])
