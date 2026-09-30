#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(Rtsne)
  library(ggplot2)
  library(ggrepel)
  library(mclust)
})

options(stringsAsFactors = FALSE)

seed <- 20260722L
npcs <- 3L
resolution <- 1.5
perplexity <- 30

input_file <- file.path(
  "outputs", "e_type_qc", "pca_shift_sum_final18",
  "ShiftSumNorm10000_Log1p_Zscore_Final18_PCA_input.csv"
)
output_dir <- file.path(
  "outputs", "e_type_qc",
  "fixed_npc3_res1.5_merged_k5_tsne"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

input_data <- read.csv(
  input_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
id_column <- "MSN_unique_ID"
feature_matrix <- as.matrix(
  input_data[, setdiff(names(input_data), id_column), drop = FALSE]
)
storage.mode(feature_matrix) <- "double"
rownames(feature_matrix) <- input_data[[id_column]]

if (
  nrow(feature_matrix) != 494L ||
    ncol(feature_matrix) != 18L ||
    anyNA(feature_matrix) ||
    any(!is.finite(feature_matrix))
) {
  stop("Unexpected final E-feature matrix.")
}

pca_fit <- prcomp(feature_matrix, center = FALSE, scale. = FALSE)
pcs <- pca_fit$x[, seq_len(npcs), drop = FALSE]
rownames(pcs) <- input_data[[id_column]]

snn <- Seurat::FindNeighbors(
  pcs,
  k.param = 20L,
  compute.SNN = TRUE,
  prune.SNN = 1 / 15,
  verbose = FALSE
)[["snn"]]
seurat_result <- Seurat::FindClusters(
  snn,
  algorithm = 1,
  resolution = resolution,
  random.seed = seed,
  n.start = 20,
  n.iter = 20,
  verbose = FALSE
)
raw_cluster <- as.integer(
  as.factor(seurat_result[, ncol(seurat_result)])
)
names(raw_cluster) <- rownames(pcs)

if (!identical(sort(unique(raw_cluster)), 1:11)) {
  stop(
    "Expected 11 raw Seurat clusters at NPC3/resolution 1.5; observed: ",
    paste(sort(unique(raw_cluster)), collapse = ", ")
  )
}

merge_map <- c(
  "1" = "S-1",
  "3" = "S-1",
  "9" = "S-1",
  "2" = "S-2",
  "5" = "S-2",
  "6" = "S-2",
  "11" = "S-2",
  "4" = "S-3",
  "8" = "S-4",
  "10" = "S-4",
  "7" = "S-5"
)
merged_levels <- paste0("S-", 1:5)
merged_cluster <- factor(
  unname(merge_map[as.character(raw_cluster)]),
  levels = merged_levels
)
if (anyNA(merged_cluster)) {
  stop("At least one raw Seurat cluster was not assigned by the merge map.")
}

set.seed(seed)
tsne_fit <- Rtsne::Rtsne(
  pcs,
  dims = 2,
  perplexity = perplexity,
  theta = 0.5,
  check_duplicates = FALSE,
  pca = FALSE,
  normalize = TRUE,
  max_iter = 1000,
  verbose = FALSE
)

plot_data <- data.frame(
  MSN_unique_ID = rownames(pcs),
  tSNE_1 = tsne_fit$Y[, 1],
  tSNE_2 = tsne_fit$Y[, 2],
  Raw_Seurat_cluster = paste0("S", raw_cluster),
  Merged_Seurat_cluster = merged_cluster,
  stringsAsFactors = FALSE
)

cluster_sizes <- table(plot_data$Merged_Seurat_cluster)
legend_labels <- paste0(
  merged_levels,
  " (n=",
  as.integer(cluster_sizes[merged_levels]),
  ")"
)
names(legend_labels) <- merged_levels

cluster_colours <- c(
  "S-1" = "#0072B2",
  "S-2" = "#E69F00",
  "S-3" = "#009E73",
  "S-4" = "#CC79A7",
  "S-5" = "#D55E00"
)

centres <- aggregate(
  cbind(tSNE_1, tSNE_2) ~ Merged_Seurat_cluster,
  data = plot_data,
  FUN = median
)

tsne_plot <- ggplot(
  plot_data,
  aes(
    x = tSNE_1,
    y = tSNE_2,
    fill = Merged_Seurat_cluster
  )
) +
  geom_point(
    shape = 21,
    size = 1.25,
    stroke = 0.18,
    colour = "white",
    alpha = 0.84
  ) +
  geom_label(
    data = centres,
    aes(
      x = tSNE_1,
      y = tSNE_2,
      label = Merged_Seurat_cluster,
      colour = Merged_Seurat_cluster
    ),
    inherit.aes = FALSE,
    fill = scales::alpha("white", 0.78),
    linewidth = 0.18,
    label.padding = grid::unit(0.08, "lines"),
    size = 2.2,
    fontface = "bold",
    show.legend = FALSE
  ) +
  scale_fill_manual(
    values = cluster_colours,
    breaks = merged_levels,
    labels = legend_labels,
    name = NULL
  ) +
  scale_colour_manual(
    values = cluster_colours,
    breaks = merged_levels,
    guide = "none"
  ) +
  coord_equal() +
  labs(
    title = "Merged E-type clusters",
    subtitle = "NPCs=3 | Seurat resolution=1.5 | t-SNE perplexity=30",
    x = "t-SNE 1",
    y = "t-SNE 2"
  ) +
  theme_classic(base_size = 7, base_family = "sans") +
  theme(
    axis.line = element_line(colour = "black", linewidth = 0.5),
    axis.ticks = element_line(colour = "black", linewidth = 0.5),
    axis.ticks.length = grid::unit(1.5, "mm"),
    axis.text = element_text(colour = "black", size = 6.2),
    axis.title = element_text(colour = "black", size = 7),
    plot.title = element_text(face = "bold", size = 8.5),
    plot.subtitle = element_text(colour = "grey25", size = 6.4),
    legend.position = "right",
    legend.text = element_text(size = 6.2),
    legend.key.height = grid::unit(3.2, "mm"),
    legend.key.width = grid::unit(3.2, "mm"),
    plot.margin = margin(4, 4, 4, 4)
  )

png_file <- file.path(
  output_dir,
  "Fixed_NPC3_Res1.5_Merged_K5_tSNE_W3_H3.png"
)
pdf_file <- file.path(
  output_dir,
  "Fixed_NPC3_Res1.5_Merged_K5_tSNE_W3_H3.pdf"
)
ggsave(
  png_file,
  tsne_plot,
  width = 3,
  height = 3,
  units = "in",
  dpi = 600,
  bg = "white"
)
ggsave(
  pdf_file,
  tsne_plot,
  width = 3,
  height = 3,
  units = "in",
  device = cairo_pdf,
  bg = "white"
)

hc_k5 <- cutree(
  hclust(dist(pcs), method = "ward.D2"),
  k = 5L
)
merged_integer <- as.integer(merged_cluster)
ari_vs_hc <- mclust::adjustedRandIndex(merged_integer, hc_k5)

normalized_mutual_information <- function(a, b) {
  a <- as.integer(as.factor(a))
  b <- as.integer(as.factor(b))
  contingency <- table(a, b)
  pxy <- contingency / sum(contingency)
  px <- rowSums(pxy)
  py <- colSums(pxy)
  nonzero <- which(pxy > 0, arr.ind = TRUE)
  mutual_information <- sum(
    pxy[nonzero] *
      log(
        pxy[nonzero] /
          (px[nonzero[, 1L]] * py[nonzero[, 2L]])
      )
  )
  entropy_x <- -sum(px[px > 0] * log(px[px > 0]))
  entropy_y <- -sum(py[py > 0] * log(py[py > 0]))
  mutual_information / sqrt(entropy_x * entropy_y)
}
nmi_vs_hc <- normalized_mutual_information(merged_integer, hc_k5)

merged_to_hc <- c(
  "S-1" = 4L,
  "S-2" = 2L,
  "S-3" = 1L,
  "S-4" = 3L,
  "S-5" = 5L
)
matched_hc_label <- unname(
  merged_to_hc[as.character(merged_cluster)]
)
matched_accuracy <- mean(matched_hc_label == hc_k5)

plot_data$HC_NPC3_K5 <- hc_k5[plot_data$MSN_unique_ID]
plot_data$Merged_cluster_matched_HC_label <- matched_hc_label
write.csv(
  plot_data,
  file.path(
    output_dir,
    "Fixed_NPC3_Res1.5_Merged_K5_tSNE_Coordinates_Assignments.csv"
  ),
  row.names = FALSE
)

merge_table <- data.frame(
  Merged_Seurat_cluster = merged_levels,
  Raw_Seurat_clusters = c(
    "S1; S3; S9",
    "S2; S5; S6; S11",
    "S4",
    "S8; S10",
    "S7"
  ),
  Cell_count = as.integer(cluster_sizes[merged_levels]),
  Matched_HC_NPC3_K5_cluster =
    as.integer(merged_to_hc[merged_levels]),
  stringsAsFactors = FALSE
)
write.csv(
  merge_table,
  file.path(
    output_dir,
    "Fixed_NPC3_Res1.5_Seurat_Merge_Map.csv"
  ),
  row.names = FALSE
)

confusion <- table(
  Merged_Seurat_cluster = merged_cluster,
  HC_NPC3_K5 = factor(hc_k5, levels = 1:5)
)
confusion_output <- data.frame(
  Merged_Seurat_cluster = rownames(confusion),
  as.data.frame.matrix(confusion),
  check.names = FALSE
)
names(confusion_output)[-1] <- paste0("HC", 1:5)
write.csv(
  confusion_output,
  file.path(
    output_dir,
    "Fixed_NPC3_Res1.5_Merged_K5_vs_HC_K5_Confusion.csv"
  ),
  row.names = FALSE
)

saveRDS(
  list(
    input_file = normalizePath(input_file, winslash = "/"),
    pca_fit = pca_fit,
    pcs = pcs,
    raw_seurat_cluster = raw_cluster,
    merge_map = merge_map,
    merged_seurat_cluster = merged_cluster,
    tsne_fit = tsne_fit,
    plot_data = plot_data,
    merge_table = merge_table,
    hc_k5 = hc_k5,
    ari_vs_hc = ari_vs_hc,
    nmi_vs_hc = nmi_vs_hc,
    matched_accuracy = matched_accuracy
  ),
  file.path(
    output_dir,
    "Fixed_NPC3_Res1.5_Merged_K5_tSNE_Working_Object.rds"
  )
)

run_log <- c(
  "Fixed NPC3/resolution 1.5 merged Seurat K5 t-SNE",
  sprintf("Cells: %d", nrow(feature_matrix)),
  sprintf("Final E features: %d", ncol(feature_matrix)),
  "Preprocessing: shift -> sum normalization x10,000 -> log1p -> z-score",
  "PCA: first 3 PCs",
  "Seurat: k.param=20, prune.SNN=1/15, Louvain algorithm=1",
  sprintf("Resolution: %.1f", resolution),
  sprintf("Random seed: %d", seed),
  "Merge map:",
  "  S-1 = S1 + S3 + S9",
  "  S-2 = S2 + S5 + S6 + S11",
  "  S-3 = S4",
  "  S-4 = S8 + S10",
  "  S-5 = S7",
  sprintf(
    "Merged cluster sizes: %s",
    paste0(
      merged_levels,
      "=",
      as.integer(cluster_sizes[merged_levels]),
      collapse = "; "
    )
  ),
  sprintf(
    "t-SNE: Rtsne on first 3 PCs, perplexity=%d, theta=0.5, max_iter=1000",
    perplexity
  ),
  sprintf("Merged K5 vs HC K5 ARI: %.6f", ari_vs_hc),
  sprintf("Merged K5 vs HC K5 NMI: %.6f", nmi_vs_hc),
  sprintf(
    "Matched-label agreement vs HC K5: %.2f%%",
    100 * matched_accuracy
  ),
  sprintf(
    "PNG: %s",
    normalizePath(png_file, winslash = "/", mustWork = TRUE)
  ),
  sprintf(
    "PDF: %s",
    normalizePath(pdf_file, winslash = "/", mustWork = TRUE)
  ),
  sprintf(
    "Output directory: %s",
    normalizePath(output_dir, winslash = "/", mustWork = TRUE)
  )
)
writeLines(
  run_log,
  file.path(
    output_dir,
    "Fixed_NPC3_Res1.5_Merged_K5_tSNE_run_log.txt"
  ),
  useBytes = TRUE
)
cat(paste(run_log, collapse = "\n"), "\n")
