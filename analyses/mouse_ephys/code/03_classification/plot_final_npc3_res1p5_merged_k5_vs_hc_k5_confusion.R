#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(mclust)
  library(clue)
})

options(stringsAsFactors = FALSE)

input_file <- file.path(
  "outputs", "e_type_qc",
  "fixed_npc3_res1.5_merged_k5_tsne",
  "Fixed_NPC3_Res1.5_Merged_K5_tSNE_Coordinates_Assignments.csv"
)
output_dir <- file.path(
  "outputs", "e_type_qc",
  "final_npc3_res1p5_merged_k5_vs_hc_k5_confusion"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

data <- read.csv(
  input_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

stopifnot(
  nrow(data) == 494L,
  length(unique(data$MSN_unique_ID)) == 494L,
  all(c(
    "Merged_Seurat_cluster",
    "HC_NPC3_K5"
  ) %in% names(data))
)

seurat_levels <- paste0("S-", 1:5)
hc_levels <- paste0("HC", 1:5)
seurat_labels <- factor(
  data$Merged_Seurat_cluster,
  levels = seurat_levels
)
hc_labels <- factor(
  paste0("HC", data$HC_NPC3_K5),
  levels = hc_levels
)

contingency_raw <- table(
  Seurat = seurat_labels,
  HC = hc_labels
)

hungarian_assignment <- as.integer(
  clue::solve_LSAP(contingency_raw, maximum = TRUE)
)
matched_count <- sum(
  contingency_raw[
    cbind(seq_len(nrow(contingency_raw)), hungarian_assignment)
  ]
)
matched_percent <- 100 * matched_count / sum(contingency_raw)

# Reorder HC columns to put the optimal one-to-one matches on the diagonal,
# while retaining the original HC identifiers in the axis labels.
hc_order <- colnames(contingency_raw)[hungarian_assignment]
contingency <- contingency_raw[, hc_order, drop = FALSE]
row_fraction <- prop.table(contingency, margin = 1L)

normalized_mutual_information <- function(a, b) {
  contingency_nmi <- table(a, b)
  pxy <- contingency_nmi / sum(contingency_nmi)
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

ari <- mclust::adjustedRandIndex(seurat_labels, hc_labels)
nmi <- normalized_mutual_information(seurat_labels, hc_labels)
seurat_to_hc_purity <- sum(apply(contingency_raw, 1L, max)) /
  sum(contingency_raw)
hc_to_seurat_purity <- sum(apply(contingency_raw, 2L, max)) /
  sum(contingency_raw)

long <- as.data.frame(contingency)
names(long) <- c(
  "Seurat_cluster",
  "HC_cluster",
  "Cell_count"
)
long$Row_fraction <- mapply(
  function(seurat_cluster, hc_cluster) {
    row_fraction[seurat_cluster, hc_cluster]
  },
  long$Seurat_cluster,
  long$HC_cluster
)
long$Matched_pair <- mapply(
  function(seurat_cluster, hc_cluster) {
    row_index <- match(seurat_cluster, seurat_levels)
    hc_cluster == hc_order[row_index]
  },
  as.character(long$Seurat_cluster),
  as.character(long$HC_cluster)
)
long$Label <- sprintf(
  "%d\n%.0f%%",
  long$Cell_count,
  100 * long$Row_fraction
)
long$Text_colour <- ifelse(
  long$Row_fraction >= 0.50,
  "white",
  "black"
)
long$Seurat_display <- factor(
  long$Seurat_cluster,
  levels = rev(seurat_levels)
)
long$HC_display <- factor(
  long$HC_cluster,
  levels = hc_order
)

plot_object <- ggplot(
  long,
  aes(
    x = HC_display,
    y = Seurat_display,
    fill = Row_fraction
  )
) +
  geom_tile(
    colour = "grey84",
    linewidth = 0.32
  ) +
  geom_tile(
    data = long[long$Matched_pair, , drop = FALSE],
    fill = NA,
    colour = "black",
    linewidth = 0.80
  ) +
  geom_text(
    aes(
      label = Label,
      colour = Text_colour
    ),
    size = 2.25,
    lineheight = 0.84,
    family = "Arial"
  ) +
  scale_colour_identity() +
  scale_fill_gradientn(
    colours = c(
      "#FFFFFF",
      "#DCEAF4",
      "#74ADD1",
      "#2166AC",
      "#053061"
    ),
    limits = c(0, 1),
    breaks = c(0, 0.5, 1),
    labels = scales::percent_format(accuracy = 1),
    name = "Within-Seurat\nrow fraction"
  ) +
  guides(
    fill = guide_colorbar(
      title.position = "left",
      title.hjust = 1,
      label.position = "bottom",
      barwidth = grid::unit(36, "mm"),
      barheight = grid::unit(2.5, "mm"),
      ticks = TRUE,
      frame.colour = "grey55",
      frame.linewidth = 0.30
    )
  ) +
  scale_x_discrete(position = "top") +
  coord_fixed() +
  labs(
    title = "NPC3 Ward.D2 HC K=5 vs Seurat res=1.5 merged K=5",
    subtitle = sprintf(
      paste0(
        "N=494 | consensus %d (%.1f%%) | ARI %.3f | NMI %.3f\n",
        "HC columns reordered by Hungarian matching; original HC IDs retained"
      ),
      matched_count,
      matched_percent,
      ari,
      nmi
    ),
    x = "Ward.D2 HC cluster",
    y = "Merged Seurat cluster"
  ) +
  theme_classic(
    base_size = 7,
    base_family = "Arial"
  ) +
  theme(
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_text(
      size = 6.5,
      colour = "black"
    ),
    axis.title = element_text(
      size = 6.5,
      colour = "black"
    ),
    plot.title = element_text(
      face = "bold",
      hjust = 0.5,
      size = 8.2,
      margin = margin(b = 1.5)
    ),
    plot.subtitle = element_text(
      hjust = 0.5,
      size = 6.2,
      colour = "grey25",
      lineheight = 0.95,
      margin = margin(b = 4)
    ),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.title = element_text(size = 6.0),
    legend.text = element_text(size = 5.8),
    panel.border = element_rect(
      colour = "grey65",
      fill = NA,
      linewidth = 0.35
    ),
    plot.margin = margin(5, 5, 5, 5)
  )

stem <- "NPC3_Final18_Seurat_res1p5_MergedK5_vs_WardD2_HC_K5_Confusion"
png_file <- file.path(output_dir, paste0(stem, ".png"))
pdf_file <- file.path(output_dir, paste0(stem, ".pdf"))
ggsave(
  png_file,
  plot_object,
  width = 3.25,
  height = 3.05,
  units = "in",
  dpi = 600,
  bg = "white"
)
ggsave(
  pdf_file,
  plot_object,
  width = 3.25,
  height = 3.05,
  units = "in",
  device = cairo_pdf,
  bg = "white"
)

counts_output <- data.frame(
  Seurat_cluster = rownames(contingency),
  as.data.frame.matrix(contingency),
  check.names = FALSE
)
row_percent_output <- data.frame(
  Seurat_cluster = rownames(row_fraction),
  100 * as.data.frame.matrix(row_fraction),
  check.names = FALSE
)
mapping_output <- data.frame(
  Seurat_cluster = seurat_levels,
  Matched_HC_cluster = hc_order,
  Matched_cell_count = contingency[
    cbind(seq_len(nrow(contingency)), seq_len(ncol(contingency)))
  ],
  Seurat_cluster_size = rowSums(contingency),
  Within_Seurat_match_percent = 100 *
    contingency[
      cbind(seq_len(nrow(contingency)), seq_len(ncol(contingency)))
    ] /
    rowSums(contingency),
  stringsAsFactors = FALSE
)
metrics_output <- data.frame(
  N = nrow(data),
  Feature_set = "Final18",
  NPCs = 3L,
  Seurat_resolution = 1.5,
  Seurat_raw_K = 11L,
  Seurat_merged_K = 5L,
  HC_K = 5L,
  Consensus_n = matched_count,
  Consensus_percent = matched_percent,
  ARI = ari,
  NMI = nmi,
  Seurat_to_HC_purity = seurat_to_hc_purity,
  HC_to_Seurat_purity = hc_to_seurat_purity,
  stringsAsFactors = FALSE
)

counts_file <- file.path(output_dir, paste0(stem, "_Counts.csv"))
row_percent_file <- file.path(
  output_dir,
  paste0(stem, "_RowPercent.csv")
)
mapping_file <- file.path(
  output_dir,
  paste0(stem, "_Hungarian_Mapping.csv")
)
metrics_file <- file.path(
  output_dir,
  paste0(stem, "_Metrics.csv")
)
write.csv(counts_output, counts_file, row.names = FALSE)
write.csv(row_percent_output, row_percent_file, row.names = FALSE)
write.csv(mapping_output, mapping_file, row.names = FALSE)
write.csv(metrics_output, metrics_file, row.names = FALSE)

summary_lines <- c(
  "Final NPC3/res=1.5 merged Seurat K5 versus Ward.D2 HC K5",
  "Feature set: Final18; high-skew metrics excluded",
  "Cells: 494",
  sprintf(
    "Consensus: %d/494 (%.2f%%)",
    matched_count,
    matched_percent
  ),
  sprintf("ARI: %.6f", ari),
  sprintf("NMI: %.6f", nmi),
  sprintf("Seurat-to-HC purity: %.6f", seurat_to_hc_purity),
  sprintf("HC-to-Seurat purity: %.6f", hc_to_seurat_purity),
  paste(
    "Hungarian mapping:",
    paste(
      paste0(
        mapping_output$Seurat_cluster,
        "->",
        mapping_output$Matched_HC_cluster
      ),
      collapse = "; "
    )
  ),
  paste(
    "PNG:",
    normalizePath(png_file, winslash = "/", mustWork = TRUE)
  ),
  paste(
    "PDF:",
    normalizePath(pdf_file, winslash = "/", mustWork = TRUE)
  )
)
log_file <- file.path(output_dir, paste0(stem, "_run_log.txt"))
writeLines(summary_lines, log_file, useBytes = TRUE)
cat(paste(summary_lines, collapse = "\n"), "\n")
