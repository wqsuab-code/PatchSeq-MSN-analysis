#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

root <- file.path(
  "C:/Users/53461/Documents/Codex/2026-07-13/zh/outputs",
  "05_branch_specific_subtype_mapping"
)
input_root <- file.path(
  "C:/Users/53461/Documents/Codex/2026-07-13/zh/work/05_input",
  "05_Codex_branch_specific_mapping_input"
)
coordinate_dir <- file.path(root, "10_six_mapping_figures")
out_dir <- file.path(root, "74_global_neuron_embedding_RPCA_Root_stability")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

ref <- fread(file.path(coordinate_dir, "global_reference_umap_coordinates.csv"))
query_coord <- fread(file.path(coordinate_dir, "global_query_mapped_umap_coordinates.csv"))
rpca <- fread(file.path(
  input_root,
  "02_global_mapping_results",
  "03B_final_RPCA_CellType_and_Root_combined.csv"
))

stopifnot(nrow(ref) == 21359L, nrow(query_coord) == 641L, nrow(rpca) == 641L)

setkey(query_coord, cell_id)
setkey(rpca, cell_id)
query_plot <- rpca[query_coord]
stopifnot(nrow(query_plot) == 641L, !anyNA(query_plot$UMAP_1), !anyNA(query_plot$UMAP_2))

query_plot[, RPCA_Root_status := fifelse(
  Final_Root_stability_from_CellType %in% c("Stable_MSN", "Stable_IN"),
  "Consensus",
  "MSN-IN disputed"
)]
query_plot[, RPCA_Root_status := factor(
  RPCA_Root_status,
  levels = c("Consensus", "MSN-IN disputed")
)]
setorder(query_plot, RPCA_Root_status)

status_counts <- query_plot[, .N, by = RPCA_Root_status]
stopifnot(
  status_counts[RPCA_Root_status == "Consensus", N] == 622L,
  status_counts[RPCA_Root_status == "MSN-IN disputed", N] == 19L
)

legend_labels <- c(
  Consensus = "Consensus (n=622)",
  `MSN-IN disputed` = "MSN-IN disputed (n=19)"
)

p <- ggplot() +
  geom_point(
    data = ref,
    aes(x = UMAP_1, y = UMAP_2, color = label),
    shape = 16,
    size = 0.28,
    alpha = 0.78
  ) +
  geom_point(
    data = query_plot,
    aes(
      x = UMAP_1,
      y = UMAP_2,
      shape = RPCA_Root_status,
      fill = RPCA_Root_status
    ),
    size = 0.92,
    stroke = 0.38,
    alpha = 1,
    color = "#111111"
  ) +
  scale_color_manual(
    values = c(D1 = "#FF1616", D2 = "#1027E8", IN = "#00D728"),
    breaks = c("D1", "D2", "IN"),
    drop = FALSE
  ) +
  scale_shape_manual(
    values = c(Consensus = 16, `MSN-IN disputed` = 21),
    labels = legend_labels,
    drop = FALSE
  ) +
  scale_fill_manual(
    values = c(Consensus = "#111111", `MSN-IN disputed` = "#FFFFFF"),
    labels = legend_labels,
    drop = FALSE
  ) +
  coord_fixed() +
  labs(
    title = "Global CellType B",
    x = "UMAP_1",
    y = "UMAP_2",
    color = "Reference",
    shape = "Query",
    fill = "Query"
  ) +
  theme_classic(base_size = 8) +
  theme(
    plot.title = element_text(size = 9.5, face = "plain", hjust = 0),
    axis.title = element_text(size = 8.5, color = "black"),
    axis.text = element_text(size = 7.2, color = "black"),
    axis.line = element_line(linewidth = 0.42, color = "black"),
    axis.ticks = element_line(linewidth = 0.38, color = "black"),
    legend.position = "right",
    legend.box = "vertical",
    legend.direction = "vertical",
    legend.justification = "center",
    legend.title = element_text(size = 7.4, face = "bold", color = "black"),
    legend.text = element_text(size = 7.2, color = "black"),
    legend.key.width = unit(8, "pt"),
    legend.key.height = unit(8, "pt"),
    legend.spacing.x = unit(3, "pt"),
    legend.margin = margin(1, 0, 0, 0, unit = "pt"),
    plot.margin = margin(2, 2, 2, 2, unit = "pt")
  ) +
  guides(
    color = guide_legend(order = 1, override.aes = list(size = 1.8, alpha = 1)),
    shape = guide_legend(order = 2, override.aes = list(size = 2.0, alpha = 1)),
    fill = "none"
  )

png_file <- file.path(out_dir, "Figure_GlobalCellType_B_Query_RPCA_Root_consensus_disputed.png")
pdf_file <- file.path(out_dir, "Figure_GlobalCellType_B_Query_RPCA_Root_consensus_disputed.pdf")
data_file <- file.path(out_dir, "GlobalCellType_B_Query_RPCA_Root_consensus_disputed_plot_data.csv")
count_file <- file.path(out_dir, "GlobalCellType_B_Query_RPCA_Root_consensus_disputed_counts.csv")

ggsave(png_file, p, width = 5, height = 4, dpi = 600, bg = "white")
ggsave(pdf_file, p, width = 5, height = 4, device = cairo_pdf, bg = "white")
fwrite(query_plot, data_file)
fwrite(status_counts, count_file)

writeLines(c(
  "Global CellType B overlay using the original RPCA cross-nPC Root-stability definition.",
  "",
  "Reference: 21,359 neurons colored by D1, D2 and IN CellType.",
  "Solid black Query: Stable_MSN or Stable_IN across the included RPCA PC settings (588 + 34 = 622).",
  "Hollow Query: Fluctuating_Root between MSN and IN across RPCA PC settings (n=19).",
  "This plot does not use Pearson labels to define point shape."
), file.path(out_dir, "README_GlobalCellType_B_Query_RPCA_Root_stability.txt"))

message("PNG: ", png_file)
message("PDF: ", pdf_file)
