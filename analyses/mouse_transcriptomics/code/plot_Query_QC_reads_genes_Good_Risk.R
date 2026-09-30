#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(SeuratObject)
  library(ggplot2)
  library(scales)
})

input_object <- file.path(
  "C:/Users/53461/Documents/Codex/2026-07-13/zh/work/05_input",
  "05_Codex_branch_specific_mapping_input",
  "00_original_inputs",
  "query_geneBody_symbol_raw_QC_no39cells_noO9273_Batch2_noO9268_Batch2_n641.rds"
)
pearson_qc_file <- "C:/Users/53461/Downloads/Top15_risk_x040_x050_QC_table.csv"
npc5_qc_file <- "C:/Users/53461/Downloads/Top15CleanMarker_MaxPC_NPC5_QC_table.csv"
out_dir <- file.path("outputs", "T-type_QC_RiskCells")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

query <- readRDS(input_object)
meta <- query[[]]
meta$Query_cell <- rownames(meta)

pearson_qc <- read.csv(pearson_qc_file, check.names = FALSE, stringsAsFactors = FALSE)
npc5_qc <- read.csv(npc5_qc_file, check.names = FALSE, stringsAsFactors = FALSE)

qc <- merge(
  npc5_qc[, c("Query_cell", "risk_class")],
  pearson_qc[, c("Query_cell", "risk_class")],
  by = "Query_cell", suffixes = c("_NPC5", "_Pearson"), all = FALSE
)
qc$QC_level <- ifelse(
  qc$risk_class_NPC5 == "retained" & qc$risk_class_Pearson == "retained",
  "Good", "Risk"
)
qc$QC_level <- factor(qc$QC_level, levels = c("Good", "Risk"))

plot_data <- merge(
  meta[, c("Query_cell", "nCount_RNA", "nFeature_RNA")],
  qc[, c("Query_cell", "risk_class_NPC5", "risk_class_Pearson", "QC_level")],
  by = "Query_cell",
  all = FALSE,
  sort = FALSE
)
plot_data <- plot_data[order(plot_data$QC_level), , drop = FALSE]

stopifnot(
  nrow(plot_data) == 641L,
  sum(plot_data$QC_level == "Good") == 622L,
  sum(plot_data$QC_level == "Risk") == 19L,
  all(is.finite(plot_data$nCount_RNA)),
  all(is.finite(plot_data$nFeature_RNA))
)

legend_labels <- c(Good = "Good (n=622)", Risk = "Risk (n=19)")

p <- ggplot(
  plot_data,
  aes(
    x = nCount_RNA,
    y = nFeature_RNA,
    color = QC_level
  )
) +
  geom_point(size = 1, alpha = 0.72, shape = 16) +
  scale_x_log10(
    breaks = 10 ^ seq(5, 7),
    labels = label_math(10^.x),
    expand = expansion(mult = c(0.04, 0.06))
  ) +
  scale_y_continuous(
    breaks = seq(5000, 17500, by = 2500),
    labels = function(x) format(x / 1000, trim = TRUE, scientific = FALSE),
    expand = expansion(mult = c(0.03, 0.05))
  ) +
  scale_color_manual(
    values = c(Good = "#1F1F1F", Risk = "#D62728"),
    labels = legend_labels,
    drop = FALSE
  ) +
  labs(
    x = "Number of reads",
    y = "Number of genes detected (Thousand)",
    color = NULL
  ) +
  theme_classic(base_size = 8) +
  theme(
    axis.title = element_text(size = 7.2, color = "black"),
    axis.text = element_text(size = 6.2, color = "black"),
    axis.line = element_line(linewidth = 0.3528, color = "black"),
    axis.ticks = element_line(linewidth = 0.3528, color = "black"),
    legend.position = "inside",
    legend.position.inside = c(0.03, 0.98),
    legend.direction = "horizontal",
    legend.justification = c(0, 1),
    legend.text = element_text(size = 5.4),
    legend.key.width = unit(5, "pt"),
    legend.key.height = unit(5, "pt"),
    legend.spacing.x = unit(2, "pt"),
    legend.margin = margin(0, 0, 1, 0, unit = "pt"),
    legend.background = element_rect(fill = "white", color = NA),
    aspect.ratio = 1,
    plot.margin = margin(2, 7, 2, 2, unit = "pt")
  ) +
  guides(color = guide_legend(override.aes = list(size = 1.8, alpha = 1)))

png_file <- file.path(out_dir, "Figure_Query_QC_reads_genes_Good_Risk.png")
pdf_file <- file.path(out_dir, "Figure_Query_QC_reads_genes_Good_Risk.pdf")
data_file <- file.path(out_dir, "Figure_Query_QC_reads_genes_Good_Risk_plot_data.csv")

ggsave(png_file, p, width = 2, height = 2, units = "in", dpi = 600, bg = "white")
ggsave(pdf_file, p, width = 2, height = 2, units = "in", device = cairo_pdf, bg = "white")
write.csv(plot_data, data_file, row.names = FALSE)

message("PNG: ", normalizePath(png_file))
message("PDF: ", normalizePath(pdf_file))
