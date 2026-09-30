#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
})

matrix_dir <- file.path(
  "C:/Users/53461/Downloads/Patchseq_Step4_for_local",
  "Patchseq_Step4_for_local/Step4_Symbol_conversion/matrices"
)
exon_matrix_file <- file.path(matrix_dir,
  "exon_only_counts_matrix_symbol_sum.tsv")
gene_body_matrix_file <- file.path(matrix_dir,
  "gene_body_exon_intron_counts_matrix_symbol_sum.tsv")
pearson_qc_file <- "C:/Users/53461/Downloads/Top15_risk_x040_x050_QC_table.csv"
npc5_qc_file <- "C:/Users/53461/Downloads/Top15CleanMarker_MaxPC_NPC5_QC_table.csv"
out_dir <- file.path("outputs", "T-type_QC_RiskCells")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

pearson_qc <- read.csv(pearson_qc_file, check.names = FALSE,
  stringsAsFactors = FALSE)
npc5_qc <- read.csv(npc5_qc_file, check.names = FALSE,
  stringsAsFactors = FALSE)

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

selected_columns <- c("GeneSymbol", qc$Query_cell)
exon_dt <- fread(exon_matrix_file, select = selected_columns,
  showProgress = FALSE)
gene_body_dt <- fread(gene_body_matrix_file, select = selected_columns,
  showProgress = FALSE)

stopifnot(
  identical(names(exon_dt), names(gene_body_dt)),
  identical(exon_dt$GeneSymbol, gene_body_dt$GeneSymbol)
)

exon_mat <- as.matrix(exon_dt[, -1])
gene_body_mat <- as.matrix(gene_body_dt[, -1])

counts <- data.frame(
  Sample_ID = colnames(exon_mat),
  Sum_exon_counts = colSums(exon_mat),
  Sum_intron_counts = colSums(pmax(gene_body_mat - exon_mat, 0)),
  stringsAsFactors = FALSE
)

rm(exon_dt, gene_body_dt, exon_mat, gene_body_mat)
gc(verbose = FALSE)

plot_data <- merge(
  counts[, c("Sample_ID", "Sum_exon_counts", "Sum_intron_counts")],
  qc[, c("Query_cell", "risk_class_NPC5", "risk_class_Pearson", "QC_level")],
  by.x = "Sample_ID", by.y = "Query_cell", all = FALSE, sort = FALSE
)
plot_data <- plot_data[order(plot_data$QC_level), , drop = FALSE]

stopifnot(
  nrow(plot_data) == 641L,
  sum(plot_data$QC_level == "Good") == 622L,
  sum(plot_data$QC_level == "Risk") == 19L,
  all(is.finite(plot_data$Sum_exon_counts)),
  all(is.finite(plot_data$Sum_intron_counts)),
  all(plot_data$Sum_exon_counts > 0),
  all(plot_data$Sum_intron_counts > 0)
)

legend_labels <- c(Good = "Good (n=622)", Risk = "Risk (n=19)")
common_limits <- c(1e3, 1e8)

p <- ggplot(plot_data,
    aes(x = Sum_exon_counts, y = Sum_intron_counts, color = QC_level)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
    linewidth = 0.3528, color = "#6B7280") +
  geom_point(size = 1, alpha = 0.72, shape = 16) +
  scale_x_log10(
    limits = common_limits,
    breaks = 10 ^ seq(3, 8),
    labels = trans_format("log10", math_format(10^.x)),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  scale_y_log10(
    limits = common_limits,
    breaks = 10 ^ seq(3, 8),
    labels = trans_format("log10", math_format(10^.x)),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  scale_color_manual(
    values = c(Good = "#1F1F1F", Risk = "#D62728"),
    labels = legend_labels, drop = FALSE
  ) +
  labs(
    x = "Sum of exon counts",
    y = "Sum of intron counts",
    color = NULL
  ) +
  theme_classic(base_size = 8) +
  theme(
    axis.title = element_text(size = 7.2, color = "black"),
    axis.text = element_text(size = 5.6, color = "black"),
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
    plot.margin = margin(3, 7, 3, 3, unit = "pt")
  ) +
  guides(color = guide_legend(override.aes = list(size = 1.8, alpha = 1)))

png_file <- file.path(out_dir,
  "Figure_Query_QC_exon_vs_intron_counts_Good_Risk.png")
pdf_file <- file.path(out_dir,
  "Figure_Query_QC_exon_vs_intron_counts_Good_Risk.pdf")
data_file <- file.path(out_dir,
  "Figure_Query_QC_exon_vs_intron_counts_Good_Risk_plot_data.csv")

ggsave(png_file, p, width = 2, height = 2, units = "in", dpi = 600,
  bg = "white")
ggsave(pdf_file, p, width = 2, height = 2, units = "in",
  device = cairo_pdf, bg = "white")
write.csv(plot_data, data_file, row.names = FALSE)

message("Matched rows: ", nrow(plot_data))
message("Exon range: ", paste(range(plot_data$Sum_exon_counts), collapse = " - "))
message("Intron range: ", paste(range(plot_data$Sum_intron_counts), collapse = " - "))
message("PNG: ", normalizePath(png_file))
message("PDF: ", normalizePath(pdf_file))
